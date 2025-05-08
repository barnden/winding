/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Simulator.h"
#include "Config.h"

#include <Eigen/Sparse>
#include <print>

using MatrixSd = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;

Simulator::Simulator(
    ParametricSurface const& f,
    std::vector<Vec2> const& init_path)
    : m_size(201)
    , m_surface(f)
    , m_spring_constant(1'000'000)
    , m_damping_coefficient(200)
    , m_timestep(0.001)
    , m_t(0.)
    , m_epsilon(0.001)
{
    if (!init_path.empty()) {
        m_size = init_path.size();
        m_position_initial = init_path;

        return;
    }

    m_position_initial = std::vector<Vec2>(m_size, Vec2::Zero());

    Vec2 p0(0., -1.);
    Vec2 p1(9.43347686131086, 0.9965750445589878);
    Vec2 dp = (p1 - p0) / (m_size - 1.);

    for (int i = 0; i < m_size; i++)
        m_position_initial[i] = p0 + (double)i * dp;
}

void Simulator::simulate(int num_iterations)
{
    int i;
    for (i = 0; (i < num_iterations && i < 10) || (i < num_iterations && !stop()); i++) {
        m_t += m_timestep;
        step();
    }

    std::println("[Simulator] {} iterations simulated.", i);
}

OffSurface::OffSurface(ParametricSurface const& f, std::vector<Vec2> const& init_path)
    : Simulator(f, init_path)
{
    m_position = m_position_initial;
    m_velocity = std::vector<Vec2>(size(), Vec2::Zero());
    m_touching = std::vector(size(), true);
    m_l = std::vector(size(), 1);
    m_r = std::vector(size(), 1);
}

void OffSurface::step()
{
    std::vector<int> touching;

    for (int i = 0; i < size(); i++) {
        if (m_touching[i])
            touching.push_back(i);
    }

    const int N = touching.size();

    MatrixSd JT(2 * N, 3 * N);
    MatrixSd dJ(3 * N, 2 * N);
    MatrixSd Mh2K(3 * N, 3 * N);
    MatrixSd mat;

    Vector forces(3 * N);
    Vector varr(2 * N);
    Vector barr;

    auto JT_ele = std::vector(6 * N, Triplet());
    auto dJ_ele = std::vector(6 * N, Triplet());
    auto Mh2K_ele = std::vector(9 * N - 12, Triplet());

    for (auto&& [i, pidx] : enumerate(touching)) {
        Vec2 const& p = m_position[pidx];
        Vec2 const& v = m_velocity[pidx];

        varr.segment<2>(2 * i) = v;

        Vec3 fu = surface().f_u(p);
        Vec3 fv = surface().f_v(p);
        Vec3 fuu = surface().f_uu(p);
        Vec3 fuv = surface().f_uv(p);
        Vec3 fvv = surface().f_vv(p);

        Vec3 dj_u = v.x() * fuu + v.y() * fuv;
        Vec3 dj_v = v.x() * fuv + v.y() * fvv;

        auto idx = 6 * i;
        JT_ele[idx + 0] = Triplet(2 * i + 0, 3 * i + 0, fu.x());
        JT_ele[idx + 1] = Triplet(2 * i + 0, 3 * i + 1, fu.y());
        JT_ele[idx + 2] = Triplet(2 * i + 0, 3 * i + 2, fu.z());
        JT_ele[idx + 3] = Triplet(2 * i + 1, 3 * i + 0, fv.x());
        JT_ele[idx + 4] = Triplet(2 * i + 1, 3 * i + 1, fv.y());
        JT_ele[idx + 5] = Triplet(2 * i + 1, 3 * i + 2, fv.z());

        dJ_ele[idx + 0] = Triplet(3 * i + 0, 2 * i + 0, dj_u.x());
        dJ_ele[idx + 2] = Triplet(3 * i + 1, 2 * i + 0, dj_u.y());
        dJ_ele[idx + 4] = Triplet(3 * i + 2, 2 * i + 0, dj_u.z());
        dJ_ele[idx + 1] = Triplet(3 * i + 0, 2 * i + 1, dj_v.x());
        dJ_ele[idx + 3] = Triplet(3 * i + 1, 2 * i + 1, dj_v.y());
        dJ_ele[idx + 5] = Triplet(3 * i + 2, 2 * i + 1, dj_v.z());

        if (i == 0) {
            Mh2K_ele[0] = Triplet(0, 0, 1);
            Mh2K_ele[1] = Triplet(1, 1, 1);
            Mh2K_ele[2] = Triplet(2, 2, 1);

            forces.segment<3>(0) = Vec3::Zero();
            continue;
        }

        if (i == N - 1) {
            Mh2K_ele[9 * N - 15] = Triplet(3 * i + 0, 3 * i + 0, 1);
            Mh2K_ele[9 * N - 14] = Triplet(3 * i + 1, 3 * i + 1, 1);
            Mh2K_ele[9 * N - 13] = Triplet(3 * i + 2, 3 * i + 2, 1);

            forces.segment<3>(3 * i) = Vec3::Zero();
            continue;
        }

        idx = 9 * i;

        int const& L = m_l[pidx];
        int const& R = m_r[pidx];

        double tmp = 1. + ((1. / L + 1. / R) * m_spring_constant * m_timestep + m_damping_coefficient) * m_timestep;
        Mh2K_ele[idx - 6] = Triplet(3 * i + 0, 3 * i + 0, tmp);
        Mh2K_ele[idx - 3] = Triplet(3 * i + 1, 3 * i + 1, tmp);
        Mh2K_ele[idx + 0] = Triplet(3 * i + 2, 3 * i + 2, tmp);

        tmp = -m_spring_constant * m_timestep * m_timestep / L;
        Mh2K_ele[idx - 5] = Triplet(3 * i + 0, 3 * i - 3, tmp);
        Mh2K_ele[idx - 2] = Triplet(3 * i + 1, 3 * i - 2, tmp);
        Mh2K_ele[idx + 1] = Triplet(3 * i + 2, 3 * i - 1, tmp);

        tmp = -m_spring_constant * m_timestep * m_timestep / R;
        Mh2K_ele[idx - 4] = Triplet(3 * i + 0, 3 * i + 3, tmp);
        Mh2K_ele[idx - 1] = Triplet(3 * i + 1, 3 * i + 4, tmp);
        Mh2K_ele[idx + 2] = Triplet(3 * i + 2, 3 * i + 5, tmp);

        // clang-format off
        Vec3 force = m_spring_constant
                   * (1. / L * (surface().f(m_position[pidx - L]) - surface().f(p))
                   +  1. / R * (surface().f(m_position[pidx + R]) - surface().f(p)));
        // clang-format on

        forces.segment<3>(3 * i) = force;
    }

    JT.setFromTriplets(JT_ele.begin(), JT_ele.end());
    dJ.setFromTriplets(dJ_ele.begin(), dJ_ele.end());
    Mh2K.setFromTriplets(Mh2K_ele.begin(), Mh2K_ele.end());

    mat = JT * Mh2K * JT.transpose();

    // World space acceleration and velocity
    Eigen::VectorXd A = forces - dJ * varr;
    Eigen::VectorXd V = JT.transpose() * varr;

    barr = (JT * V) + (m_timestep * JT * A);
    Eigen::SparseLU<MatrixSd, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(mat);
    solver.factorize(mat);

    if (solver.info() != Eigen::Success) {
        std::println(
            "Failed to factorize matrix (info: {}) (matrix dimensions: {} x {})",
            static_cast<int>(solver.info()), mat.rows(), mat.cols());

        exit(EXIT_FAILURE);
    }

    varr = solver.solve(barr);

    for (int i = 1; i < N - 1; i++) {
        int const& idx = touching[i];

        Vec2 velocity = varr.segment<2>(2 * i);

        if (Config::friction_coefficient > 0.) {
            // Equation 4.19 in Li's disseration
            Vec3 acceleration = A.segment<3>(3 * i);
            Vec3 velocity = V.segment<3>(3 * i);
            auto a_dot_n = std::abs(acceleration.dot(surface().normal(m_position[idx])));
            auto friction = (m_timestep * Config::friction_coefficient * a_dot_n) / velocity.norm();

            velocity *= std::max(1. - friction, 0.);
        }

        m_velocity[idx] = velocity;
        m_position[idx] += m_velocity[idx] * m_timestep;
    }

    lifting();
    mapping();
    landing();
}

bool OffSurface::stop()
{
    auto stopping_criteria =
        [&](auto i) {
            Eigen::Matrix<double, 3, 2> mat;
            mat << surface().f_u(m_position[i]), surface().f_v(m_position[i]);

            return (mat * m_velocity[i]).norm() > m_epsilon;
        };

    for (int i = 0; i < size(); i++)
        if (m_touching[i] && stopping_criteria(i))
            return false;

    return true;
}

void OffSurface::mapping()
{
    int cl = 0;
    int cr = 0;

    while (cl != size() - 1) {
        int len = m_r[cl];
        cr = cl + len;
        Vec3 p0 = surface().f(m_position[cl]);
        Vec3 p1 = surface().f(m_position[cr]);

        if (Config::use_ray_shoot_mapping) {
            Vec3 n0 = surface().normal(m_position[cl]);
            Vec3 n1 = surface().normal(m_position[cr]);

            Vec3 d = (p1 - p0).normalized();
            Vec3 eta = (n0 + n1).normalized();
            Vec3 ray_direction = (eta.cross(d)).cross(d).normalized();

            for (auto i = 1; i < len; i++) {
                auto s = ((double)i) / len;
                Vec3 world = s * p1 + (1. - s) * p0;

                Vec3 p = world;
                auto t = 0.;
                for (auto j = 0; j < Config::sphere_tracing_iterations; j++) {
                    auto h = surface().sdf(p);

                    if (std::abs(h) < 1e-3)
                        break;

                    t += h;
                    p = t * ray_direction + world;
                }

                Vec2 parametric = surface().closest_point(p);
                m_position[cl + i] = parametric;
            }
        } else {
            for (int i = 1; i < m_r[cl]; i++) {
                auto t = ((double)i) / len;
                Vec3 world = t * p1 + (1. - t) * p0;

                if (Config::use_bvh) {
                    m_position[cl + i] = surface().closest_point(world);
                } else {
                    Vec2 guess = m_position[cl + i];
                    m_position[cl + i] = surface().closest_point(world, guess);
                }
            }
        }

        cl = cr;
    }
}

void OffSurface::lifting()
{
    int cl = 0;
    int cm = m_r[0];

    Vec3 pl = surface().f(m_position[0]);
    Vec3 pm = surface().f(m_position[cm]);

    while (cm < size() - 1) {
        int cr = cm + m_r[cm];
        Vec3 pr = surface().f(m_position[cr]);

        if ((pl + pr - 2. * pm).dot(surface().normal(m_position[cm])) > 0.) { // lift
            m_touching[cm] = false;
            m_l[cm] = 1;
            m_r[cm] = 1;
            m_l[cr] = cr - cl;
            m_r[cl] = cr - cl;
            cm = cr;
            pm = pr;

            continue;
        }

        cl = cm;
        cm = cr;
        pl = pm;
        pm = pr;
    }
}

void OffSurface::landing()
{
    int cl = 0;
    int cr = m_r[0];

    Vec3 pl = surface().f(m_position[0]);
    Vec3 pr = surface().f(m_position[cr]);

    while (cl != size() - 1) {
        int cm = cl + 1;
        Vec3 pm;

        for (; cm < cr; cm++) {
            double mu = (double)(cm - cl) / (cr - cl);
            pm = surface().f(m_position[cm]);

            if ((pm - (mu * pr + (1. - mu) * pl)).dot(surface().normal(m_position[cm])) > 0.) {
                m_l[cm] = cm - cl;
                m_r[cl] = cm - cl;

                m_l[cr] = cr - cm;
                m_r[cm] = cr - cm;

                m_touching[cm] = true;
                m_velocity[cm] = mu * m_velocity[cr] + (1. - mu) * m_velocity[cl];

                break;
            }
        }

        if (cm == cr) {
            if (cr == size() - 1)
                return;

            cl = cr;
            cr = cl + m_r[cl];
            pl = pr;
            pr = surface().f(m_position[cr]);

            continue;
        }

        pl = pm;
        cl = cm;
    }
}
