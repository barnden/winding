#include "Simulator.h"

#include <iostream>
#include <ranges>

#include <Eigen/Sparse>
using MatrixSd = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;

Simulator::Simulator(
    ParametricSurface const& f,
    std::vector<Vec2> const& init_path)
    : m_size(201)
    , m_surface(f)
    , m_ksp(1'000'000)
    , m_kdp(200)
    , m_kpf(10'000)
    , m_kbs(0.001)
    , m_dt(0.001)
    , m_t(0.)
    , m_v_eps(0.001)
{
    if (init_path.size()) {
        m_size = init_path.size();
        m_p_init = init_path;

        return;
    }

    m_p_init = decltype(m_p_init)(m_size, Vec2::Zero());

    Vec2 p0(0., -1.);
    Vec2 p1(9.43347686131086, 0.9965750445589878);
    Vec2 dp = (p1 - p0) / (m_size - 1.);

    for (int i = 0; i < m_size; i++)
        m_p_init[i] = p0 + (double)i * dp;
}

Simulator::~Simulator() {};

void Simulator::simulate(int k)
{
    int i = 0;
    for (; (i < k && i < 10) || (i < k && !stop()); i++) {
        m_t += m_dt;
        step();
    }

    std::cout << i << " iterations simulated.\n";
}

OffSurface::OffSurface(ParametricSurface const& f, std::vector<Vec2> const& init_path)
    : Simulator(f, init_path)
{
    m_p = m_p_init;
    m_v = decltype(m_v)(size(), Vec2::Zero());
    m_touching = decltype(m_touching)(size(), true);
    m_l = decltype(m_l)(size(), 1);
    m_r = decltype(m_r)(size(), 1);
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

    std::vector<Triplet> JT_ele(6 * N, Triplet());
    std::vector<Triplet> dJ_ele(6 * N, Triplet());
    std::vector<Triplet> Mh2K_ele(9 * N - 12);

    for (auto&& [i, pidx] : enumerate(touching)) {
        Vec2 const& p = m_p[pidx];
        Vec2 const& v = m_v[pidx];

        varr(2 * i) = v.x();
        varr(2 * i + 1) = v.y();

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

            forces(0) = 0;
            forces(1) = 0;
            forces(2) = 0;

            continue;
        }

        if (i == N - 1) {
            Mh2K_ele[9 * N - 15] = Triplet(3 * i + 0, 3 * i + 0, 1);
            Mh2K_ele[9 * N - 14] = Triplet(3 * i + 1, 3 * i + 1, 1);
            Mh2K_ele[9 * N - 13] = Triplet(3 * i + 2, 3 * i + 2, 1);

            forces(3 * i + 0) = 0;
            forces(3 * i + 1) = 0;
            forces(3 * i + 2) = 0;
            continue;
        }

        idx = 9 * i;

        int const& L = m_l[pidx];
        int const& R = m_r[pidx];

        double tmp = 1. + (1. / L + 1. / R) * m_ksp * m_dt * m_dt + m_kdp * m_dt;
        Mh2K_ele[idx - 6] = Triplet(3 * i + 0, 3 * i + 0, tmp);
        Mh2K_ele[idx - 3] = Triplet(3 * i + 1, 3 * i + 1, tmp);
        Mh2K_ele[idx + 0] = Triplet(3 * i + 2, 3 * i + 2, tmp);

        tmp = -m_ksp * m_dt * m_dt / L;
        Mh2K_ele[idx - 5] = Triplet(3 * i + 0, 3 * i - 3, tmp);
        Mh2K_ele[idx - 2] = Triplet(3 * i + 1, 3 * i - 2, tmp);
        Mh2K_ele[idx + 1] = Triplet(3 * i + 2, 3 * i - 1, tmp);

        tmp = -m_ksp * m_dt * m_dt / R;
        Mh2K_ele[idx - 4] = Triplet(3 * i + 0, 3 * i + 3, tmp);
        Mh2K_ele[idx - 1] = Triplet(3 * i + 1, 3 * i + 4, tmp);
        Mh2K_ele[idx + 2] = Triplet(3 * i + 2, 3 * i + 5, tmp);

        // clang-format off
        Vec3 force = m_ksp
                   * (1. / L * (surface().f(m_p[pidx - L]) - surface().f(p))
                   +  1. / R * (surface().f(m_p[pidx + R]) - surface().f(p)));
        // clang-format on

        forces(3 * i + 0) = force.x();
        forces(3 * i + 1) = force.y();
        forces(3 * i + 2) = force.z();
    }

    JT.setFromTriplets(JT_ele.begin(), JT_ele.end());
    dJ.setFromTriplets(dJ_ele.begin(), dJ_ele.end());
    Mh2K.setFromTriplets(Mh2K_ele.begin(), Mh2K_ele.end());

    mat = JT * Mh2K * JT.transpose();
    barr = JT * (JT.transpose() * varr) + m_dt * (JT * (forces - dJ * varr));
    Eigen::SparseLU<MatrixSd, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(mat);
    solver.factorize(mat);
    varr = solver.solve(barr);

    for (int i = 1; i < N - 1; i++) {
        int const& idx = touching[i];

        m_v[idx].x() = varr(2 * i);
        m_v[idx].y() = varr(2 * i + 1);

        m_p[idx].x() += m_v[idx].x() * m_dt;
        m_p[idx].y() += m_v[idx].y() * m_dt;
    }

    lifting();
    mapping();
    landing();
}

bool OffSurface::stop()
{
    for (int i = 0; i < size(); i++)
        if (m_touching[i] && (m_v[i].x() * surface().f_u(m_p[i]) + m_v[i].y() * surface().f_v(m_p[i])).norm() > m_v_eps)
            return false;

    return true;
}

void OffSurface::mapping()
{
    int cl = 0;
    int cr = m_r[0];

    while (cl != size() - 1) {
        int len = m_r[cl];
        cr = cl + len;
        Vec3 p0 = surface().f(m_p[cl]);
        Vec3 p1 = surface().f(m_p[cr]);

        for (int i = 1; i < m_r[cl]; i++)
            m_p[cl + i] = surface().closest_point(((double)i / len) * p1 + ((double)(len - i) / len) * p0, m_p[cl + i]);

        cl = cr;
    }
}

void OffSurface::lifting()
{
    int cl = 0;
    int cm = m_r[0];

    Vec3 pl = surface().f(m_p[0]);
    Vec3 pm = surface().f(m_p[cm]);

    while (cm < size() - 1) {
        int cr = cm + m_r[cm];
        Vec3 pr = surface().f(m_p[cr]);

        if ((pl + pr - 2. * pm).dot(surface().normal(m_p[cm])) > 0.) { // lift
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

    Vec3 pl = surface().f(m_p[0]);
    Vec3 pr = surface().f(m_p[cr]);

    while (cl != size() - 1) {
        int cm = cl + 1;
        Vec3 pm;
        double mu;

        for (; cm < cr; cm++) {
            pm = surface().f(m_p[cm]);
            mu = (double)(cm - cl) / (cr - cl);

            if ((pm - (mu * pr + (1. - mu) * pl)).dot(surface().normal(m_p[cm])) > 0.) {
                m_l[cm] = cm - cl;
                m_r[cl] = cm - cl;

                m_l[cr] = cr - cm;
                m_r[cm] = cr - cm;

                m_touching[cm] = true;
                m_v[cm] = mu * m_v[cr] + (1. - mu) * m_v[cl];

                break;
            }
        }

        if (cm == cr) {
            if (cr == size() - 1)
                return;

            cl = cr;
            cr = cl + m_r[cl];
            pl = pr;
            pr = surface().f(m_p[cr]);

            continue;
        }

        pl = pm;
        cl = cm;
    }
}
