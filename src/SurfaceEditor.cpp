/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#if USE_QPMAD
#    include <qpmad/solver.h>
#else
#    include "epigraph.hpp"
#endif

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>

#include <chrono>

#include "Composer.h"
#include "Config.h"
#include "SurfaceEditor.h"

namespace fs = std::filesystem;

void write_winding_path(size_t step, ParametricSurface const& surface, std::vector<Composer::LocalFrame> const& path)
{
    auto outpath = std::format("{}/{}/path/step-{}.txt", Config::out_directory, Config::experiment, step);
    auto ostream = std::ofstream(outpath);

    for (auto&& [i, frame] : enumerate(path)) {
        Vec3 shadow = surface.f(frame.shadow);
        ostream << std::format("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n", i, frame.position.x(), frame.position.y(), frame.position.z(), frame.distance, frame.path_distance, frame.angle, shadow.x(), shadow.y(), shadow.z(), frame.direction.x(), frame.direction.y(), frame.direction.z(), frame.normal.x(), frame.normal.y(), frame.normal.z());
    }

    ostream.flush();
}

void write_maxquad(size_t step, Composer::Quad const& quad)
{
    auto outpath = std::format("{}/{}/max_quad/step-{}.txt", Config::out_directory, Config::experiment, step);
    auto ostream = std::ofstream(outpath);

    ostream << quad.p0.x() << ", " << quad.p0.y() << ", " << quad.p0.z() << '\n'
            << quad.p1.x() << ", " << quad.p1.y() << ", " << quad.p1.z() << '\n'
            << quad.p2.x() << ", " << quad.p2.y() << ", " << quad.p2.z() << '\n'
            << quad.p3.x() << ", " << quad.p3.y() << ", " << quad.p3.z() << '\n';
}

void SurfaceEditor::step(
    double timestep,
    double spring_cosntant,
    double damping_coefficient,
    double epsilon)
{
    std::cout << "[Editor] Begin step " << (m_step + 1) << '\n';

    fs::copy_file(std::format("{}/{}.txt", Config::data_directory, Config::stem), std::format("{}/{}/spline/step-0.txt", Config::out_directory, Config::experiment), fs::copy_options::overwrite_existing);

    auto composer = Composer(m_surface, m_num_revolutions, m_num_paths, m_num_particles);
    auto path = composer.simulate(timestep, spring_cosntant, damping_coefficient, epsilon, m_step);
    auto max_quad = composer.max_quad;
    auto in_surface = std::vector<std::pair<Vec2, Vec3>> {};
    auto off_surface = std::vector<std::pair<Vec2, Vec3>> {};

    write_winding_path(m_step, m_surface, path);

    bool do_push_out = (Config::alternate_push_pull && m_step % 2 == 0) || (!Config::alternate_push_pull && Config::push_out);

    auto start_build = std::chrono::steady_clock::now();
    if (do_push_out) {
        for (auto&& [i, frame] : enumerate(path)) {
            if (frame.distance < 1e-3)
                continue;

            off_surface.push_back(std::make_pair(frame.shadow, frame.position));
        }
    }

    // write_maxquad(m_step, composer.max_quad);
    bool do_pull_in = (Config::alternate_push_pull && m_step % 2 == 1) || (!Config::alternate_push_pull && Config::pull_in);

    // std::cout << "status: " << (do_push_out ? "push out\t" : "") << (do_pull_in ? "pull in\t" : "") << '\n';

    if (do_pull_in) {
        std::sort(composer.quads.begin(), composer.quads.end(), [](Composer::Quad a, Composer::Quad b) { return a.area < b.area; });

        auto N = composer.quads.size();
        // auto median = N % 2 ? composer.quads[N / 2].area : (composer.quads[N / 2].area + composer.quads[(N / 2) - 1].area) / 2.;
        auto q1 = N % 2 ? composer.quads[N / 4].area : (composer.quads[N / 4].area + composer.quads[(N / 4) - 1].area) / 2.;
        auto q3 = N % 2 ? composer.quads[(3 * N) / 4].area : (composer.quads[(3 * N) / 4].area + composer.quads[((3 * N) / 4) - 1].area) / 2.;
        auto iqr = 2. * (q3 - q1) + q3;
        // std::cout << std::format("iqr: {}, {}, {}\n", iqr, composer.quads[N-1].area, composer.quads[N-10].area);

        auto has_outliers = std::any_of(composer.quads.begin() + (3 * N) / 4, composer.quads.end(), [&iqr](Composer::Quad a) { return a.area >= iqr; });

        if (has_outliers) {
            auto constexpr num_particles = 35;

            // auto it = std::lower_bound(composer.quads.begin(),
            //                            composer.quads.end(),
            //                            composer.quads.back().area - 1e-2,
            //                            [](Composer::Quad const& a, double& v) { return a.area < v; });

            auto it = std::next(composer.quads.begin(), (3 * N) / 4);

            // Use logistic function mapping m_step -> [epsilon, 1.]
            auto constexpr epsilon = 0.90;
            auto constexpr k = 0.125;
            auto threshold = epsilon + 2. * (1. - epsilon) * (1. / (1. + std::exp(-m_step * k / 2.)) - 0.5);
            auto lower_bound = composer.quads.back().area * threshold;
            for (; it != composer.quads.end(); it = std::next(it)) {
                if (it->area >= lower_bound)
                    break;
            }

            for (; it != composer.quads.end(); it = std::next(it)) {
                auto const& quad = *it;

                Vec3 m0 = (quad.p1 + quad.p0) / 2.;
                Vec3 m1 = (quad.p3 + quad.p2) / 2.;
                Vec3 m2 = (quad.p0 + quad.p3) / 2.;
                Vec3 m3 = (quad.p1 + quad.p2) / 2.;

                Vec3 rd0 = Vec3::Zero();
                Vec3 rd1 = Vec3::Zero();

                Vec3 d0 = (m1 - m0).normalized();
                Vec3 d1 = (m3 - m2).normalized();

                for (auto i = 0; i < num_particles; i++) {
                    double t = ((double)i) / (num_particles - 1.);
                    Vec3 w0 = m0 + t * d0;
                    Vec3 w1 = m2 + t * d1;

                    Vec2 uv0 = Vec2::Zero();
                    Vec2 uv1 = Vec2::Zero();
                    uv0 = m_surface.closest_point(w0);
                    uv1 = m_surface.closest_point(w1);

                    in_surface.push_back(std::make_pair(uv0, w0));
                    in_surface.push_back(std::make_pair(uv1, w1));
                }
            }
        }
    }

    // std::cout << "[debug] finish " << (do_push_out ? "push-out" : "pull-in") << " path\n";
    auto m = m_surface.m_nv * m_surface.m_nu;
    auto k = off_surface.size() + in_surface.size();

    if (k == 0) {
        // std::cout << "[Editor] No path particles are off-surface.\n";
        return;
    } else {
        // std::cout << std::format("[Editor] Iteration {} found {} optimizable particles.\n", m_step + 1, k);
    }

    Eigen::MatrixXd M(k, 3 * m);
    Eigen::VectorXd b(k);

    for (auto&& [i, pair] : enumerate(off_surface)) {
        auto const& [shadow, world] = pair;

        Vec3 normal = m_surface.normal(shadow);
        M.row(i) = normal.transpose() * m_surface.jacobian(shadow);
        b(i) = (world - m_surface.f(shadow)).dot(normal);
    }

    for (auto&& [i, pair] : enumerate(in_surface)) {
        auto const& [shadow, world] = pair;

        Vec3 normal = -m_surface.normal(shadow);
        M.row(i + off_surface.size()) = normal.transpose() * m_surface.jacobian(shadow);
        b(i + off_surface.size()) = (world - m_surface.f(shadow)).dot(normal);
    }
    auto end_build = std::chrono::steady_clock::now();

    auto duration_build = std::chrono::duration_cast<std::chrono::microseconds>(end_build - start_build);

    std::cout << "[debug] Build time: " << duration_build.count() << " microseconds" << std::endl;

    auto start = std::chrono::steady_clock::now();

#if USE_QPMAD
    auto solver = qpmad::Solver();
    Eigen::VectorXd deltaC(3 * m);
    Eigen::MatrixXd static H = Eigen::MatrixXd::Identity(3 * m, 3 * m);
    Eigen::VectorXd static h = Eigen::VectorXd::Zero(3 * m);
    Eigen::VectorXd Aub = (1e15 * Eigen::VectorXd::Ones(k));
    auto status = solver.solve(deltaC, H, h, M, b, Aub);

    if (status != qpmad::Solver::OK) {
        std::cerr << "[Editor] Failed QP solve\n";
        exit(EXIT_FAILURE);
    }
#else
    auto qp = cvx::OptimizationProblem {};
    auto dC = qp.addVariable("deltaC", 3 * m);

    cvx::MatrixX MdC = cvx::par(M) * dC;
    auto constraint = cvx::greaterThan(MdC, cvx::par(b));
    qp.addConstraint(constraint);
    qp.addCostTerm(dC.transpose() * dC);

    auto solver = cvx::osqp::OSQPSolver(qp);
    solver.solve();

    if (solver.getExitCode() != 0) {
        std::cerr << "[Editor] Failed QP solve\n";
        exit(EXIT_FAILURE);
    }

    Eigen::VectorXd deltaC = cvx::eval(dC);

#endif

    auto end = std::chrono::steady_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::println("[Timing] Solved system in {} us", duration.count());
    std::println("[Editor] Solved least-norm QP (deltaC L2: {})", deltaC.norm());

    auto augmented_control_points = std::vector<std::vector<Vec3>>(m_surface.m_nv, std::vector<Vec3>(m_surface.m_nu, Vec3::Zero()));
    {
        auto index = 0;
        for (auto i = 0; i < m_surface.m_nv; i++) {
            for (auto j = 0; j < m_surface.m_nu; j++) {
                augmented_control_points[i][j] = m_surface.points()[i][j] + deltaC.segment<3>(3 * index);
                index++;
            }
        }
    }

    {
        auto outpath = std::format("{}/{}/spline/step-{}.txt", Config::out_directory, Config::experiment, m_step + 1);
        auto ostream = std::ofstream(outpath);

        ostream << m_surface.m_nv << ' ' << m_surface.m_nu << '\n';

        for (auto i = 0; i < m_surface.m_nv; i++) {
            for (auto j = 0; j < m_surface.m_nu; j++) {
                ostream << std::format("{} {} {}\n",
                                       augmented_control_points[i][j].x(),
                                       augmented_control_points[i][j].y(),
                                       augmented_control_points[i][j].z());
            }
        }
        ostream.close();
    }

    auto modified_surface = CubicBSpline(m_surface.m_nv, m_surface.m_nu, std::move(augmented_control_points));

    std::cout << "[Editor] Hausdorff: " << hausdorff_distance(m_surface, modified_surface) << '\n';
    std::cout << "[Editor] pct. diff: " << percent_difference(m_surface, modified_surface) << '\n';

    std::swap(m_surface, modified_surface);

    m_step++;
}
