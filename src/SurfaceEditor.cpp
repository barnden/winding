/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#if USE_QPMAD
#    include <qpmad/solver.h>
#else
#    include "epigraph.hpp"
#endif

#include <fstream>
#include <iostream>

#include "SurfaceEditor.h"

void SurfaceEditor::step(double dt, double ksp, double kdp, double eps)
{
    std::cout << "[Editor] Begin step " << (m_step + 1) << '\n';

    auto composer = Composer(m_options, m_surface, m_num_revolutions, m_num_paths, m_num_particles);
    auto path = composer.simulate(dt, ksp, kdp, eps);
    auto max_quad = composer.max_quad;
    auto in_surface = std::vector<std::pair<Vec2, Vec3>> {};
    auto off_surface = std::vector<std::pair<Vec2, Vec3>> {};
    {
        auto ostream = std::ofstream(m_options.out_path + "/" + m_options.experiment + "/path/step-" + std::to_string(m_step) + ".txt");
        for (auto&& [i, frame] : enumerate(path)) {
            bool do_push = m_options.push_out && frame.distance > 0. && (!m_options.pull_in || m_step % 2 == 1);

            if (do_push) {
                Vec2 uv = frame.shadow;
                // Vec2 uv = m_surface.closest_point(frame.position, uv, 250);
                off_surface.push_back(std::make_pair(uv, frame.position));
            }

            ostream << i << ", " << frame.position.x() << ", " << frame.position.y() << ", " << frame.position.z() << ", " << frame.distance << ", " << frame.path_distance << ", " << frame.angle << '\n';
        }
    }

    {
        auto const N = 50;
        auto ostream = std::ofstream(m_options.out_path + "/" + m_options.experiment + "/max_quad/step-" + std::to_string(m_step) + ".txt");

        Vec3 p0 = max_quad.segment<3>(0);
        Vec3 p1 = max_quad.segment<3>(3);
        Vec3 p2 = max_quad.segment<3>(6);
        Vec3 p3 = max_quad.segment<3>(9);

        ostream << p0.x() << ", " << p0.y() << ", " << p0.z() << '\n'
                << p1.x() << ", " << p1.y() << ", " << p1.z() << '\n'
                << p2.x() << ", " << p2.y() << ", " << p2.z() << '\n'
                << p3.x() << ", " << p3.y() << ", " << p3.z() << '\n';

        Vec3 q0 = (p1 + p0) / 2.;
        Vec3 q1 = (p3 + p2) / 2.;
        Vec3 q2 = (p0 + p3) / 2.;
        Vec3 q3 = (p1 + p2) / 2.;

        bool do_pull = m_options.pull_in && (!m_options.push_out || m_step % 2 == 0);

        if (do_pull) {
            for (auto i = 0; i < N; i++) {
                double t = ((double)i) / (N - 1.);

                Vec3 w0 = q0 + t * (q1 - q0);
                Vec3 w1 = q2 + t * (q3 - q2);
                Vec2 uv0 = m_surface.closest_point(w0);
                Vec2 uv1 = m_surface.closest_point(w1);

                if (m_step % 2 == 0) {
                    in_surface.push_back(std::make_pair(uv0, w0));
                    in_surface.push_back(std::make_pair(uv1, w1));
                }
            }
        }
    }
    std::cout << "[debug] finish pull-in path\n";
    auto m = m_surface.m_nv * m_surface.m_nu;
    auto k = off_surface.size() + in_surface.size();

    if (k == 0) {
        std::cout << "[Editor] No path particles are off-surface." << '\n';
        return;
    } else {
        std::cout << "[Editor] Iteration " << (m_step + 1) << " found " << k << " off-surface particles.\n";
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
    std::cout << "[Editor] Solved least-norm QP, solution L2 norm " << deltaC.norm() << '\n';
    auto augmented_control_points = std::vector<std::vector<Vec3>>(m_surface.m_nv, std::vector<Vec3>(m_surface.m_nu, Vec3::Zero()));
    auto ostream = std::ofstream(m_options.out_path + "/" + m_options.experiment + "/spline/step-" + std::to_string(m_step + 1) + ".txt");

    ostream << m_surface.m_nv << ' ' << m_surface.m_nu << '\n';

    auto index = 0;
    for (auto i = 0; i < m_surface.m_nv; i++) {
        for (auto j = 0; j < m_surface.m_nu; j++) {
            augmented_control_points[i][j] = m_surface.points()[i][j] + deltaC.segment<3>(3 * index);
            // augmented_control_points[i][j + 2 * m_surface.m_nu] = (augmented_control_points[i][j + m_surface.m_nu] = augmented_control_points[i][j]);
            index++;

            ostream << augmented_control_points[i][j].x() << ' '
                    << augmented_control_points[i][j].y() << ' '
                    << augmented_control_points[i][j].z() << '\n';
        }
    }
    ostream.close();

    m_surface = CubicBSpline(m_surface.m_nv, m_surface.m_nu, std::move(augmented_control_points));
    m_step++;
}

/**
 * void SurfaceEditor::step(double dt, double ksp, double kdp, double eps)
{
    std::cout << "[Editor] Begin step " << (m_step + 1) << '\n';

    auto composer = Composer(m_surface, m_num_revolutions, m_num_paths, m_num_particles);
    auto path = composer.simulate(dt, ksp, kdp, eps);
    auto max_quad = composer.max_quad;
    auto off_surface = std::vector<std::pair<Vec2, Vec3>> {};
    {
        auto ostream = std::ofstream(m_file_stem + "-path-" + std::to_string(m_step) + ".txt");
        for (auto&& [i, frame] : enumerate(path)) {
            if (frame.distance > 0.) {
                Vec2 uv = frame.shadow;
                // Vec2 uv = m_surface.closest_point(frame.position, uv, 250);
                off_surface.push_back(std::make_pair(uv, frame.position));
            }

            ostream << i << ", " << frame.position.x() << ", " << frame.position.y() << ", " << frame.position.z() << ", " << frame.distance << ", " << frame.path_distance << ", " << frame.angle << '\n';
        }
        ostream.close();
    }

    auto m = m_surface.m_nv * m_surface.m_nu;
    auto k = off_surface.size();

    if (k == 0) {
        std::cout << "[Editor] No path particles are off-surface." << '\n';
        return;
    } else {
        std::cout << "[Editor] Iteration " << (m_step + 1) << " found " << k << " off-surface particles.\n";
    }

    Eigen::MatrixXd M(k, 3 * m);
    Eigen::VectorXd b(k);

    for (auto&& [i, pair] : enumerate(off_surface)) {
        auto const& [shadow, world] = pair;

        Vec3 normal = m_surface.normal(shadow);
        M.row(i) = normal.transpose() * m_surface.jacobian(shadow);
        b(i) = (world - m_surface.f(shadow)).dot(normal);
    }


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
    std::cout << "[Editor] Solved least-norm QP, solution L2 norm " << deltaC.norm() << '\n';
    auto augmented_control_points = std::vector<std::vector<Vec3>>(m_surface.m_nv, std::vector<Vec3>(m_surface.m_nu, Vec3::Zero()));
    auto ostream = std::ofstream(m_file_stem + "-spline-step-" + std::to_string(m_step + 1) + ".txt");

    ostream << m_surface.m_nv << ' ' << m_surface.m_nu << '\n';

    auto index = 0;
    for (auto i = 0; i < m_surface.m_nv; i++) {
        for (auto j = 0; j < m_surface.m_nu; j++) {
            augmented_control_points[i][j] = m_surface.points()[i][j] + deltaC.segment<3>(3 * index);
            // augmented_control_points[i][j + 2 * m_surface.m_nu] = (augmented_control_points[i][j + m_surface.m_nu] = augmented_control_points[i][j]);
            index++;

            ostream << augmented_control_points[i][j].x() << ' '
                    << augmented_control_points[i][j].y() << ' '
                    << augmented_control_points[i][j].z() << '\n';
        }
    }
    ostream.close();

    m_surface = CubicBSpline(m_surface.m_nv, m_surface.m_nu, std::move(augmented_control_points));
    m_step++;
}
 */
