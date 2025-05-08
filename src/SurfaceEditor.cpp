#if USE_QPMAD
#    include <qpmad/solver.h>
#else
#    include "epigraph.hpp"
#endif

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>

#include "Composer.h"
#include "Config.h"
#include "SurfaceEditor.h"
#include "Timer.h"

namespace fs = std::filesystem;

void write_winding_path(size_t step, ParametricSurface const& surface, std::vector<LocalFrame> const& path)
{
    auto outpath = std::format("{}/{}/path/step-{}.txt", Config::out_directory, Config::experiment, step);
    auto ostream = std::ofstream(outpath);

    for (auto&& [i, frame] : enumerate(path)) {
        Vec3 shadow = surface.f(frame.shadow);
        std::println(
            ostream,
            "{}, {}, {}, {}, {}, {}, {}, {}",
            i, frame.position, frame.distance, frame.path_distance, frame.angle, shadow, frame.direction, frame.normal);
    }

    ostream.flush();
}

void write_control_points(size_t step, std::vector<std::vector<Vec3>> const& control_points)
{
    auto outpath = std::format("{}/{}/spline/step-{}.txt", Config::out_directory, Config::experiment, step + 1);
    auto ostream = std::ofstream(outpath);

    std::println(ostream, "{} {}", control_points.front().size(), control_points.size());

    for (auto i = 0uz; i < control_points.size(); i++) {
        for (auto j = 0uz; j < control_points.front().size(); j++) {
            std::println(ostream, "{}", control_points[i][j]);
        }
    }
    ostream.close();
}

auto surface_optimization(
    int num_control_points,
    int num_optimizable_particles,
    Eigen::MatrixXd M,
    Eigen::VectorXd b) -> Eigen::VectorXd
{
    auto solver_timer = Timer("Solver");
#if USE_QPMAD
    auto solver = qpmad::Solver();
    Eigen::VectorXd deltaC(3 * num_control_points);
    Eigen::MatrixXd static H = Eigen::MatrixXd::Identity(3 * num_control_points, 3 * num_control_points);
    Eigen::VectorXd static h = Eigen::VectorXd::Zero(3 * num_control_points);
    Eigen::VectorXd Aub = Eigen::VectorXd::Ones(num_optimizable_particles);

    for (auto i = 0; i < num_optimizable_particles; i++)
        Aub[i] = Infinity;

    // Solving: min (deltaC^T * H * deltaC + deltaC^T * h) subject to (b < M * deltaC < Aub)
    // Set H = Identity, h = <0>, Aub = <Infinity> to get a least-norm QP.
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

    return deltaC;
}

auto get_off_surface_particles(std::vector<LocalFrame>& path) -> std::vector<std::pair<Vec2, Vec3>>
{
    auto particles = std::vector<std::pair<Vec2, Vec3>> {};

    for (auto&& [i, frame] : enumerate(path)) {
        if (frame.distance < 1e-3)
            continue;

        particles.push_back(std::make_pair(frame.shadow, frame.position));
    }

    return particles;
}

auto get_in_surface_particles(
    ParametricSurface& surface,
    int step,
    std::vector<Quad>& quads) -> std::vector<std::pair<Vec2, Vec3>>
{
    auto particles = std::vector<std::pair<Vec2, Vec3>> {};
    std::sort(quads.begin(), quads.end(), [](Quad a, Quad b) { return a.area < b.area; });

    auto N = quads.size();
    auto q1 = N % 2 ? quads[N / 4].area : (quads[N / 4].area + quads[(N / 4) - 1].area) / 2.;
    auto q3 = N % 2 ? quads[(3 * N) / 4].area : (quads[(3 * N) / 4].area + quads[((3 * N) / 4) - 1].area) / 2.;
    auto iqr = q3 - q1;

    auto has_outliers = std::any_of(quads.begin() + (3 * N) / 4, quads.end(),
                                    [&iqr, &q3](Quad a) {
                                        // NOTE: See the NOTE on 'outlier_threshold' in Config.pp
                                        return a.area >= Config::outlier_threshold * iqr + q3;
                                    });

    if (!has_outliers)
        return particles;

    // We use IQR to see if outliers exist, however, we do not use all outliers to pull-in the surface.
    // In the paper, we talk about setting some user-defined epsilon to look at only the quads within epsilon
    // of the surface area of the largest quad.

    // Not explained in the paper:
    // Here, I define epsilon as a logistic function mapping the current optimizer step to some value in [alpha, 1.].
    // This is monotonically increasing the epsilon value used, i.e. it becomes harder for anything but the largest
    // quadrilateral to be optimized as time goes on.
    // This was an ad hoc solution based on my gut feeling, ideally I think we should just up the IQR threshold.
    auto constexpr alpha = 0.90;
    auto constexpr k = 0.125;
    auto threshold = alpha + 2. * (1. - alpha) * (1. / (1. + std::exp(-step * k / 2.)) - 0.5);
    auto lower_bound = quads.back().area * threshold;

    auto it = std::next(quads.begin(), (3 * N) / 4);
    for (; it != quads.end(); it = std::next(it)) {
        if (it->area >= lower_bound)
            break;
    }

    for (; it != quads.end(); it = std::next(it)) {
        auto const& quad = *it;

        Vec3 m0 = (quad.p1 + quad.p0) / 2.;
        Vec3 m1 = (quad.p3 + quad.p2) / 2.;
        Vec3 m2 = (quad.p0 + quad.p3) / 2.;
        Vec3 m3 = (quad.p1 + quad.p2) / 2.;

        Vec3 rd0 = Vec3::Zero();
        Vec3 rd1 = Vec3::Zero();

        Vec3 d0 = (m1 - m0).normalized();
        Vec3 d1 = (m3 - m2).normalized();

        // Regularly sample particles along the sides of the midpoint quadrilateral
        // NOTE: See NOTE on 'pull_in_particles_per_path' in Config.cpp
        for (auto i = 0; i < Config::pull_in_particles_per_path; i++) {
            double t = ((double)i) / (Config::pull_in_particles_per_path - 1.);
            Vec3 w0 = m0 + t * d0;
            Vec3 w1 = m2 + t * d1;

            Vec2 uv0 = surface.closest_point(w0);
            Vec2 uv1 = surface.closest_point(w1);

            particles.push_back(std::make_pair(uv0, w0));
            particles.push_back(std::make_pair(uv1, w1));
        }
    }

    return particles;
}

void SurfaceEditor::step(
    double timestep,
    double spring_constant,
    double damping_coefficient,
    double epsilon)
{
    std::println("[Editor] Begin step {}", m_step + 1);

    fs::copy_file(
        std::format("{}/{}.txt", Config::data_directory, Config::stem),
        std::format("{}/{}/spline/step-0.txt", Config::out_directory, Config::experiment),
        fs::copy_options::overwrite_existing);

    auto composer = Composer(m_surface, m_num_revolutions, m_num_paths, m_num_particles);
    auto path = composer.simulate(timestep, spring_constant, damping_coefficient, epsilon, m_step);
    auto optimizable_particles = std::vector<std::pair<Vec2, Vec3>> {};

    write_winding_path(m_step, m_surface, path);

    bool do_push_out = (Config::alternate_push_pull && m_step % 2 == 0) || (!Config::alternate_push_pull && Config::push_out);

    auto matrix_timer = Timer("Matrix", false);
    if (do_push_out)
        optimizable_particles = get_off_surface_particles(path);

    bool do_pull_in = (Config::alternate_push_pull && m_step % 2 == 1) || (!Config::alternate_push_pull && Config::pull_in);

    if (do_pull_in)
        optimizable_particles = get_in_surface_particles(m_surface, m_step, composer.quads);

    // Can't RAII with anon scope here since we need defs for multiple variables below
    matrix_timer.print();

    auto m = m_surface.m_nv * m_surface.m_nu;
    auto k = optimizable_particles.size();

    std::println("[Editor] Iteration {} found {} optimizable particles ({}).", m_step + 1, k, do_push_out ? "push-out" : "pull-in");

    Eigen::MatrixXd M(k, 3 * m);
    Eigen::VectorXd b(k);

    for (auto&& [i, pair] : enumerate(optimizable_particles)) {
        auto const& [shadow, world] = pair;

        Vec3 normal = m_surface.normal(shadow);

        // Flip normal if pull-in
        normal *= do_pull_in ? -1. : 1.;

        M.row(i) = normal.transpose() * m_surface.jacobian(shadow);
        b(i) = (world - m_surface.f(shadow)).dot(normal);
    }

    auto deltaC = surface_optimization(m, k, M, b);

    std::println("[Editor] Solved least-norm QP (deltaC L2: {})", deltaC.norm());

    auto augmented_control_points = std::vector(m_surface.m_nv, std::vector<Vec3>(m_surface.m_nu, Vec3::Zero()));

    {
        auto index = 0;
        for (auto i = 0; i < m_surface.m_nv; i++) {
            for (auto j = 0; j < m_surface.m_nu; j++) {
                augmented_control_points[i][j] = m_surface.points()[i][j] + deltaC.segment<3>(3 * index);
                index++;
            }
        }
    }

    write_control_points(m_step, augmented_control_points);

    auto modified_surface = CubicBSpline(m_surface.m_nv, m_surface.m_nu, std::move(augmented_control_points));

    std::println("[Editor] Hausdorff: {}", hausdorff_distance(m_surface, modified_surface));
    std::println("[Editor] pct. diff: {}", percent_difference(m_surface, modified_surface));

    std::swap(m_surface, modified_surface);

    m_step++;
}
