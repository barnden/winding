#pragma once

#include "utils.h"

#include "Simulator.h"
#include "Surfaces/Surface.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <numeric>

#define SURFACE_OFFSET 0.001
#define SPEED_DEFINED_BY_LENGTH
#define ANGLE_SPEED_CONST 2.0

class Composer2 {
    std::vector<std::vector<Vec2>> m_initial_paths;
    std::vector<int> m_winding_order;

    std::shared_ptr<Options> m_options;
    ParametricSurface const& m_surface;

    int m_num_paths;
    int m_num_particles_per_path;
    double m_l_turn;

public:
    struct LocalFrame {
        Vec3 position;
        Vec3 normal;
        Vec3 direction;
        Vec2 shadow;

        double l;
        int type;

        double distance;
        double path_distance;
        double angle;
    };
    Eigen::Matrix<double, 12, 1> max_quad;

    Composer2(std::shared_ptr<Options> const& options, ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
        : m_options(options)
        , m_surface(surface)
        , m_num_paths(num_paths)
        , m_num_particles_per_path(num_particles)
        , m_l_turn(0.3)
    {
        auto angle = (surface.v_max() - surface.v_min()) / (surface.u_max() - surface.u_min()) / num_revolutions;
        auto du = (surface.u_max() - surface.u_min()) / num_paths;
        auto ddv = (surface.v_max() - surface.v_min()) / (num_particles + 1);
        auto ddu = ddv / angle;

        m_initial_paths = std::vector<std::vector<Vec2>>(2 * num_paths, std::vector<Vec2>(num_particles));

        for (auto i = 0; i < num_paths; i++) {
            Vec2 p(surface.u_min() + i * du, surface.v_min());
            Vec2 q(surface.u_min() + i * du, surface.v_min());

            Vec3 cur = surface.f(p) + SURFACE_OFFSET * surface.normal(p);

            for (auto j = 0; j < num_particles; j++) {
                cur = surface.f(p) + SURFACE_OFFSET * surface.normal(p);
                m_initial_paths[2 * i][num_particles - 1 - j] = p;
                m_initial_paths[2 * i + 1][j] = q;

                p += Vec2(ddu, ddv);
                q += Vec2(-ddu, ddv);
            }
        }
    }

    void generate_winding_order()
    {
        m_winding_order = decltype(m_winding_order)(m_initial_paths.size(), 0);
        std::iota(std::begin(m_winding_order), std::end(m_winding_order), 0);
    }

    auto intersect(
        Vec2 up1,
        Vec2 up2,
        Vec2 down1,
        Vec2 down2,
        Vec2& intersection,
        double& t_up,
        double& t_down) -> bool
    {
        auto u = m_surface.u_max() - m_surface.u_min();

        while (up1.x() - down1.x() > (u / 2.)) {
            up1.x() -= u;
            up2.x() -= u;
        }

        while (down1.x() - up1.x() > (u / 2.)) {
            up1.x() += u;
            up2.x() += u;
        }

        Vec2 up_min = up1.cwiseMin(up2);
        Vec2 down_max = down1.cwiseMax(down2);

        if (up_min.x() > down_max.x())
            return false;

        Vec2 up_max = up1.cwiseMax(up2);
        Vec2 down_min = down1.cwiseMin(down2);

        if (up_max.x() < down_min.x())
            return false;

        if (up_min.y() > down_max.y())
            return false;

        if (up_max.y() < down_min.y())
            return false;

        Vec2 delta_down = down2 - down1;
        Vec2 delta_up = up1 - up2;
        auto dT = delta_down.x() * delta_up.y() - delta_up.x() * delta_down.y();

        if (std::abs(dT) < 1e-20)
            return false;

        Vec2 delta_p = up1 - down1;
        t_down = (delta_up.y() * delta_p.x() - delta_up.x() * delta_p.y()) / dT;

        if (t_down < 0. || t_down >= 1.)
            return false;

        t_up = (delta_down.x() * delta_p.y() - delta_down.y() * delta_p.x()) / dT;
        if (t_up < 0. || t_up >= 1.)
            return false;

        intersection = (1. - t_up) * up1 + t_up * up2;

        return true;
    }

    struct Intersection {
        double t_up;
        double t_down;

        Vec3 world;
        Vec2 parametric;

        struct {
            double angle;
            double path_distance;
        } score;
    };

    auto intersect_paths(std::vector<Vec2> const& up, std::vector<Vec2> const& down) -> std::vector<Intersection>
    {
        std::vector<Intersection> result;

        auto i = 0uz;
        auto j = down.size() - 1uz;

        double t_up;
        double t_down;

        while (i < up.size() - 1uz && j > 0uz) {
            Vec2 intersection;

            if (intersect(up[i], up[i + 1], down[j - 1], down[j], intersection, t_up, t_down)) {
                result.emplace_back(
                    t_up + i,
                    t_down + j - 1.,
                    m_surface.f(intersection),
                    intersection);
            }

            if (up[i + 1].y() < down[j - 1].y()) {
                i++;
                continue;
            }

            if (up[i + 1].y() > down[j - 1].y()) {
                j--;
                continue;
            }

            i++;
            j--;
        }

        return result;
    }

    auto simulate(
        double dt = 0.05,
        double ksp = 1'000'000,
        double kdp = 200.,
        double eps = 0.001,
        int step = -1)
    {
        generate_winding_order();

        auto paths = std::vector<std::vector<Vec2>>(m_initial_paths.size());
        auto orders = std::vector<std::vector<int>>(m_initial_paths.size());

        for (auto&& [j, i] : enumerate(m_winding_order)) {
            // FIXME: We can parallelize this
            auto simulator = OffSurface(m_options, m_surface, m_initial_paths[i]);

            simulator.ksp() = ksp;
            simulator.kdp() = kdp;
            simulator.dt() = dt;
            simulator.eps() = eps;

            simulator.simulate(1000);
            simulator.mapping();

            paths[i] = simulator.p();
            orders[i] = simulator.r();
        }

        std::vector<std::vector<Intersection>> path_intersections(m_initial_paths.size() / 2, std::vector<Intersection>(0));
        std::vector<std::vector<size_t>> up_list(m_initial_paths.size() / 2);
        std::vector<std::vector<size_t>> down_list(m_initial_paths.size() / 2);

        for (auto i = 0uz; i < m_initial_paths.size() / 2; i++) {
            for (auto j = 0uz; j < m_initial_paths.size() / 2; j++) {
                auto result = intersect_paths(paths[2 * i + 1], paths[2 * j]);
                path_intersections[i].insert(path_intersections[i].end(), result.begin(), result.end());
            }
        }

        for (auto i = 0uz; i < m_initial_paths.size() / 2; i++) {
            auto& intersections = path_intersections[i];

            up_list[i] = std::vector<size_t>(intersections.size());
            down_list[i] = std::vector<size_t>(intersections.size());

            std::iota(std::begin(up_list[i]), std::end(up_list[i]), 0uz);
            std::iota(std::begin(down_list[i]), std::end(down_list[i]), 0uz);

            auto const projection = [&intersections](std::size_t i) -> Intersection& { return intersections[i]; };

            std::ranges::sort(
                up_list[i],
                [](auto const& a, auto const& b) -> bool { return a.t_up < b.t_up; },
                projection);

            std::ranges::sort(
                down_list[i],
                [](auto const& a, auto const& b) -> bool { return a.t_down < b.t_down; },
                projection);
        }

        namespace views = std::views;

        // Compute path scores: angles and distances between paths
        for (auto&& [intersections, up_path, down_path] : views::zip(path_intersections, up_list, down_list)) {
            for (auto&& [idx, intersection] : enumerate(intersections) | views::take(intersections.size() - 1) | views::drop(1)) {
                Vec3 const& world = intersections[idx].world;

                Vec3 const& prev_up = intersections[up_path[idx - 1]].world;
                Vec3 const& next_up = intersections[up_path[idx + 1]].world;

                Vec3 const& prev_down = intersections[down_path[idx - 1]].world;
                Vec3 const& next_down = intersections[down_path[idx + 1]].world;

                Vec3 v1 = (prev_up - next_up).normalized();
                Vec3 v2 = (prev_down - next_down).normalized();

                double distance = 0.5 * ((prev_up - world).norm() + (prev_down - world).norm() + (next_up - world).norm() + (next_down - world).norm());

                intersections[idx].score = {
                    .angle = std::abs(v1.dot(v2)),
                    .path_distance = distance
                };
            }
        }
        {
            std::ofstream ofs(m_options->out_path + "/" + m_options->experiment + "/obj/step-" + std::to_string(step) + ".obj");
            std::stringstream obj_verts {};
            std::stringstream obj_faces {};

            int vertexCount = 0;
            for (auto&& [intersections, up_path, down_path] : views::zip(path_intersections, up_list, down_list)) {
                for (auto i = 0uz; i < intersections.size() - 1; i++) {
                    auto v2 = std::find(std::begin(down_path), std::end(down_path), up_path[i + 1]);
                    if (std::next(v2) == std::end(down_path))
                        continue;

                    auto i2 = *std::next(v2);

                    auto v3 = std::find(std::begin(up_path), std::end(up_path), i2);
                    if (v3 == std::begin(up_path))
                        continue;

                    auto i3 = *std::prev(v3);

                    auto v4 = std::find(std::begin(down_path), std::end(down_path), i3);
                    if (v4 == std::begin(down_path) || *std::prev(v4) != i)
                        continue;

                    Vec3 p0 = intersections[up_path[i]].world;
                    Vec3 p1 = intersections[up_path[i + 1]].world;
                    Vec3 p2 = intersections[i2].world;
                    Vec3 p3 = intersections[i3].world;

                    obj_verts << "v " << p0.transpose() << '\n'
                              << "v " << p1.transpose() << '\n'
                              << "v " << p2.transpose() << '\n'
                              << "v " << p3.transpose() << '\n';

                    obj_faces << "f " << vertexCount++ << " " << vertexCount++ << " " << vertexCount++ << " " << vertexCount++ << '\n';
                }
            }

            ofs << obj_verts.str() << obj_faces.str();
        }

        std::cout << "wrote quadmesh\n";
        exit(0);
        // for (auto&& path : up_list) {
        //     for (auto i = 0uz; i < path.size(); i++) {
        //         Vec3 p0 = path_intersections[path[i]];
        //     }
        // }

        return std::vector<LocalFrame> {};
    }
};