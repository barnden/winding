/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once
#ifndef COMPOSER_H
#    define COMPOSER_H

#    include "Config.h"
#    include "Simulator.h"
#    include "Surfaces/Surface.h"
#    include "utils.h"

#    include <format>
#    include <fstream>
#    include <iostream>
#    include <vector>

#    include <chrono>

#    define SF_OFF 0.001
#    define ANGLE_SPEED_CONST 2.0

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;
    ParametricSurface const& m_surface;

    int m_num_paths;
    int m_num_particles;
    double m_l_turn;

    struct OrderNode {
        int path;
        double u;

        OrderNode* next;
        OrderNode* prev;
    };

    struct IntersectionNode {
        Vec2 parametric = Vec2::Zero();
        Vec3 point = Vec3::Zero();

        double angle_score = 0.;
        double dist_score = 0.;
        double t_up = 0.;
        double t_down = 0.;

        IntersectionNode* prev_up = nullptr;
        IntersectionNode* next_up = nullptr;

        IntersectionNode* prev_down = nullptr;
        IntersectionNode* next_down = nullptr;

        bool is_end = false;
    };

    class IntersectionListUp {
    public:
        IntersectionNode *head, *rear;
        IntersectionListUp()
        {
            head = new IntersectionNode();
            rear = new IntersectionNode();
            head->next_up = rear;
            rear->prev_up = head;
            head->t_up = -1.0E100;
            rear->t_up = 1.0E100;
            head->is_end = true;
            rear->is_end = true;
        }
        void insert(IntersectionNode* node)
        {
            IntersectionNode* p = head;
            while (p->next_up != rear && p->next_up->t_up < node->t_up) {
                p = p->next_up;
            }
            node->prev_up = p;
            node->next_up = p->next_up;
            p->next_up->prev_up = node;
            p->next_up = node;
        }
        void interpolate(double t, double& dist, double& angle) const
        {
            IntersectionNode* p = head->next_up;

            double lambda = (t - p->t_up) / (p->next_up->t_up - p->t_up);
            dist = lambda * p->next_up->dist_score + (1 - lambda) * p->dist_score;
            angle = lambda * p->next_up->angle_score + (1 - lambda) * p->angle_score;
        }
    };

    class IntersectionListDown {
    public:
        IntersectionNode *head, *rear;
        IntersectionListDown()
        {
            head = new IntersectionNode();
            rear = new IntersectionNode();
            head->next_down = rear;
            rear->prev_down = head;
            head->t_down = -1.0E100;
            rear->t_down = 1.0E100;
            head->is_end = true;
            rear->is_end = true;
        }
        void insert(IntersectionNode* node)
        {
            IntersectionNode* p = head;
            while (p->next_down != rear && p->next_down->t_down < node->t_down) {
                p = p->next_down;
            }
            node->prev_down = p;
            node->next_down = p->next_down;
            p->next_down->prev_down = node;
            p->next_down = node;
        }
        void interpolate(double t, double& dist, double& angle) const
        {
            IntersectionNode* p = head->next_down;

            double lambda = (t - p->t_down) / (p->next_down->t_down - p->t_down);
            dist = lambda * p->next_down->dist_score + (1 - lambda) * p->dist_score;
            angle = lambda * p->next_down->angle_score + (1 - lambda) * p->angle_score;
        }
    };

public:
    struct LocalFrame {
        Vec3 position;
        Vec3 normal;
        Vec3 direction;
        Vec2 shadow;

        double distance;
        double path_distance;
        double angle;
    };

    Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
        : m_surface(surface)
        , m_num_paths(num_paths)
        , m_num_particles(num_particles)
        , m_l_turn(0.3)
    {
        auto angle = (surface.v_max() - surface.v_min()) / (surface.u_max() - surface.u_min()) / num_revolutions;
        auto du = (surface.u_max() - surface.u_min()) / num_paths;
        auto ddv = (surface.v_max() - surface.v_min()) / (num_particles + 1);
        auto ddu = ddv / angle;

        m_initial = std::vector<std::vector<Vec2>>(2 * num_paths, std::vector<Vec2>(num_particles));

        for (auto i = 0; i < num_paths; i++) {
            Vec2 p(surface.u_min() + i * du, surface.v_min());
            Vec2 q(surface.u_min() + i * du, surface.v_min());

            Vec3 cur = surface.f(p) + SF_OFF * surface.normal(p);

            for (auto j = 0; j < num_particles; j++) {
                cur = surface.f(p) + SF_OFF * surface.normal(p);
                m_initial[2 * i][num_particles - 1 - j] = p;
                m_initial[2 * i + 1][j] = q;

                p += Vec2(ddu, ddv);
                q += Vec2(-ddu, ddv);
            }
        }
    }

    void generate_winding_order();

    auto intersect(
        Vec2 up1, Vec2 up2,
        Vec2 down1, Vec2 down2,
        Vec2& intersection,
        double& t_up, double& t_down) -> bool;

    auto intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down) -> std::vector<IntersectionNode*>;

    decltype(auto) score(IntersectionNode* p)
    {
        if (p->prev_up->is_end || p->prev_down->is_end || p->next_up->is_end || p->next_down->is_end)
            return 0.;

        Vec3 v1 = (p->prev_up->point - p->next_up->point).normalized();
        Vec3 v2 = (p->prev_down->point - p->next_down->point).normalized();
        p->angle_score = abs(v1.dot(v2));
        p->dist_score = ((p->prev_up->point - p->point).norm()
                         + (p->prev_down->point - p->point).norm()
                         + (p->next_up->point - p->point).norm()
                         + (p->next_down->point - p->point).norm())
                        / 2.;

        return p->dist_score;
    }

    void score_ends(IntersectionNode* p)
    {
        if (p->is_end)
            return;

        double sum_dist = 0.;
        double sum_angle = 0.;
        int ct = 0;

        if (!p->prev_up->is_end) {
            sum_dist += p->prev_up->dist_score;
            sum_angle += p->prev_up->angle_score;

            ct++;
        }

        if (!p->prev_down->is_end) {
            sum_dist += p->prev_down->dist_score;
            sum_angle += p->prev_down->angle_score;

            ct++;
        }

        if (!p->next_up->is_end) {
            sum_dist += p->next_up->dist_score;
            sum_angle += p->next_up->angle_score;

            ct++;
        }

        if (!p->next_down->is_end) {
            sum_dist += p->next_down->dist_score;
            sum_angle += p->next_down->angle_score;

            ct++;
        }

        if (ct != 0) {
            p->dist_score = sum_dist / ct;
            p->angle_score = sum_angle / ct;
        }
    }

    void interpolate(
        LocalFrame& frame,
        std::vector<IntersectionListUp> const& up,
        std::vector<IntersectionListDown> const& down,
        int i,
        double t)
    {
        if (i % 2) {
            up[i / 2].interpolate(t, frame.path_distance, frame.angle);
        } else {
            down[i / 2].interpolate(t, frame.path_distance, frame.angle);
        }
    }

    double max_area = -INFINITY;

    struct Quad {
        Vec3 p0;
        Vec3 p1;
        Vec3 p2;
        Vec3 p3;

        double area;

        Quad()
            : p0(Vec3::Zero())
            , p1(Vec3::Zero())
            , p2(Vec3::Zero())
            , p3(Vec3::Zero())
            , area(0.) { };

        Quad(Vec3 const& p0, Vec3 const& p1, Vec3 const& p2, Vec3 const& p3)
            : p0(p0)
            , p1(p1)
            , p2(p2)
            , p3(p3)
        {
            area = 0.5 * ((p0 - p1).cross(p0 - p2).norm() + (p3 - p1).cross(p3 - p2).norm());
        }
    };
    std::vector<Quad> quads;
    Quad max_quad;

    decltype(auto) simulate(
        double timestep = 0.05,
        double spring_constant = 1'000'000.,
        double damping_coefficient = 200.,
        double epsilon = 0.001,
        int step = 0)
    {
        generate_winding_order();
        quads = std::vector<Quad>();

        auto paths = std::vector<std::vector<Vec2>>(m_initial.size());
        auto orders = std::vector<std::vector<int>>(m_initial.size());

        auto start = std::chrono::steady_clock::now();

        for (auto&& [j, i] : enumerate(m_winding_order)) {
            auto simulator = OffSurface(m_surface, m_initial[i]);

            simulator.spring_constant() = spring_constant;
            simulator.damping_coefficient() = damping_coefficient;
            simulator.timestep() = timestep;
            simulator.epsilon() = epsilon;

            simulator.simulate(1000);
            simulator.mapping();

            paths[i] = simulator.p();
            orders[i] = simulator.r();
        }

        auto end = std::chrono::steady_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;

        auto start_meshing = std::chrono::steady_clock::now();
        std::vector<IntersectionListUp> up_list(m_initial.size() / 2);
        std::vector<IntersectionListDown> down_list(m_initial.size() / 2);

        for (auto i = 0u; i < m_initial.size() / 2; i++) {
            for (auto j = 0u; j < m_initial.size() / 2; j++) {
                auto result = intersect(paths[2 * i + 1], paths[2 * j]);

                for (auto&& p : result) {
                    up_list[i].insert(p);
                    down_list[j].insert(p);
                }
            }
        }

        for (auto&& ls : up_list) {
            for (auto p = ls.head->next_up; p != ls.rear; p = p->next_up)
                score(p);
        }

        for (auto&& ls : up_list) {
            score_ends(ls.head->next_up);
            score_ends(ls.rear->prev_up);
        }

        for (auto&& ls : down_list) {
            score_ends(ls.head->next_down);
            score_ends(ls.rear->prev_down);
        }

        {
            auto outpath = std::format("{}/{}/mesh/step-{}.obj", Config::out_directory, Config::experiment, step);
            auto ostream = std::ofstream(outpath);
            for (auto& ls : up_list) {
                auto p = ls.head;
                while (p != ls.rear) {
                    if (p->next_up && p->next_up->next_down && p->next_up->next_down->prev_up && p->next_up->next_down->prev_up->prev_down == p) {
                        auto quad = Quad { p->point,
                                           p->next_up->point,
                                           p->next_up->next_down->point,
                                           p->next_up->next_down->prev_up->point };

                        quads.push_back(quad);
                        if (quad.area > max_area) {
                            max_area = quad.area;
                            max_quad = quad;
                        }

                        ostream << "v " << quad.p0.transpose() << '\n'
                                << "v " << quad.p1.transpose() << '\n'
                                << "v " << quad.p2.transpose() << '\n'
                                << "v " << quad.p3.transpose() << '\n';
                    }
                    p = p->next_up;
                }
            }

            std::cout << "max quad area: " << max_area << '\n';

            int cnt = 1;
            for (auto& ls : up_list) {
                auto p = ls.head;
                while (p != ls.rear) {
                    if (p->next_up && p->next_up->next_down && p->next_up->next_down->prev_up && p->next_up->next_down->prev_up->prev_down == p) {
                        ostream << "f " << cnt << " " << cnt + 1 << " " << cnt + 2 << " " << cnt + 3 << std::endl;
                        cnt += 4;
                    }
                    p = p->next_up;
                }
            }
        }

        auto end_meshing = std::chrono::steady_clock::now();
        auto duration_meshing = std::chrono::duration_cast<std::chrono::microseconds>(end_meshing - start_meshing);
        std::cout << "Meshing time: " << duration_meshing.count() << " microseconds" << std::endl;

        auto motion = std::vector<LocalFrame> {};
        motion.reserve(2 * m_num_paths * m_num_particles);
        for (auto&& [j, i] : enumerate(m_winding_order)) {
            // FIXME: In most cases order[idx] is equal to 1; so a lot of this computation is actually wastefu

            auto const& path = paths[i];
            auto const& order = orders[i];
            for (auto idx = 0uz; idx < path.size();) {
                auto const N = order[idx];
                size_t i0 = idx;
                size_t i1 = std::clamp(i0 + N, 0uz, path.size() - 1uz);

                Vec3 n0 = m_surface.normal(path[i0]);
                Vec3 n1 = m_surface.normal(path[i1]);

                Vec3 p0 = m_surface.f(path[i0]) + SF_OFF * n0;
                Vec3 p1 = m_surface.f(path[i1]) + SF_OFF * n1;

                Vec3 direction = (p1 - p0).normalized();

                for (auto k = 0; k < N; k++) {
                    auto t = ((double)k) / ((double)N);
                    Vec3 world = (1. - t) * p0 + t * p1;
                    Vec3 normal = ((1. - t) * n0 + t * n1).normalized();
                    Vec2 parametric = path[idx];
                    double distance = 0.;

                    if (N != 1) {
                        Vec3 surface = m_surface.f(parametric);
                        distance = (world - surface).norm();
                    }

                    auto frame = LocalFrame {
                        world,
                        normal,
                        direction,
                        parametric,
                        distance,
                        0., 0.
                    };

                    interpolate(frame, up_list, down_list, i, idx);

                    motion.push_back(frame);
                    idx++;
                }
            }
        }

        return motion;
    }

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
    decltype(auto) surface() const { return m_surface; }
};

#endif