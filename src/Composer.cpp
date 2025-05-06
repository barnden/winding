/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <chrono>
#include <format>
#include <fstream>
#include <numbers>
#include <numeric>

#include "Composer.h"
#include "Config.h"
#include "Simulator.h"

using std::numbers::pi;

IntersectionList::IntersectionList()
{
    head = new IntersectionNode();
    rear = new IntersectionNode();

    head->is_end = true;
    rear->is_end = true;
}

void IntersectionList::insert(IntersectionNode* node, IntersectionNode* IntersectionNode::* prev, IntersectionNode* IntersectionNode::* next, double IntersectionNode::* value) const
{
    IntersectionNode* p = head;
    while (p->*next != rear && (p->*next)->*value < node->*value) {
        p = p->*next;
    }

    node->*prev = p;
    node->*next = p->*next;

    (p->*next)->*prev = node;
    p->*next = node;
}

void IntersectionList::interpolate(double t, double& dist, double& angle, IntersectionNode* IntersectionNode::* next, double IntersectionNode::* value) const
{
    IntersectionNode* p = head->*next;

    double lambda = (t - p->*value) / ((p->*next)->*value - p->*value);
    dist = lambda * (p->*next)->dist_score + (1. - lambda) * p->dist_score;
    angle = lambda * (p->*next)->angle_score + (1. - lambda) * p->angle_score;
}

Composer::Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
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

auto Composer::intersect(
    Vec2 up1,
    Vec2 up2,
    Vec2 down1,
    Vec2 down2,
    Vec2& intersection,
    double& t_up,
    double& t_down) -> bool
{
    auto u = (m_surface.u_max() - m_surface.u_min());
    auto const wrap = [&u](auto x) { return u * std::floor(x / u + 0.5); };

    auto a1 = up1;
    auto a2 = up2;
    auto b1 = down1;
    auto b2 = down2;

    a2.x() -= wrap(a2.x() - a1.x());
    b2.x() -= wrap(b2.x() - b1.x());

    auto ma = (a1.x() + a2.x()) / 2.;
    auto mb = (b1.x() + b2.x()) / 2.;

    b1.x() -= wrap(mb - ma);
    b2.x() -= wrap(mb - ma);

    Vec2 deltaA = a1 - a2;
    Vec2 deltaB = b2 - b1;
    Vec2 deltaS = a1 - b1;

    auto dA = deltaA.x() * deltaB.y() - deltaA.y() * deltaB.x();

    if (std::abs(dA) < 1e-4)
        return false;

    auto dAu = deltaA.x() * deltaS.y() - deltaA.y() * deltaS.x();
    t_down = dAu / dA;

    if (t_down < 0. || t_down >= 1.)
        return false;

    auto dAt = deltaS.x() * deltaB.y() - deltaS.y() * deltaB.x();
    t_up = dAt / dA;

    if (t_up < 0. || t_up >= 1.)
        return false;

    up1.x() -= u * std::floor(up1.x() / u);
    up2.x() -= u * std::floor(up2.x() / u);

    intersection = (1. - t_up) * a1 + t_up * a2;

    return true;
}

auto Composer::intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down) -> std::vector<IntersectionNode*>
{
    std::vector<IntersectionNode*> result;

    auto i = 0uz;
    auto j = down.size() - 1uz;
    while (i < up.size() - 1 && j > 0) {
        Vec2 intersection;
        double t_up;
        double t_down;

        if (intersect(up[i], up[i + 1], down[j - 1], down[j], intersection, t_up, t_down)) {
            auto* p = new IntersectionNode();

            p->t_up = t_up + i;
            p->t_down = t_down + j - 1.;
            p->parametric = intersection;
            p->point = m_surface.f(intersection);

            result.push_back(p);
        }

        if (up[i + 1].y() < down[j - 1].y()) {
            i++;
        } else if (up[i + 1].y() > down[j - 1].y()) {
            j--;
        } else {
            i++;
            j--;
        }
    }

    return result;
}

void Composer::score_ends(IntersectionNode* p)
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

auto Composer::score(IntersectionNode* p) -> double
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

auto Composer::simulate(
    double timestep,
    double spring_constant,
    double damping_coefficient,
    double epsilon,
    int step)
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

    std::println("[Timing] Execution time: {} us", duration.count());

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

                    std::println(ostream, "v {}\nv {}\nv {}\nv {}", quad.p0, quad.p1, quad.p2, quad.p3);
                }
                p = p->next_up;
            }
        }

        std::println("[Debug] Max quad area: {}", max_area);

        int cnt = 1;
        for (auto& ls : up_list) {
            auto p = ls.head;
            while (p != ls.rear) {
                if (p->next_up && p->next_up->next_down && p->next_up->next_down->prev_up && p->next_up->next_down->prev_up->prev_down == p) {
                    std::println(ostream, "f {} {} {} {}", cnt, cnt + 1, cnt + 2, cnt + 3);
                    cnt += 4;
                }
                p = p->next_up;
            }
        }
    }

    auto end_meshing = std::chrono::steady_clock::now();
    auto duration_meshing = std::chrono::duration_cast<std::chrono::microseconds>(end_meshing - start_meshing);
    std::println("[Timing] Meshing: {} us", duration_meshing.count());

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

void Composer::generate_winding_order()
{
    if (!Config::use_winding_order) {
        m_winding_order.resize(m_initial.size());
        std::iota(m_winding_order.begin(), m_winding_order.end(), 0);

        return;
    }

    OrderNode* nodes = new OrderNode[m_initial.size()];
    OrderNode* z_pos = nodes;
    OrderNode* z_neg = nodes + 1;
    OrderNode* min_p = nodes;

    z_pos->prev = z_pos;
    z_neg->prev = z_neg;
    z_pos->next = z_pos;
    z_neg->next = z_neg;

    nodes[0].path = 0;
    nodes[1].path = 1;
    nodes[0].u = fmod(m_initial[0][0].x(), 2. * pi);

    for (auto i = 2uz; i < m_initial.size(); i++) {
        OrderNode& cur = nodes[i];
        cur.path = i;
        cur.u = fmod(m_initial[i][0].x(), 2. * pi);

        if (cur.u < 0.)
            cur.u += 2. * pi;

        if (i % 2 == 0) {
            cur.next = z_pos->next;
            cur.prev = z_pos;
            z_pos->next = nodes + i;
            cur.next->prev = nodes + i;

            if (min_p->u > cur.u)
                min_p = nodes + i;
        } else {
            cur.next = z_neg->next;
            cur.prev = z_neg;
            z_neg->next = nodes + i;
            cur.next->prev = nodes + i;
        }
    }

    z_pos = min_p;
    bool next_pos = true;
    OrderNode* next = nodes;
    double nextu;

    while (next) {
        m_winding_order.push_back(next->path);
        nextu = fmod(m_initial[next->path].back().x() - 0.1, 2. * pi);

        if (nextu < 0.)
            nextu += 2. * pi;

        if (next_pos) {
            next_pos = !next_pos;

            if (z_pos->next == z_pos) {
                next = z_neg;
                continue;
            }

            next->prev->next = next->next;
            next->next->prev = next->prev;

            if (next == z_pos)
                z_pos = z_pos->prev;

            if (nextu < z_neg->u || nextu >= z_neg->next->u) {
                next = z_neg->next;
            } else {
                next = z_neg->next->next;

                while (nextu < next->u)
                    next = next->next;
            }

            continue;
        }

        next_pos = !next_pos;
        if (z_neg->next == z_neg)
            break;

        next->prev->next = next->next;
        next->next->prev = next->prev;

        if (next == z_neg)
            z_neg = z_neg->prev;

        if (nextu < z_pos->u || nextu >= z_pos->next->u) {
            next = z_pos->next;
        } else {
            next = z_pos->next->next;

            while (nextu < next->u)
                next = next->next;
        }
    }

    delete[] nodes;
}


void Composer::interpolate(
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
