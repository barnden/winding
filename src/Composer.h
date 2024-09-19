#pragma once
#ifndef COMPOSER_H
#    define COMPOSER_H

#    include "Surface.h"
#    include <vector>

#    define SF_OFF 0.001

#    include "utils.h"

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;
    ParametricSurface const& m_surface;

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

public:
    Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
        : m_surface(surface)
    {
        auto angle = (surface.m_vMax - surface.m_vMin) / (surface.m_uMax - surface.m_uMin) / num_revolutions;
        auto du = (surface.m_uMax - surface.m_uMin) / num_paths;
        auto ddv = (surface.m_vMax - surface.m_vMin) / (num_particles + 1);
        auto ddu = ddv / angle;

        m_initial = std::vector<std::vector<Vec2>>(2 * num_paths, std::vector<Vec2>(num_particles));

        for (auto i = 0; i < num_paths; i++) {
            Vec2 p(surface.m_uMin + i * du, surface.m_vMin);
            Vec2 q(surface.m_uMin + i * du, surface.m_vMin);

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

    void generate_winding_order()
    {
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
        nodes[0].u = fmod(m_initial[0][0].x(), 2. * PI);

        for (auto i = 2uz; i < m_initial.size(); i++) {
            OrderNode& cur = nodes[i];
            cur.path = i;
            cur.u = fmod(m_initial[i][0].x(), 2. * PI);

            if (cur.u < 0.)
                cur.u += 2. * PI;

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
            nextu = fmod(m_initial[next->path].back().x() - 0.1, 2. * PI);

            if (nextu < 0.)
                nextu += 2. * PI;

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

    decltype(auto) simulate(
        double num_revolutions = 0.1,
        int num_paths = 1,
        int num_particles = 200,
        double dt = 0.05,
        double ksp = 1'000'000.,
        double kdp = 200.,
        double eps = 0.001)
    {
        generate_winding_order();

        auto paths = std::vector<std::vector<Vec2>>(m_initial.size());
        auto order = std::vector<std::vector<int>>(m_initial.size());

        for (auto&& [j, i] : enumerate(m_winding_order)) {
            auto simulator = OffSurface(m_surface, m_initial[i]);

            simulator.ksp() = ksp;
            simulator.kdp() = kdp;
            simulator.dt() = dt;
            simulator.eps() = eps;

            simulator.simulate(1000);
            simulator.mapping();

            // auto const& result = simulator.p();
            // auto pT = std::vector<Eigen::VectorXd>(result.size());
            // for (auto&& [i, p] : enumerate(result))
            //     pT[i] = Eigen::VectorXd { { (double)simulator.r()[i], p.x(), p.y() } };
            // parr[i] = pT;

            paths[i] = simulator.p();
            order[i] = simulator.r();

            // auto progress = 100. * (j + 1.) / m_winding_order.size();
            // std::cout << round(progress) << "%\n\n";
        }

        for (auto i = 0; i < m_initial.size() / 2; i++) {
            for (auto j = 0; j < m_initial.size() / 2; j++) { }
        }
    }

    decltype(auto) intersect(
        Vec2 up1,
        Vec2 up2,
        Vec2 down1,
        Vec2 down2,
        Vec2& intersection,
        double& t_up,
        double& t_down)
    {
        double half_u = (m_surface.m_uMax - m_surface.m_uMin) / 2.;

        while (up1.x() - down1.x() > half_u) {
            up1.x() -= 2. * half_u;
            up2.x() -= 2. * half_u;
        }

        while (down1.x() - up1.x() > half_u) {
            up1.x() += 2. * half_u;
            up2.x() += 2. * half_u;
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
        double dT = delta_down.x() * delta_up.y() - delta_up.x() * delta_down.y();

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
    }

    decltype(auto) intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down)
    {
        std::vector<IntersectionNode*> result;

        int i = 0;
        int j = down.size() - 1;
        double t_up;
        double t_down;
        while (i < up.size() - 1 && j > 0) {
            Vec2 intersection;
            if (intersect(up[i], up[i + 1], down[j - 1], down[j], intersection, t_up, t_down)) {
                auto *p = new IntersectionNode();

                p->t_up = t_up + i;
                p->t_down = t_down + j - 1.;
                p->parametric = intersection;
                p->point = m_surface.f(intersection);

                result.push_back(p);
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

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
};

#endif