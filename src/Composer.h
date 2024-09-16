#pragma once
#ifndef COMPOSER_H
#    define COMPOSER_H

#    include "Surface.h"
#    include <vector>

#    define SF_OFF 0.001

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;

    struct OrderNode {
        int path;
        double u;

        OrderNode* next;
        OrderNode* prev;
    };

public:
    Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
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

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
};

#endif