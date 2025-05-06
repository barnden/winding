/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <numbers>
#include <numeric>

#include "Composer.h"

using std::numbers::pi;

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

void Composer::generate_winding_order()
{
    if (Config::use_winding_order) {
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
    } else {
        m_winding_order.resize(m_initial.size());
        std::iota(m_winding_order.begin(), m_winding_order.end(), 0);
    }
}