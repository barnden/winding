/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include "Composer.h"
#include <numeric>

auto Composer::intersect(
    Vec2 up1,
    Vec2 up2,
    Vec2 down1,
    Vec2 down2,
    Vec2& intersection,
    double& t_up,
    double& t_down) -> bool
{

    double u = (m_surface.u_max() - m_surface.u_min());

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

    return true;
}

auto Composer::intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down) -> std::vector<IntersectionNode*>
{
    std::vector<IntersectionNode*> result;

    auto i = 0uz;
    auto j = down.size() - 1uz;
    double t_up;
    double t_down;
    while (i < up.size() - 1 && j > 0) {
        Vec2 intersection;
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

void Composer::generate_winding_order()
{
#if 1
    m_winding_order.resize(m_initial.size());
    std::iota(m_winding_order.begin(), m_winding_order.end(), 0);
#else

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
#endif
}