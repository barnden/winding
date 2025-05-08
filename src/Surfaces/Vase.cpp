/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include <numbers>

#include "Surfaces/Surface.h"
#include "utils.h"

using std::numbers::pi;

Vase::Vase()
    : ParametricSurface(0., 2. * pi, -1., 1.) {};

Vec3 Vase::f(Vec2 const& p) const
{
    return {
        (0.5 - 0.4 * p.y() * (1. + p.y()) * (1. - p.y())) * cos(p.x()),
        (0.5 + 0.4 * p.y() * (1. + p.y()) * (1. - p.y())) * sin(p.x()),
        p.y()
    };
}

Vec3 Vase::f_u(Vec2 const& p) const
{
    return {
        -(0.5 - 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * sin(p.x()),
        (0.5 + 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * cos(p.x()),
        0
    };
}

Vec3 Vase::f_v(Vec2 const& p) const
{
    return {
        (-0.4 + 1.2 * p.y() * p.y()) * cos(p.x()),
        (0.4 - 1.2 * p.y() * p.y()) * sin(p.x()),
        1.
    };
}

Vec3 Vase::f_uu(Vec2 const& p) const
{
    return {
        -(0.5 - 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * cos(p.x()),
        -(0.5 + 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * sin(p.x()),
        0.
    };
}

Vec3 Vase::f_uv(Vec2 const& p) const
{
    return {
        -(-0.4 + 1.2 * p.y() * p.y()) * sin(p.x()),
        (0.4 - 1.2 * p.y() * p.y()) * cos(p.x()),
        0.
    };
}

Vec3 Vase::f_vv(Vec2 const& p) const
{
    return {
        2.4 * p.y() * cos(p.x()),
        -2.4 * p.y() * sin(p.x()),
        0.
    };
}