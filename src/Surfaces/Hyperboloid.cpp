/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Surfaces/Surface.h"
#include "utils.h"

Hyperboloid::Hyperboloid()
    : ParametricSurface(0., 2. * PI, 0., 1.) {};

Vec3 Hyperboloid::f(Vec2 const& p) const
{
    return { (p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
             (p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
             p.y() };
}

Vec3 Hyperboloid::f_u(Vec2 const& p) const
{
    return {
        -(p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
        (p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
        0.
    };
}

Vec3 Hyperboloid::f_v(Vec2 const& p) const
{
    return {
        (2. * p.y() - 1.) * cos(p.x()),
        (2. * p.y() - 1.) * sin(p.x()),
        1.
    };
}

Vec3 Hyperboloid::f_uu(Vec2 const& p) const
{
    return {
        -(p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
        -(p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
        0.
    };
}

Vec3 Hyperboloid::f_uv(Vec2 const& p) const
{
    return {
        -(2. * p.y() - 1.) * sin(p.x()),
        (2. * p.y() - 1.) * cos(p.x()),
        0.
    };
}

Vec3 Hyperboloid::f_vv(Vec2 const& p) const
{
    return {
        2. * cos(p.x()),
        2. * sin(p.x()),
        0.
    };
}