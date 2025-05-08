/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include <numbers>

#include "Surfaces/Surface.h"
#include "utils.h"

using std::numbers::pi;

Torus::Torus(double r1, double r2)
    : ParametricSurface(0., 2. * pi, 0., 2. * pi - 0.001)
    , m_r1(r1)
    , m_r2(r2) {};

Vec3 Torus::f(Vec2 const& p) const
{
    return {
        (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        (m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        -m_r2 * sin(p.x())
    };
}

Vec3 Torus::f_u(Vec2 const& p) const
{
    return {
        -m_r2 * sin(p.x()) * cos(p.y()),
        -m_r2 * sin(p.x()) * sin(p.y()),
        -m_r2 * cos(p.x())
    };
}

Vec3 Torus::f_v(Vec2 const& p) const
{
    return {
        -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        0
    };
}

Vec3 Torus::f_uu(Vec2 const& p) const
{
    return {
        -m_r2 * cos(p.x()) * cos(p.y()),
        -m_r2 * cos(p.x()) * sin(p.y()),
        m_r2 * sin(p.x())
    };
}

Vec3 Torus::f_uv(Vec2 const& p) const
{
    return {
        m_r2 * sin(p.x()) * sin(p.y()),
        -m_r2 * sin(p.x()) * cos(p.y()),
        0
    };
}

Vec3 Torus::f_vv(Vec2 const& p) const
{
    return {
        -(m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        0
    };
}