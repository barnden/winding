/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Surfaces/Surface.h"
#include "utils.h"

Spring::Spring(double r1, double r2, double kh)
    : m_r1(r1)
    , m_r2(r2)
    , m_kh(kh)
{
    m_uMin = 0.;
    m_uMax = 2. * PI;
    m_vMin = 0.;
    m_vMax = 4. * PI;
}

Vec3 Spring::f(Vec2 const& p) const
{
    return {
        (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        (m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        -m_r2 * sin(p.x()) + m_kh * p.y()
    };
}

Vec3 Spring::f_u(Vec2 const& p) const
{
    return {
        -m_r2 * sin(p.x()) * cos(p.y()),
        -m_r2 * sin(p.x()) * sin(p.y()),
        -m_r2 * cos(p.x())
    };
}

Vec3 Spring::f_v(Vec2 const& p) const
{
    return {
        -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        m_kh
    };
}

Vec3 Spring::f_uu(Vec2 const& p) const
{
    return {
        -m_r2 * cos(p.x()) * cos(p.y()),
        -m_r2 * cos(p.x()) * sin(p.y()),
        m_r2 * sin(p.x())
    };
}

Vec3 Spring::f_uv(Vec2 const& p) const
{
    return {
        m_r2 * sin(p.x()) * sin(p.y()),
        -m_r2 * sin(p.x()) * cos(p.y()),
        0
    };
}

Vec3 Spring::f_vv(Vec2 const& p) const
{
    return {
        -(m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
        -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
        0
    };
}