/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Surfaces/Surface.h"
#include <iostream>
#include <ranges>

ParametricSurface::ParametricSurface(Options const& options)
    : m_options(options)
    , m_epsilon(1e-5) {};

Vec3 ParametricSurface::f_u(Vec2 const& p) const { return nf_u(p); }

Vec3 ParametricSurface::f_v(Vec2 const& p) const { return nf_v(p); }

Vec3 ParametricSurface::f_uu(Vec2 const& p) const { return nf_uu(p); }

Vec3 ParametricSurface::f_uv(Vec2 const& p) const { return nf_uv(p); }

Vec3 ParametricSurface::f_vv(Vec2 const& p) const { return nf_vv(p); }

Vec3 ParametricSurface::normal(Vec2 const& p) const
{
    return f_u(p).cross(f_v(p)).normalized();
}

Vec3 ParametricSurface::nf_u(Vec2 const& p) const
{
    Vec2 const delta(m_epsilon, 0.);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_v(Vec2 const& p) const
{
    Vec2 const delta(0., m_epsilon);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_uu(Vec2 const& p) const
{
    Vec2 const delta(2. * m_epsilon, 0.);

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
    // return ((f(p + Vec2(2.0 * m_epsilon, 0)) - f(p)) - (f(p) - f(p + Vec2(-2.0 * m_epsilon, 0)))) / 4.0 / m_epsilon / m_epsilon;
}

Vec3 ParametricSurface::nf_vv(Vec2 const& p) const
{
    Vec2 const delta(0., 2. * m_epsilon);

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
    // return ((f(p + Vec2(0, 2.0 * m_epsilon)) - f(p)) - (f(p) - f(p + Vec2(0, -2.0 * m_epsilon)))) / 4.0 / m_epsilon / m_epsilon;
}

Vec3 ParametricSurface::nf_uv(Vec2 const& p) const
{
    // clang-format off
    return (
        f(p + Vec2( m_epsilon,  m_epsilon))
      + f(p + Vec2(-m_epsilon, -m_epsilon))
      - f(p + Vec2(-m_epsilon,  m_epsilon))
      - f(p + Vec2( m_epsilon, -m_epsilon))
    ) / (4. * m_epsilon * m_epsilon);
    // clang-format on
}

Vec2 ParametricSurface::rescale(Vec2 const& p) const
{
    return Vec2(
        p.x() / (2. * PI) * (m_uMax - m_uMin) + m_uMin,
        p.y() * (m_uMax - m_uMin) + m_vMin);
}

void ParametricSurface::generate_search_grid(int nu, int nv)
{
    m_grid2D = decltype(m_grid2D)((nu + 1) * (nv + 1), Vec2::Zero());
    m_grid3D = decltype(m_grid3D)((nu + 1) * (nv + 1), Vec3::Zero());

    double du = (m_uMax - m_uMin) / nu;
    double dv = (m_vMax - m_vMin) / nv;

    int ct = 0;
    for (int i = 0; i <= nu; i++) {
        for (int j = 0; j <= nv; j++) {
            Vec2 p(m_uMin + du * i, m_vMin + dv * j);

            m_grid2D[ct] = p;
            m_grid3D[ct] = f(p);

            ct++;
        }
    }
}

Vec2 ParametricSurface::closest_point(Vec3 const& p, size_t max_iterations) const
{
    double cur;
    double mi = (p - m_grid3D[0]).norm();
    int id = 0;

    for (auto&& [i, q] : enumerate(m_grid3D)) {
        cur = (p - q).norm();

        if (cur < mi) {
            mi = cur;
            id = i;
        }
    }

    return closest_point(p, m_grid2D[id], max_iterations);
}

Vec2 ParametricSurface::closest_point(Vec3 const& p, Vec2 const& guess, size_t max_iterations) const
{
    Vec2 xk = guess;

    static constexpr auto lambda = 0.1;

    for (auto i = 0uz; i < max_iterations; i++) {
        Vec3 fp = f(xk) - p;
        Vec3 fu = f_u(xk);
        Vec3 fv = f_v(xk);
        Vec3 fuu = f_uu(xk);
        Vec3 fuv = f_uv(xk);
        Vec3 fvv = f_vv(xk);

        Vec2 Jf(fu.dot(fp), fv.dot(fp));
        Vec2 Jg(
            0.,
            lambda * (1. / (m_vMax - xk.y()) - 1. / (xk.y() - m_vMin)));

        // if (!m_options.newton_log_barrier)
            Jg.y() = 0.;

        Vec2 J = 2. * Jf + Jg;

        Eigen::Matrix2d Hf = Eigen::Matrix2d::Zero();
        Hf << fuu.dot(fp) + fu.dot(fu),
            fuv.dot(fp) + fu.dot(fv),
            fuv.dot(fp) + fu.dot(fv),
            fvv.dot(fp) + fv.dot(fv);

        Eigen::Matrix2d Hg = Eigen::Matrix2d::Zero();

        // if (m_options.newton_log_barrier)
            Hg(1, 1) = lambda * (1. / std::pow(xk.y() - m_vMin, 2.) + 1. / std::pow(m_vMax - xk.y(), 2.));

        Eigen::Matrix2d H = 2. * Hf + Hg;

        Vec2 dx = H.inverse() * J;

        if (dx.norm() <= 1e-5)
            break;

        xk -= dx;
    }

    return xk;
}
