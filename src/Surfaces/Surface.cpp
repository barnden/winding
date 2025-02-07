/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Surfaces/Surface.h"
#include "Config.h"

#include <iostream>
#include <random>
#include <ranges>

ParametricSurface::ParametricSurface(double epsilon)
    : m_epsilon(epsilon) { };

ParametricSurface::ParametricSurface(double u_min, double u_max, double v_min, double v_max, double epsilon)
    : m_uMin(u_min)
    , m_uMax(u_max)
    , m_vMin(v_min)
    , m_vMax(v_max)
    , m_epsilon(epsilon) { };

Vec3 ParametricSurface::f_u(Vec2 const& p) const { return nf_u(p); }

Vec3 ParametricSurface::f_v(Vec2 const& p) const { return nf_v(p); }

Vec3 ParametricSurface::f_uu(Vec2 const& p) const { return nf_uu(p); }

Vec3 ParametricSurface::f_uv(Vec2 const& p) const { return nf_uv(p); }

Vec3 ParametricSurface::f_vv(Vec2 const& p) const { return nf_vv(p); }

template <typename T>
auto sgn(T val) -> int
{
    return (T(0) < val) - (val < T(0));
}

auto ParametricSurface::sdf(Vec3 const& p) const -> double
{
    Vec2 P_parametric = closest_point(p);
    Vec3 v = p - f(P_parametric);
    auto sign = sgn(normal(P_parametric).dot(v));

    return sign * v.norm();
}

Vec3 ParametricSurface::normal(Vec2 const& p) const
{
    return f_u(p).cross(f_v(p)).normalized();
}

Vec3 ParametricSurface::nf_u(Vec2 const& p) const
{
    static Vec2 const delta(m_epsilon, 0.);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_v(Vec2 const& p) const
{
    static Vec2 const delta(0., m_epsilon);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_uu(Vec2 const& p) const
{
    static Vec2 const delta(2. * m_epsilon, 0.);

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
}

Vec3 ParametricSurface::nf_vv(Vec2 const& p) const
{
    static Vec2 const delta(0., 2. * m_epsilon);

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
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

void ParametricSurface::generate_search_grid(int nu, int nv) const
{
    if (!m_grid2D.empty())
        return;

    m_grid2D.clear();
    m_grid2D.resize((nu + 1) * (nv + 1));

    m_grid3D.clear();
    m_grid3D.resize((nu + 1) * (nv + 1));

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

    m_bvh = BVH(&m_grid3D, std::max(nu, nv));
}

Vec2 ParametricSurface::closest_point(Vec3 const& p, size_t max_iterations) const
{
    generate_search_grid(1024, 1024);
    auto id = m_bvh.closest_point(p);

    return closest_point(p, m_grid2D[id], max_iterations);
}

Vec2 ParametricSurface::closest_point(Vec3 const& p, Vec2 const& guess, size_t max_iterations) const
{
    Vec2 xk = guess;

    for (auto i = 0uz; i < max_iterations; i++) {
        Vec3 fp = f(xk) - p;
        Vec3 fu = f_u(xk);
        Vec3 fv = f_v(xk);
        Vec3 fuu = f_uu(xk);
        Vec3 fuv = f_uv(xk);
        Vec3 fvv = f_vv(xk);

        Vec2 J(fu.dot(fp), fv.dot(fp));

        Eigen::Matrix2d H = Eigen::Matrix2d::Zero();
        H << fuu.dot(fp) + fu.dot(fu),
            fuv.dot(fp) + fu.dot(fv),
            fuv.dot(fp) + fu.dot(fv),
            fvv.dot(fp) + fv.dot(fv);

        Vec2 dx = H.inverse() * J;

        if (dx.norm() <= 1e-5)
            break;

        xk -= dx;
    }

    return xk;
}

auto hausdorff_distance(ParametricSurface const& surfaceA, ParametricSurface const& surfaceB) -> double
{
    static auto rd = std::random_device();
    static auto generator = std::mt19937(rd());
    static auto distribution = std::uniform_real_distribution(0., 100.);

    surfaceA.generate_search_grid(1024, 1024);
    surfaceB.generate_search_grid(1024, 1024);

    auto hausdorff = -std::numeric_limits<double>::infinity();
    auto threshold = (2. / std::log2(surfaceA.m_grid3D.size()));

    for (auto&& point : surfaceA.m_grid3D) {
        if (distribution(generator) > threshold)
            continue;

        auto d_point = (point - surfaceB.f(surfaceB.closest_point(point))).norm();

        hausdorff = std::max(hausdorff, d_point);
    }

    return hausdorff;
}