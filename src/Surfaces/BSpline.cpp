/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include <fstream>

#include "BSpline.h"
#include "utils.h"
#include <iostream>

std::tuple<double, double, int, int> CubicBSpline::get_uv(Vec2 const& p) const
{
    Vec2 uv = p;

    uv.x() = std::fmod(uv.x(), m_uMax);
    uv.y() = std::fmod(uv.y(), m_vMax);

    while (uv.x() < m_uMin)
        uv.x() += m_uMax;

    while (uv.y() < m_vMin)
        uv.y() += m_vMax;

    auto iu = (int)uv.x();
    auto iv = (int)uv.y();

    if (iv == m_nv - 1)
        iv--;

    if (iv == 1)
        iv++;

    uv.x() -= iu;
    uv.y() -= iv;

    iu %= m_nu;
    iv %= m_nv;

    if (iv == 0)
        iv = 2;

    return std::make_tuple(uv.x(), uv.y(), iu, iv);
}

Eigen::RowVector3d CubicBSpline::get(int v, int u) const
{
    u %= m_nu;
    if (u < 0)
        u += m_nu;

    return m_points[v][u];
}

template <typename UFunc, typename FFunc>
[[gnu::flatten]] Vec3 CubicBSpline::interpolate(UFunc&& get_u, FFunc&& get_v, Vec2 const& p) const
{
    auto [u, v, iu, iv] = get_uv(p);

    Eigen::RowVector4d U = get_u(u) * m_basis;
    Eigen::RowVector4d V = get_v(v) * m_basis;

    // clang-format off
    Eigen::MatrixXd P = (Eigen::MatrixXd(4, 3) <<
                                V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu - 2), get(iv - 1, iu - 2), get(iv, iu - 2), get(iv + 1, iu - 2)).finished(),
                                V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu - 1), get(iv - 1, iu - 1), get(iv, iu - 1), get(iv + 1, iu - 1)).finished(),
                                V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu + 0), get(iv - 1, iu + 0), get(iv, iu + 0), get(iv + 1, iu + 0)).finished(),
                                V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu + 1), get(iv - 1, iu + 1), get(iv, iu + 1), get(iv + 1, iu + 1)).finished())
                            .finished();
    // clang-format on

    return U * P;
}

void CubicBSpline::read(std::string const& file)
{
    auto stream = std::ifstream(file);

    if (stream.fail()) {
        throw "Failed to open file " + file + '\n';
    }

    stream >> m_nu >> m_nv;

    m_points = std::vector<std::vector<Vec3>>(m_nv, std::vector<Vec3>(m_nu));

    m_uMin = 0.;
    m_uMax = m_nu;
    m_vMin = 2.;
    m_vMax = m_nv - 1.;

    for (auto i = 0; i < m_nv; i++) {
        for (auto j = 0; j < m_nu; j++) {
            stream >> m_points[i][j].x() >> m_points[i][j].y() >> m_points[i][j].z();
            // m_points[i][j + 2 * m_nu] = (m_points[i][j + m_nu] = m_points[i][j]);
        }
    }

    stream.close();
}

[[nodiscard]] Vec3 CubicBSpline::f(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[nodiscard]] Vec3 CubicBSpline::f_u(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[nodiscard]] Vec3 CubicBSpline::f_v(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
        p);
}

[[nodiscard]] Vec3 CubicBSpline::f_uv(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
        [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
        p);
}

[[nodiscard]] Vec3 CubicBSpline::f_uu(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(6. * u, 2., 0., 0.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[nodiscard]] Vec3 CubicBSpline::f_vv(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(6. * v, 2., 0., 0.); },
        p);
}

Eigen::MatrixXd CubicBSpline::jacobian(Vec2 const& p) const
{
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3 * m_nv * m_nu);

    auto const get_index = [&](int v, int u) {
        u %= m_nu;
        if (u < 0)
            u += m_nu;

        return 3 * (v * m_nu + u);
    };
    auto const [u, v, iu, iv] = get_uv(p);

    Eigen::RowVector4d U = Eigen::RowVector4d(u * u * u, u * u, u, 1.) * m_basis;
    Eigen::RowVector4d V = Eigen::RowVector4d(v * v * v, v * v, v, 1.) * m_basis;

    for (auto l = 0; l < 4; l++)
        for (auto m = 0; m < 4; m++)
            J.block<3, 3>(0, get_index(iv + (l - 2), iu + (m - 2))).diagonal().array() = V(l) * U(m);

    return J;
}