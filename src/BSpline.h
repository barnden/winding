/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include <fstream>
#include <iostream>

#include "Surface.h"

class CubicBSpline : public ParametricSurface {
    std::vector<std::vector<Vec3>> m_points;
    Eigen::Matrix4d m_basis;

    [[nodiscard]] inline std::tuple<double, double, int, int> get_uv(Vec2 const& p) const
    {
        Vec2 uv = p;

        while (uv.x() < 0.)
            uv.x() += m_uMax;

        auto iu = (int)uv.x();
        auto iv = (int)uv.y();

        if (iv == m_nv - 1)
            iv--;

        if (iv == 1)
            iv++;

        uv.x() -= iu;
        uv.y() -= iv;

        iu %= m_nu;

        return std::make_tuple(uv.x(), uv.y(), iu, iv);
    }

    [[nodiscard]] inline Eigen::RowVector3d get(int v, int u) const
    {
        u %= m_nu;
        if (u < 0)
            u += m_nu;

        return m_points[v][u];
    }

    template <typename UFunc, typename FFunc>
    [[gnu::flatten]] Vec3 interpolate(UFunc&& get_u, FFunc&& get_v, Vec2 const& p) const
    {
        auto [u, v, iu, iv] = get_uv(p);

        Eigen::RowVector4d U = get_u(u) * m_basis;
        Eigen::RowVector4d V = get_v(v) * m_basis;

        // clang-format off
        Eigen::MatrixXd P = (Eigen::MatrixXd(4, 3)
                                 << V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu - 2), get(iv - 1, iu - 2), get(iv, iu - 2), get(iv + 1, iu - 2)).finished(),
                                    V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu - 1), get(iv - 1, iu - 1), get(iv, iu - 1), get(iv + 1, iu - 1)).finished(),
                                    V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu + 0), get(iv - 1, iu + 0), get(iv, iu + 0), get(iv + 1, iu + 0)).finished(),
                                    V * (Eigen::MatrixXd(4, 3) << get(iv - 2, iu + 1), get(iv - 1, iu + 1), get(iv, iu + 1), get(iv + 1, iu + 1)).finished())
                                .finished();
        // clang-format on

        return U * P;
    }

public:
    int m_nu;
    int m_nv;
    explicit CubicBSpline(std::string const& file)
    {
        if (file.empty())
            return;

        m_basis = 1. / 6. * (Eigen::Matrix4d() << -1., 3., -3., 1., 3., -6., 3., 0., -3., 0., 3., 0., 1., 4., 1., 0.).finished();
        read(file);
    }

    CubicBSpline(int nv, int nu, decltype(m_points)&& data)
        : m_points(std::move(data))
        , m_nu(nu)
        , m_nv(nv)
    {
        m_uMin = 0.;
        m_uMax = nu;
        m_vMin = 2.;
        m_vMax = nv - 1.;
    }

    void read(std::string const& file)
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

    [[nodiscard]] Vec3 f(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
            [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
            p);
    }

    [[nodiscard]] Vec3 f_u(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
            [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
            p);
    }

    [[nodiscard]] Vec3 f_v(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
            [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
            p);
    }

    [[nodiscard]] Vec3 f_uv(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
            [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
            p);
    }

    [[nodiscard]] Vec3 f_uu(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(6. * u, 2., 0., 0.); },
            [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
            p);
    }

    [[nodiscard]] Vec3 f_vv(Vec2 const& p) const override
    {
        return interpolate(
            [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
            [](double v) { return Eigen::RowVector4d(6. * v, 2., 0., 0.); },
            p);
    }

    [[nodiscard]] Eigen::MatrixXd jacobian(Vec2 const& p) const
    {
        Eigen::Matrix3d static const I = Eigen::Matrix3d::Identity();
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3 * m_nv * m_nu);

        auto const [u, v, iu, iv] = get_uv(p);

        Eigen::RowVector4d U = Eigen::RowVector4d(u * u * u, u * u, u, 1.) * m_basis;
        Eigen::RowVector4d V = Eigen::RowVector4d(v * v * v, v * v, v, 1.) * m_basis;

        J.block<3, 3>(0, 3 * (iv - 2) * m_nu + (iu - 2)) = V(0) * U(0) * I;
        J.block<3, 3>(0, 3 * (iv - 2) * m_nu + (iu - 1)) = V(0) * U(1) * I;
        J.block<3, 3>(0, 3 * (iv - 2) * m_nu + (iu + 0)) = V(0) * U(2) * I;
        J.block<3, 3>(0, 3 * (iv - 2) * m_nu + (iu + 1)) = V(0) * U(3) * I;

        J.block<3, 3>(0, 3 * (iv - 1) * m_nu + (iu - 2)) = V(1) * U(0) * I;
        J.block<3, 3>(0, 3 * (iv - 1) * m_nu + (iu - 1)) = V(1) * U(1) * I;
        J.block<3, 3>(0, 3 * (iv - 1) * m_nu + (iu + 0)) = V(1) * U(2) * I;
        J.block<3, 3>(0, 3 * (iv - 1) * m_nu + (iu + 1)) = V(1) * U(3) * I;

        J.block<3, 3>(0, 3 * (iv + 0) * m_nu + (iu - 2)) = V(2) * U(0) * I;
        J.block<3, 3>(0, 3 * (iv + 0) * m_nu + (iu - 1)) = V(2) * U(1) * I;
        J.block<3, 3>(0, 3 * (iv + 0) * m_nu + (iu + 0)) = V(2) * U(2) * I;
        J.block<3, 3>(0, 3 * (iv + 0) * m_nu + (iu + 1)) = V(2) * U(3) * I;

        J.block<3, 3>(0, 3 * (iv + 1) * m_nu + (iu - 2)) = V(3) * U(0) * I;
        J.block<3, 3>(0, 3 * (iv + 1) * m_nu + (iu - 1)) = V(3) * U(1) * I;
        J.block<3, 3>(0, 3 * (iv + 1) * m_nu + (iu + 0)) = V(3) * U(2) * I;
        J.block<3, 3>(0, 3 * (iv + 1) * m_nu + (iu + 1)) = V(3) * U(3) * I;

        return J;
    }

    auto const& points() const { return m_points; }

    std::vector<float> get_control_polygon();
    std::vector<float> get_control_points();
};