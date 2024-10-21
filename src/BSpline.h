/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include "Surface.h"

class CubicBSpline : public ParametricSurface {
    std::vector<std::vector<Vec3>> m_points;
    Eigen::Matrix4d m_basis;

    [[nodiscard]] inline std::tuple<double, double, int, int> get_uv(Vec2 const& p) const;
    [[nodiscard]] inline Eigen::RowVector3d get(int v, int u) const;

    template <typename UFunc, typename FFunc>
    [[gnu::flatten]] Vec3 interpolate(UFunc&& get_u, FFunc&& get_v, Vec2 const& p) const;

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
        m_basis = 1. / 6. * (Eigen::Matrix4d() << -1., 3., -3., 1., 3., -6., 3., 0., -3., 0., 3., 0., 1., 4., 1., 0.).finished();

        m_uMin = 0.;
        m_uMax = nu;
        m_vMin = 2.;
        m_vMax = nv - 1.;
    }

    void read(std::string const& file);

    [[nodiscard]] Vec3 f(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_u(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_v(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_uv(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_uu(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_vv(Vec2 const& p) const override;

    [[nodiscard]] Eigen::MatrixXd jacobian(Vec2 const& p) const;

    auto const& points() const { return m_points; }

    std::vector<float> get_control_polygon();
    std::vector<float> get_control_points();
};