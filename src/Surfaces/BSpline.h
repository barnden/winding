/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include <iostream>
#include <string>
#include "Surfaces/Surface.h"
#include "utils.h"
#include "Config.h"

class CubicBSpline : public ParametricSurface {
    std::vector<std::vector<Vec3>> m_points;

    [[nodiscard]] inline auto get_uv(Vec2 const& p) const -> std::tuple<double, double, int, int>;
    [[nodiscard]] inline auto get(int v, int u) const -> Eigen::RowVector3d;

    template <typename UFunc, typename FFunc>
    [[gnu::flatten]] Vec3 interpolate(UFunc&& get_u, FFunc&& get_v, Vec2 const& p) const;

public:
    int m_nu;
    int m_nv;

    [[nodiscard]] static auto basis() -> Eigen::Matrix4d
    {
        return 1. / 6. * (Eigen::Matrix4d() << -1., 3., -3., 1., 3., -6., 3., 0., -3., 0., 3., 0., 1., 4., 1., 0.).finished();
    }

    explicit CubicBSpline()
        : ParametricSurface(1e-5)
    {
        auto file = std::format("{}/{}.txt", Config::data_directory, Config::stem);
        std::cout << "reading: " << file << '\n';
        read(file);
    }

    CubicBSpline(int nv, int nu, decltype(m_points)&& data)
        : ParametricSurface(0., nu, 2., nv - 1.)
        , m_points(std::move(data))
        , m_nu(nu)
        , m_nv(nv)
    {}

    void read(std::string const& file);

    [[nodiscard]] Vec3 f(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_u(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_v(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_uv(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_uu(Vec2 const& p) const override;
    [[nodiscard]] Vec3 f_vv(Vec2 const& p) const override;

    [[nodiscard]] Eigen::MatrixXd jacobian(Vec2 const& p) const;

    auto const& points() const { return m_points; }

    [[nodiscard]] auto get_control_polygon() -> std::vector<float>;
    [[nodiscard]] auto get_control_points() -> std::vector<float>;
};