/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once
#ifndef SURFACE_H
#    define SURFACE_H

#    include <Eigen/Dense>
#    include <memory>
#    include <vector>

#    include "utils.h"

class ParametricSurface {
    Options const& m_options;

public:
    ParametricSurface() = delete;
    ParametricSurface(Options const& options);

    ParametricSurface(ParametricSurface const& other)
        : m_options(other.m_options)
    {
        m_grid2D = other.m_grid2D;
        m_grid3D = other.m_grid3D;

        m_uMax = other.m_uMax;
        m_uMin = other.m_uMin;
        m_vMax = other.m_vMax;
        m_vMin = other.m_vMin;

        m_epsilon = other.m_epsilon;
    }

    [[nodiscard]] virtual Vec3 f(Vec2 const& p) const = 0;
    [[nodiscard]] virtual Vec3 f_u(Vec2 const& p) const;
    [[nodiscard]] virtual Vec3 f_v(Vec2 const& p) const;
    [[nodiscard]] virtual Vec3 f_uu(Vec2 const& p) const;
    [[nodiscard]] virtual Vec3 f_uv(Vec2 const& p) const;
    [[nodiscard]] virtual Vec3 f_vv(Vec2 const& p) const;
    [[nodiscard]] virtual Vec3 normal(Vec2 const& p) const;

    [[nodiscard]] Vec3 nf_u(Vec2 const& p) const;
    [[nodiscard]] Vec3 nf_v(Vec2 const& p) const;
    [[nodiscard]] Vec3 nf_uu(Vec2 const& p) const;
    [[nodiscard]] Vec3 nf_uv(Vec2 const& p) const;
    [[nodiscard]] Vec3 nf_vv(Vec2 const& p) const;

    [[nodiscard]] Vec2 rescale(Vec2 const& p) const;

    void generate_search_grid(int nu, int nv);

    [[nodiscard]] virtual Vec2 closest_point(Vec3 const& p, size_t max_iterations = 1000) const;
    [[nodiscard]] Vec2 closest_point(Vec3 const& p, Vec2 const& guess, size_t max_iterations = 1000) const;

    mutable std::vector<Vec2> m_grid2D;
    mutable std::vector<Vec3> m_grid3D;

    double m_uMax;
    double m_uMin;
    double m_vMax;
    double m_vMin;

    double m_epsilon;

    std::unique_ptr<float[]> get_display_data(int& len_out, int nc_res = 500, int nh_res = 500);
};

class Hyperboloid : public ParametricSurface {
public:
    Hyperboloid(Options const& options);

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
};

class Torus : public ParametricSurface {
    double m_r1;
    double m_r2;

public:
    Torus(Options const& options, double r1 = 1., double r2 = 0.3);

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
};

class Spring : public ParametricSurface {
    double m_r1;
    double m_r2;
    double m_kh;

public:
    Spring(Options const& options, double r1 = 1., double r2 = 0.3, double kh = 0.15);

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
};

class Vase : public ParametricSurface {
public:
    Vase(Options const& options);

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
};

class TrefoilKnot : public ParametricSurface {
public:
    TrefoilKnot(Options const& options);
    virtual Vec3 f(Vec2 const& p) const final;
};

#endif
