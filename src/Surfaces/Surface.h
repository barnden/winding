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
public:
    ParametricSurface();

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

    [[nodiscard]] Vec2 closest_point(Vec3 const& p) const;
    [[nodiscard]] Vec2 closest_point(Vec3 const& p, Vec2 const& guess, int max_iterations = 1000) const;

    std::vector<Vec2> m_grid2D;
    std::vector<Vec3> m_grid3D;

    double m_uMax;
    double m_uMin;
    double m_vMax;
    double m_vMin;

    double m_epsilon;

    std::unique_ptr<float[]> get_display_data(int& len_out, int nc_res = 500, int nh_res = 500);
};

class Hyperboloid : public ParametricSurface {
public:
    Hyperboloid();

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
    Torus(double r1 = 1., double r2 = 0.3);

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
    Spring(double r1 = 1., double r2 = 0.3, double kh = 0.15);

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
};

class Vase : public ParametricSurface {
public:
    Vase();

    virtual Vec3 f(Vec2 const& p) const final;
    virtual Vec3 f_u(Vec2 const& p) const final;
    virtual Vec3 f_v(Vec2 const& p) const final;
    virtual Vec3 f_uu(Vec2 const& p) const final;
    virtual Vec3 f_uv(Vec2 const& p) const final;
    virtual Vec3 f_vv(Vec2 const& p) const final;
};

class TrefoilKnot : public ParametricSurface {
public:
    TrefoilKnot();
    virtual Vec3 f(Vec2 const& p) const final;
};

#endif
