/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include <Eigen/Dense>
#include <vector>

#include "BVH/BVH.h"
#include "utils.h"

class ParametricSurface {
    mutable BVH m_bvh;
    void generate_search_grid(int nu, int nv) const;

protected:
    double m_uMin;
    double m_uMax;
    double m_vMin;
    double m_vMax;

    double m_epsilon;

    mutable std::vector<Vec2> m_grid2D;
    mutable std::vector<Vec3> m_grid3D;

public:
    ParametricSurface() = delete;
    ParametricSurface(double epsilon = 1e-5);
    ParametricSurface(double u_min, double u_max, double v_min, double v_max, double epsilon = 1e-5);

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

    [[nodiscard]] Vec2 closest_point(Vec3 const& p, size_t max_iterations = 1000) const;
    [[nodiscard]] Vec2 closest_point(Vec3 const& p, Vec2 const& guess, size_t max_iterations = 1000) const;

    double u_max() const { return m_uMax; }
    double u_min() const { return m_uMin; }
    double v_max() const { return m_vMax; }
    double v_min() const { return m_vMin; }
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
