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

    virtual Vec3 f(Vec2 const& p) const = 0;
    virtual Vec3 f_u(Vec2 const& p) const;
    virtual Vec3 f_v(Vec2 const& p) const;
    virtual Vec3 f_uu(Vec2 const& p) const;
    virtual Vec3 f_uv(Vec2 const& p) const;
    virtual Vec3 f_vv(Vec2 const& p) const;
    virtual Vec3 normal(Vec2 const& p) const;

    Vec3 nf_u(Vec2 const& p) const;
    Vec3 nf_v(Vec2 const& p) const;
    Vec3 nf_uu(Vec2 const& p) const;
    Vec3 nf_uv(Vec2 const& p) const;
    Vec3 nf_vv(Vec2 const& p) const;

    Vec2 rescale(Vec2 const& p) const;

    void generate_search_grid(int nu, int nv);

    Vec2 closest_point(Vec3 const& p) const;
    Vec2 closest_point(Vec3 const& p, Vec2 const& guess) const;

    std::vector<double> get_jn_ref(Vec2 const& p, double l) const;

    std::vector<Vec2> m_grid2D;
    std::vector<Vec3> m_grid3D;

    double m_uMax, m_uMin, m_vMax, m_vMin;
    double m_epsilon;

    std::unique_ptr<float[]> get_display_data(int& len_out, int nc_res = 500, int nh_res = 500);
};

#endif