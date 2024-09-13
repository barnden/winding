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

    double m_uMax;
    double m_uMin;
    double m_vMax;
    double m_vMin;

    double m_epsilon;

    std::unique_ptr<float[]> get_display_data(int& len_out, int nc_res = 500, int nh_res = 500);
};

class H1 : public ParametricSurface {
public:
    H1()
    {
        m_uMin = 0.;
        m_uMax = 2. * PI;

        m_vMin = 0.;
        m_vMax = 1.;
    }

    virtual Vec3 f(Vec2 const& p) const
    {
        return {
            (p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
            (p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
            p.y()
        };
    }

    virtual Vec3 f_u(Vec2 const& p) const
    {
        return {
            -(p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
            (p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
            0.
        };
    }

    virtual Vec3 f_v(Vec2 const& p) const
    {
        return {
            (2. * p.y() - 1.) * cos(p.x()),
            (2. * p.y() - 1.) * sin(p.x()),
            1.
        };
    }

    virtual Vec3 f_uu(Vec2 const& p) const
    {
        return {
            -(p.y() * p.y() - p.y() + 0.4) * cos(p.x()),
            -(p.y() * p.y() - p.y() + 0.4) * sin(p.x()),
            0.
        };
    }

    virtual Vec3 f_uv(Vec2 const& p) const
    {
        return {
            -(2. * p.y() - 1.) * sin(p.x()),
            (2. * p.y() - 1.) * cos(p.x()),
            0.
        };
    }

    virtual Vec3 f_vv(Vec2 const& p) const
    {
        return {
            2. * cos(p.x()),
            2. * sin(p.x()),
            0.
        };
    }
};

class Torus : public ParametricSurface {
    double r1 = 1;
    double r2 = 0.3;

public:
    Torus()
    {
        m_uMin = 0;
        m_uMax = 2 * PI;
        m_vMin = 0;
        m_vMax = 2 * PI - 0.001;
    }
    virtual Vec3 f(Vec2 const& p) const
    {
        return { (r1 + r2 * cos(p.x())) * cos(p.y()), (r1 + r2 * cos(p.x())) * sin(p.y()), -r2 * sin(p.x()) };
    }

    virtual Vec3 f_u(Vec2 const& p) const
    {
        return { -r2 * sin(p.x()) * cos(p.y()), -r2 * sin(p.x()) * sin(p.y()), -r2 * cos(p.x()) };
    }

    virtual Vec3 f_v(Vec2 const& p) const
    {
        return { -(r1 + r2 * cos(p.x())) * sin(p.y()), (r1 + r2 * cos(p.x())) * cos(p.y()), 0 };
    }

    virtual Vec3 f_uu(Vec2 const& p) const
    {
        return { -r2 * cos(p.x()) * cos(p.y()), -r2 * cos(p.x()) * sin(p.y()), r2 * sin(p.x()) };
    }

    virtual Vec3 f_uv(Vec2 const& p) const
    {
        return { r2 * sin(p.x()) * sin(p.y()), -r2 * sin(p.x()) * cos(p.y()), 0 };
    }

    virtual Vec3 f_vv(Vec2 const& p) const
    {
        return { -(r1 + r2 * cos(p.x())) * cos(p.y()), -(r1 + r2 * cos(p.x())) * sin(p.y()), 0 };
    }
};

class Spring : public ParametricSurface {
    double m_r1;
    double m_r2;
    double m_kh;

public:
    Spring()
    {
        m_uMin = 0.;
        m_uMax = 2. * PI;
        m_vMin = 0.;
        m_vMax = 4. * PI;

        m_r1 = 1.;
        m_r2 = 0.3;
        m_kh = 0.15;
    }

    virtual Vec3 f(Vec2 const& p) const
    {
        return {
            (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
            (m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
            -m_r2 * sin(p.x()) + m_kh * p.y()
        };
    }

    virtual Vec3 f_u(Vec2 const& p) const
    {
        return {
            -m_r2 * sin(p.x()) * cos(p.y()),
            -m_r2 * sin(p.x()) * sin(p.y()),
            -m_r2 * cos(p.x())
        };
    }

    virtual Vec3 f_v(Vec2 const& p) const
    {
        return {
            -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
            (m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
            m_kh
        };
    }

    virtual Vec3 f_uu(Vec2 const& p) const
    {
        return {
            -m_r2 * cos(p.x()) * cos(p.y()),
            -m_r2 * cos(p.x()) * sin(p.y()),
            m_r2 * sin(p.x())
        };
    }

    virtual Vec3 f_uv(Vec2 const& p) const
    {
        return {
            m_r2 * sin(p.x()) * sin(p.y()),
            -m_r2 * sin(p.x()) * cos(p.y()),
            0
        };
    }

    virtual Vec3 f_vv(Vec2 const& p) const
    {
        return {
            -(m_r1 + m_r2 * cos(p.x())) * cos(p.y()),
            -(m_r1 + m_r2 * cos(p.x())) * sin(p.y()),
            0
        };
    }
};

class Vase : public ParametricSurface {
public:
    Vase()
    {
        m_uMin = 0.;
        m_uMax = 2. * PI;
        m_vMin = -1.;
        m_vMax = 1.;
    }

    virtual Vec3 f(Vec2 const& p) const
    {
        return {
            (0.5 - 0.4 * p.y() * (1. + p.y()) * (1. - p.y())) * cos(p.x()),
            (0.5 + 0.4 * p.y() * (1. + p.y()) * (1. - p.y())) * sin(p.x()),
            p.y()
        };
    }

    virtual Vec3 f_u(Vec2 const& p) const
    {
        return {
            -(0.5 - 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * sin(p.x()),
            (0.5 + 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * cos(p.x()),
            0
        };
    }

    virtual Vec3 f_v(Vec2 const& p) const
    {
        return {
            (-0.4 + 1.2 * p.y() * p.y()) * cos(p.x()),
            (0.4 - 1.2 * p.y() * p.y()) * sin(p.x()),
            1.
        };
    }

    virtual Vec3 f_uu(Vec2 const& p) const
    {
        return {
            -(0.5 - 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * cos(p.x()),
            -(0.5 + 0.4 * p.y() * (1 + p.y()) * (1 - p.y())) * sin(p.x()),
            0.
        };
    }

    virtual Vec3 f_uv(Vec2 const& p) const
    {
        return {
            -(-0.4 + 1.2 * p.y() * p.y()) * sin(p.x()),
            (0.4 - 1.2 * p.y() * p.y()) * cos(p.x()),
            0.
        };
    }

    virtual Vec3 f_vv(Vec2 const& p) const
    {
        return {
            2.4 * p.y() * cos(p.x()),
            -2.4 * p.y() * sin(p.x()),
            0.
        };
    }
};

#endif
