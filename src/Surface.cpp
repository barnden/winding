#include "Surface.h"
#include <ranges>

ParametricSurface::ParametricSurface()
    : m_epsilon(1e-5) {};

Vec3 ParametricSurface::f_u(Vec2 const& p) const { return nf_u(p); }

Vec3 ParametricSurface::f_v(Vec2 const& p) const { return nf_v(p); }

Vec3 ParametricSurface::f_uu(Vec2 const& p) const { return nf_uu(p); }

Vec3 ParametricSurface::f_uv(Vec2 const& p) const { return nf_uv(p); }

Vec3 ParametricSurface::f_vv(Vec2 const& p) const { return nf_vv(p); }

Vec3 ParametricSurface::normal(Vec2 const& p) const
{
    return f_u(p).cross(f_v(p)).normalized();
}

Vec3 ParametricSurface::nf_u(Vec2 const& p) const
{
    Vec2 const delta(m_epsilon, 0.);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_v(Vec2 const& p) const
{
    Vec2 const delta(0., m_epsilon);

    return (f(p + delta) - f(p - delta)) / (2. * m_epsilon);
}

Vec3 ParametricSurface::nf_uu(Vec2 const& p) const
{
    Vec2 const delta(2. * m_epsilon, 0.);

    return (f(p + delta) - f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
}

Vec3 ParametricSurface::nf_vv(Vec2 const& p) const
{
    Vec2 const delta(0., 2. * m_epsilon);

    return (f(p + delta) - f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
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

void ParametricSurface::generate_search_grid(int nu, int nv)
{
    m_grid2D = decltype(m_grid2D)((nu + 1) * (nv + 1), Vec2::Zero());
    m_grid3D = decltype(m_grid3D)((nu + 1) * (nv + 1), Vec3::Zero());

    double du = (m_uMax - m_uMin) / nu;
    double dv = (m_vMax - m_vMin) / nv;

    int ct = 0;
    for (int i = 0; i <= nu; i++) {
        for (int j = 0; j <= nv; j++) {
            Vec2 p(m_uMin + du * i, m_vMin + dv * j);

            m_grid2D[ct] = p;
            m_grid3D[ct] = f(p);
        }
    }
}

Vec2 ParametricSurface::closest_point(Vec3 const& p) const
{
    double cur;
    double mi = (p - m_grid3D[0]).norm();
    int id = 0;

    for (auto&& [i, q] : enumerate(m_grid3D)) {
        cur = (p - q).norm();

        if (cur < mi) {
            mi = cur;
            id = i;
        }
    }

    return closest_point(p, m_grid2D[id]);
}

Vec2 ParametricSurface::closest_point(Vec3 const& p, Vec2 const& guess) const
{
    Vec2 xk = guess;
    Vec2 dx = Vec2(1., 1.);

    double h11;
    double h12;
    double h22;

    Vec3 f_ = Vec3::Zero();
    Vec3 f_u_ = Vec3::Zero();
    Vec3 f_v_ = Vec3::Zero();
    Vec3 f_uu_ = Vec3::Zero();
    Vec3 f_uv_ = Vec3::Zero();
    Vec3 f_vv_ = Vec3::Zero();

    for (int i = 0; (i < 1000) && (dx.norm() > 1e-5); i++) {
        f_ = f(xk);
        f_u_ = f_u(xk);
        f_v_ = f_v(xk);
        f_uu_ = f_uu(xk);
        f_uv_ = f_uv(xk);
        f_vv_ = f_vv(xk);

        h11 = f_uu_.dot(f_ - p) + f_u_.dot(f_u_);
        h12 = f_uv_.dot(f_ - p) + f_u_.dot(f_v_);
        h22 = f_vv_.dot(f_ - p) + f_v_.dot(f_v_);

        dx.x() = h22 * f_u_.dot(f_ - p) - h12 * f_v_.dot(f_ - p);
        dx.y() = -h12 * f_u_.dot(f_ - p) - h11 * f_v_.dot(f_ - p);

        dx /= (h11 * h22 - h12 * h12);

        xk -= 0.5 * dx;
    }

    return xk;
}

std::vector<double> ParametricSurface::get_jn_ref(Vec2 const& p, double l) const
{
    Vec3 fn0 = normal(p);
    Vec3 fnx = normal(closest_point(f(p) + Vec3(m_epsilon, 0., 0.), p));
    Vec3 fny = normal(closest_point(f(p) + Vec3(0., m_epsilon, 0.), p));
    Vec3 fnz = normal(closest_point(f(p) + Vec3(0., 0., m_epsilon), p));

    fnx = (fnx - fn0) / m_epsilon;
    fny = (fny - fn0) / m_epsilon;
    fnz = (fnz - fn0) / m_epsilon;

    return {
        fnx.x(), fny.x(), fnz.x(),
        fnx.y(), fny.y(), fnz.y(),
        fnx.z(), fny.z(), fnz.z()
    };
}

// std::unique_ptr<float[]> ParametricSurface::get_display_data(int& len_out, int nc_res, int nh_res)
// {
//     len_out = nc_res * nh_res * 4 * 6;

//     auto const result = std::make_unique<float[]>(len_out);

//     // TBD
// }