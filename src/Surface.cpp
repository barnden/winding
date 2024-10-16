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

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
    // return ((f(p + Vec2(2.0 * m_epsilon, 0)) - f(p)) - (f(p) - f(p + Vec2(-2.0 * m_epsilon, 0)))) / 4.0 / m_epsilon / m_epsilon;
}

Vec3 ParametricSurface::nf_vv(Vec2 const& p) const
{
    Vec2 const delta(0., 2. * m_epsilon);

    return (f(p + delta) + f(p - delta) - 2. * f(p)) / (4. * m_epsilon * m_epsilon);
    // return ((f(p + Vec2(0, 2.0 * m_epsilon)) - f(p)) - (f(p) - f(p + Vec2(0, -2.0 * m_epsilon)))) / 4.0 / m_epsilon / m_epsilon;
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

Vec2 ParametricSurface::closest_point(Vec3 const& p, Vec2 const& guess, int max_iterations) const
{
    Vec2 xk = guess;

    for (int i = 0; i < max_iterations; i++) {
        Vec3 fp = f(xk) - p;
        Vec3 fu = f_u(xk);
        Vec3 fv = f_v(xk);
        Vec3 fuu = f_uu(xk);
        Vec3 fuv = f_uv(xk);
        Vec3 fvv = f_vv(xk);

        double h11 = fuu.dot(fp) + fu.dot(fu);
        double h12 = fuv.dot(fp) + fu.dot(fv);
        double h22 = fvv.dot(fp) + fv.dot(fv);

        Vec2 dx(h22 * fu.dot(fp) - h12 * fv.dot(fp), -h12 * fu.dot(fp) + h11 * fv.dot(fp));
        dx /= (h11 * h22 - h12 * h12);
        xk -= 0.5 * dx;

        if (dx.norm() <= 1e-5)
            break;
    }

    return xk;
}
