#include <fstream>

#include "BSpline.h"
#include "utils.h"
#include <iostream>

[[gnu::always_inline]] std::tuple<double, double, int, int> CubicBSpline::get_uv(Vec2 const& p) const
{
    Vec2 uv = p;

    uv.x() = std::fmod(uv.x(), m_u_max);

    if (uv.x() < 0.)
        uv.x() = std::fmod(uv.x() + m_u_max, m_u_max);

    // Get the Bezier patch
    auto iu = (int)uv.x();
    auto iv = (int)uv.y();

    iv = std::clamp(iv, 2, m_nv - 2);

    // Get the local u/v coordinates on the Bezier patch
    uv.x() -= iu;
    uv.y() -= iv;

    iu %= m_nu;

    return std::make_tuple(uv.x(), uv.y(), iu, iv);
}

[[gnu::always_inline]] Eigen::RowVector3d CubicBSpline::get(int v, int u) const
{
    u %= m_nu;

    if (u < 0)
        u += m_nu;

    return m_points[v][u];
}

template <typename UFunc, typename FFunc>
[[gnu::flatten]] Vec3 CubicBSpline::interpolate(UFunc&& get_u, FFunc&& get_v, Vec2 const& p) const
{
    using Mat4x3 = Eigen::Matrix<double, 4, 3>;
    auto [u, v, iu, iv] = get_uv(p);

    Eigen::RowVector4d U = get_u(u) * basis();
    Eigen::RowVector4d V = get_v(v) * basis();

    // clang-format off
    Mat4x3 P = (Mat4x3{} <<
                V * (Mat4x3{} << get(iv - 2, iu - 2), get(iv - 1, iu - 2), get(iv, iu - 2), get(iv + 1, iu - 2)).finished(),
                V * (Mat4x3{} << get(iv - 2, iu - 1), get(iv - 1, iu - 1), get(iv, iu - 1), get(iv + 1, iu - 1)).finished(),
                V * (Mat4x3{} << get(iv - 2, iu + 0), get(iv - 1, iu + 0), get(iv, iu + 0), get(iv + 1, iu + 0)).finished(),
                V * (Mat4x3{} << get(iv - 2, iu + 1), get(iv - 1, iu + 1), get(iv, iu + 1), get(iv + 1, iu + 1)).finished())
            .finished();
    // clang-format on

    return U * P;
}

void CubicBSpline::read(std::string const& file)
{
    auto stream = std::ifstream(file);

    if (stream.fail()) {
        throw "Failed to open file " + file + '\n';
    }

    stream >> m_nu >> m_nv;

    m_points = std::vector(m_nv, std::vector<Vec3>(m_nu, Vec3::Zero()));

    m_u_min = 0.;
    m_u_max = m_nu;
    m_v_min = 2.;
    m_v_max = m_nv - 1.;

    for (auto i = 0; i < m_nv; i++) {
        for (auto j = 0; j < m_nu; j++) {
            stream >> m_points[i][j].x() >> m_points[i][j].y() >> m_points[i][j].z();
        }
    }

    stream.close();
}

[[gnu::flatten, gnu::hot]] Vec3 CubicBSpline::f(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[gnu::flatten]] Vec3 CubicBSpline::f_u(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[gnu::flatten]] Vec3 CubicBSpline::f_v(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
        p);
}

[[gnu::flatten]] Vec3 CubicBSpline::f_uv(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(3. * u * u, 2. * u, 1., 0.); },
        [](double v) { return Eigen::RowVector4d(3. * v * v, 2. * v, 1., 0.); },
        p);
}

[[gnu::flatten]] Vec3 CubicBSpline::f_uu(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(6. * u, 2., 0., 0.); },
        [](double v) { return Eigen::RowVector4d(v * v * v, v * v, v, 1.); },
        p);
}

[[gnu::flatten]] Vec3 CubicBSpline::f_vv(Vec2 const& p) const
{
    return interpolate(
        [](double u) { return Eigen::RowVector4d(u * u * u, u * u, u, 1.); },
        [](double v) { return Eigen::RowVector4d(6. * v, 2., 0., 0.); },
        p);
}

[[gnu::flatten]] Eigen::MatrixXd CubicBSpline::jacobian(Vec2 const& p) const
{
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, 3 * m_nv * m_nu);

    auto const get_index = [&](int v, int u) {
        u %= m_nu;

        if (u < 0)
            u += m_nu;

        return 3 * (v * m_nu + u);
    };
    auto const [u, v, iu, iv] = get_uv(p);

    Eigen::RowVector4d U = Eigen::RowVector4d(u * u * u, u * u, u, 1.) * basis();
    Eigen::RowVector4d V = Eigen::RowVector4d(v * v * v, v * v, v, 1.) * basis();

    for (auto l = 0; l < 4; l++)
        for (auto m = 0; m < 4; m++)
            J.block<3, 3>(0, get_index(iv + (l - 2), iu + (m - 2))).diagonal().array() = V(l) * U(m);

    return J;
}