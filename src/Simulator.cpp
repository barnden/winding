#include "Simulator.h"

#include <ranges>
#include <iostream>

#include <Eigen/Sparse>
using MatrixSd = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;
using Triplet = Eigen::Triplet<double>;

Simulator::Simulator(
    ParametricSurface const& f,
    std::vector<Vec2> const& init_path)
    : m_size(201)
    , m_surface(f)
    , m_ksp(1'000'000)
    , m_kdp(200)
    , m_kpf(10'000)
    , m_kbs(0.001)
    , m_dt(0.001)
    , m_t(0.)
    , m_v_eps(0.001)
{
    if (init_path.size()) {
        m_size = init_path.size();
        m_p_init = init_path;

        return;
    }

    m_p_init = decltype(m_p_init)(m_size, Vec2::Zero());

    Vec2 p0(0., -1.);
    Vec2 p1(9.43347686131086, 0.9965750445589878);
    Vec2 dp = (p1 - p0) / (m_size - 1.);

    for (int i = 0; i < m_size; i++)
        m_p_init[i] = p0 + (double)i * dp;
}

Simulator::~Simulator() {};

void Simulator::simulate(int k)
{
    int i = 0;
    for (; (i < k && i < 10) || (i < k && !stop()); i++) {
        m_t += m_dt;
        step();
    }

    std::cout << i << " iterations simulated.\n";
}

OffSurface::OffSurface(ParametricSurface const& f, std::vector<Vec2> const& init_path)
    : Simulator(f, init_path)
{
    m_p = m_p_init;
    m_v = decltype(m_v)(size(), Vec2::Zero());
    m_touching = decltype(m_touching)(size(), true);
    m_l = decltype(m_l)(size(), 1);
    m_r = decltype(m_r)(size(), 1);
}

void OffSurface::step()
{
    std::vector<int> touching;

    for (int i = 0; i < size(); i++) {
        if (m_touching[i])
            touching.push_back(i);
    }

    const int N = touching.size();

    MatrixSd JT(2 * N, 3 * N);
    MatrixSd dJ(3 * N, 2 * N);
    MatrixSd Mh2K(3 * N, 3 * N);
    MatrixSd mat;

    Vector farr(3 * N);
    Vector varr(2 * N);
    Vector barr;

    std::vector<Triplet> JT_ele(6 * N, Triplet());
    std::vector<Triplet> dJ_ele(6 * N, Triplet());
    std::vector<Triplet> Mh2K_ele(9 * N - 12);

    for (auto&& [i, pidx] : enumerate(touching)) {
        varr(2 * i) = m_v[pidx].x();
        varr(2 * i + 1) = m_v[pidx].y();

        Vec3 fu = surface().f_u(m_p[pidx]);
        Vec3 fv = surface().f_v(m_p[pidx]);
        Vec3 fuu = surface().f_uu(m_p[pidx]);
        Vec3 fuv = surface().f_uv(m_p[pidx]);
    }
}

bool OffSurface::stop()
{
    // for (int i = 0; i < size(); i++)
    //     if (m_touching[i] && )

    return true;
}