#include "Simulator.h"
#include "Surface.h"
#include "utils.h"

#include <Eigen/Dense>
#include <iostream>

#define SF_OFF 0.001

struct OrderNode {
    int path_id;
    OrderNode* pnext;
    OrderNode* ppre;
    double u;
};

void generate_winding_order(std::vector<int>& winding_order, std::vector<std::vector<Vec2>> const& init_sv)
{
    OrderNode* nodes = new OrderNode[init_sv.size()];
    OrderNode* z_pos = nodes;
    OrderNode* z_neg = nodes + 1;
    OrderNode* min_p = nodes;

    z_pos->ppre = z_pos;
    z_neg->ppre = z_neg;
    z_pos->pnext = z_pos;
    z_neg->pnext = z_neg;

    nodes[0].path_id = 0;
    nodes[1].path_id = 1;

    nodes[0].u = fmod(init_sv[0][0].x(), 2. * PI);

    if (nodes[0].u < 0.)
        nodes[0].u += 2. * PI;

    nodes[1].u = init_sv[1][0].x();

    for (auto i = 2uz; i < init_sv.size(); i++) {
        nodes[i].path_id = i;
        nodes[i].u = fmod(init_sv[i][0].x(), 2. * PI);

        if (nodes[i].u < 0.)
            nodes[i].u += 2. * PI;

        if (i % 2) {
            nodes[i].pnext = z_neg->pnext;
            nodes[1].ppre = z_neg;
            z_neg->pnext = nodes + i;
            nodes[i].pnext->ppre = nodes + i;

            continue;
        }

        nodes[i].pnext = z_pos->pnext;
        nodes[i].ppre = z_pos;
        z_pos->pnext = nodes + i;
        nodes[i].pnext->ppre = nodes + i;

        if (min_p->u > nodes[i].u)
            min_p = nodes + i;
    }

    z_pos = min_p;
    auto nextPos = true;
    auto* pnext = nodes;

    while (pnext) {
        winding_order.push_back(pnext->path_id);
        auto nextu = fmod(init_sv[pnext->path_id].back().x() - 0.1, 2. * PI);

        if (nextu < 0)
            nextu += 2. * PI;

        if (nextPos) {
            nextPos = false;

            if (z_pos->pnext == z_pos) {
                pnext = z_neg;

                continue;
            }

            pnext->ppre->pnext = pnext->pnext;
            pnext->pnext->ppre = pnext->ppre;

            if (pnext == z_pos)
                z_pos = z_pos->ppre;

            if (nextu < z_neg->u || nextu >= z_neg->pnext->u) {
                pnext = z_neg->pnext;

                continue;
            }

            pnext = z_neg->pnext->pnext;

            while (nextu < pnext->u)
                pnext = pnext->pnext;

            continue;
        }

        nextPos = true;

        if (z_neg->pnext == z_neg)
            break;

        pnext->ppre->pnext = pnext->pnext;
        pnext->pnext->ppre = pnext->ppre;

        if (pnext == z_neg)
            z_neg = z_neg->ppre;

        if (nextu < z_pos->u || nextu >= z_pos->pnext->u) {
            pnext = z_pos->pnext;

            continue;
        }

        pnext = z_pos->pnext->pnext;

        while (nextu < pnext->u)
            pnext = pnext->pnext;
    }

    delete[] nodes;
}

void init(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles, std::vector<std::vector<Vec2>>& init_sv)
{
    init_sv = std::vector<std::vector<Vec2>>(2 * num_paths, std::vector<Vec2>(num_particles));

    auto angle = (surface.m_vMax - surface.m_vMin) / (surface.m_uMax - surface.m_uMin) / num_revolutions;
    auto du = (surface.m_uMax - surface.m_uMin) / num_paths;
    auto ddv = (surface.m_vMax - surface.m_vMin) / (num_particles + 1);
    auto ddu = ddv / angle;

    for (auto i = 0; i < num_paths; i++) {
        Vec2 p(surface.m_uMin + i * du, surface.m_vMin);
        Vec2 q(surface.m_uMin + i * du, surface.m_vMin);

        Vec3 cur = surface.f(p) + SF_OFF * surface.normal(p);

        for (auto j = 0; j < num_particles; j++) {
            cur = surface.f(p) + SF_OFF * surface.normal(p);
            init_sv[2 * i][num_particles - 1 - j] = p;
            init_sv[2 * i + 1][j] = q;

            p += Vec2(ddu, ddv);
            q += Vec2(-ddu, ddv);
        }
    }
}

int main()
{
    auto num_revolutions = 1.;
    auto num_paths = 1;
    auto num_particles = 201;

    auto dt = 0.05;
    auto ksp = 1'000'000.;
    auto kdp = 200.;
    auto eps = 0.001;

    auto surface = H1();
    auto winding_order = std::vector<int> {};
    auto init_sv = std::vector<std::vector<Vec2>> {};

    init(surface, num_revolutions, num_paths, num_particles, init_sv);
    generate_winding_order(winding_order, init_sv);

    auto parr = std::vector<std::vector<Vec2>>(init_sv.size());
    auto rarr = std::vector<std::vector<int>>(init_sv.size());

    for (auto&& [j, i] : enumerate(winding_order)) {
        auto simulator = OffSurface(surface, init_sv[i]);

        simulator.ksp() = ksp;
        simulator.kdp() = kdp;
        simulator.dt() = dt;
        simulator.eps() = eps;

        simulator.simulate(1000);
        simulator.mapping();

        parr[i] = simulator.p();
        rarr[i] = simulator.r();

        auto progress = 100. * (j + 1.) / winding_order.size();
        std::cout << round(progress) << "%\n\n";
    }

    // for (auto&& [i, p, r] : std::views::zip(std::views::iota(0), parr[1], rarr[1])) {
    //     std::cout << i << ": (" << r << ") " << p.x() << ", " << p.y() << '\n';
    // }
}
