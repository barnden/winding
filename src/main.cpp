#include "Composer.h"
#include "Simulator.h"
#include "Surface.h"
#include "utils.h"

#include <Eigen/Dense>
#include <iostream>

int main()
{
#ifdef DEBUG_FPE_TRAP
    enable_floating_point_exceptions();
#endif
    auto num_revolutions = 10.;
    auto num_paths = 1;
    auto num_particles = 200;

    auto dt = 0.05;
    auto ksp = 1'000'000.;
    auto kdp = 200.;
    auto eps = 0.001;

    auto surface = TrefoilKnot();
    auto composer = Composer(surface, num_revolutions, num_paths, num_particles);
    composer.generate_winding_order();

    auto parr = std::vector<std::vector<Vec2>>(composer.initial_path().size());
    auto rarr = std::vector<std::vector<int>>(composer.initial_path().size());

    for (auto&& [j, i] : enumerate(composer.winding_order())) {
        auto simulator = OffSurface(surface, composer.initial_path()[i]);

        simulator.ksp() = ksp;
        simulator.kdp() = kdp;
        simulator.dt() = dt;
        simulator.eps() = eps;

        simulator.simulate(1000);
        simulator.mapping();

        parr[i] = simulator.p();
        rarr[i] = simulator.r();

        auto progress = 100. * (j + 1.) / composer.winding_order().size();
        std::cout << round(progress) << "%\n\n";
    }

    // for (auto&& [j, i] : enumerate(composer.winding_order())) {
    //     auto &p = parr[i];
    //     auto &r = rarr[i];
    // }

    for (auto&& [i, p, r] : std::views::zip(std::views::iota(0), parr[2], rarr[2])) {
        std::cout << i << ": (" << r << ") " << p.x() << ", " << p.y() << '\n';
    }
}
