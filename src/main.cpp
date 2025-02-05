#include "Composer.h"
#include "Simulator.h"
#include "SurfaceEditor.h"
#include "Surfaces/Surface.h"
#include "utils.h"
#include "Config.h"

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>


#ifdef STANDALONE
int main(int argc, char* argv[])
{
#    ifdef DEBUG_FPE_TRAP
    enable_floating_point_exceptions();
#    endif

    Config::argparse(argc, argv);

    auto num_revolutions = .6;
    auto num_paths = 32;
    auto num_particles = 200;

    auto dt = 0.05;
    auto ksp = 1'000'000.;
    auto kdp = 200.;
    auto eps = 0.001;

    auto surface = CubicBSpline();
    auto editor = SurfaceEditor(surface, num_revolutions, num_paths, num_particles);

    for (auto i = 0; i < 16; i++) {
        editor.step(dt, ksp, kdp, eps);
    }
}
#endif