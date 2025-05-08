#include <Eigen/Dense>

#include "Config.h"
#include "SurfaceEditor.h"

#ifdef STANDALONE
int main(int argc, char* argv[])
{
#    ifdef DEBUG_FPE_TRAP
    enable_floating_point_exceptions();
#    endif

    Config::argparse(argc, argv);

    auto num_revolutions = .8;
    auto num_paths = 32;
    auto num_particles = 200;

    auto timestep = 0.025;
    auto spring_constant = 1'000'000.;
    auto damping_coefficient = 200.;
    auto epsilon = 0.001;

    auto surface = CubicBSpline();
    auto editor = SurfaceEditor(surface, num_revolutions, num_paths, num_particles);

    for (auto i = 0; i < 32; i++) {
        editor.step(timestep, spring_constant, damping_coefficient, epsilon);
    }
}
#endif