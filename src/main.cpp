#include "Composer.h"
#include "Simulator.h"
#include "SurfaceEditor.h"
#include "Surfaces/Surface.h"
#include "utils.h"

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>

#include <filesystem>

namespace fs = std::filesystem;

#ifdef STANDALONE
int main(int argc, char* argv[])
{
#    ifdef DEBUG_FPE_TRAP
    enable_floating_point_exceptions();
#    endif

    auto options = std::make_shared<Options>();
    if (argc > 1)
        options->data_path = std::string(argv[1]);

    if (argc > 2)
        options->file_stem = std::string(argv[2]);

    if (argc > 3)
        options->out_path = std::string(argv[3]);

    if (argc > 3)
        options->experiment = std::string(argv[4]);

    if (argc > 4) {
        std::string flags { argv[5] };

        if (flags.contains('o')) {
            options->push_out = true;
            std::cout << "options->push_out = true\n";
        }

        if (flags.contains('i')) {
            options->pull_in = true;
            std::cout << "options->pull_in = true\n";
        }

        if (flags.contains('b')) {
            options->use_bvh = true;
            std::cout << "options->use_bvh = true\n";
        }

        if (!options->push_out && !options->pull_in) {
            std::cout << "no flag for push-out or pull-in; defaulting to push-out\n";
            options->push_out = true;
        }
    }

    {
        auto create_if_not_exists = [](auto const& path) {
            if (!fs::is_directory(path) || fs::exists(path))
                fs::create_directory(path);
        };

        create_if_not_exists(options->out_path);
        create_if_not_exists(options->out_path + "/" + options->experiment);
        create_if_not_exists(options->out_path + "/" + options->experiment + "/spline");
        create_if_not_exists(options->out_path + "/" + options->experiment + "/path");
        create_if_not_exists(options->out_path + "/" + options->experiment + "/max_quad");

        fs::copy_file(options->data_path + "/" + options->file_stem + ".txt", options->out_path + "/" + options->experiment + "/spline/step-0.txt", fs::copy_options::overwrite_existing);
    }

    auto num_revolutions = .6;
    auto num_paths = 32;
    auto num_particles = 200;

    auto dt = 0.05;
    auto ksp = 1'000'000.;
    auto kdp = 200.;
    auto eps = 0.001;

    auto surface = CubicBSpline(options);
    auto editor = SurfaceEditor(surface, options, num_revolutions, num_paths, num_particles);

    for (auto i = 0; i < 16; i++) {
        editor.step(dt, ksp, kdp, eps);
    }
}
#endif