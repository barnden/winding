#include "Composer.h"
#include "Simulator.h"
#include "SurfaceEditor.h"
#include "Surfaces/Surface.h"
#include "utils.h"

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <string_view>

#include <filesystem>

namespace fs = std::filesystem;
using namespace std::string_view_literals;

#ifdef STANDALONE
int main(int argc, char* argv[])
{
#    ifdef DEBUG_FPE_TRAP
    enable_floating_point_exceptions();
#    endif

    Options options;
    if (argc > 1)
        options.data_path = std::string(argv[1]);

    if (argc > 2)
        options.file_stem = std::string(argv[2]);

    if (argc > 3)
        options.out_path = std::string(argv[3]);

    if (argc > 3)
        options.experiment = std::string(argv[4]);

    {
        auto create_directory = [](auto const& path) {
            if (!fs::is_directory(path) || !fs::exists(path))
                return fs::create_directory(path);

            return true;
        };

        for (std::string directory : { "", "spline", "path", "max_quad", "obj" }) {
            auto path = options->out_path + "/" + options->experiment + "/" + directory;

            if (!create_directory(path))
                throw std::runtime_error("fail create_directory(\"" + path + "\")");
        }

        fs::copy_file(options.data_path + "/" + options.file_stem + ".txt", options.out_path + "/" + options.experiment + "/spline/step-0.txt", fs::copy_options::overwrite_existing);
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