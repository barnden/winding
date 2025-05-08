#include <filesystem>
#include <format>
#include <print>

#include "Config.h"

namespace fs = std::filesystem;
using namespace std::literals;
namespace Config {
std::string_view data_directory = "reference"sv;
std::string_view out_directory = "experiments"sv;
std::string_view experiment = "experiment"sv;
std::string_view stem = "vase"sv;

bool push_out = true;
bool pull_in = false;
bool alternate_push_pull = false;

bool use_bvh = true;
bool use_ray_shoot_mapping = true;

// TODO: Add argparser for the following options
int sphere_tracing_iterations = 5;
double friction_coefficient = 0.;

// NOTE: The original code had a way of generating the order for winding paths
// (see Composer.cpp), however, I have no idea how it works. It does _not_ affect
// anything else other than the ordering of the paths.
bool use_winding_order = false;

// NOTE: As stated in the paper, outlier threshold of 4 times IQR was a reasonable
// guess for most shapes. This number should be tuned for each surface.
double outlier_threshold = 4.;

// NOTE: This is the number of particles we create along each edge of the midpoint
// quadrilateral in each quad. Ideally this should be adaptive to the side length.
int pull_in_particles_per_path = 35;
}

void Config::argparse(int argc, char* argv[])
{
    if (argc != 6)
        throw std::format("[Config] Missing positional arguments, expected 6 got {}\n", argc);

    data_directory = std::string_view(argv[1]);
    stem = std::string_view(argv[2]);
    out_directory = std::string_view(argv[3]);
    experiment = std::string_view(argv[4]);

    auto flags = std::string_view(argv[5]);

    push_out = flags.contains('o');
    pull_in = flags.contains('i');
    alternate_push_pull = flags.contains('a');
    use_bvh = !flags.contains('B');
    use_ray_shoot_mapping = !flags.contains('R');

    if (push_out && pull_in && !alternate_push_pull) {
        std::println("[Config] 'push_out' and 'pull_in' specified without 'alternate_push_pull', defaulting to 'push_out'");
        pull_in = false;
    }

    if (!push_out && !pull_in) {
        std::println("[Config] Neither 'push_out' nor 'pull_in' were specified, defaulting to 'push_out'");
        push_out = true;
    }

    if (!fs::exists(data_directory))
        throw std::format("[Config] Data directory \"{}\" does not exist\n", data_directory);

    fs::create_directories(std::format("{}/{}", out_directory, experiment));

    for (auto&& directory : { "spline"sv, "path"sv, "mesh"sv }) {
        auto path = std::format("{}/{}/{}", out_directory, experiment, directory);
        fs::create_directory(path);
    }
}