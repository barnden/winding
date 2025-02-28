/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Config.h"
#include <filesystem>
#include <format>
#include <iostream>

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

// TODO: Add argparser for this
int sphere_tracing_iterations = 5;
double friction_coefficient = .0;
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
        std::cout << "[Config] 'push_out' and 'pull_in' specified without 'alternate_push_pull', defaulting to 'push_out'\n";
        pull_in = false;
    }

    if (!push_out && !pull_in) {
        std::cout << "[Config] Neither 'push_out' nor 'pull_in' were specified, defaulting to 'push_out'\n";
        push_out = true;
    }

    if (!fs::exists(data_directory))
        throw std::format("[Config] Data directory \"{}\" does not exist\n", data_directory);

    fs::create_directories(std::format("{}/{}", out_directory, experiment));

    for (auto&& directory : { "spline"sv, "path"sv, "max_quad"sv, "mesh"sv }) {
        auto path = std::format("{}/{}/{}", out_directory, experiment, directory);
        fs::create_directory(path);
    }
}