#pragma once

#include <string_view>

namespace Config {
extern std::string_view data_directory;
extern std::string_view out_directory;
extern std::string_view experiment;
extern std::string_view stem;

extern bool push_out;
extern bool pull_in;
extern bool alternate_push_pull;

extern bool use_bvh;
extern bool use_ray_shoot_mapping;
extern bool use_winding_order;
extern int sphere_tracing_iterations;
extern double friction_coefficient;
extern double outlier_threshold;
extern int pull_in_particles_per_path;

void argparse(int argc, char* argv[]);
};