#pragma once
#ifndef UTILS_H
#    define UTILS_H

#    include <Eigen/Dense>

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;

constexpr float PI = 3.14159265359;
#define enumerate(v) std::views::zip(std::views::iota(0), v)
#endif