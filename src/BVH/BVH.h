#pragma once

#include "utils.h"

#include <Eigen/Dense>
#include <memory>
#include <optional>

using AABB = Eigen::AlignedBox3d;

struct BVHNode {
    AABB volume;

    std::shared_ptr<BVHNode> left;
    std::shared_ptr<BVHNode> right;

    std::optional<std::vector<size_t>> data = std::nullopt;
};

class BVH {
    std::shared_ptr<BVHNode> m_root;
    std::vector<Vec3> const* m_points;

public:
    BVH()
        : m_root(nullptr)
        , m_points(nullptr) {};
    BVH(std::vector<Vec3> const* points, size_t num_points_per_leaf = 8);

    size_t closest_point(Vec3 const& x) const;
};