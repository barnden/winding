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
    std::vector<Vec3> const* m_normals;

public:
    BVH()
        : m_root(nullptr)
        , m_points(nullptr)
        , m_normals(nullptr) { };
    BVH(std::vector<Vec3> const* points, std::vector<Vec3> const* normals, size_t num_points_per_leaf = 8);

    template<typename AABBMetric, typename PointMetric>
    auto closest_point(AABBMetric&& aabb_metric, PointMetric&& point_metric, Vec3 const& x, Vec3 const& v=Vec3::Zero()) const -> size_t;
    auto closest_point_ray(Vec3 const& o, Vec3 const& d) const -> size_t;

    auto closest_point(Vec3 const& x) const -> size_t;
    auto closest_point_local(Vec3 const& x, Vec3 const& v) const -> size_t;
};