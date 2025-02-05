#include "BVH.h"
#include "BVH/Morton.h"
#include <iostream>
#include <stack>
#include <unordered_map>
#include <unordered_set>

[[nodiscard]] inline auto closest_point(AABB const& volume, Vec3 const& p) -> Vec3
{
    return p.cwiseMax(volume.min()).cwiseMin(volume.max());
}

inline auto create_morton_curve(std::vector<Vec3> const* points) -> std::vector<std::pair<Morton::Code, size_t>>
{
    static constexpr auto MORTON_DOMAIN = 1048576ULL;

    auto curve = std::vector<std::pair<Morton::Code, size_t>>();
    auto body_volume = AABB();
    body_volume.setEmpty();

    curve.reserve(points->size());

    for (auto&& point : *points)
        body_volume.extend(point);

    Vec3 size = (body_volume.max() - body_volume.min()).array() + 1.;
    Vec3 scale = MORTON_DOMAIN * size.cwiseInverse();

    for (auto&& [i, point] : enumerate(*points)) {
        Eigen::Matrix<Morton::Code, 1, 3> const p = (point.array() * scale.array()).cast<Morton::Code>();
        auto const code = Morton::Encode(p.x(), p.y(), p.z());

        curve.push_back({ code, i });
    }

    std::sort(curve.begin(), curve.end(),
              [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });

    return curve;
}

BVH::BVH(std::vector<Vec3> const* points, std::vector<Vec3> const* normals, size_t num_points_per_leaf)
    : m_points(points)
    , m_normals(normals)
{
    auto const curve = create_morton_curve(m_points);
    auto leaves = std::deque<std::shared_ptr<BVHNode>>();

    for (auto i = 0uz; i < curve.size(); i += num_points_per_leaf) {
        AABB volume {};
        Vec3 normal = Vec3::Zero();
        auto theta = 0.;

        auto points = std::vector<size_t>();

        for (auto j = i; j < i + num_points_per_leaf && j < curve.size(); j++) {
            [[maybe_unused]] auto&& [code, idx] = curve[j];

            points.push_back(idx);
            volume.extend(m_points->at(idx));

            normal += m_normals->at(idx);
        }
        normal.normalize();

        for (auto&& idx : points) {
            auto const& n = m_normals->at(idx);
            theta = std::max(theta, std::acos(normal.dot(n)));
        }

        leaves.push_back(std::make_shared<BVHNode>(volume, nullptr, nullptr, std::move(points)));
    }

    while (leaves.size() > 1) {
        auto left = leaves.front();
        leaves.pop_front();

        auto right = leaves.front();
        leaves.pop_front();

        auto volume = AABB();
        volume.setEmpty();

        volume.extend(left->volume);
        volume.extend(right->volume);

        leaves.push_back(std::make_shared<BVHNode>(volume, left, right, std::nullopt));
    }

    m_root = leaves.front();
}

template <typename AABBMetric, typename PointMetric>
auto BVH::closest_point(
    AABBMetric&& aabb_metric,
    PointMetric&& point_metric,
    Vec3 const& x,
    Vec3 const& v) const -> size_t
{
    auto nodes = std::deque<std::shared_ptr<BVHNode>> { m_root };

    auto alpha = std::numeric_limits<double>::infinity();
    auto closest_point = 0uz;

    while (!nodes.empty()) {
        auto const& node = nodes.front();
        nodes.pop_front();

        if (node->data.has_value()) {
            auto const d_leaf = aabb_metric(node->volume, x, v);

            if (d_leaf > alpha)
                continue;

            for (auto&& [i, idx] : enumerate(node->data.value())) {
                auto const& p = m_points->at(idx);
                auto d_point = point_metric(p, x, v);

                if (d_point < alpha) {
                    alpha = d_point;
                    closest_point = idx;
                }
            }

            continue;
        }
        auto const d_left = aabb_metric(node->left->volume, x, v);
        auto const d_right = aabb_metric(node->right->volume, x, v);

        if (d_left < d_right) {
            nodes.push_back(node->left);

            if (d_right < alpha)
                nodes.push_back(node->right);
        } else {
            nodes.push_back(node->right);

            if (d_left < alpha)
                nodes.push_back(node->left);
        }
    }

    return closest_point;
}

[[gnu::hot]] auto BVH::closest_point(Vec3 const& x) const -> size_t
{
    return BVH::closest_point(
        [](AABB const& volume, Vec3 const& p, [[maybe_unused]] Vec3 const& v) {
            Vec3 x = ::closest_point(volume, p);
            return (x - p).squaredNorm();
        },
        [](Vec3 const& p, Vec3 const& x, [[maybe_unused]] Vec3 const& v) { return (p - x).squaredNorm(); },
        x);
}

auto ray_box_intersection(AABB const& aabb, Vec3 const& o, Vec3 const& d, double t = std::numeric_limits<double>::infinity()) -> bool
{
    Vec3 t1 = (aabb.min() - o).array() * d.cwiseInverse().array();
    Vec3 t2 = (aabb.max() - o).array() * d.cwiseInverse().array();

    double tmin = std::max(std::max(std::min(t1.x(), t2.x()), std::min(t1.y(), t2.y())), std::min(t1.z(), t2.z()));
    double tmax = std::min(std::min(std::max(t1.x(), t2.x()), std::max(t1.y(), t2.y())), std::max(t1.z(), t2.z()));

    return tmax >= std::max(0., tmin) && tmin < t;
}

[[gnu::hot]] auto BVH::closest_point_ray(Vec3 const& o, Vec3 const& d) const -> size_t
{
    auto nodes = std::deque<std::shared_ptr<BVHNode>> { m_root };

    auto alpha = std::numeric_limits<double>::infinity();
    auto beta = std::numeric_limits<double>::infinity();
    auto idx = 0uz;

    for (auto&& [i, point] : enumerate(*m_points)) {
        Vec3 w_proj = (point - o).dot(d) * d + o;
        double d_point = (w_proj - point).norm();

        if (d_point < beta) {
            beta = d_point;
            idx = i;
        }
    }

    return idx;

    // auto alpha = std::numeric_limits<double>::infinity();
    // auto beta = std::numeric_limits<double>::infinity();
    // auto closest_point = 0uz;

    // while (!nodes.empty()) {
    //     auto const& node = nodes.front();
    //     nodes.pop_front();

    //     // if (!ray_box_intersection(node->volume, o, d, alpha))
    //     //     continue;

    //     if (node->data.has_value()) {
    //         for (auto&& [i, idx] : enumerate(node->data.value())) {
    //             auto const& p = m_points->at(idx);
    //             double d_proj = (p - o).dot(d);
    //             double d_point = (p - o).norm();

    //             d_proj = (d_proj < 0.) ? std::numeric_limits<double>::infinity() : d_proj;

    //             if (d_proj < alpha && d_point < beta) {
    //                 alpha = d_proj;
    //                 beta = d_point;
    //                 closest_point = idx;
    //             }
    //         }

    //         continue;
    //     }

    //     nodes.push_back(node->left);
    //     nodes.push_back(node->right);
    // }

    // return closest_point;
}

[[gnu::hot]] auto BVH::closest_point_local(Vec3 const& x, Vec3 const& v) const -> size_t
{
    return BVH::closest_point(
        [](AABB const& volume, Vec3 const& p, Vec3 const& v) {
            Vec3 x = ::closest_point(volume, p);
            return (x - p).norm();
        },
        [](Vec3 const& p, Vec3 const& x, Vec3 const& v) {
            Vec3 const w = (p - x).dot(v) * v + x;
            return (w - p).norm();
        },
        x,
        v);
}