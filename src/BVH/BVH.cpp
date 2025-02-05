/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "BVH.h"
#include "BVH/Morton.h"
#include <iostream>
#include <stack>
#include <unordered_map>
#include <unordered_set>

inline decltype(auto) create_morton_curve(std::vector<Vec3> const* points)
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

BVH::BVH(std::vector<Vec3> const* points, size_t num_points_per_leaf)
    : m_points(points)
{
    auto const curve = create_morton_curve(m_points);
    auto leaves = std::deque<std::shared_ptr<BVHNode>>();

    for (auto i = 0uz; i < curve.size(); i += num_points_per_leaf) {
        AABB volume {};
        auto points = std::vector<size_t>();

        for (auto j = i; j < i + num_points_per_leaf && j < curve.size(); j++) {
            [[maybe_unused]] auto&& [code, idx] = curve[j];

            points.push_back(idx);
            volume.extend(m_points->at(idx));
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

double distance(AABB const& volume, Vec3 const& p)
{
    Vec3 x = p.cwiseMax(volume.min()).cwiseMin(volume.max());
    return (x - p).squaredNorm();
}

[[gnu::hot]] size_t BVH::closest_point(Vec3 const& x) const
{
    auto nodes = std::deque<std::shared_ptr<BVHNode>> { m_root };

    auto alpha = std::numeric_limits<double>::infinity();
    auto closest_point = 0uz;

    while (!nodes.empty()) {
        auto const& node = nodes.front();
        nodes.pop_front();

        if (node->data.has_value()) {
            auto const d_leaf = distance(node->volume, x);

            if (d_leaf > alpha)
                continue;

            for (auto&& [i, idx] : enumerate(node->data.value())) {
                auto const& p = m_points->at(idx);
                auto const d_point = (p - x).squaredNorm();

                if (d_point < alpha) {
                    alpha = d_point;
                    closest_point = idx;
                }
            }

            continue;
        }
        auto const d_left = distance(node->left->volume, x);
        auto const d_right = distance(node->right->volume, x);

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