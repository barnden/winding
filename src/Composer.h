/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once
#ifndef COMPOSER_H
#    define COMPOSER_H

#    include "Surfaces/Surface.h"
#    include "utils.h"

#    include <print>
#    include <vector>

#    define SF_OFF            0.001
#    define ANGLE_SPEED_CONST 2.0

static constexpr double Infinity = std::numeric_limits<double>::infinity();

struct LocalFrame {
    Vec3 position;
    Vec3 normal;
    Vec3 direction;
    Vec2 shadow;

    double distance;
    double path_distance;
    double angle;
};

struct IntersectionNode {
    Vec2 parametric = Vec2::Zero();
    Vec3 point = Vec3::Zero();

    double angle_score = 0.;
    double dist_score = 0.;
    double t_up = 0.;
    double t_down = 0.;

    IntersectionNode* prev_up = nullptr;
    IntersectionNode* next_up = nullptr;

    IntersectionNode* prev_down = nullptr;
    IntersectionNode* next_down = nullptr;

    bool is_end = false;
};

class IntersectionList {
public:
    IntersectionNode *head, *rear;

    IntersectionList();

    void insert(IntersectionNode* node, IntersectionNode* IntersectionNode::* prev, IntersectionNode* IntersectionNode::* next, double IntersectionNode::* value) const;
    void interpolate(double t, double& dist, double& angle, IntersectionNode* IntersectionNode::* next, double IntersectionNode::* value) const;
};

class IntersectionListUp : public IntersectionList {
public:
    IntersectionListUp()
        : IntersectionList()
    {
        head->t_up = -Infinity;
        rear->t_up = Infinity;

        head->next_up = rear;
        rear->prev_up = head;
    }

    void insert(IntersectionNode* node) const
    {
        IntersectionList::insert(node, &IntersectionNode::prev_up, &IntersectionNode::next_up, &IntersectionNode::t_up);
    }

    void interpolate(double t, double& dist, double& angle) const
    {
        IntersectionList::interpolate(t, dist, angle, &IntersectionNode::next_up, &IntersectionNode::t_up);
    }
};

class IntersectionListDown : public IntersectionList {
public:
    IntersectionListDown()
        : IntersectionList()
    {
        head->t_down = -Infinity;
        rear->t_down = Infinity;

        head->next_down = rear;
        rear->prev_down = head;
    }

    void insert(IntersectionNode* node) const
    {
        IntersectionList::insert(node, &IntersectionNode::prev_down, &IntersectionNode::next_down, &IntersectionNode::t_down);
    }

    void interpolate(double t, double& dist, double& angle) const
    {
        IntersectionList::interpolate(t, dist, angle, &IntersectionNode::next_down, &IntersectionNode::t_down);
    }
};

struct Quad {
    Vec3 p0;
    Vec3 p1;
    Vec3 p2;
    Vec3 p3;

    double area;

    Quad()
        : p0(Vec3::Zero())
        , p1(Vec3::Zero())
        , p2(Vec3::Zero())
        , p3(Vec3::Zero())
        , area(0.) { };

    Quad(Vec3 const& p0, Vec3 const& p1, Vec3 const& p2, Vec3 const& p3)
        : p0(p0)
        , p1(p1)
        , p2(p2)
        , p3(p3)
    {
        area = 0.5 * ((p0 - p1).cross(p0 - p2).norm() + (p3 - p1).cross(p3 - p2).norm());
    }
};

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;
    ParametricSurface const& m_surface;

    int m_num_paths;
    int m_num_particles;
    double m_l_turn;

    struct OrderNode {
        int path;
        double u;

        OrderNode* next;
        OrderNode* prev;
    };

public:
    Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles);

    void generate_winding_order();

    auto intersect(
        Vec2 up1, Vec2 up2,
        Vec2 down1, Vec2 down2,
        Vec2& intersection,
        double& t_up, double& t_down) -> bool;

    auto intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down) -> std::vector<IntersectionNode*>;

    auto score(IntersectionNode* p) -> double;
    void score_ends(IntersectionNode* p);

    void interpolate(
        LocalFrame& frame,
        std::vector<IntersectionListUp> const& up,
        std::vector<IntersectionListDown> const& down,
        int i,
        double t);

    double max_area = -Infinity;

    std::vector<Quad> quads;
    Quad max_quad;

    auto simulate(
        double timestep = 0.05,
        double spring_constant = 1'000'000.,
        double damping_coefficient = 200.,
        double epsilon = 0.001,
        int step = 0);

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
    decltype(auto) surface() const { return m_surface; }
};

#endif