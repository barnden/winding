#pragma once

#include <memory>
#include <print>
#include <vector>

#include "Surfaces/Surface.h"
#include "ThreadPool.h"
#include "utils.h"

#define SF_OFF            0.001
#define ANGLE_SPEED_CONST 2.0

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
    using ptr_t = std::shared_ptr<IntersectionNode>;

    Vec2 parametric = Vec2::Zero();
    Vec3 point = Vec3::Zero();

    double angle_score = 0.;
    double dist_score = 0.;
    double t_up = 0.;
    double t_down = 0.;

    ptr_t prev_up = nullptr;
    ptr_t next_up = nullptr;

    ptr_t prev_down = nullptr;
    ptr_t next_down = nullptr;

    bool is_end = false;
};

class IntersectionList {
public:
    using node_t = IntersectionNode;
    using nodeptr_t = IntersectionNode::ptr_t;

    nodeptr_t head;
    nodeptr_t rear;

    IntersectionList();

    void insert(
        nodeptr_t node,
        nodeptr_t node_t::* prev,
        nodeptr_t node_t::* next,
        double node_t::* value) const;

    void interpolate(
        double t,
        double& dist,
        double& angle,
        nodeptr_t node_t::* next,
        double node_t::* value) const;
};

class IntersectionListUp : public IntersectionList {
public:
    IntersectionListUp();

    void insert(IntersectionNode::ptr_t node) const;
    void interpolate(double t, double& dist, double& angle) const;
};

class IntersectionListDown : public IntersectionList {
public:
    IntersectionListDown();

    void insert(IntersectionNode::ptr_t node) const;
    void interpolate(double t, double& dist, double& angle) const;
};

struct Quad {
    Vec3 p0 = Vec3::Zero();
    Vec3 p1 = Vec3::Zero();
    Vec3 p2 = Vec3::Zero();
    Vec3 p3 = Vec3::Zero();

    double area = 0.;

    Quad() = default;

    Quad(Vec3 const& p0, Vec3 const& p1, Vec3 const& p2, Vec3 const& p3)
        : p0(p0)
        , p1(p1)
        , p2(p2)
        , p3(p3)
    {
        area = 0.5 * ((p0 - p1).cross(p0 - p2).norm() + (p3 - p1).cross(p3 - p2).norm());
    }
};

struct OrderNode {
    int path;
    double u;

    OrderNode* next;
    OrderNode* prev;
};

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;
    ParametricSurface const& m_surface;

    int m_num_paths;
    int m_num_particles;

    ThreadPool<std::function<void()>> m_thread_pool;

public:
    Composer(ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles);

    void generate_winding_order();

    auto intersect(
        Vec2 up1, Vec2 up2,
        Vec2 down1, Vec2 down2,
        Vec2& intersection,
        double& t_up, double& t_down) -> bool;

    auto intersect(
        std::vector<Vec2> const& up,
        std::vector<Vec2> const& down) -> std::vector<IntersectionNode::ptr_t>;

    auto score(IntersectionNode::ptr_t p) -> double;
    void score_ends(IntersectionNode::ptr_t p);

    void interpolate(
        LocalFrame& frame,
        std::vector<IntersectionListUp> const& up,
        std::vector<IntersectionListDown> const& down,
        int i,
        double t);

    std::vector<Quad> quads;

    auto simulate(
        double timestep = 0.05,
        double spring_constant = 1'000'000.,
        double damping_coefficient = 200.,
        double epsilon = 0.001,
        int step = 0) -> std::vector<LocalFrame>;

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
    decltype(auto) surface() const { return m_surface; }
};
