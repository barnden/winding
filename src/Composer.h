#pragma once
#ifndef COMPOSER_H
#    define COMPOSER_H

#    include "Simulator.h"
#    include "Surfaces/Surface.h"
#    include <fstream>
#    include <iostream>
#    include <vector>

#    define SF_OFF 0.001
// #define SPEED_DEFINED_BY_LENGTH
#    define ANGLE_SPEED_CONST 2.0

#    include "utils.h"

class Composer {
    std::vector<std::vector<Vec2>> m_initial;
    std::vector<int> m_winding_order;

    std::shared_ptr<Options> m_options;
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

    class IntersectionListUp {
    public:
        IntersectionNode *head, *rear;
        IntersectionListUp()
        {
            head = new IntersectionNode();
            rear = new IntersectionNode();
            head->next_up = rear;
            rear->prev_up = head;
            head->t_up = -1.0E100;
            rear->t_up = 1.0E100;
            head->is_end = true;
            rear->is_end = true;
        }
        void insert(IntersectionNode* node)
        {
            IntersectionNode* p = head;
            while (p->next_up != rear && p->next_up->t_up < node->t_up) {
                p = p->next_up;
            }
            node->prev_up = p;
            node->next_up = p->next_up;
            p->next_up->prev_up = node;
            p->next_up = node;
        }
        void interpolate(double t, double& dist, double& angle) const
        {
            IntersectionNode* p = head->next_up;

            // if (p->next_up == nullptr) {
            //     dist = p->dist_score;
            //     angle = p->angle_score;

            //     return;
            // }

            // while (p->next_up->next_up != rear && t > p->next_up->t_up) {
            //     if (p->next_up == nullptr)
            //         break;
            //     p = p->next_up;
            // }
            double lambda = (t - p->t_up) / (p->next_up->t_up - p->t_up);
            dist = lambda * p->next_up->dist_score + (1 - lambda) * p->dist_score;
            angle = lambda * p->next_up->angle_score + (1 - lambda) * p->angle_score;
        }
    };

    class IntersectionListDown {
    public:
        IntersectionNode *head, *rear;
        IntersectionListDown()
        {
            head = new IntersectionNode();
            rear = new IntersectionNode();
            head->next_down = rear;
            rear->prev_down = head;
            head->t_down = -1.0E100;
            rear->t_down = 1.0E100;
            head->is_end = true;
            rear->is_end = true;
        }
        void insert(IntersectionNode* node)
        {
            IntersectionNode* p = head;
            while (p->next_down != rear && p->next_down->t_down < node->t_down) {
                p = p->next_down;
            }
            node->prev_down = p;
            node->next_down = p->next_down;
            p->next_down->prev_down = node;
            p->next_down = node;
        }
        void interpolate(double t, double& dist, double& angle) const
        {
            IntersectionNode* p = head->next_down;

            // if (p->next_down == nullptr) {
            //     dist = p->dist_score;
            //     angle = p->angle_score;

            //     return;
            // }
            // while (p->next_down->next_down != rear && t > p->next_down->t_down) {
            //     if (p->next_down == nullptr)
            //         break;
            //     p = p->next_down;
            // }
            double lambda = (t - p->t_down) / (p->next_down->t_down - p->t_down);
            dist = lambda * p->next_down->dist_score + (1 - lambda) * p->dist_score;
            angle = lambda * p->next_down->angle_score + (1 - lambda) * p->angle_score;
        }
    };

public:
    struct LocalFrame {
        Vec3 position;
        Vec3 normal;
        Vec3 direction;
        Vec2 shadow;

        double l;
        int type;

        double distance;
        double path_distance;
        double angle;
    };

    Composer(std::shared_ptr<Options> const& options, ParametricSurface const& surface, double num_revolutions, int num_paths, int num_particles)
        : m_options(options)
        , m_surface(surface)
        , m_num_paths(num_paths)
        , m_num_particles(num_particles)
        , m_l_turn(0.3)
    {
        auto angle = (surface.v_max() - surface.v_min()) / (surface.u_max() - surface.u_min()) / num_revolutions;
        auto du = (surface.u_max() - surface.u_min()) / num_paths;
        auto ddv = (surface.v_max() - surface.v_min()) / (num_particles + 1);
        auto ddu = ddv / angle;

        m_initial = std::vector<std::vector<Vec2>>(2 * num_paths, std::vector<Vec2>(num_particles));

        for (auto i = 0; i < num_paths; i++) {
            Vec2 p(surface.u_min() + i * du, surface.v_min());
            Vec2 q(surface.u_min() + i * du, surface.v_min());

            Vec3 cur = surface.f(p) + SF_OFF * surface.normal(p);

            for (auto j = 0; j < num_particles; j++) {
                cur = surface.f(p) + SF_OFF * surface.normal(p);
                m_initial[2 * i][num_particles - 1 - j] = p;
                m_initial[2 * i + 1][j] = q;

                p += Vec2(ddu, ddv);
                q += Vec2(-ddu, ddv);
            }
        }
    }

    void generate_winding_order();
    bool intersect(
        Vec2 up1,
        Vec2 up2,
        Vec2 down1,
        Vec2 down2,
        Vec2& intersection,
        double& t_up,
        double& t_down);
    std::vector<Composer::IntersectionNode*> intersect(std::vector<Vec2> const& up, std::vector<Vec2> const& down);

    decltype(auto) score(IntersectionNode* p)
    {
        if (p->prev_up->is_end || p->prev_down->is_end || p->next_up->is_end || p->next_down->is_end)
            return 0.;

        Vec3 v1 = (p->prev_up->point - p->next_up->point).normalized();
        Vec3 v2 = (p->prev_down->point - p->next_down->point).normalized();
        p->angle_score = abs(v1.dot(v2));
        p->dist_score = ((p->prev_up->point - p->point).norm()
                         + (p->prev_down->point - p->point).norm()
                         + (p->next_up->point - p->point).norm()
                         + (p->next_down->point - p->point).norm())
                        / 2.;

        return p->dist_score;
    }

    void score_ends(IntersectionNode* p)
    {
        if (p->is_end)
            return;

        double sum_dist = 0.;
        double sum_angle = 0.;
        int ct = 0;

        if (!p->prev_up->is_end) {
            sum_dist += p->prev_up->dist_score;
            sum_angle += p->prev_up->angle_score;

            ct++;
        }

        if (!p->prev_down->is_end) {
            sum_dist += p->prev_down->dist_score;
            sum_angle += p->prev_down->angle_score;

            ct++;
        }

        if (!p->next_up->is_end) {
            sum_dist += p->next_up->dist_score;
            sum_angle += p->next_up->angle_score;

            ct++;
        }

        if (!p->next_down->is_end) {
            sum_dist += p->next_down->dist_score;
            sum_angle += p->next_down->angle_score;

            ct++;
        }

        if (ct != 0) {
            p->dist_score = sum_dist / ct;
            p->angle_score = sum_angle / ct;
        }
    }

    void interpolate(
        LocalFrame& frame,
        std::vector<IntersectionListUp> const& up,
        std::vector<IntersectionListDown> const& down,
        int i,
        double t)
    {
        if (i % 2) {
            up[i / 2].interpolate(t, frame.path_distance, frame.angle);
        } else {
            down[i / 2].interpolate(t, frame.path_distance, frame.angle);
        }
    }

    double max_area = -INFINITY;
    Eigen::Matrix<double, 12, 1> max_quad;

    decltype(auto) simulate(
        double dt = 0.05,
        double ksp = 1'000'000.,
        double kdp = 200.,
        double eps = 0.001,
        int step = 1234567)
    {
        generate_winding_order();

        auto paths = std::vector<std::vector<Vec2>>(m_initial.size());
        auto orders = std::vector<std::vector<int>>(m_initial.size());

        for (auto&& [j, i] : enumerate(m_winding_order)) {
            auto simulator = OffSurface(m_options, m_surface, m_initial[i]);

            simulator.ksp() = ksp;
            simulator.kdp() = kdp;
            simulator.dt() = dt;
            simulator.eps() = eps;

            simulator.simulate(1000);
            simulator.mapping();

            paths[i] = simulator.p();
            orders[i] = simulator.r();
        }

        std::vector<IntersectionListUp> up_list(m_initial.size() / 2);
        std::vector<IntersectionListDown> down_list(m_initial.size() / 2);

        for (auto i = 0u; i < m_initial.size() / 2; i++) {
            for (auto j = 0u; j < m_initial.size() / 2; j++) {
                auto result = intersect(paths[2 * i + 1], paths[2 * j]);

                for (auto&& p : result) {
                    up_list[i].insert(p);
                    down_list[j].insert(p);
                }
            }
        }

        for (auto&& ls : up_list) {
            for (auto p = ls.head->next_up; p != ls.rear; p = p->next_up)
                score(p);
        }

        for (auto&& ls : up_list) {
            score_ends(ls.head->next_up);
            score_ends(ls.rear->prev_up);
        }

        for (auto&& ls : down_list) {
            score_ends(ls.head->next_down);
            score_ends(ls.rear->prev_down);
        }
        // std::cout << "begin quadmesh\n";
        std::ofstream ofs(m_options->out_path + "/" + m_options->experiment + "/obj/step-" + std::to_string(step) + ".obj");
        for (auto& ls : up_list) {
            auto p = ls.head;
            while (p != ls.rear) {
                if (p->next_up && p->next_up->next_down && p->next_up->next_down->prev_up && p->next_up->next_down->prev_up->prev_down == p) {
                    Vec3 p0 = p->point;
                    Vec3 p1 = p->next_up->point;
                    Vec3 p2 = p->next_up->next_down->point;
                    Vec3 p3 = p->next_up->next_down->prev_up->point;

                    double area = 0.5 * ((p1 - p3).cross(p0 - p1).norm() + (p2 - p3).cross(p2 - p1).norm());
                    if (area > max_area) {
                        max_area = area;
                        max_quad.segment<3>(0) = p0;
                        max_quad.segment<3>(3) = p1;
                        max_quad.segment<3>(6) = p2;
                        max_quad.segment<3>(9) = p3;
                    }

                    ofs << "v " << p0.transpose() << '\n';
                    ofs << "v " << p1.transpose() << '\n';
                    ofs << "v " << p2.transpose() << '\n';
                    ofs << "v " << p3.transpose() << '\n';
                }
                p = p->next_up;
            }
        }

        std::cout << "max_area: " << max_area << '\n';
        std::cout << "max_quad: " << max_quad.transpose() << '\n';

        int cnt = 1;
        for (auto& ls : up_list) {
            auto p = ls.head;
            while (p != ls.rear) {
                if (p->next_up && p->next_up->next_down && p->next_up->next_down->prev_up && p->next_up->next_down->prev_up->prev_down == p) {
                    ofs << "f " << cnt << " " << cnt + 1 << " " << cnt + 2 << " " << cnt + 3 << std::endl;
                    cnt += 4;
                }
                p = p->next_up;
            }
        }

        std::cout << "wrote quadmesh\n";
        // exit(EXIT_SUCCESS);

        auto motion = std::vector<LocalFrame>(2 * m_num_paths * (m_num_particles + 2) - 2);
        auto ct = 0;
        auto l = 0.; // what is this?

        for (auto&& [j, i] : enumerate(m_winding_order)) {
            auto const& path = paths[i];
            auto const& order = orders[i];

            auto cur = 0;
            auto type = 1;
            Vec3 pcur = m_surface.f(path[cur]);
            Vec3 pnext = Vec3::Zero();
            Vec3 ncur = m_surface.normal(path[cur]);
            Vec3 nnext = Vec3::Zero();

            while (cur < m_num_particles - 1) {
                pnext = m_surface.f(path[cur + order[cur]]);
                nnext = m_surface.normal(path[cur + order[cur]]);

                motion[ct].position = pcur + SF_OFF * ncur;
                motion[ct].normal = ncur;
                motion[ct].direction = (pnext - pcur).normalized();
                motion[ct].shadow = path[cur];
                motion[ct].l = l;
                motion[ct].type = type;
                motion[ct].distance = 0.;
                interpolate(motion[ct], up_list, down_list, i, cur);

                ct++;

                Vec3 n1 = Vec3(ncur.x(), ncur.y(), 0.).normalized();
                Vec3 n2 = Vec3(nnext.x(), nnext.y(), 0.).normalized();

                double dl = acos(std::clamp(n1.dot(n2), -1., 1.) / (order[cur] * ANGLE_SPEED_CONST));
                l += dl;

                for (int j = 1; j < order[cur]; j++) {
#    define LERP(a, b, start, end) (((double)(a - b) * start + (double)(b) * end) / (double)(a))
                    motion[ct].position = LERP(order[cur], j, pcur, pnext) + SF_OFF * ncur;
                    motion[ct].normal = LERP(order[cur], j, ncur, nnext).normalized();
                    motion[ct].direction = (pnext - pcur).normalized();
                    motion[ct].shadow = path[cur + j];
                    motion[ct].l = l;
                    motion[ct].type = type;
                    motion[ct].distance = (motion[ct].position - m_surface.f(path[cur + j])).norm();
                    interpolate(motion[ct], up_list, down_list, i, cur + j);

                    ct++;
                    l += dl;
                }

                cur += order[cur];
                pcur = pnext;
                ncur = nnext;
            }

            motion[ct].position = pcur + SF_OFF * ncur;
            motion[ct].normal = ncur;
            motion[ct].direction = motion[ct - 1].direction;
            motion[ct].shadow = path[cur];
            motion[ct].l = l;
            motion[ct].type = type;
            motion[ct].distance = 0.;
            interpolate(motion[ct], up_list, down_list, i, m_num_particles - 1);

            ct++;
            l += m_l_turn;

            if (j != (int)(m_winding_order.size() - 1)) {
                pnext = m_surface.f(m_initial[m_winding_order[j + 1]][0]);
                nnext = m_surface.normal(m_initial[m_winding_order[j + 1]][0]);

                motion[ct].normal = motion[ct - 1].direction;
                motion[ct].position = pcur + SF_OFF * motion[ct].normal;
                motion[ct].direction = (pnext - pcur).normalized();
                motion[ct].shadow = motion[ct - 1].shadow;
                motion[ct].l = l;
                motion[ct].type = -1;
                motion[ct].distance = 0.;
                interpolate(motion[ct], up_list, down_list, i, -1);

                ct++;
                Vec3 n1 = Vec3(ncur.x(), ncur.y(), 0.).normalized();
                Vec3 n2 = Vec3(nnext.x(), nnext.y(), 0.).normalized();

                l += acos(std::clamp(n1.dot(n2), -1., 1.) / (5. * ANGLE_SPEED_CONST));

                motion[ct].normal = motion[ct - 1].normal;
                motion[ct].position = pnext + SF_OFF * motion[ct].normal;
                motion[ct].direction = motion[ct - 1].direction;
                motion[ct].shadow = m_initial[m_winding_order[j + 1]][0];
                motion[ct].l = l;
                motion[ct].type = -1;
                motion[ct].distance = 0.;
                interpolate(motion[ct], up_list, down_list, i, -1);
                ct++;
                l += m_l_turn;
            }
        }

        return motion;
    }

    decltype(auto) initial_path() const { return m_initial; }
    decltype(auto) winding_order() const { return m_winding_order; }
    decltype(auto) surface() const { return m_surface; }
};

#endif