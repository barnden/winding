#pragma once

#include "Surfaces/Surface.h"
#include "utils.h"

class Simulator {
    int m_size;
    ParametricSurface const& m_surface;

protected:
    double m_spring_constant;
    double m_damping_coefficient;
    double m_timestep;
    double m_t;
    double m_epsilon;

    std::vector<Vec2> m_position_initial;

public:
    Simulator(ParametricSurface const& f, std::vector<Vec2> const& init_path);
    ~Simulator() = default;

    virtual void step() = 0;
    virtual bool stop() = 0;

    int simulate(int num_iterations);

    decltype(auto) size() const { return m_size; }
    decltype(auto) surface() { return m_surface; }

    auto& spring_constant() { return m_spring_constant; }
    auto& damping_coefficient() { return m_damping_coefficient; }
    auto& timestep() { return m_timestep; }
    auto& epsilon() { return m_epsilon; }
};

class OffSurface : public Simulator {
    std::vector<Vec2> m_position;
    std::vector<Vec2> m_velocity;

    std::vector<bool> m_touching;
    std::vector<int> m_l;
    std::vector<int> m_r;

public:
    OffSurface(ParametricSurface const& f, std::vector<Vec2> const& init_path);
    virtual void step() override;
    virtual bool stop() override;

    void mapping();
    void lifting();
    void landing();

    auto& p() { return m_position; }
    auto& v() { return m_velocity; }
    auto& touching() { return m_touching; }
    auto& r() { return m_r; }
    auto& l() { return m_l; }
};
