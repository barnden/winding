#pragma once
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Surfaces/Surface.h"
#include "utils.h"

class Simulator {
    int m_size;
    ParametricSurface const& m_surface;

protected:
    double m_ksp;
    double m_kdp;
    double m_kpf;
    double m_kbs;
    double m_dt;
    double m_t;
    double m_v_eps;

    std::vector<Vec2> m_position_initial;

public:
    Simulator(ParametricSurface const& f, std::vector<Vec2> const& init_path);
    ~Simulator();

    virtual void step() = 0;
    virtual bool stop() = 0;

    void simulate(int k);

    decltype(auto) size() const { return m_size; }
    decltype(auto) surface() { return m_surface; }

    auto& ksp() { return m_ksp; }
    auto& kdp() { return m_kdp; }
    auto& kpf() { return m_kpf; }
    auto& kbs() { return m_kbs; }
    auto& dt() { return m_dt; }
    auto& eps() { return m_v_eps; }
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

#endif
