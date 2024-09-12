#pragma once
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Surface.h"
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

    std::vector<Vec2> m_p_init;

public:
    Simulator(ParametricSurface const& f, std::vector<Vec2> const& init_path);
    ~Simulator();

    virtual void step() = 0;
    virtual bool stop() = 0;

    void simulate(int k);

    decltype(auto) size() { return m_size; }
    decltype(auto) surface() { return m_surface; }
};

class OffSurface : public Simulator {
    std::vector<Vec2> m_p;
    std::vector<Vec2> m_v;

    std::vector<bool> m_touching;
    std::vector<int> m_l;
    std::vector<int> m_r;

public:
    OffSurface(ParametricSurface const& f, std::vector<Vec2> const& init_path);
    virtual void step();
    virtual bool stop();

    void mapping();
    void lifting();
    void landing();
};

#endif