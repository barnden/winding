/*
 * Copyright (c) 2024-2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include "Surfaces/BSpline.h"

class SurfaceEditor {
    CubicBSpline m_surface;

    double m_num_revolutions;
    int m_num_paths;
    int m_num_particles;
    int m_step;

public:
    SurfaceEditor(CubicBSpline& reference_surface, double num_revolutions, int num_paths, int num_particles)
        : m_surface(reference_surface)
        , m_num_revolutions(num_revolutions)
        , m_num_paths(num_paths)
        , m_num_particles(num_particles)
        , m_step(0)
    {
        std::println("[Editor] Using {}.\n", USE_QPMAD ? "qpmad" : "Epigraph");
    }

    void step(double timestep, double spring_constant, double damping_coefficient, double epsilon);
};