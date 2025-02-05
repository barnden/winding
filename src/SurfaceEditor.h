/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include "Composer.h"
#include "Surfaces/BSpline.h"
#include <string>

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
#if USE_QPMAD
        std::cout << "[Editor] Using qpmad.\n";
#else
        std::cout << "[Editor] Using Epigraph.\n";
#endif
    }

    void step(double dt, double ksp, double kdp, double eps);
};