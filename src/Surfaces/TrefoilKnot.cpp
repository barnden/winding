/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include "Surfaces/Surface.h"
#include "utils.h"

TrefoilKnot::TrefoilKnot()
    : ParametricSurface(0., 2. * PI, 0., 2. * PI - 0.001) {};

Vec3 TrefoilKnot::f(Vec2 const& p) const
{
    Vec3 e1(0., 0., 1.);
    Vec3 p0(
        sin(p.y()) + 2. * sin(2. * p.y()),
        cos(p.y()) - 2. * cos(2. * p.y()),
        -1. * sin(3. * p.y()));
    Vec3 t(
        cos(p.y()) + 4. * cos(2. * p.y()),
        -sin(p.y()) + 4. * sin(2. * p.y()),
        -3. * cos(3. * p.y()));

    t.normalize();

    e1 -= e1.dot(t) * t;
    e1.normalize();

    Vec3 e2 = t.cross(e1);

    return p0 + 0.4 * cos(p.x()) * e1 + 0.4 * sin(p.x()) * e2;
}