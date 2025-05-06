/*
 * Copyright (c) 2025, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#include <format>
#include "utils.h"

auto operator<<(std::ostream& ostream, Vec2 const& vec) -> std::ostream& {
    return ostream << std::format("{}, {}", vec.x(), vec.y());
}
auto operator<<(std::ostream& ostream, Vec3 const& vec) -> std::ostream& {
    return ostream << std::format("{}, {}, {}", vec.x(), vec.y(), vec.z());
}