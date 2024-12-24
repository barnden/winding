/*
 * Copyright (c) 2024, Brandon G. Nguyen <brandon@nguyen.vc>
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */
#pragma once

#include <cstddef>
#include <cstdint>

#ifdef __BMI2__
#    include <immintrin.h>
#endif

namespace Morton {
using Code = uint64_t;
#ifdef __BMI2__
uint64_t constexpr MORTON_MASK = 0x1249249249249249;
Code Encode(uint32_t x, uint32_t y, uint32_t z)
{
    return _pdep_u64(x, MORTON_MASK) | _pdep_u64(y, MORTON_MASK << 1) | _pdep_u64(z, MORTON_MASK << 2);
}
#else
Code _expand(uint64_t n)
{
    n = (n | n << 32) & 0x001f00000000ffff;
    n = (n | n << 16) & 0x001f0000ff0000ff;
    n = (n | n << 8) & 0x100f00f00f00f00f;
    n = (n | n << 4) & 0x10c30c30c30c30c3;
    n = (n | n << 2) & 0x1249249249249249;

    return n;
}

Code Encode(uint64_t x, uint64_t y, uint64_t z)
{
    return (_expand(x) << 2) | (_expand(y) << 1) | (_expand(z) << 0);
}
#endif
};
