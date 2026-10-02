//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT

#ifndef fika_octree_keys_hpp
#define fika_octree_keys_hpp

// Internal: Morton keys, bounding cubes and key sorting for the octrees of the fast multipole
// methods.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "fika_point3d.hpp"

namespace fika::detail {

inline constexpr int source_depth = 21;  // Morton bits per axis

// Spreads the low 21 bits of v to every third bit.
inline auto spread_bits(std::uint64_t v) -> std::uint64_t {
  v &= 0x1fffff;
  v = (v | v << 32) & 0x1f00000000ffff;
  v = (v | v << 16) & 0x1f0000ff0000ff;
  v = (v | v << 8) & 0x100f00f00f00f00f;
  v = (v | v << 4) & 0x10c30c30c30c30c3;
  v = (v | v << 2) & 0x1249249249249249;
  return v;
}

// Interleaves x, y, z (x most significant within each octal digit).
inline auto morton(std::uint64_t x, std::uint64_t y, std::uint64_t z) -> std::uint64_t {
  return spread_bits(x) << 2 | spread_bits(y) << 1 | spread_bits(z);
}

inline auto compact_bits(std::uint64_t v) -> std::uint64_t {
  v &= 0x1249249249249249;
  v = (v ^ (v >> 2)) & 0x10c30c30c30c30c3;
  v = (v ^ (v >> 4)) & 0x100f00f00f00f00f;
  v = (v ^ (v >> 8)) & 0x1f0000ff0000ff;
  v = (v ^ (v >> 16)) & 0x1f00000000ffff;
  v = (v ^ (v >> 32)) & 0x1fffff;
  return v;
}

// Cube enclosing the points: centre and half edge (slightly padded, positive).
inline auto bounding_cube(std::span<const Point3D<double>> points)
    -> std::pair<Point3D<double>, double> {
  Point3D<double> low{std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity()};
  Point3D<double> high{-low.x, -low.y, -low.z};
  for (const auto& point : points) {
    low = {std::min(low.x, point.x), std::min(low.y, point.y), std::min(low.z, point.z)};
    high = {std::max(high.x, point.x), std::max(high.y, point.y), std::max(high.z, point.z)};
  }
  const Point3D<double> centre{0.5 * (low.x + high.x), 0.5 * (low.y + high.y),
                               0.5 * (low.z + high.z)};
  const double extent = std::max({high.x - low.x, high.y - low.y, high.z - low.z});
  return {centre, 0.5 * extent * (1.0 + 1e-9) + 1e-9};
}

// Stable LSD radix sort of (key, index) pairs by the low 3 * source_depth bits of the key, in
// 16-bit digits: pairs made in index order end sorted by (key, index), as std::sort would.
inline void radix_sort(std::vector<std::pair<std::uint64_t, std::uint32_t>>& items) {
  constexpr int digit_bits = 16;
  constexpr std::size_t buckets = std::size_t{1} << digit_bits;
  std::vector<std::pair<std::uint64_t, std::uint32_t>> buffer(items.size());
  std::vector<std::size_t> offsets(buckets + 1);
  for (int shift = 0; shift < 3 * source_depth; shift += digit_bits) {
    std::fill(offsets.begin(), offsets.end(), 0);
    for (const auto& item : items) {
      ++offsets[((item.first >> shift) & (buckets - 1)) + 1];
    }
    if (offsets[((items.front().first >> shift) & (buckets - 1)) + 1] == items.size()) {
      continue;  // one digit value: nothing moves
    }
    for (std::size_t b = 0; b < buckets; ++b) {
      offsets[b + 1] += offsets[b];
    }
    for (const auto& item : items) {
      buffer[offsets[(item.first >> shift) & (buckets - 1)]++] = item;
    }
    items.swap(buffer);
  }
}

}  // namespace fika::detail

#endif  // fika_octree_keys_hpp
