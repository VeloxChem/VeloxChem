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

#ifndef fika_solid_harmonics_hpp
#define fika_solid_harmonics_hpp

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

#include "fika_point3d.hpp"

namespace fika {

/// Highest order of the solid harmonics: twice the highest angular momentum of a shell, plus the
/// rank of the highest MM multipole (quadrupoles).
inline constexpr int max_solid_harmonic_order = 18;

/// Racah-normalized real solid harmonics S_{L,M}(R), M = -L..L, of a set of separation vectors,
/// with R = |R| and R^2, stored as rows over (L, M) and columns over separations.
///
/// S_{0,0} = 1, S_{1,M} = (y, z, x), S_{L,0} = R^L P_L(cos theta) and sum_M S_{L,M}^2 = R^(2L);
/// higher orders follow Helgaker's recursion (Molecular Electronic-Structure Theory, 6.4).
/// Order L is computed only for the first counts[L] separations, so higher orders can cover
/// fewer (the closest) separations. Storage is reused across calls.
class SolidHarmonics {
 public:
  /// counts[L]: separations needing order L, non-increasing, counts[0] <= separations.size();
  /// orders 0..counts.size() - 1 (at most max_solid_harmonic_order). R and R^2 are computed for
  /// counts[0] separations.
  void compute(std::span<const Point3D<double>> separations, std::span<const std::size_t> counts);

  auto max_order() const noexcept -> int { return static_cast<int>(counts_.size()) - 1; }

  /// Separations covered by order l.
  auto count(int l) const noexcept -> std::size_t {
    assert(l >= 0 && l <= max_order());
    return counts_[static_cast<std::size_t>(l)];
  }

  auto distances() const noexcept -> std::span<const double> { return row(0, counts_[0]); }

  auto distances_squared() const noexcept -> std::span<const double> { return row(1, counts_[0]); }

  /// S_{l,m} for the first count(l) separations.
  auto values(int l, int m) const noexcept -> std::span<const double> {
    assert(l >= 0 && l <= max_order() && m >= -l && m <= l);
    return row(static_cast<std::size_t>(2 + l * l + m + l), count(l));
  }

 private:
  auto row(std::size_t index, std::size_t size) const noexcept -> std::span<const double> {
    return std::span(rows_).subspan(index * stride_, size);
  }

  std::vector<std::size_t> counts_{0};
  std::size_t stride_ = 0;
  std::vector<double> rows_;  // (2 + (max_order + 1)^2) rows of stride_ values
};

}  // namespace fika

#endif  // fika_solid_harmonics_hpp
