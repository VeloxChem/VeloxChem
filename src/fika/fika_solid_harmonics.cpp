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

#include "fika_solid_harmonics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <utility>

namespace fika {

namespace {

/// Row of S_{l,m} in the buffer (after the R and R^2 rows).
constexpr auto row_index(int l, int m) noexcept -> std::size_t {
  return static_cast<std::size_t>(2 + l * l + m + l);
}

/// Order L + 1 from orders L and L - 1 (L >= 1) for the first n columns:
///   S_{L+1, L+1}  = c (x S_{L,L} - y S_{L,-L}),  S_{L+1,-L-1} = c (y S_{L,L} + x S_{L,-L}),
///   S_{L+1, M}    = a_M z S_{L,M} - b_M R^2 S_{L-1,M}  (|M| <= L),
/// with c = sqrt((2L+1)/(2L+2)), a_M = (2L+1)/d_M, b_M = sqrt((L+M)(L-M))/d_M and
/// d_M = sqrt((L+M+1)(L-M+1)).
template <int L>
void recursion_step(double* rows, std::size_t stride, std::size_t n) {
  const double* x = rows + row_index(1, 1) * stride;
  const double* y = rows + row_index(1, -1) * stride;
  const double* z = rows + row_index(1, 0) * stride;
  const double* r2 = rows + stride;

  const double c = std::sqrt((2.0 * L + 1.0) / (2.0 * L + 2.0));
  const double* top = rows + row_index(L, L) * stride;
  const double* bottom = rows + row_index(L, -L) * stride;
  double* next_top = rows + row_index(L + 1, L + 1) * stride;
  double* next_bottom = rows + row_index(L + 1, -L - 1) * stride;
  for (std::size_t j = 0; j < n; ++j) {
    next_top[j] = c * (x[j] * top[j] - y[j] * bottom[j]);
    next_bottom[j] = c * (y[j] * top[j] + x[j] * bottom[j]);
  }

  for (int m = -L; m <= L; ++m) {
    const double d = std::sqrt(static_cast<double>((L + m + 1) * (L - m + 1)));
    const double a = (2.0 * L + 1.0) / d;
    const double b = std::sqrt(static_cast<double>((L + m) * (L - m))) / d;
    const double* current = rows + row_index(L, m) * stride;
    double* next = rows + row_index(L + 1, m) * stride;
    if (m == L || m == -L) {
      for (std::size_t j = 0; j < n; ++j) {
        next[j] = a * z[j] * current[j];
      }
    } else {
      const double* previous = rows + row_index(L - 1, m) * stride;
      for (std::size_t j = 0; j < n; ++j) {
        next[j] = a * z[j] * current[j] - b * r2[j] * previous[j];
      }
    }
  }
}

using Step = void (*)(double*, std::size_t, std::size_t);

/// recursion_step<L> for L = 1..max_solid_harmonic_order - 1 (index L - 1).
constexpr auto steps = []<int... Ls>(std::integer_sequence<int, Ls...>) {
  return std::array<Step, sizeof...(Ls)>{&recursion_step<Ls + 1>...};
}(std::make_integer_sequence<int, max_solid_harmonic_order - 1>{});

}  // namespace

void SolidHarmonics::compute(std::span<const Point3D<double>> separations,
                             std::span<const std::size_t> counts) {
  assert(!counts.empty() && counts.size() <= max_solid_harmonic_order + 1);
  assert(counts[0] <= separations.size());
  assert(std::ranges::is_sorted(counts, std::greater<>{}));

  counts_.assign(counts.begin(), counts.end());
  const int top_order = max_order();
  const std::size_t n = counts[0];
  stride_ = (n + 7) / 8 * 8;
  const auto row_count = static_cast<std::size_t>(2 + (top_order + 1) * (top_order + 1));
  if (rows_.size() < row_count * stride_) {
    rows_.resize(row_count * stride_);
  }
  double* rows = rows_.data();

  double* r = rows;
  double* r2 = rows + stride_;
  double* s00 = rows + row_index(0, 0) * stride_;
  for (std::size_t j = 0; j < n; ++j) {
    const Point3D<double>& v = separations[j];
    r2[j] = v.x * v.x + v.y * v.y + v.z * v.z;
    r[j] = std::sqrt(r2[j]);
    s00[j] = 1.0;
  }
  if (top_order >= 1) {
    double* y = rows + row_index(1, -1) * stride_;
    double* z = rows + row_index(1, 0) * stride_;
    double* x = rows + row_index(1, 1) * stride_;
    for (std::size_t j = 0; j < counts[1]; ++j) {
      y[j] = separations[j].y;
      z[j] = separations[j].z;
      x[j] = separations[j].x;
    }
  }
  for (int l = 1; l < top_order; ++l) {
    steps[static_cast<std::size_t>(l - 1)](rows, stride_, counts[static_cast<std::size_t>(l + 1)]);
  }
}

}  // namespace fika
