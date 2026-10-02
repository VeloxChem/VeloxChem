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

#include "fika_gaunt.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <numbers>
#include <vector>

#include "fika_gaussian_normalization.hpp"
#include "fika_point3d.hpp"
#include "fika_solid_harmonics.hpp"

namespace fika {

namespace {

constexpr int table_size = max_angular_momentum + 1;
constexpr int max_order = 2 * max_angular_momentum;  // largest L of a shell pair
static_assert(max_order + max_multipole_rank <= max_solid_harmonic_order);

constexpr double gaunt_cutoff = 1e-13;

/// Gauss-Legendre nodes and weights on [-1, 1].
void gauss_legendre(int n, std::vector<double>& nodes, std::vector<double>& weights) {
  nodes.resize(static_cast<std::size_t>(n));
  weights.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double x = std::cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
    double derivative = 0.0;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double p0 = 1.0;
      double p1 = x;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      derivative = n * (x * p1 - p0) / (x * x - 1.0);
      const double step = p1 / derivative;
      x -= step;
      if (std::abs(step) < 1e-16) {
        break;
      }
    }
    nodes[static_cast<std::size_t>(i)] = x;
    weights[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - x * x) * derivative * derivative);
  }
}

/// Selection rules for M given m and m' (Section 4.4 of the overlap notes): |M| is
/// ||m| - |m'|| or |m| + |m'|; two cosine-type or two sine-type factors give M >= 0, one of each
/// gives M < 0.
auto candidate_ms(int m, int m_prime) -> std::array<int, 2> {
  const bool cosine = m >= 0;
  const bool cosine_prime = m_prime >= 0;
  const int difference = std::abs(std::abs(m) - std::abs(m_prime));
  const int sum = std::abs(m) + std::abs(m_prime);
  const int sign = cosine == cosine_prime ? 1 : -1;
  return {sign * difference, sign * sum};
}

/// Coefficients of one pair (l, l'). The products X_lm X_l'm' X_LM are polynomials of degree
/// l + l' + L <= 2 (l + l') on the sphere: l + l' + 1 Gauss-Legendre points in cos(theta) and
/// 2 (l + l') + 2 equally spaced points in phi integrate them exactly.
auto build_pair(int l, int l_prime) -> std::vector<GauntEntry> {
  const int order = l + l_prime;
  const int theta_points = order + 1;
  const int phi_points = 2 * order + 2;
  std::vector<double> nodes, weights;
  gauss_legendre(theta_points, nodes, weights);
  std::vector<Point3D<double>> points;
  std::vector<double> point_weights;
  for (int i = 0; i < theta_points; ++i) {
    const double cos_theta = nodes[static_cast<std::size_t>(i)];
    const double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    for (int j = 0; j < phi_points; ++j) {
      const double phi = 2.0 * std::numbers::pi * j / phi_points;
      points.push_back({sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta});
      point_weights.push_back(weights[static_cast<std::size_t>(i)] * 2.0 * std::numbers::pi /
                              phi_points);
    }
  }
  SolidHarmonics harmonics;  // on unit vectors: the angular functions X_LM
  const std::vector<std::size_t> counts(static_cast<std::size_t>(order) + 1, points.size());
  harmonics.compute(points, counts);

  const auto integral = [&](int m, int m_prime, int big_l, int big_m) {
    const auto a = harmonics.values(l, m);
    const auto b = harmonics.values(l_prime, m_prime);
    const auto c = harmonics.values(big_l, big_m);
    double sum = 0.0;
    for (std::size_t p = 0; p < points.size(); ++p) {
      sum += point_weights[p] * a[p] * b[p] * c[p];
    }
    return sum;
  };

  std::vector<GauntEntry> entries;
  for (int m = -l; m <= l; ++m) {
    for (int m_prime = -l_prime; m_prime <= l_prime; ++m_prime) {
      const auto ms = candidate_ms(m, m_prime);
      for (int big_l = std::abs(l - l_prime); big_l <= order; big_l += 2) {
        for (std::size_t i = 0; i < ms.size(); ++i) {
          const int big_m = ms[i];
          const bool duplicate = i == 1 && ms[1] == ms[0];
          const bool excluded = (big_m == 0 && (m >= 0) != (m_prime >= 0));
          if (duplicate || excluded || std::abs(big_m) > big_l) {
            continue;
          }
          const double value =
              (2.0 * big_l + 1.0) / (4.0 * std::numbers::pi) * integral(m, m_prime, big_l, big_m);
          if (std::abs(value) > gaunt_cutoff) {
            entries.push_back({m, m_prime, big_l, big_m, value});
          }
        }
      }
    }
  }
  return entries;
}

/// One pair's coefficients, built on first request.
struct Slot {
  std::once_flag built;
  std::vector<GauntEntry> entries;
};

}  // namespace

auto gaunt_coefficients(int l, int l_prime) -> std::span<const GauntEntry> {
  assert(l >= 0 && l <= max_angular_momentum && l_prime >= 0 && l_prime <= max_angular_momentum);
  static std::array<std::array<Slot, table_size>, table_size> table;
  Slot& slot = table[static_cast<std::size_t>(l)][static_cast<std::size_t>(l_prime)];
  std::call_once(slot.built, [&] { slot.entries = build_pair(l, l_prime); });
  return slot.entries;
}

auto multipole_gaunt_coefficients(int rank, int l) -> std::span<const GauntEntry> {
  assert(rank >= 1 && rank <= max_multipole_rank && l >= 0 && l <= max_order);
  static std::array<std::array<Slot, max_order + 1>, max_multipole_rank> table;
  Slot& slot = table[static_cast<std::size_t>(rank - 1)][static_cast<std::size_t>(l)];
  std::call_once(slot.built, [&] { slot.entries = build_pair(rank, l); });
  return slot.entries;
}

}  // namespace fika
