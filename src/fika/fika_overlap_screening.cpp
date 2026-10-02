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

#include "fika_overlap_screening.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace fika::detail {

void sort_by_distance(std::span<const AtomPair> pairs, std::span<const Point3D<double>> coordinates,
                      SortedAtomPairs& sorted, double max_distance_squared) {
  sorted.keys.clear();
  for (std::size_t i = 0; i < pairs.size(); ++i) {
    const Point3D<double>& a = coordinates[pairs[i].bra];
    const Point3D<double>& b = coordinates[pairs[i].ket];
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    const double distance_squared = dx * dx + dy * dy + dz * dz;
    if (distance_squared <= max_distance_squared) {
      sorted.keys.push_back({distance_squared, i});
    }
  }
  const std::size_t n = sorted.keys.size();
  // Keys are unique (index included), so std::sort gives the same order as a stable sort.
  std::ranges::sort(sorted.keys,
                    [](const SortedAtomPairs::SortKey& x, const SortedAtomPairs::SortKey& y) {
                      return x.distance_squared < y.distance_squared ||
                             (x.distance_squared == y.distance_squared && x.index < y.index);
                    });

  sorted.pairs.resize(n);
  sorted.separations.resize(n);
  sorted.distances_squared.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const AtomPair pair = pairs[sorted.keys[i].index];
    const Point3D<double>& a = coordinates[pair.bra];
    const Point3D<double>& b = coordinates[pair.ket];
    sorted.pairs[i] = pair;
    sorted.separations[i] = {a.x - b.x, a.y - b.y, a.z - b.z};
    sorted.distances_squared[i] = sorted.keys[i].distance_squared;
  }
}

auto ShellPairBound::operator()(double distance_squared) const noexcept -> double {
  double power = 1.0;  // R^(l_a + l_b), used only when R > 1
  if (angular_momentum > 0 && distance_squared > 1.0) {
    const double r = std::sqrt(distance_squared);
    for (int i = 0; i < angular_momentum; ++i) {
      power *= r;
    }
  }
  return prefactor * power * (constant + quadratic * distance_squared) *
         std::exp(-exponent * distance_squared);
}

auto shell_pair_bound(const BasisShell& a, const BasisShell& b) -> ShellPairBound {
  const double alpha = min_exponent(a);
  const double beta = min_exponent(b);
  const double p = alpha + beta;
  const double q = std::numbers::pi / p;
  return {max_coefficient_sum(a) * max_coefficient_sum(b) * q * std::sqrt(q), alpha * beta / p,
          angular_momentum(a) + angular_momentum(b)};
}

auto kinetic_shell_pair_bound(const BasisShell& a, const BasisShell& b) -> ShellPairBound {
  ShellPairBound bound = shell_pair_bound(a, b);
  const double alpha = max_exponent(a);
  const double beta = max_exponent(b);
  const double mu = alpha * beta / (alpha + beta);
  bound.constant = 3.0 * mu;
  bound.quadratic = 2.0 * mu * mu;
  return bound;
}

auto nuclear_attraction_shell_pair_bound(const BasisShell& a, const BasisShell& b,
                                         double charge_sum) -> ShellPairBound {
  const double alpha = min_exponent(a);
  const double beta = min_exponent(b);
  const double p = alpha + beta;
  const int l = angular_momentum(a) + angular_momentum(b);
  return {charge_sum * max_coefficient_sum(a) * max_coefficient_sum(b) * (l + 1) * 2.0 *
              std::numbers::pi / p,
          alpha * beta / p, l};
}

auto dipole_potential_shell_pair_bound(const BasisShell& a, const BasisShell& b, double dipole_sum)
    -> ShellPairBound {
  const double alpha = min_exponent(a);
  const double beta = min_exponent(b);
  const double p = alpha + beta;
  const int l = angular_momentum(a) + angular_momentum(b);
  return {dipole_sum * max_coefficient_sum(a) * max_coefficient_sum(b) * (l + 1) * 2.0 *
              std::numbers::pi / p * 2.0 * std::sqrt(p) * dipole_boys_bound,
          alpha * beta / p, l};
}

auto traceless_quadrupole_norm(const Quadrupole& quadrupole) -> double {
  const double mean = (quadrupole(0, 0) + quadrupole(1, 1) + quadrupole(2, 2)) / 3.0;
  const double xx = quadrupole(0, 0) - mean, yy = quadrupole(1, 1) - mean,
               zz = quadrupole(2, 2) - mean;
  const double xy = quadrupole(0, 1), xz = quadrupole(0, 2), yz = quadrupole(1, 2);
  // Eigenvalues of the traceless Theta: 2 s cos(phi + 2 pi k / 3), s^2 = tr(Theta^2) / 6,
  // cos(3 phi) = det(Theta) / (2 s^3).
  const double squares = xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + xz * xz + yz * yz);
  const double s = std::sqrt(squares / 6.0);
  if (s == 0.0) {
    return 0.0;
  }
  const double determinant =
      xx * (yy * zz - yz * yz) - xy * (xy * zz - yz * xz) + xz * (xy * yz - yy * xz);
  const double phi = std::acos(std::clamp(determinant / (2.0 * s * s * s), -1.0, 1.0)) / 3.0;
  double largest = 0.0;
  for (int k = 0; k < 3; ++k) {
    largest =
        std::max(largest, std::abs(2.0 * s * std::cos(phi + 2.0 * std::numbers::pi * k / 3.0)));
  }
  return largest;
}

auto quadrupole_potential_shell_pair_bound(const BasisShell& a, const BasisShell& b,
                                           double quadrupole_sum) -> ShellPairBound {
  const double alpha = min_exponent(a);
  const double beta = min_exponent(b);
  const double p = alpha + beta;
  const int l = angular_momentum(a) + angular_momentum(b);
  return {quadrupole_sum * max_coefficient_sum(a) * max_coefficient_sum(b) * (l + 1) * 4.0 *
              std::numbers::pi * quadrupole_boys_bound,
          alpha * beta / p, l};
}

namespace {

/// Largest x = R^2 >= 0 at which R^l (c0 + c1 x) exp(-mu x) is stationary: the positive root of
/// -2 mu c1 x^2 + (l c1 + 2 c1 - 2 mu c0) x + l c0 = 0 (zero if there is none). The function
/// increases below it and decreases beyond it.
auto peak_distance_squared(int l, double mu, double c0, double c1) -> double {
  if (c1 == 0.0) {
    return static_cast<double>(l) / (2.0 * mu);
  }
  const double b = l * c1 + 2.0 * c1 - 2.0 * mu * c0;
  const double root = (b + std::sqrt(b * b + 8.0 * mu * c1 * l * c0)) / (4.0 * mu * c1);
  return std::max(0.0, root);
}

/// Largest r in [low, high] with f(r) >= threshold, for f decreasing on [low, high] with
/// f(low) >= threshold > f(high).
template <typename F>
auto bisect(F&& f, double low, double high, double threshold) -> double {
  for (int iteration = 0; iteration < 200 && high - low > 1e-12 * high; ++iteration) {
    const double middle = 0.5 * (low + high);
    (f(middle) >= threshold ? low : high) = middle;
  }
  return low;
}

}  // namespace

auto ShellPairBound::cutoff_distance_squared(double threshold) const -> double {
  if (threshold <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  // For R > 1 the bound R^L (c0 + c1 R^2) exp(-mu R^2) peaks at R* (peak_distance_squared) and
  // decreases beyond max(1, R*); for R <= 1 the factor R^L is absent and the bound peaks at
  // R0 <= 1 (the L = 0 peak), decreasing beyond it.
  const double tail_start = std::max(
      1.0, std::sqrt(peak_distance_squared(angular_momentum, exponent, constant, quadratic)));
  const auto bound_at = [this](double r) { return (*this)(r * r); };
  if (bound_at(tail_start) < threshold) {
    // Only the region R < 1 can reach the threshold.
    if (quadratic == 0.0) {  // decreasing from R = 0: closed form
      if (prefactor * constant < threshold) {
        return -1.0;
      }
      return std::min(1.0, std::log(prefactor * constant / threshold) / exponent);
    }
    const double peak =
        std::min(1.0, std::sqrt(peak_distance_squared(0, exponent, constant, quadratic)));
    if (bound_at(peak) < threshold) {
      return -1.0;
    }
    const double r = bisect(bound_at, peak, 1.0, threshold);
    return r * r;
  }
  // Bisection on the decreasing tail between a point above and a point below the threshold.
  double low = tail_start;
  double high = 2.0 * tail_start;
  while (bound_at(high) >= threshold) {
    low = high;
    high *= 2.0;
  }
  const double r = bisect(bound_at, low, high, threshold);
  return r * r;
}

auto significant_pair_count(double cutoff_distance_squared,
                            std::span<const double> distances_squared) -> std::size_t {
  const auto end = std::ranges::upper_bound(distances_squared, cutoff_distance_squared);
  return static_cast<std::size_t>(end - distances_squared.begin());
}

void counts_per_order(std::span<const std::size_t> counts, std::span<const int> orders,
                      std::vector<std::size_t>& result) {
  const int max_order = orders.empty() ? -1 : *std::ranges::max_element(orders);
  result.assign(static_cast<std::size_t>(max_order + 1), 0);
  for (std::size_t i = 0; i < counts.size(); ++i) {
    // Raise only the top order, then propagate downwards once.
    auto& top = result[static_cast<std::size_t>(orders[i])];
    top = std::max(top, counts[i]);
  }
  for (int l = max_order - 1; l >= 0; --l) {
    result[static_cast<std::size_t>(l)] =
        std::max(result[static_cast<std::size_t>(l)], result[static_cast<std::size_t>(l + 1)]);
  }
}

}  // namespace fika::detail
