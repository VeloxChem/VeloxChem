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

#include "fika_basis_shell.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "fika_gaussian_normalization.hpp"

namespace fika {

namespace {

[[noreturn]] void fail(const std::string& reason) {
  throw std::invalid_argument("fika::BasisShell: " + reason);
}

void check_angular_momentum(int l) {
  if (l < 0 || l > max_angular_momentum) {
    fail("angular momentum " + std::to_string(l) + " outside 0.." +
         std::to_string(max_angular_momentum));
  }
}

void check_exponent(double exponent) {
  if (!std::isfinite(exponent) || exponent <= 0.0) {
    fail("exponent " + std::to_string(exponent) + " is not positive and finite");
  }
}

void check_input(int l, const std::vector<double>& exponents,
                 const std::vector<std::vector<double>>& contractions) {
  check_angular_momentum(l);
  if (exponents.empty()) {
    fail("no primitive exponents");
  }
  if (exponents.size() > std::numeric_limits<std::uint16_t>::max()) {
    fail("too many primitives");
  }
  for (const double exponent : exponents) {
    check_exponent(exponent);
  }
  if (contractions.empty()) {
    fail("no contracted functions");
  }
  for (const auto& contraction : contractions) {
    if (contraction.size() != exponents.size()) {
      fail(std::to_string(contraction.size()) + " coefficients for " +
           std::to_string(exponents.size()) + " exponents");
    }
    if (!std::ranges::all_of(contraction, [](double c) { return std::isfinite(c); })) {
      fail("coefficient is not finite");
    }
    if (std::ranges::all_of(contraction, [](double c) { return c == 0.0; })) {
      fail("contracted function with only zero coefficients");
    }
  }
}

/// A contracted function: nonzero effective coefficients and their primitive indices.
struct Contraction {
  std::vector<std::uint16_t> indices;
  std::vector<double> coefficients;
  double coefficient_sum = 0.0;
};

/// Overlaps of normalized primitives, s_ij = N_i N_j S_LL(a_i, a_j), row-major n x n.
template <int L>
auto normalized_primitive_overlaps(std::span<const double> exponents) -> std::vector<double> {
  const std::size_t n = exponents.size();
  std::vector<double> norms(n);
  for (std::size_t i = 0; i < n; ++i) {
    norms[i] = primitive_normalization<L>(exponents[i]);
  }
  std::vector<double> overlaps(n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const double s = norms[i] * norms[j] * primitive_overlap<L>(exponents[i], exponents[j]);
      overlaps[i * n + j] = s;
      overlaps[j * n + i] = s;
    }
  }
  return overlaps;
}

/// Validates basis-file data, sorts the exponents decreasingly (keeping coefficients with them),
/// drops zero coefficients and normalizes every contracted function: effective coefficients are
/// the file coefficients times N_l(a) / sqrt(S_kk).
auto prepare(int l, std::vector<double>& exponents,
             const std::vector<std::vector<double>>& contractions) -> std::vector<Contraction> {
  check_input(l, exponents, contractions);
  const std::size_t n = exponents.size();
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::ranges::stable_sort(
      order, [&](std::size_t a, std::size_t b) { return exponents[a] > exponents[b]; });
  std::vector<double> sorted(n);
  for (std::size_t i = 0; i < n; ++i) {
    sorted[i] = exponents[order[i]];
  }
  if (std::ranges::adjacent_find(sorted) != sorted.end()) {
    fail("duplicate exponents");
  }
  exponents = std::move(sorted);

  return dispatch_angular_momentum(l, [&]<int L>(std::integral_constant<int, L>) {
    const std::vector<double> overlaps = normalized_primitive_overlaps<L>(exponents);
    std::vector<Contraction> result;
    result.reserve(contractions.size());
    for (const auto& raw : contractions) {
      Contraction contraction;
      for (std::size_t i = 0; i < n; ++i) {
        if (const double c = raw[order[i]]; c != 0.0) {
          contraction.indices.push_back(static_cast<std::uint16_t>(i));
          contraction.coefficients.push_back(c);
        }
      }
      double self_overlap = 0.0;
      for (std::size_t p = 0; p < contraction.indices.size(); ++p) {
        for (std::size_t q = 0; q < contraction.indices.size(); ++q) {
          self_overlap += contraction.coefficients[p] * contraction.coefficients[q] *
                          overlaps[contraction.indices[p] * n + contraction.indices[q]];
        }
      }
      const double scale = 1.0 / std::sqrt(self_overlap);
      for (std::size_t p = 0; p < contraction.indices.size(); ++p) {
        contraction.coefficients[p] *=
            scale * primitive_normalization<L>(exponents[contraction.indices[p]]);
        contraction.coefficient_sum += std::abs(contraction.coefficients[p]);
      }
      result.push_back(std::move(contraction));
    }
    return result;
  });
}

auto nonzero_count(const std::vector<double>& coefficients) -> std::size_t {
  return static_cast<std::size_t>(
      std::ranges::count_if(coefficients, [](double c) { return c != 0.0; }));
}

}  // namespace

UncontractedShell::UncontractedShell(int angular_momentum, double exponent)
    : angular_momentum_(angular_momentum), exponent_(exponent) {
  check_angular_momentum(angular_momentum);
  check_exponent(exponent);
  coefficient_ = dispatch_angular_momentum(
      angular_momentum,
      [&]<int L>(std::integral_constant<int, L>) { return primitive_normalization<L>(exponent); });
}

SegmentedShell::SegmentedShell(int angular_momentum, std::vector<double> exponents,
                               const std::vector<double>& coefficients)
    : angular_momentum_(angular_momentum) {
  if (nonzero_count(coefficients) < 2) {
    fail("a segmented shell needs at least two nonzero coefficients");
  }
  Contraction contraction = std::move(prepare(angular_momentum, exponents, {coefficients}).front());
  for (const std::uint16_t index : contraction.indices) {
    exponents_.push_back(exponents[index]);
  }
  coefficients_ = std::move(contraction.coefficients);
  coefficient_sum_ = contraction.coefficient_sum;
}

GeneralShell::GeneralShell(int angular_momentum, std::vector<double> exponents,
                           const std::vector<std::vector<double>>& contractions)
    : angular_momentum_(angular_momentum) {
  const std::vector<Contraction> prepared = prepare(angular_momentum, exponents, contractions);
  exponents_ = std::move(exponents);
  offsets_.reserve(prepared.size() + 1);
  for (const Contraction& contraction : prepared) {
    primitive_indices_.insert(primitive_indices_.end(), contraction.indices.begin(),
                              contraction.indices.end());
    coefficients_.insert(coefficients_.end(), contraction.coefficients.begin(),
                         contraction.coefficients.end());
    offsets_.push_back(coefficients_.size());
    max_coefficient_sum_ = std::max(max_coefficient_sum_, contraction.coefficient_sum);
  }
}

auto make_basis_shell(int angular_momentum, std::vector<double> exponents,
                      const std::vector<std::vector<double>>& contractions) -> BasisShell {
  check_input(angular_momentum, exponents, contractions);
  if (contractions.size() > 1) {
    return GeneralShell(angular_momentum, std::move(exponents), contractions);
  }
  const std::vector<double>& coefficients = contractions.front();
  if (nonzero_count(coefficients) == 1) {
    const auto nonzero = std::ranges::find_if(coefficients, [](double c) { return c != 0.0; });
    return UncontractedShell(angular_momentum,
                             exponents[static_cast<std::size_t>(nonzero - coefficients.begin())]);
  }
  return SegmentedShell(angular_momentum, std::move(exponents), coefficients);
}

auto kind(const BasisShell& shell) noexcept -> ShellKind {
  return static_cast<ShellKind>(shell.index());
}

auto angular_momentum(const BasisShell& shell) noexcept -> int {
  return std::visit([](const auto& s) { return s.angular_momentum(); }, shell);
}

auto contraction_count(const BasisShell& shell) noexcept -> std::size_t {
  return std::visit([](const auto& s) { return s.contraction_count(); }, shell);
}

auto function_count(const BasisShell& shell) noexcept -> std::size_t {
  return static_cast<std::size_t>(2 * angular_momentum(shell) + 1) * contraction_count(shell);
}

auto exponents(const BasisShell& shell) noexcept -> std::span<const double> {
  return std::visit([](const auto& s) { return s.exponents(); }, shell);
}

auto max_coefficient_sum(const BasisShell& shell) noexcept -> double {
  return std::visit(
      [](const auto& s) -> double {
        if constexpr (std::is_same_v<std::decay_t<decltype(s)>, UncontractedShell>) {
          return s.coefficient();
        } else {
          return s.max_coefficient_sum();
        }
      },
      shell);
}

auto min_exponent(const BasisShell& shell) noexcept -> double {
  return exponents(shell).back();
}

auto max_exponent(const BasisShell& shell) noexcept -> double {
  return exponents(shell).front();
}

}  // namespace fika
