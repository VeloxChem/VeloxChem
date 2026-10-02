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

#ifndef fika_basis_shell_hpp
#define fika_basis_shell_hpp

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

namespace fika {

// Shells of real solid-harmonic Gaussians sharing one angular momentum. Coefficients are
// effective: they include the primitive normalization N_l(a) and the normalization of each
// contracted function. Exponents are stored in decreasing order.

/// A single normalized primitive: one contracted function of one primitive.
class UncontractedShell {
 public:
  /// Throws std::invalid_argument for invalid angular momentum or exponent.
  UncontractedShell(int angular_momentum, double exponent);

  auto angular_momentum() const noexcept -> int { return angular_momentum_; }
  auto primitive_count() const noexcept -> std::size_t { return 1; }
  auto contraction_count() const noexcept -> std::size_t { return 1; }
  auto exponent() const noexcept -> double { return exponent_; }
  /// Effective coefficient N_l(exponent); also the shell's sum of |coefficient|.
  auto coefficient() const noexcept -> double { return coefficient_; }

  auto exponents() const noexcept -> std::span<const double> { return {&exponent_, 1}; }
  auto coefficients() const noexcept -> std::span<const double> { return {&coefficient_, 1}; }

 private:
  int angular_momentum_;
  double exponent_;
  double coefficient_;
};

/// One contracted function of several primitives.
class SegmentedShell {
 public:
  /// coefficients[i]: coefficient of normalized primitive i (exponents[i]) as listed in
  /// basis-set files; zero coefficients are dropped. Throws std::invalid_argument for invalid
  /// input or fewer than two nonzero coefficients.
  SegmentedShell(int angular_momentum, std::vector<double> exponents,
                 const std::vector<double>& coefficients);

  auto angular_momentum() const noexcept -> int { return angular_momentum_; }
  auto primitive_count() const noexcept -> std::size_t { return exponents_.size(); }
  auto contraction_count() const noexcept -> std::size_t { return 1; }
  auto exponents() const noexcept -> std::span<const double> { return exponents_; }
  auto coefficients() const noexcept -> std::span<const double> { return coefficients_; }
  auto max_coefficient_sum() const noexcept -> double { return coefficient_sum_; }

 private:
  int angular_momentum_;
  std::vector<double> exponents_;
  std::vector<double> coefficients_;
  double coefficient_sum_;
};

/// Several contracted functions sharing primitives (general contraction), stored sparsely.
class GeneralShell {
 public:
  /// contractions[k][i]: coefficient of normalized primitive i (exponents[i]) in contracted
  /// function k, as listed in basis-set files; zero coefficients are allowed and dropped.
  /// Throws std::invalid_argument for invalid angular momentum, exponents or coefficients.
  GeneralShell(int angular_momentum, std::vector<double> exponents,
               const std::vector<std::vector<double>>& contractions);

  auto angular_momentum() const noexcept -> int { return angular_momentum_; }
  auto primitive_count() const noexcept -> std::size_t { return exponents_.size(); }
  auto contraction_count() const noexcept -> std::size_t { return offsets_.size() - 1; }

  auto exponents() const noexcept -> std::span<const double> { return exponents_; }

  /// Nonzero effective coefficients of contracted function k.
  auto coefficients(std::size_t k) const noexcept -> std::span<const double> {
    assert(k < contraction_count());
    return std::span(coefficients_).subspan(offsets_[k], offsets_[k + 1] - offsets_[k]);
  }

  /// Primitive indices (into exponents()) matching coefficients(k), increasing.
  auto primitive_indices(std::size_t k) const noexcept -> std::span<const std::uint16_t> {
    assert(k < contraction_count());
    return std::span(primitive_indices_).subspan(offsets_[k], offsets_[k + 1] - offsets_[k]);
  }

  /// Largest sum of |coefficient| over the contracted functions.
  auto max_coefficient_sum() const noexcept -> double { return max_coefficient_sum_; }

 private:
  int angular_momentum_;
  std::vector<double> exponents_;
  std::vector<std::size_t> offsets_{0};           // contraction_count() + 1 entries
  std::vector<std::uint16_t> primitive_indices_;  // nonzeros, grouped by contracted function
  std::vector<double> coefficients_;              // effective coefficients, same layout
  double max_coefficient_sum_ = 0.0;
};

using BasisShell = std::variant<UncontractedShell, SegmentedShell, GeneralShell>;

enum class ShellKind { uncontracted, segmented, general };

/// Shell from basis-file data (contractions[k][i] as for GeneralShell): one contracted function
/// with one nonzero coefficient is uncontracted, one contracted function is segmented, several
/// are a general contraction (kept as given). Throws std::invalid_argument for invalid input.
auto make_basis_shell(int angular_momentum, std::vector<double> exponents,
                      const std::vector<std::vector<double>>& contractions) -> BasisShell;

auto kind(const BasisShell& shell) noexcept -> ShellKind;
auto angular_momentum(const BasisShell& shell) noexcept -> int;
auto contraction_count(const BasisShell& shell) noexcept -> std::size_t;
/// Basis functions: (2l + 1) times contracted functions.
auto function_count(const BasisShell& shell) noexcept -> std::size_t;
auto exponents(const BasisShell& shell) noexcept -> std::span<const double>;
auto max_coefficient_sum(const BasisShell& shell) noexcept -> double;
auto min_exponent(const BasisShell& shell) noexcept -> double;
auto max_exponent(const BasisShell& shell) noexcept -> double;

/// Calls f(primitive index, effective coefficient) for the nonzero coefficients of contracted
/// function k; the index refers to the shell's exponents().
template <typename F>
void for_each_coefficient(const UncontractedShell& shell, [[maybe_unused]] std::size_t k, F&& f) {
  assert(k == 0);
  f(std::size_t{0}, shell.coefficient());
}

template <typename F>
void for_each_coefficient(const SegmentedShell& shell, [[maybe_unused]] std::size_t k, F&& f) {
  assert(k == 0);
  const auto coefficients = shell.coefficients();
  for (std::size_t i = 0; i < coefficients.size(); ++i) {
    f(i, coefficients[i]);
  }
}

template <typename F>
void for_each_coefficient(const GeneralShell& shell, std::size_t k, F&& f) {
  const auto coefficients = shell.coefficients(k);
  const auto indices = shell.primitive_indices(k);
  for (std::size_t p = 0; p < coefficients.size(); ++p) {
    f(std::size_t{indices[p]}, coefficients[p]);
  }
}

}  // namespace fika

#endif  // fika_basis_shell_hpp
