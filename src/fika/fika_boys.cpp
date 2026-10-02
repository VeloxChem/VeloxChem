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

#include "fika_boys.hpp"

#include <cassert>
#include <cmath>
#include <numbers>

#include "fika_boys_coefficients.hpp"

namespace fika::detail {

namespace {

auto evaluate(const boys_minimax::Rational& r, double x) -> double {
  double numerator = 0.0;
  for (auto c = r.numerator.rbegin(); c != r.numerator.rend(); ++c) {
    numerator = numerator * x + *c;
  }
  double denominator = 0.0;
  for (auto c = r.denominator.rbegin(); c != r.denominator.rend(); ++c) {
    denominator = denominator * x + *c;
  }
  return numerator / denominator;
}

}  // namespace

void boys(int k_max, double x, std::span<double> out) {
  assert(k_max >= 0 && k_max <= max_boys_order && x >= 0.0);
  assert(out.size() > static_cast<std::size_t>(k_max));
  const auto k = static_cast<std::size_t>(k_max);
  if (x < boys_minimax::x0) {
    // F_k_max from its approximation, then F_l = (2x F_(l+1) + exp(-x)) / (2l + 1).
    out[k] = evaluate(boys_minimax::region_a[k], x);
    const double exponential = std::exp(-x);
    for (std::size_t l = k; l-- > 0;) {
      out[l] = (2.0 * x * out[l + 1] + exponential) / static_cast<double>(2 * l + 1);
    }
    return;
  }
  // F_0, then F_(l+1) = ((2l + 1) F_l - exp(-x)) / (2x), stable for x >= x0.
  const double root = std::sqrt(x);
  out[0] = x < boys_minimax::x1 ? evaluate(boys_minimax::region_b, x)
                                : std::sqrt(std::numbers::pi) * std::erf(root) / (2.0 * root);
  const double exponential = std::exp(-x);
  const double inverse = 1.0 / (2.0 * x);
  for (std::size_t l = 0; l < k; ++l) {
    out[l + 1] = (static_cast<double>(2 * l + 1) * out[l] - exponential) * inverse;
  }
}

}  // namespace fika::detail
