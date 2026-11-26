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

#ifndef fika_gaussian_normalization_hpp
#define fika_gaussian_normalization_hpp

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace fika {

/// Highest angular momentum supported by the Gaussian kernels (l = 0..8).
inline constexpr int max_angular_momentum = 8;

template <int L>
concept supported_angular_momentum = L >= 0 && L <= max_angular_momentum;

/// Calls f(std::integral_constant<int, L>{}) with L equal to the runtime angular momentum l,
/// so kernels templated on L can be selected once per shell. Throws std::invalid_argument for l
/// outside 0..max_angular_momentum.
template <typename F>
auto dispatch_angular_momentum(int l, F&& f) -> decltype(auto) {
  switch (l) {
    case 0:
      return f(std::integral_constant<int, 0>{});
    case 1:
      return f(std::integral_constant<int, 1>{});
    case 2:
      return f(std::integral_constant<int, 2>{});
    case 3:
      return f(std::integral_constant<int, 3>{});
    case 4:
      return f(std::integral_constant<int, 4>{});
    case 5:
      return f(std::integral_constant<int, 5>{});
    case 6:
      return f(std::integral_constant<int, 6>{});
    case 7:
      return f(std::integral_constant<int, 7>{});
    case 8:
      return f(std::integral_constant<int, 8>{});
    default:
      throw std::invalid_argument("fika: unsupported angular momentum " + std::to_string(l));
  }
}
static_assert(max_angular_momentum == 8, "extend dispatch_angular_momentum");

/// x^N by repeated squaring, unrolled at compile time.
template <int N>
  requires(N >= 0)
constexpr auto integer_power(double x) noexcept -> double {
  if constexpr (N == 0) {
    return 1.0;
  } else if constexpr (N == 1) {
    return x;
  } else if constexpr (N % 2 == 0) {
    const double half = integer_power<N / 2>(x);
    return half * half;
  } else {
    return x * integer_power<N - 1>(x);
  }
}

namespace detail {

/// n!! with (-1)!! = 0!! = 1.
constexpr auto double_factorial(int n) noexcept -> double {
  double result = 1.0;
  for (int k = n; k > 1; k -= 2) {
    result *= k;
  }
  return result;
}

/// pi^(3/2).
inline constexpr double pi_three_halves = 5.5683279968317078452848179821188357;

/// 4^L / (2L-1)!!
template <int L>
inline constexpr double normalization_constant =
    integer_power<L>(4.0) / double_factorial(2 * L - 1);

/// (2L-1)!! 2^-L pi^(3/2)
template <int L>
inline constexpr double overlap_constant =
    double_factorial(2 * L - 1) / integer_power<L>(2.0) * pi_three_halves;

}  // namespace detail

/// Normalization factor of a primitive real solid-harmonic Gaussian with exponent a:
/// N_L(a) = sqrt((4a)^L / (2L-1)!!) (2a/pi)^(3/4), evaluated as sqrt(c_L a^L f sqrt(f)), f = 2a/pi.
template <int L>
  requires supported_angular_momentum<L>
inline auto primitive_normalization(double a) noexcept -> double {
  const double f = 2.0 * std::numbers::inv_pi * a;
  return std::sqrt(detail::normalization_constant<L> * integer_power<L>(a) * f * std::sqrt(f));
}

/// Same-centre overlap of unnormalized primitives with equal L and m:
/// S_LL(a, b) = (2L-1)!! / (2p)^L (pi/p)^(3/2), evaluated as K_L / (p^(L+1) sqrt(p)), p = a + b.
/// Dividing once at the end avoids amplifying the rounding of 1/p by the power L + 1.
template <int L>
  requires supported_angular_momentum<L>
inline auto primitive_overlap(double a, double b) noexcept -> double {
  const double p = a + b;
  return detail::overlap_constant<L> / (integer_power<L + 1>(p) * std::sqrt(p));
}

}  // namespace fika

#endif  // fika_gaussian_normalization_hpp
