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

#ifndef fika_gaunt_hpp
#define fika_gaunt_hpp

#include <cstddef>
#include <span>

namespace fika {

/// Real Gaunt coefficient C^{LM}_{lm,l'm'} (Racah normalization): the coefficient of
/// r^(2k) S_{LM}(r), k = (l + l' - L) / 2, in the product S_{lm}(r) S_{l'm'}(r).
struct GauntEntry {
  int m;
  int m_prime;
  int big_l;
  int big_m;
  double value;
};

/// Nonzero real Gaunt coefficients for l, l' <= max_angular_momentum, all L (|l - l'|..l + l',
/// same parity), ordered by (m, m') and then L. Each pair (l, l') is computed on its first
/// request (thread-safe) by exact quadrature over the sphere with fika's solid harmonics (so
/// conventions match everywhere).
auto gaunt_coefficients(int l, int l_prime) -> std::span<const GauntEntry>;

/// Highest rank of an MM multipole source (quadrupoles).
inline constexpr int max_multipole_rank = 2;

/// Gaunt coefficients C^{LM}_{km,lm'} of a multipole source of rank k (1..max_multipole_rank)
/// with a harmonic of degree l <= 2 max_angular_momentum (the orders of a shell pair), in the
/// layout of gaunt_coefficients (m of rank k, m' of degree l). Built on first request.
auto multipole_gaunt_coefficients(int rank, int l) -> std::span<const GauntEntry>;

}  // namespace fika

#endif  // fika_gaunt_hpp
