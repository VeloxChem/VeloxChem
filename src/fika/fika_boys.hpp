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

#ifndef fika_boys_hpp
#define fika_boys_hpp

// Internal: Boys functions F_k(x) = int_0^1 t^(2k) exp(-x t^2) dt.

#include <span>

namespace fika::detail {

/// Highest Boys order boys() evaluates.
inline constexpr int max_boys_order = 32;

/// F_0(x)..F_k_max(x) into out[0..k_max] (0 <= k_max <= max_boys_order, x >= 0). Method of
/// Vikhamar-Sandberg and Repisky (arXiv:2512.10059): on [0, x0) the rational minimax
/// approximation of F_k_max and downward recursion, on [x0, x1) that of F_0 and upward recursion.
/// Beyond x1 (where the paper drops exp(-x) and only bounds the absolute error) F_0 =
/// sqrt(pi) erf(sqrt(x)) / (2 sqrt(x)) and the exact upward recursion keep the relative error
/// small for every order. Absolute error <= 5e-14; relative error <= ~3e-13 for k <= 16.
void boys(int k_max, double x, std::span<double> out);

}  // namespace fika::detail

#endif  // fika_boys_hpp
