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
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef newints_BoysFunction_hpp
#define newints_BoysFunction_hpp

#include <array>
#include <cmath>
#include <cstddef>

#include "BoysMinimaxCoefficients.hpp"

namespace newints {

namespace boys_detail {

/// @brief Evaluates a polynomial by Horner's scheme from coefficients stored
///        ascending in degree (constant term first).
/// @param coeffs The polynomial coefficients (c[0] + c[1] x + ... + c[M-1] x^{M-1}).
/// @param x The argument.
/// @return The polynomial value at x.
template <std::size_t M>
inline auto
horner(const std::array<double, M>& coeffs, const double x) -> double
{
    double acc = coeffs[M - 1];
    for (std::size_t i = M - 1; i-- > 0;)
    {
        acc = acc * x + coeffs[i];
    }
    return acc;
}

}  // namespace boys_detail

/// @brief Evaluates the Boys functions F_0(x), ..., F_{N-1}(x) at a single argument.
///
/// Implements the three-region minimax scheme of R. Vikhamar-Sandberg and
/// M. Repisky, arXiv:2512.10059v3 (2025): the non-negative axis is split at
/// region_b_start (x0) and region_c_start (x1) into A = [0, x0), B = [x0, x1)
/// and C = [x1, inf). Region A evaluates the top order F_{N-1} from its rational
/// minimax and recurses downward (stable); region B evaluates F_0 from its
/// rational minimax and recurses upward; region C uses the asymptotic F_0 and
/// recurses upward. The maximum absolute error is bounded by ~5e-14.
///
/// @tparam N The number of Boys function values to compute (orders 0 .. N-1).
/// @param x The Boys function argument (x >= 0).
/// @param values On return, values[i] = F_i(x) for i = 0, ..., N-1.
template <std::size_t N>
auto
boys_function(const double x, std::array<double, N>& values) -> void
{
    static_assert(N >= 1, "boys_function requires at least one value (N >= 1).");
    static_assert(N <= static_cast<std::size_t>(boys_detail::max_order) + 1,
                  "boys_function supports orders up to F32 (N <= 33).");

    constexpr int k = static_cast<int>(N) - 1;

    if (x < boys_detail::region_b_start)
    {
        // Region A: top order from its rational minimax, then downward recursion.
        const double ex = 0.5 * std::exp(-x);

        double F = boys_detail::horner(boys_detail::region_a<k>::num, x) /
                   boys_detail::horner(boys_detail::region_a<k>::den, x);
        values[k] = F;

        for (int l = k - 1; l >= 0; --l)
        {
            F         = (x * F + ex) / (static_cast<double>(l) + 0.5);
            values[l] = F;
        }
    }
    else if (x < boys_detail::region_c_start)
    {
        // Region B: F_0 from its rational minimax, then upward recursion.
        const double ex = 0.5 * std::exp(-x);

        double F  = boys_detail::horner(boys_detail::region_b_num, x) /
                   boys_detail::horner(boys_detail::region_b_den, x);
        values[0] = F;

        for (int l = 1; l <= k; ++l)
        {
            F         = ((static_cast<double>(l) - 0.5) * F - ex) / x;
            values[l] = F;
        }
    }
    else
    {
        // Region C: asymptotic F_0, then upward recursion (e^{-x} negligible).
        double F  = boys_detail::half_sqrt_pi / std::sqrt(x);
        values[0] = F;

        for (int l = 1; l <= k; ++l)
        {
            F         = (static_cast<double>(l) - 0.5) * F / x;
            values[l] = F;
        }
    }
}

}  // namespace newints

#endif /* newints_BoysFunction_hpp */
