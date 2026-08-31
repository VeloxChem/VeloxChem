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



#include "SimdHarmonicsRecD.hpp"

#include <cmath>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_d_solid_harmonics(const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_d_solid_harmonics: Harmonics of angular momentum one must have 3 rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_d_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(p_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_d_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(5, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = std::sqrt(0.75);
    const auto fact_xy = std::sqrt(3.0);
    const auto fz_m1 = std::sqrt(3.0);
    const auto fz_0 = 1.5;
    const auto fr_0 = 0.5;
    const auto fz_p1 = std::sqrt(3.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    auto *ps_m2 = matrix.data(0);
    auto *ps_m1 = matrix.data(1);
    auto *ps_0 = matrix.data(2);
    auto *ps_p1 = matrix.data(3);
    auto *ps_p2 = matrix.data(4);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, ps_m2, ps_m1, ps_0, ps_p1, ps_p2 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        ps_p2[i] = fact * (x * x - y * y);

        ps_m2[i] = fact_xy * x * y;

        ps_m1[i] = fz_m1 * z * y;

        ps_0[i] = fz_0 * z * z - fr_0 * r_2;

        ps_p1[i] = fz_p1 * z * x;
    }

    return matrix;
}

}  // namespace simdfunc
