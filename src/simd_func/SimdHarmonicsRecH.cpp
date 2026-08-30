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



#include "SimdHarmonicsRecH.hpp"

#include <cmath>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_h_solid_harmonics(const CSimdMatrix &g_harmonics, const CSimdMatrix &f_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(g_harmonics.number_of_rows() == 9,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum four must have 9 rows"));

    errors::assertMsgCritical(f_harmonics.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum three must have 7 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_h_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(g_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_h_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(11, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = std::sqrt(9.0 / 10.0);
    const auto fz_m4 = 3.0;
    const auto fz_m3 = 2.25;
    const auto fr_m3 = std::sqrt(0.4375);
    const auto fz_m2 = std::sqrt(27.0 / 7.0);
    const auto fr_m2 = std::sqrt(4.0 / 7.0);
    const auto fz_m1 = std::sqrt(3.375);
    const auto fr_m1 = std::sqrt(0.625);
    const auto fz_0 = 9.0 / 5.0;
    const auto fr_0 = 4.0 / 5.0;
    const auto fz_p1 = std::sqrt(3.375);
    const auto fr_p1 = std::sqrt(0.625);
    const auto fz_p2 = std::sqrt(27.0 / 7.0);
    const auto fr_p2 = std::sqrt(4.0 / 7.0);
    const auto fz_p3 = 2.25;
    const auto fr_p3 = std::sqrt(0.4375);
    const auto fz_p4 = 3.0;

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m4 = g_harmonics.data(0);
    const auto *pt_m3 = g_harmonics.data(1);
    const auto *pt_m2 = g_harmonics.data(2);
    const auto *pt_m1 = g_harmonics.data(3);
    const auto *pt_0 = g_harmonics.data(4);
    const auto *pt_p1 = g_harmonics.data(5);
    const auto *pt_p2 = g_harmonics.data(6);
    const auto *pt_p3 = g_harmonics.data(7);
    const auto *pt_p4 = g_harmonics.data(8);

    const auto *pu_m3 = f_harmonics.data(0);
    const auto *pu_m2 = f_harmonics.data(1);
    const auto *pu_m1 = f_harmonics.data(2);
    const auto *pu_0 = f_harmonics.data(3);
    const auto *pu_p1 = f_harmonics.data(4);
    const auto *pu_p2 = f_harmonics.data(5);
    const auto *pu_p3 = f_harmonics.data(6);

    auto *ps_m5 = matrix.data(0);
    auto *ps_m4 = matrix.data(1);
    auto *ps_m3 = matrix.data(2);
    auto *ps_m2 = matrix.data(3);
    auto *ps_m1 = matrix.data(4);
    auto *ps_0 = matrix.data(5);
    auto *ps_p1 = matrix.data(6);
    auto *ps_p2 = matrix.data(7);
    auto *ps_p3 = matrix.data(8);
    auto *ps_p4 = matrix.data(9);
    auto *ps_p5 = matrix.data(10);

    // NOTE: every value entering the recursion is loaded before any result is
    // stored, so that the rows of the previous orders are read once and the
    // stores cannot force the loads to be repeated.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m4, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, ps_m5, ps_m4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4, ps_p5 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m4 = pt_m4[i];
        const auto t_m3 = pt_m3[i];
        const auto t_m2 = pt_m2[i];
        const auto t_m1 = pt_m1[i];
        const auto t_0 = pt_0[i];
        const auto t_p1 = pt_p1[i];
        const auto t_p2 = pt_p2[i];
        const auto t_p3 = pt_p3[i];
        const auto t_p4 = pt_p4[i];

        const auto u_m3 = pu_m3[i];
        const auto u_m2 = pu_m2[i];
        const auto u_m1 = pu_m1[i];
        const auto u_0 = pu_0[i];
        const auto u_p1 = pu_p1[i];
        const auto u_p2 = pu_p2[i];
        const auto u_p3 = pu_p3[i];

        ps_p5[i] = fact * (x * t_p4 - y * t_m4);

        ps_m5[i] = fact * (y * t_p4 + x * t_m4);

        ps_m4[i] = fz_m4 * z * t_m4;

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4;
    }

    return matrix;
}

}  // namespace simdfunc
