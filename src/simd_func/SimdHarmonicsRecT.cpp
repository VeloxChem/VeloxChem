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



#include "SimdHarmonicsRecT.hpp"

#include <cmath>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

namespace simdfunc {  // simdfunc namespace

auto
make_t_solid_harmonics(const CSimdMatrix &r_harmonics, const CSimdMatrix &q_harmonics, const CSimdMatrix &p_harmonics, const CSimdMatrix &coordinates) -> CSimdMatrix
{
    errors::assertMsgCritical(r_harmonics.number_of_rows() == 27,
                              std::string("SimdHarmonics.make_t_solid_harmonics: Harmonics of angular momentum thirteen must have 27 rows"));

    errors::assertMsgCritical(q_harmonics.number_of_rows() == 25,
                              std::string("SimdHarmonics.make_t_solid_harmonics: Harmonics of angular momentum twelve must have 25 rows"));

    errors::assertMsgCritical(p_harmonics.number_of_rows() == 3,
                              std::string("SimdHarmonics.make_t_solid_harmonics: Harmonics of angular momentum one must have three rows"));

    errors::assertMsgCritical(coordinates.number_of_rows() == 7,
                              std::string("SimdHarmonics.make_t_solid_harmonics: Coordinates must have seven rows"));

    errors::assertMsgCritical(r_harmonics.number_of_columns() == coordinates.number_of_columns(),
                              std::string("SimdHarmonics.make_t_solid_harmonics: Harmonics and coordinates must have the same number of columns"));

    const auto ncols = coordinates.number_of_columns();

    auto matrix = CSimdMatrix(29, ncols);

    // NOTE: the factors of the recursion depend on the angular momentum alone, so
    // they are formed once for the whole matrix instead of once for every atom pair.

    const auto fact = std::sqrt(27.0 / 28.0);
    const auto fz_m13 = std::sqrt(27.0);
    const auto fz_m12 = std::sqrt(729.0 / 52.0);
    const auto fr_m12 = std::sqrt(25.0 / 52.0);
    const auto fz_m11 = std::sqrt(243.0 / 25.0);
    const auto fr_m11 = 4.0 / 5.0;
    const auto fz_m10 = std::sqrt(7.59375);
    const auto fr_m10 = std::sqrt(0.71875);
    const auto fz_m9 = std::sqrt(729.0 / 115.0);
    const auto fr_m9 = std::sqrt(88.0 / 115.0);
    const auto fz_m8 = std::sqrt(243.0 / 44.0);
    const auto fr_m8 = std::sqrt(35.0 / 44.0);
    const auto fz_m7 = std::sqrt(243.0 / 49.0);
    const auto fr_m7 = std::sqrt(40.0 / 49.0);
    const auto fz_m6 = std::sqrt(729.0 / 160.0);
    const auto fr_m6 = std::sqrt(133.0 / 160.0);
    const auto fz_m5 = std::sqrt(81.0 / 19.0);
    const auto fr_m5 = std::sqrt(16.0 / 19.0);
    const auto fz_m4 = std::sqrt(81.0 / 20.0);
    const auto fr_m4 = std::sqrt(17.0 / 20.0);
    const auto fz_m3 = std::sqrt(729.0 / 187.0);
    const auto fr_m3 = std::sqrt(160.0 / 187.0);
    const auto fz_m2 = std::sqrt(3.796875);
    const auto fr_m2 = std::sqrt(0.859375);
    const auto fz_m1 = std::sqrt(243.0 / 65.0);
    const auto fr_m1 = std::sqrt(56.0 / 65.0);
    const auto fz_0 = 27.0 / 14.0;
    const auto fr_0 = 13.0 / 14.0;
    const auto fz_p1 = std::sqrt(243.0 / 65.0);
    const auto fr_p1 = std::sqrt(56.0 / 65.0);
    const auto fz_p2 = std::sqrt(3.796875);
    const auto fr_p2 = std::sqrt(0.859375);
    const auto fz_p3 = std::sqrt(729.0 / 187.0);
    const auto fr_p3 = std::sqrt(160.0 / 187.0);
    const auto fz_p4 = std::sqrt(81.0 / 20.0);
    const auto fr_p4 = std::sqrt(17.0 / 20.0);
    const auto fz_p5 = std::sqrt(81.0 / 19.0);
    const auto fr_p5 = std::sqrt(16.0 / 19.0);
    const auto fz_p6 = std::sqrt(729.0 / 160.0);
    const auto fr_p6 = std::sqrt(133.0 / 160.0);
    const auto fz_p7 = std::sqrt(243.0 / 49.0);
    const auto fr_p7 = std::sqrt(40.0 / 49.0);
    const auto fz_p8 = std::sqrt(243.0 / 44.0);
    const auto fr_p8 = std::sqrt(35.0 / 44.0);
    const auto fz_p9 = std::sqrt(729.0 / 115.0);
    const auto fr_p9 = std::sqrt(88.0 / 115.0);
    const auto fz_p10 = std::sqrt(7.59375);
    const auto fr_p10 = std::sqrt(0.71875);
    const auto fz_p11 = std::sqrt(243.0 / 25.0);
    const auto fr_p11 = 4.0 / 5.0;
    const auto fz_p12 = std::sqrt(729.0 / 52.0);
    const auto fr_p12 = std::sqrt(25.0 / 52.0);
    const auto fz_p13 = std::sqrt(27.0);

    // NOTE: the pointers to the rows are taken once, as the accessor of a row is
    // bounds checked and would otherwise be called for every atom pair.

    const auto *pr_y = p_harmonics.data(0);
    const auto *pr_z = p_harmonics.data(1);
    const auto *pr_x = p_harmonics.data(2);

    const auto *pr_2 = coordinates.data(6);

    const auto *pt_m13 = r_harmonics.data(0);
    const auto *pt_m12 = r_harmonics.data(1);
    const auto *pt_m11 = r_harmonics.data(2);
    const auto *pt_m10 = r_harmonics.data(3);
    const auto *pt_m9 = r_harmonics.data(4);
    const auto *pt_m8 = r_harmonics.data(5);
    const auto *pt_m7 = r_harmonics.data(6);
    const auto *pt_m6 = r_harmonics.data(7);
    const auto *pt_m5 = r_harmonics.data(8);
    const auto *pt_m4 = r_harmonics.data(9);
    const auto *pt_m3 = r_harmonics.data(10);
    const auto *pt_m2 = r_harmonics.data(11);
    const auto *pt_m1 = r_harmonics.data(12);
    const auto *pt_0 = r_harmonics.data(13);
    const auto *pt_p1 = r_harmonics.data(14);
    const auto *pt_p2 = r_harmonics.data(15);
    const auto *pt_p3 = r_harmonics.data(16);
    const auto *pt_p4 = r_harmonics.data(17);
    const auto *pt_p5 = r_harmonics.data(18);
    const auto *pt_p6 = r_harmonics.data(19);
    const auto *pt_p7 = r_harmonics.data(20);
    const auto *pt_p8 = r_harmonics.data(21);
    const auto *pt_p9 = r_harmonics.data(22);
    const auto *pt_p10 = r_harmonics.data(23);
    const auto *pt_p11 = r_harmonics.data(24);
    const auto *pt_p12 = r_harmonics.data(25);
    const auto *pt_p13 = r_harmonics.data(26);

    const auto *pu_m12 = q_harmonics.data(0);
    const auto *pu_m11 = q_harmonics.data(1);
    const auto *pu_m10 = q_harmonics.data(2);
    const auto *pu_m9 = q_harmonics.data(3);
    const auto *pu_m8 = q_harmonics.data(4);
    const auto *pu_m7 = q_harmonics.data(5);
    const auto *pu_m6 = q_harmonics.data(6);
    const auto *pu_m5 = q_harmonics.data(7);
    const auto *pu_m4 = q_harmonics.data(8);
    const auto *pu_m3 = q_harmonics.data(9);
    const auto *pu_m2 = q_harmonics.data(10);
    const auto *pu_m1 = q_harmonics.data(11);
    const auto *pu_0 = q_harmonics.data(12);
    const auto *pu_p1 = q_harmonics.data(13);
    const auto *pu_p2 = q_harmonics.data(14);
    const auto *pu_p3 = q_harmonics.data(15);
    const auto *pu_p4 = q_harmonics.data(16);
    const auto *pu_p5 = q_harmonics.data(17);
    const auto *pu_p6 = q_harmonics.data(18);
    const auto *pu_p7 = q_harmonics.data(19);
    const auto *pu_p8 = q_harmonics.data(20);
    const auto *pu_p9 = q_harmonics.data(21);
    const auto *pu_p10 = q_harmonics.data(22);
    const auto *pu_p11 = q_harmonics.data(23);
    const auto *pu_p12 = q_harmonics.data(24);

    auto *ps_m14 = matrix.data(0);
    auto *ps_m13 = matrix.data(1);
    auto *ps_m12 = matrix.data(2);
    auto *ps_m11 = matrix.data(3);
    auto *ps_m10 = matrix.data(4);
    auto *ps_m9 = matrix.data(5);
    auto *ps_m8 = matrix.data(6);
    auto *ps_m7 = matrix.data(7);
    auto *ps_m6 = matrix.data(8);
    auto *ps_m5 = matrix.data(9);
    auto *ps_m4 = matrix.data(10);
    auto *ps_m3 = matrix.data(11);
    auto *ps_m2 = matrix.data(12);
    auto *ps_m1 = matrix.data(13);
    auto *ps_0 = matrix.data(14);
    auto *ps_p1 = matrix.data(15);
    auto *ps_p2 = matrix.data(16);
    auto *ps_p3 = matrix.data(17);
    auto *ps_p4 = matrix.data(18);
    auto *ps_p5 = matrix.data(19);
    auto *ps_p6 = matrix.data(20);
    auto *ps_p7 = matrix.data(21);
    auto *ps_p8 = matrix.data(22);
    auto *ps_p9 = matrix.data(23);
    auto *ps_p10 = matrix.data(24);
    auto *ps_p11 = matrix.data(25);
    auto *ps_p12 = matrix.data(26);
    auto *ps_p13 = matrix.data(27);
    auto *ps_p14 = matrix.data(28);

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 29 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_x, pr_y, pr_z, pr_2, pt_m13, pt_m12, pt_m11, pt_m10, pt_m9, pt_m8, pt_m7, pt_m6, pt_m5, pt_m4, pt_p13, pu_m12, pu_m11, pu_m10, pu_m9, pu_m8, pu_m7, pu_m6, pu_m5, pu_m4, ps_m14, ps_m13, ps_m12, ps_m11, ps_m10, ps_m9, ps_m8, ps_m7, ps_m6, ps_m5, ps_m4, ps_p14 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto x = pr_x[i];
        const auto y = pr_y[i];
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_m13 = pt_m13[i];
        const auto t_m12 = pt_m12[i];
        const auto t_m11 = pt_m11[i];
        const auto t_m10 = pt_m10[i];
        const auto t_m9 = pt_m9[i];
        const auto t_m8 = pt_m8[i];
        const auto t_m7 = pt_m7[i];
        const auto t_m6 = pt_m6[i];
        const auto t_m5 = pt_m5[i];
        const auto t_m4 = pt_m4[i];
        const auto t_p13 = pt_p13[i];

        const auto u_m12 = pu_m12[i];
        const auto u_m11 = pu_m11[i];
        const auto u_m10 = pu_m10[i];
        const auto u_m9 = pu_m9[i];
        const auto u_m8 = pu_m8[i];
        const auto u_m7 = pu_m7[i];
        const auto u_m6 = pu_m6[i];
        const auto u_m5 = pu_m5[i];
        const auto u_m4 = pu_m4[i];

        ps_p14[i] = fact * (x * t_p13 - y * t_m13);

        ps_m14[i] = fact * (y * t_p13 + x * t_m13);

        ps_m13[i] = fz_m13 * z * t_m13;

        ps_m12[i] = fz_m12 * z * t_m12 - fr_m12 * r_2 * u_m12;

        ps_m11[i] = fz_m11 * z * t_m11 - fr_m11 * r_2 * u_m11;

        ps_m10[i] = fz_m10 * z * t_m10 - fr_m10 * r_2 * u_m10;

        ps_m9[i] = fz_m9 * z * t_m9 - fr_m9 * r_2 * u_m9;

        ps_m8[i] = fz_m8 * z * t_m8 - fr_m8 * r_2 * u_m8;

        ps_m7[i] = fz_m7 * z * t_m7 - fr_m7 * r_2 * u_m7;

        ps_m6[i] = fz_m6 * z * t_m6 - fr_m6 * r_2 * u_m6;

        ps_m5[i] = fz_m5 * z * t_m5 - fr_m5 * r_2 * u_m5;

        ps_m4[i] = fz_m4 * z * t_m4 - fr_m4 * r_2 * u_m4;
    }

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 29 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_m3, pt_m2, pt_m1, pt_0, pt_p1, pt_p2, pt_p3, pt_p4, pu_m3, pu_m2, pu_m1, pu_0, pu_p1, pu_p2, pu_p3, pu_p4, ps_m3, ps_m2, ps_m1, ps_0, ps_p1, ps_p2, ps_p3, ps_p4 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

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
        const auto u_p4 = pu_p4[i];

        ps_m3[i] = fz_m3 * z * t_m3 - fr_m3 * r_2 * u_m3;

        ps_m2[i] = fz_m2 * z * t_m2 - fr_m2 * r_2 * u_m2;

        ps_m1[i] = fz_m1 * z * t_m1 - fr_m1 * r_2 * u_m1;

        ps_0[i] = fz_0 * z * t_0 - fr_0 * r_2 * u_0;

        ps_p1[i] = fz_p1 * z * t_p1 - fr_p1 * r_2 * u_p1;

        ps_p2[i] = fz_p2 * z * t_p2 - fr_p2 * r_2 * u_p2;

        ps_p3[i] = fz_p3 * z * t_p3 - fr_p3 * r_2 * u_p3;

        ps_p4[i] = fz_p4 * z * t_p4 - fr_p4 * r_2 * u_p4;
    }

    // NOTE: the rows are formed in three loops, as the vectorizer runs out of
    // registers with all 29 of them in one. Only the components and the squared
    // distance are loaded by more than one loop, every other value once.

#pragma omp simd aligned(pr_z, pr_2, pt_p5, pt_p6, pt_p7, pt_p8, pt_p9, pt_p10, pt_p11, pt_p12, pt_p13, pu_p5, pu_p6, pu_p7, pu_p8, pu_p9, pu_p10, pu_p11, pu_p12, ps_p5, ps_p6, ps_p7, ps_p8, ps_p9, ps_p10, ps_p11, ps_p12, ps_p13 : simd::cache_line_size())
    for (size_t i = 0; i < ncols; i++)
    {
        const auto z = pr_z[i];

        const auto r_2 = pr_2[i];

        const auto t_p5 = pt_p5[i];
        const auto t_p6 = pt_p6[i];
        const auto t_p7 = pt_p7[i];
        const auto t_p8 = pt_p8[i];
        const auto t_p9 = pt_p9[i];
        const auto t_p10 = pt_p10[i];
        const auto t_p11 = pt_p11[i];
        const auto t_p12 = pt_p12[i];
        const auto t_p13 = pt_p13[i];

        const auto u_p5 = pu_p5[i];
        const auto u_p6 = pu_p6[i];
        const auto u_p7 = pu_p7[i];
        const auto u_p8 = pu_p8[i];
        const auto u_p9 = pu_p9[i];
        const auto u_p10 = pu_p10[i];
        const auto u_p11 = pu_p11[i];
        const auto u_p12 = pu_p12[i];

        ps_p5[i] = fz_p5 * z * t_p5 - fr_p5 * r_2 * u_p5;

        ps_p6[i] = fz_p6 * z * t_p6 - fr_p6 * r_2 * u_p6;

        ps_p7[i] = fz_p7 * z * t_p7 - fr_p7 * r_2 * u_p7;

        ps_p8[i] = fz_p8 * z * t_p8 - fr_p8 * r_2 * u_p8;

        ps_p9[i] = fz_p9 * z * t_p9 - fr_p9 * r_2 * u_p9;

        ps_p10[i] = fz_p10 * z * t_p10 - fr_p10 * r_2 * u_p10;

        ps_p11[i] = fz_p11 * z * t_p11 - fr_p11 * r_2 * u_p11;

        ps_p12[i] = fz_p12 * z * t_p12 - fr_p12 * r_2 * u_p12;

        ps_p13[i] = fz_p13 * z * t_p13;
    }

    return matrix;
}

}  // namespace simdfunc
