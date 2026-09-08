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


#include "SimdTransferGD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gd_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t gp, const size_t hp,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(105.0);
    const auto f_1 = 0.25 * std::sqrt(35.0);
    const auto f_2 = 0.5 * std::sqrt(35.0);
    const auto f_3 = 0.25 * std::sqrt(105.0);
    const auto f_4 = 0.75 * std::sqrt(210.0);
    const auto f_5 = 0.25 * std::sqrt(210.0);
    const auto f_6 = 0.375 * std::sqrt(70.0);
    const auto f_7 = 0.75 * std::sqrt(70.0);
    const auto f_8 = 0.125 * std::sqrt(70.0);
    const auto f_9 = 0.25 * std::sqrt(70.0);
    const auto f_10 = 0.375 * std::sqrt(210.0);
    const auto f_11 = 0.125 * std::sqrt(210.0);
    const auto f_12 = 0.5 * std::sqrt(15.0);
    const auto f_13 = 3.0 * std::sqrt(15.0);
    const auto f_14 = 0.25 * std::sqrt(5.0);
    const auto f_15 = 0.5 * std::sqrt(5.0);
    const auto f_16 = 1.5 * std::sqrt(5.0);
    const auto f_17 = 3.0 * std::sqrt(5.0);
    const auto f_18 = 0.25 * std::sqrt(15.0);
    const auto f_19 = 1.5 * std::sqrt(15.0);
    const auto f_20 = 0.75 * std::sqrt(30.0);
    const auto f_21 = std::sqrt(30.0);
    const auto f_22 = 0.375 * std::sqrt(10.0);
    const auto f_23 = 0.75 * std::sqrt(10.0);
    const auto f_24 = 0.5 * std::sqrt(10.0);
    const auto f_25 = std::sqrt(10.0);
    const auto f_26 = 0.375 * std::sqrt(30.0);
    const auto f_27 = 0.5 * std::sqrt(30.0);
    const auto f_28 = 0.375 * std::sqrt(3.0);
    const auto f_29 = 0.75 * std::sqrt(3.0);
    const auto f_30 = 3.0 * std::sqrt(3.0);
    const auto f_31 = std::sqrt(3.0);
    const auto f_32 = 0.1875 * std::sqrt(3.0);
    const auto f_33 = 1.5 * std::sqrt(3.0);
    const auto f_34 = 0.5 * std::sqrt(3.0);
    const auto f_35 = 0.125 * std::sqrt(5.0);
    const auto f_36 = 0.75 * std::sqrt(5.0);
    const auto f_37 = 0.125 * std::sqrt(15.0);
    const auto f_38 = 0.75 * std::sqrt(15.0);
    const auto f_39 = 0.125 * std::sqrt(105.0);
    const auto f_40 = 0.75 * std::sqrt(105.0);
    const auto f_41 = 0.0625 * std::sqrt(35.0);
    const auto f_42 = 0.125 * std::sqrt(35.0);
    const auto f_43 = 0.375 * std::sqrt(35.0);
    const auto f_44 = 0.75 * std::sqrt(35.0);
    const auto f_45 = 0.0625 * std::sqrt(105.0);
    const auto f_46 = 0.375 * std::sqrt(105.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;
    auto *g_35 = values + 35 * nvalues;
    auto *g_36 = values + 36 * nvalues;
    auto *g_37 = values + 37 * nvalues;
    auto *g_38 = values + 38 * nvalues;
    auto *g_39 = values + 39 * nvalues;
    auto *g_40 = values + 40 * nvalues;
    auto *g_41 = values + 41 * nvalues;
    auto *g_42 = values + 42 * nvalues;
    auto *g_43 = values + 43 * nvalues;
    auto *g_44 = values + 44 * nvalues;

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_62 = buffer.data(hp + 62);

#pragma omp simd aligned(ab_x, ab_y, gp_4, gp_5, gp_19, gp_20, hp_4, hp_11, hp_19, \
                         hp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_x[k] * gp_4[k]
                 - f_0 * ab_x[k] * gp_19[k]
                 + f_0 * hp_4[k]
                 - f_0 * hp_19[k];

        g_1[k] = f_0 * ab_y[k] * gp_5[k]
                 - f_0 * ab_y[k] * gp_20[k]
                 + f_0 * hp_11[k]
                 - f_0 * hp_32[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_3, gp_4, gp_5, gp_18, gp_19, gp_20, hp_3, hp_10, \
                         hp_14, hp_18, hp_31, hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = -f_1 * ab_x[k] * gp_3[k]
                 - f_1 * ab_y[k] * gp_4[k]
                 + f_2 * ab_z[k] * gp_5[k]
                 + f_1 * ab_x[k] * gp_18[k]
                 + f_1 * ab_y[k] * gp_19[k]
                 - f_2 * ab_z[k] * gp_20[k]
                 - f_1 * hp_3[k]
                 - f_1 * hp_10[k]
                 + f_2 * hp_14[k]
                 + f_1 * hp_18[k]
                 + f_1 * hp_31[k]
                 - f_2 * hp_35[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_3, gp_4, gp_5, gp_18, gp_19, gp_20, hp_3, hp_5, hp_10, \
                         hp_18, hp_20, hp_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_0 * ab_x[k] * gp_5[k]
                 - f_0 * ab_x[k] * gp_20[k]
                 + f_0 * hp_5[k]
                 - f_0 * hp_20[k];

        g_4[k] = f_3 * ab_x[k] * gp_3[k]
                 - f_3 * ab_y[k] * gp_4[k]
                 - f_3 * ab_x[k] * gp_18[k]
                 + f_3 * ab_y[k] * gp_19[k]
                 + f_3 * hp_3[k]
                 - f_3 * hp_10[k]
                 - f_3 * hp_18[k]
                 + f_3 * hp_31[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_13, gp_14, gp_34, gp_35, hp_13, hp_23, hp_34, \
                         hp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_4 * ab_x[k] * gp_13[k]
                 - f_5 * ab_x[k] * gp_34[k]
                 + f_4 * hp_13[k]
                 - f_5 * hp_34[k];

        g_6[k] = f_4 * ab_y[k] * gp_14[k]
                 - f_5 * ab_y[k] * gp_35[k]
                 + f_4 * hp_23[k]
                 - f_5 * hp_50[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_12, gp_13, gp_14, gp_33, gp_34, gp_35, hp_12, \
                         hp_22, hp_26, hp_33, hp_49, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_6 * ab_x[k] * gp_12[k]
                 - f_6 * ab_y[k] * gp_13[k]
                 + f_7 * ab_z[k] * gp_14[k]
                 + f_8 * ab_x[k] * gp_33[k]
                 + f_8 * ab_y[k] * gp_34[k]
                 - f_9 * ab_z[k] * gp_35[k]
                 - f_6 * hp_12[k]
                 - f_6 * hp_22[k]
                 + f_7 * hp_26[k]
                 + f_8 * hp_33[k]
                 + f_8 * hp_49[k]
                 - f_9 * hp_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_12, gp_13, gp_14, gp_33, gp_34, gp_35, hp_12, hp_14, \
                         hp_22, hp_33, hp_35, hp_49 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_4 * ab_x[k] * gp_14[k]
                 - f_5 * ab_x[k] * gp_35[k]
                 + f_4 * hp_14[k]
                 - f_5 * hp_35[k];

        g_9[k] = f_10 * ab_x[k] * gp_12[k]
                 - f_10 * ab_y[k] * gp_13[k]
                 - f_11 * ab_x[k] * gp_33[k]
                 + f_11 * ab_y[k] * gp_34[k]
                 + f_10 * hp_12[k]
                 - f_10 * hp_22[k]
                 - f_11 * hp_33[k]
                 + f_11 * hp_49[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_4, gp_5, gp_19, gp_20, gp_25, gp_26, hp_4, hp_11, \
                         hp_19, hp_25, hp_32, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_12 * ab_x[k] * gp_4[k]
                  - f_12 * ab_x[k] * gp_19[k]
                  + f_13 * ab_x[k] * gp_25[k]
                  - f_12 * hp_4[k]
                  - f_12 * hp_19[k]
                  + f_13 * hp_25[k];

        g_11[k] = -f_12 * ab_y[k] * gp_5[k]
                  - f_12 * ab_y[k] * gp_20[k]
                  + f_13 * ab_y[k] * gp_26[k]
                  - f_12 * hp_11[k]
                  - f_12 * hp_32[k]
                  + f_13 * hp_38[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_3, gp_4, gp_5, gp_18, gp_19, gp_20, gp_24, \
                         gp_25, gp_26, hp_3, hp_10, hp_14, hp_18, hp_24, hp_31, hp_35, hp_37, \
                         hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_14 * ab_x[k] * gp_3[k]
                  + f_14 * ab_y[k] * gp_4[k]
                  - f_15 * ab_z[k] * gp_5[k]
                  + f_14 * ab_x[k] * gp_18[k]
                  + f_14 * ab_y[k] * gp_19[k]
                  - f_15 * ab_z[k] * gp_20[k]
                  - f_16 * ab_x[k] * gp_24[k]
                  - f_16 * ab_y[k] * gp_25[k]
                  + f_17 * ab_z[k] * gp_26[k]
                  + f_14 * hp_3[k]
                  + f_14 * hp_10[k]
                  - f_15 * hp_14[k]
                  + f_14 * hp_18[k]
                  - f_16 * hp_24[k]
                  + f_14 * hp_31[k]
                  - f_15 * hp_35[k]
                  - f_16 * hp_37[k]
                  + f_17 * hp_41[k];
    }

#pragma omp simd aligned(ab_x, gp_5, gp_20, gp_26, hp_5, hp_20, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_12 * ab_x[k] * gp_5[k]
                  - f_12 * ab_x[k] * gp_20[k]
                  + f_13 * ab_x[k] * gp_26[k]
                  - f_12 * hp_5[k]
                  - f_12 * hp_20[k]
                  + f_13 * hp_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_3, gp_4, gp_18, gp_19, gp_24, gp_25, hp_3, hp_10, \
                         hp_18, hp_24, hp_31, hp_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = -f_18 * ab_x[k] * gp_3[k]
                  + f_18 * ab_y[k] * gp_4[k]
                  - f_18 * ab_x[k] * gp_18[k]
                  + f_18 * ab_y[k] * gp_19[k]
                  + f_19 * ab_x[k] * gp_24[k]
                  - f_19 * ab_y[k] * gp_25[k]
                  - f_18 * hp_3[k]
                  + f_18 * hp_10[k]
                  - f_18 * hp_18[k]
                  + f_19 * hp_24[k]
                  + f_18 * hp_31[k]
                  - f_19 * hp_37[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_13, gp_14, gp_34, gp_35, gp_40, gp_41, hp_13, hp_23, \
                         hp_34, hp_40, hp_50, hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_20 * ab_x[k] * gp_13[k]
                  - f_20 * ab_x[k] * gp_34[k]
                  + f_21 * ab_x[k] * gp_40[k]
                  - f_20 * hp_13[k]
                  - f_20 * hp_34[k]
                  + f_21 * hp_40[k];

        g_16[k] = -f_20 * ab_y[k] * gp_14[k]
                  - f_20 * ab_y[k] * gp_35[k]
                  + f_21 * ab_y[k] * gp_41[k]
                  - f_20 * hp_23[k]
                  - f_20 * hp_50[k]
                  + f_21 * hp_56[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_12, gp_13, gp_14, gp_33, gp_34, gp_35, gp_39, \
                         gp_40, gp_41, hp_12, hp_22, hp_26, hp_33, hp_39, hp_49, hp_53, hp_55, \
                         hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = f_22 * ab_x[k] * gp_12[k]
                  + f_22 * ab_y[k] * gp_13[k]
                  - f_23 * ab_z[k] * gp_14[k]
                  + f_22 * ab_x[k] * gp_33[k]
                  + f_22 * ab_y[k] * gp_34[k]
                  - f_23 * ab_z[k] * gp_35[k]
                  - f_24 * ab_x[k] * gp_39[k]
                  - f_24 * ab_y[k] * gp_40[k]
                  + f_25 * ab_z[k] * gp_41[k]
                  + f_22 * hp_12[k]
                  + f_22 * hp_22[k]
                  - f_23 * hp_26[k]
                  + f_22 * hp_33[k]
                  - f_24 * hp_39[k]
                  + f_22 * hp_49[k]
                  - f_23 * hp_53[k]
                  - f_24 * hp_55[k]
                  + f_25 * hp_59[k];
    }

#pragma omp simd aligned(ab_x, gp_14, gp_35, gp_41, hp_14, hp_35, \
                         hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_20 * ab_x[k] * gp_14[k]
                  - f_20 * ab_x[k] * gp_35[k]
                  + f_21 * ab_x[k] * gp_41[k]
                  - f_20 * hp_14[k]
                  - f_20 * hp_35[k]
                  + f_21 * hp_41[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_12, gp_13, gp_33, gp_34, gp_39, gp_40, hp_12, hp_22, \
                         hp_33, hp_39, hp_49, hp_55 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_26 * ab_x[k] * gp_12[k]
                  + f_26 * ab_y[k] * gp_13[k]
                  - f_26 * ab_x[k] * gp_33[k]
                  + f_26 * ab_y[k] * gp_34[k]
                  + f_27 * ab_x[k] * gp_39[k]
                  - f_27 * ab_y[k] * gp_40[k]
                  - f_26 * hp_12[k]
                  + f_26 * hp_22[k]
                  - f_26 * hp_33[k]
                  + f_27 * hp_39[k]
                  + f_26 * hp_49[k]
                  - f_27 * hp_55[k];
    }

#pragma omp simd aligned(ab_x, gp_1, gp_10, gp_16, gp_31, gp_37, gp_43, hp_1, hp_10, hp_16, \
                         hp_31, hp_37, hp_43 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = f_28 * ab_x[k] * gp_1[k]
                  + f_29 * ab_x[k] * gp_10[k]
                  - f_30 * ab_x[k] * gp_16[k]
                  + f_28 * ab_x[k] * gp_31[k]
                  - f_30 * ab_x[k] * gp_37[k]
                  + f_31 * ab_x[k] * gp_43[k]
                  + f_28 * hp_1[k]
                  + f_29 * hp_10[k]
                  - f_30 * hp_16[k]
                  + f_28 * hp_31[k]
                  - f_30 * hp_37[k]
                  + f_31 * hp_43[k];
    }

#pragma omp simd aligned(ab_y, gp_2, gp_11, gp_17, gp_32, gp_38, gp_44, hp_5, hp_20, hp_26, \
                         hp_47, hp_53, hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = f_28 * ab_y[k] * gp_2[k]
                  + f_29 * ab_y[k] * gp_11[k]
                  - f_30 * ab_y[k] * gp_17[k]
                  + f_28 * ab_y[k] * gp_32[k]
                  - f_30 * ab_y[k] * gp_38[k]
                  + f_31 * ab_y[k] * gp_44[k]
                  + f_28 * hp_5[k]
                  + f_29 * hp_20[k]
                  - f_30 * hp_26[k]
                  + f_28 * hp_47[k]
                  - f_30 * hp_53[k]
                  + f_31 * hp_59[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_0, gp_1, gp_2, gp_9, gp_10, gp_11, gp_15, gp_16, \
                         gp_17, gp_30, gp_31, gp_32, gp_36, gp_37, gp_38, gp_42, gp_43, gp_44, \
                         hp_0, hp_4, hp_8, hp_9, hp_15, hp_19, hp_23, hp_25, hp_29, hp_30, \
                         hp_36, hp_42, hp_46, hp_50, hp_52, hp_56, hp_58, \
                         hp_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -0.1875 * ab_x[k] * gp_0[k]
                  - 0.1875 * ab_y[k] * gp_1[k]
                  + 0.375 * ab_z[k] * gp_2[k]
                  - 0.375 * ab_x[k] * gp_9[k]
                  - 0.375 * ab_y[k] * gp_10[k]
                  + 0.75 * ab_z[k] * gp_11[k]
                  + 1.5 * ab_x[k] * gp_15[k]
                  + 1.5 * ab_y[k] * gp_16[k]
                  - 3.0 * ab_z[k] * gp_17[k]
                  - 0.1875 * ab_x[k] * gp_30[k]
                  - 0.1875 * ab_y[k] * gp_31[k]
                  + 0.375 * ab_z[k] * gp_32[k]
                  + 1.5 * ab_x[k] * gp_36[k]
                  + 1.5 * ab_y[k] * gp_37[k]
                  - 3.0 * ab_z[k] * gp_38[k]
                  - 0.5 * ab_x[k] * gp_42[k]
                  - 0.5 * ab_y[k] * gp_43[k]
                  + ab_z[k] * gp_44[k]
                  - 0.1875 * hp_0[k]
                  - 0.1875 * hp_4[k]
                  + 0.375 * hp_8[k]
                  - 0.375 * hp_9[k]
                  + 1.5 * hp_15[k]
                  - 0.375 * hp_19[k]
                  + 0.75 * hp_23[k]
                  + 1.5 * hp_25[k]
                  - 3.0 * hp_29[k]
                  - 0.1875 * hp_30[k]
                  + 1.5 * hp_36[k]
                  - 0.5 * hp_42[k]
                  - 0.1875 * hp_46[k]
                  + 0.375 * hp_50[k]
                  + 1.5 * hp_52[k]
                  - 3.0 * hp_56[k]
                  - 0.5 * hp_58[k]
                  + hp_62[k];
    }

#pragma omp simd aligned(ab_x, gp_2, gp_11, gp_17, gp_32, gp_38, gp_44, hp_2, hp_11, hp_17, \
                         hp_32, hp_38, hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_28 * ab_x[k] * gp_2[k]
                  + f_29 * ab_x[k] * gp_11[k]
                  - f_30 * ab_x[k] * gp_17[k]
                  + f_28 * ab_x[k] * gp_32[k]
                  - f_30 * ab_x[k] * gp_38[k]
                  + f_31 * ab_x[k] * gp_44[k]
                  + f_28 * hp_2[k]
                  + f_29 * hp_11[k]
                  - f_30 * hp_17[k]
                  + f_28 * hp_32[k]
                  - f_30 * hp_38[k]
                  + f_31 * hp_44[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_0, gp_1, gp_9, gp_10, gp_15, gp_16, gp_30, gp_31, \
                         gp_36, gp_37, gp_42, gp_43, hp_0, hp_4, hp_9, hp_15, hp_19, hp_25, \
                         hp_30, hp_36, hp_42, hp_46, hp_52, hp_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_32 * ab_x[k] * gp_0[k]
                  - f_32 * ab_y[k] * gp_1[k]
                  + f_28 * ab_x[k] * gp_9[k]
                  - f_28 * ab_y[k] * gp_10[k]
                  - f_33 * ab_x[k] * gp_15[k]
                  + f_33 * ab_y[k] * gp_16[k]
                  + f_32 * ab_x[k] * gp_30[k]
                  - f_32 * ab_y[k] * gp_31[k]
                  - f_33 * ab_x[k] * gp_36[k]
                  + f_33 * ab_y[k] * gp_37[k]
                  + f_34 * ab_x[k] * gp_42[k]
                  - f_34 * ab_y[k] * gp_43[k]
                  + f_32 * hp_0[k]
                  - f_32 * hp_4[k]
                  + f_28 * hp_9[k]
                  - f_33 * hp_15[k]
                  - f_28 * hp_19[k]
                  + f_33 * hp_25[k]
                  + f_32 * hp_30[k]
                  - f_33 * hp_36[k]
                  + f_34 * hp_42[k]
                  - f_32 * hp_46[k]
                  + f_33 * hp_52[k]
                  - f_34 * hp_58[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_7, gp_8, gp_22, gp_23, gp_28, gp_29, hp_7, hp_14, \
                         hp_22, hp_28, hp_35, hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_20 * ab_x[k] * gp_7[k]
                  - f_20 * ab_x[k] * gp_22[k]
                  + f_21 * ab_x[k] * gp_28[k]
                  - f_20 * hp_7[k]
                  - f_20 * hp_22[k]
                  + f_21 * hp_28[k];

        g_26[k] = -f_20 * ab_y[k] * gp_8[k]
                  - f_20 * ab_y[k] * gp_23[k]
                  + f_21 * ab_y[k] * gp_29[k]
                  - f_20 * hp_14[k]
                  - f_20 * hp_35[k]
                  + f_21 * hp_41[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_6, gp_7, gp_8, gp_21, gp_22, gp_23, gp_27, \
                         gp_28, gp_29, hp_6, hp_13, hp_17, hp_21, hp_27, hp_34, hp_38, hp_40, \
                         hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = f_22 * ab_x[k] * gp_6[k]
                  + f_22 * ab_y[k] * gp_7[k]
                  - f_23 * ab_z[k] * gp_8[k]
                  + f_22 * ab_x[k] * gp_21[k]
                  + f_22 * ab_y[k] * gp_22[k]
                  - f_23 * ab_z[k] * gp_23[k]
                  - f_24 * ab_x[k] * gp_27[k]
                  - f_24 * ab_y[k] * gp_28[k]
                  + f_25 * ab_z[k] * gp_29[k]
                  + f_22 * hp_6[k]
                  + f_22 * hp_13[k]
                  - f_23 * hp_17[k]
                  + f_22 * hp_21[k]
                  - f_24 * hp_27[k]
                  + f_22 * hp_34[k]
                  - f_23 * hp_38[k]
                  - f_24 * hp_40[k]
                  + f_25 * hp_44[k];
    }

#pragma omp simd aligned(ab_x, gp_8, gp_23, gp_29, hp_8, hp_23, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_20 * ab_x[k] * gp_8[k]
                  - f_20 * ab_x[k] * gp_23[k]
                  + f_21 * ab_x[k] * gp_29[k]
                  - f_20 * hp_8[k]
                  - f_20 * hp_23[k]
                  + f_21 * hp_29[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_6, gp_7, gp_21, gp_22, gp_27, gp_28, hp_6, hp_13, \
                         hp_21, hp_27, hp_34, hp_40 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = -f_26 * ab_x[k] * gp_6[k]
                  + f_26 * ab_y[k] * gp_7[k]
                  - f_26 * ab_x[k] * gp_21[k]
                  + f_26 * ab_y[k] * gp_22[k]
                  + f_27 * ab_x[k] * gp_27[k]
                  - f_27 * ab_y[k] * gp_28[k]
                  - f_26 * hp_6[k]
                  + f_26 * hp_13[k]
                  - f_26 * hp_21[k]
                  + f_27 * hp_27[k]
                  + f_26 * hp_34[k]
                  - f_27 * hp_40[k];
    }

#pragma omp simd aligned(ab_x, gp_1, gp_16, gp_31, gp_37, hp_1, hp_16, hp_31, \
                         hp_37 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_18 * ab_x[k] * gp_1[k]
                  + f_19 * ab_x[k] * gp_16[k]
                  + f_18 * ab_x[k] * gp_31[k]
                  - f_19 * ab_x[k] * gp_37[k]
                  - f_18 * hp_1[k]
                  + f_19 * hp_16[k]
                  + f_18 * hp_31[k]
                  - f_19 * hp_37[k];
    }

#pragma omp simd aligned(ab_y, gp_2, gp_17, gp_32, gp_38, hp_5, hp_26, hp_47, \
                         hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_18 * ab_y[k] * gp_2[k]
                  + f_19 * ab_y[k] * gp_17[k]
                  + f_18 * ab_y[k] * gp_32[k]
                  - f_19 * ab_y[k] * gp_38[k]
                  - f_18 * hp_5[k]
                  + f_19 * hp_26[k]
                  + f_18 * hp_47[k]
                  - f_19 * hp_53[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_0, gp_1, gp_2, gp_15, gp_16, gp_17, gp_30, \
                         gp_31, gp_32, gp_36, gp_37, gp_38, hp_0, hp_4, hp_8, hp_15, hp_25, \
                         hp_29, hp_30, hp_36, hp_46, hp_50, hp_52, \
                         hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = f_35 * ab_x[k] * gp_0[k]
                  + f_35 * ab_y[k] * gp_1[k]
                  - f_14 * ab_z[k] * gp_2[k]
                  - f_36 * ab_x[k] * gp_15[k]
                  - f_36 * ab_y[k] * gp_16[k]
                  + f_16 * ab_z[k] * gp_17[k]
                  - f_35 * ab_x[k] * gp_30[k]
                  - f_35 * ab_y[k] * gp_31[k]
                  + f_14 * ab_z[k] * gp_32[k]
                  + f_36 * ab_x[k] * gp_36[k]
                  + f_36 * ab_y[k] * gp_37[k]
                  - f_16 * ab_z[k] * gp_38[k]
                  + f_35 * hp_0[k]
                  + f_35 * hp_4[k]
                  - f_14 * hp_8[k]
                  - f_36 * hp_15[k]
                  - f_36 * hp_25[k]
                  + f_16 * hp_29[k]
                  - f_35 * hp_30[k]
                  + f_36 * hp_36[k]
                  - f_35 * hp_46[k]
                  + f_14 * hp_50[k]
                  + f_36 * hp_52[k]
                  - f_16 * hp_56[k];
    }

#pragma omp simd aligned(ab_x, gp_2, gp_17, gp_32, gp_38, hp_2, hp_17, hp_32, \
                         hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_18 * ab_x[k] * gp_2[k]
                  + f_19 * ab_x[k] * gp_17[k]
                  + f_18 * ab_x[k] * gp_32[k]
                  - f_19 * ab_x[k] * gp_38[k]
                  - f_18 * hp_2[k]
                  + f_19 * hp_17[k]
                  + f_18 * hp_32[k]
                  - f_19 * hp_38[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_0, gp_1, gp_15, gp_16, gp_30, gp_31, gp_36, gp_37, \
                         hp_0, hp_4, hp_15, hp_25, hp_30, hp_36, hp_46, \
                         hp_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_37 * ab_x[k] * gp_0[k]
                  + f_37 * ab_y[k] * gp_1[k]
                  + f_38 * ab_x[k] * gp_15[k]
                  - f_38 * ab_y[k] * gp_16[k]
                  + f_37 * ab_x[k] * gp_30[k]
                  - f_37 * ab_y[k] * gp_31[k]
                  - f_38 * ab_x[k] * gp_36[k]
                  + f_38 * ab_y[k] * gp_37[k]
                  - f_37 * hp_0[k]
                  + f_37 * hp_4[k]
                  + f_38 * hp_15[k]
                  - f_38 * hp_25[k]
                  + f_37 * hp_30[k]
                  - f_38 * hp_36[k]
                  - f_37 * hp_46[k]
                  + f_38 * hp_52[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_7, gp_8, gp_22, gp_23, hp_7, hp_14, hp_22, \
                         hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_5 * ab_x[k] * gp_7[k]
                  - f_4 * ab_x[k] * gp_22[k]
                  + f_5 * hp_7[k]
                  - f_4 * hp_22[k];

        g_36[k] = f_5 * ab_y[k] * gp_8[k]
                  - f_4 * ab_y[k] * gp_23[k]
                  + f_5 * hp_14[k]
                  - f_4 * hp_35[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_6, gp_7, gp_8, gp_21, gp_22, gp_23, hp_6, hp_13, \
                         hp_17, hp_21, hp_34, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_8 * ab_x[k] * gp_6[k]
                  - f_8 * ab_y[k] * gp_7[k]
                  + f_9 * ab_z[k] * gp_8[k]
                  + f_6 * ab_x[k] * gp_21[k]
                  + f_6 * ab_y[k] * gp_22[k]
                  - f_7 * ab_z[k] * gp_23[k]
                  - f_8 * hp_6[k]
                  - f_8 * hp_13[k]
                  + f_9 * hp_17[k]
                  + f_6 * hp_21[k]
                  + f_6 * hp_34[k]
                  - f_7 * hp_38[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_6, gp_7, gp_8, gp_21, gp_22, gp_23, hp_6, hp_8, hp_13, \
                         hp_21, hp_23, hp_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = f_5 * ab_x[k] * gp_8[k]
                  - f_4 * ab_x[k] * gp_23[k]
                  + f_5 * hp_8[k]
                  - f_4 * hp_23[k];

        g_39[k] = f_11 * ab_x[k] * gp_6[k]
                  - f_11 * ab_y[k] * gp_7[k]
                  - f_10 * ab_x[k] * gp_21[k]
                  + f_10 * ab_y[k] * gp_22[k]
                  + f_11 * hp_6[k]
                  - f_11 * hp_13[k]
                  - f_10 * hp_21[k]
                  + f_10 * hp_34[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_1, gp_2, gp_10, gp_11, gp_31, gp_32, hp_1, hp_5, \
                         hp_10, hp_20, hp_31, hp_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_40[k] = f_39 * ab_x[k] * gp_1[k]
                  - f_40 * ab_x[k] * gp_10[k]
                  + f_39 * ab_x[k] * gp_31[k]
                  + f_39 * hp_1[k]
                  - f_40 * hp_10[k]
                  + f_39 * hp_31[k];

        g_41[k] = f_39 * ab_y[k] * gp_2[k]
                  - f_40 * ab_y[k] * gp_11[k]
                  + f_39 * ab_y[k] * gp_32[k]
                  + f_39 * hp_5[k]
                  - f_40 * hp_20[k]
                  + f_39 * hp_47[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gp_0, gp_1, gp_2, gp_9, gp_10, gp_11, gp_30, gp_31, \
                         gp_32, hp_0, hp_4, hp_8, hp_9, hp_19, hp_23, hp_30, hp_46, \
                         hp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_42[k] = -f_41 * ab_x[k] * gp_0[k]
                  - f_41 * ab_y[k] * gp_1[k]
                  + f_42 * ab_z[k] * gp_2[k]
                  + f_43 * ab_x[k] * gp_9[k]
                  + f_43 * ab_y[k] * gp_10[k]
                  - f_44 * ab_z[k] * gp_11[k]
                  - f_41 * ab_x[k] * gp_30[k]
                  - f_41 * ab_y[k] * gp_31[k]
                  + f_42 * ab_z[k] * gp_32[k]
                  - f_41 * hp_0[k]
                  - f_41 * hp_4[k]
                  + f_42 * hp_8[k]
                  + f_43 * hp_9[k]
                  + f_43 * hp_19[k]
                  - f_44 * hp_23[k]
                  - f_41 * hp_30[k]
                  - f_41 * hp_46[k]
                  + f_42 * hp_50[k];
    }

#pragma omp simd aligned(ab_x, gp_2, gp_11, gp_32, hp_2, hp_11, hp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_43[k] = f_39 * ab_x[k] * gp_2[k]
                  - f_40 * ab_x[k] * gp_11[k]
                  + f_39 * ab_x[k] * gp_32[k]
                  + f_39 * hp_2[k]
                  - f_40 * hp_11[k]
                  + f_39 * hp_32[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gp_0, gp_1, gp_9, gp_10, gp_30, gp_31, hp_0, hp_4, hp_9, \
                         hp_19, hp_30, hp_46 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_44[k] = f_45 * ab_x[k] * gp_0[k]
                  - f_45 * ab_y[k] * gp_1[k]
                  - f_46 * ab_x[k] * gp_9[k]
                  + f_46 * ab_y[k] * gp_10[k]
                  + f_45 * ab_x[k] * gp_30[k]
                  - f_45 * ab_y[k] * gp_31[k]
                  + f_45 * hp_0[k]
                  - f_45 * hp_4[k]
                  - f_46 * hp_9[k]
                  + f_46 * hp_19[k]
                  + f_45 * hp_30[k]
                  - f_45 * hp_46[k];
    }
}

auto
compute_hrr_gd(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gp, const size_t hp, const size_t nmax) -> void
{
    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);
    auto *t_3 = buffer.data(target + 3);
    auto *t_4 = buffer.data(target + 4);
    auto *t_5 = buffer.data(target + 5);
    auto *t_6 = buffer.data(target + 6);
    auto *t_7 = buffer.data(target + 7);
    auto *t_8 = buffer.data(target + 8);
    auto *t_9 = buffer.data(target + 9);
    auto *t_10 = buffer.data(target + 10);
    auto *t_11 = buffer.data(target + 11);
    auto *t_12 = buffer.data(target + 12);
    auto *t_13 = buffer.data(target + 13);
    auto *t_14 = buffer.data(target + 14);
    auto *t_15 = buffer.data(target + 15);
    auto *t_16 = buffer.data(target + 16);
    auto *t_17 = buffer.data(target + 17);
    auto *t_18 = buffer.data(target + 18);
    auto *t_19 = buffer.data(target + 19);
    auto *t_20 = buffer.data(target + 20);
    auto *t_21 = buffer.data(target + 21);
    auto *t_22 = buffer.data(target + 22);
    auto *t_23 = buffer.data(target + 23);
    auto *t_24 = buffer.data(target + 24);
    auto *t_25 = buffer.data(target + 25);
    auto *t_26 = buffer.data(target + 26);
    auto *t_27 = buffer.data(target + 27);
    auto *t_28 = buffer.data(target + 28);
    auto *t_29 = buffer.data(target + 29);
    auto *t_30 = buffer.data(target + 30);
    auto *t_31 = buffer.data(target + 31);
    auto *t_32 = buffer.data(target + 32);
    auto *t_33 = buffer.data(target + 33);
    auto *t_34 = buffer.data(target + 34);
    auto *t_35 = buffer.data(target + 35);
    auto *t_36 = buffer.data(target + 36);
    auto *t_37 = buffer.data(target + 37);
    auto *t_38 = buffer.data(target + 38);
    auto *t_39 = buffer.data(target + 39);
    auto *t_40 = buffer.data(target + 40);
    auto *t_41 = buffer.data(target + 41);
    auto *t_42 = buffer.data(target + 42);
    auto *t_43 = buffer.data(target + 43);
    auto *t_44 = buffer.data(target + 44);
    auto *t_45 = buffer.data(target + 45);
    auto *t_46 = buffer.data(target + 46);
    auto *t_47 = buffer.data(target + 47);
    auto *t_48 = buffer.data(target + 48);
    auto *t_49 = buffer.data(target + 49);
    auto *t_50 = buffer.data(target + 50);
    auto *t_51 = buffer.data(target + 51);
    auto *t_52 = buffer.data(target + 52);
    auto *t_53 = buffer.data(target + 53);
    auto *t_54 = buffer.data(target + 54);
    auto *t_55 = buffer.data(target + 55);
    auto *t_56 = buffer.data(target + 56);
    auto *t_57 = buffer.data(target + 57);
    auto *t_58 = buffer.data(target + 58);
    auto *t_59 = buffer.data(target + 59);
    auto *t_60 = buffer.data(target + 60);
    auto *t_61 = buffer.data(target + 61);
    auto *t_62 = buffer.data(target + 62);
    auto *t_63 = buffer.data(target + 63);
    auto *t_64 = buffer.data(target + 64);
    auto *t_65 = buffer.data(target + 65);
    auto *t_66 = buffer.data(target + 66);
    auto *t_67 = buffer.data(target + 67);
    auto *t_68 = buffer.data(target + 68);
    auto *t_69 = buffer.data(target + 69);
    auto *t_70 = buffer.data(target + 70);
    auto *t_71 = buffer.data(target + 71);
    auto *t_72 = buffer.data(target + 72);
    auto *t_73 = buffer.data(target + 73);
    auto *t_74 = buffer.data(target + 74);
    auto *t_75 = buffer.data(target + 75);
    auto *t_76 = buffer.data(target + 76);
    auto *t_77 = buffer.data(target + 77);
    auto *t_78 = buffer.data(target + 78);
    auto *t_79 = buffer.data(target + 79);
    auto *t_80 = buffer.data(target + 80);
    auto *t_81 = buffer.data(target + 81);
    auto *t_82 = buffer.data(target + 82);
    auto *t_83 = buffer.data(target + 83);
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_3 = buffer.data(hp + 3);
    const auto *hp_4 = buffer.data(hp + 4);
    const auto *hp_5 = buffer.data(hp + 5);
    const auto *hp_6 = buffer.data(hp + 6);
    const auto *hp_7 = buffer.data(hp + 7);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_9 = buffer.data(hp + 9);
    const auto *hp_10 = buffer.data(hp + 10);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_12 = buffer.data(hp + 12);
    const auto *hp_13 = buffer.data(hp + 13);
    const auto *hp_14 = buffer.data(hp + 14);
    const auto *hp_15 = buffer.data(hp + 15);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_18 = buffer.data(hp + 18);
    const auto *hp_19 = buffer.data(hp + 19);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_21 = buffer.data(hp + 21);
    const auto *hp_22 = buffer.data(hp + 22);
    const auto *hp_23 = buffer.data(hp + 23);
    const auto *hp_24 = buffer.data(hp + 24);
    const auto *hp_25 = buffer.data(hp + 25);
    const auto *hp_26 = buffer.data(hp + 26);
    const auto *hp_27 = buffer.data(hp + 27);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);
    const auto *hp_30 = buffer.data(hp + 30);
    const auto *hp_31 = buffer.data(hp + 31);
    const auto *hp_32 = buffer.data(hp + 32);
    const auto *hp_33 = buffer.data(hp + 33);
    const auto *hp_34 = buffer.data(hp + 34);
    const auto *hp_35 = buffer.data(hp + 35);
    const auto *hp_36 = buffer.data(hp + 36);
    const auto *hp_37 = buffer.data(hp + 37);
    const auto *hp_38 = buffer.data(hp + 38);
    const auto *hp_39 = buffer.data(hp + 39);
    const auto *hp_40 = buffer.data(hp + 40);
    const auto *hp_41 = buffer.data(hp + 41);
    const auto *hp_42 = buffer.data(hp + 42);
    const auto *hp_43 = buffer.data(hp + 43);
    const auto *hp_44 = buffer.data(hp + 44);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_49 = buffer.data(hp + 49);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_58 = buffer.data(hp + 58);
    const auto *hp_59 = buffer.data(hp + 59);
    const auto *hp_62 = buffer.data(hp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, ab_y, gp_0, gp_1, gp_2, hp_0, hp_1, \
                         hp_2, hp_4, hp_5 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * gp_0[k]
                 + hp_0[k];

        t_1[k] = ab_x[k] * gp_1[k]
                 + hp_1[k];

        t_2[k] = ab_x[k] * gp_2[k]
                 + hp_2[k];

        t_3[k] = ab_y[k] * gp_1[k]
                 + hp_4[k];

        t_4[k] = ab_y[k] * gp_2[k]
                 + hp_5[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, ab_x, ab_z, gp_2, gp_3, gp_4, gp_5, hp_3, hp_4, \
                         hp_5, hp_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = ab_z[k] * gp_2[k]
                 + hp_8[k];

        t_6[k] = ab_x[k] * gp_3[k]
                 + hp_3[k];

        t_7[k] = ab_x[k] * gp_4[k]
                 + hp_4[k];

        t_8[k] = ab_x[k] * gp_5[k]
                 + hp_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, ab_x, ab_y, ab_z, gp_4, gp_5, gp_6, hp_6, \
                         hp_10, hp_11, hp_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_9[k] = ab_y[k] * gp_4[k]
                 + hp_10[k];

        t_10[k] = ab_y[k] * gp_5[k]
                  + hp_11[k];

        t_11[k] = ab_z[k] * gp_5[k]
                  + hp_14[k];

        t_12[k] = ab_x[k] * gp_6[k]
                  + hp_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, ab_x, ab_y, ab_z, gp_7, gp_8, hp_7, \
                         hp_8, hp_13, hp_14, hp_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_13[k] = ab_x[k] * gp_7[k]
                  + hp_7[k];

        t_14[k] = ab_x[k] * gp_8[k]
                  + hp_8[k];

        t_15[k] = ab_y[k] * gp_7[k]
                  + hp_13[k];

        t_16[k] = ab_y[k] * gp_8[k]
                  + hp_14[k];

        t_17[k] = ab_z[k] * gp_8[k]
                  + hp_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, ab_x, ab_y, gp_9, gp_10, gp_11, hp_9, \
                         hp_10, hp_11, hp_19, hp_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_18[k] = ab_x[k] * gp_9[k]
                  + hp_9[k];

        t_19[k] = ab_x[k] * gp_10[k]
                  + hp_10[k];

        t_20[k] = ab_x[k] * gp_11[k]
                  + hp_11[k];

        t_21[k] = ab_y[k] * gp_10[k]
                  + hp_19[k];

        t_22[k] = ab_y[k] * gp_11[k]
                  + hp_20[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, ab_x, ab_z, gp_11, gp_12, gp_13, gp_14, \
                         hp_12, hp_13, hp_14, hp_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_23[k] = ab_z[k] * gp_11[k]
                  + hp_23[k];

        t_24[k] = ab_x[k] * gp_12[k]
                  + hp_12[k];

        t_25[k] = ab_x[k] * gp_13[k]
                  + hp_13[k];

        t_26[k] = ab_x[k] * gp_14[k]
                  + hp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, ab_x, ab_y, ab_z, gp_13, gp_14, gp_15, hp_15, \
                         hp_22, hp_23, hp_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_y[k] * gp_13[k]
                  + hp_22[k];

        t_28[k] = ab_y[k] * gp_14[k]
                  + hp_23[k];

        t_29[k] = ab_z[k] * gp_14[k]
                  + hp_26[k];

        t_30[k] = ab_x[k] * gp_15[k]
                  + hp_15[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, ab_x, ab_y, ab_z, gp_16, gp_17, hp_16, \
                         hp_17, hp_25, hp_26, hp_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_31[k] = ab_x[k] * gp_16[k]
                  + hp_16[k];

        t_32[k] = ab_x[k] * gp_17[k]
                  + hp_17[k];

        t_33[k] = ab_y[k] * gp_16[k]
                  + hp_25[k];

        t_34[k] = ab_y[k] * gp_17[k]
                  + hp_26[k];

        t_35[k] = ab_z[k] * gp_17[k]
                  + hp_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, ab_x, ab_y, gp_18, gp_19, gp_20, hp_18, \
                         hp_19, hp_20, hp_31, hp_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_36[k] = ab_x[k] * gp_18[k]
                  + hp_18[k];

        t_37[k] = ab_x[k] * gp_19[k]
                  + hp_19[k];

        t_38[k] = ab_x[k] * gp_20[k]
                  + hp_20[k];

        t_39[k] = ab_y[k] * gp_19[k]
                  + hp_31[k];

        t_40[k] = ab_y[k] * gp_20[k]
                  + hp_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, ab_x, ab_z, gp_20, gp_21, gp_22, gp_23, \
                         hp_21, hp_22, hp_23, hp_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_41[k] = ab_z[k] * gp_20[k]
                  + hp_35[k];

        t_42[k] = ab_x[k] * gp_21[k]
                  + hp_21[k];

        t_43[k] = ab_x[k] * gp_22[k]
                  + hp_22[k];

        t_44[k] = ab_x[k] * gp_23[k]
                  + hp_23[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, ab_x, ab_y, ab_z, gp_22, gp_23, gp_24, hp_24, \
                         hp_34, hp_35, hp_38 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_45[k] = ab_y[k] * gp_22[k]
                  + hp_34[k];

        t_46[k] = ab_y[k] * gp_23[k]
                  + hp_35[k];

        t_47[k] = ab_z[k] * gp_23[k]
                  + hp_38[k];

        t_48[k] = ab_x[k] * gp_24[k]
                  + hp_24[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_x, ab_y, ab_z, gp_25, gp_26, hp_25, \
                         hp_26, hp_37, hp_38, hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = ab_x[k] * gp_25[k]
                  + hp_25[k];

        t_50[k] = ab_x[k] * gp_26[k]
                  + hp_26[k];

        t_51[k] = ab_y[k] * gp_25[k]
                  + hp_37[k];

        t_52[k] = ab_y[k] * gp_26[k]
                  + hp_38[k];

        t_53[k] = ab_z[k] * gp_26[k]
                  + hp_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, ab_x, ab_y, gp_27, gp_28, gp_29, hp_27, \
                         hp_28, hp_29, hp_40, hp_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = ab_x[k] * gp_27[k]
                  + hp_27[k];

        t_55[k] = ab_x[k] * gp_28[k]
                  + hp_28[k];

        t_56[k] = ab_x[k] * gp_29[k]
                  + hp_29[k];

        t_57[k] = ab_y[k] * gp_28[k]
                  + hp_40[k];

        t_58[k] = ab_y[k] * gp_29[k]
                  + hp_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, ab_x, ab_z, gp_29, gp_30, gp_31, gp_32, \
                         hp_30, hp_31, hp_32, hp_44 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_59[k] = ab_z[k] * gp_29[k]
                  + hp_44[k];

        t_60[k] = ab_x[k] * gp_30[k]
                  + hp_30[k];

        t_61[k] = ab_x[k] * gp_31[k]
                  + hp_31[k];

        t_62[k] = ab_x[k] * gp_32[k]
                  + hp_32[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, ab_x, ab_y, ab_z, gp_31, gp_32, gp_33, hp_33, \
                         hp_46, hp_47, hp_50 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = ab_y[k] * gp_31[k]
                  + hp_46[k];

        t_64[k] = ab_y[k] * gp_32[k]
                  + hp_47[k];

        t_65[k] = ab_z[k] * gp_32[k]
                  + hp_50[k];

        t_66[k] = ab_x[k] * gp_33[k]
                  + hp_33[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, t_71, ab_x, ab_y, ab_z, gp_34, gp_35, hp_34, \
                         hp_35, hp_49, hp_50, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_67[k] = ab_x[k] * gp_34[k]
                  + hp_34[k];

        t_68[k] = ab_x[k] * gp_35[k]
                  + hp_35[k];

        t_69[k] = ab_y[k] * gp_34[k]
                  + hp_49[k];

        t_70[k] = ab_y[k] * gp_35[k]
                  + hp_50[k];

        t_71[k] = ab_z[k] * gp_35[k]
                  + hp_53[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, ab_x, ab_y, gp_36, gp_37, gp_38, hp_36, \
                         hp_37, hp_38, hp_52, hp_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_72[k] = ab_x[k] * gp_36[k]
                  + hp_36[k];

        t_73[k] = ab_x[k] * gp_37[k]
                  + hp_37[k];

        t_74[k] = ab_x[k] * gp_38[k]
                  + hp_38[k];

        t_75[k] = ab_y[k] * gp_37[k]
                  + hp_52[k];

        t_76[k] = ab_y[k] * gp_38[k]
                  + hp_53[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, ab_x, ab_z, gp_38, gp_39, gp_40, gp_41, \
                         hp_39, hp_40, hp_41, hp_56 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_77[k] = ab_z[k] * gp_38[k]
                  + hp_56[k];

        t_78[k] = ab_x[k] * gp_39[k]
                  + hp_39[k];

        t_79[k] = ab_x[k] * gp_40[k]
                  + hp_40[k];

        t_80[k] = ab_x[k] * gp_41[k]
                  + hp_41[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, ab_x, ab_y, ab_z, gp_40, gp_41, gp_42, hp_42, \
                         hp_55, hp_56, hp_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_81[k] = ab_y[k] * gp_40[k]
                  + hp_55[k];

        t_82[k] = ab_y[k] * gp_41[k]
                  + hp_56[k];

        t_83[k] = ab_z[k] * gp_41[k]
                  + hp_59[k];

        t_84[k] = ab_x[k] * gp_42[k]
                  + hp_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, t_89, ab_x, ab_y, ab_z, gp_43, gp_44, hp_43, \
                         hp_44, hp_58, hp_59, hp_62 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_85[k] = ab_x[k] * gp_43[k]
                  + hp_43[k];

        t_86[k] = ab_x[k] * gp_44[k]
                  + hp_44[k];

        t_87[k] = ab_y[k] * gp_43[k]
                  + hp_58[k];

        t_88[k] = ab_y[k] * gp_44[k]
                  + hp_59[k];

        t_89[k] = ab_z[k] * gp_44[k]
                  + hp_62[k];
    }
}

}  // namespace simdtrf
