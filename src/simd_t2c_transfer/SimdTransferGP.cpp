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


#include "SimdTransferGP.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_gp(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t gs, const size_t hs, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, ab_x, ab_y, ab_z, gs_0, gs_1, hs_0, \
                         hs_1, hs_2, hs_3, hs_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = ab_x[k] * gs_0[k]
                 + hs_0[k];

        t_1[k] = ab_y[k] * gs_0[k]
                 + hs_1[k];

        t_2[k] = ab_z[k] * gs_0[k]
                 + hs_2[k];

        t_3[k] = ab_x[k] * gs_1[k]
                 + hs_1[k];

        t_4[k] = ab_y[k] * gs_1[k]
                 + hs_3[k];

        t_5[k] = ab_z[k] * gs_1[k]
                 + hs_4[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, ab_x, ab_y, ab_z, gs_2, gs_3, hs_2, hs_3, \
                         hs_4, hs_5, hs_6 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_6[k] = ab_x[k] * gs_2[k]
                 + hs_2[k];

        t_7[k] = ab_y[k] * gs_2[k]
                 + hs_4[k];

        t_8[k] = ab_z[k] * gs_2[k]
                 + hs_5[k];

        t_9[k] = ab_x[k] * gs_3[k]
                 + hs_3[k];

        t_10[k] = ab_y[k] * gs_3[k]
                  + hs_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, t_16, ab_x, ab_y, ab_z, gs_3, gs_4, \
                         gs_5, hs_4, hs_5, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_11[k] = ab_z[k] * gs_3[k]
                  + hs_7[k];

        t_12[k] = ab_x[k] * gs_4[k]
                  + hs_4[k];

        t_13[k] = ab_y[k] * gs_4[k]
                  + hs_7[k];

        t_14[k] = ab_z[k] * gs_4[k]
                  + hs_8[k];

        t_15[k] = ab_x[k] * gs_5[k]
                  + hs_5[k];

        t_16[k] = ab_y[k] * gs_5[k]
                  + hs_8[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, t_21, ab_x, ab_y, ab_z, gs_5, gs_6, gs_7, \
                         hs_6, hs_7, hs_9, hs_10, hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_17[k] = ab_z[k] * gs_5[k]
                  + hs_9[k];

        t_18[k] = ab_x[k] * gs_6[k]
                  + hs_6[k];

        t_19[k] = ab_y[k] * gs_6[k]
                  + hs_10[k];

        t_20[k] = ab_z[k] * gs_6[k]
                  + hs_11[k];

        t_21[k] = ab_x[k] * gs_7[k]
                  + hs_7[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, ab_x, ab_y, ab_z, gs_7, gs_8, hs_8, \
                         hs_11, hs_12, hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_22[k] = ab_y[k] * gs_7[k]
                  + hs_11[k];

        t_23[k] = ab_z[k] * gs_7[k]
                  + hs_12[k];

        t_24[k] = ab_x[k] * gs_8[k]
                  + hs_8[k];

        t_25[k] = ab_y[k] * gs_8[k]
                  + hs_12[k];

        t_26[k] = ab_z[k] * gs_8[k]
                  + hs_13[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, t_31, ab_x, ab_y, ab_z, gs_9, gs_10, hs_9, \
                         hs_10, hs_13, hs_14, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_27[k] = ab_x[k] * gs_9[k]
                  + hs_9[k];

        t_28[k] = ab_y[k] * gs_9[k]
                  + hs_13[k];

        t_29[k] = ab_z[k] * gs_9[k]
                  + hs_14[k];

        t_30[k] = ab_x[k] * gs_10[k]
                  + hs_10[k];

        t_31[k] = ab_y[k] * gs_10[k]
                  + hs_15[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, t_36, t_37, ab_x, ab_y, ab_z, gs_10, gs_11, \
                         gs_12, hs_11, hs_12, hs_16, hs_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_32[k] = ab_z[k] * gs_10[k]
                  + hs_16[k];

        t_33[k] = ab_x[k] * gs_11[k]
                  + hs_11[k];

        t_34[k] = ab_y[k] * gs_11[k]
                  + hs_16[k];

        t_35[k] = ab_z[k] * gs_11[k]
                  + hs_17[k];

        t_36[k] = ab_x[k] * gs_12[k]
                  + hs_12[k];

        t_37[k] = ab_y[k] * gs_12[k]
                  + hs_17[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, t_42, t_43, ab_x, ab_y, ab_z, gs_12, gs_13, \
                         gs_14, hs_13, hs_14, hs_18, hs_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_38[k] = ab_z[k] * gs_12[k]
                  + hs_18[k];

        t_39[k] = ab_x[k] * gs_13[k]
                  + hs_13[k];

        t_40[k] = ab_y[k] * gs_13[k]
                  + hs_18[k];

        t_41[k] = ab_z[k] * gs_13[k]
                  + hs_19[k];

        t_42[k] = ab_x[k] * gs_14[k]
                  + hs_14[k];

        t_43[k] = ab_y[k] * gs_14[k]
                  + hs_19[k];
    }

#pragma omp simd aligned(t_44, ab_z, gs_14, hs_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_44[k] = ab_z[k] * gs_14[k]
                  + hs_20[k];
    }
}

auto
compute_hrr_gp_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t gs, const size_t hs,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(35.0);
    const auto f_1 = 0.75 * std::sqrt(70.0);
    const auto f_2 = 0.25 * std::sqrt(70.0);
    const auto f_3 = 0.5 * std::sqrt(5.0);
    const auto f_4 = 3.0 * std::sqrt(5.0);
    const auto f_5 = 0.75 * std::sqrt(10.0);
    const auto f_6 = std::sqrt(10.0);
    const auto f_7 = 0.25 * std::sqrt(5.0);
    const auto f_8 = 1.5 * std::sqrt(5.0);
    const auto f_9 = 0.125 * std::sqrt(35.0);
    const auto f_10 = 0.75 * std::sqrt(35.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *gs_0 = buffer.data(gs + 0);
    const auto *gs_1 = buffer.data(gs + 1);
    const auto *gs_2 = buffer.data(gs + 2);
    const auto *gs_3 = buffer.data(gs + 3);
    const auto *gs_4 = buffer.data(gs + 4);
    const auto *gs_5 = buffer.data(gs + 5);
    const auto *gs_6 = buffer.data(gs + 6);
    const auto *gs_7 = buffer.data(gs + 7);
    const auto *gs_8 = buffer.data(gs + 8);
    const auto *gs_9 = buffer.data(gs + 9);
    const auto *gs_10 = buffer.data(gs + 10);
    const auto *gs_11 = buffer.data(gs + 11);
    const auto *gs_12 = buffer.data(gs + 12);
    const auto *gs_13 = buffer.data(gs + 13);
    const auto *gs_14 = buffer.data(gs + 14);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_1 = buffer.data(hs + 1);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);
    const auto *hs_17 = buffer.data(hs + 17);
    const auto *hs_18 = buffer.data(hs + 18);
    const auto *hs_19 = buffer.data(hs + 19);
    const auto *hs_20 = buffer.data(hs + 20);

#pragma omp simd aligned(ab_x, ab_y, ab_z, gs_1, gs_6, hs_1, hs_3, hs_4, hs_6, hs_10, \
                         hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * ab_y[k] * gs_1[k]
                 - f_0 * ab_y[k] * gs_6[k]
                 + f_0 * hs_3[k]
                 - f_0 * hs_10[k];

        g_1[k] = f_0 * ab_z[k] * gs_1[k]
                 - f_0 * ab_z[k] * gs_6[k]
                 + f_0 * hs_4[k]
                 - f_0 * hs_11[k];

        g_2[k] = f_0 * ab_x[k] * gs_1[k]
                 - f_0 * ab_x[k] * gs_6[k]
                 + f_0 * hs_1[k]
                 - f_0 * hs_6[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, gs_4, gs_11, hs_4, hs_7, hs_8, hs_11, hs_16, \
                         hs_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_1 * ab_y[k] * gs_4[k]
                 - f_2 * ab_y[k] * gs_11[k]
                 + f_1 * hs_7[k]
                 - f_2 * hs_16[k];

        g_4[k] = f_1 * ab_z[k] * gs_4[k]
                 - f_2 * ab_z[k] * gs_11[k]
                 + f_1 * hs_8[k]
                 - f_2 * hs_17[k];

        g_5[k] = f_1 * ab_x[k] * gs_4[k]
                 - f_2 * ab_x[k] * gs_11[k]
                 + f_1 * hs_4[k]
                 - f_2 * hs_11[k];
    }

#pragma omp simd aligned(ab_y, ab_z, gs_1, gs_6, gs_8, hs_3, hs_4, hs_10, hs_11, hs_12, \
                         hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_3 * ab_y[k] * gs_1[k]
                 - f_3 * ab_y[k] * gs_6[k]
                 + f_4 * ab_y[k] * gs_8[k]
                 - f_3 * hs_3[k]
                 - f_3 * hs_10[k]
                 + f_4 * hs_12[k];

        g_7[k] = -f_3 * ab_z[k] * gs_1[k]
                 - f_3 * ab_z[k] * gs_6[k]
                 + f_4 * ab_z[k] * gs_8[k]
                 - f_3 * hs_4[k]
                 - f_3 * hs_11[k]
                 + f_4 * hs_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gs_1, gs_4, gs_6, gs_8, gs_11, gs_13, hs_1, hs_6, hs_7, \
                         hs_8, hs_16, hs_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_3 * ab_x[k] * gs_1[k]
                 - f_3 * ab_x[k] * gs_6[k]
                 + f_4 * ab_x[k] * gs_8[k]
                 - f_3 * hs_1[k]
                 - f_3 * hs_6[k]
                 + f_4 * hs_8[k];

        g_9[k] = -f_5 * ab_y[k] * gs_4[k]
                 - f_5 * ab_y[k] * gs_11[k]
                 + f_6 * ab_y[k] * gs_13[k]
                 - f_5 * hs_7[k]
                 - f_5 * hs_16[k]
                 + f_6 * hs_18[k];
    }

#pragma omp simd aligned(ab_x, ab_z, gs_4, gs_11, gs_13, hs_4, hs_8, hs_11, hs_13, hs_17, \
                         hs_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_5 * ab_z[k] * gs_4[k]
                  - f_5 * ab_z[k] * gs_11[k]
                  + f_6 * ab_z[k] * gs_13[k]
                  - f_5 * hs_8[k]
                  - f_5 * hs_17[k]
                  + f_6 * hs_19[k];

        g_11[k] = -f_5 * ab_x[k] * gs_4[k]
                  - f_5 * ab_x[k] * gs_11[k]
                  + f_6 * ab_x[k] * gs_13[k]
                  - f_5 * hs_4[k]
                  - f_5 * hs_11[k]
                  + f_6 * hs_13[k];
    }

#pragma omp simd aligned(ab_y, gs_0, gs_3, gs_5, gs_10, gs_12, gs_14, hs_1, hs_6, hs_8, hs_15, \
                         hs_17, hs_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = 0.375 * ab_y[k] * gs_0[k]
                  + 0.75 * ab_y[k] * gs_3[k]
                  - 3.0 * ab_y[k] * gs_5[k]
                  + 0.375 * ab_y[k] * gs_10[k]
                  - 3.0 * ab_y[k] * gs_12[k]
                  + ab_y[k] * gs_14[k]
                  + 0.375 * hs_1[k]
                  + 0.75 * hs_6[k]
                  - 3.0 * hs_8[k]
                  + 0.375 * hs_15[k]
                  - 3.0 * hs_17[k]
                  + hs_19[k];
    }

#pragma omp simd aligned(ab_z, gs_0, gs_3, gs_5, gs_10, gs_12, gs_14, hs_2, hs_7, hs_9, hs_16, \
                         hs_18, hs_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = 0.375 * ab_z[k] * gs_0[k]
                  + 0.75 * ab_z[k] * gs_3[k]
                  - 3.0 * ab_z[k] * gs_5[k]
                  + 0.375 * ab_z[k] * gs_10[k]
                  - 3.0 * ab_z[k] * gs_12[k]
                  + ab_z[k] * gs_14[k]
                  + 0.375 * hs_2[k]
                  + 0.75 * hs_7[k]
                  - 3.0 * hs_9[k]
                  + 0.375 * hs_16[k]
                  - 3.0 * hs_18[k]
                  + hs_20[k];
    }

#pragma omp simd aligned(ab_x, gs_0, gs_3, gs_5, gs_10, gs_12, gs_14, hs_0, hs_3, hs_5, hs_10, \
                         hs_12, hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = 0.375 * ab_x[k] * gs_0[k]
                  + 0.75 * ab_x[k] * gs_3[k]
                  - 3.0 * ab_x[k] * gs_5[k]
                  + 0.375 * ab_x[k] * gs_10[k]
                  - 3.0 * ab_x[k] * gs_12[k]
                  + ab_x[k] * gs_14[k]
                  + 0.375 * hs_0[k]
                  + 0.75 * hs_3[k]
                  - 3.0 * hs_5[k]
                  + 0.375 * hs_10[k]
                  - 3.0 * hs_12[k]
                  + hs_14[k];
    }

#pragma omp simd aligned(ab_y, ab_z, gs_2, gs_7, gs_9, hs_4, hs_5, hs_11, hs_12, hs_13, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_5 * ab_y[k] * gs_2[k]
                  - f_5 * ab_y[k] * gs_7[k]
                  + f_6 * ab_y[k] * gs_9[k]
                  - f_5 * hs_4[k]
                  - f_5 * hs_11[k]
                  + f_6 * hs_13[k];

        g_16[k] = -f_5 * ab_z[k] * gs_2[k]
                  - f_5 * ab_z[k] * gs_7[k]
                  + f_6 * ab_z[k] * gs_9[k]
                  - f_5 * hs_5[k]
                  - f_5 * hs_12[k]
                  + f_6 * hs_14[k];
    }

#pragma omp simd aligned(ab_x, gs_2, gs_7, gs_9, hs_2, hs_7, hs_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_5 * ab_x[k] * gs_2[k]
                  - f_5 * ab_x[k] * gs_7[k]
                  + f_6 * ab_x[k] * gs_9[k]
                  - f_5 * hs_2[k]
                  - f_5 * hs_7[k]
                  + f_6 * hs_9[k];
    }

#pragma omp simd aligned(ab_y, ab_z, gs_0, gs_5, gs_10, gs_12, hs_1, hs_2, hs_8, hs_9, hs_15, \
                         hs_16, hs_17, hs_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_7 * ab_y[k] * gs_0[k]
                  + f_8 * ab_y[k] * gs_5[k]
                  + f_7 * ab_y[k] * gs_10[k]
                  - f_8 * ab_y[k] * gs_12[k]
                  - f_7 * hs_1[k]
                  + f_8 * hs_8[k]
                  + f_7 * hs_15[k]
                  - f_8 * hs_17[k];

        g_19[k] = -f_7 * ab_z[k] * gs_0[k]
                  + f_8 * ab_z[k] * gs_5[k]
                  + f_7 * ab_z[k] * gs_10[k]
                  - f_8 * ab_z[k] * gs_12[k]
                  - f_7 * hs_2[k]
                  + f_8 * hs_9[k]
                  + f_7 * hs_16[k]
                  - f_8 * hs_18[k];
    }

#pragma omp simd aligned(ab_x, ab_y, gs_0, gs_2, gs_5, gs_7, gs_10, gs_12, hs_0, hs_4, hs_5, \
                         hs_10, hs_11, hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_7 * ab_x[k] * gs_0[k]
                  + f_8 * ab_x[k] * gs_5[k]
                  + f_7 * ab_x[k] * gs_10[k]
                  - f_8 * ab_x[k] * gs_12[k]
                  - f_7 * hs_0[k]
                  + f_8 * hs_5[k]
                  + f_7 * hs_10[k]
                  - f_8 * hs_12[k];

        g_21[k] = f_2 * ab_y[k] * gs_2[k]
                  - f_1 * ab_y[k] * gs_7[k]
                  + f_2 * hs_4[k]
                  - f_1 * hs_11[k];
    }

#pragma omp simd aligned(ab_x, ab_z, gs_2, gs_7, hs_2, hs_5, hs_7, \
                         hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_2 * ab_z[k] * gs_2[k]
                  - f_1 * ab_z[k] * gs_7[k]
                  + f_2 * hs_5[k]
                  - f_1 * hs_12[k];

        g_23[k] = f_2 * ab_x[k] * gs_2[k]
                  - f_1 * ab_x[k] * gs_7[k]
                  + f_2 * hs_2[k]
                  - f_1 * hs_7[k];
    }

#pragma omp simd aligned(ab_y, ab_z, gs_0, gs_3, gs_10, hs_1, hs_2, hs_6, hs_7, hs_15, \
                         hs_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_9 * ab_y[k] * gs_0[k]
                  - f_10 * ab_y[k] * gs_3[k]
                  + f_9 * ab_y[k] * gs_10[k]
                  + f_9 * hs_1[k]
                  - f_10 * hs_6[k]
                  + f_9 * hs_15[k];

        g_25[k] = f_9 * ab_z[k] * gs_0[k]
                  - f_10 * ab_z[k] * gs_3[k]
                  + f_9 * ab_z[k] * gs_10[k]
                  + f_9 * hs_2[k]
                  - f_10 * hs_7[k]
                  + f_9 * hs_16[k];
    }

#pragma omp simd aligned(ab_x, gs_0, gs_3, gs_10, hs_0, hs_3, hs_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = f_9 * ab_x[k] * gs_0[k]
                  - f_10 * ab_x[k] * gs_3[k]
                  + f_9 * ab_x[k] * gs_10[k]
                  + f_9 * hs_0[k]
                  - f_10 * hs_3[k]
                  + f_9 * hs_10[k];
    }
}

}  // namespace simdtrf
