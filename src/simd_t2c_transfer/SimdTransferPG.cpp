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


#include "SimdTransferPG.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_pg_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t sg, const size_t sh,
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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

#pragma omp simd aligned(ab_y, sg_1, sg_4, sg_6, sg_8, sg_11, sh_3, sh_7, sh_10, sh_12, \
                         sh_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_y[k] * sg_1[k]
                 + f_0 * ab_y[k] * sg_6[k]
                 + f_0 * sh_3[k]
                 - f_0 * sh_10[k];

        g_1[k] = -f_1 * ab_y[k] * sg_4[k]
                 + f_2 * ab_y[k] * sg_11[k]
                 + f_1 * sh_7[k]
                 - f_2 * sh_16[k];

        g_2[k] = f_3 * ab_y[k] * sg_1[k]
                 + f_3 * ab_y[k] * sg_6[k]
                 - f_4 * ab_y[k] * sg_8[k]
                 - f_3 * sh_3[k]
                 - f_3 * sh_10[k]
                 + f_4 * sh_12[k];
    }

#pragma omp simd aligned(ab_y, sg_4, sg_11, sg_13, sh_7, sh_16, sh_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_5 * ab_y[k] * sg_4[k]
                 + f_5 * ab_y[k] * sg_11[k]
                 - f_6 * ab_y[k] * sg_13[k]
                 - f_5 * sh_7[k]
                 - f_5 * sh_16[k]
                 + f_6 * sh_18[k];
    }

#pragma omp simd aligned(ab_y, sg_0, sg_3, sg_5, sg_10, sg_12, sg_14, sh_1, sh_6, sh_8, sh_15, \
                         sh_17, sh_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -0.375 * ab_y[k] * sg_0[k]
                 - 0.75 * ab_y[k] * sg_3[k]
                 + 3.0 * ab_y[k] * sg_5[k]
                 - 0.375 * ab_y[k] * sg_10[k]
                 + 3.0 * ab_y[k] * sg_12[k]
                 - ab_y[k] * sg_14[k]
                 + 0.375 * sh_1[k]
                 + 0.75 * sh_6[k]
                 - 3.0 * sh_8[k]
                 + 0.375 * sh_15[k]
                 - 3.0 * sh_17[k]
                 + sh_19[k];
    }

#pragma omp simd aligned(ab_y, sg_2, sg_7, sg_9, sh_4, sh_11, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = f_5 * ab_y[k] * sg_2[k]
                 + f_5 * ab_y[k] * sg_7[k]
                 - f_6 * ab_y[k] * sg_9[k]
                 - f_5 * sh_4[k]
                 - f_5 * sh_11[k]
                 + f_6 * sh_13[k];
    }

#pragma omp simd aligned(ab_y, sg_0, sg_2, sg_5, sg_7, sg_10, sg_12, sh_1, sh_4, sh_8, sh_11, \
                         sh_15, sh_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = f_7 * ab_y[k] * sg_0[k]
                 - f_8 * ab_y[k] * sg_5[k]
                 - f_7 * ab_y[k] * sg_10[k]
                 + f_8 * ab_y[k] * sg_12[k]
                 - f_7 * sh_1[k]
                 + f_8 * sh_8[k]
                 + f_7 * sh_15[k]
                 - f_8 * sh_17[k];

        g_7[k] = -f_2 * ab_y[k] * sg_2[k]
                 + f_1 * ab_y[k] * sg_7[k]
                 + f_2 * sh_4[k]
                 - f_1 * sh_11[k];
    }

#pragma omp simd aligned(ab_y, ab_z, sg_0, sg_1, sg_3, sg_6, sg_10, sh_1, sh_4, sh_6, sh_11, \
                         sh_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_9 * ab_y[k] * sg_0[k]
                 + f_10 * ab_y[k] * sg_3[k]
                 - f_9 * ab_y[k] * sg_10[k]
                 + f_9 * sh_1[k]
                 - f_10 * sh_6[k]
                 + f_9 * sh_15[k];

        g_9[k] = -f_0 * ab_z[k] * sg_1[k]
                 + f_0 * ab_z[k] * sg_6[k]
                 + f_0 * sh_4[k]
                 - f_0 * sh_11[k];
    }

#pragma omp simd aligned(ab_z, sg_1, sg_4, sg_6, sg_8, sg_11, sg_13, sh_4, sh_8, sh_11, sh_13, \
                         sh_17, sh_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_1 * ab_z[k] * sg_4[k]
                  + f_2 * ab_z[k] * sg_11[k]
                  + f_1 * sh_8[k]
                  - f_2 * sh_17[k];

        g_11[k] = f_3 * ab_z[k] * sg_1[k]
                  + f_3 * ab_z[k] * sg_6[k]
                  - f_4 * ab_z[k] * sg_8[k]
                  - f_3 * sh_4[k]
                  - f_3 * sh_11[k]
                  + f_4 * sh_13[k];

        g_12[k] = f_5 * ab_z[k] * sg_4[k]
                  + f_5 * ab_z[k] * sg_11[k]
                  - f_6 * ab_z[k] * sg_13[k]
                  - f_5 * sh_8[k]
                  - f_5 * sh_17[k]
                  + f_6 * sh_19[k];
    }

#pragma omp simd aligned(ab_z, sg_0, sg_3, sg_5, sg_10, sg_12, sg_14, sh_2, sh_7, sh_9, sh_16, \
                         sh_18, sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -0.375 * ab_z[k] * sg_0[k]
                  - 0.75 * ab_z[k] * sg_3[k]
                  + 3.0 * ab_z[k] * sg_5[k]
                  - 0.375 * ab_z[k] * sg_10[k]
                  + 3.0 * ab_z[k] * sg_12[k]
                  - ab_z[k] * sg_14[k]
                  + 0.375 * sh_2[k]
                  + 0.75 * sh_7[k]
                  - 3.0 * sh_9[k]
                  + 0.375 * sh_16[k]
                  - 3.0 * sh_18[k]
                  + sh_20[k];
    }

#pragma omp simd aligned(ab_z, sg_2, sg_7, sg_9, sh_5, sh_12, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_5 * ab_z[k] * sg_2[k]
                  + f_5 * ab_z[k] * sg_7[k]
                  - f_6 * ab_z[k] * sg_9[k]
                  - f_5 * sh_5[k]
                  - f_5 * sh_12[k]
                  + f_6 * sh_14[k];
    }

#pragma omp simd aligned(ab_z, sg_0, sg_2, sg_5, sg_7, sg_10, sg_12, sh_2, sh_5, sh_9, sh_12, \
                         sh_16, sh_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_7 * ab_z[k] * sg_0[k]
                  - f_8 * ab_z[k] * sg_5[k]
                  - f_7 * ab_z[k] * sg_10[k]
                  + f_8 * ab_z[k] * sg_12[k]
                  - f_7 * sh_2[k]
                  + f_8 * sh_9[k]
                  + f_7 * sh_16[k]
                  - f_8 * sh_18[k];

        g_16[k] = -f_2 * ab_z[k] * sg_2[k]
                  + f_1 * ab_z[k] * sg_7[k]
                  + f_2 * sh_5[k]
                  - f_1 * sh_12[k];
    }

#pragma omp simd aligned(ab_x, ab_z, sg_0, sg_1, sg_3, sg_6, sg_10, sh_1, sh_2, sh_6, sh_7, \
                         sh_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_9 * ab_z[k] * sg_0[k]
                  + f_10 * ab_z[k] * sg_3[k]
                  - f_9 * ab_z[k] * sg_10[k]
                  + f_9 * sh_2[k]
                  - f_10 * sh_7[k]
                  + f_9 * sh_16[k];

        g_18[k] = -f_0 * ab_x[k] * sg_1[k]
                  + f_0 * ab_x[k] * sg_6[k]
                  + f_0 * sh_1[k]
                  - f_0 * sh_6[k];
    }

#pragma omp simd aligned(ab_x, sg_1, sg_4, sg_6, sg_8, sg_11, sg_13, sh_1, sh_4, sh_6, sh_8, \
                         sh_11, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = -f_1 * ab_x[k] * sg_4[k]
                  + f_2 * ab_x[k] * sg_11[k]
                  + f_1 * sh_4[k]
                  - f_2 * sh_11[k];

        g_20[k] = f_3 * ab_x[k] * sg_1[k]
                  + f_3 * ab_x[k] * sg_6[k]
                  - f_4 * ab_x[k] * sg_8[k]
                  - f_3 * sh_1[k]
                  - f_3 * sh_6[k]
                  + f_4 * sh_8[k];

        g_21[k] = f_5 * ab_x[k] * sg_4[k]
                  + f_5 * ab_x[k] * sg_11[k]
                  - f_6 * ab_x[k] * sg_13[k]
                  - f_5 * sh_4[k]
                  - f_5 * sh_11[k]
                  + f_6 * sh_13[k];
    }

#pragma omp simd aligned(ab_x, sg_0, sg_3, sg_5, sg_10, sg_12, sg_14, sh_0, sh_3, sh_5, sh_10, \
                         sh_12, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -0.375 * ab_x[k] * sg_0[k]
                  - 0.75 * ab_x[k] * sg_3[k]
                  + 3.0 * ab_x[k] * sg_5[k]
                  - 0.375 * ab_x[k] * sg_10[k]
                  + 3.0 * ab_x[k] * sg_12[k]
                  - ab_x[k] * sg_14[k]
                  + 0.375 * sh_0[k]
                  + 0.75 * sh_3[k]
                  - 3.0 * sh_5[k]
                  + 0.375 * sh_10[k]
                  - 3.0 * sh_12[k]
                  + sh_14[k];
    }

#pragma omp simd aligned(ab_x, sg_2, sg_7, sg_9, sh_2, sh_7, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_5 * ab_x[k] * sg_2[k]
                  + f_5 * ab_x[k] * sg_7[k]
                  - f_6 * ab_x[k] * sg_9[k]
                  - f_5 * sh_2[k]
                  - f_5 * sh_7[k]
                  + f_6 * sh_9[k];
    }

#pragma omp simd aligned(ab_x, sg_0, sg_2, sg_5, sg_7, sg_10, sg_12, sh_0, sh_2, sh_5, sh_7, \
                         sh_10, sh_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_7 * ab_x[k] * sg_0[k]
                  - f_8 * ab_x[k] * sg_5[k]
                  - f_7 * ab_x[k] * sg_10[k]
                  + f_8 * ab_x[k] * sg_12[k]
                  - f_7 * sh_0[k]
                  + f_8 * sh_5[k]
                  + f_7 * sh_10[k]
                  - f_8 * sh_12[k];

        g_25[k] = -f_2 * ab_x[k] * sg_2[k]
                  + f_1 * ab_x[k] * sg_7[k]
                  + f_2 * sh_2[k]
                  - f_1 * sh_7[k];
    }

#pragma omp simd aligned(ab_x, sg_0, sg_3, sg_10, sh_0, sh_3, sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_9 * ab_x[k] * sg_0[k]
                  + f_10 * ab_x[k] * sg_3[k]
                  - f_9 * ab_x[k] * sg_10[k]
                  + f_9 * sh_0[k]
                  - f_10 * sh_3[k]
                  + f_9 * sh_10[k];
    }
}

auto
compute_hrr_pg(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sg, const size_t sh, const size_t nmax) -> void
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

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_1 = buffer.data(sg + 1);
    const auto *sg_2 = buffer.data(sg + 2);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_4 = buffer.data(sg + 4);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_7 = buffer.data(sg + 7);
    const auto *sg_8 = buffer.data(sg + 8);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_11 = buffer.data(sg + 11);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_1 = buffer.data(sh + 1);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_4 = buffer.data(sh + 4);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_7 = buffer.data(sh + 7);
    const auto *sh_8 = buffer.data(sh + 8);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_11 = buffer.data(sh + 11);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_13 = buffer.data(sh + 13);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sg_0, sg_1, sg_2, sg_3, sg_4, sh_0, \
                         sh_1, sh_2, sh_3, sh_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sg_0[k]
                 + sh_0[k];

        t_1[k] = -ab_x[k] * sg_1[k]
                 + sh_1[k];

        t_2[k] = -ab_x[k] * sg_2[k]
                 + sh_2[k];

        t_3[k] = -ab_x[k] * sg_3[k]
                 + sh_3[k];

        t_4[k] = -ab_x[k] * sg_4[k]
                 + sh_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sg_5, sg_6, sg_7, sg_8, sg_9, sh_5, \
                         sh_6, sh_7, sh_8, sh_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sg_5[k]
                 + sh_5[k];

        t_6[k] = -ab_x[k] * sg_6[k]
                 + sh_6[k];

        t_7[k] = -ab_x[k] * sg_7[k]
                 + sh_7[k];

        t_8[k] = -ab_x[k] * sg_8[k]
                 + sh_8[k];

        t_9[k] = -ab_x[k] * sg_9[k]
                 + sh_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, sg_10, sg_11, sg_12, sg_13, \
                         sg_14, sh_10, sh_11, sh_12, sh_13, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * sg_10[k]
                  + sh_10[k];

        t_11[k] = -ab_x[k] * sg_11[k]
                  + sh_11[k];

        t_12[k] = -ab_x[k] * sg_12[k]
                  + sh_12[k];

        t_13[k] = -ab_x[k] * sg_13[k]
                  + sh_13[k];

        t_14[k] = -ab_x[k] * sg_14[k]
                  + sh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_y, sg_0, sg_1, sg_2, sg_3, sg_4, \
                         sh_1, sh_3, sh_4, sh_6, sh_7 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_y[k] * sg_0[k]
                  + sh_1[k];

        t_16[k] = -ab_y[k] * sg_1[k]
                  + sh_3[k];

        t_17[k] = -ab_y[k] * sg_2[k]
                  + sh_4[k];

        t_18[k] = -ab_y[k] * sg_3[k]
                  + sh_6[k];

        t_19[k] = -ab_y[k] * sg_4[k]
                  + sh_7[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_y, sg_5, sg_6, sg_7, sg_8, sg_9, \
                         sh_8, sh_10, sh_11, sh_12, sh_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_y[k] * sg_5[k]
                  + sh_8[k];

        t_21[k] = -ab_y[k] * sg_6[k]
                  + sh_10[k];

        t_22[k] = -ab_y[k] * sg_7[k]
                  + sh_11[k];

        t_23[k] = -ab_y[k] * sg_8[k]
                  + sh_12[k];

        t_24[k] = -ab_y[k] * sg_9[k]
                  + sh_13[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_y, sg_10, sg_11, sg_12, sg_13, \
                         sg_14, sh_15, sh_16, sh_17, sh_18, sh_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_y[k] * sg_10[k]
                  + sh_15[k];

        t_26[k] = -ab_y[k] * sg_11[k]
                  + sh_16[k];

        t_27[k] = -ab_y[k] * sg_12[k]
                  + sh_17[k];

        t_28[k] = -ab_y[k] * sg_13[k]
                  + sh_18[k];

        t_29[k] = -ab_y[k] * sg_14[k]
                  + sh_19[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, ab_z, sg_0, sg_1, sg_2, sg_3, sg_4, \
                         sh_2, sh_4, sh_5, sh_7, sh_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_30[k] = -ab_z[k] * sg_0[k]
                  + sh_2[k];

        t_31[k] = -ab_z[k] * sg_1[k]
                  + sh_4[k];

        t_32[k] = -ab_z[k] * sg_2[k]
                  + sh_5[k];

        t_33[k] = -ab_z[k] * sg_3[k]
                  + sh_7[k];

        t_34[k] = -ab_z[k] * sg_4[k]
                  + sh_8[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, ab_z, sg_5, sg_6, sg_7, sg_8, sg_9, \
                         sh_9, sh_11, sh_12, sh_13, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_35[k] = -ab_z[k] * sg_5[k]
                  + sh_9[k];

        t_36[k] = -ab_z[k] * sg_6[k]
                  + sh_11[k];

        t_37[k] = -ab_z[k] * sg_7[k]
                  + sh_12[k];

        t_38[k] = -ab_z[k] * sg_8[k]
                  + sh_13[k];

        t_39[k] = -ab_z[k] * sg_9[k]
                  + sh_14[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, ab_z, sg_10, sg_11, sg_12, sg_13, \
                         sg_14, sh_16, sh_17, sh_18, sh_19, sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_40[k] = -ab_z[k] * sg_10[k]
                  + sh_16[k];

        t_41[k] = -ab_z[k] * sg_11[k]
                  + sh_17[k];

        t_42[k] = -ab_z[k] * sg_12[k]
                  + sh_18[k];

        t_43[k] = -ab_z[k] * sg_13[k]
                  + sh_19[k];

        t_44[k] = -ab_z[k] * sg_14[k]
                  + sh_20[k];
    }
}

}  // namespace simdtrf
