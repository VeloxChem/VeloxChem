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


#include "SimdTransferPH.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_ph_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t sh, const size_t si,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.9375 * std::sqrt(14.0);
    const auto f_1 = 1.875 * std::sqrt(14.0);
    const auto f_2 = 0.1875 * std::sqrt(14.0);
    const auto f_3 = 1.5 * std::sqrt(35.0);
    const auto f_4 = 0.1875 * std::sqrt(70.0);
    const auto f_5 = 0.125 * std::sqrt(70.0);
    const auto f_6 = 1.5 * std::sqrt(70.0);
    const auto f_7 = 0.0625 * std::sqrt(70.0);
    const auto f_8 = 0.5 * std::sqrt(70.0);
    const auto f_9 = 0.5 * std::sqrt(105.0);
    const auto f_10 = std::sqrt(105.0);
    const auto f_11 = 0.125 * std::sqrt(15.0);
    const auto f_12 = 0.25 * std::sqrt(15.0);
    const auto f_13 = 1.5 * std::sqrt(15.0);
    const auto f_14 = std::sqrt(15.0);
    const auto f_15 = 0.25 * std::sqrt(105.0);
    const auto f_16 = 0.375 * std::sqrt(35.0);
    const auto f_17 = 2.25 * std::sqrt(35.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

#pragma omp simd aligned(ab_y, sh_1, sh_4, sh_6, sh_11, sh_15, si_3, si_7, si_10, si_16, \
                         si_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_y[k] * sh_1[k]
                 + f_1 * ab_y[k] * sh_6[k]
                 - f_2 * ab_y[k] * sh_15[k]
                 + f_0 * si_3[k]
                 - f_1 * si_10[k]
                 + f_2 * si_21[k];

        g_1[k] = -f_3 * ab_y[k] * sh_4[k]
                 + f_3 * ab_y[k] * sh_11[k]
                 + f_3 * si_7[k]
                 - f_3 * si_16[k];
    }

#pragma omp simd aligned(ab_y, sh_1, sh_6, sh_8, sh_15, sh_17, si_3, si_10, si_12, si_21, \
                         si_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_4 * ab_y[k] * sh_1[k]
                 + f_5 * ab_y[k] * sh_6[k]
                 - f_6 * ab_y[k] * sh_8[k]
                 - f_7 * ab_y[k] * sh_15[k]
                 + f_8 * ab_y[k] * sh_17[k]
                 - f_4 * si_3[k]
                 - f_5 * si_10[k]
                 + f_6 * si_12[k]
                 + f_7 * si_21[k]
                 - f_8 * si_23[k];
    }

#pragma omp simd aligned(ab_y, sh_4, sh_11, sh_13, si_7, si_16, si_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_9 * ab_y[k] * sh_4[k]
                 + f_9 * ab_y[k] * sh_11[k]
                 - f_10 * ab_y[k] * sh_13[k]
                 - f_9 * si_7[k]
                 - f_9 * si_16[k]
                 + f_10 * si_18[k];
    }

#pragma omp simd aligned(ab_y, sh_1, sh_6, sh_8, sh_15, sh_17, sh_19, si_3, si_10, si_12, \
                         si_21, si_23, si_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_11 * ab_y[k] * sh_1[k]
                 - f_12 * ab_y[k] * sh_6[k]
                 + f_13 * ab_y[k] * sh_8[k]
                 - f_11 * ab_y[k] * sh_15[k]
                 + f_13 * ab_y[k] * sh_17[k]
                 - f_14 * ab_y[k] * sh_19[k]
                 + f_11 * si_3[k]
                 + f_12 * si_10[k]
                 - f_13 * si_12[k]
                 + f_11 * si_21[k]
                 - f_13 * si_23[k]
                 + f_14 * si_25[k];
    }

#pragma omp simd aligned(ab_y, sh_2, sh_7, sh_9, sh_16, sh_18, sh_20, si_4, si_11, si_13, \
                         si_22, si_24, si_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -1.875 * ab_y[k] * sh_2[k]
                 - 3.75 * ab_y[k] * sh_7[k]
                 + 5.0 * ab_y[k] * sh_9[k]
                 - 1.875 * ab_y[k] * sh_16[k]
                 + 5.0 * ab_y[k] * sh_18[k]
                 - ab_y[k] * sh_20[k]
                 + 1.875 * si_4[k]
                 + 3.75 * si_11[k]
                 - 5.0 * si_13[k]
                 + 1.875 * si_22[k]
                 - 5.0 * si_24[k]
                 + si_26[k];
    }

#pragma omp simd aligned(ab_y, sh_0, sh_3, sh_5, sh_10, sh_12, sh_14, si_1, si_6, si_8, si_15, \
                         si_17, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = -f_11 * ab_y[k] * sh_0[k]
                 - f_12 * ab_y[k] * sh_3[k]
                 + f_13 * ab_y[k] * sh_5[k]
                 - f_11 * ab_y[k] * sh_10[k]
                 + f_13 * ab_y[k] * sh_12[k]
                 - f_14 * ab_y[k] * sh_14[k]
                 + f_11 * si_1[k]
                 + f_12 * si_6[k]
                 - f_13 * si_8[k]
                 + f_11 * si_15[k]
                 - f_13 * si_17[k]
                 + f_14 * si_19[k];
    }

#pragma omp simd aligned(ab_y, sh_2, sh_9, sh_16, sh_18, si_4, si_13, si_22, \
                         si_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_15 * ab_y[k] * sh_2[k]
                 - f_9 * ab_y[k] * sh_9[k]
                 - f_15 * ab_y[k] * sh_16[k]
                 + f_9 * ab_y[k] * sh_18[k]
                 - f_15 * si_4[k]
                 + f_9 * si_13[k]
                 + f_15 * si_22[k]
                 - f_9 * si_24[k];
    }

#pragma omp simd aligned(ab_y, sh_0, sh_3, sh_5, sh_10, sh_12, si_1, si_6, si_8, si_15, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = f_7 * ab_y[k] * sh_0[k]
                 - f_5 * ab_y[k] * sh_3[k]
                 - f_8 * ab_y[k] * sh_5[k]
                 - f_4 * ab_y[k] * sh_10[k]
                 + f_6 * ab_y[k] * sh_12[k]
                 - f_7 * si_1[k]
                 + f_5 * si_6[k]
                 + f_8 * si_8[k]
                 + f_4 * si_15[k]
                 - f_6 * si_17[k];
    }

#pragma omp simd aligned(ab_y, sh_0, sh_2, sh_3, sh_7, sh_10, sh_16, si_1, si_4, si_6, si_11, \
                         si_15, si_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = -f_16 * ab_y[k] * sh_2[k]
                 + f_17 * ab_y[k] * sh_7[k]
                 - f_16 * ab_y[k] * sh_16[k]
                 + f_16 * si_4[k]
                 - f_17 * si_11[k]
                 + f_16 * si_22[k];

        g_10[k] = -f_2 * ab_y[k] * sh_0[k]
                  + f_1 * ab_y[k] * sh_3[k]
                  - f_0 * ab_y[k] * sh_10[k]
                  + f_2 * si_1[k]
                  - f_1 * si_6[k]
                  + f_0 * si_15[k];
    }

#pragma omp simd aligned(ab_z, sh_1, sh_4, sh_6, sh_11, sh_15, si_4, si_8, si_11, si_17, \
                         si_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_0 * ab_z[k] * sh_1[k]
                  + f_1 * ab_z[k] * sh_6[k]
                  - f_2 * ab_z[k] * sh_15[k]
                  + f_0 * si_4[k]
                  - f_1 * si_11[k]
                  + f_2 * si_22[k];

        g_12[k] = -f_3 * ab_z[k] * sh_4[k]
                  + f_3 * ab_z[k] * sh_11[k]
                  + f_3 * si_8[k]
                  - f_3 * si_17[k];
    }

#pragma omp simd aligned(ab_z, sh_1, sh_6, sh_8, sh_15, sh_17, si_4, si_11, si_13, si_22, \
                         si_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_4 * ab_z[k] * sh_1[k]
                  + f_5 * ab_z[k] * sh_6[k]
                  - f_6 * ab_z[k] * sh_8[k]
                  - f_7 * ab_z[k] * sh_15[k]
                  + f_8 * ab_z[k] * sh_17[k]
                  - f_4 * si_4[k]
                  - f_5 * si_11[k]
                  + f_6 * si_13[k]
                  + f_7 * si_22[k]
                  - f_8 * si_24[k];
    }

#pragma omp simd aligned(ab_z, sh_4, sh_11, sh_13, si_8, si_17, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_9 * ab_z[k] * sh_4[k]
                  + f_9 * ab_z[k] * sh_11[k]
                  - f_10 * ab_z[k] * sh_13[k]
                  - f_9 * si_8[k]
                  - f_9 * si_17[k]
                  + f_10 * si_19[k];
    }

#pragma omp simd aligned(ab_z, sh_1, sh_6, sh_8, sh_15, sh_17, sh_19, si_4, si_11, si_13, \
                         si_22, si_24, si_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_11 * ab_z[k] * sh_1[k]
                  - f_12 * ab_z[k] * sh_6[k]
                  + f_13 * ab_z[k] * sh_8[k]
                  - f_11 * ab_z[k] * sh_15[k]
                  + f_13 * ab_z[k] * sh_17[k]
                  - f_14 * ab_z[k] * sh_19[k]
                  + f_11 * si_4[k]
                  + f_12 * si_11[k]
                  - f_13 * si_13[k]
                  + f_11 * si_22[k]
                  - f_13 * si_24[k]
                  + f_14 * si_26[k];
    }

#pragma omp simd aligned(ab_z, sh_2, sh_7, sh_9, sh_16, sh_18, sh_20, si_5, si_12, si_14, \
                         si_23, si_25, si_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = -1.875 * ab_z[k] * sh_2[k]
                  - 3.75 * ab_z[k] * sh_7[k]
                  + 5.0 * ab_z[k] * sh_9[k]
                  - 1.875 * ab_z[k] * sh_16[k]
                  + 5.0 * ab_z[k] * sh_18[k]
                  - ab_z[k] * sh_20[k]
                  + 1.875 * si_5[k]
                  + 3.75 * si_12[k]
                  - 5.0 * si_14[k]
                  + 1.875 * si_23[k]
                  - 5.0 * si_25[k]
                  + si_27[k];
    }

#pragma omp simd aligned(ab_z, sh_0, sh_3, sh_5, sh_10, sh_12, sh_14, si_2, si_7, si_9, si_16, \
                         si_18, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_11 * ab_z[k] * sh_0[k]
                  - f_12 * ab_z[k] * sh_3[k]
                  + f_13 * ab_z[k] * sh_5[k]
                  - f_11 * ab_z[k] * sh_10[k]
                  + f_13 * ab_z[k] * sh_12[k]
                  - f_14 * ab_z[k] * sh_14[k]
                  + f_11 * si_2[k]
                  + f_12 * si_7[k]
                  - f_13 * si_9[k]
                  + f_11 * si_16[k]
                  - f_13 * si_18[k]
                  + f_14 * si_20[k];
    }

#pragma omp simd aligned(ab_z, sh_2, sh_9, sh_16, sh_18, si_5, si_14, si_23, \
                         si_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = f_15 * ab_z[k] * sh_2[k]
                  - f_9 * ab_z[k] * sh_9[k]
                  - f_15 * ab_z[k] * sh_16[k]
                  + f_9 * ab_z[k] * sh_18[k]
                  - f_15 * si_5[k]
                  + f_9 * si_14[k]
                  + f_15 * si_23[k]
                  - f_9 * si_25[k];
    }

#pragma omp simd aligned(ab_z, sh_0, sh_3, sh_5, sh_10, sh_12, si_2, si_7, si_9, si_16, \
                         si_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = f_7 * ab_z[k] * sh_0[k]
                  - f_5 * ab_z[k] * sh_3[k]
                  - f_8 * ab_z[k] * sh_5[k]
                  - f_4 * ab_z[k] * sh_10[k]
                  + f_6 * ab_z[k] * sh_12[k]
                  - f_7 * si_2[k]
                  + f_5 * si_7[k]
                  + f_8 * si_9[k]
                  + f_4 * si_16[k]
                  - f_6 * si_18[k];
    }

#pragma omp simd aligned(ab_z, sh_0, sh_2, sh_3, sh_7, sh_10, sh_16, si_2, si_5, si_7, si_12, \
                         si_16, si_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_16 * ab_z[k] * sh_2[k]
                  + f_17 * ab_z[k] * sh_7[k]
                  - f_16 * ab_z[k] * sh_16[k]
                  + f_16 * si_5[k]
                  - f_17 * si_12[k]
                  + f_16 * si_23[k];

        g_21[k] = -f_2 * ab_z[k] * sh_0[k]
                  + f_1 * ab_z[k] * sh_3[k]
                  - f_0 * ab_z[k] * sh_10[k]
                  + f_2 * si_2[k]
                  - f_1 * si_7[k]
                  + f_0 * si_16[k];
    }

#pragma omp simd aligned(ab_x, sh_1, sh_4, sh_6, sh_11, sh_15, si_1, si_4, si_6, si_11, \
                         si_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = -f_0 * ab_x[k] * sh_1[k]
                  + f_1 * ab_x[k] * sh_6[k]
                  - f_2 * ab_x[k] * sh_15[k]
                  + f_0 * si_1[k]
                  - f_1 * si_6[k]
                  + f_2 * si_15[k];

        g_23[k] = -f_3 * ab_x[k] * sh_4[k]
                  + f_3 * ab_x[k] * sh_11[k]
                  + f_3 * si_4[k]
                  - f_3 * si_11[k];
    }

#pragma omp simd aligned(ab_x, sh_1, sh_6, sh_8, sh_15, sh_17, si_1, si_6, si_8, si_15, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = f_4 * ab_x[k] * sh_1[k]
                  + f_5 * ab_x[k] * sh_6[k]
                  - f_6 * ab_x[k] * sh_8[k]
                  - f_7 * ab_x[k] * sh_15[k]
                  + f_8 * ab_x[k] * sh_17[k]
                  - f_4 * si_1[k]
                  - f_5 * si_6[k]
                  + f_6 * si_8[k]
                  + f_7 * si_15[k]
                  - f_8 * si_17[k];
    }

#pragma omp simd aligned(ab_x, sh_4, sh_11, sh_13, si_4, si_11, si_13 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_9 * ab_x[k] * sh_4[k]
                  + f_9 * ab_x[k] * sh_11[k]
                  - f_10 * ab_x[k] * sh_13[k]
                  - f_9 * si_4[k]
                  - f_9 * si_11[k]
                  + f_10 * si_13[k];
    }

#pragma omp simd aligned(ab_x, sh_1, sh_6, sh_8, sh_15, sh_17, sh_19, si_1, si_6, si_8, si_15, \
                         si_17, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_11 * ab_x[k] * sh_1[k]
                  - f_12 * ab_x[k] * sh_6[k]
                  + f_13 * ab_x[k] * sh_8[k]
                  - f_11 * ab_x[k] * sh_15[k]
                  + f_13 * ab_x[k] * sh_17[k]
                  - f_14 * ab_x[k] * sh_19[k]
                  + f_11 * si_1[k]
                  + f_12 * si_6[k]
                  - f_13 * si_8[k]
                  + f_11 * si_15[k]
                  - f_13 * si_17[k]
                  + f_14 * si_19[k];
    }

#pragma omp simd aligned(ab_x, sh_2, sh_7, sh_9, sh_16, sh_18, sh_20, si_2, si_7, si_9, si_16, \
                         si_18, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_27[k] = -1.875 * ab_x[k] * sh_2[k]
                  - 3.75 * ab_x[k] * sh_7[k]
                  + 5.0 * ab_x[k] * sh_9[k]
                  - 1.875 * ab_x[k] * sh_16[k]
                  + 5.0 * ab_x[k] * sh_18[k]
                  - ab_x[k] * sh_20[k]
                  + 1.875 * si_2[k]
                  + 3.75 * si_7[k]
                  - 5.0 * si_9[k]
                  + 1.875 * si_16[k]
                  - 5.0 * si_18[k]
                  + si_20[k];
    }

#pragma omp simd aligned(ab_x, sh_0, sh_3, sh_5, sh_10, sh_12, sh_14, si_0, si_3, si_5, si_10, \
                         si_12, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = -f_11 * ab_x[k] * sh_0[k]
                  - f_12 * ab_x[k] * sh_3[k]
                  + f_13 * ab_x[k] * sh_5[k]
                  - f_11 * ab_x[k] * sh_10[k]
                  + f_13 * ab_x[k] * sh_12[k]
                  - f_14 * ab_x[k] * sh_14[k]
                  + f_11 * si_0[k]
                  + f_12 * si_3[k]
                  - f_13 * si_5[k]
                  + f_11 * si_10[k]
                  - f_13 * si_12[k]
                  + f_14 * si_14[k];
    }

#pragma omp simd aligned(ab_x, sh_2, sh_9, sh_16, sh_18, si_2, si_9, si_16, \
                         si_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_15 * ab_x[k] * sh_2[k]
                  - f_9 * ab_x[k] * sh_9[k]
                  - f_15 * ab_x[k] * sh_16[k]
                  + f_9 * ab_x[k] * sh_18[k]
                  - f_15 * si_2[k]
                  + f_9 * si_9[k]
                  + f_15 * si_16[k]
                  - f_9 * si_18[k];
    }

#pragma omp simd aligned(ab_x, sh_0, sh_3, sh_5, sh_10, sh_12, si_0, si_3, si_5, si_10, \
                         si_12 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = f_7 * ab_x[k] * sh_0[k]
                  - f_5 * ab_x[k] * sh_3[k]
                  - f_8 * ab_x[k] * sh_5[k]
                  - f_4 * ab_x[k] * sh_10[k]
                  + f_6 * ab_x[k] * sh_12[k]
                  - f_7 * si_0[k]
                  + f_5 * si_3[k]
                  + f_8 * si_5[k]
                  + f_4 * si_10[k]
                  - f_6 * si_12[k];
    }

#pragma omp simd aligned(ab_x, sh_0, sh_2, sh_3, sh_7, sh_10, sh_16, si_0, si_2, si_3, si_7, \
                         si_10, si_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_16 * ab_x[k] * sh_2[k]
                  + f_17 * ab_x[k] * sh_7[k]
                  - f_16 * ab_x[k] * sh_16[k]
                  + f_16 * si_2[k]
                  - f_17 * si_7[k]
                  + f_16 * si_16[k];

        g_32[k] = -f_2 * ab_x[k] * sh_0[k]
                  + f_1 * ab_x[k] * sh_3[k]
                  - f_0 * ab_x[k] * sh_10[k]
                  + f_2 * si_0[k]
                  - f_1 * si_3[k]
                  + f_0 * si_10[k];
    }
}

auto
compute_hrr_ph(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t sh, const size_t si, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);
    const auto *si_3 = buffer.data(si + 3);
    const auto *si_4 = buffer.data(si + 4);
    const auto *si_5 = buffer.data(si + 5);
    const auto *si_6 = buffer.data(si + 6);
    const auto *si_7 = buffer.data(si + 7);
    const auto *si_8 = buffer.data(si + 8);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_10 = buffer.data(si + 10);
    const auto *si_11 = buffer.data(si + 11);
    const auto *si_12 = buffer.data(si + 12);
    const auto *si_13 = buffer.data(si + 13);
    const auto *si_14 = buffer.data(si + 14);
    const auto *si_15 = buffer.data(si + 15);
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);
    const auto *si_21 = buffer.data(si + 21);
    const auto *si_22 = buffer.data(si + 22);
    const auto *si_23 = buffer.data(si + 23);
    const auto *si_24 = buffer.data(si + 24);
    const auto *si_25 = buffer.data(si + 25);
    const auto *si_26 = buffer.data(si + 26);
    const auto *si_27 = buffer.data(si + 27);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, sh_0, sh_1, sh_2, sh_3, sh_4, si_0, \
                         si_1, si_2, si_3, si_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * sh_0[k]
                 + si_0[k];

        t_1[k] = -ab_x[k] * sh_1[k]
                 + si_1[k];

        t_2[k] = -ab_x[k] * sh_2[k]
                 + si_2[k];

        t_3[k] = -ab_x[k] * sh_3[k]
                 + si_3[k];

        t_4[k] = -ab_x[k] * sh_4[k]
                 + si_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, sh_5, sh_6, sh_7, sh_8, sh_9, si_5, \
                         si_6, si_7, si_8, si_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * sh_5[k]
                 + si_5[k];

        t_6[k] = -ab_x[k] * sh_6[k]
                 + si_6[k];

        t_7[k] = -ab_x[k] * sh_7[k]
                 + si_7[k];

        t_8[k] = -ab_x[k] * sh_8[k]
                 + si_8[k];

        t_9[k] = -ab_x[k] * sh_9[k]
                 + si_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, sh_10, sh_11, sh_12, sh_13, \
                         sh_14, si_10, si_11, si_12, si_13, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * sh_10[k]
                  + si_10[k];

        t_11[k] = -ab_x[k] * sh_11[k]
                  + si_11[k];

        t_12[k] = -ab_x[k] * sh_12[k]
                  + si_12[k];

        t_13[k] = -ab_x[k] * sh_13[k]
                  + si_13[k];

        t_14[k] = -ab_x[k] * sh_14[k]
                  + si_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, sh_15, sh_16, sh_17, sh_18, \
                         sh_19, si_15, si_16, si_17, si_18, si_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * sh_15[k]
                  + si_15[k];

        t_16[k] = -ab_x[k] * sh_16[k]
                  + si_16[k];

        t_17[k] = -ab_x[k] * sh_17[k]
                  + si_17[k];

        t_18[k] = -ab_x[k] * sh_18[k]
                  + si_18[k];

        t_19[k] = -ab_x[k] * sh_19[k]
                  + si_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, ab_x, ab_y, sh_0, sh_1, sh_2, sh_20, si_1, \
                         si_3, si_4, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * sh_20[k]
                  + si_20[k];

        t_21[k] = -ab_y[k] * sh_0[k]
                  + si_1[k];

        t_22[k] = -ab_y[k] * sh_1[k]
                  + si_3[k];

        t_23[k] = -ab_y[k] * sh_2[k]
                  + si_4[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, ab_y, sh_3, sh_4, sh_5, sh_6, sh_7, \
                         si_6, si_7, si_8, si_10, si_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_24[k] = -ab_y[k] * sh_3[k]
                  + si_6[k];

        t_25[k] = -ab_y[k] * sh_4[k]
                  + si_7[k];

        t_26[k] = -ab_y[k] * sh_5[k]
                  + si_8[k];

        t_27[k] = -ab_y[k] * sh_6[k]
                  + si_10[k];

        t_28[k] = -ab_y[k] * sh_7[k]
                  + si_11[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_y, sh_8, sh_9, sh_10, sh_11, sh_12, \
                         si_12, si_13, si_15, si_16, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_29[k] = -ab_y[k] * sh_8[k]
                  + si_12[k];

        t_30[k] = -ab_y[k] * sh_9[k]
                  + si_13[k];

        t_31[k] = -ab_y[k] * sh_10[k]
                  + si_15[k];

        t_32[k] = -ab_y[k] * sh_11[k]
                  + si_16[k];

        t_33[k] = -ab_y[k] * sh_12[k]
                  + si_17[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, ab_y, sh_13, sh_14, sh_15, sh_16, \
                         sh_17, si_18, si_19, si_21, si_22, si_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_34[k] = -ab_y[k] * sh_13[k]
                  + si_18[k];

        t_35[k] = -ab_y[k] * sh_14[k]
                  + si_19[k];

        t_36[k] = -ab_y[k] * sh_15[k]
                  + si_21[k];

        t_37[k] = -ab_y[k] * sh_16[k]
                  + si_22[k];

        t_38[k] = -ab_y[k] * sh_17[k]
                  + si_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, ab_y, ab_z, sh_0, sh_18, sh_19, sh_20, si_2, \
                         si_24, si_25, si_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_39[k] = -ab_y[k] * sh_18[k]
                  + si_24[k];

        t_40[k] = -ab_y[k] * sh_19[k]
                  + si_25[k];

        t_41[k] = -ab_y[k] * sh_20[k]
                  + si_26[k];

        t_42[k] = -ab_z[k] * sh_0[k]
                  + si_2[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, ab_z, sh_1, sh_2, sh_3, sh_4, sh_5, \
                         si_4, si_5, si_7, si_8, si_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_43[k] = -ab_z[k] * sh_1[k]
                  + si_4[k];

        t_44[k] = -ab_z[k] * sh_2[k]
                  + si_5[k];

        t_45[k] = -ab_z[k] * sh_3[k]
                  + si_7[k];

        t_46[k] = -ab_z[k] * sh_4[k]
                  + si_8[k];

        t_47[k] = -ab_z[k] * sh_5[k]
                  + si_9[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, t_52, ab_z, sh_6, sh_7, sh_8, sh_9, sh_10, \
                         si_11, si_12, si_13, si_14, si_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_48[k] = -ab_z[k] * sh_6[k]
                  + si_11[k];

        t_49[k] = -ab_z[k] * sh_7[k]
                  + si_12[k];

        t_50[k] = -ab_z[k] * sh_8[k]
                  + si_13[k];

        t_51[k] = -ab_z[k] * sh_9[k]
                  + si_14[k];

        t_52[k] = -ab_z[k] * sh_10[k]
                  + si_16[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, ab_z, sh_11, sh_12, sh_13, sh_14, \
                         sh_15, si_17, si_18, si_19, si_20, si_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_53[k] = -ab_z[k] * sh_11[k]
                  + si_17[k];

        t_54[k] = -ab_z[k] * sh_12[k]
                  + si_18[k];

        t_55[k] = -ab_z[k] * sh_13[k]
                  + si_19[k];

        t_56[k] = -ab_z[k] * sh_14[k]
                  + si_20[k];

        t_57[k] = -ab_z[k] * sh_15[k]
                  + si_22[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, ab_z, sh_16, sh_17, sh_18, sh_19, \
                         sh_20, si_23, si_24, si_25, si_26, si_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_58[k] = -ab_z[k] * sh_16[k]
                  + si_23[k];

        t_59[k] = -ab_z[k] * sh_17[k]
                  + si_24[k];

        t_60[k] = -ab_z[k] * sh_18[k]
                  + si_25[k];

        t_61[k] = -ab_z[k] * sh_19[k]
                  + si_26[k];

        t_62[k] = -ab_z[k] * sh_20[k]
                  + si_27[k];
    }
}

}  // namespace simdovl
