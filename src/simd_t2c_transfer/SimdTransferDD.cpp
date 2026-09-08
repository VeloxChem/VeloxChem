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


#include "SimdTransferDD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_dd_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t pd, const size_t pf,
                   const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 0.25 * std::sqrt(3.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_1 = buffer.data(pd + 1);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_4 = buffer.data(pd + 4);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_13 = buffer.data(pd + 13);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_1 = buffer.data(pf + 1);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_4 = buffer.data(pf + 4);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_17 = buffer.data(pf + 17);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

#pragma omp simd aligned(ab_x, pd_6, pd_7, pd_9, pd_10, pd_11, pf_10, pf_11, pf_13, pf_14, \
                         pf_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -3.0 * ab_x[k] * pd_7[k]
                 + 3.0 * pf_11[k];

        g_1[k] = -3.0 * ab_x[k] * pd_10[k]
                 + 3.0 * pf_14[k];

        g_2[k] = f_0 * ab_x[k] * pd_6[k]
                 + f_0 * ab_x[k] * pd_9[k]
                 - f_1 * ab_x[k] * pd_11[k]
                 - f_0 * pf_10[k]
                 - f_0 * pf_13[k]
                 + f_1 * pf_15[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_6, pd_8, pd_9, pd_13, pd_16, pf_10, pf_12, pf_13, \
                         pf_23, pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -3.0 * ab_x[k] * pd_8[k]
                 + 3.0 * pf_12[k];

        g_4[k] = -1.5 * ab_x[k] * pd_6[k]
                 + 1.5 * ab_x[k] * pd_9[k]
                 + 1.5 * pf_10[k]
                 - 1.5 * pf_13[k];

        g_5[k] = -3.0 * ab_y[k] * pd_13[k]
                 + 3.0 * pf_23[k];

        g_6[k] = -3.0 * ab_y[k] * pd_16[k]
                 + 3.0 * pf_27[k];
    }

#pragma omp simd aligned(ab_y, pd_12, pd_14, pd_15, pd_17, pf_21, pf_24, pf_26, \
                         pf_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_0 * ab_y[k] * pd_12[k]
                 + f_0 * ab_y[k] * pd_15[k]
                 - f_1 * ab_y[k] * pd_17[k]
                 - f_0 * pf_21[k]
                 - f_0 * pf_26[k]
                 + f_1 * pf_28[k];

        g_8[k] = -3.0 * ab_y[k] * pd_14[k]
                 + 3.0 * pf_24[k];

        g_9[k] = -1.5 * ab_y[k] * pd_12[k]
                 + 1.5 * ab_y[k] * pd_15[k]
                 + 1.5 * pf_21[k]
                 - 1.5 * pf_26[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_1, pd_7, pd_13, pf_1, pf_13, \
                         pf_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_0 * ab_x[k] * pd_1[k]
                  + f_0 * ab_y[k] * pd_7[k]
                  - f_1 * ab_z[k] * pd_13[k]
                  - f_0 * pf_1[k]
                  - f_0 * pf_13[k]
                  + f_1 * pf_24[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_4, pd_10, pd_16, pf_4, pf_17, \
                         pf_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = f_0 * ab_x[k] * pd_4[k]
                  + f_0 * ab_y[k] * pd_10[k]
                  - f_1 * ab_z[k] * pd_16[k]
                  - f_0 * pf_4[k]
                  - f_0 * pf_17[k]
                  + f_1 * pf_28[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_0, pd_3, pd_5, pd_6, pd_9, pd_11, pd_12, pd_15, \
                         pd_17, pf_0, pf_3, pf_5, pf_11, pf_16, pf_18, pf_22, pf_27, \
                         pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -0.25 * ab_x[k] * pd_0[k]
                  - 0.25 * ab_x[k] * pd_3[k]
                  + 0.5 * ab_x[k] * pd_5[k]
                  - 0.25 * ab_y[k] * pd_6[k]
                  - 0.25 * ab_y[k] * pd_9[k]
                  + 0.5 * ab_y[k] * pd_11[k]
                  + 0.5 * ab_z[k] * pd_12[k]
                  + 0.5 * ab_z[k] * pd_15[k]
                  - ab_z[k] * pd_17[k]
                  + 0.25 * pf_0[k]
                  + 0.25 * pf_3[k]
                  - 0.5 * pf_5[k]
                  + 0.25 * pf_11[k]
                  + 0.25 * pf_16[k]
                  - 0.5 * pf_18[k]
                  - 0.5 * pf_22[k]
                  - 0.5 * pf_27[k]
                  + pf_29[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_2, pd_8, pd_14, pf_2, pf_14, \
                         pf_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_0 * ab_x[k] * pd_2[k]
                  + f_0 * ab_y[k] * pd_8[k]
                  - f_1 * ab_z[k] * pd_14[k]
                  - f_0 * pf_2[k]
                  - f_0 * pf_14[k]
                  + f_1 * pf_25[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_0, pd_3, pd_6, pd_9, pd_12, pd_15, pf_0, pf_3, \
                         pf_11, pf_16, pf_22, pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_2 * ab_x[k] * pd_0[k]
                  - f_2 * ab_x[k] * pd_3[k]
                  + f_2 * ab_y[k] * pd_6[k]
                  - f_2 * ab_y[k] * pd_9[k]
                  - f_0 * ab_z[k] * pd_12[k]
                  + f_0 * ab_z[k] * pd_15[k]
                  - f_2 * pf_0[k]
                  + f_2 * pf_3[k]
                  - f_2 * pf_11[k]
                  + f_2 * pf_16[k]
                  + f_0 * pf_22[k]
                  - f_0 * pf_27[k];
    }

#pragma omp simd aligned(ab_x, pd_12, pd_13, pd_15, pd_16, pd_17, pf_20, pf_21, pf_23, pf_24, \
                         pf_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -3.0 * ab_x[k] * pd_13[k]
                  + 3.0 * pf_21[k];

        g_16[k] = -3.0 * ab_x[k] * pd_16[k]
                  + 3.0 * pf_24[k];

        g_17[k] = f_0 * ab_x[k] * pd_12[k]
                  + f_0 * ab_x[k] * pd_15[k]
                  - f_1 * ab_x[k] * pd_17[k]
                  - f_0 * pf_20[k]
                  - f_0 * pf_23[k]
                  + f_1 * pf_25[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_1, pd_7, pd_12, pd_14, pd_15, pf_1, pf_13, pf_20, \
                         pf_22, pf_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -3.0 * ab_x[k] * pd_14[k]
                  + 3.0 * pf_22[k];

        g_19[k] = -1.5 * ab_x[k] * pd_12[k]
                  + 1.5 * ab_x[k] * pd_15[k]
                  + 1.5 * pf_20[k]
                  - 1.5 * pf_23[k];

        g_20[k] = -1.5 * ab_x[k] * pd_1[k]
                  + 1.5 * ab_y[k] * pd_7[k]
                  + 1.5 * pf_1[k]
                  - 1.5 * pf_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_4, pd_10, pf_4, pf_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -1.5 * ab_x[k] * pd_4[k]
                  + 1.5 * ab_y[k] * pd_10[k]
                  + 1.5 * pf_4[k]
                  - 1.5 * pf_17[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_0, pd_3, pd_5, pd_6, pd_9, pd_11, pf_0, pf_3, pf_5, \
                         pf_11, pf_16, pf_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_2 * ab_x[k] * pd_0[k]
                  + f_2 * ab_x[k] * pd_3[k]
                  - f_0 * ab_x[k] * pd_5[k]
                  - f_2 * ab_y[k] * pd_6[k]
                  - f_2 * ab_y[k] * pd_9[k]
                  + f_0 * ab_y[k] * pd_11[k]
                  - f_2 * pf_0[k]
                  - f_2 * pf_3[k]
                  + f_0 * pf_5[k]
                  + f_2 * pf_11[k]
                  + f_2 * pf_16[k]
                  - f_0 * pf_18[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_0, pd_2, pd_3, pd_6, pd_8, pd_9, pf_0, pf_2, pf_3, \
                         pf_11, pf_14, pf_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = -1.5 * ab_x[k] * pd_2[k]
                  + 1.5 * ab_y[k] * pd_8[k]
                  + 1.5 * pf_2[k]
                  - 1.5 * pf_14[k];

        g_24[k] = -0.75 * ab_x[k] * pd_0[k]
                  + 0.75 * ab_x[k] * pd_3[k]
                  + 0.75 * ab_y[k] * pd_6[k]
                  - 0.75 * ab_y[k] * pd_9[k]
                  + 0.75 * pf_0[k]
                  - 0.75 * pf_3[k]
                  - 0.75 * pf_11[k]
                  + 0.75 * pf_16[k];
    }
}

auto
compute_hrr_dd_sph_tri(double *values, const size_t nvalues, CSimdMatrix &buffer,
                       const CSimdMatrix &coordinates, const size_t pd, const size_t pf,
                       const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.5 * std::sqrt(3.0);
    const auto f_1 = std::sqrt(3.0);
    const auto f_2 = 0.25 * std::sqrt(3.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

    const auto *pd_0 = buffer.data(pd + 0);
    const auto *pd_2 = buffer.data(pd + 2);
    const auto *pd_3 = buffer.data(pd + 3);
    const auto *pd_5 = buffer.data(pd + 5);
    const auto *pd_6 = buffer.data(pd + 6);
    const auto *pd_7 = buffer.data(pd + 7);
    const auto *pd_8 = buffer.data(pd + 8);
    const auto *pd_9 = buffer.data(pd + 9);
    const auto *pd_10 = buffer.data(pd + 10);
    const auto *pd_11 = buffer.data(pd + 11);
    const auto *pd_12 = buffer.data(pd + 12);
    const auto *pd_14 = buffer.data(pd + 14);
    const auto *pd_15 = buffer.data(pd + 15);
    const auto *pd_16 = buffer.data(pd + 16);
    const auto *pd_17 = buffer.data(pd + 17);

    const auto *pf_0 = buffer.data(pf + 0);
    const auto *pf_2 = buffer.data(pf + 2);
    const auto *pf_3 = buffer.data(pf + 3);
    const auto *pf_5 = buffer.data(pf + 5);
    const auto *pf_10 = buffer.data(pf + 10);
    const auto *pf_11 = buffer.data(pf + 11);
    const auto *pf_12 = buffer.data(pf + 12);
    const auto *pf_13 = buffer.data(pf + 13);
    const auto *pf_14 = buffer.data(pf + 14);
    const auto *pf_15 = buffer.data(pf + 15);
    const auto *pf_16 = buffer.data(pf + 16);
    const auto *pf_18 = buffer.data(pf + 18);
    const auto *pf_20 = buffer.data(pf + 20);
    const auto *pf_21 = buffer.data(pf + 21);
    const auto *pf_22 = buffer.data(pf + 22);
    const auto *pf_23 = buffer.data(pf + 23);
    const auto *pf_24 = buffer.data(pf + 24);
    const auto *pf_25 = buffer.data(pf + 25);
    const auto *pf_26 = buffer.data(pf + 26);
    const auto *pf_27 = buffer.data(pf + 27);
    const auto *pf_28 = buffer.data(pf + 28);
    const auto *pf_29 = buffer.data(pf + 29);

#pragma omp simd aligned(ab_x, pd_6, pd_7, pd_9, pd_10, pd_11, pf_10, pf_11, pf_13, pf_14, \
                         pf_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -3.0 * ab_x[k] * pd_7[k]
                 + 3.0 * pf_11[k];

        g_1[k] = -3.0 * ab_x[k] * pd_10[k]
                 + 3.0 * pf_14[k];
        g_5[k] = g_1[k];

        g_2[k] = f_0 * ab_x[k] * pd_6[k]
                 + f_0 * ab_x[k] * pd_9[k]
                 - f_1 * ab_x[k] * pd_11[k]
                 - f_0 * pf_10[k]
                 - f_0 * pf_13[k]
                 + f_1 * pf_15[k];
        g_10[k] = g_2[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_6, pd_8, pd_9, pd_16, pf_10, pf_12, pf_13, \
                         pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = -3.0 * ab_x[k] * pd_8[k]
                 + 3.0 * pf_12[k];
        g_15[k] = g_3[k];

        g_4[k] = -1.5 * ab_x[k] * pd_6[k]
                 + 1.5 * ab_x[k] * pd_9[k]
                 + 1.5 * pf_10[k]
                 - 1.5 * pf_13[k];
        g_20[k] = g_4[k];

        g_6[k] = -3.0 * ab_y[k] * pd_16[k]
                 + 3.0 * pf_27[k];
    }

#pragma omp simd aligned(ab_y, pd_12, pd_14, pd_15, pd_17, pf_21, pf_24, pf_26, \
                         pf_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = f_0 * ab_y[k] * pd_12[k]
                 + f_0 * ab_y[k] * pd_15[k]
                 - f_1 * ab_y[k] * pd_17[k]
                 - f_0 * pf_21[k]
                 - f_0 * pf_26[k]
                 + f_1 * pf_28[k];
        g_11[k] = g_7[k];

        g_8[k] = -3.0 * ab_y[k] * pd_14[k]
                 + 3.0 * pf_24[k];
        g_16[k] = g_8[k];

        g_9[k] = -1.5 * ab_y[k] * pd_12[k]
                 + 1.5 * ab_y[k] * pd_15[k]
                 + 1.5 * pf_21[k]
                 - 1.5 * pf_26[k];
        g_21[k] = g_9[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_0, pd_3, pd_5, pd_6, pd_9, pd_11, pd_12, pd_15, \
                         pd_17, pf_0, pf_3, pf_5, pf_11, pf_16, pf_18, pf_22, pf_27, \
                         pf_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -0.25 * ab_x[k] * pd_0[k]
                  - 0.25 * ab_x[k] * pd_3[k]
                  + 0.5 * ab_x[k] * pd_5[k]
                  - 0.25 * ab_y[k] * pd_6[k]
                  - 0.25 * ab_y[k] * pd_9[k]
                  + 0.5 * ab_y[k] * pd_11[k]
                  + 0.5 * ab_z[k] * pd_12[k]
                  + 0.5 * ab_z[k] * pd_15[k]
                  - ab_z[k] * pd_17[k]
                  + 0.25 * pf_0[k]
                  + 0.25 * pf_3[k]
                  - 0.5 * pf_5[k]
                  + 0.25 * pf_11[k]
                  + 0.25 * pf_16[k]
                  - 0.5 * pf_18[k]
                  - 0.5 * pf_22[k]
                  - 0.5 * pf_27[k]
                  + pf_29[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_2, pd_8, pd_14, pf_2, pf_14, \
                         pf_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = f_0 * ab_x[k] * pd_2[k]
                  + f_0 * ab_y[k] * pd_8[k]
                  - f_1 * ab_z[k] * pd_14[k]
                  - f_0 * pf_2[k]
                  - f_0 * pf_14[k]
                  + f_1 * pf_25[k];
        g_17[k] = g_13[k];
    }

#pragma omp simd aligned(ab_x, ab_y, ab_z, pd_0, pd_3, pd_6, pd_9, pd_12, pd_15, pf_0, pf_3, \
                         pf_11, pf_16, pf_22, pf_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_14[k] = f_2 * ab_x[k] * pd_0[k]
                  - f_2 * ab_x[k] * pd_3[k]
                  + f_2 * ab_y[k] * pd_6[k]
                  - f_2 * ab_y[k] * pd_9[k]
                  - f_0 * ab_z[k] * pd_12[k]
                  + f_0 * ab_z[k] * pd_15[k]
                  - f_2 * pf_0[k]
                  + f_2 * pf_3[k]
                  - f_2 * pf_11[k]
                  + f_2 * pf_16[k]
                  + f_0 * pf_22[k]
                  - f_0 * pf_27[k];
        g_22[k] = g_14[k];
    }

#pragma omp simd aligned(ab_x, pd_12, pd_14, pd_15, pf_20, pf_22, \
                         pf_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -3.0 * ab_x[k] * pd_14[k]
                  + 3.0 * pf_22[k];

        g_19[k] = -1.5 * ab_x[k] * pd_12[k]
                  + 1.5 * ab_x[k] * pd_15[k]
                  + 1.5 * pf_20[k]
                  - 1.5 * pf_23[k];
        g_23[k] = g_19[k];
    }

#pragma omp simd aligned(ab_x, ab_y, pd_0, pd_3, pd_6, pd_9, pf_0, pf_3, pf_11, \
                         pf_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -0.75 * ab_x[k] * pd_0[k]
                  + 0.75 * ab_x[k] * pd_3[k]
                  + 0.75 * ab_y[k] * pd_6[k]
                  - 0.75 * ab_y[k] * pd_9[k]
                  + 0.75 * pf_0[k]
                  - 0.75 * pf_3[k]
                  - 0.75 * pf_11[k]
                  + 0.75 * pf_16[k];
    }
}

}  // namespace simdtrf
