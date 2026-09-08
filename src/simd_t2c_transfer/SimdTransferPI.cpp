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


#include "SimdTransferPI.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_hrr_pi_sph(double *values, const size_t nvalues, CSimdMatrix &buffer,
                   const CSimdMatrix &coordinates, const size_t si, const size_t sk,
                   const size_t nmax) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.1875 * std::sqrt(462.0);
    const auto f_1 = 0.625 * std::sqrt(462.0);
    const auto f_2 = 0.9375 * std::sqrt(154.0);
    const auto f_3 = 1.875 * std::sqrt(154.0);
    const auto f_4 = 0.1875 * std::sqrt(154.0);
    const auto f_5 = 0.75 * std::sqrt(7.0);
    const auto f_6 = 7.5 * std::sqrt(7.0);
    const auto f_7 = 0.5625 * std::sqrt(210.0);
    const auto f_8 = 0.375 * std::sqrt(210.0);
    const auto f_9 = 1.5 * std::sqrt(210.0);
    const auto f_10 = 0.1875 * std::sqrt(210.0);
    const auto f_11 = 0.5 * std::sqrt(210.0);
    const auto f_12 = 0.0625 * std::sqrt(210.0);
    const auto f_13 = 0.125 * std::sqrt(210.0);
    const auto f_14 = std::sqrt(210.0);
    const auto f_15 = 0.625 * std::sqrt(21.0);
    const auto f_16 = 1.25 * std::sqrt(21.0);
    const auto f_17 = 2.5 * std::sqrt(21.0);
    const auto f_18 = std::sqrt(21.0);
    const auto f_19 = 0.03125 * std::sqrt(210.0);
    const auto f_20 = 0.1875 * std::sqrt(7.0);
    const auto f_21 = 0.9375 * std::sqrt(7.0);
    const auto f_22 = 1.875 * std::sqrt(7.0);
    const auto f_23 = 11.25 * std::sqrt(7.0);
    const auto f_24 = 0.03125 * std::sqrt(462.0);
    const auto f_25 = 0.46875 * std::sqrt(462.0);

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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_13 = buffer.data(sk + 13);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

#pragma omp simd aligned(ab_y, si_1, si_4, si_6, si_11, si_15, si_22, sk_3, sk_7, sk_10, \
                         sk_16, sk_21, sk_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = -f_0 * ab_y[k] * si_1[k]
                 + f_1 * ab_y[k] * si_6[k]
                 - f_0 * ab_y[k] * si_15[k]
                 + f_0 * sk_3[k]
                 - f_1 * sk_10[k]
                 + f_0 * sk_21[k];

        g_1[k] = -f_2 * ab_y[k] * si_4[k]
                 + f_3 * ab_y[k] * si_11[k]
                 - f_4 * ab_y[k] * si_22[k]
                 + f_2 * sk_7[k]
                 - f_3 * sk_16[k]
                 + f_4 * sk_29[k];
    }

#pragma omp simd aligned(ab_y, si_1, si_8, si_15, si_17, sk_3, sk_12, sk_21, \
                         sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_2[k] = f_5 * ab_y[k] * si_1[k]
                 - f_6 * ab_y[k] * si_8[k]
                 - f_5 * ab_y[k] * si_15[k]
                 + f_6 * ab_y[k] * si_17[k]
                 - f_5 * sk_3[k]
                 + f_6 * sk_12[k]
                 + f_5 * sk_21[k]
                 - f_6 * sk_23[k];
    }

#pragma omp simd aligned(ab_y, si_4, si_11, si_13, si_22, si_24, sk_7, sk_16, sk_18, sk_29, \
                         sk_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_3[k] = f_7 * ab_y[k] * si_4[k]
                 + f_8 * ab_y[k] * si_11[k]
                 - f_9 * ab_y[k] * si_13[k]
                 - f_10 * ab_y[k] * si_22[k]
                 + f_11 * ab_y[k] * si_24[k]
                 - f_7 * sk_7[k]
                 - f_8 * sk_16[k]
                 + f_9 * sk_18[k]
                 + f_10 * sk_29[k]
                 - f_11 * sk_31[k];
    }

#pragma omp simd aligned(ab_y, si_1, si_6, si_8, si_15, si_17, si_19, sk_3, sk_10, sk_12, \
                         sk_21, sk_23, sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = -f_12 * ab_y[k] * si_1[k]
                 - f_13 * ab_y[k] * si_6[k]
                 + f_14 * ab_y[k] * si_8[k]
                 - f_12 * ab_y[k] * si_15[k]
                 + f_14 * ab_y[k] * si_17[k]
                 - f_14 * ab_y[k] * si_19[k]
                 + f_12 * sk_3[k]
                 + f_13 * sk_10[k]
                 - f_14 * sk_12[k]
                 + f_12 * sk_21[k]
                 - f_14 * sk_23[k]
                 + f_14 * sk_25[k];
    }

#pragma omp simd aligned(ab_y, si_4, si_11, si_13, si_22, si_24, si_26, sk_7, sk_16, sk_18, \
                         sk_29, sk_31, sk_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_5[k] = -f_15 * ab_y[k] * si_4[k]
                 - f_16 * ab_y[k] * si_11[k]
                 + f_17 * ab_y[k] * si_13[k]
                 - f_15 * ab_y[k] * si_22[k]
                 + f_17 * ab_y[k] * si_24[k]
                 - f_18 * ab_y[k] * si_26[k]
                 + f_15 * sk_7[k]
                 + f_16 * sk_16[k]
                 - f_17 * sk_18[k]
                 + f_15 * sk_29[k]
                 - f_17 * sk_31[k]
                 + f_18 * sk_33[k];
    }

#pragma omp simd aligned(ab_y, si_0, si_3, si_5, si_10, si_12, si_14, si_21, si_23, si_25, \
                         si_27, sk_1, sk_6, sk_8, sk_15, sk_17, sk_19, sk_28, sk_30, sk_32, \
                         sk_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_6[k] = 0.3125 * ab_y[k] * si_0[k]
                 + 0.9375 * ab_y[k] * si_3[k]
                 - 5.625 * ab_y[k] * si_5[k]
                 + 0.9375 * ab_y[k] * si_10[k]
                 - 11.25 * ab_y[k] * si_12[k]
                 + 7.5 * ab_y[k] * si_14[k]
                 + 0.3125 * ab_y[k] * si_21[k]
                 - 5.625 * ab_y[k] * si_23[k]
                 + 7.5 * ab_y[k] * si_25[k]
                 - ab_y[k] * si_27[k]
                 - 0.3125 * sk_1[k]
                 - 0.9375 * sk_6[k]
                 + 5.625 * sk_8[k]
                 - 0.9375 * sk_15[k]
                 + 11.25 * sk_17[k]
                 - 7.5 * sk_19[k]
                 - 0.3125 * sk_28[k]
                 + 5.625 * sk_30[k]
                 - 7.5 * sk_32[k]
                 + sk_34[k];
    }

#pragma omp simd aligned(ab_y, si_2, si_7, si_9, si_16, si_18, si_20, sk_4, sk_11, sk_13, \
                         sk_22, sk_24, sk_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_7[k] = -f_15 * ab_y[k] * si_2[k]
                 - f_16 * ab_y[k] * si_7[k]
                 + f_17 * ab_y[k] * si_9[k]
                 - f_15 * ab_y[k] * si_16[k]
                 + f_17 * ab_y[k] * si_18[k]
                 - f_18 * ab_y[k] * si_20[k]
                 + f_15 * sk_4[k]
                 + f_16 * sk_11[k]
                 - f_17 * sk_13[k]
                 + f_15 * sk_22[k]
                 - f_17 * sk_24[k]
                 + f_18 * sk_26[k];
    }

#pragma omp simd aligned(ab_y, si_0, si_3, si_5, si_10, si_14, si_21, si_23, si_25, sk_1, \
                         sk_6, sk_8, sk_15, sk_19, sk_28, sk_30, \
                         sk_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_8[k] = -f_19 * ab_y[k] * si_0[k]
                 - f_19 * ab_y[k] * si_3[k]
                 + f_11 * ab_y[k] * si_5[k]
                 + f_19 * ab_y[k] * si_10[k]
                 - f_11 * ab_y[k] * si_14[k]
                 + f_19 * ab_y[k] * si_21[k]
                 - f_11 * ab_y[k] * si_23[k]
                 + f_11 * ab_y[k] * si_25[k]
                 + f_19 * sk_1[k]
                 + f_19 * sk_6[k]
                 - f_11 * sk_8[k]
                 - f_19 * sk_15[k]
                 + f_11 * sk_19[k]
                 - f_19 * sk_28[k]
                 + f_11 * sk_30[k]
                 - f_11 * sk_32[k];
    }

#pragma omp simd aligned(ab_y, si_2, si_7, si_9, si_16, si_18, sk_4, sk_11, sk_13, sk_22, \
                         sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_9[k] = f_10 * ab_y[k] * si_2[k]
                 - f_8 * ab_y[k] * si_7[k]
                 - f_11 * ab_y[k] * si_9[k]
                 - f_7 * ab_y[k] * si_16[k]
                 + f_9 * ab_y[k] * si_18[k]
                 - f_10 * sk_4[k]
                 + f_8 * sk_11[k]
                 + f_11 * sk_13[k]
                 + f_7 * sk_22[k]
                 - f_9 * sk_24[k];
    }

#pragma omp simd aligned(ab_y, si_0, si_3, si_5, si_10, si_12, si_21, si_23, sk_1, sk_6, sk_8, \
                         sk_15, sk_17, sk_28, sk_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = f_20 * ab_y[k] * si_0[k]
                  - f_21 * ab_y[k] * si_3[k]
                  - f_22 * ab_y[k] * si_5[k]
                  - f_21 * ab_y[k] * si_10[k]
                  + f_23 * ab_y[k] * si_12[k]
                  + f_20 * ab_y[k] * si_21[k]
                  - f_22 * ab_y[k] * si_23[k]
                  - f_20 * sk_1[k]
                  + f_21 * sk_6[k]
                  + f_22 * sk_8[k]
                  + f_21 * sk_15[k]
                  - f_23 * sk_17[k]
                  - f_20 * sk_28[k]
                  + f_22 * sk_30[k];
    }

#pragma omp simd aligned(ab_y, si_2, si_7, si_16, sk_4, sk_11, sk_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_11[k] = -f_4 * ab_y[k] * si_2[k]
                  + f_3 * ab_y[k] * si_7[k]
                  - f_2 * ab_y[k] * si_16[k]
                  + f_4 * sk_4[k]
                  - f_3 * sk_11[k]
                  + f_2 * sk_22[k];
    }

#pragma omp simd aligned(ab_y, si_0, si_3, si_10, si_21, sk_1, sk_6, sk_15, \
                         sk_28 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = -f_24 * ab_y[k] * si_0[k]
                  + f_25 * ab_y[k] * si_3[k]
                  - f_25 * ab_y[k] * si_10[k]
                  + f_24 * ab_y[k] * si_21[k]
                  + f_24 * sk_1[k]
                  - f_25 * sk_6[k]
                  + f_25 * sk_15[k]
                  - f_24 * sk_28[k];
    }

#pragma omp simd aligned(ab_z, si_1, si_4, si_6, si_11, si_15, si_22, sk_4, sk_8, sk_11, \
                         sk_17, sk_22, sk_30 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_13[k] = -f_0 * ab_z[k] * si_1[k]
                  + f_1 * ab_z[k] * si_6[k]
                  - f_0 * ab_z[k] * si_15[k]
                  + f_0 * sk_4[k]
                  - f_1 * sk_11[k]
                  + f_0 * sk_22[k];

        g_14[k] = -f_2 * ab_z[k] * si_4[k]
                  + f_3 * ab_z[k] * si_11[k]
                  - f_4 * ab_z[k] * si_22[k]
                  + f_2 * sk_8[k]
                  - f_3 * sk_17[k]
                  + f_4 * sk_30[k];
    }

#pragma omp simd aligned(ab_z, si_1, si_8, si_15, si_17, sk_4, sk_13, sk_22, \
                         sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = f_5 * ab_z[k] * si_1[k]
                  - f_6 * ab_z[k] * si_8[k]
                  - f_5 * ab_z[k] * si_15[k]
                  + f_6 * ab_z[k] * si_17[k]
                  - f_5 * sk_4[k]
                  + f_6 * sk_13[k]
                  + f_5 * sk_22[k]
                  - f_6 * sk_24[k];
    }

#pragma omp simd aligned(ab_z, si_4, si_11, si_13, si_22, si_24, sk_8, sk_17, sk_19, sk_30, \
                         sk_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_16[k] = f_7 * ab_z[k] * si_4[k]
                  + f_8 * ab_z[k] * si_11[k]
                  - f_9 * ab_z[k] * si_13[k]
                  - f_10 * ab_z[k] * si_22[k]
                  + f_11 * ab_z[k] * si_24[k]
                  - f_7 * sk_8[k]
                  - f_8 * sk_17[k]
                  + f_9 * sk_19[k]
                  + f_10 * sk_30[k]
                  - f_11 * sk_32[k];
    }

#pragma omp simd aligned(ab_z, si_1, si_6, si_8, si_15, si_17, si_19, sk_4, sk_11, sk_13, \
                         sk_22, sk_24, sk_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = -f_12 * ab_z[k] * si_1[k]
                  - f_13 * ab_z[k] * si_6[k]
                  + f_14 * ab_z[k] * si_8[k]
                  - f_12 * ab_z[k] * si_15[k]
                  + f_14 * ab_z[k] * si_17[k]
                  - f_14 * ab_z[k] * si_19[k]
                  + f_12 * sk_4[k]
                  + f_13 * sk_11[k]
                  - f_14 * sk_13[k]
                  + f_12 * sk_22[k]
                  - f_14 * sk_24[k]
                  + f_14 * sk_26[k];
    }

#pragma omp simd aligned(ab_z, si_4, si_11, si_13, si_22, si_24, si_26, sk_8, sk_17, sk_19, \
                         sk_30, sk_32, sk_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_18[k] = -f_15 * ab_z[k] * si_4[k]
                  - f_16 * ab_z[k] * si_11[k]
                  + f_17 * ab_z[k] * si_13[k]
                  - f_15 * ab_z[k] * si_22[k]
                  + f_17 * ab_z[k] * si_24[k]
                  - f_18 * ab_z[k] * si_26[k]
                  + f_15 * sk_8[k]
                  + f_16 * sk_17[k]
                  - f_17 * sk_19[k]
                  + f_15 * sk_30[k]
                  - f_17 * sk_32[k]
                  + f_18 * sk_34[k];
    }

#pragma omp simd aligned(ab_z, si_0, si_3, si_5, si_10, si_12, si_14, si_21, si_23, si_25, \
                         si_27, sk_2, sk_7, sk_9, sk_16, sk_18, sk_20, sk_29, sk_31, sk_33, \
                         sk_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_19[k] = 0.3125 * ab_z[k] * si_0[k]
                  + 0.9375 * ab_z[k] * si_3[k]
                  - 5.625 * ab_z[k] * si_5[k]
                  + 0.9375 * ab_z[k] * si_10[k]
                  - 11.25 * ab_z[k] * si_12[k]
                  + 7.5 * ab_z[k] * si_14[k]
                  + 0.3125 * ab_z[k] * si_21[k]
                  - 5.625 * ab_z[k] * si_23[k]
                  + 7.5 * ab_z[k] * si_25[k]
                  - ab_z[k] * si_27[k]
                  - 0.3125 * sk_2[k]
                  - 0.9375 * sk_7[k]
                  + 5.625 * sk_9[k]
                  - 0.9375 * sk_16[k]
                  + 11.25 * sk_18[k]
                  - 7.5 * sk_20[k]
                  - 0.3125 * sk_29[k]
                  + 5.625 * sk_31[k]
                  - 7.5 * sk_33[k]
                  + sk_35[k];
    }

#pragma omp simd aligned(ab_z, si_2, si_7, si_9, si_16, si_18, si_20, sk_5, sk_12, sk_14, \
                         sk_23, sk_25, sk_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_15 * ab_z[k] * si_2[k]
                  - f_16 * ab_z[k] * si_7[k]
                  + f_17 * ab_z[k] * si_9[k]
                  - f_15 * ab_z[k] * si_16[k]
                  + f_17 * ab_z[k] * si_18[k]
                  - f_18 * ab_z[k] * si_20[k]
                  + f_15 * sk_5[k]
                  + f_16 * sk_12[k]
                  - f_17 * sk_14[k]
                  + f_15 * sk_23[k]
                  - f_17 * sk_25[k]
                  + f_18 * sk_27[k];
    }

#pragma omp simd aligned(ab_z, si_0, si_3, si_5, si_10, si_14, si_21, si_23, si_25, sk_2, \
                         sk_7, sk_9, sk_16, sk_20, sk_29, sk_31, \
                         sk_33 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_21[k] = -f_19 * ab_z[k] * si_0[k]
                  - f_19 * ab_z[k] * si_3[k]
                  + f_11 * ab_z[k] * si_5[k]
                  + f_19 * ab_z[k] * si_10[k]
                  - f_11 * ab_z[k] * si_14[k]
                  + f_19 * ab_z[k] * si_21[k]
                  - f_11 * ab_z[k] * si_23[k]
                  + f_11 * ab_z[k] * si_25[k]
                  + f_19 * sk_2[k]
                  + f_19 * sk_7[k]
                  - f_11 * sk_9[k]
                  - f_19 * sk_16[k]
                  + f_11 * sk_20[k]
                  - f_19 * sk_29[k]
                  + f_11 * sk_31[k]
                  - f_11 * sk_33[k];
    }

#pragma omp simd aligned(ab_z, si_2, si_7, si_9, si_16, si_18, sk_5, sk_12, sk_14, sk_23, \
                         sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_10 * ab_z[k] * si_2[k]
                  - f_8 * ab_z[k] * si_7[k]
                  - f_11 * ab_z[k] * si_9[k]
                  - f_7 * ab_z[k] * si_16[k]
                  + f_9 * ab_z[k] * si_18[k]
                  - f_10 * sk_5[k]
                  + f_8 * sk_12[k]
                  + f_11 * sk_14[k]
                  + f_7 * sk_23[k]
                  - f_9 * sk_25[k];
    }

#pragma omp simd aligned(ab_z, si_0, si_3, si_5, si_10, si_12, si_21, si_23, sk_2, sk_7, sk_9, \
                         sk_16, sk_18, sk_29, sk_31 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_23[k] = f_20 * ab_z[k] * si_0[k]
                  - f_21 * ab_z[k] * si_3[k]
                  - f_22 * ab_z[k] * si_5[k]
                  - f_21 * ab_z[k] * si_10[k]
                  + f_23 * ab_z[k] * si_12[k]
                  + f_20 * ab_z[k] * si_21[k]
                  - f_22 * ab_z[k] * si_23[k]
                  - f_20 * sk_2[k]
                  + f_21 * sk_7[k]
                  + f_22 * sk_9[k]
                  + f_21 * sk_16[k]
                  - f_23 * sk_18[k]
                  - f_20 * sk_29[k]
                  + f_22 * sk_31[k];
    }

#pragma omp simd aligned(ab_z, si_2, si_7, si_16, sk_5, sk_12, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_24[k] = -f_4 * ab_z[k] * si_2[k]
                  + f_3 * ab_z[k] * si_7[k]
                  - f_2 * ab_z[k] * si_16[k]
                  + f_4 * sk_5[k]
                  - f_3 * sk_12[k]
                  + f_2 * sk_23[k];
    }

#pragma omp simd aligned(ab_z, si_0, si_3, si_10, si_21, sk_2, sk_7, sk_16, \
                         sk_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = -f_24 * ab_z[k] * si_0[k]
                  + f_25 * ab_z[k] * si_3[k]
                  - f_25 * ab_z[k] * si_10[k]
                  + f_24 * ab_z[k] * si_21[k]
                  + f_24 * sk_2[k]
                  - f_25 * sk_7[k]
                  + f_25 * sk_16[k]
                  - f_24 * sk_29[k];
    }

#pragma omp simd aligned(ab_x, si_1, si_4, si_6, si_11, si_15, si_22, sk_1, sk_4, sk_6, sk_11, \
                         sk_15, sk_22 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_26[k] = -f_0 * ab_x[k] * si_1[k]
                  + f_1 * ab_x[k] * si_6[k]
                  - f_0 * ab_x[k] * si_15[k]
                  + f_0 * sk_1[k]
                  - f_1 * sk_6[k]
                  + f_0 * sk_15[k];

        g_27[k] = -f_2 * ab_x[k] * si_4[k]
                  + f_3 * ab_x[k] * si_11[k]
                  - f_4 * ab_x[k] * si_22[k]
                  + f_2 * sk_4[k]
                  - f_3 * sk_11[k]
                  + f_4 * sk_22[k];
    }

#pragma omp simd aligned(ab_x, si_1, si_8, si_15, si_17, sk_1, sk_8, sk_15, \
                         sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_28[k] = f_5 * ab_x[k] * si_1[k]
                  - f_6 * ab_x[k] * si_8[k]
                  - f_5 * ab_x[k] * si_15[k]
                  + f_6 * ab_x[k] * si_17[k]
                  - f_5 * sk_1[k]
                  + f_6 * sk_8[k]
                  + f_5 * sk_15[k]
                  - f_6 * sk_17[k];
    }

#pragma omp simd aligned(ab_x, si_4, si_11, si_13, si_22, si_24, sk_4, sk_11, sk_13, sk_22, \
                         sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_7 * ab_x[k] * si_4[k]
                  + f_8 * ab_x[k] * si_11[k]
                  - f_9 * ab_x[k] * si_13[k]
                  - f_10 * ab_x[k] * si_22[k]
                  + f_11 * ab_x[k] * si_24[k]
                  - f_7 * sk_4[k]
                  - f_8 * sk_11[k]
                  + f_9 * sk_13[k]
                  + f_10 * sk_22[k]
                  - f_11 * sk_24[k];
    }

#pragma omp simd aligned(ab_x, si_1, si_6, si_8, si_15, si_17, si_19, sk_1, sk_6, sk_8, sk_15, \
                         sk_17, sk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_30[k] = -f_12 * ab_x[k] * si_1[k]
                  - f_13 * ab_x[k] * si_6[k]
                  + f_14 * ab_x[k] * si_8[k]
                  - f_12 * ab_x[k] * si_15[k]
                  + f_14 * ab_x[k] * si_17[k]
                  - f_14 * ab_x[k] * si_19[k]
                  + f_12 * sk_1[k]
                  + f_13 * sk_6[k]
                  - f_14 * sk_8[k]
                  + f_12 * sk_15[k]
                  - f_14 * sk_17[k]
                  + f_14 * sk_19[k];
    }

#pragma omp simd aligned(ab_x, si_4, si_11, si_13, si_22, si_24, si_26, sk_4, sk_11, sk_13, \
                         sk_22, sk_24, sk_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_31[k] = -f_15 * ab_x[k] * si_4[k]
                  - f_16 * ab_x[k] * si_11[k]
                  + f_17 * ab_x[k] * si_13[k]
                  - f_15 * ab_x[k] * si_22[k]
                  + f_17 * ab_x[k] * si_24[k]
                  - f_18 * ab_x[k] * si_26[k]
                  + f_15 * sk_4[k]
                  + f_16 * sk_11[k]
                  - f_17 * sk_13[k]
                  + f_15 * sk_22[k]
                  - f_17 * sk_24[k]
                  + f_18 * sk_26[k];
    }

#pragma omp simd aligned(ab_x, si_0, si_3, si_5, si_10, si_12, si_14, si_21, si_23, si_25, \
                         si_27, sk_0, sk_3, sk_5, sk_10, sk_12, sk_14, sk_21, sk_23, sk_25, \
                         sk_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = 0.3125 * ab_x[k] * si_0[k]
                  + 0.9375 * ab_x[k] * si_3[k]
                  - 5.625 * ab_x[k] * si_5[k]
                  + 0.9375 * ab_x[k] * si_10[k]
                  - 11.25 * ab_x[k] * si_12[k]
                  + 7.5 * ab_x[k] * si_14[k]
                  + 0.3125 * ab_x[k] * si_21[k]
                  - 5.625 * ab_x[k] * si_23[k]
                  + 7.5 * ab_x[k] * si_25[k]
                  - ab_x[k] * si_27[k]
                  - 0.3125 * sk_0[k]
                  - 0.9375 * sk_3[k]
                  + 5.625 * sk_5[k]
                  - 0.9375 * sk_10[k]
                  + 11.25 * sk_12[k]
                  - 7.5 * sk_14[k]
                  - 0.3125 * sk_21[k]
                  + 5.625 * sk_23[k]
                  - 7.5 * sk_25[k]
                  + sk_27[k];
    }

#pragma omp simd aligned(ab_x, si_2, si_7, si_9, si_16, si_18, si_20, sk_2, sk_7, sk_9, sk_16, \
                         sk_18, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_33[k] = -f_15 * ab_x[k] * si_2[k]
                  - f_16 * ab_x[k] * si_7[k]
                  + f_17 * ab_x[k] * si_9[k]
                  - f_15 * ab_x[k] * si_16[k]
                  + f_17 * ab_x[k] * si_18[k]
                  - f_18 * ab_x[k] * si_20[k]
                  + f_15 * sk_2[k]
                  + f_16 * sk_7[k]
                  - f_17 * sk_9[k]
                  + f_15 * sk_16[k]
                  - f_17 * sk_18[k]
                  + f_18 * sk_20[k];
    }

#pragma omp simd aligned(ab_x, si_0, si_3, si_5, si_10, si_14, si_21, si_23, si_25, sk_0, \
                         sk_3, sk_5, sk_10, sk_14, sk_21, sk_23, \
                         sk_25 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_34[k] = -f_19 * ab_x[k] * si_0[k]
                  - f_19 * ab_x[k] * si_3[k]
                  + f_11 * ab_x[k] * si_5[k]
                  + f_19 * ab_x[k] * si_10[k]
                  - f_11 * ab_x[k] * si_14[k]
                  + f_19 * ab_x[k] * si_21[k]
                  - f_11 * ab_x[k] * si_23[k]
                  + f_11 * ab_x[k] * si_25[k]
                  + f_19 * sk_0[k]
                  + f_19 * sk_3[k]
                  - f_11 * sk_5[k]
                  - f_19 * sk_10[k]
                  + f_11 * sk_14[k]
                  - f_19 * sk_21[k]
                  + f_11 * sk_23[k]
                  - f_11 * sk_25[k];
    }

#pragma omp simd aligned(ab_x, si_2, si_7, si_9, si_16, si_18, sk_2, sk_7, sk_9, sk_16, \
                         sk_18 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_35[k] = f_10 * ab_x[k] * si_2[k]
                  - f_8 * ab_x[k] * si_7[k]
                  - f_11 * ab_x[k] * si_9[k]
                  - f_7 * ab_x[k] * si_16[k]
                  + f_9 * ab_x[k] * si_18[k]
                  - f_10 * sk_2[k]
                  + f_8 * sk_7[k]
                  + f_11 * sk_9[k]
                  + f_7 * sk_16[k]
                  - f_9 * sk_18[k];
    }

#pragma omp simd aligned(ab_x, si_0, si_3, si_5, si_10, si_12, si_21, si_23, sk_0, sk_3, sk_5, \
                         sk_10, sk_12, sk_21, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_36[k] = f_20 * ab_x[k] * si_0[k]
                  - f_21 * ab_x[k] * si_3[k]
                  - f_22 * ab_x[k] * si_5[k]
                  - f_21 * ab_x[k] * si_10[k]
                  + f_23 * ab_x[k] * si_12[k]
                  + f_20 * ab_x[k] * si_21[k]
                  - f_22 * ab_x[k] * si_23[k]
                  - f_20 * sk_0[k]
                  + f_21 * sk_3[k]
                  + f_22 * sk_5[k]
                  + f_21 * sk_10[k]
                  - f_23 * sk_12[k]
                  - f_20 * sk_21[k]
                  + f_22 * sk_23[k];
    }

#pragma omp simd aligned(ab_x, si_2, si_7, si_16, sk_2, sk_7, sk_16 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_37[k] = -f_4 * ab_x[k] * si_2[k]
                  + f_3 * ab_x[k] * si_7[k]
                  - f_2 * ab_x[k] * si_16[k]
                  + f_4 * sk_2[k]
                  - f_3 * sk_7[k]
                  + f_2 * sk_16[k];
    }

#pragma omp simd aligned(ab_x, si_0, si_3, si_10, si_21, sk_0, sk_3, sk_10, \
                         sk_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_38[k] = -f_24 * ab_x[k] * si_0[k]
                  + f_25 * ab_x[k] * si_3[k]
                  - f_25 * ab_x[k] * si_10[k]
                  + f_24 * ab_x[k] * si_21[k]
                  + f_24 * sk_0[k]
                  - f_25 * sk_3[k]
                  + f_25 * sk_10[k]
                  - f_24 * sk_21[k];
    }
}

auto
compute_hrr_pi(CSimdMatrix &buffer, const CSimdMatrix &coordinates, const size_t target,
               const size_t si, const size_t sk, const size_t nmax) -> void
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

    const auto *ab_x = coordinates.data(6);
    const auto *ab_y = coordinates.data(7);
    const auto *ab_z = coordinates.data(8);

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

    const auto *sk_0 = buffer.data(sk + 0);
    const auto *sk_1 = buffer.data(sk + 1);
    const auto *sk_2 = buffer.data(sk + 2);
    const auto *sk_3 = buffer.data(sk + 3);
    const auto *sk_4 = buffer.data(sk + 4);
    const auto *sk_5 = buffer.data(sk + 5);
    const auto *sk_6 = buffer.data(sk + 6);
    const auto *sk_7 = buffer.data(sk + 7);
    const auto *sk_8 = buffer.data(sk + 8);
    const auto *sk_9 = buffer.data(sk + 9);
    const auto *sk_10 = buffer.data(sk + 10);
    const auto *sk_11 = buffer.data(sk + 11);
    const auto *sk_12 = buffer.data(sk + 12);
    const auto *sk_13 = buffer.data(sk + 13);
    const auto *sk_14 = buffer.data(sk + 14);
    const auto *sk_15 = buffer.data(sk + 15);
    const auto *sk_16 = buffer.data(sk + 16);
    const auto *sk_17 = buffer.data(sk + 17);
    const auto *sk_18 = buffer.data(sk + 18);
    const auto *sk_19 = buffer.data(sk + 19);
    const auto *sk_20 = buffer.data(sk + 20);
    const auto *sk_21 = buffer.data(sk + 21);
    const auto *sk_22 = buffer.data(sk + 22);
    const auto *sk_23 = buffer.data(sk + 23);
    const auto *sk_24 = buffer.data(sk + 24);
    const auto *sk_25 = buffer.data(sk + 25);
    const auto *sk_26 = buffer.data(sk + 26);
    const auto *sk_27 = buffer.data(sk + 27);
    const auto *sk_28 = buffer.data(sk + 28);
    const auto *sk_29 = buffer.data(sk + 29);
    const auto *sk_30 = buffer.data(sk + 30);
    const auto *sk_31 = buffer.data(sk + 31);
    const auto *sk_32 = buffer.data(sk + 32);
    const auto *sk_33 = buffer.data(sk + 33);
    const auto *sk_34 = buffer.data(sk + 34);
    const auto *sk_35 = buffer.data(sk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, si_0, si_1, si_2, si_3, si_4, sk_0, \
                         sk_1, sk_2, sk_3, sk_4 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_0[k] = -ab_x[k] * si_0[k]
                 + sk_0[k];

        t_1[k] = -ab_x[k] * si_1[k]
                 + sk_1[k];

        t_2[k] = -ab_x[k] * si_2[k]
                 + sk_2[k];

        t_3[k] = -ab_x[k] * si_3[k]
                 + sk_3[k];

        t_4[k] = -ab_x[k] * si_4[k]
                 + sk_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, si_5, si_6, si_7, si_8, si_9, sk_5, \
                         sk_6, sk_7, sk_8, sk_9 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_5[k] = -ab_x[k] * si_5[k]
                 + sk_5[k];

        t_6[k] = -ab_x[k] * si_6[k]
                 + sk_6[k];

        t_7[k] = -ab_x[k] * si_7[k]
                 + sk_7[k];

        t_8[k] = -ab_x[k] * si_8[k]
                 + sk_8[k];

        t_9[k] = -ab_x[k] * si_9[k]
                 + sk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, si_10, si_11, si_12, si_13, \
                         si_14, sk_10, sk_11, sk_12, sk_13, sk_14 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_10[k] = -ab_x[k] * si_10[k]
                  + sk_10[k];

        t_11[k] = -ab_x[k] * si_11[k]
                  + sk_11[k];

        t_12[k] = -ab_x[k] * si_12[k]
                  + sk_12[k];

        t_13[k] = -ab_x[k] * si_13[k]
                  + sk_13[k];

        t_14[k] = -ab_x[k] * si_14[k]
                  + sk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, si_15, si_16, si_17, si_18, \
                         si_19, sk_15, sk_16, sk_17, sk_18, sk_19 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_15[k] = -ab_x[k] * si_15[k]
                  + sk_15[k];

        t_16[k] = -ab_x[k] * si_16[k]
                  + sk_16[k];

        t_17[k] = -ab_x[k] * si_17[k]
                  + sk_17[k];

        t_18[k] = -ab_x[k] * si_18[k]
                  + sk_18[k];

        t_19[k] = -ab_x[k] * si_19[k]
                  + sk_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, si_20, si_21, si_22, si_23, \
                         si_24, sk_20, sk_21, sk_22, sk_23, sk_24 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_20[k] = -ab_x[k] * si_20[k]
                  + sk_20[k];

        t_21[k] = -ab_x[k] * si_21[k]
                  + sk_21[k];

        t_22[k] = -ab_x[k] * si_22[k]
                  + sk_22[k];

        t_23[k] = -ab_x[k] * si_23[k]
                  + sk_23[k];

        t_24[k] = -ab_x[k] * si_24[k]
                  + sk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, ab_x, ab_y, si_0, si_25, si_26, si_27, sk_1, \
                         sk_25, sk_26, sk_27 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_25[k] = -ab_x[k] * si_25[k]
                  + sk_25[k];

        t_26[k] = -ab_x[k] * si_26[k]
                  + sk_26[k];

        t_27[k] = -ab_x[k] * si_27[k]
                  + sk_27[k];

        t_28[k] = -ab_y[k] * si_0[k]
                  + sk_1[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, ab_y, si_1, si_2, si_3, si_4, si_5, \
                         sk_3, sk_4, sk_6, sk_7, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_29[k] = -ab_y[k] * si_1[k]
                  + sk_3[k];

        t_30[k] = -ab_y[k] * si_2[k]
                  + sk_4[k];

        t_31[k] = -ab_y[k] * si_3[k]
                  + sk_6[k];

        t_32[k] = -ab_y[k] * si_4[k]
                  + sk_7[k];

        t_33[k] = -ab_y[k] * si_5[k]
                  + sk_8[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, ab_y, si_6, si_7, si_8, si_9, si_10, \
                         sk_10, sk_11, sk_12, sk_13, sk_15 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_34[k] = -ab_y[k] * si_6[k]
                  + sk_10[k];

        t_35[k] = -ab_y[k] * si_7[k]
                  + sk_11[k];

        t_36[k] = -ab_y[k] * si_8[k]
                  + sk_12[k];

        t_37[k] = -ab_y[k] * si_9[k]
                  + sk_13[k];

        t_38[k] = -ab_y[k] * si_10[k]
                  + sk_15[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, ab_y, si_11, si_12, si_13, si_14, \
                         si_15, sk_16, sk_17, sk_18, sk_19, sk_21 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_39[k] = -ab_y[k] * si_11[k]
                  + sk_16[k];

        t_40[k] = -ab_y[k] * si_12[k]
                  + sk_17[k];

        t_41[k] = -ab_y[k] * si_13[k]
                  + sk_18[k];

        t_42[k] = -ab_y[k] * si_14[k]
                  + sk_19[k];

        t_43[k] = -ab_y[k] * si_15[k]
                  + sk_21[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, ab_y, si_16, si_17, si_18, si_19, \
                         si_20, sk_22, sk_23, sk_24, sk_25, sk_26 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_44[k] = -ab_y[k] * si_16[k]
                  + sk_22[k];

        t_45[k] = -ab_y[k] * si_17[k]
                  + sk_23[k];

        t_46[k] = -ab_y[k] * si_18[k]
                  + sk_24[k];

        t_47[k] = -ab_y[k] * si_19[k]
                  + sk_25[k];

        t_48[k] = -ab_y[k] * si_20[k]
                  + sk_26[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, ab_y, si_21, si_22, si_23, si_24, \
                         si_25, sk_28, sk_29, sk_30, sk_31, sk_32 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_49[k] = -ab_y[k] * si_21[k]
                  + sk_28[k];

        t_50[k] = -ab_y[k] * si_22[k]
                  + sk_29[k];

        t_51[k] = -ab_y[k] * si_23[k]
                  + sk_30[k];

        t_52[k] = -ab_y[k] * si_24[k]
                  + sk_31[k];

        t_53[k] = -ab_y[k] * si_25[k]
                  + sk_32[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, ab_y, ab_z, si_0, si_1, si_26, si_27, sk_2, \
                         sk_4, sk_33, sk_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_54[k] = -ab_y[k] * si_26[k]
                  + sk_33[k];

        t_55[k] = -ab_y[k] * si_27[k]
                  + sk_34[k];

        t_56[k] = -ab_z[k] * si_0[k]
                  + sk_2[k];

        t_57[k] = -ab_z[k] * si_1[k]
                  + sk_4[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, t_62, ab_z, si_2, si_3, si_4, si_5, si_6, \
                         sk_5, sk_7, sk_8, sk_9, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_58[k] = -ab_z[k] * si_2[k]
                  + sk_5[k];

        t_59[k] = -ab_z[k] * si_3[k]
                  + sk_7[k];

        t_60[k] = -ab_z[k] * si_4[k]
                  + sk_8[k];

        t_61[k] = -ab_z[k] * si_5[k]
                  + sk_9[k];

        t_62[k] = -ab_z[k] * si_6[k]
                  + sk_11[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, ab_z, si_7, si_8, si_9, si_10, si_11, \
                         sk_12, sk_13, sk_14, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_63[k] = -ab_z[k] * si_7[k]
                  + sk_12[k];

        t_64[k] = -ab_z[k] * si_8[k]
                  + sk_13[k];

        t_65[k] = -ab_z[k] * si_9[k]
                  + sk_14[k];

        t_66[k] = -ab_z[k] * si_10[k]
                  + sk_16[k];

        t_67[k] = -ab_z[k] * si_11[k]
                  + sk_17[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, ab_z, si_12, si_13, si_14, si_15, \
                         si_16, sk_18, sk_19, sk_20, sk_22, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_68[k] = -ab_z[k] * si_12[k]
                  + sk_18[k];

        t_69[k] = -ab_z[k] * si_13[k]
                  + sk_19[k];

        t_70[k] = -ab_z[k] * si_14[k]
                  + sk_20[k];

        t_71[k] = -ab_z[k] * si_15[k]
                  + sk_22[k];

        t_72[k] = -ab_z[k] * si_16[k]
                  + sk_23[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, ab_z, si_17, si_18, si_19, si_20, \
                         si_21, sk_24, sk_25, sk_26, sk_27, sk_29 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_73[k] = -ab_z[k] * si_17[k]
                  + sk_24[k];

        t_74[k] = -ab_z[k] * si_18[k]
                  + sk_25[k];

        t_75[k] = -ab_z[k] * si_19[k]
                  + sk_26[k];

        t_76[k] = -ab_z[k] * si_20[k]
                  + sk_27[k];

        t_77[k] = -ab_z[k] * si_21[k]
                  + sk_29[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, ab_z, si_22, si_23, si_24, si_25, \
                         si_26, sk_30, sk_31, sk_32, sk_33, sk_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_78[k] = -ab_z[k] * si_22[k]
                  + sk_30[k];

        t_79[k] = -ab_z[k] * si_23[k]
                  + sk_31[k];

        t_80[k] = -ab_z[k] * si_24[k]
                  + sk_32[k];

        t_81[k] = -ab_z[k] * si_25[k]
                  + sk_33[k];

        t_82[k] = -ab_z[k] * si_26[k]
                  + sk_34[k];
    }

#pragma omp simd aligned(t_83, ab_z, si_27, sk_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        t_83[k] = -ab_z[k] * si_27[k]
                  + sk_35[k];
    }
}

}  // namespace simdovl
