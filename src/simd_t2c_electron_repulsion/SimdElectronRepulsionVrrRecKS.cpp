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


#include "SimdElectronRepulsionVrrRecKS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_ks_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_1, hs1_0, hs1_3, \
                         is_0, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, hs0_2, hs0_3, hs1_4, hs1_5, \
                         is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_6[k] = pa_z[k] * is_3[k];

        t_7[k] = pa_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_4, hs0_5, hs0_6, hs1_6, hs1_7, \
                         hs1_8, is_5, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_10[k] = pa_z[k] * is_5[k];

        t_11[k] = f_6 * hs0_6[k]
                  - f_7 * hs1_8[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, hs0_7, hs0_8, hs1_9, \
                         hs1_10, is_6, is_7, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * is_6[k];

        t_13[k] = f_6 * hs0_7[k]
                  - f_7 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_14[k] = f_8 * hs0_8[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = pa_z[k] * is_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, hs0_9, hs0_10, hs0_12, hs1_12, \
                         hs1_13, hs1_15, is_9, is_11, is_12, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_9[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_17[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_13[k]
                  + pa_x[k] * is_12[k];

        t_18[k] = pa_y[k] * is_9[k];

        t_19[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_15[k]
                  + pa_x[k] * is_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, hs0_8, hs1_10, is_14, \
                         is_16, is_17, is_18, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * is_14[k];

        t_21[k] = pa_x[k] * is_16[k];

        t_22[k] = pa_x[k] * is_17[k];

        t_23[k] = pa_x[k] * is_18[k];

        t_24[k] = pa_x[k] * is_20[k];

        t_25[k] = f_0 * hs0_8[k]
                  - f_1 * hs1_10[k]
                  + pa_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, hs0_9, hs0_10, hs0_11, hs1_12, \
                         hs1_13, hs1_14, is_14, is_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * is_14[k];

        t_27[k] = f_2 * hs0_9[k]
                  - f_3 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_28[k] = f_4 * hs0_10[k]
                  - f_5 * hs1_13[k]
                  + pa_y[k] * is_17[k];

        t_29[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_14[k]
                  + pa_y[k] * is_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, hs0_12, hs1_15, is_19, \
                         is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_15[k]
                  + pa_y[k] * is_19[k];

        t_31[k] = pa_y[k] * is_20[k];

        t_32[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_15[k]
                  + pa_z[k] * is_20[k];
    }
}

auto
compute_prim_ks_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_1, hs1_2, is_0, \
                         is_1, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_3, hs1_4, hs1_5, is_3, \
                         is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_y, hs0_6, hs0_7, hs0_11, hs1_6, hs1_7, \
                         hs1_11, is_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_11[k]
                 + pa_x[k] * is_8[k];

        t_9[k] = f_0 * hs0_7[k]
                 - f_1 * hs1_7[k]
                 + pa_y[k] * is_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, hs0_8, hs0_9, hs0_10, hs1_8, hs1_9, hs1_10, \
                         is_10, is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_8[k]
                  + pa_y[k] * is_10[k];

        t_11[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_9[k]
                  + pa_y[k] * is_11[k];

        t_12[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_10[k]
                  + pa_y[k] * is_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pa_z, hs0_11, hs1_11, is_13, \
                         is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_11[k]
                  + pa_y[k] * is_13[k];

        t_14[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_11[k]
                  + pa_z[k] * is_14[k];
    }
}

auto
compute_prim_ks_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);
    const auto *hs1_16 = buffer.data(hs1 + 16);
    const auto *hs1_17 = buffer.data(hs1 + 17);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_3, hs1_4, is_0, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_4[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_5, hs1_8, hs1_9, is_5, \
                         is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_8[k]
                 + pa_x[k] * is_8[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_9[k]
                 + pa_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_y, hs0_6, hs0_7, hs0_11, hs1_11, hs1_12, \
                         hs1_17, is_13, is_14, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_11[k]
                 + pa_x[k] * is_13[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_12[k]
                 + pa_x[k] * is_14[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_17[k]
                 + pa_x[k] * is_15[k];

        t_9[k] = f_0 * hs0_7[k]
                 - f_1 * hs1_12[k]
                 + pa_y[k] * is_16[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, hs0_8, hs0_9, hs0_10, hs1_14, hs1_15, hs1_16, \
                         is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_14[k]
                  + pa_y[k] * is_18[k];

        t_11[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_15[k]
                  + pa_y[k] * is_19[k];

        t_12[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_16[k]
                  + pa_y[k] * is_20[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pa_z, hs0_11, hs1_17, is_21, \
                         is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_17[k]
                  + pa_y[k] * is_21[k];

        t_14[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_17[k]
                  + pa_z[k] * is_22[k];
    }
}

auto
compute_prim_ks_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);
    const auto *hs0_16 = buffer.data(hs0 + 16);
    const auto *hs0_17 = buffer.data(hs0 + 17);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_3, hs1_0, hs1_3, \
                         is_0, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_3[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, hs0_4, hs0_5, hs1_4, hs1_5, \
                         is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_4[k]
                 - f_3 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_4 * hs0_5[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_6[k] = pa_z[k] * is_3[k];

        t_7[k] = pa_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_8, hs0_9, hs0_10, hs1_6, hs1_7, \
                         hs1_8, is_5, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * hs0_8[k]
                 - f_5 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_9[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_10[k] = pa_z[k] * is_5[k];

        t_11[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_8[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, hs0_11, hs0_12, hs1_9, \
                         hs1_10, is_6, is_7, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * is_6[k];

        t_13[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_14[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = pa_z[k] * is_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, hs0_14, hs0_15, hs0_17, hs1_12, \
                         hs1_13, hs1_15, is_9, is_11, is_12, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_14[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_17[k] = f_8 * hs0_15[k]
                  - f_9 * hs1_13[k]
                  + pa_x[k] * is_12[k];

        t_18[k] = pa_y[k] * is_9[k];

        t_19[k] = f_8 * hs0_17[k]
                  - f_9 * hs1_15[k]
                  + pa_x[k] * is_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pa_y, pa_z, hs0_12, hs0_14, \
                         hs1_10, hs1_12, is_14, is_16, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * is_14[k];

        t_21[k] = pa_x[k] * is_20[k];

        t_22[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_10[k]
                  + pa_y[k] * is_14[k];

        t_23[k] = pa_z[k] * is_14[k];

        t_24[k] = f_2 * hs0_14[k]
                  - f_3 * hs1_12[k]
                  + pa_y[k] * is_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, hs0_15, hs0_16, hs0_17, hs1_13, hs1_14, \
                         hs1_15, is_17, is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hs0_15[k]
                  - f_5 * hs1_13[k]
                  + pa_y[k] * is_17[k];

        t_26[k] = f_6 * hs0_16[k]
                  - f_7 * hs1_14[k]
                  + pa_y[k] * is_18[k];

        t_27[k] = f_8 * hs0_17[k]
                  - f_9 * hs1_15[k]
                  + pa_y[k] * is_19[k];

        t_28[k] = pa_y[k] * is_20[k];
    }

#pragma omp simd aligned(t_29, pa_z, hs0_17, hs1_15, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * hs0_17[k]
                  - f_1 * hs1_15[k]
                  + pa_z[k] * is_20[k];
    }
}

auto
compute_prim_ks_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);
    const auto *hs1_16 = buffer.data(hs1 + 16);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_3, hs1_4, is_0, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_4[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_5, hs1_7, hs1_8, is_5, \
                         is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_8[k]
                 + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, pa_y, hs0_6, hs0_7, hs0_11, hs1_10, hs1_11, \
                         hs1_16, is_11, is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_10[k]
                 + pa_x[k] * is_11[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_11[k]
                 + pa_x[k] * is_12[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_16[k]
                 + pa_x[k] * is_13[k];

        t_9[k] = f_0 * hs0_7[k]
                 - f_1 * hs1_11[k]
                 + pa_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, hs0_8, hs0_9, hs0_10, hs1_13, hs1_14, hs1_15, \
                         is_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_13[k]
                  + pa_y[k] * is_16[k];

        t_11[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_14[k]
                  + pa_y[k] * is_17[k];

        t_12[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_15[k]
                  + pa_y[k] * is_18[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_y, pa_z, hs0_11, hs1_16, is_19, \
                         is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_16[k]
                  + pa_y[k] * is_19[k];

        t_14[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_16[k]
                  + pa_z[k] * is_20[k];
    }
}

auto
compute_prim_ks_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_13 = buffer.data(hs0 + 13);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);
    const auto *hs0_16 = buffer.data(hs0 + 16);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_3, hs1_0, hs1_2, \
                         is_0, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_3[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, hs0_4, hs0_5, hs0_7, hs1_3, hs1_4, \
                         hs1_5, is_2, is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_4[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_5[k] = f_4 * hs0_5[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_6[k] = pa_z[k] * is_2[k];

        t_7[k] = f_4 * hs0_7[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_8, hs0_9, hs0_10, hs1_6, hs1_7, \
                         hs1_8, is_4, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * hs0_8[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = pa_z[k] * is_4[k];

        t_10[k] = f_6 * hs0_9[k]
                  - f_7 * hs1_7[k]
                  + pa_x[k] * is_7[k];

        t_11[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_8[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_z, hs0_11, hs0_13, hs0_14, hs1_9, \
                         hs1_11, hs1_12, is_6, is_9, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_13[k] = pa_z[k] * is_6[k];

        t_14[k] = f_8 * hs0_13[k]
                  - f_9 * hs1_11[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = f_8 * hs0_14[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, pa_z, hs0_11, hs0_13, hs0_16, \
                         hs1_9, hs1_11, hs1_14, is_12, is_13, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_16[k]
                  - f_9 * hs1_14[k]
                  + pa_x[k] * is_12[k];

        t_17[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_18[k] = pa_z[k] * is_13[k];

        t_19[k] = f_2 * hs0_13[k]
                  - f_3 * hs1_11[k]
                  + pa_y[k] * is_15[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pa_y, pa_z, hs0_14, hs0_15, hs0_16, hs1_12, \
                         hs1_13, hs1_14, is_16, is_17, is_18, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * hs0_14[k]
                  - f_5 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_21[k] = f_6 * hs0_15[k]
                  - f_7 * hs1_13[k]
                  + pa_y[k] * is_17[k];

        t_22[k] = f_8 * hs0_16[k]
                  - f_9 * hs1_14[k]
                  + pa_y[k] * is_18[k];

        t_23[k] = f_0 * hs0_16[k]
                  - f_1 * hs1_14[k]
                  + pa_z[k] * is_19[k];
    }
}

auto
compute_prim_ks_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_z, hs0_0, hs0_1, hs0_2, hs1_0, hs1_1, \
                         hs1_2, is_0, is_1, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_z[k] * is_0[k];

        t_2[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];

        t_3[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, hs0_3, hs0_4, hs0_5, hs1_3, hs1_4, hs1_5, is_3, \
                         is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_5[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_6[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, hs0_6, hs0_7, hs0_8, hs1_6, hs1_7, hs1_8, is_6, \
                         is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_8[k] = f_6 * hs0_7[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_9[k] = f_8 * hs0_8[k]
                 - f_9 * hs1_8[k]
                 + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, hs0_9, hs0_10, hs0_12, hs1_9, hs1_10, \
                         hs1_12, is_9, is_10, is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_8 * hs0_9[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_11[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];

        t_12[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_13[k] = pa_x[k] * is_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_x, pa_y, pa_z, hs0_8, hs1_8, \
                         is_12, is_13, is_14, is_15, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * is_13[k];

        t_15[k] = pa_x[k] * is_14[k];

        t_16[k] = pa_x[k] * is_15[k];

        t_17[k] = pa_x[k] * is_17[k];

        t_18[k] = f_0 * hs0_8[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_12[k];

        t_19[k] = pa_z[k] * is_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, hs0_9, hs0_10, hs0_11, hs1_9, hs1_10, hs1_11, \
                         is_13, is_14, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * hs0_9[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_21[k] = f_4 * hs0_10[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_14[k];

        t_22[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, hs0_12, hs1_12, is_16, \
                         is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_24[k] = pa_y[k] * is_17[k];

        t_25[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_17[k];
    }
}

auto
compute_prim_ks_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_1, hs1_0, hs1_1, \
                         is_0, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, hs0_2, hs0_3, hs1_2, hs1_3, \
                         is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_5[k];

        t_6[k] = pa_z[k] * is_3[k];

        t_7[k] = pa_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_4, hs0_5, hs0_6, hs1_4, hs1_5, \
                         hs1_6, is_5, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_7[k];

        t_10[k] = pa_z[k] * is_5[k];

        t_11[k] = f_6 * hs0_6[k]
                  - f_7 * hs1_6[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, hs0_7, hs0_8, hs1_7, hs1_8, \
                         is_6, is_7, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * is_6[k];

        t_13[k] = f_6 * hs0_7[k]
                  - f_7 * hs1_7[k]
                  + pa_x[k] * is_9[k];

        t_14[k] = f_8 * hs0_8[k]
                  - f_9 * hs1_8[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = pa_z[k] * is_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, hs0_9, hs0_10, hs0_12, hs1_9, \
                         hs1_10, hs1_12, is_9, is_11, is_12, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_9[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_11[k];

        t_17[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_12[k];

        t_18[k] = pa_y[k] * is_9[k];

        t_19[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pa_x, pa_y, hs0_8, hs1_8, is_14, \
                         is_16, is_17, is_18, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * is_14[k];

        t_21[k] = pa_x[k] * is_16[k];

        t_22[k] = pa_x[k] * is_17[k];

        t_23[k] = pa_x[k] * is_18[k];

        t_24[k] = pa_x[k] * is_20[k];

        t_25[k] = f_0 * hs0_8[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, hs0_9, hs0_10, hs0_11, hs1_9, \
                         hs1_10, hs1_11, is_14, is_16, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pa_z[k] * is_14[k];

        t_27[k] = f_2 * hs0_9[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_16[k];

        t_28[k] = f_4 * hs0_10[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_17[k];

        t_29[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pa_z, hs0_12, hs1_12, is_19, \
                         is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_19[k];

        t_31[k] = pa_y[k] * is_20[k];

        t_32[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_20[k];
    }
}

auto
compute_prim_ks_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_1, hs1_0, hs1_1, \
                         is_0, is_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, hs0_2, hs0_3, hs0_4, hs1_2, hs1_3, hs1_4, is_2, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];

        t_5[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_6[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, hs0_5, hs0_6, hs0_7, hs1_5, hs1_6, hs1_7, is_5, \
                         is_6, is_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_8[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_7[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, hs0_8, hs0_9, hs0_10, hs1_8, hs1_9, hs1_10, \
                         is_8, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_8 * hs0_8[k]
                  - f_9 * hs1_8[k]
                  + pa_x[k] * is_8[k];

        t_11[k] = f_8 * hs0_9[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_12[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_x, hs0_12, hs1_12, is_11, \
                         is_12, is_13, is_14, is_15, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_14[k] = pa_x[k] * is_12[k];

        t_15[k] = pa_x[k] * is_13[k];

        t_16[k] = pa_x[k] * is_14[k];

        t_17[k] = pa_x[k] * is_15[k];

        t_18[k] = pa_x[k] * is_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, hs0_8, hs0_9, hs0_10, hs1_8, \
                         hs1_9, hs1_10, is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * hs0_8[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_12[k];

        t_20[k] = pa_z[k] * is_12[k];

        t_21[k] = f_2 * hs0_9[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_22[k] = f_4 * hs0_10[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, hs0_11, hs0_12, hs1_11, hs1_12, \
                         is_15, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_15[k];

        t_24[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_25[k] = pa_y[k] * is_17[k];

        t_26[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_17[k];
    }
}

auto
compute_prim_ks_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t hs0, const size_t hs1, const size_t is,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);
    const auto *is_23 = buffer.data(is + 23);
    const auto *is_24 = buffer.data(is + 24);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_1, hs1_2, is_0, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_3[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_3, hs1_4, hs1_5, is_5, \
                         is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_5[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_8[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, hs0_6, hs0_7, hs0_11, hs1_7, hs1_8, hs1_12, \
                         is_13, is_14, is_17, is_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_13[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_8[k]
                 + pa_x[k] * is_14[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_12[k]
                 + pa_x[k] * is_17[k];

        t_9[k] = pa_x[k] * is_18[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, hs0_7, hs0_8, hs0_9, hs1_8, \
                         hs1_9, hs1_10, is_18, is_20, is_21, is_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * is_24[k];

        t_11[k] = f_0 * hs0_7[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_18[k];

        t_12[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_20[k];

        t_13[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_21[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, hs0_10, hs0_11, hs1_11, hs1_12, \
                         is_22, is_23, is_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_22[k];

        t_15[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_23[k];

        t_16[k] = pa_y[k] * is_24[k];

        t_17[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_24[k];
    }
}

auto
compute_prim_ks_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_1, hs1_0, hs1_3, \
                         is_0, is_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, hs0_2, hs0_3, hs1_4, hs1_5, \
                         is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_6[k] = pa_z[k] * is_3[k];

        t_7[k] = pa_y[k] * is_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_4, hs0_5, hs0_6, hs1_6, hs1_7, \
                         hs1_8, is_5, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_10[k] = pa_z[k] * is_5[k];

        t_11[k] = f_6 * hs0_6[k]
                  - f_7 * hs1_8[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_y, pa_z, hs0_7, hs0_8, hs1_9, \
                         hs1_10, is_6, is_7, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * is_6[k];

        t_13[k] = f_6 * hs0_7[k]
                  - f_7 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_14[k] = f_8 * hs0_8[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = pa_z[k] * is_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_x, pa_y, hs0_9, hs0_10, hs0_12, hs1_12, \
                         hs1_13, hs1_15, is_9, is_11, is_12, is_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_9[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_17[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_13[k]
                  + pa_x[k] * is_12[k];

        t_18[k] = pa_y[k] * is_9[k];

        t_19[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_15[k]
                  + pa_x[k] * is_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pa_x, pa_y, pa_z, hs0_8, hs0_9, hs1_10, \
                         hs1_12, is_14, is_16, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_x[k] * is_14[k];

        t_21[k] = pa_x[k] * is_20[k];

        t_22[k] = f_0 * hs0_8[k]
                  - f_1 * hs1_10[k]
                  + pa_y[k] * is_14[k];

        t_23[k] = pa_z[k] * is_14[k];

        t_24[k] = f_2 * hs0_9[k]
                  - f_3 * hs1_12[k]
                  + pa_y[k] * is_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pa_y, hs0_10, hs0_11, hs0_12, hs1_13, hs1_14, \
                         hs1_15, is_17, is_18, is_19, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hs0_10[k]
                  - f_5 * hs1_13[k]
                  + pa_y[k] * is_17[k];

        t_26[k] = f_6 * hs0_11[k]
                  - f_7 * hs1_14[k]
                  + pa_y[k] * is_18[k];

        t_27[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_15[k]
                  + pa_y[k] * is_19[k];

        t_28[k] = pa_y[k] * is_20[k];
    }

#pragma omp simd aligned(t_29, pa_z, hs0_12, hs1_15, is_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * hs0_12[k]
                  - f_1 * hs1_15[k]
                  + pa_z[k] * is_20[k];
    }
}

auto
compute_prim_ks_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_13 = buffer.data(hs0 + 13);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_3, hs1_0, hs1_1, \
                         is_0, is_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_3[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, hs0_4, hs0_5, hs0_6, hs1_2, hs1_3, hs1_4, is_2, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_4[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];

        t_5[k] = f_4 * hs0_5[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_6[k] = f_4 * hs0_6[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, hs0_7, hs0_8, hs0_9, hs1_5, hs1_6, hs1_7, is_5, \
                         is_6, is_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * hs0_7[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];

        t_8[k] = f_6 * hs0_8[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = f_6 * hs0_9[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_x, hs0_10, hs0_12, hs0_13, hs1_8, hs1_9, hs1_10, \
                         is_8, is_9, is_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_8[k]
                  + pa_x[k] * is_8[k];

        t_11[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_12[k] = f_8 * hs0_13[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pa_x, hs0_15, hs1_12, is_11, \
                         is_12, is_13, is_14, is_15, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * hs0_15[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_14[k] = pa_x[k] * is_12[k];

        t_15[k] = pa_x[k] * is_13[k];

        t_16[k] = pa_x[k] * is_14[k];

        t_17[k] = pa_x[k] * is_15[k];

        t_18[k] = pa_x[k] * is_17[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, hs0_10, hs0_12, hs0_13, hs1_8, \
                         hs1_9, hs1_10, is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * hs0_10[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_12[k];

        t_20[k] = pa_z[k] * is_12[k];

        t_21[k] = f_2 * hs0_12[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_22[k] = f_4 * hs0_13[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_14[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pa_z, hs0_14, hs0_15, hs1_11, hs1_12, \
                         is_15, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * hs0_14[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_15[k];

        t_24[k] = f_8 * hs0_15[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_25[k] = pa_y[k] * is_17[k];

        t_26[k] = f_0 * hs0_15[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_17[k];
    }
}

auto
compute_prim_ks_electron_repulsion_12(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);
    const auto *hs1_15 = buffer.data(hs1 + 15);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);
    const auto *is_20 = buffer.data(is + 20);
    const auto *is_21 = buffer.data(is + 21);
    const auto *is_22 = buffer.data(is + 22);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_2, hs1_3, is_0, \
                         is_3, is_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_3[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_4, hs1_6, hs1_7, is_5, \
                         is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_5[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_6[k]
                 + pa_x[k] * is_7[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, hs0_6, hs0_7, hs0_11, hs1_9, hs1_10, \
                         hs1_15, is_11, is_12, is_15, is_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_9[k]
                 + pa_x[k] * is_11[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_10[k]
                 + pa_x[k] * is_12[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_15[k]
                 + pa_x[k] * is_15[k];

        t_9[k] = pa_x[k] * is_16[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, hs0_7, hs0_8, hs0_9, hs1_10, \
                         hs1_12, hs1_13, is_16, is_18, is_19, is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * is_22[k];

        t_11[k] = f_0 * hs0_7[k]
                  - f_1 * hs1_10[k]
                  + pa_y[k] * is_16[k];

        t_12[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_12[k]
                  + pa_y[k] * is_18[k];

        t_13[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_13[k]
                  + pa_y[k] * is_19[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, hs0_10, hs0_11, hs1_14, hs1_15, \
                         is_20, is_21, is_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_14[k]
                  + pa_y[k] * is_20[k];

        t_15[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_15[k]
                  + pa_y[k] * is_21[k];

        t_16[k] = pa_y[k] * is_22[k];

        t_17[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_15[k]
                  + pa_z[k] * is_22[k];
    }
}

auto
compute_prim_ks_electron_repulsion_13(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_13 = buffer.data(hs0 + 13);
    const auto *hs0_14 = buffer.data(hs0 + 14);
    const auto *hs0_15 = buffer.data(hs0 + 15);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);
    const auto *hs1_13 = buffer.data(hs1 + 13);
    const auto *hs1_14 = buffer.data(hs1 + 14);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);
    const auto *is_18 = buffer.data(is + 18);
    const auto *is_19 = buffer.data(is + 19);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, hs0_0, hs0_2, hs1_0, hs1_2, \
                         is_0, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_y[k] * is_0[k];

        t_2[k] = pa_z[k] * is_0[k];

        t_3[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, hs0_3, hs0_4, hs0_6, hs1_3, hs1_4, \
                         hs1_5, is_2, is_3, is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hs0_3[k]
                 - f_3 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_5[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_6[k] = pa_z[k] * is_2[k];

        t_7[k] = f_4 * hs0_6[k]
                 - f_5 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, hs0_7, hs0_8, hs0_9, hs1_6, hs1_7, \
                         hs1_8, is_4, is_6, is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * hs0_7[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_9[k] = pa_z[k] * is_4[k];

        t_10[k] = f_6 * hs0_8[k]
                  - f_7 * hs1_7[k]
                  + pa_x[k] * is_7[k];

        t_11[k] = f_6 * hs0_9[k]
                  - f_7 * hs1_8[k]
                  + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_x, pa_z, hs0_10, hs0_12, hs0_13, hs1_9, \
                         hs1_11, hs1_12, is_6, is_9, is_10, is_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_8 * hs0_10[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_13[k] = pa_z[k] * is_6[k];

        t_14[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_11[k]
                  + pa_x[k] * is_10[k];

        t_15[k] = f_8 * hs0_13[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_x, pa_y, pa_z, hs0_10, hs0_15, \
                         hs1_9, hs1_14, is_12, is_13, is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * hs0_15[k]
                  - f_9 * hs1_14[k]
                  + pa_x[k] * is_12[k];

        t_17[k] = pa_x[k] * is_13[k];

        t_18[k] = pa_x[k] * is_19[k];

        t_19[k] = f_0 * hs0_10[k]
                  - f_1 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_20[k] = pa_z[k] * is_13[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pa_y, hs0_12, hs0_13, hs0_14, hs1_11, hs1_12, \
                         hs1_13, is_15, is_16, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * hs0_12[k]
                  - f_3 * hs1_11[k]
                  + pa_y[k] * is_15[k];

        t_22[k] = f_4 * hs0_13[k]
                  - f_5 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_23[k] = f_6 * hs0_14[k]
                  - f_7 * hs1_13[k]
                  + pa_y[k] * is_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_y, pa_z, hs0_15, hs1_14, is_18, \
                         is_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_8 * hs0_15[k]
                  - f_9 * hs1_14[k]
                  + pa_y[k] * is_18[k];

        t_25[k] = pa_y[k] * is_19[k];

        t_26[k] = f_0 * hs0_15[k]
                  - f_1 * hs1_14[k]
                  + pa_z[k] * is_19[k];
    }
}

auto
compute_prim_ks_electron_repulsion_14(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_11 = buffer.data(hs0 + 11);
    const auto *hs0_12 = buffer.data(hs0 + 12);
    const auto *hs0_13 = buffer.data(hs0 + 13);
    const auto *hs0_14 = buffer.data(hs0 + 14);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);
    const auto *hs1_12 = buffer.data(hs1 + 12);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);
    const auto *is_15 = buffer.data(is + 15);
    const auto *is_16 = buffer.data(is + 16);
    const auto *is_17 = buffer.data(is + 17);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_z, hs0_0, hs0_2, hs0_3, hs1_0, hs1_1, \
                         hs1_2, is_0, is_1, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = pa_z[k] * is_0[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];

        t_3[k] = f_2 * hs0_3[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pa_x, hs0_4, hs0_5, hs0_6, hs1_3, hs1_4, hs1_5, is_3, \
                         is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_5[k] = f_4 * hs0_5[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pa_x, hs0_7, hs0_8, hs0_9, hs1_6, hs1_7, hs1_8, is_6, \
                         is_7, is_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * hs0_7[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_8[k] = f_6 * hs0_8[k]
                 - f_7 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_9[k] = f_8 * hs0_9[k]
                 - f_9 * hs1_8[k]
                 + pa_x[k] * is_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, hs0_11, hs0_12, hs0_14, hs1_9, hs1_10, \
                         hs1_12, is_9, is_10, is_11, is_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_9[k]
                  + pa_x[k] * is_9[k];

        t_11[k] = f_8 * hs0_12[k]
                  - f_9 * hs1_10[k]
                  + pa_x[k] * is_10[k];

        t_12[k] = f_8 * hs0_14[k]
                  - f_9 * hs1_12[k]
                  + pa_x[k] * is_11[k];

        t_13[k] = pa_x[k] * is_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pa_x, pa_y, pa_z, hs0_9, hs1_8, \
                         is_12, is_13, is_14, is_15, is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pa_x[k] * is_13[k];

        t_15[k] = pa_x[k] * is_14[k];

        t_16[k] = pa_x[k] * is_15[k];

        t_17[k] = pa_x[k] * is_17[k];

        t_18[k] = f_0 * hs0_9[k]
                  - f_1 * hs1_8[k]
                  + pa_y[k] * is_12[k];

        t_19[k] = pa_z[k] * is_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, hs0_11, hs0_12, hs0_13, hs1_9, hs1_10, \
                         hs1_11, is_13, is_14, is_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * hs0_11[k]
                  - f_3 * hs1_9[k]
                  + pa_y[k] * is_13[k];

        t_21[k] = f_4 * hs0_12[k]
                  - f_5 * hs1_10[k]
                  + pa_y[k] * is_14[k];

        t_22[k] = f_6 * hs0_13[k]
                  - f_7 * hs1_11[k]
                  + pa_y[k] * is_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pa_z, hs0_14, hs1_12, is_16, \
                         is_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_8 * hs0_14[k]
                  - f_9 * hs1_12[k]
                  + pa_y[k] * is_16[k];

        t_24[k] = pa_y[k] * is_17[k];

        t_25[k] = f_0 * hs0_14[k]
                  - f_1 * hs1_12[k]
                  + pa_z[k] * is_17[k];
    }
}

auto
compute_prim_ks_electron_repulsion_15(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                      const size_t hs0, const size_t hs1, const size_t is,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / alpha;
    const auto f_1 = 3.0 * beta / (alpha * p);
    const auto f_2 = 2.0 / alpha;
    const auto f_3 = 2.0 * beta / (alpha * p);
    const auto f_4 = 1.5 / alpha;
    const auto f_5 = 1.5 * beta / (alpha * p);
    const auto f_6 = 1.0 / alpha;
    const auto f_7 = beta / (alpha * p);
    const auto f_8 = 0.5 / alpha;
    const auto f_9 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *hs0_0 = buffer.data(hs0 + 0);
    const auto *hs0_1 = buffer.data(hs0 + 1);
    const auto *hs0_2 = buffer.data(hs0 + 2);
    const auto *hs0_3 = buffer.data(hs0 + 3);
    const auto *hs0_4 = buffer.data(hs0 + 4);
    const auto *hs0_5 = buffer.data(hs0 + 5);
    const auto *hs0_6 = buffer.data(hs0 + 6);
    const auto *hs0_7 = buffer.data(hs0 + 7);
    const auto *hs0_8 = buffer.data(hs0 + 8);
    const auto *hs0_9 = buffer.data(hs0 + 9);
    const auto *hs0_10 = buffer.data(hs0 + 10);
    const auto *hs0_11 = buffer.data(hs0 + 11);

    const auto *hs1_0 = buffer.data(hs1 + 0);
    const auto *hs1_1 = buffer.data(hs1 + 1);
    const auto *hs1_2 = buffer.data(hs1 + 2);
    const auto *hs1_3 = buffer.data(hs1 + 3);
    const auto *hs1_4 = buffer.data(hs1 + 4);
    const auto *hs1_5 = buffer.data(hs1 + 5);
    const auto *hs1_6 = buffer.data(hs1 + 6);
    const auto *hs1_7 = buffer.data(hs1 + 7);
    const auto *hs1_8 = buffer.data(hs1 + 8);
    const auto *hs1_9 = buffer.data(hs1 + 9);
    const auto *hs1_10 = buffer.data(hs1 + 10);
    const auto *hs1_11 = buffer.data(hs1 + 11);

    const auto *is_0 = buffer.data(is + 0);
    const auto *is_1 = buffer.data(is + 1);
    const auto *is_2 = buffer.data(is + 2);
    const auto *is_3 = buffer.data(is + 3);
    const auto *is_4 = buffer.data(is + 4);
    const auto *is_5 = buffer.data(is + 5);
    const auto *is_6 = buffer.data(is + 6);
    const auto *is_7 = buffer.data(is + 7);
    const auto *is_8 = buffer.data(is + 8);
    const auto *is_9 = buffer.data(is + 9);
    const auto *is_10 = buffer.data(is + 10);
    const auto *is_11 = buffer.data(is + 11);
    const auto *is_12 = buffer.data(is + 12);
    const auto *is_13 = buffer.data(is + 13);
    const auto *is_14 = buffer.data(is + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, hs0_0, hs0_1, hs0_2, hs1_0, hs1_1, hs1_2, is_0, \
                         is_1, is_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hs0_0[k]
                 - f_1 * hs1_0[k]
                 + pa_x[k] * is_0[k];

        t_1[k] = f_2 * hs0_1[k]
                 - f_3 * hs1_1[k]
                 + pa_x[k] * is_1[k];

        t_2[k] = f_2 * hs0_2[k]
                 - f_3 * hs1_2[k]
                 + pa_x[k] * is_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, hs0_3, hs0_4, hs0_5, hs1_3, hs1_4, hs1_5, is_3, \
                         is_4, is_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * hs0_3[k]
                 - f_5 * hs1_3[k]
                 + pa_x[k] * is_3[k];

        t_4[k] = f_4 * hs0_4[k]
                 - f_5 * hs1_4[k]
                 + pa_x[k] * is_4[k];

        t_5[k] = f_6 * hs0_5[k]
                 - f_7 * hs1_5[k]
                 + pa_x[k] * is_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_x, hs0_6, hs0_7, hs0_11, hs1_6, hs1_7, hs1_11, \
                         is_6, is_7, is_8, is_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hs0_6[k]
                 - f_7 * hs1_6[k]
                 + pa_x[k] * is_6[k];

        t_7[k] = f_8 * hs0_7[k]
                 - f_9 * hs1_7[k]
                 + pa_x[k] * is_7[k];

        t_8[k] = f_8 * hs0_11[k]
                 - f_9 * hs1_11[k]
                 + pa_x[k] * is_8[k];

        t_9[k] = pa_x[k] * is_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_x, pa_y, hs0_7, hs0_8, hs0_9, hs1_7, \
                         hs1_8, hs1_9, is_9, is_10, is_11, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_x[k] * is_14[k];

        t_11[k] = f_0 * hs0_7[k]
                  - f_1 * hs1_7[k]
                  + pa_y[k] * is_9[k];

        t_12[k] = f_2 * hs0_8[k]
                  - f_3 * hs1_8[k]
                  + pa_y[k] * is_10[k];

        t_13[k] = f_4 * hs0_9[k]
                  - f_5 * hs1_9[k]
                  + pa_y[k] * is_11[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_y, pa_z, hs0_10, hs0_11, hs1_10, hs1_11, \
                         is_12, is_13, is_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * hs0_10[k]
                  - f_7 * hs1_10[k]
                  + pa_y[k] * is_12[k];

        t_15[k] = f_8 * hs0_11[k]
                  - f_9 * hs1_11[k]
                  + pa_y[k] * is_13[k];

        t_16[k] = pa_y[k] * is_14[k];

        t_17[k] = f_0 * hs0_11[k]
                  - f_1 * hs1_11[k]
                  + pa_z[k] * is_14[k];
    }
}

}  // namespace simdt2ceri
