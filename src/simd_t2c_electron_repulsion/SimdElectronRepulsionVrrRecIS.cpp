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


#include "SimdElectronRepulsionVrrRecIS.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_is_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_1 = buffer.data(gs0 + 1);
    const auto *gs0_2 = buffer.data(gs0 + 2);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_8 = buffer.data(gs0 + 8);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);
    const auto *gs1_11 = buffer.data(gs1 + 11);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, gs0_0, gs0_1, gs1_0, gs1_3, \
                         hs_0, hs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = pa_y[k] * hs_0[k];

        t_2[k] = pa_z[k] * hs_0[k];

        t_3[k] = f_2 * gs0_1[k]
                 - f_3 * gs1_3[k]
                 + pa_x[k] * hs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gs0_2, gs0_3, gs1_4, gs1_5, \
                         hs_3, hs_4, hs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gs0_2[k]
                 - f_3 * gs1_4[k]
                 + pa_x[k] * hs_4[k];

        t_5[k] = f_4 * gs0_3[k]
                 - f_5 * gs1_5[k]
                 + pa_x[k] * hs_5[k];

        t_6[k] = pa_z[k] * hs_3[k];

        t_7[k] = pa_y[k] * hs_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, gs0_4, gs0_5, gs0_6, gs1_6, gs1_7, \
                         gs1_9, hs_5, hs_6, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * gs0_4[k]
                 - f_5 * gs1_6[k]
                 + pa_x[k] * hs_6[k];

        t_9[k] = f_6 * gs0_5[k]
                 - f_7 * gs1_7[k]
                 + pa_x[k] * hs_7[k];

        t_10[k] = pa_z[k] * hs_5[k];

        t_11[k] = f_6 * gs0_6[k]
                  - f_7 * gs1_9[k]
                  + pa_x[k] * hs_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_x, pa_y, gs0_8, gs1_11, hs_6, \
                         hs_9, hs_10, hs_12, hs_13, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * hs_6[k];

        t_13[k] = f_6 * gs0_8[k]
                  - f_7 * gs1_11[k]
                  + pa_x[k] * hs_9[k];

        t_14[k] = pa_x[k] * hs_10[k];

        t_15[k] = pa_x[k] * hs_12[k];

        t_16[k] = pa_x[k] * hs_13[k];

        t_17[k] = pa_x[k] * hs_15[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pa_y, pa_z, gs0_5, gs0_6, gs0_7, gs1_7, \
                         gs1_9, gs1_10, hs_10, hs_12, hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * gs0_5[k]
                  - f_1 * gs1_7[k]
                  + pa_y[k] * hs_10[k];

        t_19[k] = pa_z[k] * hs_10[k];

        t_20[k] = f_2 * gs0_6[k]
                  - f_3 * gs1_9[k]
                  + pa_y[k] * hs_12[k];

        t_21[k] = f_4 * gs0_7[k]
                  - f_5 * gs1_10[k]
                  + pa_y[k] * hs_13[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pa_y, pa_z, gs0_8, gs1_11, hs_14, \
                         hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_6 * gs0_8[k]
                  - f_7 * gs1_11[k]
                  + pa_y[k] * hs_14[k];

        t_23[k] = pa_y[k] * hs_15[k];

        t_24[k] = f_0 * gs0_8[k]
                  - f_1 * gs1_11[k]
                  + pa_z[k] * hs_15[k];
    }
}

auto
compute_prim_is_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_1 = buffer.data(gs0 + 1);
    const auto *gs0_2 = buffer.data(gs0 + 2);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_8 = buffer.data(gs0 + 8);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_1 = buffer.data(gs1 + 1);
    const auto *gs1_2 = buffer.data(gs1 + 2);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_8 = buffer.data(gs1 + 8);

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

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, gs0_0, gs0_1, gs0_2, gs1_0, gs1_1, gs1_2, hs_0, \
                         hs_1, hs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = f_2 * gs0_1[k]
                 - f_3 * gs1_1[k]
                 + pa_x[k] * hs_1[k];

        t_2[k] = f_2 * gs0_2[k]
                 - f_3 * gs1_2[k]
                 + pa_x[k] * hs_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, gs0_3, gs0_4, gs0_5, gs1_3, gs1_4, gs1_5, hs_3, \
                         hs_4, hs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * gs0_3[k]
                 - f_5 * gs1_3[k]
                 + pa_x[k] * hs_3[k];

        t_4[k] = f_4 * gs0_4[k]
                 - f_5 * gs1_4[k]
                 + pa_x[k] * hs_4[k];

        t_5[k] = f_6 * gs0_5[k]
                 - f_7 * gs1_5[k]
                 + pa_x[k] * hs_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, gs0_5, gs0_6, gs0_8, gs1_5, gs1_6, gs1_8, \
                         hs_6, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gs0_8[k]
                 - f_7 * gs1_8[k]
                 + pa_x[k] * hs_6[k];

        t_7[k] = f_0 * gs0_5[k]
                 - f_1 * gs1_5[k]
                 + pa_y[k] * hs_7[k];

        t_8[k] = f_2 * gs0_6[k]
                 - f_3 * gs1_6[k]
                 + pa_y[k] * hs_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, gs0_7, gs0_8, gs1_7, gs1_8, hs_9, hs_10, \
                         hs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * gs0_7[k]
                 - f_5 * gs1_7[k]
                 + pa_y[k] * hs_9[k];

        t_10[k] = f_6 * gs0_8[k]
                  - f_7 * gs1_8[k]
                  + pa_y[k] * hs_10[k];

        t_11[k] = f_0 * gs0_8[k]
                  - f_1 * gs1_8[k]
                  + pa_z[k] * hs_11[k];
    }
}

auto
compute_prim_is_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_1 = buffer.data(gs0 + 1);
    const auto *gs0_2 = buffer.data(gs0 + 2);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_8 = buffer.data(gs0 + 8);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);
    const auto *gs1_11 = buffer.data(gs1 + 11);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);
    const auto *hs_16 = buffer.data(hs + 16);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, gs0_0, gs0_1, gs0_2, gs1_0, gs1_3, gs1_4, hs_0, \
                         hs_3, hs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = f_2 * gs0_1[k]
                 - f_3 * gs1_3[k]
                 + pa_x[k] * hs_3[k];

        t_2[k] = f_2 * gs0_2[k]
                 - f_3 * gs1_4[k]
                 + pa_x[k] * hs_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, gs0_3, gs0_4, gs0_5, gs1_5, gs1_6, gs1_7, hs_5, \
                         hs_8, hs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * gs0_3[k]
                 - f_5 * gs1_5[k]
                 + pa_x[k] * hs_5[k];

        t_4[k] = f_4 * gs0_4[k]
                 - f_5 * gs1_6[k]
                 + pa_x[k] * hs_8[k];

        t_5[k] = f_6 * gs0_5[k]
                 - f_7 * gs1_7[k]
                 + pa_x[k] * hs_9[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, gs0_5, gs0_6, gs0_8, gs1_7, gs1_9, gs1_11, \
                         hs_10, hs_11, hs_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gs0_8[k]
                 - f_7 * gs1_11[k]
                 + pa_x[k] * hs_10[k];

        t_7[k] = f_0 * gs0_5[k]
                 - f_1 * gs1_7[k]
                 + pa_y[k] * hs_11[k];

        t_8[k] = f_2 * gs0_6[k]
                 - f_3 * gs1_9[k]
                 + pa_y[k] * hs_13[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, gs0_7, gs0_8, gs1_10, gs1_11, hs_14, \
                         hs_15, hs_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * gs0_7[k]
                 - f_5 * gs1_10[k]
                 + pa_y[k] * hs_14[k];

        t_10[k] = f_6 * gs0_8[k]
                  - f_7 * gs1_11[k]
                  + pa_y[k] * hs_15[k];

        t_11[k] = f_0 * gs0_8[k]
                  - f_1 * gs1_11[k]
                  + pa_z[k] * hs_16[k];
    }
}

auto
compute_prim_is_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_9 = buffer.data(gs0 + 9);
    const auto *gs0_10 = buffer.data(gs0 + 10);
    const auto *gs0_11 = buffer.data(gs0 + 11);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);
    const auto *gs1_11 = buffer.data(gs1 + 11);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, gs0_0, gs0_3, gs1_0, gs1_3, \
                         hs_0, hs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = pa_y[k] * hs_0[k];

        t_2[k] = pa_z[k] * hs_0[k];

        t_3[k] = f_2 * gs0_3[k]
                 - f_3 * gs1_3[k]
                 + pa_x[k] * hs_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_y, pa_z, gs0_4, gs0_5, gs1_4, gs1_5, \
                         hs_3, hs_4, hs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gs0_4[k]
                 - f_3 * gs1_4[k]
                 + pa_x[k] * hs_4[k];

        t_5[k] = f_4 * gs0_5[k]
                 - f_5 * gs1_5[k]
                 + pa_x[k] * hs_5[k];

        t_6[k] = pa_z[k] * hs_3[k];

        t_7[k] = pa_y[k] * hs_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, gs0_6, gs0_7, gs0_9, gs1_6, gs1_7, \
                         gs1_9, hs_5, hs_6, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * gs0_6[k]
                 - f_5 * gs1_6[k]
                 + pa_x[k] * hs_6[k];

        t_9[k] = f_6 * gs0_7[k]
                 - f_7 * gs1_7[k]
                 + pa_x[k] * hs_7[k];

        t_10[k] = pa_z[k] * hs_5[k];

        t_11[k] = f_6 * gs0_9[k]
                  - f_7 * gs1_9[k]
                  + pa_x[k] * hs_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pa_x, pa_y, gs0_7, gs0_11, gs1_7, \
                         gs1_11, hs_6, hs_9, hs_10, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_y[k] * hs_6[k];

        t_13[k] = f_6 * gs0_11[k]
                  - f_7 * gs1_11[k]
                  + pa_x[k] * hs_9[k];

        t_14[k] = pa_x[k] * hs_10[k];

        t_15[k] = pa_x[k] * hs_15[k];

        t_16[k] = f_0 * gs0_7[k]
                  - f_1 * gs1_7[k]
                  + pa_y[k] * hs_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pa_y, pa_z, gs0_9, gs0_10, gs0_11, gs1_9, \
                         gs1_10, gs1_11, hs_10, hs_12, hs_13, hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = pa_z[k] * hs_10[k];

        t_18[k] = f_2 * gs0_9[k]
                  - f_3 * gs1_9[k]
                  + pa_y[k] * hs_12[k];

        t_19[k] = f_4 * gs0_10[k]
                  - f_5 * gs1_10[k]
                  + pa_y[k] * hs_13[k];

        t_20[k] = f_6 * gs0_11[k]
                  - f_7 * gs1_11[k]
                  + pa_y[k] * hs_14[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_y, pa_z, gs0_11, gs1_11, \
                         hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pa_y[k] * hs_15[k];

        t_22[k] = f_0 * gs0_11[k]
                  - f_1 * gs1_11[k]
                  + pa_z[k] * hs_15[k];
    }
}

auto
compute_prim_is_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_1 = buffer.data(gs0 + 1);
    const auto *gs0_2 = buffer.data(gs0 + 2);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_8 = buffer.data(gs0 + 8);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_7 = buffer.data(gs1 + 7);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);
    const auto *gs1_11 = buffer.data(gs1 + 11);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_10 = buffer.data(hs + 10);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);
    const auto *hs_15 = buffer.data(hs + 15);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, gs0_0, gs0_1, gs0_2, gs1_0, gs1_3, gs1_4, hs_0, \
                         hs_3, hs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = f_2 * gs0_1[k]
                 - f_3 * gs1_3[k]
                 + pa_x[k] * hs_3[k];

        t_2[k] = f_2 * gs0_2[k]
                 - f_3 * gs1_4[k]
                 + pa_x[k] * hs_4[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pa_x, gs0_3, gs0_4, gs0_5, gs1_5, gs1_6, gs1_7, hs_5, \
                         hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * gs0_3[k]
                 - f_5 * gs1_5[k]
                 + pa_x[k] * hs_5[k];

        t_4[k] = f_4 * gs0_4[k]
                 - f_5 * gs1_6[k]
                 + pa_x[k] * hs_7[k];

        t_5[k] = f_6 * gs0_5[k]
                 - f_7 * gs1_7[k]
                 + pa_x[k] * hs_8[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pa_x, pa_y, gs0_5, gs0_6, gs0_8, gs1_7, gs1_9, gs1_11, \
                         hs_9, hs_10, hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * gs0_8[k]
                 - f_7 * gs1_11[k]
                 + pa_x[k] * hs_9[k];

        t_7[k] = f_0 * gs0_5[k]
                 - f_1 * gs1_7[k]
                 + pa_y[k] * hs_10[k];

        t_8[k] = f_2 * gs0_6[k]
                 - f_3 * gs1_9[k]
                 + pa_y[k] * hs_12[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pa_y, pa_z, gs0_7, gs0_8, gs1_10, gs1_11, hs_13, \
                         hs_14, hs_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_4 * gs0_7[k]
                 - f_5 * gs1_10[k]
                 + pa_y[k] * hs_13[k];

        t_10[k] = f_6 * gs0_8[k]
                  - f_7 * gs1_11[k]
                  + pa_y[k] * hs_14[k];

        t_11[k] = f_0 * gs0_8[k]
                  - f_1 * gs1_11[k]
                  + pa_z[k] * hs_15[k];
    }
}

auto
compute_prim_is_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t gs0, const size_t gs1, const size_t hs,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / alpha;
    const auto f_1 = 2.5 * beta / (alpha * p);
    const auto f_2 = 1.5 / alpha;
    const auto f_3 = 1.5 * beta / (alpha * p);
    const auto f_4 = 1.0 / alpha;
    const auto f_5 = beta / (alpha * p);
    const auto f_6 = 0.5 / alpha;
    const auto f_7 = 0.5 * beta / (alpha * p);

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

    const auto *gs0_0 = buffer.data(gs0 + 0);
    const auto *gs0_3 = buffer.data(gs0 + 3);
    const auto *gs0_4 = buffer.data(gs0 + 4);
    const auto *gs0_5 = buffer.data(gs0 + 5);
    const auto *gs0_6 = buffer.data(gs0 + 6);
    const auto *gs0_7 = buffer.data(gs0 + 7);
    const auto *gs0_9 = buffer.data(gs0 + 9);
    const auto *gs0_10 = buffer.data(gs0 + 10);
    const auto *gs0_11 = buffer.data(gs0 + 11);

    const auto *gs1_0 = buffer.data(gs1 + 0);
    const auto *gs1_2 = buffer.data(gs1 + 2);
    const auto *gs1_3 = buffer.data(gs1 + 3);
    const auto *gs1_4 = buffer.data(gs1 + 4);
    const auto *gs1_5 = buffer.data(gs1 + 5);
    const auto *gs1_6 = buffer.data(gs1 + 6);
    const auto *gs1_8 = buffer.data(gs1 + 8);
    const auto *gs1_9 = buffer.data(gs1 + 9);
    const auto *gs1_10 = buffer.data(gs1 + 10);

    const auto *hs_0 = buffer.data(hs + 0);
    const auto *hs_2 = buffer.data(hs + 2);
    const auto *hs_3 = buffer.data(hs + 3);
    const auto *hs_4 = buffer.data(hs + 4);
    const auto *hs_5 = buffer.data(hs + 5);
    const auto *hs_6 = buffer.data(hs + 6);
    const auto *hs_7 = buffer.data(hs + 7);
    const auto *hs_8 = buffer.data(hs + 8);
    const auto *hs_9 = buffer.data(hs + 9);
    const auto *hs_11 = buffer.data(hs + 11);
    const auto *hs_12 = buffer.data(hs + 12);
    const auto *hs_13 = buffer.data(hs + 13);
    const auto *hs_14 = buffer.data(hs + 14);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pa_x, pa_y, pa_z, gs0_0, gs0_3, gs1_0, gs1_2, \
                         hs_0, hs_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gs0_0[k]
                 - f_1 * gs1_0[k]
                 + pa_x[k] * hs_0[k];

        t_1[k] = pa_y[k] * hs_0[k];

        t_2[k] = pa_z[k] * hs_0[k];

        t_3[k] = f_2 * gs0_3[k]
                 - f_3 * gs1_2[k]
                 + pa_x[k] * hs_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pa_x, pa_z, gs0_4, gs0_5, gs0_6, gs1_3, gs1_4, \
                         gs1_5, hs_2, hs_3, hs_4, hs_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gs0_4[k]
                 - f_3 * gs1_3[k]
                 + pa_x[k] * hs_3[k];

        t_5[k] = f_4 * gs0_5[k]
                 - f_5 * gs1_4[k]
                 + pa_x[k] * hs_4[k];

        t_6[k] = pa_z[k] * hs_2[k];

        t_7[k] = f_4 * gs0_6[k]
                 - f_5 * gs1_5[k]
                 + pa_x[k] * hs_5[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pa_x, pa_z, gs0_7, gs0_9, gs0_11, gs1_6, gs1_8, \
                         gs1_10, hs_4, hs_6, hs_7, hs_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * gs0_7[k]
                 - f_7 * gs1_6[k]
                 + pa_x[k] * hs_6[k];

        t_9[k] = pa_z[k] * hs_4[k];

        t_10[k] = f_6 * gs0_9[k]
                  - f_7 * gs1_8[k]
                  + pa_x[k] * hs_7[k];

        t_11[k] = f_6 * gs0_11[k]
                  - f_7 * gs1_10[k]
                  + pa_x[k] * hs_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pa_y, pa_z, gs0_7, gs0_9, gs0_10, gs1_6, \
                         gs1_8, gs1_9, hs_9, hs_11, hs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * gs0_7[k]
                  - f_1 * gs1_6[k]
                  + pa_y[k] * hs_9[k];

        t_13[k] = pa_z[k] * hs_9[k];

        t_14[k] = f_2 * gs0_9[k]
                  - f_3 * gs1_8[k]
                  + pa_y[k] * hs_11[k];

        t_15[k] = f_4 * gs0_10[k]
                  - f_5 * gs1_9[k]
                  + pa_y[k] * hs_12[k];
    }

#pragma omp simd aligned(t_16, t_17, pa_y, pa_z, gs0_11, gs1_10, hs_13, \
                         hs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * gs0_11[k]
                  - f_7 * gs1_10[k]
                  + pa_y[k] * hs_13[k];

        t_17[k] = f_0 * gs0_11[k]
                  - f_1 * gs1_10[k]
                  + pa_z[k] * hs_14[k];
    }
}

}  // namespace simdt2ceri
