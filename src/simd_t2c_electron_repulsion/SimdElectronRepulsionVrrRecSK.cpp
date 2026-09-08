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


#include "SimdElectronRepulsionVrrRecSK.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_sk_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);
    const auto *sh1_12 = buffer.data(sh1 + 12);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_z, sh0_0, sh0_1, sh0_2, sh1_0, sh1_1, \
                         sh1_2, si_0, si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = pb_z[k] * si_0[k];

        t_2[k] = f_2 * sh0_1[k]
                 - f_3 * sh1_1[k]
                 + pb_x[k] * si_1[k];

        t_3[k] = f_2 * sh0_2[k]
                 - f_3 * sh1_2[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, sh0_3, sh0_4, sh0_5, sh1_3, sh1_4, sh1_5, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * sh0_3[k]
                 - f_5 * sh1_3[k]
                 + pb_x[k] * si_3[k];

        t_5[k] = f_4 * sh0_4[k]
                 - f_5 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_6[k] = f_6 * sh0_5[k]
                 - f_7 * sh1_5[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, sh0_6, sh0_7, sh0_8, sh1_6, sh1_7, sh1_8, si_6, \
                         si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_6 * sh0_6[k]
                 - f_7 * sh1_6[k]
                 + pb_x[k] * si_6[k];

        t_8[k] = f_6 * sh0_7[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_9[k] = f_8 * sh0_8[k]
                 - f_9 * sh1_8[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_x, sh0_9, sh0_10, sh0_12, sh1_9, sh1_10, \
                         sh1_12, si_9, si_10, si_11, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_8 * sh0_9[k]
                  - f_9 * sh1_9[k]
                  + pb_x[k] * si_9[k];

        t_11[k] = f_8 * sh0_10[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_10[k];

        t_12[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_11[k];

        t_13[k] = pb_x[k] * si_12[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, t_18, t_19, pb_x, pb_y, pb_z, sh0_8, sh1_8, \
                         si_12, si_13, si_14, si_15, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = pb_x[k] * si_13[k];

        t_15[k] = pb_x[k] * si_14[k];

        t_16[k] = pb_x[k] * si_15[k];

        t_17[k] = pb_x[k] * si_17[k];

        t_18[k] = f_0 * sh0_8[k]
                  - f_1 * sh1_8[k]
                  + pb_y[k] * si_12[k];

        t_19[k] = pb_z[k] * si_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_y, sh0_9, sh0_10, sh0_11, sh1_9, sh1_10, sh1_11, \
                         si_13, si_14, si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * sh0_9[k]
                  - f_3 * sh1_9[k]
                  + pb_y[k] * si_13[k];

        t_21[k] = f_4 * sh0_10[k]
                  - f_5 * sh1_10[k]
                  + pb_y[k] * si_14[k];

        t_22[k] = f_6 * sh0_11[k]
                  - f_7 * sh1_11[k]
                  + pb_y[k] * si_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, pb_z, sh0_12, sh1_12, si_16, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_y[k] * si_16[k];

        t_24[k] = pb_y[k] * si_17[k];

        t_25[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_12[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_12 = buffer.data(sh1 + 12);
    const auto *sh1_13 = buffer.data(sh1 + 13);
    const auto *sh1_14 = buffer.data(sh1 + 14);
    const auto *sh1_15 = buffer.data(sh1 + 15);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_1, sh0_2, sh1_0, sh1_3, sh1_4, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_1[k]
                 - f_3 * sh1_3[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_2[k]
                 - f_3 * sh1_4[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_3, sh0_4, sh0_5, sh1_5, sh1_6, sh1_7, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_3[k]
                 - f_5 * sh1_5[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_4[k]
                 - f_5 * sh1_6[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_5[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, sh0_6, sh0_7, sh0_8, sh1_8, sh1_9, sh1_10, si_6, \
                         si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_6[k]
                 - f_7 * sh1_8[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_6 * sh0_7[k]
                 - f_7 * sh1_9[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_8[k]
                 - f_9 * sh1_10[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, sh0_9, sh0_10, sh0_12, sh1_12, sh1_13, sh1_15, \
                         si_9, si_10, si_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * sh0_9[k]
                 - f_9 * sh1_12[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_8 * sh0_10[k]
                  - f_9 * sh1_13[k]
                  + pb_x[k] * si_10[k];

        t_11[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_15[k]
                  + pb_x[k] * si_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, sh0_8, sh0_9, sh0_10, sh1_10, sh1_12, sh1_13, \
                         si_12, si_13, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sh0_8[k]
                  - f_1 * sh1_10[k]
                  + pb_y[k] * si_12[k];

        t_13[k] = f_2 * sh0_9[k]
                  - f_3 * sh1_12[k]
                  + pb_y[k] * si_13[k];

        t_14[k] = f_4 * sh0_10[k]
                  - f_5 * sh1_13[k]
                  + pb_y[k] * si_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, sh0_11, sh0_12, sh1_14, sh1_15, si_15, \
                         si_16, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * sh0_11[k]
                  - f_7 * sh1_14[k]
                  + pb_y[k] * si_15[k];

        t_16[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_15[k]
                  + pb_y[k] * si_16[k];

        t_17[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_15[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);
    const auto *sh1_12 = buffer.data(sh1 + 12);

    const auto *si_0 = buffer.data(si + 0);
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
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sh0_0, sh0_1, sh1_0, sh1_1, \
                         si_0, si_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = pb_y[k] * si_0[k];

        t_2[k] = pb_z[k] * si_0[k];

        t_3[k] = f_2 * sh0_1[k]
                 - f_3 * sh1_1[k]
                 + pb_x[k] * si_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sh0_2, sh0_3, sh1_2, sh1_3, \
                         si_3, si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sh0_2[k]
                 - f_3 * sh1_2[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_4 * sh0_3[k]
                 - f_5 * sh1_3[k]
                 + pb_x[k] * si_5[k];

        t_6[k] = pb_z[k] * si_3[k];

        t_7[k] = pb_y[k] * si_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sh0_4, sh0_5, sh0_6, sh1_4, sh1_5, \
                         sh1_6, si_5, si_6, si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sh0_4[k]
                 - f_5 * sh1_4[k]
                 + pb_x[k] * si_6[k];

        t_9[k] = f_6 * sh0_5[k]
                 - f_7 * sh1_5[k]
                 + pb_x[k] * si_7[k];

        t_10[k] = pb_z[k] * si_5[k];

        t_11[k] = f_6 * sh0_6[k]
                  - f_7 * sh1_6[k]
                  + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sh0_7, sh0_8, sh1_7, sh1_8, \
                         si_6, si_7, si_9, si_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_y[k] * si_6[k];

        t_13[k] = f_6 * sh0_7[k]
                  - f_7 * sh1_7[k]
                  + pb_x[k] * si_9[k];

        t_14[k] = f_8 * sh0_8[k]
                  - f_9 * sh1_8[k]
                  + pb_x[k] * si_10[k];

        t_15[k] = pb_z[k] * si_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, sh0_9, sh0_10, sh0_12, sh1_9, \
                         sh1_10, sh1_12, si_9, si_11, si_12, si_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * sh0_9[k]
                  - f_9 * sh1_9[k]
                  + pb_x[k] * si_11[k];

        t_17[k] = f_8 * sh0_10[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_12[k];

        t_18[k] = pb_y[k] * si_9[k];

        t_19[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, t_25, pb_x, pb_y, sh0_8, sh1_8, si_14, \
                         si_16, si_17, si_18, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_x[k] * si_14[k];

        t_21[k] = pb_x[k] * si_16[k];

        t_22[k] = pb_x[k] * si_17[k];

        t_23[k] = pb_x[k] * si_18[k];

        t_24[k] = pb_x[k] * si_20[k];

        t_25[k] = f_0 * sh0_8[k]
                  - f_1 * sh1_8[k]
                  + pb_y[k] * si_14[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_y, pb_z, sh0_9, sh0_10, sh0_11, sh1_9, \
                         sh1_10, sh1_11, si_14, si_16, si_17, si_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * si_14[k];

        t_27[k] = f_2 * sh0_9[k]
                  - f_3 * sh1_9[k]
                  + pb_y[k] * si_16[k];

        t_28[k] = f_4 * sh0_10[k]
                  - f_5 * sh1_10[k]
                  + pb_y[k] * si_17[k];

        t_29[k] = f_6 * sh0_11[k]
                  - f_7 * sh1_11[k]
                  + pb_y[k] * si_18[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pb_y, pb_z, sh0_12, sh1_12, si_19, \
                         si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_y[k] * si_19[k];

        t_31[k] = pb_y[k] * si_20[k];

        t_32[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_12[k]
                  + pb_z[k] * si_20[k];
    }
}

auto
compute_prim_sk_electron_repulsion_3(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);
    const auto *sh1_12 = buffer.data(sh1 + 12);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_1, sh0_2, sh1_0, sh1_1, sh1_2, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_1[k]
                 - f_3 * sh1_1[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_2[k]
                 - f_3 * sh1_2[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_3, sh0_4, sh0_5, sh1_3, sh1_4, sh1_5, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_3[k]
                 - f_5 * sh1_3[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_4[k]
                 - f_5 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_5[k]
                 - f_7 * sh1_5[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, sh0_6, sh0_7, sh0_8, sh1_6, sh1_7, sh1_8, si_6, \
                         si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_6[k]
                 - f_7 * sh1_6[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_6 * sh0_7[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_8[k]
                 - f_9 * sh1_8[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, sh0_9, sh0_10, sh0_12, sh1_9, sh1_10, \
                         sh1_12, si_9, si_10, si_11, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * sh0_9[k]
                 - f_9 * sh1_9[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_8 * sh0_10[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_10[k];

        t_11[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_11[k];

        t_12[k] = pb_x[k] * si_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, sh0_8, sh1_8, si_12, si_13, \
                         si_14, si_15, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * si_13[k];

        t_14[k] = pb_x[k] * si_14[k];

        t_15[k] = pb_x[k] * si_15[k];

        t_16[k] = pb_x[k] * si_17[k];

        t_17[k] = f_0 * sh0_8[k]
                  - f_1 * sh1_8[k]
                  + pb_y[k] * si_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, sh0_9, sh0_10, sh0_11, sh1_9, sh1_10, sh1_11, \
                         si_13, si_14, si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sh0_9[k]
                  - f_3 * sh1_9[k]
                  + pb_y[k] * si_13[k];

        t_19[k] = f_4 * sh0_10[k]
                  - f_5 * sh1_10[k]
                  + pb_y[k] * si_14[k];

        t_20[k] = f_6 * sh0_11[k]
                  - f_7 * sh1_11[k]
                  + pb_y[k] * si_15[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, sh0_12, sh1_12, si_16, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_12[k]
                  + pb_y[k] * si_16[k];

        t_22[k] = pb_y[k] * si_17[k];

        t_23[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_12[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_4(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);
    const auto *sh0_14 = buffer.data(sh0 + 14);
    const auto *sh0_15 = buffer.data(sh0 + 15);
    const auto *sh0_16 = buffer.data(sh0 + 16);
    const auto *sh0_17 = buffer.data(sh0 + 17);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_12 = buffer.data(sh1 + 12);
    const auto *sh1_13 = buffer.data(sh1 + 13);
    const auto *sh1_14 = buffer.data(sh1 + 14);
    const auto *sh1_15 = buffer.data(sh1 + 15);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_3, sh0_4, sh1_0, sh1_3, sh1_4, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_3[k]
                 - f_3 * sh1_3[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_4[k]
                 - f_3 * sh1_4[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_5, sh0_8, sh0_9, sh1_5, sh1_6, sh1_7, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_5[k]
                 - f_5 * sh1_5[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_8[k]
                 - f_5 * sh1_6[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_9[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, sh0_10, sh0_11, sh0_12, sh1_8, sh1_9, sh1_10, \
                         si_6, si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_10[k]
                 - f_7 * sh1_8[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_6 * sh0_11[k]
                 - f_7 * sh1_9[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_12[k]
                 - f_9 * sh1_10[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, sh0_14, sh0_15, sh0_17, sh1_12, sh1_13, \
                         sh1_15, si_9, si_10, si_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * sh0_14[k]
                 - f_9 * sh1_12[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_8 * sh0_15[k]
                  - f_9 * sh1_13[k]
                  + pb_x[k] * si_10[k];

        t_11[k] = f_8 * sh0_17[k]
                  - f_9 * sh1_15[k]
                  + pb_x[k] * si_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, sh0_12, sh0_14, sh0_15, sh1_10, sh1_12, \
                         sh1_13, si_12, si_13, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_10[k]
                  + pb_y[k] * si_12[k];

        t_13[k] = f_2 * sh0_14[k]
                  - f_3 * sh1_12[k]
                  + pb_y[k] * si_13[k];

        t_14[k] = f_4 * sh0_15[k]
                  - f_5 * sh1_13[k]
                  + pb_y[k] * si_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, sh0_16, sh0_17, sh1_14, sh1_15, si_15, \
                         si_16, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * sh0_16[k]
                  - f_7 * sh1_14[k]
                  + pb_y[k] * si_15[k];

        t_16[k] = f_8 * sh0_17[k]
                  - f_9 * sh1_15[k]
                  + pb_y[k] * si_16[k];

        t_17[k] = f_0 * sh0_17[k]
                  - f_1 * sh1_15[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_5(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_12 = buffer.data(sh1 + 12);
    const auto *sh1_13 = buffer.data(sh1 + 13);
    const auto *sh1_14 = buffer.data(sh1 + 14);
    const auto *sh1_15 = buffer.data(sh1 + 15);

    const auto *si_0 = buffer.data(si + 0);
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
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sh0_0, sh0_1, sh1_0, sh1_3, \
                         si_0, si_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = pb_y[k] * si_0[k];

        t_2[k] = pb_z[k] * si_0[k];

        t_3[k] = f_2 * sh0_1[k]
                 - f_3 * sh1_3[k]
                 + pb_x[k] * si_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_y, pb_z, sh0_2, sh0_3, sh1_4, sh1_5, \
                         si_3, si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sh0_2[k]
                 - f_3 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_4 * sh0_3[k]
                 - f_5 * sh1_5[k]
                 + pb_x[k] * si_5[k];

        t_6[k] = pb_z[k] * si_3[k];

        t_7[k] = pb_y[k] * si_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sh0_4, sh0_5, sh0_6, sh1_6, sh1_7, \
                         sh1_8, si_5, si_6, si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * sh0_4[k]
                 - f_5 * sh1_6[k]
                 + pb_x[k] * si_6[k];

        t_9[k] = f_6 * sh0_5[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_10[k] = pb_z[k] * si_5[k];

        t_11[k] = f_6 * sh0_6[k]
                  - f_7 * sh1_8[k]
                  + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, sh0_7, sh0_8, sh1_9, \
                         sh1_10, si_6, si_7, si_9, si_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_y[k] * si_6[k];

        t_13[k] = f_6 * sh0_7[k]
                  - f_7 * sh1_9[k]
                  + pb_x[k] * si_9[k];

        t_14[k] = f_8 * sh0_8[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_10[k];

        t_15[k] = pb_z[k] * si_7[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, sh0_9, sh0_10, sh0_12, sh1_12, \
                         sh1_13, sh1_15, si_9, si_11, si_12, si_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * sh0_9[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_11[k];

        t_17[k] = f_8 * sh0_10[k]
                  - f_9 * sh1_13[k]
                  + pb_x[k] * si_12[k];

        t_18[k] = pb_y[k] * si_9[k];

        t_19[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_15[k]
                  + pb_x[k] * si_13[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_x, pb_y, pb_z, sh0_8, sh0_9, sh1_10, \
                         sh1_12, si_14, si_16, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_x[k] * si_14[k];

        t_21[k] = pb_x[k] * si_20[k];

        t_22[k] = f_0 * sh0_8[k]
                  - f_1 * sh1_10[k]
                  + pb_y[k] * si_14[k];

        t_23[k] = pb_z[k] * si_14[k];

        t_24[k] = f_2 * sh0_9[k]
                  - f_3 * sh1_12[k]
                  + pb_y[k] * si_16[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_y, sh0_10, sh0_11, sh0_12, sh1_13, sh1_14, \
                         sh1_15, si_17, si_18, si_19, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * sh0_10[k]
                  - f_5 * sh1_13[k]
                  + pb_y[k] * si_17[k];

        t_26[k] = f_6 * sh0_11[k]
                  - f_7 * sh1_14[k]
                  + pb_y[k] * si_18[k];

        t_27[k] = f_8 * sh0_12[k]
                  - f_9 * sh1_15[k]
                  + pb_y[k] * si_19[k];

        t_28[k] = pb_y[k] * si_20[k];
    }

#pragma omp simd aligned(t_29, pb_z, sh0_12, sh1_15, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * sh0_12[k]
                  - f_1 * sh1_15[k]
                  + pb_z[k] * si_20[k];
    }
}

auto
compute_prim_sk_electron_repulsion_6(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_12 = buffer.data(sh0 + 12);
    const auto *sh0_13 = buffer.data(sh0 + 13);
    const auto *sh0_14 = buffer.data(sh0 + 14);
    const auto *sh0_15 = buffer.data(sh0 + 15);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);
    const auto *sh1_12 = buffer.data(sh1 + 12);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_3, sh0_4, sh1_0, sh1_1, sh1_2, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_3[k]
                 - f_3 * sh1_1[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_4[k]
                 - f_3 * sh1_2[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_5, sh0_6, sh0_7, sh1_3, sh1_4, sh1_5, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_5[k]
                 - f_5 * sh1_3[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_6[k]
                 - f_5 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_7[k]
                 - f_7 * sh1_5[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, sh0_8, sh0_9, sh0_10, sh1_6, sh1_7, sh1_8, si_6, \
                         si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_8[k]
                 - f_7 * sh1_6[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_6 * sh0_9[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_10[k]
                 - f_9 * sh1_8[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_x, sh0_12, sh0_13, sh0_15, sh1_9, sh1_10, \
                         sh1_12, si_9, si_10, si_11, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * sh0_12[k]
                 - f_9 * sh1_9[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_8 * sh0_13[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_10[k];

        t_11[k] = f_8 * sh0_15[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_11[k];

        t_12[k] = pb_x[k] * si_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, pb_x, pb_y, sh0_10, sh1_8, si_12, \
                         si_13, si_14, si_15, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_x[k] * si_13[k];

        t_14[k] = pb_x[k] * si_14[k];

        t_15[k] = pb_x[k] * si_15[k];

        t_16[k] = pb_x[k] * si_17[k];

        t_17[k] = f_0 * sh0_10[k]
                  - f_1 * sh1_8[k]
                  + pb_y[k] * si_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, sh0_12, sh0_13, sh0_14, sh1_9, sh1_10, \
                         sh1_11, si_13, si_14, si_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * sh0_12[k]
                  - f_3 * sh1_9[k]
                  + pb_y[k] * si_13[k];

        t_19[k] = f_4 * sh0_13[k]
                  - f_5 * sh1_10[k]
                  + pb_y[k] * si_14[k];

        t_20[k] = f_6 * sh0_14[k]
                  - f_7 * sh1_11[k]
                  + pb_y[k] * si_15[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, sh0_15, sh1_12, si_16, \
                         si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * sh0_15[k]
                  - f_9 * sh1_12[k]
                  + pb_y[k] * si_16[k];

        t_22[k] = pb_y[k] * si_17[k];

        t_23[k] = f_0 * sh0_15[k]
                  - f_1 * sh1_12[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_7(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_13 = buffer.data(sh0 + 13);
    const auto *sh0_14 = buffer.data(sh0 + 14);
    const auto *sh0_15 = buffer.data(sh0 + 15);
    const auto *sh0_16 = buffer.data(sh0 + 16);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_12 = buffer.data(sh1 + 12);
    const auto *sh1_13 = buffer.data(sh1 + 13);
    const auto *sh1_14 = buffer.data(sh1 + 14);
    const auto *sh1_15 = buffer.data(sh1 + 15);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_3, sh0_4, sh1_0, sh1_3, sh1_4, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_3[k]
                 - f_3 * sh1_3[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_4[k]
                 - f_3 * sh1_4[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_5, sh0_7, sh0_8, sh1_5, sh1_6, sh1_7, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_5[k]
                 - f_5 * sh1_5[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_7[k]
                 - f_5 * sh1_6[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_8[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, sh0_9, sh0_10, sh0_11, sh1_8, sh1_9, sh1_10, \
                         si_6, si_7, si_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_9[k]
                 - f_7 * sh1_8[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_6 * sh0_10[k]
                 - f_7 * sh1_9[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_11[k]
                 - f_9 * sh1_10[k]
                 + pb_x[k] * si_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, sh0_13, sh0_14, sh0_16, sh1_12, sh1_13, \
                         sh1_15, si_9, si_10, si_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * sh0_13[k]
                 - f_9 * sh1_12[k]
                 + pb_x[k] * si_9[k];

        t_10[k] = f_8 * sh0_14[k]
                  - f_9 * sh1_13[k]
                  + pb_x[k] * si_10[k];

        t_11[k] = f_8 * sh0_16[k]
                  - f_9 * sh1_15[k]
                  + pb_x[k] * si_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_y, sh0_11, sh0_13, sh0_14, sh1_10, sh1_12, \
                         sh1_13, si_12, si_13, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_0 * sh0_11[k]
                  - f_1 * sh1_10[k]
                  + pb_y[k] * si_12[k];

        t_13[k] = f_2 * sh0_13[k]
                  - f_3 * sh1_12[k]
                  + pb_y[k] * si_13[k];

        t_14[k] = f_4 * sh0_14[k]
                  - f_5 * sh1_13[k]
                  + pb_y[k] * si_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_y, pb_z, sh0_15, sh0_16, sh1_14, sh1_15, si_15, \
                         si_16, si_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * sh0_15[k]
                  - f_7 * sh1_14[k]
                  + pb_y[k] * si_15[k];

        t_16[k] = f_8 * sh0_16[k]
                  - f_9 * sh1_15[k]
                  + pb_y[k] * si_16[k];

        t_17[k] = f_0 * sh0_16[k]
                  - f_1 * sh1_15[k]
                  + pb_z[k] * si_17[k];
    }
}

auto
compute_prim_sk_electron_repulsion_8(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_10 = buffer.data(sh0 + 10);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_13 = buffer.data(sh0 + 13);
    const auto *sh0_14 = buffer.data(sh0 + 14);
    const auto *sh0_15 = buffer.data(sh0 + 15);
    const auto *sh0_16 = buffer.data(sh0 + 16);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_12 = buffer.data(sh1 + 12);
    const auto *sh1_13 = buffer.data(sh1 + 13);
    const auto *sh1_14 = buffer.data(sh1 + 14);
    const auto *sh1_15 = buffer.data(sh1 + 15);

    const auto *si_0 = buffer.data(si + 0);
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
    const auto *si_16 = buffer.data(si + 16);
    const auto *si_17 = buffer.data(si + 17);
    const auto *si_18 = buffer.data(si + 18);
    const auto *si_19 = buffer.data(si + 19);
    const auto *si_20 = buffer.data(si + 20);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, sh0_0, sh0_3, sh1_0, sh1_3, \
                         si_0, si_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = pb_y[k] * si_0[k];

        t_2[k] = pb_z[k] * si_0[k];

        t_3[k] = f_2 * sh0_3[k]
                 - f_3 * sh1_3[k]
                 + pb_x[k] * si_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pb_z, sh0_4, sh0_5, sh0_7, sh1_4, sh1_5, \
                         sh1_6, si_3, si_4, si_5, si_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * sh0_4[k]
                 - f_3 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_4 * sh0_5[k]
                 - f_5 * sh1_5[k]
                 + pb_x[k] * si_5[k];

        t_6[k] = pb_z[k] * si_3[k];

        t_7[k] = f_4 * sh0_7[k]
                 - f_5 * sh1_6[k]
                 + pb_x[k] * si_6[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pb_z, sh0_8, sh0_9, sh0_10, sh1_7, sh1_8, \
                         sh1_9, si_5, si_7, si_8, si_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_6 * sh0_8[k]
                 - f_7 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_9[k] = pb_z[k] * si_5[k];

        t_10[k] = f_6 * sh0_9[k]
                  - f_7 * sh1_8[k]
                  + pb_x[k] * si_8[k];

        t_11[k] = f_6 * sh0_10[k]
                  - f_7 * sh1_9[k]
                  + pb_x[k] * si_9[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_z, sh0_11, sh0_13, sh0_14, sh1_10, \
                         sh1_12, sh1_13, si_7, si_10, si_11, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_8 * sh0_11[k]
                  - f_9 * sh1_10[k]
                  + pb_x[k] * si_10[k];

        t_13[k] = pb_z[k] * si_7[k];

        t_14[k] = f_8 * sh0_13[k]
                  - f_9 * sh1_12[k]
                  + pb_x[k] * si_11[k];

        t_15[k] = f_8 * sh0_14[k]
                  - f_9 * sh1_13[k]
                  + pb_x[k] * si_12[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pb_y, pb_z, sh0_11, sh0_13, sh0_16, \
                         sh1_10, sh1_12, sh1_15, si_13, si_14, si_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_8 * sh0_16[k]
                  - f_9 * sh1_15[k]
                  + pb_x[k] * si_13[k];

        t_17[k] = f_0 * sh0_11[k]
                  - f_1 * sh1_10[k]
                  + pb_y[k] * si_14[k];

        t_18[k] = pb_z[k] * si_14[k];

        t_19[k] = f_2 * sh0_13[k]
                  - f_3 * sh1_12[k]
                  + pb_y[k] * si_16[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_y, pb_z, sh0_14, sh0_15, sh0_16, sh1_13, \
                         sh1_14, sh1_15, si_17, si_18, si_19, si_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_4 * sh0_14[k]
                  - f_5 * sh1_13[k]
                  + pb_y[k] * si_17[k];

        t_21[k] = f_6 * sh0_15[k]
                  - f_7 * sh1_14[k]
                  + pb_y[k] * si_18[k];

        t_22[k] = f_8 * sh0_16[k]
                  - f_9 * sh1_15[k]
                  + pb_y[k] * si_19[k];

        t_23[k] = f_0 * sh0_16[k]
                  - f_1 * sh1_15[k]
                  + pb_z[k] * si_20[k];
    }
}

auto
compute_prim_sk_electron_repulsion_9(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t sh0, const size_t sh1, const size_t si,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_1 = buffer.data(si + 1);
    const auto *si_2 = buffer.data(si + 2);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sh0_0, sh0_1, sh0_2, sh1_0, sh1_1, \
                         sh1_2, si_0, si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_0 * sh0_1[k]
                 - f_1 * sh1_1[k]
                 + pb_y[k] * si_1[k];

        t_2[k] = f_0 * sh0_2[k]
                 - f_1 * sh1_2[k]
                 + pb_z[k] * si_2[k];
    }
}

auto
compute_prim_sk_electron_repulsion_10(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sh0, const size_t sh1, const size_t si,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);

    auto *t_0 = buffer.data(target + 0);
    auto *t_1 = buffer.data(target + 1);
    auto *t_2 = buffer.data(target + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_1 = buffer.data(sh0 + 1);
    const auto *sh0_2 = buffer.data(sh0 + 2);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_14 = buffer.data(sh1 + 14);

    const auto *si_0 = buffer.data(si + 0);
    const auto *si_9 = buffer.data(si + 9);
    const auto *si_14 = buffer.data(si + 14);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sh0_0, sh0_1, sh0_2, sh1_0, sh1_9, \
                         sh1_14, si_0, si_9, si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_0 * sh0_1[k]
                 - f_1 * sh1_9[k]
                 + pb_y[k] * si_9[k];

        t_2[k] = f_0 * sh0_2[k]
                 - f_1 * sh1_14[k]
                 + pb_z[k] * si_14[k];
    }
}

auto
compute_prim_sk_electron_repulsion_11(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                      const size_t sh0, const size_t sh1, const size_t si,
                                      const size_t ncols, const double alpha, const double beta,
                                      const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / beta;
    const auto f_1 = 3.0 * alpha / (beta * p);
    const auto f_2 = 2.0 / beta;
    const auto f_3 = 2.0 * alpha / (beta * p);
    const auto f_4 = 1.5 / beta;
    const auto f_5 = 1.5 * alpha / (beta * p);
    const auto f_6 = 1.0 / beta;
    const auto f_7 = alpha / (beta * p);
    const auto f_8 = 0.5 / beta;
    const auto f_9 = 0.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sh0_0 = buffer.data(sh0 + 0);
    const auto *sh0_3 = buffer.data(sh0 + 3);
    const auto *sh0_4 = buffer.data(sh0 + 4);
    const auto *sh0_5 = buffer.data(sh0 + 5);
    const auto *sh0_6 = buffer.data(sh0 + 6);
    const auto *sh0_7 = buffer.data(sh0 + 7);
    const auto *sh0_8 = buffer.data(sh0 + 8);
    const auto *sh0_9 = buffer.data(sh0 + 9);
    const auto *sh0_11 = buffer.data(sh0 + 11);
    const auto *sh0_12 = buffer.data(sh0 + 12);
    const auto *sh0_13 = buffer.data(sh0 + 13);
    const auto *sh0_14 = buffer.data(sh0 + 14);

    const auto *sh1_0 = buffer.data(sh1 + 0);
    const auto *sh1_1 = buffer.data(sh1 + 1);
    const auto *sh1_2 = buffer.data(sh1 + 2);
    const auto *sh1_3 = buffer.data(sh1 + 3);
    const auto *sh1_4 = buffer.data(sh1 + 4);
    const auto *sh1_5 = buffer.data(sh1 + 5);
    const auto *sh1_6 = buffer.data(sh1 + 6);
    const auto *sh1_7 = buffer.data(sh1 + 7);
    const auto *sh1_8 = buffer.data(sh1 + 8);
    const auto *sh1_9 = buffer.data(sh1 + 9);
    const auto *sh1_10 = buffer.data(sh1 + 10);
    const auto *sh1_11 = buffer.data(sh1 + 11);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, sh0_0, sh0_3, sh0_4, sh1_0, sh1_1, sh1_2, si_0, \
                         si_1, si_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sh0_0[k]
                 - f_1 * sh1_0[k]
                 + pb_x[k] * si_0[k];

        t_1[k] = f_2 * sh0_3[k]
                 - f_3 * sh1_1[k]
                 + pb_x[k] * si_1[k];

        t_2[k] = f_2 * sh0_4[k]
                 - f_3 * sh1_2[k]
                 + pb_x[k] * si_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, sh0_5, sh0_6, sh0_7, sh1_3, sh1_4, sh1_5, si_3, \
                         si_4, si_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * sh0_5[k]
                 - f_5 * sh1_3[k]
                 + pb_x[k] * si_3[k];

        t_4[k] = f_4 * sh0_6[k]
                 - f_5 * sh1_4[k]
                 + pb_x[k] * si_4[k];

        t_5[k] = f_6 * sh0_7[k]
                 - f_7 * sh1_5[k]
                 + pb_x[k] * si_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_x, pb_y, sh0_8, sh0_9, sh0_14, sh1_6, sh1_7, \
                         sh1_11, si_6, si_7, si_8, si_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * sh0_8[k]
                 - f_7 * sh1_6[k]
                 + pb_x[k] * si_6[k];

        t_7[k] = f_8 * sh0_9[k]
                 - f_9 * sh1_7[k]
                 + pb_x[k] * si_7[k];

        t_8[k] = f_8 * sh0_14[k]
                 - f_9 * sh1_11[k]
                 + pb_x[k] * si_8[k];

        t_9[k] = f_0 * sh0_9[k]
                 - f_1 * sh1_7[k]
                 + pb_y[k] * si_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, sh0_11, sh0_12, sh0_13, sh1_8, sh1_9, sh1_10, \
                         si_10, si_11, si_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * sh0_11[k]
                  - f_3 * sh1_8[k]
                  + pb_y[k] * si_10[k];

        t_11[k] = f_4 * sh0_12[k]
                  - f_5 * sh1_9[k]
                  + pb_y[k] * si_11[k];

        t_12[k] = f_6 * sh0_13[k]
                  - f_7 * sh1_10[k]
                  + pb_y[k] * si_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_y, pb_z, sh0_14, sh1_11, si_13, \
                         si_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_8 * sh0_14[k]
                  - f_9 * sh1_11[k]
                  + pb_y[k] * si_13[k];

        t_14[k] = f_0 * sh0_14[k]
                  - f_1 * sh1_11[k]
                  + pb_z[k] * si_14[k];
    }
}

}  // namespace simdt2ceri
