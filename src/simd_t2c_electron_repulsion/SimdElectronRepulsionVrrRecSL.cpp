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


#include "SimdElectronRepulsionVrrRecSL.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_sl_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t si0, const size_t si1, const size_t sk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / beta;
    const auto f_1 = 3.5 * alpha / (beta * p);
    const auto f_2 = 2.5 / beta;
    const auto f_3 = 2.5 * alpha / (beta * p);
    const auto f_4 = 2.0 / beta;
    const auto f_5 = 2.0 * alpha / (beta * p);
    const auto f_6 = 1.5 / beta;
    const auto f_7 = 1.5 * alpha / (beta * p);
    const auto f_8 = 1.0 / beta;
    const auto f_9 = alpha / (beta * p);
    const auto f_10 = 0.5 / beta;
    const auto f_11 = 0.5 * alpha / (beta * p);

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

    const auto *si0_0 = buffer.data(si0 + 0);
    const auto *si0_1 = buffer.data(si0 + 1);
    const auto *si0_2 = buffer.data(si0 + 2);
    const auto *si0_3 = buffer.data(si0 + 3);
    const auto *si0_4 = buffer.data(si0 + 4);
    const auto *si0_5 = buffer.data(si0 + 5);
    const auto *si0_6 = buffer.data(si0 + 6);
    const auto *si0_7 = buffer.data(si0 + 7);
    const auto *si0_8 = buffer.data(si0 + 8);
    const auto *si0_9 = buffer.data(si0 + 9);
    const auto *si0_10 = buffer.data(si0 + 10);
    const auto *si0_11 = buffer.data(si0 + 11);
    const auto *si0_12 = buffer.data(si0 + 12);
    const auto *si0_13 = buffer.data(si0 + 13);
    const auto *si0_14 = buffer.data(si0 + 14);
    const auto *si0_15 = buffer.data(si0 + 15);
    const auto *si0_16 = buffer.data(si0 + 16);
    const auto *si0_17 = buffer.data(si0 + 17);

    const auto *si1_0 = buffer.data(si1 + 0);
    const auto *si1_3 = buffer.data(si1 + 3);
    const auto *si1_4 = buffer.data(si1 + 4);
    const auto *si1_5 = buffer.data(si1 + 5);
    const auto *si1_6 = buffer.data(si1 + 6);
    const auto *si1_7 = buffer.data(si1 + 7);
    const auto *si1_8 = buffer.data(si1 + 8);
    const auto *si1_9 = buffer.data(si1 + 9);
    const auto *si1_10 = buffer.data(si1 + 10);
    const auto *si1_11 = buffer.data(si1 + 11);
    const auto *si1_12 = buffer.data(si1 + 12);
    const auto *si1_13 = buffer.data(si1 + 13);
    const auto *si1_14 = buffer.data(si1 + 14);
    const auto *si1_16 = buffer.data(si1 + 16);
    const auto *si1_17 = buffer.data(si1 + 17);
    const auto *si1_18 = buffer.data(si1 + 18);
    const auto *si1_19 = buffer.data(si1 + 19);
    const auto *si1_20 = buffer.data(si1 + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, si0_0, si0_1, si0_2, si1_0, si1_3, si1_4, sk_0, \
                         sk_1, sk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si0_0[k]
                 - f_1 * si1_0[k]
                 + pb_x[k] * sk_0[k];

        t_1[k] = f_2 * si0_1[k]
                 - f_3 * si1_3[k]
                 + pb_x[k] * sk_1[k];

        t_2[k] = f_2 * si0_2[k]
                 - f_3 * si1_4[k]
                 + pb_x[k] * sk_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, si0_3, si0_4, si0_5, si1_5, si1_6, si1_7, sk_3, \
                         sk_4, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * si0_3[k]
                 - f_5 * si1_5[k]
                 + pb_x[k] * sk_3[k];

        t_4[k] = f_4 * si0_4[k]
                 - f_5 * si1_6[k]
                 + pb_x[k] * sk_4[k];

        t_5[k] = f_6 * si0_5[k]
                 - f_7 * si1_7[k]
                 + pb_x[k] * sk_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, si0_6, si0_7, si0_8, si1_8, si1_9, si1_10, sk_6, \
                         sk_7, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * si0_6[k]
                 - f_7 * si1_8[k]
                 + pb_x[k] * sk_6[k];

        t_7[k] = f_6 * si0_7[k]
                 - f_7 * si1_9[k]
                 + pb_x[k] * sk_7[k];

        t_8[k] = f_8 * si0_8[k]
                 - f_9 * si1_10[k]
                 + pb_x[k] * sk_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, si0_9, si0_10, si0_11, si1_11, si1_12, si1_13, \
                         sk_9, sk_10, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * si0_9[k]
                 - f_9 * si1_11[k]
                 + pb_x[k] * sk_9[k];

        t_10[k] = f_8 * si0_10[k]
                  - f_9 * si1_12[k]
                  + pb_x[k] * sk_10[k];

        t_11[k] = f_8 * si0_11[k]
                  - f_9 * si1_13[k]
                  + pb_x[k] * sk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, si0_12, si0_13, si0_14, si1_14, si1_16, \
                         si1_17, sk_12, sk_13, sk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_10 * si0_12[k]
                  - f_11 * si1_14[k]
                  + pb_x[k] * sk_12[k];

        t_13[k] = f_10 * si0_13[k]
                  - f_11 * si1_16[k]
                  + pb_x[k] * sk_13[k];

        t_14[k] = f_10 * si0_14[k]
                  - f_11 * si1_17[k]
                  + pb_x[k] * sk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, si0_12, si0_15, si0_17, si1_14, si1_18, \
                         si1_20, sk_15, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * si0_15[k]
                  - f_11 * si1_18[k]
                  + pb_x[k] * sk_15[k];

        t_16[k] = f_10 * si0_17[k]
                  - f_11 * si1_20[k]
                  + pb_x[k] * sk_16[k];

        t_17[k] = f_0 * si0_12[k]
                  - f_1 * si1_14[k]
                  + pb_y[k] * sk_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, si0_13, si0_14, si0_15, si1_16, si1_17, \
                         si1_18, sk_18, sk_19, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * si0_13[k]
                  - f_3 * si1_16[k]
                  + pb_y[k] * sk_18[k];

        t_19[k] = f_4 * si0_14[k]
                  - f_5 * si1_17[k]
                  + pb_y[k] * sk_19[k];

        t_20[k] = f_6 * si0_15[k]
                  - f_7 * si1_18[k]
                  + pb_y[k] * sk_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, si0_16, si0_17, si1_19, si1_20, sk_21, \
                         sk_22, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * si0_16[k]
                  - f_9 * si1_19[k]
                  + pb_y[k] * sk_21[k];

        t_22[k] = f_10 * si0_17[k]
                  - f_11 * si1_20[k]
                  + pb_y[k] * sk_22[k];

        t_23[k] = f_0 * si0_17[k]
                  - f_1 * si1_20[k]
                  + pb_z[k] * sk_23[k];
    }
}

auto
compute_prim_sl_electron_repulsion_1(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t si0, const size_t si1, const size_t sk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / beta;
    const auto f_1 = 3.5 * alpha / (beta * p);
    const auto f_2 = 2.5 / beta;
    const auto f_3 = 2.5 * alpha / (beta * p);
    const auto f_4 = 2.0 / beta;
    const auto f_5 = 2.0 * alpha / (beta * p);
    const auto f_6 = 1.5 / beta;
    const auto f_7 = 1.5 * alpha / (beta * p);
    const auto f_8 = 1.0 / beta;
    const auto f_9 = alpha / (beta * p);
    const auto f_10 = 0.5 / beta;
    const auto f_11 = 0.5 * alpha / (beta * p);

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

    const auto *si0_0 = buffer.data(si0 + 0);
    const auto *si0_3 = buffer.data(si0 + 3);
    const auto *si0_4 = buffer.data(si0 + 4);
    const auto *si0_5 = buffer.data(si0 + 5);
    const auto *si0_8 = buffer.data(si0 + 8);
    const auto *si0_9 = buffer.data(si0 + 9);
    const auto *si0_11 = buffer.data(si0 + 11);
    const auto *si0_13 = buffer.data(si0 + 13);
    const auto *si0_14 = buffer.data(si0 + 14);
    const auto *si0_15 = buffer.data(si0 + 15);
    const auto *si0_16 = buffer.data(si0 + 16);
    const auto *si0_17 = buffer.data(si0 + 17);
    const auto *si0_18 = buffer.data(si0 + 18);
    const auto *si0_20 = buffer.data(si0 + 20);
    const auto *si0_21 = buffer.data(si0 + 21);
    const auto *si0_22 = buffer.data(si0 + 22);
    const auto *si0_23 = buffer.data(si0 + 23);
    const auto *si0_24 = buffer.data(si0 + 24);

    const auto *si1_0 = buffer.data(si1 + 0);
    const auto *si1_3 = buffer.data(si1 + 3);
    const auto *si1_4 = buffer.data(si1 + 4);
    const auto *si1_5 = buffer.data(si1 + 5);
    const auto *si1_6 = buffer.data(si1 + 6);
    const auto *si1_7 = buffer.data(si1 + 7);
    const auto *si1_8 = buffer.data(si1 + 8);
    const auto *si1_9 = buffer.data(si1 + 9);
    const auto *si1_10 = buffer.data(si1 + 10);
    const auto *si1_11 = buffer.data(si1 + 11);
    const auto *si1_12 = buffer.data(si1 + 12);
    const auto *si1_13 = buffer.data(si1 + 13);
    const auto *si1_14 = buffer.data(si1 + 14);
    const auto *si1_16 = buffer.data(si1 + 16);
    const auto *si1_17 = buffer.data(si1 + 17);
    const auto *si1_18 = buffer.data(si1 + 18);
    const auto *si1_19 = buffer.data(si1 + 19);
    const auto *si1_20 = buffer.data(si1 + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, si0_0, si0_3, si0_4, si1_0, si1_3, si1_4, sk_0, \
                         sk_1, sk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si0_0[k]
                 - f_1 * si1_0[k]
                 + pb_x[k] * sk_0[k];

        t_1[k] = f_2 * si0_3[k]
                 - f_3 * si1_3[k]
                 + pb_x[k] * sk_1[k];

        t_2[k] = f_2 * si0_4[k]
                 - f_3 * si1_4[k]
                 + pb_x[k] * sk_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, si0_5, si0_8, si0_9, si1_5, si1_6, si1_7, sk_3, \
                         sk_4, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * si0_5[k]
                 - f_5 * si1_5[k]
                 + pb_x[k] * sk_3[k];

        t_4[k] = f_4 * si0_8[k]
                 - f_5 * si1_6[k]
                 + pb_x[k] * sk_4[k];

        t_5[k] = f_6 * si0_9[k]
                 - f_7 * si1_7[k]
                 + pb_x[k] * sk_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, si0_11, si0_13, si0_14, si1_8, si1_9, si1_10, \
                         sk_6, sk_7, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * si0_11[k]
                 - f_7 * si1_8[k]
                 + pb_x[k] * sk_6[k];

        t_7[k] = f_6 * si0_13[k]
                 - f_7 * si1_9[k]
                 + pb_x[k] * sk_7[k];

        t_8[k] = f_8 * si0_14[k]
                 - f_9 * si1_10[k]
                 + pb_x[k] * sk_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, si0_15, si0_16, si0_17, si1_11, si1_12, \
                         si1_13, sk_9, sk_10, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * si0_15[k]
                 - f_9 * si1_11[k]
                 + pb_x[k] * sk_9[k];

        t_10[k] = f_8 * si0_16[k]
                  - f_9 * si1_12[k]
                  + pb_x[k] * sk_10[k];

        t_11[k] = f_8 * si0_17[k]
                  - f_9 * si1_13[k]
                  + pb_x[k] * sk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, si0_18, si0_20, si0_21, si1_14, si1_16, \
                         si1_17, sk_12, sk_13, sk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_10 * si0_18[k]
                  - f_11 * si1_14[k]
                  + pb_x[k] * sk_12[k];

        t_13[k] = f_10 * si0_20[k]
                  - f_11 * si1_16[k]
                  + pb_x[k] * sk_13[k];

        t_14[k] = f_10 * si0_21[k]
                  - f_11 * si1_17[k]
                  + pb_x[k] * sk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, si0_18, si0_22, si0_24, si1_14, si1_18, \
                         si1_20, sk_15, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * si0_22[k]
                  - f_11 * si1_18[k]
                  + pb_x[k] * sk_15[k];

        t_16[k] = f_10 * si0_24[k]
                  - f_11 * si1_20[k]
                  + pb_x[k] * sk_16[k];

        t_17[k] = f_0 * si0_18[k]
                  - f_1 * si1_14[k]
                  + pb_y[k] * sk_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, si0_20, si0_21, si0_22, si1_16, si1_17, \
                         si1_18, sk_18, sk_19, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * si0_20[k]
                  - f_3 * si1_16[k]
                  + pb_y[k] * sk_18[k];

        t_19[k] = f_4 * si0_21[k]
                  - f_5 * si1_17[k]
                  + pb_y[k] * sk_19[k];

        t_20[k] = f_6 * si0_22[k]
                  - f_7 * si1_18[k]
                  + pb_y[k] * sk_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, si0_23, si0_24, si1_19, si1_20, sk_21, \
                         sk_22, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * si0_23[k]
                  - f_9 * si1_19[k]
                  + pb_y[k] * sk_21[k];

        t_22[k] = f_10 * si0_24[k]
                  - f_11 * si1_20[k]
                  + pb_y[k] * sk_22[k];

        t_23[k] = f_0 * si0_24[k]
                  - f_1 * si1_20[k]
                  + pb_z[k] * sk_23[k];
    }
}

auto
compute_prim_sl_electron_repulsion_2(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                     const size_t si0, const size_t si1, const size_t sk,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / beta;
    const auto f_1 = 3.5 * alpha / (beta * p);
    const auto f_2 = 2.5 / beta;
    const auto f_3 = 2.5 * alpha / (beta * p);
    const auto f_4 = 2.0 / beta;
    const auto f_5 = 2.0 * alpha / (beta * p);
    const auto f_6 = 1.5 / beta;
    const auto f_7 = 1.5 * alpha / (beta * p);
    const auto f_8 = 1.0 / beta;
    const auto f_9 = alpha / (beta * p);
    const auto f_10 = 0.5 / beta;
    const auto f_11 = 0.5 * alpha / (beta * p);

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

    const auto *si0_0 = buffer.data(si0 + 0);
    const auto *si0_3 = buffer.data(si0 + 3);
    const auto *si0_4 = buffer.data(si0 + 4);
    const auto *si0_5 = buffer.data(si0 + 5);
    const auto *si0_7 = buffer.data(si0 + 7);
    const auto *si0_8 = buffer.data(si0 + 8);
    const auto *si0_10 = buffer.data(si0 + 10);
    const auto *si0_11 = buffer.data(si0 + 11);
    const auto *si0_12 = buffer.data(si0 + 12);
    const auto *si0_13 = buffer.data(si0 + 13);
    const auto *si0_14 = buffer.data(si0 + 14);
    const auto *si0_15 = buffer.data(si0 + 15);
    const auto *si0_16 = buffer.data(si0 + 16);
    const auto *si0_18 = buffer.data(si0 + 18);
    const auto *si0_19 = buffer.data(si0 + 19);
    const auto *si0_20 = buffer.data(si0 + 20);
    const auto *si0_21 = buffer.data(si0 + 21);
    const auto *si0_22 = buffer.data(si0 + 22);

    const auto *si1_0 = buffer.data(si1 + 0);
    const auto *si1_3 = buffer.data(si1 + 3);
    const auto *si1_4 = buffer.data(si1 + 4);
    const auto *si1_5 = buffer.data(si1 + 5);
    const auto *si1_6 = buffer.data(si1 + 6);
    const auto *si1_7 = buffer.data(si1 + 7);
    const auto *si1_8 = buffer.data(si1 + 8);
    const auto *si1_9 = buffer.data(si1 + 9);
    const auto *si1_10 = buffer.data(si1 + 10);
    const auto *si1_11 = buffer.data(si1 + 11);
    const auto *si1_12 = buffer.data(si1 + 12);
    const auto *si1_13 = buffer.data(si1 + 13);
    const auto *si1_14 = buffer.data(si1 + 14);
    const auto *si1_16 = buffer.data(si1 + 16);
    const auto *si1_17 = buffer.data(si1 + 17);
    const auto *si1_18 = buffer.data(si1 + 18);
    const auto *si1_19 = buffer.data(si1 + 19);
    const auto *si1_20 = buffer.data(si1 + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, si0_0, si0_3, si0_4, si1_0, si1_3, si1_4, sk_0, \
                         sk_1, sk_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * si0_0[k]
                 - f_1 * si1_0[k]
                 + pb_x[k] * sk_0[k];

        t_1[k] = f_2 * si0_3[k]
                 - f_3 * si1_3[k]
                 + pb_x[k] * sk_1[k];

        t_2[k] = f_2 * si0_4[k]
                 - f_3 * si1_4[k]
                 + pb_x[k] * sk_2[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, si0_5, si0_7, si0_8, si1_5, si1_6, si1_7, sk_3, \
                         sk_4, sk_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_4 * si0_5[k]
                 - f_5 * si1_5[k]
                 + pb_x[k] * sk_3[k];

        t_4[k] = f_4 * si0_7[k]
                 - f_5 * si1_6[k]
                 + pb_x[k] * sk_4[k];

        t_5[k] = f_6 * si0_8[k]
                 - f_7 * si1_7[k]
                 + pb_x[k] * sk_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, si0_10, si0_11, si0_12, si1_8, si1_9, si1_10, \
                         sk_6, sk_7, sk_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * si0_10[k]
                 - f_7 * si1_8[k]
                 + pb_x[k] * sk_6[k];

        t_7[k] = f_6 * si0_11[k]
                 - f_7 * si1_9[k]
                 + pb_x[k] * sk_7[k];

        t_8[k] = f_8 * si0_12[k]
                 - f_9 * si1_10[k]
                 + pb_x[k] * sk_8[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, si0_13, si0_14, si0_15, si1_11, si1_12, \
                         si1_13, sk_9, sk_10, sk_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_8 * si0_13[k]
                 - f_9 * si1_11[k]
                 + pb_x[k] * sk_9[k];

        t_10[k] = f_8 * si0_14[k]
                  - f_9 * si1_12[k]
                  + pb_x[k] * sk_10[k];

        t_11[k] = f_8 * si0_15[k]
                  - f_9 * si1_13[k]
                  + pb_x[k] * sk_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, si0_16, si0_18, si0_19, si1_14, si1_16, \
                         si1_17, sk_12, sk_13, sk_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_10 * si0_16[k]
                  - f_11 * si1_14[k]
                  + pb_x[k] * sk_12[k];

        t_13[k] = f_10 * si0_18[k]
                  - f_11 * si1_16[k]
                  + pb_x[k] * sk_13[k];

        t_14[k] = f_10 * si0_19[k]
                  - f_11 * si1_17[k]
                  + pb_x[k] * sk_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pb_y, si0_16, si0_20, si0_22, si1_14, si1_18, \
                         si1_20, sk_15, sk_16, sk_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_10 * si0_20[k]
                  - f_11 * si1_18[k]
                  + pb_x[k] * sk_15[k];

        t_16[k] = f_10 * si0_22[k]
                  - f_11 * si1_20[k]
                  + pb_x[k] * sk_16[k];

        t_17[k] = f_0 * si0_16[k]
                  - f_1 * si1_14[k]
                  + pb_y[k] * sk_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pb_y, si0_18, si0_19, si0_20, si1_16, si1_17, \
                         si1_18, sk_18, sk_19, sk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_2 * si0_18[k]
                  - f_3 * si1_16[k]
                  + pb_y[k] * sk_18[k];

        t_19[k] = f_4 * si0_19[k]
                  - f_5 * si1_17[k]
                  + pb_y[k] * sk_19[k];

        t_20[k] = f_6 * si0_20[k]
                  - f_7 * si1_18[k]
                  + pb_y[k] * sk_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, pb_y, pb_z, si0_21, si0_22, si1_19, si1_20, sk_21, \
                         sk_22, sk_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * si0_21[k]
                  - f_9 * si1_19[k]
                  + pb_y[k] * sk_21[k];

        t_22[k] = f_10 * si0_22[k]
                  - f_11 * si1_20[k]
                  + pb_y[k] * sk_22[k];

        t_23[k] = f_0 * si0_22[k]
                  - f_1 * si1_20[k]
                  + pb_z[k] * sk_23[k];
    }
}

}  // namespace simdt2ceri
