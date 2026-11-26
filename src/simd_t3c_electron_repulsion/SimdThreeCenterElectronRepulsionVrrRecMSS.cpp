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


#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_mss_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t kss0, const size_t kss1,
                                                   const size_t lss0, const size_t lss1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / p;
    const auto f_1 = 4.0 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 3.0 / p;
    const auto f_4 = 3.0 * gamma / (p * q);
    const auto f_5 = 2.5 / p;
    const auto f_6 = 2.5 * gamma / (p * q);
    const auto f_7 = 2.0 / p;
    const auto f_8 = 2.0 * gamma / (p * q);
    const auto f_9 = 1.5 / p;
    const auto f_10 = 1.5 * gamma / (p * q);
    const auto f_11 = 1.0 / p;
    const auto f_12 = gamma / (p * q);
    const auto f_13 = 0.5 / p;
    const auto f_14 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *kss0_0 = buffer.data(kss0 + 0);
    const auto *kss0_3 = buffer.data(kss0 + 3);
    const auto *kss0_5 = buffer.data(kss0 + 5);
    const auto *kss0_6 = buffer.data(kss0 + 6);
    const auto *kss0_9 = buffer.data(kss0 + 9);
    const auto *kss0_10 = buffer.data(kss0 + 10);
    const auto *kss0_12 = buffer.data(kss0 + 12);
    const auto *kss0_14 = buffer.data(kss0 + 14);
    const auto *kss0_15 = buffer.data(kss0 + 15);
    const auto *kss0_17 = buffer.data(kss0 + 17);
    const auto *kss0_18 = buffer.data(kss0 + 18);
    const auto *kss0_20 = buffer.data(kss0 + 20);
    const auto *kss0_21 = buffer.data(kss0 + 21);
    const auto *kss0_23 = buffer.data(kss0 + 23);
    const auto *kss0_24 = buffer.data(kss0 + 24);
    const auto *kss0_25 = buffer.data(kss0 + 25);
    const auto *kss0_27 = buffer.data(kss0 + 27);
    const auto *kss0_28 = buffer.data(kss0 + 28);
    const auto *kss0_30 = buffer.data(kss0 + 30);
    const auto *kss0_31 = buffer.data(kss0 + 31);
    const auto *kss0_32 = buffer.data(kss0 + 32);
    const auto *kss0_33 = buffer.data(kss0 + 33);
    const auto *kss0_34 = buffer.data(kss0 + 34);
    const auto *kss0_35 = buffer.data(kss0 + 35);

    const auto *kss1_0 = buffer.data(kss1 + 0);
    const auto *kss1_3 = buffer.data(kss1 + 3);
    const auto *kss1_5 = buffer.data(kss1 + 5);
    const auto *kss1_6 = buffer.data(kss1 + 6);
    const auto *kss1_9 = buffer.data(kss1 + 9);
    const auto *kss1_10 = buffer.data(kss1 + 10);
    const auto *kss1_12 = buffer.data(kss1 + 12);
    const auto *kss1_14 = buffer.data(kss1 + 14);
    const auto *kss1_15 = buffer.data(kss1 + 15);
    const auto *kss1_17 = buffer.data(kss1 + 17);
    const auto *kss1_18 = buffer.data(kss1 + 18);
    const auto *kss1_20 = buffer.data(kss1 + 20);
    const auto *kss1_21 = buffer.data(kss1 + 21);
    const auto *kss1_23 = buffer.data(kss1 + 23);
    const auto *kss1_24 = buffer.data(kss1 + 24);
    const auto *kss1_25 = buffer.data(kss1 + 25);
    const auto *kss1_27 = buffer.data(kss1 + 27);
    const auto *kss1_28 = buffer.data(kss1 + 28);
    const auto *kss1_30 = buffer.data(kss1 + 30);
    const auto *kss1_31 = buffer.data(kss1 + 31);
    const auto *kss1_32 = buffer.data(kss1 + 32);
    const auto *kss1_33 = buffer.data(kss1 + 33);
    const auto *kss1_34 = buffer.data(kss1 + 34);
    const auto *kss1_35 = buffer.data(kss1 + 35);

    const auto *lss0_0 = buffer.data(lss0 + 0);
    const auto *lss0_2 = buffer.data(lss0 + 2);
    const auto *lss0_3 = buffer.data(lss0 + 3);
    const auto *lss0_5 = buffer.data(lss0 + 5);
    const auto *lss0_6 = buffer.data(lss0 + 6);
    const auto *lss0_9 = buffer.data(lss0 + 9);
    const auto *lss0_10 = buffer.data(lss0 + 10);
    const auto *lss0_12 = buffer.data(lss0 + 12);
    const auto *lss0_14 = buffer.data(lss0 + 14);
    const auto *lss0_15 = buffer.data(lss0 + 15);
    const auto *lss0_17 = buffer.data(lss0 + 17);
    const auto *lss0_18 = buffer.data(lss0 + 18);
    const auto *lss0_20 = buffer.data(lss0 + 20);
    const auto *lss0_21 = buffer.data(lss0 + 21);
    const auto *lss0_23 = buffer.data(lss0 + 23);
    const auto *lss0_24 = buffer.data(lss0 + 24);
    const auto *lss0_25 = buffer.data(lss0 + 25);
    const auto *lss0_27 = buffer.data(lss0 + 27);
    const auto *lss0_28 = buffer.data(lss0 + 28);
    const auto *lss0_30 = buffer.data(lss0 + 30);
    const auto *lss0_31 = buffer.data(lss0 + 31);
    const auto *lss0_32 = buffer.data(lss0 + 32);
    const auto *lss0_33 = buffer.data(lss0 + 33);
    const auto *lss0_35 = buffer.data(lss0 + 35);
    const auto *lss0_36 = buffer.data(lss0 + 36);
    const auto *lss0_37 = buffer.data(lss0 + 37);
    const auto *lss0_38 = buffer.data(lss0 + 38);
    const auto *lss0_39 = buffer.data(lss0 + 39);
    const auto *lss0_40 = buffer.data(lss0 + 40);
    const auto *lss0_41 = buffer.data(lss0 + 41);
    const auto *lss0_42 = buffer.data(lss0 + 42);
    const auto *lss0_43 = buffer.data(lss0 + 43);
    const auto *lss0_44 = buffer.data(lss0 + 44);

    const auto *lss1_0 = buffer.data(lss1 + 0);
    const auto *lss1_2 = buffer.data(lss1 + 2);
    const auto *lss1_3 = buffer.data(lss1 + 3);
    const auto *lss1_5 = buffer.data(lss1 + 5);
    const auto *lss1_6 = buffer.data(lss1 + 6);
    const auto *lss1_9 = buffer.data(lss1 + 9);
    const auto *lss1_10 = buffer.data(lss1 + 10);
    const auto *lss1_12 = buffer.data(lss1 + 12);
    const auto *lss1_14 = buffer.data(lss1 + 14);
    const auto *lss1_15 = buffer.data(lss1 + 15);
    const auto *lss1_17 = buffer.data(lss1 + 17);
    const auto *lss1_18 = buffer.data(lss1 + 18);
    const auto *lss1_20 = buffer.data(lss1 + 20);
    const auto *lss1_21 = buffer.data(lss1 + 21);
    const auto *lss1_23 = buffer.data(lss1 + 23);
    const auto *lss1_24 = buffer.data(lss1 + 24);
    const auto *lss1_25 = buffer.data(lss1 + 25);
    const auto *lss1_27 = buffer.data(lss1 + 27);
    const auto *lss1_28 = buffer.data(lss1 + 28);
    const auto *lss1_30 = buffer.data(lss1 + 30);
    const auto *lss1_31 = buffer.data(lss1 + 31);
    const auto *lss1_32 = buffer.data(lss1 + 32);
    const auto *lss1_33 = buffer.data(lss1 + 33);
    const auto *lss1_35 = buffer.data(lss1 + 35);
    const auto *lss1_36 = buffer.data(lss1 + 36);
    const auto *lss1_37 = buffer.data(lss1 + 37);
    const auto *lss1_38 = buffer.data(lss1 + 38);
    const auto *lss1_39 = buffer.data(lss1 + 39);
    const auto *lss1_40 = buffer.data(lss1 + 40);
    const auto *lss1_41 = buffer.data(lss1 + 41);
    const auto *lss1_42 = buffer.data(lss1 + 42);
    const auto *lss1_43 = buffer.data(lss1 + 43);
    const auto *lss1_44 = buffer.data(lss1 + 44);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, kss0_0, kss1_0, \
                         lss0_0, lss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * kss0_0[k]
                 - f_1 * kss1_0[k]
                 + pa_x[k] * lss0_0[k]
                 - f_2 * pc_x[k] * lss1_0[k];

        t_1[k] = pa_y[k] * lss0_0[k]
                 - f_2 * pc_y[k] * lss1_0[k];

        t_2[k] = pa_z[k] * lss0_0[k]
                 - f_2 * pc_z[k] * lss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, kss0_3, kss1_3, lss0_2, lss0_3, \
                         lss1_2, lss1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * kss0_3[k]
                 - f_4 * kss1_3[k]
                 + pa_x[k] * lss0_3[k]
                 - f_2 * pc_x[k] * lss1_3[k];

        t_4[k] = pa_y[k] * lss0_2[k]
                 - f_2 * pc_y[k] * lss1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, kss0_5, kss0_6, kss1_5, kss1_6, lss0_5, lss0_6, \
                         lss1_5, lss1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * kss0_5[k]
                 - f_4 * kss1_5[k]
                 + pa_x[k] * lss0_5[k]
                 - f_2 * pc_x[k] * lss1_5[k];

        t_6[k] = f_5 * kss0_6[k]
                 - f_6 * kss1_6[k]
                 + pa_x[k] * lss0_6[k]
                 - f_2 * pc_x[k] * lss1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, lss0_3, lss0_5, lss1_3, \
                         lss1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * lss0_3[k]
                 - f_2 * pc_z[k] * lss1_3[k];

        t_8[k] = pa_y[k] * lss0_5[k]
                 - f_2 * pc_y[k] * lss1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, kss0_9, kss0_10, kss1_9, kss1_10, lss0_9, \
                         lss0_10, lss1_9, lss1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * kss0_9[k]
                 - f_6 * kss1_9[k]
                 + pa_x[k] * lss0_9[k]
                 - f_2 * pc_x[k] * lss1_9[k];

        t_10[k] = f_7 * kss0_10[k]
                  - f_8 * kss1_10[k]
                  + pa_x[k] * lss0_10[k]
                  - f_2 * pc_x[k] * lss1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, kss0_12, kss1_12, lss0_6, \
                         lss0_12, lss1_6, lss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * lss0_6[k]
                  - f_2 * pc_z[k] * lss1_6[k];

        t_12[k] = f_7 * kss0_12[k]
                  - f_8 * kss1_12[k]
                  + pa_x[k] * lss0_12[k]
                  - f_2 * pc_x[k] * lss1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, kss0_14, kss1_14, lss0_9, \
                         lss0_14, lss1_9, lss1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * lss0_9[k]
                  - f_2 * pc_y[k] * lss1_9[k];

        t_14[k] = f_7 * kss0_14[k]
                  - f_8 * kss1_14[k]
                  + pa_x[k] * lss0_14[k]
                  - f_2 * pc_x[k] * lss1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, kss0_15, kss1_15, lss0_10, \
                         lss0_15, lss1_10, lss1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * kss0_15[k]
                  - f_10 * kss1_15[k]
                  + pa_x[k] * lss0_15[k]
                  - f_2 * pc_x[k] * lss1_15[k];

        t_16[k] = pa_z[k] * lss0_10[k]
                  - f_2 * pc_z[k] * lss1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, kss0_17, kss0_18, kss1_17, kss1_18, lss0_17, \
                         lss0_18, lss1_17, lss1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * kss0_17[k]
                  - f_10 * kss1_17[k]
                  + pa_x[k] * lss0_17[k]
                  - f_2 * pc_x[k] * lss1_17[k];

        t_18[k] = f_9 * kss0_18[k]
                  - f_10 * kss1_18[k]
                  + pa_x[k] * lss0_18[k]
                  - f_2 * pc_x[k] * lss1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, pc_x, pc_y, kss0_20, kss1_20, lss0_14, \
                         lss0_20, lss1_14, lss1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * lss0_14[k]
                  - f_2 * pc_y[k] * lss1_14[k];

        t_20[k] = f_9 * kss0_20[k]
                  - f_10 * kss1_20[k]
                  + pa_x[k] * lss0_20[k]
                  - f_2 * pc_x[k] * lss1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_z, pc_x, pc_z, kss0_21, kss1_21, lss0_15, \
                         lss0_21, lss1_15, lss1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * kss0_21[k]
                  - f_12 * kss1_21[k]
                  + pa_x[k] * lss0_21[k]
                  - f_2 * pc_x[k] * lss1_21[k];

        t_22[k] = pa_z[k] * lss0_15[k]
                  - f_2 * pc_z[k] * lss1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pc_x, kss0_23, kss0_24, kss1_23, kss1_24, lss0_23, \
                         lss0_24, lss1_23, lss1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * kss0_23[k]
                  - f_12 * kss1_23[k]
                  + pa_x[k] * lss0_23[k]
                  - f_2 * pc_x[k] * lss1_23[k];

        t_24[k] = f_11 * kss0_24[k]
                  - f_12 * kss1_24[k]
                  + pa_x[k] * lss0_24[k]
                  - f_2 * pc_x[k] * lss1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pa_y, pc_x, pc_y, kss0_25, kss1_25, lss0_20, \
                         lss0_25, lss1_20, lss1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * kss0_25[k]
                  - f_12 * kss1_25[k]
                  + pa_x[k] * lss0_25[k]
                  - f_2 * pc_x[k] * lss1_25[k];

        t_26[k] = pa_y[k] * lss0_20[k]
                  - f_2 * pc_y[k] * lss1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pc_x, kss0_27, kss0_28, kss1_27, kss1_28, lss0_27, \
                         lss0_28, lss1_27, lss1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * kss0_27[k]
                  - f_12 * kss1_27[k]
                  + pa_x[k] * lss0_27[k]
                  - f_2 * pc_x[k] * lss1_27[k];

        t_28[k] = f_13 * kss0_28[k]
                  - f_14 * kss1_28[k]
                  + pa_x[k] * lss0_28[k]
                  - f_2 * pc_x[k] * lss1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_x, pa_z, pc_x, pc_z, kss0_30, kss1_30, lss0_21, \
                         lss0_30, lss1_21, lss1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * lss0_21[k]
                  - f_2 * pc_z[k] * lss1_21[k];

        t_30[k] = f_13 * kss0_30[k]
                  - f_14 * kss1_30[k]
                  + pa_x[k] * lss0_30[k]
                  - f_2 * pc_x[k] * lss1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pc_x, kss0_31, kss0_32, kss1_31, kss1_32, lss0_31, \
                         lss0_32, lss1_31, lss1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * kss0_31[k]
                  - f_14 * kss1_31[k]
                  + pa_x[k] * lss0_31[k]
                  - f_2 * pc_x[k] * lss1_31[k];

        t_32[k] = f_13 * kss0_32[k]
                  - f_14 * kss1_32[k]
                  + pa_x[k] * lss0_32[k]
                  - f_2 * pc_x[k] * lss1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pa_x, pa_y, pc_x, pc_y, kss0_33, kss1_33, lss0_27, \
                         lss0_33, lss1_27, lss1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * kss0_33[k]
                  - f_14 * kss1_33[k]
                  + pa_x[k] * lss0_33[k]
                  - f_2 * pc_x[k] * lss1_33[k];

        t_34[k] = pa_y[k] * lss0_27[k]
                  - f_2 * pc_y[k] * lss1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pc_x, kss0_35, kss1_35, lss0_35, \
                         lss0_36, lss0_37, lss0_38, lss1_35, lss1_36, lss1_37, \
                         lss1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * kss0_35[k]
                  - f_14 * kss1_35[k]
                  + pa_x[k] * lss0_35[k]
                  - f_2 * pc_x[k] * lss1_35[k];

        t_36[k] = pa_x[k] * lss0_36[k]
                  - f_2 * pc_x[k] * lss1_36[k];

        t_37[k] = pa_x[k] * lss0_37[k]
                  - f_2 * pc_x[k] * lss1_37[k];

        t_38[k] = pa_x[k] * lss0_38[k]
                  - f_2 * pc_x[k] * lss1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pc_x, lss0_39, lss0_40, lss0_41, \
                         lss0_42, lss1_39, lss1_40, lss1_41, lss1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pa_x[k] * lss0_39[k]
                  - f_2 * pc_x[k] * lss1_39[k];

        t_40[k] = pa_x[k] * lss0_40[k]
                  - f_2 * pc_x[k] * lss1_40[k];

        t_41[k] = pa_x[k] * lss0_41[k]
                  - f_2 * pc_x[k] * lss1_41[k];

        t_42[k] = pa_x[k] * lss0_42[k]
                  - f_2 * pc_x[k] * lss1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pc_x, pc_y, kss0_28, kss1_28, lss0_36, \
                         lss0_43, lss0_44, lss1_36, lss1_43, lss1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_x[k] * lss0_43[k]
                  - f_2 * pc_x[k] * lss1_43[k];

        t_44[k] = pa_x[k] * lss0_44[k]
                  - f_2 * pc_x[k] * lss1_44[k];

        t_45[k] = f_0 * kss0_28[k]
                  - f_1 * kss1_28[k]
                  + pa_y[k] * lss0_36[k]
                  - f_2 * pc_y[k] * lss1_36[k];
    }

#pragma omp simd aligned(t_46, t_47, pa_y, pa_z, pc_y, pc_z, kss0_30, kss1_30, lss0_36, \
                         lss0_38, lss1_36, lss1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * lss0_36[k]
                  - f_2 * pc_z[k] * lss1_36[k];

        t_47[k] = f_3 * kss0_30[k]
                  - f_4 * kss1_30[k]
                  + pa_y[k] * lss0_38[k]
                  - f_2 * pc_y[k] * lss1_38[k];
    }

#pragma omp simd aligned(t_48, t_49, pa_y, pc_y, kss0_31, kss0_32, kss1_31, kss1_32, lss0_39, \
                         lss0_40, lss1_39, lss1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * kss0_31[k]
                  - f_6 * kss1_31[k]
                  + pa_y[k] * lss0_39[k]
                  - f_2 * pc_y[k] * lss1_39[k];

        t_49[k] = f_7 * kss0_32[k]
                  - f_8 * kss1_32[k]
                  + pa_y[k] * lss0_40[k]
                  - f_2 * pc_y[k] * lss1_40[k];
    }

#pragma omp simd aligned(t_50, t_51, pa_y, pc_y, kss0_33, kss0_34, kss1_33, kss1_34, lss0_41, \
                         lss0_42, lss1_41, lss1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * kss0_33[k]
                  - f_10 * kss1_33[k]
                  + pa_y[k] * lss0_41[k]
                  - f_2 * pc_y[k] * lss1_41[k];

        t_51[k] = f_11 * kss0_34[k]
                  - f_12 * kss1_34[k]
                  + pa_y[k] * lss0_42[k]
                  - f_2 * pc_y[k] * lss1_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pa_y, pa_z, pc_y, pc_z, kss0_35, kss1_35, lss0_43, \
                         lss0_44, lss1_43, lss1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_13 * kss0_35[k]
                  - f_14 * kss1_35[k]
                  + pa_y[k] * lss0_43[k]
                  - f_2 * pc_y[k] * lss1_43[k];

        t_53[k] = pa_y[k] * lss0_44[k]
                  - f_2 * pc_y[k] * lss1_44[k];

        t_54[k] = f_0 * kss0_35[k]
                  - f_1 * kss1_35[k]
                  + pa_z[k] * lss0_44[k]
                  - f_2 * pc_z[k] * lss1_44[k];
    }
}

}  // namespace simdt3ceri
