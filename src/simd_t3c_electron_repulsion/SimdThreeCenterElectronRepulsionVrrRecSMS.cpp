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


#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sms_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sks0, const size_t sks1,
                                                   const size_t sls0, const size_t sls1,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sks0_0 = buffer.data(sks0 + 0);
    const auto *sks0_3 = buffer.data(sks0 + 3);
    const auto *sks0_5 = buffer.data(sks0 + 5);
    const auto *sks0_6 = buffer.data(sks0 + 6);
    const auto *sks0_9 = buffer.data(sks0 + 9);
    const auto *sks0_10 = buffer.data(sks0 + 10);
    const auto *sks0_12 = buffer.data(sks0 + 12);
    const auto *sks0_14 = buffer.data(sks0 + 14);
    const auto *sks0_15 = buffer.data(sks0 + 15);
    const auto *sks0_17 = buffer.data(sks0 + 17);
    const auto *sks0_18 = buffer.data(sks0 + 18);
    const auto *sks0_20 = buffer.data(sks0 + 20);
    const auto *sks0_21 = buffer.data(sks0 + 21);
    const auto *sks0_23 = buffer.data(sks0 + 23);
    const auto *sks0_24 = buffer.data(sks0 + 24);
    const auto *sks0_25 = buffer.data(sks0 + 25);
    const auto *sks0_27 = buffer.data(sks0 + 27);
    const auto *sks0_28 = buffer.data(sks0 + 28);
    const auto *sks0_30 = buffer.data(sks0 + 30);
    const auto *sks0_31 = buffer.data(sks0 + 31);
    const auto *sks0_32 = buffer.data(sks0 + 32);
    const auto *sks0_33 = buffer.data(sks0 + 33);
    const auto *sks0_34 = buffer.data(sks0 + 34);
    const auto *sks0_35 = buffer.data(sks0 + 35);

    const auto *sks1_0 = buffer.data(sks1 + 0);
    const auto *sks1_3 = buffer.data(sks1 + 3);
    const auto *sks1_5 = buffer.data(sks1 + 5);
    const auto *sks1_6 = buffer.data(sks1 + 6);
    const auto *sks1_9 = buffer.data(sks1 + 9);
    const auto *sks1_10 = buffer.data(sks1 + 10);
    const auto *sks1_12 = buffer.data(sks1 + 12);
    const auto *sks1_14 = buffer.data(sks1 + 14);
    const auto *sks1_15 = buffer.data(sks1 + 15);
    const auto *sks1_17 = buffer.data(sks1 + 17);
    const auto *sks1_18 = buffer.data(sks1 + 18);
    const auto *sks1_20 = buffer.data(sks1 + 20);
    const auto *sks1_21 = buffer.data(sks1 + 21);
    const auto *sks1_23 = buffer.data(sks1 + 23);
    const auto *sks1_24 = buffer.data(sks1 + 24);
    const auto *sks1_25 = buffer.data(sks1 + 25);
    const auto *sks1_27 = buffer.data(sks1 + 27);
    const auto *sks1_28 = buffer.data(sks1 + 28);
    const auto *sks1_30 = buffer.data(sks1 + 30);
    const auto *sks1_31 = buffer.data(sks1 + 31);
    const auto *sks1_32 = buffer.data(sks1 + 32);
    const auto *sks1_33 = buffer.data(sks1 + 33);
    const auto *sks1_34 = buffer.data(sks1 + 34);
    const auto *sks1_35 = buffer.data(sks1 + 35);

    const auto *sls0_0 = buffer.data(sls0 + 0);
    const auto *sls0_2 = buffer.data(sls0 + 2);
    const auto *sls0_3 = buffer.data(sls0 + 3);
    const auto *sls0_5 = buffer.data(sls0 + 5);
    const auto *sls0_6 = buffer.data(sls0 + 6);
    const auto *sls0_9 = buffer.data(sls0 + 9);
    const auto *sls0_10 = buffer.data(sls0 + 10);
    const auto *sls0_12 = buffer.data(sls0 + 12);
    const auto *sls0_14 = buffer.data(sls0 + 14);
    const auto *sls0_15 = buffer.data(sls0 + 15);
    const auto *sls0_17 = buffer.data(sls0 + 17);
    const auto *sls0_18 = buffer.data(sls0 + 18);
    const auto *sls0_20 = buffer.data(sls0 + 20);
    const auto *sls0_21 = buffer.data(sls0 + 21);
    const auto *sls0_23 = buffer.data(sls0 + 23);
    const auto *sls0_24 = buffer.data(sls0 + 24);
    const auto *sls0_25 = buffer.data(sls0 + 25);
    const auto *sls0_27 = buffer.data(sls0 + 27);
    const auto *sls0_28 = buffer.data(sls0 + 28);
    const auto *sls0_30 = buffer.data(sls0 + 30);
    const auto *sls0_31 = buffer.data(sls0 + 31);
    const auto *sls0_32 = buffer.data(sls0 + 32);
    const auto *sls0_33 = buffer.data(sls0 + 33);
    const auto *sls0_35 = buffer.data(sls0 + 35);
    const auto *sls0_36 = buffer.data(sls0 + 36);
    const auto *sls0_37 = buffer.data(sls0 + 37);
    const auto *sls0_38 = buffer.data(sls0 + 38);
    const auto *sls0_39 = buffer.data(sls0 + 39);
    const auto *sls0_40 = buffer.data(sls0 + 40);
    const auto *sls0_41 = buffer.data(sls0 + 41);
    const auto *sls0_42 = buffer.data(sls0 + 42);
    const auto *sls0_43 = buffer.data(sls0 + 43);
    const auto *sls0_44 = buffer.data(sls0 + 44);

    const auto *sls1_0 = buffer.data(sls1 + 0);
    const auto *sls1_2 = buffer.data(sls1 + 2);
    const auto *sls1_3 = buffer.data(sls1 + 3);
    const auto *sls1_5 = buffer.data(sls1 + 5);
    const auto *sls1_6 = buffer.data(sls1 + 6);
    const auto *sls1_9 = buffer.data(sls1 + 9);
    const auto *sls1_10 = buffer.data(sls1 + 10);
    const auto *sls1_12 = buffer.data(sls1 + 12);
    const auto *sls1_14 = buffer.data(sls1 + 14);
    const auto *sls1_15 = buffer.data(sls1 + 15);
    const auto *sls1_17 = buffer.data(sls1 + 17);
    const auto *sls1_18 = buffer.data(sls1 + 18);
    const auto *sls1_20 = buffer.data(sls1 + 20);
    const auto *sls1_21 = buffer.data(sls1 + 21);
    const auto *sls1_23 = buffer.data(sls1 + 23);
    const auto *sls1_24 = buffer.data(sls1 + 24);
    const auto *sls1_25 = buffer.data(sls1 + 25);
    const auto *sls1_27 = buffer.data(sls1 + 27);
    const auto *sls1_28 = buffer.data(sls1 + 28);
    const auto *sls1_30 = buffer.data(sls1 + 30);
    const auto *sls1_31 = buffer.data(sls1 + 31);
    const auto *sls1_32 = buffer.data(sls1 + 32);
    const auto *sls1_33 = buffer.data(sls1 + 33);
    const auto *sls1_35 = buffer.data(sls1 + 35);
    const auto *sls1_36 = buffer.data(sls1 + 36);
    const auto *sls1_37 = buffer.data(sls1 + 37);
    const auto *sls1_38 = buffer.data(sls1 + 38);
    const auto *sls1_39 = buffer.data(sls1 + 39);
    const auto *sls1_40 = buffer.data(sls1 + 40);
    const auto *sls1_41 = buffer.data(sls1 + 41);
    const auto *sls1_42 = buffer.data(sls1 + 42);
    const auto *sls1_43 = buffer.data(sls1 + 43);
    const auto *sls1_44 = buffer.data(sls1 + 44);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sks0_0, sks1_0, \
                         sls0_0, sls1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sks0_0[k]
                 - f_1 * sks1_0[k]
                 + pb_x[k] * sls0_0[k]
                 - f_2 * pc_x[k] * sls1_0[k];

        t_1[k] = pb_y[k] * sls0_0[k]
                 - f_2 * pc_y[k] * sls1_0[k];

        t_2[k] = pb_z[k] * sls0_0[k]
                 - f_2 * pc_z[k] * sls1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sks0_3, sks1_3, sls0_2, sls0_3, \
                         sls1_2, sls1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sks0_3[k]
                 - f_4 * sks1_3[k]
                 + pb_x[k] * sls0_3[k]
                 - f_2 * pc_x[k] * sls1_3[k];

        t_4[k] = pb_y[k] * sls0_2[k]
                 - f_2 * pc_y[k] * sls1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sks0_5, sks0_6, sks1_5, sks1_6, sls0_5, sls0_6, \
                         sls1_5, sls1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sks0_5[k]
                 - f_4 * sks1_5[k]
                 + pb_x[k] * sls0_5[k]
                 - f_2 * pc_x[k] * sls1_5[k];

        t_6[k] = f_5 * sks0_6[k]
                 - f_6 * sks1_6[k]
                 + pb_x[k] * sls0_6[k]
                 - f_2 * pc_x[k] * sls1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sls0_3, sls0_5, sls1_3, \
                         sls1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sls0_3[k]
                 - f_2 * pc_z[k] * sls1_3[k];

        t_8[k] = pb_y[k] * sls0_5[k]
                 - f_2 * pc_y[k] * sls1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sks0_9, sks0_10, sks1_9, sks1_10, sls0_9, \
                         sls0_10, sls1_9, sls1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sks0_9[k]
                 - f_6 * sks1_9[k]
                 + pb_x[k] * sls0_9[k]
                 - f_2 * pc_x[k] * sls1_9[k];

        t_10[k] = f_7 * sks0_10[k]
                  - f_8 * sks1_10[k]
                  + pb_x[k] * sls0_10[k]
                  - f_2 * pc_x[k] * sls1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sks0_12, sks1_12, sls0_6, \
                         sls0_12, sls1_6, sls1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sls0_6[k]
                  - f_2 * pc_z[k] * sls1_6[k];

        t_12[k] = f_7 * sks0_12[k]
                  - f_8 * sks1_12[k]
                  + pb_x[k] * sls0_12[k]
                  - f_2 * pc_x[k] * sls1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sks0_14, sks1_14, sls0_9, \
                         sls0_14, sls1_9, sls1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sls0_9[k]
                  - f_2 * pc_y[k] * sls1_9[k];

        t_14[k] = f_7 * sks0_14[k]
                  - f_8 * sks1_14[k]
                  + pb_x[k] * sls0_14[k]
                  - f_2 * pc_x[k] * sls1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sks0_15, sks1_15, sls0_10, \
                         sls0_15, sls1_10, sls1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * sks0_15[k]
                  - f_10 * sks1_15[k]
                  + pb_x[k] * sls0_15[k]
                  - f_2 * pc_x[k] * sls1_15[k];

        t_16[k] = pb_z[k] * sls0_10[k]
                  - f_2 * pc_z[k] * sls1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sks0_17, sks0_18, sks1_17, sks1_18, sls0_17, \
                         sls0_18, sls1_17, sls1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * sks0_17[k]
                  - f_10 * sks1_17[k]
                  + pb_x[k] * sls0_17[k]
                  - f_2 * pc_x[k] * sls1_17[k];

        t_18[k] = f_9 * sks0_18[k]
                  - f_10 * sks1_18[k]
                  + pb_x[k] * sls0_18[k]
                  - f_2 * pc_x[k] * sls1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sks0_20, sks1_20, sls0_14, \
                         sls0_20, sls1_14, sls1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sls0_14[k]
                  - f_2 * pc_y[k] * sls1_14[k];

        t_20[k] = f_9 * sks0_20[k]
                  - f_10 * sks1_20[k]
                  + pb_x[k] * sls0_20[k]
                  - f_2 * pc_x[k] * sls1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sks0_21, sks1_21, sls0_15, \
                         sls0_21, sls1_15, sls1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * sks0_21[k]
                  - f_12 * sks1_21[k]
                  + pb_x[k] * sls0_21[k]
                  - f_2 * pc_x[k] * sls1_21[k];

        t_22[k] = pb_z[k] * sls0_15[k]
                  - f_2 * pc_z[k] * sls1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sks0_23, sks0_24, sks1_23, sks1_24, sls0_23, \
                         sls0_24, sls1_23, sls1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * sks0_23[k]
                  - f_12 * sks1_23[k]
                  + pb_x[k] * sls0_23[k]
                  - f_2 * pc_x[k] * sls1_23[k];

        t_24[k] = f_11 * sks0_24[k]
                  - f_12 * sks1_24[k]
                  + pb_x[k] * sls0_24[k]
                  - f_2 * pc_x[k] * sls1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sks0_25, sks1_25, sls0_20, \
                         sls0_25, sls1_20, sls1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * sks0_25[k]
                  - f_12 * sks1_25[k]
                  + pb_x[k] * sls0_25[k]
                  - f_2 * pc_x[k] * sls1_25[k];

        t_26[k] = pb_y[k] * sls0_20[k]
                  - f_2 * pc_y[k] * sls1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sks0_27, sks0_28, sks1_27, sks1_28, sls0_27, \
                         sls0_28, sls1_27, sls1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * sks0_27[k]
                  - f_12 * sks1_27[k]
                  + pb_x[k] * sls0_27[k]
                  - f_2 * pc_x[k] * sls1_27[k];

        t_28[k] = f_13 * sks0_28[k]
                  - f_14 * sks1_28[k]
                  + pb_x[k] * sls0_28[k]
                  - f_2 * pc_x[k] * sls1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sks0_30, sks1_30, sls0_21, \
                         sls0_30, sls1_21, sls1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * sls0_21[k]
                  - f_2 * pc_z[k] * sls1_21[k];

        t_30[k] = f_13 * sks0_30[k]
                  - f_14 * sks1_30[k]
                  + pb_x[k] * sls0_30[k]
                  - f_2 * pc_x[k] * sls1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sks0_31, sks0_32, sks1_31, sks1_32, sls0_31, \
                         sls0_32, sls1_31, sls1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * sks0_31[k]
                  - f_14 * sks1_31[k]
                  + pb_x[k] * sls0_31[k]
                  - f_2 * pc_x[k] * sls1_31[k];

        t_32[k] = f_13 * sks0_32[k]
                  - f_14 * sks1_32[k]
                  + pb_x[k] * sls0_32[k]
                  - f_2 * pc_x[k] * sls1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sks0_33, sks1_33, sls0_27, \
                         sls0_33, sls1_27, sls1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * sks0_33[k]
                  - f_14 * sks1_33[k]
                  + pb_x[k] * sls0_33[k]
                  - f_2 * pc_x[k] * sls1_33[k];

        t_34[k] = pb_y[k] * sls0_27[k]
                  - f_2 * pc_y[k] * sls1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pc_x, sks0_35, sks1_35, sls0_35, \
                         sls0_36, sls0_37, sls0_38, sls1_35, sls1_36, sls1_37, \
                         sls1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * sks0_35[k]
                  - f_14 * sks1_35[k]
                  + pb_x[k] * sls0_35[k]
                  - f_2 * pc_x[k] * sls1_35[k];

        t_36[k] = pb_x[k] * sls0_36[k]
                  - f_2 * pc_x[k] * sls1_36[k];

        t_37[k] = pb_x[k] * sls0_37[k]
                  - f_2 * pc_x[k] * sls1_37[k];

        t_38[k] = pb_x[k] * sls0_38[k]
                  - f_2 * pc_x[k] * sls1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pb_x, pc_x, sls0_39, sls0_40, sls0_41, \
                         sls0_42, sls1_39, sls1_40, sls1_41, sls1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_x[k] * sls0_39[k]
                  - f_2 * pc_x[k] * sls1_39[k];

        t_40[k] = pb_x[k] * sls0_40[k]
                  - f_2 * pc_x[k] * sls1_40[k];

        t_41[k] = pb_x[k] * sls0_41[k]
                  - f_2 * pc_x[k] * sls1_41[k];

        t_42[k] = pb_x[k] * sls0_42[k]
                  - f_2 * pc_x[k] * sls1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, pc_x, pc_y, sks0_28, sks1_28, sls0_36, \
                         sls0_43, sls0_44, sls1_36, sls1_43, sls1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_x[k] * sls0_43[k]
                  - f_2 * pc_x[k] * sls1_43[k];

        t_44[k] = pb_x[k] * sls0_44[k]
                  - f_2 * pc_x[k] * sls1_44[k];

        t_45[k] = f_0 * sks0_28[k]
                  - f_1 * sks1_28[k]
                  + pb_y[k] * sls0_36[k]
                  - f_2 * pc_y[k] * sls1_36[k];
    }

#pragma omp simd aligned(t_46, t_47, pb_y, pb_z, pc_y, pc_z, sks0_30, sks1_30, sls0_36, \
                         sls0_38, sls1_36, sls1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * sls0_36[k]
                  - f_2 * pc_z[k] * sls1_36[k];

        t_47[k] = f_3 * sks0_30[k]
                  - f_4 * sks1_30[k]
                  + pb_y[k] * sls0_38[k]
                  - f_2 * pc_y[k] * sls1_38[k];
    }

#pragma omp simd aligned(t_48, t_49, pb_y, pc_y, sks0_31, sks0_32, sks1_31, sks1_32, sls0_39, \
                         sls0_40, sls1_39, sls1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * sks0_31[k]
                  - f_6 * sks1_31[k]
                  + pb_y[k] * sls0_39[k]
                  - f_2 * pc_y[k] * sls1_39[k];

        t_49[k] = f_7 * sks0_32[k]
                  - f_8 * sks1_32[k]
                  + pb_y[k] * sls0_40[k]
                  - f_2 * pc_y[k] * sls1_40[k];
    }

#pragma omp simd aligned(t_50, t_51, pb_y, pc_y, sks0_33, sks0_34, sks1_33, sks1_34, sls0_41, \
                         sls0_42, sls1_41, sls1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * sks0_33[k]
                  - f_10 * sks1_33[k]
                  + pb_y[k] * sls0_41[k]
                  - f_2 * pc_y[k] * sls1_41[k];

        t_51[k] = f_11 * sks0_34[k]
                  - f_12 * sks1_34[k]
                  + pb_y[k] * sls0_42[k]
                  - f_2 * pc_y[k] * sls1_42[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_y, pb_z, pc_y, pc_z, sks0_35, sks1_35, sls0_43, \
                         sls0_44, sls1_43, sls1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_13 * sks0_35[k]
                  - f_14 * sks1_35[k]
                  + pb_y[k] * sls0_43[k]
                  - f_2 * pc_y[k] * sls1_43[k];

        t_53[k] = pb_y[k] * sls0_44[k]
                  - f_2 * pc_y[k] * sls1_44[k];

        t_54[k] = f_0 * sks0_35[k]
                  - f_1 * sks1_35[k]
                  + pb_z[k] * sls0_44[k]
                  - f_2 * pc_z[k] * sls1_44[k];
    }
}

}  // namespace simdt3ceri
