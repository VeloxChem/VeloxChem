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


#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sns_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sls0, const size_t sls1,
                                                   const size_t sms0, const size_t sms1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 4.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 3.5 / p;
    const auto f_4 = 3.5 * gamma / (p * q);
    const auto f_5 = 3.0 / p;
    const auto f_6 = 3.0 * gamma / (p * q);
    const auto f_7 = 2.5 / p;
    const auto f_8 = 2.5 * gamma / (p * q);
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * gamma / (p * q);
    const auto f_11 = 1.5 / p;
    const auto f_12 = 1.5 * gamma / (p * q);
    const auto f_13 = 1.0 / p;
    const auto f_14 = gamma / (p * q);
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 * gamma / (p * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sls0_0 = buffer.data(sls0 + 0);
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
    const auto *sls0_38 = buffer.data(sls0 + 38);
    const auto *sls0_39 = buffer.data(sls0 + 39);
    const auto *sls0_40 = buffer.data(sls0 + 40);
    const auto *sls0_41 = buffer.data(sls0 + 41);
    const auto *sls0_42 = buffer.data(sls0 + 42);
    const auto *sls0_43 = buffer.data(sls0 + 43);
    const auto *sls0_44 = buffer.data(sls0 + 44);

    const auto *sls1_0 = buffer.data(sls1 + 0);
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
    const auto *sls1_38 = buffer.data(sls1 + 38);
    const auto *sls1_39 = buffer.data(sls1 + 39);
    const auto *sls1_40 = buffer.data(sls1 + 40);
    const auto *sls1_41 = buffer.data(sls1 + 41);
    const auto *sls1_42 = buffer.data(sls1 + 42);
    const auto *sls1_43 = buffer.data(sls1 + 43);
    const auto *sls1_44 = buffer.data(sls1 + 44);

    const auto *sms0_0 = buffer.data(sms0 + 0);
    const auto *sms0_2 = buffer.data(sms0 + 2);
    const auto *sms0_3 = buffer.data(sms0 + 3);
    const auto *sms0_5 = buffer.data(sms0 + 5);
    const auto *sms0_6 = buffer.data(sms0 + 6);
    const auto *sms0_9 = buffer.data(sms0 + 9);
    const auto *sms0_10 = buffer.data(sms0 + 10);
    const auto *sms0_12 = buffer.data(sms0 + 12);
    const auto *sms0_14 = buffer.data(sms0 + 14);
    const auto *sms0_15 = buffer.data(sms0 + 15);
    const auto *sms0_17 = buffer.data(sms0 + 17);
    const auto *sms0_18 = buffer.data(sms0 + 18);
    const auto *sms0_20 = buffer.data(sms0 + 20);
    const auto *sms0_21 = buffer.data(sms0 + 21);
    const auto *sms0_23 = buffer.data(sms0 + 23);
    const auto *sms0_24 = buffer.data(sms0 + 24);
    const auto *sms0_25 = buffer.data(sms0 + 25);
    const auto *sms0_27 = buffer.data(sms0 + 27);
    const auto *sms0_28 = buffer.data(sms0 + 28);
    const auto *sms0_30 = buffer.data(sms0 + 30);
    const auto *sms0_31 = buffer.data(sms0 + 31);
    const auto *sms0_32 = buffer.data(sms0 + 32);
    const auto *sms0_33 = buffer.data(sms0 + 33);
    const auto *sms0_35 = buffer.data(sms0 + 35);
    const auto *sms0_36 = buffer.data(sms0 + 36);
    const auto *sms0_38 = buffer.data(sms0 + 38);
    const auto *sms0_39 = buffer.data(sms0 + 39);
    const auto *sms0_40 = buffer.data(sms0 + 40);
    const auto *sms0_41 = buffer.data(sms0 + 41);
    const auto *sms0_42 = buffer.data(sms0 + 42);
    const auto *sms0_44 = buffer.data(sms0 + 44);
    const auto *sms0_45 = buffer.data(sms0 + 45);
    const auto *sms0_46 = buffer.data(sms0 + 46);
    const auto *sms0_47 = buffer.data(sms0 + 47);
    const auto *sms0_48 = buffer.data(sms0 + 48);
    const auto *sms0_49 = buffer.data(sms0 + 49);
    const auto *sms0_50 = buffer.data(sms0 + 50);
    const auto *sms0_51 = buffer.data(sms0 + 51);
    const auto *sms0_52 = buffer.data(sms0 + 52);
    const auto *sms0_53 = buffer.data(sms0 + 53);
    const auto *sms0_54 = buffer.data(sms0 + 54);

    const auto *sms1_0 = buffer.data(sms1 + 0);
    const auto *sms1_2 = buffer.data(sms1 + 2);
    const auto *sms1_3 = buffer.data(sms1 + 3);
    const auto *sms1_5 = buffer.data(sms1 + 5);
    const auto *sms1_6 = buffer.data(sms1 + 6);
    const auto *sms1_9 = buffer.data(sms1 + 9);
    const auto *sms1_10 = buffer.data(sms1 + 10);
    const auto *sms1_12 = buffer.data(sms1 + 12);
    const auto *sms1_14 = buffer.data(sms1 + 14);
    const auto *sms1_15 = buffer.data(sms1 + 15);
    const auto *sms1_17 = buffer.data(sms1 + 17);
    const auto *sms1_18 = buffer.data(sms1 + 18);
    const auto *sms1_20 = buffer.data(sms1 + 20);
    const auto *sms1_21 = buffer.data(sms1 + 21);
    const auto *sms1_23 = buffer.data(sms1 + 23);
    const auto *sms1_24 = buffer.data(sms1 + 24);
    const auto *sms1_25 = buffer.data(sms1 + 25);
    const auto *sms1_27 = buffer.data(sms1 + 27);
    const auto *sms1_28 = buffer.data(sms1 + 28);
    const auto *sms1_30 = buffer.data(sms1 + 30);
    const auto *sms1_31 = buffer.data(sms1 + 31);
    const auto *sms1_32 = buffer.data(sms1 + 32);
    const auto *sms1_33 = buffer.data(sms1 + 33);
    const auto *sms1_35 = buffer.data(sms1 + 35);
    const auto *sms1_36 = buffer.data(sms1 + 36);
    const auto *sms1_38 = buffer.data(sms1 + 38);
    const auto *sms1_39 = buffer.data(sms1 + 39);
    const auto *sms1_40 = buffer.data(sms1 + 40);
    const auto *sms1_41 = buffer.data(sms1 + 41);
    const auto *sms1_42 = buffer.data(sms1 + 42);
    const auto *sms1_44 = buffer.data(sms1 + 44);
    const auto *sms1_45 = buffer.data(sms1 + 45);
    const auto *sms1_46 = buffer.data(sms1 + 46);
    const auto *sms1_47 = buffer.data(sms1 + 47);
    const auto *sms1_48 = buffer.data(sms1 + 48);
    const auto *sms1_49 = buffer.data(sms1 + 49);
    const auto *sms1_50 = buffer.data(sms1 + 50);
    const auto *sms1_51 = buffer.data(sms1 + 51);
    const auto *sms1_52 = buffer.data(sms1 + 52);
    const auto *sms1_53 = buffer.data(sms1 + 53);
    const auto *sms1_54 = buffer.data(sms1 + 54);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sls0_0, sls1_0, \
                         sms0_0, sms1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sls0_0[k]
                 - f_1 * sls1_0[k]
                 + pb_x[k] * sms0_0[k]
                 - f_2 * pc_x[k] * sms1_0[k];

        t_1[k] = pb_y[k] * sms0_0[k]
                 - f_2 * pc_y[k] * sms1_0[k];

        t_2[k] = pb_z[k] * sms0_0[k]
                 - f_2 * pc_z[k] * sms1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sls0_3, sls1_3, sms0_2, sms0_3, \
                         sms1_2, sms1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sls0_3[k]
                 - f_4 * sls1_3[k]
                 + pb_x[k] * sms0_3[k]
                 - f_2 * pc_x[k] * sms1_3[k];

        t_4[k] = pb_y[k] * sms0_2[k]
                 - f_2 * pc_y[k] * sms1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sls0_5, sls0_6, sls1_5, sls1_6, sms0_5, sms0_6, \
                         sms1_5, sms1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sls0_5[k]
                 - f_4 * sls1_5[k]
                 + pb_x[k] * sms0_5[k]
                 - f_2 * pc_x[k] * sms1_5[k];

        t_6[k] = f_5 * sls0_6[k]
                 - f_6 * sls1_6[k]
                 + pb_x[k] * sms0_6[k]
                 - f_2 * pc_x[k] * sms1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sms0_3, sms0_5, sms1_3, \
                         sms1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sms0_3[k]
                 - f_2 * pc_z[k] * sms1_3[k];

        t_8[k] = pb_y[k] * sms0_5[k]
                 - f_2 * pc_y[k] * sms1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sls0_9, sls0_10, sls1_9, sls1_10, sms0_9, \
                         sms0_10, sms1_9, sms1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sls0_9[k]
                 - f_6 * sls1_9[k]
                 + pb_x[k] * sms0_9[k]
                 - f_2 * pc_x[k] * sms1_9[k];

        t_10[k] = f_7 * sls0_10[k]
                  - f_8 * sls1_10[k]
                  + pb_x[k] * sms0_10[k]
                  - f_2 * pc_x[k] * sms1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sls0_12, sls1_12, sms0_6, \
                         sms0_12, sms1_6, sms1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sms0_6[k]
                  - f_2 * pc_z[k] * sms1_6[k];

        t_12[k] = f_7 * sls0_12[k]
                  - f_8 * sls1_12[k]
                  + pb_x[k] * sms0_12[k]
                  - f_2 * pc_x[k] * sms1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sls0_14, sls1_14, sms0_9, \
                         sms0_14, sms1_9, sms1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sms0_9[k]
                  - f_2 * pc_y[k] * sms1_9[k];

        t_14[k] = f_7 * sls0_14[k]
                  - f_8 * sls1_14[k]
                  + pb_x[k] * sms0_14[k]
                  - f_2 * pc_x[k] * sms1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sls0_15, sls1_15, sms0_10, \
                         sms0_15, sms1_10, sms1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * sls0_15[k]
                  - f_10 * sls1_15[k]
                  + pb_x[k] * sms0_15[k]
                  - f_2 * pc_x[k] * sms1_15[k];

        t_16[k] = pb_z[k] * sms0_10[k]
                  - f_2 * pc_z[k] * sms1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sls0_17, sls0_18, sls1_17, sls1_18, sms0_17, \
                         sms0_18, sms1_17, sms1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * sls0_17[k]
                  - f_10 * sls1_17[k]
                  + pb_x[k] * sms0_17[k]
                  - f_2 * pc_x[k] * sms1_17[k];

        t_18[k] = f_9 * sls0_18[k]
                  - f_10 * sls1_18[k]
                  + pb_x[k] * sms0_18[k]
                  - f_2 * pc_x[k] * sms1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sls0_20, sls1_20, sms0_14, \
                         sms0_20, sms1_14, sms1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sms0_14[k]
                  - f_2 * pc_y[k] * sms1_14[k];

        t_20[k] = f_9 * sls0_20[k]
                  - f_10 * sls1_20[k]
                  + pb_x[k] * sms0_20[k]
                  - f_2 * pc_x[k] * sms1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sls0_21, sls1_21, sms0_15, \
                         sms0_21, sms1_15, sms1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * sls0_21[k]
                  - f_12 * sls1_21[k]
                  + pb_x[k] * sms0_21[k]
                  - f_2 * pc_x[k] * sms1_21[k];

        t_22[k] = pb_z[k] * sms0_15[k]
                  - f_2 * pc_z[k] * sms1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sls0_23, sls0_24, sls1_23, sls1_24, sms0_23, \
                         sms0_24, sms1_23, sms1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * sls0_23[k]
                  - f_12 * sls1_23[k]
                  + pb_x[k] * sms0_23[k]
                  - f_2 * pc_x[k] * sms1_23[k];

        t_24[k] = f_11 * sls0_24[k]
                  - f_12 * sls1_24[k]
                  + pb_x[k] * sms0_24[k]
                  - f_2 * pc_x[k] * sms1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sls0_25, sls1_25, sms0_20, \
                         sms0_25, sms1_20, sms1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * sls0_25[k]
                  - f_12 * sls1_25[k]
                  + pb_x[k] * sms0_25[k]
                  - f_2 * pc_x[k] * sms1_25[k];

        t_26[k] = pb_y[k] * sms0_20[k]
                  - f_2 * pc_y[k] * sms1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sls0_27, sls0_28, sls1_27, sls1_28, sms0_27, \
                         sms0_28, sms1_27, sms1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * sls0_27[k]
                  - f_12 * sls1_27[k]
                  + pb_x[k] * sms0_27[k]
                  - f_2 * pc_x[k] * sms1_27[k];

        t_28[k] = f_13 * sls0_28[k]
                  - f_14 * sls1_28[k]
                  + pb_x[k] * sms0_28[k]
                  - f_2 * pc_x[k] * sms1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sls0_30, sls1_30, sms0_21, \
                         sms0_30, sms1_21, sms1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * sms0_21[k]
                  - f_2 * pc_z[k] * sms1_21[k];

        t_30[k] = f_13 * sls0_30[k]
                  - f_14 * sls1_30[k]
                  + pb_x[k] * sms0_30[k]
                  - f_2 * pc_x[k] * sms1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sls0_31, sls0_32, sls1_31, sls1_32, sms0_31, \
                         sms0_32, sms1_31, sms1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * sls0_31[k]
                  - f_14 * sls1_31[k]
                  + pb_x[k] * sms0_31[k]
                  - f_2 * pc_x[k] * sms1_31[k];

        t_32[k] = f_13 * sls0_32[k]
                  - f_14 * sls1_32[k]
                  + pb_x[k] * sms0_32[k]
                  - f_2 * pc_x[k] * sms1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sls0_33, sls1_33, sms0_27, \
                         sms0_33, sms1_27, sms1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * sls0_33[k]
                  - f_14 * sls1_33[k]
                  + pb_x[k] * sms0_33[k]
                  - f_2 * pc_x[k] * sms1_33[k];

        t_34[k] = pb_y[k] * sms0_27[k]
                  - f_2 * pc_y[k] * sms1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_x, pc_x, sls0_35, sls0_36, sls1_35, sls1_36, sms0_35, \
                         sms0_36, sms1_35, sms1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * sls0_35[k]
                  - f_14 * sls1_35[k]
                  + pb_x[k] * sms0_35[k]
                  - f_2 * pc_x[k] * sms1_35[k];

        t_36[k] = f_15 * sls0_36[k]
                  - f_16 * sls1_36[k]
                  + pb_x[k] * sms0_36[k]
                  - f_2 * pc_x[k] * sms1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_x, pb_z, pc_x, pc_z, sls0_38, sls1_38, sms0_28, \
                         sms0_38, sms1_28, sms1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_z[k] * sms0_28[k]
                  - f_2 * pc_z[k] * sms1_28[k];

        t_38[k] = f_15 * sls0_38[k]
                  - f_16 * sls1_38[k]
                  + pb_x[k] * sms0_38[k]
                  - f_2 * pc_x[k] * sms1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pb_x, pc_x, sls0_39, sls0_40, sls1_39, sls1_40, sms0_39, \
                         sms0_40, sms1_39, sms1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_15 * sls0_39[k]
                  - f_16 * sls1_39[k]
                  + pb_x[k] * sms0_39[k]
                  - f_2 * pc_x[k] * sms1_39[k];

        t_40[k] = f_15 * sls0_40[k]
                  - f_16 * sls1_40[k]
                  + pb_x[k] * sms0_40[k]
                  - f_2 * pc_x[k] * sms1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pb_x, pc_x, sls0_41, sls0_42, sls1_41, sls1_42, sms0_41, \
                         sms0_42, sms1_41, sms1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_15 * sls0_41[k]
                  - f_16 * sls1_41[k]
                  + pb_x[k] * sms0_41[k]
                  - f_2 * pc_x[k] * sms1_41[k];

        t_42[k] = f_15 * sls0_42[k]
                  - f_16 * sls1_42[k]
                  + pb_x[k] * sms0_42[k]
                  - f_2 * pc_x[k] * sms1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pb_x, pb_y, pc_x, pc_y, sls0_44, sls1_44, sms0_35, \
                         sms0_44, sms0_45, sms1_35, sms1_44, sms1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * sms0_35[k]
                  - f_2 * pc_y[k] * sms1_35[k];

        t_44[k] = f_15 * sls0_44[k]
                  - f_16 * sls1_44[k]
                  + pb_x[k] * sms0_44[k]
                  - f_2 * pc_x[k] * sms1_44[k];

        t_45[k] = pb_x[k] * sms0_45[k]
                  - f_2 * pc_x[k] * sms1_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_x, pc_x, sms0_46, sms0_47, sms0_48, \
                         sms0_49, sms1_46, sms1_47, sms1_48, sms1_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_x[k] * sms0_46[k]
                  - f_2 * pc_x[k] * sms1_46[k];

        t_47[k] = pb_x[k] * sms0_47[k]
                  - f_2 * pc_x[k] * sms1_47[k];

        t_48[k] = pb_x[k] * sms0_48[k]
                  - f_2 * pc_x[k] * sms1_48[k];

        t_49[k] = pb_x[k] * sms0_49[k]
                  - f_2 * pc_x[k] * sms1_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pc_x, sms0_50, sms0_51, sms0_52, \
                         sms0_53, sms1_50, sms1_51, sms1_52, sms1_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_x[k] * sms0_50[k]
                  - f_2 * pc_x[k] * sms1_50[k];

        t_51[k] = pb_x[k] * sms0_51[k]
                  - f_2 * pc_x[k] * sms1_51[k];

        t_52[k] = pb_x[k] * sms0_52[k]
                  - f_2 * pc_x[k] * sms1_52[k];

        t_53[k] = pb_x[k] * sms0_53[k]
                  - f_2 * pc_x[k] * sms1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sls0_36, \
                         sls1_36, sms0_45, sms0_54, sms1_45, sms1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pb_x[k] * sms0_54[k]
                  - f_2 * pc_x[k] * sms1_54[k];

        t_55[k] = f_0 * sls0_36[k]
                  - f_1 * sls1_36[k]
                  + pb_y[k] * sms0_45[k]
                  - f_2 * pc_y[k] * sms1_45[k];

        t_56[k] = pb_z[k] * sms0_45[k]
                  - f_2 * pc_z[k] * sms1_45[k];
    }

#pragma omp simd aligned(t_57, t_58, pb_y, pc_y, sls0_38, sls0_39, sls1_38, sls1_39, sms0_47, \
                         sms0_48, sms1_47, sms1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_3 * sls0_38[k]
                  - f_4 * sls1_38[k]
                  + pb_y[k] * sms0_47[k]
                  - f_2 * pc_y[k] * sms1_47[k];

        t_58[k] = f_5 * sls0_39[k]
                  - f_6 * sls1_39[k]
                  + pb_y[k] * sms0_48[k]
                  - f_2 * pc_y[k] * sms1_48[k];
    }

#pragma omp simd aligned(t_59, t_60, pb_y, pc_y, sls0_40, sls0_41, sls1_40, sls1_41, sms0_49, \
                         sms0_50, sms1_49, sms1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * sls0_40[k]
                  - f_8 * sls1_40[k]
                  + pb_y[k] * sms0_49[k]
                  - f_2 * pc_y[k] * sms1_49[k];

        t_60[k] = f_9 * sls0_41[k]
                  - f_10 * sls1_41[k]
                  + pb_y[k] * sms0_50[k]
                  - f_2 * pc_y[k] * sms1_50[k];
    }

#pragma omp simd aligned(t_61, t_62, pb_y, pc_y, sls0_42, sls0_43, sls1_42, sls1_43, sms0_51, \
                         sms0_52, sms1_51, sms1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * sls0_42[k]
                  - f_12 * sls1_42[k]
                  + pb_y[k] * sms0_51[k]
                  - f_2 * pc_y[k] * sms1_51[k];

        t_62[k] = f_13 * sls0_43[k]
                  - f_14 * sls1_43[k]
                  + pb_y[k] * sms0_52[k]
                  - f_2 * pc_y[k] * sms1_52[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pb_y, pb_z, pc_y, pc_z, sls0_44, sls1_44, sms0_53, \
                         sms0_54, sms1_53, sms1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_15 * sls0_44[k]
                  - f_16 * sls1_44[k]
                  + pb_y[k] * sms0_53[k]
                  - f_2 * pc_y[k] * sms1_53[k];

        t_64[k] = pb_y[k] * sms0_54[k]
                  - f_2 * pc_y[k] * sms1_54[k];

        t_65[k] = f_0 * sls0_44[k]
                  - f_1 * sls1_44[k]
                  + pb_z[k] * sms0_54[k]
                  - f_2 * pc_z[k] * sms1_54[k];
    }
}

}  // namespace simdt3ceri
