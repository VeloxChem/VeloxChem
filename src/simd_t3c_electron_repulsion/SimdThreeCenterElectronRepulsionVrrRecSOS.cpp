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


#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sos_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sms0, const size_t sms1,
                                                   const size_t sns0, const size_t sns1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / p;
    const auto f_1 = 5.0 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 4.0 / p;
    const auto f_4 = 4.0 * gamma / (p * q);
    const auto f_5 = 3.5 / p;
    const auto f_6 = 3.5 * gamma / (p * q);
    const auto f_7 = 3.0 / p;
    const auto f_8 = 3.0 * gamma / (p * q);
    const auto f_9 = 2.5 / p;
    const auto f_10 = 2.5 * gamma / (p * q);
    const auto f_11 = 2.0 / p;
    const auto f_12 = 2.0 * gamma / (p * q);
    const auto f_13 = 1.5 / p;
    const auto f_14 = 1.5 * gamma / (p * q);
    const auto f_15 = 1.0 / p;
    const auto f_16 = gamma / (p * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sms0_0 = buffer.data(sms0 + 0);
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
    const auto *sms0_47 = buffer.data(sms0 + 47);
    const auto *sms0_48 = buffer.data(sms0 + 48);
    const auto *sms0_49 = buffer.data(sms0 + 49);
    const auto *sms0_50 = buffer.data(sms0 + 50);
    const auto *sms0_51 = buffer.data(sms0 + 51);
    const auto *sms0_52 = buffer.data(sms0 + 52);
    const auto *sms0_53 = buffer.data(sms0 + 53);
    const auto *sms0_54 = buffer.data(sms0 + 54);

    const auto *sms1_0 = buffer.data(sms1 + 0);
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
    const auto *sms1_47 = buffer.data(sms1 + 47);
    const auto *sms1_48 = buffer.data(sms1 + 48);
    const auto *sms1_49 = buffer.data(sms1 + 49);
    const auto *sms1_50 = buffer.data(sms1 + 50);
    const auto *sms1_51 = buffer.data(sms1 + 51);
    const auto *sms1_52 = buffer.data(sms1 + 52);
    const auto *sms1_53 = buffer.data(sms1 + 53);
    const auto *sms1_54 = buffer.data(sms1 + 54);

    const auto *sns0_0 = buffer.data(sns0 + 0);
    const auto *sns0_2 = buffer.data(sns0 + 2);
    const auto *sns0_3 = buffer.data(sns0 + 3);
    const auto *sns0_5 = buffer.data(sns0 + 5);
    const auto *sns0_6 = buffer.data(sns0 + 6);
    const auto *sns0_9 = buffer.data(sns0 + 9);
    const auto *sns0_10 = buffer.data(sns0 + 10);
    const auto *sns0_12 = buffer.data(sns0 + 12);
    const auto *sns0_14 = buffer.data(sns0 + 14);
    const auto *sns0_15 = buffer.data(sns0 + 15);
    const auto *sns0_17 = buffer.data(sns0 + 17);
    const auto *sns0_18 = buffer.data(sns0 + 18);
    const auto *sns0_20 = buffer.data(sns0 + 20);
    const auto *sns0_21 = buffer.data(sns0 + 21);
    const auto *sns0_23 = buffer.data(sns0 + 23);
    const auto *sns0_24 = buffer.data(sns0 + 24);
    const auto *sns0_25 = buffer.data(sns0 + 25);
    const auto *sns0_27 = buffer.data(sns0 + 27);
    const auto *sns0_28 = buffer.data(sns0 + 28);
    const auto *sns0_30 = buffer.data(sns0 + 30);
    const auto *sns0_31 = buffer.data(sns0 + 31);
    const auto *sns0_32 = buffer.data(sns0 + 32);
    const auto *sns0_33 = buffer.data(sns0 + 33);
    const auto *sns0_35 = buffer.data(sns0 + 35);
    const auto *sns0_36 = buffer.data(sns0 + 36);
    const auto *sns0_38 = buffer.data(sns0 + 38);
    const auto *sns0_39 = buffer.data(sns0 + 39);
    const auto *sns0_40 = buffer.data(sns0 + 40);
    const auto *sns0_41 = buffer.data(sns0 + 41);
    const auto *sns0_42 = buffer.data(sns0 + 42);
    const auto *sns0_44 = buffer.data(sns0 + 44);
    const auto *sns0_45 = buffer.data(sns0 + 45);
    const auto *sns0_47 = buffer.data(sns0 + 47);
    const auto *sns0_48 = buffer.data(sns0 + 48);
    const auto *sns0_49 = buffer.data(sns0 + 49);
    const auto *sns0_50 = buffer.data(sns0 + 50);
    const auto *sns0_51 = buffer.data(sns0 + 51);
    const auto *sns0_52 = buffer.data(sns0 + 52);
    const auto *sns0_54 = buffer.data(sns0 + 54);
    const auto *sns0_55 = buffer.data(sns0 + 55);
    const auto *sns0_56 = buffer.data(sns0 + 56);
    const auto *sns0_57 = buffer.data(sns0 + 57);
    const auto *sns0_58 = buffer.data(sns0 + 58);
    const auto *sns0_59 = buffer.data(sns0 + 59);
    const auto *sns0_60 = buffer.data(sns0 + 60);
    const auto *sns0_61 = buffer.data(sns0 + 61);
    const auto *sns0_62 = buffer.data(sns0 + 62);
    const auto *sns0_63 = buffer.data(sns0 + 63);
    const auto *sns0_64 = buffer.data(sns0 + 64);
    const auto *sns0_65 = buffer.data(sns0 + 65);

    const auto *sns1_0 = buffer.data(sns1 + 0);
    const auto *sns1_2 = buffer.data(sns1 + 2);
    const auto *sns1_3 = buffer.data(sns1 + 3);
    const auto *sns1_5 = buffer.data(sns1 + 5);
    const auto *sns1_6 = buffer.data(sns1 + 6);
    const auto *sns1_9 = buffer.data(sns1 + 9);
    const auto *sns1_10 = buffer.data(sns1 + 10);
    const auto *sns1_12 = buffer.data(sns1 + 12);
    const auto *sns1_14 = buffer.data(sns1 + 14);
    const auto *sns1_15 = buffer.data(sns1 + 15);
    const auto *sns1_17 = buffer.data(sns1 + 17);
    const auto *sns1_18 = buffer.data(sns1 + 18);
    const auto *sns1_20 = buffer.data(sns1 + 20);
    const auto *sns1_21 = buffer.data(sns1 + 21);
    const auto *sns1_23 = buffer.data(sns1 + 23);
    const auto *sns1_24 = buffer.data(sns1 + 24);
    const auto *sns1_25 = buffer.data(sns1 + 25);
    const auto *sns1_27 = buffer.data(sns1 + 27);
    const auto *sns1_28 = buffer.data(sns1 + 28);
    const auto *sns1_30 = buffer.data(sns1 + 30);
    const auto *sns1_31 = buffer.data(sns1 + 31);
    const auto *sns1_32 = buffer.data(sns1 + 32);
    const auto *sns1_33 = buffer.data(sns1 + 33);
    const auto *sns1_35 = buffer.data(sns1 + 35);
    const auto *sns1_36 = buffer.data(sns1 + 36);
    const auto *sns1_38 = buffer.data(sns1 + 38);
    const auto *sns1_39 = buffer.data(sns1 + 39);
    const auto *sns1_40 = buffer.data(sns1 + 40);
    const auto *sns1_41 = buffer.data(sns1 + 41);
    const auto *sns1_42 = buffer.data(sns1 + 42);
    const auto *sns1_44 = buffer.data(sns1 + 44);
    const auto *sns1_45 = buffer.data(sns1 + 45);
    const auto *sns1_47 = buffer.data(sns1 + 47);
    const auto *sns1_48 = buffer.data(sns1 + 48);
    const auto *sns1_49 = buffer.data(sns1 + 49);
    const auto *sns1_50 = buffer.data(sns1 + 50);
    const auto *sns1_51 = buffer.data(sns1 + 51);
    const auto *sns1_52 = buffer.data(sns1 + 52);
    const auto *sns1_54 = buffer.data(sns1 + 54);
    const auto *sns1_55 = buffer.data(sns1 + 55);
    const auto *sns1_56 = buffer.data(sns1 + 56);
    const auto *sns1_57 = buffer.data(sns1 + 57);
    const auto *sns1_58 = buffer.data(sns1 + 58);
    const auto *sns1_59 = buffer.data(sns1 + 59);
    const auto *sns1_60 = buffer.data(sns1 + 60);
    const auto *sns1_61 = buffer.data(sns1 + 61);
    const auto *sns1_62 = buffer.data(sns1 + 62);
    const auto *sns1_63 = buffer.data(sns1 + 63);
    const auto *sns1_64 = buffer.data(sns1 + 64);
    const auto *sns1_65 = buffer.data(sns1 + 65);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, sms0_0, sms1_0, \
                         sns0_0, sns1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sms0_0[k]
                 - f_1 * sms1_0[k]
                 + pb_x[k] * sns0_0[k]
                 - f_2 * pc_x[k] * sns1_0[k];

        t_1[k] = pb_y[k] * sns0_0[k]
                 - f_2 * pc_y[k] * sns1_0[k];

        t_2[k] = pb_z[k] * sns0_0[k]
                 - f_2 * pc_z[k] * sns1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pb_x, pb_y, pc_x, pc_y, sms0_3, sms1_3, sns0_2, sns0_3, \
                         sns1_2, sns1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * sms0_3[k]
                 - f_4 * sms1_3[k]
                 + pb_x[k] * sns0_3[k]
                 - f_2 * pc_x[k] * sns1_3[k];

        t_4[k] = pb_y[k] * sns0_2[k]
                 - f_2 * pc_y[k] * sns1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pb_x, pc_x, sms0_5, sms0_6, sms1_5, sms1_6, sns0_5, sns0_6, \
                         sns1_5, sns1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * sms0_5[k]
                 - f_4 * sms1_5[k]
                 + pb_x[k] * sns0_5[k]
                 - f_2 * pc_x[k] * sns1_5[k];

        t_6[k] = f_5 * sms0_6[k]
                 - f_6 * sms1_6[k]
                 + pb_x[k] * sns0_6[k]
                 - f_2 * pc_x[k] * sns1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pb_y, pb_z, pc_y, pc_z, sns0_3, sns0_5, sns1_3, \
                         sns1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pb_z[k] * sns0_3[k]
                 - f_2 * pc_z[k] * sns1_3[k];

        t_8[k] = pb_y[k] * sns0_5[k]
                 - f_2 * pc_y[k] * sns1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pb_x, pc_x, sms0_9, sms0_10, sms1_9, sms1_10, sns0_9, \
                         sns0_10, sns1_9, sns1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * sms0_9[k]
                 - f_6 * sms1_9[k]
                 + pb_x[k] * sns0_9[k]
                 - f_2 * pc_x[k] * sns1_9[k];

        t_10[k] = f_7 * sms0_10[k]
                  - f_8 * sms1_10[k]
                  + pb_x[k] * sns0_10[k]
                  - f_2 * pc_x[k] * sns1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pb_x, pb_z, pc_x, pc_z, sms0_12, sms1_12, sns0_6, \
                         sns0_12, sns1_6, sns1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * sns0_6[k]
                  - f_2 * pc_z[k] * sns1_6[k];

        t_12[k] = f_7 * sms0_12[k]
                  - f_8 * sms1_12[k]
                  + pb_x[k] * sns0_12[k]
                  - f_2 * pc_x[k] * sns1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pb_x, pb_y, pc_x, pc_y, sms0_14, sms1_14, sns0_9, \
                         sns0_14, sns1_9, sns1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pb_y[k] * sns0_9[k]
                  - f_2 * pc_y[k] * sns1_9[k];

        t_14[k] = f_7 * sms0_14[k]
                  - f_8 * sms1_14[k]
                  + pb_x[k] * sns0_14[k]
                  - f_2 * pc_x[k] * sns1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pb_x, pb_z, pc_x, pc_z, sms0_15, sms1_15, sns0_10, \
                         sns0_15, sns1_10, sns1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * sms0_15[k]
                  - f_10 * sms1_15[k]
                  + pb_x[k] * sns0_15[k]
                  - f_2 * pc_x[k] * sns1_15[k];

        t_16[k] = pb_z[k] * sns0_10[k]
                  - f_2 * pc_z[k] * sns1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pb_x, pc_x, sms0_17, sms0_18, sms1_17, sms1_18, sns0_17, \
                         sns0_18, sns1_17, sns1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * sms0_17[k]
                  - f_10 * sms1_17[k]
                  + pb_x[k] * sns0_17[k]
                  - f_2 * pc_x[k] * sns1_17[k];

        t_18[k] = f_9 * sms0_18[k]
                  - f_10 * sms1_18[k]
                  + pb_x[k] * sns0_18[k]
                  - f_2 * pc_x[k] * sns1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pb_x, pb_y, pc_x, pc_y, sms0_20, sms1_20, sns0_14, \
                         sns0_20, sns1_14, sns1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pb_y[k] * sns0_14[k]
                  - f_2 * pc_y[k] * sns1_14[k];

        t_20[k] = f_9 * sms0_20[k]
                  - f_10 * sms1_20[k]
                  + pb_x[k] * sns0_20[k]
                  - f_2 * pc_x[k] * sns1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pb_x, pb_z, pc_x, pc_z, sms0_21, sms1_21, sns0_15, \
                         sns0_21, sns1_15, sns1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * sms0_21[k]
                  - f_12 * sms1_21[k]
                  + pb_x[k] * sns0_21[k]
                  - f_2 * pc_x[k] * sns1_21[k];

        t_22[k] = pb_z[k] * sns0_15[k]
                  - f_2 * pc_z[k] * sns1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pb_x, pc_x, sms0_23, sms0_24, sms1_23, sms1_24, sns0_23, \
                         sns0_24, sns1_23, sns1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * sms0_23[k]
                  - f_12 * sms1_23[k]
                  + pb_x[k] * sns0_23[k]
                  - f_2 * pc_x[k] * sns1_23[k];

        t_24[k] = f_11 * sms0_24[k]
                  - f_12 * sms1_24[k]
                  + pb_x[k] * sns0_24[k]
                  - f_2 * pc_x[k] * sns1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pb_x, pb_y, pc_x, pc_y, sms0_25, sms1_25, sns0_20, \
                         sns0_25, sns1_20, sns1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * sms0_25[k]
                  - f_12 * sms1_25[k]
                  + pb_x[k] * sns0_25[k]
                  - f_2 * pc_x[k] * sns1_25[k];

        t_26[k] = pb_y[k] * sns0_20[k]
                  - f_2 * pc_y[k] * sns1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pb_x, pc_x, sms0_27, sms0_28, sms1_27, sms1_28, sns0_27, \
                         sns0_28, sns1_27, sns1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * sms0_27[k]
                  - f_12 * sms1_27[k]
                  + pb_x[k] * sns0_27[k]
                  - f_2 * pc_x[k] * sns1_27[k];

        t_28[k] = f_13 * sms0_28[k]
                  - f_14 * sms1_28[k]
                  + pb_x[k] * sns0_28[k]
                  - f_2 * pc_x[k] * sns1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pb_x, pb_z, pc_x, pc_z, sms0_30, sms1_30, sns0_21, \
                         sns0_30, sns1_21, sns1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_z[k] * sns0_21[k]
                  - f_2 * pc_z[k] * sns1_21[k];

        t_30[k] = f_13 * sms0_30[k]
                  - f_14 * sms1_30[k]
                  + pb_x[k] * sns0_30[k]
                  - f_2 * pc_x[k] * sns1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pb_x, pc_x, sms0_31, sms0_32, sms1_31, sms1_32, sns0_31, \
                         sns0_32, sns1_31, sns1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * sms0_31[k]
                  - f_14 * sms1_31[k]
                  + pb_x[k] * sns0_31[k]
                  - f_2 * pc_x[k] * sns1_31[k];

        t_32[k] = f_13 * sms0_32[k]
                  - f_14 * sms1_32[k]
                  + pb_x[k] * sns0_32[k]
                  - f_2 * pc_x[k] * sns1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pb_x, pb_y, pc_x, pc_y, sms0_33, sms1_33, sns0_27, \
                         sns0_33, sns1_27, sns1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * sms0_33[k]
                  - f_14 * sms1_33[k]
                  + pb_x[k] * sns0_33[k]
                  - f_2 * pc_x[k] * sns1_33[k];

        t_34[k] = pb_y[k] * sns0_27[k]
                  - f_2 * pc_y[k] * sns1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pb_x, pc_x, sms0_35, sms0_36, sms1_35, sms1_36, sns0_35, \
                         sns0_36, sns1_35, sns1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * sms0_35[k]
                  - f_14 * sms1_35[k]
                  + pb_x[k] * sns0_35[k]
                  - f_2 * pc_x[k] * sns1_35[k];

        t_36[k] = f_15 * sms0_36[k]
                  - f_16 * sms1_36[k]
                  + pb_x[k] * sns0_36[k]
                  - f_2 * pc_x[k] * sns1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pb_x, pb_z, pc_x, pc_z, sms0_38, sms1_38, sns0_28, \
                         sns0_38, sns1_28, sns1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_z[k] * sns0_28[k]
                  - f_2 * pc_z[k] * sns1_28[k];

        t_38[k] = f_15 * sms0_38[k]
                  - f_16 * sms1_38[k]
                  + pb_x[k] * sns0_38[k]
                  - f_2 * pc_x[k] * sns1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pb_x, pc_x, sms0_39, sms0_40, sms1_39, sms1_40, sns0_39, \
                         sns0_40, sns1_39, sns1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_15 * sms0_39[k]
                  - f_16 * sms1_39[k]
                  + pb_x[k] * sns0_39[k]
                  - f_2 * pc_x[k] * sns1_39[k];

        t_40[k] = f_15 * sms0_40[k]
                  - f_16 * sms1_40[k]
                  + pb_x[k] * sns0_40[k]
                  - f_2 * pc_x[k] * sns1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pb_x, pc_x, sms0_41, sms0_42, sms1_41, sms1_42, sns0_41, \
                         sns0_42, sns1_41, sns1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_15 * sms0_41[k]
                  - f_16 * sms1_41[k]
                  + pb_x[k] * sns0_41[k]
                  - f_2 * pc_x[k] * sns1_41[k];

        t_42[k] = f_15 * sms0_42[k]
                  - f_16 * sms1_42[k]
                  + pb_x[k] * sns0_42[k]
                  - f_2 * pc_x[k] * sns1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pb_x, pb_y, pc_x, pc_y, sms0_44, sms1_44, sns0_35, \
                         sns0_44, sns1_35, sns1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pb_y[k] * sns0_35[k]
                  - f_2 * pc_y[k] * sns1_35[k];

        t_44[k] = f_15 * sms0_44[k]
                  - f_16 * sms1_44[k]
                  + pb_x[k] * sns0_44[k]
                  - f_2 * pc_x[k] * sns1_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pb_x, pb_z, pc_x, pc_z, sms0_45, sms1_45, sns0_36, \
                         sns0_45, sns1_36, sns1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_17 * sms0_45[k]
                  - f_18 * sms1_45[k]
                  + pb_x[k] * sns0_45[k]
                  - f_2 * pc_x[k] * sns1_45[k];

        t_46[k] = pb_z[k] * sns0_36[k]
                  - f_2 * pc_z[k] * sns1_36[k];
    }

#pragma omp simd aligned(t_47, t_48, pb_x, pc_x, sms0_47, sms0_48, sms1_47, sms1_48, sns0_47, \
                         sns0_48, sns1_47, sns1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_17 * sms0_47[k]
                  - f_18 * sms1_47[k]
                  + pb_x[k] * sns0_47[k]
                  - f_2 * pc_x[k] * sns1_47[k];

        t_48[k] = f_17 * sms0_48[k]
                  - f_18 * sms1_48[k]
                  + pb_x[k] * sns0_48[k]
                  - f_2 * pc_x[k] * sns1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pb_x, pc_x, sms0_49, sms0_50, sms1_49, sms1_50, sns0_49, \
                         sns0_50, sns1_49, sns1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_17 * sms0_49[k]
                  - f_18 * sms1_49[k]
                  + pb_x[k] * sns0_49[k]
                  - f_2 * pc_x[k] * sns1_49[k];

        t_50[k] = f_17 * sms0_50[k]
                  - f_18 * sms1_50[k]
                  + pb_x[k] * sns0_50[k]
                  - f_2 * pc_x[k] * sns1_50[k];
    }

#pragma omp simd aligned(t_51, t_52, pb_x, pc_x, sms0_51, sms0_52, sms1_51, sms1_52, sns0_51, \
                         sns0_52, sns1_51, sns1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_17 * sms0_51[k]
                  - f_18 * sms1_51[k]
                  + pb_x[k] * sns0_51[k]
                  - f_2 * pc_x[k] * sns1_51[k];

        t_52[k] = f_17 * sms0_52[k]
                  - f_18 * sms1_52[k]
                  + pb_x[k] * sns0_52[k]
                  - f_2 * pc_x[k] * sns1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, pc_x, pc_y, sms0_54, sms1_54, sns0_44, \
                         sns0_54, sns0_55, sns1_44, sns1_54, sns1_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * sns0_44[k]
                  - f_2 * pc_y[k] * sns1_44[k];

        t_54[k] = f_17 * sms0_54[k]
                  - f_18 * sms1_54[k]
                  + pb_x[k] * sns0_54[k]
                  - f_2 * pc_x[k] * sns1_54[k];

        t_55[k] = pb_x[k] * sns0_55[k]
                  - f_2 * pc_x[k] * sns1_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pb_x, pc_x, sns0_56, sns0_57, sns0_58, \
                         sns0_59, sns1_56, sns1_57, sns1_58, sns1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_x[k] * sns0_56[k]
                  - f_2 * pc_x[k] * sns1_56[k];

        t_57[k] = pb_x[k] * sns0_57[k]
                  - f_2 * pc_x[k] * sns1_57[k];

        t_58[k] = pb_x[k] * sns0_58[k]
                  - f_2 * pc_x[k] * sns1_58[k];

        t_59[k] = pb_x[k] * sns0_59[k]
                  - f_2 * pc_x[k] * sns1_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_x, pc_x, sns0_60, sns0_61, sns0_62, \
                         sns0_63, sns1_60, sns1_61, sns1_62, sns1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pb_x[k] * sns0_60[k]
                  - f_2 * pc_x[k] * sns1_60[k];

        t_61[k] = pb_x[k] * sns0_61[k]
                  - f_2 * pc_x[k] * sns1_61[k];

        t_62[k] = pb_x[k] * sns0_62[k]
                  - f_2 * pc_x[k] * sns1_62[k];

        t_63[k] = pb_x[k] * sns0_63[k]
                  - f_2 * pc_x[k] * sns1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pc_x, pc_y, sms0_45, sms1_45, sns0_55, \
                         sns0_64, sns0_65, sns1_55, sns1_64, sns1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = pb_x[k] * sns0_64[k]
                  - f_2 * pc_x[k] * sns1_64[k];

        t_65[k] = pb_x[k] * sns0_65[k]
                  - f_2 * pc_x[k] * sns1_65[k];

        t_66[k] = f_0 * sms0_45[k]
                  - f_1 * sms1_45[k]
                  + pb_y[k] * sns0_55[k]
                  - f_2 * pc_y[k] * sns1_55[k];
    }

#pragma omp simd aligned(t_67, t_68, pb_y, pb_z, pc_y, pc_z, sms0_47, sms1_47, sns0_55, \
                         sns0_57, sns1_55, sns1_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = pb_z[k] * sns0_55[k]
                  - f_2 * pc_z[k] * sns1_55[k];

        t_68[k] = f_3 * sms0_47[k]
                  - f_4 * sms1_47[k]
                  + pb_y[k] * sns0_57[k]
                  - f_2 * pc_y[k] * sns1_57[k];
    }

#pragma omp simd aligned(t_69, t_70, pb_y, pc_y, sms0_48, sms0_49, sms1_48, sms1_49, sns0_58, \
                         sns0_59, sns1_58, sns1_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_5 * sms0_48[k]
                  - f_6 * sms1_48[k]
                  + pb_y[k] * sns0_58[k]
                  - f_2 * pc_y[k] * sns1_58[k];

        t_70[k] = f_7 * sms0_49[k]
                  - f_8 * sms1_49[k]
                  + pb_y[k] * sns0_59[k]
                  - f_2 * pc_y[k] * sns1_59[k];
    }

#pragma omp simd aligned(t_71, t_72, pb_y, pc_y, sms0_50, sms0_51, sms1_50, sms1_51, sns0_60, \
                         sns0_61, sns1_60, sns1_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_9 * sms0_50[k]
                  - f_10 * sms1_50[k]
                  + pb_y[k] * sns0_60[k]
                  - f_2 * pc_y[k] * sns1_60[k];

        t_72[k] = f_11 * sms0_51[k]
                  - f_12 * sms1_51[k]
                  + pb_y[k] * sns0_61[k]
                  - f_2 * pc_y[k] * sns1_61[k];
    }

#pragma omp simd aligned(t_73, t_74, pb_y, pc_y, sms0_52, sms0_53, sms1_52, sms1_53, sns0_62, \
                         sns0_63, sns1_62, sns1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_13 * sms0_52[k]
                  - f_14 * sms1_52[k]
                  + pb_y[k] * sns0_62[k]
                  - f_2 * pc_y[k] * sns1_62[k];

        t_74[k] = f_15 * sms0_53[k]
                  - f_16 * sms1_53[k]
                  + pb_y[k] * sns0_63[k]
                  - f_2 * pc_y[k] * sns1_63[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pb_y, pb_z, pc_y, pc_z, sms0_54, sms1_54, sns0_64, \
                         sns0_65, sns1_64, sns1_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_17 * sms0_54[k]
                  - f_18 * sms1_54[k]
                  + pb_y[k] * sns0_64[k]
                  - f_2 * pc_y[k] * sns1_64[k];

        t_76[k] = pb_y[k] * sns0_65[k]
                  - f_2 * pc_y[k] * sns1_65[k];

        t_77[k] = f_0 * sms0_54[k]
                  - f_1 * sms1_54[k]
                  + pb_z[k] * sns0_65[k]
                  - f_2 * pc_z[k] * sns1_65[k];
    }
}

}  // namespace simdt3ceri
