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


#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_smd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sld0,
                                                          const size_t slp, const size_t sld1,
                                                          const size_t sms0, const size_t sms1,
                                                          const size_t smp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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
    auto *t_84 = buffer.data(target + 84);
    auto *t_85 = buffer.data(target + 85);
    auto *t_86 = buffer.data(target + 86);
    auto *t_87 = buffer.data(target + 87);
    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);
    auto *t_91 = buffer.data(target + 91);
    auto *t_92 = buffer.data(target + 92);
    auto *t_93 = buffer.data(target + 93);
    auto *t_94 = buffer.data(target + 94);
    auto *t_95 = buffer.data(target + 95);
    auto *t_96 = buffer.data(target + 96);
    auto *t_97 = buffer.data(target + 97);
    auto *t_98 = buffer.data(target + 98);
    auto *t_99 = buffer.data(target + 99);
    auto *t_100 = buffer.data(target + 100);
    auto *t_101 = buffer.data(target + 101);
    auto *t_102 = buffer.data(target + 102);
    auto *t_103 = buffer.data(target + 103);
    auto *t_104 = buffer.data(target + 104);
    auto *t_105 = buffer.data(target + 105);
    auto *t_106 = buffer.data(target + 106);
    auto *t_107 = buffer.data(target + 107);
    auto *t_108 = buffer.data(target + 108);
    auto *t_109 = buffer.data(target + 109);
    auto *t_110 = buffer.data(target + 110);
    auto *t_111 = buffer.data(target + 111);
    auto *t_112 = buffer.data(target + 112);
    auto *t_113 = buffer.data(target + 113);
    auto *t_114 = buffer.data(target + 114);
    auto *t_115 = buffer.data(target + 115);
    auto *t_116 = buffer.data(target + 116);
    auto *t_117 = buffer.data(target + 117);
    auto *t_118 = buffer.data(target + 118);
    auto *t_119 = buffer.data(target + 119);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sld0_0 = buffer.data(sld0 + 0);
    const auto *sld0_3 = buffer.data(sld0 + 3);
    const auto *sld0_5 = buffer.data(sld0 + 5);
    const auto *sld0_9 = buffer.data(sld0 + 9);
    const auto *sld0_12 = buffer.data(sld0 + 12);
    const auto *sld0_17 = buffer.data(sld0 + 17);
    const auto *sld0_18 = buffer.data(sld0 + 18);
    const auto *sld0_21 = buffer.data(sld0 + 21);
    const auto *sld0_30 = buffer.data(sld0 + 30);
    const auto *sld0_35 = buffer.data(sld0 + 35);
    const auto *sld0_36 = buffer.data(sld0 + 36);
    const auto *sld0_39 = buffer.data(sld0 + 39);
    const auto *sld0_54 = buffer.data(sld0 + 54);
    const auto *sld0_59 = buffer.data(sld0 + 59);
    const auto *sld0_60 = buffer.data(sld0 + 60);
    const auto *sld0_63 = buffer.data(sld0 + 63);
    const auto *sld0_84 = buffer.data(sld0 + 84);
    const auto *sld0_89 = buffer.data(sld0 + 89);

    const auto *slp_0 = buffer.data(slp + 0);
    const auto *slp_1 = buffer.data(slp + 1);
    const auto *slp_2 = buffer.data(slp + 2);
    const auto *slp_4 = buffer.data(slp + 4);
    const auto *slp_5 = buffer.data(slp + 5);
    const auto *slp_7 = buffer.data(slp + 7);
    const auto *slp_8 = buffer.data(slp + 8);
    const auto *slp_9 = buffer.data(slp + 9);
    const auto *slp_10 = buffer.data(slp + 10);
    const auto *slp_11 = buffer.data(slp + 11);
    const auto *slp_13 = buffer.data(slp + 13);
    const auto *slp_14 = buffer.data(slp + 14);
    const auto *slp_15 = buffer.data(slp + 15);
    const auto *slp_16 = buffer.data(slp + 16);
    const auto *slp_17 = buffer.data(slp + 17);
    const auto *slp_18 = buffer.data(slp + 18);
    const auto *slp_19 = buffer.data(slp + 19);
    const auto *slp_20 = buffer.data(slp + 20);
    const auto *slp_22 = buffer.data(slp + 22);
    const auto *slp_23 = buffer.data(slp + 23);
    const auto *slp_25 = buffer.data(slp + 25);
    const auto *slp_26 = buffer.data(slp + 26);
    const auto *slp_27 = buffer.data(slp + 27);
    const auto *slp_28 = buffer.data(slp + 28);
    const auto *slp_29 = buffer.data(slp + 29);
    const auto *slp_30 = buffer.data(slp + 30);
    const auto *slp_31 = buffer.data(slp + 31);
    const auto *slp_32 = buffer.data(slp + 32);
    const auto *slp_34 = buffer.data(slp + 34);
    const auto *slp_35 = buffer.data(slp + 35);
    const auto *slp_36 = buffer.data(slp + 36);
    const auto *slp_37 = buffer.data(slp + 37);
    const auto *slp_38 = buffer.data(slp + 38);
    const auto *slp_40 = buffer.data(slp + 40);
    const auto *slp_41 = buffer.data(slp + 41);
    const auto *slp_42 = buffer.data(slp + 42);
    const auto *slp_43 = buffer.data(slp + 43);
    const auto *slp_44 = buffer.data(slp + 44);
    const auto *slp_45 = buffer.data(slp + 45);
    const auto *slp_46 = buffer.data(slp + 46);
    const auto *slp_47 = buffer.data(slp + 47);
    const auto *slp_49 = buffer.data(slp + 49);
    const auto *slp_50 = buffer.data(slp + 50);
    const auto *slp_51 = buffer.data(slp + 51);
    const auto *slp_52 = buffer.data(slp + 52);
    const auto *slp_53 = buffer.data(slp + 53);
    const auto *slp_54 = buffer.data(slp + 54);
    const auto *slp_55 = buffer.data(slp + 55);
    const auto *slp_56 = buffer.data(slp + 56);
    const auto *slp_58 = buffer.data(slp + 58);
    const auto *slp_59 = buffer.data(slp + 59);

    const auto *sld1_0 = buffer.data(sld1 + 0);
    const auto *sld1_3 = buffer.data(sld1 + 3);
    const auto *sld1_5 = buffer.data(sld1 + 5);
    const auto *sld1_9 = buffer.data(sld1 + 9);
    const auto *sld1_12 = buffer.data(sld1 + 12);
    const auto *sld1_17 = buffer.data(sld1 + 17);
    const auto *sld1_18 = buffer.data(sld1 + 18);
    const auto *sld1_21 = buffer.data(sld1 + 21);
    const auto *sld1_30 = buffer.data(sld1 + 30);
    const auto *sld1_35 = buffer.data(sld1 + 35);
    const auto *sld1_36 = buffer.data(sld1 + 36);
    const auto *sld1_39 = buffer.data(sld1 + 39);
    const auto *sld1_54 = buffer.data(sld1 + 54);
    const auto *sld1_59 = buffer.data(sld1 + 59);
    const auto *sld1_60 = buffer.data(sld1 + 60);
    const auto *sld1_63 = buffer.data(sld1 + 63);
    const auto *sld1_84 = buffer.data(sld1 + 84);
    const auto *sld1_89 = buffer.data(sld1 + 89);

    const auto *sms0_0 = buffer.data(sms0 + 0);
    const auto *sms0_1 = buffer.data(sms0 + 1);
    const auto *sms0_2 = buffer.data(sms0 + 2);
    const auto *sms0_3 = buffer.data(sms0 + 3);
    const auto *sms0_5 = buffer.data(sms0 + 5);
    const auto *sms0_6 = buffer.data(sms0 + 6);
    const auto *sms0_7 = buffer.data(sms0 + 7);
    const auto *sms0_8 = buffer.data(sms0 + 8);
    const auto *sms0_9 = buffer.data(sms0 + 9);
    const auto *sms0_10 = buffer.data(sms0 + 10);
    const auto *sms0_11 = buffer.data(sms0 + 11);
    const auto *sms0_12 = buffer.data(sms0 + 12);
    const auto *sms0_13 = buffer.data(sms0 + 13);
    const auto *sms0_14 = buffer.data(sms0 + 14);
    const auto *sms0_15 = buffer.data(sms0 + 15);
    const auto *sms0_16 = buffer.data(sms0 + 16);
    const auto *sms0_17 = buffer.data(sms0 + 17);
    const auto *sms0_18 = buffer.data(sms0 + 18);
    const auto *sms0_19 = buffer.data(sms0 + 19);

    const auto *sms1_0 = buffer.data(sms1 + 0);
    const auto *sms1_1 = buffer.data(sms1 + 1);
    const auto *sms1_2 = buffer.data(sms1 + 2);
    const auto *sms1_3 = buffer.data(sms1 + 3);
    const auto *sms1_5 = buffer.data(sms1 + 5);
    const auto *sms1_6 = buffer.data(sms1 + 6);
    const auto *sms1_7 = buffer.data(sms1 + 7);
    const auto *sms1_8 = buffer.data(sms1 + 8);
    const auto *sms1_9 = buffer.data(sms1 + 9);
    const auto *sms1_10 = buffer.data(sms1 + 10);
    const auto *sms1_11 = buffer.data(sms1 + 11);
    const auto *sms1_12 = buffer.data(sms1 + 12);
    const auto *sms1_13 = buffer.data(sms1 + 13);
    const auto *sms1_14 = buffer.data(sms1 + 14);
    const auto *sms1_15 = buffer.data(sms1 + 15);
    const auto *sms1_16 = buffer.data(sms1 + 16);
    const auto *sms1_17 = buffer.data(sms1 + 17);
    const auto *sms1_18 = buffer.data(sms1 + 18);
    const auto *sms1_19 = buffer.data(sms1 + 19);

    const auto *smp_0 = buffer.data(smp + 0);
    const auto *smp_1 = buffer.data(smp + 1);
    const auto *smp_2 = buffer.data(smp + 2);
    const auto *smp_4 = buffer.data(smp + 4);
    const auto *smp_5 = buffer.data(smp + 5);
    const auto *smp_7 = buffer.data(smp + 7);
    const auto *smp_8 = buffer.data(smp + 8);
    const auto *smp_9 = buffer.data(smp + 9);
    const auto *smp_10 = buffer.data(smp + 10);
    const auto *smp_11 = buffer.data(smp + 11);
    const auto *smp_13 = buffer.data(smp + 13);
    const auto *smp_14 = buffer.data(smp + 14);
    const auto *smp_15 = buffer.data(smp + 15);
    const auto *smp_16 = buffer.data(smp + 16);
    const auto *smp_17 = buffer.data(smp + 17);
    const auto *smp_18 = buffer.data(smp + 18);
    const auto *smp_19 = buffer.data(smp + 19);
    const auto *smp_20 = buffer.data(smp + 20);
    const auto *smp_22 = buffer.data(smp + 22);
    const auto *smp_23 = buffer.data(smp + 23);
    const auto *smp_25 = buffer.data(smp + 25);
    const auto *smp_26 = buffer.data(smp + 26);
    const auto *smp_27 = buffer.data(smp + 27);
    const auto *smp_28 = buffer.data(smp + 28);
    const auto *smp_29 = buffer.data(smp + 29);
    const auto *smp_30 = buffer.data(smp + 30);
    const auto *smp_31 = buffer.data(smp + 31);
    const auto *smp_32 = buffer.data(smp + 32);
    const auto *smp_34 = buffer.data(smp + 34);
    const auto *smp_35 = buffer.data(smp + 35);
    const auto *smp_36 = buffer.data(smp + 36);
    const auto *smp_37 = buffer.data(smp + 37);
    const auto *smp_38 = buffer.data(smp + 38);
    const auto *smp_40 = buffer.data(smp + 40);
    const auto *smp_41 = buffer.data(smp + 41);
    const auto *smp_42 = buffer.data(smp + 42);
    const auto *smp_43 = buffer.data(smp + 43);
    const auto *smp_44 = buffer.data(smp + 44);
    const auto *smp_45 = buffer.data(smp + 45);
    const auto *smp_46 = buffer.data(smp + 46);
    const auto *smp_47 = buffer.data(smp + 47);
    const auto *smp_49 = buffer.data(smp + 49);
    const auto *smp_50 = buffer.data(smp + 50);
    const auto *smp_51 = buffer.data(smp + 51);
    const auto *smp_52 = buffer.data(smp + 52);
    const auto *smp_53 = buffer.data(smp + 53);
    const auto *smp_54 = buffer.data(smp + 54);
    const auto *smp_55 = buffer.data(smp + 55);
    const auto *smp_56 = buffer.data(smp + 56);
    const auto *smp_58 = buffer.data(smp + 58);
    const auto *smp_59 = buffer.data(smp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, slp_0, slp_1, slp_2, sms0_0, \
                         sms1_0, smp_0, smp_1, smp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * slp_0[k]
                 + f_1 * sms0_0[k]
                 - f_2 * sms1_0[k]
                 + f_3 * pc_x[k] * smp_0[k];

        t_1[k] = f_0 * slp_1[k]
                 + f_3 * pc_x[k] * smp_1[k];

        t_2[k] = f_0 * slp_2[k]
                 + f_3 * pc_x[k] * smp_2[k];

        t_3[k] = f_1 * sms0_0[k]
                 - f_2 * sms1_0[k]
                 + f_3 * pc_y[k] * smp_1[k];

        t_4[k] = f_3 * pc_y[k] * smp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, sld0_0, slp_4, sld1_0, sms0_0, \
                         sms1_0, smp_2, smp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sms0_0[k]
                 - f_2 * sms1_0[k]
                 + f_3 * pc_z[k] * smp_2[k];

        t_6[k] = pb_y[k] * sld0_0[k]
                 - f_4 * pc_y[k] * sld1_0[k];

        t_7[k] = f_5 * slp_4[k]
                 + f_3 * pc_x[k] * smp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, sld0_5, slp_1, slp_2, slp_5, \
                         sld1_5, sms0_1, sms1_1, smp_4, smp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * slp_5[k]
                 + f_3 * pc_x[k] * smp_5[k];

        t_9[k] = f_6 * slp_1[k]
                 + f_1 * sms0_1[k]
                 - f_2 * sms1_1[k]
                 + f_3 * pc_y[k] * smp_4[k];

        t_10[k] = f_6 * slp_2[k]
                  + f_3 * pc_y[k] * smp_5[k];

        t_11[k] = pb_y[k] * sld0_5[k]
                  - f_4 * pc_y[k] * sld1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, sld0_0, sld0_3, slp_7, \
                         slp_8, sld1_0, sld1_3, smp_7, smp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * sld0_0[k]
                  - f_4 * pc_z[k] * sld1_0[k];

        t_13[k] = f_5 * slp_7[k]
                  + f_3 * pc_x[k] * smp_7[k];

        t_14[k] = f_5 * slp_8[k]
                  + f_3 * pc_x[k] * smp_8[k];

        t_15[k] = pb_z[k] * sld0_3[k]
                  - f_4 * pc_z[k] * sld1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, slp_2, slp_9, sms0_2, sms0_3, \
                         sms1_2, sms1_3, smp_8, smp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * smp_8[k];

        t_17[k] = f_6 * slp_2[k]
                  + f_1 * sms0_2[k]
                  - f_2 * sms1_2[k]
                  + f_3 * pc_z[k] * smp_8[k];

        t_18[k] = f_7 * slp_9[k]
                  + f_1 * sms0_3[k]
                  - f_2 * sms1_3[k]
                  + f_3 * pc_x[k] * smp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, slp_4, slp_5, slp_10, \
                         slp_11, sms0_3, sms1_3, smp_10, smp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * slp_10[k]
                  + f_3 * pc_x[k] * smp_10[k];

        t_20[k] = f_7 * slp_11[k]
                  + f_3 * pc_x[k] * smp_11[k];

        t_21[k] = f_8 * slp_4[k]
                  + f_1 * sms0_3[k]
                  - f_2 * sms1_3[k]
                  + f_3 * pc_y[k] * smp_10[k];

        t_22[k] = f_8 * slp_5[k]
                  + f_3 * pc_y[k] * smp_11[k];

        t_23[k] = f_1 * sms0_3[k]
                  - f_2 * sms1_3[k]
                  + f_3 * pc_z[k] * smp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, sld0_12, slp_13, slp_14, sld1_12, \
                         smp_13, smp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sld0_12[k]
                  - f_4 * pc_y[k] * sld1_12[k];

        t_25[k] = f_7 * slp_13[k]
                  + f_3 * pc_x[k] * smp_13[k];

        t_26[k] = f_7 * slp_14[k]
                  + f_3 * pc_x[k] * smp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, sld0_9, sld0_17, slp_8, \
                         sld1_9, sld1_17, smp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * sld0_9[k]
                  - f_4 * pc_z[k] * sld1_9[k];

        t_28[k] = f_6 * slp_8[k]
                  + f_3 * pc_y[k] * smp_14[k];

        t_29[k] = pb_y[k] * sld0_17[k]
                  - f_4 * pc_y[k] * sld1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, slp_15, slp_16, slp_17, \
                         sms0_5, sms1_5, smp_15, smp_16, smp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * slp_15[k]
                  + f_1 * sms0_5[k]
                  - f_2 * sms1_5[k]
                  + f_3 * pc_x[k] * smp_15[k];

        t_31[k] = f_7 * slp_16[k]
                  + f_3 * pc_x[k] * smp_16[k];

        t_32[k] = f_7 * slp_17[k]
                  + f_3 * pc_x[k] * smp_17[k];

        t_33[k] = f_1 * sms0_5[k]
                  - f_2 * sms1_5[k]
                  + f_3 * pc_y[k] * smp_16[k];

        t_34[k] = f_3 * pc_y[k] * smp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, slp_8, slp_18, slp_19, sms0_5, sms0_6, \
                         sms1_5, sms1_6, smp_17, smp_18, smp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * slp_8[k]
                  + f_1 * sms0_5[k]
                  - f_2 * sms1_5[k]
                  + f_3 * pc_z[k] * smp_17[k];

        t_36[k] = f_9 * slp_18[k]
                  + f_1 * sms0_6[k]
                  - f_2 * sms1_6[k]
                  + f_3 * pc_x[k] * smp_18[k];

        t_37[k] = f_9 * slp_19[k]
                  + f_3 * pc_x[k] * smp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, slp_10, slp_11, slp_20, \
                         sms0_6, sms1_6, smp_19, smp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * slp_20[k]
                  + f_3 * pc_x[k] * smp_20[k];

        t_39[k] = f_10 * slp_10[k]
                  + f_1 * sms0_6[k]
                  - f_2 * sms1_6[k]
                  + f_3 * pc_y[k] * smp_19[k];

        t_40[k] = f_10 * slp_11[k]
                  + f_3 * pc_y[k] * smp_20[k];

        t_41[k] = f_1 * sms0_6[k]
                  - f_2 * sms1_6[k]
                  + f_3 * pc_z[k] * smp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, sld0_18, sld0_21, slp_22, \
                         slp_23, sld1_18, sld1_21, smp_22, smp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sld0_18[k]
                  - f_4 * pc_z[k] * sld1_18[k];

        t_43[k] = f_9 * slp_22[k]
                  + f_3 * pc_x[k] * smp_22[k];

        t_44[k] = f_9 * slp_23[k]
                  + f_3 * pc_x[k] * smp_23[k];

        t_45[k] = pb_z[k] * sld0_21[k]
                  - f_4 * pc_z[k] * sld1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, sld0_30, slp_11, slp_14, sld1_30, \
                         sms0_7, sms1_7, smp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * slp_14[k]
                  + f_3 * pc_y[k] * smp_23[k];

        t_47[k] = f_6 * slp_11[k]
                  + f_1 * sms0_7[k]
                  - f_2 * sms1_7[k]
                  + f_3 * pc_z[k] * smp_23[k];

        t_48[k] = pb_y[k] * sld0_30[k]
                  - f_4 * pc_y[k] * sld1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, slp_16, slp_17, slp_25, slp_26, \
                         sms0_8, sms1_8, smp_25, smp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * slp_25[k]
                  + f_3 * pc_x[k] * smp_25[k];

        t_50[k] = f_9 * slp_26[k]
                  + f_3 * pc_x[k] * smp_26[k];

        t_51[k] = f_6 * slp_16[k]
                  + f_1 * sms0_8[k]
                  - f_2 * sms1_8[k]
                  + f_3 * pc_y[k] * smp_25[k];

        t_52[k] = f_6 * slp_17[k]
                  + f_3 * pc_y[k] * smp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, sld0_35, slp_27, slp_28, sld1_35, \
                         sms0_9, sms1_9, smp_27, smp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * sld0_35[k]
                  - f_4 * pc_y[k] * sld1_35[k];

        t_54[k] = f_9 * slp_27[k]
                  + f_1 * sms0_9[k]
                  - f_2 * sms1_9[k]
                  + f_3 * pc_x[k] * smp_27[k];

        t_55[k] = f_9 * slp_28[k]
                  + f_3 * pc_x[k] * smp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, slp_17, slp_29, sms0_9, \
                         sms1_9, smp_28, smp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * slp_29[k]
                  + f_3 * pc_x[k] * smp_29[k];

        t_57[k] = f_1 * sms0_9[k]
                  - f_2 * sms1_9[k]
                  + f_3 * pc_y[k] * smp_28[k];

        t_58[k] = f_3 * pc_y[k] * smp_29[k];

        t_59[k] = f_10 * slp_17[k]
                  + f_1 * sms0_9[k]
                  - f_2 * sms1_9[k]
                  + f_3 * pc_z[k] * smp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, slp_19, slp_30, slp_31, slp_32, \
                         sms0_10, sms1_10, smp_30, smp_31, smp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * slp_30[k]
                  + f_1 * sms0_10[k]
                  - f_2 * sms1_10[k]
                  + f_3 * pc_x[k] * smp_30[k];

        t_61[k] = f_11 * slp_31[k]
                  + f_3 * pc_x[k] * smp_31[k];

        t_62[k] = f_11 * slp_32[k]
                  + f_3 * pc_x[k] * smp_32[k];

        t_63[k] = f_12 * slp_19[k]
                  + f_1 * sms0_10[k]
                  - f_2 * sms1_10[k]
                  + f_3 * pc_y[k] * smp_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, sld0_36, slp_20, \
                         slp_34, sld1_36, sms0_10, sms1_10, smp_32, \
                         smp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_12 * slp_20[k]
                  + f_3 * pc_y[k] * smp_32[k];

        t_65[k] = f_1 * sms0_10[k]
                  - f_2 * sms1_10[k]
                  + f_3 * pc_z[k] * smp_32[k];

        t_66[k] = pb_z[k] * sld0_36[k]
                  - f_4 * pc_z[k] * sld1_36[k];

        t_67[k] = f_11 * slp_34[k]
                  + f_3 * pc_x[k] * smp_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, sld0_39, slp_20, \
                         slp_23, slp_35, sld1_39, sms0_11, sms1_11, \
                         smp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * slp_35[k]
                  + f_3 * pc_x[k] * smp_35[k];

        t_69[k] = pb_z[k] * sld0_39[k]
                  - f_4 * pc_z[k] * sld1_39[k];

        t_70[k] = f_10 * slp_23[k]
                  + f_3 * pc_y[k] * smp_35[k];

        t_71[k] = f_6 * slp_20[k]
                  + f_1 * sms0_11[k]
                  - f_2 * sms1_11[k]
                  + f_3 * pc_z[k] * smp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, slp_25, slp_36, slp_37, slp_38, \
                         sms0_12, sms1_12, smp_36, smp_37, smp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * slp_36[k]
                  + f_1 * sms0_12[k]
                  - f_2 * sms1_12[k]
                  + f_3 * pc_x[k] * smp_36[k];

        t_73[k] = f_11 * slp_37[k]
                  + f_3 * pc_x[k] * smp_37[k];

        t_74[k] = f_11 * slp_38[k]
                  + f_3 * pc_x[k] * smp_38[k];

        t_75[k] = f_8 * slp_25[k]
                  + f_1 * sms0_12[k]
                  - f_2 * sms1_12[k]
                  + f_3 * pc_y[k] * smp_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, sld0_54, slp_23, slp_26, sld1_54, \
                         sms0_12, sms1_12, smp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * slp_26[k]
                  + f_3 * pc_y[k] * smp_38[k];

        t_77[k] = f_8 * slp_23[k]
                  + f_1 * sms0_12[k]
                  - f_2 * sms1_12[k]
                  + f_3 * pc_z[k] * smp_38[k];

        t_78[k] = pb_y[k] * sld0_54[k]
                  - f_4 * pc_y[k] * sld1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, slp_28, slp_29, slp_40, slp_41, \
                         sms0_13, sms1_13, smp_40, smp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * slp_40[k]
                  + f_3 * pc_x[k] * smp_40[k];

        t_80[k] = f_11 * slp_41[k]
                  + f_3 * pc_x[k] * smp_41[k];

        t_81[k] = f_6 * slp_28[k]
                  + f_1 * sms0_13[k]
                  - f_2 * sms1_13[k]
                  + f_3 * pc_y[k] * smp_40[k];

        t_82[k] = f_6 * slp_29[k]
                  + f_3 * pc_y[k] * smp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, sld0_59, slp_42, slp_43, sld1_59, \
                         sms0_14, sms1_14, smp_42, smp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * sld0_59[k]
                  - f_4 * pc_y[k] * sld1_59[k];

        t_84[k] = f_11 * slp_42[k]
                  + f_1 * sms0_14[k]
                  - f_2 * sms1_14[k]
                  + f_3 * pc_x[k] * smp_42[k];

        t_85[k] = f_11 * slp_43[k]
                  + f_3 * pc_x[k] * smp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, slp_29, slp_44, sms0_14, \
                         sms1_14, smp_43, smp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * slp_44[k]
                  + f_3 * pc_x[k] * smp_44[k];

        t_87[k] = f_1 * sms0_14[k]
                  - f_2 * sms1_14[k]
                  + f_3 * pc_y[k] * smp_43[k];

        t_88[k] = f_3 * pc_y[k] * smp_44[k];

        t_89[k] = f_12 * slp_29[k]
                  + f_1 * sms0_14[k]
                  - f_2 * sms1_14[k]
                  + f_3 * pc_z[k] * smp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, slp_31, slp_45, slp_46, slp_47, \
                         sms0_15, sms1_15, smp_45, smp_46, smp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_12 * slp_45[k]
                  + f_1 * sms0_15[k]
                  - f_2 * sms1_15[k]
                  + f_3 * pc_x[k] * smp_45[k];

        t_91[k] = f_12 * slp_46[k]
                  + f_3 * pc_x[k] * smp_46[k];

        t_92[k] = f_12 * slp_47[k]
                  + f_3 * pc_x[k] * smp_47[k];

        t_93[k] = f_11 * slp_31[k]
                  + f_1 * sms0_15[k]
                  - f_2 * sms1_15[k]
                  + f_3 * pc_y[k] * smp_46[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_z, pc_x, pc_y, pc_z, sld0_60, slp_32, \
                         slp_49, sld1_60, sms0_15, sms1_15, smp_47, \
                         smp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_11 * slp_32[k]
                  + f_3 * pc_y[k] * smp_47[k];

        t_95[k] = f_1 * sms0_15[k]
                  - f_2 * sms1_15[k]
                  + f_3 * pc_z[k] * smp_47[k];

        t_96[k] = pb_z[k] * sld0_60[k]
                  - f_4 * pc_z[k] * sld1_60[k];

        t_97[k] = f_12 * slp_49[k]
                  + f_3 * pc_x[k] * smp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_y, pc_z, sld0_63, slp_32, \
                         slp_35, slp_50, sld1_63, sms0_16, sms1_16, \
                         smp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_12 * slp_50[k]
                  + f_3 * pc_x[k] * smp_50[k];

        t_99[k] = pb_z[k] * sld0_63[k]
                  - f_4 * pc_z[k] * sld1_63[k];

        t_100[k] = f_12 * slp_35[k]
                   + f_3 * pc_y[k] * smp_50[k];

        t_101[k] = f_6 * slp_32[k]
                   + f_1 * sms0_16[k]
                   - f_2 * sms1_16[k]
                   + f_3 * pc_z[k] * smp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, slp_37, slp_51, slp_52, \
                         slp_53, sms0_17, sms1_17, smp_51, smp_52, \
                         smp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * slp_51[k]
                   + f_1 * sms0_17[k]
                   - f_2 * sms1_17[k]
                   + f_3 * pc_x[k] * smp_51[k];

        t_103[k] = f_12 * slp_52[k]
                   + f_3 * pc_x[k] * smp_52[k];

        t_104[k] = f_12 * slp_53[k]
                   + f_3 * pc_x[k] * smp_53[k];

        t_105[k] = f_10 * slp_37[k]
                   + f_1 * sms0_17[k]
                   - f_2 * sms1_17[k]
                   + f_3 * pc_y[k] * smp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, slp_35, slp_38, slp_54, \
                         sms0_17, sms0_18, sms1_17, sms1_18, smp_53, \
                         smp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * slp_38[k]
                   + f_3 * pc_y[k] * smp_53[k];

        t_107[k] = f_8 * slp_35[k]
                   + f_1 * sms0_17[k]
                   - f_2 * sms1_17[k]
                   + f_3 * pc_z[k] * smp_53[k];

        t_108[k] = f_12 * slp_54[k]
                   + f_1 * sms0_18[k]
                   - f_2 * sms1_18[k]
                   + f_3 * pc_x[k] * smp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, slp_40, slp_41, slp_55, \
                         slp_56, sms0_18, sms1_18, smp_55, smp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_12 * slp_55[k]
                   + f_3 * pc_x[k] * smp_55[k];

        t_110[k] = f_12 * slp_56[k]
                   + f_3 * pc_x[k] * smp_56[k];

        t_111[k] = f_8 * slp_40[k]
                   + f_1 * sms0_18[k]
                   - f_2 * sms1_18[k]
                   + f_3 * pc_y[k] * smp_55[k];

        t_112[k] = f_8 * slp_41[k]
                   + f_3 * pc_y[k] * smp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pc_x, pc_y, pc_z, sld0_84, slp_38, slp_58, \
                         sld1_84, sms0_18, sms1_18, smp_56, smp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * slp_38[k]
                   + f_1 * sms0_18[k]
                   - f_2 * sms1_18[k]
                   + f_3 * pc_z[k] * smp_56[k];

        t_114[k] = pb_y[k] * sld0_84[k]
                   - f_4 * pc_y[k] * sld1_84[k];

        t_115[k] = f_12 * slp_58[k]
                   + f_3 * pc_x[k] * smp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_x, pc_y, sld0_89, slp_43, \
                         slp_44, slp_59, sld1_89, sms0_19, sms1_19, smp_58, \
                         smp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_12 * slp_59[k]
                   + f_3 * pc_x[k] * smp_59[k];

        t_117[k] = f_6 * slp_43[k]
                   + f_1 * sms0_19[k]
                   - f_2 * sms1_19[k]
                   + f_3 * pc_y[k] * smp_58[k];

        t_118[k] = f_6 * slp_44[k]
                   + f_3 * pc_y[k] * smp_59[k];

        t_119[k] = pb_y[k] * sld0_89[k]
                   - f_4 * pc_y[k] * sld1_89[k];
    }
}

static auto
compute_prim_smd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sld0,
                                                          const size_t slp, const size_t sld1,
                                                          const size_t sms0, const size_t sms1,
                                                          const size_t smp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

    auto *t_120 = buffer.data(target + 120);
    auto *t_121 = buffer.data(target + 121);
    auto *t_122 = buffer.data(target + 122);
    auto *t_123 = buffer.data(target + 123);
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);
    auto *t_136 = buffer.data(target + 136);
    auto *t_137 = buffer.data(target + 137);
    auto *t_138 = buffer.data(target + 138);
    auto *t_139 = buffer.data(target + 139);
    auto *t_140 = buffer.data(target + 140);
    auto *t_141 = buffer.data(target + 141);
    auto *t_142 = buffer.data(target + 142);
    auto *t_143 = buffer.data(target + 143);
    auto *t_144 = buffer.data(target + 144);
    auto *t_145 = buffer.data(target + 145);
    auto *t_146 = buffer.data(target + 146);
    auto *t_147 = buffer.data(target + 147);
    auto *t_148 = buffer.data(target + 148);
    auto *t_149 = buffer.data(target + 149);
    auto *t_150 = buffer.data(target + 150);
    auto *t_151 = buffer.data(target + 151);
    auto *t_152 = buffer.data(target + 152);
    auto *t_153 = buffer.data(target + 153);
    auto *t_154 = buffer.data(target + 154);
    auto *t_155 = buffer.data(target + 155);
    auto *t_156 = buffer.data(target + 156);
    auto *t_157 = buffer.data(target + 157);
    auto *t_158 = buffer.data(target + 158);
    auto *t_159 = buffer.data(target + 159);
    auto *t_160 = buffer.data(target + 160);
    auto *t_161 = buffer.data(target + 161);
    auto *t_162 = buffer.data(target + 162);
    auto *t_163 = buffer.data(target + 163);
    auto *t_164 = buffer.data(target + 164);
    auto *t_165 = buffer.data(target + 165);
    auto *t_166 = buffer.data(target + 166);
    auto *t_167 = buffer.data(target + 167);
    auto *t_168 = buffer.data(target + 168);
    auto *t_169 = buffer.data(target + 169);
    auto *t_170 = buffer.data(target + 170);
    auto *t_171 = buffer.data(target + 171);
    auto *t_172 = buffer.data(target + 172);
    auto *t_173 = buffer.data(target + 173);
    auto *t_174 = buffer.data(target + 174);
    auto *t_175 = buffer.data(target + 175);
    auto *t_176 = buffer.data(target + 176);
    auto *t_177 = buffer.data(target + 177);
    auto *t_178 = buffer.data(target + 178);
    auto *t_179 = buffer.data(target + 179);
    auto *t_180 = buffer.data(target + 180);
    auto *t_181 = buffer.data(target + 181);
    auto *t_182 = buffer.data(target + 182);
    auto *t_183 = buffer.data(target + 183);
    auto *t_184 = buffer.data(target + 184);
    auto *t_185 = buffer.data(target + 185);
    auto *t_186 = buffer.data(target + 186);
    auto *t_187 = buffer.data(target + 187);
    auto *t_188 = buffer.data(target + 188);
    auto *t_189 = buffer.data(target + 189);
    auto *t_190 = buffer.data(target + 190);
    auto *t_191 = buffer.data(target + 191);
    auto *t_192 = buffer.data(target + 192);
    auto *t_193 = buffer.data(target + 193);
    auto *t_194 = buffer.data(target + 194);
    auto *t_195 = buffer.data(target + 195);
    auto *t_196 = buffer.data(target + 196);
    auto *t_197 = buffer.data(target + 197);
    auto *t_198 = buffer.data(target + 198);
    auto *t_199 = buffer.data(target + 199);
    auto *t_200 = buffer.data(target + 200);
    auto *t_201 = buffer.data(target + 201);
    auto *t_202 = buffer.data(target + 202);
    auto *t_203 = buffer.data(target + 203);
    auto *t_204 = buffer.data(target + 204);
    auto *t_205 = buffer.data(target + 205);
    auto *t_206 = buffer.data(target + 206);
    auto *t_207 = buffer.data(target + 207);
    auto *t_208 = buffer.data(target + 208);
    auto *t_209 = buffer.data(target + 209);
    auto *t_210 = buffer.data(target + 210);
    auto *t_211 = buffer.data(target + 211);
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);
    auto *t_217 = buffer.data(target + 217);
    auto *t_218 = buffer.data(target + 218);
    auto *t_219 = buffer.data(target + 219);
    auto *t_220 = buffer.data(target + 220);
    auto *t_221 = buffer.data(target + 221);
    auto *t_222 = buffer.data(target + 222);
    auto *t_223 = buffer.data(target + 223);
    auto *t_224 = buffer.data(target + 224);
    auto *t_225 = buffer.data(target + 225);
    auto *t_226 = buffer.data(target + 226);
    auto *t_227 = buffer.data(target + 227);
    auto *t_228 = buffer.data(target + 228);
    auto *t_229 = buffer.data(target + 229);
    auto *t_230 = buffer.data(target + 230);
    auto *t_231 = buffer.data(target + 231);
    auto *t_232 = buffer.data(target + 232);
    auto *t_233 = buffer.data(target + 233);
    auto *t_234 = buffer.data(target + 234);
    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sld0_90 = buffer.data(sld0 + 90);
    const auto *sld0_93 = buffer.data(sld0 + 93);
    const auto *sld0_120 = buffer.data(sld0 + 120);
    const auto *sld0_125 = buffer.data(sld0 + 125);
    const auto *sld0_126 = buffer.data(sld0 + 126);
    const auto *sld0_129 = buffer.data(sld0 + 129);
    const auto *sld0_162 = buffer.data(sld0 + 162);
    const auto *sld0_167 = buffer.data(sld0 + 167);
    const auto *sld0_168 = buffer.data(sld0 + 168);
    const auto *sld0_216 = buffer.data(sld0 + 216);
    const auto *sld0_219 = buffer.data(sld0 + 219);
    const auto *sld0_221 = buffer.data(sld0 + 221);
    const auto *sld0_225 = buffer.data(sld0 + 225);
    const auto *sld0_227 = buffer.data(sld0 + 227);
    const auto *sld0_228 = buffer.data(sld0 + 228);
    const auto *sld0_231 = buffer.data(sld0 + 231);
    const auto *sld0_233 = buffer.data(sld0 + 233);
    const auto *sld0_234 = buffer.data(sld0 + 234);
    const auto *sld0_237 = buffer.data(sld0 + 237);

    const auto *slp_44 = buffer.data(slp + 44);
    const auto *slp_46 = buffer.data(slp + 46);
    const auto *slp_47 = buffer.data(slp + 47);
    const auto *slp_50 = buffer.data(slp + 50);
    const auto *slp_52 = buffer.data(slp + 52);
    const auto *slp_53 = buffer.data(slp + 53);
    const auto *slp_55 = buffer.data(slp + 55);
    const auto *slp_56 = buffer.data(slp + 56);
    const auto *slp_58 = buffer.data(slp + 58);
    const auto *slp_59 = buffer.data(slp + 59);
    const auto *slp_60 = buffer.data(slp + 60);
    const auto *slp_61 = buffer.data(slp + 61);
    const auto *slp_62 = buffer.data(slp + 62);
    const auto *slp_63 = buffer.data(slp + 63);
    const auto *slp_64 = buffer.data(slp + 64);
    const auto *slp_65 = buffer.data(slp + 65);
    const auto *slp_67 = buffer.data(slp + 67);
    const auto *slp_68 = buffer.data(slp + 68);
    const auto *slp_69 = buffer.data(slp + 69);
    const auto *slp_70 = buffer.data(slp + 70);
    const auto *slp_71 = buffer.data(slp + 71);
    const auto *slp_72 = buffer.data(slp + 72);
    const auto *slp_73 = buffer.data(slp + 73);
    const auto *slp_74 = buffer.data(slp + 74);
    const auto *slp_75 = buffer.data(slp + 75);
    const auto *slp_76 = buffer.data(slp + 76);
    const auto *slp_77 = buffer.data(slp + 77);
    const auto *slp_79 = buffer.data(slp + 79);
    const auto *slp_80 = buffer.data(slp + 80);
    const auto *slp_81 = buffer.data(slp + 81);
    const auto *slp_82 = buffer.data(slp + 82);
    const auto *slp_83 = buffer.data(slp + 83);
    const auto *slp_84 = buffer.data(slp + 84);
    const auto *slp_85 = buffer.data(slp + 85);
    const auto *slp_86 = buffer.data(slp + 86);
    const auto *slp_88 = buffer.data(slp + 88);
    const auto *slp_89 = buffer.data(slp + 89);
    const auto *slp_90 = buffer.data(slp + 90);
    const auto *slp_91 = buffer.data(slp + 91);
    const auto *slp_92 = buffer.data(slp + 92);
    const auto *slp_93 = buffer.data(slp + 93);
    const auto *slp_94 = buffer.data(slp + 94);
    const auto *slp_95 = buffer.data(slp + 95);
    const auto *slp_96 = buffer.data(slp + 96);
    const auto *slp_97 = buffer.data(slp + 97);
    const auto *slp_98 = buffer.data(slp + 98);
    const auto *slp_99 = buffer.data(slp + 99);
    const auto *slp_100 = buffer.data(slp + 100);
    const auto *slp_101 = buffer.data(slp + 101);
    const auto *slp_103 = buffer.data(slp + 103);
    const auto *slp_104 = buffer.data(slp + 104);
    const auto *slp_105 = buffer.data(slp + 105);
    const auto *slp_106 = buffer.data(slp + 106);
    const auto *slp_107 = buffer.data(slp + 107);
    const auto *slp_108 = buffer.data(slp + 108);
    const auto *slp_109 = buffer.data(slp + 109);
    const auto *slp_110 = buffer.data(slp + 110);
    const auto *slp_112 = buffer.data(slp + 112);
    const auto *slp_113 = buffer.data(slp + 113);
    const auto *slp_114 = buffer.data(slp + 114);
    const auto *slp_115 = buffer.data(slp + 115);
    const auto *slp_116 = buffer.data(slp + 116);
    const auto *slp_117 = buffer.data(slp + 117);
    const auto *slp_118 = buffer.data(slp + 118);
    const auto *slp_119 = buffer.data(slp + 119);

    const auto *sld1_90 = buffer.data(sld1 + 90);
    const auto *sld1_93 = buffer.data(sld1 + 93);
    const auto *sld1_120 = buffer.data(sld1 + 120);
    const auto *sld1_125 = buffer.data(sld1 + 125);
    const auto *sld1_126 = buffer.data(sld1 + 126);
    const auto *sld1_129 = buffer.data(sld1 + 129);
    const auto *sld1_162 = buffer.data(sld1 + 162);
    const auto *sld1_167 = buffer.data(sld1 + 167);
    const auto *sld1_168 = buffer.data(sld1 + 168);
    const auto *sld1_216 = buffer.data(sld1 + 216);
    const auto *sld1_219 = buffer.data(sld1 + 219);
    const auto *sld1_221 = buffer.data(sld1 + 221);
    const auto *sld1_225 = buffer.data(sld1 + 225);
    const auto *sld1_227 = buffer.data(sld1 + 227);
    const auto *sld1_228 = buffer.data(sld1 + 228);
    const auto *sld1_231 = buffer.data(sld1 + 231);
    const auto *sld1_233 = buffer.data(sld1 + 233);
    const auto *sld1_234 = buffer.data(sld1 + 234);
    const auto *sld1_237 = buffer.data(sld1 + 237);

    const auto *sms0_20 = buffer.data(sms0 + 20);
    const auto *sms0_21 = buffer.data(sms0 + 21);
    const auto *sms0_22 = buffer.data(sms0 + 22);
    const auto *sms0_23 = buffer.data(sms0 + 23);
    const auto *sms0_24 = buffer.data(sms0 + 24);
    const auto *sms0_25 = buffer.data(sms0 + 25);
    const auto *sms0_26 = buffer.data(sms0 + 26);
    const auto *sms0_27 = buffer.data(sms0 + 27);
    const auto *sms0_28 = buffer.data(sms0 + 28);
    const auto *sms0_29 = buffer.data(sms0 + 29);
    const auto *sms0_30 = buffer.data(sms0 + 30);
    const auto *sms0_31 = buffer.data(sms0 + 31);
    const auto *sms0_32 = buffer.data(sms0 + 32);
    const auto *sms0_33 = buffer.data(sms0 + 33);
    const auto *sms0_34 = buffer.data(sms0 + 34);
    const auto *sms0_35 = buffer.data(sms0 + 35);

    const auto *sms1_20 = buffer.data(sms1 + 20);
    const auto *sms1_21 = buffer.data(sms1 + 21);
    const auto *sms1_22 = buffer.data(sms1 + 22);
    const auto *sms1_23 = buffer.data(sms1 + 23);
    const auto *sms1_24 = buffer.data(sms1 + 24);
    const auto *sms1_25 = buffer.data(sms1 + 25);
    const auto *sms1_26 = buffer.data(sms1 + 26);
    const auto *sms1_27 = buffer.data(sms1 + 27);
    const auto *sms1_28 = buffer.data(sms1 + 28);
    const auto *sms1_29 = buffer.data(sms1 + 29);
    const auto *sms1_30 = buffer.data(sms1 + 30);
    const auto *sms1_31 = buffer.data(sms1 + 31);
    const auto *sms1_32 = buffer.data(sms1 + 32);
    const auto *sms1_33 = buffer.data(sms1 + 33);
    const auto *sms1_34 = buffer.data(sms1 + 34);
    const auto *sms1_35 = buffer.data(sms1 + 35);

    const auto *smp_60 = buffer.data(smp + 60);
    const auto *smp_61 = buffer.data(smp + 61);
    const auto *smp_62 = buffer.data(smp + 62);
    const auto *smp_63 = buffer.data(smp + 63);
    const auto *smp_64 = buffer.data(smp + 64);
    const auto *smp_65 = buffer.data(smp + 65);
    const auto *smp_67 = buffer.data(smp + 67);
    const auto *smp_68 = buffer.data(smp + 68);
    const auto *smp_69 = buffer.data(smp + 69);
    const auto *smp_70 = buffer.data(smp + 70);
    const auto *smp_71 = buffer.data(smp + 71);
    const auto *smp_72 = buffer.data(smp + 72);
    const auto *smp_73 = buffer.data(smp + 73);
    const auto *smp_74 = buffer.data(smp + 74);
    const auto *smp_75 = buffer.data(smp + 75);
    const auto *smp_76 = buffer.data(smp + 76);
    const auto *smp_77 = buffer.data(smp + 77);
    const auto *smp_79 = buffer.data(smp + 79);
    const auto *smp_80 = buffer.data(smp + 80);
    const auto *smp_81 = buffer.data(smp + 81);
    const auto *smp_82 = buffer.data(smp + 82);
    const auto *smp_83 = buffer.data(smp + 83);
    const auto *smp_84 = buffer.data(smp + 84);
    const auto *smp_85 = buffer.data(smp + 85);
    const auto *smp_86 = buffer.data(smp + 86);
    const auto *smp_88 = buffer.data(smp + 88);
    const auto *smp_89 = buffer.data(smp + 89);
    const auto *smp_90 = buffer.data(smp + 90);
    const auto *smp_91 = buffer.data(smp + 91);
    const auto *smp_92 = buffer.data(smp + 92);
    const auto *smp_93 = buffer.data(smp + 93);
    const auto *smp_94 = buffer.data(smp + 94);
    const auto *smp_95 = buffer.data(smp + 95);
    const auto *smp_96 = buffer.data(smp + 96);
    const auto *smp_97 = buffer.data(smp + 97);
    const auto *smp_98 = buffer.data(smp + 98);
    const auto *smp_99 = buffer.data(smp + 99);
    const auto *smp_100 = buffer.data(smp + 100);
    const auto *smp_101 = buffer.data(smp + 101);
    const auto *smp_103 = buffer.data(smp + 103);
    const auto *smp_104 = buffer.data(smp + 104);
    const auto *smp_105 = buffer.data(smp + 105);
    const auto *smp_106 = buffer.data(smp + 106);
    const auto *smp_107 = buffer.data(smp + 107);
    const auto *smp_109 = buffer.data(smp + 109);
    const auto *smp_110 = buffer.data(smp + 110);
    const auto *smp_112 = buffer.data(smp + 112);
    const auto *smp_113 = buffer.data(smp + 113);
    const auto *smp_115 = buffer.data(smp + 115);
    const auto *smp_116 = buffer.data(smp + 116);
    const auto *smp_118 = buffer.data(smp + 118);
    const auto *smp_119 = buffer.data(smp + 119);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, slp_60, slp_61, \
                         slp_62, sms0_20, sms1_20, smp_60, smp_61, \
                         smp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_12 * slp_60[k]
                   + f_1 * sms0_20[k]
                   - f_2 * sms1_20[k]
                   + f_3 * pc_x[k] * smp_60[k];

        t_121[k] = f_12 * slp_61[k]
                   + f_3 * pc_x[k] * smp_61[k];

        t_122[k] = f_12 * slp_62[k]
                   + f_3 * pc_x[k] * smp_62[k];

        t_123[k] = f_1 * sms0_20[k]
                   - f_2 * sms1_20[k]
                   + f_3 * pc_y[k] * smp_61[k];

        t_124[k] = f_3 * pc_y[k] * smp_62[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_z, slp_44, slp_63, slp_64, sms0_20, \
                         sms0_21, sms1_20, sms1_21, smp_62, smp_63, \
                         smp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_11 * slp_44[k]
                   + f_1 * sms0_20[k]
                   - f_2 * sms1_20[k]
                   + f_3 * pc_z[k] * smp_62[k];

        t_126[k] = f_10 * slp_63[k]
                   + f_1 * sms0_21[k]
                   - f_2 * sms1_21[k]
                   + f_3 * pc_x[k] * smp_63[k];

        t_127[k] = f_10 * slp_64[k]
                   + f_3 * pc_x[k] * smp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, slp_46, slp_47, slp_65, \
                         sms0_21, sms1_21, smp_64, smp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_10 * slp_65[k]
                   + f_3 * pc_x[k] * smp_65[k];

        t_129[k] = f_9 * slp_46[k]
                   + f_1 * sms0_21[k]
                   - f_2 * sms1_21[k]
                   + f_3 * pc_y[k] * smp_64[k];

        t_130[k] = f_9 * slp_47[k]
                   + f_3 * pc_y[k] * smp_65[k];

        t_131[k] = f_1 * sms0_21[k]
                   - f_2 * sms1_21[k]
                   + f_3 * pc_z[k] * smp_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_z, pc_x, pc_z, sld0_90, sld0_93, \
                         slp_67, slp_68, sld1_90, sld1_93, smp_67, \
                         smp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * sld0_90[k]
                   - f_4 * pc_z[k] * sld1_90[k];

        t_133[k] = f_10 * slp_67[k]
                   + f_3 * pc_x[k] * smp_67[k];

        t_134[k] = f_10 * slp_68[k]
                   + f_3 * pc_x[k] * smp_68[k];

        t_135[k] = pb_z[k] * sld0_93[k]
                   - f_4 * pc_z[k] * sld1_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_x, pc_y, pc_z, slp_47, slp_50, slp_69, \
                         sms0_22, sms0_23, sms1_22, sms1_23, smp_68, \
                         smp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * slp_50[k]
                   + f_3 * pc_y[k] * smp_68[k];

        t_137[k] = f_6 * slp_47[k]
                   + f_1 * sms0_22[k]
                   - f_2 * sms1_22[k]
                   + f_3 * pc_z[k] * smp_68[k];

        t_138[k] = f_10 * slp_69[k]
                   + f_1 * sms0_23[k]
                   - f_2 * sms1_23[k]
                   + f_3 * pc_x[k] * smp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, slp_52, slp_53, slp_70, \
                         slp_71, sms0_23, sms1_23, smp_70, smp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_10 * slp_70[k]
                   + f_3 * pc_x[k] * smp_70[k];

        t_140[k] = f_10 * slp_71[k]
                   + f_3 * pc_x[k] * smp_71[k];

        t_141[k] = f_12 * slp_52[k]
                   + f_1 * sms0_23[k]
                   - f_2 * sms1_23[k]
                   + f_3 * pc_y[k] * smp_70[k];

        t_142[k] = f_12 * slp_53[k]
                   + f_3 * pc_y[k] * smp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_z, slp_50, slp_72, slp_73, sms0_23, \
                         sms0_24, sms1_23, sms1_24, smp_71, smp_72, \
                         smp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * slp_50[k]
                   + f_1 * sms0_23[k]
                   - f_2 * sms1_23[k]
                   + f_3 * pc_z[k] * smp_71[k];

        t_144[k] = f_10 * slp_72[k]
                   + f_1 * sms0_24[k]
                   - f_2 * sms1_24[k]
                   + f_3 * pc_x[k] * smp_72[k];

        t_145[k] = f_10 * slp_73[k]
                   + f_3 * pc_x[k] * smp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, slp_53, slp_55, slp_56, \
                         slp_74, sms0_24, sms1_24, smp_73, smp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_10 * slp_74[k]
                   + f_3 * pc_x[k] * smp_74[k];

        t_147[k] = f_10 * slp_55[k]
                   + f_1 * sms0_24[k]
                   - f_2 * sms1_24[k]
                   + f_3 * pc_y[k] * smp_73[k];

        t_148[k] = f_10 * slp_56[k]
                   + f_3 * pc_y[k] * smp_74[k];

        t_149[k] = f_10 * slp_53[k]
                   + f_1 * sms0_24[k]
                   - f_2 * sms1_24[k]
                   + f_3 * pc_z[k] * smp_74[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, slp_58, slp_75, slp_76, \
                         slp_77, sms0_25, sms1_25, smp_75, smp_76, \
                         smp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_10 * slp_75[k]
                   + f_1 * sms0_25[k]
                   - f_2 * sms1_25[k]
                   + f_3 * pc_x[k] * smp_75[k];

        t_151[k] = f_10 * slp_76[k]
                   + f_3 * pc_x[k] * smp_76[k];

        t_152[k] = f_10 * slp_77[k]
                   + f_3 * pc_x[k] * smp_77[k];

        t_153[k] = f_8 * slp_58[k]
                   + f_1 * sms0_25[k]
                   - f_2 * sms1_25[k]
                   + f_3 * pc_y[k] * smp_76[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pc_y, pc_z, sld0_120, slp_56, slp_59, \
                         sld1_120, sms0_25, sms1_25, smp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_8 * slp_59[k]
                   + f_3 * pc_y[k] * smp_77[k];

        t_155[k] = f_12 * slp_56[k]
                   + f_1 * sms0_25[k]
                   - f_2 * sms1_25[k]
                   + f_3 * pc_z[k] * smp_77[k];

        t_156[k] = pb_y[k] * sld0_120[k]
                   - f_4 * pc_y[k] * sld1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, slp_61, slp_62, slp_79, \
                         slp_80, sms0_26, sms1_26, smp_79, smp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_10 * slp_79[k]
                   + f_3 * pc_x[k] * smp_79[k];

        t_158[k] = f_10 * slp_80[k]
                   + f_3 * pc_x[k] * smp_80[k];

        t_159[k] = f_6 * slp_61[k]
                   + f_1 * sms0_26[k]
                   - f_2 * sms1_26[k]
                   + f_3 * pc_y[k] * smp_79[k];

        t_160[k] = f_6 * slp_62[k]
                   + f_3 * pc_y[k] * smp_80[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_y, pc_x, pc_y, sld0_125, slp_81, slp_82, \
                         sld1_125, sms0_27, sms1_27, smp_81, smp_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * sld0_125[k]
                   - f_4 * pc_y[k] * sld1_125[k];

        t_162[k] = f_10 * slp_81[k]
                   + f_1 * sms0_27[k]
                   - f_2 * sms1_27[k]
                   + f_3 * pc_x[k] * smp_81[k];

        t_163[k] = f_10 * slp_82[k]
                   + f_3 * pc_x[k] * smp_82[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, slp_62, slp_83, \
                         sms0_27, sms1_27, smp_82, smp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * slp_83[k]
                   + f_3 * pc_x[k] * smp_83[k];

        t_165[k] = f_1 * sms0_27[k]
                   - f_2 * sms1_27[k]
                   + f_3 * pc_y[k] * smp_82[k];

        t_166[k] = f_3 * pc_y[k] * smp_83[k];

        t_167[k] = f_9 * slp_62[k]
                   + f_1 * sms0_27[k]
                   - f_2 * sms1_27[k]
                   + f_3 * pc_z[k] * smp_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, slp_64, slp_84, slp_85, \
                         slp_86, sms0_28, sms1_28, smp_84, smp_85, \
                         smp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_8 * slp_84[k]
                   + f_1 * sms0_28[k]
                   - f_2 * sms1_28[k]
                   + f_3 * pc_x[k] * smp_84[k];

        t_169[k] = f_8 * slp_85[k]
                   + f_3 * pc_x[k] * smp_85[k];

        t_170[k] = f_8 * slp_86[k]
                   + f_3 * pc_x[k] * smp_86[k];

        t_171[k] = f_7 * slp_64[k]
                   + f_1 * sms0_28[k]
                   - f_2 * sms1_28[k]
                   + f_3 * pc_y[k] * smp_85[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, pc_x, pc_y, pc_z, sld0_126, slp_65, \
                         slp_88, sld1_126, sms0_28, sms1_28, smp_86, \
                         smp_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_7 * slp_65[k]
                   + f_3 * pc_y[k] * smp_86[k];

        t_173[k] = f_1 * sms0_28[k]
                   - f_2 * sms1_28[k]
                   + f_3 * pc_z[k] * smp_86[k];

        t_174[k] = pb_z[k] * sld0_126[k]
                   - f_4 * pc_z[k] * sld1_126[k];

        t_175[k] = f_8 * slp_88[k]
                   + f_3 * pc_x[k] * smp_88[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_z, pc_x, pc_y, pc_z, sld0_129, slp_65, \
                         slp_68, slp_89, sld1_129, sms0_29, sms1_29, \
                         smp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_8 * slp_89[k]
                   + f_3 * pc_x[k] * smp_89[k];

        t_177[k] = pb_z[k] * sld0_129[k]
                   - f_4 * pc_z[k] * sld1_129[k];

        t_178[k] = f_9 * slp_68[k]
                   + f_3 * pc_y[k] * smp_89[k];

        t_179[k] = f_6 * slp_65[k]
                   + f_1 * sms0_29[k]
                   - f_2 * sms1_29[k]
                   + f_3 * pc_z[k] * smp_89[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, slp_70, slp_90, slp_91, \
                         slp_92, sms0_30, sms1_30, smp_90, smp_91, \
                         smp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_8 * slp_90[k]
                   + f_1 * sms0_30[k]
                   - f_2 * sms1_30[k]
                   + f_3 * pc_x[k] * smp_90[k];

        t_181[k] = f_8 * slp_91[k]
                   + f_3 * pc_x[k] * smp_91[k];

        t_182[k] = f_8 * slp_92[k]
                   + f_3 * pc_x[k] * smp_92[k];

        t_183[k] = f_11 * slp_70[k]
                   + f_1 * sms0_30[k]
                   - f_2 * sms1_30[k]
                   + f_3 * pc_y[k] * smp_91[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, pc_z, slp_68, slp_71, slp_93, \
                         sms0_30, sms0_31, sms1_30, sms1_31, smp_92, \
                         smp_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_11 * slp_71[k]
                   + f_3 * pc_y[k] * smp_92[k];

        t_185[k] = f_8 * slp_68[k]
                   + f_1 * sms0_30[k]
                   - f_2 * sms1_30[k]
                   + f_3 * pc_z[k] * smp_92[k];

        t_186[k] = f_8 * slp_93[k]
                   + f_1 * sms0_31[k]
                   - f_2 * sms1_31[k]
                   + f_3 * pc_x[k] * smp_93[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pc_x, pc_y, slp_73, slp_74, slp_94, \
                         slp_95, sms0_31, sms1_31, smp_94, smp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_8 * slp_94[k]
                   + f_3 * pc_x[k] * smp_94[k];

        t_188[k] = f_8 * slp_95[k]
                   + f_3 * pc_x[k] * smp_95[k];

        t_189[k] = f_12 * slp_73[k]
                   + f_1 * sms0_31[k]
                   - f_2 * sms1_31[k]
                   + f_3 * pc_y[k] * smp_94[k];

        t_190[k] = f_12 * slp_74[k]
                   + f_3 * pc_y[k] * smp_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_x, pc_z, slp_71, slp_96, slp_97, sms0_31, \
                         sms0_32, sms1_31, sms1_32, smp_95, smp_96, \
                         smp_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_10 * slp_71[k]
                   + f_1 * sms0_31[k]
                   - f_2 * sms1_31[k]
                   + f_3 * pc_z[k] * smp_95[k];

        t_192[k] = f_8 * slp_96[k]
                   + f_1 * sms0_32[k]
                   - f_2 * sms1_32[k]
                   + f_3 * pc_x[k] * smp_96[k];

        t_193[k] = f_8 * slp_97[k]
                   + f_3 * pc_x[k] * smp_97[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, slp_74, slp_76, slp_77, \
                         slp_98, sms0_32, sms1_32, smp_97, smp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_8 * slp_98[k]
                   + f_3 * pc_x[k] * smp_98[k];

        t_195[k] = f_10 * slp_76[k]
                   + f_1 * sms0_32[k]
                   - f_2 * sms1_32[k]
                   + f_3 * pc_y[k] * smp_97[k];

        t_196[k] = f_10 * slp_77[k]
                   + f_3 * pc_y[k] * smp_98[k];

        t_197[k] = f_12 * slp_74[k]
                   + f_1 * sms0_32[k]
                   - f_2 * sms1_32[k]
                   + f_3 * pc_z[k] * smp_98[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, slp_79, slp_99, slp_100, \
                         slp_101, sms0_33, sms1_33, smp_99, smp_100, \
                         smp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_8 * slp_99[k]
                   + f_1 * sms0_33[k]
                   - f_2 * sms1_33[k]
                   + f_3 * pc_x[k] * smp_99[k];

        t_199[k] = f_8 * slp_100[k]
                   + f_3 * pc_x[k] * smp_100[k];

        t_200[k] = f_8 * slp_101[k]
                   + f_3 * pc_x[k] * smp_101[k];

        t_201[k] = f_8 * slp_79[k]
                   + f_1 * sms0_33[k]
                   - f_2 * sms1_33[k]
                   + f_3 * pc_y[k] * smp_100[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_y, pc_y, pc_z, sld0_162, slp_77, slp_80, \
                         sld1_162, sms0_33, sms1_33, smp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_8 * slp_80[k]
                   + f_3 * pc_y[k] * smp_101[k];

        t_203[k] = f_11 * slp_77[k]
                   + f_1 * sms0_33[k]
                   - f_2 * sms1_33[k]
                   + f_3 * pc_z[k] * smp_101[k];

        t_204[k] = pb_y[k] * sld0_162[k]
                   - f_4 * pc_y[k] * sld1_162[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, slp_82, slp_83, slp_103, \
                         slp_104, sms0_34, sms1_34, smp_103, smp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_8 * slp_103[k]
                   + f_3 * pc_x[k] * smp_103[k];

        t_206[k] = f_8 * slp_104[k]
                   + f_3 * pc_x[k] * smp_104[k];

        t_207[k] = f_6 * slp_82[k]
                   + f_1 * sms0_34[k]
                   - f_2 * sms1_34[k]
                   + f_3 * pc_y[k] * smp_103[k];

        t_208[k] = f_6 * slp_83[k]
                   + f_3 * pc_y[k] * smp_104[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_y, pc_x, pc_y, sld0_167, slp_105, slp_106, \
                         sld1_167, sms0_35, sms1_35, smp_105, smp_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pb_y[k] * sld0_167[k]
                   - f_4 * pc_y[k] * sld1_167[k];

        t_210[k] = f_8 * slp_105[k]
                   + f_1 * sms0_35[k]
                   - f_2 * sms1_35[k]
                   + f_3 * pc_x[k] * smp_105[k];

        t_211[k] = f_8 * slp_106[k]
                   + f_3 * pc_x[k] * smp_106[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, slp_83, slp_107, \
                         sms0_35, sms1_35, smp_106, smp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * slp_107[k]
                   + f_3 * pc_x[k] * smp_107[k];

        t_213[k] = f_1 * sms0_35[k]
                   - f_2 * sms1_35[k]
                   + f_3 * pc_y[k] * smp_106[k];

        t_214[k] = f_3 * pc_y[k] * smp_107[k];

        t_215[k] = f_7 * slp_83[k]
                   + f_1 * sms0_35[k]
                   - f_2 * sms1_35[k]
                   + f_3 * pc_z[k] * smp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pb_x, pc_x, sld0_216, sld0_219, slp_108, \
                         slp_109, slp_110, sld1_216, sld1_219, smp_109, \
                         smp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pb_x[k] * sld0_216[k]
                   + f_8 * slp_108[k]
                   - f_4 * pc_x[k] * sld1_216[k];

        t_217[k] = f_6 * slp_109[k]
                   + f_3 * pc_x[k] * smp_109[k];

        t_218[k] = f_6 * slp_110[k]
                   + f_3 * pc_x[k] * smp_110[k];

        t_219[k] = pb_x[k] * sld0_219[k]
                   - f_4 * pc_x[k] * sld1_219[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pb_x, pb_z, pc_x, pc_y, pc_z, sld0_168, \
                         sld0_221, slp_86, sld1_168, sld1_221, \
                         smp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_5 * slp_86[k]
                   + f_3 * pc_y[k] * smp_110[k];

        t_221[k] = pb_x[k] * sld0_221[k]
                   - f_4 * pc_x[k] * sld1_221[k];

        t_222[k] = pb_z[k] * sld0_168[k]
                   - f_4 * pc_z[k] * sld1_168[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_x, pc_x, pc_y, sld0_225, slp_89, \
                         slp_112, slp_113, sld1_225, smp_112, smp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_6 * slp_112[k]
                   + f_3 * pc_x[k] * smp_112[k];

        t_224[k] = f_6 * slp_113[k]
                   + f_3 * pc_x[k] * smp_113[k];

        t_225[k] = pb_x[k] * sld0_225[k]
                   - f_4 * pc_x[k] * sld1_225[k];

        t_226[k] = f_7 * slp_89[k]
                   + f_3 * pc_y[k] * smp_113[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, t_230, pb_x, pc_x, sld0_227, sld0_228, slp_114, \
                         slp_115, slp_116, sld1_227, sld1_228, smp_115, \
                         smp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = pb_x[k] * sld0_227[k]
                   - f_4 * pc_x[k] * sld1_227[k];

        t_228[k] = pb_x[k] * sld0_228[k]
                   + f_8 * slp_114[k]
                   - f_4 * pc_x[k] * sld1_228[k];

        t_229[k] = f_6 * slp_115[k]
                   + f_3 * pc_x[k] * smp_115[k];

        t_230[k] = f_6 * slp_116[k]
                   + f_3 * pc_x[k] * smp_116[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, pc_x, pc_y, sld0_231, sld0_233, \
                         sld0_234, slp_92, slp_117, sld1_231, sld1_233, sld1_234, \
                         smp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pb_x[k] * sld0_231[k]
                   - f_4 * pc_x[k] * sld1_231[k];

        t_232[k] = f_9 * slp_92[k]
                   + f_3 * pc_y[k] * smp_116[k];

        t_233[k] = pb_x[k] * sld0_233[k]
                   - f_4 * pc_x[k] * sld1_233[k];

        t_234[k] = pb_x[k] * sld0_234[k]
                   + f_8 * slp_117[k]
                   - f_4 * pc_x[k] * sld1_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pb_x, pc_x, pc_y, sld0_237, slp_95, \
                         slp_118, slp_119, sld1_237, smp_118, smp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_6 * slp_118[k]
                   + f_3 * pc_x[k] * smp_118[k];

        t_236[k] = f_6 * slp_119[k]
                   + f_3 * pc_x[k] * smp_119[k];

        t_237[k] = pb_x[k] * sld0_237[k]
                   - f_4 * pc_x[k] * sld1_237[k];

        t_238[k] = f_11 * slp_95[k]
                   + f_3 * pc_y[k] * smp_119[k];
    }
}

static auto
compute_prim_smd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sld0,
                                                          const size_t slp, const size_t sld1,
                                                          const size_t sms0, const size_t sms1,
                                                          const size_t smp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

    auto *t_239 = buffer.data(target + 239);
    auto *t_240 = buffer.data(target + 240);
    auto *t_241 = buffer.data(target + 241);
    auto *t_242 = buffer.data(target + 242);
    auto *t_243 = buffer.data(target + 243);
    auto *t_244 = buffer.data(target + 244);
    auto *t_245 = buffer.data(target + 245);
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);
    auto *t_256 = buffer.data(target + 256);
    auto *t_257 = buffer.data(target + 257);
    auto *t_258 = buffer.data(target + 258);
    auto *t_259 = buffer.data(target + 259);
    auto *t_260 = buffer.data(target + 260);
    auto *t_261 = buffer.data(target + 261);
    auto *t_262 = buffer.data(target + 262);
    auto *t_263 = buffer.data(target + 263);
    auto *t_264 = buffer.data(target + 264);
    auto *t_265 = buffer.data(target + 265);
    auto *t_266 = buffer.data(target + 266);
    auto *t_267 = buffer.data(target + 267);
    auto *t_268 = buffer.data(target + 268);
    auto *t_269 = buffer.data(target + 269);
    auto *t_270 = buffer.data(target + 270);
    auto *t_271 = buffer.data(target + 271);
    auto *t_272 = buffer.data(target + 272);
    auto *t_273 = buffer.data(target + 273);
    auto *t_274 = buffer.data(target + 274);
    auto *t_275 = buffer.data(target + 275);
    auto *t_276 = buffer.data(target + 276);
    auto *t_277 = buffer.data(target + 277);
    auto *t_278 = buffer.data(target + 278);
    auto *t_279 = buffer.data(target + 279);
    auto *t_280 = buffer.data(target + 280);
    auto *t_281 = buffer.data(target + 281);
    auto *t_282 = buffer.data(target + 282);
    auto *t_283 = buffer.data(target + 283);
    auto *t_284 = buffer.data(target + 284);
    auto *t_285 = buffer.data(target + 285);
    auto *t_286 = buffer.data(target + 286);
    auto *t_287 = buffer.data(target + 287);
    auto *t_288 = buffer.data(target + 288);
    auto *t_289 = buffer.data(target + 289);
    auto *t_290 = buffer.data(target + 290);
    auto *t_291 = buffer.data(target + 291);
    auto *t_292 = buffer.data(target + 292);
    auto *t_293 = buffer.data(target + 293);
    auto *t_294 = buffer.data(target + 294);
    auto *t_295 = buffer.data(target + 295);
    auto *t_296 = buffer.data(target + 296);
    auto *t_297 = buffer.data(target + 297);
    auto *t_298 = buffer.data(target + 298);
    auto *t_299 = buffer.data(target + 299);
    auto *t_300 = buffer.data(target + 300);
    auto *t_301 = buffer.data(target + 301);
    auto *t_302 = buffer.data(target + 302);
    auto *t_303 = buffer.data(target + 303);
    auto *t_304 = buffer.data(target + 304);
    auto *t_305 = buffer.data(target + 305);
    auto *t_306 = buffer.data(target + 306);
    auto *t_307 = buffer.data(target + 307);
    auto *t_308 = buffer.data(target + 308);
    auto *t_309 = buffer.data(target + 309);
    auto *t_310 = buffer.data(target + 310);
    auto *t_311 = buffer.data(target + 311);
    auto *t_312 = buffer.data(target + 312);
    auto *t_313 = buffer.data(target + 313);
    auto *t_314 = buffer.data(target + 314);
    auto *t_315 = buffer.data(target + 315);
    auto *t_316 = buffer.data(target + 316);
    auto *t_317 = buffer.data(target + 317);
    auto *t_318 = buffer.data(target + 318);
    auto *t_319 = buffer.data(target + 319);
    auto *t_320 = buffer.data(target + 320);
    auto *t_321 = buffer.data(target + 321);
    auto *t_322 = buffer.data(target + 322);
    auto *t_323 = buffer.data(target + 323);
    auto *t_324 = buffer.data(target + 324);
    auto *t_325 = buffer.data(target + 325);
    auto *t_326 = buffer.data(target + 326);
    auto *t_327 = buffer.data(target + 327);
    auto *t_328 = buffer.data(target + 328);
    auto *t_329 = buffer.data(target + 329);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sld0_210 = buffer.data(sld0 + 210);
    const auto *sld0_216 = buffer.data(sld0 + 216);
    const auto *sld0_219 = buffer.data(sld0 + 219);
    const auto *sld0_239 = buffer.data(sld0 + 239);
    const auto *sld0_240 = buffer.data(sld0 + 240);
    const auto *sld0_243 = buffer.data(sld0 + 243);
    const auto *sld0_245 = buffer.data(sld0 + 245);
    const auto *sld0_246 = buffer.data(sld0 + 246);
    const auto *sld0_249 = buffer.data(sld0 + 249);
    const auto *sld0_251 = buffer.data(sld0 + 251);
    const auto *sld0_252 = buffer.data(sld0 + 252);
    const auto *sld0_255 = buffer.data(sld0 + 255);
    const auto *sld0_257 = buffer.data(sld0 + 257);
    const auto *sld0_261 = buffer.data(sld0 + 261);
    const auto *sld0_263 = buffer.data(sld0 + 263);
    const auto *sld0_264 = buffer.data(sld0 + 264);
    const auto *sld0_267 = buffer.data(sld0 + 267);
    const auto *sld0_269 = buffer.data(sld0 + 269);

    const auto *slp_98 = buffer.data(slp + 98);
    const auto *slp_101 = buffer.data(slp + 101);
    const auto *slp_104 = buffer.data(slp + 104);
    const auto *slp_107 = buffer.data(slp + 107);
    const auto *slp_109 = buffer.data(slp + 109);
    const auto *slp_110 = buffer.data(slp + 110);
    const auto *slp_113 = buffer.data(slp + 113);
    const auto *slp_115 = buffer.data(slp + 115);
    const auto *slp_116 = buffer.data(slp + 116);
    const auto *slp_118 = buffer.data(slp + 118);
    const auto *slp_119 = buffer.data(slp + 119);
    const auto *slp_120 = buffer.data(slp + 120);
    const auto *slp_121 = buffer.data(slp + 121);
    const auto *slp_122 = buffer.data(slp + 122);
    const auto *slp_123 = buffer.data(slp + 123);
    const auto *slp_124 = buffer.data(slp + 124);
    const auto *slp_125 = buffer.data(slp + 125);
    const auto *slp_126 = buffer.data(slp + 126);
    const auto *slp_127 = buffer.data(slp + 127);
    const auto *slp_128 = buffer.data(slp + 128);
    const auto *slp_130 = buffer.data(slp + 130);
    const auto *slp_131 = buffer.data(slp + 131);
    const auto *slp_132 = buffer.data(slp + 132);
    const auto *slp_133 = buffer.data(slp + 133);
    const auto *slp_134 = buffer.data(slp + 134);

    const auto *sld1_210 = buffer.data(sld1 + 210);
    const auto *sld1_216 = buffer.data(sld1 + 216);
    const auto *sld1_219 = buffer.data(sld1 + 219);
    const auto *sld1_239 = buffer.data(sld1 + 239);
    const auto *sld1_240 = buffer.data(sld1 + 240);
    const auto *sld1_243 = buffer.data(sld1 + 243);
    const auto *sld1_245 = buffer.data(sld1 + 245);
    const auto *sld1_246 = buffer.data(sld1 + 246);
    const auto *sld1_249 = buffer.data(sld1 + 249);
    const auto *sld1_251 = buffer.data(sld1 + 251);
    const auto *sld1_252 = buffer.data(sld1 + 252);
    const auto *sld1_255 = buffer.data(sld1 + 255);
    const auto *sld1_257 = buffer.data(sld1 + 257);
    const auto *sld1_261 = buffer.data(sld1 + 261);
    const auto *sld1_263 = buffer.data(sld1 + 263);
    const auto *sld1_264 = buffer.data(sld1 + 264);
    const auto *sld1_267 = buffer.data(sld1 + 267);
    const auto *sld1_269 = buffer.data(sld1 + 269);

    const auto *sms0_45 = buffer.data(sms0 + 45);
    const auto *sms0_46 = buffer.data(sms0 + 46);
    const auto *sms0_47 = buffer.data(sms0 + 47);
    const auto *sms0_48 = buffer.data(sms0 + 48);
    const auto *sms0_49 = buffer.data(sms0 + 49);
    const auto *sms0_50 = buffer.data(sms0 + 50);
    const auto *sms0_51 = buffer.data(sms0 + 51);
    const auto *sms0_52 = buffer.data(sms0 + 52);
    const auto *sms0_54 = buffer.data(sms0 + 54);

    const auto *sms1_45 = buffer.data(sms1 + 45);
    const auto *sms1_46 = buffer.data(sms1 + 46);
    const auto *sms1_47 = buffer.data(sms1 + 47);
    const auto *sms1_48 = buffer.data(sms1 + 48);
    const auto *sms1_49 = buffer.data(sms1 + 49);
    const auto *sms1_50 = buffer.data(sms1 + 50);
    const auto *sms1_51 = buffer.data(sms1 + 51);
    const auto *sms1_52 = buffer.data(sms1 + 52);
    const auto *sms1_54 = buffer.data(sms1 + 54);

    const auto *smp_121 = buffer.data(smp + 121);
    const auto *smp_122 = buffer.data(smp + 122);
    const auto *smp_124 = buffer.data(smp + 124);
    const auto *smp_125 = buffer.data(smp + 125);
    const auto *smp_127 = buffer.data(smp + 127);
    const auto *smp_128 = buffer.data(smp + 128);
    const auto *smp_130 = buffer.data(smp + 130);
    const auto *smp_131 = buffer.data(smp + 131);
    const auto *smp_133 = buffer.data(smp + 133);
    const auto *smp_134 = buffer.data(smp + 134);
    const auto *smp_135 = buffer.data(smp + 135);
    const auto *smp_136 = buffer.data(smp + 136);
    const auto *smp_137 = buffer.data(smp + 137);
    const auto *smp_139 = buffer.data(smp + 139);
    const auto *smp_140 = buffer.data(smp + 140);
    const auto *smp_141 = buffer.data(smp + 141);
    const auto *smp_142 = buffer.data(smp + 142);
    const auto *smp_143 = buffer.data(smp + 143);
    const auto *smp_144 = buffer.data(smp + 144);
    const auto *smp_145 = buffer.data(smp + 145);
    const auto *smp_146 = buffer.data(smp + 146);
    const auto *smp_147 = buffer.data(smp + 147);
    const auto *smp_148 = buffer.data(smp + 148);
    const auto *smp_149 = buffer.data(smp + 149);
    const auto *smp_150 = buffer.data(smp + 150);
    const auto *smp_151 = buffer.data(smp + 151);
    const auto *smp_152 = buffer.data(smp + 152);
    const auto *smp_153 = buffer.data(smp + 153);
    const auto *smp_154 = buffer.data(smp + 154);
    const auto *smp_155 = buffer.data(smp + 155);
    const auto *smp_156 = buffer.data(smp + 156);
    const auto *smp_157 = buffer.data(smp + 157);
    const auto *smp_158 = buffer.data(smp + 158);
    const auto *smp_160 = buffer.data(smp + 160);
    const auto *smp_161 = buffer.data(smp + 161);
    const auto *smp_162 = buffer.data(smp + 162);
    const auto *smp_163 = buffer.data(smp + 163);
    const auto *smp_164 = buffer.data(smp + 164);

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pb_x, pc_x, sld0_239, sld0_240, slp_120, \
                         slp_121, slp_122, sld1_239, sld1_240, smp_121, \
                         smp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = pb_x[k] * sld0_239[k]
                   - f_4 * pc_x[k] * sld1_239[k];

        t_240[k] = pb_x[k] * sld0_240[k]
                   + f_8 * slp_120[k]
                   - f_4 * pc_x[k] * sld1_240[k];

        t_241[k] = f_6 * slp_121[k]
                   + f_3 * pc_x[k] * smp_121[k];

        t_242[k] = f_6 * slp_122[k]
                   + f_3 * pc_x[k] * smp_122[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_x, pc_x, pc_y, sld0_243, sld0_245, \
                         sld0_246, slp_98, slp_123, sld1_243, sld1_245, sld1_246, \
                         smp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pb_x[k] * sld0_243[k]
                   - f_4 * pc_x[k] * sld1_243[k];

        t_244[k] = f_12 * slp_98[k]
                   + f_3 * pc_y[k] * smp_122[k];

        t_245[k] = pb_x[k] * sld0_245[k]
                   - f_4 * pc_x[k] * sld1_245[k];

        t_246[k] = pb_x[k] * sld0_246[k]
                   + f_8 * slp_123[k]
                   - f_4 * pc_x[k] * sld1_246[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_x, pc_x, pc_y, sld0_249, slp_101, \
                         slp_124, slp_125, sld1_249, smp_124, smp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_6 * slp_124[k]
                   + f_3 * pc_x[k] * smp_124[k];

        t_248[k] = f_6 * slp_125[k]
                   + f_3 * pc_x[k] * smp_125[k];

        t_249[k] = pb_x[k] * sld0_249[k]
                   - f_4 * pc_x[k] * sld1_249[k];

        t_250[k] = f_10 * slp_101[k]
                   + f_3 * pc_y[k] * smp_125[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_x, pc_x, sld0_251, sld0_252, slp_126, \
                         slp_127, slp_128, sld1_251, sld1_252, smp_127, \
                         smp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_x[k] * sld0_251[k]
                   - f_4 * pc_x[k] * sld1_251[k];

        t_252[k] = pb_x[k] * sld0_252[k]
                   + f_8 * slp_126[k]
                   - f_4 * pc_x[k] * sld1_252[k];

        t_253[k] = f_6 * slp_127[k]
                   + f_3 * pc_x[k] * smp_127[k];

        t_254[k] = f_6 * slp_128[k]
                   + f_3 * pc_x[k] * smp_128[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pb_x, pb_y, pc_x, pc_y, sld0_210, \
                         sld0_255, sld0_257, slp_104, sld1_210, sld1_255, sld1_257, \
                         smp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pb_x[k] * sld0_255[k]
                   - f_4 * pc_x[k] * sld1_255[k];

        t_256[k] = f_8 * slp_104[k]
                   + f_3 * pc_y[k] * smp_128[k];

        t_257[k] = pb_x[k] * sld0_257[k]
                   - f_4 * pc_x[k] * sld1_257[k];

        t_258[k] = pb_y[k] * sld0_210[k]
                   - f_4 * pc_y[k] * sld1_210[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_x, pc_x, pc_y, sld0_261, slp_107, \
                         slp_130, slp_131, sld1_261, smp_130, smp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_6 * slp_130[k]
                   + f_3 * pc_x[k] * smp_130[k];

        t_260[k] = f_6 * slp_131[k]
                   + f_3 * pc_x[k] * smp_131[k];

        t_261[k] = pb_x[k] * sld0_261[k]
                   - f_4 * pc_x[k] * sld1_261[k];

        t_262[k] = f_6 * slp_107[k]
                   + f_3 * pc_y[k] * smp_131[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pb_x, pc_x, sld0_263, sld0_264, slp_132, \
                         slp_133, slp_134, sld1_263, sld1_264, smp_133, \
                         smp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = pb_x[k] * sld0_263[k]
                   - f_4 * pc_x[k] * sld1_263[k];

        t_264[k] = pb_x[k] * sld0_264[k]
                   + f_8 * slp_132[k]
                   - f_4 * pc_x[k] * sld1_264[k];

        t_265[k] = f_6 * slp_133[k]
                   + f_3 * pc_x[k] * smp_133[k];

        t_266[k] = f_6 * slp_134[k]
                   + f_3 * pc_x[k] * smp_134[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pb_x, pc_x, pc_y, sld0_267, sld0_269, \
                         sld1_267, sld1_269, sms0_45, sms1_45, smp_134, \
                         smp_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pb_x[k] * sld0_267[k]
                   - f_4 * pc_x[k] * sld1_267[k];

        t_268[k] = f_3 * pc_y[k] * smp_134[k];

        t_269[k] = pb_x[k] * sld0_269[k]
                   - f_4 * pc_x[k] * sld1_269[k];

        t_270[k] = f_1 * sms0_45[k]
                   - f_2 * sms1_45[k]
                   + f_3 * pc_x[k] * smp_135[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, slp_109, \
                         slp_110, sms0_45, sms1_45, smp_136, smp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_3 * pc_x[k] * smp_136[k];

        t_272[k] = f_3 * pc_x[k] * smp_137[k];

        t_273[k] = f_0 * slp_109[k]
                   + f_1 * sms0_45[k]
                   - f_2 * sms1_45[k]
                   + f_3 * pc_y[k] * smp_136[k];

        t_274[k] = f_0 * slp_110[k]
                   + f_3 * pc_y[k] * smp_137[k];

        t_275[k] = f_1 * sms0_45[k]
                   - f_2 * sms1_45[k]
                   + f_3 * pc_z[k] * smp_137[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pb_z, pc_x, pc_y, pc_z, sld0_216, \
                         sld0_219, slp_113, sld1_216, sld1_219, smp_139, \
                         smp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_z[k] * sld0_216[k]
                   - f_4 * pc_z[k] * sld1_216[k];

        t_277[k] = f_3 * pc_x[k] * smp_139[k];

        t_278[k] = f_3 * pc_x[k] * smp_140[k];

        t_279[k] = pb_z[k] * sld0_219[k]
                   - f_4 * pc_z[k] * sld1_219[k];

        t_280[k] = f_5 * slp_113[k]
                   + f_3 * pc_y[k] * smp_140[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pc_x, pc_z, slp_110, sms0_46, sms0_47, \
                         sms1_46, sms1_47, smp_140, smp_141, smp_142, \
                         smp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_6 * slp_110[k]
                   + f_1 * sms0_46[k]
                   - f_2 * sms1_46[k]
                   + f_3 * pc_z[k] * smp_140[k];

        t_282[k] = f_1 * sms0_47[k]
                   - f_2 * sms1_47[k]
                   + f_3 * pc_x[k] * smp_141[k];

        t_283[k] = f_3 * pc_x[k] * smp_142[k];

        t_284[k] = f_3 * pc_x[k] * smp_143[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_y, pc_z, slp_113, slp_115, slp_116, sms0_47, \
                         sms1_47, smp_142, smp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_7 * slp_115[k]
                   + f_1 * sms0_47[k]
                   - f_2 * sms1_47[k]
                   + f_3 * pc_y[k] * smp_142[k];

        t_286[k] = f_7 * slp_116[k]
                   + f_3 * pc_y[k] * smp_143[k];

        t_287[k] = f_8 * slp_113[k]
                   + f_1 * sms0_47[k]
                   - f_2 * sms1_47[k]
                   + f_3 * pc_z[k] * smp_143[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pc_x, pc_y, slp_118, slp_119, \
                         sms0_48, sms1_48, smp_144, smp_145, smp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_1 * sms0_48[k]
                   - f_2 * sms1_48[k]
                   + f_3 * pc_x[k] * smp_144[k];

        t_289[k] = f_3 * pc_x[k] * smp_145[k];

        t_290[k] = f_3 * pc_x[k] * smp_146[k];

        t_291[k] = f_9 * slp_118[k]
                   + f_1 * sms0_48[k]
                   - f_2 * sms1_48[k]
                   + f_3 * pc_y[k] * smp_145[k];

        t_292[k] = f_9 * slp_119[k]
                   + f_3 * pc_y[k] * smp_146[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_z, slp_116, sms0_48, sms0_49, \
                         sms1_48, sms1_49, smp_146, smp_147, smp_148, \
                         smp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_10 * slp_116[k]
                   + f_1 * sms0_48[k]
                   - f_2 * sms1_48[k]
                   + f_3 * pc_z[k] * smp_146[k];

        t_294[k] = f_1 * sms0_49[k]
                   - f_2 * sms1_49[k]
                   + f_3 * pc_x[k] * smp_147[k];

        t_295[k] = f_3 * pc_x[k] * smp_148[k];

        t_296[k] = f_3 * pc_x[k] * smp_149[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, slp_119, slp_121, slp_122, sms0_49, \
                         sms1_49, smp_148, smp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_11 * slp_121[k]
                   + f_1 * sms0_49[k]
                   - f_2 * sms1_49[k]
                   + f_3 * pc_y[k] * smp_148[k];

        t_298[k] = f_11 * slp_122[k]
                   + f_3 * pc_y[k] * smp_149[k];

        t_299[k] = f_12 * slp_119[k]
                   + f_1 * sms0_49[k]
                   - f_2 * sms1_49[k]
                   + f_3 * pc_z[k] * smp_149[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, pc_y, slp_124, slp_125, \
                         sms0_50, sms1_50, smp_150, smp_151, smp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * sms0_50[k]
                   - f_2 * sms1_50[k]
                   + f_3 * pc_x[k] * smp_150[k];

        t_301[k] = f_3 * pc_x[k] * smp_151[k];

        t_302[k] = f_3 * pc_x[k] * smp_152[k];

        t_303[k] = f_12 * slp_124[k]
                   + f_1 * sms0_50[k]
                   - f_2 * sms1_50[k]
                   + f_3 * pc_y[k] * smp_151[k];

        t_304[k] = f_12 * slp_125[k]
                   + f_3 * pc_y[k] * smp_152[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_z, slp_122, sms0_50, sms0_51, \
                         sms1_50, sms1_51, smp_152, smp_153, smp_154, \
                         smp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_11 * slp_122[k]
                   + f_1 * sms0_50[k]
                   - f_2 * sms1_50[k]
                   + f_3 * pc_z[k] * smp_152[k];

        t_306[k] = f_1 * sms0_51[k]
                   - f_2 * sms1_51[k]
                   + f_3 * pc_x[k] * smp_153[k];

        t_307[k] = f_3 * pc_x[k] * smp_154[k];

        t_308[k] = f_3 * pc_x[k] * smp_155[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_y, pc_z, slp_125, slp_127, slp_128, sms0_51, \
                         sms1_51, smp_154, smp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_10 * slp_127[k]
                   + f_1 * sms0_51[k]
                   - f_2 * sms1_51[k]
                   + f_3 * pc_y[k] * smp_154[k];

        t_310[k] = f_10 * slp_128[k]
                   + f_3 * pc_y[k] * smp_155[k];

        t_311[k] = f_9 * slp_125[k]
                   + f_1 * sms0_51[k]
                   - f_2 * sms1_51[k]
                   + f_3 * pc_z[k] * smp_155[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, pc_x, pc_y, slp_130, slp_131, \
                         sms0_52, sms1_52, smp_156, smp_157, smp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_1 * sms0_52[k]
                   - f_2 * sms1_52[k]
                   + f_3 * pc_x[k] * smp_156[k];

        t_313[k] = f_3 * pc_x[k] * smp_157[k];

        t_314[k] = f_3 * pc_x[k] * smp_158[k];

        t_315[k] = f_8 * slp_130[k]
                   + f_1 * sms0_52[k]
                   - f_2 * sms1_52[k]
                   + f_3 * pc_y[k] * smp_157[k];

        t_316[k] = f_8 * slp_131[k]
                   + f_3 * pc_y[k] * smp_158[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pb_y, pc_x, pc_y, pc_z, sld0_264, \
                         slp_128, sld1_264, sms0_52, sms1_52, smp_158, smp_160, \
                         smp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_7 * slp_128[k]
                   + f_1 * sms0_52[k]
                   - f_2 * sms1_52[k]
                   + f_3 * pc_z[k] * smp_158[k];

        t_318[k] = pb_y[k] * sld0_264[k]
                   - f_4 * pc_y[k] * sld1_264[k];

        t_319[k] = f_3 * pc_x[k] * smp_160[k];

        t_320[k] = f_3 * pc_x[k] * smp_161[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pb_y, pc_y, sld0_267, sld0_269, slp_133, \
                         slp_134, sld1_267, sld1_269, smp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = pb_y[k] * sld0_267[k]
                   + f_8 * slp_133[k]
                   - f_4 * pc_y[k] * sld1_267[k];

        t_322[k] = f_6 * slp_134[k]
                   + f_3 * pc_y[k] * smp_161[k];

        t_323[k] = pb_y[k] * sld0_269[k]
                   - f_4 * pc_y[k] * sld1_269[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pc_x, pc_y, pc_z, slp_134, \
                         sms0_54, sms1_54, smp_162, smp_163, smp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * sms0_54[k]
                   - f_2 * sms1_54[k]
                   + f_3 * pc_x[k] * smp_162[k];

        t_325[k] = f_3 * pc_x[k] * smp_163[k];

        t_326[k] = f_3 * pc_x[k] * smp_164[k];

        t_327[k] = f_1 * sms0_54[k]
                   - f_2 * sms1_54[k]
                   + f_3 * pc_y[k] * smp_163[k];

        t_328[k] = f_3 * pc_y[k] * smp_164[k];

        t_329[k] = f_0 * slp_134[k]
                   + f_1 * sms0_54[k]
                   - f_2 * sms1_54[k]
                   + f_3 * pc_z[k] * smp_164[k];
    }
}

auto
compute_prim_smd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sld0, const size_t slp,
                                                   const size_t sld1, const size_t sms0,
                                                   const size_t sms1, const size_t smp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_smd_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sld0, slp,
                                                              sld1, sms0, sms1, smp, ncols,
                                                              gamma, p, q);

    compute_prim_smd_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sld0, slp,
                                                              sld1, sms0, sms1, smp, ncols,
                                                              gamma, p, q);

    compute_prim_smd_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sld0, slp,
                                                              sld1, sms0, sms1, smp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
