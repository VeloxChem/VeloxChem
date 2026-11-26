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


#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_snd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smd0,
                                                          const size_t smp, const size_t smd1,
                                                          const size_t sns0, const size_t sns1,
                                                          const size_t snp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *smd0_0 = buffer.data(smd0 + 0);
    const auto *smd0_3 = buffer.data(smd0 + 3);
    const auto *smd0_5 = buffer.data(smd0 + 5);
    const auto *smd0_9 = buffer.data(smd0 + 9);
    const auto *smd0_12 = buffer.data(smd0 + 12);
    const auto *smd0_17 = buffer.data(smd0 + 17);
    const auto *smd0_18 = buffer.data(smd0 + 18);
    const auto *smd0_21 = buffer.data(smd0 + 21);
    const auto *smd0_30 = buffer.data(smd0 + 30);
    const auto *smd0_35 = buffer.data(smd0 + 35);
    const auto *smd0_36 = buffer.data(smd0 + 36);
    const auto *smd0_39 = buffer.data(smd0 + 39);
    const auto *smd0_54 = buffer.data(smd0 + 54);
    const auto *smd0_59 = buffer.data(smd0 + 59);
    const auto *smd0_60 = buffer.data(smd0 + 60);
    const auto *smd0_63 = buffer.data(smd0 + 63);
    const auto *smd0_84 = buffer.data(smd0 + 84);
    const auto *smd0_89 = buffer.data(smd0 + 89);

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

    const auto *smd1_0 = buffer.data(smd1 + 0);
    const auto *smd1_3 = buffer.data(smd1 + 3);
    const auto *smd1_5 = buffer.data(smd1 + 5);
    const auto *smd1_9 = buffer.data(smd1 + 9);
    const auto *smd1_12 = buffer.data(smd1 + 12);
    const auto *smd1_17 = buffer.data(smd1 + 17);
    const auto *smd1_18 = buffer.data(smd1 + 18);
    const auto *smd1_21 = buffer.data(smd1 + 21);
    const auto *smd1_30 = buffer.data(smd1 + 30);
    const auto *smd1_35 = buffer.data(smd1 + 35);
    const auto *smd1_36 = buffer.data(smd1 + 36);
    const auto *smd1_39 = buffer.data(smd1 + 39);
    const auto *smd1_54 = buffer.data(smd1 + 54);
    const auto *smd1_59 = buffer.data(smd1 + 59);
    const auto *smd1_60 = buffer.data(smd1 + 60);
    const auto *smd1_63 = buffer.data(smd1 + 63);
    const auto *smd1_84 = buffer.data(smd1 + 84);
    const auto *smd1_89 = buffer.data(smd1 + 89);

    const auto *sns0_0 = buffer.data(sns0 + 0);
    const auto *sns0_1 = buffer.data(sns0 + 1);
    const auto *sns0_2 = buffer.data(sns0 + 2);
    const auto *sns0_3 = buffer.data(sns0 + 3);
    const auto *sns0_5 = buffer.data(sns0 + 5);
    const auto *sns0_6 = buffer.data(sns0 + 6);
    const auto *sns0_7 = buffer.data(sns0 + 7);
    const auto *sns0_8 = buffer.data(sns0 + 8);
    const auto *sns0_9 = buffer.data(sns0 + 9);
    const auto *sns0_10 = buffer.data(sns0 + 10);
    const auto *sns0_11 = buffer.data(sns0 + 11);
    const auto *sns0_12 = buffer.data(sns0 + 12);
    const auto *sns0_13 = buffer.data(sns0 + 13);
    const auto *sns0_14 = buffer.data(sns0 + 14);
    const auto *sns0_15 = buffer.data(sns0 + 15);
    const auto *sns0_16 = buffer.data(sns0 + 16);
    const auto *sns0_17 = buffer.data(sns0 + 17);
    const auto *sns0_18 = buffer.data(sns0 + 18);
    const auto *sns0_19 = buffer.data(sns0 + 19);

    const auto *sns1_0 = buffer.data(sns1 + 0);
    const auto *sns1_1 = buffer.data(sns1 + 1);
    const auto *sns1_2 = buffer.data(sns1 + 2);
    const auto *sns1_3 = buffer.data(sns1 + 3);
    const auto *sns1_5 = buffer.data(sns1 + 5);
    const auto *sns1_6 = buffer.data(sns1 + 6);
    const auto *sns1_7 = buffer.data(sns1 + 7);
    const auto *sns1_8 = buffer.data(sns1 + 8);
    const auto *sns1_9 = buffer.data(sns1 + 9);
    const auto *sns1_10 = buffer.data(sns1 + 10);
    const auto *sns1_11 = buffer.data(sns1 + 11);
    const auto *sns1_12 = buffer.data(sns1 + 12);
    const auto *sns1_13 = buffer.data(sns1 + 13);
    const auto *sns1_14 = buffer.data(sns1 + 14);
    const auto *sns1_15 = buffer.data(sns1 + 15);
    const auto *sns1_16 = buffer.data(sns1 + 16);
    const auto *sns1_17 = buffer.data(sns1 + 17);
    const auto *sns1_18 = buffer.data(sns1 + 18);
    const auto *sns1_19 = buffer.data(sns1 + 19);

    const auto *snp_0 = buffer.data(snp + 0);
    const auto *snp_1 = buffer.data(snp + 1);
    const auto *snp_2 = buffer.data(snp + 2);
    const auto *snp_4 = buffer.data(snp + 4);
    const auto *snp_5 = buffer.data(snp + 5);
    const auto *snp_7 = buffer.data(snp + 7);
    const auto *snp_8 = buffer.data(snp + 8);
    const auto *snp_9 = buffer.data(snp + 9);
    const auto *snp_10 = buffer.data(snp + 10);
    const auto *snp_11 = buffer.data(snp + 11);
    const auto *snp_13 = buffer.data(snp + 13);
    const auto *snp_14 = buffer.data(snp + 14);
    const auto *snp_15 = buffer.data(snp + 15);
    const auto *snp_16 = buffer.data(snp + 16);
    const auto *snp_17 = buffer.data(snp + 17);
    const auto *snp_18 = buffer.data(snp + 18);
    const auto *snp_19 = buffer.data(snp + 19);
    const auto *snp_20 = buffer.data(snp + 20);
    const auto *snp_22 = buffer.data(snp + 22);
    const auto *snp_23 = buffer.data(snp + 23);
    const auto *snp_25 = buffer.data(snp + 25);
    const auto *snp_26 = buffer.data(snp + 26);
    const auto *snp_27 = buffer.data(snp + 27);
    const auto *snp_28 = buffer.data(snp + 28);
    const auto *snp_29 = buffer.data(snp + 29);
    const auto *snp_30 = buffer.data(snp + 30);
    const auto *snp_31 = buffer.data(snp + 31);
    const auto *snp_32 = buffer.data(snp + 32);
    const auto *snp_34 = buffer.data(snp + 34);
    const auto *snp_35 = buffer.data(snp + 35);
    const auto *snp_36 = buffer.data(snp + 36);
    const auto *snp_37 = buffer.data(snp + 37);
    const auto *snp_38 = buffer.data(snp + 38);
    const auto *snp_40 = buffer.data(snp + 40);
    const auto *snp_41 = buffer.data(snp + 41);
    const auto *snp_42 = buffer.data(snp + 42);
    const auto *snp_43 = buffer.data(snp + 43);
    const auto *snp_44 = buffer.data(snp + 44);
    const auto *snp_45 = buffer.data(snp + 45);
    const auto *snp_46 = buffer.data(snp + 46);
    const auto *snp_47 = buffer.data(snp + 47);
    const auto *snp_49 = buffer.data(snp + 49);
    const auto *snp_50 = buffer.data(snp + 50);
    const auto *snp_51 = buffer.data(snp + 51);
    const auto *snp_52 = buffer.data(snp + 52);
    const auto *snp_53 = buffer.data(snp + 53);
    const auto *snp_54 = buffer.data(snp + 54);
    const auto *snp_55 = buffer.data(snp + 55);
    const auto *snp_56 = buffer.data(snp + 56);
    const auto *snp_58 = buffer.data(snp + 58);
    const auto *snp_59 = buffer.data(snp + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, smp_0, smp_1, smp_2, sns0_0, \
                         sns1_0, snp_0, snp_1, snp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smp_0[k]
                 + f_1 * sns0_0[k]
                 - f_2 * sns1_0[k]
                 + f_3 * pc_x[k] * snp_0[k];

        t_1[k] = f_0 * smp_1[k]
                 + f_3 * pc_x[k] * snp_1[k];

        t_2[k] = f_0 * smp_2[k]
                 + f_3 * pc_x[k] * snp_2[k];

        t_3[k] = f_1 * sns0_0[k]
                 - f_2 * sns1_0[k]
                 + f_3 * pc_y[k] * snp_1[k];

        t_4[k] = f_3 * pc_y[k] * snp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, smd0_0, smp_4, smd1_0, sns0_0, \
                         sns1_0, snp_2, snp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sns0_0[k]
                 - f_2 * sns1_0[k]
                 + f_3 * pc_z[k] * snp_2[k];

        t_6[k] = pb_y[k] * smd0_0[k]
                 - f_4 * pc_y[k] * smd1_0[k];

        t_7[k] = f_5 * smp_4[k]
                 + f_3 * pc_x[k] * snp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, smd0_5, smp_1, smp_2, smp_5, \
                         smd1_5, sns0_1, sns1_1, snp_4, snp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * smp_5[k]
                 + f_3 * pc_x[k] * snp_5[k];

        t_9[k] = f_6 * smp_1[k]
                 + f_1 * sns0_1[k]
                 - f_2 * sns1_1[k]
                 + f_3 * pc_y[k] * snp_4[k];

        t_10[k] = f_6 * smp_2[k]
                  + f_3 * pc_y[k] * snp_5[k];

        t_11[k] = pb_y[k] * smd0_5[k]
                  - f_4 * pc_y[k] * smd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, smd0_0, smd0_3, smp_7, \
                         smp_8, smd1_0, smd1_3, snp_7, snp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * smd0_0[k]
                  - f_4 * pc_z[k] * smd1_0[k];

        t_13[k] = f_5 * smp_7[k]
                  + f_3 * pc_x[k] * snp_7[k];

        t_14[k] = f_5 * smp_8[k]
                  + f_3 * pc_x[k] * snp_8[k];

        t_15[k] = pb_z[k] * smd0_3[k]
                  - f_4 * pc_z[k] * smd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, smp_2, smp_9, sns0_2, sns0_3, \
                         sns1_2, sns1_3, snp_8, snp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * snp_8[k];

        t_17[k] = f_6 * smp_2[k]
                  + f_1 * sns0_2[k]
                  - f_2 * sns1_2[k]
                  + f_3 * pc_z[k] * snp_8[k];

        t_18[k] = f_7 * smp_9[k]
                  + f_1 * sns0_3[k]
                  - f_2 * sns1_3[k]
                  + f_3 * pc_x[k] * snp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, smp_4, smp_5, smp_10, \
                         smp_11, sns0_3, sns1_3, snp_10, snp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * smp_10[k]
                  + f_3 * pc_x[k] * snp_10[k];

        t_20[k] = f_7 * smp_11[k]
                  + f_3 * pc_x[k] * snp_11[k];

        t_21[k] = f_8 * smp_4[k]
                  + f_1 * sns0_3[k]
                  - f_2 * sns1_3[k]
                  + f_3 * pc_y[k] * snp_10[k];

        t_22[k] = f_8 * smp_5[k]
                  + f_3 * pc_y[k] * snp_11[k];

        t_23[k] = f_1 * sns0_3[k]
                  - f_2 * sns1_3[k]
                  + f_3 * pc_z[k] * snp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, smd0_12, smp_13, smp_14, smd1_12, \
                         snp_13, snp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * smd0_12[k]
                  - f_4 * pc_y[k] * smd1_12[k];

        t_25[k] = f_7 * smp_13[k]
                  + f_3 * pc_x[k] * snp_13[k];

        t_26[k] = f_7 * smp_14[k]
                  + f_3 * pc_x[k] * snp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, smd0_9, smd0_17, smp_8, \
                         smd1_9, smd1_17, snp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * smd0_9[k]
                  - f_4 * pc_z[k] * smd1_9[k];

        t_28[k] = f_6 * smp_8[k]
                  + f_3 * pc_y[k] * snp_14[k];

        t_29[k] = pb_y[k] * smd0_17[k]
                  - f_4 * pc_y[k] * smd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, smp_15, smp_16, smp_17, \
                         sns0_5, sns1_5, snp_15, snp_16, snp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * smp_15[k]
                  + f_1 * sns0_5[k]
                  - f_2 * sns1_5[k]
                  + f_3 * pc_x[k] * snp_15[k];

        t_31[k] = f_7 * smp_16[k]
                  + f_3 * pc_x[k] * snp_16[k];

        t_32[k] = f_7 * smp_17[k]
                  + f_3 * pc_x[k] * snp_17[k];

        t_33[k] = f_1 * sns0_5[k]
                  - f_2 * sns1_5[k]
                  + f_3 * pc_y[k] * snp_16[k];

        t_34[k] = f_3 * pc_y[k] * snp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, smp_8, smp_18, smp_19, sns0_5, sns0_6, \
                         sns1_5, sns1_6, snp_17, snp_18, snp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * smp_8[k]
                  + f_1 * sns0_5[k]
                  - f_2 * sns1_5[k]
                  + f_3 * pc_z[k] * snp_17[k];

        t_36[k] = f_9 * smp_18[k]
                  + f_1 * sns0_6[k]
                  - f_2 * sns1_6[k]
                  + f_3 * pc_x[k] * snp_18[k];

        t_37[k] = f_9 * smp_19[k]
                  + f_3 * pc_x[k] * snp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, smp_10, smp_11, smp_20, \
                         sns0_6, sns1_6, snp_19, snp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * smp_20[k]
                  + f_3 * pc_x[k] * snp_20[k];

        t_39[k] = f_10 * smp_10[k]
                  + f_1 * sns0_6[k]
                  - f_2 * sns1_6[k]
                  + f_3 * pc_y[k] * snp_19[k];

        t_40[k] = f_10 * smp_11[k]
                  + f_3 * pc_y[k] * snp_20[k];

        t_41[k] = f_1 * sns0_6[k]
                  - f_2 * sns1_6[k]
                  + f_3 * pc_z[k] * snp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, smd0_18, smd0_21, smp_22, \
                         smp_23, smd1_18, smd1_21, snp_22, snp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * smd0_18[k]
                  - f_4 * pc_z[k] * smd1_18[k];

        t_43[k] = f_9 * smp_22[k]
                  + f_3 * pc_x[k] * snp_22[k];

        t_44[k] = f_9 * smp_23[k]
                  + f_3 * pc_x[k] * snp_23[k];

        t_45[k] = pb_z[k] * smd0_21[k]
                  - f_4 * pc_z[k] * smd1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, smd0_30, smp_11, smp_14, smd1_30, \
                         sns0_7, sns1_7, snp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * smp_14[k]
                  + f_3 * pc_y[k] * snp_23[k];

        t_47[k] = f_6 * smp_11[k]
                  + f_1 * sns0_7[k]
                  - f_2 * sns1_7[k]
                  + f_3 * pc_z[k] * snp_23[k];

        t_48[k] = pb_y[k] * smd0_30[k]
                  - f_4 * pc_y[k] * smd1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, smp_16, smp_17, smp_25, smp_26, \
                         sns0_8, sns1_8, snp_25, snp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * smp_25[k]
                  + f_3 * pc_x[k] * snp_25[k];

        t_50[k] = f_9 * smp_26[k]
                  + f_3 * pc_x[k] * snp_26[k];

        t_51[k] = f_6 * smp_16[k]
                  + f_1 * sns0_8[k]
                  - f_2 * sns1_8[k]
                  + f_3 * pc_y[k] * snp_25[k];

        t_52[k] = f_6 * smp_17[k]
                  + f_3 * pc_y[k] * snp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, smd0_35, smp_27, smp_28, smd1_35, \
                         sns0_9, sns1_9, snp_27, snp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * smd0_35[k]
                  - f_4 * pc_y[k] * smd1_35[k];

        t_54[k] = f_9 * smp_27[k]
                  + f_1 * sns0_9[k]
                  - f_2 * sns1_9[k]
                  + f_3 * pc_x[k] * snp_27[k];

        t_55[k] = f_9 * smp_28[k]
                  + f_3 * pc_x[k] * snp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, smp_17, smp_29, sns0_9, \
                         sns1_9, snp_28, snp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * smp_29[k]
                  + f_3 * pc_x[k] * snp_29[k];

        t_57[k] = f_1 * sns0_9[k]
                  - f_2 * sns1_9[k]
                  + f_3 * pc_y[k] * snp_28[k];

        t_58[k] = f_3 * pc_y[k] * snp_29[k];

        t_59[k] = f_10 * smp_17[k]
                  + f_1 * sns0_9[k]
                  - f_2 * sns1_9[k]
                  + f_3 * pc_z[k] * snp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, smp_19, smp_30, smp_31, smp_32, \
                         sns0_10, sns1_10, snp_30, snp_31, snp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * smp_30[k]
                  + f_1 * sns0_10[k]
                  - f_2 * sns1_10[k]
                  + f_3 * pc_x[k] * snp_30[k];

        t_61[k] = f_11 * smp_31[k]
                  + f_3 * pc_x[k] * snp_31[k];

        t_62[k] = f_11 * smp_32[k]
                  + f_3 * pc_x[k] * snp_32[k];

        t_63[k] = f_12 * smp_19[k]
                  + f_1 * sns0_10[k]
                  - f_2 * sns1_10[k]
                  + f_3 * pc_y[k] * snp_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, smd0_36, smp_20, \
                         smp_34, smd1_36, sns0_10, sns1_10, snp_32, \
                         snp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_12 * smp_20[k]
                  + f_3 * pc_y[k] * snp_32[k];

        t_65[k] = f_1 * sns0_10[k]
                  - f_2 * sns1_10[k]
                  + f_3 * pc_z[k] * snp_32[k];

        t_66[k] = pb_z[k] * smd0_36[k]
                  - f_4 * pc_z[k] * smd1_36[k];

        t_67[k] = f_11 * smp_34[k]
                  + f_3 * pc_x[k] * snp_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, smd0_39, smp_20, \
                         smp_23, smp_35, smd1_39, sns0_11, sns1_11, \
                         snp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * smp_35[k]
                  + f_3 * pc_x[k] * snp_35[k];

        t_69[k] = pb_z[k] * smd0_39[k]
                  - f_4 * pc_z[k] * smd1_39[k];

        t_70[k] = f_10 * smp_23[k]
                  + f_3 * pc_y[k] * snp_35[k];

        t_71[k] = f_6 * smp_20[k]
                  + f_1 * sns0_11[k]
                  - f_2 * sns1_11[k]
                  + f_3 * pc_z[k] * snp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, smp_25, smp_36, smp_37, smp_38, \
                         sns0_12, sns1_12, snp_36, snp_37, snp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_11 * smp_36[k]
                  + f_1 * sns0_12[k]
                  - f_2 * sns1_12[k]
                  + f_3 * pc_x[k] * snp_36[k];

        t_73[k] = f_11 * smp_37[k]
                  + f_3 * pc_x[k] * snp_37[k];

        t_74[k] = f_11 * smp_38[k]
                  + f_3 * pc_x[k] * snp_38[k];

        t_75[k] = f_8 * smp_25[k]
                  + f_1 * sns0_12[k]
                  - f_2 * sns1_12[k]
                  + f_3 * pc_y[k] * snp_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, smd0_54, smp_23, smp_26, smd1_54, \
                         sns0_12, sns1_12, snp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * smp_26[k]
                  + f_3 * pc_y[k] * snp_38[k];

        t_77[k] = f_8 * smp_23[k]
                  + f_1 * sns0_12[k]
                  - f_2 * sns1_12[k]
                  + f_3 * pc_z[k] * snp_38[k];

        t_78[k] = pb_y[k] * smd0_54[k]
                  - f_4 * pc_y[k] * smd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, smp_28, smp_29, smp_40, smp_41, \
                         sns0_13, sns1_13, snp_40, snp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_11 * smp_40[k]
                  + f_3 * pc_x[k] * snp_40[k];

        t_80[k] = f_11 * smp_41[k]
                  + f_3 * pc_x[k] * snp_41[k];

        t_81[k] = f_6 * smp_28[k]
                  + f_1 * sns0_13[k]
                  - f_2 * sns1_13[k]
                  + f_3 * pc_y[k] * snp_40[k];

        t_82[k] = f_6 * smp_29[k]
                  + f_3 * pc_y[k] * snp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, smd0_59, smp_42, smp_43, smd1_59, \
                         sns0_14, sns1_14, snp_42, snp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * smd0_59[k]
                  - f_4 * pc_y[k] * smd1_59[k];

        t_84[k] = f_11 * smp_42[k]
                  + f_1 * sns0_14[k]
                  - f_2 * sns1_14[k]
                  + f_3 * pc_x[k] * snp_42[k];

        t_85[k] = f_11 * smp_43[k]
                  + f_3 * pc_x[k] * snp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, smp_29, smp_44, sns0_14, \
                         sns1_14, snp_43, snp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * smp_44[k]
                  + f_3 * pc_x[k] * snp_44[k];

        t_87[k] = f_1 * sns0_14[k]
                  - f_2 * sns1_14[k]
                  + f_3 * pc_y[k] * snp_43[k];

        t_88[k] = f_3 * pc_y[k] * snp_44[k];

        t_89[k] = f_12 * smp_29[k]
                  + f_1 * sns0_14[k]
                  - f_2 * sns1_14[k]
                  + f_3 * pc_z[k] * snp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, smp_31, smp_45, smp_46, smp_47, \
                         sns0_15, sns1_15, snp_45, snp_46, snp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_13 * smp_45[k]
                  + f_1 * sns0_15[k]
                  - f_2 * sns1_15[k]
                  + f_3 * pc_x[k] * snp_45[k];

        t_91[k] = f_13 * smp_46[k]
                  + f_3 * pc_x[k] * snp_46[k];

        t_92[k] = f_13 * smp_47[k]
                  + f_3 * pc_x[k] * snp_47[k];

        t_93[k] = f_13 * smp_31[k]
                  + f_1 * sns0_15[k]
                  - f_2 * sns1_15[k]
                  + f_3 * pc_y[k] * snp_46[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pb_z, pc_x, pc_y, pc_z, smd0_60, smp_32, \
                         smp_49, smd1_60, sns0_15, sns1_15, snp_47, \
                         snp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_13 * smp_32[k]
                  + f_3 * pc_y[k] * snp_47[k];

        t_95[k] = f_1 * sns0_15[k]
                  - f_2 * sns1_15[k]
                  + f_3 * pc_z[k] * snp_47[k];

        t_96[k] = pb_z[k] * smd0_60[k]
                  - f_4 * pc_z[k] * smd1_60[k];

        t_97[k] = f_13 * smp_49[k]
                  + f_3 * pc_x[k] * snp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pb_z, pc_x, pc_y, pc_z, smd0_63, smp_32, \
                         smp_35, smp_50, smd1_63, sns0_16, sns1_16, \
                         snp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_13 * smp_50[k]
                  + f_3 * pc_x[k] * snp_50[k];

        t_99[k] = pb_z[k] * smd0_63[k]
                  - f_4 * pc_z[k] * smd1_63[k];

        t_100[k] = f_12 * smp_35[k]
                   + f_3 * pc_y[k] * snp_50[k];

        t_101[k] = f_6 * smp_32[k]
                   + f_1 * sns0_16[k]
                   - f_2 * sns1_16[k]
                   + f_3 * pc_z[k] * snp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, smp_37, smp_51, smp_52, \
                         smp_53, sns0_17, sns1_17, snp_51, snp_52, \
                         snp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_13 * smp_51[k]
                   + f_1 * sns0_17[k]
                   - f_2 * sns1_17[k]
                   + f_3 * pc_x[k] * snp_51[k];

        t_103[k] = f_13 * smp_52[k]
                   + f_3 * pc_x[k] * snp_52[k];

        t_104[k] = f_13 * smp_53[k]
                   + f_3 * pc_x[k] * snp_53[k];

        t_105[k] = f_10 * smp_37[k]
                   + f_1 * sns0_17[k]
                   - f_2 * sns1_17[k]
                   + f_3 * pc_y[k] * snp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, smp_35, smp_38, smp_54, \
                         sns0_17, sns0_18, sns1_17, sns1_18, snp_53, \
                         snp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * smp_38[k]
                   + f_3 * pc_y[k] * snp_53[k];

        t_107[k] = f_8 * smp_35[k]
                   + f_1 * sns0_17[k]
                   - f_2 * sns1_17[k]
                   + f_3 * pc_z[k] * snp_53[k];

        t_108[k] = f_13 * smp_54[k]
                   + f_1 * sns0_18[k]
                   - f_2 * sns1_18[k]
                   + f_3 * pc_x[k] * snp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, smp_40, smp_41, smp_55, \
                         smp_56, sns0_18, sns1_18, snp_55, snp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_13 * smp_55[k]
                   + f_3 * pc_x[k] * snp_55[k];

        t_110[k] = f_13 * smp_56[k]
                   + f_3 * pc_x[k] * snp_56[k];

        t_111[k] = f_8 * smp_40[k]
                   + f_1 * sns0_18[k]
                   - f_2 * sns1_18[k]
                   + f_3 * pc_y[k] * snp_55[k];

        t_112[k] = f_8 * smp_41[k]
                   + f_3 * pc_y[k] * snp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pb_y, pc_x, pc_y, pc_z, smd0_84, smp_38, smp_58, \
                         smd1_84, sns0_18, sns1_18, snp_56, snp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * smp_38[k]
                   + f_1 * sns0_18[k]
                   - f_2 * sns1_18[k]
                   + f_3 * pc_z[k] * snp_56[k];

        t_114[k] = pb_y[k] * smd0_84[k]
                   - f_4 * pc_y[k] * smd1_84[k];

        t_115[k] = f_13 * smp_58[k]
                   + f_3 * pc_x[k] * snp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_x, pc_y, smd0_89, smp_43, \
                         smp_44, smp_59, smd1_89, sns0_19, sns1_19, snp_58, \
                         snp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_13 * smp_59[k]
                   + f_3 * pc_x[k] * snp_59[k];

        t_117[k] = f_6 * smp_43[k]
                   + f_1 * sns0_19[k]
                   - f_2 * sns1_19[k]
                   + f_3 * pc_y[k] * snp_58[k];

        t_118[k] = f_6 * smp_44[k]
                   + f_3 * pc_y[k] * snp_59[k];

        t_119[k] = pb_y[k] * smd0_89[k]
                   - f_4 * pc_y[k] * smd1_89[k];
    }
}

static auto
compute_prim_snd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smd0,
                                                          const size_t smp, const size_t smd1,
                                                          const size_t sns0, const size_t sns1,
                                                          const size_t snp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smd0_90 = buffer.data(smd0 + 90);
    const auto *smd0_93 = buffer.data(smd0 + 93);
    const auto *smd0_120 = buffer.data(smd0 + 120);
    const auto *smd0_125 = buffer.data(smd0 + 125);
    const auto *smd0_126 = buffer.data(smd0 + 126);
    const auto *smd0_129 = buffer.data(smd0 + 129);
    const auto *smd0_162 = buffer.data(smd0 + 162);
    const auto *smd0_167 = buffer.data(smd0 + 167);
    const auto *smd0_168 = buffer.data(smd0 + 168);
    const auto *smd0_171 = buffer.data(smd0 + 171);

    const auto *smp_44 = buffer.data(smp + 44);
    const auto *smp_46 = buffer.data(smp + 46);
    const auto *smp_47 = buffer.data(smp + 47);
    const auto *smp_50 = buffer.data(smp + 50);
    const auto *smp_52 = buffer.data(smp + 52);
    const auto *smp_53 = buffer.data(smp + 53);
    const auto *smp_55 = buffer.data(smp + 55);
    const auto *smp_56 = buffer.data(smp + 56);
    const auto *smp_58 = buffer.data(smp + 58);
    const auto *smp_59 = buffer.data(smp + 59);
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
    const auto *smp_108 = buffer.data(smp + 108);
    const auto *smp_109 = buffer.data(smp + 109);
    const auto *smp_110 = buffer.data(smp + 110);
    const auto *smp_112 = buffer.data(smp + 112);
    const auto *smp_113 = buffer.data(smp + 113);
    const auto *smp_114 = buffer.data(smp + 114);
    const auto *smp_115 = buffer.data(smp + 115);
    const auto *smp_116 = buffer.data(smp + 116);
    const auto *smp_117 = buffer.data(smp + 117);

    const auto *smd1_90 = buffer.data(smd1 + 90);
    const auto *smd1_93 = buffer.data(smd1 + 93);
    const auto *smd1_120 = buffer.data(smd1 + 120);
    const auto *smd1_125 = buffer.data(smd1 + 125);
    const auto *smd1_126 = buffer.data(smd1 + 126);
    const auto *smd1_129 = buffer.data(smd1 + 129);
    const auto *smd1_162 = buffer.data(smd1 + 162);
    const auto *smd1_167 = buffer.data(smd1 + 167);
    const auto *smd1_168 = buffer.data(smd1 + 168);
    const auto *smd1_171 = buffer.data(smd1 + 171);

    const auto *sns0_20 = buffer.data(sns0 + 20);
    const auto *sns0_21 = buffer.data(sns0 + 21);
    const auto *sns0_22 = buffer.data(sns0 + 22);
    const auto *sns0_23 = buffer.data(sns0 + 23);
    const auto *sns0_24 = buffer.data(sns0 + 24);
    const auto *sns0_25 = buffer.data(sns0 + 25);
    const auto *sns0_26 = buffer.data(sns0 + 26);
    const auto *sns0_27 = buffer.data(sns0 + 27);
    const auto *sns0_28 = buffer.data(sns0 + 28);
    const auto *sns0_29 = buffer.data(sns0 + 29);
    const auto *sns0_30 = buffer.data(sns0 + 30);
    const auto *sns0_31 = buffer.data(sns0 + 31);
    const auto *sns0_32 = buffer.data(sns0 + 32);
    const auto *sns0_33 = buffer.data(sns0 + 33);
    const auto *sns0_34 = buffer.data(sns0 + 34);
    const auto *sns0_35 = buffer.data(sns0 + 35);
    const auto *sns0_36 = buffer.data(sns0 + 36);
    const auto *sns0_37 = buffer.data(sns0 + 37);
    const auto *sns0_38 = buffer.data(sns0 + 38);
    const auto *sns0_39 = buffer.data(sns0 + 39);

    const auto *sns1_20 = buffer.data(sns1 + 20);
    const auto *sns1_21 = buffer.data(sns1 + 21);
    const auto *sns1_22 = buffer.data(sns1 + 22);
    const auto *sns1_23 = buffer.data(sns1 + 23);
    const auto *sns1_24 = buffer.data(sns1 + 24);
    const auto *sns1_25 = buffer.data(sns1 + 25);
    const auto *sns1_26 = buffer.data(sns1 + 26);
    const auto *sns1_27 = buffer.data(sns1 + 27);
    const auto *sns1_28 = buffer.data(sns1 + 28);
    const auto *sns1_29 = buffer.data(sns1 + 29);
    const auto *sns1_30 = buffer.data(sns1 + 30);
    const auto *sns1_31 = buffer.data(sns1 + 31);
    const auto *sns1_32 = buffer.data(sns1 + 32);
    const auto *sns1_33 = buffer.data(sns1 + 33);
    const auto *sns1_34 = buffer.data(sns1 + 34);
    const auto *sns1_35 = buffer.data(sns1 + 35);
    const auto *sns1_36 = buffer.data(sns1 + 36);
    const auto *sns1_37 = buffer.data(sns1 + 37);
    const auto *sns1_38 = buffer.data(sns1 + 38);
    const auto *sns1_39 = buffer.data(sns1 + 39);

    const auto *snp_60 = buffer.data(snp + 60);
    const auto *snp_61 = buffer.data(snp + 61);
    const auto *snp_62 = buffer.data(snp + 62);
    const auto *snp_63 = buffer.data(snp + 63);
    const auto *snp_64 = buffer.data(snp + 64);
    const auto *snp_65 = buffer.data(snp + 65);
    const auto *snp_67 = buffer.data(snp + 67);
    const auto *snp_68 = buffer.data(snp + 68);
    const auto *snp_69 = buffer.data(snp + 69);
    const auto *snp_70 = buffer.data(snp + 70);
    const auto *snp_71 = buffer.data(snp + 71);
    const auto *snp_72 = buffer.data(snp + 72);
    const auto *snp_73 = buffer.data(snp + 73);
    const auto *snp_74 = buffer.data(snp + 74);
    const auto *snp_75 = buffer.data(snp + 75);
    const auto *snp_76 = buffer.data(snp + 76);
    const auto *snp_77 = buffer.data(snp + 77);
    const auto *snp_79 = buffer.data(snp + 79);
    const auto *snp_80 = buffer.data(snp + 80);
    const auto *snp_81 = buffer.data(snp + 81);
    const auto *snp_82 = buffer.data(snp + 82);
    const auto *snp_83 = buffer.data(snp + 83);
    const auto *snp_84 = buffer.data(snp + 84);
    const auto *snp_85 = buffer.data(snp + 85);
    const auto *snp_86 = buffer.data(snp + 86);
    const auto *snp_88 = buffer.data(snp + 88);
    const auto *snp_89 = buffer.data(snp + 89);
    const auto *snp_90 = buffer.data(snp + 90);
    const auto *snp_91 = buffer.data(snp + 91);
    const auto *snp_92 = buffer.data(snp + 92);
    const auto *snp_93 = buffer.data(snp + 93);
    const auto *snp_94 = buffer.data(snp + 94);
    const auto *snp_95 = buffer.data(snp + 95);
    const auto *snp_96 = buffer.data(snp + 96);
    const auto *snp_97 = buffer.data(snp + 97);
    const auto *snp_98 = buffer.data(snp + 98);
    const auto *snp_99 = buffer.data(snp + 99);
    const auto *snp_100 = buffer.data(snp + 100);
    const auto *snp_101 = buffer.data(snp + 101);
    const auto *snp_103 = buffer.data(snp + 103);
    const auto *snp_104 = buffer.data(snp + 104);
    const auto *snp_105 = buffer.data(snp + 105);
    const auto *snp_106 = buffer.data(snp + 106);
    const auto *snp_107 = buffer.data(snp + 107);
    const auto *snp_108 = buffer.data(snp + 108);
    const auto *snp_109 = buffer.data(snp + 109);
    const auto *snp_110 = buffer.data(snp + 110);
    const auto *snp_112 = buffer.data(snp + 112);
    const auto *snp_113 = buffer.data(snp + 113);
    const auto *snp_114 = buffer.data(snp + 114);
    const auto *snp_115 = buffer.data(snp + 115);
    const auto *snp_116 = buffer.data(snp + 116);
    const auto *snp_117 = buffer.data(snp + 117);

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, smp_60, smp_61, \
                         smp_62, sns0_20, sns1_20, snp_60, snp_61, \
                         snp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_13 * smp_60[k]
                   + f_1 * sns0_20[k]
                   - f_2 * sns1_20[k]
                   + f_3 * pc_x[k] * snp_60[k];

        t_121[k] = f_13 * smp_61[k]
                   + f_3 * pc_x[k] * snp_61[k];

        t_122[k] = f_13 * smp_62[k]
                   + f_3 * pc_x[k] * snp_62[k];

        t_123[k] = f_1 * sns0_20[k]
                   - f_2 * sns1_20[k]
                   + f_3 * pc_y[k] * snp_61[k];

        t_124[k] = f_3 * pc_y[k] * snp_62[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, pc_x, pc_z, smp_44, smp_63, smp_64, sns0_20, \
                         sns0_21, sns1_20, sns1_21, snp_62, snp_63, \
                         snp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * smp_44[k]
                   + f_1 * sns0_20[k]
                   - f_2 * sns1_20[k]
                   + f_3 * pc_z[k] * snp_62[k];

        t_126[k] = f_12 * smp_63[k]
                   + f_1 * sns0_21[k]
                   - f_2 * sns1_21[k]
                   + f_3 * pc_x[k] * snp_63[k];

        t_127[k] = f_12 * smp_64[k]
                   + f_3 * pc_x[k] * snp_64[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, smp_46, smp_47, smp_65, \
                         sns0_21, sns1_21, snp_64, snp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_12 * smp_65[k]
                   + f_3 * pc_x[k] * snp_65[k];

        t_129[k] = f_11 * smp_46[k]
                   + f_1 * sns0_21[k]
                   - f_2 * sns1_21[k]
                   + f_3 * pc_y[k] * snp_64[k];

        t_130[k] = f_11 * smp_47[k]
                   + f_3 * pc_y[k] * snp_65[k];

        t_131[k] = f_1 * sns0_21[k]
                   - f_2 * sns1_21[k]
                   + f_3 * pc_z[k] * snp_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_z, pc_x, pc_z, smd0_90, smd0_93, \
                         smp_67, smp_68, smd1_90, smd1_93, snp_67, \
                         snp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * smd0_90[k]
                   - f_4 * pc_z[k] * smd1_90[k];

        t_133[k] = f_12 * smp_67[k]
                   + f_3 * pc_x[k] * snp_67[k];

        t_134[k] = f_12 * smp_68[k]
                   + f_3 * pc_x[k] * snp_68[k];

        t_135[k] = pb_z[k] * smd0_93[k]
                   - f_4 * pc_z[k] * smd1_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_x, pc_y, pc_z, smp_47, smp_50, smp_69, \
                         sns0_22, sns0_23, sns1_22, sns1_23, snp_68, \
                         snp_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_13 * smp_50[k]
                   + f_3 * pc_y[k] * snp_68[k];

        t_137[k] = f_6 * smp_47[k]
                   + f_1 * sns0_22[k]
                   - f_2 * sns1_22[k]
                   + f_3 * pc_z[k] * snp_68[k];

        t_138[k] = f_12 * smp_69[k]
                   + f_1 * sns0_23[k]
                   - f_2 * sns1_23[k]
                   + f_3 * pc_x[k] * snp_69[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, smp_52, smp_53, smp_70, \
                         smp_71, sns0_23, sns1_23, snp_70, snp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_12 * smp_70[k]
                   + f_3 * pc_x[k] * snp_70[k];

        t_140[k] = f_12 * smp_71[k]
                   + f_3 * pc_x[k] * snp_71[k];

        t_141[k] = f_12 * smp_52[k]
                   + f_1 * sns0_23[k]
                   - f_2 * sns1_23[k]
                   + f_3 * pc_y[k] * snp_70[k];

        t_142[k] = f_12 * smp_53[k]
                   + f_3 * pc_y[k] * snp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_z, smp_50, smp_72, smp_73, sns0_23, \
                         sns0_24, sns1_23, sns1_24, snp_71, snp_72, \
                         snp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * smp_50[k]
                   + f_1 * sns0_23[k]
                   - f_2 * sns1_23[k]
                   + f_3 * pc_z[k] * snp_71[k];

        t_144[k] = f_12 * smp_72[k]
                   + f_1 * sns0_24[k]
                   - f_2 * sns1_24[k]
                   + f_3 * pc_x[k] * snp_72[k];

        t_145[k] = f_12 * smp_73[k]
                   + f_3 * pc_x[k] * snp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, smp_53, smp_55, smp_56, \
                         smp_74, sns0_24, sns1_24, snp_73, snp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_12 * smp_74[k]
                   + f_3 * pc_x[k] * snp_74[k];

        t_147[k] = f_10 * smp_55[k]
                   + f_1 * sns0_24[k]
                   - f_2 * sns1_24[k]
                   + f_3 * pc_y[k] * snp_73[k];

        t_148[k] = f_10 * smp_56[k]
                   + f_3 * pc_y[k] * snp_74[k];

        t_149[k] = f_10 * smp_53[k]
                   + f_1 * sns0_24[k]
                   - f_2 * sns1_24[k]
                   + f_3 * pc_z[k] * snp_74[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, smp_58, smp_75, smp_76, \
                         smp_77, sns0_25, sns1_25, snp_75, snp_76, \
                         snp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_12 * smp_75[k]
                   + f_1 * sns0_25[k]
                   - f_2 * sns1_25[k]
                   + f_3 * pc_x[k] * snp_75[k];

        t_151[k] = f_12 * smp_76[k]
                   + f_3 * pc_x[k] * snp_76[k];

        t_152[k] = f_12 * smp_77[k]
                   + f_3 * pc_x[k] * snp_77[k];

        t_153[k] = f_8 * smp_58[k]
                   + f_1 * sns0_25[k]
                   - f_2 * sns1_25[k]
                   + f_3 * pc_y[k] * snp_76[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_y, pc_y, pc_z, smd0_120, smp_56, smp_59, \
                         smd1_120, sns0_25, sns1_25, snp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_8 * smp_59[k]
                   + f_3 * pc_y[k] * snp_77[k];

        t_155[k] = f_12 * smp_56[k]
                   + f_1 * sns0_25[k]
                   - f_2 * sns1_25[k]
                   + f_3 * pc_z[k] * snp_77[k];

        t_156[k] = pb_y[k] * smd0_120[k]
                   - f_4 * pc_y[k] * smd1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, smp_61, smp_62, smp_79, \
                         smp_80, sns0_26, sns1_26, snp_79, snp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_12 * smp_79[k]
                   + f_3 * pc_x[k] * snp_79[k];

        t_158[k] = f_12 * smp_80[k]
                   + f_3 * pc_x[k] * snp_80[k];

        t_159[k] = f_6 * smp_61[k]
                   + f_1 * sns0_26[k]
                   - f_2 * sns1_26[k]
                   + f_3 * pc_y[k] * snp_79[k];

        t_160[k] = f_6 * smp_62[k]
                   + f_3 * pc_y[k] * snp_80[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_y, pc_x, pc_y, smd0_125, smp_81, smp_82, \
                         smd1_125, sns0_27, sns1_27, snp_81, snp_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_y[k] * smd0_125[k]
                   - f_4 * pc_y[k] * smd1_125[k];

        t_162[k] = f_12 * smp_81[k]
                   + f_1 * sns0_27[k]
                   - f_2 * sns1_27[k]
                   + f_3 * pc_x[k] * snp_81[k];

        t_163[k] = f_12 * smp_82[k]
                   + f_3 * pc_x[k] * snp_82[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, smp_62, smp_83, \
                         sns0_27, sns1_27, snp_82, snp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_12 * smp_83[k]
                   + f_3 * pc_x[k] * snp_83[k];

        t_165[k] = f_1 * sns0_27[k]
                   - f_2 * sns1_27[k]
                   + f_3 * pc_y[k] * snp_82[k];

        t_166[k] = f_3 * pc_y[k] * snp_83[k];

        t_167[k] = f_11 * smp_62[k]
                   + f_1 * sns0_27[k]
                   - f_2 * sns1_27[k]
                   + f_3 * pc_z[k] * snp_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, smp_64, smp_84, smp_85, \
                         smp_86, sns0_28, sns1_28, snp_84, snp_85, \
                         snp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_10 * smp_84[k]
                   + f_1 * sns0_28[k]
                   - f_2 * sns1_28[k]
                   + f_3 * pc_x[k] * snp_84[k];

        t_169[k] = f_10 * smp_85[k]
                   + f_3 * pc_x[k] * snp_85[k];

        t_170[k] = f_10 * smp_86[k]
                   + f_3 * pc_x[k] * snp_86[k];

        t_171[k] = f_9 * smp_64[k]
                   + f_1 * sns0_28[k]
                   - f_2 * sns1_28[k]
                   + f_3 * pc_y[k] * snp_85[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pb_z, pc_x, pc_y, pc_z, smd0_126, smp_65, \
                         smp_88, smd1_126, sns0_28, sns1_28, snp_86, \
                         snp_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_9 * smp_65[k]
                   + f_3 * pc_y[k] * snp_86[k];

        t_173[k] = f_1 * sns0_28[k]
                   - f_2 * sns1_28[k]
                   + f_3 * pc_z[k] * snp_86[k];

        t_174[k] = pb_z[k] * smd0_126[k]
                   - f_4 * pc_z[k] * smd1_126[k];

        t_175[k] = f_10 * smp_88[k]
                   + f_3 * pc_x[k] * snp_88[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pb_z, pc_x, pc_y, pc_z, smd0_129, smp_65, \
                         smp_68, smp_89, smd1_129, sns0_29, sns1_29, \
                         snp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_10 * smp_89[k]
                   + f_3 * pc_x[k] * snp_89[k];

        t_177[k] = pb_z[k] * smd0_129[k]
                   - f_4 * pc_z[k] * smd1_129[k];

        t_178[k] = f_11 * smp_68[k]
                   + f_3 * pc_y[k] * snp_89[k];

        t_179[k] = f_6 * smp_65[k]
                   + f_1 * sns0_29[k]
                   - f_2 * sns1_29[k]
                   + f_3 * pc_z[k] * snp_89[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, smp_70, smp_90, smp_91, \
                         smp_92, sns0_30, sns1_30, snp_90, snp_91, \
                         snp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_10 * smp_90[k]
                   + f_1 * sns0_30[k]
                   - f_2 * sns1_30[k]
                   + f_3 * pc_x[k] * snp_90[k];

        t_181[k] = f_10 * smp_91[k]
                   + f_3 * pc_x[k] * snp_91[k];

        t_182[k] = f_10 * smp_92[k]
                   + f_3 * pc_x[k] * snp_92[k];

        t_183[k] = f_13 * smp_70[k]
                   + f_1 * sns0_30[k]
                   - f_2 * sns1_30[k]
                   + f_3 * pc_y[k] * snp_91[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pc_x, pc_y, pc_z, smp_68, smp_71, smp_93, \
                         sns0_30, sns0_31, sns1_30, sns1_31, snp_92, \
                         snp_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_13 * smp_71[k]
                   + f_3 * pc_y[k] * snp_92[k];

        t_185[k] = f_8 * smp_68[k]
                   + f_1 * sns0_30[k]
                   - f_2 * sns1_30[k]
                   + f_3 * pc_z[k] * snp_92[k];

        t_186[k] = f_10 * smp_93[k]
                   + f_1 * sns0_31[k]
                   - f_2 * sns1_31[k]
                   + f_3 * pc_x[k] * snp_93[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pc_x, pc_y, smp_73, smp_74, smp_94, \
                         smp_95, sns0_31, sns1_31, snp_94, snp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_10 * smp_94[k]
                   + f_3 * pc_x[k] * snp_94[k];

        t_188[k] = f_10 * smp_95[k]
                   + f_3 * pc_x[k] * snp_95[k];

        t_189[k] = f_12 * smp_73[k]
                   + f_1 * sns0_31[k]
                   - f_2 * sns1_31[k]
                   + f_3 * pc_y[k] * snp_94[k];

        t_190[k] = f_12 * smp_74[k]
                   + f_3 * pc_y[k] * snp_95[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pc_x, pc_z, smp_71, smp_96, smp_97, sns0_31, \
                         sns0_32, sns1_31, sns1_32, snp_95, snp_96, \
                         snp_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_10 * smp_71[k]
                   + f_1 * sns0_31[k]
                   - f_2 * sns1_31[k]
                   + f_3 * pc_z[k] * snp_95[k];

        t_192[k] = f_10 * smp_96[k]
                   + f_1 * sns0_32[k]
                   - f_2 * sns1_32[k]
                   + f_3 * pc_x[k] * snp_96[k];

        t_193[k] = f_10 * smp_97[k]
                   + f_3 * pc_x[k] * snp_97[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, pc_z, smp_74, smp_76, smp_77, \
                         smp_98, sns0_32, sns1_32, snp_97, snp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_10 * smp_98[k]
                   + f_3 * pc_x[k] * snp_98[k];

        t_195[k] = f_10 * smp_76[k]
                   + f_1 * sns0_32[k]
                   - f_2 * sns1_32[k]
                   + f_3 * pc_y[k] * snp_97[k];

        t_196[k] = f_10 * smp_77[k]
                   + f_3 * pc_y[k] * snp_98[k];

        t_197[k] = f_12 * smp_74[k]
                   + f_1 * sns0_32[k]
                   - f_2 * sns1_32[k]
                   + f_3 * pc_z[k] * snp_98[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, pc_y, smp_79, smp_99, smp_100, \
                         smp_101, sns0_33, sns1_33, snp_99, snp_100, \
                         snp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_10 * smp_99[k]
                   + f_1 * sns0_33[k]
                   - f_2 * sns1_33[k]
                   + f_3 * pc_x[k] * snp_99[k];

        t_199[k] = f_10 * smp_100[k]
                   + f_3 * pc_x[k] * snp_100[k];

        t_200[k] = f_10 * smp_101[k]
                   + f_3 * pc_x[k] * snp_101[k];

        t_201[k] = f_8 * smp_79[k]
                   + f_1 * sns0_33[k]
                   - f_2 * sns1_33[k]
                   + f_3 * pc_y[k] * snp_100[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_y, pc_y, pc_z, smd0_162, smp_77, smp_80, \
                         smd1_162, sns0_33, sns1_33, snp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_8 * smp_80[k]
                   + f_3 * pc_y[k] * snp_101[k];

        t_203[k] = f_13 * smp_77[k]
                   + f_1 * sns0_33[k]
                   - f_2 * sns1_33[k]
                   + f_3 * pc_z[k] * snp_101[k];

        t_204[k] = pb_y[k] * smd0_162[k]
                   - f_4 * pc_y[k] * smd1_162[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, smp_82, smp_83, smp_103, \
                         smp_104, sns0_34, sns1_34, snp_103, snp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_10 * smp_103[k]
                   + f_3 * pc_x[k] * snp_103[k];

        t_206[k] = f_10 * smp_104[k]
                   + f_3 * pc_x[k] * snp_104[k];

        t_207[k] = f_6 * smp_82[k]
                   + f_1 * sns0_34[k]
                   - f_2 * sns1_34[k]
                   + f_3 * pc_y[k] * snp_103[k];

        t_208[k] = f_6 * smp_83[k]
                   + f_3 * pc_y[k] * snp_104[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pb_y, pc_x, pc_y, smd0_167, smp_105, smp_106, \
                         smd1_167, sns0_35, sns1_35, snp_105, snp_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pb_y[k] * smd0_167[k]
                   - f_4 * pc_y[k] * smd1_167[k];

        t_210[k] = f_10 * smp_105[k]
                   + f_1 * sns0_35[k]
                   - f_2 * sns1_35[k]
                   + f_3 * pc_x[k] * snp_105[k];

        t_211[k] = f_10 * smp_106[k]
                   + f_3 * pc_x[k] * snp_106[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, smp_83, smp_107, \
                         sns0_35, sns1_35, snp_106, snp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_10 * smp_107[k]
                   + f_3 * pc_x[k] * snp_107[k];

        t_213[k] = f_1 * sns0_35[k]
                   - f_2 * sns1_35[k]
                   + f_3 * pc_y[k] * snp_106[k];

        t_214[k] = f_3 * pc_y[k] * snp_107[k];

        t_215[k] = f_9 * smp_83[k]
                   + f_1 * sns0_35[k]
                   - f_2 * sns1_35[k]
                   + f_3 * pc_z[k] * snp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, smp_85, smp_108, smp_109, \
                         smp_110, sns0_36, sns1_36, snp_108, snp_109, \
                         snp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_8 * smp_108[k]
                   + f_1 * sns0_36[k]
                   - f_2 * sns1_36[k]
                   + f_3 * pc_x[k] * snp_108[k];

        t_217[k] = f_8 * smp_109[k]
                   + f_3 * pc_x[k] * snp_109[k];

        t_218[k] = f_8 * smp_110[k]
                   + f_3 * pc_x[k] * snp_110[k];

        t_219[k] = f_7 * smp_85[k]
                   + f_1 * sns0_36[k]
                   - f_2 * sns1_36[k]
                   + f_3 * pc_y[k] * snp_109[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_z, pc_x, pc_y, pc_z, smd0_168, smp_86, \
                         smp_112, smd1_168, sns0_36, sns1_36, snp_110, \
                         snp_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_7 * smp_86[k]
                   + f_3 * pc_y[k] * snp_110[k];

        t_221[k] = f_1 * sns0_36[k]
                   - f_2 * sns1_36[k]
                   + f_3 * pc_z[k] * snp_110[k];

        t_222[k] = pb_z[k] * smd0_168[k]
                   - f_4 * pc_z[k] * smd1_168[k];

        t_223[k] = f_8 * smp_112[k]
                   + f_3 * pc_x[k] * snp_112[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pb_z, pc_x, pc_y, pc_z, smd0_171, smp_86, \
                         smp_89, smp_113, smd1_171, sns0_37, sns1_37, \
                         snp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_8 * smp_113[k]
                   + f_3 * pc_x[k] * snp_113[k];

        t_225[k] = pb_z[k] * smd0_171[k]
                   - f_4 * pc_z[k] * smd1_171[k];

        t_226[k] = f_9 * smp_89[k]
                   + f_3 * pc_y[k] * snp_113[k];

        t_227[k] = f_6 * smp_86[k]
                   + f_1 * sns0_37[k]
                   - f_2 * sns1_37[k]
                   + f_3 * pc_z[k] * snp_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, smp_91, smp_114, smp_115, \
                         smp_116, sns0_38, sns1_38, snp_114, snp_115, \
                         snp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_8 * smp_114[k]
                   + f_1 * sns0_38[k]
                   - f_2 * sns1_38[k]
                   + f_3 * pc_x[k] * snp_114[k];

        t_229[k] = f_8 * smp_115[k]
                   + f_3 * pc_x[k] * snp_115[k];

        t_230[k] = f_8 * smp_116[k]
                   + f_3 * pc_x[k] * snp_116[k];

        t_231[k] = f_11 * smp_91[k]
                   + f_1 * sns0_38[k]
                   - f_2 * sns1_38[k]
                   + f_3 * pc_y[k] * snp_115[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pc_x, pc_y, pc_z, smp_89, smp_92, smp_117, \
                         sns0_38, sns0_39, sns1_38, sns1_39, snp_116, \
                         snp_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_11 * smp_92[k]
                   + f_3 * pc_y[k] * snp_116[k];

        t_233[k] = f_8 * smp_89[k]
                   + f_1 * sns0_38[k]
                   - f_2 * sns1_38[k]
                   + f_3 * pc_z[k] * snp_116[k];

        t_234[k] = f_8 * smp_117[k]
                   + f_1 * sns0_39[k]
                   - f_2 * sns1_39[k]
                   + f_3 * pc_x[k] * snp_117[k];
    }
}

static auto
compute_prim_snd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smd0,
                                                          const size_t smp, const size_t smd1,
                                                          const size_t sns0, const size_t sns1,
                                                          const size_t snp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

    auto *t_235 = buffer.data(target + 235);
    auto *t_236 = buffer.data(target + 236);
    auto *t_237 = buffer.data(target + 237);
    auto *t_238 = buffer.data(target + 238);
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
    auto *t_330 = buffer.data(target + 330);
    auto *t_331 = buffer.data(target + 331);
    auto *t_332 = buffer.data(target + 332);
    auto *t_333 = buffer.data(target + 333);
    auto *t_334 = buffer.data(target + 334);
    auto *t_335 = buffer.data(target + 335);
    auto *t_336 = buffer.data(target + 336);
    auto *t_337 = buffer.data(target + 337);
    auto *t_338 = buffer.data(target + 338);
    auto *t_339 = buffer.data(target + 339);
    auto *t_340 = buffer.data(target + 340);
    auto *t_341 = buffer.data(target + 341);
    auto *t_342 = buffer.data(target + 342);
    auto *t_343 = buffer.data(target + 343);
    auto *t_344 = buffer.data(target + 344);
    auto *t_345 = buffer.data(target + 345);
    auto *t_346 = buffer.data(target + 346);
    auto *t_347 = buffer.data(target + 347);
    auto *t_348 = buffer.data(target + 348);
    auto *t_349 = buffer.data(target + 349);
    auto *t_350 = buffer.data(target + 350);
    auto *t_351 = buffer.data(target + 351);
    auto *t_352 = buffer.data(target + 352);
    auto *t_353 = buffer.data(target + 353);
    auto *t_354 = buffer.data(target + 354);
    auto *t_355 = buffer.data(target + 355);
    auto *t_356 = buffer.data(target + 356);
    auto *t_357 = buffer.data(target + 357);
    auto *t_358 = buffer.data(target + 358);
    auto *t_359 = buffer.data(target + 359);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smd0_210 = buffer.data(smd0 + 210);
    const auto *smd0_215 = buffer.data(smd0 + 215);
    const auto *smd0_216 = buffer.data(smd0 + 216);
    const auto *smd0_264 = buffer.data(smd0 + 264);
    const auto *smd0_270 = buffer.data(smd0 + 270);
    const auto *smd0_273 = buffer.data(smd0 + 273);
    const auto *smd0_275 = buffer.data(smd0 + 275);
    const auto *smd0_279 = buffer.data(smd0 + 279);
    const auto *smd0_281 = buffer.data(smd0 + 281);
    const auto *smd0_282 = buffer.data(smd0 + 282);
    const auto *smd0_285 = buffer.data(smd0 + 285);
    const auto *smd0_287 = buffer.data(smd0 + 287);
    const auto *smd0_288 = buffer.data(smd0 + 288);
    const auto *smd0_291 = buffer.data(smd0 + 291);
    const auto *smd0_293 = buffer.data(smd0 + 293);
    const auto *smd0_294 = buffer.data(smd0 + 294);
    const auto *smd0_297 = buffer.data(smd0 + 297);
    const auto *smd0_299 = buffer.data(smd0 + 299);
    const auto *smd0_300 = buffer.data(smd0 + 300);
    const auto *smd0_303 = buffer.data(smd0 + 303);
    const auto *smd0_305 = buffer.data(smd0 + 305);
    const auto *smd0_306 = buffer.data(smd0 + 306);
    const auto *smd0_309 = buffer.data(smd0 + 309);
    const auto *smd0_311 = buffer.data(smd0 + 311);
    const auto *smd0_312 = buffer.data(smd0 + 312);
    const auto *smd0_315 = buffer.data(smd0 + 315);
    const auto *smd0_317 = buffer.data(smd0 + 317);
    const auto *smd0_321 = buffer.data(smd0 + 321);
    const auto *smd0_323 = buffer.data(smd0 + 323);
    const auto *smd0_324 = buffer.data(smd0 + 324);
    const auto *smd0_327 = buffer.data(smd0 + 327);
    const auto *smd0_329 = buffer.data(smd0 + 329);

    const auto *smp_92 = buffer.data(smp + 92);
    const auto *smp_94 = buffer.data(smp + 94);
    const auto *smp_95 = buffer.data(smp + 95);
    const auto *smp_97 = buffer.data(smp + 97);
    const auto *smp_98 = buffer.data(smp + 98);
    const auto *smp_100 = buffer.data(smp + 100);
    const auto *smp_101 = buffer.data(smp + 101);
    const auto *smp_103 = buffer.data(smp + 103);
    const auto *smp_104 = buffer.data(smp + 104);
    const auto *smp_106 = buffer.data(smp + 106);
    const auto *smp_107 = buffer.data(smp + 107);
    const auto *smp_110 = buffer.data(smp + 110);
    const auto *smp_113 = buffer.data(smp + 113);
    const auto *smp_116 = buffer.data(smp + 116);
    const auto *smp_118 = buffer.data(smp + 118);
    const auto *smp_119 = buffer.data(smp + 119);
    const auto *smp_120 = buffer.data(smp + 120);
    const auto *smp_121 = buffer.data(smp + 121);
    const auto *smp_122 = buffer.data(smp + 122);
    const auto *smp_123 = buffer.data(smp + 123);
    const auto *smp_124 = buffer.data(smp + 124);
    const auto *smp_125 = buffer.data(smp + 125);
    const auto *smp_126 = buffer.data(smp + 126);
    const auto *smp_127 = buffer.data(smp + 127);
    const auto *smp_128 = buffer.data(smp + 128);
    const auto *smp_130 = buffer.data(smp + 130);
    const auto *smp_131 = buffer.data(smp + 131);
    const auto *smp_132 = buffer.data(smp + 132);
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

    const auto *smd1_210 = buffer.data(smd1 + 210);
    const auto *smd1_215 = buffer.data(smd1 + 215);
    const auto *smd1_216 = buffer.data(smd1 + 216);
    const auto *smd1_264 = buffer.data(smd1 + 264);
    const auto *smd1_270 = buffer.data(smd1 + 270);
    const auto *smd1_273 = buffer.data(smd1 + 273);
    const auto *smd1_275 = buffer.data(smd1 + 275);
    const auto *smd1_279 = buffer.data(smd1 + 279);
    const auto *smd1_281 = buffer.data(smd1 + 281);
    const auto *smd1_282 = buffer.data(smd1 + 282);
    const auto *smd1_285 = buffer.data(smd1 + 285);
    const auto *smd1_287 = buffer.data(smd1 + 287);
    const auto *smd1_288 = buffer.data(smd1 + 288);
    const auto *smd1_291 = buffer.data(smd1 + 291);
    const auto *smd1_293 = buffer.data(smd1 + 293);
    const auto *smd1_294 = buffer.data(smd1 + 294);
    const auto *smd1_297 = buffer.data(smd1 + 297);
    const auto *smd1_299 = buffer.data(smd1 + 299);
    const auto *smd1_300 = buffer.data(smd1 + 300);
    const auto *smd1_303 = buffer.data(smd1 + 303);
    const auto *smd1_305 = buffer.data(smd1 + 305);
    const auto *smd1_306 = buffer.data(smd1 + 306);
    const auto *smd1_309 = buffer.data(smd1 + 309);
    const auto *smd1_311 = buffer.data(smd1 + 311);
    const auto *smd1_312 = buffer.data(smd1 + 312);
    const auto *smd1_315 = buffer.data(smd1 + 315);
    const auto *smd1_317 = buffer.data(smd1 + 317);
    const auto *smd1_321 = buffer.data(smd1 + 321);
    const auto *smd1_323 = buffer.data(smd1 + 323);
    const auto *smd1_324 = buffer.data(smd1 + 324);
    const auto *smd1_327 = buffer.data(smd1 + 327);
    const auto *smd1_329 = buffer.data(smd1 + 329);

    const auto *sns0_39 = buffer.data(sns0 + 39);
    const auto *sns0_40 = buffer.data(sns0 + 40);
    const auto *sns0_41 = buffer.data(sns0 + 41);
    const auto *sns0_42 = buffer.data(sns0 + 42);
    const auto *sns0_43 = buffer.data(sns0 + 43);
    const auto *sns0_44 = buffer.data(sns0 + 44);
    const auto *sns0_55 = buffer.data(sns0 + 55);
    const auto *sns0_56 = buffer.data(sns0 + 56);
    const auto *sns0_57 = buffer.data(sns0 + 57);
    const auto *sns0_58 = buffer.data(sns0 + 58);
    const auto *sns0_59 = buffer.data(sns0 + 59);

    const auto *sns1_39 = buffer.data(sns1 + 39);
    const auto *sns1_40 = buffer.data(sns1 + 40);
    const auto *sns1_41 = buffer.data(sns1 + 41);
    const auto *sns1_42 = buffer.data(sns1 + 42);
    const auto *sns1_43 = buffer.data(sns1 + 43);
    const auto *sns1_44 = buffer.data(sns1 + 44);
    const auto *sns1_55 = buffer.data(sns1 + 55);
    const auto *sns1_56 = buffer.data(sns1 + 56);
    const auto *sns1_57 = buffer.data(sns1 + 57);
    const auto *sns1_58 = buffer.data(sns1 + 58);
    const auto *sns1_59 = buffer.data(sns1 + 59);

    const auto *snp_118 = buffer.data(snp + 118);
    const auto *snp_119 = buffer.data(snp + 119);
    const auto *snp_120 = buffer.data(snp + 120);
    const auto *snp_121 = buffer.data(snp + 121);
    const auto *snp_122 = buffer.data(snp + 122);
    const auto *snp_123 = buffer.data(snp + 123);
    const auto *snp_124 = buffer.data(snp + 124);
    const auto *snp_125 = buffer.data(snp + 125);
    const auto *snp_126 = buffer.data(snp + 126);
    const auto *snp_127 = buffer.data(snp + 127);
    const auto *snp_128 = buffer.data(snp + 128);
    const auto *snp_130 = buffer.data(snp + 130);
    const auto *snp_131 = buffer.data(snp + 131);
    const auto *snp_132 = buffer.data(snp + 132);
    const auto *snp_133 = buffer.data(snp + 133);
    const auto *snp_134 = buffer.data(snp + 134);
    const auto *snp_136 = buffer.data(snp + 136);
    const auto *snp_137 = buffer.data(snp + 137);
    const auto *snp_139 = buffer.data(snp + 139);
    const auto *snp_140 = buffer.data(snp + 140);
    const auto *snp_142 = buffer.data(snp + 142);
    const auto *snp_143 = buffer.data(snp + 143);
    const auto *snp_145 = buffer.data(snp + 145);
    const auto *snp_146 = buffer.data(snp + 146);
    const auto *snp_148 = buffer.data(snp + 148);
    const auto *snp_149 = buffer.data(snp + 149);
    const auto *snp_151 = buffer.data(snp + 151);
    const auto *snp_152 = buffer.data(snp + 152);
    const auto *snp_154 = buffer.data(snp + 154);
    const auto *snp_155 = buffer.data(snp + 155);
    const auto *snp_157 = buffer.data(snp + 157);
    const auto *snp_158 = buffer.data(snp + 158);
    const auto *snp_160 = buffer.data(snp + 160);
    const auto *snp_161 = buffer.data(snp + 161);
    const auto *snp_163 = buffer.data(snp + 163);
    const auto *snp_164 = buffer.data(snp + 164);
    const auto *snp_165 = buffer.data(snp + 165);
    const auto *snp_166 = buffer.data(snp + 166);
    const auto *snp_167 = buffer.data(snp + 167);
    const auto *snp_169 = buffer.data(snp + 169);
    const auto *snp_170 = buffer.data(snp + 170);
    const auto *snp_171 = buffer.data(snp + 171);
    const auto *snp_172 = buffer.data(snp + 172);
    const auto *snp_173 = buffer.data(snp + 173);
    const auto *snp_174 = buffer.data(snp + 174);
    const auto *snp_175 = buffer.data(snp + 175);
    const auto *snp_176 = buffer.data(snp + 176);
    const auto *snp_177 = buffer.data(snp + 177);
    const auto *snp_178 = buffer.data(snp + 178);
    const auto *snp_179 = buffer.data(snp + 179);

#pragma omp simd aligned(t_235, t_236, t_237, t_238, pc_x, pc_y, smp_94, smp_95, smp_118, \
                         smp_119, sns0_39, sns1_39, snp_118, snp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_8 * smp_118[k]
                   + f_3 * pc_x[k] * snp_118[k];

        t_236[k] = f_8 * smp_119[k]
                   + f_3 * pc_x[k] * snp_119[k];

        t_237[k] = f_13 * smp_94[k]
                   + f_1 * sns0_39[k]
                   - f_2 * sns1_39[k]
                   + f_3 * pc_y[k] * snp_118[k];

        t_238[k] = f_13 * smp_95[k]
                   + f_3 * pc_y[k] * snp_119[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, smp_92, smp_120, smp_121, sns0_39, \
                         sns0_40, sns1_39, sns1_40, snp_119, snp_120, \
                         snp_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_10 * smp_92[k]
                   + f_1 * sns0_39[k]
                   - f_2 * sns1_39[k]
                   + f_3 * pc_z[k] * snp_119[k];

        t_240[k] = f_8 * smp_120[k]
                   + f_1 * sns0_40[k]
                   - f_2 * sns1_40[k]
                   + f_3 * pc_x[k] * snp_120[k];

        t_241[k] = f_8 * smp_121[k]
                   + f_3 * pc_x[k] * snp_121[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pc_x, pc_y, pc_z, smp_95, smp_97, smp_98, \
                         smp_122, sns0_40, sns1_40, snp_121, snp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_8 * smp_122[k]
                   + f_3 * pc_x[k] * snp_122[k];

        t_243[k] = f_12 * smp_97[k]
                   + f_1 * sns0_40[k]
                   - f_2 * sns1_40[k]
                   + f_3 * pc_y[k] * snp_121[k];

        t_244[k] = f_12 * smp_98[k]
                   + f_3 * pc_y[k] * snp_122[k];

        t_245[k] = f_12 * smp_95[k]
                   + f_1 * sns0_40[k]
                   - f_2 * sns1_40[k]
                   + f_3 * pc_z[k] * snp_122[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, pc_y, smp_100, smp_123, smp_124, \
                         smp_125, sns0_41, sns1_41, snp_123, snp_124, \
                         snp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_8 * smp_123[k]
                   + f_1 * sns0_41[k]
                   - f_2 * sns1_41[k]
                   + f_3 * pc_x[k] * snp_123[k];

        t_247[k] = f_8 * smp_124[k]
                   + f_3 * pc_x[k] * snp_124[k];

        t_248[k] = f_8 * smp_125[k]
                   + f_3 * pc_x[k] * snp_125[k];

        t_249[k] = f_10 * smp_100[k]
                   + f_1 * sns0_41[k]
                   - f_2 * sns1_41[k]
                   + f_3 * pc_y[k] * snp_124[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pc_x, pc_y, pc_z, smp_98, smp_101, smp_126, \
                         sns0_41, sns0_42, sns1_41, sns1_42, snp_125, \
                         snp_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_10 * smp_101[k]
                   + f_3 * pc_y[k] * snp_125[k];

        t_251[k] = f_13 * smp_98[k]
                   + f_1 * sns0_41[k]
                   - f_2 * sns1_41[k]
                   + f_3 * pc_z[k] * snp_125[k];

        t_252[k] = f_8 * smp_126[k]
                   + f_1 * sns0_42[k]
                   - f_2 * sns1_42[k]
                   + f_3 * pc_x[k] * snp_126[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pc_x, pc_y, smp_103, smp_104, smp_127, \
                         smp_128, sns0_42, sns1_42, snp_127, snp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_8 * smp_127[k]
                   + f_3 * pc_x[k] * snp_127[k];

        t_254[k] = f_8 * smp_128[k]
                   + f_3 * pc_x[k] * snp_128[k];

        t_255[k] = f_8 * smp_103[k]
                   + f_1 * sns0_42[k]
                   - f_2 * sns1_42[k]
                   + f_3 * pc_y[k] * snp_127[k];

        t_256[k] = f_8 * smp_104[k]
                   + f_3 * pc_y[k] * snp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_y, pc_x, pc_y, pc_z, smd0_210, smp_101, \
                         smp_130, smd1_210, sns0_42, sns1_42, snp_128, \
                         snp_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_11 * smp_101[k]
                   + f_1 * sns0_42[k]
                   - f_2 * sns1_42[k]
                   + f_3 * pc_z[k] * snp_128[k];

        t_258[k] = pb_y[k] * smd0_210[k]
                   - f_4 * pc_y[k] * smd1_210[k];

        t_259[k] = f_8 * smp_130[k]
                   + f_3 * pc_x[k] * snp_130[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pb_y, pc_x, pc_y, smd0_215, smp_106, \
                         smp_107, smp_131, smd1_215, sns0_43, sns1_43, snp_130, \
                         snp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_8 * smp_131[k]
                   + f_3 * pc_x[k] * snp_131[k];

        t_261[k] = f_6 * smp_106[k]
                   + f_1 * sns0_43[k]
                   - f_2 * sns1_43[k]
                   + f_3 * pc_y[k] * snp_130[k];

        t_262[k] = f_6 * smp_107[k]
                   + f_3 * pc_y[k] * snp_131[k];

        t_263[k] = pb_y[k] * smd0_215[k]
                   - f_4 * pc_y[k] * smd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, pc_x, pc_y, smp_132, smp_133, \
                         smp_134, sns0_44, sns1_44, snp_132, snp_133, \
                         snp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * smp_132[k]
                   + f_1 * sns0_44[k]
                   - f_2 * sns1_44[k]
                   + f_3 * pc_x[k] * snp_132[k];

        t_265[k] = f_8 * smp_133[k]
                   + f_3 * pc_x[k] * snp_133[k];

        t_266[k] = f_8 * smp_134[k]
                   + f_3 * pc_x[k] * snp_134[k];

        t_267[k] = f_1 * sns0_44[k]
                   - f_2 * sns1_44[k]
                   + f_3 * pc_y[k] * snp_133[k];

        t_268[k] = f_3 * pc_y[k] * snp_134[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pb_x, pc_x, pc_z, smd0_270, smp_107, smp_135, \
                         smp_136, smd1_270, sns0_44, sns1_44, snp_134, \
                         snp_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_7 * smp_107[k]
                   + f_1 * sns0_44[k]
                   - f_2 * sns1_44[k]
                   + f_3 * pc_z[k] * snp_134[k];

        t_270[k] = pb_x[k] * smd0_270[k]
                   + f_8 * smp_135[k]
                   - f_4 * pc_x[k] * smd1_270[k];

        t_271[k] = f_6 * smp_136[k]
                   + f_3 * pc_x[k] * snp_136[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pb_x, pc_x, pc_y, smd0_273, smd0_275, \
                         smp_110, smp_137, smd1_273, smd1_275, \
                         snp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_6 * smp_137[k]
                   + f_3 * pc_x[k] * snp_137[k];

        t_273[k] = pb_x[k] * smd0_273[k]
                   - f_4 * pc_x[k] * smd1_273[k];

        t_274[k] = f_5 * smp_110[k]
                   + f_3 * pc_y[k] * snp_137[k];

        t_275[k] = pb_x[k] * smd0_275[k]
                   - f_4 * pc_x[k] * smd1_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pb_x, pb_z, pc_x, pc_z, smd0_216, \
                         smd0_279, smp_139, smp_140, smd1_216, smd1_279, snp_139, \
                         snp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_z[k] * smd0_216[k]
                   - f_4 * pc_z[k] * smd1_216[k];

        t_277[k] = f_6 * smp_139[k]
                   + f_3 * pc_x[k] * snp_139[k];

        t_278[k] = f_6 * smp_140[k]
                   + f_3 * pc_x[k] * snp_140[k];

        t_279[k] = pb_x[k] * smd0_279[k]
                   - f_4 * pc_x[k] * smd1_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pb_x, pc_x, pc_y, smd0_281, smd0_282, \
                         smp_113, smp_141, smp_142, smd1_281, smd1_282, snp_140, \
                         snp_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_7 * smp_113[k]
                   + f_3 * pc_y[k] * snp_140[k];

        t_281[k] = pb_x[k] * smd0_281[k]
                   - f_4 * pc_x[k] * smd1_281[k];

        t_282[k] = pb_x[k] * smd0_282[k]
                   + f_8 * smp_141[k]
                   - f_4 * pc_x[k] * smd1_282[k];

        t_283[k] = f_6 * smp_142[k]
                   + f_3 * pc_x[k] * snp_142[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pb_x, pc_x, pc_y, smd0_285, smd0_287, \
                         smp_116, smp_143, smd1_285, smd1_287, \
                         snp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_6 * smp_143[k]
                   + f_3 * pc_x[k] * snp_143[k];

        t_285[k] = pb_x[k] * smd0_285[k]
                   - f_4 * pc_x[k] * smd1_285[k];

        t_286[k] = f_9 * smp_116[k]
                   + f_3 * pc_y[k] * snp_143[k];

        t_287[k] = pb_x[k] * smd0_287[k]
                   - f_4 * pc_x[k] * smd1_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pb_x, pc_x, smd0_288, smd0_291, smp_144, \
                         smp_145, smp_146, smd1_288, smd1_291, snp_145, \
                         snp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pb_x[k] * smd0_288[k]
                   + f_8 * smp_144[k]
                   - f_4 * pc_x[k] * smd1_288[k];

        t_289[k] = f_6 * smp_145[k]
                   + f_3 * pc_x[k] * snp_145[k];

        t_290[k] = f_6 * smp_146[k]
                   + f_3 * pc_x[k] * snp_146[k];

        t_291[k] = pb_x[k] * smd0_291[k]
                   - f_4 * pc_x[k] * smd1_291[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pb_x, pc_x, pc_y, smd0_293, smd0_294, \
                         smp_119, smp_147, smp_148, smd1_293, smd1_294, snp_146, \
                         snp_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_11 * smp_119[k]
                   + f_3 * pc_y[k] * snp_146[k];

        t_293[k] = pb_x[k] * smd0_293[k]
                   - f_4 * pc_x[k] * smd1_293[k];

        t_294[k] = pb_x[k] * smd0_294[k]
                   + f_8 * smp_147[k]
                   - f_4 * pc_x[k] * smd1_294[k];

        t_295[k] = f_6 * smp_148[k]
                   + f_3 * pc_x[k] * snp_148[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pb_x, pc_x, pc_y, smd0_297, smd0_299, \
                         smp_122, smp_149, smd1_297, smd1_299, \
                         snp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_6 * smp_149[k]
                   + f_3 * pc_x[k] * snp_149[k];

        t_297[k] = pb_x[k] * smd0_297[k]
                   - f_4 * pc_x[k] * smd1_297[k];

        t_298[k] = f_13 * smp_122[k]
                   + f_3 * pc_y[k] * snp_149[k];

        t_299[k] = pb_x[k] * smd0_299[k]
                   - f_4 * pc_x[k] * smd1_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pb_x, pc_x, smd0_300, smd0_303, smp_150, \
                         smp_151, smp_152, smd1_300, smd1_303, snp_151, \
                         snp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pb_x[k] * smd0_300[k]
                   + f_8 * smp_150[k]
                   - f_4 * pc_x[k] * smd1_300[k];

        t_301[k] = f_6 * smp_151[k]
                   + f_3 * pc_x[k] * snp_151[k];

        t_302[k] = f_6 * smp_152[k]
                   + f_3 * pc_x[k] * snp_152[k];

        t_303[k] = pb_x[k] * smd0_303[k]
                   - f_4 * pc_x[k] * smd1_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pb_x, pc_x, pc_y, smd0_305, smd0_306, \
                         smp_125, smp_153, smp_154, smd1_305, smd1_306, snp_152, \
                         snp_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_12 * smp_125[k]
                   + f_3 * pc_y[k] * snp_152[k];

        t_305[k] = pb_x[k] * smd0_305[k]
                   - f_4 * pc_x[k] * smd1_305[k];

        t_306[k] = pb_x[k] * smd0_306[k]
                   + f_8 * smp_153[k]
                   - f_4 * pc_x[k] * smd1_306[k];

        t_307[k] = f_6 * smp_154[k]
                   + f_3 * pc_x[k] * snp_154[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_x, pc_x, pc_y, smd0_309, smd0_311, \
                         smp_128, smp_155, smd1_309, smd1_311, \
                         snp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_6 * smp_155[k]
                   + f_3 * pc_x[k] * snp_155[k];

        t_309[k] = pb_x[k] * smd0_309[k]
                   - f_4 * pc_x[k] * smd1_309[k];

        t_310[k] = f_10 * smp_128[k]
                   + f_3 * pc_y[k] * snp_155[k];

        t_311[k] = pb_x[k] * smd0_311[k]
                   - f_4 * pc_x[k] * smd1_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pb_x, pc_x, smd0_312, smd0_315, smp_156, \
                         smp_157, smp_158, smd1_312, smd1_315, snp_157, \
                         snp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pb_x[k] * smd0_312[k]
                   + f_8 * smp_156[k]
                   - f_4 * pc_x[k] * smd1_312[k];

        t_313[k] = f_6 * smp_157[k]
                   + f_3 * pc_x[k] * snp_157[k];

        t_314[k] = f_6 * smp_158[k]
                   + f_3 * pc_x[k] * snp_158[k];

        t_315[k] = pb_x[k] * smd0_315[k]
                   - f_4 * pc_x[k] * smd1_315[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_x, pb_y, pc_x, pc_y, smd0_264, \
                         smd0_317, smp_131, smp_160, smd1_264, smd1_317, snp_158, \
                         snp_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_8 * smp_131[k]
                   + f_3 * pc_y[k] * snp_158[k];

        t_317[k] = pb_x[k] * smd0_317[k]
                   - f_4 * pc_x[k] * smd1_317[k];

        t_318[k] = pb_y[k] * smd0_264[k]
                   - f_4 * pc_y[k] * smd1_264[k];

        t_319[k] = f_6 * smp_160[k]
                   + f_3 * pc_x[k] * snp_160[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pb_x, pc_x, pc_y, smd0_321, smd0_323, \
                         smp_134, smp_161, smd1_321, smd1_323, \
                         snp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_6 * smp_161[k]
                   + f_3 * pc_x[k] * snp_161[k];

        t_321[k] = pb_x[k] * smd0_321[k]
                   - f_4 * pc_x[k] * smd1_321[k];

        t_322[k] = f_6 * smp_134[k]
                   + f_3 * pc_y[k] * snp_161[k];

        t_323[k] = pb_x[k] * smd0_323[k]
                   - f_4 * pc_x[k] * smd1_323[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pb_x, pc_x, smd0_324, smd0_327, smp_162, \
                         smp_163, smp_164, smd1_324, smd1_327, snp_163, \
                         snp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = pb_x[k] * smd0_324[k]
                   + f_8 * smp_162[k]
                   - f_4 * pc_x[k] * smd1_324[k];

        t_325[k] = f_6 * smp_163[k]
                   + f_3 * pc_x[k] * snp_163[k];

        t_326[k] = f_6 * smp_164[k]
                   + f_3 * pc_x[k] * snp_164[k];

        t_327[k] = pb_x[k] * smd0_327[k]
                   - f_4 * pc_x[k] * smd1_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, t_332, pb_x, pc_x, pc_y, smd0_329, \
                         smd1_329, sns0_55, sns1_55, snp_164, snp_165, snp_166, \
                         snp_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_3 * pc_y[k] * snp_164[k];

        t_329[k] = pb_x[k] * smd0_329[k]
                   - f_4 * pc_x[k] * smd1_329[k];

        t_330[k] = f_1 * sns0_55[k]
                   - f_2 * sns1_55[k]
                   + f_3 * pc_x[k] * snp_165[k];

        t_331[k] = f_3 * pc_x[k] * snp_166[k];

        t_332[k] = f_3 * pc_x[k] * snp_167[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pb_z, pc_y, pc_z, smd0_270, smp_136, \
                         smp_137, smd1_270, sns0_55, sns1_55, snp_166, \
                         snp_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_0 * smp_136[k]
                   + f_1 * sns0_55[k]
                   - f_2 * sns1_55[k]
                   + f_3 * pc_y[k] * snp_166[k];

        t_334[k] = f_0 * smp_137[k]
                   + f_3 * pc_y[k] * snp_167[k];

        t_335[k] = f_1 * sns0_55[k]
                   - f_2 * sns1_55[k]
                   + f_3 * pc_z[k] * snp_167[k];

        t_336[k] = pb_z[k] * smd0_270[k]
                   - f_4 * pc_z[k] * smd1_270[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pb_z, pc_x, pc_y, pc_z, smd0_273, \
                         smp_140, smd1_273, snp_169, snp_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_3 * pc_x[k] * snp_169[k];

        t_338[k] = f_3 * pc_x[k] * snp_170[k];

        t_339[k] = pb_z[k] * smd0_273[k]
                   - f_4 * pc_z[k] * smd1_273[k];

        t_340[k] = f_5 * smp_140[k]
                   + f_3 * pc_y[k] * snp_170[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_z, smp_137, sns0_56, sns0_57, \
                         sns1_56, sns1_57, snp_170, snp_171, snp_172, \
                         snp_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_6 * smp_137[k]
                   + f_1 * sns0_56[k]
                   - f_2 * sns1_56[k]
                   + f_3 * pc_z[k] * snp_170[k];

        t_342[k] = f_1 * sns0_57[k]
                   - f_2 * sns1_57[k]
                   + f_3 * pc_x[k] * snp_171[k];

        t_343[k] = f_3 * pc_x[k] * snp_172[k];

        t_344[k] = f_3 * pc_x[k] * snp_173[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_y, pc_z, smp_140, smp_142, smp_143, sns0_57, \
                         sns1_57, snp_172, snp_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_7 * smp_142[k]
                   + f_1 * sns0_57[k]
                   - f_2 * sns1_57[k]
                   + f_3 * pc_y[k] * snp_172[k];

        t_346[k] = f_7 * smp_143[k]
                   + f_3 * pc_y[k] * snp_173[k];

        t_347[k] = f_8 * smp_140[k]
                   + f_1 * sns0_57[k]
                   - f_2 * sns1_57[k]
                   + f_3 * pc_z[k] * snp_173[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, t_352, pc_x, pc_y, smp_145, smp_146, \
                         sns0_58, sns1_58, snp_174, snp_175, snp_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_1 * sns0_58[k]
                   - f_2 * sns1_58[k]
                   + f_3 * pc_x[k] * snp_174[k];

        t_349[k] = f_3 * pc_x[k] * snp_175[k];

        t_350[k] = f_3 * pc_x[k] * snp_176[k];

        t_351[k] = f_9 * smp_145[k]
                   + f_1 * sns0_58[k]
                   - f_2 * sns1_58[k]
                   + f_3 * pc_y[k] * snp_175[k];

        t_352[k] = f_9 * smp_146[k]
                   + f_3 * pc_y[k] * snp_176[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pc_x, pc_z, smp_143, sns0_58, sns0_59, \
                         sns1_58, sns1_59, snp_176, snp_177, snp_178, \
                         snp_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_10 * smp_143[k]
                   + f_1 * sns0_58[k]
                   - f_2 * sns1_58[k]
                   + f_3 * pc_z[k] * snp_176[k];

        t_354[k] = f_1 * sns0_59[k]
                   - f_2 * sns1_59[k]
                   + f_3 * pc_x[k] * snp_177[k];

        t_355[k] = f_3 * pc_x[k] * snp_178[k];

        t_356[k] = f_3 * pc_x[k] * snp_179[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, smp_146, smp_148, smp_149, sns0_59, \
                         sns1_59, snp_178, snp_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_11 * smp_148[k]
                   + f_1 * sns0_59[k]
                   - f_2 * sns1_59[k]
                   + f_3 * pc_y[k] * snp_178[k];

        t_358[k] = f_11 * smp_149[k]
                   + f_3 * pc_y[k] * snp_179[k];

        t_359[k] = f_12 * smp_146[k]
                   + f_1 * sns0_59[k]
                   - f_2 * sns1_59[k]
                   + f_3 * pc_z[k] * snp_179[k];
    }
}

static auto
compute_prim_snd_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smd0,
                                                          const size_t smp, const size_t smd1,
                                                          const size_t sns0, const size_t sns1,
                                                          const size_t snp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);
    auto *t_387 = buffer.data(target + 387);
    auto *t_388 = buffer.data(target + 388);
    auto *t_389 = buffer.data(target + 389);
    auto *t_390 = buffer.data(target + 390);
    auto *t_391 = buffer.data(target + 391);
    auto *t_392 = buffer.data(target + 392);
    auto *t_393 = buffer.data(target + 393);
    auto *t_394 = buffer.data(target + 394);
    auto *t_395 = buffer.data(target + 395);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smd0_324 = buffer.data(smd0 + 324);
    const auto *smd0_327 = buffer.data(smd0 + 327);
    const auto *smd0_329 = buffer.data(smd0 + 329);

    const auto *smp_149 = buffer.data(smp + 149);
    const auto *smp_151 = buffer.data(smp + 151);
    const auto *smp_152 = buffer.data(smp + 152);
    const auto *smp_154 = buffer.data(smp + 154);
    const auto *smp_155 = buffer.data(smp + 155);
    const auto *smp_157 = buffer.data(smp + 157);
    const auto *smp_158 = buffer.data(smp + 158);
    const auto *smp_160 = buffer.data(smp + 160);
    const auto *smp_161 = buffer.data(smp + 161);
    const auto *smp_163 = buffer.data(smp + 163);
    const auto *smp_164 = buffer.data(smp + 164);

    const auto *smd1_324 = buffer.data(smd1 + 324);
    const auto *smd1_327 = buffer.data(smd1 + 327);
    const auto *smd1_329 = buffer.data(smd1 + 329);

    const auto *sns0_60 = buffer.data(sns0 + 60);
    const auto *sns0_61 = buffer.data(sns0 + 61);
    const auto *sns0_62 = buffer.data(sns0 + 62);
    const auto *sns0_63 = buffer.data(sns0 + 63);
    const auto *sns0_65 = buffer.data(sns0 + 65);

    const auto *sns1_60 = buffer.data(sns1 + 60);
    const auto *sns1_61 = buffer.data(sns1 + 61);
    const auto *sns1_62 = buffer.data(sns1 + 62);
    const auto *sns1_63 = buffer.data(sns1 + 63);
    const auto *sns1_65 = buffer.data(sns1 + 65);

    const auto *snp_180 = buffer.data(snp + 180);
    const auto *snp_181 = buffer.data(snp + 181);
    const auto *snp_182 = buffer.data(snp + 182);
    const auto *snp_183 = buffer.data(snp + 183);
    const auto *snp_184 = buffer.data(snp + 184);
    const auto *snp_185 = buffer.data(snp + 185);
    const auto *snp_186 = buffer.data(snp + 186);
    const auto *snp_187 = buffer.data(snp + 187);
    const auto *snp_188 = buffer.data(snp + 188);
    const auto *snp_189 = buffer.data(snp + 189);
    const auto *snp_190 = buffer.data(snp + 190);
    const auto *snp_191 = buffer.data(snp + 191);
    const auto *snp_193 = buffer.data(snp + 193);
    const auto *snp_194 = buffer.data(snp + 194);
    const auto *snp_195 = buffer.data(snp + 195);
    const auto *snp_196 = buffer.data(snp + 196);
    const auto *snp_197 = buffer.data(snp + 197);

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pc_x, pc_y, smp_151, smp_152, \
                         sns0_60, sns1_60, snp_180, snp_181, snp_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * sns0_60[k]
                   - f_2 * sns1_60[k]
                   + f_3 * pc_x[k] * snp_180[k];

        t_361[k] = f_3 * pc_x[k] * snp_181[k];

        t_362[k] = f_3 * pc_x[k] * snp_182[k];

        t_363[k] = f_13 * smp_151[k]
                   + f_1 * sns0_60[k]
                   - f_2 * sns1_60[k]
                   + f_3 * pc_y[k] * snp_181[k];

        t_364[k] = f_13 * smp_152[k]
                   + f_3 * pc_y[k] * snp_182[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_z, smp_149, sns0_60, sns0_61, \
                         sns1_60, sns1_61, snp_182, snp_183, snp_184, \
                         snp_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_13 * smp_149[k]
                   + f_1 * sns0_60[k]
                   - f_2 * sns1_60[k]
                   + f_3 * pc_z[k] * snp_182[k];

        t_366[k] = f_1 * sns0_61[k]
                   - f_2 * sns1_61[k]
                   + f_3 * pc_x[k] * snp_183[k];

        t_367[k] = f_3 * pc_x[k] * snp_184[k];

        t_368[k] = f_3 * pc_x[k] * snp_185[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pc_y, pc_z, smp_152, smp_154, smp_155, sns0_61, \
                         sns1_61, snp_184, snp_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_12 * smp_154[k]
                   + f_1 * sns0_61[k]
                   - f_2 * sns1_61[k]
                   + f_3 * pc_y[k] * snp_184[k];

        t_370[k] = f_12 * smp_155[k]
                   + f_3 * pc_y[k] * snp_185[k];

        t_371[k] = f_11 * smp_152[k]
                   + f_1 * sns0_61[k]
                   - f_2 * sns1_61[k]
                   + f_3 * pc_z[k] * snp_185[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, t_376, pc_x, pc_y, smp_157, smp_158, \
                         sns0_62, sns1_62, snp_186, snp_187, snp_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_1 * sns0_62[k]
                   - f_2 * sns1_62[k]
                   + f_3 * pc_x[k] * snp_186[k];

        t_373[k] = f_3 * pc_x[k] * snp_187[k];

        t_374[k] = f_3 * pc_x[k] * snp_188[k];

        t_375[k] = f_10 * smp_157[k]
                   + f_1 * sns0_62[k]
                   - f_2 * sns1_62[k]
                   + f_3 * pc_y[k] * snp_187[k];

        t_376[k] = f_10 * smp_158[k]
                   + f_3 * pc_y[k] * snp_188[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pc_x, pc_z, smp_155, sns0_62, sns0_63, \
                         sns1_62, sns1_63, snp_188, snp_189, snp_190, \
                         snp_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_9 * smp_155[k]
                   + f_1 * sns0_62[k]
                   - f_2 * sns1_62[k]
                   + f_3 * pc_z[k] * snp_188[k];

        t_378[k] = f_1 * sns0_63[k]
                   - f_2 * sns1_63[k]
                   + f_3 * pc_x[k] * snp_189[k];

        t_379[k] = f_3 * pc_x[k] * snp_190[k];

        t_380[k] = f_3 * pc_x[k] * snp_191[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pb_y, pc_y, pc_z, smd0_324, smp_158, \
                         smp_160, smp_161, smd1_324, sns0_63, sns1_63, snp_190, \
                         snp_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_8 * smp_160[k]
                   + f_1 * sns0_63[k]
                   - f_2 * sns1_63[k]
                   + f_3 * pc_y[k] * snp_190[k];

        t_382[k] = f_8 * smp_161[k]
                   + f_3 * pc_y[k] * snp_191[k];

        t_383[k] = f_7 * smp_158[k]
                   + f_1 * sns0_63[k]
                   - f_2 * sns1_63[k]
                   + f_3 * pc_z[k] * snp_191[k];

        t_384[k] = pb_y[k] * smd0_324[k]
                   - f_4 * pc_y[k] * smd1_324[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, t_389, pb_y, pc_x, pc_y, smd0_327, \
                         smd0_329, smp_163, smp_164, smd1_327, smd1_329, snp_193, \
                         snp_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_3 * pc_x[k] * snp_193[k];

        t_386[k] = f_3 * pc_x[k] * snp_194[k];

        t_387[k] = pb_y[k] * smd0_327[k]
                   + f_8 * smp_163[k]
                   - f_4 * pc_y[k] * smd1_327[k];

        t_388[k] = f_6 * smp_164[k]
                   + f_3 * pc_y[k] * snp_194[k];

        t_389[k] = pb_y[k] * smd0_329[k]
                   - f_4 * pc_y[k] * smd1_329[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, t_394, t_395, pc_x, pc_y, pc_z, smp_164, \
                         sns0_65, sns1_65, snp_195, snp_196, snp_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_1 * sns0_65[k]
                   - f_2 * sns1_65[k]
                   + f_3 * pc_x[k] * snp_195[k];

        t_391[k] = f_3 * pc_x[k] * snp_196[k];

        t_392[k] = f_3 * pc_x[k] * snp_197[k];

        t_393[k] = f_1 * sns0_65[k]
                   - f_2 * sns1_65[k]
                   + f_3 * pc_y[k] * snp_196[k];

        t_394[k] = f_3 * pc_y[k] * snp_197[k];

        t_395[k] = f_0 * smp_164[k]
                   + f_1 * sns0_65[k]
                   - f_2 * sns1_65[k]
                   + f_3 * pc_z[k] * snp_197[k];
    }
}

auto
compute_prim_snd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smd0, const size_t smp,
                                                   const size_t smd1, const size_t sns0,
                                                   const size_t sns1, const size_t snp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_snd_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smd0, smp,
                                                              smd1, sns0, sns1, snp, ncols,
                                                              gamma, p, q);

    compute_prim_snd_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smd0, smp,
                                                              smd1, sns0, sns1, snp, ncols,
                                                              gamma, p, q);

    compute_prim_snd_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smd0, smp,
                                                              smd1, sns0, sns1, snp, ncols,
                                                              gamma, p, q);

    compute_prim_snd_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smd0, smp,
                                                              smd1, sns0, sns1, snp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
