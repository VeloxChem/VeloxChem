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


#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_0 = buffer.data(smf0 + 0);
    const auto *smf0_6 = buffer.data(smf0 + 6);
    const auto *smf0_9 = buffer.data(smf0 + 9);
    const auto *smf0_16 = buffer.data(smf0 + 16);
    const auto *smf0_20 = buffer.data(smf0 + 20);
    const auto *smf0_29 = buffer.data(smf0 + 29);
    const auto *smf0_30 = buffer.data(smf0 + 30);
    const auto *smf0_36 = buffer.data(smf0 + 36);
    const auto *smf0_50 = buffer.data(smf0 + 50);
    const auto *smf0_59 = buffer.data(smf0 + 59);
    const auto *smf0_60 = buffer.data(smf0 + 60);
    const auto *smf0_66 = buffer.data(smf0 + 66);

    const auto *smd_0 = buffer.data(smd + 0);
    const auto *smd_3 = buffer.data(smd + 3);
    const auto *smd_4 = buffer.data(smd + 4);
    const auto *smd_5 = buffer.data(smd + 5);
    const auto *smd_6 = buffer.data(smd + 6);
    const auto *smd_9 = buffer.data(smd + 9);
    const auto *smd_10 = buffer.data(smd + 10);
    const auto *smd_11 = buffer.data(smd + 11);
    const auto *smd_12 = buffer.data(smd + 12);
    const auto *smd_15 = buffer.data(smd + 15);
    const auto *smd_16 = buffer.data(smd + 16);
    const auto *smd_17 = buffer.data(smd + 17);
    const auto *smd_18 = buffer.data(smd + 18);
    const auto *smd_21 = buffer.data(smd + 21);
    const auto *smd_22 = buffer.data(smd + 22);
    const auto *smd_23 = buffer.data(smd + 23);
    const auto *smd_24 = buffer.data(smd + 24);
    const auto *smd_27 = buffer.data(smd + 27);
    const auto *smd_28 = buffer.data(smd + 28);
    const auto *smd_29 = buffer.data(smd + 29);
    const auto *smd_30 = buffer.data(smd + 30);
    const auto *smd_33 = buffer.data(smd + 33);
    const auto *smd_34 = buffer.data(smd + 34);
    const auto *smd_35 = buffer.data(smd + 35);
    const auto *smd_36 = buffer.data(smd + 36);
    const auto *smd_39 = buffer.data(smd + 39);
    const auto *smd_40 = buffer.data(smd + 40);
    const auto *smd_41 = buffer.data(smd + 41);
    const auto *smd_42 = buffer.data(smd + 42);
    const auto *smd_45 = buffer.data(smd + 45);
    const auto *smd_46 = buffer.data(smd + 46);
    const auto *smd_47 = buffer.data(smd + 47);
    const auto *smd_48 = buffer.data(smd + 48);
    const auto *smd_51 = buffer.data(smd + 51);
    const auto *smd_52 = buffer.data(smd + 52);
    const auto *smd_53 = buffer.data(smd + 53);
    const auto *smd_54 = buffer.data(smd + 54);
    const auto *smd_57 = buffer.data(smd + 57);
    const auto *smd_58 = buffer.data(smd + 58);
    const auto *smd_59 = buffer.data(smd + 59);
    const auto *smd_60 = buffer.data(smd + 60);
    const auto *smd_63 = buffer.data(smd + 63);
    const auto *smd_64 = buffer.data(smd + 64);
    const auto *smd_65 = buffer.data(smd + 65);
    const auto *smd_69 = buffer.data(smd + 69);
    const auto *smd_70 = buffer.data(smd + 70);
    const auto *smd_71 = buffer.data(smd + 71);
    const auto *smd_72 = buffer.data(smd + 72);
    const auto *smd_75 = buffer.data(smd + 75);
    const auto *smd_76 = buffer.data(smd + 76);
    const auto *smd_77 = buffer.data(smd + 77);

    const auto *smf1_0 = buffer.data(smf1 + 0);
    const auto *smf1_6 = buffer.data(smf1 + 6);
    const auto *smf1_9 = buffer.data(smf1 + 9);
    const auto *smf1_16 = buffer.data(smf1 + 16);
    const auto *smf1_20 = buffer.data(smf1 + 20);
    const auto *smf1_29 = buffer.data(smf1 + 29);
    const auto *smf1_30 = buffer.data(smf1 + 30);
    const auto *smf1_36 = buffer.data(smf1 + 36);
    const auto *smf1_50 = buffer.data(smf1 + 50);
    const auto *smf1_59 = buffer.data(smf1 + 59);
    const auto *smf1_60 = buffer.data(smf1 + 60);
    const auto *smf1_66 = buffer.data(smf1 + 66);

    const auto *snp0_0 = buffer.data(snp0 + 0);
    const auto *snp0_1 = buffer.data(snp0 + 1);
    const auto *snp0_2 = buffer.data(snp0 + 2);
    const auto *snp0_4 = buffer.data(snp0 + 4);
    const auto *snp0_8 = buffer.data(snp0 + 8);
    const auto *snp0_9 = buffer.data(snp0 + 9);
    const auto *snp0_10 = buffer.data(snp0 + 10);
    const auto *snp0_11 = buffer.data(snp0 + 11);
    const auto *snp0_15 = buffer.data(snp0 + 15);
    const auto *snp0_16 = buffer.data(snp0 + 16);
    const auto *snp0_17 = buffer.data(snp0 + 17);
    const auto *snp0_18 = buffer.data(snp0 + 18);
    const auto *snp0_19 = buffer.data(snp0 + 19);
    const auto *snp0_20 = buffer.data(snp0 + 20);
    const auto *snp0_23 = buffer.data(snp0 + 23);
    const auto *snp0_25 = buffer.data(snp0 + 25);
    const auto *snp0_27 = buffer.data(snp0 + 27);
    const auto *snp0_28 = buffer.data(snp0 + 28);
    const auto *snp0_29 = buffer.data(snp0 + 29);
    const auto *snp0_30 = buffer.data(snp0 + 30);
    const auto *snp0_31 = buffer.data(snp0 + 31);
    const auto *snp0_32 = buffer.data(snp0 + 32);
    const auto *snp0_35 = buffer.data(snp0 + 35);
    const auto *snp0_36 = buffer.data(snp0 + 36);
    const auto *snp0_37 = buffer.data(snp0 + 37);
    const auto *snp0_38 = buffer.data(snp0 + 38);

    const auto *snp1_0 = buffer.data(snp1 + 0);
    const auto *snp1_1 = buffer.data(snp1 + 1);
    const auto *snp1_2 = buffer.data(snp1 + 2);
    const auto *snp1_4 = buffer.data(snp1 + 4);
    const auto *snp1_8 = buffer.data(snp1 + 8);
    const auto *snp1_9 = buffer.data(snp1 + 9);
    const auto *snp1_10 = buffer.data(snp1 + 10);
    const auto *snp1_11 = buffer.data(snp1 + 11);
    const auto *snp1_15 = buffer.data(snp1 + 15);
    const auto *snp1_16 = buffer.data(snp1 + 16);
    const auto *snp1_17 = buffer.data(snp1 + 17);
    const auto *snp1_18 = buffer.data(snp1 + 18);
    const auto *snp1_19 = buffer.data(snp1 + 19);
    const auto *snp1_20 = buffer.data(snp1 + 20);
    const auto *snp1_23 = buffer.data(snp1 + 23);
    const auto *snp1_25 = buffer.data(snp1 + 25);
    const auto *snp1_27 = buffer.data(snp1 + 27);
    const auto *snp1_28 = buffer.data(snp1 + 28);
    const auto *snp1_29 = buffer.data(snp1 + 29);
    const auto *snp1_30 = buffer.data(snp1 + 30);
    const auto *snp1_31 = buffer.data(snp1 + 31);
    const auto *snp1_32 = buffer.data(snp1 + 32);
    const auto *snp1_35 = buffer.data(snp1 + 35);
    const auto *snp1_36 = buffer.data(snp1 + 36);
    const auto *snp1_37 = buffer.data(snp1 + 37);
    const auto *snp1_38 = buffer.data(snp1 + 38);

    const auto *snd_0 = buffer.data(snd + 0);
    const auto *snd_3 = buffer.data(snd + 3);
    const auto *snd_4 = buffer.data(snd + 4);
    const auto *snd_5 = buffer.data(snd + 5);
    const auto *snd_6 = buffer.data(snd + 6);
    const auto *snd_9 = buffer.data(snd + 9);
    const auto *snd_10 = buffer.data(snd + 10);
    const auto *snd_11 = buffer.data(snd + 11);
    const auto *snd_12 = buffer.data(snd + 12);
    const auto *snd_15 = buffer.data(snd + 15);
    const auto *snd_16 = buffer.data(snd + 16);
    const auto *snd_17 = buffer.data(snd + 17);
    const auto *snd_18 = buffer.data(snd + 18);
    const auto *snd_21 = buffer.data(snd + 21);
    const auto *snd_22 = buffer.data(snd + 22);
    const auto *snd_23 = buffer.data(snd + 23);
    const auto *snd_24 = buffer.data(snd + 24);
    const auto *snd_27 = buffer.data(snd + 27);
    const auto *snd_28 = buffer.data(snd + 28);
    const auto *snd_29 = buffer.data(snd + 29);
    const auto *snd_30 = buffer.data(snd + 30);
    const auto *snd_33 = buffer.data(snd + 33);
    const auto *snd_34 = buffer.data(snd + 34);
    const auto *snd_35 = buffer.data(snd + 35);
    const auto *snd_36 = buffer.data(snd + 36);
    const auto *snd_39 = buffer.data(snd + 39);
    const auto *snd_40 = buffer.data(snd + 40);
    const auto *snd_41 = buffer.data(snd + 41);
    const auto *snd_42 = buffer.data(snd + 42);
    const auto *snd_45 = buffer.data(snd + 45);
    const auto *snd_46 = buffer.data(snd + 46);
    const auto *snd_47 = buffer.data(snd + 47);
    const auto *snd_48 = buffer.data(snd + 48);
    const auto *snd_51 = buffer.data(snd + 51);
    const auto *snd_52 = buffer.data(snd + 52);
    const auto *snd_53 = buffer.data(snd + 53);
    const auto *snd_54 = buffer.data(snd + 54);
    const auto *snd_57 = buffer.data(snd + 57);
    const auto *snd_58 = buffer.data(snd + 58);
    const auto *snd_59 = buffer.data(snd + 59);
    const auto *snd_60 = buffer.data(snd + 60);
    const auto *snd_63 = buffer.data(snd + 63);
    const auto *snd_64 = buffer.data(snd + 64);
    const auto *snd_65 = buffer.data(snd + 65);
    const auto *snd_66 = buffer.data(snd + 66);
    const auto *snd_69 = buffer.data(snd + 69);
    const auto *snd_70 = buffer.data(snd + 70);
    const auto *snd_71 = buffer.data(snd + 71);
    const auto *snd_72 = buffer.data(snd + 72);
    const auto *snd_75 = buffer.data(snd + 75);
    const auto *snd_76 = buffer.data(snd + 76);
    const auto *snd_77 = buffer.data(snd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, smd_0, smd_3, smd_4, \
                         snp0_0, snp1_0, snd_0, snd_3, snd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * smd_0[k]
                 + f_1 * snp0_0[k]
                 - f_2 * snp1_0[k]
                 + f_3 * pc_x[k] * snd_0[k];

        t_1[k] = f_3 * pc_y[k] * snd_0[k];

        t_2[k] = f_3 * pc_z[k] * snd_0[k];

        t_3[k] = f_0 * smd_3[k]
                 + f_3 * pc_x[k] * snd_3[k];

        t_4[k] = f_0 * smd_4[k]
                 + f_3 * pc_x[k] * snd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, smd_5, snp0_1, snp0_2, \
                         snp1_1, snp1_2, snd_3, snd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * smd_5[k]
                 + f_3 * pc_x[k] * snd_5[k];

        t_6[k] = f_1 * snp0_1[k]
                 - f_2 * snp1_1[k]
                 + f_3 * pc_y[k] * snd_3[k];

        t_7[k] = f_3 * pc_z[k] * snd_3[k];

        t_8[k] = f_3 * pc_y[k] * snd_5[k];

        t_9[k] = f_1 * snp0_2[k]
                 - f_2 * snp1_2[k]
                 + f_3 * pc_z[k] * snd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, smf0_0, smd_0, smd_9, \
                         smf1_0, snd_6, snd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * smf0_0[k]
                  - f_4 * pc_y[k] * smf1_0[k];

        t_11[k] = f_5 * smd_0[k]
                  + f_3 * pc_y[k] * snd_6[k];

        t_12[k] = f_3 * pc_z[k] * snd_6[k];

        t_13[k] = f_6 * smd_9[k]
                  + f_3 * pc_x[k] * snd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, smd_3, smd_10, smd_11, \
                         snp0_4, snp1_4, snd_9, snd_10, snd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * smd_10[k]
                  + f_3 * pc_x[k] * snd_10[k];

        t_15[k] = f_6 * smd_11[k]
                  + f_3 * pc_x[k] * snd_11[k];

        t_16[k] = f_5 * smd_3[k]
                  + f_1 * snp0_4[k]
                  - f_2 * snp1_4[k]
                  + f_3 * pc_y[k] * snd_9[k];

        t_17[k] = f_3 * pc_z[k] * snd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, smf0_0, smf0_9, \
                         smd_5, smf1_0, smf1_9, snd_11, snd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * smd_5[k]
                  + f_3 * pc_y[k] * snd_11[k];

        t_19[k] = pb_y[k] * smf0_9[k]
                  - f_4 * pc_y[k] * smf1_9[k];

        t_20[k] = pb_z[k] * smf0_0[k]
                  - f_4 * pc_z[k] * smf1_0[k];

        t_21[k] = f_3 * pc_y[k] * snd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, smd_0, smd_15, smd_16, smd_17, \
                         snd_12, snd_15, snd_16, snd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * smd_0[k]
                  + f_3 * pc_z[k] * snd_12[k];

        t_23[k] = f_6 * smd_15[k]
                  + f_3 * pc_x[k] * snd_15[k];

        t_24[k] = f_6 * smd_16[k]
                  + f_3 * pc_x[k] * snd_16[k];

        t_25[k] = f_6 * smd_17[k]
                  + f_3 * pc_x[k] * snd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, smf0_6, smd_3, smd_5, \
                         smf1_6, snp0_8, snp1_8, snd_15, snd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * smf0_6[k]
                  - f_4 * pc_z[k] * smf1_6[k];

        t_27[k] = f_5 * smd_3[k]
                  + f_3 * pc_z[k] * snd_15[k];

        t_28[k] = f_3 * pc_y[k] * snd_17[k];

        t_29[k] = f_5 * smd_5[k]
                  + f_1 * snp0_8[k]
                  - f_2 * snp1_8[k]
                  + f_3 * pc_z[k] * snd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, smd_6, smd_18, smd_21, \
                         snp0_9, snp1_9, snd_18, snd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * smd_18[k]
                  + f_1 * snp0_9[k]
                  - f_2 * snp1_9[k]
                  + f_3 * pc_x[k] * snd_18[k];

        t_31[k] = f_8 * smd_6[k]
                  + f_3 * pc_y[k] * snd_18[k];

        t_32[k] = f_3 * pc_z[k] * snd_18[k];

        t_33[k] = f_7 * smd_21[k]
                  + f_3 * pc_x[k] * snd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, smd_9, smd_22, smd_23, \
                         snp0_10, snp1_10, snd_21, snd_22, snd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * smd_22[k]
                  + f_3 * pc_x[k] * snd_22[k];

        t_35[k] = f_7 * smd_23[k]
                  + f_3 * pc_x[k] * snd_23[k];

        t_36[k] = f_8 * smd_9[k]
                  + f_1 * snp0_10[k]
                  - f_2 * snp1_10[k]
                  + f_3 * pc_y[k] * snd_21[k];

        t_37[k] = f_3 * pc_z[k] * snd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, smf0_20, smd_11, smd_12, \
                         smf1_20, snp0_11, snp1_11, snd_23, snd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * smd_11[k]
                  + f_3 * pc_y[k] * snd_23[k];

        t_39[k] = f_1 * snp0_11[k]
                  - f_2 * snp1_11[k]
                  + f_3 * pc_z[k] * snd_23[k];

        t_40[k] = pb_y[k] * smf0_20[k]
                  - f_4 * pc_y[k] * smf1_20[k];

        t_41[k] = f_5 * smd_12[k]
                  + f_3 * pc_y[k] * snd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, smd_6, smd_27, smd_28, smd_29, \
                         snd_24, snd_27, snd_28, snd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * smd_6[k]
                  + f_3 * pc_z[k] * snd_24[k];

        t_43[k] = f_7 * smd_27[k]
                  + f_3 * pc_x[k] * snd_27[k];

        t_44[k] = f_7 * smd_28[k]
                  + f_3 * pc_x[k] * snd_28[k];

        t_45[k] = f_7 * smd_29[k]
                  + f_3 * pc_x[k] * snd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, smf0_16, smf0_29, \
                         smd_9, smd_17, smf1_16, smf1_29, snd_27, \
                         snd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * smf0_16[k]
                  - f_4 * pc_z[k] * smf1_16[k];

        t_47[k] = f_5 * smd_9[k]
                  + f_3 * pc_z[k] * snd_27[k];

        t_48[k] = f_5 * smd_17[k]
                  + f_3 * pc_y[k] * snd_29[k];

        t_49[k] = pb_y[k] * smf0_29[k]
                  - f_4 * pc_y[k] * smf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, smd_12, smd_30, smd_33, \
                         snp0_15, snp1_15, snd_30, snd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * smd_30[k]
                  + f_1 * snp0_15[k]
                  - f_2 * snp1_15[k]
                  + f_3 * pc_x[k] * snd_30[k];

        t_51[k] = f_3 * pc_y[k] * snd_30[k];

        t_52[k] = f_8 * smd_12[k]
                  + f_3 * pc_z[k] * snd_30[k];

        t_53[k] = f_7 * smd_33[k]
                  + f_3 * pc_x[k] * snd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, smd_15, smd_34, \
                         smd_35, snp0_16, snp1_16, snd_33, snd_34, \
                         snd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * smd_34[k]
                  + f_3 * pc_x[k] * snd_34[k];

        t_55[k] = f_7 * smd_35[k]
                  + f_3 * pc_x[k] * snd_35[k];

        t_56[k] = f_1 * snp0_16[k]
                  - f_2 * snp1_16[k]
                  + f_3 * pc_y[k] * snd_33[k];

        t_57[k] = f_8 * smd_15[k]
                  + f_3 * pc_z[k] * snd_33[k];

        t_58[k] = f_3 * pc_y[k] * snd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, smd_17, smd_18, smd_36, \
                         snp0_17, snp0_18, snp1_17, snp1_18, snd_35, \
                         snd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * smd_17[k]
                  + f_1 * snp0_17[k]
                  - f_2 * snp1_17[k]
                  + f_3 * pc_z[k] * snd_35[k];

        t_60[k] = f_9 * smd_36[k]
                  + f_1 * snp0_18[k]
                  - f_2 * snp1_18[k]
                  + f_3 * pc_x[k] * snd_36[k];

        t_61[k] = f_10 * smd_18[k]
                  + f_3 * pc_y[k] * snd_36[k];

        t_62[k] = f_3 * pc_z[k] * snd_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, smd_21, smd_39, smd_40, smd_41, \
                         snp0_19, snp1_19, snd_39, snd_40, snd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * smd_39[k]
                  + f_3 * pc_x[k] * snd_39[k];

        t_64[k] = f_9 * smd_40[k]
                  + f_3 * pc_x[k] * snd_40[k];

        t_65[k] = f_9 * smd_41[k]
                  + f_3 * pc_x[k] * snd_41[k];

        t_66[k] = f_10 * smd_21[k]
                  + f_1 * snp0_19[k]
                  - f_2 * snp1_19[k]
                  + f_3 * pc_y[k] * snd_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, smf0_30, smd_23, smf1_30, \
                         snp0_20, snp1_20, snd_39, snd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * snd_39[k];

        t_68[k] = f_10 * smd_23[k]
                  + f_3 * pc_y[k] * snd_41[k];

        t_69[k] = f_1 * snp0_20[k]
                  - f_2 * snp1_20[k]
                  + f_3 * pc_z[k] * snd_41[k];

        t_70[k] = pb_z[k] * smf0_30[k]
                  - f_4 * pc_z[k] * smf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, smd_18, smd_24, smd_45, \
                         smd_46, snd_42, snd_45, snd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * smd_24[k]
                  + f_3 * pc_y[k] * snd_42[k];

        t_72[k] = f_5 * smd_18[k]
                  + f_3 * pc_z[k] * snd_42[k];

        t_73[k] = f_9 * smd_45[k]
                  + f_3 * pc_x[k] * snd_45[k];

        t_74[k] = f_9 * smd_46[k]
                  + f_3 * pc_x[k] * snd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, smf0_36, smd_21, \
                         smd_29, smd_47, smf1_36, snd_45, snd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * smd_47[k]
                  + f_3 * pc_x[k] * snd_47[k];

        t_76[k] = pb_z[k] * smf0_36[k]
                  - f_4 * pc_z[k] * smf1_36[k];

        t_77[k] = f_5 * smd_21[k]
                  + f_3 * pc_z[k] * snd_45[k];

        t_78[k] = f_8 * smd_29[k]
                  + f_3 * pc_y[k] * snd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, smf0_50, smd_23, smd_24, \
                         smd_30, smf1_50, snp0_23, snp1_23, snd_47, \
                         snd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * smd_23[k]
                  + f_1 * snp0_23[k]
                  - f_2 * snp1_23[k]
                  + f_3 * pc_z[k] * snd_47[k];

        t_80[k] = pb_y[k] * smf0_50[k]
                  - f_4 * pc_y[k] * smf1_50[k];

        t_81[k] = f_5 * smd_30[k]
                  + f_3 * pc_y[k] * snd_48[k];

        t_82[k] = f_8 * smd_24[k]
                  + f_3 * pc_z[k] * snd_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, smd_33, smd_51, smd_52, smd_53, \
                         snp0_25, snp1_25, snd_51, snd_52, snd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * smd_51[k]
                  + f_3 * pc_x[k] * snd_51[k];

        t_84[k] = f_9 * smd_52[k]
                  + f_3 * pc_x[k] * snd_52[k];

        t_85[k] = f_9 * smd_53[k]
                  + f_3 * pc_x[k] * snd_53[k];

        t_86[k] = f_5 * smd_33[k]
                  + f_1 * snp0_25[k]
                  - f_2 * snp1_25[k]
                  + f_3 * pc_y[k] * snd_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, smf0_59, smd_27, smd_35, smf1_59, \
                         snd_51, snd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * smd_27[k]
                  + f_3 * pc_z[k] * snd_51[k];

        t_88[k] = f_5 * smd_35[k]
                  + f_3 * pc_y[k] * snd_53[k];

        t_89[k] = pb_y[k] * smf0_59[k]
                  - f_4 * pc_y[k] * smf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, smd_30, smd_54, smd_57, \
                         snp0_27, snp1_27, snd_54, snd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * smd_54[k]
                  + f_1 * snp0_27[k]
                  - f_2 * snp1_27[k]
                  + f_3 * pc_x[k] * snd_54[k];

        t_91[k] = f_3 * pc_y[k] * snd_54[k];

        t_92[k] = f_10 * smd_30[k]
                  + f_3 * pc_z[k] * snd_54[k];

        t_93[k] = f_9 * smd_57[k]
                  + f_3 * pc_x[k] * snd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, smd_33, smd_58, \
                         smd_59, snp0_28, snp1_28, snd_57, snd_58, \
                         snd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * smd_58[k]
                  + f_3 * pc_x[k] * snd_58[k];

        t_95[k] = f_9 * smd_59[k]
                  + f_3 * pc_x[k] * snd_59[k];

        t_96[k] = f_1 * snp0_28[k]
                  - f_2 * snp1_28[k]
                  + f_3 * pc_y[k] * snd_57[k];

        t_97[k] = f_10 * smd_33[k]
                  + f_3 * pc_z[k] * snd_57[k];

        t_98[k] = f_3 * pc_y[k] * snd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, smd_35, smd_36, smd_60, \
                         snp0_29, snp0_30, snp1_29, snp1_30, snd_59, \
                         snd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * smd_35[k]
                  + f_1 * snp0_29[k]
                  - f_2 * snp1_29[k]
                  + f_3 * pc_z[k] * snd_59[k];

        t_100[k] = f_11 * smd_60[k]
                   + f_1 * snp0_30[k]
                   - f_2 * snp1_30[k]
                   + f_3 * pc_x[k] * snd_60[k];

        t_101[k] = f_12 * smd_36[k]
                   + f_3 * pc_y[k] * snd_60[k];

        t_102[k] = f_3 * pc_z[k] * snd_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, smd_39, smd_63, smd_64, \
                         smd_65, snp0_31, snp1_31, snd_63, snd_64, \
                         snd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * smd_63[k]
                   + f_3 * pc_x[k] * snd_63[k];

        t_104[k] = f_11 * smd_64[k]
                   + f_3 * pc_x[k] * snd_64[k];

        t_105[k] = f_11 * smd_65[k]
                   + f_3 * pc_x[k] * snd_65[k];

        t_106[k] = f_12 * smd_39[k]
                   + f_1 * snp0_31[k]
                   - f_2 * snp1_31[k]
                   + f_3 * pc_y[k] * snd_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, smf0_60, smd_41, \
                         smf1_60, snp0_32, snp1_32, snd_63, snd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * snd_63[k];

        t_108[k] = f_12 * smd_41[k]
                   + f_3 * pc_y[k] * snd_65[k];

        t_109[k] = f_1 * snp0_32[k]
                   - f_2 * snp1_32[k]
                   + f_3 * pc_z[k] * snd_65[k];

        t_110[k] = pb_z[k] * smf0_60[k]
                   - f_4 * pc_z[k] * smf1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, smd_36, smd_42, smd_69, \
                         smd_70, snd_66, snd_69, snd_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * smd_42[k]
                   + f_3 * pc_y[k] * snd_66[k];

        t_112[k] = f_5 * smd_36[k]
                   + f_3 * pc_z[k] * snd_66[k];

        t_113[k] = f_11 * smd_69[k]
                   + f_3 * pc_x[k] * snd_69[k];

        t_114[k] = f_11 * smd_70[k]
                   + f_3 * pc_x[k] * snd_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, smf0_66, smd_39, \
                         smd_47, smd_71, smf1_66, snd_69, snd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_11 * smd_71[k]
                   + f_3 * pc_x[k] * snd_71[k];

        t_116[k] = pb_z[k] * smf0_66[k]
                   - f_4 * pc_z[k] * smf1_66[k];

        t_117[k] = f_5 * smd_39[k]
                   + f_3 * pc_z[k] * snd_69[k];

        t_118[k] = f_10 * smd_47[k]
                   + f_3 * pc_y[k] * snd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, smd_41, smd_48, smd_72, \
                         snp0_35, snp0_36, snp1_35, snp1_36, snd_71, \
                         snd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * smd_41[k]
                   + f_1 * snp0_35[k]
                   - f_2 * snp1_35[k]
                   + f_3 * pc_z[k] * snd_71[k];

        t_120[k] = f_11 * smd_72[k]
                   + f_1 * snp0_36[k]
                   - f_2 * snp1_36[k]
                   + f_3 * pc_x[k] * snd_72[k];

        t_121[k] = f_8 * smd_48[k]
                   + f_3 * pc_y[k] * snd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, smd_42, smd_75, smd_76, \
                         smd_77, snd_72, snd_75, snd_76, snd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * smd_42[k]
                   + f_3 * pc_z[k] * snd_72[k];

        t_123[k] = f_11 * smd_75[k]
                   + f_3 * pc_x[k] * snd_75[k];

        t_124[k] = f_11 * smd_76[k]
                   + f_3 * pc_x[k] * snd_76[k];

        t_125[k] = f_11 * smd_77[k]
                   + f_3 * pc_x[k] * snd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, smd_45, smd_47, smd_51, \
                         smd_53, snp0_37, snp0_38, snp1_37, snp1_38, snd_75, \
                         snd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * smd_51[k]
                   + f_1 * snp0_37[k]
                   - f_2 * snp1_37[k]
                   + f_3 * pc_y[k] * snd_75[k];

        t_127[k] = f_8 * smd_45[k]
                   + f_3 * pc_z[k] * snd_75[k];

        t_128[k] = f_8 * smd_53[k]
                   + f_3 * pc_y[k] * snd_77[k];

        t_129[k] = f_8 * smd_47[k]
                   + f_1 * snp0_38[k]
                   - f_2 * snp1_38[k]
                   + f_3 * pc_z[k] * snd_77[k];
    }
}

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_90 = buffer.data(smf0 + 90);
    const auto *smf0_99 = buffer.data(smf0 + 99);
    const auto *smf0_100 = buffer.data(smf0 + 100);
    const auto *smf0_106 = buffer.data(smf0 + 106);
    const auto *smf0_140 = buffer.data(smf0 + 140);
    const auto *smf0_149 = buffer.data(smf0 + 149);
    const auto *smf0_150 = buffer.data(smf0 + 150);
    const auto *smf0_156 = buffer.data(smf0 + 156);

    const auto *smd_48 = buffer.data(smd + 48);
    const auto *smd_51 = buffer.data(smd + 51);
    const auto *smd_54 = buffer.data(smd + 54);
    const auto *smd_57 = buffer.data(smd + 57);
    const auto *smd_59 = buffer.data(smd + 59);
    const auto *smd_60 = buffer.data(smd + 60);
    const auto *smd_63 = buffer.data(smd + 63);
    const auto *smd_65 = buffer.data(smd + 65);
    const auto *smd_66 = buffer.data(smd + 66);
    const auto *smd_69 = buffer.data(smd + 69);
    const auto *smd_71 = buffer.data(smd + 71);
    const auto *smd_72 = buffer.data(smd + 72);
    const auto *smd_75 = buffer.data(smd + 75);
    const auto *smd_77 = buffer.data(smd + 77);
    const auto *smd_78 = buffer.data(smd + 78);
    const auto *smd_81 = buffer.data(smd + 81);
    const auto *smd_82 = buffer.data(smd + 82);
    const auto *smd_83 = buffer.data(smd + 83);
    const auto *smd_84 = buffer.data(smd + 84);
    const auto *smd_87 = buffer.data(smd + 87);
    const auto *smd_88 = buffer.data(smd + 88);
    const auto *smd_89 = buffer.data(smd + 89);
    const auto *smd_90 = buffer.data(smd + 90);
    const auto *smd_93 = buffer.data(smd + 93);
    const auto *smd_94 = buffer.data(smd + 94);
    const auto *smd_95 = buffer.data(smd + 95);
    const auto *smd_96 = buffer.data(smd + 96);
    const auto *smd_99 = buffer.data(smd + 99);
    const auto *smd_100 = buffer.data(smd + 100);
    const auto *smd_101 = buffer.data(smd + 101);
    const auto *smd_102 = buffer.data(smd + 102);
    const auto *smd_105 = buffer.data(smd + 105);
    const auto *smd_106 = buffer.data(smd + 106);
    const auto *smd_107 = buffer.data(smd + 107);
    const auto *smd_108 = buffer.data(smd + 108);
    const auto *smd_111 = buffer.data(smd + 111);
    const auto *smd_112 = buffer.data(smd + 112);
    const auto *smd_113 = buffer.data(smd + 113);
    const auto *smd_114 = buffer.data(smd + 114);
    const auto *smd_117 = buffer.data(smd + 117);
    const auto *smd_118 = buffer.data(smd + 118);
    const auto *smd_119 = buffer.data(smd + 119);
    const auto *smd_120 = buffer.data(smd + 120);
    const auto *smd_123 = buffer.data(smd + 123);
    const auto *smd_124 = buffer.data(smd + 124);
    const auto *smd_125 = buffer.data(smd + 125);
    const auto *smd_126 = buffer.data(smd + 126);
    const auto *smd_129 = buffer.data(smd + 129);
    const auto *smd_130 = buffer.data(smd + 130);
    const auto *smd_131 = buffer.data(smd + 131);
    const auto *smd_135 = buffer.data(smd + 135);
    const auto *smd_136 = buffer.data(smd + 136);
    const auto *smd_137 = buffer.data(smd + 137);
    const auto *smd_138 = buffer.data(smd + 138);
    const auto *smd_141 = buffer.data(smd + 141);
    const auto *smd_142 = buffer.data(smd + 142);
    const auto *smd_143 = buffer.data(smd + 143);
    const auto *smd_144 = buffer.data(smd + 144);
    const auto *smd_147 = buffer.data(smd + 147);
    const auto *smd_148 = buffer.data(smd + 148);
    const auto *smd_149 = buffer.data(smd + 149);
    const auto *smd_150 = buffer.data(smd + 150);
    const auto *smd_153 = buffer.data(smd + 153);
    const auto *smd_154 = buffer.data(smd + 154);

    const auto *smf1_90 = buffer.data(smf1 + 90);
    const auto *smf1_99 = buffer.data(smf1 + 99);
    const auto *smf1_100 = buffer.data(smf1 + 100);
    const auto *smf1_106 = buffer.data(smf1 + 106);
    const auto *smf1_140 = buffer.data(smf1 + 140);
    const auto *smf1_149 = buffer.data(smf1 + 149);
    const auto *smf1_150 = buffer.data(smf1 + 150);
    const auto *smf1_156 = buffer.data(smf1 + 156);

    const auto *snp0_40 = buffer.data(snp0 + 40);
    const auto *snp0_42 = buffer.data(snp0 + 42);
    const auto *snp0_43 = buffer.data(snp0 + 43);
    const auto *snp0_44 = buffer.data(snp0 + 44);
    const auto *snp0_45 = buffer.data(snp0 + 45);
    const auto *snp0_46 = buffer.data(snp0 + 46);
    const auto *snp0_47 = buffer.data(snp0 + 47);
    const auto *snp0_50 = buffer.data(snp0 + 50);
    const auto *snp0_51 = buffer.data(snp0 + 51);
    const auto *snp0_52 = buffer.data(snp0 + 52);
    const auto *snp0_53 = buffer.data(snp0 + 53);
    const auto *snp0_54 = buffer.data(snp0 + 54);
    const auto *snp0_55 = buffer.data(snp0 + 55);
    const auto *snp0_56 = buffer.data(snp0 + 56);
    const auto *snp0_58 = buffer.data(snp0 + 58);
    const auto *snp0_60 = buffer.data(snp0 + 60);
    const auto *snp0_61 = buffer.data(snp0 + 61);
    const auto *snp0_62 = buffer.data(snp0 + 62);
    const auto *snp0_63 = buffer.data(snp0 + 63);
    const auto *snp0_64 = buffer.data(snp0 + 64);
    const auto *snp0_65 = buffer.data(snp0 + 65);
    const auto *snp0_68 = buffer.data(snp0 + 68);
    const auto *snp0_69 = buffer.data(snp0 + 69);
    const auto *snp0_70 = buffer.data(snp0 + 70);
    const auto *snp0_71 = buffer.data(snp0 + 71);
    const auto *snp0_72 = buffer.data(snp0 + 72);
    const auto *snp0_73 = buffer.data(snp0 + 73);
    const auto *snp0_74 = buffer.data(snp0 + 74);
    const auto *snp0_75 = buffer.data(snp0 + 75);

    const auto *snp1_40 = buffer.data(snp1 + 40);
    const auto *snp1_42 = buffer.data(snp1 + 42);
    const auto *snp1_43 = buffer.data(snp1 + 43);
    const auto *snp1_44 = buffer.data(snp1 + 44);
    const auto *snp1_45 = buffer.data(snp1 + 45);
    const auto *snp1_46 = buffer.data(snp1 + 46);
    const auto *snp1_47 = buffer.data(snp1 + 47);
    const auto *snp1_50 = buffer.data(snp1 + 50);
    const auto *snp1_51 = buffer.data(snp1 + 51);
    const auto *snp1_52 = buffer.data(snp1 + 52);
    const auto *snp1_53 = buffer.data(snp1 + 53);
    const auto *snp1_54 = buffer.data(snp1 + 54);
    const auto *snp1_55 = buffer.data(snp1 + 55);
    const auto *snp1_56 = buffer.data(snp1 + 56);
    const auto *snp1_58 = buffer.data(snp1 + 58);
    const auto *snp1_60 = buffer.data(snp1 + 60);
    const auto *snp1_61 = buffer.data(snp1 + 61);
    const auto *snp1_62 = buffer.data(snp1 + 62);
    const auto *snp1_63 = buffer.data(snp1 + 63);
    const auto *snp1_64 = buffer.data(snp1 + 64);
    const auto *snp1_65 = buffer.data(snp1 + 65);
    const auto *snp1_68 = buffer.data(snp1 + 68);
    const auto *snp1_69 = buffer.data(snp1 + 69);
    const auto *snp1_70 = buffer.data(snp1 + 70);
    const auto *snp1_71 = buffer.data(snp1 + 71);
    const auto *snp1_72 = buffer.data(snp1 + 72);
    const auto *snp1_73 = buffer.data(snp1 + 73);
    const auto *snp1_74 = buffer.data(snp1 + 74);
    const auto *snp1_75 = buffer.data(snp1 + 75);

    const auto *snd_78 = buffer.data(snd + 78);
    const auto *snd_81 = buffer.data(snd + 81);
    const auto *snd_82 = buffer.data(snd + 82);
    const auto *snd_83 = buffer.data(snd + 83);
    const auto *snd_84 = buffer.data(snd + 84);
    const auto *snd_87 = buffer.data(snd + 87);
    const auto *snd_88 = buffer.data(snd + 88);
    const auto *snd_89 = buffer.data(snd + 89);
    const auto *snd_90 = buffer.data(snd + 90);
    const auto *snd_93 = buffer.data(snd + 93);
    const auto *snd_94 = buffer.data(snd + 94);
    const auto *snd_95 = buffer.data(snd + 95);
    const auto *snd_96 = buffer.data(snd + 96);
    const auto *snd_99 = buffer.data(snd + 99);
    const auto *snd_100 = buffer.data(snd + 100);
    const auto *snd_101 = buffer.data(snd + 101);
    const auto *snd_102 = buffer.data(snd + 102);
    const auto *snd_105 = buffer.data(snd + 105);
    const auto *snd_106 = buffer.data(snd + 106);
    const auto *snd_107 = buffer.data(snd + 107);
    const auto *snd_108 = buffer.data(snd + 108);
    const auto *snd_111 = buffer.data(snd + 111);
    const auto *snd_112 = buffer.data(snd + 112);
    const auto *snd_113 = buffer.data(snd + 113);
    const auto *snd_114 = buffer.data(snd + 114);
    const auto *snd_117 = buffer.data(snd + 117);
    const auto *snd_118 = buffer.data(snd + 118);
    const auto *snd_119 = buffer.data(snd + 119);
    const auto *snd_120 = buffer.data(snd + 120);
    const auto *snd_123 = buffer.data(snd + 123);
    const auto *snd_124 = buffer.data(snd + 124);
    const auto *snd_125 = buffer.data(snd + 125);
    const auto *snd_126 = buffer.data(snd + 126);
    const auto *snd_129 = buffer.data(snd + 129);
    const auto *snd_130 = buffer.data(snd + 130);
    const auto *snd_131 = buffer.data(snd + 131);
    const auto *snd_132 = buffer.data(snd + 132);
    const auto *snd_135 = buffer.data(snd + 135);
    const auto *snd_136 = buffer.data(snd + 136);
    const auto *snd_137 = buffer.data(snd + 137);
    const auto *snd_138 = buffer.data(snd + 138);
    const auto *snd_141 = buffer.data(snd + 141);
    const auto *snd_142 = buffer.data(snd + 142);
    const auto *snd_143 = buffer.data(snd + 143);
    const auto *snd_144 = buffer.data(snd + 144);
    const auto *snd_147 = buffer.data(snd + 147);
    const auto *snd_148 = buffer.data(snd + 148);
    const auto *snd_149 = buffer.data(snd + 149);
    const auto *snd_150 = buffer.data(snd + 150);
    const auto *snd_153 = buffer.data(snd + 153);
    const auto *snd_154 = buffer.data(snd + 154);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, smf0_90, smd_48, \
                         smd_54, smd_81, smf1_90, snd_78, snd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * smf0_90[k]
                   - f_4 * pc_y[k] * smf1_90[k];

        t_131[k] = f_5 * smd_54[k]
                   + f_3 * pc_y[k] * snd_78[k];

        t_132[k] = f_10 * smd_48[k]
                   + f_3 * pc_z[k] * snd_78[k];

        t_133[k] = f_11 * smd_81[k]
                   + f_3 * pc_x[k] * snd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, smd_51, smd_57, smd_82, \
                         smd_83, snp0_40, snp1_40, snd_81, snd_82, \
                         snd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * smd_82[k]
                   + f_3 * pc_x[k] * snd_82[k];

        t_135[k] = f_11 * smd_83[k]
                   + f_3 * pc_x[k] * snd_83[k];

        t_136[k] = f_5 * smd_57[k]
                   + f_1 * snp0_40[k]
                   - f_2 * snp1_40[k]
                   + f_3 * pc_y[k] * snd_81[k];

        t_137[k] = f_10 * smd_51[k]
                   + f_3 * pc_z[k] * snd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, smf0_99, smd_59, \
                         smd_84, smf1_99, snp0_42, snp1_42, snd_83, \
                         snd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * smd_59[k]
                   + f_3 * pc_y[k] * snd_83[k];

        t_139[k] = pb_y[k] * smf0_99[k]
                   - f_4 * pc_y[k] * smf1_99[k];

        t_140[k] = f_11 * smd_84[k]
                   + f_1 * snp0_42[k]
                   - f_2 * snp1_42[k]
                   + f_3 * pc_x[k] * snd_84[k];

        t_141[k] = f_3 * pc_y[k] * snd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, smd_54, smd_87, smd_88, \
                         smd_89, snd_84, snd_87, snd_88, snd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_12 * smd_54[k]
                   + f_3 * pc_z[k] * snd_84[k];

        t_143[k] = f_11 * smd_87[k]
                   + f_3 * pc_x[k] * snd_87[k];

        t_144[k] = f_11 * smd_88[k]
                   + f_3 * pc_x[k] * snd_88[k];

        t_145[k] = f_11 * smd_89[k]
                   + f_3 * pc_x[k] * snd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, smd_57, smd_59, snp0_43, \
                         snp0_44, snp1_43, snp1_44, snd_87, snd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * snp0_43[k]
                   - f_2 * snp1_43[k]
                   + f_3 * pc_y[k] * snd_87[k];

        t_147[k] = f_12 * smd_57[k]
                   + f_3 * pc_z[k] * snd_87[k];

        t_148[k] = f_3 * pc_y[k] * snd_89[k];

        t_149[k] = f_12 * smd_59[k]
                   + f_1 * snp0_44[k]
                   - f_2 * snp1_44[k]
                   + f_3 * pc_z[k] * snd_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, smd_60, smd_90, smd_93, \
                         snp0_45, snp1_45, snd_90, snd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_13 * smd_90[k]
                   + f_1 * snp0_45[k]
                   - f_2 * snp1_45[k]
                   + f_3 * pc_x[k] * snd_90[k];

        t_151[k] = f_13 * smd_60[k]
                   + f_3 * pc_y[k] * snd_90[k];

        t_152[k] = f_3 * pc_z[k] * snd_90[k];

        t_153[k] = f_13 * smd_93[k]
                   + f_3 * pc_x[k] * snd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, smd_63, smd_94, smd_95, \
                         snp0_46, snp1_46, snd_93, snd_94, snd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_13 * smd_94[k]
                   + f_3 * pc_x[k] * snd_94[k];

        t_155[k] = f_13 * smd_95[k]
                   + f_3 * pc_x[k] * snd_95[k];

        t_156[k] = f_13 * smd_63[k]
                   + f_1 * snp0_46[k]
                   - f_2 * snp1_46[k]
                   + f_3 * pc_y[k] * snd_93[k];

        t_157[k] = f_3 * pc_z[k] * snd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_z, pc_y, pc_z, smf0_100, smd_65, \
                         smd_66, smf1_100, snp0_47, snp1_47, snd_95, \
                         snd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_13 * smd_65[k]
                   + f_3 * pc_y[k] * snd_95[k];

        t_159[k] = f_1 * snp0_47[k]
                   - f_2 * snp1_47[k]
                   + f_3 * pc_z[k] * snd_95[k];

        t_160[k] = pb_z[k] * smf0_100[k]
                   - f_4 * pc_z[k] * smf1_100[k];

        t_161[k] = f_12 * smd_66[k]
                   + f_3 * pc_y[k] * snd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, smd_60, smd_99, smd_100, \
                         smd_101, snd_96, snd_99, snd_100, snd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * smd_60[k]
                   + f_3 * pc_z[k] * snd_96[k];

        t_163[k] = f_13 * smd_99[k]
                   + f_3 * pc_x[k] * snd_99[k];

        t_164[k] = f_13 * smd_100[k]
                   + f_3 * pc_x[k] * snd_100[k];

        t_165[k] = f_13 * smd_101[k]
                   + f_3 * pc_x[k] * snd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, smf0_106, smd_63, \
                         smd_65, smd_71, smf1_106, snp0_50, snp1_50, snd_99, \
                         snd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_z[k] * smf0_106[k]
                   - f_4 * pc_z[k] * smf1_106[k];

        t_167[k] = f_5 * smd_63[k]
                   + f_3 * pc_z[k] * snd_99[k];

        t_168[k] = f_12 * smd_71[k]
                   + f_3 * pc_y[k] * snd_101[k];

        t_169[k] = f_5 * smd_65[k]
                   + f_1 * snp0_50[k]
                   - f_2 * snp1_50[k]
                   + f_3 * pc_z[k] * snd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, smd_66, smd_72, \
                         smd_102, smd_105, snp0_51, snp1_51, snd_102, \
                         snd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * smd_102[k]
                   + f_1 * snp0_51[k]
                   - f_2 * snp1_51[k]
                   + f_3 * pc_x[k] * snd_102[k];

        t_171[k] = f_10 * smd_72[k]
                   + f_3 * pc_y[k] * snd_102[k];

        t_172[k] = f_8 * smd_66[k]
                   + f_3 * pc_z[k] * snd_102[k];

        t_173[k] = f_13 * smd_105[k]
                   + f_3 * pc_x[k] * snd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, smd_69, smd_75, \
                         smd_106, smd_107, snp0_52, snp1_52, snd_105, snd_106, \
                         snd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * smd_106[k]
                   + f_3 * pc_x[k] * snd_106[k];

        t_175[k] = f_13 * smd_107[k]
                   + f_3 * pc_x[k] * snd_107[k];

        t_176[k] = f_10 * smd_75[k]
                   + f_1 * snp0_52[k]
                   - f_2 * snp1_52[k]
                   + f_3 * pc_y[k] * snd_105[k];

        t_177[k] = f_8 * smd_69[k]
                   + f_3 * pc_z[k] * snd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, smd_71, smd_77, smd_108, \
                         snp0_53, snp0_54, snp1_53, snp1_54, snd_107, \
                         snd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * smd_77[k]
                   + f_3 * pc_y[k] * snd_107[k];

        t_179[k] = f_8 * smd_71[k]
                   + f_1 * snp0_53[k]
                   - f_2 * snp1_53[k]
                   + f_3 * pc_z[k] * snd_107[k];

        t_180[k] = f_13 * smd_108[k]
                   + f_1 * snp0_54[k]
                   - f_2 * snp1_54[k]
                   + f_3 * pc_x[k] * snd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, smd_72, smd_78, \
                         smd_111, smd_112, snd_108, snd_111, snd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * smd_78[k]
                   + f_3 * pc_y[k] * snd_108[k];

        t_182[k] = f_10 * smd_72[k]
                   + f_3 * pc_z[k] * snd_108[k];

        t_183[k] = f_13 * smd_111[k]
                   + f_3 * pc_x[k] * snd_111[k];

        t_184[k] = f_13 * smd_112[k]
                   + f_3 * pc_x[k] * snd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, smd_75, smd_81, smd_83, \
                         smd_113, snp0_55, snp1_55, snd_111, snd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_13 * smd_113[k]
                   + f_3 * pc_x[k] * snd_113[k];

        t_186[k] = f_8 * smd_81[k]
                   + f_1 * snp0_55[k]
                   - f_2 * snp1_55[k]
                   + f_3 * pc_y[k] * snd_111[k];

        t_187[k] = f_10 * smd_75[k]
                   + f_3 * pc_z[k] * snd_111[k];

        t_188[k] = f_8 * smd_83[k]
                   + f_3 * pc_y[k] * snd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pc_y, pc_z, smf0_140, smd_77, \
                         smd_78, smd_84, smf1_140, snp0_56, snp1_56, snd_113, \
                         snd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * smd_77[k]
                   + f_1 * snp0_56[k]
                   - f_2 * snp1_56[k]
                   + f_3 * pc_z[k] * snd_113[k];

        t_190[k] = pb_y[k] * smf0_140[k]
                   - f_4 * pc_y[k] * smf1_140[k];

        t_191[k] = f_5 * smd_84[k]
                   + f_3 * pc_y[k] * snd_114[k];

        t_192[k] = f_12 * smd_78[k]
                   + f_3 * pc_z[k] * snd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, smd_87, smd_117, smd_118, \
                         smd_119, snp0_58, snp1_58, snd_117, snd_118, \
                         snd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_13 * smd_117[k]
                   + f_3 * pc_x[k] * snd_117[k];

        t_194[k] = f_13 * smd_118[k]
                   + f_3 * pc_x[k] * snd_118[k];

        t_195[k] = f_13 * smd_119[k]
                   + f_3 * pc_x[k] * snd_119[k];

        t_196[k] = f_5 * smd_87[k]
                   + f_1 * snp0_58[k]
                   - f_2 * snp1_58[k]
                   + f_3 * pc_y[k] * snd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, pc_y, pc_z, smf0_149, smd_81, smd_89, \
                         smf1_149, snd_117, snd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * smd_81[k]
                   + f_3 * pc_z[k] * snd_117[k];

        t_198[k] = f_5 * smd_89[k]
                   + f_3 * pc_y[k] * snd_119[k];

        t_199[k] = pb_y[k] * smf0_149[k]
                   - f_4 * pc_y[k] * smf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, smd_84, smd_120, \
                         smd_123, snp0_60, snp1_60, snd_120, snd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_13 * smd_120[k]
                   + f_1 * snp0_60[k]
                   - f_2 * snp1_60[k]
                   + f_3 * pc_x[k] * snd_120[k];

        t_201[k] = f_3 * pc_y[k] * snd_120[k];

        t_202[k] = f_13 * smd_84[k]
                   + f_3 * pc_z[k] * snd_120[k];

        t_203[k] = f_13 * smd_123[k]
                   + f_3 * pc_x[k] * snd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pc_x, pc_y, pc_z, smd_87, smd_124, \
                         smd_125, snp0_61, snp1_61, snd_123, snd_124, \
                         snd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_13 * smd_124[k]
                   + f_3 * pc_x[k] * snd_124[k];

        t_205[k] = f_13 * smd_125[k]
                   + f_3 * pc_x[k] * snd_125[k];

        t_206[k] = f_1 * snp0_61[k]
                   - f_2 * snp1_61[k]
                   + f_3 * pc_y[k] * snd_123[k];

        t_207[k] = f_13 * smd_87[k]
                   + f_3 * pc_z[k] * snd_123[k];

        t_208[k] = f_3 * pc_y[k] * snd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, smd_89, smd_90, \
                         smd_126, snp0_62, snp0_63, snp1_62, snp1_63, snd_125, \
                         snd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_13 * smd_89[k]
                   + f_1 * snp0_62[k]
                   - f_2 * snp1_62[k]
                   + f_3 * pc_z[k] * snd_125[k];

        t_210[k] = f_12 * smd_126[k]
                   + f_1 * snp0_63[k]
                   - f_2 * snp1_63[k]
                   + f_3 * pc_x[k] * snd_126[k];

        t_211[k] = f_11 * smd_90[k]
                   + f_3 * pc_y[k] * snd_126[k];

        t_212[k] = f_3 * pc_z[k] * snd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, smd_93, smd_129, smd_130, \
                         smd_131, snp0_64, snp1_64, snd_129, snd_130, \
                         snd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_12 * smd_129[k]
                   + f_3 * pc_x[k] * snd_129[k];

        t_214[k] = f_12 * smd_130[k]
                   + f_3 * pc_x[k] * snd_130[k];

        t_215[k] = f_12 * smd_131[k]
                   + f_3 * pc_x[k] * snd_131[k];

        t_216[k] = f_11 * smd_93[k]
                   + f_1 * snp0_64[k]
                   - f_2 * snp1_64[k]
                   + f_3 * pc_y[k] * snd_129[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_z, pc_y, pc_z, smf0_150, smd_95, \
                         smf1_150, snp0_65, snp1_65, snd_129, snd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * pc_z[k] * snd_129[k];

        t_218[k] = f_11 * smd_95[k]
                   + f_3 * pc_y[k] * snd_131[k];

        t_219[k] = f_1 * snp0_65[k]
                   - f_2 * snp1_65[k]
                   + f_3 * pc_z[k] * snd_131[k];

        t_220[k] = pb_z[k] * smf0_150[k]
                   - f_4 * pc_z[k] * smf1_150[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_y, pc_z, smd_90, smd_96, \
                         smd_135, smd_136, snd_132, snd_135, snd_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_13 * smd_96[k]
                   + f_3 * pc_y[k] * snd_132[k];

        t_222[k] = f_5 * smd_90[k]
                   + f_3 * pc_z[k] * snd_132[k];

        t_223[k] = f_12 * smd_135[k]
                   + f_3 * pc_x[k] * snd_135[k];

        t_224[k] = f_12 * smd_136[k]
                   + f_3 * pc_x[k] * snd_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_z, pc_x, pc_y, pc_z, smf0_156, smd_93, \
                         smd_101, smd_137, smf1_156, snd_135, snd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_12 * smd_137[k]
                   + f_3 * pc_x[k] * snd_137[k];

        t_226[k] = pb_z[k] * smf0_156[k]
                   - f_4 * pc_z[k] * smf1_156[k];

        t_227[k] = f_5 * smd_93[k]
                   + f_3 * pc_z[k] * snd_135[k];

        t_228[k] = f_13 * smd_101[k]
                   + f_3 * pc_y[k] * snd_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, pc_z, smd_95, smd_102, smd_138, \
                         snp0_68, snp0_69, snp1_68, snp1_69, snd_137, \
                         snd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_5 * smd_95[k]
                   + f_1 * snp0_68[k]
                   - f_2 * snp1_68[k]
                   + f_3 * pc_z[k] * snd_137[k];

        t_230[k] = f_12 * smd_138[k]
                   + f_1 * snp0_69[k]
                   - f_2 * snp1_69[k]
                   + f_3 * pc_x[k] * snd_138[k];

        t_231[k] = f_12 * smd_102[k]
                   + f_3 * pc_y[k] * snd_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, smd_96, smd_141, smd_142, \
                         smd_143, snd_138, snd_141, snd_142, snd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * smd_96[k]
                   + f_3 * pc_z[k] * snd_138[k];

        t_233[k] = f_12 * smd_141[k]
                   + f_3 * pc_x[k] * snd_141[k];

        t_234[k] = f_12 * smd_142[k]
                   + f_3 * pc_x[k] * snd_142[k];

        t_235[k] = f_12 * smd_143[k]
                   + f_3 * pc_x[k] * snd_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_y, pc_z, smd_99, smd_101, smd_105, \
                         smd_107, snp0_70, snp0_71, snp1_70, snp1_71, snd_141, \
                         snd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_12 * smd_105[k]
                   + f_1 * snp0_70[k]
                   - f_2 * snp1_70[k]
                   + f_3 * pc_y[k] * snd_141[k];

        t_237[k] = f_8 * smd_99[k]
                   + f_3 * pc_z[k] * snd_141[k];

        t_238[k] = f_12 * smd_107[k]
                   + f_3 * pc_y[k] * snd_143[k];

        t_239[k] = f_8 * smd_101[k]
                   + f_1 * snp0_71[k]
                   - f_2 * snp1_71[k]
                   + f_3 * pc_z[k] * snd_143[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, pc_y, pc_z, smd_102, smd_108, \
                         smd_144, smd_147, snp0_72, snp1_72, snd_144, \
                         snd_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_12 * smd_144[k]
                   + f_1 * snp0_72[k]
                   - f_2 * snp1_72[k]
                   + f_3 * pc_x[k] * snd_144[k];

        t_241[k] = f_10 * smd_108[k]
                   + f_3 * pc_y[k] * snd_144[k];

        t_242[k] = f_10 * smd_102[k]
                   + f_3 * pc_z[k] * snd_144[k];

        t_243[k] = f_12 * smd_147[k]
                   + f_3 * pc_x[k] * snd_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pc_x, pc_y, pc_z, smd_105, smd_111, \
                         smd_148, smd_149, snp0_73, snp1_73, snd_147, snd_148, \
                         snd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_12 * smd_148[k]
                   + f_3 * pc_x[k] * snd_148[k];

        t_245[k] = f_12 * smd_149[k]
                   + f_3 * pc_x[k] * snd_149[k];

        t_246[k] = f_10 * smd_111[k]
                   + f_1 * snp0_73[k]
                   - f_2 * snp1_73[k]
                   + f_3 * pc_y[k] * snd_147[k];

        t_247[k] = f_10 * smd_105[k]
                   + f_3 * pc_z[k] * snd_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, pc_y, pc_z, smd_107, smd_113, smd_150, \
                         snp0_74, snp0_75, snp1_74, snp1_75, snd_149, \
                         snd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_10 * smd_113[k]
                   + f_3 * pc_y[k] * snd_149[k];

        t_249[k] = f_10 * smd_107[k]
                   + f_1 * snp0_74[k]
                   - f_2 * snp1_74[k]
                   + f_3 * pc_z[k] * snd_149[k];

        t_250[k] = f_12 * smd_150[k]
                   + f_1 * snp0_75[k]
                   - f_2 * snp1_75[k]
                   + f_3 * pc_x[k] * snd_150[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, pc_z, smd_108, smd_114, \
                         smd_153, smd_154, snd_150, snd_153, snd_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * smd_114[k]
                   + f_3 * pc_y[k] * snd_150[k];

        t_252[k] = f_12 * smd_108[k]
                   + f_3 * pc_z[k] * snd_150[k];

        t_253[k] = f_12 * smd_153[k]
                   + f_3 * pc_x[k] * snd_153[k];

        t_254[k] = f_12 * smd_154[k]
                   + f_3 * pc_x[k] * snd_154[k];
    }
}

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_200 = buffer.data(smf0 + 200);
    const auto *smf0_209 = buffer.data(smf0 + 209);
    const auto *smf0_210 = buffer.data(smf0 + 210);
    const auto *smf0_216 = buffer.data(smf0 + 216);
    const auto *smf0_270 = buffer.data(smf0 + 270);
    const auto *smf0_279 = buffer.data(smf0 + 279);
    const auto *smf0_280 = buffer.data(smf0 + 280);
    const auto *smf0_286 = buffer.data(smf0 + 286);

    const auto *smd_111 = buffer.data(smd + 111);
    const auto *smd_113 = buffer.data(smd + 113);
    const auto *smd_114 = buffer.data(smd + 114);
    const auto *smd_117 = buffer.data(smd + 117);
    const auto *smd_119 = buffer.data(smd + 119);
    const auto *smd_120 = buffer.data(smd + 120);
    const auto *smd_123 = buffer.data(smd + 123);
    const auto *smd_125 = buffer.data(smd + 125);
    const auto *smd_126 = buffer.data(smd + 126);
    const auto *smd_129 = buffer.data(smd + 129);
    const auto *smd_131 = buffer.data(smd + 131);
    const auto *smd_132 = buffer.data(smd + 132);
    const auto *smd_135 = buffer.data(smd + 135);
    const auto *smd_137 = buffer.data(smd + 137);
    const auto *smd_138 = buffer.data(smd + 138);
    const auto *smd_141 = buffer.data(smd + 141);
    const auto *smd_143 = buffer.data(smd + 143);
    const auto *smd_144 = buffer.data(smd + 144);
    const auto *smd_147 = buffer.data(smd + 147);
    const auto *smd_149 = buffer.data(smd + 149);
    const auto *smd_150 = buffer.data(smd + 150);
    const auto *smd_153 = buffer.data(smd + 153);
    const auto *smd_155 = buffer.data(smd + 155);
    const auto *smd_156 = buffer.data(smd + 156);
    const auto *smd_159 = buffer.data(smd + 159);
    const auto *smd_160 = buffer.data(smd + 160);
    const auto *smd_161 = buffer.data(smd + 161);
    const auto *smd_162 = buffer.data(smd + 162);
    const auto *smd_165 = buffer.data(smd + 165);
    const auto *smd_166 = buffer.data(smd + 166);
    const auto *smd_167 = buffer.data(smd + 167);
    const auto *smd_168 = buffer.data(smd + 168);
    const auto *smd_171 = buffer.data(smd + 171);
    const auto *smd_172 = buffer.data(smd + 172);
    const auto *smd_173 = buffer.data(smd + 173);
    const auto *smd_174 = buffer.data(smd + 174);
    const auto *smd_177 = buffer.data(smd + 177);
    const auto *smd_178 = buffer.data(smd + 178);
    const auto *smd_179 = buffer.data(smd + 179);
    const auto *smd_180 = buffer.data(smd + 180);
    const auto *smd_183 = buffer.data(smd + 183);
    const auto *smd_184 = buffer.data(smd + 184);
    const auto *smd_185 = buffer.data(smd + 185);
    const auto *smd_186 = buffer.data(smd + 186);
    const auto *smd_189 = buffer.data(smd + 189);
    const auto *smd_190 = buffer.data(smd + 190);
    const auto *smd_191 = buffer.data(smd + 191);
    const auto *smd_192 = buffer.data(smd + 192);
    const auto *smd_195 = buffer.data(smd + 195);
    const auto *smd_196 = buffer.data(smd + 196);
    const auto *smd_197 = buffer.data(smd + 197);
    const auto *smd_198 = buffer.data(smd + 198);
    const auto *smd_201 = buffer.data(smd + 201);
    const auto *smd_202 = buffer.data(smd + 202);
    const auto *smd_203 = buffer.data(smd + 203);
    const auto *smd_207 = buffer.data(smd + 207);
    const auto *smd_208 = buffer.data(smd + 208);
    const auto *smd_209 = buffer.data(smd + 209);
    const auto *smd_210 = buffer.data(smd + 210);
    const auto *smd_213 = buffer.data(smd + 213);
    const auto *smd_214 = buffer.data(smd + 214);
    const auto *smd_215 = buffer.data(smd + 215);
    const auto *smd_216 = buffer.data(smd + 216);
    const auto *smd_219 = buffer.data(smd + 219);
    const auto *smd_220 = buffer.data(smd + 220);
    const auto *smd_221 = buffer.data(smd + 221);
    const auto *smd_225 = buffer.data(smd + 225);
    const auto *smd_226 = buffer.data(smd + 226);
    const auto *smd_227 = buffer.data(smd + 227);

    const auto *smf1_200 = buffer.data(smf1 + 200);
    const auto *smf1_209 = buffer.data(smf1 + 209);
    const auto *smf1_210 = buffer.data(smf1 + 210);
    const auto *smf1_216 = buffer.data(smf1 + 216);
    const auto *smf1_270 = buffer.data(smf1 + 270);
    const auto *smf1_279 = buffer.data(smf1 + 279);
    const auto *smf1_280 = buffer.data(smf1 + 280);
    const auto *smf1_286 = buffer.data(smf1 + 286);

    const auto *snp0_76 = buffer.data(snp0 + 76);
    const auto *snp0_77 = buffer.data(snp0 + 77);
    const auto *snp0_79 = buffer.data(snp0 + 79);
    const auto *snp0_81 = buffer.data(snp0 + 81);
    const auto *snp0_82 = buffer.data(snp0 + 82);
    const auto *snp0_83 = buffer.data(snp0 + 83);
    const auto *snp0_84 = buffer.data(snp0 + 84);
    const auto *snp0_85 = buffer.data(snp0 + 85);
    const auto *snp0_86 = buffer.data(snp0 + 86);
    const auto *snp0_89 = buffer.data(snp0 + 89);
    const auto *snp0_90 = buffer.data(snp0 + 90);
    const auto *snp0_91 = buffer.data(snp0 + 91);
    const auto *snp0_92 = buffer.data(snp0 + 92);
    const auto *snp0_93 = buffer.data(snp0 + 93);
    const auto *snp0_94 = buffer.data(snp0 + 94);
    const auto *snp0_95 = buffer.data(snp0 + 95);
    const auto *snp0_96 = buffer.data(snp0 + 96);
    const auto *snp0_97 = buffer.data(snp0 + 97);
    const auto *snp0_98 = buffer.data(snp0 + 98);
    const auto *snp0_99 = buffer.data(snp0 + 99);
    const auto *snp0_100 = buffer.data(snp0 + 100);
    const auto *snp0_101 = buffer.data(snp0 + 101);
    const auto *snp0_103 = buffer.data(snp0 + 103);
    const auto *snp0_105 = buffer.data(snp0 + 105);
    const auto *snp0_106 = buffer.data(snp0 + 106);
    const auto *snp0_107 = buffer.data(snp0 + 107);
    const auto *snp0_108 = buffer.data(snp0 + 108);
    const auto *snp0_109 = buffer.data(snp0 + 109);
    const auto *snp0_110 = buffer.data(snp0 + 110);
    const auto *snp0_113 = buffer.data(snp0 + 113);

    const auto *snp1_76 = buffer.data(snp1 + 76);
    const auto *snp1_77 = buffer.data(snp1 + 77);
    const auto *snp1_79 = buffer.data(snp1 + 79);
    const auto *snp1_81 = buffer.data(snp1 + 81);
    const auto *snp1_82 = buffer.data(snp1 + 82);
    const auto *snp1_83 = buffer.data(snp1 + 83);
    const auto *snp1_84 = buffer.data(snp1 + 84);
    const auto *snp1_85 = buffer.data(snp1 + 85);
    const auto *snp1_86 = buffer.data(snp1 + 86);
    const auto *snp1_89 = buffer.data(snp1 + 89);
    const auto *snp1_90 = buffer.data(snp1 + 90);
    const auto *snp1_91 = buffer.data(snp1 + 91);
    const auto *snp1_92 = buffer.data(snp1 + 92);
    const auto *snp1_93 = buffer.data(snp1 + 93);
    const auto *snp1_94 = buffer.data(snp1 + 94);
    const auto *snp1_95 = buffer.data(snp1 + 95);
    const auto *snp1_96 = buffer.data(snp1 + 96);
    const auto *snp1_97 = buffer.data(snp1 + 97);
    const auto *snp1_98 = buffer.data(snp1 + 98);
    const auto *snp1_99 = buffer.data(snp1 + 99);
    const auto *snp1_100 = buffer.data(snp1 + 100);
    const auto *snp1_101 = buffer.data(snp1 + 101);
    const auto *snp1_103 = buffer.data(snp1 + 103);
    const auto *snp1_105 = buffer.data(snp1 + 105);
    const auto *snp1_106 = buffer.data(snp1 + 106);
    const auto *snp1_107 = buffer.data(snp1 + 107);
    const auto *snp1_108 = buffer.data(snp1 + 108);
    const auto *snp1_109 = buffer.data(snp1 + 109);
    const auto *snp1_110 = buffer.data(snp1 + 110);
    const auto *snp1_113 = buffer.data(snp1 + 113);

    const auto *snd_153 = buffer.data(snd + 153);
    const auto *snd_155 = buffer.data(snd + 155);
    const auto *snd_156 = buffer.data(snd + 156);
    const auto *snd_159 = buffer.data(snd + 159);
    const auto *snd_160 = buffer.data(snd + 160);
    const auto *snd_161 = buffer.data(snd + 161);
    const auto *snd_162 = buffer.data(snd + 162);
    const auto *snd_165 = buffer.data(snd + 165);
    const auto *snd_166 = buffer.data(snd + 166);
    const auto *snd_167 = buffer.data(snd + 167);
    const auto *snd_168 = buffer.data(snd + 168);
    const auto *snd_171 = buffer.data(snd + 171);
    const auto *snd_172 = buffer.data(snd + 172);
    const auto *snd_173 = buffer.data(snd + 173);
    const auto *snd_174 = buffer.data(snd + 174);
    const auto *snd_177 = buffer.data(snd + 177);
    const auto *snd_178 = buffer.data(snd + 178);
    const auto *snd_179 = buffer.data(snd + 179);
    const auto *snd_180 = buffer.data(snd + 180);
    const auto *snd_183 = buffer.data(snd + 183);
    const auto *snd_184 = buffer.data(snd + 184);
    const auto *snd_185 = buffer.data(snd + 185);
    const auto *snd_186 = buffer.data(snd + 186);
    const auto *snd_189 = buffer.data(snd + 189);
    const auto *snd_190 = buffer.data(snd + 190);
    const auto *snd_191 = buffer.data(snd + 191);
    const auto *snd_192 = buffer.data(snd + 192);
    const auto *snd_195 = buffer.data(snd + 195);
    const auto *snd_196 = buffer.data(snd + 196);
    const auto *snd_197 = buffer.data(snd + 197);
    const auto *snd_198 = buffer.data(snd + 198);
    const auto *snd_201 = buffer.data(snd + 201);
    const auto *snd_202 = buffer.data(snd + 202);
    const auto *snd_203 = buffer.data(snd + 203);
    const auto *snd_204 = buffer.data(snd + 204);
    const auto *snd_207 = buffer.data(snd + 207);
    const auto *snd_208 = buffer.data(snd + 208);
    const auto *snd_209 = buffer.data(snd + 209);
    const auto *snd_210 = buffer.data(snd + 210);
    const auto *snd_213 = buffer.data(snd + 213);
    const auto *snd_214 = buffer.data(snd + 214);
    const auto *snd_215 = buffer.data(snd + 215);
    const auto *snd_216 = buffer.data(snd + 216);
    const auto *snd_219 = buffer.data(snd + 219);
    const auto *snd_220 = buffer.data(snd + 220);
    const auto *snd_221 = buffer.data(snd + 221);
    const auto *snd_222 = buffer.data(snd + 222);
    const auto *snd_225 = buffer.data(snd + 225);
    const auto *snd_226 = buffer.data(snd + 226);
    const auto *snd_227 = buffer.data(snd + 227);

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, smd_111, smd_117, \
                         smd_119, smd_155, snp0_76, snp1_76, snd_153, \
                         snd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_12 * smd_155[k]
                   + f_3 * pc_x[k] * snd_155[k];

        t_256[k] = f_8 * smd_117[k]
                   + f_1 * snp0_76[k]
                   - f_2 * snp1_76[k]
                   + f_3 * pc_y[k] * snd_153[k];

        t_257[k] = f_12 * smd_111[k]
                   + f_3 * pc_z[k] * snd_153[k];

        t_258[k] = f_8 * smd_119[k]
                   + f_3 * pc_y[k] * snd_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_y, pc_y, pc_z, smf0_200, smd_113, \
                         smd_114, smd_120, smf1_200, snp0_77, snp1_77, snd_155, \
                         snd_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_12 * smd_113[k]
                   + f_1 * snp0_77[k]
                   - f_2 * snp1_77[k]
                   + f_3 * pc_z[k] * snd_155[k];

        t_260[k] = pb_y[k] * smf0_200[k]
                   - f_4 * pc_y[k] * smf1_200[k];

        t_261[k] = f_5 * smd_120[k]
                   + f_3 * pc_y[k] * snd_156[k];

        t_262[k] = f_13 * smd_114[k]
                   + f_3 * pc_z[k] * snd_156[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_x, pc_y, smd_123, smd_159, smd_160, \
                         smd_161, snp0_79, snp1_79, snd_159, snd_160, \
                         snd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_12 * smd_159[k]
                   + f_3 * pc_x[k] * snd_159[k];

        t_264[k] = f_12 * smd_160[k]
                   + f_3 * pc_x[k] * snd_160[k];

        t_265[k] = f_12 * smd_161[k]
                   + f_3 * pc_x[k] * snd_161[k];

        t_266[k] = f_5 * smd_123[k]
                   + f_1 * snp0_79[k]
                   - f_2 * snp1_79[k]
                   + f_3 * pc_y[k] * snd_159[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_y, pc_y, pc_z, smf0_209, smd_117, smd_125, \
                         smf1_209, snd_159, snd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_13 * smd_117[k]
                   + f_3 * pc_z[k] * snd_159[k];

        t_268[k] = f_5 * smd_125[k]
                   + f_3 * pc_y[k] * snd_161[k];

        t_269[k] = pb_y[k] * smf0_209[k]
                   - f_4 * pc_y[k] * smf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, pc_z, smd_120, smd_162, \
                         smd_165, snp0_81, snp1_81, snd_162, snd_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_12 * smd_162[k]
                   + f_1 * snp0_81[k]
                   - f_2 * snp1_81[k]
                   + f_3 * pc_x[k] * snd_162[k];

        t_271[k] = f_3 * pc_y[k] * snd_162[k];

        t_272[k] = f_11 * smd_120[k]
                   + f_3 * pc_z[k] * snd_162[k];

        t_273[k] = f_12 * smd_165[k]
                   + f_3 * pc_x[k] * snd_165[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_y, pc_z, smd_123, \
                         smd_166, smd_167, snp0_82, snp1_82, snd_165, snd_166, \
                         snd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_12 * smd_166[k]
                   + f_3 * pc_x[k] * snd_166[k];

        t_275[k] = f_12 * smd_167[k]
                   + f_3 * pc_x[k] * snd_167[k];

        t_276[k] = f_1 * snp0_82[k]
                   - f_2 * snp1_82[k]
                   + f_3 * pc_y[k] * snd_165[k];

        t_277[k] = f_11 * smd_123[k]
                   + f_3 * pc_z[k] * snd_165[k];

        t_278[k] = f_3 * pc_y[k] * snd_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pc_x, pc_y, pc_z, smd_125, smd_126, \
                         smd_168, snp0_83, snp0_84, snp1_83, snp1_84, snd_167, \
                         snd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_11 * smd_125[k]
                   + f_1 * snp0_83[k]
                   - f_2 * snp1_83[k]
                   + f_3 * pc_z[k] * snd_167[k];

        t_280[k] = f_10 * smd_168[k]
                   + f_1 * snp0_84[k]
                   - f_2 * snp1_84[k]
                   + f_3 * pc_x[k] * snd_168[k];

        t_281[k] = f_9 * smd_126[k]
                   + f_3 * pc_y[k] * snd_168[k];

        t_282[k] = f_3 * pc_z[k] * snd_168[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pc_x, pc_y, smd_129, smd_171, smd_172, \
                         smd_173, snp0_85, snp1_85, snd_171, snd_172, \
                         snd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_10 * smd_171[k]
                   + f_3 * pc_x[k] * snd_171[k];

        t_284[k] = f_10 * smd_172[k]
                   + f_3 * pc_x[k] * snd_172[k];

        t_285[k] = f_10 * smd_173[k]
                   + f_3 * pc_x[k] * snd_173[k];

        t_286[k] = f_9 * smd_129[k]
                   + f_1 * snp0_85[k]
                   - f_2 * snp1_85[k]
                   + f_3 * pc_y[k] * snd_171[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_z, pc_y, pc_z, smf0_210, smd_131, \
                         smf1_210, snp0_86, snp1_86, snd_171, snd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * snd_171[k];

        t_288[k] = f_9 * smd_131[k]
                   + f_3 * pc_y[k] * snd_173[k];

        t_289[k] = f_1 * snp0_86[k]
                   - f_2 * snp1_86[k]
                   + f_3 * pc_z[k] * snd_173[k];

        t_290[k] = pb_z[k] * smf0_210[k]
                   - f_4 * pc_z[k] * smf1_210[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, pc_y, pc_z, smd_126, smd_132, \
                         smd_177, smd_178, snd_174, snd_177, snd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_11 * smd_132[k]
                   + f_3 * pc_y[k] * snd_174[k];

        t_292[k] = f_5 * smd_126[k]
                   + f_3 * pc_z[k] * snd_174[k];

        t_293[k] = f_10 * smd_177[k]
                   + f_3 * pc_x[k] * snd_177[k];

        t_294[k] = f_10 * smd_178[k]
                   + f_3 * pc_x[k] * snd_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pb_z, pc_x, pc_y, pc_z, smf0_216, \
                         smd_129, smd_137, smd_179, smf1_216, snd_177, \
                         snd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_10 * smd_179[k]
                   + f_3 * pc_x[k] * snd_179[k];

        t_296[k] = pb_z[k] * smf0_216[k]
                   - f_4 * pc_z[k] * smf1_216[k];

        t_297[k] = f_5 * smd_129[k]
                   + f_3 * pc_z[k] * snd_177[k];

        t_298[k] = f_11 * smd_137[k]
                   + f_3 * pc_y[k] * snd_179[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_y, pc_z, smd_131, smd_138, smd_180, \
                         snp0_89, snp0_90, snp1_89, snp1_90, snd_179, \
                         snd_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_5 * smd_131[k]
                   + f_1 * snp0_89[k]
                   - f_2 * snp1_89[k]
                   + f_3 * pc_z[k] * snd_179[k];

        t_300[k] = f_10 * smd_180[k]
                   + f_1 * snp0_90[k]
                   - f_2 * snp1_90[k]
                   + f_3 * pc_x[k] * snd_180[k];

        t_301[k] = f_13 * smd_138[k]
                   + f_3 * pc_y[k] * snd_180[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_z, smd_132, smd_183, smd_184, \
                         smd_185, snd_180, snd_183, snd_184, snd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * smd_132[k]
                   + f_3 * pc_z[k] * snd_180[k];

        t_303[k] = f_10 * smd_183[k]
                   + f_3 * pc_x[k] * snd_183[k];

        t_304[k] = f_10 * smd_184[k]
                   + f_3 * pc_x[k] * snd_184[k];

        t_305[k] = f_10 * smd_185[k]
                   + f_3 * pc_x[k] * snd_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_y, pc_z, smd_135, smd_137, smd_141, \
                         smd_143, snp0_91, snp0_92, snp1_91, snp1_92, snd_183, \
                         snd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_13 * smd_141[k]
                   + f_1 * snp0_91[k]
                   - f_2 * snp1_91[k]
                   + f_3 * pc_y[k] * snd_183[k];

        t_307[k] = f_8 * smd_135[k]
                   + f_3 * pc_z[k] * snd_183[k];

        t_308[k] = f_13 * smd_143[k]
                   + f_3 * pc_y[k] * snd_185[k];

        t_309[k] = f_8 * smd_137[k]
                   + f_1 * snp0_92[k]
                   - f_2 * snp1_92[k]
                   + f_3 * pc_z[k] * snd_185[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_x, pc_y, pc_z, smd_138, smd_144, \
                         smd_186, smd_189, snp0_93, snp1_93, snd_186, \
                         snd_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_10 * smd_186[k]
                   + f_1 * snp0_93[k]
                   - f_2 * snp1_93[k]
                   + f_3 * pc_x[k] * snd_186[k];

        t_311[k] = f_12 * smd_144[k]
                   + f_3 * pc_y[k] * snd_186[k];

        t_312[k] = f_10 * smd_138[k]
                   + f_3 * pc_z[k] * snd_186[k];

        t_313[k] = f_10 * smd_189[k]
                   + f_3 * pc_x[k] * snd_189[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, smd_141, smd_147, \
                         smd_190, smd_191, snp0_94, snp1_94, snd_189, snd_190, \
                         snd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_10 * smd_190[k]
                   + f_3 * pc_x[k] * snd_190[k];

        t_315[k] = f_10 * smd_191[k]
                   + f_3 * pc_x[k] * snd_191[k];

        t_316[k] = f_12 * smd_147[k]
                   + f_1 * snp0_94[k]
                   - f_2 * snp1_94[k]
                   + f_3 * pc_y[k] * snd_189[k];

        t_317[k] = f_10 * smd_141[k]
                   + f_3 * pc_z[k] * snd_189[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_x, pc_y, pc_z, smd_143, smd_149, smd_192, \
                         snp0_95, snp0_96, snp1_95, snp1_96, snd_191, \
                         snd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_12 * smd_149[k]
                   + f_3 * pc_y[k] * snd_191[k];

        t_319[k] = f_10 * smd_143[k]
                   + f_1 * snp0_95[k]
                   - f_2 * snp1_95[k]
                   + f_3 * pc_z[k] * snd_191[k];

        t_320[k] = f_10 * smd_192[k]
                   + f_1 * snp0_96[k]
                   - f_2 * snp1_96[k]
                   + f_3 * pc_x[k] * snd_192[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pc_y, pc_z, smd_144, smd_150, \
                         smd_195, smd_196, snd_192, snd_195, snd_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_10 * smd_150[k]
                   + f_3 * pc_y[k] * snd_192[k];

        t_322[k] = f_12 * smd_144[k]
                   + f_3 * pc_z[k] * snd_192[k];

        t_323[k] = f_10 * smd_195[k]
                   + f_3 * pc_x[k] * snd_195[k];

        t_324[k] = f_10 * smd_196[k]
                   + f_3 * pc_x[k] * snd_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, smd_147, smd_153, \
                         smd_155, smd_197, snp0_97, snp1_97, snd_195, \
                         snd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_10 * smd_197[k]
                   + f_3 * pc_x[k] * snd_197[k];

        t_326[k] = f_10 * smd_153[k]
                   + f_1 * snp0_97[k]
                   - f_2 * snp1_97[k]
                   + f_3 * pc_y[k] * snd_195[k];

        t_327[k] = f_12 * smd_147[k]
                   + f_3 * pc_z[k] * snd_195[k];

        t_328[k] = f_10 * smd_155[k]
                   + f_3 * pc_y[k] * snd_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pc_x, pc_y, pc_z, smd_149, smd_156, smd_198, \
                         snp0_98, snp0_99, snp1_98, snp1_99, snd_197, \
                         snd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_12 * smd_149[k]
                   + f_1 * snp0_98[k]
                   - f_2 * snp1_98[k]
                   + f_3 * pc_z[k] * snd_197[k];

        t_330[k] = f_10 * smd_198[k]
                   + f_1 * snp0_99[k]
                   - f_2 * snp1_99[k]
                   + f_3 * pc_x[k] * snd_198[k];

        t_331[k] = f_8 * smd_156[k]
                   + f_3 * pc_y[k] * snd_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pc_x, pc_z, smd_150, smd_201, smd_202, \
                         smd_203, snd_198, snd_201, snd_202, snd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_13 * smd_150[k]
                   + f_3 * pc_z[k] * snd_198[k];

        t_333[k] = f_10 * smd_201[k]
                   + f_3 * pc_x[k] * snd_201[k];

        t_334[k] = f_10 * smd_202[k]
                   + f_3 * pc_x[k] * snd_202[k];

        t_335[k] = f_10 * smd_203[k]
                   + f_3 * pc_x[k] * snd_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_y, pc_z, smd_153, smd_155, smd_159, \
                         smd_161, snp0_100, snp0_101, snp1_100, snp1_101, snd_201, \
                         snd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_8 * smd_159[k]
                   + f_1 * snp0_100[k]
                   - f_2 * snp1_100[k]
                   + f_3 * pc_y[k] * snd_201[k];

        t_337[k] = f_13 * smd_153[k]
                   + f_3 * pc_z[k] * snd_201[k];

        t_338[k] = f_8 * smd_161[k]
                   + f_3 * pc_y[k] * snd_203[k];

        t_339[k] = f_13 * smd_155[k]
                   + f_1 * snp0_101[k]
                   - f_2 * snp1_101[k]
                   + f_3 * pc_z[k] * snd_203[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_y, pc_x, pc_y, pc_z, smf0_270, \
                         smd_156, smd_162, smd_207, smf1_270, snd_204, \
                         snd_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pb_y[k] * smf0_270[k]
                   - f_4 * pc_y[k] * smf1_270[k];

        t_341[k] = f_5 * smd_162[k]
                   + f_3 * pc_y[k] * snd_204[k];

        t_342[k] = f_11 * smd_156[k]
                   + f_3 * pc_z[k] * snd_204[k];

        t_343[k] = f_10 * smd_207[k]
                   + f_3 * pc_x[k] * snd_207[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, pc_z, smd_159, smd_165, \
                         smd_208, smd_209, snp0_103, snp1_103, snd_207, snd_208, \
                         snd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_10 * smd_208[k]
                   + f_3 * pc_x[k] * snd_208[k];

        t_345[k] = f_10 * smd_209[k]
                   + f_3 * pc_x[k] * snd_209[k];

        t_346[k] = f_5 * smd_165[k]
                   + f_1 * snp0_103[k]
                   - f_2 * snp1_103[k]
                   + f_3 * pc_y[k] * snd_207[k];

        t_347[k] = f_11 * smd_159[k]
                   + f_3 * pc_z[k] * snd_207[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_y, pc_x, pc_y, smf0_279, smd_167, \
                         smd_210, smf1_279, snp0_105, snp1_105, snd_209, \
                         snd_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_5 * smd_167[k]
                   + f_3 * pc_y[k] * snd_209[k];

        t_349[k] = pb_y[k] * smf0_279[k]
                   - f_4 * pc_y[k] * smf1_279[k];

        t_350[k] = f_10 * smd_210[k]
                   + f_1 * snp0_105[k]
                   - f_2 * snp1_105[k]
                   + f_3 * pc_x[k] * snd_210[k];

        t_351[k] = f_3 * pc_y[k] * snd_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_x, pc_z, smd_162, smd_213, smd_214, \
                         smd_215, snd_210, snd_213, snd_214, snd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_9 * smd_162[k]
                   + f_3 * pc_z[k] * snd_210[k];

        t_353[k] = f_10 * smd_213[k]
                   + f_3 * pc_x[k] * snd_213[k];

        t_354[k] = f_10 * smd_214[k]
                   + f_3 * pc_x[k] * snd_214[k];

        t_355[k] = f_10 * smd_215[k]
                   + f_3 * pc_x[k] * snd_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, smd_165, smd_167, snp0_106, \
                         snp0_107, snp1_106, snp1_107, snd_213, \
                         snd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_1 * snp0_106[k]
                   - f_2 * snp1_106[k]
                   + f_3 * pc_y[k] * snd_213[k];

        t_357[k] = f_9 * smd_165[k]
                   + f_3 * pc_z[k] * snd_213[k];

        t_358[k] = f_3 * pc_y[k] * snd_215[k];

        t_359[k] = f_9 * smd_167[k]
                   + f_1 * snp0_107[k]
                   - f_2 * snp1_107[k]
                   + f_3 * pc_z[k] * snd_215[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, smd_168, smd_216, \
                         smd_219, snp0_108, snp1_108, snd_216, \
                         snd_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_8 * smd_216[k]
                   + f_1 * snp0_108[k]
                   - f_2 * snp1_108[k]
                   + f_3 * pc_x[k] * snd_216[k];

        t_361[k] = f_7 * smd_168[k]
                   + f_3 * pc_y[k] * snd_216[k];

        t_362[k] = f_3 * pc_z[k] * snd_216[k];

        t_363[k] = f_8 * smd_219[k]
                   + f_3 * pc_x[k] * snd_219[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, smd_171, smd_220, \
                         smd_221, snp0_109, snp1_109, snd_219, snd_220, \
                         snd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_8 * smd_220[k]
                   + f_3 * pc_x[k] * snd_220[k];

        t_365[k] = f_8 * smd_221[k]
                   + f_3 * pc_x[k] * snd_221[k];

        t_366[k] = f_7 * smd_171[k]
                   + f_1 * snp0_109[k]
                   - f_2 * snp1_109[k]
                   + f_3 * pc_y[k] * snd_219[k];

        t_367[k] = f_3 * pc_z[k] * snd_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pb_z, pc_y, pc_z, smf0_280, smd_173, \
                         smd_174, smf1_280, snp0_110, snp1_110, snd_221, \
                         snd_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_7 * smd_173[k]
                   + f_3 * pc_y[k] * snd_221[k];

        t_369[k] = f_1 * snp0_110[k]
                   - f_2 * snp1_110[k]
                   + f_3 * pc_z[k] * snd_221[k];

        t_370[k] = pb_z[k] * smf0_280[k]
                   - f_4 * pc_z[k] * smf1_280[k];

        t_371[k] = f_9 * smd_174[k]
                   + f_3 * pc_y[k] * snd_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, pc_z, smd_168, smd_225, smd_226, \
                         smd_227, snd_222, snd_225, snd_226, snd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_5 * smd_168[k]
                   + f_3 * pc_z[k] * snd_222[k];

        t_373[k] = f_8 * smd_225[k]
                   + f_3 * pc_x[k] * snd_225[k];

        t_374[k] = f_8 * smd_226[k]
                   + f_3 * pc_x[k] * snd_226[k];

        t_375[k] = f_8 * smd_227[k]
                   + f_3 * pc_x[k] * snd_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pb_z, pc_y, pc_z, smf0_286, smd_171, \
                         smd_173, smd_179, smf1_286, snp0_113, snp1_113, snd_225, \
                         snd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pb_z[k] * smf0_286[k]
                   - f_4 * pc_z[k] * smf1_286[k];

        t_377[k] = f_5 * smd_171[k]
                   + f_3 * pc_z[k] * snd_225[k];

        t_378[k] = f_9 * smd_179[k]
                   + f_3 * pc_y[k] * snd_227[k];

        t_379[k] = f_5 * smd_173[k]
                   + f_1 * snp0_113[k]
                   - f_2 * snp1_113[k]
                   + f_3 * pc_z[k] * snd_227[k];
    }
}

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

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
    auto *t_396 = buffer.data(target + 396);
    auto *t_397 = buffer.data(target + 397);
    auto *t_398 = buffer.data(target + 398);
    auto *t_399 = buffer.data(target + 399);
    auto *t_400 = buffer.data(target + 400);
    auto *t_401 = buffer.data(target + 401);
    auto *t_402 = buffer.data(target + 402);
    auto *t_403 = buffer.data(target + 403);
    auto *t_404 = buffer.data(target + 404);
    auto *t_405 = buffer.data(target + 405);
    auto *t_406 = buffer.data(target + 406);
    auto *t_407 = buffer.data(target + 407);
    auto *t_408 = buffer.data(target + 408);
    auto *t_409 = buffer.data(target + 409);
    auto *t_410 = buffer.data(target + 410);
    auto *t_411 = buffer.data(target + 411);
    auto *t_412 = buffer.data(target + 412);
    auto *t_413 = buffer.data(target + 413);
    auto *t_414 = buffer.data(target + 414);
    auto *t_415 = buffer.data(target + 415);
    auto *t_416 = buffer.data(target + 416);
    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);
    auto *t_420 = buffer.data(target + 420);
    auto *t_421 = buffer.data(target + 421);
    auto *t_422 = buffer.data(target + 422);
    auto *t_423 = buffer.data(target + 423);
    auto *t_424 = buffer.data(target + 424);
    auto *t_425 = buffer.data(target + 425);
    auto *t_426 = buffer.data(target + 426);
    auto *t_427 = buffer.data(target + 427);
    auto *t_428 = buffer.data(target + 428);
    auto *t_429 = buffer.data(target + 429);
    auto *t_430 = buffer.data(target + 430);
    auto *t_431 = buffer.data(target + 431);
    auto *t_432 = buffer.data(target + 432);
    auto *t_433 = buffer.data(target + 433);
    auto *t_434 = buffer.data(target + 434);
    auto *t_435 = buffer.data(target + 435);
    auto *t_436 = buffer.data(target + 436);
    auto *t_437 = buffer.data(target + 437);
    auto *t_438 = buffer.data(target + 438);
    auto *t_439 = buffer.data(target + 439);
    auto *t_440 = buffer.data(target + 440);
    auto *t_441 = buffer.data(target + 441);
    auto *t_442 = buffer.data(target + 442);
    auto *t_443 = buffer.data(target + 443);
    auto *t_444 = buffer.data(target + 444);
    auto *t_445 = buffer.data(target + 445);
    auto *t_446 = buffer.data(target + 446);
    auto *t_447 = buffer.data(target + 447);
    auto *t_448 = buffer.data(target + 448);
    auto *t_449 = buffer.data(target + 449);
    auto *t_450 = buffer.data(target + 450);
    auto *t_451 = buffer.data(target + 451);
    auto *t_452 = buffer.data(target + 452);
    auto *t_453 = buffer.data(target + 453);
    auto *t_454 = buffer.data(target + 454);
    auto *t_455 = buffer.data(target + 455);
    auto *t_456 = buffer.data(target + 456);
    auto *t_457 = buffer.data(target + 457);
    auto *t_458 = buffer.data(target + 458);
    auto *t_459 = buffer.data(target + 459);
    auto *t_460 = buffer.data(target + 460);
    auto *t_461 = buffer.data(target + 461);
    auto *t_462 = buffer.data(target + 462);
    auto *t_463 = buffer.data(target + 463);
    auto *t_464 = buffer.data(target + 464);
    auto *t_465 = buffer.data(target + 465);
    auto *t_466 = buffer.data(target + 466);
    auto *t_467 = buffer.data(target + 467);
    auto *t_468 = buffer.data(target + 468);
    auto *t_469 = buffer.data(target + 469);
    auto *t_470 = buffer.data(target + 470);
    auto *t_471 = buffer.data(target + 471);
    auto *t_472 = buffer.data(target + 472);
    auto *t_473 = buffer.data(target + 473);
    auto *t_474 = buffer.data(target + 474);
    auto *t_475 = buffer.data(target + 475);
    auto *t_476 = buffer.data(target + 476);
    auto *t_477 = buffer.data(target + 477);
    auto *t_478 = buffer.data(target + 478);
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);
    auto *t_484 = buffer.data(target + 484);
    auto *t_485 = buffer.data(target + 485);
    auto *t_486 = buffer.data(target + 486);
    auto *t_487 = buffer.data(target + 487);
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);
    auto *t_500 = buffer.data(target + 500);
    auto *t_501 = buffer.data(target + 501);
    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_350 = buffer.data(smf0 + 350);
    const auto *smf0_359 = buffer.data(smf0 + 359);
    const auto *smf0_360 = buffer.data(smf0 + 360);
    const auto *smf0_450 = buffer.data(smf0 + 450);
    const auto *smf0_456 = buffer.data(smf0 + 456);
    const auto *smf0_459 = buffer.data(smf0 + 459);
    const auto *smf0_466 = buffer.data(smf0 + 466);
    const auto *smf0_469 = buffer.data(smf0 + 469);
    const auto *smf0_470 = buffer.data(smf0 + 470);
    const auto *smf0_476 = buffer.data(smf0 + 476);
    const auto *smf0_479 = buffer.data(smf0 + 479);
    const auto *smf0_480 = buffer.data(smf0 + 480);
    const auto *smf0_486 = buffer.data(smf0 + 486);
    const auto *smf0_489 = buffer.data(smf0 + 489);
    const auto *smf0_490 = buffer.data(smf0 + 490);
    const auto *smf0_496 = buffer.data(smf0 + 496);
    const auto *smf0_499 = buffer.data(smf0 + 499);
    const auto *smf0_500 = buffer.data(smf0 + 500);

    const auto *smd_174 = buffer.data(smd + 174);
    const auto *smd_177 = buffer.data(smd + 177);
    const auto *smd_179 = buffer.data(smd + 179);
    const auto *smd_180 = buffer.data(smd + 180);
    const auto *smd_183 = buffer.data(smd + 183);
    const auto *smd_185 = buffer.data(smd + 185);
    const auto *smd_186 = buffer.data(smd + 186);
    const auto *smd_189 = buffer.data(smd + 189);
    const auto *smd_191 = buffer.data(smd + 191);
    const auto *smd_192 = buffer.data(smd + 192);
    const auto *smd_195 = buffer.data(smd + 195);
    const auto *smd_197 = buffer.data(smd + 197);
    const auto *smd_198 = buffer.data(smd + 198);
    const auto *smd_201 = buffer.data(smd + 201);
    const auto *smd_203 = buffer.data(smd + 203);
    const auto *smd_204 = buffer.data(smd + 204);
    const auto *smd_207 = buffer.data(smd + 207);
    const auto *smd_209 = buffer.data(smd + 209);
    const auto *smd_210 = buffer.data(smd + 210);
    const auto *smd_213 = buffer.data(smd + 213);
    const auto *smd_215 = buffer.data(smd + 215);
    const auto *smd_216 = buffer.data(smd + 216);
    const auto *smd_219 = buffer.data(smd + 219);
    const auto *smd_221 = buffer.data(smd + 221);
    const auto *smd_222 = buffer.data(smd + 222);
    const auto *smd_225 = buffer.data(smd + 225);
    const auto *smd_227 = buffer.data(smd + 227);
    const auto *smd_228 = buffer.data(smd + 228);
    const auto *smd_231 = buffer.data(smd + 231);
    const auto *smd_232 = buffer.data(smd + 232);
    const auto *smd_233 = buffer.data(smd + 233);
    const auto *smd_234 = buffer.data(smd + 234);
    const auto *smd_237 = buffer.data(smd + 237);
    const auto *smd_238 = buffer.data(smd + 238);
    const auto *smd_239 = buffer.data(smd + 239);
    const auto *smd_240 = buffer.data(smd + 240);
    const auto *smd_243 = buffer.data(smd + 243);
    const auto *smd_244 = buffer.data(smd + 244);
    const auto *smd_245 = buffer.data(smd + 245);
    const auto *smd_246 = buffer.data(smd + 246);
    const auto *smd_249 = buffer.data(smd + 249);
    const auto *smd_250 = buffer.data(smd + 250);
    const auto *smd_251 = buffer.data(smd + 251);
    const auto *smd_252 = buffer.data(smd + 252);
    const auto *smd_255 = buffer.data(smd + 255);
    const auto *smd_256 = buffer.data(smd + 256);
    const auto *smd_257 = buffer.data(smd + 257);
    const auto *smd_261 = buffer.data(smd + 261);
    const auto *smd_262 = buffer.data(smd + 262);
    const auto *smd_263 = buffer.data(smd + 263);
    const auto *smd_264 = buffer.data(smd + 264);
    const auto *smd_267 = buffer.data(smd + 267);
    const auto *smd_268 = buffer.data(smd + 268);
    const auto *smd_269 = buffer.data(smd + 269);
    const auto *smd_270 = buffer.data(smd + 270);
    const auto *smd_273 = buffer.data(smd + 273);
    const auto *smd_274 = buffer.data(smd + 274);
    const auto *smd_275 = buffer.data(smd + 275);
    const auto *smd_279 = buffer.data(smd + 279);
    const auto *smd_280 = buffer.data(smd + 280);
    const auto *smd_281 = buffer.data(smd + 281);
    const auto *smd_282 = buffer.data(smd + 282);
    const auto *smd_285 = buffer.data(smd + 285);
    const auto *smd_286 = buffer.data(smd + 286);
    const auto *smd_287 = buffer.data(smd + 287);
    const auto *smd_288 = buffer.data(smd + 288);
    const auto *smd_291 = buffer.data(smd + 291);
    const auto *smd_292 = buffer.data(smd + 292);
    const auto *smd_293 = buffer.data(smd + 293);
    const auto *smd_294 = buffer.data(smd + 294);
    const auto *smd_297 = buffer.data(smd + 297);
    const auto *smd_298 = buffer.data(smd + 298);
    const auto *smd_299 = buffer.data(smd + 299);
    const auto *smd_300 = buffer.data(smd + 300);
    const auto *smd_303 = buffer.data(smd + 303);

    const auto *smf1_350 = buffer.data(smf1 + 350);
    const auto *smf1_359 = buffer.data(smf1 + 359);
    const auto *smf1_360 = buffer.data(smf1 + 360);
    const auto *smf1_450 = buffer.data(smf1 + 450);
    const auto *smf1_456 = buffer.data(smf1 + 456);
    const auto *smf1_459 = buffer.data(smf1 + 459);
    const auto *smf1_466 = buffer.data(smf1 + 466);
    const auto *smf1_469 = buffer.data(smf1 + 469);
    const auto *smf1_470 = buffer.data(smf1 + 470);
    const auto *smf1_476 = buffer.data(smf1 + 476);
    const auto *smf1_479 = buffer.data(smf1 + 479);
    const auto *smf1_480 = buffer.data(smf1 + 480);
    const auto *smf1_486 = buffer.data(smf1 + 486);
    const auto *smf1_489 = buffer.data(smf1 + 489);
    const auto *smf1_490 = buffer.data(smf1 + 490);
    const auto *smf1_496 = buffer.data(smf1 + 496);
    const auto *smf1_499 = buffer.data(smf1 + 499);
    const auto *smf1_500 = buffer.data(smf1 + 500);

    const auto *snp0_114 = buffer.data(snp0 + 114);
    const auto *snp0_115 = buffer.data(snp0 + 115);
    const auto *snp0_116 = buffer.data(snp0 + 116);
    const auto *snp0_117 = buffer.data(snp0 + 117);
    const auto *snp0_118 = buffer.data(snp0 + 118);
    const auto *snp0_119 = buffer.data(snp0 + 119);
    const auto *snp0_120 = buffer.data(snp0 + 120);
    const auto *snp0_121 = buffer.data(snp0 + 121);
    const auto *snp0_122 = buffer.data(snp0 + 122);
    const auto *snp0_123 = buffer.data(snp0 + 123);
    const auto *snp0_124 = buffer.data(snp0 + 124);
    const auto *snp0_125 = buffer.data(snp0 + 125);
    const auto *snp0_126 = buffer.data(snp0 + 126);
    const auto *snp0_127 = buffer.data(snp0 + 127);
    const auto *snp0_128 = buffer.data(snp0 + 128);
    const auto *snp0_130 = buffer.data(snp0 + 130);
    const auto *snp0_132 = buffer.data(snp0 + 132);
    const auto *snp0_133 = buffer.data(snp0 + 133);
    const auto *snp0_134 = buffer.data(snp0 + 134);

    const auto *snp1_114 = buffer.data(snp1 + 114);
    const auto *snp1_115 = buffer.data(snp1 + 115);
    const auto *snp1_116 = buffer.data(snp1 + 116);
    const auto *snp1_117 = buffer.data(snp1 + 117);
    const auto *snp1_118 = buffer.data(snp1 + 118);
    const auto *snp1_119 = buffer.data(snp1 + 119);
    const auto *snp1_120 = buffer.data(snp1 + 120);
    const auto *snp1_121 = buffer.data(snp1 + 121);
    const auto *snp1_122 = buffer.data(snp1 + 122);
    const auto *snp1_123 = buffer.data(snp1 + 123);
    const auto *snp1_124 = buffer.data(snp1 + 124);
    const auto *snp1_125 = buffer.data(snp1 + 125);
    const auto *snp1_126 = buffer.data(snp1 + 126);
    const auto *snp1_127 = buffer.data(snp1 + 127);
    const auto *snp1_128 = buffer.data(snp1 + 128);
    const auto *snp1_130 = buffer.data(snp1 + 130);
    const auto *snp1_132 = buffer.data(snp1 + 132);
    const auto *snp1_133 = buffer.data(snp1 + 133);
    const auto *snp1_134 = buffer.data(snp1 + 134);

    const auto *snd_228 = buffer.data(snd + 228);
    const auto *snd_231 = buffer.data(snd + 231);
    const auto *snd_232 = buffer.data(snd + 232);
    const auto *snd_233 = buffer.data(snd + 233);
    const auto *snd_234 = buffer.data(snd + 234);
    const auto *snd_237 = buffer.data(snd + 237);
    const auto *snd_238 = buffer.data(snd + 238);
    const auto *snd_239 = buffer.data(snd + 239);
    const auto *snd_240 = buffer.data(snd + 240);
    const auto *snd_243 = buffer.data(snd + 243);
    const auto *snd_244 = buffer.data(snd + 244);
    const auto *snd_245 = buffer.data(snd + 245);
    const auto *snd_246 = buffer.data(snd + 246);
    const auto *snd_249 = buffer.data(snd + 249);
    const auto *snd_250 = buffer.data(snd + 250);
    const auto *snd_251 = buffer.data(snd + 251);
    const auto *snd_252 = buffer.data(snd + 252);
    const auto *snd_255 = buffer.data(snd + 255);
    const auto *snd_256 = buffer.data(snd + 256);
    const auto *snd_257 = buffer.data(snd + 257);
    const auto *snd_258 = buffer.data(snd + 258);
    const auto *snd_261 = buffer.data(snd + 261);
    const auto *snd_262 = buffer.data(snd + 262);
    const auto *snd_263 = buffer.data(snd + 263);
    const auto *snd_264 = buffer.data(snd + 264);
    const auto *snd_267 = buffer.data(snd + 267);
    const auto *snd_268 = buffer.data(snd + 268);
    const auto *snd_269 = buffer.data(snd + 269);
    const auto *snd_270 = buffer.data(snd + 270);
    const auto *snd_273 = buffer.data(snd + 273);
    const auto *snd_274 = buffer.data(snd + 274);
    const auto *snd_275 = buffer.data(snd + 275);
    const auto *snd_276 = buffer.data(snd + 276);
    const auto *snd_279 = buffer.data(snd + 279);
    const auto *snd_280 = buffer.data(snd + 280);
    const auto *snd_281 = buffer.data(snd + 281);
    const auto *snd_282 = buffer.data(snd + 282);
    const auto *snd_285 = buffer.data(snd + 285);
    const auto *snd_286 = buffer.data(snd + 286);
    const auto *snd_287 = buffer.data(snd + 287);
    const auto *snd_288 = buffer.data(snd + 288);
    const auto *snd_291 = buffer.data(snd + 291);
    const auto *snd_292 = buffer.data(snd + 292);
    const auto *snd_293 = buffer.data(snd + 293);
    const auto *snd_294 = buffer.data(snd + 294);
    const auto *snd_297 = buffer.data(snd + 297);
    const auto *snd_298 = buffer.data(snd + 298);
    const auto *snd_299 = buffer.data(snd + 299);
    const auto *snd_300 = buffer.data(snd + 300);
    const auto *snd_303 = buffer.data(snd + 303);

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, smd_174, smd_180, \
                         smd_228, smd_231, snp0_114, snp1_114, snd_228, \
                         snd_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_8 * smd_228[k]
                   + f_1 * snp0_114[k]
                   - f_2 * snp1_114[k]
                   + f_3 * pc_x[k] * snd_228[k];

        t_381[k] = f_11 * smd_180[k]
                   + f_3 * pc_y[k] * snd_228[k];

        t_382[k] = f_8 * smd_174[k]
                   + f_3 * pc_z[k] * snd_228[k];

        t_383[k] = f_8 * smd_231[k]
                   + f_3 * pc_x[k] * snd_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pc_x, pc_y, pc_z, smd_177, smd_183, \
                         smd_232, smd_233, snp0_115, snp1_115, snd_231, snd_232, \
                         snd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_8 * smd_232[k]
                   + f_3 * pc_x[k] * snd_232[k];

        t_385[k] = f_8 * smd_233[k]
                   + f_3 * pc_x[k] * snd_233[k];

        t_386[k] = f_11 * smd_183[k]
                   + f_1 * snp0_115[k]
                   - f_2 * snp1_115[k]
                   + f_3 * pc_y[k] * snd_231[k];

        t_387[k] = f_8 * smd_177[k]
                   + f_3 * pc_z[k] * snd_231[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_x, pc_y, pc_z, smd_179, smd_185, smd_234, \
                         snp0_116, snp0_117, snp1_116, snp1_117, snd_233, \
                         snd_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * smd_185[k]
                   + f_3 * pc_y[k] * snd_233[k];

        t_389[k] = f_8 * smd_179[k]
                   + f_1 * snp0_116[k]
                   - f_2 * snp1_116[k]
                   + f_3 * pc_z[k] * snd_233[k];

        t_390[k] = f_8 * smd_234[k]
                   + f_1 * snp0_117[k]
                   - f_2 * snp1_117[k]
                   + f_3 * pc_x[k] * snd_234[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, smd_180, smd_186, \
                         smd_237, smd_238, snd_234, snd_237, snd_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_13 * smd_186[k]
                   + f_3 * pc_y[k] * snd_234[k];

        t_392[k] = f_10 * smd_180[k]
                   + f_3 * pc_z[k] * snd_234[k];

        t_393[k] = f_8 * smd_237[k]
                   + f_3 * pc_x[k] * snd_237[k];

        t_394[k] = f_8 * smd_238[k]
                   + f_3 * pc_x[k] * snd_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, smd_183, smd_189, \
                         smd_191, smd_239, snp0_118, snp1_118, snd_237, \
                         snd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_8 * smd_239[k]
                   + f_3 * pc_x[k] * snd_239[k];

        t_396[k] = f_13 * smd_189[k]
                   + f_1 * snp0_118[k]
                   - f_2 * snp1_118[k]
                   + f_3 * pc_y[k] * snd_237[k];

        t_397[k] = f_10 * smd_183[k]
                   + f_3 * pc_z[k] * snd_237[k];

        t_398[k] = f_13 * smd_191[k]
                   + f_3 * pc_y[k] * snd_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, smd_185, smd_192, smd_240, \
                         snp0_119, snp0_120, snp1_119, snp1_120, snd_239, \
                         snd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_10 * smd_185[k]
                   + f_1 * snp0_119[k]
                   - f_2 * snp1_119[k]
                   + f_3 * pc_z[k] * snd_239[k];

        t_400[k] = f_8 * smd_240[k]
                   + f_1 * snp0_120[k]
                   - f_2 * snp1_120[k]
                   + f_3 * pc_x[k] * snd_240[k];

        t_401[k] = f_12 * smd_192[k]
                   + f_3 * pc_y[k] * snd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, smd_186, smd_243, smd_244, \
                         smd_245, snd_240, snd_243, snd_244, snd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_12 * smd_186[k]
                   + f_3 * pc_z[k] * snd_240[k];

        t_403[k] = f_8 * smd_243[k]
                   + f_3 * pc_x[k] * snd_243[k];

        t_404[k] = f_8 * smd_244[k]
                   + f_3 * pc_x[k] * snd_244[k];

        t_405[k] = f_8 * smd_245[k]
                   + f_3 * pc_x[k] * snd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pc_y, pc_z, smd_189, smd_191, smd_195, \
                         smd_197, snp0_121, snp0_122, snp1_121, snp1_122, snd_243, \
                         snd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_12 * smd_195[k]
                   + f_1 * snp0_121[k]
                   - f_2 * snp1_121[k]
                   + f_3 * pc_y[k] * snd_243[k];

        t_407[k] = f_12 * smd_189[k]
                   + f_3 * pc_z[k] * snd_243[k];

        t_408[k] = f_12 * smd_197[k]
                   + f_3 * pc_y[k] * snd_245[k];

        t_409[k] = f_12 * smd_191[k]
                   + f_1 * snp0_122[k]
                   - f_2 * snp1_122[k]
                   + f_3 * pc_z[k] * snd_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pc_z, smd_192, smd_198, \
                         smd_246, smd_249, snp0_123, snp1_123, snd_246, \
                         snd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_8 * smd_246[k]
                   + f_1 * snp0_123[k]
                   - f_2 * snp1_123[k]
                   + f_3 * pc_x[k] * snd_246[k];

        t_411[k] = f_10 * smd_198[k]
                   + f_3 * pc_y[k] * snd_246[k];

        t_412[k] = f_13 * smd_192[k]
                   + f_3 * pc_z[k] * snd_246[k];

        t_413[k] = f_8 * smd_249[k]
                   + f_3 * pc_x[k] * snd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, smd_195, smd_201, \
                         smd_250, smd_251, snp0_124, snp1_124, snd_249, snd_250, \
                         snd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_8 * smd_250[k]
                   + f_3 * pc_x[k] * snd_250[k];

        t_415[k] = f_8 * smd_251[k]
                   + f_3 * pc_x[k] * snd_251[k];

        t_416[k] = f_10 * smd_201[k]
                   + f_1 * snp0_124[k]
                   - f_2 * snp1_124[k]
                   + f_3 * pc_y[k] * snd_249[k];

        t_417[k] = f_13 * smd_195[k]
                   + f_3 * pc_z[k] * snd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pc_x, pc_y, pc_z, smd_197, smd_203, smd_252, \
                         snp0_125, snp0_126, snp1_125, snp1_126, snd_251, \
                         snd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_10 * smd_203[k]
                   + f_3 * pc_y[k] * snd_251[k];

        t_419[k] = f_13 * smd_197[k]
                   + f_1 * snp0_125[k]
                   - f_2 * snp1_125[k]
                   + f_3 * pc_z[k] * snd_251[k];

        t_420[k] = f_8 * smd_252[k]
                   + f_1 * snp0_126[k]
                   - f_2 * snp1_126[k]
                   + f_3 * pc_x[k] * snd_252[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, smd_198, smd_204, \
                         smd_255, smd_256, snd_252, snd_255, snd_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_8 * smd_204[k]
                   + f_3 * pc_y[k] * snd_252[k];

        t_422[k] = f_11 * smd_198[k]
                   + f_3 * pc_z[k] * snd_252[k];

        t_423[k] = f_8 * smd_255[k]
                   + f_3 * pc_x[k] * snd_255[k];

        t_424[k] = f_8 * smd_256[k]
                   + f_3 * pc_x[k] * snd_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, smd_201, smd_207, \
                         smd_209, smd_257, snp0_127, snp1_127, snd_255, \
                         snd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_8 * smd_257[k]
                   + f_3 * pc_x[k] * snd_257[k];

        t_426[k] = f_8 * smd_207[k]
                   + f_1 * snp0_127[k]
                   - f_2 * snp1_127[k]
                   + f_3 * pc_y[k] * snd_255[k];

        t_427[k] = f_11 * smd_201[k]
                   + f_3 * pc_z[k] * snd_255[k];

        t_428[k] = f_8 * smd_209[k]
                   + f_3 * pc_y[k] * snd_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pb_y, pc_y, pc_z, smf0_350, smd_203, \
                         smd_204, smd_210, smf1_350, snp0_128, snp1_128, snd_257, \
                         snd_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_11 * smd_203[k]
                   + f_1 * snp0_128[k]
                   - f_2 * snp1_128[k]
                   + f_3 * pc_z[k] * snd_257[k];

        t_430[k] = pb_y[k] * smf0_350[k]
                   - f_4 * pc_y[k] * smf1_350[k];

        t_431[k] = f_5 * smd_210[k]
                   + f_3 * pc_y[k] * snd_258[k];

        t_432[k] = f_9 * smd_204[k]
                   + f_3 * pc_z[k] * snd_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, smd_213, smd_261, smd_262, \
                         smd_263, snp0_130, snp1_130, snd_261, snd_262, \
                         snd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_8 * smd_261[k]
                   + f_3 * pc_x[k] * snd_261[k];

        t_434[k] = f_8 * smd_262[k]
                   + f_3 * pc_x[k] * snd_262[k];

        t_435[k] = f_8 * smd_263[k]
                   + f_3 * pc_x[k] * snd_263[k];

        t_436[k] = f_5 * smd_213[k]
                   + f_1 * snp0_130[k]
                   - f_2 * snp1_130[k]
                   + f_3 * pc_y[k] * snd_261[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pb_y, pc_y, pc_z, smf0_359, smd_207, smd_215, \
                         smf1_359, snd_261, snd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_9 * smd_207[k]
                   + f_3 * pc_z[k] * snd_261[k];

        t_438[k] = f_5 * smd_215[k]
                   + f_3 * pc_y[k] * snd_263[k];

        t_439[k] = pb_y[k] * smf0_359[k]
                   - f_4 * pc_y[k] * smf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, pc_y, pc_z, smd_210, smd_264, \
                         smd_267, snp0_132, snp1_132, snd_264, \
                         snd_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_8 * smd_264[k]
                   + f_1 * snp0_132[k]
                   - f_2 * snp1_132[k]
                   + f_3 * pc_x[k] * snd_264[k];

        t_441[k] = f_3 * pc_y[k] * snd_264[k];

        t_442[k] = f_7 * smd_210[k]
                   + f_3 * pc_z[k] * snd_264[k];

        t_443[k] = f_8 * smd_267[k]
                   + f_3 * pc_x[k] * snd_267[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, pc_x, pc_y, pc_z, smd_213, \
                         smd_268, smd_269, snp0_133, snp1_133, snd_267, snd_268, \
                         snd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_8 * smd_268[k]
                   + f_3 * pc_x[k] * snd_268[k];

        t_445[k] = f_8 * smd_269[k]
                   + f_3 * pc_x[k] * snd_269[k];

        t_446[k] = f_1 * snp0_133[k]
                   - f_2 * snp1_133[k]
                   + f_3 * pc_y[k] * snd_267[k];

        t_447[k] = f_7 * smd_213[k]
                   + f_3 * pc_z[k] * snd_267[k];

        t_448[k] = f_3 * pc_y[k] * snd_269[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pb_x, pc_x, pc_y, pc_z, smf0_450, smd_215, \
                         smd_216, smd_270, smf1_450, snp0_134, snp1_134, snd_269, \
                         snd_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_7 * smd_215[k]
                   + f_1 * snp0_134[k]
                   - f_2 * snp1_134[k]
                   + f_3 * pc_z[k] * snd_269[k];

        t_450[k] = pb_x[k] * smf0_450[k]
                   + f_10 * smd_270[k]
                   - f_4 * pc_x[k] * smf1_450[k];

        t_451[k] = f_6 * smd_216[k]
                   + f_3 * pc_y[k] * snd_270[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, t_455, pc_x, pc_z, smd_273, smd_274, smd_275, \
                         snd_270, snd_273, snd_274, snd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_3 * pc_z[k] * snd_270[k];

        t_453[k] = f_5 * smd_273[k]
                   + f_3 * pc_x[k] * snd_273[k];

        t_454[k] = f_5 * smd_274[k]
                   + f_3 * pc_x[k] * snd_274[k];

        t_455[k] = f_5 * smd_275[k]
                   + f_3 * pc_x[k] * snd_275[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, t_459, pb_x, pc_x, pc_y, pc_z, smf0_456, \
                         smf0_459, smd_221, smf1_456, smf1_459, snd_273, \
                         snd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = pb_x[k] * smf0_456[k]
                   - f_4 * pc_x[k] * smf1_456[k];

        t_457[k] = f_3 * pc_z[k] * snd_273[k];

        t_458[k] = f_6 * smd_221[k]
                   + f_3 * pc_y[k] * snd_275[k];

        t_459[k] = pb_x[k] * smf0_459[k]
                   - f_4 * pc_x[k] * smf1_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pb_z, pc_x, pc_y, pc_z, smf0_360, \
                         smd_216, smd_222, smd_279, smf1_360, snd_276, \
                         snd_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pb_z[k] * smf0_360[k]
                   - f_4 * pc_z[k] * smf1_360[k];

        t_461[k] = f_7 * smd_222[k]
                   + f_3 * pc_y[k] * snd_276[k];

        t_462[k] = f_5 * smd_216[k]
                   + f_3 * pc_z[k] * snd_276[k];

        t_463[k] = f_5 * smd_279[k]
                   + f_3 * pc_x[k] * snd_279[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, t_467, pb_x, pc_x, pc_z, smf0_466, smd_219, \
                         smd_280, smd_281, smf1_466, snd_279, snd_280, \
                         snd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_5 * smd_280[k]
                   + f_3 * pc_x[k] * snd_280[k];

        t_465[k] = f_5 * smd_281[k]
                   + f_3 * pc_x[k] * snd_281[k];

        t_466[k] = pb_x[k] * smf0_466[k]
                   - f_4 * pc_x[k] * smf1_466[k];

        t_467[k] = f_5 * smd_219[k]
                   + f_3 * pc_z[k] * snd_279[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, t_471, pb_x, pc_x, pc_y, smf0_469, smf0_470, \
                         smd_227, smd_228, smd_282, smf1_469, smf1_470, snd_281, \
                         snd_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_7 * smd_227[k]
                   + f_3 * pc_y[k] * snd_281[k];

        t_469[k] = pb_x[k] * smf0_469[k]
                   - f_4 * pc_x[k] * smf1_469[k];

        t_470[k] = pb_x[k] * smf0_470[k]
                   + f_10 * smd_282[k]
                   - f_4 * pc_x[k] * smf1_470[k];

        t_471[k] = f_9 * smd_228[k]
                   + f_3 * pc_y[k] * snd_282[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, pc_x, pc_z, smd_222, smd_285, smd_286, \
                         smd_287, snd_282, snd_285, snd_286, snd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_8 * smd_222[k]
                   + f_3 * pc_z[k] * snd_282[k];

        t_473[k] = f_5 * smd_285[k]
                   + f_3 * pc_x[k] * snd_285[k];

        t_474[k] = f_5 * smd_286[k]
                   + f_3 * pc_x[k] * snd_286[k];

        t_475[k] = f_5 * smd_287[k]
                   + f_3 * pc_x[k] * snd_287[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, pb_x, pc_x, pc_y, pc_z, smf0_476, \
                         smf0_479, smd_225, smd_233, smf1_476, smf1_479, snd_285, \
                         snd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = pb_x[k] * smf0_476[k]
                   - f_4 * pc_x[k] * smf1_476[k];

        t_477[k] = f_8 * smd_225[k]
                   + f_3 * pc_z[k] * snd_285[k];

        t_478[k] = f_9 * smd_233[k]
                   + f_3 * pc_y[k] * snd_287[k];

        t_479[k] = pb_x[k] * smf0_479[k]
                   - f_4 * pc_x[k] * smf1_479[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pb_x, pc_x, pc_y, pc_z, smf0_480, \
                         smd_228, smd_234, smd_288, smd_291, smf1_480, snd_288, \
                         snd_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = pb_x[k] * smf0_480[k]
                   + f_10 * smd_288[k]
                   - f_4 * pc_x[k] * smf1_480[k];

        t_481[k] = f_11 * smd_234[k]
                   + f_3 * pc_y[k] * snd_288[k];

        t_482[k] = f_10 * smd_228[k]
                   + f_3 * pc_z[k] * snd_288[k];

        t_483[k] = f_5 * smd_291[k]
                   + f_3 * pc_x[k] * snd_291[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pb_x, pc_x, pc_z, smf0_486, smd_231, \
                         smd_292, smd_293, smf1_486, snd_291, snd_292, \
                         snd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_5 * smd_292[k]
                   + f_3 * pc_x[k] * snd_292[k];

        t_485[k] = f_5 * smd_293[k]
                   + f_3 * pc_x[k] * snd_293[k];

        t_486[k] = pb_x[k] * smf0_486[k]
                   - f_4 * pc_x[k] * smf1_486[k];

        t_487[k] = f_10 * smd_231[k]
                   + f_3 * pc_z[k] * snd_291[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pb_x, pc_x, pc_y, smf0_489, smf0_490, \
                         smd_239, smd_240, smd_294, smf1_489, smf1_490, snd_293, \
                         snd_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_11 * smd_239[k]
                   + f_3 * pc_y[k] * snd_293[k];

        t_489[k] = pb_x[k] * smf0_489[k]
                   - f_4 * pc_x[k] * smf1_489[k];

        t_490[k] = pb_x[k] * smf0_490[k]
                   + f_10 * smd_294[k]
                   - f_4 * pc_x[k] * smf1_490[k];

        t_491[k] = f_13 * smd_240[k]
                   + f_3 * pc_y[k] * snd_294[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, pc_z, smd_234, smd_297, smd_298, \
                         smd_299, snd_294, snd_297, snd_298, snd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_12 * smd_234[k]
                   + f_3 * pc_z[k] * snd_294[k];

        t_493[k] = f_5 * smd_297[k]
                   + f_3 * pc_x[k] * snd_297[k];

        t_494[k] = f_5 * smd_298[k]
                   + f_3 * pc_x[k] * snd_298[k];

        t_495[k] = f_5 * smd_299[k]
                   + f_3 * pc_x[k] * snd_299[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_x, pc_x, pc_y, pc_z, smf0_496, \
                         smf0_499, smd_237, smd_245, smf1_496, smf1_499, snd_297, \
                         snd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = pb_x[k] * smf0_496[k]
                   - f_4 * pc_x[k] * smf1_496[k];

        t_497[k] = f_12 * smd_237[k]
                   + f_3 * pc_z[k] * snd_297[k];

        t_498[k] = f_13 * smd_245[k]
                   + f_3 * pc_y[k] * snd_299[k];

        t_499[k] = pb_x[k] * smf0_499[k]
                   - f_4 * pc_x[k] * smf1_499[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pb_x, pc_x, pc_y, pc_z, smf0_500, \
                         smd_240, smd_246, smd_300, smd_303, smf1_500, snd_300, \
                         snd_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pb_x[k] * smf0_500[k]
                   + f_10 * smd_300[k]
                   - f_4 * pc_x[k] * smf1_500[k];

        t_501[k] = f_12 * smd_246[k]
                   + f_3 * pc_y[k] * snd_300[k];

        t_502[k] = f_13 * smd_240[k]
                   + f_3 * pc_z[k] * snd_300[k];

        t_503[k] = f_5 * smd_303[k]
                   + f_3 * pc_x[k] * snd_303[k];
    }
}

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.0 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 2.5 / q;

    auto *t_504 = buffer.data(target + 504);
    auto *t_505 = buffer.data(target + 505);
    auto *t_506 = buffer.data(target + 506);
    auto *t_507 = buffer.data(target + 507);
    auto *t_508 = buffer.data(target + 508);
    auto *t_509 = buffer.data(target + 509);
    auto *t_510 = buffer.data(target + 510);
    auto *t_511 = buffer.data(target + 511);
    auto *t_512 = buffer.data(target + 512);
    auto *t_513 = buffer.data(target + 513);
    auto *t_514 = buffer.data(target + 514);
    auto *t_515 = buffer.data(target + 515);
    auto *t_516 = buffer.data(target + 516);
    auto *t_517 = buffer.data(target + 517);
    auto *t_518 = buffer.data(target + 518);
    auto *t_519 = buffer.data(target + 519);
    auto *t_520 = buffer.data(target + 520);
    auto *t_521 = buffer.data(target + 521);
    auto *t_522 = buffer.data(target + 522);
    auto *t_523 = buffer.data(target + 523);
    auto *t_524 = buffer.data(target + 524);
    auto *t_525 = buffer.data(target + 525);
    auto *t_526 = buffer.data(target + 526);
    auto *t_527 = buffer.data(target + 527);
    auto *t_528 = buffer.data(target + 528);
    auto *t_529 = buffer.data(target + 529);
    auto *t_530 = buffer.data(target + 530);
    auto *t_531 = buffer.data(target + 531);
    auto *t_532 = buffer.data(target + 532);
    auto *t_533 = buffer.data(target + 533);
    auto *t_534 = buffer.data(target + 534);
    auto *t_535 = buffer.data(target + 535);
    auto *t_536 = buffer.data(target + 536);
    auto *t_537 = buffer.data(target + 537);
    auto *t_538 = buffer.data(target + 538);
    auto *t_539 = buffer.data(target + 539);
    auto *t_540 = buffer.data(target + 540);
    auto *t_541 = buffer.data(target + 541);
    auto *t_542 = buffer.data(target + 542);
    auto *t_543 = buffer.data(target + 543);
    auto *t_544 = buffer.data(target + 544);
    auto *t_545 = buffer.data(target + 545);
    auto *t_546 = buffer.data(target + 546);
    auto *t_547 = buffer.data(target + 547);
    auto *t_548 = buffer.data(target + 548);
    auto *t_549 = buffer.data(target + 549);
    auto *t_550 = buffer.data(target + 550);
    auto *t_551 = buffer.data(target + 551);
    auto *t_552 = buffer.data(target + 552);
    auto *t_553 = buffer.data(target + 553);
    auto *t_554 = buffer.data(target + 554);
    auto *t_555 = buffer.data(target + 555);
    auto *t_556 = buffer.data(target + 556);
    auto *t_557 = buffer.data(target + 557);
    auto *t_558 = buffer.data(target + 558);
    auto *t_559 = buffer.data(target + 559);
    auto *t_560 = buffer.data(target + 560);
    auto *t_561 = buffer.data(target + 561);
    auto *t_562 = buffer.data(target + 562);
    auto *t_563 = buffer.data(target + 563);
    auto *t_564 = buffer.data(target + 564);
    auto *t_565 = buffer.data(target + 565);
    auto *t_566 = buffer.data(target + 566);
    auto *t_567 = buffer.data(target + 567);
    auto *t_568 = buffer.data(target + 568);
    auto *t_569 = buffer.data(target + 569);
    auto *t_570 = buffer.data(target + 570);
    auto *t_571 = buffer.data(target + 571);
    auto *t_572 = buffer.data(target + 572);
    auto *t_573 = buffer.data(target + 573);
    auto *t_574 = buffer.data(target + 574);
    auto *t_575 = buffer.data(target + 575);
    auto *t_576 = buffer.data(target + 576);
    auto *t_577 = buffer.data(target + 577);
    auto *t_578 = buffer.data(target + 578);
    auto *t_579 = buffer.data(target + 579);
    auto *t_580 = buffer.data(target + 580);
    auto *t_581 = buffer.data(target + 581);
    auto *t_582 = buffer.data(target + 582);
    auto *t_583 = buffer.data(target + 583);
    auto *t_584 = buffer.data(target + 584);
    auto *t_585 = buffer.data(target + 585);
    auto *t_586 = buffer.data(target + 586);
    auto *t_587 = buffer.data(target + 587);
    auto *t_588 = buffer.data(target + 588);
    auto *t_589 = buffer.data(target + 589);
    auto *t_590 = buffer.data(target + 590);
    auto *t_591 = buffer.data(target + 591);
    auto *t_592 = buffer.data(target + 592);
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);
    auto *t_601 = buffer.data(target + 601);
    auto *t_602 = buffer.data(target + 602);
    auto *t_603 = buffer.data(target + 603);
    auto *t_604 = buffer.data(target + 604);
    auto *t_605 = buffer.data(target + 605);
    auto *t_606 = buffer.data(target + 606);
    auto *t_607 = buffer.data(target + 607);
    auto *t_608 = buffer.data(target + 608);
    auto *t_609 = buffer.data(target + 609);
    auto *t_610 = buffer.data(target + 610);
    auto *t_611 = buffer.data(target + 611);
    auto *t_612 = buffer.data(target + 612);
    auto *t_613 = buffer.data(target + 613);
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);
    auto *t_618 = buffer.data(target + 618);
    auto *t_619 = buffer.data(target + 619);
    auto *t_620 = buffer.data(target + 620);
    auto *t_621 = buffer.data(target + 621);
    auto *t_622 = buffer.data(target + 622);
    auto *t_623 = buffer.data(target + 623);
    auto *t_624 = buffer.data(target + 624);
    auto *t_625 = buffer.data(target + 625);
    auto *t_626 = buffer.data(target + 626);
    auto *t_627 = buffer.data(target + 627);
    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_440 = buffer.data(smf0 + 440);
    const auto *smf0_450 = buffer.data(smf0 + 450);
    const auto *smf0_456 = buffer.data(smf0 + 456);
    const auto *smf0_506 = buffer.data(smf0 + 506);
    const auto *smf0_509 = buffer.data(smf0 + 509);
    const auto *smf0_510 = buffer.data(smf0 + 510);
    const auto *smf0_516 = buffer.data(smf0 + 516);
    const auto *smf0_519 = buffer.data(smf0 + 519);
    const auto *smf0_520 = buffer.data(smf0 + 520);
    const auto *smf0_526 = buffer.data(smf0 + 526);
    const auto *smf0_529 = buffer.data(smf0 + 529);
    const auto *smf0_536 = buffer.data(smf0 + 536);
    const auto *smf0_539 = buffer.data(smf0 + 539);
    const auto *smf0_540 = buffer.data(smf0 + 540);
    const auto *smf0_546 = buffer.data(smf0 + 546);
    const auto *smf0_549 = buffer.data(smf0 + 549);

    const auto *smd_243 = buffer.data(smd + 243);
    const auto *smd_246 = buffer.data(smd + 246);
    const auto *smd_249 = buffer.data(smd + 249);
    const auto *smd_251 = buffer.data(smd + 251);
    const auto *smd_252 = buffer.data(smd + 252);
    const auto *smd_255 = buffer.data(smd + 255);
    const auto *smd_257 = buffer.data(smd + 257);
    const auto *smd_258 = buffer.data(smd + 258);
    const auto *smd_261 = buffer.data(smd + 261);
    const auto *smd_263 = buffer.data(smd + 263);
    const auto *smd_264 = buffer.data(smd + 264);
    const auto *smd_267 = buffer.data(smd + 267);
    const auto *smd_269 = buffer.data(smd + 269);
    const auto *smd_270 = buffer.data(smd + 270);
    const auto *smd_273 = buffer.data(smd + 273);
    const auto *smd_275 = buffer.data(smd + 275);
    const auto *smd_276 = buffer.data(smd + 276);
    const auto *smd_279 = buffer.data(smd + 279);
    const auto *smd_281 = buffer.data(smd + 281);
    const auto *smd_282 = buffer.data(smd + 282);
    const auto *smd_285 = buffer.data(smd + 285);
    const auto *smd_287 = buffer.data(smd + 287);
    const auto *smd_288 = buffer.data(smd + 288);
    const auto *smd_291 = buffer.data(smd + 291);
    const auto *smd_293 = buffer.data(smd + 293);
    const auto *smd_294 = buffer.data(smd + 294);
    const auto *smd_297 = buffer.data(smd + 297);
    const auto *smd_299 = buffer.data(smd + 299);
    const auto *smd_300 = buffer.data(smd + 300);
    const auto *smd_303 = buffer.data(smd + 303);
    const auto *smd_304 = buffer.data(smd + 304);
    const auto *smd_305 = buffer.data(smd + 305);
    const auto *smd_306 = buffer.data(smd + 306);
    const auto *smd_309 = buffer.data(smd + 309);
    const auto *smd_310 = buffer.data(smd + 310);
    const auto *smd_311 = buffer.data(smd + 311);
    const auto *smd_312 = buffer.data(smd + 312);
    const auto *smd_315 = buffer.data(smd + 315);
    const auto *smd_316 = buffer.data(smd + 316);
    const auto *smd_317 = buffer.data(smd + 317);
    const auto *smd_318 = buffer.data(smd + 318);
    const auto *smd_321 = buffer.data(smd + 321);
    const auto *smd_322 = buffer.data(smd + 322);
    const auto *smd_323 = buffer.data(smd + 323);
    const auto *smd_324 = buffer.data(smd + 324);
    const auto *smd_327 = buffer.data(smd + 327);
    const auto *smd_328 = buffer.data(smd + 328);
    const auto *smd_329 = buffer.data(smd + 329);

    const auto *smf1_440 = buffer.data(smf1 + 440);
    const auto *smf1_450 = buffer.data(smf1 + 450);
    const auto *smf1_456 = buffer.data(smf1 + 456);
    const auto *smf1_506 = buffer.data(smf1 + 506);
    const auto *smf1_509 = buffer.data(smf1 + 509);
    const auto *smf1_510 = buffer.data(smf1 + 510);
    const auto *smf1_516 = buffer.data(smf1 + 516);
    const auto *smf1_519 = buffer.data(smf1 + 519);
    const auto *smf1_520 = buffer.data(smf1 + 520);
    const auto *smf1_526 = buffer.data(smf1 + 526);
    const auto *smf1_529 = buffer.data(smf1 + 529);
    const auto *smf1_536 = buffer.data(smf1 + 536);
    const auto *smf1_539 = buffer.data(smf1 + 539);
    const auto *smf1_540 = buffer.data(smf1 + 540);
    const auto *smf1_546 = buffer.data(smf1 + 546);
    const auto *smf1_549 = buffer.data(smf1 + 549);

    const auto *snp0_165 = buffer.data(snp0 + 165);
    const auto *snp0_166 = buffer.data(snp0 + 166);
    const auto *snp0_167 = buffer.data(snp0 + 167);
    const auto *snp0_170 = buffer.data(snp0 + 170);
    const auto *snp0_171 = buffer.data(snp0 + 171);
    const auto *snp0_172 = buffer.data(snp0 + 172);
    const auto *snp0_173 = buffer.data(snp0 + 173);
    const auto *snp0_174 = buffer.data(snp0 + 174);
    const auto *snp0_175 = buffer.data(snp0 + 175);
    const auto *snp0_176 = buffer.data(snp0 + 176);
    const auto *snp0_177 = buffer.data(snp0 + 177);
    const auto *snp0_178 = buffer.data(snp0 + 178);
    const auto *snp0_179 = buffer.data(snp0 + 179);
    const auto *snp0_180 = buffer.data(snp0 + 180);
    const auto *snp0_181 = buffer.data(snp0 + 181);
    const auto *snp0_182 = buffer.data(snp0 + 182);
    const auto *snp0_183 = buffer.data(snp0 + 183);
    const auto *snp0_184 = buffer.data(snp0 + 184);
    const auto *snp0_185 = buffer.data(snp0 + 185);
    const auto *snp0_186 = buffer.data(snp0 + 186);
    const auto *snp0_187 = buffer.data(snp0 + 187);
    const auto *snp0_188 = buffer.data(snp0 + 188);
    const auto *snp0_189 = buffer.data(snp0 + 189);
    const auto *snp0_190 = buffer.data(snp0 + 190);

    const auto *snp1_165 = buffer.data(snp1 + 165);
    const auto *snp1_166 = buffer.data(snp1 + 166);
    const auto *snp1_167 = buffer.data(snp1 + 167);
    const auto *snp1_170 = buffer.data(snp1 + 170);
    const auto *snp1_171 = buffer.data(snp1 + 171);
    const auto *snp1_172 = buffer.data(snp1 + 172);
    const auto *snp1_173 = buffer.data(snp1 + 173);
    const auto *snp1_174 = buffer.data(snp1 + 174);
    const auto *snp1_175 = buffer.data(snp1 + 175);
    const auto *snp1_176 = buffer.data(snp1 + 176);
    const auto *snp1_177 = buffer.data(snp1 + 177);
    const auto *snp1_178 = buffer.data(snp1 + 178);
    const auto *snp1_179 = buffer.data(snp1 + 179);
    const auto *snp1_180 = buffer.data(snp1 + 180);
    const auto *snp1_181 = buffer.data(snp1 + 181);
    const auto *snp1_182 = buffer.data(snp1 + 182);
    const auto *snp1_183 = buffer.data(snp1 + 183);
    const auto *snp1_184 = buffer.data(snp1 + 184);
    const auto *snp1_185 = buffer.data(snp1 + 185);
    const auto *snp1_186 = buffer.data(snp1 + 186);
    const auto *snp1_187 = buffer.data(snp1 + 187);
    const auto *snp1_188 = buffer.data(snp1 + 188);
    const auto *snp1_189 = buffer.data(snp1 + 189);
    const auto *snp1_190 = buffer.data(snp1 + 190);

    const auto *snd_303 = buffer.data(snd + 303);
    const auto *snd_304 = buffer.data(snd + 304);
    const auto *snd_305 = buffer.data(snd + 305);
    const auto *snd_306 = buffer.data(snd + 306);
    const auto *snd_309 = buffer.data(snd + 309);
    const auto *snd_310 = buffer.data(snd + 310);
    const auto *snd_311 = buffer.data(snd + 311);
    const auto *snd_312 = buffer.data(snd + 312);
    const auto *snd_315 = buffer.data(snd + 315);
    const auto *snd_316 = buffer.data(snd + 316);
    const auto *snd_317 = buffer.data(snd + 317);
    const auto *snd_318 = buffer.data(snd + 318);
    const auto *snd_321 = buffer.data(snd + 321);
    const auto *snd_322 = buffer.data(snd + 322);
    const auto *snd_323 = buffer.data(snd + 323);
    const auto *snd_324 = buffer.data(snd + 324);
    const auto *snd_327 = buffer.data(snd + 327);
    const auto *snd_328 = buffer.data(snd + 328);
    const auto *snd_329 = buffer.data(snd + 329);
    const auto *snd_330 = buffer.data(snd + 330);
    const auto *snd_333 = buffer.data(snd + 333);
    const auto *snd_334 = buffer.data(snd + 334);
    const auto *snd_335 = buffer.data(snd + 335);
    const auto *snd_336 = buffer.data(snd + 336);
    const auto *snd_339 = buffer.data(snd + 339);
    const auto *snd_340 = buffer.data(snd + 340);
    const auto *snd_341 = buffer.data(snd + 341);
    const auto *snd_342 = buffer.data(snd + 342);
    const auto *snd_345 = buffer.data(snd + 345);
    const auto *snd_346 = buffer.data(snd + 346);
    const auto *snd_347 = buffer.data(snd + 347);
    const auto *snd_348 = buffer.data(snd + 348);
    const auto *snd_351 = buffer.data(snd + 351);
    const auto *snd_352 = buffer.data(snd + 352);
    const auto *snd_353 = buffer.data(snd + 353);
    const auto *snd_354 = buffer.data(snd + 354);
    const auto *snd_357 = buffer.data(snd + 357);
    const auto *snd_358 = buffer.data(snd + 358);
    const auto *snd_359 = buffer.data(snd + 359);
    const auto *snd_360 = buffer.data(snd + 360);
    const auto *snd_363 = buffer.data(snd + 363);
    const auto *snd_364 = buffer.data(snd + 364);
    const auto *snd_365 = buffer.data(snd + 365);
    const auto *snd_366 = buffer.data(snd + 366);
    const auto *snd_369 = buffer.data(snd + 369);
    const auto *snd_370 = buffer.data(snd + 370);
    const auto *snd_371 = buffer.data(snd + 371);
    const auto *snd_372 = buffer.data(snd + 372);
    const auto *snd_375 = buffer.data(snd + 375);
    const auto *snd_376 = buffer.data(snd + 376);
    const auto *snd_377 = buffer.data(snd + 377);
    const auto *snd_378 = buffer.data(snd + 378);
    const auto *snd_381 = buffer.data(snd + 381);
    const auto *snd_382 = buffer.data(snd + 382);
    const auto *snd_383 = buffer.data(snd + 383);

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pb_x, pc_x, pc_z, smf0_506, smd_243, \
                         smd_304, smd_305, smf1_506, snd_303, snd_304, \
                         snd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_5 * smd_304[k]
                   + f_3 * pc_x[k] * snd_304[k];

        t_505[k] = f_5 * smd_305[k]
                   + f_3 * pc_x[k] * snd_305[k];

        t_506[k] = pb_x[k] * smf0_506[k]
                   - f_4 * pc_x[k] * smf1_506[k];

        t_507[k] = f_13 * smd_243[k]
                   + f_3 * pc_z[k] * snd_303[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pb_x, pc_x, pc_y, smf0_509, smf0_510, \
                         smd_251, smd_252, smd_306, smf1_509, smf1_510, snd_305, \
                         snd_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_12 * smd_251[k]
                   + f_3 * pc_y[k] * snd_305[k];

        t_509[k] = pb_x[k] * smf0_509[k]
                   - f_4 * pc_x[k] * smf1_509[k];

        t_510[k] = pb_x[k] * smf0_510[k]
                   + f_10 * smd_306[k]
                   - f_4 * pc_x[k] * smf1_510[k];

        t_511[k] = f_10 * smd_252[k]
                   + f_3 * pc_y[k] * snd_306[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_z, smd_246, smd_309, smd_310, \
                         smd_311, snd_306, snd_309, snd_310, snd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * smd_246[k]
                   + f_3 * pc_z[k] * snd_306[k];

        t_513[k] = f_5 * smd_309[k]
                   + f_3 * pc_x[k] * snd_309[k];

        t_514[k] = f_5 * smd_310[k]
                   + f_3 * pc_x[k] * snd_310[k];

        t_515[k] = f_5 * smd_311[k]
                   + f_3 * pc_x[k] * snd_311[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pb_x, pc_x, pc_y, pc_z, smf0_516, \
                         smf0_519, smd_249, smd_257, smf1_516, smf1_519, snd_309, \
                         snd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = pb_x[k] * smf0_516[k]
                   - f_4 * pc_x[k] * smf1_516[k];

        t_517[k] = f_11 * smd_249[k]
                   + f_3 * pc_z[k] * snd_309[k];

        t_518[k] = f_10 * smd_257[k]
                   + f_3 * pc_y[k] * snd_311[k];

        t_519[k] = pb_x[k] * smf0_519[k]
                   - f_4 * pc_x[k] * smf1_519[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, t_523, pb_x, pc_x, pc_y, pc_z, smf0_520, \
                         smd_252, smd_258, smd_312, smd_315, smf1_520, snd_312, \
                         snd_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = pb_x[k] * smf0_520[k]
                   + f_10 * smd_312[k]
                   - f_4 * pc_x[k] * smf1_520[k];

        t_521[k] = f_8 * smd_258[k]
                   + f_3 * pc_y[k] * snd_312[k];

        t_522[k] = f_9 * smd_252[k]
                   + f_3 * pc_z[k] * snd_312[k];

        t_523[k] = f_5 * smd_315[k]
                   + f_3 * pc_x[k] * snd_315[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, t_527, pb_x, pc_x, pc_z, smf0_526, smd_255, \
                         smd_316, smd_317, smf1_526, snd_315, snd_316, \
                         snd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_5 * smd_316[k]
                   + f_3 * pc_x[k] * snd_316[k];

        t_525[k] = f_5 * smd_317[k]
                   + f_3 * pc_x[k] * snd_317[k];

        t_526[k] = pb_x[k] * smf0_526[k]
                   - f_4 * pc_x[k] * smf1_526[k];

        t_527[k] = f_9 * smd_255[k]
                   + f_3 * pc_z[k] * snd_315[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, t_531, pb_x, pb_y, pc_x, pc_y, smf0_440, \
                         smf0_529, smd_263, smd_264, smf1_440, smf1_529, snd_317, \
                         snd_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = f_8 * smd_263[k]
                   + f_3 * pc_y[k] * snd_317[k];

        t_529[k] = pb_x[k] * smf0_529[k]
                   - f_4 * pc_x[k] * smf1_529[k];

        t_530[k] = pb_y[k] * smf0_440[k]
                   - f_4 * pc_y[k] * smf1_440[k];

        t_531[k] = f_5 * smd_264[k]
                   + f_3 * pc_y[k] * snd_318[k];
    }

#pragma omp simd aligned(t_532, t_533, t_534, t_535, pc_x, pc_z, smd_258, smd_321, smd_322, \
                         smd_323, snd_318, snd_321, snd_322, snd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_532[k] = f_7 * smd_258[k]
                   + f_3 * pc_z[k] * snd_318[k];

        t_533[k] = f_5 * smd_321[k]
                   + f_3 * pc_x[k] * snd_321[k];

        t_534[k] = f_5 * smd_322[k]
                   + f_3 * pc_x[k] * snd_322[k];

        t_535[k] = f_5 * smd_323[k]
                   + f_3 * pc_x[k] * snd_323[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, t_539, pb_x, pc_x, pc_y, pc_z, smf0_536, \
                         smf0_539, smd_261, smd_269, smf1_536, smf1_539, snd_321, \
                         snd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = pb_x[k] * smf0_536[k]
                   - f_4 * pc_x[k] * smf1_536[k];

        t_537[k] = f_7 * smd_261[k]
                   + f_3 * pc_z[k] * snd_321[k];

        t_538[k] = f_5 * smd_269[k]
                   + f_3 * pc_y[k] * snd_323[k];

        t_539[k] = pb_x[k] * smf0_539[k]
                   - f_4 * pc_x[k] * smf1_539[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pb_x, pc_x, pc_y, pc_z, smf0_540, \
                         smd_264, smd_324, smd_327, smf1_540, snd_324, \
                         snd_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = pb_x[k] * smf0_540[k]
                   + f_10 * smd_324[k]
                   - f_4 * pc_x[k] * smf1_540[k];

        t_541[k] = f_3 * pc_y[k] * snd_324[k];

        t_542[k] = f_6 * smd_264[k]
                   + f_3 * pc_z[k] * snd_324[k];

        t_543[k] = f_5 * smd_327[k]
                   + f_3 * pc_x[k] * snd_327[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pb_x, pc_x, pc_z, smf0_546, smd_267, \
                         smd_328, smd_329, smf1_546, snd_327, snd_328, \
                         snd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_5 * smd_328[k]
                   + f_3 * pc_x[k] * snd_328[k];

        t_545[k] = f_5 * smd_329[k]
                   + f_3 * pc_x[k] * snd_329[k];

        t_546[k] = pb_x[k] * smf0_546[k]
                   - f_4 * pc_x[k] * smf1_546[k];

        t_547[k] = f_6 * smd_267[k]
                   + f_3 * pc_z[k] * snd_327[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pb_x, pc_x, pc_y, pc_z, smf0_549, \
                         smd_270, smf1_549, snp0_165, snp1_165, snd_329, \
                         snd_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_3 * pc_y[k] * snd_329[k];

        t_549[k] = pb_x[k] * smf0_549[k]
                   - f_4 * pc_x[k] * smf1_549[k];

        t_550[k] = f_1 * snp0_165[k]
                   - f_2 * snp1_165[k]
                   + f_3 * pc_x[k] * snd_330[k];

        t_551[k] = f_0 * smd_270[k]
                   + f_3 * pc_y[k] * snd_330[k];

        t_552[k] = f_3 * pc_z[k] * snd_330[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, t_557, t_558, pc_x, pc_y, pc_z, smd_273, \
                         smd_275, snp0_166, snp1_166, snd_333, snd_334, \
                         snd_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_3 * pc_x[k] * snd_333[k];

        t_554[k] = f_3 * pc_x[k] * snd_334[k];

        t_555[k] = f_3 * pc_x[k] * snd_335[k];

        t_556[k] = f_0 * smd_273[k]
                   + f_1 * snp0_166[k]
                   - f_2 * snp1_166[k]
                   + f_3 * pc_y[k] * snd_333[k];

        t_557[k] = f_3 * pc_z[k] * snd_333[k];

        t_558[k] = f_0 * smd_275[k]
                   + f_3 * pc_y[k] * snd_335[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pb_z, pc_y, pc_z, smf0_450, smd_270, \
                         smd_276, smf1_450, snp0_167, snp1_167, snd_335, \
                         snd_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_1 * snp0_167[k]
                   - f_2 * snp1_167[k]
                   + f_3 * pc_z[k] * snd_335[k];

        t_560[k] = pb_z[k] * smf0_450[k]
                   - f_4 * pc_z[k] * smf1_450[k];

        t_561[k] = f_6 * smd_276[k]
                   + f_3 * pc_y[k] * snd_336[k];

        t_562[k] = f_5 * smd_270[k]
                   + f_3 * pc_z[k] * snd_336[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, t_567, pb_z, pc_x, pc_z, smf0_456, \
                         smd_273, smf1_456, snd_339, snd_340, snd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_3 * pc_x[k] * snd_339[k];

        t_564[k] = f_3 * pc_x[k] * snd_340[k];

        t_565[k] = f_3 * pc_x[k] * snd_341[k];

        t_566[k] = pb_z[k] * smf0_456[k]
                   - f_4 * pc_z[k] * smf1_456[k];

        t_567[k] = f_5 * smd_273[k]
                   + f_3 * pc_z[k] * snd_339[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, smd_275, smd_281, \
                         smd_282, snp0_170, snp0_171, snp1_170, snp1_171, snd_341, \
                         snd_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_6 * smd_281[k]
                   + f_3 * pc_y[k] * snd_341[k];

        t_569[k] = f_5 * smd_275[k]
                   + f_1 * snp0_170[k]
                   - f_2 * snp1_170[k]
                   + f_3 * pc_z[k] * snd_341[k];

        t_570[k] = f_1 * snp0_171[k]
                   - f_2 * snp1_171[k]
                   + f_3 * pc_x[k] * snd_342[k];

        t_571[k] = f_7 * smd_282[k]
                   + f_3 * pc_y[k] * snd_342[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, pc_x, pc_y, pc_z, smd_276, \
                         smd_285, snp0_172, snp1_172, snd_342, snd_345, snd_346, \
                         snd_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_8 * smd_276[k]
                   + f_3 * pc_z[k] * snd_342[k];

        t_573[k] = f_3 * pc_x[k] * snd_345[k];

        t_574[k] = f_3 * pc_x[k] * snd_346[k];

        t_575[k] = f_3 * pc_x[k] * snd_347[k];

        t_576[k] = f_7 * smd_285[k]
                   + f_1 * snp0_172[k]
                   - f_2 * snp1_172[k]
                   + f_3 * pc_y[k] * snd_345[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, pc_y, pc_z, smd_279, smd_281, smd_287, snp0_173, \
                         snp1_173, snd_345, snd_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_8 * smd_279[k]
                   + f_3 * pc_z[k] * snd_345[k];

        t_578[k] = f_7 * smd_287[k]
                   + f_3 * pc_y[k] * snd_347[k];

        t_579[k] = f_8 * smd_281[k]
                   + f_1 * snp0_173[k]
                   - f_2 * snp1_173[k]
                   + f_3 * pc_z[k] * snd_347[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, t_584, pc_x, pc_y, pc_z, smd_282, \
                         smd_288, snp0_174, snp1_174, snd_348, snd_351, \
                         snd_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_1 * snp0_174[k]
                   - f_2 * snp1_174[k]
                   + f_3 * pc_x[k] * snd_348[k];

        t_581[k] = f_9 * smd_288[k]
                   + f_3 * pc_y[k] * snd_348[k];

        t_582[k] = f_10 * smd_282[k]
                   + f_3 * pc_z[k] * snd_348[k];

        t_583[k] = f_3 * pc_x[k] * snd_351[k];

        t_584[k] = f_3 * pc_x[k] * snd_352[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pc_x, pc_y, pc_z, smd_285, smd_291, \
                         smd_293, snp0_175, snp1_175, snd_351, \
                         snd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_3 * pc_x[k] * snd_353[k];

        t_586[k] = f_9 * smd_291[k]
                   + f_1 * snp0_175[k]
                   - f_2 * snp1_175[k]
                   + f_3 * pc_y[k] * snd_351[k];

        t_587[k] = f_10 * smd_285[k]
                   + f_3 * pc_z[k] * snd_351[k];

        t_588[k] = f_9 * smd_293[k]
                   + f_3 * pc_y[k] * snd_353[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pc_x, pc_y, pc_z, smd_287, smd_288, \
                         smd_294, snp0_176, snp0_177, snp1_176, snp1_177, snd_353, \
                         snd_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_10 * smd_287[k]
                   + f_1 * snp0_176[k]
                   - f_2 * snp1_176[k]
                   + f_3 * pc_z[k] * snd_353[k];

        t_590[k] = f_1 * snp0_177[k]
                   - f_2 * snp1_177[k]
                   + f_3 * pc_x[k] * snd_354[k];

        t_591[k] = f_11 * smd_294[k]
                   + f_3 * pc_y[k] * snd_354[k];

        t_592[k] = f_12 * smd_288[k]
                   + f_3 * pc_z[k] * snd_354[k];
    }

#pragma omp simd aligned(t_593, t_594, t_595, t_596, t_597, pc_x, pc_y, pc_z, smd_291, \
                         smd_297, snp0_178, snp1_178, snd_357, snd_358, \
                         snd_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = f_3 * pc_x[k] * snd_357[k];

        t_594[k] = f_3 * pc_x[k] * snd_358[k];

        t_595[k] = f_3 * pc_x[k] * snd_359[k];

        t_596[k] = f_11 * smd_297[k]
                   + f_1 * snp0_178[k]
                   - f_2 * snp1_178[k]
                   + f_3 * pc_y[k] * snd_357[k];

        t_597[k] = f_12 * smd_291[k]
                   + f_3 * pc_z[k] * snd_357[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pc_x, pc_y, pc_z, smd_293, smd_299, \
                         smd_300, snp0_179, snp0_180, snp1_179, snp1_180, snd_359, \
                         snd_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_11 * smd_299[k]
                   + f_3 * pc_y[k] * snd_359[k];

        t_599[k] = f_12 * smd_293[k]
                   + f_1 * snp0_179[k]
                   - f_2 * snp1_179[k]
                   + f_3 * pc_z[k] * snd_359[k];

        t_600[k] = f_1 * snp0_180[k]
                   - f_2 * snp1_180[k]
                   + f_3 * pc_x[k] * snd_360[k];

        t_601[k] = f_13 * smd_300[k]
                   + f_3 * pc_y[k] * snd_360[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, pc_x, pc_y, pc_z, smd_294, \
                         smd_303, snp0_181, snp1_181, snd_360, snd_363, snd_364, \
                         snd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_13 * smd_294[k]
                   + f_3 * pc_z[k] * snd_360[k];

        t_603[k] = f_3 * pc_x[k] * snd_363[k];

        t_604[k] = f_3 * pc_x[k] * snd_364[k];

        t_605[k] = f_3 * pc_x[k] * snd_365[k];

        t_606[k] = f_13 * smd_303[k]
                   + f_1 * snp0_181[k]
                   - f_2 * snp1_181[k]
                   + f_3 * pc_y[k] * snd_363[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pc_y, pc_z, smd_297, smd_299, smd_305, snp0_182, \
                         snp1_182, snd_363, snd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_13 * smd_297[k]
                   + f_3 * pc_z[k] * snd_363[k];

        t_608[k] = f_13 * smd_305[k]
                   + f_3 * pc_y[k] * snd_365[k];

        t_609[k] = f_13 * smd_299[k]
                   + f_1 * snp0_182[k]
                   - f_2 * snp1_182[k]
                   + f_3 * pc_z[k] * snd_365[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, t_614, pc_x, pc_y, pc_z, smd_300, \
                         smd_306, snp0_183, snp1_183, snd_366, snd_369, \
                         snd_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_1 * snp0_183[k]
                   - f_2 * snp1_183[k]
                   + f_3 * pc_x[k] * snd_366[k];

        t_611[k] = f_12 * smd_306[k]
                   + f_3 * pc_y[k] * snd_366[k];

        t_612[k] = f_11 * smd_300[k]
                   + f_3 * pc_z[k] * snd_366[k];

        t_613[k] = f_3 * pc_x[k] * snd_369[k];

        t_614[k] = f_3 * pc_x[k] * snd_370[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, t_618, pc_x, pc_y, pc_z, smd_303, smd_309, \
                         smd_311, snp0_184, snp1_184, snd_369, \
                         snd_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_3 * pc_x[k] * snd_371[k];

        t_616[k] = f_12 * smd_309[k]
                   + f_1 * snp0_184[k]
                   - f_2 * snp1_184[k]
                   + f_3 * pc_y[k] * snd_369[k];

        t_617[k] = f_11 * smd_303[k]
                   + f_3 * pc_z[k] * snd_369[k];

        t_618[k] = f_12 * smd_311[k]
                   + f_3 * pc_y[k] * snd_371[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, smd_305, smd_306, \
                         smd_312, snp0_185, snp0_186, snp1_185, snp1_186, snd_371, \
                         snd_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_11 * smd_305[k]
                   + f_1 * snp0_185[k]
                   - f_2 * snp1_185[k]
                   + f_3 * pc_z[k] * snd_371[k];

        t_620[k] = f_1 * snp0_186[k]
                   - f_2 * snp1_186[k]
                   + f_3 * pc_x[k] * snd_372[k];

        t_621[k] = f_10 * smd_312[k]
                   + f_3 * pc_y[k] * snd_372[k];

        t_622[k] = f_9 * smd_306[k]
                   + f_3 * pc_z[k] * snd_372[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, t_626, t_627, pc_x, pc_y, pc_z, smd_309, \
                         smd_315, snp0_187, snp1_187, snd_375, snd_376, \
                         snd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_3 * pc_x[k] * snd_375[k];

        t_624[k] = f_3 * pc_x[k] * snd_376[k];

        t_625[k] = f_3 * pc_x[k] * snd_377[k];

        t_626[k] = f_10 * smd_315[k]
                   + f_1 * snp0_187[k]
                   - f_2 * snp1_187[k]
                   + f_3 * pc_y[k] * snd_375[k];

        t_627[k] = f_9 * smd_309[k]
                   + f_3 * pc_z[k] * snd_375[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pc_x, pc_y, pc_z, smd_311, smd_317, \
                         smd_318, snp0_188, snp0_189, snp1_188, snp1_189, snd_377, \
                         snd_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_10 * smd_317[k]
                   + f_3 * pc_y[k] * snd_377[k];

        t_629[k] = f_9 * smd_311[k]
                   + f_1 * snp0_188[k]
                   - f_2 * snp1_188[k]
                   + f_3 * pc_z[k] * snd_377[k];

        t_630[k] = f_1 * snp0_189[k]
                   - f_2 * snp1_189[k]
                   + f_3 * pc_x[k] * snd_378[k];

        t_631[k] = f_8 * smd_318[k]
                   + f_3 * pc_y[k] * snd_378[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, t_636, pc_x, pc_y, pc_z, smd_312, \
                         smd_321, snp0_190, snp1_190, snd_378, snd_381, snd_382, \
                         snd_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_7 * smd_312[k]
                   + f_3 * pc_z[k] * snd_378[k];

        t_633[k] = f_3 * pc_x[k] * snd_381[k];

        t_634[k] = f_3 * pc_x[k] * snd_382[k];

        t_635[k] = f_3 * pc_x[k] * snd_383[k];

        t_636[k] = f_8 * smd_321[k]
                   + f_1 * snp0_190[k]
                   - f_2 * snp1_190[k]
                   + f_3 * pc_y[k] * snd_381[k];
    }
}

static auto
compute_prim_snf_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t smf0,
                                                          const size_t smd, const size_t smf1,
                                                          const size_t snp0, const size_t snp1,
                                                          const size_t snd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 4.5 / q;
    const auto f_7 = 4.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.5 / q;

    auto *t_637 = buffer.data(target + 637);
    auto *t_638 = buffer.data(target + 638);
    auto *t_639 = buffer.data(target + 639);
    auto *t_640 = buffer.data(target + 640);
    auto *t_641 = buffer.data(target + 641);
    auto *t_642 = buffer.data(target + 642);
    auto *t_643 = buffer.data(target + 643);
    auto *t_644 = buffer.data(target + 644);
    auto *t_645 = buffer.data(target + 645);
    auto *t_646 = buffer.data(target + 646);
    auto *t_647 = buffer.data(target + 647);
    auto *t_648 = buffer.data(target + 648);
    auto *t_649 = buffer.data(target + 649);
    auto *t_650 = buffer.data(target + 650);
    auto *t_651 = buffer.data(target + 651);
    auto *t_652 = buffer.data(target + 652);
    auto *t_653 = buffer.data(target + 653);
    auto *t_654 = buffer.data(target + 654);
    auto *t_655 = buffer.data(target + 655);
    auto *t_656 = buffer.data(target + 656);
    auto *t_657 = buffer.data(target + 657);
    auto *t_658 = buffer.data(target + 658);
    auto *t_659 = buffer.data(target + 659);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *smf0_540 = buffer.data(smf0 + 540);
    const auto *smf0_546 = buffer.data(smf0 + 546);
    const auto *smf0_549 = buffer.data(smf0 + 549);

    const auto *smd_315 = buffer.data(smd + 315);
    const auto *smd_317 = buffer.data(smd + 317);
    const auto *smd_318 = buffer.data(smd + 318);
    const auto *smd_321 = buffer.data(smd + 321);
    const auto *smd_323 = buffer.data(smd + 323);
    const auto *smd_324 = buffer.data(smd + 324);
    const auto *smd_327 = buffer.data(smd + 327);
    const auto *smd_329 = buffer.data(smd + 329);

    const auto *smf1_540 = buffer.data(smf1 + 540);
    const auto *smf1_546 = buffer.data(smf1 + 546);
    const auto *smf1_549 = buffer.data(smf1 + 549);

    const auto *snp0_191 = buffer.data(snp0 + 191);
    const auto *snp0_195 = buffer.data(snp0 + 195);
    const auto *snp0_196 = buffer.data(snp0 + 196);
    const auto *snp0_197 = buffer.data(snp0 + 197);

    const auto *snp1_191 = buffer.data(snp1 + 191);
    const auto *snp1_195 = buffer.data(snp1 + 195);
    const auto *snp1_196 = buffer.data(snp1 + 196);
    const auto *snp1_197 = buffer.data(snp1 + 197);

    const auto *snd_381 = buffer.data(snd + 381);
    const auto *snd_383 = buffer.data(snd + 383);
    const auto *snd_384 = buffer.data(snd + 384);
    const auto *snd_387 = buffer.data(snd + 387);
    const auto *snd_388 = buffer.data(snd + 388);
    const auto *snd_389 = buffer.data(snd + 389);
    const auto *snd_390 = buffer.data(snd + 390);
    const auto *snd_393 = buffer.data(snd + 393);
    const auto *snd_394 = buffer.data(snd + 394);
    const auto *snd_395 = buffer.data(snd + 395);

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pb_y, pc_y, pc_z, smf0_540, smd_315, \
                         smd_317, smd_323, smf1_540, snp0_191, snp1_191, snd_381, \
                         snd_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_7 * smd_315[k]
                   + f_3 * pc_z[k] * snd_381[k];

        t_638[k] = f_8 * smd_323[k]
                   + f_3 * pc_y[k] * snd_383[k];

        t_639[k] = f_7 * smd_317[k]
                   + f_1 * snp0_191[k]
                   - f_2 * snp1_191[k]
                   + f_3 * pc_z[k] * snd_383[k];

        t_640[k] = pb_y[k] * smf0_540[k]
                   - f_4 * pc_y[k] * smf1_540[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, t_645, pc_x, pc_y, pc_z, smd_318, \
                         smd_324, snd_384, snd_387, snd_388, snd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_5 * smd_324[k]
                   + f_3 * pc_y[k] * snd_384[k];

        t_642[k] = f_6 * smd_318[k]
                   + f_3 * pc_z[k] * snd_384[k];

        t_643[k] = f_3 * pc_x[k] * snd_387[k];

        t_644[k] = f_3 * pc_x[k] * snd_388[k];

        t_645[k] = f_3 * pc_x[k] * snd_389[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pb_y, pc_y, pc_z, smf0_546, smf0_549, \
                         smd_321, smd_327, smd_329, smf1_546, smf1_549, snd_387, \
                         snd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = pb_y[k] * smf0_546[k]
                   + f_10 * smd_327[k]
                   - f_4 * pc_y[k] * smf1_546[k];

        t_647[k] = f_6 * smd_321[k]
                   + f_3 * pc_z[k] * snd_387[k];

        t_648[k] = f_5 * smd_329[k]
                   + f_3 * pc_y[k] * snd_389[k];

        t_649[k] = pb_y[k] * smf0_549[k]
                   - f_4 * pc_y[k] * smf1_549[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, t_654, t_655, pc_x, pc_y, pc_z, smd_324, \
                         snp0_195, snp1_195, snd_390, snd_393, snd_394, \
                         snd_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_1 * snp0_195[k]
                   - f_2 * snp1_195[k]
                   + f_3 * pc_x[k] * snd_390[k];

        t_651[k] = f_3 * pc_y[k] * snd_390[k];

        t_652[k] = f_0 * smd_324[k]
                   + f_3 * pc_z[k] * snd_390[k];

        t_653[k] = f_3 * pc_x[k] * snd_393[k];

        t_654[k] = f_3 * pc_x[k] * snd_394[k];

        t_655[k] = f_3 * pc_x[k] * snd_395[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pc_y, pc_z, smd_327, smd_329, snp0_196, \
                         snp0_197, snp1_196, snp1_197, snd_393, \
                         snd_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_1 * snp0_196[k]
                   - f_2 * snp1_196[k]
                   + f_3 * pc_y[k] * snd_393[k];

        t_657[k] = f_0 * smd_327[k]
                   + f_3 * pc_z[k] * snd_393[k];

        t_658[k] = f_3 * pc_y[k] * snd_395[k];

        t_659[k] = f_0 * smd_329[k]
                   + f_1 * snp0_197[k]
                   - f_2 * snp1_197[k]
                   + f_3 * pc_z[k] * snd_395[k];
    }
}

auto
compute_prim_snf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t smf0, const size_t smd,
                                                   const size_t smf1, const size_t snp0,
                                                   const size_t snp1, const size_t snd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_snf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);

    compute_prim_snf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);

    compute_prim_snf_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);

    compute_prim_snf_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);

    compute_prim_snf_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);

    compute_prim_snf_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, smf0, smd,
                                                              smf1, snp0, snp1, snd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
