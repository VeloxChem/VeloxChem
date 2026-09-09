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


#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
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

    const auto *snf0_0 = buffer.data(snf0 + 0);
    const auto *snf0_6 = buffer.data(snf0 + 6);
    const auto *snf0_9 = buffer.data(snf0 + 9);
    const auto *snf0_16 = buffer.data(snf0 + 16);
    const auto *snf0_20 = buffer.data(snf0 + 20);
    const auto *snf0_29 = buffer.data(snf0 + 29);
    const auto *snf0_30 = buffer.data(snf0 + 30);
    const auto *snf0_36 = buffer.data(snf0 + 36);
    const auto *snf0_50 = buffer.data(snf0 + 50);
    const auto *snf0_59 = buffer.data(snf0 + 59);
    const auto *snf0_60 = buffer.data(snf0 + 60);
    const auto *snf0_66 = buffer.data(snf0 + 66);

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
    const auto *snd_69 = buffer.data(snd + 69);
    const auto *snd_70 = buffer.data(snd + 70);
    const auto *snd_71 = buffer.data(snd + 71);
    const auto *snd_72 = buffer.data(snd + 72);
    const auto *snd_75 = buffer.data(snd + 75);
    const auto *snd_76 = buffer.data(snd + 76);
    const auto *snd_77 = buffer.data(snd + 77);

    const auto *snf1_0 = buffer.data(snf1 + 0);
    const auto *snf1_6 = buffer.data(snf1 + 6);
    const auto *snf1_9 = buffer.data(snf1 + 9);
    const auto *snf1_16 = buffer.data(snf1 + 16);
    const auto *snf1_20 = buffer.data(snf1 + 20);
    const auto *snf1_29 = buffer.data(snf1 + 29);
    const auto *snf1_30 = buffer.data(snf1 + 30);
    const auto *snf1_36 = buffer.data(snf1 + 36);
    const auto *snf1_50 = buffer.data(snf1 + 50);
    const auto *snf1_59 = buffer.data(snf1 + 59);
    const auto *snf1_60 = buffer.data(snf1 + 60);
    const auto *snf1_66 = buffer.data(snf1 + 66);

    const auto *sop0_0 = buffer.data(sop0 + 0);
    const auto *sop0_1 = buffer.data(sop0 + 1);
    const auto *sop0_2 = buffer.data(sop0 + 2);
    const auto *sop0_4 = buffer.data(sop0 + 4);
    const auto *sop0_8 = buffer.data(sop0 + 8);
    const auto *sop0_9 = buffer.data(sop0 + 9);
    const auto *sop0_10 = buffer.data(sop0 + 10);
    const auto *sop0_11 = buffer.data(sop0 + 11);
    const auto *sop0_15 = buffer.data(sop0 + 15);
    const auto *sop0_16 = buffer.data(sop0 + 16);
    const auto *sop0_17 = buffer.data(sop0 + 17);
    const auto *sop0_18 = buffer.data(sop0 + 18);
    const auto *sop0_19 = buffer.data(sop0 + 19);
    const auto *sop0_20 = buffer.data(sop0 + 20);
    const auto *sop0_23 = buffer.data(sop0 + 23);
    const auto *sop0_25 = buffer.data(sop0 + 25);
    const auto *sop0_27 = buffer.data(sop0 + 27);
    const auto *sop0_28 = buffer.data(sop0 + 28);
    const auto *sop0_29 = buffer.data(sop0 + 29);
    const auto *sop0_30 = buffer.data(sop0 + 30);
    const auto *sop0_31 = buffer.data(sop0 + 31);
    const auto *sop0_32 = buffer.data(sop0 + 32);
    const auto *sop0_35 = buffer.data(sop0 + 35);
    const auto *sop0_36 = buffer.data(sop0 + 36);
    const auto *sop0_37 = buffer.data(sop0 + 37);
    const auto *sop0_38 = buffer.data(sop0 + 38);

    const auto *sop1_0 = buffer.data(sop1 + 0);
    const auto *sop1_1 = buffer.data(sop1 + 1);
    const auto *sop1_2 = buffer.data(sop1 + 2);
    const auto *sop1_4 = buffer.data(sop1 + 4);
    const auto *sop1_8 = buffer.data(sop1 + 8);
    const auto *sop1_9 = buffer.data(sop1 + 9);
    const auto *sop1_10 = buffer.data(sop1 + 10);
    const auto *sop1_11 = buffer.data(sop1 + 11);
    const auto *sop1_15 = buffer.data(sop1 + 15);
    const auto *sop1_16 = buffer.data(sop1 + 16);
    const auto *sop1_17 = buffer.data(sop1 + 17);
    const auto *sop1_18 = buffer.data(sop1 + 18);
    const auto *sop1_19 = buffer.data(sop1 + 19);
    const auto *sop1_20 = buffer.data(sop1 + 20);
    const auto *sop1_23 = buffer.data(sop1 + 23);
    const auto *sop1_25 = buffer.data(sop1 + 25);
    const auto *sop1_27 = buffer.data(sop1 + 27);
    const auto *sop1_28 = buffer.data(sop1 + 28);
    const auto *sop1_29 = buffer.data(sop1 + 29);
    const auto *sop1_30 = buffer.data(sop1 + 30);
    const auto *sop1_31 = buffer.data(sop1 + 31);
    const auto *sop1_32 = buffer.data(sop1 + 32);
    const auto *sop1_35 = buffer.data(sop1 + 35);
    const auto *sop1_36 = buffer.data(sop1 + 36);
    const auto *sop1_37 = buffer.data(sop1 + 37);
    const auto *sop1_38 = buffer.data(sop1 + 38);

    const auto *sod_0 = buffer.data(sod + 0);
    const auto *sod_3 = buffer.data(sod + 3);
    const auto *sod_4 = buffer.data(sod + 4);
    const auto *sod_5 = buffer.data(sod + 5);
    const auto *sod_6 = buffer.data(sod + 6);
    const auto *sod_9 = buffer.data(sod + 9);
    const auto *sod_10 = buffer.data(sod + 10);
    const auto *sod_11 = buffer.data(sod + 11);
    const auto *sod_12 = buffer.data(sod + 12);
    const auto *sod_15 = buffer.data(sod + 15);
    const auto *sod_16 = buffer.data(sod + 16);
    const auto *sod_17 = buffer.data(sod + 17);
    const auto *sod_18 = buffer.data(sod + 18);
    const auto *sod_21 = buffer.data(sod + 21);
    const auto *sod_22 = buffer.data(sod + 22);
    const auto *sod_23 = buffer.data(sod + 23);
    const auto *sod_24 = buffer.data(sod + 24);
    const auto *sod_27 = buffer.data(sod + 27);
    const auto *sod_28 = buffer.data(sod + 28);
    const auto *sod_29 = buffer.data(sod + 29);
    const auto *sod_30 = buffer.data(sod + 30);
    const auto *sod_33 = buffer.data(sod + 33);
    const auto *sod_34 = buffer.data(sod + 34);
    const auto *sod_35 = buffer.data(sod + 35);
    const auto *sod_36 = buffer.data(sod + 36);
    const auto *sod_39 = buffer.data(sod + 39);
    const auto *sod_40 = buffer.data(sod + 40);
    const auto *sod_41 = buffer.data(sod + 41);
    const auto *sod_42 = buffer.data(sod + 42);
    const auto *sod_45 = buffer.data(sod + 45);
    const auto *sod_46 = buffer.data(sod + 46);
    const auto *sod_47 = buffer.data(sod + 47);
    const auto *sod_48 = buffer.data(sod + 48);
    const auto *sod_51 = buffer.data(sod + 51);
    const auto *sod_52 = buffer.data(sod + 52);
    const auto *sod_53 = buffer.data(sod + 53);
    const auto *sod_54 = buffer.data(sod + 54);
    const auto *sod_57 = buffer.data(sod + 57);
    const auto *sod_58 = buffer.data(sod + 58);
    const auto *sod_59 = buffer.data(sod + 59);
    const auto *sod_60 = buffer.data(sod + 60);
    const auto *sod_63 = buffer.data(sod + 63);
    const auto *sod_64 = buffer.data(sod + 64);
    const auto *sod_65 = buffer.data(sod + 65);
    const auto *sod_66 = buffer.data(sod + 66);
    const auto *sod_69 = buffer.data(sod + 69);
    const auto *sod_70 = buffer.data(sod + 70);
    const auto *sod_71 = buffer.data(sod + 71);
    const auto *sod_72 = buffer.data(sod + 72);
    const auto *sod_75 = buffer.data(sod + 75);
    const auto *sod_76 = buffer.data(sod + 76);
    const auto *sod_77 = buffer.data(sod + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, snd_0, snd_3, snd_4, \
                         sop0_0, sop1_0, sod_0, sod_3, sod_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * snd_0[k]
                 + f_1 * sop0_0[k]
                 - f_2 * sop1_0[k]
                 + f_3 * pc_x[k] * sod_0[k];

        t_1[k] = f_3 * pc_y[k] * sod_0[k];

        t_2[k] = f_3 * pc_z[k] * sod_0[k];

        t_3[k] = f_0 * snd_3[k]
                 + f_3 * pc_x[k] * sod_3[k];

        t_4[k] = f_0 * snd_4[k]
                 + f_3 * pc_x[k] * sod_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, snd_5, sop0_1, sop0_2, \
                         sop1_1, sop1_2, sod_3, sod_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * snd_5[k]
                 + f_3 * pc_x[k] * sod_5[k];

        t_6[k] = f_1 * sop0_1[k]
                 - f_2 * sop1_1[k]
                 + f_3 * pc_y[k] * sod_3[k];

        t_7[k] = f_3 * pc_z[k] * sod_3[k];

        t_8[k] = f_3 * pc_y[k] * sod_5[k];

        t_9[k] = f_1 * sop0_2[k]
                 - f_2 * sop1_2[k]
                 + f_3 * pc_z[k] * sod_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, snf0_0, snd_0, snd_9, \
                         snf1_0, sod_6, sod_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * snf0_0[k]
                  - f_4 * pc_y[k] * snf1_0[k];

        t_11[k] = f_5 * snd_0[k]
                  + f_3 * pc_y[k] * sod_6[k];

        t_12[k] = f_3 * pc_z[k] * sod_6[k];

        t_13[k] = f_6 * snd_9[k]
                  + f_3 * pc_x[k] * sod_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, snd_3, snd_10, snd_11, \
                         sop0_4, sop1_4, sod_9, sod_10, sod_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * snd_10[k]
                  + f_3 * pc_x[k] * sod_10[k];

        t_15[k] = f_6 * snd_11[k]
                  + f_3 * pc_x[k] * sod_11[k];

        t_16[k] = f_5 * snd_3[k]
                  + f_1 * sop0_4[k]
                  - f_2 * sop1_4[k]
                  + f_3 * pc_y[k] * sod_9[k];

        t_17[k] = f_3 * pc_z[k] * sod_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, snf0_0, snf0_9, \
                         snd_5, snf1_0, snf1_9, sod_11, sod_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * snd_5[k]
                  + f_3 * pc_y[k] * sod_11[k];

        t_19[k] = pb_y[k] * snf0_9[k]
                  - f_4 * pc_y[k] * snf1_9[k];

        t_20[k] = pb_z[k] * snf0_0[k]
                  - f_4 * pc_z[k] * snf1_0[k];

        t_21[k] = f_3 * pc_y[k] * sod_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, snd_0, snd_15, snd_16, snd_17, \
                         sod_12, sod_15, sod_16, sod_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * snd_0[k]
                  + f_3 * pc_z[k] * sod_12[k];

        t_23[k] = f_6 * snd_15[k]
                  + f_3 * pc_x[k] * sod_15[k];

        t_24[k] = f_6 * snd_16[k]
                  + f_3 * pc_x[k] * sod_16[k];

        t_25[k] = f_6 * snd_17[k]
                  + f_3 * pc_x[k] * sod_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, snf0_6, snd_3, snd_5, \
                         snf1_6, sop0_8, sop1_8, sod_15, sod_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * snf0_6[k]
                  - f_4 * pc_z[k] * snf1_6[k];

        t_27[k] = f_5 * snd_3[k]
                  + f_3 * pc_z[k] * sod_15[k];

        t_28[k] = f_3 * pc_y[k] * sod_17[k];

        t_29[k] = f_5 * snd_5[k]
                  + f_1 * sop0_8[k]
                  - f_2 * sop1_8[k]
                  + f_3 * pc_z[k] * sod_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, snd_6, snd_18, snd_21, \
                         sop0_9, sop1_9, sod_18, sod_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * snd_18[k]
                  + f_1 * sop0_9[k]
                  - f_2 * sop1_9[k]
                  + f_3 * pc_x[k] * sod_18[k];

        t_31[k] = f_8 * snd_6[k]
                  + f_3 * pc_y[k] * sod_18[k];

        t_32[k] = f_3 * pc_z[k] * sod_18[k];

        t_33[k] = f_7 * snd_21[k]
                  + f_3 * pc_x[k] * sod_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, snd_9, snd_22, snd_23, \
                         sop0_10, sop1_10, sod_21, sod_22, sod_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * snd_22[k]
                  + f_3 * pc_x[k] * sod_22[k];

        t_35[k] = f_7 * snd_23[k]
                  + f_3 * pc_x[k] * sod_23[k];

        t_36[k] = f_8 * snd_9[k]
                  + f_1 * sop0_10[k]
                  - f_2 * sop1_10[k]
                  + f_3 * pc_y[k] * sod_21[k];

        t_37[k] = f_3 * pc_z[k] * sod_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, snf0_20, snd_11, snd_12, \
                         snf1_20, sop0_11, sop1_11, sod_23, sod_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * snd_11[k]
                  + f_3 * pc_y[k] * sod_23[k];

        t_39[k] = f_1 * sop0_11[k]
                  - f_2 * sop1_11[k]
                  + f_3 * pc_z[k] * sod_23[k];

        t_40[k] = pb_y[k] * snf0_20[k]
                  - f_4 * pc_y[k] * snf1_20[k];

        t_41[k] = f_5 * snd_12[k]
                  + f_3 * pc_y[k] * sod_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, snd_6, snd_27, snd_28, snd_29, \
                         sod_24, sod_27, sod_28, sod_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * snd_6[k]
                  + f_3 * pc_z[k] * sod_24[k];

        t_43[k] = f_7 * snd_27[k]
                  + f_3 * pc_x[k] * sod_27[k];

        t_44[k] = f_7 * snd_28[k]
                  + f_3 * pc_x[k] * sod_28[k];

        t_45[k] = f_7 * snd_29[k]
                  + f_3 * pc_x[k] * sod_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, snf0_16, snf0_29, \
                         snd_9, snd_17, snf1_16, snf1_29, sod_27, \
                         sod_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * snf0_16[k]
                  - f_4 * pc_z[k] * snf1_16[k];

        t_47[k] = f_5 * snd_9[k]
                  + f_3 * pc_z[k] * sod_27[k];

        t_48[k] = f_5 * snd_17[k]
                  + f_3 * pc_y[k] * sod_29[k];

        t_49[k] = pb_y[k] * snf0_29[k]
                  - f_4 * pc_y[k] * snf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, snd_12, snd_30, snd_33, \
                         sop0_15, sop1_15, sod_30, sod_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * snd_30[k]
                  + f_1 * sop0_15[k]
                  - f_2 * sop1_15[k]
                  + f_3 * pc_x[k] * sod_30[k];

        t_51[k] = f_3 * pc_y[k] * sod_30[k];

        t_52[k] = f_8 * snd_12[k]
                  + f_3 * pc_z[k] * sod_30[k];

        t_53[k] = f_7 * snd_33[k]
                  + f_3 * pc_x[k] * sod_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, snd_15, snd_34, \
                         snd_35, sop0_16, sop1_16, sod_33, sod_34, \
                         sod_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * snd_34[k]
                  + f_3 * pc_x[k] * sod_34[k];

        t_55[k] = f_7 * snd_35[k]
                  + f_3 * pc_x[k] * sod_35[k];

        t_56[k] = f_1 * sop0_16[k]
                  - f_2 * sop1_16[k]
                  + f_3 * pc_y[k] * sod_33[k];

        t_57[k] = f_8 * snd_15[k]
                  + f_3 * pc_z[k] * sod_33[k];

        t_58[k] = f_3 * pc_y[k] * sod_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, snd_17, snd_18, snd_36, \
                         sop0_17, sop0_18, sop1_17, sop1_18, sod_35, \
                         sod_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * snd_17[k]
                  + f_1 * sop0_17[k]
                  - f_2 * sop1_17[k]
                  + f_3 * pc_z[k] * sod_35[k];

        t_60[k] = f_9 * snd_36[k]
                  + f_1 * sop0_18[k]
                  - f_2 * sop1_18[k]
                  + f_3 * pc_x[k] * sod_36[k];

        t_61[k] = f_10 * snd_18[k]
                  + f_3 * pc_y[k] * sod_36[k];

        t_62[k] = f_3 * pc_z[k] * sod_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, snd_21, snd_39, snd_40, snd_41, \
                         sop0_19, sop1_19, sod_39, sod_40, sod_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * snd_39[k]
                  + f_3 * pc_x[k] * sod_39[k];

        t_64[k] = f_9 * snd_40[k]
                  + f_3 * pc_x[k] * sod_40[k];

        t_65[k] = f_9 * snd_41[k]
                  + f_3 * pc_x[k] * sod_41[k];

        t_66[k] = f_10 * snd_21[k]
                  + f_1 * sop0_19[k]
                  - f_2 * sop1_19[k]
                  + f_3 * pc_y[k] * sod_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, snf0_30, snd_23, snf1_30, \
                         sop0_20, sop1_20, sod_39, sod_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * sod_39[k];

        t_68[k] = f_10 * snd_23[k]
                  + f_3 * pc_y[k] * sod_41[k];

        t_69[k] = f_1 * sop0_20[k]
                  - f_2 * sop1_20[k]
                  + f_3 * pc_z[k] * sod_41[k];

        t_70[k] = pb_z[k] * snf0_30[k]
                  - f_4 * pc_z[k] * snf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, snd_18, snd_24, snd_45, \
                         snd_46, sod_42, sod_45, sod_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * snd_24[k]
                  + f_3 * pc_y[k] * sod_42[k];

        t_72[k] = f_5 * snd_18[k]
                  + f_3 * pc_z[k] * sod_42[k];

        t_73[k] = f_9 * snd_45[k]
                  + f_3 * pc_x[k] * sod_45[k];

        t_74[k] = f_9 * snd_46[k]
                  + f_3 * pc_x[k] * sod_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, snf0_36, snd_21, \
                         snd_29, snd_47, snf1_36, sod_45, sod_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_9 * snd_47[k]
                  + f_3 * pc_x[k] * sod_47[k];

        t_76[k] = pb_z[k] * snf0_36[k]
                  - f_4 * pc_z[k] * snf1_36[k];

        t_77[k] = f_5 * snd_21[k]
                  + f_3 * pc_z[k] * sod_45[k];

        t_78[k] = f_8 * snd_29[k]
                  + f_3 * pc_y[k] * sod_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, snf0_50, snd_23, snd_24, \
                         snd_30, snf1_50, sop0_23, sop1_23, sod_47, \
                         sod_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * snd_23[k]
                  + f_1 * sop0_23[k]
                  - f_2 * sop1_23[k]
                  + f_3 * pc_z[k] * sod_47[k];

        t_80[k] = pb_y[k] * snf0_50[k]
                  - f_4 * pc_y[k] * snf1_50[k];

        t_81[k] = f_5 * snd_30[k]
                  + f_3 * pc_y[k] * sod_48[k];

        t_82[k] = f_8 * snd_24[k]
                  + f_3 * pc_z[k] * sod_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, snd_33, snd_51, snd_52, snd_53, \
                         sop0_25, sop1_25, sod_51, sod_52, sod_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * snd_51[k]
                  + f_3 * pc_x[k] * sod_51[k];

        t_84[k] = f_9 * snd_52[k]
                  + f_3 * pc_x[k] * sod_52[k];

        t_85[k] = f_9 * snd_53[k]
                  + f_3 * pc_x[k] * sod_53[k];

        t_86[k] = f_5 * snd_33[k]
                  + f_1 * sop0_25[k]
                  - f_2 * sop1_25[k]
                  + f_3 * pc_y[k] * sod_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, snf0_59, snd_27, snd_35, snf1_59, \
                         sod_51, sod_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * snd_27[k]
                  + f_3 * pc_z[k] * sod_51[k];

        t_88[k] = f_5 * snd_35[k]
                  + f_3 * pc_y[k] * sod_53[k];

        t_89[k] = pb_y[k] * snf0_59[k]
                  - f_4 * pc_y[k] * snf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, snd_30, snd_54, snd_57, \
                         sop0_27, sop1_27, sod_54, sod_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_9 * snd_54[k]
                  + f_1 * sop0_27[k]
                  - f_2 * sop1_27[k]
                  + f_3 * pc_x[k] * sod_54[k];

        t_91[k] = f_3 * pc_y[k] * sod_54[k];

        t_92[k] = f_10 * snd_30[k]
                  + f_3 * pc_z[k] * sod_54[k];

        t_93[k] = f_9 * snd_57[k]
                  + f_3 * pc_x[k] * sod_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, snd_33, snd_58, \
                         snd_59, sop0_28, sop1_28, sod_57, sod_58, \
                         sod_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_9 * snd_58[k]
                  + f_3 * pc_x[k] * sod_58[k];

        t_95[k] = f_9 * snd_59[k]
                  + f_3 * pc_x[k] * sod_59[k];

        t_96[k] = f_1 * sop0_28[k]
                  - f_2 * sop1_28[k]
                  + f_3 * pc_y[k] * sod_57[k];

        t_97[k] = f_10 * snd_33[k]
                  + f_3 * pc_z[k] * sod_57[k];

        t_98[k] = f_3 * pc_y[k] * sod_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_y, pc_z, snd_35, snd_36, snd_60, \
                         sop0_29, sop0_30, sop1_29, sop1_30, sod_59, \
                         sod_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_10 * snd_35[k]
                  + f_1 * sop0_29[k]
                  - f_2 * sop1_29[k]
                  + f_3 * pc_z[k] * sod_59[k];

        t_100[k] = f_11 * snd_60[k]
                   + f_1 * sop0_30[k]
                   - f_2 * sop1_30[k]
                   + f_3 * pc_x[k] * sod_60[k];

        t_101[k] = f_12 * snd_36[k]
                   + f_3 * pc_y[k] * sod_60[k];

        t_102[k] = f_3 * pc_z[k] * sod_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pc_x, pc_y, snd_39, snd_63, snd_64, \
                         snd_65, sop0_31, sop1_31, sod_63, sod_64, \
                         sod_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_11 * snd_63[k]
                   + f_3 * pc_x[k] * sod_63[k];

        t_104[k] = f_11 * snd_64[k]
                   + f_3 * pc_x[k] * sod_64[k];

        t_105[k] = f_11 * snd_65[k]
                   + f_3 * pc_x[k] * sod_65[k];

        t_106[k] = f_12 * snd_39[k]
                   + f_1 * sop0_31[k]
                   - f_2 * sop1_31[k]
                   + f_3 * pc_y[k] * sod_63[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pb_z, pc_y, pc_z, snf0_60, snd_41, \
                         snf1_60, sop0_32, sop1_32, sod_63, sod_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * pc_z[k] * sod_63[k];

        t_108[k] = f_12 * snd_41[k]
                   + f_3 * pc_y[k] * sod_65[k];

        t_109[k] = f_1 * sop0_32[k]
                   - f_2 * sop1_32[k]
                   + f_3 * pc_z[k] * sod_65[k];

        t_110[k] = pb_z[k] * snf0_60[k]
                   - f_4 * pc_z[k] * snf1_60[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, snd_36, snd_42, snd_69, \
                         snd_70, sod_66, sod_69, sod_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * snd_42[k]
                   + f_3 * pc_y[k] * sod_66[k];

        t_112[k] = f_5 * snd_36[k]
                   + f_3 * pc_z[k] * sod_66[k];

        t_113[k] = f_11 * snd_69[k]
                   + f_3 * pc_x[k] * sod_69[k];

        t_114[k] = f_11 * snd_70[k]
                   + f_3 * pc_x[k] * sod_70[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pb_z, pc_x, pc_y, pc_z, snf0_66, snd_39, \
                         snd_47, snd_71, snf1_66, sod_69, sod_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_11 * snd_71[k]
                   + f_3 * pc_x[k] * sod_71[k];

        t_116[k] = pb_z[k] * snf0_66[k]
                   - f_4 * pc_z[k] * snf1_66[k];

        t_117[k] = f_5 * snd_39[k]
                   + f_3 * pc_z[k] * sod_69[k];

        t_118[k] = f_10 * snd_47[k]
                   + f_3 * pc_y[k] * sod_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_y, pc_z, snd_41, snd_48, snd_72, \
                         sop0_35, sop0_36, sop1_35, sop1_36, sod_71, \
                         sod_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * snd_41[k]
                   + f_1 * sop0_35[k]
                   - f_2 * sop1_35[k]
                   + f_3 * pc_z[k] * sod_71[k];

        t_120[k] = f_11 * snd_72[k]
                   + f_1 * sop0_36[k]
                   - f_2 * sop1_36[k]
                   + f_3 * pc_x[k] * sod_72[k];

        t_121[k] = f_8 * snd_48[k]
                   + f_3 * pc_y[k] * sod_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, snd_42, snd_75, snd_76, \
                         snd_77, sod_72, sod_75, sod_76, sod_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * snd_42[k]
                   + f_3 * pc_z[k] * sod_72[k];

        t_123[k] = f_11 * snd_75[k]
                   + f_3 * pc_x[k] * sod_75[k];

        t_124[k] = f_11 * snd_76[k]
                   + f_3 * pc_x[k] * sod_76[k];

        t_125[k] = f_11 * snd_77[k]
                   + f_3 * pc_x[k] * sod_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_y, pc_z, snd_45, snd_47, snd_51, \
                         snd_53, sop0_37, sop0_38, sop1_37, sop1_38, sod_75, \
                         sod_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_8 * snd_51[k]
                   + f_1 * sop0_37[k]
                   - f_2 * sop1_37[k]
                   + f_3 * pc_y[k] * sod_75[k];

        t_127[k] = f_8 * snd_45[k]
                   + f_3 * pc_z[k] * sod_75[k];

        t_128[k] = f_8 * snd_53[k]
                   + f_3 * pc_y[k] * sod_77[k];

        t_129[k] = f_8 * snd_47[k]
                   + f_1 * sop0_38[k]
                   - f_2 * sop1_38[k]
                   + f_3 * pc_z[k] * sod_77[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
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
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *snf0_90 = buffer.data(snf0 + 90);
    const auto *snf0_99 = buffer.data(snf0 + 99);
    const auto *snf0_100 = buffer.data(snf0 + 100);
    const auto *snf0_106 = buffer.data(snf0 + 106);
    const auto *snf0_140 = buffer.data(snf0 + 140);
    const auto *snf0_149 = buffer.data(snf0 + 149);
    const auto *snf0_150 = buffer.data(snf0 + 150);
    const auto *snf0_156 = buffer.data(snf0 + 156);

    const auto *snd_48 = buffer.data(snd + 48);
    const auto *snd_51 = buffer.data(snd + 51);
    const auto *snd_54 = buffer.data(snd + 54);
    const auto *snd_57 = buffer.data(snd + 57);
    const auto *snd_59 = buffer.data(snd + 59);
    const auto *snd_60 = buffer.data(snd + 60);
    const auto *snd_63 = buffer.data(snd + 63);
    const auto *snd_65 = buffer.data(snd + 65);
    const auto *snd_66 = buffer.data(snd + 66);
    const auto *snd_69 = buffer.data(snd + 69);
    const auto *snd_71 = buffer.data(snd + 71);
    const auto *snd_72 = buffer.data(snd + 72);
    const auto *snd_75 = buffer.data(snd + 75);
    const auto *snd_77 = buffer.data(snd + 77);
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

    const auto *snf1_90 = buffer.data(snf1 + 90);
    const auto *snf1_99 = buffer.data(snf1 + 99);
    const auto *snf1_100 = buffer.data(snf1 + 100);
    const auto *snf1_106 = buffer.data(snf1 + 106);
    const auto *snf1_140 = buffer.data(snf1 + 140);
    const auto *snf1_149 = buffer.data(snf1 + 149);
    const auto *snf1_150 = buffer.data(snf1 + 150);
    const auto *snf1_156 = buffer.data(snf1 + 156);

    const auto *sop0_40 = buffer.data(sop0 + 40);
    const auto *sop0_42 = buffer.data(sop0 + 42);
    const auto *sop0_43 = buffer.data(sop0 + 43);
    const auto *sop0_44 = buffer.data(sop0 + 44);
    const auto *sop0_45 = buffer.data(sop0 + 45);
    const auto *sop0_46 = buffer.data(sop0 + 46);
    const auto *sop0_47 = buffer.data(sop0 + 47);
    const auto *sop0_50 = buffer.data(sop0 + 50);
    const auto *sop0_51 = buffer.data(sop0 + 51);
    const auto *sop0_52 = buffer.data(sop0 + 52);
    const auto *sop0_53 = buffer.data(sop0 + 53);
    const auto *sop0_54 = buffer.data(sop0 + 54);
    const auto *sop0_55 = buffer.data(sop0 + 55);
    const auto *sop0_56 = buffer.data(sop0 + 56);
    const auto *sop0_58 = buffer.data(sop0 + 58);
    const auto *sop0_60 = buffer.data(sop0 + 60);
    const auto *sop0_61 = buffer.data(sop0 + 61);
    const auto *sop0_62 = buffer.data(sop0 + 62);
    const auto *sop0_63 = buffer.data(sop0 + 63);
    const auto *sop0_64 = buffer.data(sop0 + 64);
    const auto *sop0_65 = buffer.data(sop0 + 65);
    const auto *sop0_68 = buffer.data(sop0 + 68);
    const auto *sop0_69 = buffer.data(sop0 + 69);
    const auto *sop0_70 = buffer.data(sop0 + 70);
    const auto *sop0_71 = buffer.data(sop0 + 71);
    const auto *sop0_72 = buffer.data(sop0 + 72);
    const auto *sop0_73 = buffer.data(sop0 + 73);
    const auto *sop0_74 = buffer.data(sop0 + 74);
    const auto *sop0_75 = buffer.data(sop0 + 75);

    const auto *sop1_40 = buffer.data(sop1 + 40);
    const auto *sop1_42 = buffer.data(sop1 + 42);
    const auto *sop1_43 = buffer.data(sop1 + 43);
    const auto *sop1_44 = buffer.data(sop1 + 44);
    const auto *sop1_45 = buffer.data(sop1 + 45);
    const auto *sop1_46 = buffer.data(sop1 + 46);
    const auto *sop1_47 = buffer.data(sop1 + 47);
    const auto *sop1_50 = buffer.data(sop1 + 50);
    const auto *sop1_51 = buffer.data(sop1 + 51);
    const auto *sop1_52 = buffer.data(sop1 + 52);
    const auto *sop1_53 = buffer.data(sop1 + 53);
    const auto *sop1_54 = buffer.data(sop1 + 54);
    const auto *sop1_55 = buffer.data(sop1 + 55);
    const auto *sop1_56 = buffer.data(sop1 + 56);
    const auto *sop1_58 = buffer.data(sop1 + 58);
    const auto *sop1_60 = buffer.data(sop1 + 60);
    const auto *sop1_61 = buffer.data(sop1 + 61);
    const auto *sop1_62 = buffer.data(sop1 + 62);
    const auto *sop1_63 = buffer.data(sop1 + 63);
    const auto *sop1_64 = buffer.data(sop1 + 64);
    const auto *sop1_65 = buffer.data(sop1 + 65);
    const auto *sop1_68 = buffer.data(sop1 + 68);
    const auto *sop1_69 = buffer.data(sop1 + 69);
    const auto *sop1_70 = buffer.data(sop1 + 70);
    const auto *sop1_71 = buffer.data(sop1 + 71);
    const auto *sop1_72 = buffer.data(sop1 + 72);
    const auto *sop1_73 = buffer.data(sop1 + 73);
    const auto *sop1_74 = buffer.data(sop1 + 74);
    const auto *sop1_75 = buffer.data(sop1 + 75);

    const auto *sod_78 = buffer.data(sod + 78);
    const auto *sod_81 = buffer.data(sod + 81);
    const auto *sod_82 = buffer.data(sod + 82);
    const auto *sod_83 = buffer.data(sod + 83);
    const auto *sod_84 = buffer.data(sod + 84);
    const auto *sod_87 = buffer.data(sod + 87);
    const auto *sod_88 = buffer.data(sod + 88);
    const auto *sod_89 = buffer.data(sod + 89);
    const auto *sod_90 = buffer.data(sod + 90);
    const auto *sod_93 = buffer.data(sod + 93);
    const auto *sod_94 = buffer.data(sod + 94);
    const auto *sod_95 = buffer.data(sod + 95);
    const auto *sod_96 = buffer.data(sod + 96);
    const auto *sod_99 = buffer.data(sod + 99);
    const auto *sod_100 = buffer.data(sod + 100);
    const auto *sod_101 = buffer.data(sod + 101);
    const auto *sod_102 = buffer.data(sod + 102);
    const auto *sod_105 = buffer.data(sod + 105);
    const auto *sod_106 = buffer.data(sod + 106);
    const auto *sod_107 = buffer.data(sod + 107);
    const auto *sod_108 = buffer.data(sod + 108);
    const auto *sod_111 = buffer.data(sod + 111);
    const auto *sod_112 = buffer.data(sod + 112);
    const auto *sod_113 = buffer.data(sod + 113);
    const auto *sod_114 = buffer.data(sod + 114);
    const auto *sod_117 = buffer.data(sod + 117);
    const auto *sod_118 = buffer.data(sod + 118);
    const auto *sod_119 = buffer.data(sod + 119);
    const auto *sod_120 = buffer.data(sod + 120);
    const auto *sod_123 = buffer.data(sod + 123);
    const auto *sod_124 = buffer.data(sod + 124);
    const auto *sod_125 = buffer.data(sod + 125);
    const auto *sod_126 = buffer.data(sod + 126);
    const auto *sod_129 = buffer.data(sod + 129);
    const auto *sod_130 = buffer.data(sod + 130);
    const auto *sod_131 = buffer.data(sod + 131);
    const auto *sod_132 = buffer.data(sod + 132);
    const auto *sod_135 = buffer.data(sod + 135);
    const auto *sod_136 = buffer.data(sod + 136);
    const auto *sod_137 = buffer.data(sod + 137);
    const auto *sod_138 = buffer.data(sod + 138);
    const auto *sod_141 = buffer.data(sod + 141);
    const auto *sod_142 = buffer.data(sod + 142);
    const auto *sod_143 = buffer.data(sod + 143);
    const auto *sod_144 = buffer.data(sod + 144);
    const auto *sod_147 = buffer.data(sod + 147);
    const auto *sod_148 = buffer.data(sod + 148);
    const auto *sod_149 = buffer.data(sod + 149);
    const auto *sod_150 = buffer.data(sod + 150);
    const auto *sod_153 = buffer.data(sod + 153);
    const auto *sod_154 = buffer.data(sod + 154);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, snf0_90, snd_48, \
                         snd_54, snd_81, snf1_90, sod_78, sod_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * snf0_90[k]
                   - f_4 * pc_y[k] * snf1_90[k];

        t_131[k] = f_5 * snd_54[k]
                   + f_3 * pc_y[k] * sod_78[k];

        t_132[k] = f_10 * snd_48[k]
                   + f_3 * pc_z[k] * sod_78[k];

        t_133[k] = f_11 * snd_81[k]
                   + f_3 * pc_x[k] * sod_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, snd_51, snd_57, snd_82, \
                         snd_83, sop0_40, sop1_40, sod_81, sod_82, \
                         sod_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * snd_82[k]
                   + f_3 * pc_x[k] * sod_82[k];

        t_135[k] = f_11 * snd_83[k]
                   + f_3 * pc_x[k] * sod_83[k];

        t_136[k] = f_5 * snd_57[k]
                   + f_1 * sop0_40[k]
                   - f_2 * sop1_40[k]
                   + f_3 * pc_y[k] * sod_81[k];

        t_137[k] = f_10 * snd_51[k]
                   + f_3 * pc_z[k] * sod_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_y, pc_x, pc_y, snf0_99, snd_59, \
                         snd_84, snf1_99, sop0_42, sop1_42, sod_83, \
                         sod_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * snd_59[k]
                   + f_3 * pc_y[k] * sod_83[k];

        t_139[k] = pb_y[k] * snf0_99[k]
                   - f_4 * pc_y[k] * snf1_99[k];

        t_140[k] = f_11 * snd_84[k]
                   + f_1 * sop0_42[k]
                   - f_2 * sop1_42[k]
                   + f_3 * pc_x[k] * sod_84[k];

        t_141[k] = f_3 * pc_y[k] * sod_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, snd_54, snd_87, snd_88, \
                         snd_89, sod_84, sod_87, sod_88, sod_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_12 * snd_54[k]
                   + f_3 * pc_z[k] * sod_84[k];

        t_143[k] = f_11 * snd_87[k]
                   + f_3 * pc_x[k] * sod_87[k];

        t_144[k] = f_11 * snd_88[k]
                   + f_3 * pc_x[k] * sod_88[k];

        t_145[k] = f_11 * snd_89[k]
                   + f_3 * pc_x[k] * sod_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pc_y, pc_z, snd_57, snd_59, sop0_43, \
                         sop0_44, sop1_43, sop1_44, sod_87, sod_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * sop0_43[k]
                   - f_2 * sop1_43[k]
                   + f_3 * pc_y[k] * sod_87[k];

        t_147[k] = f_12 * snd_57[k]
                   + f_3 * pc_z[k] * sod_87[k];

        t_148[k] = f_3 * pc_y[k] * sod_89[k];

        t_149[k] = f_12 * snd_59[k]
                   + f_1 * sop0_44[k]
                   - f_2 * sop1_44[k]
                   + f_3 * pc_z[k] * sod_89[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pc_x, pc_y, pc_z, snd_60, snd_90, snd_93, \
                         sop0_45, sop1_45, sod_90, sod_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_13 * snd_90[k]
                   + f_1 * sop0_45[k]
                   - f_2 * sop1_45[k]
                   + f_3 * pc_x[k] * sod_90[k];

        t_151[k] = f_14 * snd_60[k]
                   + f_3 * pc_y[k] * sod_90[k];

        t_152[k] = f_3 * pc_z[k] * sod_90[k];

        t_153[k] = f_13 * snd_93[k]
                   + f_3 * pc_x[k] * sod_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, snd_63, snd_94, snd_95, \
                         sop0_46, sop1_46, sod_93, sod_94, sod_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_13 * snd_94[k]
                   + f_3 * pc_x[k] * sod_94[k];

        t_155[k] = f_13 * snd_95[k]
                   + f_3 * pc_x[k] * sod_95[k];

        t_156[k] = f_14 * snd_63[k]
                   + f_1 * sop0_46[k]
                   - f_2 * sop1_46[k]
                   + f_3 * pc_y[k] * sod_93[k];

        t_157[k] = f_3 * pc_z[k] * sod_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pb_z, pc_y, pc_z, snf0_100, snd_65, \
                         snd_66, snf1_100, sop0_47, sop1_47, sod_95, \
                         sod_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_14 * snd_65[k]
                   + f_3 * pc_y[k] * sod_95[k];

        t_159[k] = f_1 * sop0_47[k]
                   - f_2 * sop1_47[k]
                   + f_3 * pc_z[k] * sod_95[k];

        t_160[k] = pb_z[k] * snf0_100[k]
                   - f_4 * pc_z[k] * snf1_100[k];

        t_161[k] = f_12 * snd_66[k]
                   + f_3 * pc_y[k] * sod_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, snd_60, snd_99, snd_100, \
                         snd_101, sod_96, sod_99, sod_100, sod_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * snd_60[k]
                   + f_3 * pc_z[k] * sod_96[k];

        t_163[k] = f_13 * snd_99[k]
                   + f_3 * pc_x[k] * sod_99[k];

        t_164[k] = f_13 * snd_100[k]
                   + f_3 * pc_x[k] * sod_100[k];

        t_165[k] = f_13 * snd_101[k]
                   + f_3 * pc_x[k] * sod_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, snf0_106, snd_63, \
                         snd_65, snd_71, snf1_106, sop0_50, sop1_50, sod_99, \
                         sod_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_z[k] * snf0_106[k]
                   - f_4 * pc_z[k] * snf1_106[k];

        t_167[k] = f_5 * snd_63[k]
                   + f_3 * pc_z[k] * sod_99[k];

        t_168[k] = f_12 * snd_71[k]
                   + f_3 * pc_y[k] * sod_101[k];

        t_169[k] = f_5 * snd_65[k]
                   + f_1 * sop0_50[k]
                   - f_2 * sop1_50[k]
                   + f_3 * pc_z[k] * sod_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, snd_66, snd_72, \
                         snd_102, snd_105, sop0_51, sop1_51, sod_102, \
                         sod_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_13 * snd_102[k]
                   + f_1 * sop0_51[k]
                   - f_2 * sop1_51[k]
                   + f_3 * pc_x[k] * sod_102[k];

        t_171[k] = f_10 * snd_72[k]
                   + f_3 * pc_y[k] * sod_102[k];

        t_172[k] = f_8 * snd_66[k]
                   + f_3 * pc_z[k] * sod_102[k];

        t_173[k] = f_13 * snd_105[k]
                   + f_3 * pc_x[k] * sod_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, snd_69, snd_75, \
                         snd_106, snd_107, sop0_52, sop1_52, sod_105, sod_106, \
                         sod_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_13 * snd_106[k]
                   + f_3 * pc_x[k] * sod_106[k];

        t_175[k] = f_13 * snd_107[k]
                   + f_3 * pc_x[k] * sod_107[k];

        t_176[k] = f_10 * snd_75[k]
                   + f_1 * sop0_52[k]
                   - f_2 * sop1_52[k]
                   + f_3 * pc_y[k] * sod_105[k];

        t_177[k] = f_8 * snd_69[k]
                   + f_3 * pc_z[k] * sod_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, snd_71, snd_77, snd_108, \
                         sop0_53, sop0_54, sop1_53, sop1_54, sod_107, \
                         sod_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * snd_77[k]
                   + f_3 * pc_y[k] * sod_107[k];

        t_179[k] = f_8 * snd_71[k]
                   + f_1 * sop0_53[k]
                   - f_2 * sop1_53[k]
                   + f_3 * pc_z[k] * sod_107[k];

        t_180[k] = f_13 * snd_108[k]
                   + f_1 * sop0_54[k]
                   - f_2 * sop1_54[k]
                   + f_3 * pc_x[k] * sod_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, snd_72, snd_78, \
                         snd_111, snd_112, sod_108, sod_111, sod_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * snd_78[k]
                   + f_3 * pc_y[k] * sod_108[k];

        t_182[k] = f_10 * snd_72[k]
                   + f_3 * pc_z[k] * sod_108[k];

        t_183[k] = f_13 * snd_111[k]
                   + f_3 * pc_x[k] * sod_111[k];

        t_184[k] = f_13 * snd_112[k]
                   + f_3 * pc_x[k] * sod_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, snd_75, snd_81, snd_83, \
                         snd_113, sop0_55, sop1_55, sod_111, sod_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_13 * snd_113[k]
                   + f_3 * pc_x[k] * sod_113[k];

        t_186[k] = f_8 * snd_81[k]
                   + f_1 * sop0_55[k]
                   - f_2 * sop1_55[k]
                   + f_3 * pc_y[k] * sod_111[k];

        t_187[k] = f_10 * snd_75[k]
                   + f_3 * pc_z[k] * sod_111[k];

        t_188[k] = f_8 * snd_83[k]
                   + f_3 * pc_y[k] * sod_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_y, pc_y, pc_z, snf0_140, snd_77, \
                         snd_78, snd_84, snf1_140, sop0_56, sop1_56, sod_113, \
                         sod_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * snd_77[k]
                   + f_1 * sop0_56[k]
                   - f_2 * sop1_56[k]
                   + f_3 * pc_z[k] * sod_113[k];

        t_190[k] = pb_y[k] * snf0_140[k]
                   - f_4 * pc_y[k] * snf1_140[k];

        t_191[k] = f_5 * snd_84[k]
                   + f_3 * pc_y[k] * sod_114[k];

        t_192[k] = f_12 * snd_78[k]
                   + f_3 * pc_z[k] * sod_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, snd_87, snd_117, snd_118, \
                         snd_119, sop0_58, sop1_58, sod_117, sod_118, \
                         sod_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_13 * snd_117[k]
                   + f_3 * pc_x[k] * sod_117[k];

        t_194[k] = f_13 * snd_118[k]
                   + f_3 * pc_x[k] * sod_118[k];

        t_195[k] = f_13 * snd_119[k]
                   + f_3 * pc_x[k] * sod_119[k];

        t_196[k] = f_5 * snd_87[k]
                   + f_1 * sop0_58[k]
                   - f_2 * sop1_58[k]
                   + f_3 * pc_y[k] * sod_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_y, pc_y, pc_z, snf0_149, snd_81, snd_89, \
                         snf1_149, sod_117, sod_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * snd_81[k]
                   + f_3 * pc_z[k] * sod_117[k];

        t_198[k] = f_5 * snd_89[k]
                   + f_3 * pc_y[k] * sod_119[k];

        t_199[k] = pb_y[k] * snf0_149[k]
                   - f_4 * pc_y[k] * snf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, snd_84, snd_120, \
                         snd_123, sop0_60, sop1_60, sod_120, sod_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_13 * snd_120[k]
                   + f_1 * sop0_60[k]
                   - f_2 * sop1_60[k]
                   + f_3 * pc_x[k] * sod_120[k];

        t_201[k] = f_3 * pc_y[k] * sod_120[k];

        t_202[k] = f_14 * snd_84[k]
                   + f_3 * pc_z[k] * sod_120[k];

        t_203[k] = f_13 * snd_123[k]
                   + f_3 * pc_x[k] * sod_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, t_208, pc_x, pc_y, pc_z, snd_87, snd_124, \
                         snd_125, sop0_61, sop1_61, sod_123, sod_124, \
                         sod_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_13 * snd_124[k]
                   + f_3 * pc_x[k] * sod_124[k];

        t_205[k] = f_13 * snd_125[k]
                   + f_3 * pc_x[k] * sod_125[k];

        t_206[k] = f_1 * sop0_61[k]
                   - f_2 * sop1_61[k]
                   + f_3 * pc_y[k] * sod_123[k];

        t_207[k] = f_14 * snd_87[k]
                   + f_3 * pc_z[k] * sod_123[k];

        t_208[k] = f_3 * pc_y[k] * sod_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, snd_89, snd_90, \
                         snd_126, sop0_62, sop0_63, sop1_62, sop1_63, sod_125, \
                         sod_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_14 * snd_89[k]
                   + f_1 * sop0_62[k]
                   - f_2 * sop1_62[k]
                   + f_3 * pc_z[k] * sod_125[k];

        t_210[k] = f_14 * snd_126[k]
                   + f_1 * sop0_63[k]
                   - f_2 * sop1_63[k]
                   + f_3 * pc_x[k] * sod_126[k];

        t_211[k] = f_13 * snd_90[k]
                   + f_3 * pc_y[k] * sod_126[k];

        t_212[k] = f_3 * pc_z[k] * sod_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, snd_93, snd_129, snd_130, \
                         snd_131, sop0_64, sop1_64, sod_129, sod_130, \
                         sod_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_14 * snd_129[k]
                   + f_3 * pc_x[k] * sod_129[k];

        t_214[k] = f_14 * snd_130[k]
                   + f_3 * pc_x[k] * sod_130[k];

        t_215[k] = f_14 * snd_131[k]
                   + f_3 * pc_x[k] * sod_131[k];

        t_216[k] = f_13 * snd_93[k]
                   + f_1 * sop0_64[k]
                   - f_2 * sop1_64[k]
                   + f_3 * pc_y[k] * sod_129[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pb_z, pc_y, pc_z, snf0_150, snd_95, \
                         snf1_150, sop0_65, sop1_65, sod_129, sod_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_3 * pc_z[k] * sod_129[k];

        t_218[k] = f_13 * snd_95[k]
                   + f_3 * pc_y[k] * sod_131[k];

        t_219[k] = f_1 * sop0_65[k]
                   - f_2 * sop1_65[k]
                   + f_3 * pc_z[k] * sod_131[k];

        t_220[k] = pb_z[k] * snf0_150[k]
                   - f_4 * pc_z[k] * snf1_150[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pc_x, pc_y, pc_z, snd_90, snd_96, \
                         snd_135, snd_136, sod_132, sod_135, sod_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_14 * snd_96[k]
                   + f_3 * pc_y[k] * sod_132[k];

        t_222[k] = f_5 * snd_90[k]
                   + f_3 * pc_z[k] * sod_132[k];

        t_223[k] = f_14 * snd_135[k]
                   + f_3 * pc_x[k] * sod_135[k];

        t_224[k] = f_14 * snd_136[k]
                   + f_3 * pc_x[k] * sod_136[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pb_z, pc_x, pc_y, pc_z, snf0_156, snd_93, \
                         snd_101, snd_137, snf1_156, sod_135, sod_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_14 * snd_137[k]
                   + f_3 * pc_x[k] * sod_137[k];

        t_226[k] = pb_z[k] * snf0_156[k]
                   - f_4 * pc_z[k] * snf1_156[k];

        t_227[k] = f_5 * snd_93[k]
                   + f_3 * pc_z[k] * sod_135[k];

        t_228[k] = f_14 * snd_101[k]
                   + f_3 * pc_y[k] * sod_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, pc_z, snd_95, snd_102, snd_138, \
                         sop0_68, sop0_69, sop1_68, sop1_69, sod_137, \
                         sod_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_5 * snd_95[k]
                   + f_1 * sop0_68[k]
                   - f_2 * sop1_68[k]
                   + f_3 * pc_z[k] * sod_137[k];

        t_230[k] = f_14 * snd_138[k]
                   + f_1 * sop0_69[k]
                   - f_2 * sop1_69[k]
                   + f_3 * pc_x[k] * sod_138[k];

        t_231[k] = f_12 * snd_102[k]
                   + f_3 * pc_y[k] * sod_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_z, snd_96, snd_141, snd_142, \
                         snd_143, sod_138, sod_141, sod_142, sod_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * snd_96[k]
                   + f_3 * pc_z[k] * sod_138[k];

        t_233[k] = f_14 * snd_141[k]
                   + f_3 * pc_x[k] * sod_141[k];

        t_234[k] = f_14 * snd_142[k]
                   + f_3 * pc_x[k] * sod_142[k];

        t_235[k] = f_14 * snd_143[k]
                   + f_3 * pc_x[k] * sod_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_y, pc_z, snd_99, snd_101, snd_105, \
                         snd_107, sop0_70, sop0_71, sop1_70, sop1_71, sod_141, \
                         sod_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_12 * snd_105[k]
                   + f_1 * sop0_70[k]
                   - f_2 * sop1_70[k]
                   + f_3 * pc_y[k] * sod_141[k];

        t_237[k] = f_8 * snd_99[k]
                   + f_3 * pc_z[k] * sod_141[k];

        t_238[k] = f_12 * snd_107[k]
                   + f_3 * pc_y[k] * sod_143[k];

        t_239[k] = f_8 * snd_101[k]
                   + f_1 * sop0_71[k]
                   - f_2 * sop1_71[k]
                   + f_3 * pc_z[k] * sod_143[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, pc_y, pc_z, snd_102, snd_108, \
                         snd_144, snd_147, sop0_72, sop1_72, sod_144, \
                         sod_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_14 * snd_144[k]
                   + f_1 * sop0_72[k]
                   - f_2 * sop1_72[k]
                   + f_3 * pc_x[k] * sod_144[k];

        t_241[k] = f_10 * snd_108[k]
                   + f_3 * pc_y[k] * sod_144[k];

        t_242[k] = f_10 * snd_102[k]
                   + f_3 * pc_z[k] * sod_144[k];

        t_243[k] = f_14 * snd_147[k]
                   + f_3 * pc_x[k] * sod_147[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pc_x, pc_y, pc_z, snd_105, snd_111, \
                         snd_148, snd_149, sop0_73, sop1_73, sod_147, sod_148, \
                         sod_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_14 * snd_148[k]
                   + f_3 * pc_x[k] * sod_148[k];

        t_245[k] = f_14 * snd_149[k]
                   + f_3 * pc_x[k] * sod_149[k];

        t_246[k] = f_10 * snd_111[k]
                   + f_1 * sop0_73[k]
                   - f_2 * sop1_73[k]
                   + f_3 * pc_y[k] * sod_147[k];

        t_247[k] = f_10 * snd_105[k]
                   + f_3 * pc_z[k] * sod_147[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, pc_y, pc_z, snd_107, snd_113, snd_150, \
                         sop0_74, sop0_75, sop1_74, sop1_75, sod_149, \
                         sod_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_10 * snd_113[k]
                   + f_3 * pc_y[k] * sod_149[k];

        t_249[k] = f_10 * snd_107[k]
                   + f_1 * sop0_74[k]
                   - f_2 * sop1_74[k]
                   + f_3 * pc_z[k] * sod_149[k];

        t_250[k] = f_14 * snd_150[k]
                   + f_1 * sop0_75[k]
                   - f_2 * sop1_75[k]
                   + f_3 * pc_x[k] * sod_150[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, pc_z, snd_108, snd_114, \
                         snd_153, snd_154, sod_150, sod_153, sod_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_8 * snd_114[k]
                   + f_3 * pc_y[k] * sod_150[k];

        t_252[k] = f_12 * snd_108[k]
                   + f_3 * pc_z[k] * sod_150[k];

        t_253[k] = f_14 * snd_153[k]
                   + f_3 * pc_x[k] * sod_153[k];

        t_254[k] = f_14 * snd_154[k]
                   + f_3 * pc_x[k] * sod_154[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
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
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *snf0_200 = buffer.data(snf0 + 200);
    const auto *snf0_209 = buffer.data(snf0 + 209);
    const auto *snf0_210 = buffer.data(snf0 + 210);
    const auto *snf0_216 = buffer.data(snf0 + 216);
    const auto *snf0_270 = buffer.data(snf0 + 270);
    const auto *snf0_279 = buffer.data(snf0 + 279);
    const auto *snf0_280 = buffer.data(snf0 + 280);
    const auto *snf0_286 = buffer.data(snf0 + 286);

    const auto *snd_111 = buffer.data(snd + 111);
    const auto *snd_113 = buffer.data(snd + 113);
    const auto *snd_114 = buffer.data(snd + 114);
    const auto *snd_117 = buffer.data(snd + 117);
    const auto *snd_119 = buffer.data(snd + 119);
    const auto *snd_120 = buffer.data(snd + 120);
    const auto *snd_123 = buffer.data(snd + 123);
    const auto *snd_125 = buffer.data(snd + 125);
    const auto *snd_126 = buffer.data(snd + 126);
    const auto *snd_129 = buffer.data(snd + 129);
    const auto *snd_131 = buffer.data(snd + 131);
    const auto *snd_132 = buffer.data(snd + 132);
    const auto *snd_135 = buffer.data(snd + 135);
    const auto *snd_137 = buffer.data(snd + 137);
    const auto *snd_138 = buffer.data(snd + 138);
    const auto *snd_141 = buffer.data(snd + 141);
    const auto *snd_143 = buffer.data(snd + 143);
    const auto *snd_144 = buffer.data(snd + 144);
    const auto *snd_147 = buffer.data(snd + 147);
    const auto *snd_149 = buffer.data(snd + 149);
    const auto *snd_150 = buffer.data(snd + 150);
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
    const auto *snd_225 = buffer.data(snd + 225);
    const auto *snd_226 = buffer.data(snd + 226);
    const auto *snd_227 = buffer.data(snd + 227);

    const auto *snf1_200 = buffer.data(snf1 + 200);
    const auto *snf1_209 = buffer.data(snf1 + 209);
    const auto *snf1_210 = buffer.data(snf1 + 210);
    const auto *snf1_216 = buffer.data(snf1 + 216);
    const auto *snf1_270 = buffer.data(snf1 + 270);
    const auto *snf1_279 = buffer.data(snf1 + 279);
    const auto *snf1_280 = buffer.data(snf1 + 280);
    const auto *snf1_286 = buffer.data(snf1 + 286);

    const auto *sop0_76 = buffer.data(sop0 + 76);
    const auto *sop0_77 = buffer.data(sop0 + 77);
    const auto *sop0_79 = buffer.data(sop0 + 79);
    const auto *sop0_81 = buffer.data(sop0 + 81);
    const auto *sop0_82 = buffer.data(sop0 + 82);
    const auto *sop0_83 = buffer.data(sop0 + 83);
    const auto *sop0_84 = buffer.data(sop0 + 84);
    const auto *sop0_85 = buffer.data(sop0 + 85);
    const auto *sop0_86 = buffer.data(sop0 + 86);
    const auto *sop0_89 = buffer.data(sop0 + 89);
    const auto *sop0_90 = buffer.data(sop0 + 90);
    const auto *sop0_91 = buffer.data(sop0 + 91);
    const auto *sop0_92 = buffer.data(sop0 + 92);
    const auto *sop0_93 = buffer.data(sop0 + 93);
    const auto *sop0_94 = buffer.data(sop0 + 94);
    const auto *sop0_95 = buffer.data(sop0 + 95);
    const auto *sop0_96 = buffer.data(sop0 + 96);
    const auto *sop0_97 = buffer.data(sop0 + 97);
    const auto *sop0_98 = buffer.data(sop0 + 98);
    const auto *sop0_99 = buffer.data(sop0 + 99);
    const auto *sop0_100 = buffer.data(sop0 + 100);
    const auto *sop0_101 = buffer.data(sop0 + 101);
    const auto *sop0_103 = buffer.data(sop0 + 103);
    const auto *sop0_105 = buffer.data(sop0 + 105);
    const auto *sop0_106 = buffer.data(sop0 + 106);
    const auto *sop0_107 = buffer.data(sop0 + 107);
    const auto *sop0_108 = buffer.data(sop0 + 108);
    const auto *sop0_109 = buffer.data(sop0 + 109);
    const auto *sop0_110 = buffer.data(sop0 + 110);
    const auto *sop0_113 = buffer.data(sop0 + 113);

    const auto *sop1_76 = buffer.data(sop1 + 76);
    const auto *sop1_77 = buffer.data(sop1 + 77);
    const auto *sop1_79 = buffer.data(sop1 + 79);
    const auto *sop1_81 = buffer.data(sop1 + 81);
    const auto *sop1_82 = buffer.data(sop1 + 82);
    const auto *sop1_83 = buffer.data(sop1 + 83);
    const auto *sop1_84 = buffer.data(sop1 + 84);
    const auto *sop1_85 = buffer.data(sop1 + 85);
    const auto *sop1_86 = buffer.data(sop1 + 86);
    const auto *sop1_89 = buffer.data(sop1 + 89);
    const auto *sop1_90 = buffer.data(sop1 + 90);
    const auto *sop1_91 = buffer.data(sop1 + 91);
    const auto *sop1_92 = buffer.data(sop1 + 92);
    const auto *sop1_93 = buffer.data(sop1 + 93);
    const auto *sop1_94 = buffer.data(sop1 + 94);
    const auto *sop1_95 = buffer.data(sop1 + 95);
    const auto *sop1_96 = buffer.data(sop1 + 96);
    const auto *sop1_97 = buffer.data(sop1 + 97);
    const auto *sop1_98 = buffer.data(sop1 + 98);
    const auto *sop1_99 = buffer.data(sop1 + 99);
    const auto *sop1_100 = buffer.data(sop1 + 100);
    const auto *sop1_101 = buffer.data(sop1 + 101);
    const auto *sop1_103 = buffer.data(sop1 + 103);
    const auto *sop1_105 = buffer.data(sop1 + 105);
    const auto *sop1_106 = buffer.data(sop1 + 106);
    const auto *sop1_107 = buffer.data(sop1 + 107);
    const auto *sop1_108 = buffer.data(sop1 + 108);
    const auto *sop1_109 = buffer.data(sop1 + 109);
    const auto *sop1_110 = buffer.data(sop1 + 110);
    const auto *sop1_113 = buffer.data(sop1 + 113);

    const auto *sod_153 = buffer.data(sod + 153);
    const auto *sod_155 = buffer.data(sod + 155);
    const auto *sod_156 = buffer.data(sod + 156);
    const auto *sod_159 = buffer.data(sod + 159);
    const auto *sod_160 = buffer.data(sod + 160);
    const auto *sod_161 = buffer.data(sod + 161);
    const auto *sod_162 = buffer.data(sod + 162);
    const auto *sod_165 = buffer.data(sod + 165);
    const auto *sod_166 = buffer.data(sod + 166);
    const auto *sod_167 = buffer.data(sod + 167);
    const auto *sod_168 = buffer.data(sod + 168);
    const auto *sod_171 = buffer.data(sod + 171);
    const auto *sod_172 = buffer.data(sod + 172);
    const auto *sod_173 = buffer.data(sod + 173);
    const auto *sod_174 = buffer.data(sod + 174);
    const auto *sod_177 = buffer.data(sod + 177);
    const auto *sod_178 = buffer.data(sod + 178);
    const auto *sod_179 = buffer.data(sod + 179);
    const auto *sod_180 = buffer.data(sod + 180);
    const auto *sod_183 = buffer.data(sod + 183);
    const auto *sod_184 = buffer.data(sod + 184);
    const auto *sod_185 = buffer.data(sod + 185);
    const auto *sod_186 = buffer.data(sod + 186);
    const auto *sod_189 = buffer.data(sod + 189);
    const auto *sod_190 = buffer.data(sod + 190);
    const auto *sod_191 = buffer.data(sod + 191);
    const auto *sod_192 = buffer.data(sod + 192);
    const auto *sod_195 = buffer.data(sod + 195);
    const auto *sod_196 = buffer.data(sod + 196);
    const auto *sod_197 = buffer.data(sod + 197);
    const auto *sod_198 = buffer.data(sod + 198);
    const auto *sod_201 = buffer.data(sod + 201);
    const auto *sod_202 = buffer.data(sod + 202);
    const auto *sod_203 = buffer.data(sod + 203);
    const auto *sod_204 = buffer.data(sod + 204);
    const auto *sod_207 = buffer.data(sod + 207);
    const auto *sod_208 = buffer.data(sod + 208);
    const auto *sod_209 = buffer.data(sod + 209);
    const auto *sod_210 = buffer.data(sod + 210);
    const auto *sod_213 = buffer.data(sod + 213);
    const auto *sod_214 = buffer.data(sod + 214);
    const auto *sod_215 = buffer.data(sod + 215);
    const auto *sod_216 = buffer.data(sod + 216);
    const auto *sod_219 = buffer.data(sod + 219);
    const auto *sod_220 = buffer.data(sod + 220);
    const auto *sod_221 = buffer.data(sod + 221);
    const auto *sod_222 = buffer.data(sod + 222);
    const auto *sod_225 = buffer.data(sod + 225);
    const auto *sod_226 = buffer.data(sod + 226);
    const auto *sod_227 = buffer.data(sod + 227);

#pragma omp simd aligned(t_255, t_256, t_257, t_258, pc_x, pc_y, pc_z, snd_111, snd_117, \
                         snd_119, snd_155, sop0_76, sop1_76, sod_153, \
                         sod_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_14 * snd_155[k]
                   + f_3 * pc_x[k] * sod_155[k];

        t_256[k] = f_8 * snd_117[k]
                   + f_1 * sop0_76[k]
                   - f_2 * sop1_76[k]
                   + f_3 * pc_y[k] * sod_153[k];

        t_257[k] = f_12 * snd_111[k]
                   + f_3 * pc_z[k] * sod_153[k];

        t_258[k] = f_8 * snd_119[k]
                   + f_3 * pc_y[k] * sod_155[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, pb_y, pc_y, pc_z, snf0_200, snd_113, \
                         snd_114, snd_120, snf1_200, sop0_77, sop1_77, sod_155, \
                         sod_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_12 * snd_113[k]
                   + f_1 * sop0_77[k]
                   - f_2 * sop1_77[k]
                   + f_3 * pc_z[k] * sod_155[k];

        t_260[k] = pb_y[k] * snf0_200[k]
                   - f_4 * pc_y[k] * snf1_200[k];

        t_261[k] = f_5 * snd_120[k]
                   + f_3 * pc_y[k] * sod_156[k];

        t_262[k] = f_14 * snd_114[k]
                   + f_3 * pc_z[k] * sod_156[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, t_266, pc_x, pc_y, snd_123, snd_159, snd_160, \
                         snd_161, sop0_79, sop1_79, sod_159, sod_160, \
                         sod_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_14 * snd_159[k]
                   + f_3 * pc_x[k] * sod_159[k];

        t_264[k] = f_14 * snd_160[k]
                   + f_3 * pc_x[k] * sod_160[k];

        t_265[k] = f_14 * snd_161[k]
                   + f_3 * pc_x[k] * sod_161[k];

        t_266[k] = f_5 * snd_123[k]
                   + f_1 * sop0_79[k]
                   - f_2 * sop1_79[k]
                   + f_3 * pc_y[k] * sod_159[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_y, pc_y, pc_z, snf0_209, snd_117, snd_125, \
                         snf1_209, sod_159, sod_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * snd_117[k]
                   + f_3 * pc_z[k] * sod_159[k];

        t_268[k] = f_5 * snd_125[k]
                   + f_3 * pc_y[k] * sod_161[k];

        t_269[k] = pb_y[k] * snf0_209[k]
                   - f_4 * pc_y[k] * snf1_209[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, pc_z, snd_120, snd_162, \
                         snd_165, sop0_81, sop1_81, sod_162, sod_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_14 * snd_162[k]
                   + f_1 * sop0_81[k]
                   - f_2 * sop1_81[k]
                   + f_3 * pc_x[k] * sod_162[k];

        t_271[k] = f_3 * pc_y[k] * sod_162[k];

        t_272[k] = f_13 * snd_120[k]
                   + f_3 * pc_z[k] * sod_162[k];

        t_273[k] = f_14 * snd_165[k]
                   + f_3 * pc_x[k] * sod_165[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_y, pc_z, snd_123, \
                         snd_166, snd_167, sop0_82, sop1_82, sod_165, sod_166, \
                         sod_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_14 * snd_166[k]
                   + f_3 * pc_x[k] * sod_166[k];

        t_275[k] = f_14 * snd_167[k]
                   + f_3 * pc_x[k] * sod_167[k];

        t_276[k] = f_1 * sop0_82[k]
                   - f_2 * sop1_82[k]
                   + f_3 * pc_y[k] * sod_165[k];

        t_277[k] = f_13 * snd_123[k]
                   + f_3 * pc_z[k] * sod_165[k];

        t_278[k] = f_3 * pc_y[k] * sod_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pc_x, pc_y, pc_z, snd_125, snd_126, \
                         snd_168, sop0_83, sop0_84, sop1_83, sop1_84, sod_167, \
                         sod_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_13 * snd_125[k]
                   + f_1 * sop0_83[k]
                   - f_2 * sop1_83[k]
                   + f_3 * pc_z[k] * sod_167[k];

        t_280[k] = f_12 * snd_168[k]
                   + f_1 * sop0_84[k]
                   - f_2 * sop1_84[k]
                   + f_3 * pc_x[k] * sod_168[k];

        t_281[k] = f_11 * snd_126[k]
                   + f_3 * pc_y[k] * sod_168[k];

        t_282[k] = f_3 * pc_z[k] * sod_168[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pc_x, pc_y, snd_129, snd_171, snd_172, \
                         snd_173, sop0_85, sop1_85, sod_171, sod_172, \
                         sod_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_12 * snd_171[k]
                   + f_3 * pc_x[k] * sod_171[k];

        t_284[k] = f_12 * snd_172[k]
                   + f_3 * pc_x[k] * sod_172[k];

        t_285[k] = f_12 * snd_173[k]
                   + f_3 * pc_x[k] * sod_173[k];

        t_286[k] = f_11 * snd_129[k]
                   + f_1 * sop0_85[k]
                   - f_2 * sop1_85[k]
                   + f_3 * pc_y[k] * sod_171[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pb_z, pc_y, pc_z, snf0_210, snd_131, \
                         snf1_210, sop0_86, sop1_86, sod_171, sod_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_3 * pc_z[k] * sod_171[k];

        t_288[k] = f_11 * snd_131[k]
                   + f_3 * pc_y[k] * sod_173[k];

        t_289[k] = f_1 * sop0_86[k]
                   - f_2 * sop1_86[k]
                   + f_3 * pc_z[k] * sod_173[k];

        t_290[k] = pb_z[k] * snf0_210[k]
                   - f_4 * pc_z[k] * snf1_210[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, pc_y, pc_z, snd_126, snd_132, \
                         snd_177, snd_178, sod_174, sod_177, sod_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_13 * snd_132[k]
                   + f_3 * pc_y[k] * sod_174[k];

        t_292[k] = f_5 * snd_126[k]
                   + f_3 * pc_z[k] * sod_174[k];

        t_293[k] = f_12 * snd_177[k]
                   + f_3 * pc_x[k] * sod_177[k];

        t_294[k] = f_12 * snd_178[k]
                   + f_3 * pc_x[k] * sod_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pb_z, pc_x, pc_y, pc_z, snf0_216, \
                         snd_129, snd_137, snd_179, snf1_216, sod_177, \
                         sod_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_12 * snd_179[k]
                   + f_3 * pc_x[k] * sod_179[k];

        t_296[k] = pb_z[k] * snf0_216[k]
                   - f_4 * pc_z[k] * snf1_216[k];

        t_297[k] = f_5 * snd_129[k]
                   + f_3 * pc_z[k] * sod_177[k];

        t_298[k] = f_13 * snd_137[k]
                   + f_3 * pc_y[k] * sod_179[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pc_x, pc_y, pc_z, snd_131, snd_138, snd_180, \
                         sop0_89, sop0_90, sop1_89, sop1_90, sod_179, \
                         sod_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_5 * snd_131[k]
                   + f_1 * sop0_89[k]
                   - f_2 * sop1_89[k]
                   + f_3 * pc_z[k] * sod_179[k];

        t_300[k] = f_12 * snd_180[k]
                   + f_1 * sop0_90[k]
                   - f_2 * sop1_90[k]
                   + f_3 * pc_x[k] * sod_180[k];

        t_301[k] = f_14 * snd_138[k]
                   + f_3 * pc_y[k] * sod_180[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_z, snd_132, snd_183, snd_184, \
                         snd_185, sod_180, sod_183, sod_184, sod_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_8 * snd_132[k]
                   + f_3 * pc_z[k] * sod_180[k];

        t_303[k] = f_12 * snd_183[k]
                   + f_3 * pc_x[k] * sod_183[k];

        t_304[k] = f_12 * snd_184[k]
                   + f_3 * pc_x[k] * sod_184[k];

        t_305[k] = f_12 * snd_185[k]
                   + f_3 * pc_x[k] * sod_185[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_y, pc_z, snd_135, snd_137, snd_141, \
                         snd_143, sop0_91, sop0_92, sop1_91, sop1_92, sod_183, \
                         sod_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_14 * snd_141[k]
                   + f_1 * sop0_91[k]
                   - f_2 * sop1_91[k]
                   + f_3 * pc_y[k] * sod_183[k];

        t_307[k] = f_8 * snd_135[k]
                   + f_3 * pc_z[k] * sod_183[k];

        t_308[k] = f_14 * snd_143[k]
                   + f_3 * pc_y[k] * sod_185[k];

        t_309[k] = f_8 * snd_137[k]
                   + f_1 * sop0_92[k]
                   - f_2 * sop1_92[k]
                   + f_3 * pc_z[k] * sod_185[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pc_x, pc_y, pc_z, snd_138, snd_144, \
                         snd_186, snd_189, sop0_93, sop1_93, sod_186, \
                         sod_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_12 * snd_186[k]
                   + f_1 * sop0_93[k]
                   - f_2 * sop1_93[k]
                   + f_3 * pc_x[k] * sod_186[k];

        t_311[k] = f_12 * snd_144[k]
                   + f_3 * pc_y[k] * sod_186[k];

        t_312[k] = f_10 * snd_138[k]
                   + f_3 * pc_z[k] * sod_186[k];

        t_313[k] = f_12 * snd_189[k]
                   + f_3 * pc_x[k] * sod_189[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pc_x, pc_y, pc_z, snd_141, snd_147, \
                         snd_190, snd_191, sop0_94, sop1_94, sod_189, sod_190, \
                         sod_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = f_12 * snd_190[k]
                   + f_3 * pc_x[k] * sod_190[k];

        t_315[k] = f_12 * snd_191[k]
                   + f_3 * pc_x[k] * sod_191[k];

        t_316[k] = f_12 * snd_147[k]
                   + f_1 * sop0_94[k]
                   - f_2 * sop1_94[k]
                   + f_3 * pc_y[k] * sod_189[k];

        t_317[k] = f_10 * snd_141[k]
                   + f_3 * pc_z[k] * sod_189[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_x, pc_y, pc_z, snd_143, snd_149, snd_192, \
                         sop0_95, sop0_96, sop1_95, sop1_96, sod_191, \
                         sod_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_12 * snd_149[k]
                   + f_3 * pc_y[k] * sod_191[k];

        t_319[k] = f_10 * snd_143[k]
                   + f_1 * sop0_95[k]
                   - f_2 * sop1_95[k]
                   + f_3 * pc_z[k] * sod_191[k];

        t_320[k] = f_12 * snd_192[k]
                   + f_1 * sop0_96[k]
                   - f_2 * sop1_96[k]
                   + f_3 * pc_x[k] * sod_192[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pc_y, pc_z, snd_144, snd_150, \
                         snd_195, snd_196, sod_192, sod_195, sod_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_10 * snd_150[k]
                   + f_3 * pc_y[k] * sod_192[k];

        t_322[k] = f_12 * snd_144[k]
                   + f_3 * pc_z[k] * sod_192[k];

        t_323[k] = f_12 * snd_195[k]
                   + f_3 * pc_x[k] * sod_195[k];

        t_324[k] = f_12 * snd_196[k]
                   + f_3 * pc_x[k] * sod_196[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, snd_147, snd_153, \
                         snd_155, snd_197, sop0_97, sop1_97, sod_195, \
                         sod_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_12 * snd_197[k]
                   + f_3 * pc_x[k] * sod_197[k];

        t_326[k] = f_10 * snd_153[k]
                   + f_1 * sop0_97[k]
                   - f_2 * sop1_97[k]
                   + f_3 * pc_y[k] * sod_195[k];

        t_327[k] = f_12 * snd_147[k]
                   + f_3 * pc_z[k] * sod_195[k];

        t_328[k] = f_10 * snd_155[k]
                   + f_3 * pc_y[k] * sod_197[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pc_x, pc_y, pc_z, snd_149, snd_156, snd_198, \
                         sop0_98, sop0_99, sop1_98, sop1_99, sod_197, \
                         sod_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_12 * snd_149[k]
                   + f_1 * sop0_98[k]
                   - f_2 * sop1_98[k]
                   + f_3 * pc_z[k] * sod_197[k];

        t_330[k] = f_12 * snd_198[k]
                   + f_1 * sop0_99[k]
                   - f_2 * sop1_99[k]
                   + f_3 * pc_x[k] * sod_198[k];

        t_331[k] = f_8 * snd_156[k]
                   + f_3 * pc_y[k] * sod_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pc_x, pc_z, snd_150, snd_201, snd_202, \
                         snd_203, sod_198, sod_201, sod_202, sod_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_14 * snd_150[k]
                   + f_3 * pc_z[k] * sod_198[k];

        t_333[k] = f_12 * snd_201[k]
                   + f_3 * pc_x[k] * sod_201[k];

        t_334[k] = f_12 * snd_202[k]
                   + f_3 * pc_x[k] * sod_202[k];

        t_335[k] = f_12 * snd_203[k]
                   + f_3 * pc_x[k] * sod_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_y, pc_z, snd_153, snd_155, snd_159, \
                         snd_161, sop0_100, sop0_101, sop1_100, sop1_101, sod_201, \
                         sod_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_8 * snd_159[k]
                   + f_1 * sop0_100[k]
                   - f_2 * sop1_100[k]
                   + f_3 * pc_y[k] * sod_201[k];

        t_337[k] = f_14 * snd_153[k]
                   + f_3 * pc_z[k] * sod_201[k];

        t_338[k] = f_8 * snd_161[k]
                   + f_3 * pc_y[k] * sod_203[k];

        t_339[k] = f_14 * snd_155[k]
                   + f_1 * sop0_101[k]
                   - f_2 * sop1_101[k]
                   + f_3 * pc_z[k] * sod_203[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pb_y, pc_x, pc_y, pc_z, snf0_270, \
                         snd_156, snd_162, snd_207, snf1_270, sod_204, \
                         sod_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pb_y[k] * snf0_270[k]
                   - f_4 * pc_y[k] * snf1_270[k];

        t_341[k] = f_5 * snd_162[k]
                   + f_3 * pc_y[k] * sod_204[k];

        t_342[k] = f_13 * snd_156[k]
                   + f_3 * pc_z[k] * sod_204[k];

        t_343[k] = f_12 * snd_207[k]
                   + f_3 * pc_x[k] * sod_207[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pc_x, pc_y, pc_z, snd_159, snd_165, \
                         snd_208, snd_209, sop0_103, sop1_103, sod_207, sod_208, \
                         sod_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_12 * snd_208[k]
                   + f_3 * pc_x[k] * sod_208[k];

        t_345[k] = f_12 * snd_209[k]
                   + f_3 * pc_x[k] * sod_209[k];

        t_346[k] = f_5 * snd_165[k]
                   + f_1 * sop0_103[k]
                   - f_2 * sop1_103[k]
                   + f_3 * pc_y[k] * sod_207[k];

        t_347[k] = f_13 * snd_159[k]
                   + f_3 * pc_z[k] * sod_207[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pb_y, pc_x, pc_y, snf0_279, snd_167, \
                         snd_210, snf1_279, sop0_105, sop1_105, sod_209, \
                         sod_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_5 * snd_167[k]
                   + f_3 * pc_y[k] * sod_209[k];

        t_349[k] = pb_y[k] * snf0_279[k]
                   - f_4 * pc_y[k] * snf1_279[k];

        t_350[k] = f_12 * snd_210[k]
                   + f_1 * sop0_105[k]
                   - f_2 * sop1_105[k]
                   + f_3 * pc_x[k] * sod_210[k];

        t_351[k] = f_3 * pc_y[k] * sod_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_x, pc_z, snd_162, snd_213, snd_214, \
                         snd_215, sod_210, sod_213, sod_214, sod_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_11 * snd_162[k]
                   + f_3 * pc_z[k] * sod_210[k];

        t_353[k] = f_12 * snd_213[k]
                   + f_3 * pc_x[k] * sod_213[k];

        t_354[k] = f_12 * snd_214[k]
                   + f_3 * pc_x[k] * sod_214[k];

        t_355[k] = f_12 * snd_215[k]
                   + f_3 * pc_x[k] * sod_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pc_y, pc_z, snd_165, snd_167, sop0_106, \
                         sop0_107, sop1_106, sop1_107, sod_213, \
                         sod_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_1 * sop0_106[k]
                   - f_2 * sop1_106[k]
                   + f_3 * pc_y[k] * sod_213[k];

        t_357[k] = f_11 * snd_165[k]
                   + f_3 * pc_z[k] * sod_213[k];

        t_358[k] = f_3 * pc_y[k] * sod_215[k];

        t_359[k] = f_11 * snd_167[k]
                   + f_1 * sop0_107[k]
                   - f_2 * sop1_107[k]
                   + f_3 * pc_z[k] * sod_215[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, pc_x, pc_y, pc_z, snd_168, snd_216, \
                         snd_219, sop0_108, sop1_108, sod_216, \
                         sod_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_10 * snd_216[k]
                   + f_1 * sop0_108[k]
                   - f_2 * sop1_108[k]
                   + f_3 * pc_x[k] * sod_216[k];

        t_361[k] = f_9 * snd_168[k]
                   + f_3 * pc_y[k] * sod_216[k];

        t_362[k] = f_3 * pc_z[k] * sod_216[k];

        t_363[k] = f_10 * snd_219[k]
                   + f_3 * pc_x[k] * sod_219[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, snd_171, snd_220, \
                         snd_221, sop0_109, sop1_109, sod_219, sod_220, \
                         sod_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_10 * snd_220[k]
                   + f_3 * pc_x[k] * sod_220[k];

        t_365[k] = f_10 * snd_221[k]
                   + f_3 * pc_x[k] * sod_221[k];

        t_366[k] = f_9 * snd_171[k]
                   + f_1 * sop0_109[k]
                   - f_2 * sop1_109[k]
                   + f_3 * pc_y[k] * sod_219[k];

        t_367[k] = f_3 * pc_z[k] * sod_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pb_z, pc_y, pc_z, snf0_280, snd_173, \
                         snd_174, snf1_280, sop0_110, sop1_110, sod_221, \
                         sod_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_9 * snd_173[k]
                   + f_3 * pc_y[k] * sod_221[k];

        t_369[k] = f_1 * sop0_110[k]
                   - f_2 * sop1_110[k]
                   + f_3 * pc_z[k] * sod_221[k];

        t_370[k] = pb_z[k] * snf0_280[k]
                   - f_4 * pc_z[k] * snf1_280[k];

        t_371[k] = f_11 * snd_174[k]
                   + f_3 * pc_y[k] * sod_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, pc_z, snd_168, snd_225, snd_226, \
                         snd_227, sod_222, sod_225, sod_226, sod_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_5 * snd_168[k]
                   + f_3 * pc_z[k] * sod_222[k];

        t_373[k] = f_10 * snd_225[k]
                   + f_3 * pc_x[k] * sod_225[k];

        t_374[k] = f_10 * snd_226[k]
                   + f_3 * pc_x[k] * sod_226[k];

        t_375[k] = f_10 * snd_227[k]
                   + f_3 * pc_x[k] * sod_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pb_z, pc_y, pc_z, snf0_286, snd_171, \
                         snd_173, snd_179, snf1_286, sop0_113, sop1_113, sod_225, \
                         sod_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pb_z[k] * snf0_286[k]
                   - f_4 * pc_z[k] * snf1_286[k];

        t_377[k] = f_5 * snd_171[k]
                   + f_3 * pc_z[k] * sod_225[k];

        t_378[k] = f_11 * snd_179[k]
                   + f_3 * pc_y[k] * sod_227[k];

        t_379[k] = f_5 * snd_173[k]
                   + f_1 * sop0_113[k]
                   - f_2 * sop1_113[k]
                   + f_3 * pc_z[k] * sod_227[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
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
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snf0_350 = buffer.data(snf0 + 350);
    const auto *snf0_359 = buffer.data(snf0 + 359);
    const auto *snf0_360 = buffer.data(snf0 + 360);
    const auto *snf0_366 = buffer.data(snf0 + 366);

    const auto *snd_174 = buffer.data(snd + 174);
    const auto *snd_177 = buffer.data(snd + 177);
    const auto *snd_179 = buffer.data(snd + 179);
    const auto *snd_180 = buffer.data(snd + 180);
    const auto *snd_183 = buffer.data(snd + 183);
    const auto *snd_185 = buffer.data(snd + 185);
    const auto *snd_186 = buffer.data(snd + 186);
    const auto *snd_189 = buffer.data(snd + 189);
    const auto *snd_191 = buffer.data(snd + 191);
    const auto *snd_192 = buffer.data(snd + 192);
    const auto *snd_195 = buffer.data(snd + 195);
    const auto *snd_197 = buffer.data(snd + 197);
    const auto *snd_198 = buffer.data(snd + 198);
    const auto *snd_201 = buffer.data(snd + 201);
    const auto *snd_203 = buffer.data(snd + 203);
    const auto *snd_204 = buffer.data(snd + 204);
    const auto *snd_207 = buffer.data(snd + 207);
    const auto *snd_209 = buffer.data(snd + 209);
    const auto *snd_210 = buffer.data(snd + 210);
    const auto *snd_213 = buffer.data(snd + 213);
    const auto *snd_215 = buffer.data(snd + 215);
    const auto *snd_216 = buffer.data(snd + 216);
    const auto *snd_219 = buffer.data(snd + 219);
    const auto *snd_221 = buffer.data(snd + 221);
    const auto *snd_222 = buffer.data(snd + 222);
    const auto *snd_225 = buffer.data(snd + 225);
    const auto *snd_227 = buffer.data(snd + 227);
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

    const auto *snf1_350 = buffer.data(snf1 + 350);
    const auto *snf1_359 = buffer.data(snf1 + 359);
    const auto *snf1_360 = buffer.data(snf1 + 360);
    const auto *snf1_366 = buffer.data(snf1 + 366);

    const auto *sop0_114 = buffer.data(sop0 + 114);
    const auto *sop0_115 = buffer.data(sop0 + 115);
    const auto *sop0_116 = buffer.data(sop0 + 116);
    const auto *sop0_117 = buffer.data(sop0 + 117);
    const auto *sop0_118 = buffer.data(sop0 + 118);
    const auto *sop0_119 = buffer.data(sop0 + 119);
    const auto *sop0_120 = buffer.data(sop0 + 120);
    const auto *sop0_121 = buffer.data(sop0 + 121);
    const auto *sop0_122 = buffer.data(sop0 + 122);
    const auto *sop0_123 = buffer.data(sop0 + 123);
    const auto *sop0_124 = buffer.data(sop0 + 124);
    const auto *sop0_125 = buffer.data(sop0 + 125);
    const auto *sop0_126 = buffer.data(sop0 + 126);
    const auto *sop0_127 = buffer.data(sop0 + 127);
    const auto *sop0_128 = buffer.data(sop0 + 128);
    const auto *sop0_130 = buffer.data(sop0 + 130);
    const auto *sop0_132 = buffer.data(sop0 + 132);
    const auto *sop0_133 = buffer.data(sop0 + 133);
    const auto *sop0_134 = buffer.data(sop0 + 134);
    const auto *sop0_135 = buffer.data(sop0 + 135);
    const auto *sop0_136 = buffer.data(sop0 + 136);
    const auto *sop0_137 = buffer.data(sop0 + 137);
    const auto *sop0_140 = buffer.data(sop0 + 140);
    const auto *sop0_141 = buffer.data(sop0 + 141);
    const auto *sop0_142 = buffer.data(sop0 + 142);
    const auto *sop0_143 = buffer.data(sop0 + 143);
    const auto *sop0_144 = buffer.data(sop0 + 144);
    const auto *sop0_145 = buffer.data(sop0 + 145);
    const auto *sop0_146 = buffer.data(sop0 + 146);
    const auto *sop0_147 = buffer.data(sop0 + 147);
    const auto *sop0_148 = buffer.data(sop0 + 148);
    const auto *sop0_149 = buffer.data(sop0 + 149);
    const auto *sop0_150 = buffer.data(sop0 + 150);

    const auto *sop1_114 = buffer.data(sop1 + 114);
    const auto *sop1_115 = buffer.data(sop1 + 115);
    const auto *sop1_116 = buffer.data(sop1 + 116);
    const auto *sop1_117 = buffer.data(sop1 + 117);
    const auto *sop1_118 = buffer.data(sop1 + 118);
    const auto *sop1_119 = buffer.data(sop1 + 119);
    const auto *sop1_120 = buffer.data(sop1 + 120);
    const auto *sop1_121 = buffer.data(sop1 + 121);
    const auto *sop1_122 = buffer.data(sop1 + 122);
    const auto *sop1_123 = buffer.data(sop1 + 123);
    const auto *sop1_124 = buffer.data(sop1 + 124);
    const auto *sop1_125 = buffer.data(sop1 + 125);
    const auto *sop1_126 = buffer.data(sop1 + 126);
    const auto *sop1_127 = buffer.data(sop1 + 127);
    const auto *sop1_128 = buffer.data(sop1 + 128);
    const auto *sop1_130 = buffer.data(sop1 + 130);
    const auto *sop1_132 = buffer.data(sop1 + 132);
    const auto *sop1_133 = buffer.data(sop1 + 133);
    const auto *sop1_134 = buffer.data(sop1 + 134);
    const auto *sop1_135 = buffer.data(sop1 + 135);
    const auto *sop1_136 = buffer.data(sop1 + 136);
    const auto *sop1_137 = buffer.data(sop1 + 137);
    const auto *sop1_140 = buffer.data(sop1 + 140);
    const auto *sop1_141 = buffer.data(sop1 + 141);
    const auto *sop1_142 = buffer.data(sop1 + 142);
    const auto *sop1_143 = buffer.data(sop1 + 143);
    const auto *sop1_144 = buffer.data(sop1 + 144);
    const auto *sop1_145 = buffer.data(sop1 + 145);
    const auto *sop1_146 = buffer.data(sop1 + 146);
    const auto *sop1_147 = buffer.data(sop1 + 147);
    const auto *sop1_148 = buffer.data(sop1 + 148);
    const auto *sop1_149 = buffer.data(sop1 + 149);
    const auto *sop1_150 = buffer.data(sop1 + 150);

    const auto *sod_228 = buffer.data(sod + 228);
    const auto *sod_231 = buffer.data(sod + 231);
    const auto *sod_232 = buffer.data(sod + 232);
    const auto *sod_233 = buffer.data(sod + 233);
    const auto *sod_234 = buffer.data(sod + 234);
    const auto *sod_237 = buffer.data(sod + 237);
    const auto *sod_238 = buffer.data(sod + 238);
    const auto *sod_239 = buffer.data(sod + 239);
    const auto *sod_240 = buffer.data(sod + 240);
    const auto *sod_243 = buffer.data(sod + 243);
    const auto *sod_244 = buffer.data(sod + 244);
    const auto *sod_245 = buffer.data(sod + 245);
    const auto *sod_246 = buffer.data(sod + 246);
    const auto *sod_249 = buffer.data(sod + 249);
    const auto *sod_250 = buffer.data(sod + 250);
    const auto *sod_251 = buffer.data(sod + 251);
    const auto *sod_252 = buffer.data(sod + 252);
    const auto *sod_255 = buffer.data(sod + 255);
    const auto *sod_256 = buffer.data(sod + 256);
    const auto *sod_257 = buffer.data(sod + 257);
    const auto *sod_258 = buffer.data(sod + 258);
    const auto *sod_261 = buffer.data(sod + 261);
    const auto *sod_262 = buffer.data(sod + 262);
    const auto *sod_263 = buffer.data(sod + 263);
    const auto *sod_264 = buffer.data(sod + 264);
    const auto *sod_267 = buffer.data(sod + 267);
    const auto *sod_268 = buffer.data(sod + 268);
    const auto *sod_269 = buffer.data(sod + 269);
    const auto *sod_270 = buffer.data(sod + 270);
    const auto *sod_273 = buffer.data(sod + 273);
    const auto *sod_274 = buffer.data(sod + 274);
    const auto *sod_275 = buffer.data(sod + 275);
    const auto *sod_276 = buffer.data(sod + 276);
    const auto *sod_279 = buffer.data(sod + 279);
    const auto *sod_280 = buffer.data(sod + 280);
    const auto *sod_281 = buffer.data(sod + 281);
    const auto *sod_282 = buffer.data(sod + 282);
    const auto *sod_285 = buffer.data(sod + 285);
    const auto *sod_286 = buffer.data(sod + 286);
    const auto *sod_287 = buffer.data(sod + 287);
    const auto *sod_288 = buffer.data(sod + 288);
    const auto *sod_291 = buffer.data(sod + 291);
    const auto *sod_292 = buffer.data(sod + 292);
    const auto *sod_293 = buffer.data(sod + 293);
    const auto *sod_294 = buffer.data(sod + 294);
    const auto *sod_297 = buffer.data(sod + 297);
    const auto *sod_298 = buffer.data(sod + 298);
    const auto *sod_299 = buffer.data(sod + 299);
    const auto *sod_300 = buffer.data(sod + 300);

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, snd_174, snd_180, \
                         snd_228, snd_231, sop0_114, sop1_114, sod_228, \
                         sod_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_10 * snd_228[k]
                   + f_1 * sop0_114[k]
                   - f_2 * sop1_114[k]
                   + f_3 * pc_x[k] * sod_228[k];

        t_381[k] = f_13 * snd_180[k]
                   + f_3 * pc_y[k] * sod_228[k];

        t_382[k] = f_8 * snd_174[k]
                   + f_3 * pc_z[k] * sod_228[k];

        t_383[k] = f_10 * snd_231[k]
                   + f_3 * pc_x[k] * sod_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pc_x, pc_y, pc_z, snd_177, snd_183, \
                         snd_232, snd_233, sop0_115, sop1_115, sod_231, sod_232, \
                         sod_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_10 * snd_232[k]
                   + f_3 * pc_x[k] * sod_232[k];

        t_385[k] = f_10 * snd_233[k]
                   + f_3 * pc_x[k] * sod_233[k];

        t_386[k] = f_13 * snd_183[k]
                   + f_1 * sop0_115[k]
                   - f_2 * sop1_115[k]
                   + f_3 * pc_y[k] * sod_231[k];

        t_387[k] = f_8 * snd_177[k]
                   + f_3 * pc_z[k] * sod_231[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_x, pc_y, pc_z, snd_179, snd_185, snd_234, \
                         sop0_116, sop0_117, sop1_116, sop1_117, sod_233, \
                         sod_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_13 * snd_185[k]
                   + f_3 * pc_y[k] * sod_233[k];

        t_389[k] = f_8 * snd_179[k]
                   + f_1 * sop0_116[k]
                   - f_2 * sop1_116[k]
                   + f_3 * pc_z[k] * sod_233[k];

        t_390[k] = f_10 * snd_234[k]
                   + f_1 * sop0_117[k]
                   - f_2 * sop1_117[k]
                   + f_3 * pc_x[k] * sod_234[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, snd_180, snd_186, \
                         snd_237, snd_238, sod_234, sod_237, sod_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_14 * snd_186[k]
                   + f_3 * pc_y[k] * sod_234[k];

        t_392[k] = f_10 * snd_180[k]
                   + f_3 * pc_z[k] * sod_234[k];

        t_393[k] = f_10 * snd_237[k]
                   + f_3 * pc_x[k] * sod_237[k];

        t_394[k] = f_10 * snd_238[k]
                   + f_3 * pc_x[k] * sod_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, snd_183, snd_189, \
                         snd_191, snd_239, sop0_118, sop1_118, sod_237, \
                         sod_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_10 * snd_239[k]
                   + f_3 * pc_x[k] * sod_239[k];

        t_396[k] = f_14 * snd_189[k]
                   + f_1 * sop0_118[k]
                   - f_2 * sop1_118[k]
                   + f_3 * pc_y[k] * sod_237[k];

        t_397[k] = f_10 * snd_183[k]
                   + f_3 * pc_z[k] * sod_237[k];

        t_398[k] = f_14 * snd_191[k]
                   + f_3 * pc_y[k] * sod_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, snd_185, snd_192, snd_240, \
                         sop0_119, sop0_120, sop1_119, sop1_120, sod_239, \
                         sod_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_10 * snd_185[k]
                   + f_1 * sop0_119[k]
                   - f_2 * sop1_119[k]
                   + f_3 * pc_z[k] * sod_239[k];

        t_400[k] = f_10 * snd_240[k]
                   + f_1 * sop0_120[k]
                   - f_2 * sop1_120[k]
                   + f_3 * pc_x[k] * sod_240[k];

        t_401[k] = f_12 * snd_192[k]
                   + f_3 * pc_y[k] * sod_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, snd_186, snd_243, snd_244, \
                         snd_245, sod_240, sod_243, sod_244, sod_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_12 * snd_186[k]
                   + f_3 * pc_z[k] * sod_240[k];

        t_403[k] = f_10 * snd_243[k]
                   + f_3 * pc_x[k] * sod_243[k];

        t_404[k] = f_10 * snd_244[k]
                   + f_3 * pc_x[k] * sod_244[k];

        t_405[k] = f_10 * snd_245[k]
                   + f_3 * pc_x[k] * sod_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pc_y, pc_z, snd_189, snd_191, snd_195, \
                         snd_197, sop0_121, sop0_122, sop1_121, sop1_122, sod_243, \
                         sod_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_12 * snd_195[k]
                   + f_1 * sop0_121[k]
                   - f_2 * sop1_121[k]
                   + f_3 * pc_y[k] * sod_243[k];

        t_407[k] = f_12 * snd_189[k]
                   + f_3 * pc_z[k] * sod_243[k];

        t_408[k] = f_12 * snd_197[k]
                   + f_3 * pc_y[k] * sod_245[k];

        t_409[k] = f_12 * snd_191[k]
                   + f_1 * sop0_122[k]
                   - f_2 * sop1_122[k]
                   + f_3 * pc_z[k] * sod_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pc_z, snd_192, snd_198, \
                         snd_246, snd_249, sop0_123, sop1_123, sod_246, \
                         sod_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_10 * snd_246[k]
                   + f_1 * sop0_123[k]
                   - f_2 * sop1_123[k]
                   + f_3 * pc_x[k] * sod_246[k];

        t_411[k] = f_10 * snd_198[k]
                   + f_3 * pc_y[k] * sod_246[k];

        t_412[k] = f_14 * snd_192[k]
                   + f_3 * pc_z[k] * sod_246[k];

        t_413[k] = f_10 * snd_249[k]
                   + f_3 * pc_x[k] * sod_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, snd_195, snd_201, \
                         snd_250, snd_251, sop0_124, sop1_124, sod_249, sod_250, \
                         sod_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_10 * snd_250[k]
                   + f_3 * pc_x[k] * sod_250[k];

        t_415[k] = f_10 * snd_251[k]
                   + f_3 * pc_x[k] * sod_251[k];

        t_416[k] = f_10 * snd_201[k]
                   + f_1 * sop0_124[k]
                   - f_2 * sop1_124[k]
                   + f_3 * pc_y[k] * sod_249[k];

        t_417[k] = f_14 * snd_195[k]
                   + f_3 * pc_z[k] * sod_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pc_x, pc_y, pc_z, snd_197, snd_203, snd_252, \
                         sop0_125, sop0_126, sop1_125, sop1_126, sod_251, \
                         sod_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_10 * snd_203[k]
                   + f_3 * pc_y[k] * sod_251[k];

        t_419[k] = f_14 * snd_197[k]
                   + f_1 * sop0_125[k]
                   - f_2 * sop1_125[k]
                   + f_3 * pc_z[k] * sod_251[k];

        t_420[k] = f_10 * snd_252[k]
                   + f_1 * sop0_126[k]
                   - f_2 * sop1_126[k]
                   + f_3 * pc_x[k] * sod_252[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, snd_198, snd_204, \
                         snd_255, snd_256, sod_252, sod_255, sod_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_8 * snd_204[k]
                   + f_3 * pc_y[k] * sod_252[k];

        t_422[k] = f_13 * snd_198[k]
                   + f_3 * pc_z[k] * sod_252[k];

        t_423[k] = f_10 * snd_255[k]
                   + f_3 * pc_x[k] * sod_255[k];

        t_424[k] = f_10 * snd_256[k]
                   + f_3 * pc_x[k] * sod_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, snd_201, snd_207, \
                         snd_209, snd_257, sop0_127, sop1_127, sod_255, \
                         sod_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_10 * snd_257[k]
                   + f_3 * pc_x[k] * sod_257[k];

        t_426[k] = f_8 * snd_207[k]
                   + f_1 * sop0_127[k]
                   - f_2 * sop1_127[k]
                   + f_3 * pc_y[k] * sod_255[k];

        t_427[k] = f_13 * snd_201[k]
                   + f_3 * pc_z[k] * sod_255[k];

        t_428[k] = f_8 * snd_209[k]
                   + f_3 * pc_y[k] * sod_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pb_y, pc_y, pc_z, snf0_350, snd_203, \
                         snd_204, snd_210, snf1_350, sop0_128, sop1_128, sod_257, \
                         sod_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_13 * snd_203[k]
                   + f_1 * sop0_128[k]
                   - f_2 * sop1_128[k]
                   + f_3 * pc_z[k] * sod_257[k];

        t_430[k] = pb_y[k] * snf0_350[k]
                   - f_4 * pc_y[k] * snf1_350[k];

        t_431[k] = f_5 * snd_210[k]
                   + f_3 * pc_y[k] * sod_258[k];

        t_432[k] = f_11 * snd_204[k]
                   + f_3 * pc_z[k] * sod_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, snd_213, snd_261, snd_262, \
                         snd_263, sop0_130, sop1_130, sod_261, sod_262, \
                         sod_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_10 * snd_261[k]
                   + f_3 * pc_x[k] * sod_261[k];

        t_434[k] = f_10 * snd_262[k]
                   + f_3 * pc_x[k] * sod_262[k];

        t_435[k] = f_10 * snd_263[k]
                   + f_3 * pc_x[k] * sod_263[k];

        t_436[k] = f_5 * snd_213[k]
                   + f_1 * sop0_130[k]
                   - f_2 * sop1_130[k]
                   + f_3 * pc_y[k] * sod_261[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pb_y, pc_y, pc_z, snf0_359, snd_207, snd_215, \
                         snf1_359, sod_261, sod_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_11 * snd_207[k]
                   + f_3 * pc_z[k] * sod_261[k];

        t_438[k] = f_5 * snd_215[k]
                   + f_3 * pc_y[k] * sod_263[k];

        t_439[k] = pb_y[k] * snf0_359[k]
                   - f_4 * pc_y[k] * snf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, pc_y, pc_z, snd_210, snd_264, \
                         snd_267, sop0_132, sop1_132, sod_264, \
                         sod_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_10 * snd_264[k]
                   + f_1 * sop0_132[k]
                   - f_2 * sop1_132[k]
                   + f_3 * pc_x[k] * sod_264[k];

        t_441[k] = f_3 * pc_y[k] * sod_264[k];

        t_442[k] = f_9 * snd_210[k]
                   + f_3 * pc_z[k] * sod_264[k];

        t_443[k] = f_10 * snd_267[k]
                   + f_3 * pc_x[k] * sod_267[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, t_448, pc_x, pc_y, pc_z, snd_213, \
                         snd_268, snd_269, sop0_133, sop1_133, sod_267, sod_268, \
                         sod_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_10 * snd_268[k]
                   + f_3 * pc_x[k] * sod_268[k];

        t_445[k] = f_10 * snd_269[k]
                   + f_3 * pc_x[k] * sod_269[k];

        t_446[k] = f_1 * sop0_133[k]
                   - f_2 * sop1_133[k]
                   + f_3 * pc_y[k] * sod_267[k];

        t_447[k] = f_9 * snd_213[k]
                   + f_3 * pc_z[k] * sod_267[k];

        t_448[k] = f_3 * pc_y[k] * sod_269[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, snd_215, snd_216, \
                         snd_270, sop0_134, sop0_135, sop1_134, sop1_135, sod_269, \
                         sod_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_9 * snd_215[k]
                   + f_1 * sop0_134[k]
                   - f_2 * sop1_134[k]
                   + f_3 * pc_z[k] * sod_269[k];

        t_450[k] = f_8 * snd_270[k]
                   + f_1 * sop0_135[k]
                   - f_2 * sop1_135[k]
                   + f_3 * pc_x[k] * sod_270[k];

        t_451[k] = f_7 * snd_216[k]
                   + f_3 * pc_y[k] * sod_270[k];

        t_452[k] = f_3 * pc_z[k] * sod_270[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, snd_219, snd_273, snd_274, \
                         snd_275, sop0_136, sop1_136, sod_273, sod_274, \
                         sod_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_8 * snd_273[k]
                   + f_3 * pc_x[k] * sod_273[k];

        t_454[k] = f_8 * snd_274[k]
                   + f_3 * pc_x[k] * sod_274[k];

        t_455[k] = f_8 * snd_275[k]
                   + f_3 * pc_x[k] * sod_275[k];

        t_456[k] = f_7 * snd_219[k]
                   + f_1 * sop0_136[k]
                   - f_2 * sop1_136[k]
                   + f_3 * pc_y[k] * sod_273[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pb_z, pc_y, pc_z, snf0_360, snd_221, \
                         snf1_360, sop0_137, sop1_137, sod_273, \
                         sod_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * sod_273[k];

        t_458[k] = f_7 * snd_221[k]
                   + f_3 * pc_y[k] * sod_275[k];

        t_459[k] = f_1 * sop0_137[k]
                   - f_2 * sop1_137[k]
                   + f_3 * pc_z[k] * sod_275[k];

        t_460[k] = pb_z[k] * snf0_360[k]
                   - f_4 * pc_z[k] * snf1_360[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pc_x, pc_y, pc_z, snd_216, snd_222, \
                         snd_279, snd_280, sod_276, sod_279, sod_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_9 * snd_222[k]
                   + f_3 * pc_y[k] * sod_276[k];

        t_462[k] = f_5 * snd_216[k]
                   + f_3 * pc_z[k] * sod_276[k];

        t_463[k] = f_8 * snd_279[k]
                   + f_3 * pc_x[k] * sod_279[k];

        t_464[k] = f_8 * snd_280[k]
                   + f_3 * pc_x[k] * sod_280[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pb_z, pc_x, pc_y, pc_z, snf0_366, \
                         snd_219, snd_227, snd_281, snf1_366, sod_279, \
                         sod_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_8 * snd_281[k]
                   + f_3 * pc_x[k] * sod_281[k];

        t_466[k] = pb_z[k] * snf0_366[k]
                   - f_4 * pc_z[k] * snf1_366[k];

        t_467[k] = f_5 * snd_219[k]
                   + f_3 * pc_z[k] * sod_279[k];

        t_468[k] = f_9 * snd_227[k]
                   + f_3 * pc_y[k] * sod_281[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, pc_x, pc_y, pc_z, snd_221, snd_228, snd_282, \
                         sop0_140, sop0_141, sop1_140, sop1_141, sod_281, \
                         sod_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_5 * snd_221[k]
                   + f_1 * sop0_140[k]
                   - f_2 * sop1_140[k]
                   + f_3 * pc_z[k] * sod_281[k];

        t_470[k] = f_8 * snd_282[k]
                   + f_1 * sop0_141[k]
                   - f_2 * sop1_141[k]
                   + f_3 * pc_x[k] * sod_282[k];

        t_471[k] = f_11 * snd_228[k]
                   + f_3 * pc_y[k] * sod_282[k];
    }

#pragma omp simd aligned(t_472, t_473, t_474, t_475, pc_x, pc_z, snd_222, snd_285, snd_286, \
                         snd_287, sod_282, sod_285, sod_286, sod_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_472[k] = f_8 * snd_222[k]
                   + f_3 * pc_z[k] * sod_282[k];

        t_473[k] = f_8 * snd_285[k]
                   + f_3 * pc_x[k] * sod_285[k];

        t_474[k] = f_8 * snd_286[k]
                   + f_3 * pc_x[k] * sod_286[k];

        t_475[k] = f_8 * snd_287[k]
                   + f_3 * pc_x[k] * sod_287[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, t_479, pc_y, pc_z, snd_225, snd_227, snd_231, \
                         snd_233, sop0_142, sop0_143, sop1_142, sop1_143, sod_285, \
                         sod_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_11 * snd_231[k]
                   + f_1 * sop0_142[k]
                   - f_2 * sop1_142[k]
                   + f_3 * pc_y[k] * sod_285[k];

        t_477[k] = f_8 * snd_225[k]
                   + f_3 * pc_z[k] * sod_285[k];

        t_478[k] = f_11 * snd_233[k]
                   + f_3 * pc_y[k] * sod_287[k];

        t_479[k] = f_8 * snd_227[k]
                   + f_1 * sop0_143[k]
                   - f_2 * sop1_143[k]
                   + f_3 * pc_z[k] * sod_287[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pc_x, pc_y, pc_z, snd_228, snd_234, \
                         snd_288, snd_291, sop0_144, sop1_144, sod_288, \
                         sod_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_8 * snd_288[k]
                   + f_1 * sop0_144[k]
                   - f_2 * sop1_144[k]
                   + f_3 * pc_x[k] * sod_288[k];

        t_481[k] = f_13 * snd_234[k]
                   + f_3 * pc_y[k] * sod_288[k];

        t_482[k] = f_10 * snd_228[k]
                   + f_3 * pc_z[k] * sod_288[k];

        t_483[k] = f_8 * snd_291[k]
                   + f_3 * pc_x[k] * sod_291[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, pc_x, pc_y, pc_z, snd_231, snd_237, \
                         snd_292, snd_293, sop0_145, sop1_145, sod_291, sod_292, \
                         sod_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_8 * snd_292[k]
                   + f_3 * pc_x[k] * sod_292[k];

        t_485[k] = f_8 * snd_293[k]
                   + f_3 * pc_x[k] * sod_293[k];

        t_486[k] = f_13 * snd_237[k]
                   + f_1 * sop0_145[k]
                   - f_2 * sop1_145[k]
                   + f_3 * pc_y[k] * sod_291[k];

        t_487[k] = f_10 * snd_231[k]
                   + f_3 * pc_z[k] * sod_291[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, pc_z, snd_233, snd_239, snd_294, \
                         sop0_146, sop0_147, sop1_146, sop1_147, sod_293, \
                         sod_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_13 * snd_239[k]
                   + f_3 * pc_y[k] * sod_293[k];

        t_489[k] = f_10 * snd_233[k]
                   + f_1 * sop0_146[k]
                   - f_2 * sop1_146[k]
                   + f_3 * pc_z[k] * sod_293[k];

        t_490[k] = f_8 * snd_294[k]
                   + f_1 * sop0_147[k]
                   - f_2 * sop1_147[k]
                   + f_3 * pc_x[k] * sod_294[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, pc_z, snd_234, snd_240, \
                         snd_297, snd_298, sod_294, sod_297, sod_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_14 * snd_240[k]
                   + f_3 * pc_y[k] * sod_294[k];

        t_492[k] = f_12 * snd_234[k]
                   + f_3 * pc_z[k] * sod_294[k];

        t_493[k] = f_8 * snd_297[k]
                   + f_3 * pc_x[k] * sod_297[k];

        t_494[k] = f_8 * snd_298[k]
                   + f_3 * pc_x[k] * sod_298[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, pc_z, snd_237, snd_243, \
                         snd_245, snd_299, sop0_148, sop1_148, sod_297, \
                         sod_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_8 * snd_299[k]
                   + f_3 * pc_x[k] * sod_299[k];

        t_496[k] = f_14 * snd_243[k]
                   + f_1 * sop0_148[k]
                   - f_2 * sop1_148[k]
                   + f_3 * pc_y[k] * sod_297[k];

        t_497[k] = f_12 * snd_237[k]
                   + f_3 * pc_z[k] * sod_297[k];

        t_498[k] = f_14 * snd_245[k]
                   + f_3 * pc_y[k] * sod_299[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_x, pc_y, pc_z, snd_239, snd_246, snd_300, \
                         sop0_149, sop0_150, sop1_149, sop1_150, sod_299, \
                         sod_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_12 * snd_239[k]
                   + f_1 * sop0_149[k]
                   - f_2 * sop1_149[k]
                   + f_3 * pc_z[k] * sod_299[k];

        t_500[k] = f_8 * snd_300[k]
                   + f_1 * sop0_150[k]
                   - f_2 * sop1_150[k]
                   + f_3 * pc_x[k] * sod_300[k];

        t_501[k] = f_12 * snd_246[k]
                   + f_3 * pc_y[k] * sod_300[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
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
    const auto f_6 = 5.0 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

    auto *t_502 = buffer.data(target + 502);
    auto *t_503 = buffer.data(target + 503);
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snf0_440 = buffer.data(snf0 + 440);
    const auto *snf0_449 = buffer.data(snf0 + 449);
    const auto *snf0_450 = buffer.data(snf0 + 450);
    const auto *snf0_550 = buffer.data(snf0 + 550);
    const auto *snf0_556 = buffer.data(snf0 + 556);
    const auto *snf0_559 = buffer.data(snf0 + 559);
    const auto *snf0_566 = buffer.data(snf0 + 566);
    const auto *snf0_569 = buffer.data(snf0 + 569);
    const auto *snf0_570 = buffer.data(snf0 + 570);
    const auto *snf0_576 = buffer.data(snf0 + 576);
    const auto *snf0_579 = buffer.data(snf0 + 579);
    const auto *snf0_580 = buffer.data(snf0 + 580);
    const auto *snf0_586 = buffer.data(snf0 + 586);
    const auto *snf0_589 = buffer.data(snf0 + 589);
    const auto *snf0_590 = buffer.data(snf0 + 590);
    const auto *snf0_596 = buffer.data(snf0 + 596);
    const auto *snf0_599 = buffer.data(snf0 + 599);
    const auto *snf0_600 = buffer.data(snf0 + 600);
    const auto *snf0_606 = buffer.data(snf0 + 606);
    const auto *snf0_609 = buffer.data(snf0 + 609);
    const auto *snf0_610 = buffer.data(snf0 + 610);
    const auto *snf0_616 = buffer.data(snf0 + 616);
    const auto *snf0_619 = buffer.data(snf0 + 619);
    const auto *snf0_620 = buffer.data(snf0 + 620);
    const auto *snf0_626 = buffer.data(snf0 + 626);

    const auto *snd_240 = buffer.data(snd + 240);
    const auto *snd_243 = buffer.data(snd + 243);
    const auto *snd_245 = buffer.data(snd + 245);
    const auto *snd_246 = buffer.data(snd + 246);
    const auto *snd_249 = buffer.data(snd + 249);
    const auto *snd_251 = buffer.data(snd + 251);
    const auto *snd_252 = buffer.data(snd + 252);
    const auto *snd_255 = buffer.data(snd + 255);
    const auto *snd_257 = buffer.data(snd + 257);
    const auto *snd_258 = buffer.data(snd + 258);
    const auto *snd_261 = buffer.data(snd + 261);
    const auto *snd_263 = buffer.data(snd + 263);
    const auto *snd_264 = buffer.data(snd + 264);
    const auto *snd_267 = buffer.data(snd + 267);
    const auto *snd_269 = buffer.data(snd + 269);
    const auto *snd_270 = buffer.data(snd + 270);
    const auto *snd_273 = buffer.data(snd + 273);
    const auto *snd_275 = buffer.data(snd + 275);
    const auto *snd_276 = buffer.data(snd + 276);
    const auto *snd_279 = buffer.data(snd + 279);
    const auto *snd_281 = buffer.data(snd + 281);
    const auto *snd_282 = buffer.data(snd + 282);
    const auto *snd_285 = buffer.data(snd + 285);
    const auto *snd_287 = buffer.data(snd + 287);
    const auto *snd_288 = buffer.data(snd + 288);
    const auto *snd_291 = buffer.data(snd + 291);
    const auto *snd_293 = buffer.data(snd + 293);
    const auto *snd_294 = buffer.data(snd + 294);
    const auto *snd_297 = buffer.data(snd + 297);
    const auto *snd_299 = buffer.data(snd + 299);
    const auto *snd_300 = buffer.data(snd + 300);
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

    const auto *snf1_440 = buffer.data(snf1 + 440);
    const auto *snf1_449 = buffer.data(snf1 + 449);
    const auto *snf1_450 = buffer.data(snf1 + 450);
    const auto *snf1_550 = buffer.data(snf1 + 550);
    const auto *snf1_556 = buffer.data(snf1 + 556);
    const auto *snf1_559 = buffer.data(snf1 + 559);
    const auto *snf1_566 = buffer.data(snf1 + 566);
    const auto *snf1_569 = buffer.data(snf1 + 569);
    const auto *snf1_570 = buffer.data(snf1 + 570);
    const auto *snf1_576 = buffer.data(snf1 + 576);
    const auto *snf1_579 = buffer.data(snf1 + 579);
    const auto *snf1_580 = buffer.data(snf1 + 580);
    const auto *snf1_586 = buffer.data(snf1 + 586);
    const auto *snf1_589 = buffer.data(snf1 + 589);
    const auto *snf1_590 = buffer.data(snf1 + 590);
    const auto *snf1_596 = buffer.data(snf1 + 596);
    const auto *snf1_599 = buffer.data(snf1 + 599);
    const auto *snf1_600 = buffer.data(snf1 + 600);
    const auto *snf1_606 = buffer.data(snf1 + 606);
    const auto *snf1_609 = buffer.data(snf1 + 609);
    const auto *snf1_610 = buffer.data(snf1 + 610);
    const auto *snf1_616 = buffer.data(snf1 + 616);
    const auto *snf1_619 = buffer.data(snf1 + 619);
    const auto *snf1_620 = buffer.data(snf1 + 620);
    const auto *snf1_626 = buffer.data(snf1 + 626);

    const auto *sop0_151 = buffer.data(sop0 + 151);
    const auto *sop0_152 = buffer.data(sop0 + 152);
    const auto *sop0_153 = buffer.data(sop0 + 153);
    const auto *sop0_154 = buffer.data(sop0 + 154);
    const auto *sop0_155 = buffer.data(sop0 + 155);
    const auto *sop0_156 = buffer.data(sop0 + 156);
    const auto *sop0_157 = buffer.data(sop0 + 157);
    const auto *sop0_158 = buffer.data(sop0 + 158);
    const auto *sop0_160 = buffer.data(sop0 + 160);
    const auto *sop0_162 = buffer.data(sop0 + 162);
    const auto *sop0_163 = buffer.data(sop0 + 163);
    const auto *sop0_164 = buffer.data(sop0 + 164);

    const auto *sop1_151 = buffer.data(sop1 + 151);
    const auto *sop1_152 = buffer.data(sop1 + 152);
    const auto *sop1_153 = buffer.data(sop1 + 153);
    const auto *sop1_154 = buffer.data(sop1 + 154);
    const auto *sop1_155 = buffer.data(sop1 + 155);
    const auto *sop1_156 = buffer.data(sop1 + 156);
    const auto *sop1_157 = buffer.data(sop1 + 157);
    const auto *sop1_158 = buffer.data(sop1 + 158);
    const auto *sop1_160 = buffer.data(sop1 + 160);
    const auto *sop1_162 = buffer.data(sop1 + 162);
    const auto *sop1_163 = buffer.data(sop1 + 163);
    const auto *sop1_164 = buffer.data(sop1 + 164);

    const auto *sod_300 = buffer.data(sod + 300);
    const auto *sod_303 = buffer.data(sod + 303);
    const auto *sod_304 = buffer.data(sod + 304);
    const auto *sod_305 = buffer.data(sod + 305);
    const auto *sod_306 = buffer.data(sod + 306);
    const auto *sod_309 = buffer.data(sod + 309);
    const auto *sod_310 = buffer.data(sod + 310);
    const auto *sod_311 = buffer.data(sod + 311);
    const auto *sod_312 = buffer.data(sod + 312);
    const auto *sod_315 = buffer.data(sod + 315);
    const auto *sod_316 = buffer.data(sod + 316);
    const auto *sod_317 = buffer.data(sod + 317);
    const auto *sod_318 = buffer.data(sod + 318);
    const auto *sod_321 = buffer.data(sod + 321);
    const auto *sod_322 = buffer.data(sod + 322);
    const auto *sod_323 = buffer.data(sod + 323);
    const auto *sod_324 = buffer.data(sod + 324);
    const auto *sod_327 = buffer.data(sod + 327);
    const auto *sod_328 = buffer.data(sod + 328);
    const auto *sod_329 = buffer.data(sod + 329);
    const auto *sod_330 = buffer.data(sod + 330);
    const auto *sod_333 = buffer.data(sod + 333);
    const auto *sod_334 = buffer.data(sod + 334);
    const auto *sod_335 = buffer.data(sod + 335);
    const auto *sod_336 = buffer.data(sod + 336);
    const auto *sod_339 = buffer.data(sod + 339);
    const auto *sod_340 = buffer.data(sod + 340);
    const auto *sod_341 = buffer.data(sod + 341);
    const auto *sod_342 = buffer.data(sod + 342);
    const auto *sod_345 = buffer.data(sod + 345);
    const auto *sod_346 = buffer.data(sod + 346);
    const auto *sod_347 = buffer.data(sod + 347);
    const auto *sod_348 = buffer.data(sod + 348);
    const auto *sod_351 = buffer.data(sod + 351);
    const auto *sod_352 = buffer.data(sod + 352);
    const auto *sod_353 = buffer.data(sod + 353);
    const auto *sod_354 = buffer.data(sod + 354);
    const auto *sod_357 = buffer.data(sod + 357);
    const auto *sod_358 = buffer.data(sod + 358);
    const auto *sod_359 = buffer.data(sod + 359);
    const auto *sod_360 = buffer.data(sod + 360);
    const auto *sod_363 = buffer.data(sod + 363);
    const auto *sod_364 = buffer.data(sod + 364);
    const auto *sod_365 = buffer.data(sod + 365);
    const auto *sod_366 = buffer.data(sod + 366);
    const auto *sod_369 = buffer.data(sod + 369);
    const auto *sod_370 = buffer.data(sod + 370);
    const auto *sod_371 = buffer.data(sod + 371);
    const auto *sod_372 = buffer.data(sod + 372);
    const auto *sod_375 = buffer.data(sod + 375);
    const auto *sod_376 = buffer.data(sod + 376);
    const auto *sod_377 = buffer.data(sod + 377);

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pc_x, pc_z, snd_240, snd_303, snd_304, \
                         snd_305, sod_300, sod_303, sod_304, sod_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_14 * snd_240[k]
                   + f_3 * pc_z[k] * sod_300[k];

        t_503[k] = f_8 * snd_303[k]
                   + f_3 * pc_x[k] * sod_303[k];

        t_504[k] = f_8 * snd_304[k]
                   + f_3 * pc_x[k] * sod_304[k];

        t_505[k] = f_8 * snd_305[k]
                   + f_3 * pc_x[k] * sod_305[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, t_509, pc_y, pc_z, snd_243, snd_245, snd_249, \
                         snd_251, sop0_151, sop0_152, sop1_151, sop1_152, sod_303, \
                         sod_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_12 * snd_249[k]
                   + f_1 * sop0_151[k]
                   - f_2 * sop1_151[k]
                   + f_3 * pc_y[k] * sod_303[k];

        t_507[k] = f_14 * snd_243[k]
                   + f_3 * pc_z[k] * sod_303[k];

        t_508[k] = f_12 * snd_251[k]
                   + f_3 * pc_y[k] * sod_305[k];

        t_509[k] = f_14 * snd_245[k]
                   + f_1 * sop0_152[k]
                   - f_2 * sop1_152[k]
                   + f_3 * pc_z[k] * sod_305[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, t_513, pc_x, pc_y, pc_z, snd_246, snd_252, \
                         snd_306, snd_309, sop0_153, sop1_153, sod_306, \
                         sod_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_8 * snd_306[k]
                   + f_1 * sop0_153[k]
                   - f_2 * sop1_153[k]
                   + f_3 * pc_x[k] * sod_306[k];

        t_511[k] = f_10 * snd_252[k]
                   + f_3 * pc_y[k] * sod_306[k];

        t_512[k] = f_13 * snd_246[k]
                   + f_3 * pc_z[k] * sod_306[k];

        t_513[k] = f_8 * snd_309[k]
                   + f_3 * pc_x[k] * sod_309[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, pc_x, pc_y, pc_z, snd_249, snd_255, \
                         snd_310, snd_311, sop0_154, sop1_154, sod_309, sod_310, \
                         sod_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_8 * snd_310[k]
                   + f_3 * pc_x[k] * sod_310[k];

        t_515[k] = f_8 * snd_311[k]
                   + f_3 * pc_x[k] * sod_311[k];

        t_516[k] = f_10 * snd_255[k]
                   + f_1 * sop0_154[k]
                   - f_2 * sop1_154[k]
                   + f_3 * pc_y[k] * sod_309[k];

        t_517[k] = f_13 * snd_249[k]
                   + f_3 * pc_z[k] * sod_309[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pc_x, pc_y, pc_z, snd_251, snd_257, snd_312, \
                         sop0_155, sop0_156, sop1_155, sop1_156, sod_311, \
                         sod_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_10 * snd_257[k]
                   + f_3 * pc_y[k] * sod_311[k];

        t_519[k] = f_13 * snd_251[k]
                   + f_1 * sop0_155[k]
                   - f_2 * sop1_155[k]
                   + f_3 * pc_z[k] * sod_311[k];

        t_520[k] = f_8 * snd_312[k]
                   + f_1 * sop0_156[k]
                   - f_2 * sop1_156[k]
                   + f_3 * pc_x[k] * sod_312[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, pc_z, snd_252, snd_258, \
                         snd_315, snd_316, sod_312, sod_315, sod_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_8 * snd_258[k]
                   + f_3 * pc_y[k] * sod_312[k];

        t_522[k] = f_11 * snd_252[k]
                   + f_3 * pc_z[k] * sod_312[k];

        t_523[k] = f_8 * snd_315[k]
                   + f_3 * pc_x[k] * sod_315[k];

        t_524[k] = f_8 * snd_316[k]
                   + f_3 * pc_x[k] * sod_316[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, pc_x, pc_y, pc_z, snd_255, snd_261, \
                         snd_263, snd_317, sop0_157, sop1_157, sod_315, \
                         sod_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_8 * snd_317[k]
                   + f_3 * pc_x[k] * sod_317[k];

        t_526[k] = f_8 * snd_261[k]
                   + f_1 * sop0_157[k]
                   - f_2 * sop1_157[k]
                   + f_3 * pc_y[k] * sod_315[k];

        t_527[k] = f_11 * snd_255[k]
                   + f_3 * pc_z[k] * sod_315[k];

        t_528[k] = f_8 * snd_263[k]
                   + f_3 * pc_y[k] * sod_317[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pb_y, pc_y, pc_z, snf0_440, snd_257, \
                         snd_258, snd_264, snf1_440, sop0_158, sop1_158, sod_317, \
                         sod_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_11 * snd_257[k]
                   + f_1 * sop0_158[k]
                   - f_2 * sop1_158[k]
                   + f_3 * pc_z[k] * sod_317[k];

        t_530[k] = pb_y[k] * snf0_440[k]
                   - f_4 * pc_y[k] * snf1_440[k];

        t_531[k] = f_5 * snd_264[k]
                   + f_3 * pc_y[k] * sod_318[k];

        t_532[k] = f_9 * snd_258[k]
                   + f_3 * pc_z[k] * sod_318[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, snd_267, snd_321, snd_322, \
                         snd_323, sop0_160, sop1_160, sod_321, sod_322, \
                         sod_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_8 * snd_321[k]
                   + f_3 * pc_x[k] * sod_321[k];

        t_534[k] = f_8 * snd_322[k]
                   + f_3 * pc_x[k] * sod_322[k];

        t_535[k] = f_8 * snd_323[k]
                   + f_3 * pc_x[k] * sod_323[k];

        t_536[k] = f_5 * snd_267[k]
                   + f_1 * sop0_160[k]
                   - f_2 * sop1_160[k]
                   + f_3 * pc_y[k] * sod_321[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pb_y, pc_y, pc_z, snf0_449, snd_261, snd_269, \
                         snf1_449, sod_321, sod_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_9 * snd_261[k]
                   + f_3 * pc_z[k] * sod_321[k];

        t_538[k] = f_5 * snd_269[k]
                   + f_3 * pc_y[k] * sod_323[k];

        t_539[k] = pb_y[k] * snf0_449[k]
                   - f_4 * pc_y[k] * snf1_449[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, t_543, pc_x, pc_y, pc_z, snd_264, snd_324, \
                         snd_327, sop0_162, sop1_162, sod_324, \
                         sod_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_8 * snd_324[k]
                   + f_1 * sop0_162[k]
                   - f_2 * sop1_162[k]
                   + f_3 * pc_x[k] * sod_324[k];

        t_541[k] = f_3 * pc_y[k] * sod_324[k];

        t_542[k] = f_7 * snd_264[k]
                   + f_3 * pc_z[k] * sod_324[k];

        t_543[k] = f_8 * snd_327[k]
                   + f_3 * pc_x[k] * sod_327[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, t_548, pc_x, pc_y, pc_z, snd_267, \
                         snd_328, snd_329, sop0_163, sop1_163, sod_327, sod_328, \
                         sod_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_8 * snd_328[k]
                   + f_3 * pc_x[k] * sod_328[k];

        t_545[k] = f_8 * snd_329[k]
                   + f_3 * pc_x[k] * sod_329[k];

        t_546[k] = f_1 * sop0_163[k]
                   - f_2 * sop1_163[k]
                   + f_3 * pc_y[k] * sod_327[k];

        t_547[k] = f_7 * snd_267[k]
                   + f_3 * pc_z[k] * sod_327[k];

        t_548[k] = f_3 * pc_y[k] * sod_329[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pb_x, pc_x, pc_y, pc_z, snf0_550, snd_269, \
                         snd_270, snd_330, snf1_550, sop0_164, sop1_164, sod_329, \
                         sod_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_7 * snd_269[k]
                   + f_1 * sop0_164[k]
                   - f_2 * sop1_164[k]
                   + f_3 * pc_z[k] * sod_329[k];

        t_550[k] = pb_x[k] * snf0_550[k]
                   + f_10 * snd_330[k]
                   - f_4 * pc_x[k] * snf1_550[k];

        t_551[k] = f_6 * snd_270[k]
                   + f_3 * pc_y[k] * sod_330[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pc_x, pc_z, snd_333, snd_334, snd_335, \
                         sod_330, sod_333, sod_334, sod_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = f_3 * pc_z[k] * sod_330[k];

        t_553[k] = f_5 * snd_333[k]
                   + f_3 * pc_x[k] * sod_333[k];

        t_554[k] = f_5 * snd_334[k]
                   + f_3 * pc_x[k] * sod_334[k];

        t_555[k] = f_5 * snd_335[k]
                   + f_3 * pc_x[k] * sod_335[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, pb_x, pc_x, pc_y, pc_z, snf0_556, \
                         snf0_559, snd_275, snf1_556, snf1_559, sod_333, \
                         sod_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = pb_x[k] * snf0_556[k]
                   - f_4 * pc_x[k] * snf1_556[k];

        t_557[k] = f_3 * pc_z[k] * sod_333[k];

        t_558[k] = f_6 * snd_275[k]
                   + f_3 * pc_y[k] * sod_335[k];

        t_559[k] = pb_x[k] * snf0_559[k]
                   - f_4 * pc_x[k] * snf1_559[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, pb_z, pc_x, pc_y, pc_z, snf0_450, \
                         snd_270, snd_276, snd_339, snf1_450, sod_336, \
                         sod_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = pb_z[k] * snf0_450[k]
                   - f_4 * pc_z[k] * snf1_450[k];

        t_561[k] = f_7 * snd_276[k]
                   + f_3 * pc_y[k] * sod_336[k];

        t_562[k] = f_5 * snd_270[k]
                   + f_3 * pc_z[k] * sod_336[k];

        t_563[k] = f_5 * snd_339[k]
                   + f_3 * pc_x[k] * sod_339[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, t_567, pb_x, pc_x, pc_z, snf0_566, snd_273, \
                         snd_340, snd_341, snf1_566, sod_339, sod_340, \
                         sod_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_5 * snd_340[k]
                   + f_3 * pc_x[k] * sod_340[k];

        t_565[k] = f_5 * snd_341[k]
                   + f_3 * pc_x[k] * sod_341[k];

        t_566[k] = pb_x[k] * snf0_566[k]
                   - f_4 * pc_x[k] * snf1_566[k];

        t_567[k] = f_5 * snd_273[k]
                   + f_3 * pc_z[k] * sod_339[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, t_571, pb_x, pc_x, pc_y, snf0_569, snf0_570, \
                         snd_281, snd_282, snd_342, snf1_569, snf1_570, sod_341, \
                         sod_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_7 * snd_281[k]
                   + f_3 * pc_y[k] * sod_341[k];

        t_569[k] = pb_x[k] * snf0_569[k]
                   - f_4 * pc_x[k] * snf1_569[k];

        t_570[k] = pb_x[k] * snf0_570[k]
                   + f_10 * snd_342[k]
                   - f_4 * pc_x[k] * snf1_570[k];

        t_571[k] = f_9 * snd_282[k]
                   + f_3 * pc_y[k] * sod_342[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_z, snd_276, snd_345, snd_346, \
                         snd_347, sod_342, sod_345, sod_346, sod_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_8 * snd_276[k]
                   + f_3 * pc_z[k] * sod_342[k];

        t_573[k] = f_5 * snd_345[k]
                   + f_3 * pc_x[k] * sod_345[k];

        t_574[k] = f_5 * snd_346[k]
                   + f_3 * pc_x[k] * sod_346[k];

        t_575[k] = f_5 * snd_347[k]
                   + f_3 * pc_x[k] * sod_347[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pb_x, pc_x, pc_y, pc_z, snf0_576, \
                         snf0_579, snd_279, snd_287, snf1_576, snf1_579, sod_345, \
                         sod_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = pb_x[k] * snf0_576[k]
                   - f_4 * pc_x[k] * snf1_576[k];

        t_577[k] = f_8 * snd_279[k]
                   + f_3 * pc_z[k] * sod_345[k];

        t_578[k] = f_9 * snd_287[k]
                   + f_3 * pc_y[k] * sod_347[k];

        t_579[k] = pb_x[k] * snf0_579[k]
                   - f_4 * pc_x[k] * snf1_579[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pb_x, pc_x, pc_y, pc_z, snf0_580, \
                         snd_282, snd_288, snd_348, snd_351, snf1_580, sod_348, \
                         sod_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = pb_x[k] * snf0_580[k]
                   + f_10 * snd_348[k]
                   - f_4 * pc_x[k] * snf1_580[k];

        t_581[k] = f_11 * snd_288[k]
                   + f_3 * pc_y[k] * sod_348[k];

        t_582[k] = f_10 * snd_282[k]
                   + f_3 * pc_z[k] * sod_348[k];

        t_583[k] = f_5 * snd_351[k]
                   + f_3 * pc_x[k] * sod_351[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_x, pc_x, pc_z, snf0_586, snd_285, \
                         snd_352, snd_353, snf1_586, sod_351, sod_352, \
                         sod_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_5 * snd_352[k]
                   + f_3 * pc_x[k] * sod_352[k];

        t_585[k] = f_5 * snd_353[k]
                   + f_3 * pc_x[k] * sod_353[k];

        t_586[k] = pb_x[k] * snf0_586[k]
                   - f_4 * pc_x[k] * snf1_586[k];

        t_587[k] = f_10 * snd_285[k]
                   + f_3 * pc_z[k] * sod_351[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_x, pc_x, pc_y, snf0_589, snf0_590, \
                         snd_293, snd_294, snd_354, snf1_589, snf1_590, sod_353, \
                         sod_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_11 * snd_293[k]
                   + f_3 * pc_y[k] * sod_353[k];

        t_589[k] = pb_x[k] * snf0_589[k]
                   - f_4 * pc_x[k] * snf1_589[k];

        t_590[k] = pb_x[k] * snf0_590[k]
                   + f_10 * snd_354[k]
                   - f_4 * pc_x[k] * snf1_590[k];

        t_591[k] = f_13 * snd_294[k]
                   + f_3 * pc_y[k] * sod_354[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, snd_288, snd_357, snd_358, \
                         snd_359, sod_354, sod_357, sod_358, sod_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_12 * snd_288[k]
                   + f_3 * pc_z[k] * sod_354[k];

        t_593[k] = f_5 * snd_357[k]
                   + f_3 * pc_x[k] * sod_357[k];

        t_594[k] = f_5 * snd_358[k]
                   + f_3 * pc_x[k] * sod_358[k];

        t_595[k] = f_5 * snd_359[k]
                   + f_3 * pc_x[k] * sod_359[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_x, pc_x, pc_y, pc_z, snf0_596, \
                         snf0_599, snd_291, snd_299, snf1_596, snf1_599, sod_357, \
                         sod_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = pb_x[k] * snf0_596[k]
                   - f_4 * pc_x[k] * snf1_596[k];

        t_597[k] = f_12 * snd_291[k]
                   + f_3 * pc_z[k] * sod_357[k];

        t_598[k] = f_13 * snd_299[k]
                   + f_3 * pc_y[k] * sod_359[k];

        t_599[k] = pb_x[k] * snf0_599[k]
                   - f_4 * pc_x[k] * snf1_599[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pb_x, pc_x, pc_y, pc_z, snf0_600, \
                         snd_294, snd_300, snd_360, snd_363, snf1_600, sod_360, \
                         sod_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pb_x[k] * snf0_600[k]
                   + f_10 * snd_360[k]
                   - f_4 * pc_x[k] * snf1_600[k];

        t_601[k] = f_14 * snd_300[k]
                   + f_3 * pc_y[k] * sod_360[k];

        t_602[k] = f_14 * snd_294[k]
                   + f_3 * pc_z[k] * sod_360[k];

        t_603[k] = f_5 * snd_363[k]
                   + f_3 * pc_x[k] * sod_363[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pb_x, pc_x, pc_z, snf0_606, snd_297, \
                         snd_364, snd_365, snf1_606, sod_363, sod_364, \
                         sod_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_5 * snd_364[k]
                   + f_3 * pc_x[k] * sod_364[k];

        t_605[k] = f_5 * snd_365[k]
                   + f_3 * pc_x[k] * sod_365[k];

        t_606[k] = pb_x[k] * snf0_606[k]
                   - f_4 * pc_x[k] * snf1_606[k];

        t_607[k] = f_14 * snd_297[k]
                   + f_3 * pc_z[k] * sod_363[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pb_x, pc_x, pc_y, snf0_609, snf0_610, \
                         snd_305, snd_306, snd_366, snf1_609, snf1_610, sod_365, \
                         sod_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_14 * snd_305[k]
                   + f_3 * pc_y[k] * sod_365[k];

        t_609[k] = pb_x[k] * snf0_609[k]
                   - f_4 * pc_x[k] * snf1_609[k];

        t_610[k] = pb_x[k] * snf0_610[k]
                   + f_10 * snd_366[k]
                   - f_4 * pc_x[k] * snf1_610[k];

        t_611[k] = f_12 * snd_306[k]
                   + f_3 * pc_y[k] * sod_366[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_x, pc_z, snd_300, snd_369, snd_370, \
                         snd_371, sod_366, sod_369, sod_370, sod_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_13 * snd_300[k]
                   + f_3 * pc_z[k] * sod_366[k];

        t_613[k] = f_5 * snd_369[k]
                   + f_3 * pc_x[k] * sod_369[k];

        t_614[k] = f_5 * snd_370[k]
                   + f_3 * pc_x[k] * sod_370[k];

        t_615[k] = f_5 * snd_371[k]
                   + f_3 * pc_x[k] * sod_371[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pb_x, pc_x, pc_y, pc_z, snf0_616, \
                         snf0_619, snd_303, snd_311, snf1_616, snf1_619, sod_369, \
                         sod_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pb_x[k] * snf0_616[k]
                   - f_4 * pc_x[k] * snf1_616[k];

        t_617[k] = f_13 * snd_303[k]
                   + f_3 * pc_z[k] * sod_369[k];

        t_618[k] = f_12 * snd_311[k]
                   + f_3 * pc_y[k] * sod_371[k];

        t_619[k] = pb_x[k] * snf0_619[k]
                   - f_4 * pc_x[k] * snf1_619[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, t_623, pb_x, pc_x, pc_y, pc_z, snf0_620, \
                         snd_306, snd_312, snd_372, snd_375, snf1_620, sod_372, \
                         sod_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = pb_x[k] * snf0_620[k]
                   + f_10 * snd_372[k]
                   - f_4 * pc_x[k] * snf1_620[k];

        t_621[k] = f_10 * snd_312[k]
                   + f_3 * pc_y[k] * sod_372[k];

        t_622[k] = f_11 * snd_306[k]
                   + f_3 * pc_z[k] * sod_372[k];

        t_623[k] = f_5 * snd_375[k]
                   + f_3 * pc_x[k] * sod_375[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, t_627, pb_x, pc_x, pc_z, snf0_626, snd_309, \
                         snd_376, snd_377, snf1_626, sod_375, sod_376, \
                         sod_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_5 * snd_376[k]
                   + f_3 * pc_x[k] * sod_376[k];

        t_625[k] = f_5 * snd_377[k]
                   + f_3 * pc_x[k] * sod_377[k];

        t_626[k] = pb_x[k] * snf0_626[k]
                   - f_4 * pc_x[k] * snf1_626[k];

        t_627[k] = f_11 * snd_309[k]
                   + f_3 * pc_z[k] * sod_375[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 4.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 4.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 3.5 / q;
    const auto f_12 = 2.0 / q;
    const auto f_13 = 3.0 / q;
    const auto f_14 = 2.5 / q;

    auto *t_628 = buffer.data(target + 628);
    auto *t_629 = buffer.data(target + 629);
    auto *t_630 = buffer.data(target + 630);
    auto *t_631 = buffer.data(target + 631);
    auto *t_632 = buffer.data(target + 632);
    auto *t_633 = buffer.data(target + 633);
    auto *t_634 = buffer.data(target + 634);
    auto *t_635 = buffer.data(target + 635);
    auto *t_636 = buffer.data(target + 636);
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
    auto *t_660 = buffer.data(target + 660);
    auto *t_661 = buffer.data(target + 661);
    auto *t_662 = buffer.data(target + 662);
    auto *t_663 = buffer.data(target + 663);
    auto *t_664 = buffer.data(target + 664);
    auto *t_665 = buffer.data(target + 665);
    auto *t_666 = buffer.data(target + 666);
    auto *t_667 = buffer.data(target + 667);
    auto *t_668 = buffer.data(target + 668);
    auto *t_669 = buffer.data(target + 669);
    auto *t_670 = buffer.data(target + 670);
    auto *t_671 = buffer.data(target + 671);
    auto *t_672 = buffer.data(target + 672);
    auto *t_673 = buffer.data(target + 673);
    auto *t_674 = buffer.data(target + 674);
    auto *t_675 = buffer.data(target + 675);
    auto *t_676 = buffer.data(target + 676);
    auto *t_677 = buffer.data(target + 677);
    auto *t_678 = buffer.data(target + 678);
    auto *t_679 = buffer.data(target + 679);
    auto *t_680 = buffer.data(target + 680);
    auto *t_681 = buffer.data(target + 681);
    auto *t_682 = buffer.data(target + 682);
    auto *t_683 = buffer.data(target + 683);
    auto *t_684 = buffer.data(target + 684);
    auto *t_685 = buffer.data(target + 685);
    auto *t_686 = buffer.data(target + 686);
    auto *t_687 = buffer.data(target + 687);
    auto *t_688 = buffer.data(target + 688);
    auto *t_689 = buffer.data(target + 689);
    auto *t_690 = buffer.data(target + 690);
    auto *t_691 = buffer.data(target + 691);
    auto *t_692 = buffer.data(target + 692);
    auto *t_693 = buffer.data(target + 693);
    auto *t_694 = buffer.data(target + 694);
    auto *t_695 = buffer.data(target + 695);
    auto *t_696 = buffer.data(target + 696);
    auto *t_697 = buffer.data(target + 697);
    auto *t_698 = buffer.data(target + 698);
    auto *t_699 = buffer.data(target + 699);
    auto *t_700 = buffer.data(target + 700);
    auto *t_701 = buffer.data(target + 701);
    auto *t_702 = buffer.data(target + 702);
    auto *t_703 = buffer.data(target + 703);
    auto *t_704 = buffer.data(target + 704);
    auto *t_705 = buffer.data(target + 705);
    auto *t_706 = buffer.data(target + 706);
    auto *t_707 = buffer.data(target + 707);
    auto *t_708 = buffer.data(target + 708);
    auto *t_709 = buffer.data(target + 709);
    auto *t_710 = buffer.data(target + 710);
    auto *t_711 = buffer.data(target + 711);
    auto *t_712 = buffer.data(target + 712);
    auto *t_713 = buffer.data(target + 713);
    auto *t_714 = buffer.data(target + 714);
    auto *t_715 = buffer.data(target + 715);
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);
    auto *t_720 = buffer.data(target + 720);
    auto *t_721 = buffer.data(target + 721);
    auto *t_722 = buffer.data(target + 722);
    auto *t_723 = buffer.data(target + 723);
    auto *t_724 = buffer.data(target + 724);
    auto *t_725 = buffer.data(target + 725);
    auto *t_726 = buffer.data(target + 726);
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);
    auto *t_731 = buffer.data(target + 731);
    auto *t_732 = buffer.data(target + 732);
    auto *t_733 = buffer.data(target + 733);
    auto *t_734 = buffer.data(target + 734);
    auto *t_735 = buffer.data(target + 735);
    auto *t_736 = buffer.data(target + 736);
    auto *t_737 = buffer.data(target + 737);
    auto *t_738 = buffer.data(target + 738);
    auto *t_739 = buffer.data(target + 739);
    auto *t_740 = buffer.data(target + 740);
    auto *t_741 = buffer.data(target + 741);
    auto *t_742 = buffer.data(target + 742);
    auto *t_743 = buffer.data(target + 743);
    auto *t_744 = buffer.data(target + 744);
    auto *t_745 = buffer.data(target + 745);
    auto *t_746 = buffer.data(target + 746);
    auto *t_747 = buffer.data(target + 747);
    auto *t_748 = buffer.data(target + 748);
    auto *t_749 = buffer.data(target + 749);
    auto *t_750 = buffer.data(target + 750);
    auto *t_751 = buffer.data(target + 751);
    auto *t_752 = buffer.data(target + 752);
    auto *t_753 = buffer.data(target + 753);
    auto *t_754 = buffer.data(target + 754);
    auto *t_755 = buffer.data(target + 755);
    auto *t_756 = buffer.data(target + 756);
    auto *t_757 = buffer.data(target + 757);
    auto *t_758 = buffer.data(target + 758);
    auto *t_759 = buffer.data(target + 759);
    auto *t_760 = buffer.data(target + 760);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snf0_540 = buffer.data(snf0 + 540);
    const auto *snf0_550 = buffer.data(snf0 + 550);
    const auto *snf0_556 = buffer.data(snf0 + 556);
    const auto *snf0_629 = buffer.data(snf0 + 629);
    const auto *snf0_630 = buffer.data(snf0 + 630);
    const auto *snf0_636 = buffer.data(snf0 + 636);
    const auto *snf0_639 = buffer.data(snf0 + 639);
    const auto *snf0_646 = buffer.data(snf0 + 646);
    const auto *snf0_649 = buffer.data(snf0 + 649);
    const auto *snf0_650 = buffer.data(snf0 + 650);
    const auto *snf0_656 = buffer.data(snf0 + 656);
    const auto *snf0_659 = buffer.data(snf0 + 659);

    const auto *snd_312 = buffer.data(snd + 312);
    const auto *snd_315 = buffer.data(snd + 315);
    const auto *snd_317 = buffer.data(snd + 317);
    const auto *snd_318 = buffer.data(snd + 318);
    const auto *snd_321 = buffer.data(snd + 321);
    const auto *snd_323 = buffer.data(snd + 323);
    const auto *snd_324 = buffer.data(snd + 324);
    const auto *snd_327 = buffer.data(snd + 327);
    const auto *snd_329 = buffer.data(snd + 329);
    const auto *snd_330 = buffer.data(snd + 330);
    const auto *snd_333 = buffer.data(snd + 333);
    const auto *snd_335 = buffer.data(snd + 335);
    const auto *snd_336 = buffer.data(snd + 336);
    const auto *snd_339 = buffer.data(snd + 339);
    const auto *snd_341 = buffer.data(snd + 341);
    const auto *snd_342 = buffer.data(snd + 342);
    const auto *snd_345 = buffer.data(snd + 345);
    const auto *snd_347 = buffer.data(snd + 347);
    const auto *snd_348 = buffer.data(snd + 348);
    const auto *snd_351 = buffer.data(snd + 351);
    const auto *snd_353 = buffer.data(snd + 353);
    const auto *snd_354 = buffer.data(snd + 354);
    const auto *snd_357 = buffer.data(snd + 357);
    const auto *snd_359 = buffer.data(snd + 359);
    const auto *snd_360 = buffer.data(snd + 360);
    const auto *snd_363 = buffer.data(snd + 363);
    const auto *snd_365 = buffer.data(snd + 365);
    const auto *snd_366 = buffer.data(snd + 366);
    const auto *snd_369 = buffer.data(snd + 369);
    const auto *snd_371 = buffer.data(snd + 371);
    const auto *snd_372 = buffer.data(snd + 372);
    const auto *snd_375 = buffer.data(snd + 375);
    const auto *snd_377 = buffer.data(snd + 377);
    const auto *snd_378 = buffer.data(snd + 378);
    const auto *snd_381 = buffer.data(snd + 381);
    const auto *snd_382 = buffer.data(snd + 382);
    const auto *snd_383 = buffer.data(snd + 383);
    const auto *snd_384 = buffer.data(snd + 384);
    const auto *snd_387 = buffer.data(snd + 387);
    const auto *snd_388 = buffer.data(snd + 388);
    const auto *snd_389 = buffer.data(snd + 389);
    const auto *snd_390 = buffer.data(snd + 390);
    const auto *snd_393 = buffer.data(snd + 393);
    const auto *snd_394 = buffer.data(snd + 394);
    const auto *snd_395 = buffer.data(snd + 395);

    const auto *snf1_540 = buffer.data(snf1 + 540);
    const auto *snf1_550 = buffer.data(snf1 + 550);
    const auto *snf1_556 = buffer.data(snf1 + 556);
    const auto *snf1_629 = buffer.data(snf1 + 629);
    const auto *snf1_630 = buffer.data(snf1 + 630);
    const auto *snf1_636 = buffer.data(snf1 + 636);
    const auto *snf1_639 = buffer.data(snf1 + 639);
    const auto *snf1_646 = buffer.data(snf1 + 646);
    const auto *snf1_649 = buffer.data(snf1 + 649);
    const auto *snf1_650 = buffer.data(snf1 + 650);
    const auto *snf1_656 = buffer.data(snf1 + 656);
    const auto *snf1_659 = buffer.data(snf1 + 659);

    const auto *sop0_198 = buffer.data(sop0 + 198);
    const auto *sop0_199 = buffer.data(sop0 + 199);
    const auto *sop0_200 = buffer.data(sop0 + 200);
    const auto *sop0_203 = buffer.data(sop0 + 203);
    const auto *sop0_204 = buffer.data(sop0 + 204);
    const auto *sop0_205 = buffer.data(sop0 + 205);
    const auto *sop0_206 = buffer.data(sop0 + 206);
    const auto *sop0_207 = buffer.data(sop0 + 207);
    const auto *sop0_208 = buffer.data(sop0 + 208);
    const auto *sop0_209 = buffer.data(sop0 + 209);
    const auto *sop0_210 = buffer.data(sop0 + 210);
    const auto *sop0_211 = buffer.data(sop0 + 211);
    const auto *sop0_212 = buffer.data(sop0 + 212);
    const auto *sop0_213 = buffer.data(sop0 + 213);
    const auto *sop0_214 = buffer.data(sop0 + 214);
    const auto *sop0_215 = buffer.data(sop0 + 215);
    const auto *sop0_216 = buffer.data(sop0 + 216);
    const auto *sop0_217 = buffer.data(sop0 + 217);
    const auto *sop0_218 = buffer.data(sop0 + 218);
    const auto *sop0_219 = buffer.data(sop0 + 219);
    const auto *sop0_220 = buffer.data(sop0 + 220);
    const auto *sop0_221 = buffer.data(sop0 + 221);
    const auto *sop0_222 = buffer.data(sop0 + 222);
    const auto *sop0_223 = buffer.data(sop0 + 223);
    const auto *sop0_224 = buffer.data(sop0 + 224);
    const auto *sop0_225 = buffer.data(sop0 + 225);
    const auto *sop0_226 = buffer.data(sop0 + 226);
    const auto *sop0_227 = buffer.data(sop0 + 227);

    const auto *sop1_198 = buffer.data(sop1 + 198);
    const auto *sop1_199 = buffer.data(sop1 + 199);
    const auto *sop1_200 = buffer.data(sop1 + 200);
    const auto *sop1_203 = buffer.data(sop1 + 203);
    const auto *sop1_204 = buffer.data(sop1 + 204);
    const auto *sop1_205 = buffer.data(sop1 + 205);
    const auto *sop1_206 = buffer.data(sop1 + 206);
    const auto *sop1_207 = buffer.data(sop1 + 207);
    const auto *sop1_208 = buffer.data(sop1 + 208);
    const auto *sop1_209 = buffer.data(sop1 + 209);
    const auto *sop1_210 = buffer.data(sop1 + 210);
    const auto *sop1_211 = buffer.data(sop1 + 211);
    const auto *sop1_212 = buffer.data(sop1 + 212);
    const auto *sop1_213 = buffer.data(sop1 + 213);
    const auto *sop1_214 = buffer.data(sop1 + 214);
    const auto *sop1_215 = buffer.data(sop1 + 215);
    const auto *sop1_216 = buffer.data(sop1 + 216);
    const auto *sop1_217 = buffer.data(sop1 + 217);
    const auto *sop1_218 = buffer.data(sop1 + 218);
    const auto *sop1_219 = buffer.data(sop1 + 219);
    const auto *sop1_220 = buffer.data(sop1 + 220);
    const auto *sop1_221 = buffer.data(sop1 + 221);
    const auto *sop1_222 = buffer.data(sop1 + 222);
    const auto *sop1_223 = buffer.data(sop1 + 223);
    const auto *sop1_224 = buffer.data(sop1 + 224);
    const auto *sop1_225 = buffer.data(sop1 + 225);
    const auto *sop1_226 = buffer.data(sop1 + 226);
    const auto *sop1_227 = buffer.data(sop1 + 227);

    const auto *sod_377 = buffer.data(sod + 377);
    const auto *sod_378 = buffer.data(sod + 378);
    const auto *sod_381 = buffer.data(sod + 381);
    const auto *sod_382 = buffer.data(sod + 382);
    const auto *sod_383 = buffer.data(sod + 383);
    const auto *sod_384 = buffer.data(sod + 384);
    const auto *sod_387 = buffer.data(sod + 387);
    const auto *sod_388 = buffer.data(sod + 388);
    const auto *sod_389 = buffer.data(sod + 389);
    const auto *sod_390 = buffer.data(sod + 390);
    const auto *sod_393 = buffer.data(sod + 393);
    const auto *sod_394 = buffer.data(sod + 394);
    const auto *sod_395 = buffer.data(sod + 395);
    const auto *sod_396 = buffer.data(sod + 396);
    const auto *sod_399 = buffer.data(sod + 399);
    const auto *sod_400 = buffer.data(sod + 400);
    const auto *sod_401 = buffer.data(sod + 401);
    const auto *sod_402 = buffer.data(sod + 402);
    const auto *sod_405 = buffer.data(sod + 405);
    const auto *sod_406 = buffer.data(sod + 406);
    const auto *sod_407 = buffer.data(sod + 407);
    const auto *sod_408 = buffer.data(sod + 408);
    const auto *sod_411 = buffer.data(sod + 411);
    const auto *sod_412 = buffer.data(sod + 412);
    const auto *sod_413 = buffer.data(sod + 413);
    const auto *sod_414 = buffer.data(sod + 414);
    const auto *sod_417 = buffer.data(sod + 417);
    const auto *sod_418 = buffer.data(sod + 418);
    const auto *sod_419 = buffer.data(sod + 419);
    const auto *sod_420 = buffer.data(sod + 420);
    const auto *sod_423 = buffer.data(sod + 423);
    const auto *sod_424 = buffer.data(sod + 424);
    const auto *sod_425 = buffer.data(sod + 425);
    const auto *sod_426 = buffer.data(sod + 426);
    const auto *sod_429 = buffer.data(sod + 429);
    const auto *sod_430 = buffer.data(sod + 430);
    const auto *sod_431 = buffer.data(sod + 431);
    const auto *sod_432 = buffer.data(sod + 432);
    const auto *sod_435 = buffer.data(sod + 435);
    const auto *sod_436 = buffer.data(sod + 436);
    const auto *sod_437 = buffer.data(sod + 437);
    const auto *sod_438 = buffer.data(sod + 438);
    const auto *sod_441 = buffer.data(sod + 441);
    const auto *sod_442 = buffer.data(sod + 442);
    const auto *sod_443 = buffer.data(sod + 443);
    const auto *sod_444 = buffer.data(sod + 444);
    const auto *sod_447 = buffer.data(sod + 447);
    const auto *sod_448 = buffer.data(sod + 448);
    const auto *sod_449 = buffer.data(sod + 449);
    const auto *sod_450 = buffer.data(sod + 450);
    const auto *sod_453 = buffer.data(sod + 453);
    const auto *sod_454 = buffer.data(sod + 454);
    const auto *sod_455 = buffer.data(sod + 455);

#pragma omp simd aligned(t_628, t_629, t_630, t_631, pb_x, pc_x, pc_y, snf0_629, snf0_630, \
                         snd_317, snd_318, snd_378, snf1_629, snf1_630, sod_377, \
                         sod_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_10 * snd_317[k]
                   + f_3 * pc_y[k] * sod_377[k];

        t_629[k] = pb_x[k] * snf0_629[k]
                   - f_4 * pc_x[k] * snf1_629[k];

        t_630[k] = pb_x[k] * snf0_630[k]
                   + f_10 * snd_378[k]
                   - f_4 * pc_x[k] * snf1_630[k];

        t_631[k] = f_8 * snd_318[k]
                   + f_3 * pc_y[k] * sod_378[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, t_635, pc_x, pc_z, snd_312, snd_381, snd_382, \
                         snd_383, sod_378, sod_381, sod_382, sod_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_9 * snd_312[k]
                   + f_3 * pc_z[k] * sod_378[k];

        t_633[k] = f_5 * snd_381[k]
                   + f_3 * pc_x[k] * sod_381[k];

        t_634[k] = f_5 * snd_382[k]
                   + f_3 * pc_x[k] * sod_382[k];

        t_635[k] = f_5 * snd_383[k]
                   + f_3 * pc_x[k] * sod_383[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, t_639, pb_x, pc_x, pc_y, pc_z, snf0_636, \
                         snf0_639, snd_315, snd_323, snf1_636, snf1_639, sod_381, \
                         sod_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = pb_x[k] * snf0_636[k]
                   - f_4 * pc_x[k] * snf1_636[k];

        t_637[k] = f_9 * snd_315[k]
                   + f_3 * pc_z[k] * sod_381[k];

        t_638[k] = f_8 * snd_323[k]
                   + f_3 * pc_y[k] * sod_383[k];

        t_639[k] = pb_x[k] * snf0_639[k]
                   - f_4 * pc_x[k] * snf1_639[k];
    }

#pragma omp simd aligned(t_640, t_641, t_642, t_643, pb_y, pc_x, pc_y, pc_z, snf0_540, \
                         snd_318, snd_324, snd_387, snf1_540, sod_384, \
                         sod_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_640[k] = pb_y[k] * snf0_540[k]
                   - f_4 * pc_y[k] * snf1_540[k];

        t_641[k] = f_5 * snd_324[k]
                   + f_3 * pc_y[k] * sod_384[k];

        t_642[k] = f_7 * snd_318[k]
                   + f_3 * pc_z[k] * sod_384[k];

        t_643[k] = f_5 * snd_387[k]
                   + f_3 * pc_x[k] * sod_387[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, t_647, pb_x, pc_x, pc_z, snf0_646, snd_321, \
                         snd_388, snd_389, snf1_646, sod_387, sod_388, \
                         sod_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_5 * snd_388[k]
                   + f_3 * pc_x[k] * sod_388[k];

        t_645[k] = f_5 * snd_389[k]
                   + f_3 * pc_x[k] * sod_389[k];

        t_646[k] = pb_x[k] * snf0_646[k]
                   - f_4 * pc_x[k] * snf1_646[k];

        t_647[k] = f_7 * snd_321[k]
                   + f_3 * pc_z[k] * sod_387[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, t_651, pb_x, pc_x, pc_y, snf0_649, snf0_650, \
                         snd_329, snd_390, snf1_649, snf1_650, sod_389, \
                         sod_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_5 * snd_329[k]
                   + f_3 * pc_y[k] * sod_389[k];

        t_649[k] = pb_x[k] * snf0_649[k]
                   - f_4 * pc_x[k] * snf1_649[k];

        t_650[k] = pb_x[k] * snf0_650[k]
                   + f_10 * snd_390[k]
                   - f_4 * pc_x[k] * snf1_650[k];

        t_651[k] = f_3 * pc_y[k] * sod_390[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, t_655, pc_x, pc_z, snd_324, snd_393, snd_394, \
                         snd_395, sod_390, sod_393, sod_394, sod_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_6 * snd_324[k]
                   + f_3 * pc_z[k] * sod_390[k];

        t_653[k] = f_5 * snd_393[k]
                   + f_3 * pc_x[k] * sod_393[k];

        t_654[k] = f_5 * snd_394[k]
                   + f_3 * pc_x[k] * sod_394[k];

        t_655[k] = f_5 * snd_395[k]
                   + f_3 * pc_x[k] * sod_395[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pb_x, pc_x, pc_y, pc_z, snf0_656, \
                         snf0_659, snd_327, snf1_656, snf1_659, sod_393, \
                         sod_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = pb_x[k] * snf0_656[k]
                   - f_4 * pc_x[k] * snf1_656[k];

        t_657[k] = f_6 * snd_327[k]
                   + f_3 * pc_z[k] * sod_393[k];

        t_658[k] = f_3 * pc_y[k] * sod_395[k];

        t_659[k] = pb_x[k] * snf0_659[k]
                   - f_4 * pc_x[k] * snf1_659[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, pc_x, pc_y, pc_z, snd_330, \
                         sop0_198, sop1_198, sod_396, sod_399, sod_400, \
                         sod_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_1 * sop0_198[k]
                   - f_2 * sop1_198[k]
                   + f_3 * pc_x[k] * sod_396[k];

        t_661[k] = f_0 * snd_330[k]
                   + f_3 * pc_y[k] * sod_396[k];

        t_662[k] = f_3 * pc_z[k] * sod_396[k];

        t_663[k] = f_3 * pc_x[k] * sod_399[k];

        t_664[k] = f_3 * pc_x[k] * sod_400[k];

        t_665[k] = f_3 * pc_x[k] * sod_401[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, snd_333, snd_335, sop0_199, \
                         sop0_200, sop1_199, sop1_200, sod_399, \
                         sod_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_0 * snd_333[k]
                   + f_1 * sop0_199[k]
                   - f_2 * sop1_199[k]
                   + f_3 * pc_y[k] * sod_399[k];

        t_667[k] = f_3 * pc_z[k] * sod_399[k];

        t_668[k] = f_0 * snd_335[k]
                   + f_3 * pc_y[k] * sod_401[k];

        t_669[k] = f_1 * sop0_200[k]
                   - f_2 * sop1_200[k]
                   + f_3 * pc_z[k] * sod_401[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, t_674, pb_z, pc_x, pc_y, pc_z, snf0_550, \
                         snd_330, snd_336, snf1_550, sod_402, sod_405, \
                         sod_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = pb_z[k] * snf0_550[k]
                   - f_4 * pc_z[k] * snf1_550[k];

        t_671[k] = f_6 * snd_336[k]
                   + f_3 * pc_y[k] * sod_402[k];

        t_672[k] = f_5 * snd_330[k]
                   + f_3 * pc_z[k] * sod_402[k];

        t_673[k] = f_3 * pc_x[k] * sod_405[k];

        t_674[k] = f_3 * pc_x[k] * sod_406[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, t_678, pb_z, pc_x, pc_y, pc_z, snf0_556, \
                         snd_333, snd_341, snf1_556, sod_405, sod_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = f_3 * pc_x[k] * sod_407[k];

        t_676[k] = pb_z[k] * snf0_556[k]
                   - f_4 * pc_z[k] * snf1_556[k];

        t_677[k] = f_5 * snd_333[k]
                   + f_3 * pc_z[k] * sod_405[k];

        t_678[k] = f_6 * snd_341[k]
                   + f_3 * pc_y[k] * sod_407[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, t_682, pc_x, pc_y, pc_z, snd_335, snd_336, \
                         snd_342, sop0_203, sop0_204, sop1_203, sop1_204, sod_407, \
                         sod_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_5 * snd_335[k]
                   + f_1 * sop0_203[k]
                   - f_2 * sop1_203[k]
                   + f_3 * pc_z[k] * sod_407[k];

        t_680[k] = f_1 * sop0_204[k]
                   - f_2 * sop1_204[k]
                   + f_3 * pc_x[k] * sod_408[k];

        t_681[k] = f_7 * snd_342[k]
                   + f_3 * pc_y[k] * sod_408[k];

        t_682[k] = f_8 * snd_336[k]
                   + f_3 * pc_z[k] * sod_408[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, t_686, t_687, pc_x, pc_y, pc_z, snd_339, \
                         snd_345, sop0_205, sop1_205, sod_411, sod_412, \
                         sod_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_3 * pc_x[k] * sod_411[k];

        t_684[k] = f_3 * pc_x[k] * sod_412[k];

        t_685[k] = f_3 * pc_x[k] * sod_413[k];

        t_686[k] = f_7 * snd_345[k]
                   + f_1 * sop0_205[k]
                   - f_2 * sop1_205[k]
                   + f_3 * pc_y[k] * sod_411[k];

        t_687[k] = f_8 * snd_339[k]
                   + f_3 * pc_z[k] * sod_411[k];
    }

#pragma omp simd aligned(t_688, t_689, t_690, t_691, pc_x, pc_y, pc_z, snd_341, snd_347, \
                         snd_348, sop0_206, sop0_207, sop1_206, sop1_207, sod_413, \
                         sod_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_688[k] = f_7 * snd_347[k]
                   + f_3 * pc_y[k] * sod_413[k];

        t_689[k] = f_8 * snd_341[k]
                   + f_1 * sop0_206[k]
                   - f_2 * sop1_206[k]
                   + f_3 * pc_z[k] * sod_413[k];

        t_690[k] = f_1 * sop0_207[k]
                   - f_2 * sop1_207[k]
                   + f_3 * pc_x[k] * sod_414[k];

        t_691[k] = f_9 * snd_348[k]
                   + f_3 * pc_y[k] * sod_414[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, pc_x, pc_y, pc_z, snd_342, \
                         snd_351, sop0_208, sop1_208, sod_414, sod_417, sod_418, \
                         sod_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_10 * snd_342[k]
                   + f_3 * pc_z[k] * sod_414[k];

        t_693[k] = f_3 * pc_x[k] * sod_417[k];

        t_694[k] = f_3 * pc_x[k] * sod_418[k];

        t_695[k] = f_3 * pc_x[k] * sod_419[k];

        t_696[k] = f_9 * snd_351[k]
                   + f_1 * sop0_208[k]
                   - f_2 * sop1_208[k]
                   + f_3 * pc_y[k] * sod_417[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, snd_345, snd_347, snd_353, sop0_209, \
                         sop1_209, sod_417, sod_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_10 * snd_345[k]
                   + f_3 * pc_z[k] * sod_417[k];

        t_698[k] = f_9 * snd_353[k]
                   + f_3 * pc_y[k] * sod_419[k];

        t_699[k] = f_10 * snd_347[k]
                   + f_1 * sop0_209[k]
                   - f_2 * sop1_209[k]
                   + f_3 * pc_z[k] * sod_419[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, t_704, pc_x, pc_y, pc_z, snd_348, \
                         snd_354, sop0_210, sop1_210, sod_420, sod_423, \
                         sod_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_1 * sop0_210[k]
                   - f_2 * sop1_210[k]
                   + f_3 * pc_x[k] * sod_420[k];

        t_701[k] = f_11 * snd_354[k]
                   + f_3 * pc_y[k] * sod_420[k];

        t_702[k] = f_12 * snd_348[k]
                   + f_3 * pc_z[k] * sod_420[k];

        t_703[k] = f_3 * pc_x[k] * sod_423[k];

        t_704[k] = f_3 * pc_x[k] * sod_424[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, pc_x, pc_y, pc_z, snd_351, snd_357, \
                         snd_359, sop0_211, sop1_211, sod_423, \
                         sod_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_3 * pc_x[k] * sod_425[k];

        t_706[k] = f_11 * snd_357[k]
                   + f_1 * sop0_211[k]
                   - f_2 * sop1_211[k]
                   + f_3 * pc_y[k] * sod_423[k];

        t_707[k] = f_12 * snd_351[k]
                   + f_3 * pc_z[k] * sod_423[k];

        t_708[k] = f_11 * snd_359[k]
                   + f_3 * pc_y[k] * sod_425[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, t_712, pc_x, pc_y, pc_z, snd_353, snd_354, \
                         snd_360, sop0_212, sop0_213, sop1_212, sop1_213, sod_425, \
                         sod_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_12 * snd_353[k]
                   + f_1 * sop0_212[k]
                   - f_2 * sop1_212[k]
                   + f_3 * pc_z[k] * sod_425[k];

        t_710[k] = f_1 * sop0_213[k]
                   - f_2 * sop1_213[k]
                   + f_3 * pc_x[k] * sod_426[k];

        t_711[k] = f_13 * snd_360[k]
                   + f_3 * pc_y[k] * sod_426[k];

        t_712[k] = f_14 * snd_354[k]
                   + f_3 * pc_z[k] * sod_426[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, t_716, t_717, pc_x, pc_y, pc_z, snd_357, \
                         snd_363, sop0_214, sop1_214, sod_429, sod_430, \
                         sod_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_3 * pc_x[k] * sod_429[k];

        t_714[k] = f_3 * pc_x[k] * sod_430[k];

        t_715[k] = f_3 * pc_x[k] * sod_431[k];

        t_716[k] = f_13 * snd_363[k]
                   + f_1 * sop0_214[k]
                   - f_2 * sop1_214[k]
                   + f_3 * pc_y[k] * sod_429[k];

        t_717[k] = f_14 * snd_357[k]
                   + f_3 * pc_z[k] * sod_429[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pc_x, pc_y, pc_z, snd_359, snd_365, \
                         snd_366, sop0_215, sop0_216, sop1_215, sop1_216, sod_431, \
                         sod_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_13 * snd_365[k]
                   + f_3 * pc_y[k] * sod_431[k];

        t_719[k] = f_14 * snd_359[k]
                   + f_1 * sop0_215[k]
                   - f_2 * sop1_215[k]
                   + f_3 * pc_z[k] * sod_431[k];

        t_720[k] = f_1 * sop0_216[k]
                   - f_2 * sop1_216[k]
                   + f_3 * pc_x[k] * sod_432[k];

        t_721[k] = f_14 * snd_366[k]
                   + f_3 * pc_y[k] * sod_432[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, pc_x, pc_y, pc_z, snd_360, \
                         snd_369, sop0_217, sop1_217, sod_432, sod_435, sod_436, \
                         sod_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_13 * snd_360[k]
                   + f_3 * pc_z[k] * sod_432[k];

        t_723[k] = f_3 * pc_x[k] * sod_435[k];

        t_724[k] = f_3 * pc_x[k] * sod_436[k];

        t_725[k] = f_3 * pc_x[k] * sod_437[k];

        t_726[k] = f_14 * snd_369[k]
                   + f_1 * sop0_217[k]
                   - f_2 * sop1_217[k]
                   + f_3 * pc_y[k] * sod_435[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, pc_y, pc_z, snd_363, snd_365, snd_371, sop0_218, \
                         sop1_218, sod_435, sod_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_13 * snd_363[k]
                   + f_3 * pc_z[k] * sod_435[k];

        t_728[k] = f_14 * snd_371[k]
                   + f_3 * pc_y[k] * sod_437[k];

        t_729[k] = f_13 * snd_365[k]
                   + f_1 * sop0_218[k]
                   - f_2 * sop1_218[k]
                   + f_3 * pc_z[k] * sod_437[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, t_734, pc_x, pc_y, pc_z, snd_366, \
                         snd_372, sop0_219, sop1_219, sod_438, sod_441, \
                         sod_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_1 * sop0_219[k]
                   - f_2 * sop1_219[k]
                   + f_3 * pc_x[k] * sod_438[k];

        t_731[k] = f_12 * snd_372[k]
                   + f_3 * pc_y[k] * sod_438[k];

        t_732[k] = f_11 * snd_366[k]
                   + f_3 * pc_z[k] * sod_438[k];

        t_733[k] = f_3 * pc_x[k] * sod_441[k];

        t_734[k] = f_3 * pc_x[k] * sod_442[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pc_x, pc_y, pc_z, snd_369, snd_375, \
                         snd_377, sop0_220, sop1_220, sod_441, \
                         sod_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_3 * pc_x[k] * sod_443[k];

        t_736[k] = f_12 * snd_375[k]
                   + f_1 * sop0_220[k]
                   - f_2 * sop1_220[k]
                   + f_3 * pc_y[k] * sod_441[k];

        t_737[k] = f_11 * snd_369[k]
                   + f_3 * pc_z[k] * sod_441[k];

        t_738[k] = f_12 * snd_377[k]
                   + f_3 * pc_y[k] * sod_443[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, t_742, pc_x, pc_y, pc_z, snd_371, snd_372, \
                         snd_378, sop0_221, sop0_222, sop1_221, sop1_222, sod_443, \
                         sod_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_11 * snd_371[k]
                   + f_1 * sop0_221[k]
                   - f_2 * sop1_221[k]
                   + f_3 * pc_z[k] * sod_443[k];

        t_740[k] = f_1 * sop0_222[k]
                   - f_2 * sop1_222[k]
                   + f_3 * pc_x[k] * sod_444[k];

        t_741[k] = f_10 * snd_378[k]
                   + f_3 * pc_y[k] * sod_444[k];

        t_742[k] = f_9 * snd_372[k]
                   + f_3 * pc_z[k] * sod_444[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, t_746, t_747, pc_x, pc_y, pc_z, snd_375, \
                         snd_381, sop0_223, sop1_223, sod_447, sod_448, \
                         sod_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = f_3 * pc_x[k] * sod_447[k];

        t_744[k] = f_3 * pc_x[k] * sod_448[k];

        t_745[k] = f_3 * pc_x[k] * sod_449[k];

        t_746[k] = f_10 * snd_381[k]
                   + f_1 * sop0_223[k]
                   - f_2 * sop1_223[k]
                   + f_3 * pc_y[k] * sod_447[k];

        t_747[k] = f_9 * snd_375[k]
                   + f_3 * pc_z[k] * sod_447[k];
    }

#pragma omp simd aligned(t_748, t_749, t_750, t_751, pc_x, pc_y, pc_z, snd_377, snd_383, \
                         snd_384, sop0_224, sop0_225, sop1_224, sop1_225, sod_449, \
                         sod_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_748[k] = f_10 * snd_383[k]
                   + f_3 * pc_y[k] * sod_449[k];

        t_749[k] = f_9 * snd_377[k]
                   + f_1 * sop0_224[k]
                   - f_2 * sop1_224[k]
                   + f_3 * pc_z[k] * sod_449[k];

        t_750[k] = f_1 * sop0_225[k]
                   - f_2 * sop1_225[k]
                   + f_3 * pc_x[k] * sod_450[k];

        t_751[k] = f_8 * snd_384[k]
                   + f_3 * pc_y[k] * sod_450[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, t_756, pc_x, pc_y, pc_z, snd_378, \
                         snd_387, sop0_226, sop1_226, sod_450, sod_453, sod_454, \
                         sod_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_7 * snd_378[k]
                   + f_3 * pc_z[k] * sod_450[k];

        t_753[k] = f_3 * pc_x[k] * sod_453[k];

        t_754[k] = f_3 * pc_x[k] * sod_454[k];

        t_755[k] = f_3 * pc_x[k] * sod_455[k];

        t_756[k] = f_8 * snd_387[k]
                   + f_1 * sop0_226[k]
                   - f_2 * sop1_226[k]
                   + f_3 * pc_y[k] * sod_453[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pb_y, pc_y, pc_z, snf0_650, snd_381, \
                         snd_383, snd_389, snf1_650, sop0_227, sop1_227, sod_453, \
                         sod_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_7 * snd_381[k]
                   + f_3 * pc_z[k] * sod_453[k];

        t_758[k] = f_8 * snd_389[k]
                   + f_3 * pc_y[k] * sod_455[k];

        t_759[k] = f_7 * snd_383[k]
                   + f_1 * sop0_227[k]
                   - f_2 * sop1_227[k]
                   + f_3 * pc_z[k] * sod_455[k];

        t_760[k] = pb_y[k] * snf0_650[k]
                   - f_4 * pc_y[k] * snf1_650[k];
    }
}

static auto
compute_prim_sof_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t snf0,
                                                          const size_t snd, const size_t snf1,
                                                          const size_t sop0, const size_t sop1,
                                                          const size_t sod, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_10 = 1.5 / q;

    auto *t_761 = buffer.data(target + 761);
    auto *t_762 = buffer.data(target + 762);
    auto *t_763 = buffer.data(target + 763);
    auto *t_764 = buffer.data(target + 764);
    auto *t_765 = buffer.data(target + 765);
    auto *t_766 = buffer.data(target + 766);
    auto *t_767 = buffer.data(target + 767);
    auto *t_768 = buffer.data(target + 768);
    auto *t_769 = buffer.data(target + 769);
    auto *t_770 = buffer.data(target + 770);
    auto *t_771 = buffer.data(target + 771);
    auto *t_772 = buffer.data(target + 772);
    auto *t_773 = buffer.data(target + 773);
    auto *t_774 = buffer.data(target + 774);
    auto *t_775 = buffer.data(target + 775);
    auto *t_776 = buffer.data(target + 776);
    auto *t_777 = buffer.data(target + 777);
    auto *t_778 = buffer.data(target + 778);
    auto *t_779 = buffer.data(target + 779);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *snf0_656 = buffer.data(snf0 + 656);
    const auto *snf0_659 = buffer.data(snf0 + 659);

    const auto *snd_384 = buffer.data(snd + 384);
    const auto *snd_387 = buffer.data(snd + 387);
    const auto *snd_390 = buffer.data(snd + 390);
    const auto *snd_393 = buffer.data(snd + 393);
    const auto *snd_395 = buffer.data(snd + 395);

    const auto *snf1_656 = buffer.data(snf1 + 656);
    const auto *snf1_659 = buffer.data(snf1 + 659);

    const auto *sop0_231 = buffer.data(sop0 + 231);
    const auto *sop0_232 = buffer.data(sop0 + 232);
    const auto *sop0_233 = buffer.data(sop0 + 233);

    const auto *sop1_231 = buffer.data(sop1 + 231);
    const auto *sop1_232 = buffer.data(sop1 + 232);
    const auto *sop1_233 = buffer.data(sop1 + 233);

    const auto *sod_456 = buffer.data(sod + 456);
    const auto *sod_459 = buffer.data(sod + 459);
    const auto *sod_460 = buffer.data(sod + 460);
    const auto *sod_461 = buffer.data(sod + 461);
    const auto *sod_462 = buffer.data(sod + 462);
    const auto *sod_465 = buffer.data(sod + 465);
    const auto *sod_466 = buffer.data(sod + 466);
    const auto *sod_467 = buffer.data(sod + 467);

#pragma omp simd aligned(t_761, t_762, t_763, t_764, t_765, pc_x, pc_y, pc_z, snd_384, \
                         snd_390, sod_456, sod_459, sod_460, sod_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_5 * snd_390[k]
                   + f_3 * pc_y[k] * sod_456[k];

        t_762[k] = f_6 * snd_384[k]
                   + f_3 * pc_z[k] * sod_456[k];

        t_763[k] = f_3 * pc_x[k] * sod_459[k];

        t_764[k] = f_3 * pc_x[k] * sod_460[k];

        t_765[k] = f_3 * pc_x[k] * sod_461[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pb_y, pc_y, pc_z, snf0_656, snf0_659, \
                         snd_387, snd_393, snd_395, snf1_656, snf1_659, sod_459, \
                         sod_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = pb_y[k] * snf0_656[k]
                   + f_10 * snd_393[k]
                   - f_4 * pc_y[k] * snf1_656[k];

        t_767[k] = f_6 * snd_387[k]
                   + f_3 * pc_z[k] * sod_459[k];

        t_768[k] = f_5 * snd_395[k]
                   + f_3 * pc_y[k] * sod_461[k];

        t_769[k] = pb_y[k] * snf0_659[k]
                   - f_4 * pc_y[k] * snf1_659[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, t_775, pc_x, pc_y, pc_z, snd_390, \
                         sop0_231, sop1_231, sod_462, sod_465, sod_466, \
                         sod_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_1 * sop0_231[k]
                   - f_2 * sop1_231[k]
                   + f_3 * pc_x[k] * sod_462[k];

        t_771[k] = f_3 * pc_y[k] * sod_462[k];

        t_772[k] = f_0 * snd_390[k]
                   + f_3 * pc_z[k] * sod_462[k];

        t_773[k] = f_3 * pc_x[k] * sod_465[k];

        t_774[k] = f_3 * pc_x[k] * sod_466[k];

        t_775[k] = f_3 * pc_x[k] * sod_467[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, t_779, pc_y, pc_z, snd_393, snd_395, sop0_232, \
                         sop0_233, sop1_232, sop1_233, sod_465, \
                         sod_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_1 * sop0_232[k]
                   - f_2 * sop1_232[k]
                   + f_3 * pc_y[k] * sod_465[k];

        t_777[k] = f_0 * snd_393[k]
                   + f_3 * pc_z[k] * sod_465[k];

        t_778[k] = f_3 * pc_y[k] * sod_467[k];

        t_779[k] = f_0 * snd_395[k]
                   + f_1 * sop0_233[k]
                   - f_2 * sop1_233[k]
                   + f_3 * pc_z[k] * sod_467[k];
    }
}

auto
compute_prim_sof_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t snf0, const size_t snd,
                                                   const size_t snf1, const size_t sop0,
                                                   const size_t sop1, const size_t sod,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sof_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);

    compute_prim_sof_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, snf0, snd,
                                                              snf1, sop0, sop1, sod, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
