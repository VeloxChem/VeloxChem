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


#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_0 = buffer.data(ksh0 + 0);
    const auto *ksh0_3 = buffer.data(ksh0 + 3);
    const auto *ksh0_5 = buffer.data(ksh0 + 5);
    const auto *ksh0_6 = buffer.data(ksh0 + 6);
    const auto *ksh0_9 = buffer.data(ksh0 + 9);
    const auto *ksh0_15 = buffer.data(ksh0 + 15);
    const auto *ksh0_20 = buffer.data(ksh0 + 20);
    const auto *ksh0_24 = buffer.data(ksh0 + 24);
    const auto *ksh0_27 = buffer.data(ksh0 + 27);
    const auto *ksh0_36 = buffer.data(ksh0 + 36);
    const auto *ksh0_42 = buffer.data(ksh0 + 42);
    const auto *ksh0_47 = buffer.data(ksh0 + 47);
    const auto *ksh0_51 = buffer.data(ksh0 + 51);
    const auto *ksh0_62 = buffer.data(ksh0 + 62);

    const auto *ksg_0 = buffer.data(ksg + 0);
    const auto *ksg_1 = buffer.data(ksg + 1);
    const auto *ksg_2 = buffer.data(ksg + 2);
    const auto *ksg_3 = buffer.data(ksg + 3);
    const auto *ksg_5 = buffer.data(ksg + 5);
    const auto *ksg_10 = buffer.data(ksg + 10);
    const auto *ksg_12 = buffer.data(ksg + 12);
    const auto *ksg_14 = buffer.data(ksg + 14);
    const auto *ksg_15 = buffer.data(ksg + 15);
    const auto *ksg_18 = buffer.data(ksg + 18);
    const auto *ksg_20 = buffer.data(ksg + 20);
    const auto *ksg_25 = buffer.data(ksg + 25);
    const auto *ksg_27 = buffer.data(ksg + 27);
    const auto *ksg_28 = buffer.data(ksg + 28);
    const auto *ksg_29 = buffer.data(ksg + 29);
    const auto *ksg_30 = buffer.data(ksg + 30);
    const auto *ksg_32 = buffer.data(ksg + 32);
    const auto *ksg_35 = buffer.data(ksg + 35);
    const auto *ksg_40 = buffer.data(ksg + 40);
    const auto *ksg_41 = buffer.data(ksg + 41);
    const auto *ksg_42 = buffer.data(ksg + 42);
    const auto *ksg_43 = buffer.data(ksg + 43);
    const auto *ksg_44 = buffer.data(ksg + 44);
    const auto *ksg_45 = buffer.data(ksg + 45);
    const auto *ksg_48 = buffer.data(ksg + 48);
    const auto *ksg_51 = buffer.data(ksg + 51);
    const auto *ksg_55 = buffer.data(ksg + 55);
    const auto *ksg_57 = buffer.data(ksg + 57);
    const auto *ksg_58 = buffer.data(ksg + 58);
    const auto *ksg_59 = buffer.data(ksg + 59);
    const auto *ksg_70 = buffer.data(ksg + 70);
    const auto *ksg_71 = buffer.data(ksg + 71);
    const auto *ksg_72 = buffer.data(ksg + 72);
    const auto *ksg_73 = buffer.data(ksg + 73);
    const auto *ksg_74 = buffer.data(ksg + 74);
    const auto *ksg_75 = buffer.data(ksg + 75);
    const auto *ksg_80 = buffer.data(ksg + 80);
    const auto *ksg_84 = buffer.data(ksg + 84);
    const auto *ksg_85 = buffer.data(ksg + 85);
    const auto *ksg_86 = buffer.data(ksg + 86);
    const auto *ksg_87 = buffer.data(ksg + 87);
    const auto *ksg_89 = buffer.data(ksg + 89);
    const auto *ksg_90 = buffer.data(ksg + 90);
    const auto *ksg_93 = buffer.data(ksg + 93);

    const auto *ksh1_0 = buffer.data(ksh1 + 0);
    const auto *ksh1_3 = buffer.data(ksh1 + 3);
    const auto *ksh1_5 = buffer.data(ksh1 + 5);
    const auto *ksh1_6 = buffer.data(ksh1 + 6);
    const auto *ksh1_9 = buffer.data(ksh1 + 9);
    const auto *ksh1_15 = buffer.data(ksh1 + 15);
    const auto *ksh1_20 = buffer.data(ksh1 + 20);
    const auto *ksh1_24 = buffer.data(ksh1 + 24);
    const auto *ksh1_27 = buffer.data(ksh1 + 27);
    const auto *ksh1_36 = buffer.data(ksh1 + 36);
    const auto *ksh1_42 = buffer.data(ksh1 + 42);
    const auto *ksh1_47 = buffer.data(ksh1 + 47);
    const auto *ksh1_51 = buffer.data(ksh1 + 51);
    const auto *ksh1_62 = buffer.data(ksh1 + 62);

    const auto *lsf0_0 = buffer.data(lsf0 + 0);
    const auto *lsf0_1 = buffer.data(lsf0 + 1);
    const auto *lsf0_2 = buffer.data(lsf0 + 2);
    const auto *lsf0_6 = buffer.data(lsf0 + 6);
    const auto *lsf0_8 = buffer.data(lsf0 + 8);
    const auto *lsf0_9 = buffer.data(lsf0 + 9);
    const auto *lsf0_16 = buffer.data(lsf0 + 16);
    const auto *lsf0_17 = buffer.data(lsf0 + 17);
    const auto *lsf0_22 = buffer.data(lsf0 + 22);
    const auto *lsf0_27 = buffer.data(lsf0 + 27);
    const auto *lsf0_28 = buffer.data(lsf0 + 28);
    const auto *lsf0_29 = buffer.data(lsf0 + 29);
    const auto *lsf0_30 = buffer.data(lsf0 + 30);
    const auto *lsf0_32 = buffer.data(lsf0 + 32);
    const auto *lsf0_33 = buffer.data(lsf0 + 33);
    const auto *lsf0_36 = buffer.data(lsf0 + 36);
    const auto *lsf0_37 = buffer.data(lsf0 + 37);
    const auto *lsf0_39 = buffer.data(lsf0 + 39);
    const auto *lsf0_48 = buffer.data(lsf0 + 48);
    const auto *lsf0_49 = buffer.data(lsf0 + 49);
    const auto *lsf0_50 = buffer.data(lsf0 + 50);
    const auto *lsf0_51 = buffer.data(lsf0 + 51);
    const auto *lsf0_52 = buffer.data(lsf0 + 52);
    const auto *lsf0_55 = buffer.data(lsf0 + 55);
    const auto *lsf0_56 = buffer.data(lsf0 + 56);
    const auto *lsf0_57 = buffer.data(lsf0 + 57);
    const auto *lsf0_58 = buffer.data(lsf0 + 58);
    const auto *lsf0_59 = buffer.data(lsf0 + 59);
    const auto *lsf0_60 = buffer.data(lsf0 + 60);
    const auto *lsf0_63 = buffer.data(lsf0 + 63);

    const auto *lsf1_0 = buffer.data(lsf1 + 0);
    const auto *lsf1_1 = buffer.data(lsf1 + 1);
    const auto *lsf1_2 = buffer.data(lsf1 + 2);
    const auto *lsf1_6 = buffer.data(lsf1 + 6);
    const auto *lsf1_8 = buffer.data(lsf1 + 8);
    const auto *lsf1_9 = buffer.data(lsf1 + 9);
    const auto *lsf1_16 = buffer.data(lsf1 + 16);
    const auto *lsf1_17 = buffer.data(lsf1 + 17);
    const auto *lsf1_22 = buffer.data(lsf1 + 22);
    const auto *lsf1_27 = buffer.data(lsf1 + 27);
    const auto *lsf1_28 = buffer.data(lsf1 + 28);
    const auto *lsf1_29 = buffer.data(lsf1 + 29);
    const auto *lsf1_30 = buffer.data(lsf1 + 30);
    const auto *lsf1_32 = buffer.data(lsf1 + 32);
    const auto *lsf1_33 = buffer.data(lsf1 + 33);
    const auto *lsf1_36 = buffer.data(lsf1 + 36);
    const auto *lsf1_37 = buffer.data(lsf1 + 37);
    const auto *lsf1_39 = buffer.data(lsf1 + 39);
    const auto *lsf1_48 = buffer.data(lsf1 + 48);
    const auto *lsf1_49 = buffer.data(lsf1 + 49);
    const auto *lsf1_50 = buffer.data(lsf1 + 50);
    const auto *lsf1_51 = buffer.data(lsf1 + 51);
    const auto *lsf1_52 = buffer.data(lsf1 + 52);
    const auto *lsf1_55 = buffer.data(lsf1 + 55);
    const auto *lsf1_56 = buffer.data(lsf1 + 56);
    const auto *lsf1_57 = buffer.data(lsf1 + 57);
    const auto *lsf1_58 = buffer.data(lsf1 + 58);
    const auto *lsf1_59 = buffer.data(lsf1 + 59);
    const auto *lsf1_60 = buffer.data(lsf1 + 60);
    const auto *lsf1_63 = buffer.data(lsf1 + 63);

    const auto *lsg_0 = buffer.data(lsg + 0);
    const auto *lsg_1 = buffer.data(lsg + 1);
    const auto *lsg_2 = buffer.data(lsg + 2);
    const auto *lsg_3 = buffer.data(lsg + 3);
    const auto *lsg_5 = buffer.data(lsg + 5);
    const auto *lsg_6 = buffer.data(lsg + 6);
    const auto *lsg_9 = buffer.data(lsg + 9);
    const auto *lsg_10 = buffer.data(lsg + 10);
    const auto *lsg_12 = buffer.data(lsg + 12);
    const auto *lsg_13 = buffer.data(lsg + 13);
    const auto *lsg_14 = buffer.data(lsg + 14);
    const auto *lsg_15 = buffer.data(lsg + 15);
    const auto *lsg_16 = buffer.data(lsg + 16);
    const auto *lsg_18 = buffer.data(lsg + 18);
    const auto *lsg_20 = buffer.data(lsg + 20);
    const auto *lsg_21 = buffer.data(lsg + 21);
    const auto *lsg_25 = buffer.data(lsg + 25);
    const auto *lsg_26 = buffer.data(lsg + 26);
    const auto *lsg_27 = buffer.data(lsg + 27);
    const auto *lsg_28 = buffer.data(lsg + 28);
    const auto *lsg_29 = buffer.data(lsg + 29);
    const auto *lsg_30 = buffer.data(lsg + 30);
    const auto *lsg_32 = buffer.data(lsg + 32);
    const auto *lsg_34 = buffer.data(lsg + 34);
    const auto *lsg_35 = buffer.data(lsg + 35);
    const auto *lsg_39 = buffer.data(lsg + 39);
    const auto *lsg_40 = buffer.data(lsg + 40);
    const auto *lsg_41 = buffer.data(lsg + 41);
    const auto *lsg_42 = buffer.data(lsg + 42);
    const auto *lsg_43 = buffer.data(lsg + 43);
    const auto *lsg_44 = buffer.data(lsg + 44);
    const auto *lsg_45 = buffer.data(lsg + 45);
    const auto *lsg_46 = buffer.data(lsg + 46);
    const auto *lsg_47 = buffer.data(lsg + 47);
    const auto *lsg_48 = buffer.data(lsg + 48);
    const auto *lsg_50 = buffer.data(lsg + 50);
    const auto *lsg_51 = buffer.data(lsg + 51);
    const auto *lsg_55 = buffer.data(lsg + 55);
    const auto *lsg_56 = buffer.data(lsg + 56);
    const auto *lsg_57 = buffer.data(lsg + 57);
    const auto *lsg_58 = buffer.data(lsg + 58);
    const auto *lsg_59 = buffer.data(lsg + 59);
    const auto *lsg_60 = buffer.data(lsg + 60);
    const auto *lsg_62 = buffer.data(lsg + 62);
    const auto *lsg_63 = buffer.data(lsg + 63);
    const auto *lsg_65 = buffer.data(lsg + 65);
    const auto *lsg_70 = buffer.data(lsg + 70);
    const auto *lsg_71 = buffer.data(lsg + 71);
    const auto *lsg_72 = buffer.data(lsg + 72);
    const auto *lsg_73 = buffer.data(lsg + 73);
    const auto *lsg_74 = buffer.data(lsg + 74);
    const auto *lsg_75 = buffer.data(lsg + 75);
    const auto *lsg_76 = buffer.data(lsg + 76);
    const auto *lsg_77 = buffer.data(lsg + 77);
    const auto *lsg_78 = buffer.data(lsg + 78);
    const auto *lsg_79 = buffer.data(lsg + 79);
    const auto *lsg_80 = buffer.data(lsg + 80);
    const auto *lsg_84 = buffer.data(lsg + 84);
    const auto *lsg_85 = buffer.data(lsg + 85);
    const auto *lsg_86 = buffer.data(lsg + 86);
    const auto *lsg_87 = buffer.data(lsg + 87);
    const auto *lsg_88 = buffer.data(lsg + 88);
    const auto *lsg_89 = buffer.data(lsg + 89);
    const auto *lsg_90 = buffer.data(lsg + 90);
    const auto *lsg_93 = buffer.data(lsg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksg_0, lsf0_0, \
                         lsf1_0, lsg_0, lsg_1, lsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksg_0[k]
                 + f_1 * lsf0_0[k]
                 - f_2 * lsf1_0[k]
                 + f_3 * pc_x[k] * lsg_0[k];

        t_1[k] = f_3 * pc_y[k] * lsg_0[k];

        t_2[k] = f_3 * pc_z[k] * lsg_0[k];

        t_3[k] = f_4 * lsf0_0[k]
                 - f_5 * lsf1_0[k]
                 + f_3 * pc_y[k] * lsg_1[k];

        t_4[k] = f_3 * pc_y[k] * lsg_2[k];

        t_5[k] = f_4 * lsf0_0[k]
                 - f_5 * lsf1_0[k]
                 + f_3 * pc_z[k] * lsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, ksg_10, lsf0_1, lsf0_2, \
                         lsf1_1, lsf1_2, lsg_3, lsg_5, lsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * lsf0_1[k]
                 - f_7 * lsf1_1[k]
                 + f_3 * pc_y[k] * lsg_3[k];

        t_7[k] = f_3 * pc_z[k] * lsg_3[k];

        t_8[k] = f_3 * pc_y[k] * lsg_5[k];

        t_9[k] = f_6 * lsf0_2[k]
                 - f_7 * lsf1_2[k]
                 + f_3 * pc_z[k] * lsg_5[k];

        t_10[k] = f_0 * ksg_10[k]
                  + f_3 * pc_x[k] * lsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, ksg_12, ksg_14, lsg_6, \
                         lsg_9, lsg_12, lsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * lsg_6[k];

        t_12[k] = f_0 * ksg_12[k]
                  + f_3 * pc_x[k] * lsg_12[k];

        t_13[k] = f_3 * pc_y[k] * lsg_9[k];

        t_14[k] = f_0 * ksg_14[k]
                  + f_3 * pc_x[k] * lsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, lsf0_6, lsf0_8, lsf0_9, lsf1_6, \
                         lsf1_8, lsf1_9, lsg_10, lsg_12, lsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * lsf0_6[k]
                  - f_2 * lsf1_6[k]
                  + f_3 * pc_y[k] * lsg_10[k];

        t_16[k] = f_3 * pc_z[k] * lsg_10[k];

        t_17[k] = f_6 * lsf0_8[k]
                  - f_7 * lsf1_8[k]
                  + f_3 * pc_y[k] * lsg_12[k];

        t_18[k] = f_4 * lsf0_9[k]
                  - f_5 * lsf1_9[k]
                  + f_3 * pc_y[k] * lsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, ksh0_0, ksg_0, \
                         ksh1_0, lsf0_9, lsf1_9, lsg_14, lsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * lsg_14[k];

        t_20[k] = f_1 * lsf0_9[k]
                  - f_2 * lsf1_9[k]
                  + f_3 * pc_z[k] * lsg_14[k];

        t_21[k] = pa_y[k] * ksh0_0[k]
                  - f_8 * pc_y[k] * ksh1_0[k];

        t_22[k] = f_9 * ksg_0[k]
                  + f_3 * pc_y[k] * lsg_15[k];

        t_23[k] = f_3 * pc_z[k] * lsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, ksh0_3, ksh0_5, ksh0_6, \
                         ksg_1, ksg_3, ksh1_3, ksh1_5, ksh1_6, lsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * ksh0_3[k]
                  + f_10 * ksg_1[k]
                  - f_8 * pc_y[k] * ksh1_3[k];

        t_25[k] = f_3 * pc_z[k] * lsg_16[k];

        t_26[k] = pa_y[k] * ksh0_5[k]
                  - f_8 * pc_y[k] * ksh1_5[k];

        t_27[k] = pa_y[k] * ksh0_6[k]
                  + f_11 * ksg_3[k]
                  - f_8 * pc_y[k] * ksh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, ksh0_9, ksg_5, \
                         ksg_25, ksh1_9, lsg_18, lsg_20, lsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * lsg_18[k];

        t_29[k] = f_9 * ksg_5[k]
                  + f_3 * pc_y[k] * lsg_20[k];

        t_30[k] = pa_y[k] * ksh0_9[k]
                  - f_8 * pc_y[k] * ksh1_9[k];

        t_31[k] = f_12 * ksg_25[k]
                  + f_3 * pc_x[k] * lsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, ksg_27, ksg_28, ksg_29, lsg_21, \
                         lsg_27, lsg_28, lsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * lsg_21[k];

        t_33[k] = f_12 * ksg_27[k]
                  + f_3 * pc_x[k] * lsg_27[k];

        t_34[k] = f_12 * ksg_28[k]
                  + f_3 * pc_x[k] * lsg_28[k];

        t_35[k] = f_12 * ksg_29[k]
                  + f_3 * pc_x[k] * lsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, ksg_10, lsf0_16, lsf0_17, \
                         lsf1_16, lsf1_17, lsg_25, lsg_26, lsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * ksg_10[k]
                  + f_1 * lsf0_16[k]
                  - f_2 * lsf1_16[k]
                  + f_3 * pc_y[k] * lsg_25[k];

        t_37[k] = f_3 * pc_z[k] * lsg_25[k];

        t_38[k] = f_4 * lsf0_16[k]
                  - f_5 * lsf1_16[k]
                  + f_3 * pc_z[k] * lsg_26[k];

        t_39[k] = f_6 * lsf0_17[k]
                  - f_7 * lsf1_17[k]
                  + f_3 * pc_z[k] * lsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, ksh0_0, ksh0_20, \
                         ksg_14, ksh1_0, ksh1_20, lsg_29, lsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * ksg_14[k]
                  + f_3 * pc_y[k] * lsg_29[k];

        t_41[k] = pa_y[k] * ksh0_20[k]
                  - f_8 * pc_y[k] * ksh1_20[k];

        t_42[k] = pa_z[k] * ksh0_0[k]
                  - f_8 * pc_z[k] * ksh1_0[k];

        t_43[k] = f_3 * pc_y[k] * lsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, ksh0_3, ksh0_5, ksg_0, \
                         ksg_2, ksh1_3, ksh1_5, lsg_30, lsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * ksg_0[k]
                  + f_3 * pc_z[k] * lsg_30[k];

        t_45[k] = pa_z[k] * ksh0_3[k]
                  - f_8 * pc_z[k] * ksh1_3[k];

        t_46[k] = f_3 * pc_y[k] * lsg_32[k];

        t_47[k] = pa_z[k] * ksh0_5[k]
                  + f_10 * ksg_2[k]
                  - f_8 * pc_z[k] * ksh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, ksh0_6, ksh0_9, ksg_5, \
                         ksh1_6, ksh1_9, lsf0_22, lsf1_22, lsg_34, \
                         lsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * ksh0_6[k]
                  - f_8 * pc_z[k] * ksh1_6[k];

        t_49[k] = f_4 * lsf0_22[k]
                  - f_5 * lsf1_22[k]
                  + f_3 * pc_y[k] * lsg_34[k];

        t_50[k] = f_3 * pc_y[k] * lsg_35[k];

        t_51[k] = pa_z[k] * ksh0_9[k]
                  + f_11 * ksg_5[k]
                  - f_8 * pc_z[k] * ksh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, ksg_40, ksg_41, ksg_42, \
                         ksg_44, lsg_39, lsg_40, lsg_41, lsg_42, \
                         lsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * ksg_40[k]
                  + f_3 * pc_x[k] * lsg_40[k];

        t_53[k] = f_12 * ksg_41[k]
                  + f_3 * pc_x[k] * lsg_41[k];

        t_54[k] = f_12 * ksg_42[k]
                  + f_3 * pc_x[k] * lsg_42[k];

        t_55[k] = f_3 * pc_y[k] * lsg_39[k];

        t_56[k] = f_12 * ksg_44[k]
                  + f_3 * pc_x[k] * lsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, ksh0_15, ksh1_15, lsf0_27, \
                         lsf0_28, lsf1_27, lsf1_28, lsg_41, lsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * ksh0_15[k]
                  - f_8 * pc_z[k] * ksh1_15[k];

        t_58[k] = f_13 * lsf0_27[k]
                  - f_14 * lsf1_27[k]
                  + f_3 * pc_y[k] * lsg_41[k];

        t_59[k] = f_6 * lsf0_28[k]
                  - f_7 * lsf1_28[k]
                  + f_3 * pc_y[k] * lsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, ksg_14, ksg_45, lsf0_29, \
                         lsf0_30, lsf1_29, lsf1_30, lsg_43, lsg_44, \
                         lsg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * lsf0_29[k]
                  - f_5 * lsf1_29[k]
                  + f_3 * pc_y[k] * lsg_43[k];

        t_61[k] = f_3 * pc_y[k] * lsg_44[k];

        t_62[k] = f_9 * ksg_14[k]
                  + f_1 * lsf0_29[k]
                  - f_2 * lsf1_29[k]
                  + f_3 * pc_z[k] * lsg_44[k];

        t_63[k] = f_15 * ksg_45[k]
                  + f_1 * lsf0_30[k]
                  - f_2 * lsf1_30[k]
                  + f_3 * pc_x[k] * lsg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, ksg_15, ksg_48, lsf0_33, \
                         lsf1_33, lsg_45, lsg_46, lsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * ksg_15[k]
                  + f_3 * pc_y[k] * lsg_45[k];

        t_65[k] = f_3 * pc_z[k] * lsg_45[k];

        t_66[k] = f_15 * ksg_48[k]
                  + f_6 * lsf0_33[k]
                  - f_7 * lsf1_33[k]
                  + f_3 * pc_x[k] * lsg_48[k];

        t_67[k] = f_3 * pc_z[k] * lsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, ksg_51, lsf0_30, lsf0_36, lsf1_30, \
                         lsf1_36, lsg_47, lsg_48, lsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * lsf0_30[k]
                  - f_5 * lsf1_30[k]
                  + f_3 * pc_z[k] * lsg_47[k];

        t_69[k] = f_15 * ksg_51[k]
                  + f_4 * lsf0_36[k]
                  - f_5 * lsf1_36[k]
                  + f_3 * pc_x[k] * lsg_51[k];

        t_70[k] = f_3 * pc_z[k] * lsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, ksg_20, ksg_55, lsf0_32, \
                         lsf1_32, lsg_50, lsg_51, lsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * ksg_20[k]
                  + f_3 * pc_y[k] * lsg_50[k];

        t_72[k] = f_6 * lsf0_32[k]
                  - f_7 * lsf1_32[k]
                  + f_3 * pc_z[k] * lsg_50[k];

        t_73[k] = f_15 * ksg_55[k]
                  + f_3 * pc_x[k] * lsg_55[k];

        t_74[k] = f_3 * pc_z[k] * lsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, ksg_25, ksg_57, ksg_58, ksg_59, \
                         lsf0_36, lsf1_36, lsg_55, lsg_57, lsg_58, \
                         lsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * ksg_57[k]
                  + f_3 * pc_x[k] * lsg_57[k];

        t_76[k] = f_15 * ksg_58[k]
                  + f_3 * pc_x[k] * lsg_58[k];

        t_77[k] = f_15 * ksg_59[k]
                  + f_3 * pc_x[k] * lsg_59[k];

        t_78[k] = f_10 * ksg_25[k]
                  + f_1 * lsf0_36[k]
                  - f_2 * lsf1_36[k]
                  + f_3 * pc_y[k] * lsg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, ksg_29, lsf0_36, lsf0_37, \
                         lsf1_36, lsf1_37, lsg_55, lsg_56, lsg_57, \
                         lsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * lsg_55[k];

        t_80[k] = f_4 * lsf0_36[k]
                  - f_5 * lsf1_36[k]
                  + f_3 * pc_z[k] * lsg_56[k];

        t_81[k] = f_6 * lsf0_37[k]
                  - f_7 * lsf1_37[k]
                  + f_3 * pc_z[k] * lsg_57[k];

        t_82[k] = f_10 * ksg_29[k]
                  + f_3 * pc_y[k] * lsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, ksh0_42, ksg_15, ksg_30, \
                         ksh1_42, lsf0_39, lsf1_39, lsg_59, lsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * lsf0_39[k]
                  - f_2 * lsf1_39[k]
                  + f_3 * pc_z[k] * lsg_59[k];

        t_84[k] = pa_y[k] * ksh0_42[k]
                  - f_8 * pc_y[k] * ksh1_42[k];

        t_85[k] = f_9 * ksg_30[k]
                  + f_3 * pc_y[k] * lsg_60[k];

        t_86[k] = f_9 * ksg_15[k]
                  + f_3 * pc_z[k] * lsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, ksh0_24, ksh0_27, \
                         ksh0_47, ksg_32, ksh1_24, ksh1_27, ksh1_47, \
                         lsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * ksh0_24[k]
                  - f_8 * pc_z[k] * ksh1_24[k];

        t_88[k] = f_9 * ksg_32[k]
                  + f_3 * pc_y[k] * lsg_62[k];

        t_89[k] = pa_y[k] * ksh0_47[k]
                  - f_8 * pc_y[k] * ksh1_47[k];

        t_90[k] = pa_z[k] * ksh0_27[k]
                  - f_8 * pc_z[k] * ksh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, ksh0_51, ksg_18, \
                         ksg_35, ksg_70, ksh1_51, lsg_63, lsg_65, \
                         lsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * ksg_18[k]
                  + f_3 * pc_z[k] * lsg_63[k];

        t_92[k] = f_9 * ksg_35[k]
                  + f_3 * pc_y[k] * lsg_65[k];

        t_93[k] = pa_y[k] * ksh0_51[k]
                  - f_8 * pc_y[k] * ksh1_51[k];

        t_94[k] = f_15 * ksg_70[k]
                  + f_3 * pc_x[k] * lsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, ksg_71, ksg_72, ksg_73, ksg_74, lsg_71, \
                         lsg_72, lsg_73, lsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * ksg_71[k]
                  + f_3 * pc_x[k] * lsg_71[k];

        t_96[k] = f_15 * ksg_72[k]
                  + f_3 * pc_x[k] * lsg_72[k];

        t_97[k] = f_15 * ksg_73[k]
                  + f_3 * pc_x[k] * lsg_73[k];

        t_98[k] = f_15 * ksg_74[k]
                  + f_3 * pc_x[k] * lsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, ksh0_36, ksg_25, ksg_42, \
                         ksh1_36, lsf0_48, lsf1_48, lsg_70, lsg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * ksh0_36[k]
                  - f_8 * pc_z[k] * ksh1_36[k];

        t_100[k] = f_9 * ksg_25[k]
                   + f_3 * pc_z[k] * lsg_70[k];

        t_101[k] = f_9 * ksg_42[k]
                   + f_6 * lsf0_48[k]
                   - f_7 * lsf1_48[k]
                   + f_3 * pc_y[k] * lsg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, ksh0_62, ksg_43, ksg_44, ksh1_62, \
                         lsf0_49, lsf1_49, lsg_73, lsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * ksg_43[k]
                   + f_4 * lsf0_49[k]
                   - f_5 * lsf1_49[k]
                   + f_3 * pc_y[k] * lsg_73[k];

        t_103[k] = f_9 * ksg_44[k]
                   + f_3 * pc_y[k] * lsg_74[k];

        t_104[k] = pa_y[k] * ksh0_62[k]
                   - f_8 * pc_y[k] * ksh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, ksg_30, ksg_75, \
                         lsf0_50, lsf1_50, lsg_75, lsg_76, lsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * ksg_75[k]
                   + f_1 * lsf0_50[k]
                   - f_2 * lsf1_50[k]
                   + f_3 * pc_x[k] * lsg_75[k];

        t_106[k] = f_3 * pc_y[k] * lsg_75[k];

        t_107[k] = f_10 * ksg_30[k]
                   + f_3 * pc_z[k] * lsg_75[k];

        t_108[k] = f_4 * lsf0_50[k]
                   - f_5 * lsf1_50[k]
                   + f_3 * pc_y[k] * lsg_76[k];

        t_109[k] = f_3 * pc_y[k] * lsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, ksg_80, lsf0_51, lsf0_52, \
                         lsf0_55, lsf1_51, lsf1_52, lsf1_55, lsg_78, lsg_79, \
                         lsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * ksg_80[k]
                   + f_6 * lsf0_55[k]
                   - f_7 * lsf1_55[k]
                   + f_3 * pc_x[k] * lsg_80[k];

        t_111[k] = f_6 * lsf0_51[k]
                   - f_7 * lsf1_51[k]
                   + f_3 * pc_y[k] * lsg_78[k];

        t_112[k] = f_4 * lsf0_52[k]
                   - f_5 * lsf1_52[k]
                   + f_3 * pc_y[k] * lsg_79[k];

        t_113[k] = f_3 * pc_y[k] * lsg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, ksg_84, ksg_85, ksg_86, ksg_87, \
                         lsf0_59, lsf1_59, lsg_84, lsg_85, lsg_86, \
                         lsg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * ksg_84[k]
                   + f_4 * lsf0_59[k]
                   - f_5 * lsf1_59[k]
                   + f_3 * pc_x[k] * lsg_84[k];

        t_115[k] = f_15 * ksg_85[k]
                   + f_3 * pc_x[k] * lsg_85[k];

        t_116[k] = f_15 * ksg_86[k]
                   + f_3 * pc_x[k] * lsg_86[k];

        t_117[k] = f_15 * ksg_87[k]
                   + f_3 * pc_x[k] * lsg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, ksg_89, lsf0_56, lsf0_57, \
                         lsf1_56, lsf1_57, lsg_84, lsg_85, lsg_86, \
                         lsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * lsg_84[k];

        t_119[k] = f_15 * ksg_89[k]
                   + f_3 * pc_x[k] * lsg_89[k];

        t_120[k] = f_1 * lsf0_56[k]
                   - f_2 * lsf1_56[k]
                   + f_3 * pc_y[k] * lsg_85[k];

        t_121[k] = f_13 * lsf0_57[k]
                   - f_14 * lsf1_57[k]
                   + f_3 * pc_y[k] * lsg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, ksg_44, lsf0_58, lsf0_59, \
                         lsf1_58, lsf1_59, lsg_87, lsg_88, lsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * lsf0_58[k]
                   - f_7 * lsf1_58[k]
                   + f_3 * pc_y[k] * lsg_87[k];

        t_123[k] = f_4 * lsf0_59[k]
                   - f_5 * lsf1_59[k]
                   + f_3 * pc_y[k] * lsg_88[k];

        t_124[k] = f_3 * pc_y[k] * lsg_89[k];

        t_125[k] = f_10 * ksg_44[k]
                   + f_1 * lsf0_59[k]
                   - f_2 * lsf1_59[k]
                   + f_3 * pc_z[k] * lsg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, ksg_45, ksg_90, ksg_93, \
                         lsf0_60, lsf0_63, lsf1_60, lsf1_63, lsg_90, \
                         lsg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * ksg_90[k]
                   + f_1 * lsf0_60[k]
                   - f_2 * lsf1_60[k]
                   + f_3 * pc_x[k] * lsg_90[k];

        t_127[k] = f_11 * ksg_45[k]
                   + f_3 * pc_y[k] * lsg_90[k];

        t_128[k] = f_3 * pc_z[k] * lsg_90[k];

        t_129[k] = f_16 * ksg_93[k]
                   + f_6 * lsf0_63[k]
                   - f_7 * lsf1_63[k]
                   + f_3 * pc_x[k] * lsg_93[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_63 = buffer.data(ksh0 + 63);
    const auto *ksh0_66 = buffer.data(ksh0 + 66);
    const auto *ksh0_69 = buffer.data(ksh0 + 69);
    const auto *ksh0_78 = buffer.data(ksh0 + 78);
    const auto *ksh0_105 = buffer.data(ksh0 + 105);
    const auto *ksh0_108 = buffer.data(ksh0 + 108);
    const auto *ksh0_110 = buffer.data(ksh0 + 110);
    const auto *ksh0_111 = buffer.data(ksh0 + 111);
    const auto *ksh0_114 = buffer.data(ksh0 + 114);
    const auto *ksh0_125 = buffer.data(ksh0 + 125);
    const auto *ksh0_126 = buffer.data(ksh0 + 126);
    const auto *ksh0_129 = buffer.data(ksh0 + 129);
    const auto *ksh0_132 = buffer.data(ksh0 + 132);
    const auto *ksh0_141 = buffer.data(ksh0 + 141);

    const auto *ksg_45 = buffer.data(ksg + 45);
    const auto *ksg_48 = buffer.data(ksg + 48);
    const auto *ksg_50 = buffer.data(ksg + 50);
    const auto *ksg_55 = buffer.data(ksg + 55);
    const auto *ksg_59 = buffer.data(ksg + 59);
    const auto *ksg_60 = buffer.data(ksg + 60);
    const auto *ksg_62 = buffer.data(ksg + 62);
    const auto *ksg_63 = buffer.data(ksg + 63);
    const auto *ksg_65 = buffer.data(ksg + 65);
    const auto *ksg_70 = buffer.data(ksg + 70);
    const auto *ksg_72 = buffer.data(ksg + 72);
    const auto *ksg_73 = buffer.data(ksg + 73);
    const auto *ksg_74 = buffer.data(ksg + 74);
    const auto *ksg_75 = buffer.data(ksg + 75);
    const auto *ksg_76 = buffer.data(ksg + 76);
    const auto *ksg_77 = buffer.data(ksg + 77);
    const auto *ksg_78 = buffer.data(ksg + 78);
    const auto *ksg_80 = buffer.data(ksg + 80);
    const auto *ksg_85 = buffer.data(ksg + 85);
    const auto *ksg_87 = buffer.data(ksg + 87);
    const auto *ksg_88 = buffer.data(ksg + 88);
    const auto *ksg_89 = buffer.data(ksg + 89);
    const auto *ksg_90 = buffer.data(ksg + 90);
    const auto *ksg_93 = buffer.data(ksg + 93);
    const auto *ksg_95 = buffer.data(ksg + 95);
    const auto *ksg_96 = buffer.data(ksg + 96);
    const auto *ksg_100 = buffer.data(ksg + 100);
    const auto *ksg_102 = buffer.data(ksg + 102);
    const auto *ksg_103 = buffer.data(ksg + 103);
    const auto *ksg_104 = buffer.data(ksg + 104);
    const auto *ksg_105 = buffer.data(ksg + 105);
    const auto *ksg_107 = buffer.data(ksg + 107);
    const auto *ksg_110 = buffer.data(ksg + 110);
    const auto *ksg_114 = buffer.data(ksg + 114);
    const auto *ksg_115 = buffer.data(ksg + 115);
    const auto *ksg_116 = buffer.data(ksg + 116);
    const auto *ksg_117 = buffer.data(ksg + 117);
    const auto *ksg_118 = buffer.data(ksg + 118);
    const auto *ksg_119 = buffer.data(ksg + 119);
    const auto *ksg_130 = buffer.data(ksg + 130);
    const auto *ksg_131 = buffer.data(ksg + 131);
    const auto *ksg_132 = buffer.data(ksg + 132);
    const auto *ksg_133 = buffer.data(ksg + 133);
    const auto *ksg_134 = buffer.data(ksg + 134);
    const auto *ksg_135 = buffer.data(ksg + 135);
    const auto *ksg_140 = buffer.data(ksg + 140);
    const auto *ksg_144 = buffer.data(ksg + 144);
    const auto *ksg_145 = buffer.data(ksg + 145);
    const auto *ksg_146 = buffer.data(ksg + 146);
    const auto *ksg_147 = buffer.data(ksg + 147);
    const auto *ksg_149 = buffer.data(ksg + 149);
    const auto *ksg_150 = buffer.data(ksg + 150);
    const auto *ksg_153 = buffer.data(ksg + 153);
    const auto *ksg_156 = buffer.data(ksg + 156);
    const auto *ksg_160 = buffer.data(ksg + 160);
    const auto *ksg_162 = buffer.data(ksg + 162);
    const auto *ksg_163 = buffer.data(ksg + 163);
    const auto *ksg_164 = buffer.data(ksg + 164);
    const auto *ksg_170 = buffer.data(ksg + 170);
    const auto *ksg_174 = buffer.data(ksg + 174);
    const auto *ksg_175 = buffer.data(ksg + 175);
    const auto *ksg_176 = buffer.data(ksg + 176);
    const auto *ksg_177 = buffer.data(ksg + 177);
    const auto *ksg_178 = buffer.data(ksg + 178);
    const auto *ksg_179 = buffer.data(ksg + 179);

    const auto *ksh1_63 = buffer.data(ksh1 + 63);
    const auto *ksh1_66 = buffer.data(ksh1 + 66);
    const auto *ksh1_69 = buffer.data(ksh1 + 69);
    const auto *ksh1_78 = buffer.data(ksh1 + 78);
    const auto *ksh1_105 = buffer.data(ksh1 + 105);
    const auto *ksh1_108 = buffer.data(ksh1 + 108);
    const auto *ksh1_110 = buffer.data(ksh1 + 110);
    const auto *ksh1_111 = buffer.data(ksh1 + 111);
    const auto *ksh1_114 = buffer.data(ksh1 + 114);
    const auto *ksh1_125 = buffer.data(ksh1 + 125);
    const auto *ksh1_126 = buffer.data(ksh1 + 126);
    const auto *ksh1_129 = buffer.data(ksh1 + 129);
    const auto *ksh1_132 = buffer.data(ksh1 + 132);
    const auto *ksh1_141 = buffer.data(ksh1 + 141);

    const auto *lsf0_60 = buffer.data(lsf0 + 60);
    const auto *lsf0_62 = buffer.data(lsf0 + 62);
    const auto *lsf0_66 = buffer.data(lsf0 + 66);
    const auto *lsf0_67 = buffer.data(lsf0 + 67);
    const auto *lsf0_69 = buffer.data(lsf0 + 69);
    const auto *lsf0_75 = buffer.data(lsf0 + 75);
    const auto *lsf0_78 = buffer.data(lsf0 + 78);
    const auto *lsf0_79 = buffer.data(lsf0 + 79);
    const auto *lsf0_86 = buffer.data(lsf0 + 86);
    const auto *lsf0_88 = buffer.data(lsf0 + 88);
    const auto *lsf0_89 = buffer.data(lsf0 + 89);
    const auto *lsf0_90 = buffer.data(lsf0 + 90);
    const auto *lsf0_91 = buffer.data(lsf0 + 91);
    const auto *lsf0_92 = buffer.data(lsf0 + 92);
    const auto *lsf0_95 = buffer.data(lsf0 + 95);
    const auto *lsf0_96 = buffer.data(lsf0 + 96);
    const auto *lsf0_97 = buffer.data(lsf0 + 97);
    const auto *lsf0_98 = buffer.data(lsf0 + 98);
    const auto *lsf0_99 = buffer.data(lsf0 + 99);
    const auto *lsf0_100 = buffer.data(lsf0 + 100);
    const auto *lsf0_102 = buffer.data(lsf0 + 102);
    const auto *lsf0_103 = buffer.data(lsf0 + 103);
    const auto *lsf0_106 = buffer.data(lsf0 + 106);
    const auto *lsf0_107 = buffer.data(lsf0 + 107);
    const auto *lsf0_109 = buffer.data(lsf0 + 109);
    const auto *lsf0_115 = buffer.data(lsf0 + 115);
    const auto *lsf0_118 = buffer.data(lsf0 + 118);
    const auto *lsf0_119 = buffer.data(lsf0 + 119);

    const auto *lsf1_60 = buffer.data(lsf1 + 60);
    const auto *lsf1_62 = buffer.data(lsf1 + 62);
    const auto *lsf1_66 = buffer.data(lsf1 + 66);
    const auto *lsf1_67 = buffer.data(lsf1 + 67);
    const auto *lsf1_69 = buffer.data(lsf1 + 69);
    const auto *lsf1_75 = buffer.data(lsf1 + 75);
    const auto *lsf1_78 = buffer.data(lsf1 + 78);
    const auto *lsf1_79 = buffer.data(lsf1 + 79);
    const auto *lsf1_86 = buffer.data(lsf1 + 86);
    const auto *lsf1_88 = buffer.data(lsf1 + 88);
    const auto *lsf1_89 = buffer.data(lsf1 + 89);
    const auto *lsf1_90 = buffer.data(lsf1 + 90);
    const auto *lsf1_91 = buffer.data(lsf1 + 91);
    const auto *lsf1_92 = buffer.data(lsf1 + 92);
    const auto *lsf1_95 = buffer.data(lsf1 + 95);
    const auto *lsf1_96 = buffer.data(lsf1 + 96);
    const auto *lsf1_97 = buffer.data(lsf1 + 97);
    const auto *lsf1_98 = buffer.data(lsf1 + 98);
    const auto *lsf1_99 = buffer.data(lsf1 + 99);
    const auto *lsf1_100 = buffer.data(lsf1 + 100);
    const auto *lsf1_102 = buffer.data(lsf1 + 102);
    const auto *lsf1_103 = buffer.data(lsf1 + 103);
    const auto *lsf1_106 = buffer.data(lsf1 + 106);
    const auto *lsf1_107 = buffer.data(lsf1 + 107);
    const auto *lsf1_109 = buffer.data(lsf1 + 109);
    const auto *lsf1_115 = buffer.data(lsf1 + 115);
    const auto *lsf1_118 = buffer.data(lsf1 + 118);
    const auto *lsf1_119 = buffer.data(lsf1 + 119);

    const auto *lsg_91 = buffer.data(lsg + 91);
    const auto *lsg_92 = buffer.data(lsg + 92);
    const auto *lsg_93 = buffer.data(lsg + 93);
    const auto *lsg_95 = buffer.data(lsg + 95);
    const auto *lsg_96 = buffer.data(lsg + 96);
    const auto *lsg_100 = buffer.data(lsg + 100);
    const auto *lsg_101 = buffer.data(lsg + 101);
    const auto *lsg_102 = buffer.data(lsg + 102);
    const auto *lsg_103 = buffer.data(lsg + 103);
    const auto *lsg_104 = buffer.data(lsg + 104);
    const auto *lsg_105 = buffer.data(lsg + 105);
    const auto *lsg_107 = buffer.data(lsg + 107);
    const auto *lsg_108 = buffer.data(lsg + 108);
    const auto *lsg_110 = buffer.data(lsg + 110);
    const auto *lsg_114 = buffer.data(lsg + 114);
    const auto *lsg_115 = buffer.data(lsg + 115);
    const auto *lsg_116 = buffer.data(lsg + 116);
    const auto *lsg_117 = buffer.data(lsg + 117);
    const auto *lsg_118 = buffer.data(lsg + 118);
    const auto *lsg_119 = buffer.data(lsg + 119);
    const auto *lsg_120 = buffer.data(lsg + 120);
    const auto *lsg_122 = buffer.data(lsg + 122);
    const auto *lsg_123 = buffer.data(lsg + 123);
    const auto *lsg_125 = buffer.data(lsg + 125);
    const auto *lsg_130 = buffer.data(lsg + 130);
    const auto *lsg_131 = buffer.data(lsg + 131);
    const auto *lsg_132 = buffer.data(lsg + 132);
    const auto *lsg_133 = buffer.data(lsg + 133);
    const auto *lsg_134 = buffer.data(lsg + 134);
    const auto *lsg_135 = buffer.data(lsg + 135);
    const auto *lsg_136 = buffer.data(lsg + 136);
    const auto *lsg_137 = buffer.data(lsg + 137);
    const auto *lsg_138 = buffer.data(lsg + 138);
    const auto *lsg_139 = buffer.data(lsg + 139);
    const auto *lsg_140 = buffer.data(lsg + 140);
    const auto *lsg_144 = buffer.data(lsg + 144);
    const auto *lsg_145 = buffer.data(lsg + 145);
    const auto *lsg_146 = buffer.data(lsg + 146);
    const auto *lsg_147 = buffer.data(lsg + 147);
    const auto *lsg_148 = buffer.data(lsg + 148);
    const auto *lsg_149 = buffer.data(lsg + 149);
    const auto *lsg_150 = buffer.data(lsg + 150);
    const auto *lsg_151 = buffer.data(lsg + 151);
    const auto *lsg_152 = buffer.data(lsg + 152);
    const auto *lsg_153 = buffer.data(lsg + 153);
    const auto *lsg_155 = buffer.data(lsg + 155);
    const auto *lsg_156 = buffer.data(lsg + 156);
    const auto *lsg_160 = buffer.data(lsg + 160);
    const auto *lsg_161 = buffer.data(lsg + 161);
    const auto *lsg_162 = buffer.data(lsg + 162);
    const auto *lsg_163 = buffer.data(lsg + 163);
    const auto *lsg_164 = buffer.data(lsg + 164);
    const auto *lsg_165 = buffer.data(lsg + 165);
    const auto *lsg_167 = buffer.data(lsg + 167);
    const auto *lsg_168 = buffer.data(lsg + 168);
    const auto *lsg_170 = buffer.data(lsg + 170);
    const auto *lsg_174 = buffer.data(lsg + 174);
    const auto *lsg_175 = buffer.data(lsg + 175);
    const auto *lsg_176 = buffer.data(lsg + 176);
    const auto *lsg_177 = buffer.data(lsg + 177);
    const auto *lsg_178 = buffer.data(lsg + 178);
    const auto *lsg_179 = buffer.data(lsg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, ksg_96, lsf0_60, lsf0_66, \
                         lsf1_60, lsf1_66, lsg_91, lsg_92, lsg_93, \
                         lsg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * lsg_91[k];

        t_131[k] = f_4 * lsf0_60[k]
                   - f_5 * lsf1_60[k]
                   + f_3 * pc_z[k] * lsg_92[k];

        t_132[k] = f_16 * ksg_96[k]
                   + f_4 * lsf0_66[k]
                   - f_5 * lsf1_66[k]
                   + f_3 * pc_x[k] * lsg_96[k];

        t_133[k] = f_3 * pc_z[k] * lsg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, ksg_50, ksg_100, \
                         lsf0_62, lsf1_62, lsg_95, lsg_96, lsg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * ksg_50[k]
                   + f_3 * pc_y[k] * lsg_95[k];

        t_135[k] = f_6 * lsf0_62[k]
                   - f_7 * lsf1_62[k]
                   + f_3 * pc_z[k] * lsg_95[k];

        t_136[k] = f_16 * ksg_100[k]
                   + f_3 * pc_x[k] * lsg_100[k];

        t_137[k] = f_3 * pc_z[k] * lsg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, ksg_55, ksg_102, ksg_103, \
                         ksg_104, lsf0_66, lsf1_66, lsg_100, lsg_102, lsg_103, \
                         lsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * ksg_102[k]
                   + f_3 * pc_x[k] * lsg_102[k];

        t_139[k] = f_16 * ksg_103[k]
                   + f_3 * pc_x[k] * lsg_103[k];

        t_140[k] = f_16 * ksg_104[k]
                   + f_3 * pc_x[k] * lsg_104[k];

        t_141[k] = f_11 * ksg_55[k]
                   + f_1 * lsf0_66[k]
                   - f_2 * lsf1_66[k]
                   + f_3 * pc_y[k] * lsg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, ksg_59, lsf0_66, lsf0_67, \
                         lsf1_66, lsf1_67, lsg_100, lsg_101, lsg_102, \
                         lsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * lsg_100[k];

        t_143[k] = f_4 * lsf0_66[k]
                   - f_5 * lsf1_66[k]
                   + f_3 * pc_z[k] * lsg_101[k];

        t_144[k] = f_6 * lsf0_67[k]
                   - f_7 * lsf1_67[k]
                   + f_3 * pc_z[k] * lsg_102[k];

        t_145[k] = f_11 * ksg_59[k]
                   + f_3 * pc_y[k] * lsg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, ksh0_63, ksg_45, \
                         ksg_60, ksh1_63, lsf0_69, lsf1_69, lsg_104, \
                         lsg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * lsf0_69[k]
                   - f_2 * lsf1_69[k]
                   + f_3 * pc_z[k] * lsg_104[k];

        t_147[k] = pa_z[k] * ksh0_63[k]
                   - f_8 * pc_z[k] * ksh1_63[k];

        t_148[k] = f_10 * ksg_60[k]
                   + f_3 * pc_y[k] * lsg_105[k];

        t_149[k] = f_9 * ksg_45[k]
                   + f_3 * pc_z[k] * lsg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, ksh0_66, ksg_62, \
                         ksg_110, ksh1_66, lsf0_75, lsf1_75, lsg_107, \
                         lsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * ksh0_66[k]
                   - f_8 * pc_z[k] * ksh1_66[k];

        t_151[k] = f_10 * ksg_62[k]
                   + f_3 * pc_y[k] * lsg_107[k];

        t_152[k] = f_16 * ksg_110[k]
                   + f_6 * lsf0_75[k]
                   - f_7 * lsf1_75[k]
                   + f_3 * pc_x[k] * lsg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, ksh0_69, ksg_48, ksg_65, \
                         ksh1_69, lsg_108, lsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * ksh0_69[k]
                   - f_8 * pc_z[k] * ksh1_69[k];

        t_154[k] = f_9 * ksg_48[k]
                   + f_3 * pc_z[k] * lsg_108[k];

        t_155[k] = f_10 * ksg_65[k]
                   + f_3 * pc_y[k] * lsg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, ksg_114, ksg_115, ksg_116, ksg_117, \
                         lsf0_79, lsf1_79, lsg_114, lsg_115, lsg_116, \
                         lsg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * ksg_114[k]
                   + f_4 * lsf0_79[k]
                   - f_5 * lsf1_79[k]
                   + f_3 * pc_x[k] * lsg_114[k];

        t_157[k] = f_16 * ksg_115[k]
                   + f_3 * pc_x[k] * lsg_115[k];

        t_158[k] = f_16 * ksg_116[k]
                   + f_3 * pc_x[k] * lsg_116[k];

        t_159[k] = f_16 * ksg_117[k]
                   + f_3 * pc_x[k] * lsg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, ksh0_78, ksg_55, \
                         ksg_118, ksg_119, ksh1_78, lsg_115, lsg_118, \
                         lsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * ksg_118[k]
                   + f_3 * pc_x[k] * lsg_118[k];

        t_161[k] = f_16 * ksg_119[k]
                   + f_3 * pc_x[k] * lsg_119[k];

        t_162[k] = pa_z[k] * ksh0_78[k]
                   - f_8 * pc_z[k] * ksh1_78[k];

        t_163[k] = f_9 * ksg_55[k]
                   + f_3 * pc_z[k] * lsg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, ksg_72, ksg_73, ksg_74, lsf0_78, lsf0_79, \
                         lsf1_78, lsf1_79, lsg_117, lsg_118, lsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * ksg_72[k]
                   + f_6 * lsf0_78[k]
                   - f_7 * lsf1_78[k]
                   + f_3 * pc_y[k] * lsg_117[k];

        t_165[k] = f_10 * ksg_73[k]
                   + f_4 * lsf0_79[k]
                   - f_5 * lsf1_79[k]
                   + f_3 * pc_y[k] * lsg_118[k];

        t_166[k] = f_10 * ksg_74[k]
                   + f_3 * pc_y[k] * lsg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, ksh0_105, ksg_59, \
                         ksg_60, ksg_75, ksh1_105, lsf0_79, lsf1_79, lsg_119, \
                         lsg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * ksg_59[k]
                   + f_1 * lsf0_79[k]
                   - f_2 * lsf1_79[k]
                   + f_3 * pc_z[k] * lsg_119[k];

        t_168[k] = pa_y[k] * ksh0_105[k]
                   - f_8 * pc_y[k] * ksh1_105[k];

        t_169[k] = f_9 * ksg_75[k]
                   + f_3 * pc_y[k] * lsg_120[k];

        t_170[k] = f_10 * ksg_60[k]
                   + f_3 * pc_z[k] * lsg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, ksh0_108, ksh0_110, ksh0_111, \
                         ksg_76, ksg_77, ksg_78, ksh1_108, ksh1_110, ksh1_111, \
                         lsg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * ksh0_108[k]
                   + f_10 * ksg_76[k]
                   - f_8 * pc_y[k] * ksh1_108[k];

        t_172[k] = f_9 * ksg_77[k]
                   + f_3 * pc_y[k] * lsg_122[k];

        t_173[k] = pa_y[k] * ksh0_110[k]
                   - f_8 * pc_y[k] * ksh1_110[k];

        t_174[k] = pa_y[k] * ksh0_111[k]
                   + f_11 * ksg_78[k]
                   - f_8 * pc_y[k] * ksh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, ksh0_114, ksg_63, \
                         ksg_80, ksg_130, ksh1_114, lsg_123, lsg_125, \
                         lsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * ksg_63[k]
                   + f_3 * pc_z[k] * lsg_123[k];

        t_176[k] = f_9 * ksg_80[k]
                   + f_3 * pc_y[k] * lsg_125[k];

        t_177[k] = pa_y[k] * ksh0_114[k]
                   - f_8 * pc_y[k] * ksh1_114[k];

        t_178[k] = f_16 * ksg_130[k]
                   + f_3 * pc_x[k] * lsg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, ksg_131, ksg_132, ksg_133, ksg_134, \
                         lsg_131, lsg_132, lsg_133, lsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * ksg_131[k]
                   + f_3 * pc_x[k] * lsg_131[k];

        t_180[k] = f_16 * ksg_132[k]
                   + f_3 * pc_x[k] * lsg_132[k];

        t_181[k] = f_16 * ksg_133[k]
                   + f_3 * pc_x[k] * lsg_133[k];

        t_182[k] = f_16 * ksg_134[k]
                   + f_3 * pc_x[k] * lsg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, ksg_70, ksg_85, ksg_87, lsf0_86, \
                         lsf0_88, lsf1_86, lsf1_88, lsg_130, lsg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * ksg_85[k]
                   + f_1 * lsf0_86[k]
                   - f_2 * lsf1_86[k]
                   + f_3 * pc_y[k] * lsg_130[k];

        t_184[k] = f_10 * ksg_70[k]
                   + f_3 * pc_z[k] * lsg_130[k];

        t_185[k] = f_9 * ksg_87[k]
                   + f_6 * lsf0_88[k]
                   - f_7 * lsf1_88[k]
                   + f_3 * pc_y[k] * lsg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, ksh0_125, ksg_88, ksg_89, ksh1_125, \
                         lsf0_89, lsf1_89, lsg_133, lsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * ksg_88[k]
                   + f_4 * lsf0_89[k]
                   - f_5 * lsf1_89[k]
                   + f_3 * pc_y[k] * lsg_133[k];

        t_187[k] = f_9 * ksg_89[k]
                   + f_3 * pc_y[k] * lsg_134[k];

        t_188[k] = pa_y[k] * ksh0_125[k]
                   - f_8 * pc_y[k] * ksh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, ksg_75, ksg_135, \
                         lsf0_90, lsf1_90, lsg_135, lsg_136, lsg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * ksg_135[k]
                   + f_1 * lsf0_90[k]
                   - f_2 * lsf1_90[k]
                   + f_3 * pc_x[k] * lsg_135[k];

        t_190[k] = f_3 * pc_y[k] * lsg_135[k];

        t_191[k] = f_11 * ksg_75[k]
                   + f_3 * pc_z[k] * lsg_135[k];

        t_192[k] = f_4 * lsf0_90[k]
                   - f_5 * lsf1_90[k]
                   + f_3 * pc_y[k] * lsg_136[k];

        t_193[k] = f_3 * pc_y[k] * lsg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, ksg_140, lsf0_91, lsf0_92, \
                         lsf0_95, lsf1_91, lsf1_92, lsf1_95, lsg_138, lsg_139, \
                         lsg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * ksg_140[k]
                   + f_6 * lsf0_95[k]
                   - f_7 * lsf1_95[k]
                   + f_3 * pc_x[k] * lsg_140[k];

        t_195[k] = f_6 * lsf0_91[k]
                   - f_7 * lsf1_91[k]
                   + f_3 * pc_y[k] * lsg_138[k];

        t_196[k] = f_4 * lsf0_92[k]
                   - f_5 * lsf1_92[k]
                   + f_3 * pc_y[k] * lsg_139[k];

        t_197[k] = f_3 * pc_y[k] * lsg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, ksg_144, ksg_145, ksg_146, ksg_147, \
                         lsf0_99, lsf1_99, lsg_144, lsg_145, lsg_146, \
                         lsg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * ksg_144[k]
                   + f_4 * lsf0_99[k]
                   - f_5 * lsf1_99[k]
                   + f_3 * pc_x[k] * lsg_144[k];

        t_199[k] = f_16 * ksg_145[k]
                   + f_3 * pc_x[k] * lsg_145[k];

        t_200[k] = f_16 * ksg_146[k]
                   + f_3 * pc_x[k] * lsg_146[k];

        t_201[k] = f_16 * ksg_147[k]
                   + f_3 * pc_x[k] * lsg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, ksg_149, lsf0_96, lsf0_97, \
                         lsf1_96, lsf1_97, lsg_144, lsg_145, lsg_146, \
                         lsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * lsg_144[k];

        t_203[k] = f_16 * ksg_149[k]
                   + f_3 * pc_x[k] * lsg_149[k];

        t_204[k] = f_1 * lsf0_96[k]
                   - f_2 * lsf1_96[k]
                   + f_3 * pc_y[k] * lsg_145[k];

        t_205[k] = f_13 * lsf0_97[k]
                   - f_14 * lsf1_97[k]
                   + f_3 * pc_y[k] * lsg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, ksg_89, lsf0_98, lsf0_99, \
                         lsf1_98, lsf1_99, lsg_147, lsg_148, lsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * lsf0_98[k]
                   - f_7 * lsf1_98[k]
                   + f_3 * pc_y[k] * lsg_147[k];

        t_207[k] = f_4 * lsf0_99[k]
                   - f_5 * lsf1_99[k]
                   + f_3 * pc_y[k] * lsg_148[k];

        t_208[k] = f_3 * pc_y[k] * lsg_149[k];

        t_209[k] = f_11 * ksg_89[k]
                   + f_1 * lsf0_99[k]
                   - f_2 * lsf1_99[k]
                   + f_3 * pc_z[k] * lsg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, ksg_90, ksg_150, \
                         ksg_153, lsf0_100, lsf0_103, lsf1_100, lsf1_103, lsg_150, \
                         lsg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_17 * ksg_150[k]
                   + f_1 * lsf0_100[k]
                   - f_2 * lsf1_100[k]
                   + f_3 * pc_x[k] * lsg_150[k];

        t_211[k] = f_17 * ksg_90[k]
                   + f_3 * pc_y[k] * lsg_150[k];

        t_212[k] = f_3 * pc_z[k] * lsg_150[k];

        t_213[k] = f_17 * ksg_153[k]
                   + f_6 * lsf0_103[k]
                   - f_7 * lsf1_103[k]
                   + f_3 * pc_x[k] * lsg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, ksg_156, lsf0_100, lsf0_106, \
                         lsf1_100, lsf1_106, lsg_151, lsg_152, lsg_153, \
                         lsg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * lsg_151[k];

        t_215[k] = f_4 * lsf0_100[k]
                   - f_5 * lsf1_100[k]
                   + f_3 * pc_z[k] * lsg_152[k];

        t_216[k] = f_17 * ksg_156[k]
                   + f_4 * lsf0_106[k]
                   - f_5 * lsf1_106[k]
                   + f_3 * pc_x[k] * lsg_156[k];

        t_217[k] = f_3 * pc_z[k] * lsg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, ksg_95, ksg_160, \
                         lsf0_102, lsf1_102, lsg_155, lsg_156, \
                         lsg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_17 * ksg_95[k]
                   + f_3 * pc_y[k] * lsg_155[k];

        t_219[k] = f_6 * lsf0_102[k]
                   - f_7 * lsf1_102[k]
                   + f_3 * pc_z[k] * lsg_155[k];

        t_220[k] = f_17 * ksg_160[k]
                   + f_3 * pc_x[k] * lsg_160[k];

        t_221[k] = f_3 * pc_z[k] * lsg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, ksg_100, ksg_162, ksg_163, \
                         ksg_164, lsf0_106, lsf1_106, lsg_160, lsg_162, lsg_163, \
                         lsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_17 * ksg_162[k]
                   + f_3 * pc_x[k] * lsg_162[k];

        t_223[k] = f_17 * ksg_163[k]
                   + f_3 * pc_x[k] * lsg_163[k];

        t_224[k] = f_17 * ksg_164[k]
                   + f_3 * pc_x[k] * lsg_164[k];

        t_225[k] = f_17 * ksg_100[k]
                   + f_1 * lsf0_106[k]
                   - f_2 * lsf1_106[k]
                   + f_3 * pc_y[k] * lsg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, ksg_104, lsf0_106, lsf0_107, \
                         lsf1_106, lsf1_107, lsg_160, lsg_161, lsg_162, \
                         lsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * lsg_160[k];

        t_227[k] = f_4 * lsf0_106[k]
                   - f_5 * lsf1_106[k]
                   + f_3 * pc_z[k] * lsg_161[k];

        t_228[k] = f_6 * lsf0_107[k]
                   - f_7 * lsf1_107[k]
                   + f_3 * pc_z[k] * lsg_162[k];

        t_229[k] = f_17 * ksg_104[k]
                   + f_3 * pc_y[k] * lsg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, ksh0_126, ksg_90, \
                         ksg_105, ksh1_126, lsf0_109, lsf1_109, lsg_164, \
                         lsg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * lsf0_109[k]
                   - f_2 * lsf1_109[k]
                   + f_3 * pc_z[k] * lsg_164[k];

        t_231[k] = pa_z[k] * ksh0_126[k]
                   - f_8 * pc_z[k] * ksh1_126[k];

        t_232[k] = f_11 * ksg_105[k]
                   + f_3 * pc_y[k] * lsg_165[k];

        t_233[k] = f_9 * ksg_90[k]
                   + f_3 * pc_z[k] * lsg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, ksh0_129, ksg_107, \
                         ksg_170, ksh1_129, lsf0_115, lsf1_115, lsg_167, \
                         lsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * ksh0_129[k]
                   - f_8 * pc_z[k] * ksh1_129[k];

        t_235[k] = f_11 * ksg_107[k]
                   + f_3 * pc_y[k] * lsg_167[k];

        t_236[k] = f_17 * ksg_170[k]
                   + f_6 * lsf0_115[k]
                   - f_7 * lsf1_115[k]
                   + f_3 * pc_x[k] * lsg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, ksh0_132, ksg_93, ksg_110, \
                         ksh1_132, lsg_168, lsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * ksh0_132[k]
                   - f_8 * pc_z[k] * ksh1_132[k];

        t_238[k] = f_9 * ksg_93[k]
                   + f_3 * pc_z[k] * lsg_168[k];

        t_239[k] = f_11 * ksg_110[k]
                   + f_3 * pc_y[k] * lsg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, ksg_174, ksg_175, ksg_176, ksg_177, \
                         lsf0_119, lsf1_119, lsg_174, lsg_175, lsg_176, \
                         lsg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_17 * ksg_174[k]
                   + f_4 * lsf0_119[k]
                   - f_5 * lsf1_119[k]
                   + f_3 * pc_x[k] * lsg_174[k];

        t_241[k] = f_17 * ksg_175[k]
                   + f_3 * pc_x[k] * lsg_175[k];

        t_242[k] = f_17 * ksg_176[k]
                   + f_3 * pc_x[k] * lsg_176[k];

        t_243[k] = f_17 * ksg_177[k]
                   + f_3 * pc_x[k] * lsg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, ksh0_141, ksg_100, \
                         ksg_178, ksg_179, ksh1_141, lsg_175, lsg_178, \
                         lsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_17 * ksg_178[k]
                   + f_3 * pc_x[k] * lsg_178[k];

        t_245[k] = f_17 * ksg_179[k]
                   + f_3 * pc_x[k] * lsg_179[k];

        t_246[k] = pa_z[k] * ksh0_141[k]
                   - f_8 * pc_z[k] * ksh1_141[k];

        t_247[k] = f_9 * ksg_100[k]
                   + f_3 * pc_z[k] * lsg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, ksg_117, ksg_118, ksg_119, lsf0_118, \
                         lsf0_119, lsf1_118, lsf1_119, lsg_177, lsg_178, \
                         lsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * ksg_117[k]
                   + f_6 * lsf0_118[k]
                   - f_7 * lsf1_118[k]
                   + f_3 * pc_y[k] * lsg_177[k];

        t_249[k] = f_11 * ksg_118[k]
                   + f_4 * lsf0_119[k]
                   - f_5 * lsf1_119[k]
                   + f_3 * pc_y[k] * lsg_178[k];

        t_250[k] = f_11 * ksg_119[k]
                   + f_3 * pc_y[k] * lsg_179[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

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
    auto *t_360 = buffer.data(target + 360);
    auto *t_361 = buffer.data(target + 361);
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_189 = buffer.data(ksh0 + 189);
    const auto *ksh0_192 = buffer.data(ksh0 + 192);
    const auto *ksh0_194 = buffer.data(ksh0 + 194);
    const auto *ksh0_195 = buffer.data(ksh0 + 195);
    const auto *ksh0_198 = buffer.data(ksh0 + 198);
    const auto *ksh0_209 = buffer.data(ksh0 + 209);
    const auto *ksh0_210 = buffer.data(ksh0 + 210);
    const auto *ksh0_213 = buffer.data(ksh0 + 213);
    const auto *ksh0_216 = buffer.data(ksh0 + 216);
    const auto *ksh0_225 = buffer.data(ksh0 + 225);

    const auto *ksg_104 = buffer.data(ksg + 104);
    const auto *ksg_105 = buffer.data(ksg + 105);
    const auto *ksg_108 = buffer.data(ksg + 108);
    const auto *ksg_115 = buffer.data(ksg + 115);
    const auto *ksg_119 = buffer.data(ksg + 119);
    const auto *ksg_120 = buffer.data(ksg + 120);
    const auto *ksg_122 = buffer.data(ksg + 122);
    const auto *ksg_123 = buffer.data(ksg + 123);
    const auto *ksg_125 = buffer.data(ksg + 125);
    const auto *ksg_130 = buffer.data(ksg + 130);
    const auto *ksg_132 = buffer.data(ksg + 132);
    const auto *ksg_133 = buffer.data(ksg + 133);
    const auto *ksg_134 = buffer.data(ksg + 134);
    const auto *ksg_135 = buffer.data(ksg + 135);
    const auto *ksg_136 = buffer.data(ksg + 136);
    const auto *ksg_137 = buffer.data(ksg + 137);
    const auto *ksg_138 = buffer.data(ksg + 138);
    const auto *ksg_140 = buffer.data(ksg + 140);
    const auto *ksg_145 = buffer.data(ksg + 145);
    const auto *ksg_147 = buffer.data(ksg + 147);
    const auto *ksg_148 = buffer.data(ksg + 148);
    const auto *ksg_149 = buffer.data(ksg + 149);
    const auto *ksg_150 = buffer.data(ksg + 150);
    const auto *ksg_153 = buffer.data(ksg + 153);
    const auto *ksg_155 = buffer.data(ksg + 155);
    const auto *ksg_160 = buffer.data(ksg + 160);
    const auto *ksg_164 = buffer.data(ksg + 164);
    const auto *ksg_165 = buffer.data(ksg + 165);
    const auto *ksg_167 = buffer.data(ksg + 167);
    const auto *ksg_168 = buffer.data(ksg + 168);
    const auto *ksg_170 = buffer.data(ksg + 170);
    const auto *ksg_177 = buffer.data(ksg + 177);
    const auto *ksg_178 = buffer.data(ksg + 178);
    const auto *ksg_179 = buffer.data(ksg + 179);
    const auto *ksg_180 = buffer.data(ksg + 180);
    const auto *ksg_182 = buffer.data(ksg + 182);
    const auto *ksg_183 = buffer.data(ksg + 183);
    const auto *ksg_185 = buffer.data(ksg + 185);
    const auto *ksg_186 = buffer.data(ksg + 186);
    const auto *ksg_189 = buffer.data(ksg + 189);
    const auto *ksg_190 = buffer.data(ksg + 190);
    const auto *ksg_191 = buffer.data(ksg + 191);
    const auto *ksg_192 = buffer.data(ksg + 192);
    const auto *ksg_193 = buffer.data(ksg + 193);
    const auto *ksg_194 = buffer.data(ksg + 194);
    const auto *ksg_205 = buffer.data(ksg + 205);
    const auto *ksg_206 = buffer.data(ksg + 206);
    const auto *ksg_207 = buffer.data(ksg + 207);
    const auto *ksg_208 = buffer.data(ksg + 208);
    const auto *ksg_209 = buffer.data(ksg + 209);
    const auto *ksg_210 = buffer.data(ksg + 210);
    const auto *ksg_215 = buffer.data(ksg + 215);
    const auto *ksg_219 = buffer.data(ksg + 219);
    const auto *ksg_220 = buffer.data(ksg + 220);
    const auto *ksg_221 = buffer.data(ksg + 221);
    const auto *ksg_222 = buffer.data(ksg + 222);
    const auto *ksg_224 = buffer.data(ksg + 224);
    const auto *ksg_225 = buffer.data(ksg + 225);
    const auto *ksg_228 = buffer.data(ksg + 228);
    const auto *ksg_231 = buffer.data(ksg + 231);
    const auto *ksg_235 = buffer.data(ksg + 235);
    const auto *ksg_237 = buffer.data(ksg + 237);
    const auto *ksg_238 = buffer.data(ksg + 238);
    const auto *ksg_239 = buffer.data(ksg + 239);
    const auto *ksg_245 = buffer.data(ksg + 245);
    const auto *ksg_249 = buffer.data(ksg + 249);
    const auto *ksg_250 = buffer.data(ksg + 250);
    const auto *ksg_251 = buffer.data(ksg + 251);
    const auto *ksg_252 = buffer.data(ksg + 252);
    const auto *ksg_253 = buffer.data(ksg + 253);
    const auto *ksg_254 = buffer.data(ksg + 254);
    const auto *ksg_255 = buffer.data(ksg + 255);
    const auto *ksg_258 = buffer.data(ksg + 258);
    const auto *ksg_260 = buffer.data(ksg + 260);
    const auto *ksg_261 = buffer.data(ksg + 261);
    const auto *ksg_264 = buffer.data(ksg + 264);
    const auto *ksg_265 = buffer.data(ksg + 265);
    const auto *ksg_266 = buffer.data(ksg + 266);

    const auto *ksh1_189 = buffer.data(ksh1 + 189);
    const auto *ksh1_192 = buffer.data(ksh1 + 192);
    const auto *ksh1_194 = buffer.data(ksh1 + 194);
    const auto *ksh1_195 = buffer.data(ksh1 + 195);
    const auto *ksh1_198 = buffer.data(ksh1 + 198);
    const auto *ksh1_209 = buffer.data(ksh1 + 209);
    const auto *ksh1_210 = buffer.data(ksh1 + 210);
    const auto *ksh1_213 = buffer.data(ksh1 + 213);
    const auto *ksh1_216 = buffer.data(ksh1 + 216);
    const auto *ksh1_225 = buffer.data(ksh1 + 225);

    const auto *lsf0_119 = buffer.data(lsf0 + 119);
    const auto *lsf0_120 = buffer.data(lsf0 + 120);
    const auto *lsf0_123 = buffer.data(lsf0 + 123);
    const auto *lsf0_125 = buffer.data(lsf0 + 125);
    const auto *lsf0_126 = buffer.data(lsf0 + 126);
    const auto *lsf0_128 = buffer.data(lsf0 + 128);
    const auto *lsf0_129 = buffer.data(lsf0 + 129);
    const auto *lsf0_136 = buffer.data(lsf0 + 136);
    const auto *lsf0_138 = buffer.data(lsf0 + 138);
    const auto *lsf0_139 = buffer.data(lsf0 + 139);
    const auto *lsf0_140 = buffer.data(lsf0 + 140);
    const auto *lsf0_141 = buffer.data(lsf0 + 141);
    const auto *lsf0_142 = buffer.data(lsf0 + 142);
    const auto *lsf0_145 = buffer.data(lsf0 + 145);
    const auto *lsf0_146 = buffer.data(lsf0 + 146);
    const auto *lsf0_147 = buffer.data(lsf0 + 147);
    const auto *lsf0_148 = buffer.data(lsf0 + 148);
    const auto *lsf0_149 = buffer.data(lsf0 + 149);
    const auto *lsf0_150 = buffer.data(lsf0 + 150);
    const auto *lsf0_152 = buffer.data(lsf0 + 152);
    const auto *lsf0_153 = buffer.data(lsf0 + 153);
    const auto *lsf0_156 = buffer.data(lsf0 + 156);
    const auto *lsf0_157 = buffer.data(lsf0 + 157);
    const auto *lsf0_159 = buffer.data(lsf0 + 159);
    const auto *lsf0_165 = buffer.data(lsf0 + 165);
    const auto *lsf0_168 = buffer.data(lsf0 + 168);
    const auto *lsf0_169 = buffer.data(lsf0 + 169);
    const auto *lsf0_170 = buffer.data(lsf0 + 170);
    const auto *lsf0_173 = buffer.data(lsf0 + 173);
    const auto *lsf0_175 = buffer.data(lsf0 + 175);
    const auto *lsf0_176 = buffer.data(lsf0 + 176);
    const auto *lsf0_179 = buffer.data(lsf0 + 179);

    const auto *lsf1_119 = buffer.data(lsf1 + 119);
    const auto *lsf1_120 = buffer.data(lsf1 + 120);
    const auto *lsf1_123 = buffer.data(lsf1 + 123);
    const auto *lsf1_125 = buffer.data(lsf1 + 125);
    const auto *lsf1_126 = buffer.data(lsf1 + 126);
    const auto *lsf1_128 = buffer.data(lsf1 + 128);
    const auto *lsf1_129 = buffer.data(lsf1 + 129);
    const auto *lsf1_136 = buffer.data(lsf1 + 136);
    const auto *lsf1_138 = buffer.data(lsf1 + 138);
    const auto *lsf1_139 = buffer.data(lsf1 + 139);
    const auto *lsf1_140 = buffer.data(lsf1 + 140);
    const auto *lsf1_141 = buffer.data(lsf1 + 141);
    const auto *lsf1_142 = buffer.data(lsf1 + 142);
    const auto *lsf1_145 = buffer.data(lsf1 + 145);
    const auto *lsf1_146 = buffer.data(lsf1 + 146);
    const auto *lsf1_147 = buffer.data(lsf1 + 147);
    const auto *lsf1_148 = buffer.data(lsf1 + 148);
    const auto *lsf1_149 = buffer.data(lsf1 + 149);
    const auto *lsf1_150 = buffer.data(lsf1 + 150);
    const auto *lsf1_152 = buffer.data(lsf1 + 152);
    const auto *lsf1_153 = buffer.data(lsf1 + 153);
    const auto *lsf1_156 = buffer.data(lsf1 + 156);
    const auto *lsf1_157 = buffer.data(lsf1 + 157);
    const auto *lsf1_159 = buffer.data(lsf1 + 159);
    const auto *lsf1_165 = buffer.data(lsf1 + 165);
    const auto *lsf1_168 = buffer.data(lsf1 + 168);
    const auto *lsf1_169 = buffer.data(lsf1 + 169);
    const auto *lsf1_170 = buffer.data(lsf1 + 170);
    const auto *lsf1_173 = buffer.data(lsf1 + 173);
    const auto *lsf1_175 = buffer.data(lsf1 + 175);
    const auto *lsf1_176 = buffer.data(lsf1 + 176);
    const auto *lsf1_179 = buffer.data(lsf1 + 179);

    const auto *lsg_179 = buffer.data(lsg + 179);
    const auto *lsg_180 = buffer.data(lsg + 180);
    const auto *lsg_182 = buffer.data(lsg + 182);
    const auto *lsg_183 = buffer.data(lsg + 183);
    const auto *lsg_185 = buffer.data(lsg + 185);
    const auto *lsg_186 = buffer.data(lsg + 186);
    const auto *lsg_189 = buffer.data(lsg + 189);
    const auto *lsg_190 = buffer.data(lsg + 190);
    const auto *lsg_191 = buffer.data(lsg + 191);
    const auto *lsg_192 = buffer.data(lsg + 192);
    const auto *lsg_193 = buffer.data(lsg + 193);
    const auto *lsg_194 = buffer.data(lsg + 194);
    const auto *lsg_195 = buffer.data(lsg + 195);
    const auto *lsg_197 = buffer.data(lsg + 197);
    const auto *lsg_198 = buffer.data(lsg + 198);
    const auto *lsg_200 = buffer.data(lsg + 200);
    const auto *lsg_205 = buffer.data(lsg + 205);
    const auto *lsg_206 = buffer.data(lsg + 206);
    const auto *lsg_207 = buffer.data(lsg + 207);
    const auto *lsg_208 = buffer.data(lsg + 208);
    const auto *lsg_209 = buffer.data(lsg + 209);
    const auto *lsg_210 = buffer.data(lsg + 210);
    const auto *lsg_211 = buffer.data(lsg + 211);
    const auto *lsg_212 = buffer.data(lsg + 212);
    const auto *lsg_213 = buffer.data(lsg + 213);
    const auto *lsg_214 = buffer.data(lsg + 214);
    const auto *lsg_215 = buffer.data(lsg + 215);
    const auto *lsg_219 = buffer.data(lsg + 219);
    const auto *lsg_220 = buffer.data(lsg + 220);
    const auto *lsg_221 = buffer.data(lsg + 221);
    const auto *lsg_222 = buffer.data(lsg + 222);
    const auto *lsg_223 = buffer.data(lsg + 223);
    const auto *lsg_224 = buffer.data(lsg + 224);
    const auto *lsg_225 = buffer.data(lsg + 225);
    const auto *lsg_226 = buffer.data(lsg + 226);
    const auto *lsg_227 = buffer.data(lsg + 227);
    const auto *lsg_228 = buffer.data(lsg + 228);
    const auto *lsg_230 = buffer.data(lsg + 230);
    const auto *lsg_231 = buffer.data(lsg + 231);
    const auto *lsg_235 = buffer.data(lsg + 235);
    const auto *lsg_236 = buffer.data(lsg + 236);
    const auto *lsg_237 = buffer.data(lsg + 237);
    const auto *lsg_238 = buffer.data(lsg + 238);
    const auto *lsg_239 = buffer.data(lsg + 239);
    const auto *lsg_240 = buffer.data(lsg + 240);
    const auto *lsg_242 = buffer.data(lsg + 242);
    const auto *lsg_243 = buffer.data(lsg + 243);
    const auto *lsg_245 = buffer.data(lsg + 245);
    const auto *lsg_249 = buffer.data(lsg + 249);
    const auto *lsg_250 = buffer.data(lsg + 250);
    const auto *lsg_251 = buffer.data(lsg + 251);
    const auto *lsg_252 = buffer.data(lsg + 252);
    const auto *lsg_253 = buffer.data(lsg + 253);
    const auto *lsg_254 = buffer.data(lsg + 254);
    const auto *lsg_255 = buffer.data(lsg + 255);
    const auto *lsg_257 = buffer.data(lsg + 257);
    const auto *lsg_258 = buffer.data(lsg + 258);
    const auto *lsg_260 = buffer.data(lsg + 260);
    const auto *lsg_261 = buffer.data(lsg + 261);
    const auto *lsg_264 = buffer.data(lsg + 264);
    const auto *lsg_265 = buffer.data(lsg + 265);
    const auto *lsg_266 = buffer.data(lsg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, ksg_104, ksg_120, ksg_180, \
                         lsf0_119, lsf0_120, lsf1_119, lsf1_120, lsg_179, \
                         lsg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * ksg_104[k]
                   + f_1 * lsf0_119[k]
                   - f_2 * lsf1_119[k]
                   + f_3 * pc_z[k] * lsg_179[k];

        t_252[k] = f_17 * ksg_180[k]
                   + f_1 * lsf0_120[k]
                   - f_2 * lsf1_120[k]
                   + f_3 * pc_x[k] * lsg_180[k];

        t_253[k] = f_10 * ksg_120[k]
                   + f_3 * pc_y[k] * lsg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, ksg_105, ksg_122, ksg_183, \
                         lsf0_123, lsf1_123, lsg_180, lsg_182, \
                         lsg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * ksg_105[k]
                   + f_3 * pc_z[k] * lsg_180[k];

        t_255[k] = f_17 * ksg_183[k]
                   + f_6 * lsf0_123[k]
                   - f_7 * lsf1_123[k]
                   + f_3 * pc_x[k] * lsg_183[k];

        t_256[k] = f_10 * ksg_122[k]
                   + f_3 * pc_y[k] * lsg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, ksg_108, ksg_185, ksg_186, lsf0_125, \
                         lsf0_126, lsf1_125, lsf1_126, lsg_183, lsg_185, \
                         lsg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_17 * ksg_185[k]
                   + f_6 * lsf0_125[k]
                   - f_7 * lsf1_125[k]
                   + f_3 * pc_x[k] * lsg_185[k];

        t_258[k] = f_17 * ksg_186[k]
                   + f_4 * lsf0_126[k]
                   - f_5 * lsf1_126[k]
                   + f_3 * pc_x[k] * lsg_186[k];

        t_259[k] = f_10 * ksg_108[k]
                   + f_3 * pc_z[k] * lsg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, ksg_125, ksg_189, ksg_190, \
                         ksg_191, lsf0_129, lsf1_129, lsg_185, lsg_189, lsg_190, \
                         lsg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * ksg_125[k]
                   + f_3 * pc_y[k] * lsg_185[k];

        t_261[k] = f_17 * ksg_189[k]
                   + f_4 * lsf0_129[k]
                   - f_5 * lsf1_129[k]
                   + f_3 * pc_x[k] * lsg_189[k];

        t_262[k] = f_17 * ksg_190[k]
                   + f_3 * pc_x[k] * lsg_190[k];

        t_263[k] = f_17 * ksg_191[k]
                   + f_3 * pc_x[k] * lsg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, ksg_130, ksg_192, ksg_193, \
                         ksg_194, lsf0_126, lsf1_126, lsg_190, lsg_192, lsg_193, \
                         lsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * ksg_192[k]
                   + f_3 * pc_x[k] * lsg_192[k];

        t_265[k] = f_17 * ksg_193[k]
                   + f_3 * pc_x[k] * lsg_193[k];

        t_266[k] = f_17 * ksg_194[k]
                   + f_3 * pc_x[k] * lsg_194[k];

        t_267[k] = f_10 * ksg_130[k]
                   + f_1 * lsf0_126[k]
                   - f_2 * lsf1_126[k]
                   + f_3 * pc_y[k] * lsg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, ksg_115, ksg_132, ksg_133, lsf0_128, \
                         lsf0_129, lsf1_128, lsf1_129, lsg_190, lsg_192, \
                         lsg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * ksg_115[k]
                   + f_3 * pc_z[k] * lsg_190[k];

        t_269[k] = f_10 * ksg_132[k]
                   + f_6 * lsf0_128[k]
                   - f_7 * lsf1_128[k]
                   + f_3 * pc_y[k] * lsg_192[k];

        t_270[k] = f_10 * ksg_133[k]
                   + f_4 * lsf0_129[k]
                   - f_5 * lsf1_129[k]
                   + f_3 * pc_y[k] * lsg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, ksh0_189, ksg_119, \
                         ksg_134, ksg_135, ksh1_189, lsf0_129, lsf1_129, lsg_194, \
                         lsg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * ksg_134[k]
                   + f_3 * pc_y[k] * lsg_194[k];

        t_272[k] = f_10 * ksg_119[k]
                   + f_1 * lsf0_129[k]
                   - f_2 * lsf1_129[k]
                   + f_3 * pc_z[k] * lsg_194[k];

        t_273[k] = pa_y[k] * ksh0_189[k]
                   - f_8 * pc_y[k] * ksh1_189[k];

        t_274[k] = f_9 * ksg_135[k]
                   + f_3 * pc_y[k] * lsg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, ksh0_192, ksh0_194, \
                         ksg_120, ksg_136, ksg_137, ksh1_192, ksh1_194, lsg_195, \
                         lsg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * ksg_120[k]
                   + f_3 * pc_z[k] * lsg_195[k];

        t_276[k] = pa_y[k] * ksh0_192[k]
                   + f_10 * ksg_136[k]
                   - f_8 * pc_y[k] * ksh1_192[k];

        t_277[k] = f_9 * ksg_137[k]
                   + f_3 * pc_y[k] * lsg_197[k];

        t_278[k] = pa_y[k] * ksh0_194[k]
                   - f_8 * pc_y[k] * ksh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, ksh0_195, ksh0_198, \
                         ksg_123, ksg_138, ksg_140, ksh1_195, ksh1_198, lsg_198, \
                         lsg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * ksh0_195[k]
                   + f_11 * ksg_138[k]
                   - f_8 * pc_y[k] * ksh1_195[k];

        t_280[k] = f_11 * ksg_123[k]
                   + f_3 * pc_z[k] * lsg_198[k];

        t_281[k] = f_9 * ksg_140[k]
                   + f_3 * pc_y[k] * lsg_200[k];

        t_282[k] = pa_y[k] * ksh0_198[k]
                   - f_8 * pc_y[k] * ksh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, ksg_205, ksg_206, ksg_207, \
                         ksg_208, ksg_209, lsg_205, lsg_206, lsg_207, lsg_208, \
                         lsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_17 * ksg_205[k]
                   + f_3 * pc_x[k] * lsg_205[k];

        t_284[k] = f_17 * ksg_206[k]
                   + f_3 * pc_x[k] * lsg_206[k];

        t_285[k] = f_17 * ksg_207[k]
                   + f_3 * pc_x[k] * lsg_207[k];

        t_286[k] = f_17 * ksg_208[k]
                   + f_3 * pc_x[k] * lsg_208[k];

        t_287[k] = f_17 * ksg_209[k]
                   + f_3 * pc_x[k] * lsg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, ksg_130, ksg_145, ksg_147, lsf0_136, \
                         lsf0_138, lsf1_136, lsf1_138, lsg_205, \
                         lsg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * ksg_145[k]
                   + f_1 * lsf0_136[k]
                   - f_2 * lsf1_136[k]
                   + f_3 * pc_y[k] * lsg_205[k];

        t_289[k] = f_11 * ksg_130[k]
                   + f_3 * pc_z[k] * lsg_205[k];

        t_290[k] = f_9 * ksg_147[k]
                   + f_6 * lsf0_138[k]
                   - f_7 * lsf1_138[k]
                   + f_3 * pc_y[k] * lsg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, ksh0_209, ksg_148, ksg_149, \
                         ksh1_209, lsf0_139, lsf1_139, lsg_208, \
                         lsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * ksg_148[k]
                   + f_4 * lsf0_139[k]
                   - f_5 * lsf1_139[k]
                   + f_3 * pc_y[k] * lsg_208[k];

        t_292[k] = f_9 * ksg_149[k]
                   + f_3 * pc_y[k] * lsg_209[k];

        t_293[k] = pa_y[k] * ksh0_209[k]
                   - f_8 * pc_y[k] * ksh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, ksg_135, \
                         ksg_210, lsf0_140, lsf1_140, lsg_210, lsg_211, \
                         lsg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_17 * ksg_210[k]
                   + f_1 * lsf0_140[k]
                   - f_2 * lsf1_140[k]
                   + f_3 * pc_x[k] * lsg_210[k];

        t_295[k] = f_3 * pc_y[k] * lsg_210[k];

        t_296[k] = f_17 * ksg_135[k]
                   + f_3 * pc_z[k] * lsg_210[k];

        t_297[k] = f_4 * lsf0_140[k]
                   - f_5 * lsf1_140[k]
                   + f_3 * pc_y[k] * lsg_211[k];

        t_298[k] = f_3 * pc_y[k] * lsg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, ksg_215, lsf0_141, lsf0_142, \
                         lsf0_145, lsf1_141, lsf1_142, lsf1_145, lsg_213, lsg_214, \
                         lsg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_17 * ksg_215[k]
                   + f_6 * lsf0_145[k]
                   - f_7 * lsf1_145[k]
                   + f_3 * pc_x[k] * lsg_215[k];

        t_300[k] = f_6 * lsf0_141[k]
                   - f_7 * lsf1_141[k]
                   + f_3 * pc_y[k] * lsg_213[k];

        t_301[k] = f_4 * lsf0_142[k]
                   - f_5 * lsf1_142[k]
                   + f_3 * pc_y[k] * lsg_214[k];

        t_302[k] = f_3 * pc_y[k] * lsg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, ksg_219, ksg_220, ksg_221, ksg_222, \
                         lsf0_149, lsf1_149, lsg_219, lsg_220, lsg_221, \
                         lsg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_17 * ksg_219[k]
                   + f_4 * lsf0_149[k]
                   - f_5 * lsf1_149[k]
                   + f_3 * pc_x[k] * lsg_219[k];

        t_304[k] = f_17 * ksg_220[k]
                   + f_3 * pc_x[k] * lsg_220[k];

        t_305[k] = f_17 * ksg_221[k]
                   + f_3 * pc_x[k] * lsg_221[k];

        t_306[k] = f_17 * ksg_222[k]
                   + f_3 * pc_x[k] * lsg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, ksg_224, lsf0_146, lsf0_147, \
                         lsf1_146, lsf1_147, lsg_219, lsg_220, lsg_221, \
                         lsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * lsg_219[k];

        t_308[k] = f_17 * ksg_224[k]
                   + f_3 * pc_x[k] * lsg_224[k];

        t_309[k] = f_1 * lsf0_146[k]
                   - f_2 * lsf1_146[k]
                   + f_3 * pc_y[k] * lsg_220[k];

        t_310[k] = f_13 * lsf0_147[k]
                   - f_14 * lsf1_147[k]
                   + f_3 * pc_y[k] * lsg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, ksg_149, lsf0_148, lsf0_149, \
                         lsf1_148, lsf1_149, lsg_222, lsg_223, \
                         lsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * lsf0_148[k]
                   - f_7 * lsf1_148[k]
                   + f_3 * pc_y[k] * lsg_222[k];

        t_312[k] = f_4 * lsf0_149[k]
                   - f_5 * lsf1_149[k]
                   + f_3 * pc_y[k] * lsg_223[k];

        t_313[k] = f_3 * pc_y[k] * lsg_224[k];

        t_314[k] = f_17 * ksg_149[k]
                   + f_1 * lsf0_149[k]
                   - f_2 * lsf1_149[k]
                   + f_3 * pc_z[k] * lsg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, ksg_150, ksg_225, \
                         ksg_228, lsf0_150, lsf0_153, lsf1_150, lsf1_153, lsg_225, \
                         lsg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * ksg_225[k]
                   + f_1 * lsf0_150[k]
                   - f_2 * lsf1_150[k]
                   + f_3 * pc_x[k] * lsg_225[k];

        t_316[k] = f_16 * ksg_150[k]
                   + f_3 * pc_y[k] * lsg_225[k];

        t_317[k] = f_3 * pc_z[k] * lsg_225[k];

        t_318[k] = f_11 * ksg_228[k]
                   + f_6 * lsf0_153[k]
                   - f_7 * lsf1_153[k]
                   + f_3 * pc_x[k] * lsg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, ksg_231, lsf0_150, lsf0_156, \
                         lsf1_150, lsf1_156, lsg_226, lsg_227, lsg_228, \
                         lsg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * lsg_226[k];

        t_320[k] = f_4 * lsf0_150[k]
                   - f_5 * lsf1_150[k]
                   + f_3 * pc_z[k] * lsg_227[k];

        t_321[k] = f_11 * ksg_231[k]
                   + f_4 * lsf0_156[k]
                   - f_5 * lsf1_156[k]
                   + f_3 * pc_x[k] * lsg_231[k];

        t_322[k] = f_3 * pc_z[k] * lsg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, ksg_155, ksg_235, \
                         lsf0_152, lsf1_152, lsg_230, lsg_231, \
                         lsg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_16 * ksg_155[k]
                   + f_3 * pc_y[k] * lsg_230[k];

        t_324[k] = f_6 * lsf0_152[k]
                   - f_7 * lsf1_152[k]
                   + f_3 * pc_z[k] * lsg_230[k];

        t_325[k] = f_11 * ksg_235[k]
                   + f_3 * pc_x[k] * lsg_235[k];

        t_326[k] = f_3 * pc_z[k] * lsg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, ksg_160, ksg_237, ksg_238, \
                         ksg_239, lsf0_156, lsf1_156, lsg_235, lsg_237, lsg_238, \
                         lsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_11 * ksg_237[k]
                   + f_3 * pc_x[k] * lsg_237[k];

        t_328[k] = f_11 * ksg_238[k]
                   + f_3 * pc_x[k] * lsg_238[k];

        t_329[k] = f_11 * ksg_239[k]
                   + f_3 * pc_x[k] * lsg_239[k];

        t_330[k] = f_16 * ksg_160[k]
                   + f_1 * lsf0_156[k]
                   - f_2 * lsf1_156[k]
                   + f_3 * pc_y[k] * lsg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, ksg_164, lsf0_156, lsf0_157, \
                         lsf1_156, lsf1_157, lsg_235, lsg_236, lsg_237, \
                         lsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * lsg_235[k];

        t_332[k] = f_4 * lsf0_156[k]
                   - f_5 * lsf1_156[k]
                   + f_3 * pc_z[k] * lsg_236[k];

        t_333[k] = f_6 * lsf0_157[k]
                   - f_7 * lsf1_157[k]
                   + f_3 * pc_z[k] * lsg_237[k];

        t_334[k] = f_16 * ksg_164[k]
                   + f_3 * pc_y[k] * lsg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, ksh0_210, ksg_150, \
                         ksg_165, ksh1_210, lsf0_159, lsf1_159, lsg_239, \
                         lsg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * lsf0_159[k]
                   - f_2 * lsf1_159[k]
                   + f_3 * pc_z[k] * lsg_239[k];

        t_336[k] = pa_z[k] * ksh0_210[k]
                   - f_8 * pc_z[k] * ksh1_210[k];

        t_337[k] = f_17 * ksg_165[k]
                   + f_3 * pc_y[k] * lsg_240[k];

        t_338[k] = f_9 * ksg_150[k]
                   + f_3 * pc_z[k] * lsg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, ksh0_213, ksg_167, \
                         ksg_245, ksh1_213, lsf0_165, lsf1_165, lsg_242, \
                         lsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * ksh0_213[k]
                   - f_8 * pc_z[k] * ksh1_213[k];

        t_340[k] = f_17 * ksg_167[k]
                   + f_3 * pc_y[k] * lsg_242[k];

        t_341[k] = f_11 * ksg_245[k]
                   + f_6 * lsf0_165[k]
                   - f_7 * lsf1_165[k]
                   + f_3 * pc_x[k] * lsg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, ksh0_216, ksg_153, ksg_170, \
                         ksh1_216, lsg_243, lsg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * ksh0_216[k]
                   - f_8 * pc_z[k] * ksh1_216[k];

        t_343[k] = f_9 * ksg_153[k]
                   + f_3 * pc_z[k] * lsg_243[k];

        t_344[k] = f_17 * ksg_170[k]
                   + f_3 * pc_y[k] * lsg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, ksg_249, ksg_250, ksg_251, ksg_252, \
                         lsf0_169, lsf1_169, lsg_249, lsg_250, lsg_251, \
                         lsg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_11 * ksg_249[k]
                   + f_4 * lsf0_169[k]
                   - f_5 * lsf1_169[k]
                   + f_3 * pc_x[k] * lsg_249[k];

        t_346[k] = f_11 * ksg_250[k]
                   + f_3 * pc_x[k] * lsg_250[k];

        t_347[k] = f_11 * ksg_251[k]
                   + f_3 * pc_x[k] * lsg_251[k];

        t_348[k] = f_11 * ksg_252[k]
                   + f_3 * pc_x[k] * lsg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, ksh0_225, ksg_160, \
                         ksg_253, ksg_254, ksh1_225, lsg_250, lsg_253, \
                         lsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_11 * ksg_253[k]
                   + f_3 * pc_x[k] * lsg_253[k];

        t_350[k] = f_11 * ksg_254[k]
                   + f_3 * pc_x[k] * lsg_254[k];

        t_351[k] = pa_z[k] * ksh0_225[k]
                   - f_8 * pc_z[k] * ksh1_225[k];

        t_352[k] = f_9 * ksg_160[k]
                   + f_3 * pc_z[k] * lsg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, ksg_177, ksg_178, ksg_179, lsf0_168, \
                         lsf0_169, lsf1_168, lsf1_169, lsg_252, lsg_253, \
                         lsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_17 * ksg_177[k]
                   + f_6 * lsf0_168[k]
                   - f_7 * lsf1_168[k]
                   + f_3 * pc_y[k] * lsg_252[k];

        t_354[k] = f_17 * ksg_178[k]
                   + f_4 * lsf0_169[k]
                   - f_5 * lsf1_169[k]
                   + f_3 * pc_y[k] * lsg_253[k];

        t_355[k] = f_17 * ksg_179[k]
                   + f_3 * pc_y[k] * lsg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, ksg_164, ksg_180, ksg_255, \
                         lsf0_169, lsf0_170, lsf1_169, lsf1_170, lsg_254, \
                         lsg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * ksg_164[k]
                   + f_1 * lsf0_169[k]
                   - f_2 * lsf1_169[k]
                   + f_3 * pc_z[k] * lsg_254[k];

        t_357[k] = f_11 * ksg_255[k]
                   + f_1 * lsf0_170[k]
                   - f_2 * lsf1_170[k]
                   + f_3 * pc_x[k] * lsg_255[k];

        t_358[k] = f_11 * ksg_180[k]
                   + f_3 * pc_y[k] * lsg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, ksg_165, ksg_182, ksg_258, \
                         lsf0_173, lsf1_173, lsg_255, lsg_257, \
                         lsg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * ksg_165[k]
                   + f_3 * pc_z[k] * lsg_255[k];

        t_360[k] = f_11 * ksg_258[k]
                   + f_6 * lsf0_173[k]
                   - f_7 * lsf1_173[k]
                   + f_3 * pc_x[k] * lsg_258[k];

        t_361[k] = f_11 * ksg_182[k]
                   + f_3 * pc_y[k] * lsg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, ksg_168, ksg_260, ksg_261, lsf0_175, \
                         lsf0_176, lsf1_175, lsf1_176, lsg_258, lsg_260, \
                         lsg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_11 * ksg_260[k]
                   + f_6 * lsf0_175[k]
                   - f_7 * lsf1_175[k]
                   + f_3 * pc_x[k] * lsg_260[k];

        t_363[k] = f_11 * ksg_261[k]
                   + f_4 * lsf0_176[k]
                   - f_5 * lsf1_176[k]
                   + f_3 * pc_x[k] * lsg_261[k];

        t_364[k] = f_10 * ksg_168[k]
                   + f_3 * pc_z[k] * lsg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, ksg_185, ksg_264, ksg_265, \
                         ksg_266, lsf0_179, lsf1_179, lsg_260, lsg_264, lsg_265, \
                         lsg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * ksg_185[k]
                   + f_3 * pc_y[k] * lsg_260[k];

        t_366[k] = f_11 * ksg_264[k]
                   + f_4 * lsf0_179[k]
                   - f_5 * lsf1_179[k]
                   + f_3 * pc_x[k] * lsg_264[k];

        t_367[k] = f_11 * ksg_265[k]
                   + f_3 * pc_x[k] * lsg_265[k];

        t_368[k] = f_11 * ksg_266[k]
                   + f_3 * pc_x[k] * lsg_266[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_294 = buffer.data(ksh0 + 294);
    const auto *ksh0_297 = buffer.data(ksh0 + 297);
    const auto *ksh0_299 = buffer.data(ksh0 + 299);
    const auto *ksh0_300 = buffer.data(ksh0 + 300);
    const auto *ksh0_303 = buffer.data(ksh0 + 303);
    const auto *ksh0_314 = buffer.data(ksh0 + 314);
    const auto *ksh0_315 = buffer.data(ksh0 + 315);
    const auto *ksh0_318 = buffer.data(ksh0 + 318);
    const auto *ksh0_321 = buffer.data(ksh0 + 321);
    const auto *ksh0_330 = buffer.data(ksh0 + 330);

    const auto *ksg_175 = buffer.data(ksg + 175);
    const auto *ksg_179 = buffer.data(ksg + 179);
    const auto *ksg_180 = buffer.data(ksg + 180);
    const auto *ksg_183 = buffer.data(ksg + 183);
    const auto *ksg_190 = buffer.data(ksg + 190);
    const auto *ksg_192 = buffer.data(ksg + 192);
    const auto *ksg_193 = buffer.data(ksg + 193);
    const auto *ksg_194 = buffer.data(ksg + 194);
    const auto *ksg_195 = buffer.data(ksg + 195);
    const auto *ksg_197 = buffer.data(ksg + 197);
    const auto *ksg_198 = buffer.data(ksg + 198);
    const auto *ksg_200 = buffer.data(ksg + 200);
    const auto *ksg_205 = buffer.data(ksg + 205);
    const auto *ksg_207 = buffer.data(ksg + 207);
    const auto *ksg_208 = buffer.data(ksg + 208);
    const auto *ksg_209 = buffer.data(ksg + 209);
    const auto *ksg_210 = buffer.data(ksg + 210);
    const auto *ksg_211 = buffer.data(ksg + 211);
    const auto *ksg_212 = buffer.data(ksg + 212);
    const auto *ksg_213 = buffer.data(ksg + 213);
    const auto *ksg_215 = buffer.data(ksg + 215);
    const auto *ksg_220 = buffer.data(ksg + 220);
    const auto *ksg_222 = buffer.data(ksg + 222);
    const auto *ksg_223 = buffer.data(ksg + 223);
    const auto *ksg_224 = buffer.data(ksg + 224);
    const auto *ksg_225 = buffer.data(ksg + 225);
    const auto *ksg_228 = buffer.data(ksg + 228);
    const auto *ksg_230 = buffer.data(ksg + 230);
    const auto *ksg_235 = buffer.data(ksg + 235);
    const auto *ksg_239 = buffer.data(ksg + 239);
    const auto *ksg_240 = buffer.data(ksg + 240);
    const auto *ksg_242 = buffer.data(ksg + 242);
    const auto *ksg_245 = buffer.data(ksg + 245);
    const auto *ksg_252 = buffer.data(ksg + 252);
    const auto *ksg_253 = buffer.data(ksg + 253);
    const auto *ksg_254 = buffer.data(ksg + 254);
    const auto *ksg_255 = buffer.data(ksg + 255);
    const auto *ksg_257 = buffer.data(ksg + 257);
    const auto *ksg_267 = buffer.data(ksg + 267);
    const auto *ksg_268 = buffer.data(ksg + 268);
    const auto *ksg_269 = buffer.data(ksg + 269);
    const auto *ksg_270 = buffer.data(ksg + 270);
    const auto *ksg_273 = buffer.data(ksg + 273);
    const auto *ksg_275 = buffer.data(ksg + 275);
    const auto *ksg_276 = buffer.data(ksg + 276);
    const auto *ksg_279 = buffer.data(ksg + 279);
    const auto *ksg_280 = buffer.data(ksg + 280);
    const auto *ksg_281 = buffer.data(ksg + 281);
    const auto *ksg_282 = buffer.data(ksg + 282);
    const auto *ksg_283 = buffer.data(ksg + 283);
    const auto *ksg_284 = buffer.data(ksg + 284);
    const auto *ksg_295 = buffer.data(ksg + 295);
    const auto *ksg_296 = buffer.data(ksg + 296);
    const auto *ksg_297 = buffer.data(ksg + 297);
    const auto *ksg_298 = buffer.data(ksg + 298);
    const auto *ksg_299 = buffer.data(ksg + 299);
    const auto *ksg_300 = buffer.data(ksg + 300);
    const auto *ksg_305 = buffer.data(ksg + 305);
    const auto *ksg_309 = buffer.data(ksg + 309);
    const auto *ksg_310 = buffer.data(ksg + 310);
    const auto *ksg_311 = buffer.data(ksg + 311);
    const auto *ksg_312 = buffer.data(ksg + 312);
    const auto *ksg_314 = buffer.data(ksg + 314);
    const auto *ksg_315 = buffer.data(ksg + 315);
    const auto *ksg_318 = buffer.data(ksg + 318);
    const auto *ksg_321 = buffer.data(ksg + 321);
    const auto *ksg_325 = buffer.data(ksg + 325);
    const auto *ksg_327 = buffer.data(ksg + 327);
    const auto *ksg_328 = buffer.data(ksg + 328);
    const auto *ksg_329 = buffer.data(ksg + 329);
    const auto *ksg_335 = buffer.data(ksg + 335);
    const auto *ksg_339 = buffer.data(ksg + 339);
    const auto *ksg_340 = buffer.data(ksg + 340);
    const auto *ksg_341 = buffer.data(ksg + 341);
    const auto *ksg_342 = buffer.data(ksg + 342);
    const auto *ksg_343 = buffer.data(ksg + 343);
    const auto *ksg_344 = buffer.data(ksg + 344);
    const auto *ksg_345 = buffer.data(ksg + 345);
    const auto *ksg_348 = buffer.data(ksg + 348);

    const auto *ksh1_294 = buffer.data(ksh1 + 294);
    const auto *ksh1_297 = buffer.data(ksh1 + 297);
    const auto *ksh1_299 = buffer.data(ksh1 + 299);
    const auto *ksh1_300 = buffer.data(ksh1 + 300);
    const auto *ksh1_303 = buffer.data(ksh1 + 303);
    const auto *ksh1_314 = buffer.data(ksh1 + 314);
    const auto *ksh1_315 = buffer.data(ksh1 + 315);
    const auto *ksh1_318 = buffer.data(ksh1 + 318);
    const auto *ksh1_321 = buffer.data(ksh1 + 321);
    const auto *ksh1_330 = buffer.data(ksh1 + 330);

    const auto *lsf0_176 = buffer.data(lsf0 + 176);
    const auto *lsf0_178 = buffer.data(lsf0 + 178);
    const auto *lsf0_179 = buffer.data(lsf0 + 179);
    const auto *lsf0_180 = buffer.data(lsf0 + 180);
    const auto *lsf0_183 = buffer.data(lsf0 + 183);
    const auto *lsf0_185 = buffer.data(lsf0 + 185);
    const auto *lsf0_186 = buffer.data(lsf0 + 186);
    const auto *lsf0_188 = buffer.data(lsf0 + 188);
    const auto *lsf0_189 = buffer.data(lsf0 + 189);
    const auto *lsf0_196 = buffer.data(lsf0 + 196);
    const auto *lsf0_198 = buffer.data(lsf0 + 198);
    const auto *lsf0_199 = buffer.data(lsf0 + 199);
    const auto *lsf0_200 = buffer.data(lsf0 + 200);
    const auto *lsf0_201 = buffer.data(lsf0 + 201);
    const auto *lsf0_202 = buffer.data(lsf0 + 202);
    const auto *lsf0_205 = buffer.data(lsf0 + 205);
    const auto *lsf0_206 = buffer.data(lsf0 + 206);
    const auto *lsf0_207 = buffer.data(lsf0 + 207);
    const auto *lsf0_208 = buffer.data(lsf0 + 208);
    const auto *lsf0_209 = buffer.data(lsf0 + 209);
    const auto *lsf0_210 = buffer.data(lsf0 + 210);
    const auto *lsf0_212 = buffer.data(lsf0 + 212);
    const auto *lsf0_213 = buffer.data(lsf0 + 213);
    const auto *lsf0_216 = buffer.data(lsf0 + 216);
    const auto *lsf0_217 = buffer.data(lsf0 + 217);
    const auto *lsf0_219 = buffer.data(lsf0 + 219);
    const auto *lsf0_225 = buffer.data(lsf0 + 225);
    const auto *lsf0_228 = buffer.data(lsf0 + 228);
    const auto *lsf0_229 = buffer.data(lsf0 + 229);
    const auto *lsf0_230 = buffer.data(lsf0 + 230);
    const auto *lsf0_233 = buffer.data(lsf0 + 233);

    const auto *lsf1_176 = buffer.data(lsf1 + 176);
    const auto *lsf1_178 = buffer.data(lsf1 + 178);
    const auto *lsf1_179 = buffer.data(lsf1 + 179);
    const auto *lsf1_180 = buffer.data(lsf1 + 180);
    const auto *lsf1_183 = buffer.data(lsf1 + 183);
    const auto *lsf1_185 = buffer.data(lsf1 + 185);
    const auto *lsf1_186 = buffer.data(lsf1 + 186);
    const auto *lsf1_188 = buffer.data(lsf1 + 188);
    const auto *lsf1_189 = buffer.data(lsf1 + 189);
    const auto *lsf1_196 = buffer.data(lsf1 + 196);
    const auto *lsf1_198 = buffer.data(lsf1 + 198);
    const auto *lsf1_199 = buffer.data(lsf1 + 199);
    const auto *lsf1_200 = buffer.data(lsf1 + 200);
    const auto *lsf1_201 = buffer.data(lsf1 + 201);
    const auto *lsf1_202 = buffer.data(lsf1 + 202);
    const auto *lsf1_205 = buffer.data(lsf1 + 205);
    const auto *lsf1_206 = buffer.data(lsf1 + 206);
    const auto *lsf1_207 = buffer.data(lsf1 + 207);
    const auto *lsf1_208 = buffer.data(lsf1 + 208);
    const auto *lsf1_209 = buffer.data(lsf1 + 209);
    const auto *lsf1_210 = buffer.data(lsf1 + 210);
    const auto *lsf1_212 = buffer.data(lsf1 + 212);
    const auto *lsf1_213 = buffer.data(lsf1 + 213);
    const auto *lsf1_216 = buffer.data(lsf1 + 216);
    const auto *lsf1_217 = buffer.data(lsf1 + 217);
    const auto *lsf1_219 = buffer.data(lsf1 + 219);
    const auto *lsf1_225 = buffer.data(lsf1 + 225);
    const auto *lsf1_228 = buffer.data(lsf1 + 228);
    const auto *lsf1_229 = buffer.data(lsf1 + 229);
    const auto *lsf1_230 = buffer.data(lsf1 + 230);
    const auto *lsf1_233 = buffer.data(lsf1 + 233);

    const auto *lsg_265 = buffer.data(lsg + 265);
    const auto *lsg_267 = buffer.data(lsg + 267);
    const auto *lsg_268 = buffer.data(lsg + 268);
    const auto *lsg_269 = buffer.data(lsg + 269);
    const auto *lsg_270 = buffer.data(lsg + 270);
    const auto *lsg_272 = buffer.data(lsg + 272);
    const auto *lsg_273 = buffer.data(lsg + 273);
    const auto *lsg_275 = buffer.data(lsg + 275);
    const auto *lsg_276 = buffer.data(lsg + 276);
    const auto *lsg_279 = buffer.data(lsg + 279);
    const auto *lsg_280 = buffer.data(lsg + 280);
    const auto *lsg_281 = buffer.data(lsg + 281);
    const auto *lsg_282 = buffer.data(lsg + 282);
    const auto *lsg_283 = buffer.data(lsg + 283);
    const auto *lsg_284 = buffer.data(lsg + 284);
    const auto *lsg_285 = buffer.data(lsg + 285);
    const auto *lsg_287 = buffer.data(lsg + 287);
    const auto *lsg_288 = buffer.data(lsg + 288);
    const auto *lsg_290 = buffer.data(lsg + 290);
    const auto *lsg_295 = buffer.data(lsg + 295);
    const auto *lsg_296 = buffer.data(lsg + 296);
    const auto *lsg_297 = buffer.data(lsg + 297);
    const auto *lsg_298 = buffer.data(lsg + 298);
    const auto *lsg_299 = buffer.data(lsg + 299);
    const auto *lsg_300 = buffer.data(lsg + 300);
    const auto *lsg_301 = buffer.data(lsg + 301);
    const auto *lsg_302 = buffer.data(lsg + 302);
    const auto *lsg_303 = buffer.data(lsg + 303);
    const auto *lsg_304 = buffer.data(lsg + 304);
    const auto *lsg_305 = buffer.data(lsg + 305);
    const auto *lsg_309 = buffer.data(lsg + 309);
    const auto *lsg_310 = buffer.data(lsg + 310);
    const auto *lsg_311 = buffer.data(lsg + 311);
    const auto *lsg_312 = buffer.data(lsg + 312);
    const auto *lsg_313 = buffer.data(lsg + 313);
    const auto *lsg_314 = buffer.data(lsg + 314);
    const auto *lsg_315 = buffer.data(lsg + 315);
    const auto *lsg_316 = buffer.data(lsg + 316);
    const auto *lsg_317 = buffer.data(lsg + 317);
    const auto *lsg_318 = buffer.data(lsg + 318);
    const auto *lsg_320 = buffer.data(lsg + 320);
    const auto *lsg_321 = buffer.data(lsg + 321);
    const auto *lsg_325 = buffer.data(lsg + 325);
    const auto *lsg_326 = buffer.data(lsg + 326);
    const auto *lsg_327 = buffer.data(lsg + 327);
    const auto *lsg_328 = buffer.data(lsg + 328);
    const auto *lsg_329 = buffer.data(lsg + 329);
    const auto *lsg_330 = buffer.data(lsg + 330);
    const auto *lsg_332 = buffer.data(lsg + 332);
    const auto *lsg_333 = buffer.data(lsg + 333);
    const auto *lsg_335 = buffer.data(lsg + 335);
    const auto *lsg_339 = buffer.data(lsg + 339);
    const auto *lsg_340 = buffer.data(lsg + 340);
    const auto *lsg_341 = buffer.data(lsg + 341);
    const auto *lsg_342 = buffer.data(lsg + 342);
    const auto *lsg_343 = buffer.data(lsg + 343);
    const auto *lsg_344 = buffer.data(lsg + 344);
    const auto *lsg_345 = buffer.data(lsg + 345);
    const auto *lsg_347 = buffer.data(lsg + 347);
    const auto *lsg_348 = buffer.data(lsg + 348);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, ksg_190, ksg_267, ksg_268, \
                         ksg_269, lsf0_176, lsf1_176, lsg_265, lsg_267, lsg_268, \
                         lsg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_11 * ksg_267[k]
                   + f_3 * pc_x[k] * lsg_267[k];

        t_370[k] = f_11 * ksg_268[k]
                   + f_3 * pc_x[k] * lsg_268[k];

        t_371[k] = f_11 * ksg_269[k]
                   + f_3 * pc_x[k] * lsg_269[k];

        t_372[k] = f_11 * ksg_190[k]
                   + f_1 * lsf0_176[k]
                   - f_2 * lsf1_176[k]
                   + f_3 * pc_y[k] * lsg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, ksg_175, ksg_192, ksg_193, lsf0_178, \
                         lsf0_179, lsf1_178, lsf1_179, lsg_265, lsg_267, \
                         lsg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * ksg_175[k]
                   + f_3 * pc_z[k] * lsg_265[k];

        t_374[k] = f_11 * ksg_192[k]
                   + f_6 * lsf0_178[k]
                   - f_7 * lsf1_178[k]
                   + f_3 * pc_y[k] * lsg_267[k];

        t_375[k] = f_11 * ksg_193[k]
                   + f_4 * lsf0_179[k]
                   - f_5 * lsf1_179[k]
                   + f_3 * pc_y[k] * lsg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, ksg_179, ksg_194, ksg_270, \
                         lsf0_179, lsf0_180, lsf1_179, lsf1_180, lsg_269, \
                         lsg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * ksg_194[k]
                   + f_3 * pc_y[k] * lsg_269[k];

        t_377[k] = f_10 * ksg_179[k]
                   + f_1 * lsf0_179[k]
                   - f_2 * lsf1_179[k]
                   + f_3 * pc_z[k] * lsg_269[k];

        t_378[k] = f_11 * ksg_270[k]
                   + f_1 * lsf0_180[k]
                   - f_2 * lsf1_180[k]
                   + f_3 * pc_x[k] * lsg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, ksg_180, ksg_195, \
                         ksg_197, ksg_273, lsf0_183, lsf1_183, lsg_270, lsg_272, \
                         lsg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * ksg_195[k]
                   + f_3 * pc_y[k] * lsg_270[k];

        t_380[k] = f_11 * ksg_180[k]
                   + f_3 * pc_z[k] * lsg_270[k];

        t_381[k] = f_11 * ksg_273[k]
                   + f_6 * lsf0_183[k]
                   - f_7 * lsf1_183[k]
                   + f_3 * pc_x[k] * lsg_273[k];

        t_382[k] = f_10 * ksg_197[k]
                   + f_3 * pc_y[k] * lsg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, ksg_183, ksg_275, ksg_276, lsf0_185, \
                         lsf0_186, lsf1_185, lsf1_186, lsg_273, lsg_275, \
                         lsg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_11 * ksg_275[k]
                   + f_6 * lsf0_185[k]
                   - f_7 * lsf1_185[k]
                   + f_3 * pc_x[k] * lsg_275[k];

        t_384[k] = f_11 * ksg_276[k]
                   + f_4 * lsf0_186[k]
                   - f_5 * lsf1_186[k]
                   + f_3 * pc_x[k] * lsg_276[k];

        t_385[k] = f_11 * ksg_183[k]
                   + f_3 * pc_z[k] * lsg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, ksg_200, ksg_279, ksg_280, \
                         ksg_281, lsf0_189, lsf1_189, lsg_275, lsg_279, lsg_280, \
                         lsg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * ksg_200[k]
                   + f_3 * pc_y[k] * lsg_275[k];

        t_387[k] = f_11 * ksg_279[k]
                   + f_4 * lsf0_189[k]
                   - f_5 * lsf1_189[k]
                   + f_3 * pc_x[k] * lsg_279[k];

        t_388[k] = f_11 * ksg_280[k]
                   + f_3 * pc_x[k] * lsg_280[k];

        t_389[k] = f_11 * ksg_281[k]
                   + f_3 * pc_x[k] * lsg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, ksg_205, ksg_282, ksg_283, \
                         ksg_284, lsf0_186, lsf1_186, lsg_280, lsg_282, lsg_283, \
                         lsg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_11 * ksg_282[k]
                   + f_3 * pc_x[k] * lsg_282[k];

        t_391[k] = f_11 * ksg_283[k]
                   + f_3 * pc_x[k] * lsg_283[k];

        t_392[k] = f_11 * ksg_284[k]
                   + f_3 * pc_x[k] * lsg_284[k];

        t_393[k] = f_10 * ksg_205[k]
                   + f_1 * lsf0_186[k]
                   - f_2 * lsf1_186[k]
                   + f_3 * pc_y[k] * lsg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, ksg_190, ksg_207, ksg_208, lsf0_188, \
                         lsf0_189, lsf1_188, lsf1_189, lsg_280, lsg_282, \
                         lsg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * ksg_190[k]
                   + f_3 * pc_z[k] * lsg_280[k];

        t_395[k] = f_10 * ksg_207[k]
                   + f_6 * lsf0_188[k]
                   - f_7 * lsf1_188[k]
                   + f_3 * pc_y[k] * lsg_282[k];

        t_396[k] = f_10 * ksg_208[k]
                   + f_4 * lsf0_189[k]
                   - f_5 * lsf1_189[k]
                   + f_3 * pc_y[k] * lsg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, ksh0_294, ksg_194, \
                         ksg_209, ksg_210, ksh1_294, lsf0_189, lsf1_189, lsg_284, \
                         lsg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * ksg_209[k]
                   + f_3 * pc_y[k] * lsg_284[k];

        t_398[k] = f_11 * ksg_194[k]
                   + f_1 * lsf0_189[k]
                   - f_2 * lsf1_189[k]
                   + f_3 * pc_z[k] * lsg_284[k];

        t_399[k] = pa_y[k] * ksh0_294[k]
                   - f_8 * pc_y[k] * ksh1_294[k];

        t_400[k] = f_9 * ksg_210[k]
                   + f_3 * pc_y[k] * lsg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, ksh0_297, ksh0_299, \
                         ksg_195, ksg_211, ksg_212, ksh1_297, ksh1_299, lsg_285, \
                         lsg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_17 * ksg_195[k]
                   + f_3 * pc_z[k] * lsg_285[k];

        t_402[k] = pa_y[k] * ksh0_297[k]
                   + f_10 * ksg_211[k]
                   - f_8 * pc_y[k] * ksh1_297[k];

        t_403[k] = f_9 * ksg_212[k]
                   + f_3 * pc_y[k] * lsg_287[k];

        t_404[k] = pa_y[k] * ksh0_299[k]
                   - f_8 * pc_y[k] * ksh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, ksh0_300, ksh0_303, \
                         ksg_198, ksg_213, ksg_215, ksh1_300, ksh1_303, lsg_288, \
                         lsg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * ksh0_300[k]
                   + f_11 * ksg_213[k]
                   - f_8 * pc_y[k] * ksh1_300[k];

        t_406[k] = f_17 * ksg_198[k]
                   + f_3 * pc_z[k] * lsg_288[k];

        t_407[k] = f_9 * ksg_215[k]
                   + f_3 * pc_y[k] * lsg_290[k];

        t_408[k] = pa_y[k] * ksh0_303[k]
                   - f_8 * pc_y[k] * ksh1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, ksg_295, ksg_296, ksg_297, \
                         ksg_298, ksg_299, lsg_295, lsg_296, lsg_297, lsg_298, \
                         lsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_11 * ksg_295[k]
                   + f_3 * pc_x[k] * lsg_295[k];

        t_410[k] = f_11 * ksg_296[k]
                   + f_3 * pc_x[k] * lsg_296[k];

        t_411[k] = f_11 * ksg_297[k]
                   + f_3 * pc_x[k] * lsg_297[k];

        t_412[k] = f_11 * ksg_298[k]
                   + f_3 * pc_x[k] * lsg_298[k];

        t_413[k] = f_11 * ksg_299[k]
                   + f_3 * pc_x[k] * lsg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, ksg_205, ksg_220, ksg_222, lsf0_196, \
                         lsf0_198, lsf1_196, lsf1_198, lsg_295, \
                         lsg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * ksg_220[k]
                   + f_1 * lsf0_196[k]
                   - f_2 * lsf1_196[k]
                   + f_3 * pc_y[k] * lsg_295[k];

        t_415[k] = f_17 * ksg_205[k]
                   + f_3 * pc_z[k] * lsg_295[k];

        t_416[k] = f_9 * ksg_222[k]
                   + f_6 * lsf0_198[k]
                   - f_7 * lsf1_198[k]
                   + f_3 * pc_y[k] * lsg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, ksh0_314, ksg_223, ksg_224, \
                         ksh1_314, lsf0_199, lsf1_199, lsg_298, \
                         lsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * ksg_223[k]
                   + f_4 * lsf0_199[k]
                   - f_5 * lsf1_199[k]
                   + f_3 * pc_y[k] * lsg_298[k];

        t_418[k] = f_9 * ksg_224[k]
                   + f_3 * pc_y[k] * lsg_299[k];

        t_419[k] = pa_y[k] * ksh0_314[k]
                   - f_8 * pc_y[k] * ksh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, ksg_210, \
                         ksg_300, lsf0_200, lsf1_200, lsg_300, lsg_301, \
                         lsg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_11 * ksg_300[k]
                   + f_1 * lsf0_200[k]
                   - f_2 * lsf1_200[k]
                   + f_3 * pc_x[k] * lsg_300[k];

        t_421[k] = f_3 * pc_y[k] * lsg_300[k];

        t_422[k] = f_16 * ksg_210[k]
                   + f_3 * pc_z[k] * lsg_300[k];

        t_423[k] = f_4 * lsf0_200[k]
                   - f_5 * lsf1_200[k]
                   + f_3 * pc_y[k] * lsg_301[k];

        t_424[k] = f_3 * pc_y[k] * lsg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, ksg_305, lsf0_201, lsf0_202, \
                         lsf0_205, lsf1_201, lsf1_202, lsf1_205, lsg_303, lsg_304, \
                         lsg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_11 * ksg_305[k]
                   + f_6 * lsf0_205[k]
                   - f_7 * lsf1_205[k]
                   + f_3 * pc_x[k] * lsg_305[k];

        t_426[k] = f_6 * lsf0_201[k]
                   - f_7 * lsf1_201[k]
                   + f_3 * pc_y[k] * lsg_303[k];

        t_427[k] = f_4 * lsf0_202[k]
                   - f_5 * lsf1_202[k]
                   + f_3 * pc_y[k] * lsg_304[k];

        t_428[k] = f_3 * pc_y[k] * lsg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, ksg_309, ksg_310, ksg_311, ksg_312, \
                         lsf0_209, lsf1_209, lsg_309, lsg_310, lsg_311, \
                         lsg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_11 * ksg_309[k]
                   + f_4 * lsf0_209[k]
                   - f_5 * lsf1_209[k]
                   + f_3 * pc_x[k] * lsg_309[k];

        t_430[k] = f_11 * ksg_310[k]
                   + f_3 * pc_x[k] * lsg_310[k];

        t_431[k] = f_11 * ksg_311[k]
                   + f_3 * pc_x[k] * lsg_311[k];

        t_432[k] = f_11 * ksg_312[k]
                   + f_3 * pc_x[k] * lsg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, ksg_314, lsf0_206, lsf0_207, \
                         lsf1_206, lsf1_207, lsg_309, lsg_310, lsg_311, \
                         lsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * lsg_309[k];

        t_434[k] = f_11 * ksg_314[k]
                   + f_3 * pc_x[k] * lsg_314[k];

        t_435[k] = f_1 * lsf0_206[k]
                   - f_2 * lsf1_206[k]
                   + f_3 * pc_y[k] * lsg_310[k];

        t_436[k] = f_13 * lsf0_207[k]
                   - f_14 * lsf1_207[k]
                   + f_3 * pc_y[k] * lsg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, ksg_224, lsf0_208, lsf0_209, \
                         lsf1_208, lsf1_209, lsg_312, lsg_313, \
                         lsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * lsf0_208[k]
                   - f_7 * lsf1_208[k]
                   + f_3 * pc_y[k] * lsg_312[k];

        t_438[k] = f_4 * lsf0_209[k]
                   - f_5 * lsf1_209[k]
                   + f_3 * pc_y[k] * lsg_313[k];

        t_439[k] = f_3 * pc_y[k] * lsg_314[k];

        t_440[k] = f_16 * ksg_224[k]
                   + f_1 * lsf0_209[k]
                   - f_2 * lsf1_209[k]
                   + f_3 * pc_z[k] * lsg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, ksg_225, ksg_315, \
                         ksg_318, lsf0_210, lsf0_213, lsf1_210, lsf1_213, lsg_315, \
                         lsg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_10 * ksg_315[k]
                   + f_1 * lsf0_210[k]
                   - f_2 * lsf1_210[k]
                   + f_3 * pc_x[k] * lsg_315[k];

        t_442[k] = f_15 * ksg_225[k]
                   + f_3 * pc_y[k] * lsg_315[k];

        t_443[k] = f_3 * pc_z[k] * lsg_315[k];

        t_444[k] = f_10 * ksg_318[k]
                   + f_6 * lsf0_213[k]
                   - f_7 * lsf1_213[k]
                   + f_3 * pc_x[k] * lsg_318[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_z, ksg_321, lsf0_210, lsf0_216, \
                         lsf1_210, lsf1_216, lsg_316, lsg_317, lsg_318, \
                         lsg_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * lsg_316[k];

        t_446[k] = f_4 * lsf0_210[k]
                   - f_5 * lsf1_210[k]
                   + f_3 * pc_z[k] * lsg_317[k];

        t_447[k] = f_10 * ksg_321[k]
                   + f_4 * lsf0_216[k]
                   - f_5 * lsf1_216[k]
                   + f_3 * pc_x[k] * lsg_321[k];

        t_448[k] = f_3 * pc_z[k] * lsg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, ksg_230, ksg_325, \
                         lsf0_212, lsf1_212, lsg_320, lsg_321, \
                         lsg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_15 * ksg_230[k]
                   + f_3 * pc_y[k] * lsg_320[k];

        t_450[k] = f_6 * lsf0_212[k]
                   - f_7 * lsf1_212[k]
                   + f_3 * pc_z[k] * lsg_320[k];

        t_451[k] = f_10 * ksg_325[k]
                   + f_3 * pc_x[k] * lsg_325[k];

        t_452[k] = f_3 * pc_z[k] * lsg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, ksg_235, ksg_327, ksg_328, \
                         ksg_329, lsf0_216, lsf1_216, lsg_325, lsg_327, lsg_328, \
                         lsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_10 * ksg_327[k]
                   + f_3 * pc_x[k] * lsg_327[k];

        t_454[k] = f_10 * ksg_328[k]
                   + f_3 * pc_x[k] * lsg_328[k];

        t_455[k] = f_10 * ksg_329[k]
                   + f_3 * pc_x[k] * lsg_329[k];

        t_456[k] = f_15 * ksg_235[k]
                   + f_1 * lsf0_216[k]
                   - f_2 * lsf1_216[k]
                   + f_3 * pc_y[k] * lsg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, ksg_239, lsf0_216, lsf0_217, \
                         lsf1_216, lsf1_217, lsg_325, lsg_326, lsg_327, \
                         lsg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * lsg_325[k];

        t_458[k] = f_4 * lsf0_216[k]
                   - f_5 * lsf1_216[k]
                   + f_3 * pc_z[k] * lsg_326[k];

        t_459[k] = f_6 * lsf0_217[k]
                   - f_7 * lsf1_217[k]
                   + f_3 * pc_z[k] * lsg_327[k];

        t_460[k] = f_15 * ksg_239[k]
                   + f_3 * pc_y[k] * lsg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_z, pc_y, pc_z, ksh0_315, ksg_225, \
                         ksg_240, ksh1_315, lsf0_219, lsf1_219, lsg_329, \
                         lsg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * lsf0_219[k]
                   - f_2 * lsf1_219[k]
                   + f_3 * pc_z[k] * lsg_329[k];

        t_462[k] = pa_z[k] * ksh0_315[k]
                   - f_8 * pc_z[k] * ksh1_315[k];

        t_463[k] = f_16 * ksg_240[k]
                   + f_3 * pc_y[k] * lsg_330[k];

        t_464[k] = f_9 * ksg_225[k]
                   + f_3 * pc_z[k] * lsg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_z, pc_x, pc_y, pc_z, ksh0_318, ksg_242, \
                         ksg_335, ksh1_318, lsf0_225, lsf1_225, lsg_332, \
                         lsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * ksh0_318[k]
                   - f_8 * pc_z[k] * ksh1_318[k];

        t_466[k] = f_16 * ksg_242[k]
                   + f_3 * pc_y[k] * lsg_332[k];

        t_467[k] = f_10 * ksg_335[k]
                   + f_6 * lsf0_225[k]
                   - f_7 * lsf1_225[k]
                   + f_3 * pc_x[k] * lsg_335[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, ksh0_321, ksg_228, ksg_245, \
                         ksh1_321, lsg_333, lsg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * ksh0_321[k]
                   - f_8 * pc_z[k] * ksh1_321[k];

        t_469[k] = f_9 * ksg_228[k]
                   + f_3 * pc_z[k] * lsg_333[k];

        t_470[k] = f_16 * ksg_245[k]
                   + f_3 * pc_y[k] * lsg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, ksg_339, ksg_340, ksg_341, ksg_342, \
                         lsf0_229, lsf1_229, lsg_339, lsg_340, lsg_341, \
                         lsg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_10 * ksg_339[k]
                   + f_4 * lsf0_229[k]
                   - f_5 * lsf1_229[k]
                   + f_3 * pc_x[k] * lsg_339[k];

        t_472[k] = f_10 * ksg_340[k]
                   + f_3 * pc_x[k] * lsg_340[k];

        t_473[k] = f_10 * ksg_341[k]
                   + f_3 * pc_x[k] * lsg_341[k];

        t_474[k] = f_10 * ksg_342[k]
                   + f_3 * pc_x[k] * lsg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, ksh0_330, ksg_235, \
                         ksg_343, ksg_344, ksh1_330, lsg_340, lsg_343, \
                         lsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_10 * ksg_343[k]
                   + f_3 * pc_x[k] * lsg_343[k];

        t_476[k] = f_10 * ksg_344[k]
                   + f_3 * pc_x[k] * lsg_344[k];

        t_477[k] = pa_z[k] * ksh0_330[k]
                   - f_8 * pc_z[k] * ksh1_330[k];

        t_478[k] = f_9 * ksg_235[k]
                   + f_3 * pc_z[k] * lsg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pc_y, ksg_252, ksg_253, ksg_254, lsf0_228, \
                         lsf0_229, lsf1_228, lsf1_229, lsg_342, lsg_343, \
                         lsg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_16 * ksg_252[k]
                   + f_6 * lsf0_228[k]
                   - f_7 * lsf1_228[k]
                   + f_3 * pc_y[k] * lsg_342[k];

        t_480[k] = f_16 * ksg_253[k]
                   + f_4 * lsf0_229[k]
                   - f_5 * lsf1_229[k]
                   + f_3 * pc_y[k] * lsg_343[k];

        t_481[k] = f_16 * ksg_254[k]
                   + f_3 * pc_y[k] * lsg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, ksg_239, ksg_255, ksg_345, \
                         lsf0_229, lsf0_230, lsf1_229, lsf1_230, lsg_344, \
                         lsg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * ksg_239[k]
                   + f_1 * lsf0_229[k]
                   - f_2 * lsf1_229[k]
                   + f_3 * pc_z[k] * lsg_344[k];

        t_483[k] = f_10 * ksg_345[k]
                   + f_1 * lsf0_230[k]
                   - f_2 * lsf1_230[k]
                   + f_3 * pc_x[k] * lsg_345[k];

        t_484[k] = f_17 * ksg_255[k]
                   + f_3 * pc_y[k] * lsg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, ksg_240, ksg_257, ksg_348, \
                         lsf0_233, lsf1_233, lsg_345, lsg_347, \
                         lsg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * ksg_240[k]
                   + f_3 * pc_z[k] * lsg_345[k];

        t_486[k] = f_10 * ksg_348[k]
                   + f_6 * lsf0_233[k]
                   - f_7 * lsf1_233[k]
                   + f_3 * pc_x[k] * lsg_348[k];

        t_487[k] = f_17 * ksg_257[k]
                   + f_3 * pc_y[k] * lsg_347[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_420 = buffer.data(ksh0 + 420);
    const auto *ksh0_423 = buffer.data(ksh0 + 423);
    const auto *ksh0_425 = buffer.data(ksh0 + 425);
    const auto *ksh0_426 = buffer.data(ksh0 + 426);
    const auto *ksh0_429 = buffer.data(ksh0 + 429);
    const auto *ksh0_440 = buffer.data(ksh0 + 440);
    const auto *ksh0_588 = buffer.data(ksh0 + 588);
    const auto *ksh0_591 = buffer.data(ksh0 + 591);
    const auto *ksh0_594 = buffer.data(ksh0 + 594);
    const auto *ksh0_603 = buffer.data(ksh0 + 603);

    const auto *ksg_243 = buffer.data(ksg + 243);
    const auto *ksg_250 = buffer.data(ksg + 250);
    const auto *ksg_254 = buffer.data(ksg + 254);
    const auto *ksg_255 = buffer.data(ksg + 255);
    const auto *ksg_258 = buffer.data(ksg + 258);
    const auto *ksg_260 = buffer.data(ksg + 260);
    const auto *ksg_265 = buffer.data(ksg + 265);
    const auto *ksg_267 = buffer.data(ksg + 267);
    const auto *ksg_268 = buffer.data(ksg + 268);
    const auto *ksg_269 = buffer.data(ksg + 269);
    const auto *ksg_270 = buffer.data(ksg + 270);
    const auto *ksg_272 = buffer.data(ksg + 272);
    const auto *ksg_273 = buffer.data(ksg + 273);
    const auto *ksg_275 = buffer.data(ksg + 275);
    const auto *ksg_280 = buffer.data(ksg + 280);
    const auto *ksg_282 = buffer.data(ksg + 282);
    const auto *ksg_283 = buffer.data(ksg + 283);
    const auto *ksg_284 = buffer.data(ksg + 284);
    const auto *ksg_285 = buffer.data(ksg + 285);
    const auto *ksg_287 = buffer.data(ksg + 287);
    const auto *ksg_288 = buffer.data(ksg + 288);
    const auto *ksg_290 = buffer.data(ksg + 290);
    const auto *ksg_295 = buffer.data(ksg + 295);
    const auto *ksg_297 = buffer.data(ksg + 297);
    const auto *ksg_298 = buffer.data(ksg + 298);
    const auto *ksg_299 = buffer.data(ksg + 299);
    const auto *ksg_300 = buffer.data(ksg + 300);
    const auto *ksg_301 = buffer.data(ksg + 301);
    const auto *ksg_302 = buffer.data(ksg + 302);
    const auto *ksg_303 = buffer.data(ksg + 303);
    const auto *ksg_305 = buffer.data(ksg + 305);
    const auto *ksg_310 = buffer.data(ksg + 310);
    const auto *ksg_312 = buffer.data(ksg + 312);
    const auto *ksg_313 = buffer.data(ksg + 313);
    const auto *ksg_314 = buffer.data(ksg + 314);
    const auto *ksg_315 = buffer.data(ksg + 315);
    const auto *ksg_320 = buffer.data(ksg + 320);
    const auto *ksg_350 = buffer.data(ksg + 350);
    const auto *ksg_351 = buffer.data(ksg + 351);
    const auto *ksg_354 = buffer.data(ksg + 354);
    const auto *ksg_355 = buffer.data(ksg + 355);
    const auto *ksg_356 = buffer.data(ksg + 356);
    const auto *ksg_357 = buffer.data(ksg + 357);
    const auto *ksg_358 = buffer.data(ksg + 358);
    const auto *ksg_359 = buffer.data(ksg + 359);
    const auto *ksg_360 = buffer.data(ksg + 360);
    const auto *ksg_363 = buffer.data(ksg + 363);
    const auto *ksg_365 = buffer.data(ksg + 365);
    const auto *ksg_366 = buffer.data(ksg + 366);
    const auto *ksg_369 = buffer.data(ksg + 369);
    const auto *ksg_370 = buffer.data(ksg + 370);
    const auto *ksg_371 = buffer.data(ksg + 371);
    const auto *ksg_372 = buffer.data(ksg + 372);
    const auto *ksg_373 = buffer.data(ksg + 373);
    const auto *ksg_374 = buffer.data(ksg + 374);
    const auto *ksg_375 = buffer.data(ksg + 375);
    const auto *ksg_378 = buffer.data(ksg + 378);
    const auto *ksg_380 = buffer.data(ksg + 380);
    const auto *ksg_381 = buffer.data(ksg + 381);
    const auto *ksg_384 = buffer.data(ksg + 384);
    const auto *ksg_385 = buffer.data(ksg + 385);
    const auto *ksg_386 = buffer.data(ksg + 386);
    const auto *ksg_387 = buffer.data(ksg + 387);
    const auto *ksg_388 = buffer.data(ksg + 388);
    const auto *ksg_389 = buffer.data(ksg + 389);
    const auto *ksg_400 = buffer.data(ksg + 400);
    const auto *ksg_401 = buffer.data(ksg + 401);
    const auto *ksg_402 = buffer.data(ksg + 402);
    const auto *ksg_403 = buffer.data(ksg + 403);
    const auto *ksg_404 = buffer.data(ksg + 404);
    const auto *ksg_405 = buffer.data(ksg + 405);
    const auto *ksg_410 = buffer.data(ksg + 410);
    const auto *ksg_414 = buffer.data(ksg + 414);
    const auto *ksg_415 = buffer.data(ksg + 415);
    const auto *ksg_416 = buffer.data(ksg + 416);
    const auto *ksg_417 = buffer.data(ksg + 417);
    const auto *ksg_419 = buffer.data(ksg + 419);
    const auto *ksg_420 = buffer.data(ksg + 420);
    const auto *ksg_423 = buffer.data(ksg + 423);
    const auto *ksg_426 = buffer.data(ksg + 426);
    const auto *ksg_430 = buffer.data(ksg + 430);
    const auto *ksg_432 = buffer.data(ksg + 432);
    const auto *ksg_433 = buffer.data(ksg + 433);
    const auto *ksg_434 = buffer.data(ksg + 434);

    const auto *ksh1_420 = buffer.data(ksh1 + 420);
    const auto *ksh1_423 = buffer.data(ksh1 + 423);
    const auto *ksh1_425 = buffer.data(ksh1 + 425);
    const auto *ksh1_426 = buffer.data(ksh1 + 426);
    const auto *ksh1_429 = buffer.data(ksh1 + 429);
    const auto *ksh1_440 = buffer.data(ksh1 + 440);
    const auto *ksh1_588 = buffer.data(ksh1 + 588);
    const auto *ksh1_591 = buffer.data(ksh1 + 591);
    const auto *ksh1_594 = buffer.data(ksh1 + 594);
    const auto *ksh1_603 = buffer.data(ksh1 + 603);

    const auto *lsf0_235 = buffer.data(lsf0 + 235);
    const auto *lsf0_236 = buffer.data(lsf0 + 236);
    const auto *lsf0_238 = buffer.data(lsf0 + 238);
    const auto *lsf0_239 = buffer.data(lsf0 + 239);
    const auto *lsf0_240 = buffer.data(lsf0 + 240);
    const auto *lsf0_243 = buffer.data(lsf0 + 243);
    const auto *lsf0_245 = buffer.data(lsf0 + 245);
    const auto *lsf0_246 = buffer.data(lsf0 + 246);
    const auto *lsf0_248 = buffer.data(lsf0 + 248);
    const auto *lsf0_249 = buffer.data(lsf0 + 249);
    const auto *lsf0_250 = buffer.data(lsf0 + 250);
    const auto *lsf0_253 = buffer.data(lsf0 + 253);
    const auto *lsf0_255 = buffer.data(lsf0 + 255);
    const auto *lsf0_256 = buffer.data(lsf0 + 256);
    const auto *lsf0_258 = buffer.data(lsf0 + 258);
    const auto *lsf0_259 = buffer.data(lsf0 + 259);
    const auto *lsf0_266 = buffer.data(lsf0 + 266);
    const auto *lsf0_268 = buffer.data(lsf0 + 268);
    const auto *lsf0_269 = buffer.data(lsf0 + 269);
    const auto *lsf0_270 = buffer.data(lsf0 + 270);
    const auto *lsf0_271 = buffer.data(lsf0 + 271);
    const auto *lsf0_272 = buffer.data(lsf0 + 272);
    const auto *lsf0_275 = buffer.data(lsf0 + 275);
    const auto *lsf0_276 = buffer.data(lsf0 + 276);
    const auto *lsf0_277 = buffer.data(lsf0 + 277);
    const auto *lsf0_278 = buffer.data(lsf0 + 278);
    const auto *lsf0_279 = buffer.data(lsf0 + 279);
    const auto *lsf0_280 = buffer.data(lsf0 + 280);
    const auto *lsf0_282 = buffer.data(lsf0 + 282);

    const auto *lsf1_235 = buffer.data(lsf1 + 235);
    const auto *lsf1_236 = buffer.data(lsf1 + 236);
    const auto *lsf1_238 = buffer.data(lsf1 + 238);
    const auto *lsf1_239 = buffer.data(lsf1 + 239);
    const auto *lsf1_240 = buffer.data(lsf1 + 240);
    const auto *lsf1_243 = buffer.data(lsf1 + 243);
    const auto *lsf1_245 = buffer.data(lsf1 + 245);
    const auto *lsf1_246 = buffer.data(lsf1 + 246);
    const auto *lsf1_248 = buffer.data(lsf1 + 248);
    const auto *lsf1_249 = buffer.data(lsf1 + 249);
    const auto *lsf1_250 = buffer.data(lsf1 + 250);
    const auto *lsf1_253 = buffer.data(lsf1 + 253);
    const auto *lsf1_255 = buffer.data(lsf1 + 255);
    const auto *lsf1_256 = buffer.data(lsf1 + 256);
    const auto *lsf1_258 = buffer.data(lsf1 + 258);
    const auto *lsf1_259 = buffer.data(lsf1 + 259);
    const auto *lsf1_266 = buffer.data(lsf1 + 266);
    const auto *lsf1_268 = buffer.data(lsf1 + 268);
    const auto *lsf1_269 = buffer.data(lsf1 + 269);
    const auto *lsf1_270 = buffer.data(lsf1 + 270);
    const auto *lsf1_271 = buffer.data(lsf1 + 271);
    const auto *lsf1_272 = buffer.data(lsf1 + 272);
    const auto *lsf1_275 = buffer.data(lsf1 + 275);
    const auto *lsf1_276 = buffer.data(lsf1 + 276);
    const auto *lsf1_277 = buffer.data(lsf1 + 277);
    const auto *lsf1_278 = buffer.data(lsf1 + 278);
    const auto *lsf1_279 = buffer.data(lsf1 + 279);
    const auto *lsf1_280 = buffer.data(lsf1 + 280);
    const auto *lsf1_282 = buffer.data(lsf1 + 282);

    const auto *lsg_348 = buffer.data(lsg + 348);
    const auto *lsg_350 = buffer.data(lsg + 350);
    const auto *lsg_351 = buffer.data(lsg + 351);
    const auto *lsg_354 = buffer.data(lsg + 354);
    const auto *lsg_355 = buffer.data(lsg + 355);
    const auto *lsg_356 = buffer.data(lsg + 356);
    const auto *lsg_357 = buffer.data(lsg + 357);
    const auto *lsg_358 = buffer.data(lsg + 358);
    const auto *lsg_359 = buffer.data(lsg + 359);
    const auto *lsg_360 = buffer.data(lsg + 360);
    const auto *lsg_362 = buffer.data(lsg + 362);
    const auto *lsg_363 = buffer.data(lsg + 363);
    const auto *lsg_365 = buffer.data(lsg + 365);
    const auto *lsg_366 = buffer.data(lsg + 366);
    const auto *lsg_369 = buffer.data(lsg + 369);
    const auto *lsg_370 = buffer.data(lsg + 370);
    const auto *lsg_371 = buffer.data(lsg + 371);
    const auto *lsg_372 = buffer.data(lsg + 372);
    const auto *lsg_373 = buffer.data(lsg + 373);
    const auto *lsg_374 = buffer.data(lsg + 374);
    const auto *lsg_375 = buffer.data(lsg + 375);
    const auto *lsg_377 = buffer.data(lsg + 377);
    const auto *lsg_378 = buffer.data(lsg + 378);
    const auto *lsg_380 = buffer.data(lsg + 380);
    const auto *lsg_381 = buffer.data(lsg + 381);
    const auto *lsg_384 = buffer.data(lsg + 384);
    const auto *lsg_385 = buffer.data(lsg + 385);
    const auto *lsg_386 = buffer.data(lsg + 386);
    const auto *lsg_387 = buffer.data(lsg + 387);
    const auto *lsg_388 = buffer.data(lsg + 388);
    const auto *lsg_389 = buffer.data(lsg + 389);
    const auto *lsg_390 = buffer.data(lsg + 390);
    const auto *lsg_392 = buffer.data(lsg + 392);
    const auto *lsg_393 = buffer.data(lsg + 393);
    const auto *lsg_395 = buffer.data(lsg + 395);
    const auto *lsg_400 = buffer.data(lsg + 400);
    const auto *lsg_401 = buffer.data(lsg + 401);
    const auto *lsg_402 = buffer.data(lsg + 402);
    const auto *lsg_403 = buffer.data(lsg + 403);
    const auto *lsg_404 = buffer.data(lsg + 404);
    const auto *lsg_405 = buffer.data(lsg + 405);
    const auto *lsg_406 = buffer.data(lsg + 406);
    const auto *lsg_407 = buffer.data(lsg + 407);
    const auto *lsg_408 = buffer.data(lsg + 408);
    const auto *lsg_409 = buffer.data(lsg + 409);
    const auto *lsg_410 = buffer.data(lsg + 410);
    const auto *lsg_414 = buffer.data(lsg + 414);
    const auto *lsg_415 = buffer.data(lsg + 415);
    const auto *lsg_416 = buffer.data(lsg + 416);
    const auto *lsg_417 = buffer.data(lsg + 417);
    const auto *lsg_418 = buffer.data(lsg + 418);
    const auto *lsg_419 = buffer.data(lsg + 419);
    const auto *lsg_420 = buffer.data(lsg + 420);
    const auto *lsg_421 = buffer.data(lsg + 421);
    const auto *lsg_422 = buffer.data(lsg + 422);
    const auto *lsg_423 = buffer.data(lsg + 423);
    const auto *lsg_425 = buffer.data(lsg + 425);
    const auto *lsg_426 = buffer.data(lsg + 426);
    const auto *lsg_430 = buffer.data(lsg + 430);
    const auto *lsg_432 = buffer.data(lsg + 432);
    const auto *lsg_433 = buffer.data(lsg + 433);
    const auto *lsg_434 = buffer.data(lsg + 434);

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_z, ksg_243, ksg_350, ksg_351, lsf0_235, \
                         lsf0_236, lsf1_235, lsf1_236, lsg_348, lsg_350, \
                         lsg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_10 * ksg_350[k]
                   + f_6 * lsf0_235[k]
                   - f_7 * lsf1_235[k]
                   + f_3 * pc_x[k] * lsg_350[k];

        t_489[k] = f_10 * ksg_351[k]
                   + f_4 * lsf0_236[k]
                   - f_5 * lsf1_236[k]
                   + f_3 * pc_x[k] * lsg_351[k];

        t_490[k] = f_10 * ksg_243[k]
                   + f_3 * pc_z[k] * lsg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pc_x, pc_y, ksg_260, ksg_354, ksg_355, \
                         ksg_356, lsf0_239, lsf1_239, lsg_350, lsg_354, lsg_355, \
                         lsg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_17 * ksg_260[k]
                   + f_3 * pc_y[k] * lsg_350[k];

        t_492[k] = f_10 * ksg_354[k]
                   + f_4 * lsf0_239[k]
                   - f_5 * lsf1_239[k]
                   + f_3 * pc_x[k] * lsg_354[k];

        t_493[k] = f_10 * ksg_355[k]
                   + f_3 * pc_x[k] * lsg_355[k];

        t_494[k] = f_10 * ksg_356[k]
                   + f_3 * pc_x[k] * lsg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, ksg_265, ksg_357, ksg_358, \
                         ksg_359, lsf0_236, lsf1_236, lsg_355, lsg_357, lsg_358, \
                         lsg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_10 * ksg_357[k]
                   + f_3 * pc_x[k] * lsg_357[k];

        t_496[k] = f_10 * ksg_358[k]
                   + f_3 * pc_x[k] * lsg_358[k];

        t_497[k] = f_10 * ksg_359[k]
                   + f_3 * pc_x[k] * lsg_359[k];

        t_498[k] = f_17 * ksg_265[k]
                   + f_1 * lsf0_236[k]
                   - f_2 * lsf1_236[k]
                   + f_3 * pc_y[k] * lsg_355[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, pc_z, ksg_250, ksg_267, ksg_268, lsf0_238, \
                         lsf0_239, lsf1_238, lsf1_239, lsg_355, lsg_357, \
                         lsg_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_10 * ksg_250[k]
                   + f_3 * pc_z[k] * lsg_355[k];

        t_500[k] = f_17 * ksg_267[k]
                   + f_6 * lsf0_238[k]
                   - f_7 * lsf1_238[k]
                   + f_3 * pc_y[k] * lsg_357[k];

        t_501[k] = f_17 * ksg_268[k]
                   + f_4 * lsf0_239[k]
                   - f_5 * lsf1_239[k]
                   + f_3 * pc_y[k] * lsg_358[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, ksg_254, ksg_269, ksg_360, \
                         lsf0_239, lsf0_240, lsf1_239, lsf1_240, lsg_359, \
                         lsg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_17 * ksg_269[k]
                   + f_3 * pc_y[k] * lsg_359[k];

        t_503[k] = f_10 * ksg_254[k]
                   + f_1 * lsf0_239[k]
                   - f_2 * lsf1_239[k]
                   + f_3 * pc_z[k] * lsg_359[k];

        t_504[k] = f_10 * ksg_360[k]
                   + f_1 * lsf0_240[k]
                   - f_2 * lsf1_240[k]
                   + f_3 * pc_x[k] * lsg_360[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, ksg_255, ksg_270, \
                         ksg_272, ksg_363, lsf0_243, lsf1_243, lsg_360, lsg_362, \
                         lsg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_11 * ksg_270[k]
                   + f_3 * pc_y[k] * lsg_360[k];

        t_506[k] = f_11 * ksg_255[k]
                   + f_3 * pc_z[k] * lsg_360[k];

        t_507[k] = f_10 * ksg_363[k]
                   + f_6 * lsf0_243[k]
                   - f_7 * lsf1_243[k]
                   + f_3 * pc_x[k] * lsg_363[k];

        t_508[k] = f_11 * ksg_272[k]
                   + f_3 * pc_y[k] * lsg_362[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, ksg_258, ksg_365, ksg_366, lsf0_245, \
                         lsf0_246, lsf1_245, lsf1_246, lsg_363, lsg_365, \
                         lsg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_10 * ksg_365[k]
                   + f_6 * lsf0_245[k]
                   - f_7 * lsf1_245[k]
                   + f_3 * pc_x[k] * lsg_365[k];

        t_510[k] = f_10 * ksg_366[k]
                   + f_4 * lsf0_246[k]
                   - f_5 * lsf1_246[k]
                   + f_3 * pc_x[k] * lsg_366[k];

        t_511[k] = f_11 * ksg_258[k]
                   + f_3 * pc_z[k] * lsg_363[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pc_x, pc_y, ksg_275, ksg_369, ksg_370, \
                         ksg_371, lsf0_249, lsf1_249, lsg_365, lsg_369, lsg_370, \
                         lsg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_11 * ksg_275[k]
                   + f_3 * pc_y[k] * lsg_365[k];

        t_513[k] = f_10 * ksg_369[k]
                   + f_4 * lsf0_249[k]
                   - f_5 * lsf1_249[k]
                   + f_3 * pc_x[k] * lsg_369[k];

        t_514[k] = f_10 * ksg_370[k]
                   + f_3 * pc_x[k] * lsg_370[k];

        t_515[k] = f_10 * ksg_371[k]
                   + f_3 * pc_x[k] * lsg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, pc_y, ksg_280, ksg_372, ksg_373, \
                         ksg_374, lsf0_246, lsf1_246, lsg_370, lsg_372, lsg_373, \
                         lsg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_10 * ksg_372[k]
                   + f_3 * pc_x[k] * lsg_372[k];

        t_517[k] = f_10 * ksg_373[k]
                   + f_3 * pc_x[k] * lsg_373[k];

        t_518[k] = f_10 * ksg_374[k]
                   + f_3 * pc_x[k] * lsg_374[k];

        t_519[k] = f_11 * ksg_280[k]
                   + f_1 * lsf0_246[k]
                   - f_2 * lsf1_246[k]
                   + f_3 * pc_y[k] * lsg_370[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, ksg_265, ksg_282, ksg_283, lsf0_248, \
                         lsf0_249, lsf1_248, lsf1_249, lsg_370, lsg_372, \
                         lsg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_11 * ksg_265[k]
                   + f_3 * pc_z[k] * lsg_370[k];

        t_521[k] = f_11 * ksg_282[k]
                   + f_6 * lsf0_248[k]
                   - f_7 * lsf1_248[k]
                   + f_3 * pc_y[k] * lsg_372[k];

        t_522[k] = f_11 * ksg_283[k]
                   + f_4 * lsf0_249[k]
                   - f_5 * lsf1_249[k]
                   + f_3 * pc_y[k] * lsg_373[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, ksg_269, ksg_284, ksg_375, \
                         lsf0_249, lsf0_250, lsf1_249, lsf1_250, lsg_374, \
                         lsg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_11 * ksg_284[k]
                   + f_3 * pc_y[k] * lsg_374[k];

        t_524[k] = f_11 * ksg_269[k]
                   + f_1 * lsf0_249[k]
                   - f_2 * lsf1_249[k]
                   + f_3 * pc_z[k] * lsg_374[k];

        t_525[k] = f_10 * ksg_375[k]
                   + f_1 * lsf0_250[k]
                   - f_2 * lsf1_250[k]
                   + f_3 * pc_x[k] * lsg_375[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_x, pc_y, pc_z, ksg_270, ksg_285, \
                         ksg_287, ksg_378, lsf0_253, lsf1_253, lsg_375, lsg_377, \
                         lsg_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * ksg_285[k]
                   + f_3 * pc_y[k] * lsg_375[k];

        t_527[k] = f_17 * ksg_270[k]
                   + f_3 * pc_z[k] * lsg_375[k];

        t_528[k] = f_10 * ksg_378[k]
                   + f_6 * lsf0_253[k]
                   - f_7 * lsf1_253[k]
                   + f_3 * pc_x[k] * lsg_378[k];

        t_529[k] = f_10 * ksg_287[k]
                   + f_3 * pc_y[k] * lsg_377[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, pc_z, ksg_273, ksg_380, ksg_381, lsf0_255, \
                         lsf0_256, lsf1_255, lsf1_256, lsg_378, lsg_380, \
                         lsg_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_10 * ksg_380[k]
                   + f_6 * lsf0_255[k]
                   - f_7 * lsf1_255[k]
                   + f_3 * pc_x[k] * lsg_380[k];

        t_531[k] = f_10 * ksg_381[k]
                   + f_4 * lsf0_256[k]
                   - f_5 * lsf1_256[k]
                   + f_3 * pc_x[k] * lsg_381[k];

        t_532[k] = f_17 * ksg_273[k]
                   + f_3 * pc_z[k] * lsg_378[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_y, ksg_290, ksg_384, ksg_385, \
                         ksg_386, lsf0_259, lsf1_259, lsg_380, lsg_384, lsg_385, \
                         lsg_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_10 * ksg_290[k]
                   + f_3 * pc_y[k] * lsg_380[k];

        t_534[k] = f_10 * ksg_384[k]
                   + f_4 * lsf0_259[k]
                   - f_5 * lsf1_259[k]
                   + f_3 * pc_x[k] * lsg_384[k];

        t_535[k] = f_10 * ksg_385[k]
                   + f_3 * pc_x[k] * lsg_385[k];

        t_536[k] = f_10 * ksg_386[k]
                   + f_3 * pc_x[k] * lsg_386[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, ksg_295, ksg_387, ksg_388, \
                         ksg_389, lsf0_256, lsf1_256, lsg_385, lsg_387, lsg_388, \
                         lsg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_10 * ksg_387[k]
                   + f_3 * pc_x[k] * lsg_387[k];

        t_538[k] = f_10 * ksg_388[k]
                   + f_3 * pc_x[k] * lsg_388[k];

        t_539[k] = f_10 * ksg_389[k]
                   + f_3 * pc_x[k] * lsg_389[k];

        t_540[k] = f_10 * ksg_295[k]
                   + f_1 * lsf0_256[k]
                   - f_2 * lsf1_256[k]
                   + f_3 * pc_y[k] * lsg_385[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pc_y, pc_z, ksg_280, ksg_297, ksg_298, lsf0_258, \
                         lsf0_259, lsf1_258, lsf1_259, lsg_385, lsg_387, \
                         lsg_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_17 * ksg_280[k]
                   + f_3 * pc_z[k] * lsg_385[k];

        t_542[k] = f_10 * ksg_297[k]
                   + f_6 * lsf0_258[k]
                   - f_7 * lsf1_258[k]
                   + f_3 * pc_y[k] * lsg_387[k];

        t_543[k] = f_10 * ksg_298[k]
                   + f_4 * lsf0_259[k]
                   - f_5 * lsf1_259[k]
                   + f_3 * pc_y[k] * lsg_388[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_y, pc_z, ksh0_420, ksg_284, \
                         ksg_299, ksg_300, ksh1_420, lsf0_259, lsf1_259, lsg_389, \
                         lsg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = f_10 * ksg_299[k]
                   + f_3 * pc_y[k] * lsg_389[k];

        t_545[k] = f_17 * ksg_284[k]
                   + f_1 * lsf0_259[k]
                   - f_2 * lsf1_259[k]
                   + f_3 * pc_z[k] * lsg_389[k];

        t_546[k] = pa_y[k] * ksh0_420[k]
                   - f_8 * pc_y[k] * ksh1_420[k];

        t_547[k] = f_9 * ksg_300[k]
                   + f_3 * pc_y[k] * lsg_390[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, pa_y, pc_y, pc_z, ksh0_423, ksh0_425, \
                         ksg_285, ksg_301, ksg_302, ksh1_423, ksh1_425, lsg_390, \
                         lsg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * ksg_285[k]
                   + f_3 * pc_z[k] * lsg_390[k];

        t_549[k] = pa_y[k] * ksh0_423[k]
                   + f_10 * ksg_301[k]
                   - f_8 * pc_y[k] * ksh1_423[k];

        t_550[k] = f_9 * ksg_302[k]
                   + f_3 * pc_y[k] * lsg_392[k];

        t_551[k] = pa_y[k] * ksh0_425[k]
                   - f_8 * pc_y[k] * ksh1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, t_555, pa_y, pc_y, pc_z, ksh0_426, ksh0_429, \
                         ksg_288, ksg_303, ksg_305, ksh1_426, ksh1_429, lsg_393, \
                         lsg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_y[k] * ksh0_426[k]
                   + f_11 * ksg_303[k]
                   - f_8 * pc_y[k] * ksh1_426[k];

        t_553[k] = f_16 * ksg_288[k]
                   + f_3 * pc_z[k] * lsg_393[k];

        t_554[k] = f_9 * ksg_305[k]
                   + f_3 * pc_y[k] * lsg_395[k];

        t_555[k] = pa_y[k] * ksh0_429[k]
                   - f_8 * pc_y[k] * ksh1_429[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, t_559, t_560, pc_x, ksg_400, ksg_401, ksg_402, \
                         ksg_403, ksg_404, lsg_400, lsg_401, lsg_402, lsg_403, \
                         lsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_10 * ksg_400[k]
                   + f_3 * pc_x[k] * lsg_400[k];

        t_557[k] = f_10 * ksg_401[k]
                   + f_3 * pc_x[k] * lsg_401[k];

        t_558[k] = f_10 * ksg_402[k]
                   + f_3 * pc_x[k] * lsg_402[k];

        t_559[k] = f_10 * ksg_403[k]
                   + f_3 * pc_x[k] * lsg_403[k];

        t_560[k] = f_10 * ksg_404[k]
                   + f_3 * pc_x[k] * lsg_404[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, pc_y, pc_z, ksg_295, ksg_310, ksg_312, lsf0_266, \
                         lsf0_268, lsf1_266, lsf1_268, lsg_400, \
                         lsg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_9 * ksg_310[k]
                   + f_1 * lsf0_266[k]
                   - f_2 * lsf1_266[k]
                   + f_3 * pc_y[k] * lsg_400[k];

        t_562[k] = f_16 * ksg_295[k]
                   + f_3 * pc_z[k] * lsg_400[k];

        t_563[k] = f_9 * ksg_312[k]
                   + f_6 * lsf0_268[k]
                   - f_7 * lsf1_268[k]
                   + f_3 * pc_y[k] * lsg_402[k];
    }

#pragma omp simd aligned(t_564, t_565, t_566, pa_y, pc_y, ksh0_440, ksg_313, ksg_314, \
                         ksh1_440, lsf0_269, lsf1_269, lsg_403, \
                         lsg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_564[k] = f_9 * ksg_313[k]
                   + f_4 * lsf0_269[k]
                   - f_5 * lsf1_269[k]
                   + f_3 * pc_y[k] * lsg_403[k];

        t_565[k] = f_9 * ksg_314[k]
                   + f_3 * pc_y[k] * lsg_404[k];

        t_566[k] = pa_y[k] * ksh0_440[k]
                   - f_8 * pc_y[k] * ksh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, ksg_300, \
                         ksg_405, lsf0_270, lsf1_270, lsg_405, lsg_406, \
                         lsg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_10 * ksg_405[k]
                   + f_1 * lsf0_270[k]
                   - f_2 * lsf1_270[k]
                   + f_3 * pc_x[k] * lsg_405[k];

        t_568[k] = f_3 * pc_y[k] * lsg_405[k];

        t_569[k] = f_15 * ksg_300[k]
                   + f_3 * pc_z[k] * lsg_405[k];

        t_570[k] = f_4 * lsf0_270[k]
                   - f_5 * lsf1_270[k]
                   + f_3 * pc_y[k] * lsg_406[k];

        t_571[k] = f_3 * pc_y[k] * lsg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, ksg_410, lsf0_271, lsf0_272, \
                         lsf0_275, lsf1_271, lsf1_272, lsf1_275, lsg_408, lsg_409, \
                         lsg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_10 * ksg_410[k]
                   + f_6 * lsf0_275[k]
                   - f_7 * lsf1_275[k]
                   + f_3 * pc_x[k] * lsg_410[k];

        t_573[k] = f_6 * lsf0_271[k]
                   - f_7 * lsf1_271[k]
                   + f_3 * pc_y[k] * lsg_408[k];

        t_574[k] = f_4 * lsf0_272[k]
                   - f_5 * lsf1_272[k]
                   + f_3 * pc_y[k] * lsg_409[k];

        t_575[k] = f_3 * pc_y[k] * lsg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, pc_x, ksg_414, ksg_415, ksg_416, ksg_417, \
                         lsf0_279, lsf1_279, lsg_414, lsg_415, lsg_416, \
                         lsg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_10 * ksg_414[k]
                   + f_4 * lsf0_279[k]
                   - f_5 * lsf1_279[k]
                   + f_3 * pc_x[k] * lsg_414[k];

        t_577[k] = f_10 * ksg_415[k]
                   + f_3 * pc_x[k] * lsg_415[k];

        t_578[k] = f_10 * ksg_416[k]
                   + f_3 * pc_x[k] * lsg_416[k];

        t_579[k] = f_10 * ksg_417[k]
                   + f_3 * pc_x[k] * lsg_417[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, t_583, pc_x, pc_y, ksg_419, lsf0_276, lsf0_277, \
                         lsf1_276, lsf1_277, lsg_414, lsg_415, lsg_416, \
                         lsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_3 * pc_y[k] * lsg_414[k];

        t_581[k] = f_10 * ksg_419[k]
                   + f_3 * pc_x[k] * lsg_419[k];

        t_582[k] = f_1 * lsf0_276[k]
                   - f_2 * lsf1_276[k]
                   + f_3 * pc_y[k] * lsg_415[k];

        t_583[k] = f_13 * lsf0_277[k]
                   - f_14 * lsf1_277[k]
                   + f_3 * pc_y[k] * lsg_416[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, ksg_314, lsf0_278, lsf0_279, \
                         lsf1_278, lsf1_279, lsg_417, lsg_418, \
                         lsg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * lsf0_278[k]
                   - f_7 * lsf1_278[k]
                   + f_3 * pc_y[k] * lsg_417[k];

        t_585[k] = f_4 * lsf0_279[k]
                   - f_5 * lsf1_279[k]
                   + f_3 * pc_y[k] * lsg_418[k];

        t_586[k] = f_3 * pc_y[k] * lsg_419[k];

        t_587[k] = f_15 * ksg_314[k]
                   + f_1 * lsf0_279[k]
                   - f_2 * lsf1_279[k]
                   + f_3 * pc_z[k] * lsg_419[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_x, pc_x, pc_y, pc_z, ksh0_588, \
                         ksh0_591, ksg_315, ksg_420, ksg_423, ksh1_588, ksh1_591, \
                         lsg_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_x[k] * ksh0_588[k]
                   + f_16 * ksg_420[k]
                   - f_8 * pc_x[k] * ksh1_588[k];

        t_589[k] = f_12 * ksg_315[k]
                   + f_3 * pc_y[k] * lsg_420[k];

        t_590[k] = f_3 * pc_z[k] * lsg_420[k];

        t_591[k] = pa_x[k] * ksh0_591[k]
                   + f_11 * ksg_423[k]
                   - f_8 * pc_x[k] * ksh1_591[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_x, pc_x, pc_z, ksh0_594, ksg_426, \
                         ksh1_594, lsf0_280, lsf1_280, lsg_421, lsg_422, \
                         lsg_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * lsg_421[k];

        t_593[k] = f_4 * lsf0_280[k]
                   - f_5 * lsf1_280[k]
                   + f_3 * pc_z[k] * lsg_422[k];

        t_594[k] = pa_x[k] * ksh0_594[k]
                   + f_10 * ksg_426[k]
                   - f_8 * pc_x[k] * ksh1_594[k];

        t_595[k] = f_3 * pc_z[k] * lsg_423[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, ksg_320, ksg_430, \
                         lsf0_282, lsf1_282, lsg_425, lsg_426, \
                         lsg_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_12 * ksg_320[k]
                   + f_3 * pc_y[k] * lsg_425[k];

        t_597[k] = f_6 * lsf0_282[k]
                   - f_7 * lsf1_282[k]
                   + f_3 * pc_z[k] * lsg_425[k];

        t_598[k] = f_9 * ksg_430[k]
                   + f_3 * pc_x[k] * lsg_430[k];

        t_599[k] = f_3 * pc_z[k] * lsg_426[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pa_x, pc_x, ksh0_603, ksg_432, ksg_433, \
                         ksg_434, ksh1_603, lsg_432, lsg_433, lsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_9 * ksg_432[k]
                   + f_3 * pc_x[k] * lsg_432[k];

        t_601[k] = f_9 * ksg_433[k]
                   + f_3 * pc_x[k] * lsg_433[k];

        t_602[k] = f_9 * ksg_434[k]
                   + f_3 * pc_x[k] * lsg_434[k];

        t_603[k] = pa_x[k] * ksh0_603[k]
                   - f_8 * pc_x[k] * ksh1_603[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_441 = buffer.data(ksh0 + 441);
    const auto *ksh0_444 = buffer.data(ksh0 + 444);
    const auto *ksh0_447 = buffer.data(ksh0 + 447);
    const auto *ksh0_567 = buffer.data(ksh0 + 567);
    const auto *ksh0_572 = buffer.data(ksh0 + 572);
    const auto *ksh0_576 = buffer.data(ksh0 + 576);
    const auto *ksh0_605 = buffer.data(ksh0 + 605);
    const auto *ksh0_606 = buffer.data(ksh0 + 606);
    const auto *ksh0_608 = buffer.data(ksh0 + 608);
    const auto *ksh0_614 = buffer.data(ksh0 + 614);
    const auto *ksh0_618 = buffer.data(ksh0 + 618);
    const auto *ksh0_624 = buffer.data(ksh0 + 624);
    const auto *ksh0_626 = buffer.data(ksh0 + 626);
    const auto *ksh0_627 = buffer.data(ksh0 + 627);
    const auto *ksh0_629 = buffer.data(ksh0 + 629);
    const auto *ksh0_630 = buffer.data(ksh0 + 630);
    const auto *ksh0_633 = buffer.data(ksh0 + 633);
    const auto *ksh0_635 = buffer.data(ksh0 + 635);
    const auto *ksh0_636 = buffer.data(ksh0 + 636);
    const auto *ksh0_639 = buffer.data(ksh0 + 639);
    const auto *ksh0_645 = buffer.data(ksh0 + 645);
    const auto *ksh0_647 = buffer.data(ksh0 + 647);
    const auto *ksh0_648 = buffer.data(ksh0 + 648);
    const auto *ksh0_650 = buffer.data(ksh0 + 650);
    const auto *ksh0_651 = buffer.data(ksh0 + 651);
    const auto *ksh0_654 = buffer.data(ksh0 + 654);
    const auto *ksh0_656 = buffer.data(ksh0 + 656);
    const auto *ksh0_657 = buffer.data(ksh0 + 657);
    const auto *ksh0_660 = buffer.data(ksh0 + 660);
    const auto *ksh0_666 = buffer.data(ksh0 + 666);
    const auto *ksh0_668 = buffer.data(ksh0 + 668);
    const auto *ksh0_669 = buffer.data(ksh0 + 669);
    const auto *ksh0_671 = buffer.data(ksh0 + 671);
    const auto *ksh0_672 = buffer.data(ksh0 + 672);
    const auto *ksh0_675 = buffer.data(ksh0 + 675);
    const auto *ksh0_677 = buffer.data(ksh0 + 677);
    const auto *ksh0_678 = buffer.data(ksh0 + 678);
    const auto *ksh0_681 = buffer.data(ksh0 + 681);
    const auto *ksh0_687 = buffer.data(ksh0 + 687);
    const auto *ksh0_689 = buffer.data(ksh0 + 689);
    const auto *ksh0_690 = buffer.data(ksh0 + 690);
    const auto *ksh0_692 = buffer.data(ksh0 + 692);
    const auto *ksh0_693 = buffer.data(ksh0 + 693);
    const auto *ksh0_696 = buffer.data(ksh0 + 696);
    const auto *ksh0_698 = buffer.data(ksh0 + 698);
    const auto *ksh0_699 = buffer.data(ksh0 + 699);
    const auto *ksh0_702 = buffer.data(ksh0 + 702);
    const auto *ksh0_708 = buffer.data(ksh0 + 708);
    const auto *ksh0_710 = buffer.data(ksh0 + 710);
    const auto *ksh0_711 = buffer.data(ksh0 + 711);
    const auto *ksh0_713 = buffer.data(ksh0 + 713);
    const auto *ksh0_717 = buffer.data(ksh0 + 717);
    const auto *ksh0_720 = buffer.data(ksh0 + 720);

    const auto *ksg_315 = buffer.data(ksg + 315);
    const auto *ksg_318 = buffer.data(ksg + 318);
    const auto *ksg_325 = buffer.data(ksg + 325);
    const auto *ksg_329 = buffer.data(ksg + 329);
    const auto *ksg_330 = buffer.data(ksg + 330);
    const auto *ksg_332 = buffer.data(ksg + 332);
    const auto *ksg_333 = buffer.data(ksg + 333);
    const auto *ksg_335 = buffer.data(ksg + 335);
    const auto *ksg_340 = buffer.data(ksg + 340);
    const auto *ksg_344 = buffer.data(ksg + 344);
    const auto *ksg_345 = buffer.data(ksg + 345);
    const auto *ksg_347 = buffer.data(ksg + 347);
    const auto *ksg_348 = buffer.data(ksg + 348);
    const auto *ksg_350 = buffer.data(ksg + 350);
    const auto *ksg_355 = buffer.data(ksg + 355);
    const auto *ksg_359 = buffer.data(ksg + 359);
    const auto *ksg_360 = buffer.data(ksg + 360);
    const auto *ksg_362 = buffer.data(ksg + 362);
    const auto *ksg_363 = buffer.data(ksg + 363);
    const auto *ksg_365 = buffer.data(ksg + 365);
    const auto *ksg_370 = buffer.data(ksg + 370);
    const auto *ksg_374 = buffer.data(ksg + 374);
    const auto *ksg_375 = buffer.data(ksg + 375);
    const auto *ksg_377 = buffer.data(ksg + 377);
    const auto *ksg_378 = buffer.data(ksg + 378);
    const auto *ksg_380 = buffer.data(ksg + 380);
    const auto *ksg_385 = buffer.data(ksg + 385);
    const auto *ksg_389 = buffer.data(ksg + 389);
    const auto *ksg_390 = buffer.data(ksg + 390);
    const auto *ksg_392 = buffer.data(ksg + 392);
    const auto *ksg_393 = buffer.data(ksg + 393);
    const auto *ksg_395 = buffer.data(ksg + 395);
    const auto *ksg_404 = buffer.data(ksg + 404);
    const auto *ksg_405 = buffer.data(ksg + 405);
    const auto *ksg_407 = buffer.data(ksg + 407);
    const auto *ksg_410 = buffer.data(ksg + 410);
    const auto *ksg_440 = buffer.data(ksg + 440);
    const auto *ksg_444 = buffer.data(ksg + 444);
    const auto *ksg_445 = buffer.data(ksg + 445);
    const auto *ksg_446 = buffer.data(ksg + 446);
    const auto *ksg_447 = buffer.data(ksg + 447);
    const auto *ksg_448 = buffer.data(ksg + 448);
    const auto *ksg_449 = buffer.data(ksg + 449);
    const auto *ksg_450 = buffer.data(ksg + 450);
    const auto *ksg_453 = buffer.data(ksg + 453);
    const auto *ksg_455 = buffer.data(ksg + 455);
    const auto *ksg_456 = buffer.data(ksg + 456);
    const auto *ksg_459 = buffer.data(ksg + 459);
    const auto *ksg_460 = buffer.data(ksg + 460);
    const auto *ksg_461 = buffer.data(ksg + 461);
    const auto *ksg_462 = buffer.data(ksg + 462);
    const auto *ksg_463 = buffer.data(ksg + 463);
    const auto *ksg_464 = buffer.data(ksg + 464);
    const auto *ksg_465 = buffer.data(ksg + 465);
    const auto *ksg_468 = buffer.data(ksg + 468);
    const auto *ksg_470 = buffer.data(ksg + 470);
    const auto *ksg_471 = buffer.data(ksg + 471);
    const auto *ksg_474 = buffer.data(ksg + 474);
    const auto *ksg_475 = buffer.data(ksg + 475);
    const auto *ksg_476 = buffer.data(ksg + 476);
    const auto *ksg_477 = buffer.data(ksg + 477);
    const auto *ksg_478 = buffer.data(ksg + 478);
    const auto *ksg_479 = buffer.data(ksg + 479);
    const auto *ksg_480 = buffer.data(ksg + 480);
    const auto *ksg_483 = buffer.data(ksg + 483);
    const auto *ksg_485 = buffer.data(ksg + 485);
    const auto *ksg_486 = buffer.data(ksg + 486);
    const auto *ksg_489 = buffer.data(ksg + 489);
    const auto *ksg_490 = buffer.data(ksg + 490);
    const auto *ksg_491 = buffer.data(ksg + 491);
    const auto *ksg_492 = buffer.data(ksg + 492);
    const auto *ksg_493 = buffer.data(ksg + 493);
    const auto *ksg_494 = buffer.data(ksg + 494);
    const auto *ksg_495 = buffer.data(ksg + 495);
    const auto *ksg_498 = buffer.data(ksg + 498);
    const auto *ksg_500 = buffer.data(ksg + 500);
    const auto *ksg_501 = buffer.data(ksg + 501);
    const auto *ksg_504 = buffer.data(ksg + 504);
    const auto *ksg_505 = buffer.data(ksg + 505);
    const auto *ksg_506 = buffer.data(ksg + 506);
    const auto *ksg_507 = buffer.data(ksg + 507);
    const auto *ksg_508 = buffer.data(ksg + 508);
    const auto *ksg_509 = buffer.data(ksg + 509);
    const auto *ksg_513 = buffer.data(ksg + 513);
    const auto *ksg_516 = buffer.data(ksg + 516);
    const auto *ksg_520 = buffer.data(ksg + 520);
    const auto *ksg_521 = buffer.data(ksg + 521);
    const auto *ksg_522 = buffer.data(ksg + 522);

    const auto *ksh1_441 = buffer.data(ksh1 + 441);
    const auto *ksh1_444 = buffer.data(ksh1 + 444);
    const auto *ksh1_447 = buffer.data(ksh1 + 447);
    const auto *ksh1_567 = buffer.data(ksh1 + 567);
    const auto *ksh1_572 = buffer.data(ksh1 + 572);
    const auto *ksh1_576 = buffer.data(ksh1 + 576);
    const auto *ksh1_605 = buffer.data(ksh1 + 605);
    const auto *ksh1_606 = buffer.data(ksh1 + 606);
    const auto *ksh1_608 = buffer.data(ksh1 + 608);
    const auto *ksh1_614 = buffer.data(ksh1 + 614);
    const auto *ksh1_618 = buffer.data(ksh1 + 618);
    const auto *ksh1_624 = buffer.data(ksh1 + 624);
    const auto *ksh1_626 = buffer.data(ksh1 + 626);
    const auto *ksh1_627 = buffer.data(ksh1 + 627);
    const auto *ksh1_629 = buffer.data(ksh1 + 629);
    const auto *ksh1_630 = buffer.data(ksh1 + 630);
    const auto *ksh1_633 = buffer.data(ksh1 + 633);
    const auto *ksh1_635 = buffer.data(ksh1 + 635);
    const auto *ksh1_636 = buffer.data(ksh1 + 636);
    const auto *ksh1_639 = buffer.data(ksh1 + 639);
    const auto *ksh1_645 = buffer.data(ksh1 + 645);
    const auto *ksh1_647 = buffer.data(ksh1 + 647);
    const auto *ksh1_648 = buffer.data(ksh1 + 648);
    const auto *ksh1_650 = buffer.data(ksh1 + 650);
    const auto *ksh1_651 = buffer.data(ksh1 + 651);
    const auto *ksh1_654 = buffer.data(ksh1 + 654);
    const auto *ksh1_656 = buffer.data(ksh1 + 656);
    const auto *ksh1_657 = buffer.data(ksh1 + 657);
    const auto *ksh1_660 = buffer.data(ksh1 + 660);
    const auto *ksh1_666 = buffer.data(ksh1 + 666);
    const auto *ksh1_668 = buffer.data(ksh1 + 668);
    const auto *ksh1_669 = buffer.data(ksh1 + 669);
    const auto *ksh1_671 = buffer.data(ksh1 + 671);
    const auto *ksh1_672 = buffer.data(ksh1 + 672);
    const auto *ksh1_675 = buffer.data(ksh1 + 675);
    const auto *ksh1_677 = buffer.data(ksh1 + 677);
    const auto *ksh1_678 = buffer.data(ksh1 + 678);
    const auto *ksh1_681 = buffer.data(ksh1 + 681);
    const auto *ksh1_687 = buffer.data(ksh1 + 687);
    const auto *ksh1_689 = buffer.data(ksh1 + 689);
    const auto *ksh1_690 = buffer.data(ksh1 + 690);
    const auto *ksh1_692 = buffer.data(ksh1 + 692);
    const auto *ksh1_693 = buffer.data(ksh1 + 693);
    const auto *ksh1_696 = buffer.data(ksh1 + 696);
    const auto *ksh1_698 = buffer.data(ksh1 + 698);
    const auto *ksh1_699 = buffer.data(ksh1 + 699);
    const auto *ksh1_702 = buffer.data(ksh1 + 702);
    const auto *ksh1_708 = buffer.data(ksh1 + 708);
    const auto *ksh1_710 = buffer.data(ksh1 + 710);
    const auto *ksh1_711 = buffer.data(ksh1 + 711);
    const auto *ksh1_713 = buffer.data(ksh1 + 713);
    const auto *ksh1_717 = buffer.data(ksh1 + 717);
    const auto *ksh1_720 = buffer.data(ksh1 + 720);

    const auto *lsg_430 = buffer.data(lsg + 430);
    const auto *lsg_434 = buffer.data(lsg + 434);
    const auto *lsg_435 = buffer.data(lsg + 435);
    const auto *lsg_437 = buffer.data(lsg + 437);
    const auto *lsg_438 = buffer.data(lsg + 438);
    const auto *lsg_440 = buffer.data(lsg + 440);
    const auto *lsg_445 = buffer.data(lsg + 445);
    const auto *lsg_446 = buffer.data(lsg + 446);
    const auto *lsg_447 = buffer.data(lsg + 447);
    const auto *lsg_448 = buffer.data(lsg + 448);
    const auto *lsg_449 = buffer.data(lsg + 449);
    const auto *lsg_450 = buffer.data(lsg + 450);
    const auto *lsg_452 = buffer.data(lsg + 452);
    const auto *lsg_453 = buffer.data(lsg + 453);
    const auto *lsg_455 = buffer.data(lsg + 455);
    const auto *lsg_460 = buffer.data(lsg + 460);
    const auto *lsg_461 = buffer.data(lsg + 461);
    const auto *lsg_462 = buffer.data(lsg + 462);
    const auto *lsg_463 = buffer.data(lsg + 463);
    const auto *lsg_464 = buffer.data(lsg + 464);
    const auto *lsg_465 = buffer.data(lsg + 465);
    const auto *lsg_467 = buffer.data(lsg + 467);
    const auto *lsg_468 = buffer.data(lsg + 468);
    const auto *lsg_470 = buffer.data(lsg + 470);
    const auto *lsg_475 = buffer.data(lsg + 475);
    const auto *lsg_476 = buffer.data(lsg + 476);
    const auto *lsg_477 = buffer.data(lsg + 477);
    const auto *lsg_478 = buffer.data(lsg + 478);
    const auto *lsg_479 = buffer.data(lsg + 479);
    const auto *lsg_480 = buffer.data(lsg + 480);
    const auto *lsg_482 = buffer.data(lsg + 482);
    const auto *lsg_483 = buffer.data(lsg + 483);
    const auto *lsg_485 = buffer.data(lsg + 485);
    const auto *lsg_490 = buffer.data(lsg + 490);
    const auto *lsg_491 = buffer.data(lsg + 491);
    const auto *lsg_492 = buffer.data(lsg + 492);
    const auto *lsg_493 = buffer.data(lsg + 493);
    const auto *lsg_494 = buffer.data(lsg + 494);
    const auto *lsg_495 = buffer.data(lsg + 495);
    const auto *lsg_497 = buffer.data(lsg + 497);
    const auto *lsg_498 = buffer.data(lsg + 498);
    const auto *lsg_500 = buffer.data(lsg + 500);
    const auto *lsg_505 = buffer.data(lsg + 505);
    const auto *lsg_506 = buffer.data(lsg + 506);
    const auto *lsg_507 = buffer.data(lsg + 507);
    const auto *lsg_508 = buffer.data(lsg + 508);
    const auto *lsg_509 = buffer.data(lsg + 509);
    const auto *lsg_510 = buffer.data(lsg + 510);
    const auto *lsg_512 = buffer.data(lsg + 512);
    const auto *lsg_513 = buffer.data(lsg + 513);
    const auto *lsg_515 = buffer.data(lsg + 515);
    const auto *lsg_520 = buffer.data(lsg + 520);
    const auto *lsg_521 = buffer.data(lsg + 521);
    const auto *lsg_522 = buffer.data(lsg + 522);

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_x, pc_x, pc_y, pc_z, ksh0_605, \
                         ksh0_606, ksg_329, ksh1_605, ksh1_606, lsg_430, \
                         lsg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * lsg_430[k];

        t_605[k] = pa_x[k] * ksh0_605[k]
                   - f_8 * pc_x[k] * ksh1_605[k];

        t_606[k] = pa_x[k] * ksh0_606[k]
                   - f_8 * pc_x[k] * ksh1_606[k];

        t_607[k] = f_12 * ksg_329[k]
                   + f_3 * pc_y[k] * lsg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_x, pa_z, pc_x, pc_y, pc_z, ksh0_441, \
                         ksh0_608, ksg_315, ksg_330, ksh1_441, ksh1_608, \
                         lsg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pa_x[k] * ksh0_608[k]
                   - f_8 * pc_x[k] * ksh1_608[k];

        t_609[k] = pa_z[k] * ksh0_441[k]
                   - f_8 * pc_z[k] * ksh1_441[k];

        t_610[k] = f_15 * ksg_330[k]
                   + f_3 * pc_y[k] * lsg_435[k];

        t_611[k] = f_9 * ksg_315[k]
                   + f_3 * pc_z[k] * lsg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_x, pa_z, pc_x, pc_y, pc_z, ksh0_444, \
                         ksh0_614, ksg_332, ksg_440, ksh1_444, ksh1_614, \
                         lsg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pa_z[k] * ksh0_444[k]
                   - f_8 * pc_z[k] * ksh1_444[k];

        t_613[k] = f_15 * ksg_332[k]
                   + f_3 * pc_y[k] * lsg_437[k];

        t_614[k] = pa_x[k] * ksh0_614[k]
                   + f_11 * ksg_440[k]
                   - f_8 * pc_x[k] * ksh1_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pa_z, pc_y, pc_z, ksh0_447, ksg_318, ksg_335, \
                         ksh1_447, lsg_438, lsg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = pa_z[k] * ksh0_447[k]
                   - f_8 * pc_z[k] * ksh1_447[k];

        t_616[k] = f_9 * ksg_318[k]
                   + f_3 * pc_z[k] * lsg_438[k];

        t_617[k] = f_15 * ksg_335[k]
                   + f_3 * pc_y[k] * lsg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pa_x, pc_x, ksh0_618, ksg_444, ksg_445, \
                         ksg_446, ksg_447, ksh1_618, lsg_445, lsg_446, \
                         lsg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = pa_x[k] * ksh0_618[k]
                   + f_10 * ksg_444[k]
                   - f_8 * pc_x[k] * ksh1_618[k];

        t_619[k] = f_9 * ksg_445[k]
                   + f_3 * pc_x[k] * lsg_445[k];

        t_620[k] = f_9 * ksg_446[k]
                   + f_3 * pc_x[k] * lsg_446[k];

        t_621[k] = f_9 * ksg_447[k]
                   + f_3 * pc_x[k] * lsg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_x, pc_x, pc_z, ksh0_624, ksg_325, \
                         ksg_448, ksg_449, ksh1_624, lsg_445, lsg_448, \
                         lsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_9 * ksg_448[k]
                   + f_3 * pc_x[k] * lsg_448[k];

        t_623[k] = f_9 * ksg_449[k]
                   + f_3 * pc_x[k] * lsg_449[k];

        t_624[k] = pa_x[k] * ksh0_624[k]
                   - f_8 * pc_x[k] * ksh1_624[k];

        t_625[k] = f_9 * ksg_325[k]
                   + f_3 * pc_z[k] * lsg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pa_x, pc_x, pc_y, ksh0_626, ksh0_627, \
                         ksh0_629, ksg_344, ksh1_626, ksh1_627, ksh1_629, \
                         lsg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_x[k] * ksh0_626[k]
                   - f_8 * pc_x[k] * ksh1_626[k];

        t_627[k] = pa_x[k] * ksh0_627[k]
                   - f_8 * pc_x[k] * ksh1_627[k];

        t_628[k] = f_15 * ksg_344[k]
                   + f_3 * pc_y[k] * lsg_449[k];

        t_629[k] = pa_x[k] * ksh0_629[k]
                   - f_8 * pc_x[k] * ksh1_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pa_x, pc_x, pc_y, pc_z, ksh0_630, ksg_330, \
                         ksg_345, ksg_450, ksh1_630, lsg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pa_x[k] * ksh0_630[k]
                   + f_16 * ksg_450[k]
                   - f_8 * pc_x[k] * ksh1_630[k];

        t_631[k] = f_16 * ksg_345[k]
                   + f_3 * pc_y[k] * lsg_450[k];

        t_632[k] = f_10 * ksg_330[k]
                   + f_3 * pc_z[k] * lsg_450[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pa_x, pc_x, pc_y, ksh0_633, ksh0_635, ksg_347, \
                         ksg_453, ksg_455, ksh1_633, ksh1_635, \
                         lsg_452 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = pa_x[k] * ksh0_633[k]
                   + f_11 * ksg_453[k]
                   - f_8 * pc_x[k] * ksh1_633[k];

        t_634[k] = f_16 * ksg_347[k]
                   + f_3 * pc_y[k] * lsg_452[k];

        t_635[k] = pa_x[k] * ksh0_635[k]
                   + f_11 * ksg_455[k]
                   - f_8 * pc_x[k] * ksh1_635[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pa_x, pc_x, pc_y, pc_z, ksh0_636, ksg_333, \
                         ksg_350, ksg_456, ksh1_636, lsg_453, lsg_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = pa_x[k] * ksh0_636[k]
                   + f_10 * ksg_456[k]
                   - f_8 * pc_x[k] * ksh1_636[k];

        t_637[k] = f_10 * ksg_333[k]
                   + f_3 * pc_z[k] * lsg_453[k];

        t_638[k] = f_16 * ksg_350[k]
                   + f_3 * pc_y[k] * lsg_455[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pa_x, pc_x, ksh0_639, ksg_459, ksg_460, \
                         ksg_461, ksg_462, ksh1_639, lsg_460, lsg_461, \
                         lsg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pa_x[k] * ksh0_639[k]
                   + f_10 * ksg_459[k]
                   - f_8 * pc_x[k] * ksh1_639[k];

        t_640[k] = f_9 * ksg_460[k]
                   + f_3 * pc_x[k] * lsg_460[k];

        t_641[k] = f_9 * ksg_461[k]
                   + f_3 * pc_x[k] * lsg_461[k];

        t_642[k] = f_9 * ksg_462[k]
                   + f_3 * pc_x[k] * lsg_462[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pa_x, pc_x, pc_z, ksh0_645, ksg_340, \
                         ksg_463, ksg_464, ksh1_645, lsg_460, lsg_463, \
                         lsg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_9 * ksg_463[k]
                   + f_3 * pc_x[k] * lsg_463[k];

        t_644[k] = f_9 * ksg_464[k]
                   + f_3 * pc_x[k] * lsg_464[k];

        t_645[k] = pa_x[k] * ksh0_645[k]
                   - f_8 * pc_x[k] * ksh1_645[k];

        t_646[k] = f_10 * ksg_340[k]
                   + f_3 * pc_z[k] * lsg_460[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pa_x, pc_x, pc_y, ksh0_647, ksh0_648, \
                         ksh0_650, ksg_359, ksh1_647, ksh1_648, ksh1_650, \
                         lsg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = pa_x[k] * ksh0_647[k]
                   - f_8 * pc_x[k] * ksh1_647[k];

        t_648[k] = pa_x[k] * ksh0_648[k]
                   - f_8 * pc_x[k] * ksh1_648[k];

        t_649[k] = f_16 * ksg_359[k]
                   + f_3 * pc_y[k] * lsg_464[k];

        t_650[k] = pa_x[k] * ksh0_650[k]
                   - f_8 * pc_x[k] * ksh1_650[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pa_x, pc_x, pc_y, pc_z, ksh0_651, ksg_345, \
                         ksg_360, ksg_465, ksh1_651, lsg_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = pa_x[k] * ksh0_651[k]
                   + f_16 * ksg_465[k]
                   - f_8 * pc_x[k] * ksh1_651[k];

        t_652[k] = f_17 * ksg_360[k]
                   + f_3 * pc_y[k] * lsg_465[k];

        t_653[k] = f_11 * ksg_345[k]
                   + f_3 * pc_z[k] * lsg_465[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pa_x, pc_x, pc_y, ksh0_654, ksh0_656, ksg_362, \
                         ksg_468, ksg_470, ksh1_654, ksh1_656, \
                         lsg_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = pa_x[k] * ksh0_654[k]
                   + f_11 * ksg_468[k]
                   - f_8 * pc_x[k] * ksh1_654[k];

        t_655[k] = f_17 * ksg_362[k]
                   + f_3 * pc_y[k] * lsg_467[k];

        t_656[k] = pa_x[k] * ksh0_656[k]
                   + f_11 * ksg_470[k]
                   - f_8 * pc_x[k] * ksh1_656[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_x, pc_x, pc_y, pc_z, ksh0_657, ksg_348, \
                         ksg_365, ksg_471, ksh1_657, lsg_468, lsg_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pa_x[k] * ksh0_657[k]
                   + f_10 * ksg_471[k]
                   - f_8 * pc_x[k] * ksh1_657[k];

        t_658[k] = f_11 * ksg_348[k]
                   + f_3 * pc_z[k] * lsg_468[k];

        t_659[k] = f_17 * ksg_365[k]
                   + f_3 * pc_y[k] * lsg_470[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, pa_x, pc_x, ksh0_660, ksg_474, ksg_475, \
                         ksg_476, ksg_477, ksh1_660, lsg_475, lsg_476, \
                         lsg_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = pa_x[k] * ksh0_660[k]
                   + f_10 * ksg_474[k]
                   - f_8 * pc_x[k] * ksh1_660[k];

        t_661[k] = f_9 * ksg_475[k]
                   + f_3 * pc_x[k] * lsg_475[k];

        t_662[k] = f_9 * ksg_476[k]
                   + f_3 * pc_x[k] * lsg_476[k];

        t_663[k] = f_9 * ksg_477[k]
                   + f_3 * pc_x[k] * lsg_477[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, pa_x, pc_x, pc_z, ksh0_666, ksg_355, \
                         ksg_478, ksg_479, ksh1_666, lsg_475, lsg_478, \
                         lsg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_9 * ksg_478[k]
                   + f_3 * pc_x[k] * lsg_478[k];

        t_665[k] = f_9 * ksg_479[k]
                   + f_3 * pc_x[k] * lsg_479[k];

        t_666[k] = pa_x[k] * ksh0_666[k]
                   - f_8 * pc_x[k] * ksh1_666[k];

        t_667[k] = f_11 * ksg_355[k]
                   + f_3 * pc_z[k] * lsg_475[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, t_671, pa_x, pc_x, pc_y, ksh0_668, ksh0_669, \
                         ksh0_671, ksg_374, ksh1_668, ksh1_669, ksh1_671, \
                         lsg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = pa_x[k] * ksh0_668[k]
                   - f_8 * pc_x[k] * ksh1_668[k];

        t_669[k] = pa_x[k] * ksh0_669[k]
                   - f_8 * pc_x[k] * ksh1_669[k];

        t_670[k] = f_17 * ksg_374[k]
                   + f_3 * pc_y[k] * lsg_479[k];

        t_671[k] = pa_x[k] * ksh0_671[k]
                   - f_8 * pc_x[k] * ksh1_671[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, pa_x, pc_x, pc_y, pc_z, ksh0_672, ksg_360, \
                         ksg_375, ksg_480, ksh1_672, lsg_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = pa_x[k] * ksh0_672[k]
                   + f_16 * ksg_480[k]
                   - f_8 * pc_x[k] * ksh1_672[k];

        t_673[k] = f_11 * ksg_375[k]
                   + f_3 * pc_y[k] * lsg_480[k];

        t_674[k] = f_17 * ksg_360[k]
                   + f_3 * pc_z[k] * lsg_480[k];
    }

#pragma omp simd aligned(t_675, t_676, t_677, pa_x, pc_x, pc_y, ksh0_675, ksh0_677, ksg_377, \
                         ksg_483, ksg_485, ksh1_675, ksh1_677, \
                         lsg_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_675[k] = pa_x[k] * ksh0_675[k]
                   + f_11 * ksg_483[k]
                   - f_8 * pc_x[k] * ksh1_675[k];

        t_676[k] = f_11 * ksg_377[k]
                   + f_3 * pc_y[k] * lsg_482[k];

        t_677[k] = pa_x[k] * ksh0_677[k]
                   + f_11 * ksg_485[k]
                   - f_8 * pc_x[k] * ksh1_677[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pa_x, pc_x, pc_y, pc_z, ksh0_678, ksg_363, \
                         ksg_380, ksg_486, ksh1_678, lsg_483, lsg_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = pa_x[k] * ksh0_678[k]
                   + f_10 * ksg_486[k]
                   - f_8 * pc_x[k] * ksh1_678[k];

        t_679[k] = f_17 * ksg_363[k]
                   + f_3 * pc_z[k] * lsg_483[k];

        t_680[k] = f_11 * ksg_380[k]
                   + f_3 * pc_y[k] * lsg_485[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, pa_x, pc_x, ksh0_681, ksg_489, ksg_490, \
                         ksg_491, ksg_492, ksh1_681, lsg_490, lsg_491, \
                         lsg_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = pa_x[k] * ksh0_681[k]
                   + f_10 * ksg_489[k]
                   - f_8 * pc_x[k] * ksh1_681[k];

        t_682[k] = f_9 * ksg_490[k]
                   + f_3 * pc_x[k] * lsg_490[k];

        t_683[k] = f_9 * ksg_491[k]
                   + f_3 * pc_x[k] * lsg_491[k];

        t_684[k] = f_9 * ksg_492[k]
                   + f_3 * pc_x[k] * lsg_492[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, pa_x, pc_x, pc_z, ksh0_687, ksg_370, \
                         ksg_493, ksg_494, ksh1_687, lsg_490, lsg_493, \
                         lsg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_9 * ksg_493[k]
                   + f_3 * pc_x[k] * lsg_493[k];

        t_686[k] = f_9 * ksg_494[k]
                   + f_3 * pc_x[k] * lsg_494[k];

        t_687[k] = pa_x[k] * ksh0_687[k]
                   - f_8 * pc_x[k] * ksh1_687[k];

        t_688[k] = f_17 * ksg_370[k]
                   + f_3 * pc_z[k] * lsg_490[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pa_x, pc_x, pc_y, ksh0_689, ksh0_690, \
                         ksh0_692, ksg_389, ksh1_689, ksh1_690, ksh1_692, \
                         lsg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pa_x[k] * ksh0_689[k]
                   - f_8 * pc_x[k] * ksh1_689[k];

        t_690[k] = pa_x[k] * ksh0_690[k]
                   - f_8 * pc_x[k] * ksh1_690[k];

        t_691[k] = f_11 * ksg_389[k]
                   + f_3 * pc_y[k] * lsg_494[k];

        t_692[k] = pa_x[k] * ksh0_692[k]
                   - f_8 * pc_x[k] * ksh1_692[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pa_x, pc_x, pc_y, pc_z, ksh0_693, ksg_375, \
                         ksg_390, ksg_495, ksh1_693, lsg_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pa_x[k] * ksh0_693[k]
                   + f_16 * ksg_495[k]
                   - f_8 * pc_x[k] * ksh1_693[k];

        t_694[k] = f_10 * ksg_390[k]
                   + f_3 * pc_y[k] * lsg_495[k];

        t_695[k] = f_16 * ksg_375[k]
                   + f_3 * pc_z[k] * lsg_495[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pa_x, pc_x, pc_y, ksh0_696, ksh0_698, ksg_392, \
                         ksg_498, ksg_500, ksh1_696, ksh1_698, \
                         lsg_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = pa_x[k] * ksh0_696[k]
                   + f_11 * ksg_498[k]
                   - f_8 * pc_x[k] * ksh1_696[k];

        t_697[k] = f_10 * ksg_392[k]
                   + f_3 * pc_y[k] * lsg_497[k];

        t_698[k] = pa_x[k] * ksh0_698[k]
                   + f_11 * ksg_500[k]
                   - f_8 * pc_x[k] * ksh1_698[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pa_x, pc_x, pc_y, pc_z, ksh0_699, ksg_378, \
                         ksg_395, ksg_501, ksh1_699, lsg_498, lsg_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = pa_x[k] * ksh0_699[k]
                   + f_10 * ksg_501[k]
                   - f_8 * pc_x[k] * ksh1_699[k];

        t_700[k] = f_16 * ksg_378[k]
                   + f_3 * pc_z[k] * lsg_498[k];

        t_701[k] = f_10 * ksg_395[k]
                   + f_3 * pc_y[k] * lsg_500[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pa_x, pc_x, ksh0_702, ksg_504, ksg_505, \
                         ksg_506, ksg_507, ksh1_702, lsg_505, lsg_506, \
                         lsg_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = pa_x[k] * ksh0_702[k]
                   + f_10 * ksg_504[k]
                   - f_8 * pc_x[k] * ksh1_702[k];

        t_703[k] = f_9 * ksg_505[k]
                   + f_3 * pc_x[k] * lsg_505[k];

        t_704[k] = f_9 * ksg_506[k]
                   + f_3 * pc_x[k] * lsg_506[k];

        t_705[k] = f_9 * ksg_507[k]
                   + f_3 * pc_x[k] * lsg_507[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pa_x, pc_x, pc_z, ksh0_708, ksg_385, \
                         ksg_508, ksg_509, ksh1_708, lsg_505, lsg_508, \
                         lsg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_9 * ksg_508[k]
                   + f_3 * pc_x[k] * lsg_508[k];

        t_707[k] = f_9 * ksg_509[k]
                   + f_3 * pc_x[k] * lsg_509[k];

        t_708[k] = pa_x[k] * ksh0_708[k]
                   - f_8 * pc_x[k] * ksh1_708[k];

        t_709[k] = f_16 * ksg_385[k]
                   + f_3 * pc_z[k] * lsg_505[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pa_x, pc_x, pc_y, ksh0_710, ksh0_711, \
                         ksh0_713, ksg_404, ksh1_710, ksh1_711, ksh1_713, \
                         lsg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = pa_x[k] * ksh0_710[k]
                   - f_8 * pc_x[k] * ksh1_710[k];

        t_711[k] = pa_x[k] * ksh0_711[k]
                   - f_8 * pc_x[k] * ksh1_711[k];

        t_712[k] = f_10 * ksg_404[k]
                   + f_3 * pc_y[k] * lsg_509[k];

        t_713[k] = pa_x[k] * ksh0_713[k]
                   - f_8 * pc_x[k] * ksh1_713[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pa_y, pc_y, pc_z, ksh0_567, ksg_390, ksg_405, \
                         ksh1_567, lsg_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pa_y[k] * ksh0_567[k]
                   - f_8 * pc_y[k] * ksh1_567[k];

        t_715[k] = f_9 * ksg_405[k]
                   + f_3 * pc_y[k] * lsg_510[k];

        t_716[k] = f_15 * ksg_390[k]
                   + f_3 * pc_z[k] * lsg_510[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pa_x, pa_y, pc_x, pc_y, ksh0_572, ksh0_717, \
                         ksg_407, ksg_513, ksh1_572, ksh1_717, \
                         lsg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = pa_x[k] * ksh0_717[k]
                   + f_11 * ksg_513[k]
                   - f_8 * pc_x[k] * ksh1_717[k];

        t_718[k] = f_9 * ksg_407[k]
                   + f_3 * pc_y[k] * lsg_512[k];

        t_719[k] = pa_y[k] * ksh0_572[k]
                   - f_8 * pc_y[k] * ksh1_572[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pa_x, pc_x, pc_y, pc_z, ksh0_720, ksg_393, \
                         ksg_410, ksg_516, ksh1_720, lsg_513, lsg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = pa_x[k] * ksh0_720[k]
                   + f_10 * ksg_516[k]
                   - f_8 * pc_x[k] * ksh1_720[k];

        t_721[k] = f_15 * ksg_393[k]
                   + f_3 * pc_z[k] * lsg_513[k];

        t_722[k] = f_9 * ksg_410[k]
                   + f_3 * pc_y[k] * lsg_515[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, t_726, pa_y, pc_x, pc_y, ksh0_576, ksg_520, \
                         ksg_521, ksg_522, ksh1_576, lsg_520, lsg_521, \
                         lsg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = pa_y[k] * ksh0_576[k]
                   - f_8 * pc_y[k] * ksh1_576[k];

        t_724[k] = f_9 * ksg_520[k]
                   + f_3 * pc_x[k] * lsg_520[k];

        t_725[k] = f_9 * ksg_521[k]
                   + f_3 * pc_x[k] * lsg_521[k];

        t_726[k] = f_9 * ksg_522[k]
                   + f_3 * pc_x[k] * lsg_522[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_780 = buffer.data(target + 780);
    auto *t_781 = buffer.data(target + 781);
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);
    auto *t_791 = buffer.data(target + 791);
    auto *t_792 = buffer.data(target + 792);
    auto *t_793 = buffer.data(target + 793);
    auto *t_794 = buffer.data(target + 794);
    auto *t_795 = buffer.data(target + 795);
    auto *t_796 = buffer.data(target + 796);
    auto *t_797 = buffer.data(target + 797);
    auto *t_798 = buffer.data(target + 798);
    auto *t_799 = buffer.data(target + 799);
    auto *t_800 = buffer.data(target + 800);
    auto *t_801 = buffer.data(target + 801);
    auto *t_802 = buffer.data(target + 802);
    auto *t_803 = buffer.data(target + 803);
    auto *t_804 = buffer.data(target + 804);
    auto *t_805 = buffer.data(target + 805);
    auto *t_806 = buffer.data(target + 806);
    auto *t_807 = buffer.data(target + 807);
    auto *t_808 = buffer.data(target + 808);
    auto *t_809 = buffer.data(target + 809);
    auto *t_810 = buffer.data(target + 810);
    auto *t_811 = buffer.data(target + 811);
    auto *t_812 = buffer.data(target + 812);
    auto *t_813 = buffer.data(target + 813);
    auto *t_814 = buffer.data(target + 814);
    auto *t_815 = buffer.data(target + 815);
    auto *t_816 = buffer.data(target + 816);
    auto *t_817 = buffer.data(target + 817);
    auto *t_818 = buffer.data(target + 818);
    auto *t_819 = buffer.data(target + 819);
    auto *t_820 = buffer.data(target + 820);
    auto *t_821 = buffer.data(target + 821);
    auto *t_822 = buffer.data(target + 822);
    auto *t_823 = buffer.data(target + 823);
    auto *t_824 = buffer.data(target + 824);
    auto *t_825 = buffer.data(target + 825);
    auto *t_826 = buffer.data(target + 826);
    auto *t_827 = buffer.data(target + 827);
    auto *t_828 = buffer.data(target + 828);
    auto *t_829 = buffer.data(target + 829);
    auto *t_830 = buffer.data(target + 830);
    auto *t_831 = buffer.data(target + 831);
    auto *t_832 = buffer.data(target + 832);
    auto *t_833 = buffer.data(target + 833);
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);
    auto *t_840 = buffer.data(target + 840);
    auto *t_841 = buffer.data(target + 841);
    auto *t_842 = buffer.data(target + 842);
    auto *t_843 = buffer.data(target + 843);
    auto *t_844 = buffer.data(target + 844);
    auto *t_845 = buffer.data(target + 845);
    auto *t_846 = buffer.data(target + 846);
    auto *t_847 = buffer.data(target + 847);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_588 = buffer.data(ksh0 + 588);
    const auto *ksh0_589 = buffer.data(ksh0 + 589);
    const auto *ksh0_591 = buffer.data(ksh0 + 591);
    const auto *ksh0_594 = buffer.data(ksh0 + 594);
    const auto *ksh0_603 = buffer.data(ksh0 + 603);
    const auto *ksh0_605 = buffer.data(ksh0 + 605);
    const auto *ksh0_606 = buffer.data(ksh0 + 606);
    const auto *ksh0_729 = buffer.data(ksh0 + 729);
    const auto *ksh0_731 = buffer.data(ksh0 + 731);
    const auto *ksh0_732 = buffer.data(ksh0 + 732);
    const auto *ksh0_734 = buffer.data(ksh0 + 734);
    const auto *ksh0_735 = buffer.data(ksh0 + 735);
    const auto *ksh0_740 = buffer.data(ksh0 + 740);
    const auto *ksh0_744 = buffer.data(ksh0 + 744);
    const auto *ksh0_750 = buffer.data(ksh0 + 750);
    const auto *ksh0_751 = buffer.data(ksh0 + 751);
    const auto *ksh0_752 = buffer.data(ksh0 + 752);
    const auto *ksh0_753 = buffer.data(ksh0 + 753);
    const auto *ksh0_755 = buffer.data(ksh0 + 755);

    const auto *ksg_400 = buffer.data(ksg + 400);
    const auto *ksg_405 = buffer.data(ksg + 405);
    const auto *ksg_419 = buffer.data(ksg + 419);
    const auto *ksg_430 = buffer.data(ksg + 430);
    const auto *ksg_431 = buffer.data(ksg + 431);
    const auto *ksg_432 = buffer.data(ksg + 432);
    const auto *ksg_434 = buffer.data(ksg + 434);
    const auto *ksg_445 = buffer.data(ksg + 445);
    const auto *ksg_449 = buffer.data(ksg + 449);
    const auto *ksg_460 = buffer.data(ksg + 460);
    const auto *ksg_462 = buffer.data(ksg + 462);
    const auto *ksg_463 = buffer.data(ksg + 463);
    const auto *ksg_464 = buffer.data(ksg + 464);
    const auto *ksg_475 = buffer.data(ksg + 475);
    const auto *ksg_477 = buffer.data(ksg + 477);
    const auto *ksg_478 = buffer.data(ksg + 478);
    const auto *ksg_479 = buffer.data(ksg + 479);
    const auto *ksg_523 = buffer.data(ksg + 523);
    const auto *ksg_524 = buffer.data(ksg + 524);
    const auto *ksg_525 = buffer.data(ksg + 525);
    const auto *ksg_530 = buffer.data(ksg + 530);
    const auto *ksg_534 = buffer.data(ksg + 534);
    const auto *ksg_535 = buffer.data(ksg + 535);
    const auto *ksg_536 = buffer.data(ksg + 536);
    const auto *ksg_537 = buffer.data(ksg + 537);
    const auto *ksg_539 = buffer.data(ksg + 539);

    const auto *ksh1_588 = buffer.data(ksh1 + 588);
    const auto *ksh1_589 = buffer.data(ksh1 + 589);
    const auto *ksh1_591 = buffer.data(ksh1 + 591);
    const auto *ksh1_594 = buffer.data(ksh1 + 594);
    const auto *ksh1_603 = buffer.data(ksh1 + 603);
    const auto *ksh1_605 = buffer.data(ksh1 + 605);
    const auto *ksh1_606 = buffer.data(ksh1 + 606);
    const auto *ksh1_729 = buffer.data(ksh1 + 729);
    const auto *ksh1_731 = buffer.data(ksh1 + 731);
    const auto *ksh1_732 = buffer.data(ksh1 + 732);
    const auto *ksh1_734 = buffer.data(ksh1 + 734);
    const auto *ksh1_735 = buffer.data(ksh1 + 735);
    const auto *ksh1_740 = buffer.data(ksh1 + 740);
    const auto *ksh1_744 = buffer.data(ksh1 + 744);
    const auto *ksh1_750 = buffer.data(ksh1 + 750);
    const auto *ksh1_751 = buffer.data(ksh1 + 751);
    const auto *ksh1_752 = buffer.data(ksh1 + 752);
    const auto *ksh1_753 = buffer.data(ksh1 + 753);
    const auto *ksh1_755 = buffer.data(ksh1 + 755);

    const auto *lsf0_350 = buffer.data(lsf0 + 350);
    const auto *lsf0_351 = buffer.data(lsf0 + 351);
    const auto *lsf0_352 = buffer.data(lsf0 + 352);
    const auto *lsf0_360 = buffer.data(lsf0 + 360);
    const auto *lsf0_361 = buffer.data(lsf0 + 361);
    const auto *lsf0_363 = buffer.data(lsf0 + 363);
    const auto *lsf0_365 = buffer.data(lsf0 + 365);
    const auto *lsf0_366 = buffer.data(lsf0 + 366);
    const auto *lsf0_367 = buffer.data(lsf0 + 367);
    const auto *lsf0_368 = buffer.data(lsf0 + 368);
    const auto *lsf0_369 = buffer.data(lsf0 + 369);
    const auto *lsf0_372 = buffer.data(lsf0 + 372);
    const auto *lsf0_374 = buffer.data(lsf0 + 374);
    const auto *lsf0_375 = buffer.data(lsf0 + 375);
    const auto *lsf0_377 = buffer.data(lsf0 + 377);
    const auto *lsf0_378 = buffer.data(lsf0 + 378);
    const auto *lsf0_379 = buffer.data(lsf0 + 379);
    const auto *lsf0_380 = buffer.data(lsf0 + 380);
    const auto *lsf0_381 = buffer.data(lsf0 + 381);
    const auto *lsf0_382 = buffer.data(lsf0 + 382);
    const auto *lsf0_383 = buffer.data(lsf0 + 383);
    const auto *lsf0_384 = buffer.data(lsf0 + 384);
    const auto *lsf0_385 = buffer.data(lsf0 + 385);
    const auto *lsf0_386 = buffer.data(lsf0 + 386);
    const auto *lsf0_387 = buffer.data(lsf0 + 387);
    const auto *lsf0_388 = buffer.data(lsf0 + 388);
    const auto *lsf0_389 = buffer.data(lsf0 + 389);
    const auto *lsf0_390 = buffer.data(lsf0 + 390);
    const auto *lsf0_391 = buffer.data(lsf0 + 391);
    const auto *lsf0_392 = buffer.data(lsf0 + 392);
    const auto *lsf0_393 = buffer.data(lsf0 + 393);
    const auto *lsf0_394 = buffer.data(lsf0 + 394);
    const auto *lsf0_395 = buffer.data(lsf0 + 395);
    const auto *lsf0_396 = buffer.data(lsf0 + 396);
    const auto *lsf0_397 = buffer.data(lsf0 + 397);
    const auto *lsf0_398 = buffer.data(lsf0 + 398);
    const auto *lsf0_399 = buffer.data(lsf0 + 399);
    const auto *lsf0_400 = buffer.data(lsf0 + 400);
    const auto *lsf0_401 = buffer.data(lsf0 + 401);
    const auto *lsf0_402 = buffer.data(lsf0 + 402);
    const auto *lsf0_403 = buffer.data(lsf0 + 403);
    const auto *lsf0_404 = buffer.data(lsf0 + 404);
    const auto *lsf0_405 = buffer.data(lsf0 + 405);
    const auto *lsf0_406 = buffer.data(lsf0 + 406);
    const auto *lsf0_407 = buffer.data(lsf0 + 407);

    const auto *lsf1_350 = buffer.data(lsf1 + 350);
    const auto *lsf1_351 = buffer.data(lsf1 + 351);
    const auto *lsf1_352 = buffer.data(lsf1 + 352);
    const auto *lsf1_360 = buffer.data(lsf1 + 360);
    const auto *lsf1_361 = buffer.data(lsf1 + 361);
    const auto *lsf1_363 = buffer.data(lsf1 + 363);
    const auto *lsf1_365 = buffer.data(lsf1 + 365);
    const auto *lsf1_366 = buffer.data(lsf1 + 366);
    const auto *lsf1_367 = buffer.data(lsf1 + 367);
    const auto *lsf1_368 = buffer.data(lsf1 + 368);
    const auto *lsf1_369 = buffer.data(lsf1 + 369);
    const auto *lsf1_372 = buffer.data(lsf1 + 372);
    const auto *lsf1_374 = buffer.data(lsf1 + 374);
    const auto *lsf1_375 = buffer.data(lsf1 + 375);
    const auto *lsf1_377 = buffer.data(lsf1 + 377);
    const auto *lsf1_378 = buffer.data(lsf1 + 378);
    const auto *lsf1_379 = buffer.data(lsf1 + 379);
    const auto *lsf1_380 = buffer.data(lsf1 + 380);
    const auto *lsf1_381 = buffer.data(lsf1 + 381);
    const auto *lsf1_382 = buffer.data(lsf1 + 382);
    const auto *lsf1_383 = buffer.data(lsf1 + 383);
    const auto *lsf1_384 = buffer.data(lsf1 + 384);
    const auto *lsf1_385 = buffer.data(lsf1 + 385);
    const auto *lsf1_386 = buffer.data(lsf1 + 386);
    const auto *lsf1_387 = buffer.data(lsf1 + 387);
    const auto *lsf1_388 = buffer.data(lsf1 + 388);
    const auto *lsf1_389 = buffer.data(lsf1 + 389);
    const auto *lsf1_390 = buffer.data(lsf1 + 390);
    const auto *lsf1_391 = buffer.data(lsf1 + 391);
    const auto *lsf1_392 = buffer.data(lsf1 + 392);
    const auto *lsf1_393 = buffer.data(lsf1 + 393);
    const auto *lsf1_394 = buffer.data(lsf1 + 394);
    const auto *lsf1_395 = buffer.data(lsf1 + 395);
    const auto *lsf1_396 = buffer.data(lsf1 + 396);
    const auto *lsf1_397 = buffer.data(lsf1 + 397);
    const auto *lsf1_398 = buffer.data(lsf1 + 398);
    const auto *lsf1_399 = buffer.data(lsf1 + 399);
    const auto *lsf1_400 = buffer.data(lsf1 + 400);
    const auto *lsf1_401 = buffer.data(lsf1 + 401);
    const auto *lsf1_402 = buffer.data(lsf1 + 402);
    const auto *lsf1_403 = buffer.data(lsf1 + 403);
    const auto *lsf1_404 = buffer.data(lsf1 + 404);
    const auto *lsf1_405 = buffer.data(lsf1 + 405);
    const auto *lsf1_406 = buffer.data(lsf1 + 406);
    const auto *lsf1_407 = buffer.data(lsf1 + 407);

    const auto *lsg_520 = buffer.data(lsg + 520);
    const auto *lsg_523 = buffer.data(lsg + 523);
    const auto *lsg_524 = buffer.data(lsg + 524);
    const auto *lsg_525 = buffer.data(lsg + 525);
    const auto *lsg_526 = buffer.data(lsg + 526);
    const auto *lsg_527 = buffer.data(lsg + 527);
    const auto *lsg_528 = buffer.data(lsg + 528);
    const auto *lsg_529 = buffer.data(lsg + 529);
    const auto *lsg_530 = buffer.data(lsg + 530);
    const auto *lsg_534 = buffer.data(lsg + 534);
    const auto *lsg_535 = buffer.data(lsg + 535);
    const auto *lsg_536 = buffer.data(lsg + 536);
    const auto *lsg_537 = buffer.data(lsg + 537);
    const auto *lsg_539 = buffer.data(lsg + 539);
    const auto *lsg_540 = buffer.data(lsg + 540);
    const auto *lsg_541 = buffer.data(lsg + 541);
    const auto *lsg_543 = buffer.data(lsg + 543);
    const auto *lsg_545 = buffer.data(lsg + 545);
    const auto *lsg_546 = buffer.data(lsg + 546);
    const auto *lsg_548 = buffer.data(lsg + 548);
    const auto *lsg_549 = buffer.data(lsg + 549);
    const auto *lsg_550 = buffer.data(lsg + 550);
    const auto *lsg_551 = buffer.data(lsg + 551);
    const auto *lsg_552 = buffer.data(lsg + 552);
    const auto *lsg_553 = buffer.data(lsg + 553);
    const auto *lsg_554 = buffer.data(lsg + 554);
    const auto *lsg_557 = buffer.data(lsg + 557);
    const auto *lsg_559 = buffer.data(lsg + 559);
    const auto *lsg_560 = buffer.data(lsg + 560);
    const auto *lsg_562 = buffer.data(lsg + 562);
    const auto *lsg_563 = buffer.data(lsg + 563);
    const auto *lsg_564 = buffer.data(lsg + 564);
    const auto *lsg_565 = buffer.data(lsg + 565);
    const auto *lsg_566 = buffer.data(lsg + 566);
    const auto *lsg_567 = buffer.data(lsg + 567);
    const auto *lsg_568 = buffer.data(lsg + 568);
    const auto *lsg_569 = buffer.data(lsg + 569);
    const auto *lsg_570 = buffer.data(lsg + 570);
    const auto *lsg_571 = buffer.data(lsg + 571);
    const auto *lsg_572 = buffer.data(lsg + 572);
    const auto *lsg_573 = buffer.data(lsg + 573);
    const auto *lsg_574 = buffer.data(lsg + 574);
    const auto *lsg_575 = buffer.data(lsg + 575);
    const auto *lsg_576 = buffer.data(lsg + 576);
    const auto *lsg_577 = buffer.data(lsg + 577);
    const auto *lsg_578 = buffer.data(lsg + 578);
    const auto *lsg_579 = buffer.data(lsg + 579);
    const auto *lsg_580 = buffer.data(lsg + 580);
    const auto *lsg_581 = buffer.data(lsg + 581);
    const auto *lsg_582 = buffer.data(lsg + 582);
    const auto *lsg_583 = buffer.data(lsg + 583);
    const auto *lsg_584 = buffer.data(lsg + 584);
    const auto *lsg_585 = buffer.data(lsg + 585);
    const auto *lsg_586 = buffer.data(lsg + 586);
    const auto *lsg_587 = buffer.data(lsg + 587);
    const auto *lsg_588 = buffer.data(lsg + 588);
    const auto *lsg_589 = buffer.data(lsg + 589);
    const auto *lsg_590 = buffer.data(lsg + 590);
    const auto *lsg_591 = buffer.data(lsg + 591);
    const auto *lsg_592 = buffer.data(lsg + 592);
    const auto *lsg_593 = buffer.data(lsg + 593);
    const auto *lsg_594 = buffer.data(lsg + 594);
    const auto *lsg_595 = buffer.data(lsg + 595);
    const auto *lsg_596 = buffer.data(lsg + 596);
    const auto *lsg_597 = buffer.data(lsg + 597);
    const auto *lsg_598 = buffer.data(lsg + 598);
    const auto *lsg_599 = buffer.data(lsg + 599);
    const auto *lsg_600 = buffer.data(lsg + 600);
    const auto *lsg_601 = buffer.data(lsg + 601);
    const auto *lsg_602 = buffer.data(lsg + 602);
    const auto *lsg_603 = buffer.data(lsg + 603);
    const auto *lsg_604 = buffer.data(lsg + 604);
    const auto *lsg_605 = buffer.data(lsg + 605);
    const auto *lsg_606 = buffer.data(lsg + 606);
    const auto *lsg_607 = buffer.data(lsg + 607);

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pa_x, pc_x, pc_z, ksh0_729, ksg_400, \
                         ksg_523, ksg_524, ksh1_729, lsg_520, lsg_523, \
                         lsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_9 * ksg_523[k]
                   + f_3 * pc_x[k] * lsg_523[k];

        t_728[k] = f_9 * ksg_524[k]
                   + f_3 * pc_x[k] * lsg_524[k];

        t_729[k] = pa_x[k] * ksh0_729[k]
                   - f_8 * pc_x[k] * ksh1_729[k];

        t_730[k] = f_15 * ksg_400[k]
                   + f_3 * pc_z[k] * lsg_520[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, t_734, pa_x, pc_x, pc_y, ksh0_731, ksh0_732, \
                         ksh0_734, ksg_419, ksh1_731, ksh1_732, ksh1_734, \
                         lsg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = pa_x[k] * ksh0_731[k]
                   - f_8 * pc_x[k] * ksh1_731[k];

        t_732[k] = pa_x[k] * ksh0_732[k]
                   - f_8 * pc_x[k] * ksh1_732[k];

        t_733[k] = f_9 * ksg_419[k]
                   + f_3 * pc_y[k] * lsg_524[k];

        t_734[k] = pa_x[k] * ksh0_734[k]
                   - f_8 * pc_x[k] * ksh1_734[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, pa_x, pc_x, pc_y, pc_z, ksh0_735, \
                         ksg_405, ksg_525, ksh1_735, lsf0_350, lsf1_350, lsg_525, \
                         lsg_526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = pa_x[k] * ksh0_735[k]
                   + f_16 * ksg_525[k]
                   - f_8 * pc_x[k] * ksh1_735[k];

        t_736[k] = f_3 * pc_y[k] * lsg_525[k];

        t_737[k] = f_12 * ksg_405[k]
                   + f_3 * pc_z[k] * lsg_525[k];

        t_738[k] = f_4 * lsf0_350[k]
                   - f_5 * lsf1_350[k]
                   + f_3 * pc_y[k] * lsg_526[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pa_x, pc_x, pc_y, ksh0_740, ksg_530, ksh1_740, \
                         lsf0_351, lsf1_351, lsg_527, lsg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_3 * pc_y[k] * lsg_527[k];

        t_740[k] = pa_x[k] * ksh0_740[k]
                   + f_11 * ksg_530[k]
                   - f_8 * pc_x[k] * ksh1_740[k];

        t_741[k] = f_6 * lsf0_351[k]
                   - f_7 * lsf1_351[k]
                   + f_3 * pc_y[k] * lsg_528[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pa_x, pc_x, pc_y, ksh0_744, ksg_534, \
                         ksg_535, ksh1_744, lsf0_352, lsf1_352, lsg_529, lsg_530, \
                         lsg_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_4 * lsf0_352[k]
                   - f_5 * lsf1_352[k]
                   + f_3 * pc_y[k] * lsg_529[k];

        t_743[k] = f_3 * pc_y[k] * lsg_530[k];

        t_744[k] = pa_x[k] * ksh0_744[k]
                   + f_10 * ksg_534[k]
                   - f_8 * pc_x[k] * ksh1_744[k];

        t_745[k] = f_9 * ksg_535[k]
                   + f_3 * pc_x[k] * lsg_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, ksg_536, ksg_537, ksg_539, \
                         lsg_534, lsg_536, lsg_537, lsg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_9 * ksg_536[k]
                   + f_3 * pc_x[k] * lsg_536[k];

        t_747[k] = f_9 * ksg_537[k]
                   + f_3 * pc_x[k] * lsg_537[k];

        t_748[k] = f_3 * pc_y[k] * lsg_534[k];

        t_749[k] = f_9 * ksg_539[k]
                   + f_3 * pc_x[k] * lsg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pa_x, pc_x, ksh0_750, ksh0_751, ksh0_752, \
                         ksh0_753, ksh1_750, ksh1_751, ksh1_752, \
                         ksh1_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = pa_x[k] * ksh0_750[k]
                   - f_8 * pc_x[k] * ksh1_750[k];

        t_751[k] = pa_x[k] * ksh0_751[k]
                   - f_8 * pc_x[k] * ksh1_751[k];

        t_752[k] = pa_x[k] * ksh0_752[k]
                   - f_8 * pc_x[k] * ksh1_752[k];

        t_753[k] = pa_x[k] * ksh0_753[k]
                   - f_8 * pc_x[k] * ksh1_753[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pa_x, pc_x, pc_y, ksh0_755, ksh1_755, \
                         lsf0_360, lsf0_361, lsf1_360, lsf1_361, lsg_539, lsg_540, \
                         lsg_541 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_3 * pc_y[k] * lsg_539[k];

        t_755[k] = pa_x[k] * ksh0_755[k]
                   - f_8 * pc_x[k] * ksh1_755[k];

        t_756[k] = f_1 * lsf0_360[k]
                   - f_2 * lsf1_360[k]
                   + f_3 * pc_x[k] * lsg_540[k];

        t_757[k] = f_13 * lsf0_361[k]
                   - f_14 * lsf1_361[k]
                   + f_3 * pc_x[k] * lsg_541[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, t_761, pc_x, pc_z, lsf0_363, lsf0_365, lsf1_363, \
                         lsf1_365, lsg_540, lsg_541, lsg_543, lsg_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_3 * pc_z[k] * lsg_540[k];

        t_759[k] = f_6 * lsf0_363[k]
                   - f_7 * lsf1_363[k]
                   + f_3 * pc_x[k] * lsg_543[k];

        t_760[k] = f_3 * pc_z[k] * lsg_541[k];

        t_761[k] = f_6 * lsf0_365[k]
                   - f_7 * lsf1_365[k]
                   + f_3 * pc_x[k] * lsg_545[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, t_765, pc_x, pc_z, lsf0_366, lsf0_368, lsf0_369, \
                         lsf1_366, lsf1_368, lsf1_369, lsg_543, lsg_546, lsg_548, \
                         lsg_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_4 * lsf0_366[k]
                   - f_5 * lsf1_366[k]
                   + f_3 * pc_x[k] * lsg_546[k];

        t_763[k] = f_3 * pc_z[k] * lsg_543[k];

        t_764[k] = f_4 * lsf0_368[k]
                   - f_5 * lsf1_368[k]
                   + f_3 * pc_x[k] * lsg_548[k];

        t_765[k] = f_4 * lsf0_369[k]
                   - f_5 * lsf1_369[k]
                   + f_3 * pc_x[k] * lsg_549[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, t_770, t_771, pc_x, pc_y, ksg_430, \
                         lsf0_366, lsf1_366, lsg_550, lsg_551, lsg_552, lsg_553, \
                         lsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_3 * pc_x[k] * lsg_550[k];

        t_767[k] = f_3 * pc_x[k] * lsg_551[k];

        t_768[k] = f_3 * pc_x[k] * lsg_552[k];

        t_769[k] = f_3 * pc_x[k] * lsg_553[k];

        t_770[k] = f_3 * pc_x[k] * lsg_554[k];

        t_771[k] = f_0 * ksg_430[k]
                   + f_1 * lsf0_366[k]
                   - f_2 * lsf1_366[k]
                   + f_3 * pc_y[k] * lsg_550[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, pc_y, pc_z, ksg_434, lsf0_366, lsf0_367, \
                         lsf1_366, lsf1_367, lsg_550, lsg_551, lsg_552, \
                         lsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_3 * pc_z[k] * lsg_550[k];

        t_773[k] = f_4 * lsf0_366[k]
                   - f_5 * lsf1_366[k]
                   + f_3 * pc_z[k] * lsg_551[k];

        t_774[k] = f_6 * lsf0_367[k]
                   - f_7 * lsf1_367[k]
                   + f_3 * pc_z[k] * lsg_552[k];

        t_775[k] = f_0 * ksg_434[k]
                   + f_3 * pc_y[k] * lsg_554[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pa_z, pc_z, ksh0_588, ksh0_589, ksh1_588, \
                         ksh1_589, lsf0_369, lsf1_369, lsg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_1 * lsf0_369[k]
                   - f_2 * lsf1_369[k]
                   + f_3 * pc_z[k] * lsg_554[k];

        t_777[k] = pa_z[k] * ksh0_588[k]
                   - f_8 * pc_z[k] * ksh1_588[k];

        t_778[k] = pa_z[k] * ksh0_589[k]
                   - f_8 * pc_z[k] * ksh1_589[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pa_z, pc_x, pc_z, ksh0_591, ksh1_591, lsf0_372, \
                         lsf0_374, lsf1_372, lsf1_374, lsg_557, \
                         lsg_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_13 * lsf0_372[k]
                   - f_14 * lsf1_372[k]
                   + f_3 * pc_x[k] * lsg_557[k];

        t_780[k] = pa_z[k] * ksh0_591[k]
                   - f_8 * pc_z[k] * ksh1_591[k];

        t_781[k] = f_6 * lsf0_374[k]
                   - f_7 * lsf1_374[k]
                   + f_3 * pc_x[k] * lsg_559[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pa_z, pc_x, pc_z, ksh0_594, ksh1_594, lsf0_375, \
                         lsf0_377, lsf1_375, lsf1_377, lsg_560, \
                         lsg_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_6 * lsf0_375[k]
                   - f_7 * lsf1_375[k]
                   + f_3 * pc_x[k] * lsg_560[k];

        t_783[k] = pa_z[k] * ksh0_594[k]
                   - f_8 * pc_z[k] * ksh1_594[k];

        t_784[k] = f_4 * lsf0_377[k]
                   - f_5 * lsf1_377[k]
                   + f_3 * pc_x[k] * lsg_562[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, t_789, pc_x, lsf0_378, lsf0_379, \
                         lsf1_378, lsf1_379, lsg_563, lsg_564, lsg_565, lsg_566, \
                         lsg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_4 * lsf0_378[k]
                   - f_5 * lsf1_378[k]
                   + f_3 * pc_x[k] * lsg_563[k];

        t_786[k] = f_4 * lsf0_379[k]
                   - f_5 * lsf1_379[k]
                   + f_3 * pc_x[k] * lsg_564[k];

        t_787[k] = f_3 * pc_x[k] * lsg_565[k];

        t_788[k] = f_3 * pc_x[k] * lsg_566[k];

        t_789[k] = f_3 * pc_x[k] * lsg_567[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pa_z, pc_x, pc_z, ksh0_603, ksg_430, \
                         ksh1_603, lsg_565, lsg_568, lsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_3 * pc_x[k] * lsg_568[k];

        t_791[k] = f_3 * pc_x[k] * lsg_569[k];

        t_792[k] = pa_z[k] * ksh0_603[k]
                   - f_8 * pc_z[k] * ksh1_603[k];

        t_793[k] = f_9 * ksg_430[k]
                   + f_3 * pc_z[k] * lsg_565[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, pa_z, pc_y, pc_z, ksh0_605, ksh0_606, ksg_431, \
                         ksg_432, ksg_449, ksh1_605, ksh1_606, \
                         lsg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = pa_z[k] * ksh0_605[k]
                   + f_10 * ksg_431[k]
                   - f_8 * pc_z[k] * ksh1_605[k];

        t_795[k] = pa_z[k] * ksh0_606[k]
                   + f_11 * ksg_432[k]
                   - f_8 * pc_z[k] * ksh1_606[k];

        t_796[k] = f_12 * ksg_449[k]
                   + f_3 * pc_y[k] * lsg_569[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, pc_x, pc_z, ksg_434, lsf0_379, lsf0_380, \
                         lsf0_381, lsf1_379, lsf1_380, lsf1_381, lsg_569, lsg_570, \
                         lsg_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_9 * ksg_434[k]
                   + f_1 * lsf0_379[k]
                   - f_2 * lsf1_379[k]
                   + f_3 * pc_z[k] * lsg_569[k];

        t_798[k] = f_1 * lsf0_380[k]
                   - f_2 * lsf1_380[k]
                   + f_3 * pc_x[k] * lsg_570[k];

        t_799[k] = f_13 * lsf0_381[k]
                   - f_14 * lsf1_381[k]
                   + f_3 * pc_x[k] * lsg_571[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, lsf0_382, lsf0_383, lsf0_384, lsf1_382, \
                         lsf1_383, lsf1_384, lsg_572, lsg_573, \
                         lsg_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_13 * lsf0_382[k]
                   - f_14 * lsf1_382[k]
                   + f_3 * pc_x[k] * lsg_572[k];

        t_801[k] = f_6 * lsf0_383[k]
                   - f_7 * lsf1_383[k]
                   + f_3 * pc_x[k] * lsg_573[k];

        t_802[k] = f_6 * lsf0_384[k]
                   - f_7 * lsf1_384[k]
                   + f_3 * pc_x[k] * lsg_574[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_x, lsf0_385, lsf0_386, lsf0_387, lsf1_385, \
                         lsf1_386, lsf1_387, lsg_575, lsg_576, \
                         lsg_577 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_6 * lsf0_385[k]
                   - f_7 * lsf1_385[k]
                   + f_3 * pc_x[k] * lsg_575[k];

        t_804[k] = f_4 * lsf0_386[k]
                   - f_5 * lsf1_386[k]
                   + f_3 * pc_x[k] * lsg_576[k];

        t_805[k] = f_4 * lsf0_387[k]
                   - f_5 * lsf1_387[k]
                   + f_3 * pc_x[k] * lsg_577[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, t_809, t_810, pc_x, lsf0_388, lsf0_389, \
                         lsf1_388, lsf1_389, lsg_578, lsg_579, lsg_580, lsg_581, \
                         lsg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_4 * lsf0_388[k]
                   - f_5 * lsf1_388[k]
                   + f_3 * pc_x[k] * lsg_578[k];

        t_807[k] = f_4 * lsf0_389[k]
                   - f_5 * lsf1_389[k]
                   + f_3 * pc_x[k] * lsg_579[k];

        t_808[k] = f_3 * pc_x[k] * lsg_580[k];

        t_809[k] = f_3 * pc_x[k] * lsg_581[k];

        t_810[k] = f_3 * pc_x[k] * lsg_582[k];
    }

#pragma omp simd aligned(t_811, t_812, t_813, t_814, pc_x, pc_y, pc_z, ksg_445, ksg_460, \
                         lsf0_386, lsf1_386, lsg_580, lsg_583, \
                         lsg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_811[k] = f_3 * pc_x[k] * lsg_583[k];

        t_812[k] = f_3 * pc_x[k] * lsg_584[k];

        t_813[k] = f_15 * ksg_460[k]
                   + f_1 * lsf0_386[k]
                   - f_2 * lsf1_386[k]
                   + f_3 * pc_y[k] * lsg_580[k];

        t_814[k] = f_10 * ksg_445[k]
                   + f_3 * pc_z[k] * lsg_580[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_y, ksg_462, ksg_463, ksg_464, lsf0_388, \
                         lsf0_389, lsf1_388, lsf1_389, lsg_582, lsg_583, \
                         lsg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_15 * ksg_462[k]
                   + f_6 * lsf0_388[k]
                   - f_7 * lsf1_388[k]
                   + f_3 * pc_y[k] * lsg_582[k];

        t_816[k] = f_15 * ksg_463[k]
                   + f_4 * lsf0_389[k]
                   - f_5 * lsf1_389[k]
                   + f_3 * pc_y[k] * lsg_583[k];

        t_817[k] = f_15 * ksg_464[k]
                   + f_3 * pc_y[k] * lsg_584[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_z, ksg_449, lsf0_389, lsf0_390, \
                         lsf0_391, lsf1_389, lsf1_390, lsf1_391, lsg_584, lsg_585, \
                         lsg_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_10 * ksg_449[k]
                   + f_1 * lsf0_389[k]
                   - f_2 * lsf1_389[k]
                   + f_3 * pc_z[k] * lsg_584[k];

        t_819[k] = f_1 * lsf0_390[k]
                   - f_2 * lsf1_390[k]
                   + f_3 * pc_x[k] * lsg_585[k];

        t_820[k] = f_13 * lsf0_391[k]
                   - f_14 * lsf1_391[k]
                   + f_3 * pc_x[k] * lsg_586[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, lsf0_392, lsf0_393, lsf0_394, lsf1_392, \
                         lsf1_393, lsf1_394, lsg_587, lsg_588, \
                         lsg_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_13 * lsf0_392[k]
                   - f_14 * lsf1_392[k]
                   + f_3 * pc_x[k] * lsg_587[k];

        t_822[k] = f_6 * lsf0_393[k]
                   - f_7 * lsf1_393[k]
                   + f_3 * pc_x[k] * lsg_588[k];

        t_823[k] = f_6 * lsf0_394[k]
                   - f_7 * lsf1_394[k]
                   + f_3 * pc_x[k] * lsg_589[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, lsf0_395, lsf0_396, lsf0_397, lsf1_395, \
                         lsf1_396, lsf1_397, lsg_590, lsg_591, \
                         lsg_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_6 * lsf0_395[k]
                   - f_7 * lsf1_395[k]
                   + f_3 * pc_x[k] * lsg_590[k];

        t_825[k] = f_4 * lsf0_396[k]
                   - f_5 * lsf1_396[k]
                   + f_3 * pc_x[k] * lsg_591[k];

        t_826[k] = f_4 * lsf0_397[k]
                   - f_5 * lsf1_397[k]
                   + f_3 * pc_x[k] * lsg_592[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, t_830, t_831, pc_x, lsf0_398, lsf0_399, \
                         lsf1_398, lsf1_399, lsg_593, lsg_594, lsg_595, lsg_596, \
                         lsg_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_4 * lsf0_398[k]
                   - f_5 * lsf1_398[k]
                   + f_3 * pc_x[k] * lsg_593[k];

        t_828[k] = f_4 * lsf0_399[k]
                   - f_5 * lsf1_399[k]
                   + f_3 * pc_x[k] * lsg_594[k];

        t_829[k] = f_3 * pc_x[k] * lsg_595[k];

        t_830[k] = f_3 * pc_x[k] * lsg_596[k];

        t_831[k] = f_3 * pc_x[k] * lsg_597[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pc_x, pc_y, pc_z, ksg_460, ksg_475, \
                         lsf0_396, lsf1_396, lsg_595, lsg_598, \
                         lsg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_3 * pc_x[k] * lsg_598[k];

        t_833[k] = f_3 * pc_x[k] * lsg_599[k];

        t_834[k] = f_16 * ksg_475[k]
                   + f_1 * lsf0_396[k]
                   - f_2 * lsf1_396[k]
                   + f_3 * pc_y[k] * lsg_595[k];

        t_835[k] = f_11 * ksg_460[k]
                   + f_3 * pc_z[k] * lsg_595[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, pc_y, ksg_477, ksg_478, ksg_479, lsf0_398, \
                         lsf0_399, lsf1_398, lsf1_399, lsg_597, lsg_598, \
                         lsg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * ksg_477[k]
                   + f_6 * lsf0_398[k]
                   - f_7 * lsf1_398[k]
                   + f_3 * pc_y[k] * lsg_597[k];

        t_837[k] = f_16 * ksg_478[k]
                   + f_4 * lsf0_399[k]
                   - f_5 * lsf1_399[k]
                   + f_3 * pc_y[k] * lsg_598[k];

        t_838[k] = f_16 * ksg_479[k]
                   + f_3 * pc_y[k] * lsg_599[k];
    }

#pragma omp simd aligned(t_839, t_840, t_841, pc_x, pc_z, ksg_464, lsf0_399, lsf0_400, \
                         lsf0_401, lsf1_399, lsf1_400, lsf1_401, lsg_599, lsg_600, \
                         lsg_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_839[k] = f_11 * ksg_464[k]
                   + f_1 * lsf0_399[k]
                   - f_2 * lsf1_399[k]
                   + f_3 * pc_z[k] * lsg_599[k];

        t_840[k] = f_1 * lsf0_400[k]
                   - f_2 * lsf1_400[k]
                   + f_3 * pc_x[k] * lsg_600[k];

        t_841[k] = f_13 * lsf0_401[k]
                   - f_14 * lsf1_401[k]
                   + f_3 * pc_x[k] * lsg_601[k];
    }

#pragma omp simd aligned(t_842, t_843, t_844, pc_x, lsf0_402, lsf0_403, lsf0_404, lsf1_402, \
                         lsf1_403, lsf1_404, lsg_602, lsg_603, \
                         lsg_604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_842[k] = f_13 * lsf0_402[k]
                   - f_14 * lsf1_402[k]
                   + f_3 * pc_x[k] * lsg_602[k];

        t_843[k] = f_6 * lsf0_403[k]
                   - f_7 * lsf1_403[k]
                   + f_3 * pc_x[k] * lsg_603[k];

        t_844[k] = f_6 * lsf0_404[k]
                   - f_7 * lsf1_404[k]
                   + f_3 * pc_x[k] * lsg_604[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, lsf0_405, lsf0_406, lsf0_407, lsf1_405, \
                         lsf1_406, lsf1_407, lsg_605, lsg_606, \
                         lsg_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_6 * lsf0_405[k]
                   - f_7 * lsf1_405[k]
                   + f_3 * pc_x[k] * lsg_605[k];

        t_846[k] = f_4 * lsf0_406[k]
                   - f_5 * lsf1_406[k]
                   + f_3 * pc_x[k] * lsg_606[k];

        t_847[k] = f_4 * lsf0_407[k]
                   - f_5 * lsf1_407[k]
                   + f_3 * pc_x[k] * lsg_607[k];
    }
}

static auto
compute_prim_lsh_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksh0,
                                                          const size_t ksg, const size_t ksh1,
                                                          const size_t lsf0, const size_t lsf1,
                                                          const size_t lsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 3.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;
    const auto f_17 = 2.0 / q;

    auto *t_848 = buffer.data(target + 848);
    auto *t_849 = buffer.data(target + 849);
    auto *t_850 = buffer.data(target + 850);
    auto *t_851 = buffer.data(target + 851);
    auto *t_852 = buffer.data(target + 852);
    auto *t_853 = buffer.data(target + 853);
    auto *t_854 = buffer.data(target + 854);
    auto *t_855 = buffer.data(target + 855);
    auto *t_856 = buffer.data(target + 856);
    auto *t_857 = buffer.data(target + 857);
    auto *t_858 = buffer.data(target + 858);
    auto *t_859 = buffer.data(target + 859);
    auto *t_860 = buffer.data(target + 860);
    auto *t_861 = buffer.data(target + 861);
    auto *t_862 = buffer.data(target + 862);
    auto *t_863 = buffer.data(target + 863);
    auto *t_864 = buffer.data(target + 864);
    auto *t_865 = buffer.data(target + 865);
    auto *t_866 = buffer.data(target + 866);
    auto *t_867 = buffer.data(target + 867);
    auto *t_868 = buffer.data(target + 868);
    auto *t_869 = buffer.data(target + 869);
    auto *t_870 = buffer.data(target + 870);
    auto *t_871 = buffer.data(target + 871);
    auto *t_872 = buffer.data(target + 872);
    auto *t_873 = buffer.data(target + 873);
    auto *t_874 = buffer.data(target + 874);
    auto *t_875 = buffer.data(target + 875);
    auto *t_876 = buffer.data(target + 876);
    auto *t_877 = buffer.data(target + 877);
    auto *t_878 = buffer.data(target + 878);
    auto *t_879 = buffer.data(target + 879);
    auto *t_880 = buffer.data(target + 880);
    auto *t_881 = buffer.data(target + 881);
    auto *t_882 = buffer.data(target + 882);
    auto *t_883 = buffer.data(target + 883);
    auto *t_884 = buffer.data(target + 884);
    auto *t_885 = buffer.data(target + 885);
    auto *t_886 = buffer.data(target + 886);
    auto *t_887 = buffer.data(target + 887);
    auto *t_888 = buffer.data(target + 888);
    auto *t_889 = buffer.data(target + 889);
    auto *t_890 = buffer.data(target + 890);
    auto *t_891 = buffer.data(target + 891);
    auto *t_892 = buffer.data(target + 892);
    auto *t_893 = buffer.data(target + 893);
    auto *t_894 = buffer.data(target + 894);
    auto *t_895 = buffer.data(target + 895);
    auto *t_896 = buffer.data(target + 896);
    auto *t_897 = buffer.data(target + 897);
    auto *t_898 = buffer.data(target + 898);
    auto *t_899 = buffer.data(target + 899);
    auto *t_900 = buffer.data(target + 900);
    auto *t_901 = buffer.data(target + 901);
    auto *t_902 = buffer.data(target + 902);
    auto *t_903 = buffer.data(target + 903);
    auto *t_904 = buffer.data(target + 904);
    auto *t_905 = buffer.data(target + 905);
    auto *t_906 = buffer.data(target + 906);
    auto *t_907 = buffer.data(target + 907);
    auto *t_908 = buffer.data(target + 908);
    auto *t_909 = buffer.data(target + 909);
    auto *t_910 = buffer.data(target + 910);
    auto *t_911 = buffer.data(target + 911);
    auto *t_912 = buffer.data(target + 912);
    auto *t_913 = buffer.data(target + 913);
    auto *t_914 = buffer.data(target + 914);
    auto *t_915 = buffer.data(target + 915);
    auto *t_916 = buffer.data(target + 916);
    auto *t_917 = buffer.data(target + 917);
    auto *t_918 = buffer.data(target + 918);
    auto *t_919 = buffer.data(target + 919);
    auto *t_920 = buffer.data(target + 920);
    auto *t_921 = buffer.data(target + 921);
    auto *t_922 = buffer.data(target + 922);
    auto *t_923 = buffer.data(target + 923);
    auto *t_924 = buffer.data(target + 924);
    auto *t_925 = buffer.data(target + 925);
    auto *t_926 = buffer.data(target + 926);
    auto *t_927 = buffer.data(target + 927);
    auto *t_928 = buffer.data(target + 928);
    auto *t_929 = buffer.data(target + 929);
    auto *t_930 = buffer.data(target + 930);
    auto *t_931 = buffer.data(target + 931);
    auto *t_932 = buffer.data(target + 932);
    auto *t_933 = buffer.data(target + 933);
    auto *t_934 = buffer.data(target + 934);
    auto *t_935 = buffer.data(target + 935);
    auto *t_936 = buffer.data(target + 936);
    auto *t_937 = buffer.data(target + 937);
    auto *t_938 = buffer.data(target + 938);
    auto *t_939 = buffer.data(target + 939);
    auto *t_940 = buffer.data(target + 940);
    auto *t_941 = buffer.data(target + 941);
    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh0_735 = buffer.data(ksh0 + 735);
    const auto *ksh0_737 = buffer.data(ksh0 + 737);
    const auto *ksh0_740 = buffer.data(ksh0 + 740);
    const auto *ksh0_744 = buffer.data(ksh0 + 744);
    const auto *ksh0_750 = buffer.data(ksh0 + 750);
    const auto *ksh0_752 = buffer.data(ksh0 + 752);
    const auto *ksh0_753 = buffer.data(ksh0 + 753);
    const auto *ksh0_755 = buffer.data(ksh0 + 755);

    const auto *ksg_475 = buffer.data(ksg + 475);
    const auto *ksg_479 = buffer.data(ksg + 479);
    const auto *ksg_490 = buffer.data(ksg + 490);
    const auto *ksg_492 = buffer.data(ksg + 492);
    const auto *ksg_493 = buffer.data(ksg + 493);
    const auto *ksg_494 = buffer.data(ksg + 494);
    const auto *ksg_505 = buffer.data(ksg + 505);
    const auto *ksg_507 = buffer.data(ksg + 507);
    const auto *ksg_508 = buffer.data(ksg + 508);
    const auto *ksg_509 = buffer.data(ksg + 509);
    const auto *ksg_520 = buffer.data(ksg + 520);
    const auto *ksg_522 = buffer.data(ksg + 522);
    const auto *ksg_523 = buffer.data(ksg + 523);
    const auto *ksg_524 = buffer.data(ksg + 524);
    const auto *ksg_535 = buffer.data(ksg + 535);
    const auto *ksg_537 = buffer.data(ksg + 537);
    const auto *ksg_538 = buffer.data(ksg + 538);
    const auto *ksg_539 = buffer.data(ksg + 539);

    const auto *ksh1_735 = buffer.data(ksh1 + 735);
    const auto *ksh1_737 = buffer.data(ksh1 + 737);
    const auto *ksh1_740 = buffer.data(ksh1 + 740);
    const auto *ksh1_744 = buffer.data(ksh1 + 744);
    const auto *ksh1_750 = buffer.data(ksh1 + 750);
    const auto *ksh1_752 = buffer.data(ksh1 + 752);
    const auto *ksh1_753 = buffer.data(ksh1 + 753);
    const auto *ksh1_755 = buffer.data(ksh1 + 755);

    const auto *lsf0_406 = buffer.data(lsf0 + 406);
    const auto *lsf0_408 = buffer.data(lsf0 + 408);
    const auto *lsf0_409 = buffer.data(lsf0 + 409);
    const auto *lsf0_410 = buffer.data(lsf0 + 410);
    const auto *lsf0_411 = buffer.data(lsf0 + 411);
    const auto *lsf0_412 = buffer.data(lsf0 + 412);
    const auto *lsf0_413 = buffer.data(lsf0 + 413);
    const auto *lsf0_414 = buffer.data(lsf0 + 414);
    const auto *lsf0_415 = buffer.data(lsf0 + 415);
    const auto *lsf0_416 = buffer.data(lsf0 + 416);
    const auto *lsf0_417 = buffer.data(lsf0 + 417);
    const auto *lsf0_418 = buffer.data(lsf0 + 418);
    const auto *lsf0_419 = buffer.data(lsf0 + 419);
    const auto *lsf0_420 = buffer.data(lsf0 + 420);
    const auto *lsf0_421 = buffer.data(lsf0 + 421);
    const auto *lsf0_422 = buffer.data(lsf0 + 422);
    const auto *lsf0_423 = buffer.data(lsf0 + 423);
    const auto *lsf0_424 = buffer.data(lsf0 + 424);
    const auto *lsf0_425 = buffer.data(lsf0 + 425);
    const auto *lsf0_426 = buffer.data(lsf0 + 426);
    const auto *lsf0_427 = buffer.data(lsf0 + 427);
    const auto *lsf0_428 = buffer.data(lsf0 + 428);
    const auto *lsf0_429 = buffer.data(lsf0 + 429);
    const auto *lsf0_431 = buffer.data(lsf0 + 431);
    const auto *lsf0_433 = buffer.data(lsf0 + 433);
    const auto *lsf0_434 = buffer.data(lsf0 + 434);
    const auto *lsf0_436 = buffer.data(lsf0 + 436);
    const auto *lsf0_437 = buffer.data(lsf0 + 437);
    const auto *lsf0_438 = buffer.data(lsf0 + 438);
    const auto *lsf0_440 = buffer.data(lsf0 + 440);
    const auto *lsf0_442 = buffer.data(lsf0 + 442);
    const auto *lsf0_443 = buffer.data(lsf0 + 443);
    const auto *lsf0_445 = buffer.data(lsf0 + 445);
    const auto *lsf0_446 = buffer.data(lsf0 + 446);
    const auto *lsf0_447 = buffer.data(lsf0 + 447);
    const auto *lsf0_448 = buffer.data(lsf0 + 448);
    const auto *lsf0_449 = buffer.data(lsf0 + 449);

    const auto *lsf1_406 = buffer.data(lsf1 + 406);
    const auto *lsf1_408 = buffer.data(lsf1 + 408);
    const auto *lsf1_409 = buffer.data(lsf1 + 409);
    const auto *lsf1_410 = buffer.data(lsf1 + 410);
    const auto *lsf1_411 = buffer.data(lsf1 + 411);
    const auto *lsf1_412 = buffer.data(lsf1 + 412);
    const auto *lsf1_413 = buffer.data(lsf1 + 413);
    const auto *lsf1_414 = buffer.data(lsf1 + 414);
    const auto *lsf1_415 = buffer.data(lsf1 + 415);
    const auto *lsf1_416 = buffer.data(lsf1 + 416);
    const auto *lsf1_417 = buffer.data(lsf1 + 417);
    const auto *lsf1_418 = buffer.data(lsf1 + 418);
    const auto *lsf1_419 = buffer.data(lsf1 + 419);
    const auto *lsf1_420 = buffer.data(lsf1 + 420);
    const auto *lsf1_421 = buffer.data(lsf1 + 421);
    const auto *lsf1_422 = buffer.data(lsf1 + 422);
    const auto *lsf1_423 = buffer.data(lsf1 + 423);
    const auto *lsf1_424 = buffer.data(lsf1 + 424);
    const auto *lsf1_425 = buffer.data(lsf1 + 425);
    const auto *lsf1_426 = buffer.data(lsf1 + 426);
    const auto *lsf1_427 = buffer.data(lsf1 + 427);
    const auto *lsf1_428 = buffer.data(lsf1 + 428);
    const auto *lsf1_429 = buffer.data(lsf1 + 429);
    const auto *lsf1_431 = buffer.data(lsf1 + 431);
    const auto *lsf1_433 = buffer.data(lsf1 + 433);
    const auto *lsf1_434 = buffer.data(lsf1 + 434);
    const auto *lsf1_436 = buffer.data(lsf1 + 436);
    const auto *lsf1_437 = buffer.data(lsf1 + 437);
    const auto *lsf1_438 = buffer.data(lsf1 + 438);
    const auto *lsf1_440 = buffer.data(lsf1 + 440);
    const auto *lsf1_442 = buffer.data(lsf1 + 442);
    const auto *lsf1_443 = buffer.data(lsf1 + 443);
    const auto *lsf1_445 = buffer.data(lsf1 + 445);
    const auto *lsf1_446 = buffer.data(lsf1 + 446);
    const auto *lsf1_447 = buffer.data(lsf1 + 447);
    const auto *lsf1_448 = buffer.data(lsf1 + 448);
    const auto *lsf1_449 = buffer.data(lsf1 + 449);

    const auto *lsg_608 = buffer.data(lsg + 608);
    const auto *lsg_609 = buffer.data(lsg + 609);
    const auto *lsg_610 = buffer.data(lsg + 610);
    const auto *lsg_611 = buffer.data(lsg + 611);
    const auto *lsg_612 = buffer.data(lsg + 612);
    const auto *lsg_613 = buffer.data(lsg + 613);
    const auto *lsg_614 = buffer.data(lsg + 614);
    const auto *lsg_615 = buffer.data(lsg + 615);
    const auto *lsg_616 = buffer.data(lsg + 616);
    const auto *lsg_617 = buffer.data(lsg + 617);
    const auto *lsg_618 = buffer.data(lsg + 618);
    const auto *lsg_619 = buffer.data(lsg + 619);
    const auto *lsg_620 = buffer.data(lsg + 620);
    const auto *lsg_621 = buffer.data(lsg + 621);
    const auto *lsg_622 = buffer.data(lsg + 622);
    const auto *lsg_623 = buffer.data(lsg + 623);
    const auto *lsg_624 = buffer.data(lsg + 624);
    const auto *lsg_625 = buffer.data(lsg + 625);
    const auto *lsg_626 = buffer.data(lsg + 626);
    const auto *lsg_627 = buffer.data(lsg + 627);
    const auto *lsg_628 = buffer.data(lsg + 628);
    const auto *lsg_629 = buffer.data(lsg + 629);
    const auto *lsg_630 = buffer.data(lsg + 630);
    const auto *lsg_631 = buffer.data(lsg + 631);
    const auto *lsg_632 = buffer.data(lsg + 632);
    const auto *lsg_633 = buffer.data(lsg + 633);
    const auto *lsg_634 = buffer.data(lsg + 634);
    const auto *lsg_635 = buffer.data(lsg + 635);
    const auto *lsg_636 = buffer.data(lsg + 636);
    const auto *lsg_637 = buffer.data(lsg + 637);
    const auto *lsg_638 = buffer.data(lsg + 638);
    const auto *lsg_639 = buffer.data(lsg + 639);
    const auto *lsg_640 = buffer.data(lsg + 640);
    const auto *lsg_641 = buffer.data(lsg + 641);
    const auto *lsg_642 = buffer.data(lsg + 642);
    const auto *lsg_643 = buffer.data(lsg + 643);
    const auto *lsg_644 = buffer.data(lsg + 644);
    const auto *lsg_646 = buffer.data(lsg + 646);
    const auto *lsg_648 = buffer.data(lsg + 648);
    const auto *lsg_649 = buffer.data(lsg + 649);
    const auto *lsg_651 = buffer.data(lsg + 651);
    const auto *lsg_652 = buffer.data(lsg + 652);
    const auto *lsg_653 = buffer.data(lsg + 653);
    const auto *lsg_655 = buffer.data(lsg + 655);
    const auto *lsg_656 = buffer.data(lsg + 656);
    const auto *lsg_657 = buffer.data(lsg + 657);
    const auto *lsg_658 = buffer.data(lsg + 658);
    const auto *lsg_659 = buffer.data(lsg + 659);
    const auto *lsg_660 = buffer.data(lsg + 660);
    const auto *lsg_662 = buffer.data(lsg + 662);
    const auto *lsg_663 = buffer.data(lsg + 663);
    const auto *lsg_665 = buffer.data(lsg + 665);
    const auto *lsg_666 = buffer.data(lsg + 666);
    const auto *lsg_667 = buffer.data(lsg + 667);
    const auto *lsg_669 = buffer.data(lsg + 669);
    const auto *lsg_670 = buffer.data(lsg + 670);
    const auto *lsg_671 = buffer.data(lsg + 671);
    const auto *lsg_672 = buffer.data(lsg + 672);
    const auto *lsg_673 = buffer.data(lsg + 673);
    const auto *lsg_674 = buffer.data(lsg + 674);

#pragma omp simd aligned(t_848, t_849, t_850, t_851, t_852, pc_x, lsf0_408, lsf0_409, \
                         lsf1_408, lsf1_409, lsg_608, lsg_609, lsg_610, lsg_611, \
                         lsg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_4 * lsf0_408[k]
                   - f_5 * lsf1_408[k]
                   + f_3 * pc_x[k] * lsg_608[k];

        t_849[k] = f_4 * lsf0_409[k]
                   - f_5 * lsf1_409[k]
                   + f_3 * pc_x[k] * lsg_609[k];

        t_850[k] = f_3 * pc_x[k] * lsg_610[k];

        t_851[k] = f_3 * pc_x[k] * lsg_611[k];

        t_852[k] = f_3 * pc_x[k] * lsg_612[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pc_x, pc_y, pc_z, ksg_475, ksg_490, \
                         lsf0_406, lsf1_406, lsg_610, lsg_613, \
                         lsg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_3 * pc_x[k] * lsg_613[k];

        t_854[k] = f_3 * pc_x[k] * lsg_614[k];

        t_855[k] = f_17 * ksg_490[k]
                   + f_1 * lsf0_406[k]
                   - f_2 * lsf1_406[k]
                   + f_3 * pc_y[k] * lsg_610[k];

        t_856[k] = f_17 * ksg_475[k]
                   + f_3 * pc_z[k] * lsg_610[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pc_y, ksg_492, ksg_493, ksg_494, lsf0_408, \
                         lsf0_409, lsf1_408, lsf1_409, lsg_612, lsg_613, \
                         lsg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_17 * ksg_492[k]
                   + f_6 * lsf0_408[k]
                   - f_7 * lsf1_408[k]
                   + f_3 * pc_y[k] * lsg_612[k];

        t_858[k] = f_17 * ksg_493[k]
                   + f_4 * lsf0_409[k]
                   - f_5 * lsf1_409[k]
                   + f_3 * pc_y[k] * lsg_613[k];

        t_859[k] = f_17 * ksg_494[k]
                   + f_3 * pc_y[k] * lsg_614[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pc_x, pc_z, ksg_479, lsf0_409, lsf0_410, \
                         lsf0_411, lsf1_409, lsf1_410, lsf1_411, lsg_614, lsg_615, \
                         lsg_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = f_17 * ksg_479[k]
                   + f_1 * lsf0_409[k]
                   - f_2 * lsf1_409[k]
                   + f_3 * pc_z[k] * lsg_614[k];

        t_861[k] = f_1 * lsf0_410[k]
                   - f_2 * lsf1_410[k]
                   + f_3 * pc_x[k] * lsg_615[k];

        t_862[k] = f_13 * lsf0_411[k]
                   - f_14 * lsf1_411[k]
                   + f_3 * pc_x[k] * lsg_616[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_x, lsf0_412, lsf0_413, lsf0_414, lsf1_412, \
                         lsf1_413, lsf1_414, lsg_617, lsg_618, \
                         lsg_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_13 * lsf0_412[k]
                   - f_14 * lsf1_412[k]
                   + f_3 * pc_x[k] * lsg_617[k];

        t_864[k] = f_6 * lsf0_413[k]
                   - f_7 * lsf1_413[k]
                   + f_3 * pc_x[k] * lsg_618[k];

        t_865[k] = f_6 * lsf0_414[k]
                   - f_7 * lsf1_414[k]
                   + f_3 * pc_x[k] * lsg_619[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, lsf0_415, lsf0_416, lsf0_417, lsf1_415, \
                         lsf1_416, lsf1_417, lsg_620, lsg_621, \
                         lsg_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_6 * lsf0_415[k]
                   - f_7 * lsf1_415[k]
                   + f_3 * pc_x[k] * lsg_620[k];

        t_867[k] = f_4 * lsf0_416[k]
                   - f_5 * lsf1_416[k]
                   + f_3 * pc_x[k] * lsg_621[k];

        t_868[k] = f_4 * lsf0_417[k]
                   - f_5 * lsf1_417[k]
                   + f_3 * pc_x[k] * lsg_622[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, t_872, t_873, pc_x, lsf0_418, lsf0_419, \
                         lsf1_418, lsf1_419, lsg_623, lsg_624, lsg_625, lsg_626, \
                         lsg_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_4 * lsf0_418[k]
                   - f_5 * lsf1_418[k]
                   + f_3 * pc_x[k] * lsg_623[k];

        t_870[k] = f_4 * lsf0_419[k]
                   - f_5 * lsf1_419[k]
                   + f_3 * pc_x[k] * lsg_624[k];

        t_871[k] = f_3 * pc_x[k] * lsg_625[k];

        t_872[k] = f_3 * pc_x[k] * lsg_626[k];

        t_873[k] = f_3 * pc_x[k] * lsg_627[k];
    }

#pragma omp simd aligned(t_874, t_875, t_876, t_877, pc_x, pc_y, pc_z, ksg_490, ksg_505, \
                         lsf0_416, lsf1_416, lsg_625, lsg_628, \
                         lsg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_874[k] = f_3 * pc_x[k] * lsg_628[k];

        t_875[k] = f_3 * pc_x[k] * lsg_629[k];

        t_876[k] = f_11 * ksg_505[k]
                   + f_1 * lsf0_416[k]
                   - f_2 * lsf1_416[k]
                   + f_3 * pc_y[k] * lsg_625[k];

        t_877[k] = f_16 * ksg_490[k]
                   + f_3 * pc_z[k] * lsg_625[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_y, ksg_507, ksg_508, ksg_509, lsf0_418, \
                         lsf0_419, lsf1_418, lsf1_419, lsg_627, lsg_628, \
                         lsg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_11 * ksg_507[k]
                   + f_6 * lsf0_418[k]
                   - f_7 * lsf1_418[k]
                   + f_3 * pc_y[k] * lsg_627[k];

        t_879[k] = f_11 * ksg_508[k]
                   + f_4 * lsf0_419[k]
                   - f_5 * lsf1_419[k]
                   + f_3 * pc_y[k] * lsg_628[k];

        t_880[k] = f_11 * ksg_509[k]
                   + f_3 * pc_y[k] * lsg_629[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, pc_x, pc_z, ksg_494, lsf0_419, lsf0_420, \
                         lsf0_421, lsf1_419, lsf1_420, lsf1_421, lsg_629, lsg_630, \
                         lsg_631 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_16 * ksg_494[k]
                   + f_1 * lsf0_419[k]
                   - f_2 * lsf1_419[k]
                   + f_3 * pc_z[k] * lsg_629[k];

        t_882[k] = f_1 * lsf0_420[k]
                   - f_2 * lsf1_420[k]
                   + f_3 * pc_x[k] * lsg_630[k];

        t_883[k] = f_13 * lsf0_421[k]
                   - f_14 * lsf1_421[k]
                   + f_3 * pc_x[k] * lsg_631[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, pc_x, lsf0_422, lsf0_423, lsf0_424, lsf1_422, \
                         lsf1_423, lsf1_424, lsg_632, lsg_633, \
                         lsg_634 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_13 * lsf0_422[k]
                   - f_14 * lsf1_422[k]
                   + f_3 * pc_x[k] * lsg_632[k];

        t_885[k] = f_6 * lsf0_423[k]
                   - f_7 * lsf1_423[k]
                   + f_3 * pc_x[k] * lsg_633[k];

        t_886[k] = f_6 * lsf0_424[k]
                   - f_7 * lsf1_424[k]
                   + f_3 * pc_x[k] * lsg_634[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pc_x, lsf0_425, lsf0_426, lsf0_427, lsf1_425, \
                         lsf1_426, lsf1_427, lsg_635, lsg_636, \
                         lsg_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_6 * lsf0_425[k]
                   - f_7 * lsf1_425[k]
                   + f_3 * pc_x[k] * lsg_635[k];

        t_888[k] = f_4 * lsf0_426[k]
                   - f_5 * lsf1_426[k]
                   + f_3 * pc_x[k] * lsg_636[k];

        t_889[k] = f_4 * lsf0_427[k]
                   - f_5 * lsf1_427[k]
                   + f_3 * pc_x[k] * lsg_637[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, t_894, pc_x, lsf0_428, lsf0_429, \
                         lsf1_428, lsf1_429, lsg_638, lsg_639, lsg_640, lsg_641, \
                         lsg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_4 * lsf0_428[k]
                   - f_5 * lsf1_428[k]
                   + f_3 * pc_x[k] * lsg_638[k];

        t_891[k] = f_4 * lsf0_429[k]
                   - f_5 * lsf1_429[k]
                   + f_3 * pc_x[k] * lsg_639[k];

        t_892[k] = f_3 * pc_x[k] * lsg_640[k];

        t_893[k] = f_3 * pc_x[k] * lsg_641[k];

        t_894[k] = f_3 * pc_x[k] * lsg_642[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pc_x, pc_y, pc_z, ksg_505, ksg_520, \
                         lsf0_426, lsf1_426, lsg_640, lsg_643, \
                         lsg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_3 * pc_x[k] * lsg_643[k];

        t_896[k] = f_3 * pc_x[k] * lsg_644[k];

        t_897[k] = f_10 * ksg_520[k]
                   + f_1 * lsf0_426[k]
                   - f_2 * lsf1_426[k]
                   + f_3 * pc_y[k] * lsg_640[k];

        t_898[k] = f_15 * ksg_505[k]
                   + f_3 * pc_z[k] * lsg_640[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, pc_y, ksg_522, ksg_523, ksg_524, lsf0_428, \
                         lsf0_429, lsf1_428, lsf1_429, lsg_642, lsg_643, \
                         lsg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = f_10 * ksg_522[k]
                   + f_6 * lsf0_428[k]
                   - f_7 * lsf1_428[k]
                   + f_3 * pc_y[k] * lsg_642[k];

        t_900[k] = f_10 * ksg_523[k]
                   + f_4 * lsf0_429[k]
                   - f_5 * lsf1_429[k]
                   + f_3 * pc_y[k] * lsg_643[k];

        t_901[k] = f_10 * ksg_524[k]
                   + f_3 * pc_y[k] * lsg_644[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pa_y, pc_x, pc_y, pc_z, ksh0_735, ksg_509, \
                         ksh1_735, lsf0_429, lsf0_431, lsf1_429, lsf1_431, lsg_644, \
                         lsg_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_15 * ksg_509[k]
                   + f_1 * lsf0_429[k]
                   - f_2 * lsf1_429[k]
                   + f_3 * pc_z[k] * lsg_644[k];

        t_903[k] = pa_y[k] * ksh0_735[k]
                   - f_8 * pc_y[k] * ksh1_735[k];

        t_904[k] = f_13 * lsf0_431[k]
                   - f_14 * lsf1_431[k]
                   + f_3 * pc_x[k] * lsg_646[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pa_y, pc_x, pc_y, ksh0_737, ksh1_737, lsf0_433, \
                         lsf0_434, lsf1_433, lsf1_434, lsg_648, \
                         lsg_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = pa_y[k] * ksh0_737[k]
                   - f_8 * pc_y[k] * ksh1_737[k];

        t_906[k] = f_6 * lsf0_433[k]
                   - f_7 * lsf1_433[k]
                   + f_3 * pc_x[k] * lsg_648[k];

        t_907[k] = f_6 * lsf0_434[k]
                   - f_7 * lsf1_434[k]
                   + f_3 * pc_x[k] * lsg_649[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pa_y, pc_x, pc_y, ksh0_740, ksh1_740, lsf0_436, \
                         lsf0_437, lsf1_436, lsf1_437, lsg_651, \
                         lsg_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = pa_y[k] * ksh0_740[k]
                   - f_8 * pc_y[k] * ksh1_740[k];

        t_909[k] = f_4 * lsf0_436[k]
                   - f_5 * lsf1_436[k]
                   + f_3 * pc_x[k] * lsg_651[k];

        t_910[k] = f_4 * lsf0_437[k]
                   - f_5 * lsf1_437[k]
                   + f_3 * pc_x[k] * lsg_652[k];
    }

#pragma omp simd aligned(t_911, t_912, t_913, t_914, t_915, pa_y, pc_x, pc_y, ksh0_744, \
                         ksh1_744, lsf0_438, lsf1_438, lsg_653, lsg_655, lsg_656, \
                         lsg_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_4 * lsf0_438[k]
                   - f_5 * lsf1_438[k]
                   + f_3 * pc_x[k] * lsg_653[k];

        t_912[k] = pa_y[k] * ksh0_744[k]
                   - f_8 * pc_y[k] * ksh1_744[k];

        t_913[k] = f_3 * pc_x[k] * lsg_655[k];

        t_914[k] = f_3 * pc_x[k] * lsg_656[k];

        t_915[k] = f_3 * pc_x[k] * lsg_657[k];
    }

#pragma omp simd aligned(t_916, t_917, t_918, t_919, pa_y, pc_x, pc_y, pc_z, ksh0_750, \
                         ksg_520, ksg_535, ksh1_750, lsg_655, lsg_658, \
                         lsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_916[k] = f_3 * pc_x[k] * lsg_658[k];

        t_917[k] = f_3 * pc_x[k] * lsg_659[k];

        t_918[k] = pa_y[k] * ksh0_750[k]
                   + f_16 * ksg_535[k]
                   - f_8 * pc_y[k] * ksh1_750[k];

        t_919[k] = f_12 * ksg_520[k]
                   + f_3 * pc_z[k] * lsg_655[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pa_y, pc_y, ksh0_752, ksh0_753, ksh0_755, \
                         ksg_537, ksg_538, ksg_539, ksh1_752, ksh1_753, ksh1_755, \
                         lsg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pa_y[k] * ksh0_752[k]
                   + f_11 * ksg_537[k]
                   - f_8 * pc_y[k] * ksh1_752[k];

        t_921[k] = pa_y[k] * ksh0_753[k]
                   + f_10 * ksg_538[k]
                   - f_8 * pc_y[k] * ksh1_753[k];

        t_922[k] = f_9 * ksg_539[k]
                   + f_3 * pc_y[k] * lsg_659[k];

        t_923[k] = pa_y[k] * ksh0_755[k]
                   - f_8 * pc_y[k] * ksh1_755[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, t_928, pc_x, pc_y, lsf0_440, lsf0_442, \
                         lsf0_443, lsf1_440, lsf1_442, lsf1_443, lsg_660, lsg_662, \
                         lsg_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_1 * lsf0_440[k]
                   - f_2 * lsf1_440[k]
                   + f_3 * pc_x[k] * lsg_660[k];

        t_925[k] = f_3 * pc_y[k] * lsg_660[k];

        t_926[k] = f_13 * lsf0_442[k]
                   - f_14 * lsf1_442[k]
                   + f_3 * pc_x[k] * lsg_662[k];

        t_927[k] = f_6 * lsf0_443[k]
                   - f_7 * lsf1_443[k]
                   + f_3 * pc_x[k] * lsg_663[k];

        t_928[k] = f_3 * pc_y[k] * lsg_662[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, t_932, pc_x, pc_y, lsf0_445, lsf0_446, lsf0_447, \
                         lsf1_445, lsf1_446, lsf1_447, lsg_665, lsg_666, \
                         lsg_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = f_6 * lsf0_445[k]
                   - f_7 * lsf1_445[k]
                   + f_3 * pc_x[k] * lsg_665[k];

        t_930[k] = f_4 * lsf0_446[k]
                   - f_5 * lsf1_446[k]
                   + f_3 * pc_x[k] * lsg_666[k];

        t_931[k] = f_4 * lsf0_447[k]
                   - f_5 * lsf1_447[k]
                   + f_3 * pc_x[k] * lsg_667[k];

        t_932[k] = f_3 * pc_y[k] * lsg_665[k];
    }

#pragma omp simd aligned(t_933, t_934, t_935, t_936, t_937, t_938, pc_x, lsf0_449, lsf1_449, \
                         lsg_669, lsg_670, lsg_671, lsg_672, lsg_673, \
                         lsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_933[k] = f_4 * lsf0_449[k]
                   - f_5 * lsf1_449[k]
                   + f_3 * pc_x[k] * lsg_669[k];

        t_934[k] = f_3 * pc_x[k] * lsg_670[k];

        t_935[k] = f_3 * pc_x[k] * lsg_671[k];

        t_936[k] = f_3 * pc_x[k] * lsg_672[k];

        t_937[k] = f_3 * pc_x[k] * lsg_673[k];

        t_938[k] = f_3 * pc_x[k] * lsg_674[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, pc_y, lsf0_446, lsf0_447, lsf0_448, lsf1_446, \
                         lsf1_447, lsf1_448, lsg_670, lsg_671, \
                         lsg_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_1 * lsf0_446[k]
                   - f_2 * lsf1_446[k]
                   + f_3 * pc_y[k] * lsg_670[k];

        t_940[k] = f_13 * lsf0_447[k]
                   - f_14 * lsf1_447[k]
                   + f_3 * pc_y[k] * lsg_671[k];

        t_941[k] = f_6 * lsf0_448[k]
                   - f_7 * lsf1_448[k]
                   + f_3 * pc_y[k] * lsg_672[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, pc_y, pc_z, ksg_539, lsf0_449, lsf1_449, \
                         lsg_673, lsg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_4 * lsf0_449[k]
                   - f_5 * lsf1_449[k]
                   + f_3 * pc_y[k] * lsg_673[k];

        t_943[k] = f_3 * pc_y[k] * lsg_674[k];

        t_944[k] = f_0 * ksg_539[k]
                   + f_1 * lsf0_449[k]
                   - f_2 * lsf1_449[k]
                   + f_3 * pc_z[k] * lsg_674[k];
    }
}

auto
compute_prim_lsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksh0, const size_t ksg,
                                                   const size_t ksh1, const size_t lsf0,
                                                   const size_t lsf1, const size_t lsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsg, ncols, gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);

    compute_prim_lsh_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, ksh0, ksg,
                                                              ksh1, lsf0, lsf1, lsg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
