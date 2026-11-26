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


#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_12 = 3.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *ish0_0 = buffer.data(ish0 + 0);
    const auto *ish0_3 = buffer.data(ish0 + 3);
    const auto *ish0_5 = buffer.data(ish0 + 5);
    const auto *ish0_6 = buffer.data(ish0 + 6);
    const auto *ish0_9 = buffer.data(ish0 + 9);
    const auto *ish0_15 = buffer.data(ish0 + 15);
    const auto *ish0_20 = buffer.data(ish0 + 20);
    const auto *ish0_24 = buffer.data(ish0 + 24);
    const auto *ish0_27 = buffer.data(ish0 + 27);
    const auto *ish0_36 = buffer.data(ish0 + 36);
    const auto *ish0_42 = buffer.data(ish0 + 42);
    const auto *ish0_47 = buffer.data(ish0 + 47);
    const auto *ish0_51 = buffer.data(ish0 + 51);
    const auto *ish0_62 = buffer.data(ish0 + 62);

    const auto *isg_0 = buffer.data(isg + 0);
    const auto *isg_1 = buffer.data(isg + 1);
    const auto *isg_2 = buffer.data(isg + 2);
    const auto *isg_3 = buffer.data(isg + 3);
    const auto *isg_5 = buffer.data(isg + 5);
    const auto *isg_10 = buffer.data(isg + 10);
    const auto *isg_12 = buffer.data(isg + 12);
    const auto *isg_14 = buffer.data(isg + 14);
    const auto *isg_15 = buffer.data(isg + 15);
    const auto *isg_18 = buffer.data(isg + 18);
    const auto *isg_20 = buffer.data(isg + 20);
    const auto *isg_25 = buffer.data(isg + 25);
    const auto *isg_27 = buffer.data(isg + 27);
    const auto *isg_28 = buffer.data(isg + 28);
    const auto *isg_29 = buffer.data(isg + 29);
    const auto *isg_30 = buffer.data(isg + 30);
    const auto *isg_32 = buffer.data(isg + 32);
    const auto *isg_35 = buffer.data(isg + 35);
    const auto *isg_40 = buffer.data(isg + 40);
    const auto *isg_41 = buffer.data(isg + 41);
    const auto *isg_42 = buffer.data(isg + 42);
    const auto *isg_43 = buffer.data(isg + 43);
    const auto *isg_44 = buffer.data(isg + 44);
    const auto *isg_45 = buffer.data(isg + 45);
    const auto *isg_48 = buffer.data(isg + 48);
    const auto *isg_51 = buffer.data(isg + 51);
    const auto *isg_55 = buffer.data(isg + 55);
    const auto *isg_57 = buffer.data(isg + 57);
    const auto *isg_58 = buffer.data(isg + 58);
    const auto *isg_59 = buffer.data(isg + 59);
    const auto *isg_70 = buffer.data(isg + 70);
    const auto *isg_71 = buffer.data(isg + 71);
    const auto *isg_72 = buffer.data(isg + 72);
    const auto *isg_73 = buffer.data(isg + 73);
    const auto *isg_74 = buffer.data(isg + 74);
    const auto *isg_75 = buffer.data(isg + 75);
    const auto *isg_80 = buffer.data(isg + 80);
    const auto *isg_84 = buffer.data(isg + 84);
    const auto *isg_85 = buffer.data(isg + 85);
    const auto *isg_86 = buffer.data(isg + 86);
    const auto *isg_87 = buffer.data(isg + 87);
    const auto *isg_89 = buffer.data(isg + 89);
    const auto *isg_90 = buffer.data(isg + 90);
    const auto *isg_93 = buffer.data(isg + 93);

    const auto *ish1_0 = buffer.data(ish1 + 0);
    const auto *ish1_3 = buffer.data(ish1 + 3);
    const auto *ish1_5 = buffer.data(ish1 + 5);
    const auto *ish1_6 = buffer.data(ish1 + 6);
    const auto *ish1_9 = buffer.data(ish1 + 9);
    const auto *ish1_15 = buffer.data(ish1 + 15);
    const auto *ish1_20 = buffer.data(ish1 + 20);
    const auto *ish1_24 = buffer.data(ish1 + 24);
    const auto *ish1_27 = buffer.data(ish1 + 27);
    const auto *ish1_36 = buffer.data(ish1 + 36);
    const auto *ish1_42 = buffer.data(ish1 + 42);
    const auto *ish1_47 = buffer.data(ish1 + 47);
    const auto *ish1_51 = buffer.data(ish1 + 51);
    const auto *ish1_62 = buffer.data(ish1 + 62);

    const auto *ksf0_0 = buffer.data(ksf0 + 0);
    const auto *ksf0_1 = buffer.data(ksf0 + 1);
    const auto *ksf0_2 = buffer.data(ksf0 + 2);
    const auto *ksf0_6 = buffer.data(ksf0 + 6);
    const auto *ksf0_8 = buffer.data(ksf0 + 8);
    const auto *ksf0_9 = buffer.data(ksf0 + 9);
    const auto *ksf0_16 = buffer.data(ksf0 + 16);
    const auto *ksf0_17 = buffer.data(ksf0 + 17);
    const auto *ksf0_22 = buffer.data(ksf0 + 22);
    const auto *ksf0_27 = buffer.data(ksf0 + 27);
    const auto *ksf0_28 = buffer.data(ksf0 + 28);
    const auto *ksf0_29 = buffer.data(ksf0 + 29);
    const auto *ksf0_30 = buffer.data(ksf0 + 30);
    const auto *ksf0_32 = buffer.data(ksf0 + 32);
    const auto *ksf0_33 = buffer.data(ksf0 + 33);
    const auto *ksf0_36 = buffer.data(ksf0 + 36);
    const auto *ksf0_37 = buffer.data(ksf0 + 37);
    const auto *ksf0_39 = buffer.data(ksf0 + 39);
    const auto *ksf0_48 = buffer.data(ksf0 + 48);
    const auto *ksf0_49 = buffer.data(ksf0 + 49);
    const auto *ksf0_50 = buffer.data(ksf0 + 50);
    const auto *ksf0_51 = buffer.data(ksf0 + 51);
    const auto *ksf0_52 = buffer.data(ksf0 + 52);
    const auto *ksf0_55 = buffer.data(ksf0 + 55);
    const auto *ksf0_56 = buffer.data(ksf0 + 56);
    const auto *ksf0_57 = buffer.data(ksf0 + 57);
    const auto *ksf0_58 = buffer.data(ksf0 + 58);
    const auto *ksf0_59 = buffer.data(ksf0 + 59);
    const auto *ksf0_60 = buffer.data(ksf0 + 60);
    const auto *ksf0_63 = buffer.data(ksf0 + 63);

    const auto *ksf1_0 = buffer.data(ksf1 + 0);
    const auto *ksf1_1 = buffer.data(ksf1 + 1);
    const auto *ksf1_2 = buffer.data(ksf1 + 2);
    const auto *ksf1_6 = buffer.data(ksf1 + 6);
    const auto *ksf1_8 = buffer.data(ksf1 + 8);
    const auto *ksf1_9 = buffer.data(ksf1 + 9);
    const auto *ksf1_16 = buffer.data(ksf1 + 16);
    const auto *ksf1_17 = buffer.data(ksf1 + 17);
    const auto *ksf1_22 = buffer.data(ksf1 + 22);
    const auto *ksf1_27 = buffer.data(ksf1 + 27);
    const auto *ksf1_28 = buffer.data(ksf1 + 28);
    const auto *ksf1_29 = buffer.data(ksf1 + 29);
    const auto *ksf1_30 = buffer.data(ksf1 + 30);
    const auto *ksf1_32 = buffer.data(ksf1 + 32);
    const auto *ksf1_33 = buffer.data(ksf1 + 33);
    const auto *ksf1_36 = buffer.data(ksf1 + 36);
    const auto *ksf1_37 = buffer.data(ksf1 + 37);
    const auto *ksf1_39 = buffer.data(ksf1 + 39);
    const auto *ksf1_48 = buffer.data(ksf1 + 48);
    const auto *ksf1_49 = buffer.data(ksf1 + 49);
    const auto *ksf1_50 = buffer.data(ksf1 + 50);
    const auto *ksf1_51 = buffer.data(ksf1 + 51);
    const auto *ksf1_52 = buffer.data(ksf1 + 52);
    const auto *ksf1_55 = buffer.data(ksf1 + 55);
    const auto *ksf1_56 = buffer.data(ksf1 + 56);
    const auto *ksf1_57 = buffer.data(ksf1 + 57);
    const auto *ksf1_58 = buffer.data(ksf1 + 58);
    const auto *ksf1_59 = buffer.data(ksf1 + 59);
    const auto *ksf1_60 = buffer.data(ksf1 + 60);
    const auto *ksf1_63 = buffer.data(ksf1 + 63);

    const auto *ksg_0 = buffer.data(ksg + 0);
    const auto *ksg_1 = buffer.data(ksg + 1);
    const auto *ksg_2 = buffer.data(ksg + 2);
    const auto *ksg_3 = buffer.data(ksg + 3);
    const auto *ksg_5 = buffer.data(ksg + 5);
    const auto *ksg_6 = buffer.data(ksg + 6);
    const auto *ksg_9 = buffer.data(ksg + 9);
    const auto *ksg_10 = buffer.data(ksg + 10);
    const auto *ksg_12 = buffer.data(ksg + 12);
    const auto *ksg_13 = buffer.data(ksg + 13);
    const auto *ksg_14 = buffer.data(ksg + 14);
    const auto *ksg_15 = buffer.data(ksg + 15);
    const auto *ksg_16 = buffer.data(ksg + 16);
    const auto *ksg_18 = buffer.data(ksg + 18);
    const auto *ksg_20 = buffer.data(ksg + 20);
    const auto *ksg_21 = buffer.data(ksg + 21);
    const auto *ksg_25 = buffer.data(ksg + 25);
    const auto *ksg_26 = buffer.data(ksg + 26);
    const auto *ksg_27 = buffer.data(ksg + 27);
    const auto *ksg_28 = buffer.data(ksg + 28);
    const auto *ksg_29 = buffer.data(ksg + 29);
    const auto *ksg_30 = buffer.data(ksg + 30);
    const auto *ksg_32 = buffer.data(ksg + 32);
    const auto *ksg_34 = buffer.data(ksg + 34);
    const auto *ksg_35 = buffer.data(ksg + 35);
    const auto *ksg_39 = buffer.data(ksg + 39);
    const auto *ksg_40 = buffer.data(ksg + 40);
    const auto *ksg_41 = buffer.data(ksg + 41);
    const auto *ksg_42 = buffer.data(ksg + 42);
    const auto *ksg_43 = buffer.data(ksg + 43);
    const auto *ksg_44 = buffer.data(ksg + 44);
    const auto *ksg_45 = buffer.data(ksg + 45);
    const auto *ksg_46 = buffer.data(ksg + 46);
    const auto *ksg_47 = buffer.data(ksg + 47);
    const auto *ksg_48 = buffer.data(ksg + 48);
    const auto *ksg_50 = buffer.data(ksg + 50);
    const auto *ksg_51 = buffer.data(ksg + 51);
    const auto *ksg_55 = buffer.data(ksg + 55);
    const auto *ksg_56 = buffer.data(ksg + 56);
    const auto *ksg_57 = buffer.data(ksg + 57);
    const auto *ksg_58 = buffer.data(ksg + 58);
    const auto *ksg_59 = buffer.data(ksg + 59);
    const auto *ksg_60 = buffer.data(ksg + 60);
    const auto *ksg_62 = buffer.data(ksg + 62);
    const auto *ksg_63 = buffer.data(ksg + 63);
    const auto *ksg_65 = buffer.data(ksg + 65);
    const auto *ksg_70 = buffer.data(ksg + 70);
    const auto *ksg_71 = buffer.data(ksg + 71);
    const auto *ksg_72 = buffer.data(ksg + 72);
    const auto *ksg_73 = buffer.data(ksg + 73);
    const auto *ksg_74 = buffer.data(ksg + 74);
    const auto *ksg_75 = buffer.data(ksg + 75);
    const auto *ksg_76 = buffer.data(ksg + 76);
    const auto *ksg_77 = buffer.data(ksg + 77);
    const auto *ksg_78 = buffer.data(ksg + 78);
    const auto *ksg_79 = buffer.data(ksg + 79);
    const auto *ksg_80 = buffer.data(ksg + 80);
    const auto *ksg_84 = buffer.data(ksg + 84);
    const auto *ksg_85 = buffer.data(ksg + 85);
    const auto *ksg_86 = buffer.data(ksg + 86);
    const auto *ksg_87 = buffer.data(ksg + 87);
    const auto *ksg_88 = buffer.data(ksg + 88);
    const auto *ksg_89 = buffer.data(ksg + 89);
    const auto *ksg_90 = buffer.data(ksg + 90);
    const auto *ksg_93 = buffer.data(ksg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, isg_0, ksf0_0, \
                         ksf1_0, ksg_0, ksg_1, ksg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isg_0[k]
                 + f_1 * ksf0_0[k]
                 - f_2 * ksf1_0[k]
                 + f_3 * pc_x[k] * ksg_0[k];

        t_1[k] = f_3 * pc_y[k] * ksg_0[k];

        t_2[k] = f_3 * pc_z[k] * ksg_0[k];

        t_3[k] = f_4 * ksf0_0[k]
                 - f_5 * ksf1_0[k]
                 + f_3 * pc_y[k] * ksg_1[k];

        t_4[k] = f_3 * pc_y[k] * ksg_2[k];

        t_5[k] = f_4 * ksf0_0[k]
                 - f_5 * ksf1_0[k]
                 + f_3 * pc_z[k] * ksg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, isg_10, ksf0_1, ksf0_2, \
                         ksf1_1, ksf1_2, ksg_3, ksg_5, ksg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ksf0_1[k]
                 - f_7 * ksf1_1[k]
                 + f_3 * pc_y[k] * ksg_3[k];

        t_7[k] = f_3 * pc_z[k] * ksg_3[k];

        t_8[k] = f_3 * pc_y[k] * ksg_5[k];

        t_9[k] = f_6 * ksf0_2[k]
                 - f_7 * ksf1_2[k]
                 + f_3 * pc_z[k] * ksg_5[k];

        t_10[k] = f_0 * isg_10[k]
                  + f_3 * pc_x[k] * ksg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, isg_12, isg_14, ksg_6, \
                         ksg_9, ksg_12, ksg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ksg_6[k];

        t_12[k] = f_0 * isg_12[k]
                  + f_3 * pc_x[k] * ksg_12[k];

        t_13[k] = f_3 * pc_y[k] * ksg_9[k];

        t_14[k] = f_0 * isg_14[k]
                  + f_3 * pc_x[k] * ksg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, ksf0_6, ksf0_8, ksf0_9, ksf1_6, \
                         ksf1_8, ksf1_9, ksg_10, ksg_12, ksg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * ksf0_6[k]
                  - f_2 * ksf1_6[k]
                  + f_3 * pc_y[k] * ksg_10[k];

        t_16[k] = f_3 * pc_z[k] * ksg_10[k];

        t_17[k] = f_6 * ksf0_8[k]
                  - f_7 * ksf1_8[k]
                  + f_3 * pc_y[k] * ksg_12[k];

        t_18[k] = f_4 * ksf0_9[k]
                  - f_5 * ksf1_9[k]
                  + f_3 * pc_y[k] * ksg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, ish0_0, isg_0, \
                         ish1_0, ksf0_9, ksf1_9, ksg_14, ksg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * ksg_14[k];

        t_20[k] = f_1 * ksf0_9[k]
                  - f_2 * ksf1_9[k]
                  + f_3 * pc_z[k] * ksg_14[k];

        t_21[k] = pa_y[k] * ish0_0[k]
                  - f_8 * pc_y[k] * ish1_0[k];

        t_22[k] = f_9 * isg_0[k]
                  + f_3 * pc_y[k] * ksg_15[k];

        t_23[k] = f_3 * pc_z[k] * ksg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, ish0_3, ish0_5, ish0_6, \
                         isg_1, isg_3, ish1_3, ish1_5, ish1_6, ksg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * ish0_3[k]
                  + f_10 * isg_1[k]
                  - f_8 * pc_y[k] * ish1_3[k];

        t_25[k] = f_3 * pc_z[k] * ksg_16[k];

        t_26[k] = pa_y[k] * ish0_5[k]
                  - f_8 * pc_y[k] * ish1_5[k];

        t_27[k] = pa_y[k] * ish0_6[k]
                  + f_11 * isg_3[k]
                  - f_8 * pc_y[k] * ish1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, ish0_9, isg_5, \
                         isg_25, ish1_9, ksg_18, ksg_20, ksg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * ksg_18[k];

        t_29[k] = f_9 * isg_5[k]
                  + f_3 * pc_y[k] * ksg_20[k];

        t_30[k] = pa_y[k] * ish0_9[k]
                  - f_8 * pc_y[k] * ish1_9[k];

        t_31[k] = f_12 * isg_25[k]
                  + f_3 * pc_x[k] * ksg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, isg_27, isg_28, isg_29, ksg_21, \
                         ksg_27, ksg_28, ksg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * ksg_21[k];

        t_33[k] = f_12 * isg_27[k]
                  + f_3 * pc_x[k] * ksg_27[k];

        t_34[k] = f_12 * isg_28[k]
                  + f_3 * pc_x[k] * ksg_28[k];

        t_35[k] = f_12 * isg_29[k]
                  + f_3 * pc_x[k] * ksg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, isg_10, ksf0_16, ksf0_17, \
                         ksf1_16, ksf1_17, ksg_25, ksg_26, ksg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * isg_10[k]
                  + f_1 * ksf0_16[k]
                  - f_2 * ksf1_16[k]
                  + f_3 * pc_y[k] * ksg_25[k];

        t_37[k] = f_3 * pc_z[k] * ksg_25[k];

        t_38[k] = f_4 * ksf0_16[k]
                  - f_5 * ksf1_16[k]
                  + f_3 * pc_z[k] * ksg_26[k];

        t_39[k] = f_6 * ksf0_17[k]
                  - f_7 * ksf1_17[k]
                  + f_3 * pc_z[k] * ksg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, ish0_0, ish0_20, \
                         isg_14, ish1_0, ish1_20, ksg_29, ksg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * isg_14[k]
                  + f_3 * pc_y[k] * ksg_29[k];

        t_41[k] = pa_y[k] * ish0_20[k]
                  - f_8 * pc_y[k] * ish1_20[k];

        t_42[k] = pa_z[k] * ish0_0[k]
                  - f_8 * pc_z[k] * ish1_0[k];

        t_43[k] = f_3 * pc_y[k] * ksg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, ish0_3, ish0_5, isg_0, \
                         isg_2, ish1_3, ish1_5, ksg_30, ksg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * isg_0[k]
                  + f_3 * pc_z[k] * ksg_30[k];

        t_45[k] = pa_z[k] * ish0_3[k]
                  - f_8 * pc_z[k] * ish1_3[k];

        t_46[k] = f_3 * pc_y[k] * ksg_32[k];

        t_47[k] = pa_z[k] * ish0_5[k]
                  + f_10 * isg_2[k]
                  - f_8 * pc_z[k] * ish1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, ish0_6, ish0_9, isg_5, \
                         ish1_6, ish1_9, ksf0_22, ksf1_22, ksg_34, \
                         ksg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * ish0_6[k]
                  - f_8 * pc_z[k] * ish1_6[k];

        t_49[k] = f_4 * ksf0_22[k]
                  - f_5 * ksf1_22[k]
                  + f_3 * pc_y[k] * ksg_34[k];

        t_50[k] = f_3 * pc_y[k] * ksg_35[k];

        t_51[k] = pa_z[k] * ish0_9[k]
                  + f_11 * isg_5[k]
                  - f_8 * pc_z[k] * ish1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, isg_40, isg_41, isg_42, \
                         isg_44, ksg_39, ksg_40, ksg_41, ksg_42, \
                         ksg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * isg_40[k]
                  + f_3 * pc_x[k] * ksg_40[k];

        t_53[k] = f_12 * isg_41[k]
                  + f_3 * pc_x[k] * ksg_41[k];

        t_54[k] = f_12 * isg_42[k]
                  + f_3 * pc_x[k] * ksg_42[k];

        t_55[k] = f_3 * pc_y[k] * ksg_39[k];

        t_56[k] = f_12 * isg_44[k]
                  + f_3 * pc_x[k] * ksg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, ish0_15, ish1_15, ksf0_27, \
                         ksf0_28, ksf1_27, ksf1_28, ksg_41, ksg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * ish0_15[k]
                  - f_8 * pc_z[k] * ish1_15[k];

        t_58[k] = f_13 * ksf0_27[k]
                  - f_14 * ksf1_27[k]
                  + f_3 * pc_y[k] * ksg_41[k];

        t_59[k] = f_6 * ksf0_28[k]
                  - f_7 * ksf1_28[k]
                  + f_3 * pc_y[k] * ksg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, isg_14, isg_45, ksf0_29, \
                         ksf0_30, ksf1_29, ksf1_30, ksg_43, ksg_44, \
                         ksg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * ksf0_29[k]
                  - f_5 * ksf1_29[k]
                  + f_3 * pc_y[k] * ksg_43[k];

        t_61[k] = f_3 * pc_y[k] * ksg_44[k];

        t_62[k] = f_9 * isg_14[k]
                  + f_1 * ksf0_29[k]
                  - f_2 * ksf1_29[k]
                  + f_3 * pc_z[k] * ksg_44[k];

        t_63[k] = f_15 * isg_45[k]
                  + f_1 * ksf0_30[k]
                  - f_2 * ksf1_30[k]
                  + f_3 * pc_x[k] * ksg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, isg_15, isg_48, ksf0_33, \
                         ksf1_33, ksg_45, ksg_46, ksg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * isg_15[k]
                  + f_3 * pc_y[k] * ksg_45[k];

        t_65[k] = f_3 * pc_z[k] * ksg_45[k];

        t_66[k] = f_15 * isg_48[k]
                  + f_6 * ksf0_33[k]
                  - f_7 * ksf1_33[k]
                  + f_3 * pc_x[k] * ksg_48[k];

        t_67[k] = f_3 * pc_z[k] * ksg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, isg_51, ksf0_30, ksf0_36, ksf1_30, \
                         ksf1_36, ksg_47, ksg_48, ksg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * ksf0_30[k]
                  - f_5 * ksf1_30[k]
                  + f_3 * pc_z[k] * ksg_47[k];

        t_69[k] = f_15 * isg_51[k]
                  + f_4 * ksf0_36[k]
                  - f_5 * ksf1_36[k]
                  + f_3 * pc_x[k] * ksg_51[k];

        t_70[k] = f_3 * pc_z[k] * ksg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, isg_20, isg_55, ksf0_32, \
                         ksf1_32, ksg_50, ksg_51, ksg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * isg_20[k]
                  + f_3 * pc_y[k] * ksg_50[k];

        t_72[k] = f_6 * ksf0_32[k]
                  - f_7 * ksf1_32[k]
                  + f_3 * pc_z[k] * ksg_50[k];

        t_73[k] = f_15 * isg_55[k]
                  + f_3 * pc_x[k] * ksg_55[k];

        t_74[k] = f_3 * pc_z[k] * ksg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, isg_25, isg_57, isg_58, isg_59, \
                         ksf0_36, ksf1_36, ksg_55, ksg_57, ksg_58, \
                         ksg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * isg_57[k]
                  + f_3 * pc_x[k] * ksg_57[k];

        t_76[k] = f_15 * isg_58[k]
                  + f_3 * pc_x[k] * ksg_58[k];

        t_77[k] = f_15 * isg_59[k]
                  + f_3 * pc_x[k] * ksg_59[k];

        t_78[k] = f_10 * isg_25[k]
                  + f_1 * ksf0_36[k]
                  - f_2 * ksf1_36[k]
                  + f_3 * pc_y[k] * ksg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, isg_29, ksf0_36, ksf0_37, \
                         ksf1_36, ksf1_37, ksg_55, ksg_56, ksg_57, \
                         ksg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * ksg_55[k];

        t_80[k] = f_4 * ksf0_36[k]
                  - f_5 * ksf1_36[k]
                  + f_3 * pc_z[k] * ksg_56[k];

        t_81[k] = f_6 * ksf0_37[k]
                  - f_7 * ksf1_37[k]
                  + f_3 * pc_z[k] * ksg_57[k];

        t_82[k] = f_10 * isg_29[k]
                  + f_3 * pc_y[k] * ksg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, ish0_42, isg_15, isg_30, \
                         ish1_42, ksf0_39, ksf1_39, ksg_59, ksg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * ksf0_39[k]
                  - f_2 * ksf1_39[k]
                  + f_3 * pc_z[k] * ksg_59[k];

        t_84[k] = pa_y[k] * ish0_42[k]
                  - f_8 * pc_y[k] * ish1_42[k];

        t_85[k] = f_9 * isg_30[k]
                  + f_3 * pc_y[k] * ksg_60[k];

        t_86[k] = f_9 * isg_15[k]
                  + f_3 * pc_z[k] * ksg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, ish0_24, ish0_27, \
                         ish0_47, isg_32, ish1_24, ish1_27, ish1_47, \
                         ksg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * ish0_24[k]
                  - f_8 * pc_z[k] * ish1_24[k];

        t_88[k] = f_9 * isg_32[k]
                  + f_3 * pc_y[k] * ksg_62[k];

        t_89[k] = pa_y[k] * ish0_47[k]
                  - f_8 * pc_y[k] * ish1_47[k];

        t_90[k] = pa_z[k] * ish0_27[k]
                  - f_8 * pc_z[k] * ish1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, ish0_51, isg_18, \
                         isg_35, isg_70, ish1_51, ksg_63, ksg_65, \
                         ksg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * isg_18[k]
                  + f_3 * pc_z[k] * ksg_63[k];

        t_92[k] = f_9 * isg_35[k]
                  + f_3 * pc_y[k] * ksg_65[k];

        t_93[k] = pa_y[k] * ish0_51[k]
                  - f_8 * pc_y[k] * ish1_51[k];

        t_94[k] = f_15 * isg_70[k]
                  + f_3 * pc_x[k] * ksg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, isg_71, isg_72, isg_73, isg_74, ksg_71, \
                         ksg_72, ksg_73, ksg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * isg_71[k]
                  + f_3 * pc_x[k] * ksg_71[k];

        t_96[k] = f_15 * isg_72[k]
                  + f_3 * pc_x[k] * ksg_72[k];

        t_97[k] = f_15 * isg_73[k]
                  + f_3 * pc_x[k] * ksg_73[k];

        t_98[k] = f_15 * isg_74[k]
                  + f_3 * pc_x[k] * ksg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, ish0_36, isg_25, isg_42, \
                         ish1_36, ksf0_48, ksf1_48, ksg_70, ksg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * ish0_36[k]
                  - f_8 * pc_z[k] * ish1_36[k];

        t_100[k] = f_9 * isg_25[k]
                   + f_3 * pc_z[k] * ksg_70[k];

        t_101[k] = f_9 * isg_42[k]
                   + f_6 * ksf0_48[k]
                   - f_7 * ksf1_48[k]
                   + f_3 * pc_y[k] * ksg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, ish0_62, isg_43, isg_44, ish1_62, \
                         ksf0_49, ksf1_49, ksg_73, ksg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * isg_43[k]
                   + f_4 * ksf0_49[k]
                   - f_5 * ksf1_49[k]
                   + f_3 * pc_y[k] * ksg_73[k];

        t_103[k] = f_9 * isg_44[k]
                   + f_3 * pc_y[k] * ksg_74[k];

        t_104[k] = pa_y[k] * ish0_62[k]
                   - f_8 * pc_y[k] * ish1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, isg_30, isg_75, \
                         ksf0_50, ksf1_50, ksg_75, ksg_76, ksg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * isg_75[k]
                   + f_1 * ksf0_50[k]
                   - f_2 * ksf1_50[k]
                   + f_3 * pc_x[k] * ksg_75[k];

        t_106[k] = f_3 * pc_y[k] * ksg_75[k];

        t_107[k] = f_10 * isg_30[k]
                   + f_3 * pc_z[k] * ksg_75[k];

        t_108[k] = f_4 * ksf0_50[k]
                   - f_5 * ksf1_50[k]
                   + f_3 * pc_y[k] * ksg_76[k];

        t_109[k] = f_3 * pc_y[k] * ksg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, isg_80, ksf0_51, ksf0_52, \
                         ksf0_55, ksf1_51, ksf1_52, ksf1_55, ksg_78, ksg_79, \
                         ksg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * isg_80[k]
                   + f_6 * ksf0_55[k]
                   - f_7 * ksf1_55[k]
                   + f_3 * pc_x[k] * ksg_80[k];

        t_111[k] = f_6 * ksf0_51[k]
                   - f_7 * ksf1_51[k]
                   + f_3 * pc_y[k] * ksg_78[k];

        t_112[k] = f_4 * ksf0_52[k]
                   - f_5 * ksf1_52[k]
                   + f_3 * pc_y[k] * ksg_79[k];

        t_113[k] = f_3 * pc_y[k] * ksg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, isg_84, isg_85, isg_86, isg_87, \
                         ksf0_59, ksf1_59, ksg_84, ksg_85, ksg_86, \
                         ksg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * isg_84[k]
                   + f_4 * ksf0_59[k]
                   - f_5 * ksf1_59[k]
                   + f_3 * pc_x[k] * ksg_84[k];

        t_115[k] = f_15 * isg_85[k]
                   + f_3 * pc_x[k] * ksg_85[k];

        t_116[k] = f_15 * isg_86[k]
                   + f_3 * pc_x[k] * ksg_86[k];

        t_117[k] = f_15 * isg_87[k]
                   + f_3 * pc_x[k] * ksg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, isg_89, ksf0_56, ksf0_57, \
                         ksf1_56, ksf1_57, ksg_84, ksg_85, ksg_86, \
                         ksg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * ksg_84[k];

        t_119[k] = f_15 * isg_89[k]
                   + f_3 * pc_x[k] * ksg_89[k];

        t_120[k] = f_1 * ksf0_56[k]
                   - f_2 * ksf1_56[k]
                   + f_3 * pc_y[k] * ksg_85[k];

        t_121[k] = f_13 * ksf0_57[k]
                   - f_14 * ksf1_57[k]
                   + f_3 * pc_y[k] * ksg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, isg_44, ksf0_58, ksf0_59, \
                         ksf1_58, ksf1_59, ksg_87, ksg_88, ksg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * ksf0_58[k]
                   - f_7 * ksf1_58[k]
                   + f_3 * pc_y[k] * ksg_87[k];

        t_123[k] = f_4 * ksf0_59[k]
                   - f_5 * ksf1_59[k]
                   + f_3 * pc_y[k] * ksg_88[k];

        t_124[k] = f_3 * pc_y[k] * ksg_89[k];

        t_125[k] = f_10 * isg_44[k]
                   + f_1 * ksf0_59[k]
                   - f_2 * ksf1_59[k]
                   + f_3 * pc_z[k] * ksg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, isg_45, isg_90, isg_93, \
                         ksf0_60, ksf0_63, ksf1_60, ksf1_63, ksg_90, \
                         ksg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_16 * isg_90[k]
                   + f_1 * ksf0_60[k]
                   - f_2 * ksf1_60[k]
                   + f_3 * pc_x[k] * ksg_90[k];

        t_127[k] = f_11 * isg_45[k]
                   + f_3 * pc_y[k] * ksg_90[k];

        t_128[k] = f_3 * pc_z[k] * ksg_90[k];

        t_129[k] = f_16 * isg_93[k]
                   + f_6 * ksf0_63[k]
                   - f_7 * ksf1_63[k]
                   + f_3 * pc_x[k] * ksg_93[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
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
    const auto f_16 = 2.0 / q;

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

    const auto *ish0_63 = buffer.data(ish0 + 63);
    const auto *ish0_66 = buffer.data(ish0 + 66);
    const auto *ish0_69 = buffer.data(ish0 + 69);
    const auto *ish0_78 = buffer.data(ish0 + 78);
    const auto *ish0_105 = buffer.data(ish0 + 105);
    const auto *ish0_108 = buffer.data(ish0 + 108);
    const auto *ish0_110 = buffer.data(ish0 + 110);
    const auto *ish0_111 = buffer.data(ish0 + 111);
    const auto *ish0_114 = buffer.data(ish0 + 114);
    const auto *ish0_125 = buffer.data(ish0 + 125);
    const auto *ish0_126 = buffer.data(ish0 + 126);
    const auto *ish0_129 = buffer.data(ish0 + 129);
    const auto *ish0_132 = buffer.data(ish0 + 132);
    const auto *ish0_141 = buffer.data(ish0 + 141);

    const auto *isg_45 = buffer.data(isg + 45);
    const auto *isg_48 = buffer.data(isg + 48);
    const auto *isg_50 = buffer.data(isg + 50);
    const auto *isg_55 = buffer.data(isg + 55);
    const auto *isg_59 = buffer.data(isg + 59);
    const auto *isg_60 = buffer.data(isg + 60);
    const auto *isg_62 = buffer.data(isg + 62);
    const auto *isg_63 = buffer.data(isg + 63);
    const auto *isg_65 = buffer.data(isg + 65);
    const auto *isg_70 = buffer.data(isg + 70);
    const auto *isg_72 = buffer.data(isg + 72);
    const auto *isg_73 = buffer.data(isg + 73);
    const auto *isg_74 = buffer.data(isg + 74);
    const auto *isg_75 = buffer.data(isg + 75);
    const auto *isg_76 = buffer.data(isg + 76);
    const auto *isg_77 = buffer.data(isg + 77);
    const auto *isg_78 = buffer.data(isg + 78);
    const auto *isg_80 = buffer.data(isg + 80);
    const auto *isg_85 = buffer.data(isg + 85);
    const auto *isg_87 = buffer.data(isg + 87);
    const auto *isg_88 = buffer.data(isg + 88);
    const auto *isg_89 = buffer.data(isg + 89);
    const auto *isg_90 = buffer.data(isg + 90);
    const auto *isg_93 = buffer.data(isg + 93);
    const auto *isg_95 = buffer.data(isg + 95);
    const auto *isg_96 = buffer.data(isg + 96);
    const auto *isg_100 = buffer.data(isg + 100);
    const auto *isg_102 = buffer.data(isg + 102);
    const auto *isg_103 = buffer.data(isg + 103);
    const auto *isg_104 = buffer.data(isg + 104);
    const auto *isg_105 = buffer.data(isg + 105);
    const auto *isg_107 = buffer.data(isg + 107);
    const auto *isg_110 = buffer.data(isg + 110);
    const auto *isg_114 = buffer.data(isg + 114);
    const auto *isg_115 = buffer.data(isg + 115);
    const auto *isg_116 = buffer.data(isg + 116);
    const auto *isg_117 = buffer.data(isg + 117);
    const auto *isg_118 = buffer.data(isg + 118);
    const auto *isg_119 = buffer.data(isg + 119);
    const auto *isg_130 = buffer.data(isg + 130);
    const auto *isg_131 = buffer.data(isg + 131);
    const auto *isg_132 = buffer.data(isg + 132);
    const auto *isg_133 = buffer.data(isg + 133);
    const auto *isg_134 = buffer.data(isg + 134);
    const auto *isg_135 = buffer.data(isg + 135);
    const auto *isg_140 = buffer.data(isg + 140);
    const auto *isg_144 = buffer.data(isg + 144);
    const auto *isg_145 = buffer.data(isg + 145);
    const auto *isg_146 = buffer.data(isg + 146);
    const auto *isg_147 = buffer.data(isg + 147);
    const auto *isg_149 = buffer.data(isg + 149);
    const auto *isg_150 = buffer.data(isg + 150);
    const auto *isg_153 = buffer.data(isg + 153);
    const auto *isg_156 = buffer.data(isg + 156);
    const auto *isg_160 = buffer.data(isg + 160);
    const auto *isg_162 = buffer.data(isg + 162);
    const auto *isg_163 = buffer.data(isg + 163);
    const auto *isg_164 = buffer.data(isg + 164);
    const auto *isg_170 = buffer.data(isg + 170);
    const auto *isg_174 = buffer.data(isg + 174);
    const auto *isg_175 = buffer.data(isg + 175);
    const auto *isg_176 = buffer.data(isg + 176);
    const auto *isg_177 = buffer.data(isg + 177);
    const auto *isg_178 = buffer.data(isg + 178);
    const auto *isg_179 = buffer.data(isg + 179);

    const auto *ish1_63 = buffer.data(ish1 + 63);
    const auto *ish1_66 = buffer.data(ish1 + 66);
    const auto *ish1_69 = buffer.data(ish1 + 69);
    const auto *ish1_78 = buffer.data(ish1 + 78);
    const auto *ish1_105 = buffer.data(ish1 + 105);
    const auto *ish1_108 = buffer.data(ish1 + 108);
    const auto *ish1_110 = buffer.data(ish1 + 110);
    const auto *ish1_111 = buffer.data(ish1 + 111);
    const auto *ish1_114 = buffer.data(ish1 + 114);
    const auto *ish1_125 = buffer.data(ish1 + 125);
    const auto *ish1_126 = buffer.data(ish1 + 126);
    const auto *ish1_129 = buffer.data(ish1 + 129);
    const auto *ish1_132 = buffer.data(ish1 + 132);
    const auto *ish1_141 = buffer.data(ish1 + 141);

    const auto *ksf0_60 = buffer.data(ksf0 + 60);
    const auto *ksf0_62 = buffer.data(ksf0 + 62);
    const auto *ksf0_66 = buffer.data(ksf0 + 66);
    const auto *ksf0_67 = buffer.data(ksf0 + 67);
    const auto *ksf0_69 = buffer.data(ksf0 + 69);
    const auto *ksf0_75 = buffer.data(ksf0 + 75);
    const auto *ksf0_78 = buffer.data(ksf0 + 78);
    const auto *ksf0_79 = buffer.data(ksf0 + 79);
    const auto *ksf0_86 = buffer.data(ksf0 + 86);
    const auto *ksf0_88 = buffer.data(ksf0 + 88);
    const auto *ksf0_89 = buffer.data(ksf0 + 89);
    const auto *ksf0_90 = buffer.data(ksf0 + 90);
    const auto *ksf0_91 = buffer.data(ksf0 + 91);
    const auto *ksf0_92 = buffer.data(ksf0 + 92);
    const auto *ksf0_95 = buffer.data(ksf0 + 95);
    const auto *ksf0_96 = buffer.data(ksf0 + 96);
    const auto *ksf0_97 = buffer.data(ksf0 + 97);
    const auto *ksf0_98 = buffer.data(ksf0 + 98);
    const auto *ksf0_99 = buffer.data(ksf0 + 99);
    const auto *ksf0_100 = buffer.data(ksf0 + 100);
    const auto *ksf0_102 = buffer.data(ksf0 + 102);
    const auto *ksf0_103 = buffer.data(ksf0 + 103);
    const auto *ksf0_106 = buffer.data(ksf0 + 106);
    const auto *ksf0_107 = buffer.data(ksf0 + 107);
    const auto *ksf0_109 = buffer.data(ksf0 + 109);
    const auto *ksf0_115 = buffer.data(ksf0 + 115);
    const auto *ksf0_118 = buffer.data(ksf0 + 118);
    const auto *ksf0_119 = buffer.data(ksf0 + 119);

    const auto *ksf1_60 = buffer.data(ksf1 + 60);
    const auto *ksf1_62 = buffer.data(ksf1 + 62);
    const auto *ksf1_66 = buffer.data(ksf1 + 66);
    const auto *ksf1_67 = buffer.data(ksf1 + 67);
    const auto *ksf1_69 = buffer.data(ksf1 + 69);
    const auto *ksf1_75 = buffer.data(ksf1 + 75);
    const auto *ksf1_78 = buffer.data(ksf1 + 78);
    const auto *ksf1_79 = buffer.data(ksf1 + 79);
    const auto *ksf1_86 = buffer.data(ksf1 + 86);
    const auto *ksf1_88 = buffer.data(ksf1 + 88);
    const auto *ksf1_89 = buffer.data(ksf1 + 89);
    const auto *ksf1_90 = buffer.data(ksf1 + 90);
    const auto *ksf1_91 = buffer.data(ksf1 + 91);
    const auto *ksf1_92 = buffer.data(ksf1 + 92);
    const auto *ksf1_95 = buffer.data(ksf1 + 95);
    const auto *ksf1_96 = buffer.data(ksf1 + 96);
    const auto *ksf1_97 = buffer.data(ksf1 + 97);
    const auto *ksf1_98 = buffer.data(ksf1 + 98);
    const auto *ksf1_99 = buffer.data(ksf1 + 99);
    const auto *ksf1_100 = buffer.data(ksf1 + 100);
    const auto *ksf1_102 = buffer.data(ksf1 + 102);
    const auto *ksf1_103 = buffer.data(ksf1 + 103);
    const auto *ksf1_106 = buffer.data(ksf1 + 106);
    const auto *ksf1_107 = buffer.data(ksf1 + 107);
    const auto *ksf1_109 = buffer.data(ksf1 + 109);
    const auto *ksf1_115 = buffer.data(ksf1 + 115);
    const auto *ksf1_118 = buffer.data(ksf1 + 118);
    const auto *ksf1_119 = buffer.data(ksf1 + 119);

    const auto *ksg_91 = buffer.data(ksg + 91);
    const auto *ksg_92 = buffer.data(ksg + 92);
    const auto *ksg_93 = buffer.data(ksg + 93);
    const auto *ksg_95 = buffer.data(ksg + 95);
    const auto *ksg_96 = buffer.data(ksg + 96);
    const auto *ksg_100 = buffer.data(ksg + 100);
    const auto *ksg_101 = buffer.data(ksg + 101);
    const auto *ksg_102 = buffer.data(ksg + 102);
    const auto *ksg_103 = buffer.data(ksg + 103);
    const auto *ksg_104 = buffer.data(ksg + 104);
    const auto *ksg_105 = buffer.data(ksg + 105);
    const auto *ksg_107 = buffer.data(ksg + 107);
    const auto *ksg_108 = buffer.data(ksg + 108);
    const auto *ksg_110 = buffer.data(ksg + 110);
    const auto *ksg_114 = buffer.data(ksg + 114);
    const auto *ksg_115 = buffer.data(ksg + 115);
    const auto *ksg_116 = buffer.data(ksg + 116);
    const auto *ksg_117 = buffer.data(ksg + 117);
    const auto *ksg_118 = buffer.data(ksg + 118);
    const auto *ksg_119 = buffer.data(ksg + 119);
    const auto *ksg_120 = buffer.data(ksg + 120);
    const auto *ksg_122 = buffer.data(ksg + 122);
    const auto *ksg_123 = buffer.data(ksg + 123);
    const auto *ksg_125 = buffer.data(ksg + 125);
    const auto *ksg_130 = buffer.data(ksg + 130);
    const auto *ksg_131 = buffer.data(ksg + 131);
    const auto *ksg_132 = buffer.data(ksg + 132);
    const auto *ksg_133 = buffer.data(ksg + 133);
    const auto *ksg_134 = buffer.data(ksg + 134);
    const auto *ksg_135 = buffer.data(ksg + 135);
    const auto *ksg_136 = buffer.data(ksg + 136);
    const auto *ksg_137 = buffer.data(ksg + 137);
    const auto *ksg_138 = buffer.data(ksg + 138);
    const auto *ksg_139 = buffer.data(ksg + 139);
    const auto *ksg_140 = buffer.data(ksg + 140);
    const auto *ksg_144 = buffer.data(ksg + 144);
    const auto *ksg_145 = buffer.data(ksg + 145);
    const auto *ksg_146 = buffer.data(ksg + 146);
    const auto *ksg_147 = buffer.data(ksg + 147);
    const auto *ksg_148 = buffer.data(ksg + 148);
    const auto *ksg_149 = buffer.data(ksg + 149);
    const auto *ksg_150 = buffer.data(ksg + 150);
    const auto *ksg_151 = buffer.data(ksg + 151);
    const auto *ksg_152 = buffer.data(ksg + 152);
    const auto *ksg_153 = buffer.data(ksg + 153);
    const auto *ksg_155 = buffer.data(ksg + 155);
    const auto *ksg_156 = buffer.data(ksg + 156);
    const auto *ksg_160 = buffer.data(ksg + 160);
    const auto *ksg_161 = buffer.data(ksg + 161);
    const auto *ksg_162 = buffer.data(ksg + 162);
    const auto *ksg_163 = buffer.data(ksg + 163);
    const auto *ksg_164 = buffer.data(ksg + 164);
    const auto *ksg_165 = buffer.data(ksg + 165);
    const auto *ksg_167 = buffer.data(ksg + 167);
    const auto *ksg_168 = buffer.data(ksg + 168);
    const auto *ksg_170 = buffer.data(ksg + 170);
    const auto *ksg_174 = buffer.data(ksg + 174);
    const auto *ksg_175 = buffer.data(ksg + 175);
    const auto *ksg_176 = buffer.data(ksg + 176);
    const auto *ksg_177 = buffer.data(ksg + 177);
    const auto *ksg_178 = buffer.data(ksg + 178);
    const auto *ksg_179 = buffer.data(ksg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, isg_96, ksf0_60, ksf0_66, \
                         ksf1_60, ksf1_66, ksg_91, ksg_92, ksg_93, \
                         ksg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * ksg_91[k];

        t_131[k] = f_4 * ksf0_60[k]
                   - f_5 * ksf1_60[k]
                   + f_3 * pc_z[k] * ksg_92[k];

        t_132[k] = f_16 * isg_96[k]
                   + f_4 * ksf0_66[k]
                   - f_5 * ksf1_66[k]
                   + f_3 * pc_x[k] * ksg_96[k];

        t_133[k] = f_3 * pc_z[k] * ksg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, isg_50, isg_100, \
                         ksf0_62, ksf1_62, ksg_95, ksg_96, ksg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * isg_50[k]
                   + f_3 * pc_y[k] * ksg_95[k];

        t_135[k] = f_6 * ksf0_62[k]
                   - f_7 * ksf1_62[k]
                   + f_3 * pc_z[k] * ksg_95[k];

        t_136[k] = f_16 * isg_100[k]
                   + f_3 * pc_x[k] * ksg_100[k];

        t_137[k] = f_3 * pc_z[k] * ksg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, isg_55, isg_102, isg_103, \
                         isg_104, ksf0_66, ksf1_66, ksg_100, ksg_102, ksg_103, \
                         ksg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_16 * isg_102[k]
                   + f_3 * pc_x[k] * ksg_102[k];

        t_139[k] = f_16 * isg_103[k]
                   + f_3 * pc_x[k] * ksg_103[k];

        t_140[k] = f_16 * isg_104[k]
                   + f_3 * pc_x[k] * ksg_104[k];

        t_141[k] = f_11 * isg_55[k]
                   + f_1 * ksf0_66[k]
                   - f_2 * ksf1_66[k]
                   + f_3 * pc_y[k] * ksg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, isg_59, ksf0_66, ksf0_67, \
                         ksf1_66, ksf1_67, ksg_100, ksg_101, ksg_102, \
                         ksg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * ksg_100[k];

        t_143[k] = f_4 * ksf0_66[k]
                   - f_5 * ksf1_66[k]
                   + f_3 * pc_z[k] * ksg_101[k];

        t_144[k] = f_6 * ksf0_67[k]
                   - f_7 * ksf1_67[k]
                   + f_3 * pc_z[k] * ksg_102[k];

        t_145[k] = f_11 * isg_59[k]
                   + f_3 * pc_y[k] * ksg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, ish0_63, isg_45, \
                         isg_60, ish1_63, ksf0_69, ksf1_69, ksg_104, \
                         ksg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * ksf0_69[k]
                   - f_2 * ksf1_69[k]
                   + f_3 * pc_z[k] * ksg_104[k];

        t_147[k] = pa_z[k] * ish0_63[k]
                   - f_8 * pc_z[k] * ish1_63[k];

        t_148[k] = f_10 * isg_60[k]
                   + f_3 * pc_y[k] * ksg_105[k];

        t_149[k] = f_9 * isg_45[k]
                   + f_3 * pc_z[k] * ksg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, ish0_66, isg_62, \
                         isg_110, ish1_66, ksf0_75, ksf1_75, ksg_107, \
                         ksg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * ish0_66[k]
                   - f_8 * pc_z[k] * ish1_66[k];

        t_151[k] = f_10 * isg_62[k]
                   + f_3 * pc_y[k] * ksg_107[k];

        t_152[k] = f_16 * isg_110[k]
                   + f_6 * ksf0_75[k]
                   - f_7 * ksf1_75[k]
                   + f_3 * pc_x[k] * ksg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, ish0_69, isg_48, isg_65, \
                         ish1_69, ksg_108, ksg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * ish0_69[k]
                   - f_8 * pc_z[k] * ish1_69[k];

        t_154[k] = f_9 * isg_48[k]
                   + f_3 * pc_z[k] * ksg_108[k];

        t_155[k] = f_10 * isg_65[k]
                   + f_3 * pc_y[k] * ksg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, isg_114, isg_115, isg_116, isg_117, \
                         ksf0_79, ksf1_79, ksg_114, ksg_115, ksg_116, \
                         ksg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_16 * isg_114[k]
                   + f_4 * ksf0_79[k]
                   - f_5 * ksf1_79[k]
                   + f_3 * pc_x[k] * ksg_114[k];

        t_157[k] = f_16 * isg_115[k]
                   + f_3 * pc_x[k] * ksg_115[k];

        t_158[k] = f_16 * isg_116[k]
                   + f_3 * pc_x[k] * ksg_116[k];

        t_159[k] = f_16 * isg_117[k]
                   + f_3 * pc_x[k] * ksg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, ish0_78, isg_55, \
                         isg_118, isg_119, ish1_78, ksg_115, ksg_118, \
                         ksg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_16 * isg_118[k]
                   + f_3 * pc_x[k] * ksg_118[k];

        t_161[k] = f_16 * isg_119[k]
                   + f_3 * pc_x[k] * ksg_119[k];

        t_162[k] = pa_z[k] * ish0_78[k]
                   - f_8 * pc_z[k] * ish1_78[k];

        t_163[k] = f_9 * isg_55[k]
                   + f_3 * pc_z[k] * ksg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, isg_72, isg_73, isg_74, ksf0_78, ksf0_79, \
                         ksf1_78, ksf1_79, ksg_117, ksg_118, ksg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * isg_72[k]
                   + f_6 * ksf0_78[k]
                   - f_7 * ksf1_78[k]
                   + f_3 * pc_y[k] * ksg_117[k];

        t_165[k] = f_10 * isg_73[k]
                   + f_4 * ksf0_79[k]
                   - f_5 * ksf1_79[k]
                   + f_3 * pc_y[k] * ksg_118[k];

        t_166[k] = f_10 * isg_74[k]
                   + f_3 * pc_y[k] * ksg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, ish0_105, isg_59, \
                         isg_60, isg_75, ish1_105, ksf0_79, ksf1_79, ksg_119, \
                         ksg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * isg_59[k]
                   + f_1 * ksf0_79[k]
                   - f_2 * ksf1_79[k]
                   + f_3 * pc_z[k] * ksg_119[k];

        t_168[k] = pa_y[k] * ish0_105[k]
                   - f_8 * pc_y[k] * ish1_105[k];

        t_169[k] = f_9 * isg_75[k]
                   + f_3 * pc_y[k] * ksg_120[k];

        t_170[k] = f_10 * isg_60[k]
                   + f_3 * pc_z[k] * ksg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, ish0_108, ish0_110, ish0_111, \
                         isg_76, isg_77, isg_78, ish1_108, ish1_110, ish1_111, \
                         ksg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * ish0_108[k]
                   + f_10 * isg_76[k]
                   - f_8 * pc_y[k] * ish1_108[k];

        t_172[k] = f_9 * isg_77[k]
                   + f_3 * pc_y[k] * ksg_122[k];

        t_173[k] = pa_y[k] * ish0_110[k]
                   - f_8 * pc_y[k] * ish1_110[k];

        t_174[k] = pa_y[k] * ish0_111[k]
                   + f_11 * isg_78[k]
                   - f_8 * pc_y[k] * ish1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, ish0_114, isg_63, \
                         isg_80, isg_130, ish1_114, ksg_123, ksg_125, \
                         ksg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * isg_63[k]
                   + f_3 * pc_z[k] * ksg_123[k];

        t_176[k] = f_9 * isg_80[k]
                   + f_3 * pc_y[k] * ksg_125[k];

        t_177[k] = pa_y[k] * ish0_114[k]
                   - f_8 * pc_y[k] * ish1_114[k];

        t_178[k] = f_16 * isg_130[k]
                   + f_3 * pc_x[k] * ksg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, isg_131, isg_132, isg_133, isg_134, \
                         ksg_131, ksg_132, ksg_133, ksg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_16 * isg_131[k]
                   + f_3 * pc_x[k] * ksg_131[k];

        t_180[k] = f_16 * isg_132[k]
                   + f_3 * pc_x[k] * ksg_132[k];

        t_181[k] = f_16 * isg_133[k]
                   + f_3 * pc_x[k] * ksg_133[k];

        t_182[k] = f_16 * isg_134[k]
                   + f_3 * pc_x[k] * ksg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, isg_70, isg_85, isg_87, ksf0_86, \
                         ksf0_88, ksf1_86, ksf1_88, ksg_130, ksg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * isg_85[k]
                   + f_1 * ksf0_86[k]
                   - f_2 * ksf1_86[k]
                   + f_3 * pc_y[k] * ksg_130[k];

        t_184[k] = f_10 * isg_70[k]
                   + f_3 * pc_z[k] * ksg_130[k];

        t_185[k] = f_9 * isg_87[k]
                   + f_6 * ksf0_88[k]
                   - f_7 * ksf1_88[k]
                   + f_3 * pc_y[k] * ksg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, ish0_125, isg_88, isg_89, ish1_125, \
                         ksf0_89, ksf1_89, ksg_133, ksg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * isg_88[k]
                   + f_4 * ksf0_89[k]
                   - f_5 * ksf1_89[k]
                   + f_3 * pc_y[k] * ksg_133[k];

        t_187[k] = f_9 * isg_89[k]
                   + f_3 * pc_y[k] * ksg_134[k];

        t_188[k] = pa_y[k] * ish0_125[k]
                   - f_8 * pc_y[k] * ish1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, isg_75, isg_135, \
                         ksf0_90, ksf1_90, ksg_135, ksg_136, ksg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_16 * isg_135[k]
                   + f_1 * ksf0_90[k]
                   - f_2 * ksf1_90[k]
                   + f_3 * pc_x[k] * ksg_135[k];

        t_190[k] = f_3 * pc_y[k] * ksg_135[k];

        t_191[k] = f_11 * isg_75[k]
                   + f_3 * pc_z[k] * ksg_135[k];

        t_192[k] = f_4 * ksf0_90[k]
                   - f_5 * ksf1_90[k]
                   + f_3 * pc_y[k] * ksg_136[k];

        t_193[k] = f_3 * pc_y[k] * ksg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, isg_140, ksf0_91, ksf0_92, \
                         ksf0_95, ksf1_91, ksf1_92, ksf1_95, ksg_138, ksg_139, \
                         ksg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_16 * isg_140[k]
                   + f_6 * ksf0_95[k]
                   - f_7 * ksf1_95[k]
                   + f_3 * pc_x[k] * ksg_140[k];

        t_195[k] = f_6 * ksf0_91[k]
                   - f_7 * ksf1_91[k]
                   + f_3 * pc_y[k] * ksg_138[k];

        t_196[k] = f_4 * ksf0_92[k]
                   - f_5 * ksf1_92[k]
                   + f_3 * pc_y[k] * ksg_139[k];

        t_197[k] = f_3 * pc_y[k] * ksg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, isg_144, isg_145, isg_146, isg_147, \
                         ksf0_99, ksf1_99, ksg_144, ksg_145, ksg_146, \
                         ksg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_16 * isg_144[k]
                   + f_4 * ksf0_99[k]
                   - f_5 * ksf1_99[k]
                   + f_3 * pc_x[k] * ksg_144[k];

        t_199[k] = f_16 * isg_145[k]
                   + f_3 * pc_x[k] * ksg_145[k];

        t_200[k] = f_16 * isg_146[k]
                   + f_3 * pc_x[k] * ksg_146[k];

        t_201[k] = f_16 * isg_147[k]
                   + f_3 * pc_x[k] * ksg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, isg_149, ksf0_96, ksf0_97, \
                         ksf1_96, ksf1_97, ksg_144, ksg_145, ksg_146, \
                         ksg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * ksg_144[k];

        t_203[k] = f_16 * isg_149[k]
                   + f_3 * pc_x[k] * ksg_149[k];

        t_204[k] = f_1 * ksf0_96[k]
                   - f_2 * ksf1_96[k]
                   + f_3 * pc_y[k] * ksg_145[k];

        t_205[k] = f_13 * ksf0_97[k]
                   - f_14 * ksf1_97[k]
                   + f_3 * pc_y[k] * ksg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, isg_89, ksf0_98, ksf0_99, \
                         ksf1_98, ksf1_99, ksg_147, ksg_148, ksg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * ksf0_98[k]
                   - f_7 * ksf1_98[k]
                   + f_3 * pc_y[k] * ksg_147[k];

        t_207[k] = f_4 * ksf0_99[k]
                   - f_5 * ksf1_99[k]
                   + f_3 * pc_y[k] * ksg_148[k];

        t_208[k] = f_3 * pc_y[k] * ksg_149[k];

        t_209[k] = f_11 * isg_89[k]
                   + f_1 * ksf0_99[k]
                   - f_2 * ksf1_99[k]
                   + f_3 * pc_z[k] * ksg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, isg_90, isg_150, \
                         isg_153, ksf0_100, ksf0_103, ksf1_100, ksf1_103, ksg_150, \
                         ksg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_11 * isg_150[k]
                   + f_1 * ksf0_100[k]
                   - f_2 * ksf1_100[k]
                   + f_3 * pc_x[k] * ksg_150[k];

        t_211[k] = f_16 * isg_90[k]
                   + f_3 * pc_y[k] * ksg_150[k];

        t_212[k] = f_3 * pc_z[k] * ksg_150[k];

        t_213[k] = f_11 * isg_153[k]
                   + f_6 * ksf0_103[k]
                   - f_7 * ksf1_103[k]
                   + f_3 * pc_x[k] * ksg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, isg_156, ksf0_100, ksf0_106, \
                         ksf1_100, ksf1_106, ksg_151, ksg_152, ksg_153, \
                         ksg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * ksg_151[k];

        t_215[k] = f_4 * ksf0_100[k]
                   - f_5 * ksf1_100[k]
                   + f_3 * pc_z[k] * ksg_152[k];

        t_216[k] = f_11 * isg_156[k]
                   + f_4 * ksf0_106[k]
                   - f_5 * ksf1_106[k]
                   + f_3 * pc_x[k] * ksg_156[k];

        t_217[k] = f_3 * pc_z[k] * ksg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, isg_95, isg_160, \
                         ksf0_102, ksf1_102, ksg_155, ksg_156, \
                         ksg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_16 * isg_95[k]
                   + f_3 * pc_y[k] * ksg_155[k];

        t_219[k] = f_6 * ksf0_102[k]
                   - f_7 * ksf1_102[k]
                   + f_3 * pc_z[k] * ksg_155[k];

        t_220[k] = f_11 * isg_160[k]
                   + f_3 * pc_x[k] * ksg_160[k];

        t_221[k] = f_3 * pc_z[k] * ksg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, isg_100, isg_162, isg_163, \
                         isg_164, ksf0_106, ksf1_106, ksg_160, ksg_162, ksg_163, \
                         ksg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_11 * isg_162[k]
                   + f_3 * pc_x[k] * ksg_162[k];

        t_223[k] = f_11 * isg_163[k]
                   + f_3 * pc_x[k] * ksg_163[k];

        t_224[k] = f_11 * isg_164[k]
                   + f_3 * pc_x[k] * ksg_164[k];

        t_225[k] = f_16 * isg_100[k]
                   + f_1 * ksf0_106[k]
                   - f_2 * ksf1_106[k]
                   + f_3 * pc_y[k] * ksg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, isg_104, ksf0_106, ksf0_107, \
                         ksf1_106, ksf1_107, ksg_160, ksg_161, ksg_162, \
                         ksg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * ksg_160[k];

        t_227[k] = f_4 * ksf0_106[k]
                   - f_5 * ksf1_106[k]
                   + f_3 * pc_z[k] * ksg_161[k];

        t_228[k] = f_6 * ksf0_107[k]
                   - f_7 * ksf1_107[k]
                   + f_3 * pc_z[k] * ksg_162[k];

        t_229[k] = f_16 * isg_104[k]
                   + f_3 * pc_y[k] * ksg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, ish0_126, isg_90, \
                         isg_105, ish1_126, ksf0_109, ksf1_109, ksg_164, \
                         ksg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * ksf0_109[k]
                   - f_2 * ksf1_109[k]
                   + f_3 * pc_z[k] * ksg_164[k];

        t_231[k] = pa_z[k] * ish0_126[k]
                   - f_8 * pc_z[k] * ish1_126[k];

        t_232[k] = f_11 * isg_105[k]
                   + f_3 * pc_y[k] * ksg_165[k];

        t_233[k] = f_9 * isg_90[k]
                   + f_3 * pc_z[k] * ksg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, ish0_129, isg_107, \
                         isg_170, ish1_129, ksf0_115, ksf1_115, ksg_167, \
                         ksg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * ish0_129[k]
                   - f_8 * pc_z[k] * ish1_129[k];

        t_235[k] = f_11 * isg_107[k]
                   + f_3 * pc_y[k] * ksg_167[k];

        t_236[k] = f_11 * isg_170[k]
                   + f_6 * ksf0_115[k]
                   - f_7 * ksf1_115[k]
                   + f_3 * pc_x[k] * ksg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, ish0_132, isg_93, isg_110, \
                         ish1_132, ksg_168, ksg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * ish0_132[k]
                   - f_8 * pc_z[k] * ish1_132[k];

        t_238[k] = f_9 * isg_93[k]
                   + f_3 * pc_z[k] * ksg_168[k];

        t_239[k] = f_11 * isg_110[k]
                   + f_3 * pc_y[k] * ksg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, isg_174, isg_175, isg_176, isg_177, \
                         ksf0_119, ksf1_119, ksg_174, ksg_175, ksg_176, \
                         ksg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_11 * isg_174[k]
                   + f_4 * ksf0_119[k]
                   - f_5 * ksf1_119[k]
                   + f_3 * pc_x[k] * ksg_174[k];

        t_241[k] = f_11 * isg_175[k]
                   + f_3 * pc_x[k] * ksg_175[k];

        t_242[k] = f_11 * isg_176[k]
                   + f_3 * pc_x[k] * ksg_176[k];

        t_243[k] = f_11 * isg_177[k]
                   + f_3 * pc_x[k] * ksg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, ish0_141, isg_100, \
                         isg_178, isg_179, ish1_141, ksg_175, ksg_178, \
                         ksg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_11 * isg_178[k]
                   + f_3 * pc_x[k] * ksg_178[k];

        t_245[k] = f_11 * isg_179[k]
                   + f_3 * pc_x[k] * ksg_179[k];

        t_246[k] = pa_z[k] * ish0_141[k]
                   - f_8 * pc_z[k] * ish1_141[k];

        t_247[k] = f_9 * isg_100[k]
                   + f_3 * pc_z[k] * ksg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, isg_117, isg_118, isg_119, ksf0_118, \
                         ksf0_119, ksf1_118, ksf1_119, ksg_177, ksg_178, \
                         ksg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * isg_117[k]
                   + f_6 * ksf0_118[k]
                   - f_7 * ksf1_118[k]
                   + f_3 * pc_y[k] * ksg_177[k];

        t_249[k] = f_11 * isg_118[k]
                   + f_4 * ksf0_119[k]
                   - f_5 * ksf1_119[k]
                   + f_3 * pc_y[k] * ksg_178[k];

        t_250[k] = f_11 * isg_119[k]
                   + f_3 * pc_y[k] * ksg_179[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *ish0_189 = buffer.data(ish0 + 189);
    const auto *ish0_192 = buffer.data(ish0 + 192);
    const auto *ish0_194 = buffer.data(ish0 + 194);
    const auto *ish0_195 = buffer.data(ish0 + 195);
    const auto *ish0_198 = buffer.data(ish0 + 198);
    const auto *ish0_209 = buffer.data(ish0 + 209);
    const auto *ish0_210 = buffer.data(ish0 + 210);
    const auto *ish0_213 = buffer.data(ish0 + 213);
    const auto *ish0_216 = buffer.data(ish0 + 216);
    const auto *ish0_225 = buffer.data(ish0 + 225);

    const auto *isg_104 = buffer.data(isg + 104);
    const auto *isg_105 = buffer.data(isg + 105);
    const auto *isg_108 = buffer.data(isg + 108);
    const auto *isg_115 = buffer.data(isg + 115);
    const auto *isg_119 = buffer.data(isg + 119);
    const auto *isg_120 = buffer.data(isg + 120);
    const auto *isg_122 = buffer.data(isg + 122);
    const auto *isg_123 = buffer.data(isg + 123);
    const auto *isg_125 = buffer.data(isg + 125);
    const auto *isg_130 = buffer.data(isg + 130);
    const auto *isg_132 = buffer.data(isg + 132);
    const auto *isg_133 = buffer.data(isg + 133);
    const auto *isg_134 = buffer.data(isg + 134);
    const auto *isg_135 = buffer.data(isg + 135);
    const auto *isg_136 = buffer.data(isg + 136);
    const auto *isg_137 = buffer.data(isg + 137);
    const auto *isg_138 = buffer.data(isg + 138);
    const auto *isg_140 = buffer.data(isg + 140);
    const auto *isg_145 = buffer.data(isg + 145);
    const auto *isg_147 = buffer.data(isg + 147);
    const auto *isg_148 = buffer.data(isg + 148);
    const auto *isg_149 = buffer.data(isg + 149);
    const auto *isg_150 = buffer.data(isg + 150);
    const auto *isg_153 = buffer.data(isg + 153);
    const auto *isg_155 = buffer.data(isg + 155);
    const auto *isg_160 = buffer.data(isg + 160);
    const auto *isg_164 = buffer.data(isg + 164);
    const auto *isg_165 = buffer.data(isg + 165);
    const auto *isg_167 = buffer.data(isg + 167);
    const auto *isg_168 = buffer.data(isg + 168);
    const auto *isg_170 = buffer.data(isg + 170);
    const auto *isg_177 = buffer.data(isg + 177);
    const auto *isg_178 = buffer.data(isg + 178);
    const auto *isg_179 = buffer.data(isg + 179);
    const auto *isg_180 = buffer.data(isg + 180);
    const auto *isg_182 = buffer.data(isg + 182);
    const auto *isg_183 = buffer.data(isg + 183);
    const auto *isg_185 = buffer.data(isg + 185);
    const auto *isg_186 = buffer.data(isg + 186);
    const auto *isg_189 = buffer.data(isg + 189);
    const auto *isg_190 = buffer.data(isg + 190);
    const auto *isg_191 = buffer.data(isg + 191);
    const auto *isg_192 = buffer.data(isg + 192);
    const auto *isg_193 = buffer.data(isg + 193);
    const auto *isg_194 = buffer.data(isg + 194);
    const auto *isg_205 = buffer.data(isg + 205);
    const auto *isg_206 = buffer.data(isg + 206);
    const auto *isg_207 = buffer.data(isg + 207);
    const auto *isg_208 = buffer.data(isg + 208);
    const auto *isg_209 = buffer.data(isg + 209);
    const auto *isg_210 = buffer.data(isg + 210);
    const auto *isg_215 = buffer.data(isg + 215);
    const auto *isg_219 = buffer.data(isg + 219);
    const auto *isg_220 = buffer.data(isg + 220);
    const auto *isg_221 = buffer.data(isg + 221);
    const auto *isg_222 = buffer.data(isg + 222);
    const auto *isg_224 = buffer.data(isg + 224);
    const auto *isg_225 = buffer.data(isg + 225);
    const auto *isg_228 = buffer.data(isg + 228);
    const auto *isg_231 = buffer.data(isg + 231);
    const auto *isg_235 = buffer.data(isg + 235);
    const auto *isg_237 = buffer.data(isg + 237);
    const auto *isg_238 = buffer.data(isg + 238);
    const auto *isg_239 = buffer.data(isg + 239);
    const auto *isg_245 = buffer.data(isg + 245);
    const auto *isg_249 = buffer.data(isg + 249);
    const auto *isg_250 = buffer.data(isg + 250);
    const auto *isg_251 = buffer.data(isg + 251);
    const auto *isg_252 = buffer.data(isg + 252);
    const auto *isg_253 = buffer.data(isg + 253);
    const auto *isg_254 = buffer.data(isg + 254);
    const auto *isg_255 = buffer.data(isg + 255);
    const auto *isg_258 = buffer.data(isg + 258);
    const auto *isg_260 = buffer.data(isg + 260);
    const auto *isg_261 = buffer.data(isg + 261);
    const auto *isg_264 = buffer.data(isg + 264);
    const auto *isg_265 = buffer.data(isg + 265);
    const auto *isg_266 = buffer.data(isg + 266);

    const auto *ish1_189 = buffer.data(ish1 + 189);
    const auto *ish1_192 = buffer.data(ish1 + 192);
    const auto *ish1_194 = buffer.data(ish1 + 194);
    const auto *ish1_195 = buffer.data(ish1 + 195);
    const auto *ish1_198 = buffer.data(ish1 + 198);
    const auto *ish1_209 = buffer.data(ish1 + 209);
    const auto *ish1_210 = buffer.data(ish1 + 210);
    const auto *ish1_213 = buffer.data(ish1 + 213);
    const auto *ish1_216 = buffer.data(ish1 + 216);
    const auto *ish1_225 = buffer.data(ish1 + 225);

    const auto *ksf0_119 = buffer.data(ksf0 + 119);
    const auto *ksf0_120 = buffer.data(ksf0 + 120);
    const auto *ksf0_123 = buffer.data(ksf0 + 123);
    const auto *ksf0_125 = buffer.data(ksf0 + 125);
    const auto *ksf0_126 = buffer.data(ksf0 + 126);
    const auto *ksf0_128 = buffer.data(ksf0 + 128);
    const auto *ksf0_129 = buffer.data(ksf0 + 129);
    const auto *ksf0_136 = buffer.data(ksf0 + 136);
    const auto *ksf0_138 = buffer.data(ksf0 + 138);
    const auto *ksf0_139 = buffer.data(ksf0 + 139);
    const auto *ksf0_140 = buffer.data(ksf0 + 140);
    const auto *ksf0_141 = buffer.data(ksf0 + 141);
    const auto *ksf0_142 = buffer.data(ksf0 + 142);
    const auto *ksf0_145 = buffer.data(ksf0 + 145);
    const auto *ksf0_146 = buffer.data(ksf0 + 146);
    const auto *ksf0_147 = buffer.data(ksf0 + 147);
    const auto *ksf0_148 = buffer.data(ksf0 + 148);
    const auto *ksf0_149 = buffer.data(ksf0 + 149);
    const auto *ksf0_150 = buffer.data(ksf0 + 150);
    const auto *ksf0_152 = buffer.data(ksf0 + 152);
    const auto *ksf0_153 = buffer.data(ksf0 + 153);
    const auto *ksf0_156 = buffer.data(ksf0 + 156);
    const auto *ksf0_157 = buffer.data(ksf0 + 157);
    const auto *ksf0_159 = buffer.data(ksf0 + 159);
    const auto *ksf0_165 = buffer.data(ksf0 + 165);
    const auto *ksf0_168 = buffer.data(ksf0 + 168);
    const auto *ksf0_169 = buffer.data(ksf0 + 169);
    const auto *ksf0_170 = buffer.data(ksf0 + 170);
    const auto *ksf0_173 = buffer.data(ksf0 + 173);
    const auto *ksf0_175 = buffer.data(ksf0 + 175);
    const auto *ksf0_176 = buffer.data(ksf0 + 176);
    const auto *ksf0_179 = buffer.data(ksf0 + 179);

    const auto *ksf1_119 = buffer.data(ksf1 + 119);
    const auto *ksf1_120 = buffer.data(ksf1 + 120);
    const auto *ksf1_123 = buffer.data(ksf1 + 123);
    const auto *ksf1_125 = buffer.data(ksf1 + 125);
    const auto *ksf1_126 = buffer.data(ksf1 + 126);
    const auto *ksf1_128 = buffer.data(ksf1 + 128);
    const auto *ksf1_129 = buffer.data(ksf1 + 129);
    const auto *ksf1_136 = buffer.data(ksf1 + 136);
    const auto *ksf1_138 = buffer.data(ksf1 + 138);
    const auto *ksf1_139 = buffer.data(ksf1 + 139);
    const auto *ksf1_140 = buffer.data(ksf1 + 140);
    const auto *ksf1_141 = buffer.data(ksf1 + 141);
    const auto *ksf1_142 = buffer.data(ksf1 + 142);
    const auto *ksf1_145 = buffer.data(ksf1 + 145);
    const auto *ksf1_146 = buffer.data(ksf1 + 146);
    const auto *ksf1_147 = buffer.data(ksf1 + 147);
    const auto *ksf1_148 = buffer.data(ksf1 + 148);
    const auto *ksf1_149 = buffer.data(ksf1 + 149);
    const auto *ksf1_150 = buffer.data(ksf1 + 150);
    const auto *ksf1_152 = buffer.data(ksf1 + 152);
    const auto *ksf1_153 = buffer.data(ksf1 + 153);
    const auto *ksf1_156 = buffer.data(ksf1 + 156);
    const auto *ksf1_157 = buffer.data(ksf1 + 157);
    const auto *ksf1_159 = buffer.data(ksf1 + 159);
    const auto *ksf1_165 = buffer.data(ksf1 + 165);
    const auto *ksf1_168 = buffer.data(ksf1 + 168);
    const auto *ksf1_169 = buffer.data(ksf1 + 169);
    const auto *ksf1_170 = buffer.data(ksf1 + 170);
    const auto *ksf1_173 = buffer.data(ksf1 + 173);
    const auto *ksf1_175 = buffer.data(ksf1 + 175);
    const auto *ksf1_176 = buffer.data(ksf1 + 176);
    const auto *ksf1_179 = buffer.data(ksf1 + 179);

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
    const auto *ksg_195 = buffer.data(ksg + 195);
    const auto *ksg_197 = buffer.data(ksg + 197);
    const auto *ksg_198 = buffer.data(ksg + 198);
    const auto *ksg_200 = buffer.data(ksg + 200);
    const auto *ksg_205 = buffer.data(ksg + 205);
    const auto *ksg_206 = buffer.data(ksg + 206);
    const auto *ksg_207 = buffer.data(ksg + 207);
    const auto *ksg_208 = buffer.data(ksg + 208);
    const auto *ksg_209 = buffer.data(ksg + 209);
    const auto *ksg_210 = buffer.data(ksg + 210);
    const auto *ksg_211 = buffer.data(ksg + 211);
    const auto *ksg_212 = buffer.data(ksg + 212);
    const auto *ksg_213 = buffer.data(ksg + 213);
    const auto *ksg_214 = buffer.data(ksg + 214);
    const auto *ksg_215 = buffer.data(ksg + 215);
    const auto *ksg_219 = buffer.data(ksg + 219);
    const auto *ksg_220 = buffer.data(ksg + 220);
    const auto *ksg_221 = buffer.data(ksg + 221);
    const auto *ksg_222 = buffer.data(ksg + 222);
    const auto *ksg_223 = buffer.data(ksg + 223);
    const auto *ksg_224 = buffer.data(ksg + 224);
    const auto *ksg_225 = buffer.data(ksg + 225);
    const auto *ksg_226 = buffer.data(ksg + 226);
    const auto *ksg_227 = buffer.data(ksg + 227);
    const auto *ksg_228 = buffer.data(ksg + 228);
    const auto *ksg_230 = buffer.data(ksg + 230);
    const auto *ksg_231 = buffer.data(ksg + 231);
    const auto *ksg_235 = buffer.data(ksg + 235);
    const auto *ksg_236 = buffer.data(ksg + 236);
    const auto *ksg_237 = buffer.data(ksg + 237);
    const auto *ksg_238 = buffer.data(ksg + 238);
    const auto *ksg_239 = buffer.data(ksg + 239);
    const auto *ksg_240 = buffer.data(ksg + 240);
    const auto *ksg_242 = buffer.data(ksg + 242);
    const auto *ksg_243 = buffer.data(ksg + 243);
    const auto *ksg_245 = buffer.data(ksg + 245);
    const auto *ksg_249 = buffer.data(ksg + 249);
    const auto *ksg_250 = buffer.data(ksg + 250);
    const auto *ksg_251 = buffer.data(ksg + 251);
    const auto *ksg_252 = buffer.data(ksg + 252);
    const auto *ksg_253 = buffer.data(ksg + 253);
    const auto *ksg_254 = buffer.data(ksg + 254);
    const auto *ksg_255 = buffer.data(ksg + 255);
    const auto *ksg_257 = buffer.data(ksg + 257);
    const auto *ksg_258 = buffer.data(ksg + 258);
    const auto *ksg_260 = buffer.data(ksg + 260);
    const auto *ksg_261 = buffer.data(ksg + 261);
    const auto *ksg_264 = buffer.data(ksg + 264);
    const auto *ksg_265 = buffer.data(ksg + 265);
    const auto *ksg_266 = buffer.data(ksg + 266);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, isg_104, isg_120, isg_180, \
                         ksf0_119, ksf0_120, ksf1_119, ksf1_120, ksg_179, \
                         ksg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * isg_104[k]
                   + f_1 * ksf0_119[k]
                   - f_2 * ksf1_119[k]
                   + f_3 * pc_z[k] * ksg_179[k];

        t_252[k] = f_11 * isg_180[k]
                   + f_1 * ksf0_120[k]
                   - f_2 * ksf1_120[k]
                   + f_3 * pc_x[k] * ksg_180[k];

        t_253[k] = f_10 * isg_120[k]
                   + f_3 * pc_y[k] * ksg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, isg_105, isg_122, isg_183, \
                         ksf0_123, ksf1_123, ksg_180, ksg_182, \
                         ksg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * isg_105[k]
                   + f_3 * pc_z[k] * ksg_180[k];

        t_255[k] = f_11 * isg_183[k]
                   + f_6 * ksf0_123[k]
                   - f_7 * ksf1_123[k]
                   + f_3 * pc_x[k] * ksg_183[k];

        t_256[k] = f_10 * isg_122[k]
                   + f_3 * pc_y[k] * ksg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, isg_108, isg_185, isg_186, ksf0_125, \
                         ksf0_126, ksf1_125, ksf1_126, ksg_183, ksg_185, \
                         ksg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_11 * isg_185[k]
                   + f_6 * ksf0_125[k]
                   - f_7 * ksf1_125[k]
                   + f_3 * pc_x[k] * ksg_185[k];

        t_258[k] = f_11 * isg_186[k]
                   + f_4 * ksf0_126[k]
                   - f_5 * ksf1_126[k]
                   + f_3 * pc_x[k] * ksg_186[k];

        t_259[k] = f_10 * isg_108[k]
                   + f_3 * pc_z[k] * ksg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, isg_125, isg_189, isg_190, \
                         isg_191, ksf0_129, ksf1_129, ksg_185, ksg_189, ksg_190, \
                         ksg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * isg_125[k]
                   + f_3 * pc_y[k] * ksg_185[k];

        t_261[k] = f_11 * isg_189[k]
                   + f_4 * ksf0_129[k]
                   - f_5 * ksf1_129[k]
                   + f_3 * pc_x[k] * ksg_189[k];

        t_262[k] = f_11 * isg_190[k]
                   + f_3 * pc_x[k] * ksg_190[k];

        t_263[k] = f_11 * isg_191[k]
                   + f_3 * pc_x[k] * ksg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, isg_130, isg_192, isg_193, \
                         isg_194, ksf0_126, ksf1_126, ksg_190, ksg_192, ksg_193, \
                         ksg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_11 * isg_192[k]
                   + f_3 * pc_x[k] * ksg_192[k];

        t_265[k] = f_11 * isg_193[k]
                   + f_3 * pc_x[k] * ksg_193[k];

        t_266[k] = f_11 * isg_194[k]
                   + f_3 * pc_x[k] * ksg_194[k];

        t_267[k] = f_10 * isg_130[k]
                   + f_1 * ksf0_126[k]
                   - f_2 * ksf1_126[k]
                   + f_3 * pc_y[k] * ksg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, isg_115, isg_132, isg_133, ksf0_128, \
                         ksf0_129, ksf1_128, ksf1_129, ksg_190, ksg_192, \
                         ksg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * isg_115[k]
                   + f_3 * pc_z[k] * ksg_190[k];

        t_269[k] = f_10 * isg_132[k]
                   + f_6 * ksf0_128[k]
                   - f_7 * ksf1_128[k]
                   + f_3 * pc_y[k] * ksg_192[k];

        t_270[k] = f_10 * isg_133[k]
                   + f_4 * ksf0_129[k]
                   - f_5 * ksf1_129[k]
                   + f_3 * pc_y[k] * ksg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, ish0_189, isg_119, \
                         isg_134, isg_135, ish1_189, ksf0_129, ksf1_129, ksg_194, \
                         ksg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * isg_134[k]
                   + f_3 * pc_y[k] * ksg_194[k];

        t_272[k] = f_10 * isg_119[k]
                   + f_1 * ksf0_129[k]
                   - f_2 * ksf1_129[k]
                   + f_3 * pc_z[k] * ksg_194[k];

        t_273[k] = pa_y[k] * ish0_189[k]
                   - f_8 * pc_y[k] * ish1_189[k];

        t_274[k] = f_9 * isg_135[k]
                   + f_3 * pc_y[k] * ksg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, ish0_192, ish0_194, \
                         isg_120, isg_136, isg_137, ish1_192, ish1_194, ksg_195, \
                         ksg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * isg_120[k]
                   + f_3 * pc_z[k] * ksg_195[k];

        t_276[k] = pa_y[k] * ish0_192[k]
                   + f_10 * isg_136[k]
                   - f_8 * pc_y[k] * ish1_192[k];

        t_277[k] = f_9 * isg_137[k]
                   + f_3 * pc_y[k] * ksg_197[k];

        t_278[k] = pa_y[k] * ish0_194[k]
                   - f_8 * pc_y[k] * ish1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, ish0_195, ish0_198, \
                         isg_123, isg_138, isg_140, ish1_195, ish1_198, ksg_198, \
                         ksg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * ish0_195[k]
                   + f_11 * isg_138[k]
                   - f_8 * pc_y[k] * ish1_195[k];

        t_280[k] = f_11 * isg_123[k]
                   + f_3 * pc_z[k] * ksg_198[k];

        t_281[k] = f_9 * isg_140[k]
                   + f_3 * pc_y[k] * ksg_200[k];

        t_282[k] = pa_y[k] * ish0_198[k]
                   - f_8 * pc_y[k] * ish1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, isg_205, isg_206, isg_207, \
                         isg_208, isg_209, ksg_205, ksg_206, ksg_207, ksg_208, \
                         ksg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_11 * isg_205[k]
                   + f_3 * pc_x[k] * ksg_205[k];

        t_284[k] = f_11 * isg_206[k]
                   + f_3 * pc_x[k] * ksg_206[k];

        t_285[k] = f_11 * isg_207[k]
                   + f_3 * pc_x[k] * ksg_207[k];

        t_286[k] = f_11 * isg_208[k]
                   + f_3 * pc_x[k] * ksg_208[k];

        t_287[k] = f_11 * isg_209[k]
                   + f_3 * pc_x[k] * ksg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, isg_130, isg_145, isg_147, ksf0_136, \
                         ksf0_138, ksf1_136, ksf1_138, ksg_205, \
                         ksg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * isg_145[k]
                   + f_1 * ksf0_136[k]
                   - f_2 * ksf1_136[k]
                   + f_3 * pc_y[k] * ksg_205[k];

        t_289[k] = f_11 * isg_130[k]
                   + f_3 * pc_z[k] * ksg_205[k];

        t_290[k] = f_9 * isg_147[k]
                   + f_6 * ksf0_138[k]
                   - f_7 * ksf1_138[k]
                   + f_3 * pc_y[k] * ksg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, ish0_209, isg_148, isg_149, \
                         ish1_209, ksf0_139, ksf1_139, ksg_208, \
                         ksg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * isg_148[k]
                   + f_4 * ksf0_139[k]
                   - f_5 * ksf1_139[k]
                   + f_3 * pc_y[k] * ksg_208[k];

        t_292[k] = f_9 * isg_149[k]
                   + f_3 * pc_y[k] * ksg_209[k];

        t_293[k] = pa_y[k] * ish0_209[k]
                   - f_8 * pc_y[k] * ish1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, isg_135, \
                         isg_210, ksf0_140, ksf1_140, ksg_210, ksg_211, \
                         ksg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_11 * isg_210[k]
                   + f_1 * ksf0_140[k]
                   - f_2 * ksf1_140[k]
                   + f_3 * pc_x[k] * ksg_210[k];

        t_295[k] = f_3 * pc_y[k] * ksg_210[k];

        t_296[k] = f_16 * isg_135[k]
                   + f_3 * pc_z[k] * ksg_210[k];

        t_297[k] = f_4 * ksf0_140[k]
                   - f_5 * ksf1_140[k]
                   + f_3 * pc_y[k] * ksg_211[k];

        t_298[k] = f_3 * pc_y[k] * ksg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, isg_215, ksf0_141, ksf0_142, \
                         ksf0_145, ksf1_141, ksf1_142, ksf1_145, ksg_213, ksg_214, \
                         ksg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_11 * isg_215[k]
                   + f_6 * ksf0_145[k]
                   - f_7 * ksf1_145[k]
                   + f_3 * pc_x[k] * ksg_215[k];

        t_300[k] = f_6 * ksf0_141[k]
                   - f_7 * ksf1_141[k]
                   + f_3 * pc_y[k] * ksg_213[k];

        t_301[k] = f_4 * ksf0_142[k]
                   - f_5 * ksf1_142[k]
                   + f_3 * pc_y[k] * ksg_214[k];

        t_302[k] = f_3 * pc_y[k] * ksg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, isg_219, isg_220, isg_221, isg_222, \
                         ksf0_149, ksf1_149, ksg_219, ksg_220, ksg_221, \
                         ksg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_11 * isg_219[k]
                   + f_4 * ksf0_149[k]
                   - f_5 * ksf1_149[k]
                   + f_3 * pc_x[k] * ksg_219[k];

        t_304[k] = f_11 * isg_220[k]
                   + f_3 * pc_x[k] * ksg_220[k];

        t_305[k] = f_11 * isg_221[k]
                   + f_3 * pc_x[k] * ksg_221[k];

        t_306[k] = f_11 * isg_222[k]
                   + f_3 * pc_x[k] * ksg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, isg_224, ksf0_146, ksf0_147, \
                         ksf1_146, ksf1_147, ksg_219, ksg_220, ksg_221, \
                         ksg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * ksg_219[k];

        t_308[k] = f_11 * isg_224[k]
                   + f_3 * pc_x[k] * ksg_224[k];

        t_309[k] = f_1 * ksf0_146[k]
                   - f_2 * ksf1_146[k]
                   + f_3 * pc_y[k] * ksg_220[k];

        t_310[k] = f_13 * ksf0_147[k]
                   - f_14 * ksf1_147[k]
                   + f_3 * pc_y[k] * ksg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, isg_149, ksf0_148, ksf0_149, \
                         ksf1_148, ksf1_149, ksg_222, ksg_223, \
                         ksg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * ksf0_148[k]
                   - f_7 * ksf1_148[k]
                   + f_3 * pc_y[k] * ksg_222[k];

        t_312[k] = f_4 * ksf0_149[k]
                   - f_5 * ksf1_149[k]
                   + f_3 * pc_y[k] * ksg_223[k];

        t_313[k] = f_3 * pc_y[k] * ksg_224[k];

        t_314[k] = f_16 * isg_149[k]
                   + f_1 * ksf0_149[k]
                   - f_2 * ksf1_149[k]
                   + f_3 * pc_z[k] * ksg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, isg_150, isg_225, \
                         isg_228, ksf0_150, ksf0_153, ksf1_150, ksf1_153, ksg_225, \
                         ksg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_10 * isg_225[k]
                   + f_1 * ksf0_150[k]
                   - f_2 * ksf1_150[k]
                   + f_3 * pc_x[k] * ksg_225[k];

        t_316[k] = f_15 * isg_150[k]
                   + f_3 * pc_y[k] * ksg_225[k];

        t_317[k] = f_3 * pc_z[k] * ksg_225[k];

        t_318[k] = f_10 * isg_228[k]
                   + f_6 * ksf0_153[k]
                   - f_7 * ksf1_153[k]
                   + f_3 * pc_x[k] * ksg_228[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pc_x, pc_z, isg_231, ksf0_150, ksf0_156, \
                         ksf1_150, ksf1_156, ksg_226, ksg_227, ksg_228, \
                         ksg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * ksg_226[k];

        t_320[k] = f_4 * ksf0_150[k]
                   - f_5 * ksf1_150[k]
                   + f_3 * pc_z[k] * ksg_227[k];

        t_321[k] = f_10 * isg_231[k]
                   + f_4 * ksf0_156[k]
                   - f_5 * ksf1_156[k]
                   + f_3 * pc_x[k] * ksg_231[k];

        t_322[k] = f_3 * pc_z[k] * ksg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, isg_155, isg_235, \
                         ksf0_152, ksf1_152, ksg_230, ksg_231, \
                         ksg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_15 * isg_155[k]
                   + f_3 * pc_y[k] * ksg_230[k];

        t_324[k] = f_6 * ksf0_152[k]
                   - f_7 * ksf1_152[k]
                   + f_3 * pc_z[k] * ksg_230[k];

        t_325[k] = f_10 * isg_235[k]
                   + f_3 * pc_x[k] * ksg_235[k];

        t_326[k] = f_3 * pc_z[k] * ksg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pc_x, pc_y, isg_160, isg_237, isg_238, \
                         isg_239, ksf0_156, ksf1_156, ksg_235, ksg_237, ksg_238, \
                         ksg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_10 * isg_237[k]
                   + f_3 * pc_x[k] * ksg_237[k];

        t_328[k] = f_10 * isg_238[k]
                   + f_3 * pc_x[k] * ksg_238[k];

        t_329[k] = f_10 * isg_239[k]
                   + f_3 * pc_x[k] * ksg_239[k];

        t_330[k] = f_15 * isg_160[k]
                   + f_1 * ksf0_156[k]
                   - f_2 * ksf1_156[k]
                   + f_3 * pc_y[k] * ksg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, isg_164, ksf0_156, ksf0_157, \
                         ksf1_156, ksf1_157, ksg_235, ksg_236, ksg_237, \
                         ksg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * ksg_235[k];

        t_332[k] = f_4 * ksf0_156[k]
                   - f_5 * ksf1_156[k]
                   + f_3 * pc_z[k] * ksg_236[k];

        t_333[k] = f_6 * ksf0_157[k]
                   - f_7 * ksf1_157[k]
                   + f_3 * pc_z[k] * ksg_237[k];

        t_334[k] = f_15 * isg_164[k]
                   + f_3 * pc_y[k] * ksg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_z, pc_y, pc_z, ish0_210, isg_150, \
                         isg_165, ish1_210, ksf0_159, ksf1_159, ksg_239, \
                         ksg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * ksf0_159[k]
                   - f_2 * ksf1_159[k]
                   + f_3 * pc_z[k] * ksg_239[k];

        t_336[k] = pa_z[k] * ish0_210[k]
                   - f_8 * pc_z[k] * ish1_210[k];

        t_337[k] = f_16 * isg_165[k]
                   + f_3 * pc_y[k] * ksg_240[k];

        t_338[k] = f_9 * isg_150[k]
                   + f_3 * pc_z[k] * ksg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_z, pc_x, pc_y, pc_z, ish0_213, isg_167, \
                         isg_245, ish1_213, ksf0_165, ksf1_165, ksg_242, \
                         ksg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * ish0_213[k]
                   - f_8 * pc_z[k] * ish1_213[k];

        t_340[k] = f_16 * isg_167[k]
                   + f_3 * pc_y[k] * ksg_242[k];

        t_341[k] = f_10 * isg_245[k]
                   + f_6 * ksf0_165[k]
                   - f_7 * ksf1_165[k]
                   + f_3 * pc_x[k] * ksg_245[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, ish0_216, isg_153, isg_170, \
                         ish1_216, ksg_243, ksg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * ish0_216[k]
                   - f_8 * pc_z[k] * ish1_216[k];

        t_343[k] = f_9 * isg_153[k]
                   + f_3 * pc_z[k] * ksg_243[k];

        t_344[k] = f_16 * isg_170[k]
                   + f_3 * pc_y[k] * ksg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pc_x, isg_249, isg_250, isg_251, isg_252, \
                         ksf0_169, ksf1_169, ksg_249, ksg_250, ksg_251, \
                         ksg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_10 * isg_249[k]
                   + f_4 * ksf0_169[k]
                   - f_5 * ksf1_169[k]
                   + f_3 * pc_x[k] * ksg_249[k];

        t_346[k] = f_10 * isg_250[k]
                   + f_3 * pc_x[k] * ksg_250[k];

        t_347[k] = f_10 * isg_251[k]
                   + f_3 * pc_x[k] * ksg_251[k];

        t_348[k] = f_10 * isg_252[k]
                   + f_3 * pc_x[k] * ksg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, ish0_225, isg_160, \
                         isg_253, isg_254, ish1_225, ksg_250, ksg_253, \
                         ksg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_10 * isg_253[k]
                   + f_3 * pc_x[k] * ksg_253[k];

        t_350[k] = f_10 * isg_254[k]
                   + f_3 * pc_x[k] * ksg_254[k];

        t_351[k] = pa_z[k] * ish0_225[k]
                   - f_8 * pc_z[k] * ish1_225[k];

        t_352[k] = f_9 * isg_160[k]
                   + f_3 * pc_z[k] * ksg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, isg_177, isg_178, isg_179, ksf0_168, \
                         ksf0_169, ksf1_168, ksf1_169, ksg_252, ksg_253, \
                         ksg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * isg_177[k]
                   + f_6 * ksf0_168[k]
                   - f_7 * ksf1_168[k]
                   + f_3 * pc_y[k] * ksg_252[k];

        t_354[k] = f_16 * isg_178[k]
                   + f_4 * ksf0_169[k]
                   - f_5 * ksf1_169[k]
                   + f_3 * pc_y[k] * ksg_253[k];

        t_355[k] = f_16 * isg_179[k]
                   + f_3 * pc_y[k] * ksg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, isg_164, isg_180, isg_255, \
                         ksf0_169, ksf0_170, ksf1_169, ksf1_170, ksg_254, \
                         ksg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * isg_164[k]
                   + f_1 * ksf0_169[k]
                   - f_2 * ksf1_169[k]
                   + f_3 * pc_z[k] * ksg_254[k];

        t_357[k] = f_10 * isg_255[k]
                   + f_1 * ksf0_170[k]
                   - f_2 * ksf1_170[k]
                   + f_3 * pc_x[k] * ksg_255[k];

        t_358[k] = f_11 * isg_180[k]
                   + f_3 * pc_y[k] * ksg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, isg_165, isg_182, isg_258, \
                         ksf0_173, ksf1_173, ksg_255, ksg_257, \
                         ksg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * isg_165[k]
                   + f_3 * pc_z[k] * ksg_255[k];

        t_360[k] = f_10 * isg_258[k]
                   + f_6 * ksf0_173[k]
                   - f_7 * ksf1_173[k]
                   + f_3 * pc_x[k] * ksg_258[k];

        t_361[k] = f_11 * isg_182[k]
                   + f_3 * pc_y[k] * ksg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, pc_z, isg_168, isg_260, isg_261, ksf0_175, \
                         ksf0_176, ksf1_175, ksf1_176, ksg_258, ksg_260, \
                         ksg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_10 * isg_260[k]
                   + f_6 * ksf0_175[k]
                   - f_7 * ksf1_175[k]
                   + f_3 * pc_x[k] * ksg_260[k];

        t_363[k] = f_10 * isg_261[k]
                   + f_4 * ksf0_176[k]
                   - f_5 * ksf1_176[k]
                   + f_3 * pc_x[k] * ksg_261[k];

        t_364[k] = f_10 * isg_168[k]
                   + f_3 * pc_z[k] * ksg_258[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pc_x, pc_y, isg_185, isg_264, isg_265, \
                         isg_266, ksf0_179, ksf1_179, ksg_260, ksg_264, ksg_265, \
                         ksg_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_11 * isg_185[k]
                   + f_3 * pc_y[k] * ksg_260[k];

        t_366[k] = f_10 * isg_264[k]
                   + f_4 * ksf0_179[k]
                   - f_5 * ksf1_179[k]
                   + f_3 * pc_x[k] * ksg_264[k];

        t_367[k] = f_10 * isg_265[k]
                   + f_3 * pc_x[k] * ksg_265[k];

        t_368[k] = f_10 * isg_266[k]
                   + f_3 * pc_x[k] * ksg_266[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
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
    const auto f_12 = 3.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_488 = buffer.data(target + 488);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ish0_294 = buffer.data(ish0 + 294);
    const auto *ish0_297 = buffer.data(ish0 + 297);
    const auto *ish0_299 = buffer.data(ish0 + 299);
    const auto *ish0_300 = buffer.data(ish0 + 300);
    const auto *ish0_303 = buffer.data(ish0 + 303);
    const auto *ish0_314 = buffer.data(ish0 + 314);
    const auto *ish0_315 = buffer.data(ish0 + 315);
    const auto *ish0_318 = buffer.data(ish0 + 318);
    const auto *ish0_321 = buffer.data(ish0 + 321);
    const auto *ish0_441 = buffer.data(ish0 + 441);
    const auto *ish0_444 = buffer.data(ish0 + 444);
    const auto *ish0_447 = buffer.data(ish0 + 447);
    const auto *ish0_456 = buffer.data(ish0 + 456);
    const auto *ish0_458 = buffer.data(ish0 + 458);
    const auto *ish0_459 = buffer.data(ish0 + 459);
    const auto *ish0_461 = buffer.data(ish0 + 461);
    const auto *ish0_467 = buffer.data(ish0 + 467);
    const auto *ish0_471 = buffer.data(ish0 + 471);
    const auto *ish0_477 = buffer.data(ish0 + 477);
    const auto *ish0_479 = buffer.data(ish0 + 479);
    const auto *ish0_480 = buffer.data(ish0 + 480);
    const auto *ish0_482 = buffer.data(ish0 + 482);
    const auto *ish0_483 = buffer.data(ish0 + 483);
    const auto *ish0_486 = buffer.data(ish0 + 486);
    const auto *ish0_488 = buffer.data(ish0 + 488);

    const auto *isg_175 = buffer.data(isg + 175);
    const auto *isg_179 = buffer.data(isg + 179);
    const auto *isg_180 = buffer.data(isg + 180);
    const auto *isg_183 = buffer.data(isg + 183);
    const auto *isg_190 = buffer.data(isg + 190);
    const auto *isg_192 = buffer.data(isg + 192);
    const auto *isg_193 = buffer.data(isg + 193);
    const auto *isg_194 = buffer.data(isg + 194);
    const auto *isg_195 = buffer.data(isg + 195);
    const auto *isg_197 = buffer.data(isg + 197);
    const auto *isg_198 = buffer.data(isg + 198);
    const auto *isg_200 = buffer.data(isg + 200);
    const auto *isg_205 = buffer.data(isg + 205);
    const auto *isg_207 = buffer.data(isg + 207);
    const auto *isg_208 = buffer.data(isg + 208);
    const auto *isg_209 = buffer.data(isg + 209);
    const auto *isg_210 = buffer.data(isg + 210);
    const auto *isg_211 = buffer.data(isg + 211);
    const auto *isg_212 = buffer.data(isg + 212);
    const auto *isg_213 = buffer.data(isg + 213);
    const auto *isg_215 = buffer.data(isg + 215);
    const auto *isg_220 = buffer.data(isg + 220);
    const auto *isg_222 = buffer.data(isg + 222);
    const auto *isg_223 = buffer.data(isg + 223);
    const auto *isg_224 = buffer.data(isg + 224);
    const auto *isg_225 = buffer.data(isg + 225);
    const auto *isg_228 = buffer.data(isg + 228);
    const auto *isg_230 = buffer.data(isg + 230);
    const auto *isg_235 = buffer.data(isg + 235);
    const auto *isg_239 = buffer.data(isg + 239);
    const auto *isg_240 = buffer.data(isg + 240);
    const auto *isg_242 = buffer.data(isg + 242);
    const auto *isg_245 = buffer.data(isg + 245);
    const auto *isg_254 = buffer.data(isg + 254);
    const auto *isg_255 = buffer.data(isg + 255);
    const auto *isg_257 = buffer.data(isg + 257);
    const auto *isg_267 = buffer.data(isg + 267);
    const auto *isg_268 = buffer.data(isg + 268);
    const auto *isg_269 = buffer.data(isg + 269);
    const auto *isg_270 = buffer.data(isg + 270);
    const auto *isg_273 = buffer.data(isg + 273);
    const auto *isg_275 = buffer.data(isg + 275);
    const auto *isg_276 = buffer.data(isg + 276);
    const auto *isg_279 = buffer.data(isg + 279);
    const auto *isg_280 = buffer.data(isg + 280);
    const auto *isg_281 = buffer.data(isg + 281);
    const auto *isg_282 = buffer.data(isg + 282);
    const auto *isg_283 = buffer.data(isg + 283);
    const auto *isg_284 = buffer.data(isg + 284);
    const auto *isg_295 = buffer.data(isg + 295);
    const auto *isg_296 = buffer.data(isg + 296);
    const auto *isg_297 = buffer.data(isg + 297);
    const auto *isg_298 = buffer.data(isg + 298);
    const auto *isg_299 = buffer.data(isg + 299);
    const auto *isg_300 = buffer.data(isg + 300);
    const auto *isg_305 = buffer.data(isg + 305);
    const auto *isg_309 = buffer.data(isg + 309);
    const auto *isg_310 = buffer.data(isg + 310);
    const auto *isg_311 = buffer.data(isg + 311);
    const auto *isg_312 = buffer.data(isg + 312);
    const auto *isg_314 = buffer.data(isg + 314);
    const auto *isg_315 = buffer.data(isg + 315);
    const auto *isg_318 = buffer.data(isg + 318);
    const auto *isg_321 = buffer.data(isg + 321);
    const auto *isg_325 = buffer.data(isg + 325);
    const auto *isg_327 = buffer.data(isg + 327);
    const auto *isg_328 = buffer.data(isg + 328);
    const auto *isg_329 = buffer.data(isg + 329);
    const auto *isg_335 = buffer.data(isg + 335);
    const auto *isg_339 = buffer.data(isg + 339);
    const auto *isg_340 = buffer.data(isg + 340);
    const auto *isg_341 = buffer.data(isg + 341);
    const auto *isg_342 = buffer.data(isg + 342);
    const auto *isg_343 = buffer.data(isg + 343);
    const auto *isg_344 = buffer.data(isg + 344);
    const auto *isg_345 = buffer.data(isg + 345);
    const auto *isg_348 = buffer.data(isg + 348);
    const auto *isg_350 = buffer.data(isg + 350);

    const auto *ish1_294 = buffer.data(ish1 + 294);
    const auto *ish1_297 = buffer.data(ish1 + 297);
    const auto *ish1_299 = buffer.data(ish1 + 299);
    const auto *ish1_300 = buffer.data(ish1 + 300);
    const auto *ish1_303 = buffer.data(ish1 + 303);
    const auto *ish1_314 = buffer.data(ish1 + 314);
    const auto *ish1_315 = buffer.data(ish1 + 315);
    const auto *ish1_318 = buffer.data(ish1 + 318);
    const auto *ish1_321 = buffer.data(ish1 + 321);
    const auto *ish1_441 = buffer.data(ish1 + 441);
    const auto *ish1_444 = buffer.data(ish1 + 444);
    const auto *ish1_447 = buffer.data(ish1 + 447);
    const auto *ish1_456 = buffer.data(ish1 + 456);
    const auto *ish1_458 = buffer.data(ish1 + 458);
    const auto *ish1_459 = buffer.data(ish1 + 459);
    const auto *ish1_461 = buffer.data(ish1 + 461);
    const auto *ish1_467 = buffer.data(ish1 + 467);
    const auto *ish1_471 = buffer.data(ish1 + 471);
    const auto *ish1_477 = buffer.data(ish1 + 477);
    const auto *ish1_479 = buffer.data(ish1 + 479);
    const auto *ish1_480 = buffer.data(ish1 + 480);
    const auto *ish1_482 = buffer.data(ish1 + 482);
    const auto *ish1_483 = buffer.data(ish1 + 483);
    const auto *ish1_486 = buffer.data(ish1 + 486);
    const auto *ish1_488 = buffer.data(ish1 + 488);

    const auto *ksf0_176 = buffer.data(ksf0 + 176);
    const auto *ksf0_178 = buffer.data(ksf0 + 178);
    const auto *ksf0_179 = buffer.data(ksf0 + 179);
    const auto *ksf0_180 = buffer.data(ksf0 + 180);
    const auto *ksf0_183 = buffer.data(ksf0 + 183);
    const auto *ksf0_185 = buffer.data(ksf0 + 185);
    const auto *ksf0_186 = buffer.data(ksf0 + 186);
    const auto *ksf0_188 = buffer.data(ksf0 + 188);
    const auto *ksf0_189 = buffer.data(ksf0 + 189);
    const auto *ksf0_196 = buffer.data(ksf0 + 196);
    const auto *ksf0_198 = buffer.data(ksf0 + 198);
    const auto *ksf0_199 = buffer.data(ksf0 + 199);
    const auto *ksf0_200 = buffer.data(ksf0 + 200);
    const auto *ksf0_201 = buffer.data(ksf0 + 201);
    const auto *ksf0_202 = buffer.data(ksf0 + 202);
    const auto *ksf0_205 = buffer.data(ksf0 + 205);
    const auto *ksf0_206 = buffer.data(ksf0 + 206);
    const auto *ksf0_207 = buffer.data(ksf0 + 207);
    const auto *ksf0_208 = buffer.data(ksf0 + 208);
    const auto *ksf0_209 = buffer.data(ksf0 + 209);
    const auto *ksf0_210 = buffer.data(ksf0 + 210);
    const auto *ksf0_212 = buffer.data(ksf0 + 212);

    const auto *ksf1_176 = buffer.data(ksf1 + 176);
    const auto *ksf1_178 = buffer.data(ksf1 + 178);
    const auto *ksf1_179 = buffer.data(ksf1 + 179);
    const auto *ksf1_180 = buffer.data(ksf1 + 180);
    const auto *ksf1_183 = buffer.data(ksf1 + 183);
    const auto *ksf1_185 = buffer.data(ksf1 + 185);
    const auto *ksf1_186 = buffer.data(ksf1 + 186);
    const auto *ksf1_188 = buffer.data(ksf1 + 188);
    const auto *ksf1_189 = buffer.data(ksf1 + 189);
    const auto *ksf1_196 = buffer.data(ksf1 + 196);
    const auto *ksf1_198 = buffer.data(ksf1 + 198);
    const auto *ksf1_199 = buffer.data(ksf1 + 199);
    const auto *ksf1_200 = buffer.data(ksf1 + 200);
    const auto *ksf1_201 = buffer.data(ksf1 + 201);
    const auto *ksf1_202 = buffer.data(ksf1 + 202);
    const auto *ksf1_205 = buffer.data(ksf1 + 205);
    const auto *ksf1_206 = buffer.data(ksf1 + 206);
    const auto *ksf1_207 = buffer.data(ksf1 + 207);
    const auto *ksf1_208 = buffer.data(ksf1 + 208);
    const auto *ksf1_209 = buffer.data(ksf1 + 209);
    const auto *ksf1_210 = buffer.data(ksf1 + 210);
    const auto *ksf1_212 = buffer.data(ksf1 + 212);

    const auto *ksg_265 = buffer.data(ksg + 265);
    const auto *ksg_267 = buffer.data(ksg + 267);
    const auto *ksg_268 = buffer.data(ksg + 268);
    const auto *ksg_269 = buffer.data(ksg + 269);
    const auto *ksg_270 = buffer.data(ksg + 270);
    const auto *ksg_272 = buffer.data(ksg + 272);
    const auto *ksg_273 = buffer.data(ksg + 273);
    const auto *ksg_275 = buffer.data(ksg + 275);
    const auto *ksg_276 = buffer.data(ksg + 276);
    const auto *ksg_279 = buffer.data(ksg + 279);
    const auto *ksg_280 = buffer.data(ksg + 280);
    const auto *ksg_281 = buffer.data(ksg + 281);
    const auto *ksg_282 = buffer.data(ksg + 282);
    const auto *ksg_283 = buffer.data(ksg + 283);
    const auto *ksg_284 = buffer.data(ksg + 284);
    const auto *ksg_285 = buffer.data(ksg + 285);
    const auto *ksg_287 = buffer.data(ksg + 287);
    const auto *ksg_288 = buffer.data(ksg + 288);
    const auto *ksg_290 = buffer.data(ksg + 290);
    const auto *ksg_295 = buffer.data(ksg + 295);
    const auto *ksg_296 = buffer.data(ksg + 296);
    const auto *ksg_297 = buffer.data(ksg + 297);
    const auto *ksg_298 = buffer.data(ksg + 298);
    const auto *ksg_299 = buffer.data(ksg + 299);
    const auto *ksg_300 = buffer.data(ksg + 300);
    const auto *ksg_301 = buffer.data(ksg + 301);
    const auto *ksg_302 = buffer.data(ksg + 302);
    const auto *ksg_303 = buffer.data(ksg + 303);
    const auto *ksg_304 = buffer.data(ksg + 304);
    const auto *ksg_305 = buffer.data(ksg + 305);
    const auto *ksg_309 = buffer.data(ksg + 309);
    const auto *ksg_310 = buffer.data(ksg + 310);
    const auto *ksg_311 = buffer.data(ksg + 311);
    const auto *ksg_312 = buffer.data(ksg + 312);
    const auto *ksg_313 = buffer.data(ksg + 313);
    const auto *ksg_314 = buffer.data(ksg + 314);
    const auto *ksg_315 = buffer.data(ksg + 315);
    const auto *ksg_316 = buffer.data(ksg + 316);
    const auto *ksg_317 = buffer.data(ksg + 317);
    const auto *ksg_318 = buffer.data(ksg + 318);
    const auto *ksg_320 = buffer.data(ksg + 320);
    const auto *ksg_321 = buffer.data(ksg + 321);
    const auto *ksg_325 = buffer.data(ksg + 325);
    const auto *ksg_327 = buffer.data(ksg + 327);
    const auto *ksg_328 = buffer.data(ksg + 328);
    const auto *ksg_329 = buffer.data(ksg + 329);
    const auto *ksg_330 = buffer.data(ksg + 330);
    const auto *ksg_332 = buffer.data(ksg + 332);
    const auto *ksg_333 = buffer.data(ksg + 333);
    const auto *ksg_335 = buffer.data(ksg + 335);
    const auto *ksg_340 = buffer.data(ksg + 340);
    const auto *ksg_341 = buffer.data(ksg + 341);
    const auto *ksg_342 = buffer.data(ksg + 342);
    const auto *ksg_343 = buffer.data(ksg + 343);
    const auto *ksg_344 = buffer.data(ksg + 344);
    const auto *ksg_345 = buffer.data(ksg + 345);
    const auto *ksg_347 = buffer.data(ksg + 347);

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pc_x, pc_y, isg_190, isg_267, isg_268, \
                         isg_269, ksf0_176, ksf1_176, ksg_265, ksg_267, ksg_268, \
                         ksg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_10 * isg_267[k]
                   + f_3 * pc_x[k] * ksg_267[k];

        t_370[k] = f_10 * isg_268[k]
                   + f_3 * pc_x[k] * ksg_268[k];

        t_371[k] = f_10 * isg_269[k]
                   + f_3 * pc_x[k] * ksg_269[k];

        t_372[k] = f_11 * isg_190[k]
                   + f_1 * ksf0_176[k]
                   - f_2 * ksf1_176[k]
                   + f_3 * pc_y[k] * ksg_265[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_y, pc_z, isg_175, isg_192, isg_193, ksf0_178, \
                         ksf0_179, ksf1_178, ksf1_179, ksg_265, ksg_267, \
                         ksg_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_10 * isg_175[k]
                   + f_3 * pc_z[k] * ksg_265[k];

        t_374[k] = f_11 * isg_192[k]
                   + f_6 * ksf0_178[k]
                   - f_7 * ksf1_178[k]
                   + f_3 * pc_y[k] * ksg_267[k];

        t_375[k] = f_11 * isg_193[k]
                   + f_4 * ksf0_179[k]
                   - f_5 * ksf1_179[k]
                   + f_3 * pc_y[k] * ksg_268[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pc_x, pc_y, pc_z, isg_179, isg_194, isg_270, \
                         ksf0_179, ksf0_180, ksf1_179, ksf1_180, ksg_269, \
                         ksg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_11 * isg_194[k]
                   + f_3 * pc_y[k] * ksg_269[k];

        t_377[k] = f_10 * isg_179[k]
                   + f_1 * ksf0_179[k]
                   - f_2 * ksf1_179[k]
                   + f_3 * pc_z[k] * ksg_269[k];

        t_378[k] = f_10 * isg_270[k]
                   + f_1 * ksf0_180[k]
                   - f_2 * ksf1_180[k]
                   + f_3 * pc_x[k] * ksg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, t_382, pc_x, pc_y, pc_z, isg_180, isg_195, \
                         isg_197, isg_273, ksf0_183, ksf1_183, ksg_270, ksg_272, \
                         ksg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_10 * isg_195[k]
                   + f_3 * pc_y[k] * ksg_270[k];

        t_380[k] = f_11 * isg_180[k]
                   + f_3 * pc_z[k] * ksg_270[k];

        t_381[k] = f_10 * isg_273[k]
                   + f_6 * ksf0_183[k]
                   - f_7 * ksf1_183[k]
                   + f_3 * pc_x[k] * ksg_273[k];

        t_382[k] = f_10 * isg_197[k]
                   + f_3 * pc_y[k] * ksg_272[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, pc_z, isg_183, isg_275, isg_276, ksf0_185, \
                         ksf0_186, ksf1_185, ksf1_186, ksg_273, ksg_275, \
                         ksg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_10 * isg_275[k]
                   + f_6 * ksf0_185[k]
                   - f_7 * ksf1_185[k]
                   + f_3 * pc_x[k] * ksg_275[k];

        t_384[k] = f_10 * isg_276[k]
                   + f_4 * ksf0_186[k]
                   - f_5 * ksf1_186[k]
                   + f_3 * pc_x[k] * ksg_276[k];

        t_385[k] = f_11 * isg_183[k]
                   + f_3 * pc_z[k] * ksg_273[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pc_x, pc_y, isg_200, isg_279, isg_280, \
                         isg_281, ksf0_189, ksf1_189, ksg_275, ksg_279, ksg_280, \
                         ksg_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_10 * isg_200[k]
                   + f_3 * pc_y[k] * ksg_275[k];

        t_387[k] = f_10 * isg_279[k]
                   + f_4 * ksf0_189[k]
                   - f_5 * ksf1_189[k]
                   + f_3 * pc_x[k] * ksg_279[k];

        t_388[k] = f_10 * isg_280[k]
                   + f_3 * pc_x[k] * ksg_280[k];

        t_389[k] = f_10 * isg_281[k]
                   + f_3 * pc_x[k] * ksg_281[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, pc_y, isg_205, isg_282, isg_283, \
                         isg_284, ksf0_186, ksf1_186, ksg_280, ksg_282, ksg_283, \
                         ksg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_10 * isg_282[k]
                   + f_3 * pc_x[k] * ksg_282[k];

        t_391[k] = f_10 * isg_283[k]
                   + f_3 * pc_x[k] * ksg_283[k];

        t_392[k] = f_10 * isg_284[k]
                   + f_3 * pc_x[k] * ksg_284[k];

        t_393[k] = f_10 * isg_205[k]
                   + f_1 * ksf0_186[k]
                   - f_2 * ksf1_186[k]
                   + f_3 * pc_y[k] * ksg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, isg_190, isg_207, isg_208, ksf0_188, \
                         ksf0_189, ksf1_188, ksf1_189, ksg_280, ksg_282, \
                         ksg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * isg_190[k]
                   + f_3 * pc_z[k] * ksg_280[k];

        t_395[k] = f_10 * isg_207[k]
                   + f_6 * ksf0_188[k]
                   - f_7 * ksf1_188[k]
                   + f_3 * pc_y[k] * ksg_282[k];

        t_396[k] = f_10 * isg_208[k]
                   + f_4 * ksf0_189[k]
                   - f_5 * ksf1_189[k]
                   + f_3 * pc_y[k] * ksg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pc_y, pc_z, ish0_294, isg_194, \
                         isg_209, isg_210, ish1_294, ksf0_189, ksf1_189, ksg_284, \
                         ksg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * isg_209[k]
                   + f_3 * pc_y[k] * ksg_284[k];

        t_398[k] = f_11 * isg_194[k]
                   + f_1 * ksf0_189[k]
                   - f_2 * ksf1_189[k]
                   + f_3 * pc_z[k] * ksg_284[k];

        t_399[k] = pa_y[k] * ish0_294[k]
                   - f_8 * pc_y[k] * ish1_294[k];

        t_400[k] = f_9 * isg_210[k]
                   + f_3 * pc_y[k] * ksg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_y, pc_y, pc_z, ish0_297, ish0_299, \
                         isg_195, isg_211, isg_212, ish1_297, ish1_299, ksg_285, \
                         ksg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_16 * isg_195[k]
                   + f_3 * pc_z[k] * ksg_285[k];

        t_402[k] = pa_y[k] * ish0_297[k]
                   + f_10 * isg_211[k]
                   - f_8 * pc_y[k] * ish1_297[k];

        t_403[k] = f_9 * isg_212[k]
                   + f_3 * pc_y[k] * ksg_287[k];

        t_404[k] = pa_y[k] * ish0_299[k]
                   - f_8 * pc_y[k] * ish1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pa_y, pc_y, pc_z, ish0_300, ish0_303, \
                         isg_198, isg_213, isg_215, ish1_300, ish1_303, ksg_288, \
                         ksg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_y[k] * ish0_300[k]
                   + f_11 * isg_213[k]
                   - f_8 * pc_y[k] * ish1_300[k];

        t_406[k] = f_16 * isg_198[k]
                   + f_3 * pc_z[k] * ksg_288[k];

        t_407[k] = f_9 * isg_215[k]
                   + f_3 * pc_y[k] * ksg_290[k];

        t_408[k] = pa_y[k] * ish0_303[k]
                   - f_8 * pc_y[k] * ish1_303[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, isg_295, isg_296, isg_297, \
                         isg_298, isg_299, ksg_295, ksg_296, ksg_297, ksg_298, \
                         ksg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_10 * isg_295[k]
                   + f_3 * pc_x[k] * ksg_295[k];

        t_410[k] = f_10 * isg_296[k]
                   + f_3 * pc_x[k] * ksg_296[k];

        t_411[k] = f_10 * isg_297[k]
                   + f_3 * pc_x[k] * ksg_297[k];

        t_412[k] = f_10 * isg_298[k]
                   + f_3 * pc_x[k] * ksg_298[k];

        t_413[k] = f_10 * isg_299[k]
                   + f_3 * pc_x[k] * ksg_299[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_y, pc_z, isg_205, isg_220, isg_222, ksf0_196, \
                         ksf0_198, ksf1_196, ksf1_198, ksg_295, \
                         ksg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_9 * isg_220[k]
                   + f_1 * ksf0_196[k]
                   - f_2 * ksf1_196[k]
                   + f_3 * pc_y[k] * ksg_295[k];

        t_415[k] = f_16 * isg_205[k]
                   + f_3 * pc_z[k] * ksg_295[k];

        t_416[k] = f_9 * isg_222[k]
                   + f_6 * ksf0_198[k]
                   - f_7 * ksf1_198[k]
                   + f_3 * pc_y[k] * ksg_297[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, ish0_314, isg_223, isg_224, \
                         ish1_314, ksf0_199, ksf1_199, ksg_298, \
                         ksg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_9 * isg_223[k]
                   + f_4 * ksf0_199[k]
                   - f_5 * ksf1_199[k]
                   + f_3 * pc_y[k] * ksg_298[k];

        t_418[k] = f_9 * isg_224[k]
                   + f_3 * pc_y[k] * ksg_299[k];

        t_419[k] = pa_y[k] * ish0_314[k]
                   - f_8 * pc_y[k] * ish1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, isg_210, \
                         isg_300, ksf0_200, ksf1_200, ksg_300, ksg_301, \
                         ksg_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_10 * isg_300[k]
                   + f_1 * ksf0_200[k]
                   - f_2 * ksf1_200[k]
                   + f_3 * pc_x[k] * ksg_300[k];

        t_421[k] = f_3 * pc_y[k] * ksg_300[k];

        t_422[k] = f_15 * isg_210[k]
                   + f_3 * pc_z[k] * ksg_300[k];

        t_423[k] = f_4 * ksf0_200[k]
                   - f_5 * ksf1_200[k]
                   + f_3 * pc_y[k] * ksg_301[k];

        t_424[k] = f_3 * pc_y[k] * ksg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, isg_305, ksf0_201, ksf0_202, \
                         ksf0_205, ksf1_201, ksf1_202, ksf1_205, ksg_303, ksg_304, \
                         ksg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_10 * isg_305[k]
                   + f_6 * ksf0_205[k]
                   - f_7 * ksf1_205[k]
                   + f_3 * pc_x[k] * ksg_305[k];

        t_426[k] = f_6 * ksf0_201[k]
                   - f_7 * ksf1_201[k]
                   + f_3 * pc_y[k] * ksg_303[k];

        t_427[k] = f_4 * ksf0_202[k]
                   - f_5 * ksf1_202[k]
                   + f_3 * pc_y[k] * ksg_304[k];

        t_428[k] = f_3 * pc_y[k] * ksg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pc_x, isg_309, isg_310, isg_311, isg_312, \
                         ksf0_209, ksf1_209, ksg_309, ksg_310, ksg_311, \
                         ksg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_10 * isg_309[k]
                   + f_4 * ksf0_209[k]
                   - f_5 * ksf1_209[k]
                   + f_3 * pc_x[k] * ksg_309[k];

        t_430[k] = f_10 * isg_310[k]
                   + f_3 * pc_x[k] * ksg_310[k];

        t_431[k] = f_10 * isg_311[k]
                   + f_3 * pc_x[k] * ksg_311[k];

        t_432[k] = f_10 * isg_312[k]
                   + f_3 * pc_x[k] * ksg_312[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, isg_314, ksf0_206, ksf0_207, \
                         ksf1_206, ksf1_207, ksg_309, ksg_310, ksg_311, \
                         ksg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_3 * pc_y[k] * ksg_309[k];

        t_434[k] = f_10 * isg_314[k]
                   + f_3 * pc_x[k] * ksg_314[k];

        t_435[k] = f_1 * ksf0_206[k]
                   - f_2 * ksf1_206[k]
                   + f_3 * pc_y[k] * ksg_310[k];

        t_436[k] = f_13 * ksf0_207[k]
                   - f_14 * ksf1_207[k]
                   + f_3 * pc_y[k] * ksg_311[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_y, pc_z, isg_224, ksf0_208, ksf0_209, \
                         ksf1_208, ksf1_209, ksg_312, ksg_313, \
                         ksg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_6 * ksf0_208[k]
                   - f_7 * ksf1_208[k]
                   + f_3 * pc_y[k] * ksg_312[k];

        t_438[k] = f_4 * ksf0_209[k]
                   - f_5 * ksf1_209[k]
                   + f_3 * pc_y[k] * ksg_313[k];

        t_439[k] = f_3 * pc_y[k] * ksg_314[k];

        t_440[k] = f_15 * isg_224[k]
                   + f_1 * ksf0_209[k]
                   - f_2 * ksf1_209[k]
                   + f_3 * pc_z[k] * ksg_314[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_x, pc_x, pc_y, pc_z, ish0_441, \
                         ish0_444, isg_225, isg_315, isg_318, ish1_441, ish1_444, \
                         ksg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_x[k] * ish0_441[k]
                   + f_15 * isg_315[k]
                   - f_8 * pc_x[k] * ish1_441[k];

        t_442[k] = f_12 * isg_225[k]
                   + f_3 * pc_y[k] * ksg_315[k];

        t_443[k] = f_3 * pc_z[k] * ksg_315[k];

        t_444[k] = pa_x[k] * ish0_444[k]
                   + f_11 * isg_318[k]
                   - f_8 * pc_x[k] * ish1_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_x, pc_x, pc_z, ish0_447, isg_321, \
                         ish1_447, ksf0_210, ksf1_210, ksg_316, ksg_317, \
                         ksg_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_z[k] * ksg_316[k];

        t_446[k] = f_4 * ksf0_210[k]
                   - f_5 * ksf1_210[k]
                   + f_3 * pc_z[k] * ksg_317[k];

        t_447[k] = pa_x[k] * ish0_447[k]
                   + f_10 * isg_321[k]
                   - f_8 * pc_x[k] * ish1_447[k];

        t_448[k] = f_3 * pc_z[k] * ksg_318[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, isg_230, isg_325, \
                         ksf0_212, ksf1_212, ksg_320, ksg_321, \
                         ksg_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_12 * isg_230[k]
                   + f_3 * pc_y[k] * ksg_320[k];

        t_450[k] = f_6 * ksf0_212[k]
                   - f_7 * ksf1_212[k]
                   + f_3 * pc_z[k] * ksg_320[k];

        t_451[k] = f_9 * isg_325[k]
                   + f_3 * pc_x[k] * ksg_325[k];

        t_452[k] = f_3 * pc_z[k] * ksg_321[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pa_x, pc_x, ish0_456, isg_327, isg_328, \
                         isg_329, ish1_456, ksg_327, ksg_328, ksg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_9 * isg_327[k]
                   + f_3 * pc_x[k] * ksg_327[k];

        t_454[k] = f_9 * isg_328[k]
                   + f_3 * pc_x[k] * ksg_328[k];

        t_455[k] = f_9 * isg_329[k]
                   + f_3 * pc_x[k] * ksg_329[k];

        t_456[k] = pa_x[k] * ish0_456[k]
                   - f_8 * pc_x[k] * ish1_456[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pa_x, pc_x, pc_y, pc_z, ish0_458, \
                         ish0_459, isg_239, ish1_458, ish1_459, ksg_325, \
                         ksg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * ksg_325[k];

        t_458[k] = pa_x[k] * ish0_458[k]
                   - f_8 * pc_x[k] * ish1_458[k];

        t_459[k] = pa_x[k] * ish0_459[k]
                   - f_8 * pc_x[k] * ish1_459[k];

        t_460[k] = f_12 * isg_239[k]
                   + f_3 * pc_y[k] * ksg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_x, pa_z, pc_x, pc_y, pc_z, ish0_315, \
                         ish0_461, isg_225, isg_240, ish1_315, ish1_461, \
                         ksg_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = pa_x[k] * ish0_461[k]
                   - f_8 * pc_x[k] * ish1_461[k];

        t_462[k] = pa_z[k] * ish0_315[k]
                   - f_8 * pc_z[k] * ish1_315[k];

        t_463[k] = f_15 * isg_240[k]
                   + f_3 * pc_y[k] * ksg_330[k];

        t_464[k] = f_9 * isg_225[k]
                   + f_3 * pc_z[k] * ksg_330[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pa_x, pa_z, pc_x, pc_y, pc_z, ish0_318, \
                         ish0_467, isg_242, isg_335, ish1_318, ish1_467, \
                         ksg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = pa_z[k] * ish0_318[k]
                   - f_8 * pc_z[k] * ish1_318[k];

        t_466[k] = f_15 * isg_242[k]
                   + f_3 * pc_y[k] * ksg_332[k];

        t_467[k] = pa_x[k] * ish0_467[k]
                   + f_11 * isg_335[k]
                   - f_8 * pc_x[k] * ish1_467[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pa_z, pc_y, pc_z, ish0_321, isg_228, isg_245, \
                         ish1_321, ksg_333, ksg_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = pa_z[k] * ish0_321[k]
                   - f_8 * pc_z[k] * ish1_321[k];

        t_469[k] = f_9 * isg_228[k]
                   + f_3 * pc_z[k] * ksg_333[k];

        t_470[k] = f_15 * isg_245[k]
                   + f_3 * pc_y[k] * ksg_335[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pa_x, pc_x, ish0_471, isg_339, isg_340, \
                         isg_341, isg_342, ish1_471, ksg_340, ksg_341, \
                         ksg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pa_x[k] * ish0_471[k]
                   + f_10 * isg_339[k]
                   - f_8 * pc_x[k] * ish1_471[k];

        t_472[k] = f_9 * isg_340[k]
                   + f_3 * pc_x[k] * ksg_340[k];

        t_473[k] = f_9 * isg_341[k]
                   + f_3 * pc_x[k] * ksg_341[k];

        t_474[k] = f_9 * isg_342[k]
                   + f_3 * pc_x[k] * ksg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_x, pc_x, pc_z, ish0_477, isg_235, \
                         isg_343, isg_344, ish1_477, ksg_340, ksg_343, \
                         ksg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_9 * isg_343[k]
                   + f_3 * pc_x[k] * ksg_343[k];

        t_476[k] = f_9 * isg_344[k]
                   + f_3 * pc_x[k] * ksg_344[k];

        t_477[k] = pa_x[k] * ish0_477[k]
                   - f_8 * pc_x[k] * ish1_477[k];

        t_478[k] = f_9 * isg_235[k]
                   + f_3 * pc_z[k] * ksg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pa_x, pc_x, pc_y, ish0_479, ish0_480, \
                         ish0_482, isg_254, ish1_479, ish1_480, ish1_482, \
                         ksg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pa_x[k] * ish0_479[k]
                   - f_8 * pc_x[k] * ish1_479[k];

        t_480[k] = pa_x[k] * ish0_480[k]
                   - f_8 * pc_x[k] * ish1_480[k];

        t_481[k] = f_15 * isg_254[k]
                   + f_3 * pc_y[k] * ksg_344[k];

        t_482[k] = pa_x[k] * ish0_482[k]
                   - f_8 * pc_x[k] * ish1_482[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pa_x, pc_x, pc_y, pc_z, ish0_483, isg_240, \
                         isg_255, isg_345, ish1_483, ksg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pa_x[k] * ish0_483[k]
                   + f_15 * isg_345[k]
                   - f_8 * pc_x[k] * ish1_483[k];

        t_484[k] = f_16 * isg_255[k]
                   + f_3 * pc_y[k] * ksg_345[k];

        t_485[k] = f_10 * isg_240[k]
                   + f_3 * pc_z[k] * ksg_345[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pa_x, pc_x, pc_y, ish0_486, ish0_488, isg_257, \
                         isg_348, isg_350, ish1_486, ish1_488, \
                         ksg_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = pa_x[k] * ish0_486[k]
                   + f_11 * isg_348[k]
                   - f_8 * pc_x[k] * ish1_486[k];

        t_487[k] = f_16 * isg_257[k]
                   + f_3 * pc_y[k] * ksg_347[k];

        t_488[k] = pa_x[k] * ish0_488[k]
                   + f_11 * isg_350[k]
                   - f_8 * pc_x[k] * ish1_488[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_12 = 3.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ish0_420 = buffer.data(ish0 + 420);
    const auto *ish0_425 = buffer.data(ish0 + 425);
    const auto *ish0_429 = buffer.data(ish0 + 429);
    const auto *ish0_441 = buffer.data(ish0 + 441);
    const auto *ish0_442 = buffer.data(ish0 + 442);
    const auto *ish0_444 = buffer.data(ish0 + 444);
    const auto *ish0_489 = buffer.data(ish0 + 489);
    const auto *ish0_492 = buffer.data(ish0 + 492);
    const auto *ish0_498 = buffer.data(ish0 + 498);
    const auto *ish0_500 = buffer.data(ish0 + 500);
    const auto *ish0_501 = buffer.data(ish0 + 501);
    const auto *ish0_503 = buffer.data(ish0 + 503);
    const auto *ish0_504 = buffer.data(ish0 + 504);
    const auto *ish0_507 = buffer.data(ish0 + 507);
    const auto *ish0_509 = buffer.data(ish0 + 509);
    const auto *ish0_510 = buffer.data(ish0 + 510);
    const auto *ish0_513 = buffer.data(ish0 + 513);
    const auto *ish0_519 = buffer.data(ish0 + 519);
    const auto *ish0_521 = buffer.data(ish0 + 521);
    const auto *ish0_522 = buffer.data(ish0 + 522);
    const auto *ish0_524 = buffer.data(ish0 + 524);
    const auto *ish0_525 = buffer.data(ish0 + 525);
    const auto *ish0_528 = buffer.data(ish0 + 528);
    const auto *ish0_530 = buffer.data(ish0 + 530);
    const auto *ish0_531 = buffer.data(ish0 + 531);
    const auto *ish0_534 = buffer.data(ish0 + 534);
    const auto *ish0_540 = buffer.data(ish0 + 540);
    const auto *ish0_542 = buffer.data(ish0 + 542);
    const auto *ish0_543 = buffer.data(ish0 + 543);
    const auto *ish0_545 = buffer.data(ish0 + 545);
    const auto *ish0_549 = buffer.data(ish0 + 549);
    const auto *ish0_552 = buffer.data(ish0 + 552);
    const auto *ish0_561 = buffer.data(ish0 + 561);
    const auto *ish0_563 = buffer.data(ish0 + 563);
    const auto *ish0_564 = buffer.data(ish0 + 564);
    const auto *ish0_566 = buffer.data(ish0 + 566);
    const auto *ish0_567 = buffer.data(ish0 + 567);
    const auto *ish0_572 = buffer.data(ish0 + 572);
    const auto *ish0_576 = buffer.data(ish0 + 576);
    const auto *ish0_582 = buffer.data(ish0 + 582);
    const auto *ish0_583 = buffer.data(ish0 + 583);
    const auto *ish0_584 = buffer.data(ish0 + 584);
    const auto *ish0_585 = buffer.data(ish0 + 585);
    const auto *ish0_587 = buffer.data(ish0 + 587);

    const auto *isg_243 = buffer.data(isg + 243);
    const auto *isg_250 = buffer.data(isg + 250);
    const auto *isg_255 = buffer.data(isg + 255);
    const auto *isg_258 = buffer.data(isg + 258);
    const auto *isg_260 = buffer.data(isg + 260);
    const auto *isg_265 = buffer.data(isg + 265);
    const auto *isg_269 = buffer.data(isg + 269);
    const auto *isg_270 = buffer.data(isg + 270);
    const auto *isg_272 = buffer.data(isg + 272);
    const auto *isg_273 = buffer.data(isg + 273);
    const auto *isg_275 = buffer.data(isg + 275);
    const auto *isg_280 = buffer.data(isg + 280);
    const auto *isg_284 = buffer.data(isg + 284);
    const auto *isg_285 = buffer.data(isg + 285);
    const auto *isg_287 = buffer.data(isg + 287);
    const auto *isg_288 = buffer.data(isg + 288);
    const auto *isg_290 = buffer.data(isg + 290);
    const auto *isg_295 = buffer.data(isg + 295);
    const auto *isg_299 = buffer.data(isg + 299);
    const auto *isg_300 = buffer.data(isg + 300);
    const auto *isg_302 = buffer.data(isg + 302);
    const auto *isg_305 = buffer.data(isg + 305);
    const auto *isg_314 = buffer.data(isg + 314);
    const auto *isg_325 = buffer.data(isg + 325);
    const auto *isg_329 = buffer.data(isg + 329);
    const auto *isg_351 = buffer.data(isg + 351);
    const auto *isg_354 = buffer.data(isg + 354);
    const auto *isg_355 = buffer.data(isg + 355);
    const auto *isg_356 = buffer.data(isg + 356);
    const auto *isg_357 = buffer.data(isg + 357);
    const auto *isg_358 = buffer.data(isg + 358);
    const auto *isg_359 = buffer.data(isg + 359);
    const auto *isg_360 = buffer.data(isg + 360);
    const auto *isg_363 = buffer.data(isg + 363);
    const auto *isg_365 = buffer.data(isg + 365);
    const auto *isg_366 = buffer.data(isg + 366);
    const auto *isg_369 = buffer.data(isg + 369);
    const auto *isg_370 = buffer.data(isg + 370);
    const auto *isg_371 = buffer.data(isg + 371);
    const auto *isg_372 = buffer.data(isg + 372);
    const auto *isg_373 = buffer.data(isg + 373);
    const auto *isg_374 = buffer.data(isg + 374);
    const auto *isg_375 = buffer.data(isg + 375);
    const auto *isg_378 = buffer.data(isg + 378);
    const auto *isg_380 = buffer.data(isg + 380);
    const auto *isg_381 = buffer.data(isg + 381);
    const auto *isg_384 = buffer.data(isg + 384);
    const auto *isg_385 = buffer.data(isg + 385);
    const auto *isg_386 = buffer.data(isg + 386);
    const auto *isg_387 = buffer.data(isg + 387);
    const auto *isg_388 = buffer.data(isg + 388);
    const auto *isg_389 = buffer.data(isg + 389);
    const auto *isg_393 = buffer.data(isg + 393);
    const auto *isg_396 = buffer.data(isg + 396);
    const auto *isg_400 = buffer.data(isg + 400);
    const auto *isg_401 = buffer.data(isg + 401);
    const auto *isg_402 = buffer.data(isg + 402);
    const auto *isg_403 = buffer.data(isg + 403);
    const auto *isg_404 = buffer.data(isg + 404);
    const auto *isg_405 = buffer.data(isg + 405);
    const auto *isg_410 = buffer.data(isg + 410);
    const auto *isg_414 = buffer.data(isg + 414);
    const auto *isg_415 = buffer.data(isg + 415);
    const auto *isg_416 = buffer.data(isg + 416);
    const auto *isg_417 = buffer.data(isg + 417);
    const auto *isg_419 = buffer.data(isg + 419);

    const auto *ish1_420 = buffer.data(ish1 + 420);
    const auto *ish1_425 = buffer.data(ish1 + 425);
    const auto *ish1_429 = buffer.data(ish1 + 429);
    const auto *ish1_441 = buffer.data(ish1 + 441);
    const auto *ish1_442 = buffer.data(ish1 + 442);
    const auto *ish1_444 = buffer.data(ish1 + 444);
    const auto *ish1_489 = buffer.data(ish1 + 489);
    const auto *ish1_492 = buffer.data(ish1 + 492);
    const auto *ish1_498 = buffer.data(ish1 + 498);
    const auto *ish1_500 = buffer.data(ish1 + 500);
    const auto *ish1_501 = buffer.data(ish1 + 501);
    const auto *ish1_503 = buffer.data(ish1 + 503);
    const auto *ish1_504 = buffer.data(ish1 + 504);
    const auto *ish1_507 = buffer.data(ish1 + 507);
    const auto *ish1_509 = buffer.data(ish1 + 509);
    const auto *ish1_510 = buffer.data(ish1 + 510);
    const auto *ish1_513 = buffer.data(ish1 + 513);
    const auto *ish1_519 = buffer.data(ish1 + 519);
    const auto *ish1_521 = buffer.data(ish1 + 521);
    const auto *ish1_522 = buffer.data(ish1 + 522);
    const auto *ish1_524 = buffer.data(ish1 + 524);
    const auto *ish1_525 = buffer.data(ish1 + 525);
    const auto *ish1_528 = buffer.data(ish1 + 528);
    const auto *ish1_530 = buffer.data(ish1 + 530);
    const auto *ish1_531 = buffer.data(ish1 + 531);
    const auto *ish1_534 = buffer.data(ish1 + 534);
    const auto *ish1_540 = buffer.data(ish1 + 540);
    const auto *ish1_542 = buffer.data(ish1 + 542);
    const auto *ish1_543 = buffer.data(ish1 + 543);
    const auto *ish1_545 = buffer.data(ish1 + 545);
    const auto *ish1_549 = buffer.data(ish1 + 549);
    const auto *ish1_552 = buffer.data(ish1 + 552);
    const auto *ish1_561 = buffer.data(ish1 + 561);
    const auto *ish1_563 = buffer.data(ish1 + 563);
    const auto *ish1_564 = buffer.data(ish1 + 564);
    const auto *ish1_566 = buffer.data(ish1 + 566);
    const auto *ish1_567 = buffer.data(ish1 + 567);
    const auto *ish1_572 = buffer.data(ish1 + 572);
    const auto *ish1_576 = buffer.data(ish1 + 576);
    const auto *ish1_582 = buffer.data(ish1 + 582);
    const auto *ish1_583 = buffer.data(ish1 + 583);
    const auto *ish1_584 = buffer.data(ish1 + 584);
    const auto *ish1_585 = buffer.data(ish1 + 585);
    const auto *ish1_587 = buffer.data(ish1 + 587);

    const auto *ksf0_270 = buffer.data(ksf0 + 270);
    const auto *ksf0_271 = buffer.data(ksf0 + 271);
    const auto *ksf0_272 = buffer.data(ksf0 + 272);
    const auto *ksf0_280 = buffer.data(ksf0 + 280);
    const auto *ksf0_281 = buffer.data(ksf0 + 281);
    const auto *ksf0_283 = buffer.data(ksf0 + 283);
    const auto *ksf0_285 = buffer.data(ksf0 + 285);
    const auto *ksf0_286 = buffer.data(ksf0 + 286);
    const auto *ksf0_287 = buffer.data(ksf0 + 287);
    const auto *ksf0_288 = buffer.data(ksf0 + 288);
    const auto *ksf0_289 = buffer.data(ksf0 + 289);
    const auto *ksf0_292 = buffer.data(ksf0 + 292);
    const auto *ksf0_294 = buffer.data(ksf0 + 294);

    const auto *ksf1_270 = buffer.data(ksf1 + 270);
    const auto *ksf1_271 = buffer.data(ksf1 + 271);
    const auto *ksf1_272 = buffer.data(ksf1 + 272);
    const auto *ksf1_280 = buffer.data(ksf1 + 280);
    const auto *ksf1_281 = buffer.data(ksf1 + 281);
    const auto *ksf1_283 = buffer.data(ksf1 + 283);
    const auto *ksf1_285 = buffer.data(ksf1 + 285);
    const auto *ksf1_286 = buffer.data(ksf1 + 286);
    const auto *ksf1_287 = buffer.data(ksf1 + 287);
    const auto *ksf1_288 = buffer.data(ksf1 + 288);
    const auto *ksf1_289 = buffer.data(ksf1 + 289);
    const auto *ksf1_292 = buffer.data(ksf1 + 292);
    const auto *ksf1_294 = buffer.data(ksf1 + 294);

    const auto *ksg_348 = buffer.data(ksg + 348);
    const auto *ksg_350 = buffer.data(ksg + 350);
    const auto *ksg_355 = buffer.data(ksg + 355);
    const auto *ksg_356 = buffer.data(ksg + 356);
    const auto *ksg_357 = buffer.data(ksg + 357);
    const auto *ksg_358 = buffer.data(ksg + 358);
    const auto *ksg_359 = buffer.data(ksg + 359);
    const auto *ksg_360 = buffer.data(ksg + 360);
    const auto *ksg_362 = buffer.data(ksg + 362);
    const auto *ksg_363 = buffer.data(ksg + 363);
    const auto *ksg_365 = buffer.data(ksg + 365);
    const auto *ksg_370 = buffer.data(ksg + 370);
    const auto *ksg_371 = buffer.data(ksg + 371);
    const auto *ksg_372 = buffer.data(ksg + 372);
    const auto *ksg_373 = buffer.data(ksg + 373);
    const auto *ksg_374 = buffer.data(ksg + 374);
    const auto *ksg_375 = buffer.data(ksg + 375);
    const auto *ksg_377 = buffer.data(ksg + 377);
    const auto *ksg_378 = buffer.data(ksg + 378);
    const auto *ksg_380 = buffer.data(ksg + 380);
    const auto *ksg_385 = buffer.data(ksg + 385);
    const auto *ksg_386 = buffer.data(ksg + 386);
    const auto *ksg_387 = buffer.data(ksg + 387);
    const auto *ksg_388 = buffer.data(ksg + 388);
    const auto *ksg_389 = buffer.data(ksg + 389);
    const auto *ksg_390 = buffer.data(ksg + 390);
    const auto *ksg_392 = buffer.data(ksg + 392);
    const auto *ksg_393 = buffer.data(ksg + 393);
    const auto *ksg_395 = buffer.data(ksg + 395);
    const auto *ksg_400 = buffer.data(ksg + 400);
    const auto *ksg_401 = buffer.data(ksg + 401);
    const auto *ksg_402 = buffer.data(ksg + 402);
    const auto *ksg_403 = buffer.data(ksg + 403);
    const auto *ksg_404 = buffer.data(ksg + 404);
    const auto *ksg_405 = buffer.data(ksg + 405);
    const auto *ksg_406 = buffer.data(ksg + 406);
    const auto *ksg_407 = buffer.data(ksg + 407);
    const auto *ksg_408 = buffer.data(ksg + 408);
    const auto *ksg_409 = buffer.data(ksg + 409);
    const auto *ksg_410 = buffer.data(ksg + 410);
    const auto *ksg_414 = buffer.data(ksg + 414);
    const auto *ksg_415 = buffer.data(ksg + 415);
    const auto *ksg_416 = buffer.data(ksg + 416);
    const auto *ksg_417 = buffer.data(ksg + 417);
    const auto *ksg_419 = buffer.data(ksg + 419);
    const auto *ksg_420 = buffer.data(ksg + 420);
    const auto *ksg_421 = buffer.data(ksg + 421);
    const auto *ksg_423 = buffer.data(ksg + 423);
    const auto *ksg_425 = buffer.data(ksg + 425);
    const auto *ksg_426 = buffer.data(ksg + 426);
    const auto *ksg_428 = buffer.data(ksg + 428);
    const auto *ksg_429 = buffer.data(ksg + 429);
    const auto *ksg_430 = buffer.data(ksg + 430);
    const auto *ksg_431 = buffer.data(ksg + 431);
    const auto *ksg_432 = buffer.data(ksg + 432);
    const auto *ksg_433 = buffer.data(ksg + 433);
    const auto *ksg_434 = buffer.data(ksg + 434);
    const auto *ksg_437 = buffer.data(ksg + 437);
    const auto *ksg_439 = buffer.data(ksg + 439);

#pragma omp simd aligned(t_489, t_490, t_491, pa_x, pc_x, pc_y, pc_z, ish0_489, isg_243, \
                         isg_260, isg_351, ish1_489, ksg_348, ksg_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = pa_x[k] * ish0_489[k]
                   + f_10 * isg_351[k]
                   - f_8 * pc_x[k] * ish1_489[k];

        t_490[k] = f_10 * isg_243[k]
                   + f_3 * pc_z[k] * ksg_348[k];

        t_491[k] = f_16 * isg_260[k]
                   + f_3 * pc_y[k] * ksg_350[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pa_x, pc_x, ish0_492, isg_354, isg_355, \
                         isg_356, isg_357, ish1_492, ksg_355, ksg_356, \
                         ksg_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = pa_x[k] * ish0_492[k]
                   + f_10 * isg_354[k]
                   - f_8 * pc_x[k] * ish1_492[k];

        t_493[k] = f_9 * isg_355[k]
                   + f_3 * pc_x[k] * ksg_355[k];

        t_494[k] = f_9 * isg_356[k]
                   + f_3 * pc_x[k] * ksg_356[k];

        t_495[k] = f_9 * isg_357[k]
                   + f_3 * pc_x[k] * ksg_357[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pa_x, pc_x, pc_z, ish0_498, isg_250, \
                         isg_358, isg_359, ish1_498, ksg_355, ksg_358, \
                         ksg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_9 * isg_358[k]
                   + f_3 * pc_x[k] * ksg_358[k];

        t_497[k] = f_9 * isg_359[k]
                   + f_3 * pc_x[k] * ksg_359[k];

        t_498[k] = pa_x[k] * ish0_498[k]
                   - f_8 * pc_x[k] * ish1_498[k];

        t_499[k] = f_10 * isg_250[k]
                   + f_3 * pc_z[k] * ksg_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_x, pc_x, pc_y, ish0_500, ish0_501, \
                         ish0_503, isg_269, ish1_500, ish1_501, ish1_503, \
                         ksg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pa_x[k] * ish0_500[k]
                   - f_8 * pc_x[k] * ish1_500[k];

        t_501[k] = pa_x[k] * ish0_501[k]
                   - f_8 * pc_x[k] * ish1_501[k];

        t_502[k] = f_16 * isg_269[k]
                   + f_3 * pc_y[k] * ksg_359[k];

        t_503[k] = pa_x[k] * ish0_503[k]
                   - f_8 * pc_x[k] * ish1_503[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pa_x, pc_x, pc_y, pc_z, ish0_504, isg_255, \
                         isg_270, isg_360, ish1_504, ksg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = pa_x[k] * ish0_504[k]
                   + f_15 * isg_360[k]
                   - f_8 * pc_x[k] * ish1_504[k];

        t_505[k] = f_11 * isg_270[k]
                   + f_3 * pc_y[k] * ksg_360[k];

        t_506[k] = f_11 * isg_255[k]
                   + f_3 * pc_z[k] * ksg_360[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pa_x, pc_x, pc_y, ish0_507, ish0_509, isg_272, \
                         isg_363, isg_365, ish1_507, ish1_509, \
                         ksg_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = pa_x[k] * ish0_507[k]
                   + f_11 * isg_363[k]
                   - f_8 * pc_x[k] * ish1_507[k];

        t_508[k] = f_11 * isg_272[k]
                   + f_3 * pc_y[k] * ksg_362[k];

        t_509[k] = pa_x[k] * ish0_509[k]
                   + f_11 * isg_365[k]
                   - f_8 * pc_x[k] * ish1_509[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pa_x, pc_x, pc_y, pc_z, ish0_510, isg_258, \
                         isg_275, isg_366, ish1_510, ksg_363, ksg_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = pa_x[k] * ish0_510[k]
                   + f_10 * isg_366[k]
                   - f_8 * pc_x[k] * ish1_510[k];

        t_511[k] = f_11 * isg_258[k]
                   + f_3 * pc_z[k] * ksg_363[k];

        t_512[k] = f_11 * isg_275[k]
                   + f_3 * pc_y[k] * ksg_365[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, t_516, pa_x, pc_x, ish0_513, isg_369, isg_370, \
                         isg_371, isg_372, ish1_513, ksg_370, ksg_371, \
                         ksg_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = pa_x[k] * ish0_513[k]
                   + f_10 * isg_369[k]
                   - f_8 * pc_x[k] * ish1_513[k];

        t_514[k] = f_9 * isg_370[k]
                   + f_3 * pc_x[k] * ksg_370[k];

        t_515[k] = f_9 * isg_371[k]
                   + f_3 * pc_x[k] * ksg_371[k];

        t_516[k] = f_9 * isg_372[k]
                   + f_3 * pc_x[k] * ksg_372[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pa_x, pc_x, pc_z, ish0_519, isg_265, \
                         isg_373, isg_374, ish1_519, ksg_370, ksg_373, \
                         ksg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_9 * isg_373[k]
                   + f_3 * pc_x[k] * ksg_373[k];

        t_518[k] = f_9 * isg_374[k]
                   + f_3 * pc_x[k] * ksg_374[k];

        t_519[k] = pa_x[k] * ish0_519[k]
                   - f_8 * pc_x[k] * ish1_519[k];

        t_520[k] = f_11 * isg_265[k]
                   + f_3 * pc_z[k] * ksg_370[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pa_x, pc_x, pc_y, ish0_521, ish0_522, \
                         ish0_524, isg_284, ish1_521, ish1_522, ish1_524, \
                         ksg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = pa_x[k] * ish0_521[k]
                   - f_8 * pc_x[k] * ish1_521[k];

        t_522[k] = pa_x[k] * ish0_522[k]
                   - f_8 * pc_x[k] * ish1_522[k];

        t_523[k] = f_11 * isg_284[k]
                   + f_3 * pc_y[k] * ksg_374[k];

        t_524[k] = pa_x[k] * ish0_524[k]
                   - f_8 * pc_x[k] * ish1_524[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, pa_x, pc_x, pc_y, pc_z, ish0_525, isg_270, \
                         isg_285, isg_375, ish1_525, ksg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = pa_x[k] * ish0_525[k]
                   + f_15 * isg_375[k]
                   - f_8 * pc_x[k] * ish1_525[k];

        t_526[k] = f_10 * isg_285[k]
                   + f_3 * pc_y[k] * ksg_375[k];

        t_527[k] = f_16 * isg_270[k]
                   + f_3 * pc_z[k] * ksg_375[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pa_x, pc_x, pc_y, ish0_528, ish0_530, isg_287, \
                         isg_378, isg_380, ish1_528, ish1_530, \
                         ksg_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = pa_x[k] * ish0_528[k]
                   + f_11 * isg_378[k]
                   - f_8 * pc_x[k] * ish1_528[k];

        t_529[k] = f_10 * isg_287[k]
                   + f_3 * pc_y[k] * ksg_377[k];

        t_530[k] = pa_x[k] * ish0_530[k]
                   + f_11 * isg_380[k]
                   - f_8 * pc_x[k] * ish1_530[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pa_x, pc_x, pc_y, pc_z, ish0_531, isg_273, \
                         isg_290, isg_381, ish1_531, ksg_378, ksg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = pa_x[k] * ish0_531[k]
                   + f_10 * isg_381[k]
                   - f_8 * pc_x[k] * ish1_531[k];

        t_532[k] = f_16 * isg_273[k]
                   + f_3 * pc_z[k] * ksg_378[k];

        t_533[k] = f_10 * isg_290[k]
                   + f_3 * pc_y[k] * ksg_380[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pa_x, pc_x, ish0_534, isg_384, isg_385, \
                         isg_386, isg_387, ish1_534, ksg_385, ksg_386, \
                         ksg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = pa_x[k] * ish0_534[k]
                   + f_10 * isg_384[k]
                   - f_8 * pc_x[k] * ish1_534[k];

        t_535[k] = f_9 * isg_385[k]
                   + f_3 * pc_x[k] * ksg_385[k];

        t_536[k] = f_9 * isg_386[k]
                   + f_3 * pc_x[k] * ksg_386[k];

        t_537[k] = f_9 * isg_387[k]
                   + f_3 * pc_x[k] * ksg_387[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_x, pc_x, pc_z, ish0_540, isg_280, \
                         isg_388, isg_389, ish1_540, ksg_385, ksg_388, \
                         ksg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_9 * isg_388[k]
                   + f_3 * pc_x[k] * ksg_388[k];

        t_539[k] = f_9 * isg_389[k]
                   + f_3 * pc_x[k] * ksg_389[k];

        t_540[k] = pa_x[k] * ish0_540[k]
                   - f_8 * pc_x[k] * ish1_540[k];

        t_541[k] = f_16 * isg_280[k]
                   + f_3 * pc_z[k] * ksg_385[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pa_x, pc_x, pc_y, ish0_542, ish0_543, \
                         ish0_545, isg_299, ish1_542, ish1_543, ish1_545, \
                         ksg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = pa_x[k] * ish0_542[k]
                   - f_8 * pc_x[k] * ish1_542[k];

        t_543[k] = pa_x[k] * ish0_543[k]
                   - f_8 * pc_x[k] * ish1_543[k];

        t_544[k] = f_10 * isg_299[k]
                   + f_3 * pc_y[k] * ksg_389[k];

        t_545[k] = pa_x[k] * ish0_545[k]
                   - f_8 * pc_x[k] * ish1_545[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pa_y, pc_y, pc_z, ish0_420, isg_285, isg_300, \
                         ish1_420, ksg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * ish0_420[k]
                   - f_8 * pc_y[k] * ish1_420[k];

        t_547[k] = f_9 * isg_300[k]
                   + f_3 * pc_y[k] * ksg_390[k];

        t_548[k] = f_15 * isg_285[k]
                   + f_3 * pc_z[k] * ksg_390[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_x, pa_y, pc_x, pc_y, ish0_425, ish0_549, \
                         isg_302, isg_393, ish1_425, ish1_549, \
                         ksg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pa_x[k] * ish0_549[k]
                   + f_11 * isg_393[k]
                   - f_8 * pc_x[k] * ish1_549[k];

        t_550[k] = f_9 * isg_302[k]
                   + f_3 * pc_y[k] * ksg_392[k];

        t_551[k] = pa_y[k] * ish0_425[k]
                   - f_8 * pc_y[k] * ish1_425[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_x, pc_x, pc_y, pc_z, ish0_552, isg_288, \
                         isg_305, isg_396, ish1_552, ksg_393, ksg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_x[k] * ish0_552[k]
                   + f_10 * isg_396[k]
                   - f_8 * pc_x[k] * ish1_552[k];

        t_553[k] = f_15 * isg_288[k]
                   + f_3 * pc_z[k] * ksg_393[k];

        t_554[k] = f_9 * isg_305[k]
                   + f_3 * pc_y[k] * ksg_395[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pa_y, pc_x, pc_y, ish0_429, isg_400, \
                         isg_401, isg_402, ish1_429, ksg_400, ksg_401, \
                         ksg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_y[k] * ish0_429[k]
                   - f_8 * pc_y[k] * ish1_429[k];

        t_556[k] = f_9 * isg_400[k]
                   + f_3 * pc_x[k] * ksg_400[k];

        t_557[k] = f_9 * isg_401[k]
                   + f_3 * pc_x[k] * ksg_401[k];

        t_558[k] = f_9 * isg_402[k]
                   + f_3 * pc_x[k] * ksg_402[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_x, pc_x, pc_z, ish0_561, isg_295, \
                         isg_403, isg_404, ish1_561, ksg_400, ksg_403, \
                         ksg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_9 * isg_403[k]
                   + f_3 * pc_x[k] * ksg_403[k];

        t_560[k] = f_9 * isg_404[k]
                   + f_3 * pc_x[k] * ksg_404[k];

        t_561[k] = pa_x[k] * ish0_561[k]
                   - f_8 * pc_x[k] * ish1_561[k];

        t_562[k] = f_15 * isg_295[k]
                   + f_3 * pc_z[k] * ksg_400[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pa_x, pc_x, pc_y, ish0_563, ish0_564, \
                         ish0_566, isg_314, ish1_563, ish1_564, ish1_566, \
                         ksg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pa_x[k] * ish0_563[k]
                   - f_8 * pc_x[k] * ish1_563[k];

        t_564[k] = pa_x[k] * ish0_564[k]
                   - f_8 * pc_x[k] * ish1_564[k];

        t_565[k] = f_9 * isg_314[k]
                   + f_3 * pc_y[k] * ksg_404[k];

        t_566[k] = pa_x[k] * ish0_566[k]
                   - f_8 * pc_x[k] * ish1_566[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, pa_x, pc_x, pc_y, pc_z, ish0_567, \
                         isg_300, isg_405, ish1_567, ksf0_270, ksf1_270, ksg_405, \
                         ksg_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = pa_x[k] * ish0_567[k]
                   + f_15 * isg_405[k]
                   - f_8 * pc_x[k] * ish1_567[k];

        t_568[k] = f_3 * pc_y[k] * ksg_405[k];

        t_569[k] = f_12 * isg_300[k]
                   + f_3 * pc_z[k] * ksg_405[k];

        t_570[k] = f_4 * ksf0_270[k]
                   - f_5 * ksf1_270[k]
                   + f_3 * pc_y[k] * ksg_406[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pa_x, pc_x, pc_y, ish0_572, isg_410, ish1_572, \
                         ksf0_271, ksf1_271, ksg_407, ksg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_3 * pc_y[k] * ksg_407[k];

        t_572[k] = pa_x[k] * ish0_572[k]
                   + f_11 * isg_410[k]
                   - f_8 * pc_x[k] * ish1_572[k];

        t_573[k] = f_6 * ksf0_271[k]
                   - f_7 * ksf1_271[k]
                   + f_3 * pc_y[k] * ksg_408[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_x, pc_x, pc_y, ish0_576, isg_414, \
                         isg_415, ish1_576, ksf0_272, ksf1_272, ksg_409, ksg_410, \
                         ksg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_4 * ksf0_272[k]
                   - f_5 * ksf1_272[k]
                   + f_3 * pc_y[k] * ksg_409[k];

        t_575[k] = f_3 * pc_y[k] * ksg_410[k];

        t_576[k] = pa_x[k] * ish0_576[k]
                   + f_10 * isg_414[k]
                   - f_8 * pc_x[k] * ish1_576[k];

        t_577[k] = f_9 * isg_415[k]
                   + f_3 * pc_x[k] * ksg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, pc_y, isg_416, isg_417, isg_419, \
                         ksg_414, ksg_416, ksg_417, ksg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_9 * isg_416[k]
                   + f_3 * pc_x[k] * ksg_416[k];

        t_579[k] = f_9 * isg_417[k]
                   + f_3 * pc_x[k] * ksg_417[k];

        t_580[k] = f_3 * pc_y[k] * ksg_414[k];

        t_581[k] = f_9 * isg_419[k]
                   + f_3 * pc_x[k] * ksg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pa_x, pc_x, ish0_582, ish0_583, ish0_584, \
                         ish0_585, ish1_582, ish1_583, ish1_584, \
                         ish1_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pa_x[k] * ish0_582[k]
                   - f_8 * pc_x[k] * ish1_582[k];

        t_583[k] = pa_x[k] * ish0_583[k]
                   - f_8 * pc_x[k] * ish1_583[k];

        t_584[k] = pa_x[k] * ish0_584[k]
                   - f_8 * pc_x[k] * ish1_584[k];

        t_585[k] = pa_x[k] * ish0_585[k]
                   - f_8 * pc_x[k] * ish1_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pa_x, pc_x, pc_y, ish0_587, ish1_587, \
                         ksf0_280, ksf0_281, ksf1_280, ksf1_281, ksg_419, ksg_420, \
                         ksg_421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * ksg_419[k];

        t_587[k] = pa_x[k] * ish0_587[k]
                   - f_8 * pc_x[k] * ish1_587[k];

        t_588[k] = f_1 * ksf0_280[k]
                   - f_2 * ksf1_280[k]
                   + f_3 * pc_x[k] * ksg_420[k];

        t_589[k] = f_13 * ksf0_281[k]
                   - f_14 * ksf1_281[k]
                   + f_3 * pc_x[k] * ksg_421[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, pc_z, ksf0_283, ksf0_285, ksf1_283, \
                         ksf1_285, ksg_420, ksg_421, ksg_423, ksg_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_3 * pc_z[k] * ksg_420[k];

        t_591[k] = f_6 * ksf0_283[k]
                   - f_7 * ksf1_283[k]
                   + f_3 * pc_x[k] * ksg_423[k];

        t_592[k] = f_3 * pc_z[k] * ksg_421[k];

        t_593[k] = f_6 * ksf0_285[k]
                   - f_7 * ksf1_285[k]
                   + f_3 * pc_x[k] * ksg_425[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pc_x, pc_z, ksf0_286, ksf0_288, ksf0_289, \
                         ksf1_286, ksf1_288, ksf1_289, ksg_423, ksg_426, ksg_428, \
                         ksg_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_4 * ksf0_286[k]
                   - f_5 * ksf1_286[k]
                   + f_3 * pc_x[k] * ksg_426[k];

        t_595[k] = f_3 * pc_z[k] * ksg_423[k];

        t_596[k] = f_4 * ksf0_288[k]
                   - f_5 * ksf1_288[k]
                   + f_3 * pc_x[k] * ksg_428[k];

        t_597[k] = f_4 * ksf0_289[k]
                   - f_5 * ksf1_289[k]
                   + f_3 * pc_x[k] * ksg_429[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, t_602, t_603, pc_x, pc_y, isg_325, \
                         ksf0_286, ksf1_286, ksg_430, ksg_431, ksg_432, ksg_433, \
                         ksg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_3 * pc_x[k] * ksg_430[k];

        t_599[k] = f_3 * pc_x[k] * ksg_431[k];

        t_600[k] = f_3 * pc_x[k] * ksg_432[k];

        t_601[k] = f_3 * pc_x[k] * ksg_433[k];

        t_602[k] = f_3 * pc_x[k] * ksg_434[k];

        t_603[k] = f_0 * isg_325[k]
                   + f_1 * ksf0_286[k]
                   - f_2 * ksf1_286[k]
                   + f_3 * pc_y[k] * ksg_430[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pc_y, pc_z, isg_329, ksf0_286, ksf0_287, \
                         ksf1_286, ksf1_287, ksg_430, ksg_431, ksg_432, \
                         ksg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * ksg_430[k];

        t_605[k] = f_4 * ksf0_286[k]
                   - f_5 * ksf1_286[k]
                   + f_3 * pc_z[k] * ksg_431[k];

        t_606[k] = f_6 * ksf0_287[k]
                   - f_7 * ksf1_287[k]
                   + f_3 * pc_z[k] * ksg_432[k];

        t_607[k] = f_0 * isg_329[k]
                   + f_3 * pc_y[k] * ksg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, pa_z, pc_z, ish0_441, ish0_442, ish1_441, \
                         ish1_442, ksf0_289, ksf1_289, ksg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * ksf0_289[k]
                   - f_2 * ksf1_289[k]
                   + f_3 * pc_z[k] * ksg_434[k];

        t_609[k] = pa_z[k] * ish0_441[k]
                   - f_8 * pc_z[k] * ish1_441[k];

        t_610[k] = pa_z[k] * ish0_442[k]
                   - f_8 * pc_z[k] * ish1_442[k];
    }

#pragma omp simd aligned(t_611, t_612, t_613, pa_z, pc_x, pc_z, ish0_444, ish1_444, ksf0_292, \
                         ksf0_294, ksf1_292, ksf1_294, ksg_437, \
                         ksg_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_611[k] = f_13 * ksf0_292[k]
                   - f_14 * ksf1_292[k]
                   + f_3 * pc_x[k] * ksg_437[k];

        t_612[k] = pa_z[k] * ish0_444[k]
                   - f_8 * pc_z[k] * ish1_444[k];

        t_613[k] = f_6 * ksf0_294[k]
                   - f_7 * ksf1_294[k]
                   + f_3 * pc_x[k] * ksg_439[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
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
    const auto f_12 = 3.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / q;

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
    auto *t_727 = buffer.data(target + 727);
    auto *t_728 = buffer.data(target + 728);
    auto *t_729 = buffer.data(target + 729);
    auto *t_730 = buffer.data(target + 730);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ish0_447 = buffer.data(ish0 + 447);
    const auto *ish0_456 = buffer.data(ish0 + 456);
    const auto *ish0_458 = buffer.data(ish0 + 458);
    const auto *ish0_459 = buffer.data(ish0 + 459);
    const auto *ish0_567 = buffer.data(ish0 + 567);
    const auto *ish0_569 = buffer.data(ish0 + 569);
    const auto *ish0_572 = buffer.data(ish0 + 572);
    const auto *ish0_576 = buffer.data(ish0 + 576);
    const auto *ish0_582 = buffer.data(ish0 + 582);

    const auto *isg_325 = buffer.data(isg + 325);
    const auto *isg_326 = buffer.data(isg + 326);
    const auto *isg_327 = buffer.data(isg + 327);
    const auto *isg_329 = buffer.data(isg + 329);
    const auto *isg_340 = buffer.data(isg + 340);
    const auto *isg_344 = buffer.data(isg + 344);
    const auto *isg_355 = buffer.data(isg + 355);
    const auto *isg_357 = buffer.data(isg + 357);
    const auto *isg_358 = buffer.data(isg + 358);
    const auto *isg_359 = buffer.data(isg + 359);
    const auto *isg_370 = buffer.data(isg + 370);
    const auto *isg_372 = buffer.data(isg + 372);
    const auto *isg_373 = buffer.data(isg + 373);
    const auto *isg_374 = buffer.data(isg + 374);
    const auto *isg_385 = buffer.data(isg + 385);
    const auto *isg_387 = buffer.data(isg + 387);
    const auto *isg_388 = buffer.data(isg + 388);
    const auto *isg_389 = buffer.data(isg + 389);
    const auto *isg_400 = buffer.data(isg + 400);
    const auto *isg_402 = buffer.data(isg + 402);
    const auto *isg_403 = buffer.data(isg + 403);
    const auto *isg_404 = buffer.data(isg + 404);
    const auto *isg_415 = buffer.data(isg + 415);

    const auto *ish1_447 = buffer.data(ish1 + 447);
    const auto *ish1_456 = buffer.data(ish1 + 456);
    const auto *ish1_458 = buffer.data(ish1 + 458);
    const auto *ish1_459 = buffer.data(ish1 + 459);
    const auto *ish1_567 = buffer.data(ish1 + 567);
    const auto *ish1_569 = buffer.data(ish1 + 569);
    const auto *ish1_572 = buffer.data(ish1 + 572);
    const auto *ish1_576 = buffer.data(ish1 + 576);
    const auto *ish1_582 = buffer.data(ish1 + 582);

    const auto *ksf0_295 = buffer.data(ksf0 + 295);
    const auto *ksf0_297 = buffer.data(ksf0 + 297);
    const auto *ksf0_298 = buffer.data(ksf0 + 298);
    const auto *ksf0_299 = buffer.data(ksf0 + 299);
    const auto *ksf0_300 = buffer.data(ksf0 + 300);
    const auto *ksf0_301 = buffer.data(ksf0 + 301);
    const auto *ksf0_302 = buffer.data(ksf0 + 302);
    const auto *ksf0_303 = buffer.data(ksf0 + 303);
    const auto *ksf0_304 = buffer.data(ksf0 + 304);
    const auto *ksf0_305 = buffer.data(ksf0 + 305);
    const auto *ksf0_306 = buffer.data(ksf0 + 306);
    const auto *ksf0_307 = buffer.data(ksf0 + 307);
    const auto *ksf0_308 = buffer.data(ksf0 + 308);
    const auto *ksf0_309 = buffer.data(ksf0 + 309);
    const auto *ksf0_310 = buffer.data(ksf0 + 310);
    const auto *ksf0_311 = buffer.data(ksf0 + 311);
    const auto *ksf0_312 = buffer.data(ksf0 + 312);
    const auto *ksf0_313 = buffer.data(ksf0 + 313);
    const auto *ksf0_314 = buffer.data(ksf0 + 314);
    const auto *ksf0_315 = buffer.data(ksf0 + 315);
    const auto *ksf0_316 = buffer.data(ksf0 + 316);
    const auto *ksf0_317 = buffer.data(ksf0 + 317);
    const auto *ksf0_318 = buffer.data(ksf0 + 318);
    const auto *ksf0_319 = buffer.data(ksf0 + 319);
    const auto *ksf0_320 = buffer.data(ksf0 + 320);
    const auto *ksf0_321 = buffer.data(ksf0 + 321);
    const auto *ksf0_322 = buffer.data(ksf0 + 322);
    const auto *ksf0_323 = buffer.data(ksf0 + 323);
    const auto *ksf0_324 = buffer.data(ksf0 + 324);
    const auto *ksf0_325 = buffer.data(ksf0 + 325);
    const auto *ksf0_326 = buffer.data(ksf0 + 326);
    const auto *ksf0_327 = buffer.data(ksf0 + 327);
    const auto *ksf0_328 = buffer.data(ksf0 + 328);
    const auto *ksf0_329 = buffer.data(ksf0 + 329);
    const auto *ksf0_330 = buffer.data(ksf0 + 330);
    const auto *ksf0_331 = buffer.data(ksf0 + 331);
    const auto *ksf0_332 = buffer.data(ksf0 + 332);
    const auto *ksf0_333 = buffer.data(ksf0 + 333);
    const auto *ksf0_334 = buffer.data(ksf0 + 334);
    const auto *ksf0_335 = buffer.data(ksf0 + 335);
    const auto *ksf0_336 = buffer.data(ksf0 + 336);
    const auto *ksf0_337 = buffer.data(ksf0 + 337);
    const auto *ksf0_338 = buffer.data(ksf0 + 338);
    const auto *ksf0_339 = buffer.data(ksf0 + 339);
    const auto *ksf0_341 = buffer.data(ksf0 + 341);
    const auto *ksf0_343 = buffer.data(ksf0 + 343);
    const auto *ksf0_344 = buffer.data(ksf0 + 344);
    const auto *ksf0_346 = buffer.data(ksf0 + 346);
    const auto *ksf0_347 = buffer.data(ksf0 + 347);
    const auto *ksf0_348 = buffer.data(ksf0 + 348);

    const auto *ksf1_295 = buffer.data(ksf1 + 295);
    const auto *ksf1_297 = buffer.data(ksf1 + 297);
    const auto *ksf1_298 = buffer.data(ksf1 + 298);
    const auto *ksf1_299 = buffer.data(ksf1 + 299);
    const auto *ksf1_300 = buffer.data(ksf1 + 300);
    const auto *ksf1_301 = buffer.data(ksf1 + 301);
    const auto *ksf1_302 = buffer.data(ksf1 + 302);
    const auto *ksf1_303 = buffer.data(ksf1 + 303);
    const auto *ksf1_304 = buffer.data(ksf1 + 304);
    const auto *ksf1_305 = buffer.data(ksf1 + 305);
    const auto *ksf1_306 = buffer.data(ksf1 + 306);
    const auto *ksf1_307 = buffer.data(ksf1 + 307);
    const auto *ksf1_308 = buffer.data(ksf1 + 308);
    const auto *ksf1_309 = buffer.data(ksf1 + 309);
    const auto *ksf1_310 = buffer.data(ksf1 + 310);
    const auto *ksf1_311 = buffer.data(ksf1 + 311);
    const auto *ksf1_312 = buffer.data(ksf1 + 312);
    const auto *ksf1_313 = buffer.data(ksf1 + 313);
    const auto *ksf1_314 = buffer.data(ksf1 + 314);
    const auto *ksf1_315 = buffer.data(ksf1 + 315);
    const auto *ksf1_316 = buffer.data(ksf1 + 316);
    const auto *ksf1_317 = buffer.data(ksf1 + 317);
    const auto *ksf1_318 = buffer.data(ksf1 + 318);
    const auto *ksf1_319 = buffer.data(ksf1 + 319);
    const auto *ksf1_320 = buffer.data(ksf1 + 320);
    const auto *ksf1_321 = buffer.data(ksf1 + 321);
    const auto *ksf1_322 = buffer.data(ksf1 + 322);
    const auto *ksf1_323 = buffer.data(ksf1 + 323);
    const auto *ksf1_324 = buffer.data(ksf1 + 324);
    const auto *ksf1_325 = buffer.data(ksf1 + 325);
    const auto *ksf1_326 = buffer.data(ksf1 + 326);
    const auto *ksf1_327 = buffer.data(ksf1 + 327);
    const auto *ksf1_328 = buffer.data(ksf1 + 328);
    const auto *ksf1_329 = buffer.data(ksf1 + 329);
    const auto *ksf1_330 = buffer.data(ksf1 + 330);
    const auto *ksf1_331 = buffer.data(ksf1 + 331);
    const auto *ksf1_332 = buffer.data(ksf1 + 332);
    const auto *ksf1_333 = buffer.data(ksf1 + 333);
    const auto *ksf1_334 = buffer.data(ksf1 + 334);
    const auto *ksf1_335 = buffer.data(ksf1 + 335);
    const auto *ksf1_336 = buffer.data(ksf1 + 336);
    const auto *ksf1_337 = buffer.data(ksf1 + 337);
    const auto *ksf1_338 = buffer.data(ksf1 + 338);
    const auto *ksf1_339 = buffer.data(ksf1 + 339);
    const auto *ksf1_341 = buffer.data(ksf1 + 341);
    const auto *ksf1_343 = buffer.data(ksf1 + 343);
    const auto *ksf1_344 = buffer.data(ksf1 + 344);
    const auto *ksf1_346 = buffer.data(ksf1 + 346);
    const auto *ksf1_347 = buffer.data(ksf1 + 347);
    const auto *ksf1_348 = buffer.data(ksf1 + 348);

    const auto *ksg_440 = buffer.data(ksg + 440);
    const auto *ksg_442 = buffer.data(ksg + 442);
    const auto *ksg_443 = buffer.data(ksg + 443);
    const auto *ksg_444 = buffer.data(ksg + 444);
    const auto *ksg_445 = buffer.data(ksg + 445);
    const auto *ksg_446 = buffer.data(ksg + 446);
    const auto *ksg_447 = buffer.data(ksg + 447);
    const auto *ksg_448 = buffer.data(ksg + 448);
    const auto *ksg_449 = buffer.data(ksg + 449);
    const auto *ksg_450 = buffer.data(ksg + 450);
    const auto *ksg_451 = buffer.data(ksg + 451);
    const auto *ksg_452 = buffer.data(ksg + 452);
    const auto *ksg_453 = buffer.data(ksg + 453);
    const auto *ksg_454 = buffer.data(ksg + 454);
    const auto *ksg_455 = buffer.data(ksg + 455);
    const auto *ksg_456 = buffer.data(ksg + 456);
    const auto *ksg_457 = buffer.data(ksg + 457);
    const auto *ksg_458 = buffer.data(ksg + 458);
    const auto *ksg_459 = buffer.data(ksg + 459);
    const auto *ksg_460 = buffer.data(ksg + 460);
    const auto *ksg_461 = buffer.data(ksg + 461);
    const auto *ksg_462 = buffer.data(ksg + 462);
    const auto *ksg_463 = buffer.data(ksg + 463);
    const auto *ksg_464 = buffer.data(ksg + 464);
    const auto *ksg_465 = buffer.data(ksg + 465);
    const auto *ksg_466 = buffer.data(ksg + 466);
    const auto *ksg_467 = buffer.data(ksg + 467);
    const auto *ksg_468 = buffer.data(ksg + 468);
    const auto *ksg_469 = buffer.data(ksg + 469);
    const auto *ksg_470 = buffer.data(ksg + 470);
    const auto *ksg_471 = buffer.data(ksg + 471);
    const auto *ksg_472 = buffer.data(ksg + 472);
    const auto *ksg_473 = buffer.data(ksg + 473);
    const auto *ksg_474 = buffer.data(ksg + 474);
    const auto *ksg_475 = buffer.data(ksg + 475);
    const auto *ksg_476 = buffer.data(ksg + 476);
    const auto *ksg_477 = buffer.data(ksg + 477);
    const auto *ksg_478 = buffer.data(ksg + 478);
    const auto *ksg_479 = buffer.data(ksg + 479);
    const auto *ksg_480 = buffer.data(ksg + 480);
    const auto *ksg_481 = buffer.data(ksg + 481);
    const auto *ksg_482 = buffer.data(ksg + 482);
    const auto *ksg_483 = buffer.data(ksg + 483);
    const auto *ksg_484 = buffer.data(ksg + 484);
    const auto *ksg_485 = buffer.data(ksg + 485);
    const auto *ksg_486 = buffer.data(ksg + 486);
    const auto *ksg_487 = buffer.data(ksg + 487);
    const auto *ksg_488 = buffer.data(ksg + 488);
    const auto *ksg_489 = buffer.data(ksg + 489);
    const auto *ksg_490 = buffer.data(ksg + 490);
    const auto *ksg_491 = buffer.data(ksg + 491);
    const auto *ksg_492 = buffer.data(ksg + 492);
    const auto *ksg_493 = buffer.data(ksg + 493);
    const auto *ksg_494 = buffer.data(ksg + 494);
    const auto *ksg_495 = buffer.data(ksg + 495);
    const auto *ksg_496 = buffer.data(ksg + 496);
    const auto *ksg_497 = buffer.data(ksg + 497);
    const auto *ksg_498 = buffer.data(ksg + 498);
    const auto *ksg_499 = buffer.data(ksg + 499);
    const auto *ksg_500 = buffer.data(ksg + 500);
    const auto *ksg_501 = buffer.data(ksg + 501);
    const auto *ksg_502 = buffer.data(ksg + 502);
    const auto *ksg_503 = buffer.data(ksg + 503);
    const auto *ksg_504 = buffer.data(ksg + 504);
    const auto *ksg_505 = buffer.data(ksg + 505);
    const auto *ksg_506 = buffer.data(ksg + 506);
    const auto *ksg_507 = buffer.data(ksg + 507);
    const auto *ksg_508 = buffer.data(ksg + 508);
    const auto *ksg_509 = buffer.data(ksg + 509);
    const auto *ksg_511 = buffer.data(ksg + 511);
    const auto *ksg_513 = buffer.data(ksg + 513);
    const auto *ksg_514 = buffer.data(ksg + 514);
    const auto *ksg_516 = buffer.data(ksg + 516);
    const auto *ksg_517 = buffer.data(ksg + 517);
    const auto *ksg_518 = buffer.data(ksg + 518);
    const auto *ksg_520 = buffer.data(ksg + 520);
    const auto *ksg_521 = buffer.data(ksg + 521);
    const auto *ksg_522 = buffer.data(ksg + 522);
    const auto *ksg_523 = buffer.data(ksg + 523);
    const auto *ksg_524 = buffer.data(ksg + 524);

#pragma omp simd aligned(t_614, t_615, t_616, pa_z, pc_x, pc_z, ish0_447, ish1_447, ksf0_295, \
                         ksf0_297, ksf1_295, ksf1_297, ksg_440, \
                         ksg_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_6 * ksf0_295[k]
                   - f_7 * ksf1_295[k]
                   + f_3 * pc_x[k] * ksg_440[k];

        t_615[k] = pa_z[k] * ish0_447[k]
                   - f_8 * pc_z[k] * ish1_447[k];

        t_616[k] = f_4 * ksf0_297[k]
                   - f_5 * ksf1_297[k]
                   + f_3 * pc_x[k] * ksg_442[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, t_621, pc_x, ksf0_298, ksf0_299, \
                         ksf1_298, ksf1_299, ksg_443, ksg_444, ksg_445, ksg_446, \
                         ksg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_4 * ksf0_298[k]
                   - f_5 * ksf1_298[k]
                   + f_3 * pc_x[k] * ksg_443[k];

        t_618[k] = f_4 * ksf0_299[k]
                   - f_5 * ksf1_299[k]
                   + f_3 * pc_x[k] * ksg_444[k];

        t_619[k] = f_3 * pc_x[k] * ksg_445[k];

        t_620[k] = f_3 * pc_x[k] * ksg_446[k];

        t_621[k] = f_3 * pc_x[k] * ksg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_z, pc_x, pc_z, ish0_456, isg_325, \
                         ish1_456, ksg_445, ksg_448, ksg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_3 * pc_x[k] * ksg_448[k];

        t_623[k] = f_3 * pc_x[k] * ksg_449[k];

        t_624[k] = pa_z[k] * ish0_456[k]
                   - f_8 * pc_z[k] * ish1_456[k];

        t_625[k] = f_9 * isg_325[k]
                   + f_3 * pc_z[k] * ksg_445[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_z, pc_y, pc_z, ish0_458, ish0_459, isg_326, \
                         isg_327, isg_344, ish1_458, ish1_459, \
                         ksg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_z[k] * ish0_458[k]
                   + f_10 * isg_326[k]
                   - f_8 * pc_z[k] * ish1_458[k];

        t_627[k] = pa_z[k] * ish0_459[k]
                   + f_11 * isg_327[k]
                   - f_8 * pc_z[k] * ish1_459[k];

        t_628[k] = f_12 * isg_344[k]
                   + f_3 * pc_y[k] * ksg_449[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, pc_x, pc_z, isg_329, ksf0_299, ksf0_300, \
                         ksf0_301, ksf1_299, ksf1_300, ksf1_301, ksg_449, ksg_450, \
                         ksg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_9 * isg_329[k]
                   + f_1 * ksf0_299[k]
                   - f_2 * ksf1_299[k]
                   + f_3 * pc_z[k] * ksg_449[k];

        t_630[k] = f_1 * ksf0_300[k]
                   - f_2 * ksf1_300[k]
                   + f_3 * pc_x[k] * ksg_450[k];

        t_631[k] = f_13 * ksf0_301[k]
                   - f_14 * ksf1_301[k]
                   + f_3 * pc_x[k] * ksg_451[k];
    }

#pragma omp simd aligned(t_632, t_633, t_634, pc_x, ksf0_302, ksf0_303, ksf0_304, ksf1_302, \
                         ksf1_303, ksf1_304, ksg_452, ksg_453, \
                         ksg_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_632[k] = f_13 * ksf0_302[k]
                   - f_14 * ksf1_302[k]
                   + f_3 * pc_x[k] * ksg_452[k];

        t_633[k] = f_6 * ksf0_303[k]
                   - f_7 * ksf1_303[k]
                   + f_3 * pc_x[k] * ksg_453[k];

        t_634[k] = f_6 * ksf0_304[k]
                   - f_7 * ksf1_304[k]
                   + f_3 * pc_x[k] * ksg_454[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, pc_x, ksf0_305, ksf0_306, ksf0_307, ksf1_305, \
                         ksf1_306, ksf1_307, ksg_455, ksg_456, \
                         ksg_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_6 * ksf0_305[k]
                   - f_7 * ksf1_305[k]
                   + f_3 * pc_x[k] * ksg_455[k];

        t_636[k] = f_4 * ksf0_306[k]
                   - f_5 * ksf1_306[k]
                   + f_3 * pc_x[k] * ksg_456[k];

        t_637[k] = f_4 * ksf0_307[k]
                   - f_5 * ksf1_307[k]
                   + f_3 * pc_x[k] * ksg_457[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, t_642, pc_x, ksf0_308, ksf0_309, \
                         ksf1_308, ksf1_309, ksg_458, ksg_459, ksg_460, ksg_461, \
                         ksg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_4 * ksf0_308[k]
                   - f_5 * ksf1_308[k]
                   + f_3 * pc_x[k] * ksg_458[k];

        t_639[k] = f_4 * ksf0_309[k]
                   - f_5 * ksf1_309[k]
                   + f_3 * pc_x[k] * ksg_459[k];

        t_640[k] = f_3 * pc_x[k] * ksg_460[k];

        t_641[k] = f_3 * pc_x[k] * ksg_461[k];

        t_642[k] = f_3 * pc_x[k] * ksg_462[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pc_x, pc_y, pc_z, isg_340, isg_355, \
                         ksf0_306, ksf1_306, ksg_460, ksg_463, \
                         ksg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_3 * pc_x[k] * ksg_463[k];

        t_644[k] = f_3 * pc_x[k] * ksg_464[k];

        t_645[k] = f_15 * isg_355[k]
                   + f_1 * ksf0_306[k]
                   - f_2 * ksf1_306[k]
                   + f_3 * pc_y[k] * ksg_460[k];

        t_646[k] = f_10 * isg_340[k]
                   + f_3 * pc_z[k] * ksg_460[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_y, isg_357, isg_358, isg_359, ksf0_308, \
                         ksf0_309, ksf1_308, ksf1_309, ksg_462, ksg_463, \
                         ksg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_15 * isg_357[k]
                   + f_6 * ksf0_308[k]
                   - f_7 * ksf1_308[k]
                   + f_3 * pc_y[k] * ksg_462[k];

        t_648[k] = f_15 * isg_358[k]
                   + f_4 * ksf0_309[k]
                   - f_5 * ksf1_309[k]
                   + f_3 * pc_y[k] * ksg_463[k];

        t_649[k] = f_15 * isg_359[k]
                   + f_3 * pc_y[k] * ksg_464[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_z, isg_344, ksf0_309, ksf0_310, \
                         ksf0_311, ksf1_309, ksf1_310, ksf1_311, ksg_464, ksg_465, \
                         ksg_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_10 * isg_344[k]
                   + f_1 * ksf0_309[k]
                   - f_2 * ksf1_309[k]
                   + f_3 * pc_z[k] * ksg_464[k];

        t_651[k] = f_1 * ksf0_310[k]
                   - f_2 * ksf1_310[k]
                   + f_3 * pc_x[k] * ksg_465[k];

        t_652[k] = f_13 * ksf0_311[k]
                   - f_14 * ksf1_311[k]
                   + f_3 * pc_x[k] * ksg_466[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, ksf0_312, ksf0_313, ksf0_314, ksf1_312, \
                         ksf1_313, ksf1_314, ksg_467, ksg_468, \
                         ksg_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_13 * ksf0_312[k]
                   - f_14 * ksf1_312[k]
                   + f_3 * pc_x[k] * ksg_467[k];

        t_654[k] = f_6 * ksf0_313[k]
                   - f_7 * ksf1_313[k]
                   + f_3 * pc_x[k] * ksg_468[k];

        t_655[k] = f_6 * ksf0_314[k]
                   - f_7 * ksf1_314[k]
                   + f_3 * pc_x[k] * ksg_469[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, ksf0_315, ksf0_316, ksf0_317, ksf1_315, \
                         ksf1_316, ksf1_317, ksg_470, ksg_471, \
                         ksg_472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_6 * ksf0_315[k]
                   - f_7 * ksf1_315[k]
                   + f_3 * pc_x[k] * ksg_470[k];

        t_657[k] = f_4 * ksf0_316[k]
                   - f_5 * ksf1_316[k]
                   + f_3 * pc_x[k] * ksg_471[k];

        t_658[k] = f_4 * ksf0_317[k]
                   - f_5 * ksf1_317[k]
                   + f_3 * pc_x[k] * ksg_472[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, ksf0_318, ksf0_319, \
                         ksf1_318, ksf1_319, ksg_473, ksg_474, ksg_475, ksg_476, \
                         ksg_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_4 * ksf0_318[k]
                   - f_5 * ksf1_318[k]
                   + f_3 * pc_x[k] * ksg_473[k];

        t_660[k] = f_4 * ksf0_319[k]
                   - f_5 * ksf1_319[k]
                   + f_3 * pc_x[k] * ksg_474[k];

        t_661[k] = f_3 * pc_x[k] * ksg_475[k];

        t_662[k] = f_3 * pc_x[k] * ksg_476[k];

        t_663[k] = f_3 * pc_x[k] * ksg_477[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, t_667, pc_x, pc_y, pc_z, isg_355, isg_370, \
                         ksf0_316, ksf1_316, ksg_475, ksg_478, \
                         ksg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_3 * pc_x[k] * ksg_478[k];

        t_665[k] = f_3 * pc_x[k] * ksg_479[k];

        t_666[k] = f_16 * isg_370[k]
                   + f_1 * ksf0_316[k]
                   - f_2 * ksf1_316[k]
                   + f_3 * pc_y[k] * ksg_475[k];

        t_667[k] = f_11 * isg_355[k]
                   + f_3 * pc_z[k] * ksg_475[k];
    }

#pragma omp simd aligned(t_668, t_669, t_670, pc_y, isg_372, isg_373, isg_374, ksf0_318, \
                         ksf0_319, ksf1_318, ksf1_319, ksg_477, ksg_478, \
                         ksg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_668[k] = f_16 * isg_372[k]
                   + f_6 * ksf0_318[k]
                   - f_7 * ksf1_318[k]
                   + f_3 * pc_y[k] * ksg_477[k];

        t_669[k] = f_16 * isg_373[k]
                   + f_4 * ksf0_319[k]
                   - f_5 * ksf1_319[k]
                   + f_3 * pc_y[k] * ksg_478[k];

        t_670[k] = f_16 * isg_374[k]
                   + f_3 * pc_y[k] * ksg_479[k];
    }

#pragma omp simd aligned(t_671, t_672, t_673, pc_x, pc_z, isg_359, ksf0_319, ksf0_320, \
                         ksf0_321, ksf1_319, ksf1_320, ksf1_321, ksg_479, ksg_480, \
                         ksg_481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_671[k] = f_11 * isg_359[k]
                   + f_1 * ksf0_319[k]
                   - f_2 * ksf1_319[k]
                   + f_3 * pc_z[k] * ksg_479[k];

        t_672[k] = f_1 * ksf0_320[k]
                   - f_2 * ksf1_320[k]
                   + f_3 * pc_x[k] * ksg_480[k];

        t_673[k] = f_13 * ksf0_321[k]
                   - f_14 * ksf1_321[k]
                   + f_3 * pc_x[k] * ksg_481[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pc_x, ksf0_322, ksf0_323, ksf0_324, ksf1_322, \
                         ksf1_323, ksf1_324, ksg_482, ksg_483, \
                         ksg_484 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_13 * ksf0_322[k]
                   - f_14 * ksf1_322[k]
                   + f_3 * pc_x[k] * ksg_482[k];

        t_675[k] = f_6 * ksf0_323[k]
                   - f_7 * ksf1_323[k]
                   + f_3 * pc_x[k] * ksg_483[k];

        t_676[k] = f_6 * ksf0_324[k]
                   - f_7 * ksf1_324[k]
                   + f_3 * pc_x[k] * ksg_484[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, ksf0_325, ksf0_326, ksf0_327, ksf1_325, \
                         ksf1_326, ksf1_327, ksg_485, ksg_486, \
                         ksg_487 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_6 * ksf0_325[k]
                   - f_7 * ksf1_325[k]
                   + f_3 * pc_x[k] * ksg_485[k];

        t_678[k] = f_4 * ksf0_326[k]
                   - f_5 * ksf1_326[k]
                   + f_3 * pc_x[k] * ksg_486[k];

        t_679[k] = f_4 * ksf0_327[k]
                   - f_5 * ksf1_327[k]
                   + f_3 * pc_x[k] * ksg_487[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, t_684, pc_x, ksf0_328, ksf0_329, \
                         ksf1_328, ksf1_329, ksg_488, ksg_489, ksg_490, ksg_491, \
                         ksg_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * ksf0_328[k]
                   - f_5 * ksf1_328[k]
                   + f_3 * pc_x[k] * ksg_488[k];

        t_681[k] = f_4 * ksf0_329[k]
                   - f_5 * ksf1_329[k]
                   + f_3 * pc_x[k] * ksg_489[k];

        t_682[k] = f_3 * pc_x[k] * ksg_490[k];

        t_683[k] = f_3 * pc_x[k] * ksg_491[k];

        t_684[k] = f_3 * pc_x[k] * ksg_492[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, pc_x, pc_y, pc_z, isg_370, isg_385, \
                         ksf0_326, ksf1_326, ksg_490, ksg_493, \
                         ksg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_3 * pc_x[k] * ksg_493[k];

        t_686[k] = f_3 * pc_x[k] * ksg_494[k];

        t_687[k] = f_11 * isg_385[k]
                   + f_1 * ksf0_326[k]
                   - f_2 * ksf1_326[k]
                   + f_3 * pc_y[k] * ksg_490[k];

        t_688[k] = f_16 * isg_370[k]
                   + f_3 * pc_z[k] * ksg_490[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pc_y, isg_387, isg_388, isg_389, ksf0_328, \
                         ksf0_329, ksf1_328, ksf1_329, ksg_492, ksg_493, \
                         ksg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_11 * isg_387[k]
                   + f_6 * ksf0_328[k]
                   - f_7 * ksf1_328[k]
                   + f_3 * pc_y[k] * ksg_492[k];

        t_690[k] = f_11 * isg_388[k]
                   + f_4 * ksf0_329[k]
                   - f_5 * ksf1_329[k]
                   + f_3 * pc_y[k] * ksg_493[k];

        t_691[k] = f_11 * isg_389[k]
                   + f_3 * pc_y[k] * ksg_494[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pc_x, pc_z, isg_374, ksf0_329, ksf0_330, \
                         ksf0_331, ksf1_329, ksf1_330, ksf1_331, ksg_494, ksg_495, \
                         ksg_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_16 * isg_374[k]
                   + f_1 * ksf0_329[k]
                   - f_2 * ksf1_329[k]
                   + f_3 * pc_z[k] * ksg_494[k];

        t_693[k] = f_1 * ksf0_330[k]
                   - f_2 * ksf1_330[k]
                   + f_3 * pc_x[k] * ksg_495[k];

        t_694[k] = f_13 * ksf0_331[k]
                   - f_14 * ksf1_331[k]
                   + f_3 * pc_x[k] * ksg_496[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pc_x, ksf0_332, ksf0_333, ksf0_334, ksf1_332, \
                         ksf1_333, ksf1_334, ksg_497, ksg_498, \
                         ksg_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_13 * ksf0_332[k]
                   - f_14 * ksf1_332[k]
                   + f_3 * pc_x[k] * ksg_497[k];

        t_696[k] = f_6 * ksf0_333[k]
                   - f_7 * ksf1_333[k]
                   + f_3 * pc_x[k] * ksg_498[k];

        t_697[k] = f_6 * ksf0_334[k]
                   - f_7 * ksf1_334[k]
                   + f_3 * pc_x[k] * ksg_499[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, ksf0_335, ksf0_336, ksf0_337, ksf1_335, \
                         ksf1_336, ksf1_337, ksg_500, ksg_501, \
                         ksg_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_6 * ksf0_335[k]
                   - f_7 * ksf1_335[k]
                   + f_3 * pc_x[k] * ksg_500[k];

        t_699[k] = f_4 * ksf0_336[k]
                   - f_5 * ksf1_336[k]
                   + f_3 * pc_x[k] * ksg_501[k];

        t_700[k] = f_4 * ksf0_337[k]
                   - f_5 * ksf1_337[k]
                   + f_3 * pc_x[k] * ksg_502[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, t_705, pc_x, ksf0_338, ksf0_339, \
                         ksf1_338, ksf1_339, ksg_503, ksg_504, ksg_505, ksg_506, \
                         ksg_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_4 * ksf0_338[k]
                   - f_5 * ksf1_338[k]
                   + f_3 * pc_x[k] * ksg_503[k];

        t_702[k] = f_4 * ksf0_339[k]
                   - f_5 * ksf1_339[k]
                   + f_3 * pc_x[k] * ksg_504[k];

        t_703[k] = f_3 * pc_x[k] * ksg_505[k];

        t_704[k] = f_3 * pc_x[k] * ksg_506[k];

        t_705[k] = f_3 * pc_x[k] * ksg_507[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pc_x, pc_y, pc_z, isg_385, isg_400, \
                         ksf0_336, ksf1_336, ksg_505, ksg_508, \
                         ksg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_3 * pc_x[k] * ksg_508[k];

        t_707[k] = f_3 * pc_x[k] * ksg_509[k];

        t_708[k] = f_10 * isg_400[k]
                   + f_1 * ksf0_336[k]
                   - f_2 * ksf1_336[k]
                   + f_3 * pc_y[k] * ksg_505[k];

        t_709[k] = f_15 * isg_385[k]
                   + f_3 * pc_z[k] * ksg_505[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, pc_y, isg_402, isg_403, isg_404, ksf0_338, \
                         ksf0_339, ksf1_338, ksf1_339, ksg_507, ksg_508, \
                         ksg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_10 * isg_402[k]
                   + f_6 * ksf0_338[k]
                   - f_7 * ksf1_338[k]
                   + f_3 * pc_y[k] * ksg_507[k];

        t_711[k] = f_10 * isg_403[k]
                   + f_4 * ksf0_339[k]
                   - f_5 * ksf1_339[k]
                   + f_3 * pc_y[k] * ksg_508[k];

        t_712[k] = f_10 * isg_404[k]
                   + f_3 * pc_y[k] * ksg_509[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, pa_y, pc_x, pc_y, pc_z, ish0_567, isg_389, \
                         ish1_567, ksf0_339, ksf0_341, ksf1_339, ksf1_341, ksg_509, \
                         ksg_511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_15 * isg_389[k]
                   + f_1 * ksf0_339[k]
                   - f_2 * ksf1_339[k]
                   + f_3 * pc_z[k] * ksg_509[k];

        t_714[k] = pa_y[k] * ish0_567[k]
                   - f_8 * pc_y[k] * ish1_567[k];

        t_715[k] = f_13 * ksf0_341[k]
                   - f_14 * ksf1_341[k]
                   + f_3 * pc_x[k] * ksg_511[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, pa_y, pc_x, pc_y, ish0_569, ish1_569, ksf0_343, \
                         ksf0_344, ksf1_343, ksf1_344, ksg_513, \
                         ksg_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = pa_y[k] * ish0_569[k]
                   - f_8 * pc_y[k] * ish1_569[k];

        t_717[k] = f_6 * ksf0_343[k]
                   - f_7 * ksf1_343[k]
                   + f_3 * pc_x[k] * ksg_513[k];

        t_718[k] = f_6 * ksf0_344[k]
                   - f_7 * ksf1_344[k]
                   + f_3 * pc_x[k] * ksg_514[k];
    }

#pragma omp simd aligned(t_719, t_720, t_721, pa_y, pc_x, pc_y, ish0_572, ish1_572, ksf0_346, \
                         ksf0_347, ksf1_346, ksf1_347, ksg_516, \
                         ksg_517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_719[k] = pa_y[k] * ish0_572[k]
                   - f_8 * pc_y[k] * ish1_572[k];

        t_720[k] = f_4 * ksf0_346[k]
                   - f_5 * ksf1_346[k]
                   + f_3 * pc_x[k] * ksg_516[k];

        t_721[k] = f_4 * ksf0_347[k]
                   - f_5 * ksf1_347[k]
                   + f_3 * pc_x[k] * ksg_517[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, t_726, pa_y, pc_x, pc_y, ish0_576, \
                         ish1_576, ksf0_348, ksf1_348, ksg_518, ksg_520, ksg_521, \
                         ksg_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_4 * ksf0_348[k]
                   - f_5 * ksf1_348[k]
                   + f_3 * pc_x[k] * ksg_518[k];

        t_723[k] = pa_y[k] * ish0_576[k]
                   - f_8 * pc_y[k] * ish1_576[k];

        t_724[k] = f_3 * pc_x[k] * ksg_520[k];

        t_725[k] = f_3 * pc_x[k] * ksg_521[k];

        t_726[k] = f_3 * pc_x[k] * ksg_522[k];
    }

#pragma omp simd aligned(t_727, t_728, t_729, t_730, pa_y, pc_x, pc_y, pc_z, ish0_582, \
                         isg_400, isg_415, ish1_582, ksg_520, ksg_523, \
                         ksg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_727[k] = f_3 * pc_x[k] * ksg_523[k];

        t_728[k] = f_3 * pc_x[k] * ksg_524[k];

        t_729[k] = pa_y[k] * ish0_582[k]
                   + f_15 * isg_415[k]
                   - f_8 * pc_y[k] * ish1_582[k];

        t_730[k] = f_12 * isg_400[k]
                   + f_3 * pc_z[k] * ksg_520[k];
    }
}

static auto
compute_prim_ksh_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ish0,
                                                          const size_t isg, const size_t ish1,
                                                          const size_t ksf0, const size_t ksf1,
                                                          const size_t ksg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ish0_584 = buffer.data(ish0 + 584);
    const auto *ish0_585 = buffer.data(ish0 + 585);
    const auto *ish0_587 = buffer.data(ish0 + 587);

    const auto *isg_417 = buffer.data(isg + 417);
    const auto *isg_418 = buffer.data(isg + 418);
    const auto *isg_419 = buffer.data(isg + 419);

    const auto *ish1_584 = buffer.data(ish1 + 584);
    const auto *ish1_585 = buffer.data(ish1 + 585);
    const auto *ish1_587 = buffer.data(ish1 + 587);

    const auto *ksf0_350 = buffer.data(ksf0 + 350);
    const auto *ksf0_352 = buffer.data(ksf0 + 352);
    const auto *ksf0_353 = buffer.data(ksf0 + 353);
    const auto *ksf0_355 = buffer.data(ksf0 + 355);
    const auto *ksf0_356 = buffer.data(ksf0 + 356);
    const auto *ksf0_357 = buffer.data(ksf0 + 357);
    const auto *ksf0_358 = buffer.data(ksf0 + 358);
    const auto *ksf0_359 = buffer.data(ksf0 + 359);

    const auto *ksf1_350 = buffer.data(ksf1 + 350);
    const auto *ksf1_352 = buffer.data(ksf1 + 352);
    const auto *ksf1_353 = buffer.data(ksf1 + 353);
    const auto *ksf1_355 = buffer.data(ksf1 + 355);
    const auto *ksf1_356 = buffer.data(ksf1 + 356);
    const auto *ksf1_357 = buffer.data(ksf1 + 357);
    const auto *ksf1_358 = buffer.data(ksf1 + 358);
    const auto *ksf1_359 = buffer.data(ksf1 + 359);

    const auto *ksg_524 = buffer.data(ksg + 524);
    const auto *ksg_525 = buffer.data(ksg + 525);
    const auto *ksg_527 = buffer.data(ksg + 527);
    const auto *ksg_528 = buffer.data(ksg + 528);
    const auto *ksg_530 = buffer.data(ksg + 530);
    const auto *ksg_531 = buffer.data(ksg + 531);
    const auto *ksg_532 = buffer.data(ksg + 532);
    const auto *ksg_534 = buffer.data(ksg + 534);
    const auto *ksg_535 = buffer.data(ksg + 535);
    const auto *ksg_536 = buffer.data(ksg + 536);
    const auto *ksg_537 = buffer.data(ksg + 537);
    const auto *ksg_538 = buffer.data(ksg + 538);
    const auto *ksg_539 = buffer.data(ksg + 539);

#pragma omp simd aligned(t_731, t_732, t_733, t_734, pa_y, pc_y, ish0_584, ish0_585, ish0_587, \
                         isg_417, isg_418, isg_419, ish1_584, ish1_585, ish1_587, \
                         ksg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = pa_y[k] * ish0_584[k]
                   + f_11 * isg_417[k]
                   - f_8 * pc_y[k] * ish1_584[k];

        t_732[k] = pa_y[k] * ish0_585[k]
                   + f_10 * isg_418[k]
                   - f_8 * pc_y[k] * ish1_585[k];

        t_733[k] = f_9 * isg_419[k]
                   + f_3 * pc_y[k] * ksg_524[k];

        t_734[k] = pa_y[k] * ish0_587[k]
                   - f_8 * pc_y[k] * ish1_587[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, t_738, t_739, pc_x, pc_y, ksf0_350, ksf0_352, \
                         ksf0_353, ksf1_350, ksf1_352, ksf1_353, ksg_525, ksg_527, \
                         ksg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_1 * ksf0_350[k]
                   - f_2 * ksf1_350[k]
                   + f_3 * pc_x[k] * ksg_525[k];

        t_736[k] = f_3 * pc_y[k] * ksg_525[k];

        t_737[k] = f_13 * ksf0_352[k]
                   - f_14 * ksf1_352[k]
                   + f_3 * pc_x[k] * ksg_527[k];

        t_738[k] = f_6 * ksf0_353[k]
                   - f_7 * ksf1_353[k]
                   + f_3 * pc_x[k] * ksg_528[k];

        t_739[k] = f_3 * pc_y[k] * ksg_527[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, pc_y, ksf0_355, ksf0_356, ksf0_357, \
                         ksf1_355, ksf1_356, ksf1_357, ksg_530, ksg_531, \
                         ksg_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_6 * ksf0_355[k]
                   - f_7 * ksf1_355[k]
                   + f_3 * pc_x[k] * ksg_530[k];

        t_741[k] = f_4 * ksf0_356[k]
                   - f_5 * ksf1_356[k]
                   + f_3 * pc_x[k] * ksg_531[k];

        t_742[k] = f_4 * ksf0_357[k]
                   - f_5 * ksf1_357[k]
                   + f_3 * pc_x[k] * ksg_532[k];

        t_743[k] = f_3 * pc_y[k] * ksg_530[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, t_749, pc_x, ksf0_359, ksf1_359, \
                         ksg_534, ksg_535, ksg_536, ksg_537, ksg_538, \
                         ksg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_4 * ksf0_359[k]
                   - f_5 * ksf1_359[k]
                   + f_3 * pc_x[k] * ksg_534[k];

        t_745[k] = f_3 * pc_x[k] * ksg_535[k];

        t_746[k] = f_3 * pc_x[k] * ksg_536[k];

        t_747[k] = f_3 * pc_x[k] * ksg_537[k];

        t_748[k] = f_3 * pc_x[k] * ksg_538[k];

        t_749[k] = f_3 * pc_x[k] * ksg_539[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_y, ksf0_356, ksf0_357, ksf0_358, ksf1_356, \
                         ksf1_357, ksf1_358, ksg_535, ksg_536, \
                         ksg_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_1 * ksf0_356[k]
                   - f_2 * ksf1_356[k]
                   + f_3 * pc_y[k] * ksg_535[k];

        t_751[k] = f_13 * ksf0_357[k]
                   - f_14 * ksf1_357[k]
                   + f_3 * pc_y[k] * ksg_536[k];

        t_752[k] = f_6 * ksf0_358[k]
                   - f_7 * ksf1_358[k]
                   + f_3 * pc_y[k] * ksg_537[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pc_y, pc_z, isg_419, ksf0_359, ksf1_359, \
                         ksg_538, ksg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_4 * ksf0_359[k]
                   - f_5 * ksf1_359[k]
                   + f_3 * pc_y[k] * ksg_538[k];

        t_754[k] = f_3 * pc_y[k] * ksg_539[k];

        t_755[k] = f_0 * isg_419[k]
                   + f_1 * ksf0_359[k]
                   - f_2 * ksf1_359[k]
                   + f_3 * pc_z[k] * ksg_539[k];
    }
}

auto
compute_prim_ksh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ish0, const size_t isg,
                                                   const size_t ish1, const size_t ksf0,
                                                   const size_t ksf1, const size_t ksg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);

    compute_prim_ksh_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, ish0, isg,
                                                              ish1, ksf0, ksf1, ksg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
