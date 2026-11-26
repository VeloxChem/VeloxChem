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


#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t hsf0, const size_t hsf1,
                                                          const size_t hsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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

    const auto *gsh0_0 = buffer.data(gsh0 + 0);
    const auto *gsh0_3 = buffer.data(gsh0 + 3);
    const auto *gsh0_5 = buffer.data(gsh0 + 5);
    const auto *gsh0_6 = buffer.data(gsh0 + 6);
    const auto *gsh0_9 = buffer.data(gsh0 + 9);
    const auto *gsh0_15 = buffer.data(gsh0 + 15);
    const auto *gsh0_20 = buffer.data(gsh0 + 20);
    const auto *gsh0_24 = buffer.data(gsh0 + 24);
    const auto *gsh0_27 = buffer.data(gsh0 + 27);
    const auto *gsh0_36 = buffer.data(gsh0 + 36);
    const auto *gsh0_42 = buffer.data(gsh0 + 42);
    const auto *gsh0_47 = buffer.data(gsh0 + 47);
    const auto *gsh0_51 = buffer.data(gsh0 + 51);
    const auto *gsh0_62 = buffer.data(gsh0 + 62);

    const auto *gsg_0 = buffer.data(gsg + 0);
    const auto *gsg_1 = buffer.data(gsg + 1);
    const auto *gsg_2 = buffer.data(gsg + 2);
    const auto *gsg_3 = buffer.data(gsg + 3);
    const auto *gsg_5 = buffer.data(gsg + 5);
    const auto *gsg_10 = buffer.data(gsg + 10);
    const auto *gsg_12 = buffer.data(gsg + 12);
    const auto *gsg_14 = buffer.data(gsg + 14);
    const auto *gsg_15 = buffer.data(gsg + 15);
    const auto *gsg_18 = buffer.data(gsg + 18);
    const auto *gsg_20 = buffer.data(gsg + 20);
    const auto *gsg_25 = buffer.data(gsg + 25);
    const auto *gsg_27 = buffer.data(gsg + 27);
    const auto *gsg_28 = buffer.data(gsg + 28);
    const auto *gsg_29 = buffer.data(gsg + 29);
    const auto *gsg_30 = buffer.data(gsg + 30);
    const auto *gsg_32 = buffer.data(gsg + 32);
    const auto *gsg_35 = buffer.data(gsg + 35);
    const auto *gsg_40 = buffer.data(gsg + 40);
    const auto *gsg_41 = buffer.data(gsg + 41);
    const auto *gsg_42 = buffer.data(gsg + 42);
    const auto *gsg_43 = buffer.data(gsg + 43);
    const auto *gsg_44 = buffer.data(gsg + 44);
    const auto *gsg_45 = buffer.data(gsg + 45);
    const auto *gsg_48 = buffer.data(gsg + 48);
    const auto *gsg_51 = buffer.data(gsg + 51);
    const auto *gsg_55 = buffer.data(gsg + 55);
    const auto *gsg_57 = buffer.data(gsg + 57);
    const auto *gsg_58 = buffer.data(gsg + 58);
    const auto *gsg_59 = buffer.data(gsg + 59);
    const auto *gsg_70 = buffer.data(gsg + 70);
    const auto *gsg_71 = buffer.data(gsg + 71);
    const auto *gsg_72 = buffer.data(gsg + 72);
    const auto *gsg_73 = buffer.data(gsg + 73);
    const auto *gsg_74 = buffer.data(gsg + 74);
    const auto *gsg_75 = buffer.data(gsg + 75);
    const auto *gsg_80 = buffer.data(gsg + 80);
    const auto *gsg_84 = buffer.data(gsg + 84);
    const auto *gsg_85 = buffer.data(gsg + 85);
    const auto *gsg_86 = buffer.data(gsg + 86);
    const auto *gsg_87 = buffer.data(gsg + 87);
    const auto *gsg_89 = buffer.data(gsg + 89);
    const auto *gsg_90 = buffer.data(gsg + 90);
    const auto *gsg_93 = buffer.data(gsg + 93);

    const auto *gsh1_0 = buffer.data(gsh1 + 0);
    const auto *gsh1_3 = buffer.data(gsh1 + 3);
    const auto *gsh1_5 = buffer.data(gsh1 + 5);
    const auto *gsh1_6 = buffer.data(gsh1 + 6);
    const auto *gsh1_9 = buffer.data(gsh1 + 9);
    const auto *gsh1_15 = buffer.data(gsh1 + 15);
    const auto *gsh1_20 = buffer.data(gsh1 + 20);
    const auto *gsh1_24 = buffer.data(gsh1 + 24);
    const auto *gsh1_27 = buffer.data(gsh1 + 27);
    const auto *gsh1_36 = buffer.data(gsh1 + 36);
    const auto *gsh1_42 = buffer.data(gsh1 + 42);
    const auto *gsh1_47 = buffer.data(gsh1 + 47);
    const auto *gsh1_51 = buffer.data(gsh1 + 51);
    const auto *gsh1_62 = buffer.data(gsh1 + 62);

    const auto *hsf0_0 = buffer.data(hsf0 + 0);
    const auto *hsf0_1 = buffer.data(hsf0 + 1);
    const auto *hsf0_2 = buffer.data(hsf0 + 2);
    const auto *hsf0_6 = buffer.data(hsf0 + 6);
    const auto *hsf0_8 = buffer.data(hsf0 + 8);
    const auto *hsf0_9 = buffer.data(hsf0 + 9);
    const auto *hsf0_16 = buffer.data(hsf0 + 16);
    const auto *hsf0_17 = buffer.data(hsf0 + 17);
    const auto *hsf0_22 = buffer.data(hsf0 + 22);
    const auto *hsf0_27 = buffer.data(hsf0 + 27);
    const auto *hsf0_28 = buffer.data(hsf0 + 28);
    const auto *hsf0_29 = buffer.data(hsf0 + 29);
    const auto *hsf0_30 = buffer.data(hsf0 + 30);
    const auto *hsf0_32 = buffer.data(hsf0 + 32);
    const auto *hsf0_33 = buffer.data(hsf0 + 33);
    const auto *hsf0_36 = buffer.data(hsf0 + 36);
    const auto *hsf0_37 = buffer.data(hsf0 + 37);
    const auto *hsf0_39 = buffer.data(hsf0 + 39);
    const auto *hsf0_48 = buffer.data(hsf0 + 48);
    const auto *hsf0_49 = buffer.data(hsf0 + 49);
    const auto *hsf0_50 = buffer.data(hsf0 + 50);
    const auto *hsf0_51 = buffer.data(hsf0 + 51);
    const auto *hsf0_52 = buffer.data(hsf0 + 52);
    const auto *hsf0_55 = buffer.data(hsf0 + 55);
    const auto *hsf0_56 = buffer.data(hsf0 + 56);
    const auto *hsf0_57 = buffer.data(hsf0 + 57);
    const auto *hsf0_58 = buffer.data(hsf0 + 58);
    const auto *hsf0_59 = buffer.data(hsf0 + 59);
    const auto *hsf0_60 = buffer.data(hsf0 + 60);
    const auto *hsf0_63 = buffer.data(hsf0 + 63);

    const auto *hsf1_0 = buffer.data(hsf1 + 0);
    const auto *hsf1_1 = buffer.data(hsf1 + 1);
    const auto *hsf1_2 = buffer.data(hsf1 + 2);
    const auto *hsf1_6 = buffer.data(hsf1 + 6);
    const auto *hsf1_8 = buffer.data(hsf1 + 8);
    const auto *hsf1_9 = buffer.data(hsf1 + 9);
    const auto *hsf1_16 = buffer.data(hsf1 + 16);
    const auto *hsf1_17 = buffer.data(hsf1 + 17);
    const auto *hsf1_22 = buffer.data(hsf1 + 22);
    const auto *hsf1_27 = buffer.data(hsf1 + 27);
    const auto *hsf1_28 = buffer.data(hsf1 + 28);
    const auto *hsf1_29 = buffer.data(hsf1 + 29);
    const auto *hsf1_30 = buffer.data(hsf1 + 30);
    const auto *hsf1_32 = buffer.data(hsf1 + 32);
    const auto *hsf1_33 = buffer.data(hsf1 + 33);
    const auto *hsf1_36 = buffer.data(hsf1 + 36);
    const auto *hsf1_37 = buffer.data(hsf1 + 37);
    const auto *hsf1_39 = buffer.data(hsf1 + 39);
    const auto *hsf1_48 = buffer.data(hsf1 + 48);
    const auto *hsf1_49 = buffer.data(hsf1 + 49);
    const auto *hsf1_50 = buffer.data(hsf1 + 50);
    const auto *hsf1_51 = buffer.data(hsf1 + 51);
    const auto *hsf1_52 = buffer.data(hsf1 + 52);
    const auto *hsf1_55 = buffer.data(hsf1 + 55);
    const auto *hsf1_56 = buffer.data(hsf1 + 56);
    const auto *hsf1_57 = buffer.data(hsf1 + 57);
    const auto *hsf1_58 = buffer.data(hsf1 + 58);
    const auto *hsf1_59 = buffer.data(hsf1 + 59);
    const auto *hsf1_60 = buffer.data(hsf1 + 60);
    const auto *hsf1_63 = buffer.data(hsf1 + 63);

    const auto *hsg_0 = buffer.data(hsg + 0);
    const auto *hsg_1 = buffer.data(hsg + 1);
    const auto *hsg_2 = buffer.data(hsg + 2);
    const auto *hsg_3 = buffer.data(hsg + 3);
    const auto *hsg_5 = buffer.data(hsg + 5);
    const auto *hsg_6 = buffer.data(hsg + 6);
    const auto *hsg_9 = buffer.data(hsg + 9);
    const auto *hsg_10 = buffer.data(hsg + 10);
    const auto *hsg_12 = buffer.data(hsg + 12);
    const auto *hsg_13 = buffer.data(hsg + 13);
    const auto *hsg_14 = buffer.data(hsg + 14);
    const auto *hsg_15 = buffer.data(hsg + 15);
    const auto *hsg_16 = buffer.data(hsg + 16);
    const auto *hsg_18 = buffer.data(hsg + 18);
    const auto *hsg_20 = buffer.data(hsg + 20);
    const auto *hsg_21 = buffer.data(hsg + 21);
    const auto *hsg_25 = buffer.data(hsg + 25);
    const auto *hsg_26 = buffer.data(hsg + 26);
    const auto *hsg_27 = buffer.data(hsg + 27);
    const auto *hsg_28 = buffer.data(hsg + 28);
    const auto *hsg_29 = buffer.data(hsg + 29);
    const auto *hsg_30 = buffer.data(hsg + 30);
    const auto *hsg_32 = buffer.data(hsg + 32);
    const auto *hsg_34 = buffer.data(hsg + 34);
    const auto *hsg_35 = buffer.data(hsg + 35);
    const auto *hsg_39 = buffer.data(hsg + 39);
    const auto *hsg_40 = buffer.data(hsg + 40);
    const auto *hsg_41 = buffer.data(hsg + 41);
    const auto *hsg_42 = buffer.data(hsg + 42);
    const auto *hsg_43 = buffer.data(hsg + 43);
    const auto *hsg_44 = buffer.data(hsg + 44);
    const auto *hsg_45 = buffer.data(hsg + 45);
    const auto *hsg_46 = buffer.data(hsg + 46);
    const auto *hsg_47 = buffer.data(hsg + 47);
    const auto *hsg_48 = buffer.data(hsg + 48);
    const auto *hsg_50 = buffer.data(hsg + 50);
    const auto *hsg_51 = buffer.data(hsg + 51);
    const auto *hsg_55 = buffer.data(hsg + 55);
    const auto *hsg_56 = buffer.data(hsg + 56);
    const auto *hsg_57 = buffer.data(hsg + 57);
    const auto *hsg_58 = buffer.data(hsg + 58);
    const auto *hsg_59 = buffer.data(hsg + 59);
    const auto *hsg_60 = buffer.data(hsg + 60);
    const auto *hsg_62 = buffer.data(hsg + 62);
    const auto *hsg_63 = buffer.data(hsg + 63);
    const auto *hsg_65 = buffer.data(hsg + 65);
    const auto *hsg_70 = buffer.data(hsg + 70);
    const auto *hsg_71 = buffer.data(hsg + 71);
    const auto *hsg_72 = buffer.data(hsg + 72);
    const auto *hsg_73 = buffer.data(hsg + 73);
    const auto *hsg_74 = buffer.data(hsg + 74);
    const auto *hsg_75 = buffer.data(hsg + 75);
    const auto *hsg_76 = buffer.data(hsg + 76);
    const auto *hsg_77 = buffer.data(hsg + 77);
    const auto *hsg_78 = buffer.data(hsg + 78);
    const auto *hsg_79 = buffer.data(hsg + 79);
    const auto *hsg_80 = buffer.data(hsg + 80);
    const auto *hsg_84 = buffer.data(hsg + 84);
    const auto *hsg_85 = buffer.data(hsg + 85);
    const auto *hsg_86 = buffer.data(hsg + 86);
    const auto *hsg_87 = buffer.data(hsg + 87);
    const auto *hsg_88 = buffer.data(hsg + 88);
    const auto *hsg_89 = buffer.data(hsg + 89);
    const auto *hsg_90 = buffer.data(hsg + 90);
    const auto *hsg_93 = buffer.data(hsg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsg_0, hsf0_0, \
                         hsf1_0, hsg_0, hsg_1, hsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsg_0[k]
                 + f_1 * hsf0_0[k]
                 - f_2 * hsf1_0[k]
                 + f_3 * pc_x[k] * hsg_0[k];

        t_1[k] = f_3 * pc_y[k] * hsg_0[k];

        t_2[k] = f_3 * pc_z[k] * hsg_0[k];

        t_3[k] = f_4 * hsf0_0[k]
                 - f_5 * hsf1_0[k]
                 + f_3 * pc_y[k] * hsg_1[k];

        t_4[k] = f_3 * pc_y[k] * hsg_2[k];

        t_5[k] = f_4 * hsf0_0[k]
                 - f_5 * hsf1_0[k]
                 + f_3 * pc_z[k] * hsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, gsg_10, hsf0_1, hsf0_2, \
                         hsf1_1, hsf1_2, hsg_3, hsg_5, hsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hsf0_1[k]
                 - f_7 * hsf1_1[k]
                 + f_3 * pc_y[k] * hsg_3[k];

        t_7[k] = f_3 * pc_z[k] * hsg_3[k];

        t_8[k] = f_3 * pc_y[k] * hsg_5[k];

        t_9[k] = f_6 * hsf0_2[k]
                 - f_7 * hsf1_2[k]
                 + f_3 * pc_z[k] * hsg_5[k];

        t_10[k] = f_0 * gsg_10[k]
                  + f_3 * pc_x[k] * hsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, gsg_12, gsg_14, hsg_6, \
                         hsg_9, hsg_12, hsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * hsg_6[k];

        t_12[k] = f_0 * gsg_12[k]
                  + f_3 * pc_x[k] * hsg_12[k];

        t_13[k] = f_3 * pc_y[k] * hsg_9[k];

        t_14[k] = f_0 * gsg_14[k]
                  + f_3 * pc_x[k] * hsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, hsf0_6, hsf0_8, hsf0_9, hsf1_6, \
                         hsf1_8, hsf1_9, hsg_10, hsg_12, hsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * hsf0_6[k]
                  - f_2 * hsf1_6[k]
                  + f_3 * pc_y[k] * hsg_10[k];

        t_16[k] = f_3 * pc_z[k] * hsg_10[k];

        t_17[k] = f_6 * hsf0_8[k]
                  - f_7 * hsf1_8[k]
                  + f_3 * pc_y[k] * hsg_12[k];

        t_18[k] = f_4 * hsf0_9[k]
                  - f_5 * hsf1_9[k]
                  + f_3 * pc_y[k] * hsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, gsh0_0, gsg_0, \
                         gsh1_0, hsf0_9, hsf1_9, hsg_14, hsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * hsg_14[k];

        t_20[k] = f_1 * hsf0_9[k]
                  - f_2 * hsf1_9[k]
                  + f_3 * pc_z[k] * hsg_14[k];

        t_21[k] = pa_y[k] * gsh0_0[k]
                  - f_8 * pc_y[k] * gsh1_0[k];

        t_22[k] = f_9 * gsg_0[k]
                  + f_3 * pc_y[k] * hsg_15[k];

        t_23[k] = f_3 * pc_z[k] * hsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, gsh0_3, gsh0_5, gsh0_6, \
                         gsg_1, gsg_3, gsh1_3, gsh1_5, gsh1_6, hsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * gsh0_3[k]
                  + f_10 * gsg_1[k]
                  - f_8 * pc_y[k] * gsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * hsg_16[k];

        t_26[k] = pa_y[k] * gsh0_5[k]
                  - f_8 * pc_y[k] * gsh1_5[k];

        t_27[k] = pa_y[k] * gsh0_6[k]
                  + f_11 * gsg_3[k]
                  - f_8 * pc_y[k] * gsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, gsh0_9, gsg_5, \
                         gsg_25, gsh1_9, hsg_18, hsg_20, hsg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * hsg_18[k];

        t_29[k] = f_9 * gsg_5[k]
                  + f_3 * pc_y[k] * hsg_20[k];

        t_30[k] = pa_y[k] * gsh0_9[k]
                  - f_8 * pc_y[k] * gsh1_9[k];

        t_31[k] = f_12 * gsg_25[k]
                  + f_3 * pc_x[k] * hsg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, gsg_27, gsg_28, gsg_29, hsg_21, \
                         hsg_27, hsg_28, hsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * hsg_21[k];

        t_33[k] = f_12 * gsg_27[k]
                  + f_3 * pc_x[k] * hsg_27[k];

        t_34[k] = f_12 * gsg_28[k]
                  + f_3 * pc_x[k] * hsg_28[k];

        t_35[k] = f_12 * gsg_29[k]
                  + f_3 * pc_x[k] * hsg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, gsg_10, hsf0_16, hsf0_17, \
                         hsf1_16, hsf1_17, hsg_25, hsg_26, hsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * gsg_10[k]
                  + f_1 * hsf0_16[k]
                  - f_2 * hsf1_16[k]
                  + f_3 * pc_y[k] * hsg_25[k];

        t_37[k] = f_3 * pc_z[k] * hsg_25[k];

        t_38[k] = f_4 * hsf0_16[k]
                  - f_5 * hsf1_16[k]
                  + f_3 * pc_z[k] * hsg_26[k];

        t_39[k] = f_6 * hsf0_17[k]
                  - f_7 * hsf1_17[k]
                  + f_3 * pc_z[k] * hsg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, gsh0_0, gsh0_20, \
                         gsg_14, gsh1_0, gsh1_20, hsg_29, hsg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * gsg_14[k]
                  + f_3 * pc_y[k] * hsg_29[k];

        t_41[k] = pa_y[k] * gsh0_20[k]
                  - f_8 * pc_y[k] * gsh1_20[k];

        t_42[k] = pa_z[k] * gsh0_0[k]
                  - f_8 * pc_z[k] * gsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * hsg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, gsh0_3, gsh0_5, gsg_0, \
                         gsg_2, gsh1_3, gsh1_5, hsg_30, hsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * gsg_0[k]
                  + f_3 * pc_z[k] * hsg_30[k];

        t_45[k] = pa_z[k] * gsh0_3[k]
                  - f_8 * pc_z[k] * gsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * hsg_32[k];

        t_47[k] = pa_z[k] * gsh0_5[k]
                  + f_10 * gsg_2[k]
                  - f_8 * pc_z[k] * gsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, gsh0_6, gsh0_9, gsg_5, \
                         gsh1_6, gsh1_9, hsf0_22, hsf1_22, hsg_34, \
                         hsg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * gsh0_6[k]
                  - f_8 * pc_z[k] * gsh1_6[k];

        t_49[k] = f_4 * hsf0_22[k]
                  - f_5 * hsf1_22[k]
                  + f_3 * pc_y[k] * hsg_34[k];

        t_50[k] = f_3 * pc_y[k] * hsg_35[k];

        t_51[k] = pa_z[k] * gsh0_9[k]
                  + f_11 * gsg_5[k]
                  - f_8 * pc_z[k] * gsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, gsg_40, gsg_41, gsg_42, \
                         gsg_44, hsg_39, hsg_40, hsg_41, hsg_42, \
                         hsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * gsg_40[k]
                  + f_3 * pc_x[k] * hsg_40[k];

        t_53[k] = f_12 * gsg_41[k]
                  + f_3 * pc_x[k] * hsg_41[k];

        t_54[k] = f_12 * gsg_42[k]
                  + f_3 * pc_x[k] * hsg_42[k];

        t_55[k] = f_3 * pc_y[k] * hsg_39[k];

        t_56[k] = f_12 * gsg_44[k]
                  + f_3 * pc_x[k] * hsg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, gsh0_15, gsh1_15, hsf0_27, \
                         hsf0_28, hsf1_27, hsf1_28, hsg_41, hsg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * gsh0_15[k]
                  - f_8 * pc_z[k] * gsh1_15[k];

        t_58[k] = f_13 * hsf0_27[k]
                  - f_14 * hsf1_27[k]
                  + f_3 * pc_y[k] * hsg_41[k];

        t_59[k] = f_6 * hsf0_28[k]
                  - f_7 * hsf1_28[k]
                  + f_3 * pc_y[k] * hsg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, gsg_14, gsg_45, hsf0_29, \
                         hsf0_30, hsf1_29, hsf1_30, hsg_43, hsg_44, \
                         hsg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * hsf0_29[k]
                  - f_5 * hsf1_29[k]
                  + f_3 * pc_y[k] * hsg_43[k];

        t_61[k] = f_3 * pc_y[k] * hsg_44[k];

        t_62[k] = f_9 * gsg_14[k]
                  + f_1 * hsf0_29[k]
                  - f_2 * hsf1_29[k]
                  + f_3 * pc_z[k] * hsg_44[k];

        t_63[k] = f_11 * gsg_45[k]
                  + f_1 * hsf0_30[k]
                  - f_2 * hsf1_30[k]
                  + f_3 * pc_x[k] * hsg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, gsg_15, gsg_48, hsf0_33, \
                         hsf1_33, hsg_45, hsg_46, hsg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * gsg_15[k]
                  + f_3 * pc_y[k] * hsg_45[k];

        t_65[k] = f_3 * pc_z[k] * hsg_45[k];

        t_66[k] = f_11 * gsg_48[k]
                  + f_6 * hsf0_33[k]
                  - f_7 * hsf1_33[k]
                  + f_3 * pc_x[k] * hsg_48[k];

        t_67[k] = f_3 * pc_z[k] * hsg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, gsg_51, hsf0_30, hsf0_36, hsf1_30, \
                         hsf1_36, hsg_47, hsg_48, hsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * hsf0_30[k]
                  - f_5 * hsf1_30[k]
                  + f_3 * pc_z[k] * hsg_47[k];

        t_69[k] = f_11 * gsg_51[k]
                  + f_4 * hsf0_36[k]
                  - f_5 * hsf1_36[k]
                  + f_3 * pc_x[k] * hsg_51[k];

        t_70[k] = f_3 * pc_z[k] * hsg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, gsg_20, gsg_55, hsf0_32, \
                         hsf1_32, hsg_50, hsg_51, hsg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * gsg_20[k]
                  + f_3 * pc_y[k] * hsg_50[k];

        t_72[k] = f_6 * hsf0_32[k]
                  - f_7 * hsf1_32[k]
                  + f_3 * pc_z[k] * hsg_50[k];

        t_73[k] = f_11 * gsg_55[k]
                  + f_3 * pc_x[k] * hsg_55[k];

        t_74[k] = f_3 * pc_z[k] * hsg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, gsg_25, gsg_57, gsg_58, gsg_59, \
                         hsf0_36, hsf1_36, hsg_55, hsg_57, hsg_58, \
                         hsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_11 * gsg_57[k]
                  + f_3 * pc_x[k] * hsg_57[k];

        t_76[k] = f_11 * gsg_58[k]
                  + f_3 * pc_x[k] * hsg_58[k];

        t_77[k] = f_11 * gsg_59[k]
                  + f_3 * pc_x[k] * hsg_59[k];

        t_78[k] = f_10 * gsg_25[k]
                  + f_1 * hsf0_36[k]
                  - f_2 * hsf1_36[k]
                  + f_3 * pc_y[k] * hsg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, gsg_29, hsf0_36, hsf0_37, \
                         hsf1_36, hsf1_37, hsg_55, hsg_56, hsg_57, \
                         hsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * hsg_55[k];

        t_80[k] = f_4 * hsf0_36[k]
                  - f_5 * hsf1_36[k]
                  + f_3 * pc_z[k] * hsg_56[k];

        t_81[k] = f_6 * hsf0_37[k]
                  - f_7 * hsf1_37[k]
                  + f_3 * pc_z[k] * hsg_57[k];

        t_82[k] = f_10 * gsg_29[k]
                  + f_3 * pc_y[k] * hsg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, gsh0_42, gsg_15, gsg_30, \
                         gsh1_42, hsf0_39, hsf1_39, hsg_59, hsg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * hsf0_39[k]
                  - f_2 * hsf1_39[k]
                  + f_3 * pc_z[k] * hsg_59[k];

        t_84[k] = pa_y[k] * gsh0_42[k]
                  - f_8 * pc_y[k] * gsh1_42[k];

        t_85[k] = f_9 * gsg_30[k]
                  + f_3 * pc_y[k] * hsg_60[k];

        t_86[k] = f_9 * gsg_15[k]
                  + f_3 * pc_z[k] * hsg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, gsh0_24, gsh0_27, \
                         gsh0_47, gsg_32, gsh1_24, gsh1_27, gsh1_47, \
                         hsg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * gsh0_24[k]
                  - f_8 * pc_z[k] * gsh1_24[k];

        t_88[k] = f_9 * gsg_32[k]
                  + f_3 * pc_y[k] * hsg_62[k];

        t_89[k] = pa_y[k] * gsh0_47[k]
                  - f_8 * pc_y[k] * gsh1_47[k];

        t_90[k] = pa_z[k] * gsh0_27[k]
                  - f_8 * pc_z[k] * gsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, gsh0_51, gsg_18, \
                         gsg_35, gsg_70, gsh1_51, hsg_63, hsg_65, \
                         hsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * gsg_18[k]
                  + f_3 * pc_z[k] * hsg_63[k];

        t_92[k] = f_9 * gsg_35[k]
                  + f_3 * pc_y[k] * hsg_65[k];

        t_93[k] = pa_y[k] * gsh0_51[k]
                  - f_8 * pc_y[k] * gsh1_51[k];

        t_94[k] = f_11 * gsg_70[k]
                  + f_3 * pc_x[k] * hsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, gsg_71, gsg_72, gsg_73, gsg_74, hsg_71, \
                         hsg_72, hsg_73, hsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * gsg_71[k]
                  + f_3 * pc_x[k] * hsg_71[k];

        t_96[k] = f_11 * gsg_72[k]
                  + f_3 * pc_x[k] * hsg_72[k];

        t_97[k] = f_11 * gsg_73[k]
                  + f_3 * pc_x[k] * hsg_73[k];

        t_98[k] = f_11 * gsg_74[k]
                  + f_3 * pc_x[k] * hsg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, gsh0_36, gsg_25, gsg_42, \
                         gsh1_36, hsf0_48, hsf1_48, hsg_70, hsg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * gsh0_36[k]
                  - f_8 * pc_z[k] * gsh1_36[k];

        t_100[k] = f_9 * gsg_25[k]
                   + f_3 * pc_z[k] * hsg_70[k];

        t_101[k] = f_9 * gsg_42[k]
                   + f_6 * hsf0_48[k]
                   - f_7 * hsf1_48[k]
                   + f_3 * pc_y[k] * hsg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, gsh0_62, gsg_43, gsg_44, gsh1_62, \
                         hsf0_49, hsf1_49, hsg_73, hsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * gsg_43[k]
                   + f_4 * hsf0_49[k]
                   - f_5 * hsf1_49[k]
                   + f_3 * pc_y[k] * hsg_73[k];

        t_103[k] = f_9 * gsg_44[k]
                   + f_3 * pc_y[k] * hsg_74[k];

        t_104[k] = pa_y[k] * gsh0_62[k]
                   - f_8 * pc_y[k] * gsh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, gsg_30, gsg_75, \
                         hsf0_50, hsf1_50, hsg_75, hsg_76, hsg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_11 * gsg_75[k]
                   + f_1 * hsf0_50[k]
                   - f_2 * hsf1_50[k]
                   + f_3 * pc_x[k] * hsg_75[k];

        t_106[k] = f_3 * pc_y[k] * hsg_75[k];

        t_107[k] = f_10 * gsg_30[k]
                   + f_3 * pc_z[k] * hsg_75[k];

        t_108[k] = f_4 * hsf0_50[k]
                   - f_5 * hsf1_50[k]
                   + f_3 * pc_y[k] * hsg_76[k];

        t_109[k] = f_3 * pc_y[k] * hsg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, gsg_80, hsf0_51, hsf0_52, \
                         hsf0_55, hsf1_51, hsf1_52, hsf1_55, hsg_78, hsg_79, \
                         hsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * gsg_80[k]
                   + f_6 * hsf0_55[k]
                   - f_7 * hsf1_55[k]
                   + f_3 * pc_x[k] * hsg_80[k];

        t_111[k] = f_6 * hsf0_51[k]
                   - f_7 * hsf1_51[k]
                   + f_3 * pc_y[k] * hsg_78[k];

        t_112[k] = f_4 * hsf0_52[k]
                   - f_5 * hsf1_52[k]
                   + f_3 * pc_y[k] * hsg_79[k];

        t_113[k] = f_3 * pc_y[k] * hsg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, gsg_84, gsg_85, gsg_86, gsg_87, \
                         hsf0_59, hsf1_59, hsg_84, hsg_85, hsg_86, \
                         hsg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * gsg_84[k]
                   + f_4 * hsf0_59[k]
                   - f_5 * hsf1_59[k]
                   + f_3 * pc_x[k] * hsg_84[k];

        t_115[k] = f_11 * gsg_85[k]
                   + f_3 * pc_x[k] * hsg_85[k];

        t_116[k] = f_11 * gsg_86[k]
                   + f_3 * pc_x[k] * hsg_86[k];

        t_117[k] = f_11 * gsg_87[k]
                   + f_3 * pc_x[k] * hsg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, gsg_89, hsf0_56, hsf0_57, \
                         hsf1_56, hsf1_57, hsg_84, hsg_85, hsg_86, \
                         hsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * hsg_84[k];

        t_119[k] = f_11 * gsg_89[k]
                   + f_3 * pc_x[k] * hsg_89[k];

        t_120[k] = f_1 * hsf0_56[k]
                   - f_2 * hsf1_56[k]
                   + f_3 * pc_y[k] * hsg_85[k];

        t_121[k] = f_13 * hsf0_57[k]
                   - f_14 * hsf1_57[k]
                   + f_3 * pc_y[k] * hsg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, gsg_44, hsf0_58, hsf0_59, \
                         hsf1_58, hsf1_59, hsg_87, hsg_88, hsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * hsf0_58[k]
                   - f_7 * hsf1_58[k]
                   + f_3 * pc_y[k] * hsg_87[k];

        t_123[k] = f_4 * hsf0_59[k]
                   - f_5 * hsf1_59[k]
                   + f_3 * pc_y[k] * hsg_88[k];

        t_124[k] = f_3 * pc_y[k] * hsg_89[k];

        t_125[k] = f_10 * gsg_44[k]
                   + f_1 * hsf0_59[k]
                   - f_2 * hsf1_59[k]
                   + f_3 * pc_z[k] * hsg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, gsg_45, gsg_90, gsg_93, \
                         hsf0_60, hsf0_63, hsf1_60, hsf1_63, hsg_90, \
                         hsg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_10 * gsg_90[k]
                   + f_1 * hsf0_60[k]
                   - f_2 * hsf1_60[k]
                   + f_3 * pc_x[k] * hsg_90[k];

        t_127[k] = f_11 * gsg_45[k]
                   + f_3 * pc_y[k] * hsg_90[k];

        t_128[k] = f_3 * pc_z[k] * hsg_90[k];

        t_129[k] = f_10 * gsg_93[k]
                   + f_6 * hsf0_63[k]
                   - f_7 * hsf1_63[k]
                   + f_3 * pc_x[k] * hsg_93[k];
    }
}

static auto
compute_prim_hsh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t hsf0, const size_t hsf1,
                                                          const size_t hsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsh0_63 = buffer.data(gsh0 + 63);
    const auto *gsh0_66 = buffer.data(gsh0 + 66);
    const auto *gsh0_69 = buffer.data(gsh0 + 69);
    const auto *gsh0_78 = buffer.data(gsh0 + 78);
    const auto *gsh0_105 = buffer.data(gsh0 + 105);
    const auto *gsh0_108 = buffer.data(gsh0 + 108);
    const auto *gsh0_110 = buffer.data(gsh0 + 110);
    const auto *gsh0_111 = buffer.data(gsh0 + 111);
    const auto *gsh0_114 = buffer.data(gsh0 + 114);
    const auto *gsh0_125 = buffer.data(gsh0 + 125);
    const auto *gsh0_126 = buffer.data(gsh0 + 126);
    const auto *gsh0_129 = buffer.data(gsh0 + 129);
    const auto *gsh0_132 = buffer.data(gsh0 + 132);
    const auto *gsh0_210 = buffer.data(gsh0 + 210);
    const auto *gsh0_213 = buffer.data(gsh0 + 213);
    const auto *gsh0_216 = buffer.data(gsh0 + 216);
    const auto *gsh0_225 = buffer.data(gsh0 + 225);
    const auto *gsh0_227 = buffer.data(gsh0 + 227);
    const auto *gsh0_228 = buffer.data(gsh0 + 228);
    const auto *gsh0_230 = buffer.data(gsh0 + 230);
    const auto *gsh0_236 = buffer.data(gsh0 + 236);
    const auto *gsh0_240 = buffer.data(gsh0 + 240);
    const auto *gsh0_246 = buffer.data(gsh0 + 246);
    const auto *gsh0_248 = buffer.data(gsh0 + 248);
    const auto *gsh0_249 = buffer.data(gsh0 + 249);
    const auto *gsh0_251 = buffer.data(gsh0 + 251);
    const auto *gsh0_252 = buffer.data(gsh0 + 252);

    const auto *gsg_45 = buffer.data(gsg + 45);
    const auto *gsg_48 = buffer.data(gsg + 48);
    const auto *gsg_50 = buffer.data(gsg + 50);
    const auto *gsg_55 = buffer.data(gsg + 55);
    const auto *gsg_59 = buffer.data(gsg + 59);
    const auto *gsg_60 = buffer.data(gsg + 60);
    const auto *gsg_62 = buffer.data(gsg + 62);
    const auto *gsg_63 = buffer.data(gsg + 63);
    const auto *gsg_65 = buffer.data(gsg + 65);
    const auto *gsg_70 = buffer.data(gsg + 70);
    const auto *gsg_72 = buffer.data(gsg + 72);
    const auto *gsg_73 = buffer.data(gsg + 73);
    const auto *gsg_74 = buffer.data(gsg + 74);
    const auto *gsg_75 = buffer.data(gsg + 75);
    const auto *gsg_76 = buffer.data(gsg + 76);
    const auto *gsg_77 = buffer.data(gsg + 77);
    const auto *gsg_78 = buffer.data(gsg + 78);
    const auto *gsg_80 = buffer.data(gsg + 80);
    const auto *gsg_85 = buffer.data(gsg + 85);
    const auto *gsg_87 = buffer.data(gsg + 87);
    const auto *gsg_88 = buffer.data(gsg + 88);
    const auto *gsg_89 = buffer.data(gsg + 89);
    const auto *gsg_90 = buffer.data(gsg + 90);
    const auto *gsg_93 = buffer.data(gsg + 93);
    const auto *gsg_95 = buffer.data(gsg + 95);
    const auto *gsg_96 = buffer.data(gsg + 96);
    const auto *gsg_100 = buffer.data(gsg + 100);
    const auto *gsg_102 = buffer.data(gsg + 102);
    const auto *gsg_103 = buffer.data(gsg + 103);
    const auto *gsg_104 = buffer.data(gsg + 104);
    const auto *gsg_105 = buffer.data(gsg + 105);
    const auto *gsg_107 = buffer.data(gsg + 107);
    const auto *gsg_110 = buffer.data(gsg + 110);
    const auto *gsg_114 = buffer.data(gsg + 114);
    const auto *gsg_115 = buffer.data(gsg + 115);
    const auto *gsg_116 = buffer.data(gsg + 116);
    const auto *gsg_117 = buffer.data(gsg + 117);
    const auto *gsg_118 = buffer.data(gsg + 118);
    const auto *gsg_119 = buffer.data(gsg + 119);
    const auto *gsg_120 = buffer.data(gsg + 120);
    const auto *gsg_130 = buffer.data(gsg + 130);
    const auto *gsg_131 = buffer.data(gsg + 131);
    const auto *gsg_132 = buffer.data(gsg + 132);
    const auto *gsg_133 = buffer.data(gsg + 133);
    const auto *gsg_134 = buffer.data(gsg + 134);
    const auto *gsg_135 = buffer.data(gsg + 135);
    const auto *gsg_140 = buffer.data(gsg + 140);
    const auto *gsg_144 = buffer.data(gsg + 144);
    const auto *gsg_145 = buffer.data(gsg + 145);
    const auto *gsg_146 = buffer.data(gsg + 146);
    const auto *gsg_147 = buffer.data(gsg + 147);
    const auto *gsg_149 = buffer.data(gsg + 149);
    const auto *gsg_150 = buffer.data(gsg + 150);
    const auto *gsg_153 = buffer.data(gsg + 153);
    const auto *gsg_156 = buffer.data(gsg + 156);
    const auto *gsg_160 = buffer.data(gsg + 160);
    const auto *gsg_162 = buffer.data(gsg + 162);
    const auto *gsg_163 = buffer.data(gsg + 163);
    const auto *gsg_164 = buffer.data(gsg + 164);
    const auto *gsg_170 = buffer.data(gsg + 170);
    const auto *gsg_174 = buffer.data(gsg + 174);
    const auto *gsg_175 = buffer.data(gsg + 175);
    const auto *gsg_176 = buffer.data(gsg + 176);
    const auto *gsg_177 = buffer.data(gsg + 177);
    const auto *gsg_178 = buffer.data(gsg + 178);
    const auto *gsg_179 = buffer.data(gsg + 179);
    const auto *gsg_180 = buffer.data(gsg + 180);

    const auto *gsh1_63 = buffer.data(gsh1 + 63);
    const auto *gsh1_66 = buffer.data(gsh1 + 66);
    const auto *gsh1_69 = buffer.data(gsh1 + 69);
    const auto *gsh1_78 = buffer.data(gsh1 + 78);
    const auto *gsh1_105 = buffer.data(gsh1 + 105);
    const auto *gsh1_108 = buffer.data(gsh1 + 108);
    const auto *gsh1_110 = buffer.data(gsh1 + 110);
    const auto *gsh1_111 = buffer.data(gsh1 + 111);
    const auto *gsh1_114 = buffer.data(gsh1 + 114);
    const auto *gsh1_125 = buffer.data(gsh1 + 125);
    const auto *gsh1_126 = buffer.data(gsh1 + 126);
    const auto *gsh1_129 = buffer.data(gsh1 + 129);
    const auto *gsh1_132 = buffer.data(gsh1 + 132);
    const auto *gsh1_210 = buffer.data(gsh1 + 210);
    const auto *gsh1_213 = buffer.data(gsh1 + 213);
    const auto *gsh1_216 = buffer.data(gsh1 + 216);
    const auto *gsh1_225 = buffer.data(gsh1 + 225);
    const auto *gsh1_227 = buffer.data(gsh1 + 227);
    const auto *gsh1_228 = buffer.data(gsh1 + 228);
    const auto *gsh1_230 = buffer.data(gsh1 + 230);
    const auto *gsh1_236 = buffer.data(gsh1 + 236);
    const auto *gsh1_240 = buffer.data(gsh1 + 240);
    const auto *gsh1_246 = buffer.data(gsh1 + 246);
    const auto *gsh1_248 = buffer.data(gsh1 + 248);
    const auto *gsh1_249 = buffer.data(gsh1 + 249);
    const auto *gsh1_251 = buffer.data(gsh1 + 251);
    const auto *gsh1_252 = buffer.data(gsh1 + 252);

    const auto *hsf0_60 = buffer.data(hsf0 + 60);
    const auto *hsf0_62 = buffer.data(hsf0 + 62);
    const auto *hsf0_66 = buffer.data(hsf0 + 66);
    const auto *hsf0_67 = buffer.data(hsf0 + 67);
    const auto *hsf0_69 = buffer.data(hsf0 + 69);
    const auto *hsf0_75 = buffer.data(hsf0 + 75);
    const auto *hsf0_78 = buffer.data(hsf0 + 78);
    const auto *hsf0_79 = buffer.data(hsf0 + 79);
    const auto *hsf0_86 = buffer.data(hsf0 + 86);
    const auto *hsf0_88 = buffer.data(hsf0 + 88);
    const auto *hsf0_89 = buffer.data(hsf0 + 89);
    const auto *hsf0_90 = buffer.data(hsf0 + 90);
    const auto *hsf0_91 = buffer.data(hsf0 + 91);
    const auto *hsf0_92 = buffer.data(hsf0 + 92);
    const auto *hsf0_95 = buffer.data(hsf0 + 95);
    const auto *hsf0_96 = buffer.data(hsf0 + 96);
    const auto *hsf0_97 = buffer.data(hsf0 + 97);
    const auto *hsf0_98 = buffer.data(hsf0 + 98);
    const auto *hsf0_99 = buffer.data(hsf0 + 99);
    const auto *hsf0_100 = buffer.data(hsf0 + 100);
    const auto *hsf0_102 = buffer.data(hsf0 + 102);

    const auto *hsf1_60 = buffer.data(hsf1 + 60);
    const auto *hsf1_62 = buffer.data(hsf1 + 62);
    const auto *hsf1_66 = buffer.data(hsf1 + 66);
    const auto *hsf1_67 = buffer.data(hsf1 + 67);
    const auto *hsf1_69 = buffer.data(hsf1 + 69);
    const auto *hsf1_75 = buffer.data(hsf1 + 75);
    const auto *hsf1_78 = buffer.data(hsf1 + 78);
    const auto *hsf1_79 = buffer.data(hsf1 + 79);
    const auto *hsf1_86 = buffer.data(hsf1 + 86);
    const auto *hsf1_88 = buffer.data(hsf1 + 88);
    const auto *hsf1_89 = buffer.data(hsf1 + 89);
    const auto *hsf1_90 = buffer.data(hsf1 + 90);
    const auto *hsf1_91 = buffer.data(hsf1 + 91);
    const auto *hsf1_92 = buffer.data(hsf1 + 92);
    const auto *hsf1_95 = buffer.data(hsf1 + 95);
    const auto *hsf1_96 = buffer.data(hsf1 + 96);
    const auto *hsf1_97 = buffer.data(hsf1 + 97);
    const auto *hsf1_98 = buffer.data(hsf1 + 98);
    const auto *hsf1_99 = buffer.data(hsf1 + 99);
    const auto *hsf1_100 = buffer.data(hsf1 + 100);
    const auto *hsf1_102 = buffer.data(hsf1 + 102);

    const auto *hsg_91 = buffer.data(hsg + 91);
    const auto *hsg_92 = buffer.data(hsg + 92);
    const auto *hsg_93 = buffer.data(hsg + 93);
    const auto *hsg_95 = buffer.data(hsg + 95);
    const auto *hsg_96 = buffer.data(hsg + 96);
    const auto *hsg_100 = buffer.data(hsg + 100);
    const auto *hsg_101 = buffer.data(hsg + 101);
    const auto *hsg_102 = buffer.data(hsg + 102);
    const auto *hsg_103 = buffer.data(hsg + 103);
    const auto *hsg_104 = buffer.data(hsg + 104);
    const auto *hsg_105 = buffer.data(hsg + 105);
    const auto *hsg_107 = buffer.data(hsg + 107);
    const auto *hsg_108 = buffer.data(hsg + 108);
    const auto *hsg_110 = buffer.data(hsg + 110);
    const auto *hsg_114 = buffer.data(hsg + 114);
    const auto *hsg_115 = buffer.data(hsg + 115);
    const auto *hsg_116 = buffer.data(hsg + 116);
    const auto *hsg_117 = buffer.data(hsg + 117);
    const auto *hsg_118 = buffer.data(hsg + 118);
    const auto *hsg_119 = buffer.data(hsg + 119);
    const auto *hsg_120 = buffer.data(hsg + 120);
    const auto *hsg_122 = buffer.data(hsg + 122);
    const auto *hsg_123 = buffer.data(hsg + 123);
    const auto *hsg_125 = buffer.data(hsg + 125);
    const auto *hsg_130 = buffer.data(hsg + 130);
    const auto *hsg_131 = buffer.data(hsg + 131);
    const auto *hsg_132 = buffer.data(hsg + 132);
    const auto *hsg_133 = buffer.data(hsg + 133);
    const auto *hsg_134 = buffer.data(hsg + 134);
    const auto *hsg_135 = buffer.data(hsg + 135);
    const auto *hsg_136 = buffer.data(hsg + 136);
    const auto *hsg_137 = buffer.data(hsg + 137);
    const auto *hsg_138 = buffer.data(hsg + 138);
    const auto *hsg_139 = buffer.data(hsg + 139);
    const auto *hsg_140 = buffer.data(hsg + 140);
    const auto *hsg_144 = buffer.data(hsg + 144);
    const auto *hsg_145 = buffer.data(hsg + 145);
    const auto *hsg_146 = buffer.data(hsg + 146);
    const auto *hsg_147 = buffer.data(hsg + 147);
    const auto *hsg_148 = buffer.data(hsg + 148);
    const auto *hsg_149 = buffer.data(hsg + 149);
    const auto *hsg_150 = buffer.data(hsg + 150);
    const auto *hsg_151 = buffer.data(hsg + 151);
    const auto *hsg_152 = buffer.data(hsg + 152);
    const auto *hsg_153 = buffer.data(hsg + 153);
    const auto *hsg_155 = buffer.data(hsg + 155);
    const auto *hsg_156 = buffer.data(hsg + 156);
    const auto *hsg_160 = buffer.data(hsg + 160);
    const auto *hsg_162 = buffer.data(hsg + 162);
    const auto *hsg_163 = buffer.data(hsg + 163);
    const auto *hsg_164 = buffer.data(hsg + 164);
    const auto *hsg_165 = buffer.data(hsg + 165);
    const auto *hsg_167 = buffer.data(hsg + 167);
    const auto *hsg_168 = buffer.data(hsg + 168);
    const auto *hsg_170 = buffer.data(hsg + 170);
    const auto *hsg_175 = buffer.data(hsg + 175);
    const auto *hsg_176 = buffer.data(hsg + 176);
    const auto *hsg_177 = buffer.data(hsg + 177);
    const auto *hsg_178 = buffer.data(hsg + 178);
    const auto *hsg_179 = buffer.data(hsg + 179);
    const auto *hsg_180 = buffer.data(hsg + 180);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, gsg_96, hsf0_60, hsf0_66, \
                         hsf1_60, hsf1_66, hsg_91, hsg_92, hsg_93, \
                         hsg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * hsg_91[k];

        t_131[k] = f_4 * hsf0_60[k]
                   - f_5 * hsf1_60[k]
                   + f_3 * pc_z[k] * hsg_92[k];

        t_132[k] = f_10 * gsg_96[k]
                   + f_4 * hsf0_66[k]
                   - f_5 * hsf1_66[k]
                   + f_3 * pc_x[k] * hsg_96[k];

        t_133[k] = f_3 * pc_z[k] * hsg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, gsg_50, gsg_100, \
                         hsf0_62, hsf1_62, hsg_95, hsg_96, hsg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * gsg_50[k]
                   + f_3 * pc_y[k] * hsg_95[k];

        t_135[k] = f_6 * hsf0_62[k]
                   - f_7 * hsf1_62[k]
                   + f_3 * pc_z[k] * hsg_95[k];

        t_136[k] = f_10 * gsg_100[k]
                   + f_3 * pc_x[k] * hsg_100[k];

        t_137[k] = f_3 * pc_z[k] * hsg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, gsg_55, gsg_102, gsg_103, \
                         gsg_104, hsf0_66, hsf1_66, hsg_100, hsg_102, hsg_103, \
                         hsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_10 * gsg_102[k]
                   + f_3 * pc_x[k] * hsg_102[k];

        t_139[k] = f_10 * gsg_103[k]
                   + f_3 * pc_x[k] * hsg_103[k];

        t_140[k] = f_10 * gsg_104[k]
                   + f_3 * pc_x[k] * hsg_104[k];

        t_141[k] = f_11 * gsg_55[k]
                   + f_1 * hsf0_66[k]
                   - f_2 * hsf1_66[k]
                   + f_3 * pc_y[k] * hsg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, gsg_59, hsf0_66, hsf0_67, \
                         hsf1_66, hsf1_67, hsg_100, hsg_101, hsg_102, \
                         hsg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * hsg_100[k];

        t_143[k] = f_4 * hsf0_66[k]
                   - f_5 * hsf1_66[k]
                   + f_3 * pc_z[k] * hsg_101[k];

        t_144[k] = f_6 * hsf0_67[k]
                   - f_7 * hsf1_67[k]
                   + f_3 * pc_z[k] * hsg_102[k];

        t_145[k] = f_11 * gsg_59[k]
                   + f_3 * pc_y[k] * hsg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, gsh0_63, gsg_45, \
                         gsg_60, gsh1_63, hsf0_69, hsf1_69, hsg_104, \
                         hsg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * hsf0_69[k]
                   - f_2 * hsf1_69[k]
                   + f_3 * pc_z[k] * hsg_104[k];

        t_147[k] = pa_z[k] * gsh0_63[k]
                   - f_8 * pc_z[k] * gsh1_63[k];

        t_148[k] = f_10 * gsg_60[k]
                   + f_3 * pc_y[k] * hsg_105[k];

        t_149[k] = f_9 * gsg_45[k]
                   + f_3 * pc_z[k] * hsg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, gsh0_66, gsg_62, \
                         gsg_110, gsh1_66, hsf0_75, hsf1_75, hsg_107, \
                         hsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * gsh0_66[k]
                   - f_8 * pc_z[k] * gsh1_66[k];

        t_151[k] = f_10 * gsg_62[k]
                   + f_3 * pc_y[k] * hsg_107[k];

        t_152[k] = f_10 * gsg_110[k]
                   + f_6 * hsf0_75[k]
                   - f_7 * hsf1_75[k]
                   + f_3 * pc_x[k] * hsg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, gsh0_69, gsg_48, gsg_65, \
                         gsh1_69, hsg_108, hsg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * gsh0_69[k]
                   - f_8 * pc_z[k] * gsh1_69[k];

        t_154[k] = f_9 * gsg_48[k]
                   + f_3 * pc_z[k] * hsg_108[k];

        t_155[k] = f_10 * gsg_65[k]
                   + f_3 * pc_y[k] * hsg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, gsg_114, gsg_115, gsg_116, gsg_117, \
                         hsf0_79, hsf1_79, hsg_114, hsg_115, hsg_116, \
                         hsg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_10 * gsg_114[k]
                   + f_4 * hsf0_79[k]
                   - f_5 * hsf1_79[k]
                   + f_3 * pc_x[k] * hsg_114[k];

        t_157[k] = f_10 * gsg_115[k]
                   + f_3 * pc_x[k] * hsg_115[k];

        t_158[k] = f_10 * gsg_116[k]
                   + f_3 * pc_x[k] * hsg_116[k];

        t_159[k] = f_10 * gsg_117[k]
                   + f_3 * pc_x[k] * hsg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, gsh0_78, gsg_55, \
                         gsg_118, gsg_119, gsh1_78, hsg_115, hsg_118, \
                         hsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_10 * gsg_118[k]
                   + f_3 * pc_x[k] * hsg_118[k];

        t_161[k] = f_10 * gsg_119[k]
                   + f_3 * pc_x[k] * hsg_119[k];

        t_162[k] = pa_z[k] * gsh0_78[k]
                   - f_8 * pc_z[k] * gsh1_78[k];

        t_163[k] = f_9 * gsg_55[k]
                   + f_3 * pc_z[k] * hsg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, gsg_72, gsg_73, gsg_74, hsf0_78, hsf0_79, \
                         hsf1_78, hsf1_79, hsg_117, hsg_118, hsg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * gsg_72[k]
                   + f_6 * hsf0_78[k]
                   - f_7 * hsf1_78[k]
                   + f_3 * pc_y[k] * hsg_117[k];

        t_165[k] = f_10 * gsg_73[k]
                   + f_4 * hsf0_79[k]
                   - f_5 * hsf1_79[k]
                   + f_3 * pc_y[k] * hsg_118[k];

        t_166[k] = f_10 * gsg_74[k]
                   + f_3 * pc_y[k] * hsg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, gsh0_105, gsg_59, \
                         gsg_60, gsg_75, gsh1_105, hsf0_79, hsf1_79, hsg_119, \
                         hsg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * gsg_59[k]
                   + f_1 * hsf0_79[k]
                   - f_2 * hsf1_79[k]
                   + f_3 * pc_z[k] * hsg_119[k];

        t_168[k] = pa_y[k] * gsh0_105[k]
                   - f_8 * pc_y[k] * gsh1_105[k];

        t_169[k] = f_9 * gsg_75[k]
                   + f_3 * pc_y[k] * hsg_120[k];

        t_170[k] = f_10 * gsg_60[k]
                   + f_3 * pc_z[k] * hsg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, gsh0_108, gsh0_110, gsh0_111, \
                         gsg_76, gsg_77, gsg_78, gsh1_108, gsh1_110, gsh1_111, \
                         hsg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * gsh0_108[k]
                   + f_10 * gsg_76[k]
                   - f_8 * pc_y[k] * gsh1_108[k];

        t_172[k] = f_9 * gsg_77[k]
                   + f_3 * pc_y[k] * hsg_122[k];

        t_173[k] = pa_y[k] * gsh0_110[k]
                   - f_8 * pc_y[k] * gsh1_110[k];

        t_174[k] = pa_y[k] * gsh0_111[k]
                   + f_11 * gsg_78[k]
                   - f_8 * pc_y[k] * gsh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, gsh0_114, gsg_63, \
                         gsg_80, gsg_130, gsh1_114, hsg_123, hsg_125, \
                         hsg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * gsg_63[k]
                   + f_3 * pc_z[k] * hsg_123[k];

        t_176[k] = f_9 * gsg_80[k]
                   + f_3 * pc_y[k] * hsg_125[k];

        t_177[k] = pa_y[k] * gsh0_114[k]
                   - f_8 * pc_y[k] * gsh1_114[k];

        t_178[k] = f_10 * gsg_130[k]
                   + f_3 * pc_x[k] * hsg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, gsg_131, gsg_132, gsg_133, gsg_134, \
                         hsg_131, hsg_132, hsg_133, hsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_10 * gsg_131[k]
                   + f_3 * pc_x[k] * hsg_131[k];

        t_180[k] = f_10 * gsg_132[k]
                   + f_3 * pc_x[k] * hsg_132[k];

        t_181[k] = f_10 * gsg_133[k]
                   + f_3 * pc_x[k] * hsg_133[k];

        t_182[k] = f_10 * gsg_134[k]
                   + f_3 * pc_x[k] * hsg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, gsg_70, gsg_85, gsg_87, hsf0_86, \
                         hsf0_88, hsf1_86, hsf1_88, hsg_130, hsg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * gsg_85[k]
                   + f_1 * hsf0_86[k]
                   - f_2 * hsf1_86[k]
                   + f_3 * pc_y[k] * hsg_130[k];

        t_184[k] = f_10 * gsg_70[k]
                   + f_3 * pc_z[k] * hsg_130[k];

        t_185[k] = f_9 * gsg_87[k]
                   + f_6 * hsf0_88[k]
                   - f_7 * hsf1_88[k]
                   + f_3 * pc_y[k] * hsg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, gsh0_125, gsg_88, gsg_89, gsh1_125, \
                         hsf0_89, hsf1_89, hsg_133, hsg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * gsg_88[k]
                   + f_4 * hsf0_89[k]
                   - f_5 * hsf1_89[k]
                   + f_3 * pc_y[k] * hsg_133[k];

        t_187[k] = f_9 * gsg_89[k]
                   + f_3 * pc_y[k] * hsg_134[k];

        t_188[k] = pa_y[k] * gsh0_125[k]
                   - f_8 * pc_y[k] * gsh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, gsg_75, gsg_135, \
                         hsf0_90, hsf1_90, hsg_135, hsg_136, hsg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * gsg_135[k]
                   + f_1 * hsf0_90[k]
                   - f_2 * hsf1_90[k]
                   + f_3 * pc_x[k] * hsg_135[k];

        t_190[k] = f_3 * pc_y[k] * hsg_135[k];

        t_191[k] = f_11 * gsg_75[k]
                   + f_3 * pc_z[k] * hsg_135[k];

        t_192[k] = f_4 * hsf0_90[k]
                   - f_5 * hsf1_90[k]
                   + f_3 * pc_y[k] * hsg_136[k];

        t_193[k] = f_3 * pc_y[k] * hsg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, gsg_140, hsf0_91, hsf0_92, \
                         hsf0_95, hsf1_91, hsf1_92, hsf1_95, hsg_138, hsg_139, \
                         hsg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_10 * gsg_140[k]
                   + f_6 * hsf0_95[k]
                   - f_7 * hsf1_95[k]
                   + f_3 * pc_x[k] * hsg_140[k];

        t_195[k] = f_6 * hsf0_91[k]
                   - f_7 * hsf1_91[k]
                   + f_3 * pc_y[k] * hsg_138[k];

        t_196[k] = f_4 * hsf0_92[k]
                   - f_5 * hsf1_92[k]
                   + f_3 * pc_y[k] * hsg_139[k];

        t_197[k] = f_3 * pc_y[k] * hsg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, gsg_144, gsg_145, gsg_146, gsg_147, \
                         hsf0_99, hsf1_99, hsg_144, hsg_145, hsg_146, \
                         hsg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_10 * gsg_144[k]
                   + f_4 * hsf0_99[k]
                   - f_5 * hsf1_99[k]
                   + f_3 * pc_x[k] * hsg_144[k];

        t_199[k] = f_10 * gsg_145[k]
                   + f_3 * pc_x[k] * hsg_145[k];

        t_200[k] = f_10 * gsg_146[k]
                   + f_3 * pc_x[k] * hsg_146[k];

        t_201[k] = f_10 * gsg_147[k]
                   + f_3 * pc_x[k] * hsg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, gsg_149, hsf0_96, hsf0_97, \
                         hsf1_96, hsf1_97, hsg_144, hsg_145, hsg_146, \
                         hsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * hsg_144[k];

        t_203[k] = f_10 * gsg_149[k]
                   + f_3 * pc_x[k] * hsg_149[k];

        t_204[k] = f_1 * hsf0_96[k]
                   - f_2 * hsf1_96[k]
                   + f_3 * pc_y[k] * hsg_145[k];

        t_205[k] = f_13 * hsf0_97[k]
                   - f_14 * hsf1_97[k]
                   + f_3 * pc_y[k] * hsg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, gsg_89, hsf0_98, hsf0_99, \
                         hsf1_98, hsf1_99, hsg_147, hsg_148, hsg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * hsf0_98[k]
                   - f_7 * hsf1_98[k]
                   + f_3 * pc_y[k] * hsg_147[k];

        t_207[k] = f_4 * hsf0_99[k]
                   - f_5 * hsf1_99[k]
                   + f_3 * pc_y[k] * hsg_148[k];

        t_208[k] = f_3 * pc_y[k] * hsg_149[k];

        t_209[k] = f_11 * gsg_89[k]
                   + f_1 * hsf0_99[k]
                   - f_2 * hsf1_99[k]
                   + f_3 * pc_z[k] * hsg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_x, pc_x, pc_y, pc_z, gsh0_210, \
                         gsh0_213, gsg_90, gsg_150, gsg_153, gsh1_210, gsh1_213, \
                         hsg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = pa_x[k] * gsh0_210[k]
                   + f_0 * gsg_150[k]
                   - f_8 * pc_x[k] * gsh1_210[k];

        t_211[k] = f_12 * gsg_90[k]
                   + f_3 * pc_y[k] * hsg_150[k];

        t_212[k] = f_3 * pc_z[k] * hsg_150[k];

        t_213[k] = pa_x[k] * gsh0_213[k]
                   + f_11 * gsg_153[k]
                   - f_8 * pc_x[k] * gsh1_213[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_x, pc_x, pc_z, gsh0_216, gsg_156, \
                         gsh1_216, hsf0_100, hsf1_100, hsg_151, hsg_152, \
                         hsg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * hsg_151[k];

        t_215[k] = f_4 * hsf0_100[k]
                   - f_5 * hsf1_100[k]
                   + f_3 * pc_z[k] * hsg_152[k];

        t_216[k] = pa_x[k] * gsh0_216[k]
                   + f_10 * gsg_156[k]
                   - f_8 * pc_x[k] * gsh1_216[k];

        t_217[k] = f_3 * pc_z[k] * hsg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, gsg_95, gsg_160, \
                         hsf0_102, hsf1_102, hsg_155, hsg_156, \
                         hsg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_12 * gsg_95[k]
                   + f_3 * pc_y[k] * hsg_155[k];

        t_219[k] = f_6 * hsf0_102[k]
                   - f_7 * hsf1_102[k]
                   + f_3 * pc_z[k] * hsg_155[k];

        t_220[k] = f_9 * gsg_160[k]
                   + f_3 * pc_x[k] * hsg_160[k];

        t_221[k] = f_3 * pc_z[k] * hsg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pa_x, pc_x, gsh0_225, gsg_162, gsg_163, \
                         gsg_164, gsh1_225, hsg_162, hsg_163, hsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_9 * gsg_162[k]
                   + f_3 * pc_x[k] * hsg_162[k];

        t_223[k] = f_9 * gsg_163[k]
                   + f_3 * pc_x[k] * hsg_163[k];

        t_224[k] = f_9 * gsg_164[k]
                   + f_3 * pc_x[k] * hsg_164[k];

        t_225[k] = pa_x[k] * gsh0_225[k]
                   - f_8 * pc_x[k] * gsh1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_x, pc_x, pc_y, pc_z, gsh0_227, \
                         gsh0_228, gsg_104, gsh1_227, gsh1_228, hsg_160, \
                         hsg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * hsg_160[k];

        t_227[k] = pa_x[k] * gsh0_227[k]
                   - f_8 * pc_x[k] * gsh1_227[k];

        t_228[k] = pa_x[k] * gsh0_228[k]
                   - f_8 * pc_x[k] * gsh1_228[k];

        t_229[k] = f_12 * gsg_104[k]
                   + f_3 * pc_y[k] * hsg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_x, pa_z, pc_x, pc_y, pc_z, gsh0_126, \
                         gsh0_230, gsg_90, gsg_105, gsh1_126, gsh1_230, \
                         hsg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pa_x[k] * gsh0_230[k]
                   - f_8 * pc_x[k] * gsh1_230[k];

        t_231[k] = pa_z[k] * gsh0_126[k]
                   - f_8 * pc_z[k] * gsh1_126[k];

        t_232[k] = f_11 * gsg_105[k]
                   + f_3 * pc_y[k] * hsg_165[k];

        t_233[k] = f_9 * gsg_90[k]
                   + f_3 * pc_z[k] * hsg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_x, pa_z, pc_x, pc_y, pc_z, gsh0_129, \
                         gsh0_236, gsg_107, gsg_170, gsh1_129, gsh1_236, \
                         hsg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * gsh0_129[k]
                   - f_8 * pc_z[k] * gsh1_129[k];

        t_235[k] = f_11 * gsg_107[k]
                   + f_3 * pc_y[k] * hsg_167[k];

        t_236[k] = pa_x[k] * gsh0_236[k]
                   + f_11 * gsg_170[k]
                   - f_8 * pc_x[k] * gsh1_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, gsh0_132, gsg_93, gsg_110, \
                         gsh1_132, hsg_168, hsg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * gsh0_132[k]
                   - f_8 * pc_z[k] * gsh1_132[k];

        t_238[k] = f_9 * gsg_93[k]
                   + f_3 * pc_z[k] * hsg_168[k];

        t_239[k] = f_11 * gsg_110[k]
                   + f_3 * pc_y[k] * hsg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pa_x, pc_x, gsh0_240, gsg_174, gsg_175, \
                         gsg_176, gsg_177, gsh1_240, hsg_175, hsg_176, \
                         hsg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pa_x[k] * gsh0_240[k]
                   + f_10 * gsg_174[k]
                   - f_8 * pc_x[k] * gsh1_240[k];

        t_241[k] = f_9 * gsg_175[k]
                   + f_3 * pc_x[k] * hsg_175[k];

        t_242[k] = f_9 * gsg_176[k]
                   + f_3 * pc_x[k] * hsg_176[k];

        t_243[k] = f_9 * gsg_177[k]
                   + f_3 * pc_x[k] * hsg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_x, pc_x, pc_z, gsh0_246, gsg_100, \
                         gsg_178, gsg_179, gsh1_246, hsg_175, hsg_178, \
                         hsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_9 * gsg_178[k]
                   + f_3 * pc_x[k] * hsg_178[k];

        t_245[k] = f_9 * gsg_179[k]
                   + f_3 * pc_x[k] * hsg_179[k];

        t_246[k] = pa_x[k] * gsh0_246[k]
                   - f_8 * pc_x[k] * gsh1_246[k];

        t_247[k] = f_9 * gsg_100[k]
                   + f_3 * pc_z[k] * hsg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_x, pc_x, pc_y, gsh0_248, gsh0_249, \
                         gsh0_251, gsg_119, gsh1_248, gsh1_249, gsh1_251, \
                         hsg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pa_x[k] * gsh0_248[k]
                   - f_8 * pc_x[k] * gsh1_248[k];

        t_249[k] = pa_x[k] * gsh0_249[k]
                   - f_8 * pc_x[k] * gsh1_249[k];

        t_250[k] = f_11 * gsg_119[k]
                   + f_3 * pc_y[k] * hsg_179[k];

        t_251[k] = pa_x[k] * gsh0_251[k]
                   - f_8 * pc_x[k] * gsh1_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pa_x, pc_x, pc_y, pc_z, gsh0_252, gsg_105, \
                         gsg_120, gsg_180, gsh1_252, hsg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_x[k] * gsh0_252[k]
                   + f_0 * gsg_180[k]
                   - f_8 * pc_x[k] * gsh1_252[k];

        t_253[k] = f_10 * gsg_120[k]
                   + f_3 * pc_y[k] * hsg_180[k];

        t_254[k] = f_10 * gsg_105[k]
                   + f_3 * pc_z[k] * hsg_180[k];
    }
}

static auto
compute_prim_hsh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t hsf0, const size_t hsf1,
                                                          const size_t hsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsh0_189 = buffer.data(gsh0 + 189);
    const auto *gsh0_194 = buffer.data(gsh0 + 194);
    const auto *gsh0_198 = buffer.data(gsh0 + 198);
    const auto *gsh0_210 = buffer.data(gsh0 + 210);
    const auto *gsh0_211 = buffer.data(gsh0 + 211);
    const auto *gsh0_213 = buffer.data(gsh0 + 213);
    const auto *gsh0_216 = buffer.data(gsh0 + 216);
    const auto *gsh0_225 = buffer.data(gsh0 + 225);
    const auto *gsh0_227 = buffer.data(gsh0 + 227);
    const auto *gsh0_228 = buffer.data(gsh0 + 228);
    const auto *gsh0_255 = buffer.data(gsh0 + 255);
    const auto *gsh0_257 = buffer.data(gsh0 + 257);
    const auto *gsh0_258 = buffer.data(gsh0 + 258);
    const auto *gsh0_261 = buffer.data(gsh0 + 261);
    const auto *gsh0_267 = buffer.data(gsh0 + 267);
    const auto *gsh0_269 = buffer.data(gsh0 + 269);
    const auto *gsh0_270 = buffer.data(gsh0 + 270);
    const auto *gsh0_272 = buffer.data(gsh0 + 272);
    const auto *gsh0_276 = buffer.data(gsh0 + 276);
    const auto *gsh0_279 = buffer.data(gsh0 + 279);
    const auto *gsh0_288 = buffer.data(gsh0 + 288);
    const auto *gsh0_290 = buffer.data(gsh0 + 290);
    const auto *gsh0_291 = buffer.data(gsh0 + 291);
    const auto *gsh0_293 = buffer.data(gsh0 + 293);
    const auto *gsh0_294 = buffer.data(gsh0 + 294);
    const auto *gsh0_299 = buffer.data(gsh0 + 299);
    const auto *gsh0_303 = buffer.data(gsh0 + 303);
    const auto *gsh0_309 = buffer.data(gsh0 + 309);
    const auto *gsh0_310 = buffer.data(gsh0 + 310);
    const auto *gsh0_311 = buffer.data(gsh0 + 311);
    const auto *gsh0_312 = buffer.data(gsh0 + 312);
    const auto *gsh0_314 = buffer.data(gsh0 + 314);

    const auto *gsg_108 = buffer.data(gsg + 108);
    const auto *gsg_115 = buffer.data(gsg + 115);
    const auto *gsg_120 = buffer.data(gsg + 120);
    const auto *gsg_122 = buffer.data(gsg + 122);
    const auto *gsg_123 = buffer.data(gsg + 123);
    const auto *gsg_125 = buffer.data(gsg + 125);
    const auto *gsg_130 = buffer.data(gsg + 130);
    const auto *gsg_134 = buffer.data(gsg + 134);
    const auto *gsg_135 = buffer.data(gsg + 135);
    const auto *gsg_137 = buffer.data(gsg + 137);
    const auto *gsg_140 = buffer.data(gsg + 140);
    const auto *gsg_149 = buffer.data(gsg + 149);
    const auto *gsg_160 = buffer.data(gsg + 160);
    const auto *gsg_161 = buffer.data(gsg + 161);
    const auto *gsg_162 = buffer.data(gsg + 162);
    const auto *gsg_164 = buffer.data(gsg + 164);
    const auto *gsg_175 = buffer.data(gsg + 175);
    const auto *gsg_179 = buffer.data(gsg + 179);
    const auto *gsg_183 = buffer.data(gsg + 183);
    const auto *gsg_185 = buffer.data(gsg + 185);
    const auto *gsg_186 = buffer.data(gsg + 186);
    const auto *gsg_189 = buffer.data(gsg + 189);
    const auto *gsg_190 = buffer.data(gsg + 190);
    const auto *gsg_191 = buffer.data(gsg + 191);
    const auto *gsg_192 = buffer.data(gsg + 192);
    const auto *gsg_193 = buffer.data(gsg + 193);
    const auto *gsg_194 = buffer.data(gsg + 194);
    const auto *gsg_198 = buffer.data(gsg + 198);
    const auto *gsg_201 = buffer.data(gsg + 201);
    const auto *gsg_205 = buffer.data(gsg + 205);
    const auto *gsg_206 = buffer.data(gsg + 206);
    const auto *gsg_207 = buffer.data(gsg + 207);
    const auto *gsg_208 = buffer.data(gsg + 208);
    const auto *gsg_209 = buffer.data(gsg + 209);
    const auto *gsg_210 = buffer.data(gsg + 210);
    const auto *gsg_215 = buffer.data(gsg + 215);
    const auto *gsg_219 = buffer.data(gsg + 219);
    const auto *gsg_220 = buffer.data(gsg + 220);
    const auto *gsg_221 = buffer.data(gsg + 221);
    const auto *gsg_222 = buffer.data(gsg + 222);
    const auto *gsg_224 = buffer.data(gsg + 224);

    const auto *gsh1_189 = buffer.data(gsh1 + 189);
    const auto *gsh1_194 = buffer.data(gsh1 + 194);
    const auto *gsh1_198 = buffer.data(gsh1 + 198);
    const auto *gsh1_210 = buffer.data(gsh1 + 210);
    const auto *gsh1_211 = buffer.data(gsh1 + 211);
    const auto *gsh1_213 = buffer.data(gsh1 + 213);
    const auto *gsh1_216 = buffer.data(gsh1 + 216);
    const auto *gsh1_225 = buffer.data(gsh1 + 225);
    const auto *gsh1_227 = buffer.data(gsh1 + 227);
    const auto *gsh1_228 = buffer.data(gsh1 + 228);
    const auto *gsh1_255 = buffer.data(gsh1 + 255);
    const auto *gsh1_257 = buffer.data(gsh1 + 257);
    const auto *gsh1_258 = buffer.data(gsh1 + 258);
    const auto *gsh1_261 = buffer.data(gsh1 + 261);
    const auto *gsh1_267 = buffer.data(gsh1 + 267);
    const auto *gsh1_269 = buffer.data(gsh1 + 269);
    const auto *gsh1_270 = buffer.data(gsh1 + 270);
    const auto *gsh1_272 = buffer.data(gsh1 + 272);
    const auto *gsh1_276 = buffer.data(gsh1 + 276);
    const auto *gsh1_279 = buffer.data(gsh1 + 279);
    const auto *gsh1_288 = buffer.data(gsh1 + 288);
    const auto *gsh1_290 = buffer.data(gsh1 + 290);
    const auto *gsh1_291 = buffer.data(gsh1 + 291);
    const auto *gsh1_293 = buffer.data(gsh1 + 293);
    const auto *gsh1_294 = buffer.data(gsh1 + 294);
    const auto *gsh1_299 = buffer.data(gsh1 + 299);
    const auto *gsh1_303 = buffer.data(gsh1 + 303);
    const auto *gsh1_309 = buffer.data(gsh1 + 309);
    const auto *gsh1_310 = buffer.data(gsh1 + 310);
    const auto *gsh1_311 = buffer.data(gsh1 + 311);
    const auto *gsh1_312 = buffer.data(gsh1 + 312);
    const auto *gsh1_314 = buffer.data(gsh1 + 314);

    const auto *hsf0_140 = buffer.data(hsf0 + 140);
    const auto *hsf0_141 = buffer.data(hsf0 + 141);
    const auto *hsf0_142 = buffer.data(hsf0 + 142);
    const auto *hsf0_150 = buffer.data(hsf0 + 150);
    const auto *hsf0_151 = buffer.data(hsf0 + 151);
    const auto *hsf0_153 = buffer.data(hsf0 + 153);
    const auto *hsf0_155 = buffer.data(hsf0 + 155);
    const auto *hsf0_156 = buffer.data(hsf0 + 156);
    const auto *hsf0_157 = buffer.data(hsf0 + 157);
    const auto *hsf0_158 = buffer.data(hsf0 + 158);
    const auto *hsf0_159 = buffer.data(hsf0 + 159);
    const auto *hsf0_162 = buffer.data(hsf0 + 162);
    const auto *hsf0_164 = buffer.data(hsf0 + 164);
    const auto *hsf0_165 = buffer.data(hsf0 + 165);
    const auto *hsf0_167 = buffer.data(hsf0 + 167);
    const auto *hsf0_168 = buffer.data(hsf0 + 168);
    const auto *hsf0_169 = buffer.data(hsf0 + 169);
    const auto *hsf0_170 = buffer.data(hsf0 + 170);
    const auto *hsf0_171 = buffer.data(hsf0 + 171);
    const auto *hsf0_172 = buffer.data(hsf0 + 172);
    const auto *hsf0_173 = buffer.data(hsf0 + 173);
    const auto *hsf0_174 = buffer.data(hsf0 + 174);
    const auto *hsf0_175 = buffer.data(hsf0 + 175);
    const auto *hsf0_176 = buffer.data(hsf0 + 176);
    const auto *hsf0_177 = buffer.data(hsf0 + 177);
    const auto *hsf0_178 = buffer.data(hsf0 + 178);
    const auto *hsf0_179 = buffer.data(hsf0 + 179);
    const auto *hsf0_180 = buffer.data(hsf0 + 180);
    const auto *hsf0_181 = buffer.data(hsf0 + 181);

    const auto *hsf1_140 = buffer.data(hsf1 + 140);
    const auto *hsf1_141 = buffer.data(hsf1 + 141);
    const auto *hsf1_142 = buffer.data(hsf1 + 142);
    const auto *hsf1_150 = buffer.data(hsf1 + 150);
    const auto *hsf1_151 = buffer.data(hsf1 + 151);
    const auto *hsf1_153 = buffer.data(hsf1 + 153);
    const auto *hsf1_155 = buffer.data(hsf1 + 155);
    const auto *hsf1_156 = buffer.data(hsf1 + 156);
    const auto *hsf1_157 = buffer.data(hsf1 + 157);
    const auto *hsf1_158 = buffer.data(hsf1 + 158);
    const auto *hsf1_159 = buffer.data(hsf1 + 159);
    const auto *hsf1_162 = buffer.data(hsf1 + 162);
    const auto *hsf1_164 = buffer.data(hsf1 + 164);
    const auto *hsf1_165 = buffer.data(hsf1 + 165);
    const auto *hsf1_167 = buffer.data(hsf1 + 167);
    const auto *hsf1_168 = buffer.data(hsf1 + 168);
    const auto *hsf1_169 = buffer.data(hsf1 + 169);
    const auto *hsf1_170 = buffer.data(hsf1 + 170);
    const auto *hsf1_171 = buffer.data(hsf1 + 171);
    const auto *hsf1_172 = buffer.data(hsf1 + 172);
    const auto *hsf1_173 = buffer.data(hsf1 + 173);
    const auto *hsf1_174 = buffer.data(hsf1 + 174);
    const auto *hsf1_175 = buffer.data(hsf1 + 175);
    const auto *hsf1_176 = buffer.data(hsf1 + 176);
    const auto *hsf1_177 = buffer.data(hsf1 + 177);
    const auto *hsf1_178 = buffer.data(hsf1 + 178);
    const auto *hsf1_179 = buffer.data(hsf1 + 179);
    const auto *hsf1_180 = buffer.data(hsf1 + 180);
    const auto *hsf1_181 = buffer.data(hsf1 + 181);

    const auto *hsg_182 = buffer.data(hsg + 182);
    const auto *hsg_183 = buffer.data(hsg + 183);
    const auto *hsg_185 = buffer.data(hsg + 185);
    const auto *hsg_190 = buffer.data(hsg + 190);
    const auto *hsg_191 = buffer.data(hsg + 191);
    const auto *hsg_192 = buffer.data(hsg + 192);
    const auto *hsg_193 = buffer.data(hsg + 193);
    const auto *hsg_194 = buffer.data(hsg + 194);
    const auto *hsg_195 = buffer.data(hsg + 195);
    const auto *hsg_197 = buffer.data(hsg + 197);
    const auto *hsg_198 = buffer.data(hsg + 198);
    const auto *hsg_200 = buffer.data(hsg + 200);
    const auto *hsg_205 = buffer.data(hsg + 205);
    const auto *hsg_206 = buffer.data(hsg + 206);
    const auto *hsg_207 = buffer.data(hsg + 207);
    const auto *hsg_208 = buffer.data(hsg + 208);
    const auto *hsg_209 = buffer.data(hsg + 209);
    const auto *hsg_210 = buffer.data(hsg + 210);
    const auto *hsg_211 = buffer.data(hsg + 211);
    const auto *hsg_212 = buffer.data(hsg + 212);
    const auto *hsg_213 = buffer.data(hsg + 213);
    const auto *hsg_214 = buffer.data(hsg + 214);
    const auto *hsg_215 = buffer.data(hsg + 215);
    const auto *hsg_219 = buffer.data(hsg + 219);
    const auto *hsg_220 = buffer.data(hsg + 220);
    const auto *hsg_221 = buffer.data(hsg + 221);
    const auto *hsg_222 = buffer.data(hsg + 222);
    const auto *hsg_224 = buffer.data(hsg + 224);
    const auto *hsg_225 = buffer.data(hsg + 225);
    const auto *hsg_226 = buffer.data(hsg + 226);
    const auto *hsg_228 = buffer.data(hsg + 228);
    const auto *hsg_230 = buffer.data(hsg + 230);
    const auto *hsg_231 = buffer.data(hsg + 231);
    const auto *hsg_233 = buffer.data(hsg + 233);
    const auto *hsg_234 = buffer.data(hsg + 234);
    const auto *hsg_235 = buffer.data(hsg + 235);
    const auto *hsg_236 = buffer.data(hsg + 236);
    const auto *hsg_237 = buffer.data(hsg + 237);
    const auto *hsg_238 = buffer.data(hsg + 238);
    const auto *hsg_239 = buffer.data(hsg + 239);
    const auto *hsg_242 = buffer.data(hsg + 242);
    const auto *hsg_244 = buffer.data(hsg + 244);
    const auto *hsg_245 = buffer.data(hsg + 245);
    const auto *hsg_247 = buffer.data(hsg + 247);
    const auto *hsg_248 = buffer.data(hsg + 248);
    const auto *hsg_249 = buffer.data(hsg + 249);
    const auto *hsg_250 = buffer.data(hsg + 250);
    const auto *hsg_251 = buffer.data(hsg + 251);
    const auto *hsg_252 = buffer.data(hsg + 252);
    const auto *hsg_253 = buffer.data(hsg + 253);
    const auto *hsg_254 = buffer.data(hsg + 254);
    const auto *hsg_255 = buffer.data(hsg + 255);
    const auto *hsg_256 = buffer.data(hsg + 256);
    const auto *hsg_257 = buffer.data(hsg + 257);
    const auto *hsg_258 = buffer.data(hsg + 258);
    const auto *hsg_259 = buffer.data(hsg + 259);
    const auto *hsg_260 = buffer.data(hsg + 260);
    const auto *hsg_261 = buffer.data(hsg + 261);
    const auto *hsg_262 = buffer.data(hsg + 262);
    const auto *hsg_263 = buffer.data(hsg + 263);
    const auto *hsg_264 = buffer.data(hsg + 264);
    const auto *hsg_265 = buffer.data(hsg + 265);
    const auto *hsg_266 = buffer.data(hsg + 266);
    const auto *hsg_267 = buffer.data(hsg + 267);
    const auto *hsg_268 = buffer.data(hsg + 268);
    const auto *hsg_269 = buffer.data(hsg + 269);
    const auto *hsg_270 = buffer.data(hsg + 270);
    const auto *hsg_271 = buffer.data(hsg + 271);

#pragma omp simd aligned(t_255, t_256, t_257, pa_x, pc_x, pc_y, gsh0_255, gsh0_257, gsg_122, \
                         gsg_183, gsg_185, gsh1_255, gsh1_257, \
                         hsg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pa_x[k] * gsh0_255[k]
                   + f_11 * gsg_183[k]
                   - f_8 * pc_x[k] * gsh1_255[k];

        t_256[k] = f_10 * gsg_122[k]
                   + f_3 * pc_y[k] * hsg_182[k];

        t_257[k] = pa_x[k] * gsh0_257[k]
                   + f_11 * gsg_185[k]
                   - f_8 * pc_x[k] * gsh1_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_x, pc_x, pc_y, pc_z, gsh0_258, gsg_108, \
                         gsg_125, gsg_186, gsh1_258, hsg_183, hsg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pa_x[k] * gsh0_258[k]
                   + f_10 * gsg_186[k]
                   - f_8 * pc_x[k] * gsh1_258[k];

        t_259[k] = f_10 * gsg_108[k]
                   + f_3 * pc_z[k] * hsg_183[k];

        t_260[k] = f_10 * gsg_125[k]
                   + f_3 * pc_y[k] * hsg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_x, pc_x, gsh0_261, gsg_189, gsg_190, \
                         gsg_191, gsg_192, gsh1_261, hsg_190, hsg_191, \
                         hsg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_x[k] * gsh0_261[k]
                   + f_10 * gsg_189[k]
                   - f_8 * pc_x[k] * gsh1_261[k];

        t_262[k] = f_9 * gsg_190[k]
                   + f_3 * pc_x[k] * hsg_190[k];

        t_263[k] = f_9 * gsg_191[k]
                   + f_3 * pc_x[k] * hsg_191[k];

        t_264[k] = f_9 * gsg_192[k]
                   + f_3 * pc_x[k] * hsg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pa_x, pc_x, pc_z, gsh0_267, gsg_115, \
                         gsg_193, gsg_194, gsh1_267, hsg_190, hsg_193, \
                         hsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_9 * gsg_193[k]
                   + f_3 * pc_x[k] * hsg_193[k];

        t_266[k] = f_9 * gsg_194[k]
                   + f_3 * pc_x[k] * hsg_194[k];

        t_267[k] = pa_x[k] * gsh0_267[k]
                   - f_8 * pc_x[k] * gsh1_267[k];

        t_268[k] = f_10 * gsg_115[k]
                   + f_3 * pc_z[k] * hsg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pa_x, pc_x, pc_y, gsh0_269, gsh0_270, \
                         gsh0_272, gsg_134, gsh1_269, gsh1_270, gsh1_272, \
                         hsg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pa_x[k] * gsh0_269[k]
                   - f_8 * pc_x[k] * gsh1_269[k];

        t_270[k] = pa_x[k] * gsh0_270[k]
                   - f_8 * pc_x[k] * gsh1_270[k];

        t_271[k] = f_10 * gsg_134[k]
                   + f_3 * pc_y[k] * hsg_194[k];

        t_272[k] = pa_x[k] * gsh0_272[k]
                   - f_8 * pc_x[k] * gsh1_272[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_y, pc_y, pc_z, gsh0_189, gsg_120, gsg_135, \
                         gsh1_189, hsg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pa_y[k] * gsh0_189[k]
                   - f_8 * pc_y[k] * gsh1_189[k];

        t_274[k] = f_9 * gsg_135[k]
                   + f_3 * pc_y[k] * hsg_195[k];

        t_275[k] = f_11 * gsg_120[k]
                   + f_3 * pc_z[k] * hsg_195[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pa_x, pa_y, pc_x, pc_y, gsh0_194, gsh0_276, \
                         gsg_137, gsg_198, gsh1_194, gsh1_276, \
                         hsg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_x[k] * gsh0_276[k]
                   + f_11 * gsg_198[k]
                   - f_8 * pc_x[k] * gsh1_276[k];

        t_277[k] = f_9 * gsg_137[k]
                   + f_3 * pc_y[k] * hsg_197[k];

        t_278[k] = pa_y[k] * gsh0_194[k]
                   - f_8 * pc_y[k] * gsh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_x, pc_x, pc_y, pc_z, gsh0_279, gsg_123, \
                         gsg_140, gsg_201, gsh1_279, hsg_198, hsg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_x[k] * gsh0_279[k]
                   + f_10 * gsg_201[k]
                   - f_8 * pc_x[k] * gsh1_279[k];

        t_280[k] = f_11 * gsg_123[k]
                   + f_3 * pc_z[k] * hsg_198[k];

        t_281[k] = f_9 * gsg_140[k]
                   + f_3 * pc_y[k] * hsg_200[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_y, pc_x, pc_y, gsh0_198, gsg_205, \
                         gsg_206, gsg_207, gsh1_198, hsg_205, hsg_206, \
                         hsg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pa_y[k] * gsh0_198[k]
                   - f_8 * pc_y[k] * gsh1_198[k];

        t_283[k] = f_9 * gsg_205[k]
                   + f_3 * pc_x[k] * hsg_205[k];

        t_284[k] = f_9 * gsg_206[k]
                   + f_3 * pc_x[k] * hsg_206[k];

        t_285[k] = f_9 * gsg_207[k]
                   + f_3 * pc_x[k] * hsg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_x, pc_x, pc_z, gsh0_288, gsg_130, \
                         gsg_208, gsg_209, gsh1_288, hsg_205, hsg_208, \
                         hsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_9 * gsg_208[k]
                   + f_3 * pc_x[k] * hsg_208[k];

        t_287[k] = f_9 * gsg_209[k]
                   + f_3 * pc_x[k] * hsg_209[k];

        t_288[k] = pa_x[k] * gsh0_288[k]
                   - f_8 * pc_x[k] * gsh1_288[k];

        t_289[k] = f_11 * gsg_130[k]
                   + f_3 * pc_z[k] * hsg_205[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_x, pc_x, pc_y, gsh0_290, gsh0_291, \
                         gsh0_293, gsg_149, gsh1_290, gsh1_291, gsh1_293, \
                         hsg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_x[k] * gsh0_290[k]
                   - f_8 * pc_x[k] * gsh1_290[k];

        t_291[k] = pa_x[k] * gsh0_291[k]
                   - f_8 * pc_x[k] * gsh1_291[k];

        t_292[k] = f_9 * gsg_149[k]
                   + f_3 * pc_y[k] * hsg_209[k];

        t_293[k] = pa_x[k] * gsh0_293[k]
                   - f_8 * pc_x[k] * gsh1_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_x, pc_x, pc_y, pc_z, gsh0_294, \
                         gsg_135, gsg_210, gsh1_294, hsf0_140, hsf1_140, hsg_210, \
                         hsg_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pa_x[k] * gsh0_294[k]
                   + f_0 * gsg_210[k]
                   - f_8 * pc_x[k] * gsh1_294[k];

        t_295[k] = f_3 * pc_y[k] * hsg_210[k];

        t_296[k] = f_12 * gsg_135[k]
                   + f_3 * pc_z[k] * hsg_210[k];

        t_297[k] = f_4 * hsf0_140[k]
                   - f_5 * hsf1_140[k]
                   + f_3 * pc_y[k] * hsg_211[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pa_x, pc_x, pc_y, gsh0_299, gsg_215, gsh1_299, \
                         hsf0_141, hsf1_141, hsg_212, hsg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * hsg_212[k];

        t_299[k] = pa_x[k] * gsh0_299[k]
                   + f_11 * gsg_215[k]
                   - f_8 * pc_x[k] * gsh1_299[k];

        t_300[k] = f_6 * hsf0_141[k]
                   - f_7 * hsf1_141[k]
                   + f_3 * pc_y[k] * hsg_213[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_x, pc_x, pc_y, gsh0_303, gsg_219, \
                         gsg_220, gsh1_303, hsf0_142, hsf1_142, hsg_214, hsg_215, \
                         hsg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_4 * hsf0_142[k]
                   - f_5 * hsf1_142[k]
                   + f_3 * pc_y[k] * hsg_214[k];

        t_302[k] = f_3 * pc_y[k] * hsg_215[k];

        t_303[k] = pa_x[k] * gsh0_303[k]
                   + f_10 * gsg_219[k]
                   - f_8 * pc_x[k] * gsh1_303[k];

        t_304[k] = f_9 * gsg_220[k]
                   + f_3 * pc_x[k] * hsg_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, gsg_221, gsg_222, gsg_224, \
                         hsg_219, hsg_221, hsg_222, hsg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_9 * gsg_221[k]
                   + f_3 * pc_x[k] * hsg_221[k];

        t_306[k] = f_9 * gsg_222[k]
                   + f_3 * pc_x[k] * hsg_222[k];

        t_307[k] = f_3 * pc_y[k] * hsg_219[k];

        t_308[k] = f_9 * gsg_224[k]
                   + f_3 * pc_x[k] * hsg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_x, pc_x, gsh0_309, gsh0_310, gsh0_311, \
                         gsh0_312, gsh1_309, gsh1_310, gsh1_311, \
                         gsh1_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * gsh0_309[k]
                   - f_8 * pc_x[k] * gsh1_309[k];

        t_310[k] = pa_x[k] * gsh0_310[k]
                   - f_8 * pc_x[k] * gsh1_310[k];

        t_311[k] = pa_x[k] * gsh0_311[k]
                   - f_8 * pc_x[k] * gsh1_311[k];

        t_312[k] = pa_x[k] * gsh0_312[k]
                   - f_8 * pc_x[k] * gsh1_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, pa_x, pc_x, pc_y, gsh0_314, gsh1_314, \
                         hsf0_150, hsf0_151, hsf1_150, hsf1_151, hsg_224, hsg_225, \
                         hsg_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * hsg_224[k];

        t_314[k] = pa_x[k] * gsh0_314[k]
                   - f_8 * pc_x[k] * gsh1_314[k];

        t_315[k] = f_1 * hsf0_150[k]
                   - f_2 * hsf1_150[k]
                   + f_3 * pc_x[k] * hsg_225[k];

        t_316[k] = f_13 * hsf0_151[k]
                   - f_14 * hsf1_151[k]
                   + f_3 * pc_x[k] * hsg_226[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pc_x, pc_z, hsf0_153, hsf0_155, hsf1_153, \
                         hsf1_155, hsg_225, hsg_226, hsg_228, hsg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_3 * pc_z[k] * hsg_225[k];

        t_318[k] = f_6 * hsf0_153[k]
                   - f_7 * hsf1_153[k]
                   + f_3 * pc_x[k] * hsg_228[k];

        t_319[k] = f_3 * pc_z[k] * hsg_226[k];

        t_320[k] = f_6 * hsf0_155[k]
                   - f_7 * hsf1_155[k]
                   + f_3 * pc_x[k] * hsg_230[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pc_z, hsf0_156, hsf0_158, hsf0_159, \
                         hsf1_156, hsf1_158, hsf1_159, hsg_228, hsg_231, hsg_233, \
                         hsg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_4 * hsf0_156[k]
                   - f_5 * hsf1_156[k]
                   + f_3 * pc_x[k] * hsg_231[k];

        t_322[k] = f_3 * pc_z[k] * hsg_228[k];

        t_323[k] = f_4 * hsf0_158[k]
                   - f_5 * hsf1_158[k]
                   + f_3 * pc_x[k] * hsg_233[k];

        t_324[k] = f_4 * hsf0_159[k]
                   - f_5 * hsf1_159[k]
                   + f_3 * pc_x[k] * hsg_234[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, t_330, pc_x, pc_y, gsg_160, \
                         hsf0_156, hsf1_156, hsg_235, hsg_236, hsg_237, hsg_238, \
                         hsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_3 * pc_x[k] * hsg_235[k];

        t_326[k] = f_3 * pc_x[k] * hsg_236[k];

        t_327[k] = f_3 * pc_x[k] * hsg_237[k];

        t_328[k] = f_3 * pc_x[k] * hsg_238[k];

        t_329[k] = f_3 * pc_x[k] * hsg_239[k];

        t_330[k] = f_0 * gsg_160[k]
                   + f_1 * hsf0_156[k]
                   - f_2 * hsf1_156[k]
                   + f_3 * pc_y[k] * hsg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pc_y, pc_z, gsg_164, hsf0_156, hsf0_157, \
                         hsf1_156, hsf1_157, hsg_235, hsg_236, hsg_237, \
                         hsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * hsg_235[k];

        t_332[k] = f_4 * hsf0_156[k]
                   - f_5 * hsf1_156[k]
                   + f_3 * pc_z[k] * hsg_236[k];

        t_333[k] = f_6 * hsf0_157[k]
                   - f_7 * hsf1_157[k]
                   + f_3 * pc_z[k] * hsg_237[k];

        t_334[k] = f_0 * gsg_164[k]
                   + f_3 * pc_y[k] * hsg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pa_z, pc_z, gsh0_210, gsh0_211, gsh1_210, \
                         gsh1_211, hsf0_159, hsf1_159, hsg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_1 * hsf0_159[k]
                   - f_2 * hsf1_159[k]
                   + f_3 * pc_z[k] * hsg_239[k];

        t_336[k] = pa_z[k] * gsh0_210[k]
                   - f_8 * pc_z[k] * gsh1_210[k];

        t_337[k] = pa_z[k] * gsh0_211[k]
                   - f_8 * pc_z[k] * gsh1_211[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pa_z, pc_x, pc_z, gsh0_213, gsh1_213, hsf0_162, \
                         hsf0_164, hsf1_162, hsf1_164, hsg_242, \
                         hsg_244 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_13 * hsf0_162[k]
                   - f_14 * hsf1_162[k]
                   + f_3 * pc_x[k] * hsg_242[k];

        t_339[k] = pa_z[k] * gsh0_213[k]
                   - f_8 * pc_z[k] * gsh1_213[k];

        t_340[k] = f_6 * hsf0_164[k]
                   - f_7 * hsf1_164[k]
                   + f_3 * pc_x[k] * hsg_244[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pa_z, pc_x, pc_z, gsh0_216, gsh1_216, hsf0_165, \
                         hsf0_167, hsf1_165, hsf1_167, hsg_245, \
                         hsg_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_6 * hsf0_165[k]
                   - f_7 * hsf1_165[k]
                   + f_3 * pc_x[k] * hsg_245[k];

        t_342[k] = pa_z[k] * gsh0_216[k]
                   - f_8 * pc_z[k] * gsh1_216[k];

        t_343[k] = f_4 * hsf0_167[k]
                   - f_5 * hsf1_167[k]
                   + f_3 * pc_x[k] * hsg_247[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, t_348, pc_x, hsf0_168, hsf0_169, \
                         hsf1_168, hsf1_169, hsg_248, hsg_249, hsg_250, hsg_251, \
                         hsg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_4 * hsf0_168[k]
                   - f_5 * hsf1_168[k]
                   + f_3 * pc_x[k] * hsg_248[k];

        t_345[k] = f_4 * hsf0_169[k]
                   - f_5 * hsf1_169[k]
                   + f_3 * pc_x[k] * hsg_249[k];

        t_346[k] = f_3 * pc_x[k] * hsg_250[k];

        t_347[k] = f_3 * pc_x[k] * hsg_251[k];

        t_348[k] = f_3 * pc_x[k] * hsg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_z, pc_x, pc_z, gsh0_225, gsg_160, \
                         gsh1_225, hsg_250, hsg_253, hsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_3 * pc_x[k] * hsg_253[k];

        t_350[k] = f_3 * pc_x[k] * hsg_254[k];

        t_351[k] = pa_z[k] * gsh0_225[k]
                   - f_8 * pc_z[k] * gsh1_225[k];

        t_352[k] = f_9 * gsg_160[k]
                   + f_3 * pc_z[k] * hsg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pa_z, pc_y, pc_z, gsh0_227, gsh0_228, gsg_161, \
                         gsg_162, gsg_179, gsh1_227, gsh1_228, \
                         hsg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pa_z[k] * gsh0_227[k]
                   + f_10 * gsg_161[k]
                   - f_8 * pc_z[k] * gsh1_227[k];

        t_354[k] = pa_z[k] * gsh0_228[k]
                   + f_11 * gsg_162[k]
                   - f_8 * pc_z[k] * gsh1_228[k];

        t_355[k] = f_12 * gsg_179[k]
                   + f_3 * pc_y[k] * hsg_254[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_z, gsg_164, hsf0_169, hsf0_170, \
                         hsf0_171, hsf1_169, hsf1_170, hsf1_171, hsg_254, hsg_255, \
                         hsg_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_9 * gsg_164[k]
                   + f_1 * hsf0_169[k]
                   - f_2 * hsf1_169[k]
                   + f_3 * pc_z[k] * hsg_254[k];

        t_357[k] = f_1 * hsf0_170[k]
                   - f_2 * hsf1_170[k]
                   + f_3 * pc_x[k] * hsg_255[k];

        t_358[k] = f_13 * hsf0_171[k]
                   - f_14 * hsf1_171[k]
                   + f_3 * pc_x[k] * hsg_256[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, hsf0_172, hsf0_173, hsf0_174, hsf1_172, \
                         hsf1_173, hsf1_174, hsg_257, hsg_258, \
                         hsg_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_13 * hsf0_172[k]
                   - f_14 * hsf1_172[k]
                   + f_3 * pc_x[k] * hsg_257[k];

        t_360[k] = f_6 * hsf0_173[k]
                   - f_7 * hsf1_173[k]
                   + f_3 * pc_x[k] * hsg_258[k];

        t_361[k] = f_6 * hsf0_174[k]
                   - f_7 * hsf1_174[k]
                   + f_3 * pc_x[k] * hsg_259[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pc_x, hsf0_175, hsf0_176, hsf0_177, hsf1_175, \
                         hsf1_176, hsf1_177, hsg_260, hsg_261, \
                         hsg_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_6 * hsf0_175[k]
                   - f_7 * hsf1_175[k]
                   + f_3 * pc_x[k] * hsg_260[k];

        t_363[k] = f_4 * hsf0_176[k]
                   - f_5 * hsf1_176[k]
                   + f_3 * pc_x[k] * hsg_261[k];

        t_364[k] = f_4 * hsf0_177[k]
                   - f_5 * hsf1_177[k]
                   + f_3 * pc_x[k] * hsg_262[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, pc_x, hsf0_178, hsf0_179, \
                         hsf1_178, hsf1_179, hsg_263, hsg_264, hsg_265, hsg_266, \
                         hsg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_4 * hsf0_178[k]
                   - f_5 * hsf1_178[k]
                   + f_3 * pc_x[k] * hsg_263[k];

        t_366[k] = f_4 * hsf0_179[k]
                   - f_5 * hsf1_179[k]
                   + f_3 * pc_x[k] * hsg_264[k];

        t_367[k] = f_3 * pc_x[k] * hsg_265[k];

        t_368[k] = f_3 * pc_x[k] * hsg_266[k];

        t_369[k] = f_3 * pc_x[k] * hsg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pc_x, pc_y, pc_z, gsg_175, gsg_190, \
                         hsf0_176, hsf1_176, hsg_265, hsg_268, \
                         hsg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_3 * pc_x[k] * hsg_268[k];

        t_371[k] = f_3 * pc_x[k] * hsg_269[k];

        t_372[k] = f_11 * gsg_190[k]
                   + f_1 * hsf0_176[k]
                   - f_2 * hsf1_176[k]
                   + f_3 * pc_y[k] * hsg_265[k];

        t_373[k] = f_10 * gsg_175[k]
                   + f_3 * pc_z[k] * hsg_265[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pc_y, gsg_192, gsg_193, gsg_194, hsf0_178, \
                         hsf0_179, hsf1_178, hsf1_179, hsg_267, hsg_268, \
                         hsg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_11 * gsg_192[k]
                   + f_6 * hsf0_178[k]
                   - f_7 * hsf1_178[k]
                   + f_3 * pc_y[k] * hsg_267[k];

        t_375[k] = f_11 * gsg_193[k]
                   + f_4 * hsf0_179[k]
                   - f_5 * hsf1_179[k]
                   + f_3 * pc_y[k] * hsg_268[k];

        t_376[k] = f_11 * gsg_194[k]
                   + f_3 * pc_y[k] * hsg_269[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pc_x, pc_z, gsg_179, hsf0_179, hsf0_180, \
                         hsf0_181, hsf1_179, hsf1_180, hsf1_181, hsg_269, hsg_270, \
                         hsg_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_10 * gsg_179[k]
                   + f_1 * hsf0_179[k]
                   - f_2 * hsf1_179[k]
                   + f_3 * pc_z[k] * hsg_269[k];

        t_378[k] = f_1 * hsf0_180[k]
                   - f_2 * hsf1_180[k]
                   + f_3 * pc_x[k] * hsg_270[k];

        t_379[k] = f_13 * hsf0_181[k]
                   - f_14 * hsf1_181[k]
                   + f_3 * pc_x[k] * hsg_271[k];
    }
}

static auto
compute_prim_hsh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t hsf0, const size_t hsf1,
                                                          const size_t hsg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsh0_294 = buffer.data(gsh0 + 294);
    const auto *gsh0_296 = buffer.data(gsh0 + 296);
    const auto *gsh0_299 = buffer.data(gsh0 + 299);
    const auto *gsh0_303 = buffer.data(gsh0 + 303);
    const auto *gsh0_309 = buffer.data(gsh0 + 309);
    const auto *gsh0_311 = buffer.data(gsh0 + 311);
    const auto *gsh0_312 = buffer.data(gsh0 + 312);
    const auto *gsh0_314 = buffer.data(gsh0 + 314);

    const auto *gsg_190 = buffer.data(gsg + 190);
    const auto *gsg_194 = buffer.data(gsg + 194);
    const auto *gsg_205 = buffer.data(gsg + 205);
    const auto *gsg_207 = buffer.data(gsg + 207);
    const auto *gsg_208 = buffer.data(gsg + 208);
    const auto *gsg_209 = buffer.data(gsg + 209);
    const auto *gsg_220 = buffer.data(gsg + 220);
    const auto *gsg_222 = buffer.data(gsg + 222);
    const auto *gsg_223 = buffer.data(gsg + 223);
    const auto *gsg_224 = buffer.data(gsg + 224);

    const auto *gsh1_294 = buffer.data(gsh1 + 294);
    const auto *gsh1_296 = buffer.data(gsh1 + 296);
    const auto *gsh1_299 = buffer.data(gsh1 + 299);
    const auto *gsh1_303 = buffer.data(gsh1 + 303);
    const auto *gsh1_309 = buffer.data(gsh1 + 309);
    const auto *gsh1_311 = buffer.data(gsh1 + 311);
    const auto *gsh1_312 = buffer.data(gsh1 + 312);
    const auto *gsh1_314 = buffer.data(gsh1 + 314);

    const auto *hsf0_182 = buffer.data(hsf0 + 182);
    const auto *hsf0_183 = buffer.data(hsf0 + 183);
    const auto *hsf0_184 = buffer.data(hsf0 + 184);
    const auto *hsf0_185 = buffer.data(hsf0 + 185);
    const auto *hsf0_186 = buffer.data(hsf0 + 186);
    const auto *hsf0_187 = buffer.data(hsf0 + 187);
    const auto *hsf0_188 = buffer.data(hsf0 + 188);
    const auto *hsf0_189 = buffer.data(hsf0 + 189);
    const auto *hsf0_191 = buffer.data(hsf0 + 191);
    const auto *hsf0_193 = buffer.data(hsf0 + 193);
    const auto *hsf0_194 = buffer.data(hsf0 + 194);
    const auto *hsf0_196 = buffer.data(hsf0 + 196);
    const auto *hsf0_197 = buffer.data(hsf0 + 197);
    const auto *hsf0_198 = buffer.data(hsf0 + 198);
    const auto *hsf0_200 = buffer.data(hsf0 + 200);
    const auto *hsf0_202 = buffer.data(hsf0 + 202);
    const auto *hsf0_203 = buffer.data(hsf0 + 203);
    const auto *hsf0_205 = buffer.data(hsf0 + 205);
    const auto *hsf0_206 = buffer.data(hsf0 + 206);
    const auto *hsf0_207 = buffer.data(hsf0 + 207);
    const auto *hsf0_208 = buffer.data(hsf0 + 208);
    const auto *hsf0_209 = buffer.data(hsf0 + 209);

    const auto *hsf1_182 = buffer.data(hsf1 + 182);
    const auto *hsf1_183 = buffer.data(hsf1 + 183);
    const auto *hsf1_184 = buffer.data(hsf1 + 184);
    const auto *hsf1_185 = buffer.data(hsf1 + 185);
    const auto *hsf1_186 = buffer.data(hsf1 + 186);
    const auto *hsf1_187 = buffer.data(hsf1 + 187);
    const auto *hsf1_188 = buffer.data(hsf1 + 188);
    const auto *hsf1_189 = buffer.data(hsf1 + 189);
    const auto *hsf1_191 = buffer.data(hsf1 + 191);
    const auto *hsf1_193 = buffer.data(hsf1 + 193);
    const auto *hsf1_194 = buffer.data(hsf1 + 194);
    const auto *hsf1_196 = buffer.data(hsf1 + 196);
    const auto *hsf1_197 = buffer.data(hsf1 + 197);
    const auto *hsf1_198 = buffer.data(hsf1 + 198);
    const auto *hsf1_200 = buffer.data(hsf1 + 200);
    const auto *hsf1_202 = buffer.data(hsf1 + 202);
    const auto *hsf1_203 = buffer.data(hsf1 + 203);
    const auto *hsf1_205 = buffer.data(hsf1 + 205);
    const auto *hsf1_206 = buffer.data(hsf1 + 206);
    const auto *hsf1_207 = buffer.data(hsf1 + 207);
    const auto *hsf1_208 = buffer.data(hsf1 + 208);
    const auto *hsf1_209 = buffer.data(hsf1 + 209);

    const auto *hsg_272 = buffer.data(hsg + 272);
    const auto *hsg_273 = buffer.data(hsg + 273);
    const auto *hsg_274 = buffer.data(hsg + 274);
    const auto *hsg_275 = buffer.data(hsg + 275);
    const auto *hsg_276 = buffer.data(hsg + 276);
    const auto *hsg_277 = buffer.data(hsg + 277);
    const auto *hsg_278 = buffer.data(hsg + 278);
    const auto *hsg_279 = buffer.data(hsg + 279);
    const auto *hsg_280 = buffer.data(hsg + 280);
    const auto *hsg_281 = buffer.data(hsg + 281);
    const auto *hsg_282 = buffer.data(hsg + 282);
    const auto *hsg_283 = buffer.data(hsg + 283);
    const auto *hsg_284 = buffer.data(hsg + 284);
    const auto *hsg_286 = buffer.data(hsg + 286);
    const auto *hsg_288 = buffer.data(hsg + 288);
    const auto *hsg_289 = buffer.data(hsg + 289);
    const auto *hsg_291 = buffer.data(hsg + 291);
    const auto *hsg_292 = buffer.data(hsg + 292);
    const auto *hsg_293 = buffer.data(hsg + 293);
    const auto *hsg_295 = buffer.data(hsg + 295);
    const auto *hsg_296 = buffer.data(hsg + 296);
    const auto *hsg_297 = buffer.data(hsg + 297);
    const auto *hsg_298 = buffer.data(hsg + 298);
    const auto *hsg_299 = buffer.data(hsg + 299);
    const auto *hsg_300 = buffer.data(hsg + 300);
    const auto *hsg_302 = buffer.data(hsg + 302);
    const auto *hsg_303 = buffer.data(hsg + 303);
    const auto *hsg_305 = buffer.data(hsg + 305);
    const auto *hsg_306 = buffer.data(hsg + 306);
    const auto *hsg_307 = buffer.data(hsg + 307);
    const auto *hsg_309 = buffer.data(hsg + 309);
    const auto *hsg_310 = buffer.data(hsg + 310);
    const auto *hsg_311 = buffer.data(hsg + 311);
    const auto *hsg_312 = buffer.data(hsg + 312);
    const auto *hsg_313 = buffer.data(hsg + 313);
    const auto *hsg_314 = buffer.data(hsg + 314);

#pragma omp simd aligned(t_380, t_381, t_382, pc_x, hsf0_182, hsf0_183, hsf0_184, hsf1_182, \
                         hsf1_183, hsf1_184, hsg_272, hsg_273, \
                         hsg_274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_13 * hsf0_182[k]
                   - f_14 * hsf1_182[k]
                   + f_3 * pc_x[k] * hsg_272[k];

        t_381[k] = f_6 * hsf0_183[k]
                   - f_7 * hsf1_183[k]
                   + f_3 * pc_x[k] * hsg_273[k];

        t_382[k] = f_6 * hsf0_184[k]
                   - f_7 * hsf1_184[k]
                   + f_3 * pc_x[k] * hsg_274[k];
    }

#pragma omp simd aligned(t_383, t_384, t_385, pc_x, hsf0_185, hsf0_186, hsf0_187, hsf1_185, \
                         hsf1_186, hsf1_187, hsg_275, hsg_276, \
                         hsg_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_383[k] = f_6 * hsf0_185[k]
                   - f_7 * hsf1_185[k]
                   + f_3 * pc_x[k] * hsg_275[k];

        t_384[k] = f_4 * hsf0_186[k]
                   - f_5 * hsf1_186[k]
                   + f_3 * pc_x[k] * hsg_276[k];

        t_385[k] = f_4 * hsf0_187[k]
                   - f_5 * hsf1_187[k]
                   + f_3 * pc_x[k] * hsg_277[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, t_390, pc_x, hsf0_188, hsf0_189, \
                         hsf1_188, hsf1_189, hsg_278, hsg_279, hsg_280, hsg_281, \
                         hsg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_4 * hsf0_188[k]
                   - f_5 * hsf1_188[k]
                   + f_3 * pc_x[k] * hsg_278[k];

        t_387[k] = f_4 * hsf0_189[k]
                   - f_5 * hsf1_189[k]
                   + f_3 * pc_x[k] * hsg_279[k];

        t_388[k] = f_3 * pc_x[k] * hsg_280[k];

        t_389[k] = f_3 * pc_x[k] * hsg_281[k];

        t_390[k] = f_3 * pc_x[k] * hsg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, gsg_190, gsg_205, \
                         hsf0_186, hsf1_186, hsg_280, hsg_283, \
                         hsg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_3 * pc_x[k] * hsg_283[k];

        t_392[k] = f_3 * pc_x[k] * hsg_284[k];

        t_393[k] = f_10 * gsg_205[k]
                   + f_1 * hsf0_186[k]
                   - f_2 * hsf1_186[k]
                   + f_3 * pc_y[k] * hsg_280[k];

        t_394[k] = f_11 * gsg_190[k]
                   + f_3 * pc_z[k] * hsg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_y, gsg_207, gsg_208, gsg_209, hsf0_188, \
                         hsf0_189, hsf1_188, hsf1_189, hsg_282, hsg_283, \
                         hsg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_10 * gsg_207[k]
                   + f_6 * hsf0_188[k]
                   - f_7 * hsf1_188[k]
                   + f_3 * pc_y[k] * hsg_282[k];

        t_396[k] = f_10 * gsg_208[k]
                   + f_4 * hsf0_189[k]
                   - f_5 * hsf1_189[k]
                   + f_3 * pc_y[k] * hsg_283[k];

        t_397[k] = f_10 * gsg_209[k]
                   + f_3 * pc_y[k] * hsg_284[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_y, pc_x, pc_y, pc_z, gsh0_294, gsg_194, \
                         gsh1_294, hsf0_189, hsf0_191, hsf1_189, hsf1_191, hsg_284, \
                         hsg_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_11 * gsg_194[k]
                   + f_1 * hsf0_189[k]
                   - f_2 * hsf1_189[k]
                   + f_3 * pc_z[k] * hsg_284[k];

        t_399[k] = pa_y[k] * gsh0_294[k]
                   - f_8 * pc_y[k] * gsh1_294[k];

        t_400[k] = f_13 * hsf0_191[k]
                   - f_14 * hsf1_191[k]
                   + f_3 * pc_x[k] * hsg_286[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_y, pc_x, pc_y, gsh0_296, gsh1_296, hsf0_193, \
                         hsf0_194, hsf1_193, hsf1_194, hsg_288, \
                         hsg_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = pa_y[k] * gsh0_296[k]
                   - f_8 * pc_y[k] * gsh1_296[k];

        t_402[k] = f_6 * hsf0_193[k]
                   - f_7 * hsf1_193[k]
                   + f_3 * pc_x[k] * hsg_288[k];

        t_403[k] = f_6 * hsf0_194[k]
                   - f_7 * hsf1_194[k]
                   + f_3 * pc_x[k] * hsg_289[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_y, pc_x, pc_y, gsh0_299, gsh1_299, hsf0_196, \
                         hsf0_197, hsf1_196, hsf1_197, hsg_291, \
                         hsg_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * gsh0_299[k]
                   - f_8 * pc_y[k] * gsh1_299[k];

        t_405[k] = f_4 * hsf0_196[k]
                   - f_5 * hsf1_196[k]
                   + f_3 * pc_x[k] * hsg_291[k];

        t_406[k] = f_4 * hsf0_197[k]
                   - f_5 * hsf1_197[k]
                   + f_3 * pc_x[k] * hsg_292[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pa_y, pc_x, pc_y, gsh0_303, \
                         gsh1_303, hsf0_198, hsf1_198, hsg_293, hsg_295, hsg_296, \
                         hsg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_4 * hsf0_198[k]
                   - f_5 * hsf1_198[k]
                   + f_3 * pc_x[k] * hsg_293[k];

        t_408[k] = pa_y[k] * gsh0_303[k]
                   - f_8 * pc_y[k] * gsh1_303[k];

        t_409[k] = f_3 * pc_x[k] * hsg_295[k];

        t_410[k] = f_3 * pc_x[k] * hsg_296[k];

        t_411[k] = f_3 * pc_x[k] * hsg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_y, pc_x, pc_y, pc_z, gsh0_309, \
                         gsg_205, gsg_220, gsh1_309, hsg_295, hsg_298, \
                         hsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_3 * pc_x[k] * hsg_298[k];

        t_413[k] = f_3 * pc_x[k] * hsg_299[k];

        t_414[k] = pa_y[k] * gsh0_309[k]
                   + f_0 * gsg_220[k]
                   - f_8 * pc_y[k] * gsh1_309[k];

        t_415[k] = f_12 * gsg_205[k]
                   + f_3 * pc_z[k] * hsg_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_y, pc_y, gsh0_311, gsh0_312, gsh0_314, \
                         gsg_222, gsg_223, gsg_224, gsh1_311, gsh1_312, gsh1_314, \
                         hsg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pa_y[k] * gsh0_311[k]
                   + f_11 * gsg_222[k]
                   - f_8 * pc_y[k] * gsh1_311[k];

        t_417[k] = pa_y[k] * gsh0_312[k]
                   + f_10 * gsg_223[k]
                   - f_8 * pc_y[k] * gsh1_312[k];

        t_418[k] = f_9 * gsg_224[k]
                   + f_3 * pc_y[k] * hsg_299[k];

        t_419[k] = pa_y[k] * gsh0_314[k]
                   - f_8 * pc_y[k] * gsh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, hsf0_200, hsf0_202, \
                         hsf0_203, hsf1_200, hsf1_202, hsf1_203, hsg_300, hsg_302, \
                         hsg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_1 * hsf0_200[k]
                   - f_2 * hsf1_200[k]
                   + f_3 * pc_x[k] * hsg_300[k];

        t_421[k] = f_3 * pc_y[k] * hsg_300[k];

        t_422[k] = f_13 * hsf0_202[k]
                   - f_14 * hsf1_202[k]
                   + f_3 * pc_x[k] * hsg_302[k];

        t_423[k] = f_6 * hsf0_203[k]
                   - f_7 * hsf1_203[k]
                   + f_3 * pc_x[k] * hsg_303[k];

        t_424[k] = f_3 * pc_y[k] * hsg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, hsf0_205, hsf0_206, hsf0_207, \
                         hsf1_205, hsf1_206, hsf1_207, hsg_305, hsg_306, \
                         hsg_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_6 * hsf0_205[k]
                   - f_7 * hsf1_205[k]
                   + f_3 * pc_x[k] * hsg_305[k];

        t_426[k] = f_4 * hsf0_206[k]
                   - f_5 * hsf1_206[k]
                   + f_3 * pc_x[k] * hsg_306[k];

        t_427[k] = f_4 * hsf0_207[k]
                   - f_5 * hsf1_207[k]
                   + f_3 * pc_x[k] * hsg_307[k];

        t_428[k] = f_3 * pc_y[k] * hsg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, pc_x, hsf0_209, hsf1_209, \
                         hsg_309, hsg_310, hsg_311, hsg_312, hsg_313, \
                         hsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_4 * hsf0_209[k]
                   - f_5 * hsf1_209[k]
                   + f_3 * pc_x[k] * hsg_309[k];

        t_430[k] = f_3 * pc_x[k] * hsg_310[k];

        t_431[k] = f_3 * pc_x[k] * hsg_311[k];

        t_432[k] = f_3 * pc_x[k] * hsg_312[k];

        t_433[k] = f_3 * pc_x[k] * hsg_313[k];

        t_434[k] = f_3 * pc_x[k] * hsg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pc_y, hsf0_206, hsf0_207, hsf0_208, hsf1_206, \
                         hsf1_207, hsf1_208, hsg_310, hsg_311, \
                         hsg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * hsf0_206[k]
                   - f_2 * hsf1_206[k]
                   + f_3 * pc_y[k] * hsg_310[k];

        t_436[k] = f_13 * hsf0_207[k]
                   - f_14 * hsf1_207[k]
                   + f_3 * pc_y[k] * hsg_311[k];

        t_437[k] = f_6 * hsf0_208[k]
                   - f_7 * hsf1_208[k]
                   + f_3 * pc_y[k] * hsg_312[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pc_y, pc_z, gsg_224, hsf0_209, hsf1_209, \
                         hsg_313, hsg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_4 * hsf0_209[k]
                   - f_5 * hsf1_209[k]
                   + f_3 * pc_y[k] * hsg_313[k];

        t_439[k] = f_3 * pc_y[k] * hsg_314[k];

        t_440[k] = f_0 * gsg_224[k]
                   + f_1 * hsf0_209[k]
                   - f_2 * hsf1_209[k]
                   + f_3 * pc_z[k] * hsg_314[k];
    }
}

auto
compute_prim_hsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsh0, const size_t gsg,
                                                   const size_t gsh1, const size_t hsf0,
                                                   const size_t hsf1, const size_t hsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsh_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsh0, gsg,
                                                              gsh1, hsf0, hsf1, hsg, ncols,
                                                              gamma, p, q);

    compute_prim_hsh_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsh0, gsg,
                                                              gsh1, hsf0, hsf1, hsg, ncols,
                                                              gamma, p, q);

    compute_prim_hsh_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, gsh0, gsg,
                                                              gsh1, hsf0, hsf1, hsg, ncols,
                                                              gamma, p, q);

    compute_prim_hsh_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, gsh0, gsg,
                                                              gsh1, hsf0, hsf1, hsg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
