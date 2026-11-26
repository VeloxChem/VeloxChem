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


#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ish_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsh0,
                                                          const size_t hsg, const size_t hsh1,
                                                          const size_t isf0, const size_t isf1,
                                                          const size_t isg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.0 / q;

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

    const auto *hsh0_0 = buffer.data(hsh0 + 0);
    const auto *hsh0_3 = buffer.data(hsh0 + 3);
    const auto *hsh0_5 = buffer.data(hsh0 + 5);
    const auto *hsh0_6 = buffer.data(hsh0 + 6);
    const auto *hsh0_9 = buffer.data(hsh0 + 9);
    const auto *hsh0_15 = buffer.data(hsh0 + 15);
    const auto *hsh0_20 = buffer.data(hsh0 + 20);
    const auto *hsh0_24 = buffer.data(hsh0 + 24);
    const auto *hsh0_27 = buffer.data(hsh0 + 27);
    const auto *hsh0_36 = buffer.data(hsh0 + 36);
    const auto *hsh0_42 = buffer.data(hsh0 + 42);
    const auto *hsh0_47 = buffer.data(hsh0 + 47);
    const auto *hsh0_51 = buffer.data(hsh0 + 51);
    const auto *hsh0_62 = buffer.data(hsh0 + 62);

    const auto *hsg_0 = buffer.data(hsg + 0);
    const auto *hsg_1 = buffer.data(hsg + 1);
    const auto *hsg_2 = buffer.data(hsg + 2);
    const auto *hsg_3 = buffer.data(hsg + 3);
    const auto *hsg_5 = buffer.data(hsg + 5);
    const auto *hsg_10 = buffer.data(hsg + 10);
    const auto *hsg_12 = buffer.data(hsg + 12);
    const auto *hsg_14 = buffer.data(hsg + 14);
    const auto *hsg_15 = buffer.data(hsg + 15);
    const auto *hsg_18 = buffer.data(hsg + 18);
    const auto *hsg_20 = buffer.data(hsg + 20);
    const auto *hsg_25 = buffer.data(hsg + 25);
    const auto *hsg_27 = buffer.data(hsg + 27);
    const auto *hsg_28 = buffer.data(hsg + 28);
    const auto *hsg_29 = buffer.data(hsg + 29);
    const auto *hsg_30 = buffer.data(hsg + 30);
    const auto *hsg_32 = buffer.data(hsg + 32);
    const auto *hsg_35 = buffer.data(hsg + 35);
    const auto *hsg_40 = buffer.data(hsg + 40);
    const auto *hsg_41 = buffer.data(hsg + 41);
    const auto *hsg_42 = buffer.data(hsg + 42);
    const auto *hsg_43 = buffer.data(hsg + 43);
    const auto *hsg_44 = buffer.data(hsg + 44);
    const auto *hsg_45 = buffer.data(hsg + 45);
    const auto *hsg_48 = buffer.data(hsg + 48);
    const auto *hsg_51 = buffer.data(hsg + 51);
    const auto *hsg_55 = buffer.data(hsg + 55);
    const auto *hsg_57 = buffer.data(hsg + 57);
    const auto *hsg_58 = buffer.data(hsg + 58);
    const auto *hsg_59 = buffer.data(hsg + 59);
    const auto *hsg_70 = buffer.data(hsg + 70);
    const auto *hsg_71 = buffer.data(hsg + 71);
    const auto *hsg_72 = buffer.data(hsg + 72);
    const auto *hsg_73 = buffer.data(hsg + 73);
    const auto *hsg_74 = buffer.data(hsg + 74);
    const auto *hsg_75 = buffer.data(hsg + 75);
    const auto *hsg_80 = buffer.data(hsg + 80);
    const auto *hsg_84 = buffer.data(hsg + 84);
    const auto *hsg_85 = buffer.data(hsg + 85);
    const auto *hsg_86 = buffer.data(hsg + 86);
    const auto *hsg_87 = buffer.data(hsg + 87);
    const auto *hsg_89 = buffer.data(hsg + 89);
    const auto *hsg_90 = buffer.data(hsg + 90);
    const auto *hsg_93 = buffer.data(hsg + 93);

    const auto *hsh1_0 = buffer.data(hsh1 + 0);
    const auto *hsh1_3 = buffer.data(hsh1 + 3);
    const auto *hsh1_5 = buffer.data(hsh1 + 5);
    const auto *hsh1_6 = buffer.data(hsh1 + 6);
    const auto *hsh1_9 = buffer.data(hsh1 + 9);
    const auto *hsh1_15 = buffer.data(hsh1 + 15);
    const auto *hsh1_20 = buffer.data(hsh1 + 20);
    const auto *hsh1_24 = buffer.data(hsh1 + 24);
    const auto *hsh1_27 = buffer.data(hsh1 + 27);
    const auto *hsh1_36 = buffer.data(hsh1 + 36);
    const auto *hsh1_42 = buffer.data(hsh1 + 42);
    const auto *hsh1_47 = buffer.data(hsh1 + 47);
    const auto *hsh1_51 = buffer.data(hsh1 + 51);
    const auto *hsh1_62 = buffer.data(hsh1 + 62);

    const auto *isf0_0 = buffer.data(isf0 + 0);
    const auto *isf0_1 = buffer.data(isf0 + 1);
    const auto *isf0_2 = buffer.data(isf0 + 2);
    const auto *isf0_6 = buffer.data(isf0 + 6);
    const auto *isf0_8 = buffer.data(isf0 + 8);
    const auto *isf0_9 = buffer.data(isf0 + 9);
    const auto *isf0_16 = buffer.data(isf0 + 16);
    const auto *isf0_17 = buffer.data(isf0 + 17);
    const auto *isf0_22 = buffer.data(isf0 + 22);
    const auto *isf0_27 = buffer.data(isf0 + 27);
    const auto *isf0_28 = buffer.data(isf0 + 28);
    const auto *isf0_29 = buffer.data(isf0 + 29);
    const auto *isf0_30 = buffer.data(isf0 + 30);
    const auto *isf0_32 = buffer.data(isf0 + 32);
    const auto *isf0_33 = buffer.data(isf0 + 33);
    const auto *isf0_36 = buffer.data(isf0 + 36);
    const auto *isf0_37 = buffer.data(isf0 + 37);
    const auto *isf0_39 = buffer.data(isf0 + 39);
    const auto *isf0_48 = buffer.data(isf0 + 48);
    const auto *isf0_49 = buffer.data(isf0 + 49);
    const auto *isf0_50 = buffer.data(isf0 + 50);
    const auto *isf0_51 = buffer.data(isf0 + 51);
    const auto *isf0_52 = buffer.data(isf0 + 52);
    const auto *isf0_55 = buffer.data(isf0 + 55);
    const auto *isf0_56 = buffer.data(isf0 + 56);
    const auto *isf0_57 = buffer.data(isf0 + 57);
    const auto *isf0_58 = buffer.data(isf0 + 58);
    const auto *isf0_59 = buffer.data(isf0 + 59);
    const auto *isf0_60 = buffer.data(isf0 + 60);
    const auto *isf0_63 = buffer.data(isf0 + 63);

    const auto *isf1_0 = buffer.data(isf1 + 0);
    const auto *isf1_1 = buffer.data(isf1 + 1);
    const auto *isf1_2 = buffer.data(isf1 + 2);
    const auto *isf1_6 = buffer.data(isf1 + 6);
    const auto *isf1_8 = buffer.data(isf1 + 8);
    const auto *isf1_9 = buffer.data(isf1 + 9);
    const auto *isf1_16 = buffer.data(isf1 + 16);
    const auto *isf1_17 = buffer.data(isf1 + 17);
    const auto *isf1_22 = buffer.data(isf1 + 22);
    const auto *isf1_27 = buffer.data(isf1 + 27);
    const auto *isf1_28 = buffer.data(isf1 + 28);
    const auto *isf1_29 = buffer.data(isf1 + 29);
    const auto *isf1_30 = buffer.data(isf1 + 30);
    const auto *isf1_32 = buffer.data(isf1 + 32);
    const auto *isf1_33 = buffer.data(isf1 + 33);
    const auto *isf1_36 = buffer.data(isf1 + 36);
    const auto *isf1_37 = buffer.data(isf1 + 37);
    const auto *isf1_39 = buffer.data(isf1 + 39);
    const auto *isf1_48 = buffer.data(isf1 + 48);
    const auto *isf1_49 = buffer.data(isf1 + 49);
    const auto *isf1_50 = buffer.data(isf1 + 50);
    const auto *isf1_51 = buffer.data(isf1 + 51);
    const auto *isf1_52 = buffer.data(isf1 + 52);
    const auto *isf1_55 = buffer.data(isf1 + 55);
    const auto *isf1_56 = buffer.data(isf1 + 56);
    const auto *isf1_57 = buffer.data(isf1 + 57);
    const auto *isf1_58 = buffer.data(isf1 + 58);
    const auto *isf1_59 = buffer.data(isf1 + 59);
    const auto *isf1_60 = buffer.data(isf1 + 60);
    const auto *isf1_63 = buffer.data(isf1 + 63);

    const auto *isg_0 = buffer.data(isg + 0);
    const auto *isg_1 = buffer.data(isg + 1);
    const auto *isg_2 = buffer.data(isg + 2);
    const auto *isg_3 = buffer.data(isg + 3);
    const auto *isg_5 = buffer.data(isg + 5);
    const auto *isg_6 = buffer.data(isg + 6);
    const auto *isg_9 = buffer.data(isg + 9);
    const auto *isg_10 = buffer.data(isg + 10);
    const auto *isg_12 = buffer.data(isg + 12);
    const auto *isg_13 = buffer.data(isg + 13);
    const auto *isg_14 = buffer.data(isg + 14);
    const auto *isg_15 = buffer.data(isg + 15);
    const auto *isg_16 = buffer.data(isg + 16);
    const auto *isg_18 = buffer.data(isg + 18);
    const auto *isg_20 = buffer.data(isg + 20);
    const auto *isg_21 = buffer.data(isg + 21);
    const auto *isg_25 = buffer.data(isg + 25);
    const auto *isg_26 = buffer.data(isg + 26);
    const auto *isg_27 = buffer.data(isg + 27);
    const auto *isg_28 = buffer.data(isg + 28);
    const auto *isg_29 = buffer.data(isg + 29);
    const auto *isg_30 = buffer.data(isg + 30);
    const auto *isg_32 = buffer.data(isg + 32);
    const auto *isg_34 = buffer.data(isg + 34);
    const auto *isg_35 = buffer.data(isg + 35);
    const auto *isg_39 = buffer.data(isg + 39);
    const auto *isg_40 = buffer.data(isg + 40);
    const auto *isg_41 = buffer.data(isg + 41);
    const auto *isg_42 = buffer.data(isg + 42);
    const auto *isg_43 = buffer.data(isg + 43);
    const auto *isg_44 = buffer.data(isg + 44);
    const auto *isg_45 = buffer.data(isg + 45);
    const auto *isg_46 = buffer.data(isg + 46);
    const auto *isg_47 = buffer.data(isg + 47);
    const auto *isg_48 = buffer.data(isg + 48);
    const auto *isg_50 = buffer.data(isg + 50);
    const auto *isg_51 = buffer.data(isg + 51);
    const auto *isg_55 = buffer.data(isg + 55);
    const auto *isg_56 = buffer.data(isg + 56);
    const auto *isg_57 = buffer.data(isg + 57);
    const auto *isg_58 = buffer.data(isg + 58);
    const auto *isg_59 = buffer.data(isg + 59);
    const auto *isg_60 = buffer.data(isg + 60);
    const auto *isg_62 = buffer.data(isg + 62);
    const auto *isg_63 = buffer.data(isg + 63);
    const auto *isg_65 = buffer.data(isg + 65);
    const auto *isg_70 = buffer.data(isg + 70);
    const auto *isg_71 = buffer.data(isg + 71);
    const auto *isg_72 = buffer.data(isg + 72);
    const auto *isg_73 = buffer.data(isg + 73);
    const auto *isg_74 = buffer.data(isg + 74);
    const auto *isg_75 = buffer.data(isg + 75);
    const auto *isg_76 = buffer.data(isg + 76);
    const auto *isg_77 = buffer.data(isg + 77);
    const auto *isg_78 = buffer.data(isg + 78);
    const auto *isg_79 = buffer.data(isg + 79);
    const auto *isg_80 = buffer.data(isg + 80);
    const auto *isg_84 = buffer.data(isg + 84);
    const auto *isg_85 = buffer.data(isg + 85);
    const auto *isg_86 = buffer.data(isg + 86);
    const auto *isg_87 = buffer.data(isg + 87);
    const auto *isg_88 = buffer.data(isg + 88);
    const auto *isg_89 = buffer.data(isg + 89);
    const auto *isg_90 = buffer.data(isg + 90);
    const auto *isg_93 = buffer.data(isg + 93);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsg_0, isf0_0, \
                         isf1_0, isg_0, isg_1, isg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsg_0[k]
                 + f_1 * isf0_0[k]
                 - f_2 * isf1_0[k]
                 + f_3 * pc_x[k] * isg_0[k];

        t_1[k] = f_3 * pc_y[k] * isg_0[k];

        t_2[k] = f_3 * pc_z[k] * isg_0[k];

        t_3[k] = f_4 * isf0_0[k]
                 - f_5 * isf1_0[k]
                 + f_3 * pc_y[k] * isg_1[k];

        t_4[k] = f_3 * pc_y[k] * isg_2[k];

        t_5[k] = f_4 * isf0_0[k]
                 - f_5 * isf1_0[k]
                 + f_3 * pc_z[k] * isg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, hsg_10, isf0_1, isf0_2, \
                         isf1_1, isf1_2, isg_3, isg_5, isg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * isf0_1[k]
                 - f_7 * isf1_1[k]
                 + f_3 * pc_y[k] * isg_3[k];

        t_7[k] = f_3 * pc_z[k] * isg_3[k];

        t_8[k] = f_3 * pc_y[k] * isg_5[k];

        t_9[k] = f_6 * isf0_2[k]
                 - f_7 * isf1_2[k]
                 + f_3 * pc_z[k] * isg_5[k];

        t_10[k] = f_0 * hsg_10[k]
                  + f_3 * pc_x[k] * isg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, hsg_12, hsg_14, isg_6, \
                         isg_9, isg_12, isg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * isg_6[k];

        t_12[k] = f_0 * hsg_12[k]
                  + f_3 * pc_x[k] * isg_12[k];

        t_13[k] = f_3 * pc_y[k] * isg_9[k];

        t_14[k] = f_0 * hsg_14[k]
                  + f_3 * pc_x[k] * isg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, isf0_6, isf0_8, isf0_9, isf1_6, \
                         isf1_8, isf1_9, isg_10, isg_12, isg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * isf0_6[k]
                  - f_2 * isf1_6[k]
                  + f_3 * pc_y[k] * isg_10[k];

        t_16[k] = f_3 * pc_z[k] * isg_10[k];

        t_17[k] = f_6 * isf0_8[k]
                  - f_7 * isf1_8[k]
                  + f_3 * pc_y[k] * isg_12[k];

        t_18[k] = f_4 * isf0_9[k]
                  - f_5 * isf1_9[k]
                  + f_3 * pc_y[k] * isg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, hsh0_0, hsg_0, \
                         hsh1_0, isf0_9, isf1_9, isg_14, isg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * isg_14[k];

        t_20[k] = f_1 * isf0_9[k]
                  - f_2 * isf1_9[k]
                  + f_3 * pc_z[k] * isg_14[k];

        t_21[k] = pa_y[k] * hsh0_0[k]
                  - f_8 * pc_y[k] * hsh1_0[k];

        t_22[k] = f_9 * hsg_0[k]
                  + f_3 * pc_y[k] * isg_15[k];

        t_23[k] = f_3 * pc_z[k] * isg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pa_y, pc_y, pc_z, hsh0_3, hsh0_5, hsh0_6, \
                         hsg_1, hsg_3, hsh1_3, hsh1_5, hsh1_6, isg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_y[k] * hsh0_3[k]
                  + f_10 * hsg_1[k]
                  - f_8 * pc_y[k] * hsh1_3[k];

        t_25[k] = f_3 * pc_z[k] * isg_16[k];

        t_26[k] = pa_y[k] * hsh0_5[k]
                  - f_8 * pc_y[k] * hsh1_5[k];

        t_27[k] = pa_y[k] * hsh0_6[k]
                  + f_11 * hsg_3[k]
                  - f_8 * pc_y[k] * hsh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pa_y, pc_x, pc_y, pc_z, hsh0_9, hsg_5, \
                         hsg_25, hsh1_9, isg_18, isg_20, isg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * isg_18[k];

        t_29[k] = f_9 * hsg_5[k]
                  + f_3 * pc_y[k] * isg_20[k];

        t_30[k] = pa_y[k] * hsh0_9[k]
                  - f_8 * pc_y[k] * hsh1_9[k];

        t_31[k] = f_12 * hsg_25[k]
                  + f_3 * pc_x[k] * isg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_z, hsg_27, hsg_28, hsg_29, isg_21, \
                         isg_27, isg_28, isg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_3 * pc_z[k] * isg_21[k];

        t_33[k] = f_12 * hsg_27[k]
                  + f_3 * pc_x[k] * isg_27[k];

        t_34[k] = f_12 * hsg_28[k]
                  + f_3 * pc_x[k] * isg_28[k];

        t_35[k] = f_12 * hsg_29[k]
                  + f_3 * pc_x[k] * isg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, hsg_10, isf0_16, isf0_17, \
                         isf1_16, isf1_17, isg_25, isg_26, isg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * hsg_10[k]
                  + f_1 * isf0_16[k]
                  - f_2 * isf1_16[k]
                  + f_3 * pc_y[k] * isg_25[k];

        t_37[k] = f_3 * pc_z[k] * isg_25[k];

        t_38[k] = f_4 * isf0_16[k]
                  - f_5 * isf1_16[k]
                  + f_3 * pc_z[k] * isg_26[k];

        t_39[k] = f_6 * isf0_17[k]
                  - f_7 * isf1_17[k]
                  + f_3 * pc_z[k] * isg_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pa_z, pc_y, pc_z, hsh0_0, hsh0_20, \
                         hsg_14, hsh1_0, hsh1_20, isg_29, isg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_9 * hsg_14[k]
                  + f_3 * pc_y[k] * isg_29[k];

        t_41[k] = pa_y[k] * hsh0_20[k]
                  - f_8 * pc_y[k] * hsh1_20[k];

        t_42[k] = pa_z[k] * hsh0_0[k]
                  - f_8 * pc_z[k] * hsh1_0[k];

        t_43[k] = f_3 * pc_y[k] * isg_30[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_y, pc_z, hsh0_3, hsh0_5, hsg_0, \
                         hsg_2, hsh1_3, hsh1_5, isg_30, isg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * hsg_0[k]
                  + f_3 * pc_z[k] * isg_30[k];

        t_45[k] = pa_z[k] * hsh0_3[k]
                  - f_8 * pc_z[k] * hsh1_3[k];

        t_46[k] = f_3 * pc_y[k] * isg_32[k];

        t_47[k] = pa_z[k] * hsh0_5[k]
                  + f_10 * hsg_2[k]
                  - f_8 * pc_z[k] * hsh1_5[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_z, pc_y, pc_z, hsh0_6, hsh0_9, hsg_5, \
                         hsh1_6, hsh1_9, isf0_22, isf1_22, isg_34, \
                         isg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_z[k] * hsh0_6[k]
                  - f_8 * pc_z[k] * hsh1_6[k];

        t_49[k] = f_4 * isf0_22[k]
                  - f_5 * isf1_22[k]
                  + f_3 * pc_y[k] * isg_34[k];

        t_50[k] = f_3 * pc_y[k] * isg_35[k];

        t_51[k] = pa_z[k] * hsh0_9[k]
                  + f_11 * hsg_5[k]
                  - f_8 * pc_z[k] * hsh1_9[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, t_56, pc_x, pc_y, hsg_40, hsg_41, hsg_42, \
                         hsg_44, isg_39, isg_40, isg_41, isg_42, \
                         isg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_12 * hsg_40[k]
                  + f_3 * pc_x[k] * isg_40[k];

        t_53[k] = f_12 * hsg_41[k]
                  + f_3 * pc_x[k] * isg_41[k];

        t_54[k] = f_12 * hsg_42[k]
                  + f_3 * pc_x[k] * isg_42[k];

        t_55[k] = f_3 * pc_y[k] * isg_39[k];

        t_56[k] = f_12 * hsg_44[k]
                  + f_3 * pc_x[k] * isg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pa_z, pc_y, pc_z, hsh0_15, hsh1_15, isf0_27, \
                         isf0_28, isf1_27, isf1_28, isg_41, isg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pa_z[k] * hsh0_15[k]
                  - f_8 * pc_z[k] * hsh1_15[k];

        t_58[k] = f_13 * isf0_27[k]
                  - f_14 * isf1_27[k]
                  + f_3 * pc_y[k] * isg_41[k];

        t_59[k] = f_6 * isf0_28[k]
                  - f_7 * isf1_28[k]
                  + f_3 * pc_y[k] * isg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, pc_z, hsg_14, hsg_45, isf0_29, \
                         isf0_30, isf1_29, isf1_30, isg_43, isg_44, \
                         isg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_4 * isf0_29[k]
                  - f_5 * isf1_29[k]
                  + f_3 * pc_y[k] * isg_43[k];

        t_61[k] = f_3 * pc_y[k] * isg_44[k];

        t_62[k] = f_9 * hsg_14[k]
                  + f_1 * isf0_29[k]
                  - f_2 * isf1_29[k]
                  + f_3 * pc_z[k] * isg_44[k];

        t_63[k] = f_15 * hsg_45[k]
                  + f_1 * isf0_30[k]
                  - f_2 * isf1_30[k]
                  + f_3 * pc_x[k] * isg_45[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pc_x, pc_y, pc_z, hsg_15, hsg_48, isf0_33, \
                         isf1_33, isg_45, isg_46, isg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_10 * hsg_15[k]
                  + f_3 * pc_y[k] * isg_45[k];

        t_65[k] = f_3 * pc_z[k] * isg_45[k];

        t_66[k] = f_15 * hsg_48[k]
                  + f_6 * isf0_33[k]
                  - f_7 * isf1_33[k]
                  + f_3 * pc_x[k] * isg_48[k];

        t_67[k] = f_3 * pc_z[k] * isg_46[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pc_x, pc_z, hsg_51, isf0_30, isf0_36, isf1_30, \
                         isf1_36, isg_47, isg_48, isg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_4 * isf0_30[k]
                  - f_5 * isf1_30[k]
                  + f_3 * pc_z[k] * isg_47[k];

        t_69[k] = f_15 * hsg_51[k]
                  + f_4 * isf0_36[k]
                  - f_5 * isf1_36[k]
                  + f_3 * pc_x[k] * isg_51[k];

        t_70[k] = f_3 * pc_z[k] * isg_48[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, hsg_20, hsg_55, isf0_32, \
                         isf1_32, isg_50, isg_51, isg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_10 * hsg_20[k]
                  + f_3 * pc_y[k] * isg_50[k];

        t_72[k] = f_6 * isf0_32[k]
                  - f_7 * isf1_32[k]
                  + f_3 * pc_z[k] * isg_50[k];

        t_73[k] = f_15 * hsg_55[k]
                  + f_3 * pc_x[k] * isg_55[k];

        t_74[k] = f_3 * pc_z[k] * isg_51[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, hsg_25, hsg_57, hsg_58, hsg_59, \
                         isf0_36, isf1_36, isg_55, isg_57, isg_58, \
                         isg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_15 * hsg_57[k]
                  + f_3 * pc_x[k] * isg_57[k];

        t_76[k] = f_15 * hsg_58[k]
                  + f_3 * pc_x[k] * isg_58[k];

        t_77[k] = f_15 * hsg_59[k]
                  + f_3 * pc_x[k] * isg_59[k];

        t_78[k] = f_10 * hsg_25[k]
                  + f_1 * isf0_36[k]
                  - f_2 * isf1_36[k]
                  + f_3 * pc_y[k] * isg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_y, pc_z, hsg_29, isf0_36, isf0_37, \
                         isf1_36, isf1_37, isg_55, isg_56, isg_57, \
                         isg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * isg_55[k];

        t_80[k] = f_4 * isf0_36[k]
                  - f_5 * isf1_36[k]
                  + f_3 * pc_z[k] * isg_56[k];

        t_81[k] = f_6 * isf0_37[k]
                  - f_7 * isf1_37[k]
                  + f_3 * pc_z[k] * isg_57[k];

        t_82[k] = f_10 * hsg_29[k]
                  + f_3 * pc_y[k] * isg_59[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pa_y, pc_y, pc_z, hsh0_42, hsg_15, hsg_30, \
                         hsh1_42, isf0_39, isf1_39, isg_59, isg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_1 * isf0_39[k]
                  - f_2 * isf1_39[k]
                  + f_3 * pc_z[k] * isg_59[k];

        t_84[k] = pa_y[k] * hsh0_42[k]
                  - f_8 * pc_y[k] * hsh1_42[k];

        t_85[k] = f_9 * hsg_30[k]
                  + f_3 * pc_y[k] * isg_60[k];

        t_86[k] = f_9 * hsg_15[k]
                  + f_3 * pc_z[k] * isg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, hsh0_24, hsh0_27, \
                         hsh0_47, hsg_32, hsh1_24, hsh1_27, hsh1_47, \
                         isg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pa_z[k] * hsh0_24[k]
                  - f_8 * pc_z[k] * hsh1_24[k];

        t_88[k] = f_9 * hsg_32[k]
                  + f_3 * pc_y[k] * isg_62[k];

        t_89[k] = pa_y[k] * hsh0_47[k]
                  - f_8 * pc_y[k] * hsh1_47[k];

        t_90[k] = pa_z[k] * hsh0_27[k]
                  - f_8 * pc_z[k] * hsh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, pc_z, hsh0_51, hsg_18, \
                         hsg_35, hsg_70, hsh1_51, isg_63, isg_65, \
                         isg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * hsg_18[k]
                  + f_3 * pc_z[k] * isg_63[k];

        t_92[k] = f_9 * hsg_35[k]
                  + f_3 * pc_y[k] * isg_65[k];

        t_93[k] = pa_y[k] * hsh0_51[k]
                  - f_8 * pc_y[k] * hsh1_51[k];

        t_94[k] = f_15 * hsg_70[k]
                  + f_3 * pc_x[k] * isg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, hsg_71, hsg_72, hsg_73, hsg_74, isg_71, \
                         isg_72, isg_73, isg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_15 * hsg_71[k]
                  + f_3 * pc_x[k] * isg_71[k];

        t_96[k] = f_15 * hsg_72[k]
                  + f_3 * pc_x[k] * isg_72[k];

        t_97[k] = f_15 * hsg_73[k]
                  + f_3 * pc_x[k] * isg_73[k];

        t_98[k] = f_15 * hsg_74[k]
                  + f_3 * pc_x[k] * isg_74[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, hsh0_36, hsg_25, hsg_42, \
                         hsh1_36, isf0_48, isf1_48, isg_70, isg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * hsh0_36[k]
                  - f_8 * pc_z[k] * hsh1_36[k];

        t_100[k] = f_9 * hsg_25[k]
                   + f_3 * pc_z[k] * isg_70[k];

        t_101[k] = f_9 * hsg_42[k]
                   + f_6 * isf0_48[k]
                   - f_7 * isf1_48[k]
                   + f_3 * pc_y[k] * isg_72[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, pa_y, pc_y, hsh0_62, hsg_43, hsg_44, hsh1_62, \
                         isf0_49, isf1_49, isg_73, isg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_9 * hsg_43[k]
                   + f_4 * isf0_49[k]
                   - f_5 * isf1_49[k]
                   + f_3 * pc_y[k] * isg_73[k];

        t_103[k] = f_9 * hsg_44[k]
                   + f_3 * pc_y[k] * isg_74[k];

        t_104[k] = pa_y[k] * hsh0_62[k]
                   - f_8 * pc_y[k] * hsh1_62[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, hsg_30, hsg_75, \
                         isf0_50, isf1_50, isg_75, isg_76, isg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_15 * hsg_75[k]
                   + f_1 * isf0_50[k]
                   - f_2 * isf1_50[k]
                   + f_3 * pc_x[k] * isg_75[k];

        t_106[k] = f_3 * pc_y[k] * isg_75[k];

        t_107[k] = f_10 * hsg_30[k]
                   + f_3 * pc_z[k] * isg_75[k];

        t_108[k] = f_4 * isf0_50[k]
                   - f_5 * isf1_50[k]
                   + f_3 * pc_y[k] * isg_76[k];

        t_109[k] = f_3 * pc_y[k] * isg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_y, hsg_80, isf0_51, isf0_52, \
                         isf0_55, isf1_51, isf1_52, isf1_55, isg_78, isg_79, \
                         isg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_15 * hsg_80[k]
                   + f_6 * isf0_55[k]
                   - f_7 * isf1_55[k]
                   + f_3 * pc_x[k] * isg_80[k];

        t_111[k] = f_6 * isf0_51[k]
                   - f_7 * isf1_51[k]
                   + f_3 * pc_y[k] * isg_78[k];

        t_112[k] = f_4 * isf0_52[k]
                   - f_5 * isf1_52[k]
                   + f_3 * pc_y[k] * isg_79[k];

        t_113[k] = f_3 * pc_y[k] * isg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, hsg_84, hsg_85, hsg_86, hsg_87, \
                         isf0_59, isf1_59, isg_84, isg_85, isg_86, \
                         isg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_15 * hsg_84[k]
                   + f_4 * isf0_59[k]
                   - f_5 * isf1_59[k]
                   + f_3 * pc_x[k] * isg_84[k];

        t_115[k] = f_15 * hsg_85[k]
                   + f_3 * pc_x[k] * isg_85[k];

        t_116[k] = f_15 * hsg_86[k]
                   + f_3 * pc_x[k] * isg_86[k];

        t_117[k] = f_15 * hsg_87[k]
                   + f_3 * pc_x[k] * isg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, hsg_89, isf0_56, isf0_57, \
                         isf1_56, isf1_57, isg_84, isg_85, isg_86, \
                         isg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_3 * pc_y[k] * isg_84[k];

        t_119[k] = f_15 * hsg_89[k]
                   + f_3 * pc_x[k] * isg_89[k];

        t_120[k] = f_1 * isf0_56[k]
                   - f_2 * isf1_56[k]
                   + f_3 * pc_y[k] * isg_85[k];

        t_121[k] = f_13 * isf0_57[k]
                   - f_14 * isf1_57[k]
                   + f_3 * pc_y[k] * isg_86[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, hsg_44, isf0_58, isf0_59, \
                         isf1_58, isf1_59, isg_87, isg_88, isg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_6 * isf0_58[k]
                   - f_7 * isf1_58[k]
                   + f_3 * pc_y[k] * isg_87[k];

        t_123[k] = f_4 * isf0_59[k]
                   - f_5 * isf1_59[k]
                   + f_3 * pc_y[k] * isg_88[k];

        t_124[k] = f_3 * pc_y[k] * isg_89[k];

        t_125[k] = f_10 * hsg_44[k]
                   + f_1 * isf0_59[k]
                   - f_2 * isf1_59[k]
                   + f_3 * pc_z[k] * isg_89[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, hsg_45, hsg_90, hsg_93, \
                         isf0_60, isf0_63, isf1_60, isf1_63, isg_90, \
                         isg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_11 * hsg_90[k]
                   + f_1 * isf0_60[k]
                   - f_2 * isf1_60[k]
                   + f_3 * pc_x[k] * isg_90[k];

        t_127[k] = f_11 * hsg_45[k]
                   + f_3 * pc_y[k] * isg_90[k];

        t_128[k] = f_3 * pc_z[k] * isg_90[k];

        t_129[k] = f_11 * hsg_93[k]
                   + f_6 * isf0_63[k]
                   - f_7 * isf1_63[k]
                   + f_3 * pc_x[k] * isg_93[k];
    }
}

static auto
compute_prim_ish_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsh0,
                                                          const size_t hsg, const size_t hsh1,
                                                          const size_t isf0, const size_t isf1,
                                                          const size_t isg, const size_t ncols,
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
    const auto f_15 = 2.0 / q;

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

    const auto *hsh0_63 = buffer.data(hsh0 + 63);
    const auto *hsh0_66 = buffer.data(hsh0 + 66);
    const auto *hsh0_69 = buffer.data(hsh0 + 69);
    const auto *hsh0_78 = buffer.data(hsh0 + 78);
    const auto *hsh0_105 = buffer.data(hsh0 + 105);
    const auto *hsh0_108 = buffer.data(hsh0 + 108);
    const auto *hsh0_110 = buffer.data(hsh0 + 110);
    const auto *hsh0_111 = buffer.data(hsh0 + 111);
    const auto *hsh0_114 = buffer.data(hsh0 + 114);
    const auto *hsh0_125 = buffer.data(hsh0 + 125);
    const auto *hsh0_126 = buffer.data(hsh0 + 126);
    const auto *hsh0_129 = buffer.data(hsh0 + 129);
    const auto *hsh0_132 = buffer.data(hsh0 + 132);
    const auto *hsh0_141 = buffer.data(hsh0 + 141);

    const auto *hsg_45 = buffer.data(hsg + 45);
    const auto *hsg_48 = buffer.data(hsg + 48);
    const auto *hsg_50 = buffer.data(hsg + 50);
    const auto *hsg_55 = buffer.data(hsg + 55);
    const auto *hsg_59 = buffer.data(hsg + 59);
    const auto *hsg_60 = buffer.data(hsg + 60);
    const auto *hsg_62 = buffer.data(hsg + 62);
    const auto *hsg_63 = buffer.data(hsg + 63);
    const auto *hsg_65 = buffer.data(hsg + 65);
    const auto *hsg_70 = buffer.data(hsg + 70);
    const auto *hsg_72 = buffer.data(hsg + 72);
    const auto *hsg_73 = buffer.data(hsg + 73);
    const auto *hsg_74 = buffer.data(hsg + 74);
    const auto *hsg_75 = buffer.data(hsg + 75);
    const auto *hsg_76 = buffer.data(hsg + 76);
    const auto *hsg_77 = buffer.data(hsg + 77);
    const auto *hsg_78 = buffer.data(hsg + 78);
    const auto *hsg_80 = buffer.data(hsg + 80);
    const auto *hsg_85 = buffer.data(hsg + 85);
    const auto *hsg_87 = buffer.data(hsg + 87);
    const auto *hsg_88 = buffer.data(hsg + 88);
    const auto *hsg_89 = buffer.data(hsg + 89);
    const auto *hsg_90 = buffer.data(hsg + 90);
    const auto *hsg_93 = buffer.data(hsg + 93);
    const auto *hsg_95 = buffer.data(hsg + 95);
    const auto *hsg_96 = buffer.data(hsg + 96);
    const auto *hsg_100 = buffer.data(hsg + 100);
    const auto *hsg_102 = buffer.data(hsg + 102);
    const auto *hsg_103 = buffer.data(hsg + 103);
    const auto *hsg_104 = buffer.data(hsg + 104);
    const auto *hsg_105 = buffer.data(hsg + 105);
    const auto *hsg_107 = buffer.data(hsg + 107);
    const auto *hsg_110 = buffer.data(hsg + 110);
    const auto *hsg_114 = buffer.data(hsg + 114);
    const auto *hsg_115 = buffer.data(hsg + 115);
    const auto *hsg_116 = buffer.data(hsg + 116);
    const auto *hsg_117 = buffer.data(hsg + 117);
    const auto *hsg_118 = buffer.data(hsg + 118);
    const auto *hsg_119 = buffer.data(hsg + 119);
    const auto *hsg_130 = buffer.data(hsg + 130);
    const auto *hsg_131 = buffer.data(hsg + 131);
    const auto *hsg_132 = buffer.data(hsg + 132);
    const auto *hsg_133 = buffer.data(hsg + 133);
    const auto *hsg_134 = buffer.data(hsg + 134);
    const auto *hsg_135 = buffer.data(hsg + 135);
    const auto *hsg_140 = buffer.data(hsg + 140);
    const auto *hsg_144 = buffer.data(hsg + 144);
    const auto *hsg_145 = buffer.data(hsg + 145);
    const auto *hsg_146 = buffer.data(hsg + 146);
    const auto *hsg_147 = buffer.data(hsg + 147);
    const auto *hsg_149 = buffer.data(hsg + 149);
    const auto *hsg_150 = buffer.data(hsg + 150);
    const auto *hsg_153 = buffer.data(hsg + 153);
    const auto *hsg_156 = buffer.data(hsg + 156);
    const auto *hsg_160 = buffer.data(hsg + 160);
    const auto *hsg_162 = buffer.data(hsg + 162);
    const auto *hsg_163 = buffer.data(hsg + 163);
    const auto *hsg_164 = buffer.data(hsg + 164);
    const auto *hsg_170 = buffer.data(hsg + 170);
    const auto *hsg_174 = buffer.data(hsg + 174);
    const auto *hsg_175 = buffer.data(hsg + 175);
    const auto *hsg_176 = buffer.data(hsg + 176);
    const auto *hsg_177 = buffer.data(hsg + 177);
    const auto *hsg_178 = buffer.data(hsg + 178);
    const auto *hsg_179 = buffer.data(hsg + 179);

    const auto *hsh1_63 = buffer.data(hsh1 + 63);
    const auto *hsh1_66 = buffer.data(hsh1 + 66);
    const auto *hsh1_69 = buffer.data(hsh1 + 69);
    const auto *hsh1_78 = buffer.data(hsh1 + 78);
    const auto *hsh1_105 = buffer.data(hsh1 + 105);
    const auto *hsh1_108 = buffer.data(hsh1 + 108);
    const auto *hsh1_110 = buffer.data(hsh1 + 110);
    const auto *hsh1_111 = buffer.data(hsh1 + 111);
    const auto *hsh1_114 = buffer.data(hsh1 + 114);
    const auto *hsh1_125 = buffer.data(hsh1 + 125);
    const auto *hsh1_126 = buffer.data(hsh1 + 126);
    const auto *hsh1_129 = buffer.data(hsh1 + 129);
    const auto *hsh1_132 = buffer.data(hsh1 + 132);
    const auto *hsh1_141 = buffer.data(hsh1 + 141);

    const auto *isf0_60 = buffer.data(isf0 + 60);
    const auto *isf0_62 = buffer.data(isf0 + 62);
    const auto *isf0_66 = buffer.data(isf0 + 66);
    const auto *isf0_67 = buffer.data(isf0 + 67);
    const auto *isf0_69 = buffer.data(isf0 + 69);
    const auto *isf0_75 = buffer.data(isf0 + 75);
    const auto *isf0_78 = buffer.data(isf0 + 78);
    const auto *isf0_79 = buffer.data(isf0 + 79);
    const auto *isf0_86 = buffer.data(isf0 + 86);
    const auto *isf0_88 = buffer.data(isf0 + 88);
    const auto *isf0_89 = buffer.data(isf0 + 89);
    const auto *isf0_90 = buffer.data(isf0 + 90);
    const auto *isf0_91 = buffer.data(isf0 + 91);
    const auto *isf0_92 = buffer.data(isf0 + 92);
    const auto *isf0_95 = buffer.data(isf0 + 95);
    const auto *isf0_96 = buffer.data(isf0 + 96);
    const auto *isf0_97 = buffer.data(isf0 + 97);
    const auto *isf0_98 = buffer.data(isf0 + 98);
    const auto *isf0_99 = buffer.data(isf0 + 99);
    const auto *isf0_100 = buffer.data(isf0 + 100);
    const auto *isf0_102 = buffer.data(isf0 + 102);
    const auto *isf0_103 = buffer.data(isf0 + 103);
    const auto *isf0_106 = buffer.data(isf0 + 106);
    const auto *isf0_107 = buffer.data(isf0 + 107);
    const auto *isf0_109 = buffer.data(isf0 + 109);
    const auto *isf0_115 = buffer.data(isf0 + 115);
    const auto *isf0_118 = buffer.data(isf0 + 118);
    const auto *isf0_119 = buffer.data(isf0 + 119);

    const auto *isf1_60 = buffer.data(isf1 + 60);
    const auto *isf1_62 = buffer.data(isf1 + 62);
    const auto *isf1_66 = buffer.data(isf1 + 66);
    const auto *isf1_67 = buffer.data(isf1 + 67);
    const auto *isf1_69 = buffer.data(isf1 + 69);
    const auto *isf1_75 = buffer.data(isf1 + 75);
    const auto *isf1_78 = buffer.data(isf1 + 78);
    const auto *isf1_79 = buffer.data(isf1 + 79);
    const auto *isf1_86 = buffer.data(isf1 + 86);
    const auto *isf1_88 = buffer.data(isf1 + 88);
    const auto *isf1_89 = buffer.data(isf1 + 89);
    const auto *isf1_90 = buffer.data(isf1 + 90);
    const auto *isf1_91 = buffer.data(isf1 + 91);
    const auto *isf1_92 = buffer.data(isf1 + 92);
    const auto *isf1_95 = buffer.data(isf1 + 95);
    const auto *isf1_96 = buffer.data(isf1 + 96);
    const auto *isf1_97 = buffer.data(isf1 + 97);
    const auto *isf1_98 = buffer.data(isf1 + 98);
    const auto *isf1_99 = buffer.data(isf1 + 99);
    const auto *isf1_100 = buffer.data(isf1 + 100);
    const auto *isf1_102 = buffer.data(isf1 + 102);
    const auto *isf1_103 = buffer.data(isf1 + 103);
    const auto *isf1_106 = buffer.data(isf1 + 106);
    const auto *isf1_107 = buffer.data(isf1 + 107);
    const auto *isf1_109 = buffer.data(isf1 + 109);
    const auto *isf1_115 = buffer.data(isf1 + 115);
    const auto *isf1_118 = buffer.data(isf1 + 118);
    const auto *isf1_119 = buffer.data(isf1 + 119);

    const auto *isg_91 = buffer.data(isg + 91);
    const auto *isg_92 = buffer.data(isg + 92);
    const auto *isg_93 = buffer.data(isg + 93);
    const auto *isg_95 = buffer.data(isg + 95);
    const auto *isg_96 = buffer.data(isg + 96);
    const auto *isg_100 = buffer.data(isg + 100);
    const auto *isg_101 = buffer.data(isg + 101);
    const auto *isg_102 = buffer.data(isg + 102);
    const auto *isg_103 = buffer.data(isg + 103);
    const auto *isg_104 = buffer.data(isg + 104);
    const auto *isg_105 = buffer.data(isg + 105);
    const auto *isg_107 = buffer.data(isg + 107);
    const auto *isg_108 = buffer.data(isg + 108);
    const auto *isg_110 = buffer.data(isg + 110);
    const auto *isg_114 = buffer.data(isg + 114);
    const auto *isg_115 = buffer.data(isg + 115);
    const auto *isg_116 = buffer.data(isg + 116);
    const auto *isg_117 = buffer.data(isg + 117);
    const auto *isg_118 = buffer.data(isg + 118);
    const auto *isg_119 = buffer.data(isg + 119);
    const auto *isg_120 = buffer.data(isg + 120);
    const auto *isg_122 = buffer.data(isg + 122);
    const auto *isg_123 = buffer.data(isg + 123);
    const auto *isg_125 = buffer.data(isg + 125);
    const auto *isg_130 = buffer.data(isg + 130);
    const auto *isg_131 = buffer.data(isg + 131);
    const auto *isg_132 = buffer.data(isg + 132);
    const auto *isg_133 = buffer.data(isg + 133);
    const auto *isg_134 = buffer.data(isg + 134);
    const auto *isg_135 = buffer.data(isg + 135);
    const auto *isg_136 = buffer.data(isg + 136);
    const auto *isg_137 = buffer.data(isg + 137);
    const auto *isg_138 = buffer.data(isg + 138);
    const auto *isg_139 = buffer.data(isg + 139);
    const auto *isg_140 = buffer.data(isg + 140);
    const auto *isg_144 = buffer.data(isg + 144);
    const auto *isg_145 = buffer.data(isg + 145);
    const auto *isg_146 = buffer.data(isg + 146);
    const auto *isg_147 = buffer.data(isg + 147);
    const auto *isg_148 = buffer.data(isg + 148);
    const auto *isg_149 = buffer.data(isg + 149);
    const auto *isg_150 = buffer.data(isg + 150);
    const auto *isg_151 = buffer.data(isg + 151);
    const auto *isg_152 = buffer.data(isg + 152);
    const auto *isg_153 = buffer.data(isg + 153);
    const auto *isg_155 = buffer.data(isg + 155);
    const auto *isg_156 = buffer.data(isg + 156);
    const auto *isg_160 = buffer.data(isg + 160);
    const auto *isg_161 = buffer.data(isg + 161);
    const auto *isg_162 = buffer.data(isg + 162);
    const auto *isg_163 = buffer.data(isg + 163);
    const auto *isg_164 = buffer.data(isg + 164);
    const auto *isg_165 = buffer.data(isg + 165);
    const auto *isg_167 = buffer.data(isg + 167);
    const auto *isg_168 = buffer.data(isg + 168);
    const auto *isg_170 = buffer.data(isg + 170);
    const auto *isg_174 = buffer.data(isg + 174);
    const auto *isg_175 = buffer.data(isg + 175);
    const auto *isg_176 = buffer.data(isg + 176);
    const auto *isg_177 = buffer.data(isg + 177);
    const auto *isg_178 = buffer.data(isg + 178);
    const auto *isg_179 = buffer.data(isg + 179);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_x, pc_z, hsg_96, isf0_60, isf0_66, \
                         isf1_60, isf1_66, isg_91, isg_92, isg_93, \
                         isg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_3 * pc_z[k] * isg_91[k];

        t_131[k] = f_4 * isf0_60[k]
                   - f_5 * isf1_60[k]
                   + f_3 * pc_z[k] * isg_92[k];

        t_132[k] = f_11 * hsg_96[k]
                   + f_4 * isf0_66[k]
                   - f_5 * isf1_66[k]
                   + f_3 * pc_x[k] * isg_96[k];

        t_133[k] = f_3 * pc_z[k] * isg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, hsg_50, hsg_100, \
                         isf0_62, isf1_62, isg_95, isg_96, isg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_11 * hsg_50[k]
                   + f_3 * pc_y[k] * isg_95[k];

        t_135[k] = f_6 * isf0_62[k]
                   - f_7 * isf1_62[k]
                   + f_3 * pc_z[k] * isg_95[k];

        t_136[k] = f_11 * hsg_100[k]
                   + f_3 * pc_x[k] * isg_100[k];

        t_137[k] = f_3 * pc_z[k] * isg_96[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pc_x, pc_y, hsg_55, hsg_102, hsg_103, \
                         hsg_104, isf0_66, isf1_66, isg_100, isg_102, isg_103, \
                         isg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_11 * hsg_102[k]
                   + f_3 * pc_x[k] * isg_102[k];

        t_139[k] = f_11 * hsg_103[k]
                   + f_3 * pc_x[k] * isg_103[k];

        t_140[k] = f_11 * hsg_104[k]
                   + f_3 * pc_x[k] * isg_104[k];

        t_141[k] = f_11 * hsg_55[k]
                   + f_1 * isf0_66[k]
                   - f_2 * isf1_66[k]
                   + f_3 * pc_y[k] * isg_100[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, pc_z, hsg_59, isf0_66, isf0_67, \
                         isf1_66, isf1_67, isg_100, isg_101, isg_102, \
                         isg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_3 * pc_z[k] * isg_100[k];

        t_143[k] = f_4 * isf0_66[k]
                   - f_5 * isf1_66[k]
                   + f_3 * pc_z[k] * isg_101[k];

        t_144[k] = f_6 * isf0_67[k]
                   - f_7 * isf1_67[k]
                   + f_3 * pc_z[k] * isg_102[k];

        t_145[k] = f_11 * hsg_59[k]
                   + f_3 * pc_y[k] * isg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, hsh0_63, hsg_45, \
                         hsg_60, hsh1_63, isf0_69, isf1_69, isg_104, \
                         isg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * isf0_69[k]
                   - f_2 * isf1_69[k]
                   + f_3 * pc_z[k] * isg_104[k];

        t_147[k] = pa_z[k] * hsh0_63[k]
                   - f_8 * pc_z[k] * hsh1_63[k];

        t_148[k] = f_10 * hsg_60[k]
                   + f_3 * pc_y[k] * isg_105[k];

        t_149[k] = f_9 * hsg_45[k]
                   + f_3 * pc_z[k] * isg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_z, pc_x, pc_y, pc_z, hsh0_66, hsg_62, \
                         hsg_110, hsh1_66, isf0_75, isf1_75, isg_107, \
                         isg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * hsh0_66[k]
                   - f_8 * pc_z[k] * hsh1_66[k];

        t_151[k] = f_10 * hsg_62[k]
                   + f_3 * pc_y[k] * isg_107[k];

        t_152[k] = f_11 * hsg_110[k]
                   + f_6 * isf0_75[k]
                   - f_7 * isf1_75[k]
                   + f_3 * pc_x[k] * isg_110[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pa_z, pc_y, pc_z, hsh0_69, hsg_48, hsg_65, \
                         hsh1_69, isg_108, isg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = pa_z[k] * hsh0_69[k]
                   - f_8 * pc_z[k] * hsh1_69[k];

        t_154[k] = f_9 * hsg_48[k]
                   + f_3 * pc_z[k] * isg_108[k];

        t_155[k] = f_10 * hsg_65[k]
                   + f_3 * pc_y[k] * isg_110[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, hsg_114, hsg_115, hsg_116, hsg_117, \
                         isf0_79, isf1_79, isg_114, isg_115, isg_116, \
                         isg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_11 * hsg_114[k]
                   + f_4 * isf0_79[k]
                   - f_5 * isf1_79[k]
                   + f_3 * pc_x[k] * isg_114[k];

        t_157[k] = f_11 * hsg_115[k]
                   + f_3 * pc_x[k] * isg_115[k];

        t_158[k] = f_11 * hsg_116[k]
                   + f_3 * pc_x[k] * isg_116[k];

        t_159[k] = f_11 * hsg_117[k]
                   + f_3 * pc_x[k] * isg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pc_x, pc_z, hsh0_78, hsg_55, \
                         hsg_118, hsg_119, hsh1_78, isg_115, isg_118, \
                         isg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_11 * hsg_118[k]
                   + f_3 * pc_x[k] * isg_118[k];

        t_161[k] = f_11 * hsg_119[k]
                   + f_3 * pc_x[k] * isg_119[k];

        t_162[k] = pa_z[k] * hsh0_78[k]
                   - f_8 * pc_z[k] * hsh1_78[k];

        t_163[k] = f_9 * hsg_55[k]
                   + f_3 * pc_z[k] * isg_115[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pc_y, hsg_72, hsg_73, hsg_74, isf0_78, isf0_79, \
                         isf1_78, isf1_79, isg_117, isg_118, isg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_10 * hsg_72[k]
                   + f_6 * isf0_78[k]
                   - f_7 * isf1_78[k]
                   + f_3 * pc_y[k] * isg_117[k];

        t_165[k] = f_10 * hsg_73[k]
                   + f_4 * isf0_79[k]
                   - f_5 * isf1_79[k]
                   + f_3 * pc_y[k] * isg_118[k];

        t_166[k] = f_10 * hsg_74[k]
                   + f_3 * pc_y[k] * isg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_y, pc_y, pc_z, hsh0_105, hsg_59, \
                         hsg_60, hsg_75, hsh1_105, isf0_79, isf1_79, isg_119, \
                         isg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * hsg_59[k]
                   + f_1 * isf0_79[k]
                   - f_2 * isf1_79[k]
                   + f_3 * pc_z[k] * isg_119[k];

        t_168[k] = pa_y[k] * hsh0_105[k]
                   - f_8 * pc_y[k] * hsh1_105[k];

        t_169[k] = f_9 * hsg_75[k]
                   + f_3 * pc_y[k] * isg_120[k];

        t_170[k] = f_10 * hsg_60[k]
                   + f_3 * pc_z[k] * isg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_y, pc_y, hsh0_108, hsh0_110, hsh0_111, \
                         hsg_76, hsg_77, hsg_78, hsh1_108, hsh1_110, hsh1_111, \
                         isg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_y[k] * hsh0_108[k]
                   + f_10 * hsg_76[k]
                   - f_8 * pc_y[k] * hsh1_108[k];

        t_172[k] = f_9 * hsg_77[k]
                   + f_3 * pc_y[k] * isg_122[k];

        t_173[k] = pa_y[k] * hsh0_110[k]
                   - f_8 * pc_y[k] * hsh1_110[k];

        t_174[k] = pa_y[k] * hsh0_111[k]
                   + f_11 * hsg_78[k]
                   - f_8 * pc_y[k] * hsh1_111[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_y, pc_x, pc_y, pc_z, hsh0_114, hsg_63, \
                         hsg_80, hsg_130, hsh1_114, isg_123, isg_125, \
                         isg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_10 * hsg_63[k]
                   + f_3 * pc_z[k] * isg_123[k];

        t_176[k] = f_9 * hsg_80[k]
                   + f_3 * pc_y[k] * isg_125[k];

        t_177[k] = pa_y[k] * hsh0_114[k]
                   - f_8 * pc_y[k] * hsh1_114[k];

        t_178[k] = f_11 * hsg_130[k]
                   + f_3 * pc_x[k] * isg_130[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, hsg_131, hsg_132, hsg_133, hsg_134, \
                         isg_131, isg_132, isg_133, isg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_11 * hsg_131[k]
                   + f_3 * pc_x[k] * isg_131[k];

        t_180[k] = f_11 * hsg_132[k]
                   + f_3 * pc_x[k] * isg_132[k];

        t_181[k] = f_11 * hsg_133[k]
                   + f_3 * pc_x[k] * isg_133[k];

        t_182[k] = f_11 * hsg_134[k]
                   + f_3 * pc_x[k] * isg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, hsg_70, hsg_85, hsg_87, isf0_86, \
                         isf0_88, isf1_86, isf1_88, isg_130, isg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * hsg_85[k]
                   + f_1 * isf0_86[k]
                   - f_2 * isf1_86[k]
                   + f_3 * pc_y[k] * isg_130[k];

        t_184[k] = f_10 * hsg_70[k]
                   + f_3 * pc_z[k] * isg_130[k];

        t_185[k] = f_9 * hsg_87[k]
                   + f_6 * isf0_88[k]
                   - f_7 * isf1_88[k]
                   + f_3 * pc_y[k] * isg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pa_y, pc_y, hsh0_125, hsg_88, hsg_89, hsh1_125, \
                         isf0_89, isf1_89, isg_133, isg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * hsg_88[k]
                   + f_4 * isf0_89[k]
                   - f_5 * isf1_89[k]
                   + f_3 * pc_y[k] * isg_133[k];

        t_187[k] = f_9 * hsg_89[k]
                   + f_3 * pc_y[k] * isg_134[k];

        t_188[k] = pa_y[k] * hsh0_125[k]
                   - f_8 * pc_y[k] * hsh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, hsg_75, hsg_135, \
                         isf0_90, isf1_90, isg_135, isg_136, isg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_11 * hsg_135[k]
                   + f_1 * isf0_90[k]
                   - f_2 * isf1_90[k]
                   + f_3 * pc_x[k] * isg_135[k];

        t_190[k] = f_3 * pc_y[k] * isg_135[k];

        t_191[k] = f_11 * hsg_75[k]
                   + f_3 * pc_z[k] * isg_135[k];

        t_192[k] = f_4 * isf0_90[k]
                   - f_5 * isf1_90[k]
                   + f_3 * pc_y[k] * isg_136[k];

        t_193[k] = f_3 * pc_y[k] * isg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pc_x, pc_y, hsg_140, isf0_91, isf0_92, \
                         isf0_95, isf1_91, isf1_92, isf1_95, isg_138, isg_139, \
                         isg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_11 * hsg_140[k]
                   + f_6 * isf0_95[k]
                   - f_7 * isf1_95[k]
                   + f_3 * pc_x[k] * isg_140[k];

        t_195[k] = f_6 * isf0_91[k]
                   - f_7 * isf1_91[k]
                   + f_3 * pc_y[k] * isg_138[k];

        t_196[k] = f_4 * isf0_92[k]
                   - f_5 * isf1_92[k]
                   + f_3 * pc_y[k] * isg_139[k];

        t_197[k] = f_3 * pc_y[k] * isg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pc_x, hsg_144, hsg_145, hsg_146, hsg_147, \
                         isf0_99, isf1_99, isg_144, isg_145, isg_146, \
                         isg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_11 * hsg_144[k]
                   + f_4 * isf0_99[k]
                   - f_5 * isf1_99[k]
                   + f_3 * pc_x[k] * isg_144[k];

        t_199[k] = f_11 * hsg_145[k]
                   + f_3 * pc_x[k] * isg_145[k];

        t_200[k] = f_11 * hsg_146[k]
                   + f_3 * pc_x[k] * isg_146[k];

        t_201[k] = f_11 * hsg_147[k]
                   + f_3 * pc_x[k] * isg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pc_x, pc_y, hsg_149, isf0_96, isf0_97, \
                         isf1_96, isf1_97, isg_144, isg_145, isg_146, \
                         isg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_3 * pc_y[k] * isg_144[k];

        t_203[k] = f_11 * hsg_149[k]
                   + f_3 * pc_x[k] * isg_149[k];

        t_204[k] = f_1 * isf0_96[k]
                   - f_2 * isf1_96[k]
                   + f_3 * pc_y[k] * isg_145[k];

        t_205[k] = f_13 * isf0_97[k]
                   - f_14 * isf1_97[k]
                   + f_3 * pc_y[k] * isg_146[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, hsg_89, isf0_98, isf0_99, \
                         isf1_98, isf1_99, isg_147, isg_148, isg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * isf0_98[k]
                   - f_7 * isf1_98[k]
                   + f_3 * pc_y[k] * isg_147[k];

        t_207[k] = f_4 * isf0_99[k]
                   - f_5 * isf1_99[k]
                   + f_3 * pc_y[k] * isg_148[k];

        t_208[k] = f_3 * pc_y[k] * isg_149[k];

        t_209[k] = f_11 * hsg_89[k]
                   + f_1 * isf0_99[k]
                   - f_2 * isf1_99[k]
                   + f_3 * pc_z[k] * isg_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, pc_y, pc_z, hsg_90, hsg_150, \
                         hsg_153, isf0_100, isf0_103, isf1_100, isf1_103, isg_150, \
                         isg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_10 * hsg_150[k]
                   + f_1 * isf0_100[k]
                   - f_2 * isf1_100[k]
                   + f_3 * pc_x[k] * isg_150[k];

        t_211[k] = f_15 * hsg_90[k]
                   + f_3 * pc_y[k] * isg_150[k];

        t_212[k] = f_3 * pc_z[k] * isg_150[k];

        t_213[k] = f_10 * hsg_153[k]
                   + f_6 * isf0_103[k]
                   - f_7 * isf1_103[k]
                   + f_3 * pc_x[k] * isg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, hsg_156, isf0_100, isf0_106, \
                         isf1_100, isf1_106, isg_151, isg_152, isg_153, \
                         isg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_3 * pc_z[k] * isg_151[k];

        t_215[k] = f_4 * isf0_100[k]
                   - f_5 * isf1_100[k]
                   + f_3 * pc_z[k] * isg_152[k];

        t_216[k] = f_10 * hsg_156[k]
                   + f_4 * isf0_106[k]
                   - f_5 * isf1_106[k]
                   + f_3 * pc_x[k] * isg_156[k];

        t_217[k] = f_3 * pc_z[k] * isg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, hsg_95, hsg_160, \
                         isf0_102, isf1_102, isg_155, isg_156, \
                         isg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_15 * hsg_95[k]
                   + f_3 * pc_y[k] * isg_155[k];

        t_219[k] = f_6 * isf0_102[k]
                   - f_7 * isf1_102[k]
                   + f_3 * pc_z[k] * isg_155[k];

        t_220[k] = f_10 * hsg_160[k]
                   + f_3 * pc_x[k] * isg_160[k];

        t_221[k] = f_3 * pc_z[k] * isg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, hsg_100, hsg_162, hsg_163, \
                         hsg_164, isf0_106, isf1_106, isg_160, isg_162, isg_163, \
                         isg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_10 * hsg_162[k]
                   + f_3 * pc_x[k] * isg_162[k];

        t_223[k] = f_10 * hsg_163[k]
                   + f_3 * pc_x[k] * isg_163[k];

        t_224[k] = f_10 * hsg_164[k]
                   + f_3 * pc_x[k] * isg_164[k];

        t_225[k] = f_15 * hsg_100[k]
                   + f_1 * isf0_106[k]
                   - f_2 * isf1_106[k]
                   + f_3 * pc_y[k] * isg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, hsg_104, isf0_106, isf0_107, \
                         isf1_106, isf1_107, isg_160, isg_161, isg_162, \
                         isg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * isg_160[k];

        t_227[k] = f_4 * isf0_106[k]
                   - f_5 * isf1_106[k]
                   + f_3 * pc_z[k] * isg_161[k];

        t_228[k] = f_6 * isf0_107[k]
                   - f_7 * isf1_107[k]
                   + f_3 * pc_z[k] * isg_162[k];

        t_229[k] = f_15 * hsg_104[k]
                   + f_3 * pc_y[k] * isg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pa_z, pc_y, pc_z, hsh0_126, hsg_90, \
                         hsg_105, hsh1_126, isf0_109, isf1_109, isg_164, \
                         isg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_1 * isf0_109[k]
                   - f_2 * isf1_109[k]
                   + f_3 * pc_z[k] * isg_164[k];

        t_231[k] = pa_z[k] * hsh0_126[k]
                   - f_8 * pc_z[k] * hsh1_126[k];

        t_232[k] = f_11 * hsg_105[k]
                   + f_3 * pc_y[k] * isg_165[k];

        t_233[k] = f_9 * hsg_90[k]
                   + f_3 * pc_z[k] * isg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pa_z, pc_x, pc_y, pc_z, hsh0_129, hsg_107, \
                         hsg_170, hsh1_129, isf0_115, isf1_115, isg_167, \
                         isg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pa_z[k] * hsh0_129[k]
                   - f_8 * pc_z[k] * hsh1_129[k];

        t_235[k] = f_11 * hsg_107[k]
                   + f_3 * pc_y[k] * isg_167[k];

        t_236[k] = f_10 * hsg_170[k]
                   + f_6 * isf0_115[k]
                   - f_7 * isf1_115[k]
                   + f_3 * pc_x[k] * isg_170[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_z, pc_y, pc_z, hsh0_132, hsg_93, hsg_110, \
                         hsh1_132, isg_168, isg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_z[k] * hsh0_132[k]
                   - f_8 * pc_z[k] * hsh1_132[k];

        t_238[k] = f_9 * hsg_93[k]
                   + f_3 * pc_z[k] * isg_168[k];

        t_239[k] = f_11 * hsg_110[k]
                   + f_3 * pc_y[k] * isg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pc_x, hsg_174, hsg_175, hsg_176, hsg_177, \
                         isf0_119, isf1_119, isg_174, isg_175, isg_176, \
                         isg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_10 * hsg_174[k]
                   + f_4 * isf0_119[k]
                   - f_5 * isf1_119[k]
                   + f_3 * pc_x[k] * isg_174[k];

        t_241[k] = f_10 * hsg_175[k]
                   + f_3 * pc_x[k] * isg_175[k];

        t_242[k] = f_10 * hsg_176[k]
                   + f_3 * pc_x[k] * isg_176[k];

        t_243[k] = f_10 * hsg_177[k]
                   + f_3 * pc_x[k] * isg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_z, pc_x, pc_z, hsh0_141, hsg_100, \
                         hsg_178, hsg_179, hsh1_141, isg_175, isg_178, \
                         isg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_10 * hsg_178[k]
                   + f_3 * pc_x[k] * isg_178[k];

        t_245[k] = f_10 * hsg_179[k]
                   + f_3 * pc_x[k] * isg_179[k];

        t_246[k] = pa_z[k] * hsh0_141[k]
                   - f_8 * pc_z[k] * hsh1_141[k];

        t_247[k] = f_9 * hsg_100[k]
                   + f_3 * pc_z[k] * isg_175[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, hsg_117, hsg_118, hsg_119, isf0_118, \
                         isf0_119, isf1_118, isf1_119, isg_177, isg_178, \
                         isg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * hsg_117[k]
                   + f_6 * isf0_118[k]
                   - f_7 * isf1_118[k]
                   + f_3 * pc_y[k] * isg_177[k];

        t_249[k] = f_11 * hsg_118[k]
                   + f_4 * isf0_119[k]
                   - f_5 * isf1_119[k]
                   + f_3 * pc_y[k] * isg_178[k];

        t_250[k] = f_11 * hsg_119[k]
                   + f_3 * pc_y[k] * isg_179[k];
    }
}

static auto
compute_prim_ish_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsh0,
                                                          const size_t hsg, const size_t hsh1,
                                                          const size_t isf0, const size_t isf1,
                                                          const size_t isg, const size_t ncols,
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
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.0 / q;

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
    auto *t_369 = buffer.data(target + 369);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsh0_189 = buffer.data(hsh0 + 189);
    const auto *hsh0_192 = buffer.data(hsh0 + 192);
    const auto *hsh0_194 = buffer.data(hsh0 + 194);
    const auto *hsh0_195 = buffer.data(hsh0 + 195);
    const auto *hsh0_198 = buffer.data(hsh0 + 198);
    const auto *hsh0_209 = buffer.data(hsh0 + 209);
    const auto *hsh0_210 = buffer.data(hsh0 + 210);
    const auto *hsh0_213 = buffer.data(hsh0 + 213);
    const auto *hsh0_216 = buffer.data(hsh0 + 216);
    const auto *hsh0_315 = buffer.data(hsh0 + 315);
    const auto *hsh0_318 = buffer.data(hsh0 + 318);
    const auto *hsh0_321 = buffer.data(hsh0 + 321);
    const auto *hsh0_330 = buffer.data(hsh0 + 330);
    const auto *hsh0_332 = buffer.data(hsh0 + 332);
    const auto *hsh0_333 = buffer.data(hsh0 + 333);
    const auto *hsh0_335 = buffer.data(hsh0 + 335);
    const auto *hsh0_341 = buffer.data(hsh0 + 341);
    const auto *hsh0_345 = buffer.data(hsh0 + 345);
    const auto *hsh0_351 = buffer.data(hsh0 + 351);
    const auto *hsh0_353 = buffer.data(hsh0 + 353);
    const auto *hsh0_354 = buffer.data(hsh0 + 354);
    const auto *hsh0_356 = buffer.data(hsh0 + 356);
    const auto *hsh0_357 = buffer.data(hsh0 + 357);
    const auto *hsh0_360 = buffer.data(hsh0 + 360);
    const auto *hsh0_362 = buffer.data(hsh0 + 362);
    const auto *hsh0_363 = buffer.data(hsh0 + 363);
    const auto *hsh0_366 = buffer.data(hsh0 + 366);

    const auto *hsg_104 = buffer.data(hsg + 104);
    const auto *hsg_105 = buffer.data(hsg + 105);
    const auto *hsg_108 = buffer.data(hsg + 108);
    const auto *hsg_115 = buffer.data(hsg + 115);
    const auto *hsg_119 = buffer.data(hsg + 119);
    const auto *hsg_120 = buffer.data(hsg + 120);
    const auto *hsg_122 = buffer.data(hsg + 122);
    const auto *hsg_123 = buffer.data(hsg + 123);
    const auto *hsg_125 = buffer.data(hsg + 125);
    const auto *hsg_130 = buffer.data(hsg + 130);
    const auto *hsg_132 = buffer.data(hsg + 132);
    const auto *hsg_133 = buffer.data(hsg + 133);
    const auto *hsg_134 = buffer.data(hsg + 134);
    const auto *hsg_135 = buffer.data(hsg + 135);
    const auto *hsg_136 = buffer.data(hsg + 136);
    const auto *hsg_137 = buffer.data(hsg + 137);
    const auto *hsg_138 = buffer.data(hsg + 138);
    const auto *hsg_140 = buffer.data(hsg + 140);
    const auto *hsg_145 = buffer.data(hsg + 145);
    const auto *hsg_147 = buffer.data(hsg + 147);
    const auto *hsg_148 = buffer.data(hsg + 148);
    const auto *hsg_149 = buffer.data(hsg + 149);
    const auto *hsg_150 = buffer.data(hsg + 150);
    const auto *hsg_153 = buffer.data(hsg + 153);
    const auto *hsg_155 = buffer.data(hsg + 155);
    const auto *hsg_160 = buffer.data(hsg + 160);
    const auto *hsg_164 = buffer.data(hsg + 164);
    const auto *hsg_165 = buffer.data(hsg + 165);
    const auto *hsg_167 = buffer.data(hsg + 167);
    const auto *hsg_168 = buffer.data(hsg + 168);
    const auto *hsg_170 = buffer.data(hsg + 170);
    const auto *hsg_179 = buffer.data(hsg + 179);
    const auto *hsg_180 = buffer.data(hsg + 180);
    const auto *hsg_182 = buffer.data(hsg + 182);
    const auto *hsg_183 = buffer.data(hsg + 183);
    const auto *hsg_185 = buffer.data(hsg + 185);
    const auto *hsg_186 = buffer.data(hsg + 186);
    const auto *hsg_189 = buffer.data(hsg + 189);
    const auto *hsg_190 = buffer.data(hsg + 190);
    const auto *hsg_191 = buffer.data(hsg + 191);
    const auto *hsg_192 = buffer.data(hsg + 192);
    const auto *hsg_193 = buffer.data(hsg + 193);
    const auto *hsg_194 = buffer.data(hsg + 194);
    const auto *hsg_205 = buffer.data(hsg + 205);
    const auto *hsg_206 = buffer.data(hsg + 206);
    const auto *hsg_207 = buffer.data(hsg + 207);
    const auto *hsg_208 = buffer.data(hsg + 208);
    const auto *hsg_209 = buffer.data(hsg + 209);
    const auto *hsg_210 = buffer.data(hsg + 210);
    const auto *hsg_215 = buffer.data(hsg + 215);
    const auto *hsg_219 = buffer.data(hsg + 219);
    const auto *hsg_220 = buffer.data(hsg + 220);
    const auto *hsg_221 = buffer.data(hsg + 221);
    const auto *hsg_222 = buffer.data(hsg + 222);
    const auto *hsg_224 = buffer.data(hsg + 224);
    const auto *hsg_225 = buffer.data(hsg + 225);
    const auto *hsg_228 = buffer.data(hsg + 228);
    const auto *hsg_231 = buffer.data(hsg + 231);
    const auto *hsg_235 = buffer.data(hsg + 235);
    const auto *hsg_237 = buffer.data(hsg + 237);
    const auto *hsg_238 = buffer.data(hsg + 238);
    const auto *hsg_239 = buffer.data(hsg + 239);
    const auto *hsg_245 = buffer.data(hsg + 245);
    const auto *hsg_249 = buffer.data(hsg + 249);
    const auto *hsg_250 = buffer.data(hsg + 250);
    const auto *hsg_251 = buffer.data(hsg + 251);
    const auto *hsg_252 = buffer.data(hsg + 252);
    const auto *hsg_253 = buffer.data(hsg + 253);
    const auto *hsg_254 = buffer.data(hsg + 254);
    const auto *hsg_255 = buffer.data(hsg + 255);
    const auto *hsg_258 = buffer.data(hsg + 258);
    const auto *hsg_260 = buffer.data(hsg + 260);
    const auto *hsg_261 = buffer.data(hsg + 261);
    const auto *hsg_264 = buffer.data(hsg + 264);
    const auto *hsg_265 = buffer.data(hsg + 265);
    const auto *hsg_266 = buffer.data(hsg + 266);
    const auto *hsg_267 = buffer.data(hsg + 267);

    const auto *hsh1_189 = buffer.data(hsh1 + 189);
    const auto *hsh1_192 = buffer.data(hsh1 + 192);
    const auto *hsh1_194 = buffer.data(hsh1 + 194);
    const auto *hsh1_195 = buffer.data(hsh1 + 195);
    const auto *hsh1_198 = buffer.data(hsh1 + 198);
    const auto *hsh1_209 = buffer.data(hsh1 + 209);
    const auto *hsh1_210 = buffer.data(hsh1 + 210);
    const auto *hsh1_213 = buffer.data(hsh1 + 213);
    const auto *hsh1_216 = buffer.data(hsh1 + 216);
    const auto *hsh1_315 = buffer.data(hsh1 + 315);
    const auto *hsh1_318 = buffer.data(hsh1 + 318);
    const auto *hsh1_321 = buffer.data(hsh1 + 321);
    const auto *hsh1_330 = buffer.data(hsh1 + 330);
    const auto *hsh1_332 = buffer.data(hsh1 + 332);
    const auto *hsh1_333 = buffer.data(hsh1 + 333);
    const auto *hsh1_335 = buffer.data(hsh1 + 335);
    const auto *hsh1_341 = buffer.data(hsh1 + 341);
    const auto *hsh1_345 = buffer.data(hsh1 + 345);
    const auto *hsh1_351 = buffer.data(hsh1 + 351);
    const auto *hsh1_353 = buffer.data(hsh1 + 353);
    const auto *hsh1_354 = buffer.data(hsh1 + 354);
    const auto *hsh1_356 = buffer.data(hsh1 + 356);
    const auto *hsh1_357 = buffer.data(hsh1 + 357);
    const auto *hsh1_360 = buffer.data(hsh1 + 360);
    const auto *hsh1_362 = buffer.data(hsh1 + 362);
    const auto *hsh1_363 = buffer.data(hsh1 + 363);
    const auto *hsh1_366 = buffer.data(hsh1 + 366);

    const auto *isf0_119 = buffer.data(isf0 + 119);
    const auto *isf0_120 = buffer.data(isf0 + 120);
    const auto *isf0_123 = buffer.data(isf0 + 123);
    const auto *isf0_125 = buffer.data(isf0 + 125);
    const auto *isf0_126 = buffer.data(isf0 + 126);
    const auto *isf0_128 = buffer.data(isf0 + 128);
    const auto *isf0_129 = buffer.data(isf0 + 129);
    const auto *isf0_136 = buffer.data(isf0 + 136);
    const auto *isf0_138 = buffer.data(isf0 + 138);
    const auto *isf0_139 = buffer.data(isf0 + 139);
    const auto *isf0_140 = buffer.data(isf0 + 140);
    const auto *isf0_141 = buffer.data(isf0 + 141);
    const auto *isf0_142 = buffer.data(isf0 + 142);
    const auto *isf0_145 = buffer.data(isf0 + 145);
    const auto *isf0_146 = buffer.data(isf0 + 146);
    const auto *isf0_147 = buffer.data(isf0 + 147);
    const auto *isf0_148 = buffer.data(isf0 + 148);
    const auto *isf0_149 = buffer.data(isf0 + 149);
    const auto *isf0_150 = buffer.data(isf0 + 150);
    const auto *isf0_152 = buffer.data(isf0 + 152);

    const auto *isf1_119 = buffer.data(isf1 + 119);
    const auto *isf1_120 = buffer.data(isf1 + 120);
    const auto *isf1_123 = buffer.data(isf1 + 123);
    const auto *isf1_125 = buffer.data(isf1 + 125);
    const auto *isf1_126 = buffer.data(isf1 + 126);
    const auto *isf1_128 = buffer.data(isf1 + 128);
    const auto *isf1_129 = buffer.data(isf1 + 129);
    const auto *isf1_136 = buffer.data(isf1 + 136);
    const auto *isf1_138 = buffer.data(isf1 + 138);
    const auto *isf1_139 = buffer.data(isf1 + 139);
    const auto *isf1_140 = buffer.data(isf1 + 140);
    const auto *isf1_141 = buffer.data(isf1 + 141);
    const auto *isf1_142 = buffer.data(isf1 + 142);
    const auto *isf1_145 = buffer.data(isf1 + 145);
    const auto *isf1_146 = buffer.data(isf1 + 146);
    const auto *isf1_147 = buffer.data(isf1 + 147);
    const auto *isf1_148 = buffer.data(isf1 + 148);
    const auto *isf1_149 = buffer.data(isf1 + 149);
    const auto *isf1_150 = buffer.data(isf1 + 150);
    const auto *isf1_152 = buffer.data(isf1 + 152);

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
    const auto *isg_195 = buffer.data(isg + 195);
    const auto *isg_197 = buffer.data(isg + 197);
    const auto *isg_198 = buffer.data(isg + 198);
    const auto *isg_200 = buffer.data(isg + 200);
    const auto *isg_205 = buffer.data(isg + 205);
    const auto *isg_206 = buffer.data(isg + 206);
    const auto *isg_207 = buffer.data(isg + 207);
    const auto *isg_208 = buffer.data(isg + 208);
    const auto *isg_209 = buffer.data(isg + 209);
    const auto *isg_210 = buffer.data(isg + 210);
    const auto *isg_211 = buffer.data(isg + 211);
    const auto *isg_212 = buffer.data(isg + 212);
    const auto *isg_213 = buffer.data(isg + 213);
    const auto *isg_214 = buffer.data(isg + 214);
    const auto *isg_215 = buffer.data(isg + 215);
    const auto *isg_219 = buffer.data(isg + 219);
    const auto *isg_220 = buffer.data(isg + 220);
    const auto *isg_221 = buffer.data(isg + 221);
    const auto *isg_222 = buffer.data(isg + 222);
    const auto *isg_223 = buffer.data(isg + 223);
    const auto *isg_224 = buffer.data(isg + 224);
    const auto *isg_225 = buffer.data(isg + 225);
    const auto *isg_226 = buffer.data(isg + 226);
    const auto *isg_227 = buffer.data(isg + 227);
    const auto *isg_228 = buffer.data(isg + 228);
    const auto *isg_230 = buffer.data(isg + 230);
    const auto *isg_231 = buffer.data(isg + 231);
    const auto *isg_235 = buffer.data(isg + 235);
    const auto *isg_237 = buffer.data(isg + 237);
    const auto *isg_238 = buffer.data(isg + 238);
    const auto *isg_239 = buffer.data(isg + 239);
    const auto *isg_240 = buffer.data(isg + 240);
    const auto *isg_242 = buffer.data(isg + 242);
    const auto *isg_243 = buffer.data(isg + 243);
    const auto *isg_245 = buffer.data(isg + 245);
    const auto *isg_250 = buffer.data(isg + 250);
    const auto *isg_251 = buffer.data(isg + 251);
    const auto *isg_252 = buffer.data(isg + 252);
    const auto *isg_253 = buffer.data(isg + 253);
    const auto *isg_254 = buffer.data(isg + 254);
    const auto *isg_255 = buffer.data(isg + 255);
    const auto *isg_257 = buffer.data(isg + 257);
    const auto *isg_258 = buffer.data(isg + 258);
    const auto *isg_260 = buffer.data(isg + 260);
    const auto *isg_265 = buffer.data(isg + 265);
    const auto *isg_266 = buffer.data(isg + 266);
    const auto *isg_267 = buffer.data(isg + 267);

#pragma omp simd aligned(t_251, t_252, t_253, pc_x, pc_y, pc_z, hsg_104, hsg_120, hsg_180, \
                         isf0_119, isf0_120, isf1_119, isf1_120, isg_179, \
                         isg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_9 * hsg_104[k]
                   + f_1 * isf0_119[k]
                   - f_2 * isf1_119[k]
                   + f_3 * pc_z[k] * isg_179[k];

        t_252[k] = f_10 * hsg_180[k]
                   + f_1 * isf0_120[k]
                   - f_2 * isf1_120[k]
                   + f_3 * pc_x[k] * isg_180[k];

        t_253[k] = f_10 * hsg_120[k]
                   + f_3 * pc_y[k] * isg_180[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, hsg_105, hsg_122, hsg_183, \
                         isf0_123, isf1_123, isg_180, isg_182, \
                         isg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * hsg_105[k]
                   + f_3 * pc_z[k] * isg_180[k];

        t_255[k] = f_10 * hsg_183[k]
                   + f_6 * isf0_123[k]
                   - f_7 * isf1_123[k]
                   + f_3 * pc_x[k] * isg_183[k];

        t_256[k] = f_10 * hsg_122[k]
                   + f_3 * pc_y[k] * isg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pc_x, pc_z, hsg_108, hsg_185, hsg_186, isf0_125, \
                         isf0_126, isf1_125, isf1_126, isg_183, isg_185, \
                         isg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_10 * hsg_185[k]
                   + f_6 * isf0_125[k]
                   - f_7 * isf1_125[k]
                   + f_3 * pc_x[k] * isg_185[k];

        t_258[k] = f_10 * hsg_186[k]
                   + f_4 * isf0_126[k]
                   - f_5 * isf1_126[k]
                   + f_3 * pc_x[k] * isg_186[k];

        t_259[k] = f_10 * hsg_108[k]
                   + f_3 * pc_z[k] * isg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, hsg_125, hsg_189, hsg_190, \
                         hsg_191, isf0_129, isf1_129, isg_185, isg_189, isg_190, \
                         isg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_10 * hsg_125[k]
                   + f_3 * pc_y[k] * isg_185[k];

        t_261[k] = f_10 * hsg_189[k]
                   + f_4 * isf0_129[k]
                   - f_5 * isf1_129[k]
                   + f_3 * pc_x[k] * isg_189[k];

        t_262[k] = f_10 * hsg_190[k]
                   + f_3 * pc_x[k] * isg_190[k];

        t_263[k] = f_10 * hsg_191[k]
                   + f_3 * pc_x[k] * isg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, hsg_130, hsg_192, hsg_193, \
                         hsg_194, isf0_126, isf1_126, isg_190, isg_192, isg_193, \
                         isg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * hsg_192[k]
                   + f_3 * pc_x[k] * isg_192[k];

        t_265[k] = f_10 * hsg_193[k]
                   + f_3 * pc_x[k] * isg_193[k];

        t_266[k] = f_10 * hsg_194[k]
                   + f_3 * pc_x[k] * isg_194[k];

        t_267[k] = f_10 * hsg_130[k]
                   + f_1 * isf0_126[k]
                   - f_2 * isf1_126[k]
                   + f_3 * pc_y[k] * isg_190[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_y, pc_z, hsg_115, hsg_132, hsg_133, isf0_128, \
                         isf0_129, isf1_128, isf1_129, isg_190, isg_192, \
                         isg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_10 * hsg_115[k]
                   + f_3 * pc_z[k] * isg_190[k];

        t_269[k] = f_10 * hsg_132[k]
                   + f_6 * isf0_128[k]
                   - f_7 * isf1_128[k]
                   + f_3 * pc_y[k] * isg_192[k];

        t_270[k] = f_10 * hsg_133[k]
                   + f_4 * isf0_129[k]
                   - f_5 * isf1_129[k]
                   + f_3 * pc_y[k] * isg_193[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pa_y, pc_y, pc_z, hsh0_189, hsg_119, \
                         hsg_134, hsg_135, hsh1_189, isf0_129, isf1_129, isg_194, \
                         isg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_10 * hsg_134[k]
                   + f_3 * pc_y[k] * isg_194[k];

        t_272[k] = f_10 * hsg_119[k]
                   + f_1 * isf0_129[k]
                   - f_2 * isf1_129[k]
                   + f_3 * pc_z[k] * isg_194[k];

        t_273[k] = pa_y[k] * hsh0_189[k]
                   - f_8 * pc_y[k] * hsh1_189[k];

        t_274[k] = f_9 * hsg_135[k]
                   + f_3 * pc_y[k] * isg_195[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pa_y, pc_y, pc_z, hsh0_192, hsh0_194, \
                         hsg_120, hsg_136, hsg_137, hsh1_192, hsh1_194, isg_195, \
                         isg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_11 * hsg_120[k]
                   + f_3 * pc_z[k] * isg_195[k];

        t_276[k] = pa_y[k] * hsh0_192[k]
                   + f_10 * hsg_136[k]
                   - f_8 * pc_y[k] * hsh1_192[k];

        t_277[k] = f_9 * hsg_137[k]
                   + f_3 * pc_y[k] * isg_197[k];

        t_278[k] = pa_y[k] * hsh0_194[k]
                   - f_8 * pc_y[k] * hsh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_y, pc_y, pc_z, hsh0_195, hsh0_198, \
                         hsg_123, hsg_138, hsg_140, hsh1_195, hsh1_198, isg_198, \
                         isg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pa_y[k] * hsh0_195[k]
                   + f_11 * hsg_138[k]
                   - f_8 * pc_y[k] * hsh1_195[k];

        t_280[k] = f_11 * hsg_123[k]
                   + f_3 * pc_z[k] * isg_198[k];

        t_281[k] = f_9 * hsg_140[k]
                   + f_3 * pc_y[k] * isg_200[k];

        t_282[k] = pa_y[k] * hsh0_198[k]
                   - f_8 * pc_y[k] * hsh1_198[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, t_287, pc_x, hsg_205, hsg_206, hsg_207, \
                         hsg_208, hsg_209, isg_205, isg_206, isg_207, isg_208, \
                         isg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_10 * hsg_205[k]
                   + f_3 * pc_x[k] * isg_205[k];

        t_284[k] = f_10 * hsg_206[k]
                   + f_3 * pc_x[k] * isg_206[k];

        t_285[k] = f_10 * hsg_207[k]
                   + f_3 * pc_x[k] * isg_207[k];

        t_286[k] = f_10 * hsg_208[k]
                   + f_3 * pc_x[k] * isg_208[k];

        t_287[k] = f_10 * hsg_209[k]
                   + f_3 * pc_x[k] * isg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pc_y, pc_z, hsg_130, hsg_145, hsg_147, isf0_136, \
                         isf0_138, isf1_136, isf1_138, isg_205, \
                         isg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_9 * hsg_145[k]
                   + f_1 * isf0_136[k]
                   - f_2 * isf1_136[k]
                   + f_3 * pc_y[k] * isg_205[k];

        t_289[k] = f_11 * hsg_130[k]
                   + f_3 * pc_z[k] * isg_205[k];

        t_290[k] = f_9 * hsg_147[k]
                   + f_6 * isf0_138[k]
                   - f_7 * isf1_138[k]
                   + f_3 * pc_y[k] * isg_207[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pc_y, hsh0_209, hsg_148, hsg_149, \
                         hsh1_209, isf0_139, isf1_139, isg_208, \
                         isg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_9 * hsg_148[k]
                   + f_4 * isf0_139[k]
                   - f_5 * isf1_139[k]
                   + f_3 * pc_y[k] * isg_208[k];

        t_292[k] = f_9 * hsg_149[k]
                   + f_3 * pc_y[k] * isg_209[k];

        t_293[k] = pa_y[k] * hsh0_209[k]
                   - f_8 * pc_y[k] * hsh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, hsg_135, \
                         hsg_210, isf0_140, isf1_140, isg_210, isg_211, \
                         isg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_10 * hsg_210[k]
                   + f_1 * isf0_140[k]
                   - f_2 * isf1_140[k]
                   + f_3 * pc_x[k] * isg_210[k];

        t_295[k] = f_3 * pc_y[k] * isg_210[k];

        t_296[k] = f_15 * hsg_135[k]
                   + f_3 * pc_z[k] * isg_210[k];

        t_297[k] = f_4 * isf0_140[k]
                   - f_5 * isf1_140[k]
                   + f_3 * pc_y[k] * isg_211[k];

        t_298[k] = f_3 * pc_y[k] * isg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, hsg_215, isf0_141, isf0_142, \
                         isf0_145, isf1_141, isf1_142, isf1_145, isg_213, isg_214, \
                         isg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_10 * hsg_215[k]
                   + f_6 * isf0_145[k]
                   - f_7 * isf1_145[k]
                   + f_3 * pc_x[k] * isg_215[k];

        t_300[k] = f_6 * isf0_141[k]
                   - f_7 * isf1_141[k]
                   + f_3 * pc_y[k] * isg_213[k];

        t_301[k] = f_4 * isf0_142[k]
                   - f_5 * isf1_142[k]
                   + f_3 * pc_y[k] * isg_214[k];

        t_302[k] = f_3 * pc_y[k] * isg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pc_x, hsg_219, hsg_220, hsg_221, hsg_222, \
                         isf0_149, isf1_149, isg_219, isg_220, isg_221, \
                         isg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_10 * hsg_219[k]
                   + f_4 * isf0_149[k]
                   - f_5 * isf1_149[k]
                   + f_3 * pc_x[k] * isg_219[k];

        t_304[k] = f_10 * hsg_220[k]
                   + f_3 * pc_x[k] * isg_220[k];

        t_305[k] = f_10 * hsg_221[k]
                   + f_3 * pc_x[k] * isg_221[k];

        t_306[k] = f_10 * hsg_222[k]
                   + f_3 * pc_x[k] * isg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, t_310, pc_x, pc_y, hsg_224, isf0_146, isf0_147, \
                         isf1_146, isf1_147, isg_219, isg_220, isg_221, \
                         isg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_3 * pc_y[k] * isg_219[k];

        t_308[k] = f_10 * hsg_224[k]
                   + f_3 * pc_x[k] * isg_224[k];

        t_309[k] = f_1 * isf0_146[k]
                   - f_2 * isf1_146[k]
                   + f_3 * pc_y[k] * isg_220[k];

        t_310[k] = f_13 * isf0_147[k]
                   - f_14 * isf1_147[k]
                   + f_3 * pc_y[k] * isg_221[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, hsg_149, isf0_148, isf0_149, \
                         isf1_148, isf1_149, isg_222, isg_223, \
                         isg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_6 * isf0_148[k]
                   - f_7 * isf1_148[k]
                   + f_3 * pc_y[k] * isg_222[k];

        t_312[k] = f_4 * isf0_149[k]
                   - f_5 * isf1_149[k]
                   + f_3 * pc_y[k] * isg_223[k];

        t_313[k] = f_3 * pc_y[k] * isg_224[k];

        t_314[k] = f_15 * hsg_149[k]
                   + f_1 * isf0_149[k]
                   - f_2 * isf1_149[k]
                   + f_3 * pc_z[k] * isg_224[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pa_x, pc_x, pc_y, pc_z, hsh0_315, \
                         hsh0_318, hsg_150, hsg_225, hsg_228, hsh1_315, hsh1_318, \
                         isg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pa_x[k] * hsh0_315[k]
                   + f_12 * hsg_225[k]
                   - f_8 * pc_x[k] * hsh1_315[k];

        t_316[k] = f_12 * hsg_150[k]
                   + f_3 * pc_y[k] * isg_225[k];

        t_317[k] = f_3 * pc_z[k] * isg_225[k];

        t_318[k] = pa_x[k] * hsh0_318[k]
                   + f_11 * hsg_228[k]
                   - f_8 * pc_x[k] * hsh1_318[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, t_322, pa_x, pc_x, pc_z, hsh0_321, hsg_231, \
                         hsh1_321, isf0_150, isf1_150, isg_226, isg_227, \
                         isg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_3 * pc_z[k] * isg_226[k];

        t_320[k] = f_4 * isf0_150[k]
                   - f_5 * isf1_150[k]
                   + f_3 * pc_z[k] * isg_227[k];

        t_321[k] = pa_x[k] * hsh0_321[k]
                   + f_10 * hsg_231[k]
                   - f_8 * pc_x[k] * hsh1_321[k];

        t_322[k] = f_3 * pc_z[k] * isg_228[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, t_326, pc_x, pc_y, pc_z, hsg_155, hsg_235, \
                         isf0_152, isf1_152, isg_230, isg_231, \
                         isg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_12 * hsg_155[k]
                   + f_3 * pc_y[k] * isg_230[k];

        t_324[k] = f_6 * isf0_152[k]
                   - f_7 * isf1_152[k]
                   + f_3 * pc_z[k] * isg_230[k];

        t_325[k] = f_9 * hsg_235[k]
                   + f_3 * pc_x[k] * isg_235[k];

        t_326[k] = f_3 * pc_z[k] * isg_231[k];
    }

#pragma omp simd aligned(t_327, t_328, t_329, t_330, pa_x, pc_x, hsh0_330, hsg_237, hsg_238, \
                         hsg_239, hsh1_330, isg_237, isg_238, isg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_327[k] = f_9 * hsg_237[k]
                   + f_3 * pc_x[k] * isg_237[k];

        t_328[k] = f_9 * hsg_238[k]
                   + f_3 * pc_x[k] * isg_238[k];

        t_329[k] = f_9 * hsg_239[k]
                   + f_3 * pc_x[k] * isg_239[k];

        t_330[k] = pa_x[k] * hsh0_330[k]
                   - f_8 * pc_x[k] * hsh1_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pa_x, pc_x, pc_y, pc_z, hsh0_332, \
                         hsh0_333, hsg_164, hsh1_332, hsh1_333, isg_235, \
                         isg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * isg_235[k];

        t_332[k] = pa_x[k] * hsh0_332[k]
                   - f_8 * pc_x[k] * hsh1_332[k];

        t_333[k] = pa_x[k] * hsh0_333[k]
                   - f_8 * pc_x[k] * hsh1_333[k];

        t_334[k] = f_12 * hsg_164[k]
                   + f_3 * pc_y[k] * isg_239[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_x, pa_z, pc_x, pc_y, pc_z, hsh0_210, \
                         hsh0_335, hsg_150, hsg_165, hsh1_210, hsh1_335, \
                         isg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pa_x[k] * hsh0_335[k]
                   - f_8 * pc_x[k] * hsh1_335[k];

        t_336[k] = pa_z[k] * hsh0_210[k]
                   - f_8 * pc_z[k] * hsh1_210[k];

        t_337[k] = f_15 * hsg_165[k]
                   + f_3 * pc_y[k] * isg_240[k];

        t_338[k] = f_9 * hsg_150[k]
                   + f_3 * pc_z[k] * isg_240[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pa_x, pa_z, pc_x, pc_y, pc_z, hsh0_213, \
                         hsh0_341, hsg_167, hsg_245, hsh1_213, hsh1_341, \
                         isg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = pa_z[k] * hsh0_213[k]
                   - f_8 * pc_z[k] * hsh1_213[k];

        t_340[k] = f_15 * hsg_167[k]
                   + f_3 * pc_y[k] * isg_242[k];

        t_341[k] = pa_x[k] * hsh0_341[k]
                   + f_11 * hsg_245[k]
                   - f_8 * pc_x[k] * hsh1_341[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pa_z, pc_y, pc_z, hsh0_216, hsg_153, hsg_170, \
                         hsh1_216, isg_243, isg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = pa_z[k] * hsh0_216[k]
                   - f_8 * pc_z[k] * hsh1_216[k];

        t_343[k] = f_9 * hsg_153[k]
                   + f_3 * pc_z[k] * isg_243[k];

        t_344[k] = f_15 * hsg_170[k]
                   + f_3 * pc_y[k] * isg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pa_x, pc_x, hsh0_345, hsg_249, hsg_250, \
                         hsg_251, hsg_252, hsh1_345, isg_250, isg_251, \
                         isg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pa_x[k] * hsh0_345[k]
                   + f_10 * hsg_249[k]
                   - f_8 * pc_x[k] * hsh1_345[k];

        t_346[k] = f_9 * hsg_250[k]
                   + f_3 * pc_x[k] * isg_250[k];

        t_347[k] = f_9 * hsg_251[k]
                   + f_3 * pc_x[k] * isg_251[k];

        t_348[k] = f_9 * hsg_252[k]
                   + f_3 * pc_x[k] * isg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_x, pc_x, pc_z, hsh0_351, hsg_160, \
                         hsg_253, hsg_254, hsh1_351, isg_250, isg_253, \
                         isg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_9 * hsg_253[k]
                   + f_3 * pc_x[k] * isg_253[k];

        t_350[k] = f_9 * hsg_254[k]
                   + f_3 * pc_x[k] * isg_254[k];

        t_351[k] = pa_x[k] * hsh0_351[k]
                   - f_8 * pc_x[k] * hsh1_351[k];

        t_352[k] = f_9 * hsg_160[k]
                   + f_3 * pc_z[k] * isg_250[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pa_x, pc_x, pc_y, hsh0_353, hsh0_354, \
                         hsh0_356, hsg_179, hsh1_353, hsh1_354, hsh1_356, \
                         isg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = pa_x[k] * hsh0_353[k]
                   - f_8 * pc_x[k] * hsh1_353[k];

        t_354[k] = pa_x[k] * hsh0_354[k]
                   - f_8 * pc_x[k] * hsh1_354[k];

        t_355[k] = f_15 * hsg_179[k]
                   + f_3 * pc_y[k] * isg_254[k];

        t_356[k] = pa_x[k] * hsh0_356[k]
                   - f_8 * pc_x[k] * hsh1_356[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pa_x, pc_x, pc_y, pc_z, hsh0_357, hsg_165, \
                         hsg_180, hsg_255, hsh1_357, isg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_x[k] * hsh0_357[k]
                   + f_12 * hsg_255[k]
                   - f_8 * pc_x[k] * hsh1_357[k];

        t_358[k] = f_11 * hsg_180[k]
                   + f_3 * pc_y[k] * isg_255[k];

        t_359[k] = f_10 * hsg_165[k]
                   + f_3 * pc_z[k] * isg_255[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pa_x, pc_x, pc_y, hsh0_360, hsh0_362, hsg_182, \
                         hsg_258, hsg_260, hsh1_360, hsh1_362, \
                         isg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pa_x[k] * hsh0_360[k]
                   + f_11 * hsg_258[k]
                   - f_8 * pc_x[k] * hsh1_360[k];

        t_361[k] = f_11 * hsg_182[k]
                   + f_3 * pc_y[k] * isg_257[k];

        t_362[k] = pa_x[k] * hsh0_362[k]
                   + f_11 * hsg_260[k]
                   - f_8 * pc_x[k] * hsh1_362[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_x, pc_x, pc_y, pc_z, hsh0_363, hsg_168, \
                         hsg_185, hsg_261, hsh1_363, isg_258, isg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pa_x[k] * hsh0_363[k]
                   + f_10 * hsg_261[k]
                   - f_8 * pc_x[k] * hsh1_363[k];

        t_364[k] = f_10 * hsg_168[k]
                   + f_3 * pc_z[k] * isg_258[k];

        t_365[k] = f_11 * hsg_185[k]
                   + f_3 * pc_y[k] * isg_260[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_x, pc_x, hsh0_366, hsg_264, hsg_265, \
                         hsg_266, hsg_267, hsh1_366, isg_265, isg_266, \
                         isg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = pa_x[k] * hsh0_366[k]
                   + f_10 * hsg_264[k]
                   - f_8 * pc_x[k] * hsh1_366[k];

        t_367[k] = f_9 * hsg_265[k]
                   + f_3 * pc_x[k] * isg_265[k];

        t_368[k] = f_9 * hsg_266[k]
                   + f_3 * pc_x[k] * isg_266[k];

        t_369[k] = f_9 * hsg_267[k]
                   + f_3 * pc_x[k] * isg_267[k];
    }
}

static auto
compute_prim_ish_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsh0,
                                                          const size_t hsg, const size_t hsh1,
                                                          const size_t isf0, const size_t isf1,
                                                          const size_t isg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.0 / q;

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
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsh0_294 = buffer.data(hsh0 + 294);
    const auto *hsh0_299 = buffer.data(hsh0 + 299);
    const auto *hsh0_303 = buffer.data(hsh0 + 303);
    const auto *hsh0_315 = buffer.data(hsh0 + 315);
    const auto *hsh0_316 = buffer.data(hsh0 + 316);
    const auto *hsh0_318 = buffer.data(hsh0 + 318);
    const auto *hsh0_321 = buffer.data(hsh0 + 321);
    const auto *hsh0_330 = buffer.data(hsh0 + 330);
    const auto *hsh0_332 = buffer.data(hsh0 + 332);
    const auto *hsh0_333 = buffer.data(hsh0 + 333);
    const auto *hsh0_372 = buffer.data(hsh0 + 372);
    const auto *hsh0_374 = buffer.data(hsh0 + 374);
    const auto *hsh0_375 = buffer.data(hsh0 + 375);
    const auto *hsh0_377 = buffer.data(hsh0 + 377);
    const auto *hsh0_378 = buffer.data(hsh0 + 378);
    const auto *hsh0_381 = buffer.data(hsh0 + 381);
    const auto *hsh0_383 = buffer.data(hsh0 + 383);
    const auto *hsh0_384 = buffer.data(hsh0 + 384);
    const auto *hsh0_387 = buffer.data(hsh0 + 387);
    const auto *hsh0_393 = buffer.data(hsh0 + 393);
    const auto *hsh0_395 = buffer.data(hsh0 + 395);
    const auto *hsh0_396 = buffer.data(hsh0 + 396);
    const auto *hsh0_398 = buffer.data(hsh0 + 398);
    const auto *hsh0_402 = buffer.data(hsh0 + 402);
    const auto *hsh0_405 = buffer.data(hsh0 + 405);
    const auto *hsh0_414 = buffer.data(hsh0 + 414);
    const auto *hsh0_416 = buffer.data(hsh0 + 416);
    const auto *hsh0_417 = buffer.data(hsh0 + 417);
    const auto *hsh0_419 = buffer.data(hsh0 + 419);
    const auto *hsh0_420 = buffer.data(hsh0 + 420);
    const auto *hsh0_425 = buffer.data(hsh0 + 425);
    const auto *hsh0_429 = buffer.data(hsh0 + 429);
    const auto *hsh0_435 = buffer.data(hsh0 + 435);
    const auto *hsh0_436 = buffer.data(hsh0 + 436);
    const auto *hsh0_437 = buffer.data(hsh0 + 437);
    const auto *hsh0_438 = buffer.data(hsh0 + 438);
    const auto *hsh0_440 = buffer.data(hsh0 + 440);

    const auto *hsg_175 = buffer.data(hsg + 175);
    const auto *hsg_180 = buffer.data(hsg + 180);
    const auto *hsg_183 = buffer.data(hsg + 183);
    const auto *hsg_190 = buffer.data(hsg + 190);
    const auto *hsg_194 = buffer.data(hsg + 194);
    const auto *hsg_195 = buffer.data(hsg + 195);
    const auto *hsg_197 = buffer.data(hsg + 197);
    const auto *hsg_198 = buffer.data(hsg + 198);
    const auto *hsg_200 = buffer.data(hsg + 200);
    const auto *hsg_205 = buffer.data(hsg + 205);
    const auto *hsg_209 = buffer.data(hsg + 209);
    const auto *hsg_210 = buffer.data(hsg + 210);
    const auto *hsg_212 = buffer.data(hsg + 212);
    const auto *hsg_215 = buffer.data(hsg + 215);
    const auto *hsg_224 = buffer.data(hsg + 224);
    const auto *hsg_235 = buffer.data(hsg + 235);
    const auto *hsg_236 = buffer.data(hsg + 236);
    const auto *hsg_237 = buffer.data(hsg + 237);
    const auto *hsg_239 = buffer.data(hsg + 239);
    const auto *hsg_254 = buffer.data(hsg + 254);
    const auto *hsg_268 = buffer.data(hsg + 268);
    const auto *hsg_269 = buffer.data(hsg + 269);
    const auto *hsg_270 = buffer.data(hsg + 270);
    const auto *hsg_273 = buffer.data(hsg + 273);
    const auto *hsg_275 = buffer.data(hsg + 275);
    const auto *hsg_276 = buffer.data(hsg + 276);
    const auto *hsg_279 = buffer.data(hsg + 279);
    const auto *hsg_280 = buffer.data(hsg + 280);
    const auto *hsg_281 = buffer.data(hsg + 281);
    const auto *hsg_282 = buffer.data(hsg + 282);
    const auto *hsg_283 = buffer.data(hsg + 283);
    const auto *hsg_284 = buffer.data(hsg + 284);
    const auto *hsg_288 = buffer.data(hsg + 288);
    const auto *hsg_291 = buffer.data(hsg + 291);
    const auto *hsg_295 = buffer.data(hsg + 295);
    const auto *hsg_296 = buffer.data(hsg + 296);
    const auto *hsg_297 = buffer.data(hsg + 297);
    const auto *hsg_298 = buffer.data(hsg + 298);
    const auto *hsg_299 = buffer.data(hsg + 299);
    const auto *hsg_300 = buffer.data(hsg + 300);
    const auto *hsg_305 = buffer.data(hsg + 305);
    const auto *hsg_309 = buffer.data(hsg + 309);
    const auto *hsg_310 = buffer.data(hsg + 310);
    const auto *hsg_311 = buffer.data(hsg + 311);
    const auto *hsg_312 = buffer.data(hsg + 312);
    const auto *hsg_314 = buffer.data(hsg + 314);

    const auto *hsh1_294 = buffer.data(hsh1 + 294);
    const auto *hsh1_299 = buffer.data(hsh1 + 299);
    const auto *hsh1_303 = buffer.data(hsh1 + 303);
    const auto *hsh1_315 = buffer.data(hsh1 + 315);
    const auto *hsh1_316 = buffer.data(hsh1 + 316);
    const auto *hsh1_318 = buffer.data(hsh1 + 318);
    const auto *hsh1_321 = buffer.data(hsh1 + 321);
    const auto *hsh1_330 = buffer.data(hsh1 + 330);
    const auto *hsh1_332 = buffer.data(hsh1 + 332);
    const auto *hsh1_333 = buffer.data(hsh1 + 333);
    const auto *hsh1_372 = buffer.data(hsh1 + 372);
    const auto *hsh1_374 = buffer.data(hsh1 + 374);
    const auto *hsh1_375 = buffer.data(hsh1 + 375);
    const auto *hsh1_377 = buffer.data(hsh1 + 377);
    const auto *hsh1_378 = buffer.data(hsh1 + 378);
    const auto *hsh1_381 = buffer.data(hsh1 + 381);
    const auto *hsh1_383 = buffer.data(hsh1 + 383);
    const auto *hsh1_384 = buffer.data(hsh1 + 384);
    const auto *hsh1_387 = buffer.data(hsh1 + 387);
    const auto *hsh1_393 = buffer.data(hsh1 + 393);
    const auto *hsh1_395 = buffer.data(hsh1 + 395);
    const auto *hsh1_396 = buffer.data(hsh1 + 396);
    const auto *hsh1_398 = buffer.data(hsh1 + 398);
    const auto *hsh1_402 = buffer.data(hsh1 + 402);
    const auto *hsh1_405 = buffer.data(hsh1 + 405);
    const auto *hsh1_414 = buffer.data(hsh1 + 414);
    const auto *hsh1_416 = buffer.data(hsh1 + 416);
    const auto *hsh1_417 = buffer.data(hsh1 + 417);
    const auto *hsh1_419 = buffer.data(hsh1 + 419);
    const auto *hsh1_420 = buffer.data(hsh1 + 420);
    const auto *hsh1_425 = buffer.data(hsh1 + 425);
    const auto *hsh1_429 = buffer.data(hsh1 + 429);
    const auto *hsh1_435 = buffer.data(hsh1 + 435);
    const auto *hsh1_436 = buffer.data(hsh1 + 436);
    const auto *hsh1_437 = buffer.data(hsh1 + 437);
    const auto *hsh1_438 = buffer.data(hsh1 + 438);
    const auto *hsh1_440 = buffer.data(hsh1 + 440);

    const auto *isf0_200 = buffer.data(isf0 + 200);
    const auto *isf0_201 = buffer.data(isf0 + 201);
    const auto *isf0_202 = buffer.data(isf0 + 202);
    const auto *isf0_210 = buffer.data(isf0 + 210);
    const auto *isf0_211 = buffer.data(isf0 + 211);
    const auto *isf0_213 = buffer.data(isf0 + 213);
    const auto *isf0_215 = buffer.data(isf0 + 215);
    const auto *isf0_216 = buffer.data(isf0 + 216);
    const auto *isf0_217 = buffer.data(isf0 + 217);
    const auto *isf0_218 = buffer.data(isf0 + 218);
    const auto *isf0_219 = buffer.data(isf0 + 219);
    const auto *isf0_222 = buffer.data(isf0 + 222);
    const auto *isf0_224 = buffer.data(isf0 + 224);
    const auto *isf0_225 = buffer.data(isf0 + 225);
    const auto *isf0_227 = buffer.data(isf0 + 227);
    const auto *isf0_228 = buffer.data(isf0 + 228);
    const auto *isf0_229 = buffer.data(isf0 + 229);
    const auto *isf0_230 = buffer.data(isf0 + 230);
    const auto *isf0_231 = buffer.data(isf0 + 231);
    const auto *isf0_232 = buffer.data(isf0 + 232);
    const auto *isf0_233 = buffer.data(isf0 + 233);
    const auto *isf0_234 = buffer.data(isf0 + 234);
    const auto *isf0_235 = buffer.data(isf0 + 235);
    const auto *isf0_236 = buffer.data(isf0 + 236);
    const auto *isf0_237 = buffer.data(isf0 + 237);
    const auto *isf0_238 = buffer.data(isf0 + 238);
    const auto *isf0_239 = buffer.data(isf0 + 239);

    const auto *isf1_200 = buffer.data(isf1 + 200);
    const auto *isf1_201 = buffer.data(isf1 + 201);
    const auto *isf1_202 = buffer.data(isf1 + 202);
    const auto *isf1_210 = buffer.data(isf1 + 210);
    const auto *isf1_211 = buffer.data(isf1 + 211);
    const auto *isf1_213 = buffer.data(isf1 + 213);
    const auto *isf1_215 = buffer.data(isf1 + 215);
    const auto *isf1_216 = buffer.data(isf1 + 216);
    const auto *isf1_217 = buffer.data(isf1 + 217);
    const auto *isf1_218 = buffer.data(isf1 + 218);
    const auto *isf1_219 = buffer.data(isf1 + 219);
    const auto *isf1_222 = buffer.data(isf1 + 222);
    const auto *isf1_224 = buffer.data(isf1 + 224);
    const auto *isf1_225 = buffer.data(isf1 + 225);
    const auto *isf1_227 = buffer.data(isf1 + 227);
    const auto *isf1_228 = buffer.data(isf1 + 228);
    const auto *isf1_229 = buffer.data(isf1 + 229);
    const auto *isf1_230 = buffer.data(isf1 + 230);
    const auto *isf1_231 = buffer.data(isf1 + 231);
    const auto *isf1_232 = buffer.data(isf1 + 232);
    const auto *isf1_233 = buffer.data(isf1 + 233);
    const auto *isf1_234 = buffer.data(isf1 + 234);
    const auto *isf1_235 = buffer.data(isf1 + 235);
    const auto *isf1_236 = buffer.data(isf1 + 236);
    const auto *isf1_237 = buffer.data(isf1 + 237);
    const auto *isf1_238 = buffer.data(isf1 + 238);
    const auto *isf1_239 = buffer.data(isf1 + 239);

    const auto *isg_265 = buffer.data(isg + 265);
    const auto *isg_268 = buffer.data(isg + 268);
    const auto *isg_269 = buffer.data(isg + 269);
    const auto *isg_270 = buffer.data(isg + 270);
    const auto *isg_272 = buffer.data(isg + 272);
    const auto *isg_273 = buffer.data(isg + 273);
    const auto *isg_275 = buffer.data(isg + 275);
    const auto *isg_280 = buffer.data(isg + 280);
    const auto *isg_281 = buffer.data(isg + 281);
    const auto *isg_282 = buffer.data(isg + 282);
    const auto *isg_283 = buffer.data(isg + 283);
    const auto *isg_284 = buffer.data(isg + 284);
    const auto *isg_285 = buffer.data(isg + 285);
    const auto *isg_287 = buffer.data(isg + 287);
    const auto *isg_288 = buffer.data(isg + 288);
    const auto *isg_290 = buffer.data(isg + 290);
    const auto *isg_295 = buffer.data(isg + 295);
    const auto *isg_296 = buffer.data(isg + 296);
    const auto *isg_297 = buffer.data(isg + 297);
    const auto *isg_298 = buffer.data(isg + 298);
    const auto *isg_299 = buffer.data(isg + 299);
    const auto *isg_300 = buffer.data(isg + 300);
    const auto *isg_301 = buffer.data(isg + 301);
    const auto *isg_302 = buffer.data(isg + 302);
    const auto *isg_303 = buffer.data(isg + 303);
    const auto *isg_304 = buffer.data(isg + 304);
    const auto *isg_305 = buffer.data(isg + 305);
    const auto *isg_309 = buffer.data(isg + 309);
    const auto *isg_310 = buffer.data(isg + 310);
    const auto *isg_311 = buffer.data(isg + 311);
    const auto *isg_312 = buffer.data(isg + 312);
    const auto *isg_314 = buffer.data(isg + 314);
    const auto *isg_315 = buffer.data(isg + 315);
    const auto *isg_316 = buffer.data(isg + 316);
    const auto *isg_318 = buffer.data(isg + 318);
    const auto *isg_320 = buffer.data(isg + 320);
    const auto *isg_321 = buffer.data(isg + 321);
    const auto *isg_323 = buffer.data(isg + 323);
    const auto *isg_324 = buffer.data(isg + 324);
    const auto *isg_325 = buffer.data(isg + 325);
    const auto *isg_326 = buffer.data(isg + 326);
    const auto *isg_327 = buffer.data(isg + 327);
    const auto *isg_328 = buffer.data(isg + 328);
    const auto *isg_329 = buffer.data(isg + 329);
    const auto *isg_332 = buffer.data(isg + 332);
    const auto *isg_334 = buffer.data(isg + 334);
    const auto *isg_335 = buffer.data(isg + 335);
    const auto *isg_337 = buffer.data(isg + 337);
    const auto *isg_338 = buffer.data(isg + 338);
    const auto *isg_339 = buffer.data(isg + 339);
    const auto *isg_340 = buffer.data(isg + 340);
    const auto *isg_341 = buffer.data(isg + 341);
    const auto *isg_342 = buffer.data(isg + 342);
    const auto *isg_343 = buffer.data(isg + 343);
    const auto *isg_344 = buffer.data(isg + 344);
    const auto *isg_345 = buffer.data(isg + 345);
    const auto *isg_346 = buffer.data(isg + 346);
    const auto *isg_347 = buffer.data(isg + 347);
    const auto *isg_348 = buffer.data(isg + 348);
    const auto *isg_349 = buffer.data(isg + 349);
    const auto *isg_350 = buffer.data(isg + 350);
    const auto *isg_351 = buffer.data(isg + 351);
    const auto *isg_352 = buffer.data(isg + 352);
    const auto *isg_353 = buffer.data(isg + 353);
    const auto *isg_354 = buffer.data(isg + 354);
    const auto *isg_355 = buffer.data(isg + 355);
    const auto *isg_356 = buffer.data(isg + 356);
    const auto *isg_357 = buffer.data(isg + 357);

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_x, pc_x, pc_z, hsh0_372, hsg_175, \
                         hsg_268, hsg_269, hsh1_372, isg_265, isg_268, \
                         isg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_9 * hsg_268[k]
                   + f_3 * pc_x[k] * isg_268[k];

        t_371[k] = f_9 * hsg_269[k]
                   + f_3 * pc_x[k] * isg_269[k];

        t_372[k] = pa_x[k] * hsh0_372[k]
                   - f_8 * pc_x[k] * hsh1_372[k];

        t_373[k] = f_10 * hsg_175[k]
                   + f_3 * pc_z[k] * isg_265[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pa_x, pc_x, pc_y, hsh0_374, hsh0_375, \
                         hsh0_377, hsg_194, hsh1_374, hsh1_375, hsh1_377, \
                         isg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_x[k] * hsh0_374[k]
                   - f_8 * pc_x[k] * hsh1_374[k];

        t_375[k] = pa_x[k] * hsh0_375[k]
                   - f_8 * pc_x[k] * hsh1_375[k];

        t_376[k] = f_11 * hsg_194[k]
                   + f_3 * pc_y[k] * isg_269[k];

        t_377[k] = pa_x[k] * hsh0_377[k]
                   - f_8 * pc_x[k] * hsh1_377[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pa_x, pc_x, pc_y, pc_z, hsh0_378, hsg_180, \
                         hsg_195, hsg_270, hsh1_378, isg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_x[k] * hsh0_378[k]
                   + f_12 * hsg_270[k]
                   - f_8 * pc_x[k] * hsh1_378[k];

        t_379[k] = f_10 * hsg_195[k]
                   + f_3 * pc_y[k] * isg_270[k];

        t_380[k] = f_11 * hsg_180[k]
                   + f_3 * pc_z[k] * isg_270[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, pa_x, pc_x, pc_y, hsh0_381, hsh0_383, hsg_197, \
                         hsg_273, hsg_275, hsh1_381, hsh1_383, \
                         isg_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = pa_x[k] * hsh0_381[k]
                   + f_11 * hsg_273[k]
                   - f_8 * pc_x[k] * hsh1_381[k];

        t_382[k] = f_10 * hsg_197[k]
                   + f_3 * pc_y[k] * isg_272[k];

        t_383[k] = pa_x[k] * hsh0_383[k]
                   + f_11 * hsg_275[k]
                   - f_8 * pc_x[k] * hsh1_383[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pa_x, pc_x, pc_y, pc_z, hsh0_384, hsg_183, \
                         hsg_200, hsg_276, hsh1_384, isg_273, isg_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = pa_x[k] * hsh0_384[k]
                   + f_10 * hsg_276[k]
                   - f_8 * pc_x[k] * hsh1_384[k];

        t_385[k] = f_11 * hsg_183[k]
                   + f_3 * pc_z[k] * isg_273[k];

        t_386[k] = f_10 * hsg_200[k]
                   + f_3 * pc_y[k] * isg_275[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_x, pc_x, hsh0_387, hsg_279, hsg_280, \
                         hsg_281, hsg_282, hsh1_387, isg_280, isg_281, \
                         isg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = pa_x[k] * hsh0_387[k]
                   + f_10 * hsg_279[k]
                   - f_8 * pc_x[k] * hsh1_387[k];

        t_388[k] = f_9 * hsg_280[k]
                   + f_3 * pc_x[k] * isg_280[k];

        t_389[k] = f_9 * hsg_281[k]
                   + f_3 * pc_x[k] * isg_281[k];

        t_390[k] = f_9 * hsg_282[k]
                   + f_3 * pc_x[k] * isg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_x, pc_x, pc_z, hsh0_393, hsg_190, \
                         hsg_283, hsg_284, hsh1_393, isg_280, isg_283, \
                         isg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_9 * hsg_283[k]
                   + f_3 * pc_x[k] * isg_283[k];

        t_392[k] = f_9 * hsg_284[k]
                   + f_3 * pc_x[k] * isg_284[k];

        t_393[k] = pa_x[k] * hsh0_393[k]
                   - f_8 * pc_x[k] * hsh1_393[k];

        t_394[k] = f_11 * hsg_190[k]
                   + f_3 * pc_z[k] * isg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_x, pc_x, pc_y, hsh0_395, hsh0_396, \
                         hsh0_398, hsg_209, hsh1_395, hsh1_396, hsh1_398, \
                         isg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_x[k] * hsh0_395[k]
                   - f_8 * pc_x[k] * hsh1_395[k];

        t_396[k] = pa_x[k] * hsh0_396[k]
                   - f_8 * pc_x[k] * hsh1_396[k];

        t_397[k] = f_10 * hsg_209[k]
                   + f_3 * pc_y[k] * isg_284[k];

        t_398[k] = pa_x[k] * hsh0_398[k]
                   - f_8 * pc_x[k] * hsh1_398[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pa_y, pc_y, pc_z, hsh0_294, hsg_195, hsg_210, \
                         hsh1_294, isg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = pa_y[k] * hsh0_294[k]
                   - f_8 * pc_y[k] * hsh1_294[k];

        t_400[k] = f_9 * hsg_210[k]
                   + f_3 * pc_y[k] * isg_285[k];

        t_401[k] = f_15 * hsg_195[k]
                   + f_3 * pc_z[k] * isg_285[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_x, pa_y, pc_x, pc_y, hsh0_299, hsh0_402, \
                         hsg_212, hsg_288, hsh1_299, hsh1_402, \
                         isg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_x[k] * hsh0_402[k]
                   + f_11 * hsg_288[k]
                   - f_8 * pc_x[k] * hsh1_402[k];

        t_403[k] = f_9 * hsg_212[k]
                   + f_3 * pc_y[k] * isg_287[k];

        t_404[k] = pa_y[k] * hsh0_299[k]
                   - f_8 * pc_y[k] * hsh1_299[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, pa_x, pc_x, pc_y, pc_z, hsh0_405, hsg_198, \
                         hsg_215, hsg_291, hsh1_405, isg_288, isg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_x[k] * hsh0_405[k]
                   + f_10 * hsg_291[k]
                   - f_8 * pc_x[k] * hsh1_405[k];

        t_406[k] = f_15 * hsg_198[k]
                   + f_3 * pc_z[k] * isg_288[k];

        t_407[k] = f_9 * hsg_215[k]
                   + f_3 * pc_y[k] * isg_290[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_y, pc_x, pc_y, hsh0_303, hsg_295, \
                         hsg_296, hsg_297, hsh1_303, isg_295, isg_296, \
                         isg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = pa_y[k] * hsh0_303[k]
                   - f_8 * pc_y[k] * hsh1_303[k];

        t_409[k] = f_9 * hsg_295[k]
                   + f_3 * pc_x[k] * isg_295[k];

        t_410[k] = f_9 * hsg_296[k]
                   + f_3 * pc_x[k] * isg_296[k];

        t_411[k] = f_9 * hsg_297[k]
                   + f_3 * pc_x[k] * isg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_x, pc_x, pc_z, hsh0_414, hsg_205, \
                         hsg_298, hsg_299, hsh1_414, isg_295, isg_298, \
                         isg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_9 * hsg_298[k]
                   + f_3 * pc_x[k] * isg_298[k];

        t_413[k] = f_9 * hsg_299[k]
                   + f_3 * pc_x[k] * isg_299[k];

        t_414[k] = pa_x[k] * hsh0_414[k]
                   - f_8 * pc_x[k] * hsh1_414[k];

        t_415[k] = f_15 * hsg_205[k]
                   + f_3 * pc_z[k] * isg_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_x, pc_x, pc_y, hsh0_416, hsh0_417, \
                         hsh0_419, hsg_224, hsh1_416, hsh1_417, hsh1_419, \
                         isg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pa_x[k] * hsh0_416[k]
                   - f_8 * pc_x[k] * hsh1_416[k];

        t_417[k] = pa_x[k] * hsh0_417[k]
                   - f_8 * pc_x[k] * hsh1_417[k];

        t_418[k] = f_9 * hsg_224[k]
                   + f_3 * pc_y[k] * isg_299[k];

        t_419[k] = pa_x[k] * hsh0_419[k]
                   - f_8 * pc_x[k] * hsh1_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_x, pc_x, pc_y, pc_z, hsh0_420, \
                         hsg_210, hsg_300, hsh1_420, isf0_200, isf1_200, isg_300, \
                         isg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pa_x[k] * hsh0_420[k]
                   + f_12 * hsg_300[k]
                   - f_8 * pc_x[k] * hsh1_420[k];

        t_421[k] = f_3 * pc_y[k] * isg_300[k];

        t_422[k] = f_12 * hsg_210[k]
                   + f_3 * pc_z[k] * isg_300[k];

        t_423[k] = f_4 * isf0_200[k]
                   - f_5 * isf1_200[k]
                   + f_3 * pc_y[k] * isg_301[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_x, pc_x, pc_y, hsh0_425, hsg_305, hsh1_425, \
                         isf0_201, isf1_201, isg_302, isg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_y[k] * isg_302[k];

        t_425[k] = pa_x[k] * hsh0_425[k]
                   + f_11 * hsg_305[k]
                   - f_8 * pc_x[k] * hsh1_425[k];

        t_426[k] = f_6 * isf0_201[k]
                   - f_7 * isf1_201[k]
                   + f_3 * pc_y[k] * isg_303[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pa_x, pc_x, pc_y, hsh0_429, hsg_309, \
                         hsg_310, hsh1_429, isf0_202, isf1_202, isg_304, isg_305, \
                         isg_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_4 * isf0_202[k]
                   - f_5 * isf1_202[k]
                   + f_3 * pc_y[k] * isg_304[k];

        t_428[k] = f_3 * pc_y[k] * isg_305[k];

        t_429[k] = pa_x[k] * hsh0_429[k]
                   + f_10 * hsg_309[k]
                   - f_8 * pc_x[k] * hsh1_429[k];

        t_430[k] = f_9 * hsg_310[k]
                   + f_3 * pc_x[k] * isg_310[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, pc_x, pc_y, hsg_311, hsg_312, hsg_314, \
                         isg_309, isg_311, isg_312, isg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_9 * hsg_311[k]
                   + f_3 * pc_x[k] * isg_311[k];

        t_432[k] = f_9 * hsg_312[k]
                   + f_3 * pc_x[k] * isg_312[k];

        t_433[k] = f_3 * pc_y[k] * isg_309[k];

        t_434[k] = f_9 * hsg_314[k]
                   + f_3 * pc_x[k] * isg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pa_x, pc_x, hsh0_435, hsh0_436, hsh0_437, \
                         hsh0_438, hsh1_435, hsh1_436, hsh1_437, \
                         hsh1_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = pa_x[k] * hsh0_435[k]
                   - f_8 * pc_x[k] * hsh1_435[k];

        t_436[k] = pa_x[k] * hsh0_436[k]
                   - f_8 * pc_x[k] * hsh1_436[k];

        t_437[k] = pa_x[k] * hsh0_437[k]
                   - f_8 * pc_x[k] * hsh1_437[k];

        t_438[k] = pa_x[k] * hsh0_438[k]
                   - f_8 * pc_x[k] * hsh1_438[k];
    }

#pragma omp simd aligned(t_439, t_440, t_441, t_442, pa_x, pc_x, pc_y, hsh0_440, hsh1_440, \
                         isf0_210, isf0_211, isf1_210, isf1_211, isg_314, isg_315, \
                         isg_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * isg_314[k];

        t_440[k] = pa_x[k] * hsh0_440[k]
                   - f_8 * pc_x[k] * hsh1_440[k];

        t_441[k] = f_1 * isf0_210[k]
                   - f_2 * isf1_210[k]
                   + f_3 * pc_x[k] * isg_315[k];

        t_442[k] = f_13 * isf0_211[k]
                   - f_14 * isf1_211[k]
                   + f_3 * pc_x[k] * isg_316[k];
    }

#pragma omp simd aligned(t_443, t_444, t_445, t_446, pc_x, pc_z, isf0_213, isf0_215, isf1_213, \
                         isf1_215, isg_315, isg_316, isg_318, isg_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_443[k] = f_3 * pc_z[k] * isg_315[k];

        t_444[k] = f_6 * isf0_213[k]
                   - f_7 * isf1_213[k]
                   + f_3 * pc_x[k] * isg_318[k];

        t_445[k] = f_3 * pc_z[k] * isg_316[k];

        t_446[k] = f_6 * isf0_215[k]
                   - f_7 * isf1_215[k]
                   + f_3 * pc_x[k] * isg_320[k];
    }

#pragma omp simd aligned(t_447, t_448, t_449, t_450, pc_x, pc_z, isf0_216, isf0_218, isf0_219, \
                         isf1_216, isf1_218, isf1_219, isg_318, isg_321, isg_323, \
                         isg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_447[k] = f_4 * isf0_216[k]
                   - f_5 * isf1_216[k]
                   + f_3 * pc_x[k] * isg_321[k];

        t_448[k] = f_3 * pc_z[k] * isg_318[k];

        t_449[k] = f_4 * isf0_218[k]
                   - f_5 * isf1_218[k]
                   + f_3 * pc_x[k] * isg_323[k];

        t_450[k] = f_4 * isf0_219[k]
                   - f_5 * isf1_219[k]
                   + f_3 * pc_x[k] * isg_324[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, t_455, t_456, pc_x, pc_y, hsg_235, \
                         isf0_216, isf1_216, isg_325, isg_326, isg_327, isg_328, \
                         isg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_3 * pc_x[k] * isg_325[k];

        t_452[k] = f_3 * pc_x[k] * isg_326[k];

        t_453[k] = f_3 * pc_x[k] * isg_327[k];

        t_454[k] = f_3 * pc_x[k] * isg_328[k];

        t_455[k] = f_3 * pc_x[k] * isg_329[k];

        t_456[k] = f_0 * hsg_235[k]
                   + f_1 * isf0_216[k]
                   - f_2 * isf1_216[k]
                   + f_3 * pc_y[k] * isg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pc_y, pc_z, hsg_239, isf0_216, isf0_217, \
                         isf1_216, isf1_217, isg_325, isg_326, isg_327, \
                         isg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_z[k] * isg_325[k];

        t_458[k] = f_4 * isf0_216[k]
                   - f_5 * isf1_216[k]
                   + f_3 * pc_z[k] * isg_326[k];

        t_459[k] = f_6 * isf0_217[k]
                   - f_7 * isf1_217[k]
                   + f_3 * pc_z[k] * isg_327[k];

        t_460[k] = f_0 * hsg_239[k]
                   + f_3 * pc_y[k] * isg_329[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, pa_z, pc_z, hsh0_315, hsh0_316, hsh1_315, \
                         hsh1_316, isf0_219, isf1_219, isg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_1 * isf0_219[k]
                   - f_2 * isf1_219[k]
                   + f_3 * pc_z[k] * isg_329[k];

        t_462[k] = pa_z[k] * hsh0_315[k]
                   - f_8 * pc_z[k] * hsh1_315[k];

        t_463[k] = pa_z[k] * hsh0_316[k]
                   - f_8 * pc_z[k] * hsh1_316[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pa_z, pc_x, pc_z, hsh0_318, hsh1_318, isf0_222, \
                         isf0_224, isf1_222, isf1_224, isg_332, \
                         isg_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_13 * isf0_222[k]
                   - f_14 * isf1_222[k]
                   + f_3 * pc_x[k] * isg_332[k];

        t_465[k] = pa_z[k] * hsh0_318[k]
                   - f_8 * pc_z[k] * hsh1_318[k];

        t_466[k] = f_6 * isf0_224[k]
                   - f_7 * isf1_224[k]
                   + f_3 * pc_x[k] * isg_334[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pa_z, pc_x, pc_z, hsh0_321, hsh1_321, isf0_225, \
                         isf0_227, isf1_225, isf1_227, isg_335, \
                         isg_337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_6 * isf0_225[k]
                   - f_7 * isf1_225[k]
                   + f_3 * pc_x[k] * isg_335[k];

        t_468[k] = pa_z[k] * hsh0_321[k]
                   - f_8 * pc_z[k] * hsh1_321[k];

        t_469[k] = f_4 * isf0_227[k]
                   - f_5 * isf1_227[k]
                   + f_3 * pc_x[k] * isg_337[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pc_x, isf0_228, isf0_229, \
                         isf1_228, isf1_229, isg_338, isg_339, isg_340, isg_341, \
                         isg_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_4 * isf0_228[k]
                   - f_5 * isf1_228[k]
                   + f_3 * pc_x[k] * isg_338[k];

        t_471[k] = f_4 * isf0_229[k]
                   - f_5 * isf1_229[k]
                   + f_3 * pc_x[k] * isg_339[k];

        t_472[k] = f_3 * pc_x[k] * isg_340[k];

        t_473[k] = f_3 * pc_x[k] * isg_341[k];

        t_474[k] = f_3 * pc_x[k] * isg_342[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_z, pc_x, pc_z, hsh0_330, hsg_235, \
                         hsh1_330, isg_340, isg_343, isg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_3 * pc_x[k] * isg_343[k];

        t_476[k] = f_3 * pc_x[k] * isg_344[k];

        t_477[k] = pa_z[k] * hsh0_330[k]
                   - f_8 * pc_z[k] * hsh1_330[k];

        t_478[k] = f_9 * hsg_235[k]
                   + f_3 * pc_z[k] * isg_340[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pa_z, pc_y, pc_z, hsh0_332, hsh0_333, hsg_236, \
                         hsg_237, hsg_254, hsh1_332, hsh1_333, \
                         isg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = pa_z[k] * hsh0_332[k]
                   + f_10 * hsg_236[k]
                   - f_8 * pc_z[k] * hsh1_332[k];

        t_480[k] = pa_z[k] * hsh0_333[k]
                   + f_11 * hsg_237[k]
                   - f_8 * pc_z[k] * hsh1_333[k];

        t_481[k] = f_12 * hsg_254[k]
                   + f_3 * pc_y[k] * isg_344[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_z, hsg_239, isf0_229, isf0_230, \
                         isf0_231, isf1_229, isf1_230, isf1_231, isg_344, isg_345, \
                         isg_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_9 * hsg_239[k]
                   + f_1 * isf0_229[k]
                   - f_2 * isf1_229[k]
                   + f_3 * pc_z[k] * isg_344[k];

        t_483[k] = f_1 * isf0_230[k]
                   - f_2 * isf1_230[k]
                   + f_3 * pc_x[k] * isg_345[k];

        t_484[k] = f_13 * isf0_231[k]
                   - f_14 * isf1_231[k]
                   + f_3 * pc_x[k] * isg_346[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, isf0_232, isf0_233, isf0_234, isf1_232, \
                         isf1_233, isf1_234, isg_347, isg_348, \
                         isg_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_13 * isf0_232[k]
                   - f_14 * isf1_232[k]
                   + f_3 * pc_x[k] * isg_347[k];

        t_486[k] = f_6 * isf0_233[k]
                   - f_7 * isf1_233[k]
                   + f_3 * pc_x[k] * isg_348[k];

        t_487[k] = f_6 * isf0_234[k]
                   - f_7 * isf1_234[k]
                   + f_3 * pc_x[k] * isg_349[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, isf0_235, isf0_236, isf0_237, isf1_235, \
                         isf1_236, isf1_237, isg_350, isg_351, \
                         isg_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_6 * isf0_235[k]
                   - f_7 * isf1_235[k]
                   + f_3 * pc_x[k] * isg_350[k];

        t_489[k] = f_4 * isf0_236[k]
                   - f_5 * isf1_236[k]
                   + f_3 * pc_x[k] * isg_351[k];

        t_490[k] = f_4 * isf0_237[k]
                   - f_5 * isf1_237[k]
                   + f_3 * pc_x[k] * isg_352[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, isf0_238, isf0_239, \
                         isf1_238, isf1_239, isg_353, isg_354, isg_355, isg_356, \
                         isg_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_4 * isf0_238[k]
                   - f_5 * isf1_238[k]
                   + f_3 * pc_x[k] * isg_353[k];

        t_492[k] = f_4 * isf0_239[k]
                   - f_5 * isf1_239[k]
                   + f_3 * pc_x[k] * isg_354[k];

        t_493[k] = f_3 * pc_x[k] * isg_355[k];

        t_494[k] = f_3 * pc_x[k] * isg_356[k];

        t_495[k] = f_3 * pc_x[k] * isg_357[k];
    }
}

static auto
compute_prim_ish_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsh0,
                                                          const size_t hsg, const size_t hsh1,
                                                          const size_t isf0, const size_t isf1,
                                                          const size_t isg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.5 / gamma;
    const auto f_14 = 1.5 * p / (gamma * q);
    const auto f_15 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsh0_420 = buffer.data(hsh0 + 420);
    const auto *hsh0_422 = buffer.data(hsh0 + 422);
    const auto *hsh0_425 = buffer.data(hsh0 + 425);
    const auto *hsh0_429 = buffer.data(hsh0 + 429);
    const auto *hsh0_435 = buffer.data(hsh0 + 435);
    const auto *hsh0_437 = buffer.data(hsh0 + 437);
    const auto *hsh0_438 = buffer.data(hsh0 + 438);
    const auto *hsh0_440 = buffer.data(hsh0 + 440);

    const auto *hsg_250 = buffer.data(hsg + 250);
    const auto *hsg_254 = buffer.data(hsg + 254);
    const auto *hsg_265 = buffer.data(hsg + 265);
    const auto *hsg_267 = buffer.data(hsg + 267);
    const auto *hsg_268 = buffer.data(hsg + 268);
    const auto *hsg_269 = buffer.data(hsg + 269);
    const auto *hsg_280 = buffer.data(hsg + 280);
    const auto *hsg_282 = buffer.data(hsg + 282);
    const auto *hsg_283 = buffer.data(hsg + 283);
    const auto *hsg_284 = buffer.data(hsg + 284);
    const auto *hsg_295 = buffer.data(hsg + 295);
    const auto *hsg_297 = buffer.data(hsg + 297);
    const auto *hsg_298 = buffer.data(hsg + 298);
    const auto *hsg_299 = buffer.data(hsg + 299);
    const auto *hsg_310 = buffer.data(hsg + 310);
    const auto *hsg_312 = buffer.data(hsg + 312);
    const auto *hsg_313 = buffer.data(hsg + 313);
    const auto *hsg_314 = buffer.data(hsg + 314);

    const auto *hsh1_420 = buffer.data(hsh1 + 420);
    const auto *hsh1_422 = buffer.data(hsh1 + 422);
    const auto *hsh1_425 = buffer.data(hsh1 + 425);
    const auto *hsh1_429 = buffer.data(hsh1 + 429);
    const auto *hsh1_435 = buffer.data(hsh1 + 435);
    const auto *hsh1_437 = buffer.data(hsh1 + 437);
    const auto *hsh1_438 = buffer.data(hsh1 + 438);
    const auto *hsh1_440 = buffer.data(hsh1 + 440);

    const auto *isf0_236 = buffer.data(isf0 + 236);
    const auto *isf0_238 = buffer.data(isf0 + 238);
    const auto *isf0_239 = buffer.data(isf0 + 239);
    const auto *isf0_240 = buffer.data(isf0 + 240);
    const auto *isf0_241 = buffer.data(isf0 + 241);
    const auto *isf0_242 = buffer.data(isf0 + 242);
    const auto *isf0_243 = buffer.data(isf0 + 243);
    const auto *isf0_244 = buffer.data(isf0 + 244);
    const auto *isf0_245 = buffer.data(isf0 + 245);
    const auto *isf0_246 = buffer.data(isf0 + 246);
    const auto *isf0_247 = buffer.data(isf0 + 247);
    const auto *isf0_248 = buffer.data(isf0 + 248);
    const auto *isf0_249 = buffer.data(isf0 + 249);
    const auto *isf0_250 = buffer.data(isf0 + 250);
    const auto *isf0_251 = buffer.data(isf0 + 251);
    const auto *isf0_252 = buffer.data(isf0 + 252);
    const auto *isf0_253 = buffer.data(isf0 + 253);
    const auto *isf0_254 = buffer.data(isf0 + 254);
    const auto *isf0_255 = buffer.data(isf0 + 255);
    const auto *isf0_256 = buffer.data(isf0 + 256);
    const auto *isf0_257 = buffer.data(isf0 + 257);
    const auto *isf0_258 = buffer.data(isf0 + 258);
    const auto *isf0_259 = buffer.data(isf0 + 259);
    const auto *isf0_261 = buffer.data(isf0 + 261);
    const auto *isf0_263 = buffer.data(isf0 + 263);
    const auto *isf0_264 = buffer.data(isf0 + 264);
    const auto *isf0_266 = buffer.data(isf0 + 266);
    const auto *isf0_267 = buffer.data(isf0 + 267);
    const auto *isf0_268 = buffer.data(isf0 + 268);
    const auto *isf0_270 = buffer.data(isf0 + 270);
    const auto *isf0_272 = buffer.data(isf0 + 272);
    const auto *isf0_273 = buffer.data(isf0 + 273);
    const auto *isf0_275 = buffer.data(isf0 + 275);
    const auto *isf0_276 = buffer.data(isf0 + 276);
    const auto *isf0_277 = buffer.data(isf0 + 277);
    const auto *isf0_278 = buffer.data(isf0 + 278);
    const auto *isf0_279 = buffer.data(isf0 + 279);

    const auto *isf1_236 = buffer.data(isf1 + 236);
    const auto *isf1_238 = buffer.data(isf1 + 238);
    const auto *isf1_239 = buffer.data(isf1 + 239);
    const auto *isf1_240 = buffer.data(isf1 + 240);
    const auto *isf1_241 = buffer.data(isf1 + 241);
    const auto *isf1_242 = buffer.data(isf1 + 242);
    const auto *isf1_243 = buffer.data(isf1 + 243);
    const auto *isf1_244 = buffer.data(isf1 + 244);
    const auto *isf1_245 = buffer.data(isf1 + 245);
    const auto *isf1_246 = buffer.data(isf1 + 246);
    const auto *isf1_247 = buffer.data(isf1 + 247);
    const auto *isf1_248 = buffer.data(isf1 + 248);
    const auto *isf1_249 = buffer.data(isf1 + 249);
    const auto *isf1_250 = buffer.data(isf1 + 250);
    const auto *isf1_251 = buffer.data(isf1 + 251);
    const auto *isf1_252 = buffer.data(isf1 + 252);
    const auto *isf1_253 = buffer.data(isf1 + 253);
    const auto *isf1_254 = buffer.data(isf1 + 254);
    const auto *isf1_255 = buffer.data(isf1 + 255);
    const auto *isf1_256 = buffer.data(isf1 + 256);
    const auto *isf1_257 = buffer.data(isf1 + 257);
    const auto *isf1_258 = buffer.data(isf1 + 258);
    const auto *isf1_259 = buffer.data(isf1 + 259);
    const auto *isf1_261 = buffer.data(isf1 + 261);
    const auto *isf1_263 = buffer.data(isf1 + 263);
    const auto *isf1_264 = buffer.data(isf1 + 264);
    const auto *isf1_266 = buffer.data(isf1 + 266);
    const auto *isf1_267 = buffer.data(isf1 + 267);
    const auto *isf1_268 = buffer.data(isf1 + 268);
    const auto *isf1_270 = buffer.data(isf1 + 270);
    const auto *isf1_272 = buffer.data(isf1 + 272);
    const auto *isf1_273 = buffer.data(isf1 + 273);
    const auto *isf1_275 = buffer.data(isf1 + 275);
    const auto *isf1_276 = buffer.data(isf1 + 276);
    const auto *isf1_277 = buffer.data(isf1 + 277);
    const auto *isf1_278 = buffer.data(isf1 + 278);
    const auto *isf1_279 = buffer.data(isf1 + 279);

    const auto *isg_355 = buffer.data(isg + 355);
    const auto *isg_357 = buffer.data(isg + 357);
    const auto *isg_358 = buffer.data(isg + 358);
    const auto *isg_359 = buffer.data(isg + 359);
    const auto *isg_360 = buffer.data(isg + 360);
    const auto *isg_361 = buffer.data(isg + 361);
    const auto *isg_362 = buffer.data(isg + 362);
    const auto *isg_363 = buffer.data(isg + 363);
    const auto *isg_364 = buffer.data(isg + 364);
    const auto *isg_365 = buffer.data(isg + 365);
    const auto *isg_366 = buffer.data(isg + 366);
    const auto *isg_367 = buffer.data(isg + 367);
    const auto *isg_368 = buffer.data(isg + 368);
    const auto *isg_369 = buffer.data(isg + 369);
    const auto *isg_370 = buffer.data(isg + 370);
    const auto *isg_371 = buffer.data(isg + 371);
    const auto *isg_372 = buffer.data(isg + 372);
    const auto *isg_373 = buffer.data(isg + 373);
    const auto *isg_374 = buffer.data(isg + 374);
    const auto *isg_375 = buffer.data(isg + 375);
    const auto *isg_376 = buffer.data(isg + 376);
    const auto *isg_377 = buffer.data(isg + 377);
    const auto *isg_378 = buffer.data(isg + 378);
    const auto *isg_379 = buffer.data(isg + 379);
    const auto *isg_380 = buffer.data(isg + 380);
    const auto *isg_381 = buffer.data(isg + 381);
    const auto *isg_382 = buffer.data(isg + 382);
    const auto *isg_383 = buffer.data(isg + 383);
    const auto *isg_384 = buffer.data(isg + 384);
    const auto *isg_385 = buffer.data(isg + 385);
    const auto *isg_386 = buffer.data(isg + 386);
    const auto *isg_387 = buffer.data(isg + 387);
    const auto *isg_388 = buffer.data(isg + 388);
    const auto *isg_389 = buffer.data(isg + 389);
    const auto *isg_391 = buffer.data(isg + 391);
    const auto *isg_393 = buffer.data(isg + 393);
    const auto *isg_394 = buffer.data(isg + 394);
    const auto *isg_396 = buffer.data(isg + 396);
    const auto *isg_397 = buffer.data(isg + 397);
    const auto *isg_398 = buffer.data(isg + 398);
    const auto *isg_400 = buffer.data(isg + 400);
    const auto *isg_401 = buffer.data(isg + 401);
    const auto *isg_402 = buffer.data(isg + 402);
    const auto *isg_403 = buffer.data(isg + 403);
    const auto *isg_404 = buffer.data(isg + 404);
    const auto *isg_405 = buffer.data(isg + 405);
    const auto *isg_407 = buffer.data(isg + 407);
    const auto *isg_408 = buffer.data(isg + 408);
    const auto *isg_410 = buffer.data(isg + 410);
    const auto *isg_411 = buffer.data(isg + 411);
    const auto *isg_412 = buffer.data(isg + 412);
    const auto *isg_414 = buffer.data(isg + 414);
    const auto *isg_415 = buffer.data(isg + 415);
    const auto *isg_416 = buffer.data(isg + 416);
    const auto *isg_417 = buffer.data(isg + 417);
    const auto *isg_418 = buffer.data(isg + 418);
    const auto *isg_419 = buffer.data(isg + 419);

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_x, pc_y, pc_z, hsg_250, hsg_265, \
                         isf0_236, isf1_236, isg_355, isg_358, \
                         isg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_3 * pc_x[k] * isg_358[k];

        t_497[k] = f_3 * pc_x[k] * isg_359[k];

        t_498[k] = f_15 * hsg_265[k]
                   + f_1 * isf0_236[k]
                   - f_2 * isf1_236[k]
                   + f_3 * pc_y[k] * isg_355[k];

        t_499[k] = f_10 * hsg_250[k]
                   + f_3 * pc_z[k] * isg_355[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pc_y, hsg_267, hsg_268, hsg_269, isf0_238, \
                         isf0_239, isf1_238, isf1_239, isg_357, isg_358, \
                         isg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_15 * hsg_267[k]
                   + f_6 * isf0_238[k]
                   - f_7 * isf1_238[k]
                   + f_3 * pc_y[k] * isg_357[k];

        t_501[k] = f_15 * hsg_268[k]
                   + f_4 * isf0_239[k]
                   - f_5 * isf1_239[k]
                   + f_3 * pc_y[k] * isg_358[k];

        t_502[k] = f_15 * hsg_269[k]
                   + f_3 * pc_y[k] * isg_359[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pc_x, pc_z, hsg_254, isf0_239, isf0_240, \
                         isf0_241, isf1_239, isf1_240, isf1_241, isg_359, isg_360, \
                         isg_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_10 * hsg_254[k]
                   + f_1 * isf0_239[k]
                   - f_2 * isf1_239[k]
                   + f_3 * pc_z[k] * isg_359[k];

        t_504[k] = f_1 * isf0_240[k]
                   - f_2 * isf1_240[k]
                   + f_3 * pc_x[k] * isg_360[k];

        t_505[k] = f_13 * isf0_241[k]
                   - f_14 * isf1_241[k]
                   + f_3 * pc_x[k] * isg_361[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pc_x, isf0_242, isf0_243, isf0_244, isf1_242, \
                         isf1_243, isf1_244, isg_362, isg_363, \
                         isg_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_13 * isf0_242[k]
                   - f_14 * isf1_242[k]
                   + f_3 * pc_x[k] * isg_362[k];

        t_507[k] = f_6 * isf0_243[k]
                   - f_7 * isf1_243[k]
                   + f_3 * pc_x[k] * isg_363[k];

        t_508[k] = f_6 * isf0_244[k]
                   - f_7 * isf1_244[k]
                   + f_3 * pc_x[k] * isg_364[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, isf0_245, isf0_246, isf0_247, isf1_245, \
                         isf1_246, isf1_247, isg_365, isg_366, \
                         isg_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_6 * isf0_245[k]
                   - f_7 * isf1_245[k]
                   + f_3 * pc_x[k] * isg_365[k];

        t_510[k] = f_4 * isf0_246[k]
                   - f_5 * isf1_246[k]
                   + f_3 * pc_x[k] * isg_366[k];

        t_511[k] = f_4 * isf0_247[k]
                   - f_5 * isf1_247[k]
                   + f_3 * pc_x[k] * isg_367[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, t_516, pc_x, isf0_248, isf0_249, \
                         isf1_248, isf1_249, isg_368, isg_369, isg_370, isg_371, \
                         isg_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_4 * isf0_248[k]
                   - f_5 * isf1_248[k]
                   + f_3 * pc_x[k] * isg_368[k];

        t_513[k] = f_4 * isf0_249[k]
                   - f_5 * isf1_249[k]
                   + f_3 * pc_x[k] * isg_369[k];

        t_514[k] = f_3 * pc_x[k] * isg_370[k];

        t_515[k] = f_3 * pc_x[k] * isg_371[k];

        t_516[k] = f_3 * pc_x[k] * isg_372[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, pc_x, pc_y, pc_z, hsg_265, hsg_280, \
                         isf0_246, isf1_246, isg_370, isg_373, \
                         isg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_3 * pc_x[k] * isg_373[k];

        t_518[k] = f_3 * pc_x[k] * isg_374[k];

        t_519[k] = f_11 * hsg_280[k]
                   + f_1 * isf0_246[k]
                   - f_2 * isf1_246[k]
                   + f_3 * pc_y[k] * isg_370[k];

        t_520[k] = f_11 * hsg_265[k]
                   + f_3 * pc_z[k] * isg_370[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, pc_y, hsg_282, hsg_283, hsg_284, isf0_248, \
                         isf0_249, isf1_248, isf1_249, isg_372, isg_373, \
                         isg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_11 * hsg_282[k]
                   + f_6 * isf0_248[k]
                   - f_7 * isf1_248[k]
                   + f_3 * pc_y[k] * isg_372[k];

        t_522[k] = f_11 * hsg_283[k]
                   + f_4 * isf0_249[k]
                   - f_5 * isf1_249[k]
                   + f_3 * pc_y[k] * isg_373[k];

        t_523[k] = f_11 * hsg_284[k]
                   + f_3 * pc_y[k] * isg_374[k];
    }

#pragma omp simd aligned(t_524, t_525, t_526, pc_x, pc_z, hsg_269, isf0_249, isf0_250, \
                         isf0_251, isf1_249, isf1_250, isf1_251, isg_374, isg_375, \
                         isg_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_524[k] = f_11 * hsg_269[k]
                   + f_1 * isf0_249[k]
                   - f_2 * isf1_249[k]
                   + f_3 * pc_z[k] * isg_374[k];

        t_525[k] = f_1 * isf0_250[k]
                   - f_2 * isf1_250[k]
                   + f_3 * pc_x[k] * isg_375[k];

        t_526[k] = f_13 * isf0_251[k]
                   - f_14 * isf1_251[k]
                   + f_3 * pc_x[k] * isg_376[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_x, isf0_252, isf0_253, isf0_254, isf1_252, \
                         isf1_253, isf1_254, isg_377, isg_378, \
                         isg_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_13 * isf0_252[k]
                   - f_14 * isf1_252[k]
                   + f_3 * pc_x[k] * isg_377[k];

        t_528[k] = f_6 * isf0_253[k]
                   - f_7 * isf1_253[k]
                   + f_3 * pc_x[k] * isg_378[k];

        t_529[k] = f_6 * isf0_254[k]
                   - f_7 * isf1_254[k]
                   + f_3 * pc_x[k] * isg_379[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pc_x, isf0_255, isf0_256, isf0_257, isf1_255, \
                         isf1_256, isf1_257, isg_380, isg_381, \
                         isg_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_6 * isf0_255[k]
                   - f_7 * isf1_255[k]
                   + f_3 * pc_x[k] * isg_380[k];

        t_531[k] = f_4 * isf0_256[k]
                   - f_5 * isf1_256[k]
                   + f_3 * pc_x[k] * isg_381[k];

        t_532[k] = f_4 * isf0_257[k]
                   - f_5 * isf1_257[k]
                   + f_3 * pc_x[k] * isg_382[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, t_537, pc_x, isf0_258, isf0_259, \
                         isf1_258, isf1_259, isg_383, isg_384, isg_385, isg_386, \
                         isg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_4 * isf0_258[k]
                   - f_5 * isf1_258[k]
                   + f_3 * pc_x[k] * isg_383[k];

        t_534[k] = f_4 * isf0_259[k]
                   - f_5 * isf1_259[k]
                   + f_3 * pc_x[k] * isg_384[k];

        t_535[k] = f_3 * pc_x[k] * isg_385[k];

        t_536[k] = f_3 * pc_x[k] * isg_386[k];

        t_537[k] = f_3 * pc_x[k] * isg_387[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pc_x, pc_y, pc_z, hsg_280, hsg_295, \
                         isf0_256, isf1_256, isg_385, isg_388, \
                         isg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_3 * pc_x[k] * isg_388[k];

        t_539[k] = f_3 * pc_x[k] * isg_389[k];

        t_540[k] = f_10 * hsg_295[k]
                   + f_1 * isf0_256[k]
                   - f_2 * isf1_256[k]
                   + f_3 * pc_y[k] * isg_385[k];

        t_541[k] = f_15 * hsg_280[k]
                   + f_3 * pc_z[k] * isg_385[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_y, hsg_297, hsg_298, hsg_299, isf0_258, \
                         isf0_259, isf1_258, isf1_259, isg_387, isg_388, \
                         isg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_10 * hsg_297[k]
                   + f_6 * isf0_258[k]
                   - f_7 * isf1_258[k]
                   + f_3 * pc_y[k] * isg_387[k];

        t_543[k] = f_10 * hsg_298[k]
                   + f_4 * isf0_259[k]
                   - f_5 * isf1_259[k]
                   + f_3 * pc_y[k] * isg_388[k];

        t_544[k] = f_10 * hsg_299[k]
                   + f_3 * pc_y[k] * isg_389[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pa_y, pc_x, pc_y, pc_z, hsh0_420, hsg_284, \
                         hsh1_420, isf0_259, isf0_261, isf1_259, isf1_261, isg_389, \
                         isg_391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_15 * hsg_284[k]
                   + f_1 * isf0_259[k]
                   - f_2 * isf1_259[k]
                   + f_3 * pc_z[k] * isg_389[k];

        t_546[k] = pa_y[k] * hsh0_420[k]
                   - f_8 * pc_y[k] * hsh1_420[k];

        t_547[k] = f_13 * isf0_261[k]
                   - f_14 * isf1_261[k]
                   + f_3 * pc_x[k] * isg_391[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pa_y, pc_x, pc_y, hsh0_422, hsh1_422, isf0_263, \
                         isf0_264, isf1_263, isf1_264, isg_393, \
                         isg_394 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = pa_y[k] * hsh0_422[k]
                   - f_8 * pc_y[k] * hsh1_422[k];

        t_549[k] = f_6 * isf0_263[k]
                   - f_7 * isf1_263[k]
                   + f_3 * pc_x[k] * isg_393[k];

        t_550[k] = f_6 * isf0_264[k]
                   - f_7 * isf1_264[k]
                   + f_3 * pc_x[k] * isg_394[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pa_y, pc_x, pc_y, hsh0_425, hsh1_425, isf0_266, \
                         isf0_267, isf1_266, isf1_267, isg_396, \
                         isg_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = pa_y[k] * hsh0_425[k]
                   - f_8 * pc_y[k] * hsh1_425[k];

        t_552[k] = f_4 * isf0_266[k]
                   - f_5 * isf1_266[k]
                   + f_3 * pc_x[k] * isg_396[k];

        t_553[k] = f_4 * isf0_267[k]
                   - f_5 * isf1_267[k]
                   + f_3 * pc_x[k] * isg_397[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, pa_y, pc_x, pc_y, hsh0_429, \
                         hsh1_429, isf0_268, isf1_268, isg_398, isg_400, isg_401, \
                         isg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_4 * isf0_268[k]
                   - f_5 * isf1_268[k]
                   + f_3 * pc_x[k] * isg_398[k];

        t_555[k] = pa_y[k] * hsh0_429[k]
                   - f_8 * pc_y[k] * hsh1_429[k];

        t_556[k] = f_3 * pc_x[k] * isg_400[k];

        t_557[k] = f_3 * pc_x[k] * isg_401[k];

        t_558[k] = f_3 * pc_x[k] * isg_402[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, hsh0_435, \
                         hsg_295, hsg_310, hsh1_435, isg_400, isg_403, \
                         isg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_3 * pc_x[k] * isg_403[k];

        t_560[k] = f_3 * pc_x[k] * isg_404[k];

        t_561[k] = pa_y[k] * hsh0_435[k]
                   + f_12 * hsg_310[k]
                   - f_8 * pc_y[k] * hsh1_435[k];

        t_562[k] = f_12 * hsg_295[k]
                   + f_3 * pc_z[k] * isg_400[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pa_y, pc_y, hsh0_437, hsh0_438, hsh0_440, \
                         hsg_312, hsg_313, hsg_314, hsh1_437, hsh1_438, hsh1_440, \
                         isg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pa_y[k] * hsh0_437[k]
                   + f_11 * hsg_312[k]
                   - f_8 * pc_y[k] * hsh1_437[k];

        t_564[k] = pa_y[k] * hsh0_438[k]
                   + f_10 * hsg_313[k]
                   - f_8 * pc_y[k] * hsh1_438[k];

        t_565[k] = f_9 * hsg_314[k]
                   + f_3 * pc_y[k] * isg_404[k];

        t_566[k] = pa_y[k] * hsh0_440[k]
                   - f_8 * pc_y[k] * hsh1_440[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, isf0_270, isf0_272, \
                         isf0_273, isf1_270, isf1_272, isf1_273, isg_405, isg_407, \
                         isg_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_1 * isf0_270[k]
                   - f_2 * isf1_270[k]
                   + f_3 * pc_x[k] * isg_405[k];

        t_568[k] = f_3 * pc_y[k] * isg_405[k];

        t_569[k] = f_13 * isf0_272[k]
                   - f_14 * isf1_272[k]
                   + f_3 * pc_x[k] * isg_407[k];

        t_570[k] = f_6 * isf0_273[k]
                   - f_7 * isf1_273[k]
                   + f_3 * pc_x[k] * isg_408[k];

        t_571[k] = f_3 * pc_y[k] * isg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, isf0_275, isf0_276, isf0_277, \
                         isf1_275, isf1_276, isf1_277, isg_410, isg_411, \
                         isg_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_6 * isf0_275[k]
                   - f_7 * isf1_275[k]
                   + f_3 * pc_x[k] * isg_410[k];

        t_573[k] = f_4 * isf0_276[k]
                   - f_5 * isf1_276[k]
                   + f_3 * pc_x[k] * isg_411[k];

        t_574[k] = f_4 * isf0_277[k]
                   - f_5 * isf1_277[k]
                   + f_3 * pc_x[k] * isg_412[k];

        t_575[k] = f_3 * pc_y[k] * isg_410[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, pc_x, isf0_279, isf1_279, \
                         isg_414, isg_415, isg_416, isg_417, isg_418, \
                         isg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_4 * isf0_279[k]
                   - f_5 * isf1_279[k]
                   + f_3 * pc_x[k] * isg_414[k];

        t_577[k] = f_3 * pc_x[k] * isg_415[k];

        t_578[k] = f_3 * pc_x[k] * isg_416[k];

        t_579[k] = f_3 * pc_x[k] * isg_417[k];

        t_580[k] = f_3 * pc_x[k] * isg_418[k];

        t_581[k] = f_3 * pc_x[k] * isg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, isf0_276, isf0_277, isf0_278, isf1_276, \
                         isf1_277, isf1_278, isg_415, isg_416, \
                         isg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * isf0_276[k]
                   - f_2 * isf1_276[k]
                   + f_3 * pc_y[k] * isg_415[k];

        t_583[k] = f_13 * isf0_277[k]
                   - f_14 * isf1_277[k]
                   + f_3 * pc_y[k] * isg_416[k];

        t_584[k] = f_6 * isf0_278[k]
                   - f_7 * isf1_278[k]
                   + f_3 * pc_y[k] * isg_417[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_y, pc_z, hsg_314, isf0_279, isf1_279, \
                         isg_418, isg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_4 * isf0_279[k]
                   - f_5 * isf1_279[k]
                   + f_3 * pc_y[k] * isg_418[k];

        t_586[k] = f_3 * pc_y[k] * isg_419[k];

        t_587[k] = f_0 * hsg_314[k]
                   + f_1 * isf0_279[k]
                   - f_2 * isf1_279[k]
                   + f_3 * pc_z[k] * isg_419[k];
    }
}

auto
compute_prim_ish_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsh0, const size_t hsg,
                                                   const size_t hsh1, const size_t isf0,
                                                   const size_t isf1, const size_t isg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ish_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsh0, hsg,
                                                              hsh1, isf0, isf1, isg, ncols,
                                                              gamma, p, q);

    compute_prim_ish_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsh0, hsg,
                                                              hsh1, isf0, isf1, isg, ncols,
                                                              gamma, p, q);

    compute_prim_ish_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsh0, hsg,
                                                              hsh1, isf0, isf1, isg, ncols,
                                                              gamma, p, q);

    compute_prim_ish_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, hsh0, hsg,
                                                              hsh1, isf0, isf1, isg, ncols,
                                                              gamma, p, q);

    compute_prim_ish_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, hsh0, hsg,
                                                              hsh1, isf0, isf1, isg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
