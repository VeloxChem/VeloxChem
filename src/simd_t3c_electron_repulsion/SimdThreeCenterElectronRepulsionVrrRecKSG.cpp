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


#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isg0,
                                                          const size_t isf, const size_t isg1,
                                                          const size_t ksd0, const size_t ksd1,
                                                          const size_t ksf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *isg0_0 = buffer.data(isg0 + 0);
    const auto *isg0_3 = buffer.data(isg0 + 3);
    const auto *isg0_5 = buffer.data(isg0 + 5);
    const auto *isg0_10 = buffer.data(isg0 + 10);
    const auto *isg0_14 = buffer.data(isg0 + 14);
    const auto *isg0_18 = buffer.data(isg0 + 18);
    const auto *isg0_25 = buffer.data(isg0 + 25);
    const auto *isg0_30 = buffer.data(isg0 + 30);
    const auto *isg0_35 = buffer.data(isg0 + 35);
    const auto *isg0_44 = buffer.data(isg0 + 44);
    const auto *isg0_45 = buffer.data(isg0 + 45);
    const auto *isg0_48 = buffer.data(isg0 + 48);
    const auto *isg0_55 = buffer.data(isg0 + 55);
    const auto *isg0_75 = buffer.data(isg0 + 75);
    const auto *isg0_78 = buffer.data(isg0 + 78);
    const auto *isg0_80 = buffer.data(isg0 + 80);

    const auto *isf_0 = buffer.data(isf + 0);
    const auto *isf_1 = buffer.data(isf + 1);
    const auto *isf_2 = buffer.data(isf + 2);
    const auto *isf_6 = buffer.data(isf + 6);
    const auto *isf_9 = buffer.data(isf + 9);
    const auto *isf_10 = buffer.data(isf + 10);
    const auto *isf_16 = buffer.data(isf + 16);
    const auto *isf_18 = buffer.data(isf + 18);
    const auto *isf_19 = buffer.data(isf + 19);
    const auto *isf_20 = buffer.data(isf + 20);
    const auto *isf_22 = buffer.data(isf + 22);
    const auto *isf_26 = buffer.data(isf + 26);
    const auto *isf_27 = buffer.data(isf + 27);
    const auto *isf_28 = buffer.data(isf + 28);
    const auto *isf_29 = buffer.data(isf + 29);
    const auto *isf_30 = buffer.data(isf + 30);
    const auto *isf_33 = buffer.data(isf + 33);
    const auto *isf_36 = buffer.data(isf + 36);
    const auto *isf_38 = buffer.data(isf + 38);
    const auto *isf_39 = buffer.data(isf + 39);
    const auto *isf_40 = buffer.data(isf + 40);
    const auto *isf_42 = buffer.data(isf + 42);
    const auto *isf_46 = buffer.data(isf + 46);
    const auto *isf_47 = buffer.data(isf + 47);
    const auto *isf_48 = buffer.data(isf + 48);
    const auto *isf_49 = buffer.data(isf + 49);
    const auto *isf_50 = buffer.data(isf + 50);
    const auto *isf_51 = buffer.data(isf + 51);
    const auto *isf_52 = buffer.data(isf + 52);
    const auto *isf_55 = buffer.data(isf + 55);
    const auto *isf_56 = buffer.data(isf + 56);
    const auto *isf_57 = buffer.data(isf + 57);
    const auto *isf_59 = buffer.data(isf + 59);
    const auto *isf_60 = buffer.data(isf + 60);
    const auto *isf_63 = buffer.data(isf + 63);
    const auto *isf_66 = buffer.data(isf + 66);
    const auto *isf_68 = buffer.data(isf + 68);
    const auto *isf_69 = buffer.data(isf + 69);
    const auto *isf_75 = buffer.data(isf + 75);
    const auto *isf_76 = buffer.data(isf + 76);
    const auto *isf_77 = buffer.data(isf + 77);
    const auto *isf_78 = buffer.data(isf + 78);
    const auto *isf_79 = buffer.data(isf + 79);
    const auto *isf_86 = buffer.data(isf + 86);
    const auto *isf_87 = buffer.data(isf + 87);
    const auto *isf_88 = buffer.data(isf + 88);
    const auto *isf_89 = buffer.data(isf + 89);

    const auto *isg1_0 = buffer.data(isg1 + 0);
    const auto *isg1_3 = buffer.data(isg1 + 3);
    const auto *isg1_5 = buffer.data(isg1 + 5);
    const auto *isg1_10 = buffer.data(isg1 + 10);
    const auto *isg1_14 = buffer.data(isg1 + 14);
    const auto *isg1_18 = buffer.data(isg1 + 18);
    const auto *isg1_25 = buffer.data(isg1 + 25);
    const auto *isg1_30 = buffer.data(isg1 + 30);
    const auto *isg1_35 = buffer.data(isg1 + 35);
    const auto *isg1_44 = buffer.data(isg1 + 44);
    const auto *isg1_45 = buffer.data(isg1 + 45);
    const auto *isg1_48 = buffer.data(isg1 + 48);
    const auto *isg1_55 = buffer.data(isg1 + 55);
    const auto *isg1_75 = buffer.data(isg1 + 75);
    const auto *isg1_78 = buffer.data(isg1 + 78);
    const auto *isg1_80 = buffer.data(isg1 + 80);

    const auto *ksd0_0 = buffer.data(ksd0 + 0);
    const auto *ksd0_3 = buffer.data(ksd0 + 3);
    const auto *ksd0_5 = buffer.data(ksd0 + 5);
    const auto *ksd0_9 = buffer.data(ksd0 + 9);
    const auto *ksd0_16 = buffer.data(ksd0 + 16);
    const auto *ksd0_17 = buffer.data(ksd0 + 17);
    const auto *ksd0_18 = buffer.data(ksd0 + 18);
    const auto *ksd0_21 = buffer.data(ksd0 + 21);
    const auto *ksd0_23 = buffer.data(ksd0 + 23);
    const auto *ksd0_29 = buffer.data(ksd0 + 29);
    const auto *ksd0_30 = buffer.data(ksd0 + 30);
    const auto *ksd0_33 = buffer.data(ksd0 + 33);
    const auto *ksd0_34 = buffer.data(ksd0 + 34);
    const auto *ksd0_35 = buffer.data(ksd0 + 35);
    const auto *ksd0_36 = buffer.data(ksd0 + 36);
    const auto *ksd0_39 = buffer.data(ksd0 + 39);
    const auto *ksd0_41 = buffer.data(ksd0 + 41);
    const auto *ksd0_47 = buffer.data(ksd0 + 47);

    const auto *ksd1_0 = buffer.data(ksd1 + 0);
    const auto *ksd1_3 = buffer.data(ksd1 + 3);
    const auto *ksd1_5 = buffer.data(ksd1 + 5);
    const auto *ksd1_9 = buffer.data(ksd1 + 9);
    const auto *ksd1_16 = buffer.data(ksd1 + 16);
    const auto *ksd1_17 = buffer.data(ksd1 + 17);
    const auto *ksd1_18 = buffer.data(ksd1 + 18);
    const auto *ksd1_21 = buffer.data(ksd1 + 21);
    const auto *ksd1_23 = buffer.data(ksd1 + 23);
    const auto *ksd1_29 = buffer.data(ksd1 + 29);
    const auto *ksd1_30 = buffer.data(ksd1 + 30);
    const auto *ksd1_33 = buffer.data(ksd1 + 33);
    const auto *ksd1_34 = buffer.data(ksd1 + 34);
    const auto *ksd1_35 = buffer.data(ksd1 + 35);
    const auto *ksd1_36 = buffer.data(ksd1 + 36);
    const auto *ksd1_39 = buffer.data(ksd1 + 39);
    const auto *ksd1_41 = buffer.data(ksd1 + 41);
    const auto *ksd1_47 = buffer.data(ksd1 + 47);

    const auto *ksf_0 = buffer.data(ksf + 0);
    const auto *ksf_1 = buffer.data(ksf + 1);
    const auto *ksf_2 = buffer.data(ksf + 2);
    const auto *ksf_3 = buffer.data(ksf + 3);
    const auto *ksf_5 = buffer.data(ksf + 5);
    const auto *ksf_6 = buffer.data(ksf + 6);
    const auto *ksf_8 = buffer.data(ksf + 8);
    const auto *ksf_9 = buffer.data(ksf + 9);
    const auto *ksf_10 = buffer.data(ksf + 10);
    const auto *ksf_11 = buffer.data(ksf + 11);
    const auto *ksf_13 = buffer.data(ksf + 13);
    const auto *ksf_16 = buffer.data(ksf + 16);
    const auto *ksf_17 = buffer.data(ksf + 17);
    const auto *ksf_18 = buffer.data(ksf + 18);
    const auto *ksf_19 = buffer.data(ksf + 19);
    const auto *ksf_20 = buffer.data(ksf + 20);
    const auto *ksf_22 = buffer.data(ksf + 22);
    const auto *ksf_25 = buffer.data(ksf + 25);
    const auto *ksf_26 = buffer.data(ksf + 26);
    const auto *ksf_27 = buffer.data(ksf + 27);
    const auto *ksf_28 = buffer.data(ksf + 28);
    const auto *ksf_29 = buffer.data(ksf + 29);
    const auto *ksf_30 = buffer.data(ksf + 30);
    const auto *ksf_31 = buffer.data(ksf + 31);
    const auto *ksf_32 = buffer.data(ksf + 32);
    const auto *ksf_33 = buffer.data(ksf + 33);
    const auto *ksf_36 = buffer.data(ksf + 36);
    const auto *ksf_37 = buffer.data(ksf + 37);
    const auto *ksf_38 = buffer.data(ksf + 38);
    const auto *ksf_39 = buffer.data(ksf + 39);
    const auto *ksf_40 = buffer.data(ksf + 40);
    const auto *ksf_42 = buffer.data(ksf + 42);
    const auto *ksf_46 = buffer.data(ksf + 46);
    const auto *ksf_47 = buffer.data(ksf + 47);
    const auto *ksf_48 = buffer.data(ksf + 48);
    const auto *ksf_49 = buffer.data(ksf + 49);
    const auto *ksf_50 = buffer.data(ksf + 50);
    const auto *ksf_51 = buffer.data(ksf + 51);
    const auto *ksf_52 = buffer.data(ksf + 52);
    const auto *ksf_55 = buffer.data(ksf + 55);
    const auto *ksf_56 = buffer.data(ksf + 56);
    const auto *ksf_57 = buffer.data(ksf + 57);
    const auto *ksf_58 = buffer.data(ksf + 58);
    const auto *ksf_59 = buffer.data(ksf + 59);
    const auto *ksf_60 = buffer.data(ksf + 60);
    const auto *ksf_61 = buffer.data(ksf + 61);
    const auto *ksf_62 = buffer.data(ksf + 62);
    const auto *ksf_63 = buffer.data(ksf + 63);
    const auto *ksf_66 = buffer.data(ksf + 66);
    const auto *ksf_67 = buffer.data(ksf + 67);
    const auto *ksf_68 = buffer.data(ksf + 68);
    const auto *ksf_69 = buffer.data(ksf + 69);
    const auto *ksf_70 = buffer.data(ksf + 70);
    const auto *ksf_72 = buffer.data(ksf + 72);
    const auto *ksf_75 = buffer.data(ksf + 75);
    const auto *ksf_76 = buffer.data(ksf + 76);
    const auto *ksf_77 = buffer.data(ksf + 77);
    const auto *ksf_78 = buffer.data(ksf + 78);
    const auto *ksf_79 = buffer.data(ksf + 79);
    const auto *ksf_80 = buffer.data(ksf + 80);
    const auto *ksf_82 = buffer.data(ksf + 82);
    const auto *ksf_86 = buffer.data(ksf + 86);
    const auto *ksf_87 = buffer.data(ksf + 87);
    const auto *ksf_88 = buffer.data(ksf + 88);
    const auto *ksf_89 = buffer.data(ksf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, isf_0, ksd0_0, \
                         ksd1_0, ksf_0, ksf_1, ksf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * isf_0[k]
                 + f_1 * ksd0_0[k]
                 - f_2 * ksd1_0[k]
                 + f_3 * pc_x[k] * ksf_0[k];

        t_1[k] = f_3 * pc_y[k] * ksf_0[k];

        t_2[k] = f_3 * pc_z[k] * ksf_0[k];

        t_3[k] = f_4 * ksd0_0[k]
                 - f_5 * ksd1_0[k]
                 + f_3 * pc_y[k] * ksf_1[k];

        t_4[k] = f_3 * pc_y[k] * ksf_2[k];

        t_5[k] = f_4 * ksd0_0[k]
                 - f_5 * ksd1_0[k]
                 + f_3 * pc_z[k] * ksf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, isf_6, isf_9, ksd0_3, \
                         ksd1_3, ksf_3, ksf_5, ksf_6, ksf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * isf_6[k]
                 + f_3 * pc_x[k] * ksf_6[k];

        t_7[k] = f_3 * pc_z[k] * ksf_3[k];

        t_8[k] = f_3 * pc_y[k] * ksf_5[k];

        t_9[k] = f_0 * isf_9[k]
                 + f_3 * pc_x[k] * ksf_9[k];

        t_10[k] = f_1 * ksd0_3[k]
                  - f_2 * ksd1_3[k]
                  + f_3 * pc_y[k] * ksf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, isg0_0, isg1_0, \
                         ksd0_5, ksd1_5, ksf_6, ksf_8, ksf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ksf_6[k];

        t_12[k] = f_4 * ksd0_5[k]
                  - f_5 * ksd1_5[k]
                  + f_3 * pc_y[k] * ksf_8[k];

        t_13[k] = f_3 * pc_y[k] * ksf_9[k];

        t_14[k] = f_1 * ksd0_5[k]
                  - f_2 * ksd1_5[k]
                  + f_3 * pc_z[k] * ksf_9[k];

        t_15[k] = pa_y[k] * isg0_0[k]
                  - f_6 * pc_y[k] * isg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, isg0_3, isg0_5, \
                         isf_0, isf_1, isg1_3, isg1_5, ksf_10, ksf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * isf_0[k]
                  + f_3 * pc_y[k] * ksf_10[k];

        t_17[k] = f_3 * pc_z[k] * ksf_10[k];

        t_18[k] = pa_y[k] * isg0_3[k]
                  + f_8 * isf_1[k]
                  - f_6 * pc_y[k] * isg1_3[k];

        t_19[k] = f_3 * pc_z[k] * ksf_11[k];

        t_20[k] = pa_y[k] * isg0_5[k]
                  - f_6 * pc_y[k] * isg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, isf_16, isf_18, isf_19, ksf_13, \
                         ksf_16, ksf_18, ksf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * isf_16[k]
                  + f_3 * pc_x[k] * ksf_16[k];

        t_22[k] = f_3 * pc_z[k] * ksf_13[k];

        t_23[k] = f_9 * isf_18[k]
                  + f_3 * pc_x[k] * ksf_18[k];

        t_24[k] = f_9 * isf_19[k]
                  + f_3 * pc_x[k] * ksf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, isf_6, isf_9, ksd0_9, ksd1_9, \
                         ksf_16, ksf_17, ksf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * isf_6[k]
                  + f_1 * ksd0_9[k]
                  - f_2 * ksd1_9[k]
                  + f_3 * pc_y[k] * ksf_16[k];

        t_26[k] = f_3 * pc_z[k] * ksf_16[k];

        t_27[k] = f_4 * ksd0_9[k]
                  - f_5 * ksd1_9[k]
                  + f_3 * pc_z[k] * ksf_17[k];

        t_28[k] = f_7 * isf_9[k]
                  + f_3 * pc_y[k] * ksf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, isg0_0, isg0_14, \
                         isf_0, isg1_0, isg1_14, ksf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * isg0_14[k]
                  - f_6 * pc_y[k] * isg1_14[k];

        t_30[k] = pa_z[k] * isg0_0[k]
                  - f_6 * pc_z[k] * isg1_0[k];

        t_31[k] = f_3 * pc_y[k] * ksf_20[k];

        t_32[k] = f_7 * isf_0[k]
                  + f_3 * pc_z[k] * ksf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, isg0_3, isg0_5, \
                         isf_2, isf_26, isg1_3, isg1_5, ksf_22, \
                         ksf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * isg0_3[k]
                  - f_6 * pc_z[k] * isg1_3[k];

        t_34[k] = f_3 * pc_y[k] * ksf_22[k];

        t_35[k] = pa_z[k] * isg0_5[k]
                  + f_8 * isf_2[k]
                  - f_6 * pc_z[k] * isg1_5[k];

        t_36[k] = f_9 * isf_26[k]
                  + f_3 * pc_x[k] * ksf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, isg0_10, isf_27, \
                         isf_29, isg1_10, ksf_25, ksf_27, ksf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * isf_27[k]
                  + f_3 * pc_x[k] * ksf_27[k];

        t_38[k] = f_3 * pc_y[k] * ksf_25[k];

        t_39[k] = f_9 * isf_29[k]
                  + f_3 * pc_x[k] * ksf_29[k];

        t_40[k] = pa_z[k] * isg0_10[k]
                  - f_6 * pc_z[k] * isg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, isf_9, ksd0_16, ksd0_17, ksd1_16, \
                         ksd1_17, ksf_27, ksf_28, ksf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * ksd0_16[k]
                  - f_11 * ksd1_16[k]
                  + f_3 * pc_y[k] * ksf_27[k];

        t_42[k] = f_4 * ksd0_17[k]
                  - f_5 * ksd1_17[k]
                  + f_3 * pc_y[k] * ksf_28[k];

        t_43[k] = f_3 * pc_y[k] * ksf_29[k];

        t_44[k] = f_7 * isf_9[k]
                  + f_1 * ksd0_17[k]
                  - f_2 * ksd1_17[k]
                  + f_3 * pc_z[k] * ksf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, isf_10, isf_30, isf_33, \
                         ksd0_18, ksd0_21, ksd1_18, ksd1_21, ksf_30, \
                         ksf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * isf_30[k]
                  + f_1 * ksd0_18[k]
                  - f_2 * ksd1_18[k]
                  + f_3 * pc_x[k] * ksf_30[k];

        t_46[k] = f_8 * isf_10[k]
                  + f_3 * pc_y[k] * ksf_30[k];

        t_47[k] = f_3 * pc_z[k] * ksf_30[k];

        t_48[k] = f_12 * isf_33[k]
                  + f_4 * ksd0_21[k]
                  - f_5 * ksd1_21[k]
                  + f_3 * pc_x[k] * ksf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, isf_36, isf_38, ksd0_18, \
                         ksd1_18, ksf_31, ksf_32, ksf_33, ksf_36, \
                         ksf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * ksf_31[k];

        t_50[k] = f_4 * ksd0_18[k]
                  - f_5 * ksd1_18[k]
                  + f_3 * pc_z[k] * ksf_32[k];

        t_51[k] = f_12 * isf_36[k]
                  + f_3 * pc_x[k] * ksf_36[k];

        t_52[k] = f_3 * pc_z[k] * ksf_33[k];

        t_53[k] = f_12 * isf_38[k]
                  + f_3 * pc_x[k] * ksf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, isf_16, isf_19, \
                         isf_39, ksd0_21, ksd1_21, ksf_36, ksf_37, \
                         ksf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * isf_39[k]
                  + f_3 * pc_x[k] * ksf_39[k];

        t_55[k] = f_8 * isf_16[k]
                  + f_1 * ksd0_21[k]
                  - f_2 * ksd1_21[k]
                  + f_3 * pc_y[k] * ksf_36[k];

        t_56[k] = f_3 * pc_z[k] * ksf_36[k];

        t_57[k] = f_4 * ksd0_21[k]
                  - f_5 * ksd1_21[k]
                  + f_3 * pc_z[k] * ksf_37[k];

        t_58[k] = f_8 * isf_19[k]
                  + f_3 * pc_y[k] * ksf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, isg0_30, isf_10, isf_20, \
                         isg1_30, ksd0_23, ksd1_23, ksf_39, ksf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * ksd0_23[k]
                  - f_2 * ksd1_23[k]
                  + f_3 * pc_z[k] * ksf_39[k];

        t_60[k] = pa_y[k] * isg0_30[k]
                  - f_6 * pc_y[k] * isg1_30[k];

        t_61[k] = f_7 * isf_20[k]
                  + f_3 * pc_y[k] * ksf_40[k];

        t_62[k] = f_7 * isf_10[k]
                  + f_3 * pc_z[k] * ksf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, isg0_18, isg0_35, isf_22, \
                         isg1_18, isg1_35, ksf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * isg0_18[k]
                  - f_6 * pc_z[k] * isg1_18[k];

        t_64[k] = f_7 * isf_22[k]
                  + f_3 * pc_y[k] * ksf_42[k];

        t_65[k] = pa_y[k] * isg0_35[k]
                  - f_6 * pc_y[k] * isg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, isf_46, isf_47, isf_48, isf_49, ksf_46, \
                         ksf_47, ksf_48, ksf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * isf_46[k]
                  + f_3 * pc_x[k] * ksf_46[k];

        t_67[k] = f_12 * isf_47[k]
                  + f_3 * pc_x[k] * ksf_47[k];

        t_68[k] = f_12 * isf_48[k]
                  + f_3 * pc_x[k] * ksf_48[k];

        t_69[k] = f_12 * isf_49[k]
                  + f_3 * pc_x[k] * ksf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, isg0_25, isf_16, isf_28, isg1_25, \
                         ksd0_29, ksd1_29, ksf_46, ksf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * isg0_25[k]
                  - f_6 * pc_z[k] * isg1_25[k];

        t_71[k] = f_7 * isf_16[k]
                  + f_3 * pc_z[k] * ksf_46[k];

        t_72[k] = f_7 * isf_28[k]
                  + f_4 * ksd0_29[k]
                  - f_5 * ksd1_29[k]
                  + f_3 * pc_y[k] * ksf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, isg0_44, isf_29, isf_50, \
                         isg1_44, ksd0_30, ksd1_30, ksf_49, ksf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * isf_29[k]
                  + f_3 * pc_y[k] * ksf_49[k];

        t_74[k] = pa_y[k] * isg0_44[k]
                  - f_6 * pc_y[k] * isg1_44[k];

        t_75[k] = f_12 * isf_50[k]
                  + f_1 * ksd0_30[k]
                  - f_2 * ksd1_30[k]
                  + f_3 * pc_x[k] * ksf_50[k];

        t_76[k] = f_3 * pc_y[k] * ksf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, isf_20, ksd0_30, ksd1_30, ksf_50, \
                         ksf_51, ksf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * isf_20[k]
                  + f_3 * pc_z[k] * ksf_50[k];

        t_78[k] = f_4 * ksd0_30[k]
                  - f_5 * ksd1_30[k]
                  + f_3 * pc_y[k] * ksf_51[k];

        t_79[k] = f_3 * pc_y[k] * ksf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, isf_55, isf_56, isf_57, ksd0_35, \
                         ksd1_35, ksf_55, ksf_56, ksf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * isf_55[k]
                  + f_4 * ksd0_35[k]
                  - f_5 * ksd1_35[k]
                  + f_3 * pc_x[k] * ksf_55[k];

        t_81[k] = f_12 * isf_56[k]
                  + f_3 * pc_x[k] * ksf_56[k];

        t_82[k] = f_12 * isf_57[k]
                  + f_3 * pc_x[k] * ksf_57[k];

        t_83[k] = f_3 * pc_y[k] * ksf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, isf_59, ksd0_33, ksd0_34, ksd1_33, \
                         ksd1_34, ksf_56, ksf_57, ksf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * isf_59[k]
                  + f_3 * pc_x[k] * ksf_59[k];

        t_85[k] = f_1 * ksd0_33[k]
                  - f_2 * ksd1_33[k]
                  + f_3 * pc_y[k] * ksf_56[k];

        t_86[k] = f_10 * ksd0_34[k]
                  - f_11 * ksd1_34[k]
                  + f_3 * pc_y[k] * ksf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, isf_29, isf_60, ksd0_35, \
                         ksd0_36, ksd1_35, ksd1_36, ksf_58, ksf_59, \
                         ksf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * ksd0_35[k]
                  - f_5 * ksd1_35[k]
                  + f_3 * pc_y[k] * ksf_58[k];

        t_88[k] = f_3 * pc_y[k] * ksf_59[k];

        t_89[k] = f_8 * isf_29[k]
                  + f_1 * ksd0_35[k]
                  - f_2 * ksd1_35[k]
                  + f_3 * pc_z[k] * ksf_59[k];

        t_90[k] = f_13 * isf_60[k]
                  + f_1 * ksd0_36[k]
                  - f_2 * ksd1_36[k]
                  + f_3 * pc_x[k] * ksf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, isf_30, isf_63, ksd0_39, \
                         ksd1_39, ksf_60, ksf_61, ksf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * isf_30[k]
                  + f_3 * pc_y[k] * ksf_60[k];

        t_92[k] = f_3 * pc_z[k] * ksf_60[k];

        t_93[k] = f_13 * isf_63[k]
                  + f_4 * ksd0_39[k]
                  - f_5 * ksd1_39[k]
                  + f_3 * pc_x[k] * ksf_63[k];

        t_94[k] = f_3 * pc_z[k] * ksf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, isf_66, isf_68, ksd0_36, ksd1_36, \
                         ksf_62, ksf_63, ksf_66, ksf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * ksd0_36[k]
                  - f_5 * ksd1_36[k]
                  + f_3 * pc_z[k] * ksf_62[k];

        t_96[k] = f_13 * isf_66[k]
                  + f_3 * pc_x[k] * ksf_66[k];

        t_97[k] = f_3 * pc_z[k] * ksf_63[k];

        t_98[k] = f_13 * isf_68[k]
                  + f_3 * pc_x[k] * ksf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, isf_36, isf_39, \
                         isf_69, ksd0_39, ksd1_39, ksf_66, ksf_67, \
                         ksf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * isf_69[k]
                  + f_3 * pc_x[k] * ksf_69[k];

        t_100[k] = f_14 * isf_36[k]
                   + f_1 * ksd0_39[k]
                   - f_2 * ksd1_39[k]
                   + f_3 * pc_y[k] * ksf_66[k];

        t_101[k] = f_3 * pc_z[k] * ksf_66[k];

        t_102[k] = f_4 * ksd0_39[k]
                   - f_5 * ksd1_39[k]
                   + f_3 * pc_z[k] * ksf_67[k];

        t_103[k] = f_14 * isf_39[k]
                   + f_3 * pc_y[k] * ksf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, isg0_45, isf_30, \
                         isf_40, isg1_45, ksd0_41, ksd1_41, ksf_69, \
                         ksf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * ksd0_41[k]
                   - f_2 * ksd1_41[k]
                   + f_3 * pc_z[k] * ksf_69[k];

        t_105[k] = pa_z[k] * isg0_45[k]
                   - f_6 * pc_z[k] * isg1_45[k];

        t_106[k] = f_8 * isf_40[k]
                   + f_3 * pc_y[k] * ksf_70[k];

        t_107[k] = f_7 * isf_30[k]
                   + f_3 * pc_z[k] * ksf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, isg0_48, isf_42, isf_75, \
                         isg1_48, ksd0_47, ksd1_47, ksf_72, ksf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * isg0_48[k]
                   - f_6 * pc_z[k] * isg1_48[k];

        t_109[k] = f_8 * isf_42[k]
                   + f_3 * pc_y[k] * ksf_72[k];

        t_110[k] = f_13 * isf_75[k]
                   + f_4 * ksd0_47[k]
                   - f_5 * ksd1_47[k]
                   + f_3 * pc_x[k] * ksf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, isf_76, isf_77, isf_78, isf_79, \
                         ksf_76, ksf_77, ksf_78, ksf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * isf_76[k]
                   + f_3 * pc_x[k] * ksf_76[k];

        t_112[k] = f_13 * isf_77[k]
                   + f_3 * pc_x[k] * ksf_77[k];

        t_113[k] = f_13 * isf_78[k]
                   + f_3 * pc_x[k] * ksf_78[k];

        t_114[k] = f_13 * isf_79[k]
                   + f_3 * pc_x[k] * ksf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, isg0_55, isf_36, isf_48, \
                         isg1_55, ksd0_47, ksd1_47, ksf_76, ksf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * isg0_55[k]
                   - f_6 * pc_z[k] * isg1_55[k];

        t_116[k] = f_7 * isf_36[k]
                   + f_3 * pc_z[k] * ksf_76[k];

        t_117[k] = f_8 * isf_48[k]
                   + f_4 * ksd0_47[k]
                   - f_5 * ksd1_47[k]
                   + f_3 * pc_y[k] * ksf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, isg0_75, isf_39, \
                         isf_49, isf_50, isg1_75, ksd0_47, ksd1_47, ksf_79, \
                         ksf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * isf_49[k]
                   + f_3 * pc_y[k] * ksf_79[k];

        t_119[k] = f_7 * isf_39[k]
                   + f_1 * ksd0_47[k]
                   - f_2 * ksd1_47[k]
                   + f_3 * pc_z[k] * ksf_79[k];

        t_120[k] = pa_y[k] * isg0_75[k]
                   - f_6 * pc_y[k] * isg1_75[k];

        t_121[k] = f_7 * isf_50[k]
                   + f_3 * pc_y[k] * ksf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, isg0_78, isg0_80, \
                         isf_40, isf_51, isf_52, isg1_78, isg1_80, ksf_80, \
                         ksf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * isf_40[k]
                   + f_3 * pc_z[k] * ksf_80[k];

        t_123[k] = pa_y[k] * isg0_78[k]
                   + f_8 * isf_51[k]
                   - f_6 * pc_y[k] * isg1_78[k];

        t_124[k] = f_7 * isf_52[k]
                   + f_3 * pc_y[k] * ksf_82[k];

        t_125[k] = pa_y[k] * isg0_80[k]
                   - f_6 * pc_y[k] * isg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, isf_86, isf_87, isf_88, isf_89, \
                         ksf_86, ksf_87, ksf_88, ksf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * isf_86[k]
                   + f_3 * pc_x[k] * ksf_86[k];

        t_127[k] = f_13 * isf_87[k]
                   + f_3 * pc_x[k] * ksf_87[k];

        t_128[k] = f_13 * isf_88[k]
                   + f_3 * pc_x[k] * ksf_88[k];

        t_129[k] = f_13 * isf_89[k]
                   + f_3 * pc_x[k] * ksf_89[k];
    }
}

static auto
compute_prim_ksg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isg0,
                                                          const size_t isf, const size_t isg1,
                                                          const size_t ksd0, const size_t ksd1,
                                                          const size_t ksf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isg0_89 = buffer.data(isg0 + 89);
    const auto *isg0_90 = buffer.data(isg0 + 90);
    const auto *isg0_93 = buffer.data(isg0 + 93);
    const auto *isg0_100 = buffer.data(isg0 + 100);
    const auto *isg0_135 = buffer.data(isg0 + 135);
    const auto *isg0_138 = buffer.data(isg0 + 138);
    const auto *isg0_140 = buffer.data(isg0 + 140);
    const auto *isg0_149 = buffer.data(isg0 + 149);
    const auto *isg0_150 = buffer.data(isg0 + 150);
    const auto *isg0_153 = buffer.data(isg0 + 153);
    const auto *isg0_160 = buffer.data(isg0 + 160);

    const auto *isf_46 = buffer.data(isf + 46);
    const auto *isf_50 = buffer.data(isf + 50);
    const auto *isf_56 = buffer.data(isf + 56);
    const auto *isf_58 = buffer.data(isf + 58);
    const auto *isf_59 = buffer.data(isf + 59);
    const auto *isf_60 = buffer.data(isf + 60);
    const auto *isf_66 = buffer.data(isf + 66);
    const auto *isf_69 = buffer.data(isf + 69);
    const auto *isf_70 = buffer.data(isf + 70);
    const auto *isf_72 = buffer.data(isf + 72);
    const auto *isf_76 = buffer.data(isf + 76);
    const auto *isf_78 = buffer.data(isf + 78);
    const auto *isf_79 = buffer.data(isf + 79);
    const auto *isf_80 = buffer.data(isf + 80);
    const auto *isf_82 = buffer.data(isf + 82);
    const auto *isf_86 = buffer.data(isf + 86);
    const auto *isf_88 = buffer.data(isf + 88);
    const auto *isf_89 = buffer.data(isf + 89);
    const auto *isf_90 = buffer.data(isf + 90);
    const auto *isf_91 = buffer.data(isf + 91);
    const auto *isf_92 = buffer.data(isf + 92);
    const auto *isf_95 = buffer.data(isf + 95);
    const auto *isf_96 = buffer.data(isf + 96);
    const auto *isf_97 = buffer.data(isf + 97);
    const auto *isf_98 = buffer.data(isf + 98);
    const auto *isf_99 = buffer.data(isf + 99);
    const auto *isf_100 = buffer.data(isf + 100);
    const auto *isf_103 = buffer.data(isf + 103);
    const auto *isf_106 = buffer.data(isf + 106);
    const auto *isf_108 = buffer.data(isf + 108);
    const auto *isf_109 = buffer.data(isf + 109);
    const auto *isf_110 = buffer.data(isf + 110);
    const auto *isf_112 = buffer.data(isf + 112);
    const auto *isf_115 = buffer.data(isf + 115);
    const auto *isf_116 = buffer.data(isf + 116);
    const auto *isf_117 = buffer.data(isf + 117);
    const auto *isf_118 = buffer.data(isf + 118);
    const auto *isf_119 = buffer.data(isf + 119);
    const auto *isf_120 = buffer.data(isf + 120);
    const auto *isf_123 = buffer.data(isf + 123);
    const auto *isf_125 = buffer.data(isf + 125);
    const auto *isf_126 = buffer.data(isf + 126);
    const auto *isf_127 = buffer.data(isf + 127);
    const auto *isf_128 = buffer.data(isf + 128);
    const auto *isf_129 = buffer.data(isf + 129);
    const auto *isf_136 = buffer.data(isf + 136);
    const auto *isf_137 = buffer.data(isf + 137);
    const auto *isf_138 = buffer.data(isf + 138);
    const auto *isf_139 = buffer.data(isf + 139);
    const auto *isf_140 = buffer.data(isf + 140);
    const auto *isf_145 = buffer.data(isf + 145);
    const auto *isf_146 = buffer.data(isf + 146);
    const auto *isf_147 = buffer.data(isf + 147);
    const auto *isf_149 = buffer.data(isf + 149);
    const auto *isf_150 = buffer.data(isf + 150);
    const auto *isf_153 = buffer.data(isf + 153);
    const auto *isf_156 = buffer.data(isf + 156);
    const auto *isf_158 = buffer.data(isf + 158);
    const auto *isf_159 = buffer.data(isf + 159);
    const auto *isf_165 = buffer.data(isf + 165);
    const auto *isf_166 = buffer.data(isf + 166);
    const auto *isf_167 = buffer.data(isf + 167);
    const auto *isf_168 = buffer.data(isf + 168);
    const auto *isf_169 = buffer.data(isf + 169);

    const auto *isg1_89 = buffer.data(isg1 + 89);
    const auto *isg1_90 = buffer.data(isg1 + 90);
    const auto *isg1_93 = buffer.data(isg1 + 93);
    const auto *isg1_100 = buffer.data(isg1 + 100);
    const auto *isg1_135 = buffer.data(isg1 + 135);
    const auto *isg1_138 = buffer.data(isg1 + 138);
    const auto *isg1_140 = buffer.data(isg1 + 140);
    const auto *isg1_149 = buffer.data(isg1 + 149);
    const auto *isg1_150 = buffer.data(isg1 + 150);
    const auto *isg1_153 = buffer.data(isg1 + 153);
    const auto *isg1_160 = buffer.data(isg1 + 160);

    const auto *ksd0_51 = buffer.data(ksd0 + 51);
    const auto *ksd0_53 = buffer.data(ksd0 + 53);
    const auto *ksd0_54 = buffer.data(ksd0 + 54);
    const auto *ksd0_57 = buffer.data(ksd0 + 57);
    const auto *ksd0_58 = buffer.data(ksd0 + 58);
    const auto *ksd0_59 = buffer.data(ksd0 + 59);
    const auto *ksd0_60 = buffer.data(ksd0 + 60);
    const auto *ksd0_63 = buffer.data(ksd0 + 63);
    const auto *ksd0_65 = buffer.data(ksd0 + 65);
    const auto *ksd0_71 = buffer.data(ksd0 + 71);
    const auto *ksd0_72 = buffer.data(ksd0 + 72);
    const auto *ksd0_75 = buffer.data(ksd0 + 75);
    const auto *ksd0_77 = buffer.data(ksd0 + 77);
    const auto *ksd0_81 = buffer.data(ksd0 + 81);
    const auto *ksd0_83 = buffer.data(ksd0 + 83);
    const auto *ksd0_84 = buffer.data(ksd0 + 84);
    const auto *ksd0_87 = buffer.data(ksd0 + 87);
    const auto *ksd0_88 = buffer.data(ksd0 + 88);
    const auto *ksd0_89 = buffer.data(ksd0 + 89);
    const auto *ksd0_90 = buffer.data(ksd0 + 90);
    const auto *ksd0_93 = buffer.data(ksd0 + 93);
    const auto *ksd0_95 = buffer.data(ksd0 + 95);
    const auto *ksd0_101 = buffer.data(ksd0 + 101);

    const auto *ksd1_51 = buffer.data(ksd1 + 51);
    const auto *ksd1_53 = buffer.data(ksd1 + 53);
    const auto *ksd1_54 = buffer.data(ksd1 + 54);
    const auto *ksd1_57 = buffer.data(ksd1 + 57);
    const auto *ksd1_58 = buffer.data(ksd1 + 58);
    const auto *ksd1_59 = buffer.data(ksd1 + 59);
    const auto *ksd1_60 = buffer.data(ksd1 + 60);
    const auto *ksd1_63 = buffer.data(ksd1 + 63);
    const auto *ksd1_65 = buffer.data(ksd1 + 65);
    const auto *ksd1_71 = buffer.data(ksd1 + 71);
    const auto *ksd1_72 = buffer.data(ksd1 + 72);
    const auto *ksd1_75 = buffer.data(ksd1 + 75);
    const auto *ksd1_77 = buffer.data(ksd1 + 77);
    const auto *ksd1_81 = buffer.data(ksd1 + 81);
    const auto *ksd1_83 = buffer.data(ksd1 + 83);
    const auto *ksd1_84 = buffer.data(ksd1 + 84);
    const auto *ksd1_87 = buffer.data(ksd1 + 87);
    const auto *ksd1_88 = buffer.data(ksd1 + 88);
    const auto *ksd1_89 = buffer.data(ksd1 + 89);
    const auto *ksd1_90 = buffer.data(ksd1 + 90);
    const auto *ksd1_93 = buffer.data(ksd1 + 93);
    const auto *ksd1_95 = buffer.data(ksd1 + 95);
    const auto *ksd1_101 = buffer.data(ksd1 + 101);

    const auto *ksf_86 = buffer.data(ksf + 86);
    const auto *ksf_88 = buffer.data(ksf + 88);
    const auto *ksf_89 = buffer.data(ksf + 89);
    const auto *ksf_90 = buffer.data(ksf + 90);
    const auto *ksf_91 = buffer.data(ksf + 91);
    const auto *ksf_92 = buffer.data(ksf + 92);
    const auto *ksf_95 = buffer.data(ksf + 95);
    const auto *ksf_96 = buffer.data(ksf + 96);
    const auto *ksf_97 = buffer.data(ksf + 97);
    const auto *ksf_98 = buffer.data(ksf + 98);
    const auto *ksf_99 = buffer.data(ksf + 99);
    const auto *ksf_100 = buffer.data(ksf + 100);
    const auto *ksf_101 = buffer.data(ksf + 101);
    const auto *ksf_102 = buffer.data(ksf + 102);
    const auto *ksf_103 = buffer.data(ksf + 103);
    const auto *ksf_106 = buffer.data(ksf + 106);
    const auto *ksf_107 = buffer.data(ksf + 107);
    const auto *ksf_108 = buffer.data(ksf + 108);
    const auto *ksf_109 = buffer.data(ksf + 109);
    const auto *ksf_110 = buffer.data(ksf + 110);
    const auto *ksf_112 = buffer.data(ksf + 112);
    const auto *ksf_115 = buffer.data(ksf + 115);
    const auto *ksf_116 = buffer.data(ksf + 116);
    const auto *ksf_117 = buffer.data(ksf + 117);
    const auto *ksf_118 = buffer.data(ksf + 118);
    const auto *ksf_119 = buffer.data(ksf + 119);
    const auto *ksf_120 = buffer.data(ksf + 120);
    const auto *ksf_122 = buffer.data(ksf + 122);
    const auto *ksf_123 = buffer.data(ksf + 123);
    const auto *ksf_125 = buffer.data(ksf + 125);
    const auto *ksf_126 = buffer.data(ksf + 126);
    const auto *ksf_127 = buffer.data(ksf + 127);
    const auto *ksf_128 = buffer.data(ksf + 128);
    const auto *ksf_129 = buffer.data(ksf + 129);
    const auto *ksf_130 = buffer.data(ksf + 130);
    const auto *ksf_132 = buffer.data(ksf + 132);
    const auto *ksf_136 = buffer.data(ksf + 136);
    const auto *ksf_137 = buffer.data(ksf + 137);
    const auto *ksf_138 = buffer.data(ksf + 138);
    const auto *ksf_139 = buffer.data(ksf + 139);
    const auto *ksf_140 = buffer.data(ksf + 140);
    const auto *ksf_141 = buffer.data(ksf + 141);
    const auto *ksf_142 = buffer.data(ksf + 142);
    const auto *ksf_145 = buffer.data(ksf + 145);
    const auto *ksf_146 = buffer.data(ksf + 146);
    const auto *ksf_147 = buffer.data(ksf + 147);
    const auto *ksf_148 = buffer.data(ksf + 148);
    const auto *ksf_149 = buffer.data(ksf + 149);
    const auto *ksf_150 = buffer.data(ksf + 150);
    const auto *ksf_151 = buffer.data(ksf + 151);
    const auto *ksf_152 = buffer.data(ksf + 152);
    const auto *ksf_153 = buffer.data(ksf + 153);
    const auto *ksf_156 = buffer.data(ksf + 156);
    const auto *ksf_157 = buffer.data(ksf + 157);
    const auto *ksf_158 = buffer.data(ksf + 158);
    const auto *ksf_159 = buffer.data(ksf + 159);
    const auto *ksf_160 = buffer.data(ksf + 160);
    const auto *ksf_162 = buffer.data(ksf + 162);
    const auto *ksf_165 = buffer.data(ksf + 165);
    const auto *ksf_166 = buffer.data(ksf + 166);
    const auto *ksf_167 = buffer.data(ksf + 167);
    const auto *ksf_168 = buffer.data(ksf + 168);
    const auto *ksf_169 = buffer.data(ksf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, isf_46, isf_56, isf_58, ksd0_51, \
                         ksd0_53, ksd1_51, ksd1_53, ksf_86, ksf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * isf_56[k]
                   + f_1 * ksd0_51[k]
                   - f_2 * ksd1_51[k]
                   + f_3 * pc_y[k] * ksf_86[k];

        t_131[k] = f_8 * isf_46[k]
                   + f_3 * pc_z[k] * ksf_86[k];

        t_132[k] = f_7 * isf_58[k]
                   + f_4 * ksd0_53[k]
                   - f_5 * ksd1_53[k]
                   + f_3 * pc_y[k] * ksf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, isg0_89, isf_59, \
                         isf_90, isg1_89, ksd0_54, ksd1_54, ksf_89, \
                         ksf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * isf_59[k]
                   + f_3 * pc_y[k] * ksf_89[k];

        t_134[k] = pa_y[k] * isg0_89[k]
                   - f_6 * pc_y[k] * isg1_89[k];

        t_135[k] = f_13 * isf_90[k]
                   + f_1 * ksd0_54[k]
                   - f_2 * ksd1_54[k]
                   + f_3 * pc_x[k] * ksf_90[k];

        t_136[k] = f_3 * pc_y[k] * ksf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, isf_50, ksd0_54, ksd1_54, ksf_90, \
                         ksf_91, ksf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * isf_50[k]
                   + f_3 * pc_z[k] * ksf_90[k];

        t_138[k] = f_4 * ksd0_54[k]
                   - f_5 * ksd1_54[k]
                   + f_3 * pc_y[k] * ksf_91[k];

        t_139[k] = f_3 * pc_y[k] * ksf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, isf_95, isf_96, isf_97, \
                         ksd0_59, ksd1_59, ksf_95, ksf_96, ksf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * isf_95[k]
                   + f_4 * ksd0_59[k]
                   - f_5 * ksd1_59[k]
                   + f_3 * pc_x[k] * ksf_95[k];

        t_141[k] = f_13 * isf_96[k]
                   + f_3 * pc_x[k] * ksf_96[k];

        t_142[k] = f_13 * isf_97[k]
                   + f_3 * pc_x[k] * ksf_97[k];

        t_143[k] = f_3 * pc_y[k] * ksf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, isf_99, ksd0_57, ksd0_58, ksd1_57, \
                         ksd1_58, ksf_96, ksf_97, ksf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * isf_99[k]
                   + f_3 * pc_x[k] * ksf_99[k];

        t_145[k] = f_1 * ksd0_57[k]
                   - f_2 * ksd1_57[k]
                   + f_3 * pc_y[k] * ksf_96[k];

        t_146[k] = f_10 * ksd0_58[k]
                   - f_11 * ksd1_58[k]
                   + f_3 * pc_y[k] * ksf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, isf_59, isf_100, \
                         ksd0_59, ksd0_60, ksd1_59, ksd1_60, ksf_98, ksf_99, \
                         ksf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * ksd0_59[k]
                   - f_5 * ksd1_59[k]
                   + f_3 * pc_y[k] * ksf_98[k];

        t_148[k] = f_3 * pc_y[k] * ksf_99[k];

        t_149[k] = f_14 * isf_59[k]
                   + f_1 * ksd0_59[k]
                   - f_2 * ksd1_59[k]
                   + f_3 * pc_z[k] * ksf_99[k];

        t_150[k] = f_14 * isf_100[k]
                   + f_1 * ksd0_60[k]
                   - f_2 * ksd1_60[k]
                   + f_3 * pc_x[k] * ksf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, isf_60, isf_103, \
                         ksd0_63, ksd1_63, ksf_100, ksf_101, ksf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_13 * isf_60[k]
                   + f_3 * pc_y[k] * ksf_100[k];

        t_152[k] = f_3 * pc_z[k] * ksf_100[k];

        t_153[k] = f_14 * isf_103[k]
                   + f_4 * ksd0_63[k]
                   - f_5 * ksd1_63[k]
                   + f_3 * pc_x[k] * ksf_103[k];

        t_154[k] = f_3 * pc_z[k] * ksf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, isf_106, isf_108, ksd0_60, \
                         ksd1_60, ksf_102, ksf_103, ksf_106, ksf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * ksd0_60[k]
                   - f_5 * ksd1_60[k]
                   + f_3 * pc_z[k] * ksf_102[k];

        t_156[k] = f_14 * isf_106[k]
                   + f_3 * pc_x[k] * ksf_106[k];

        t_157[k] = f_3 * pc_z[k] * ksf_103[k];

        t_158[k] = f_14 * isf_108[k]
                   + f_3 * pc_x[k] * ksf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, isf_66, isf_69, \
                         isf_109, ksd0_63, ksd1_63, ksf_106, ksf_107, \
                         ksf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_14 * isf_109[k]
                   + f_3 * pc_x[k] * ksf_109[k];

        t_160[k] = f_13 * isf_66[k]
                   + f_1 * ksd0_63[k]
                   - f_2 * ksd1_63[k]
                   + f_3 * pc_y[k] * ksf_106[k];

        t_161[k] = f_3 * pc_z[k] * ksf_106[k];

        t_162[k] = f_4 * ksd0_63[k]
                   - f_5 * ksd1_63[k]
                   + f_3 * pc_z[k] * ksf_107[k];

        t_163[k] = f_13 * isf_69[k]
                   + f_3 * pc_y[k] * ksf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, isg0_90, isf_60, \
                         isf_70, isg1_90, ksd0_65, ksd1_65, ksf_109, \
                         ksf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * ksd0_65[k]
                   - f_2 * ksd1_65[k]
                   + f_3 * pc_z[k] * ksf_109[k];

        t_165[k] = pa_z[k] * isg0_90[k]
                   - f_6 * pc_z[k] * isg1_90[k];

        t_166[k] = f_14 * isf_70[k]
                   + f_3 * pc_y[k] * ksf_110[k];

        t_167[k] = f_7 * isf_60[k]
                   + f_3 * pc_z[k] * ksf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, isg0_93, isf_72, \
                         isf_115, isg1_93, ksd0_71, ksd1_71, ksf_112, \
                         ksf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * isg0_93[k]
                   - f_6 * pc_z[k] * isg1_93[k];

        t_169[k] = f_14 * isf_72[k]
                   + f_3 * pc_y[k] * ksf_112[k];

        t_170[k] = f_14 * isf_115[k]
                   + f_4 * ksd0_71[k]
                   - f_5 * ksd1_71[k]
                   + f_3 * pc_x[k] * ksf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, isf_116, isf_117, isf_118, isf_119, \
                         ksf_116, ksf_117, ksf_118, ksf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_14 * isf_116[k]
                   + f_3 * pc_x[k] * ksf_116[k];

        t_172[k] = f_14 * isf_117[k]
                   + f_3 * pc_x[k] * ksf_117[k];

        t_173[k] = f_14 * isf_118[k]
                   + f_3 * pc_x[k] * ksf_118[k];

        t_174[k] = f_14 * isf_119[k]
                   + f_3 * pc_x[k] * ksf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, isg0_100, isf_66, isf_78, \
                         isg1_100, ksd0_71, ksd1_71, ksf_116, ksf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * isg0_100[k]
                   - f_6 * pc_z[k] * isg1_100[k];

        t_176[k] = f_7 * isf_66[k]
                   + f_3 * pc_z[k] * ksf_116[k];

        t_177[k] = f_14 * isf_78[k]
                   + f_4 * ksd0_71[k]
                   - f_5 * ksd1_71[k]
                   + f_3 * pc_y[k] * ksf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, isf_69, isf_79, isf_120, \
                         ksd0_71, ksd0_72, ksd1_71, ksd1_72, ksf_119, \
                         ksf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * isf_79[k]
                   + f_3 * pc_y[k] * ksf_119[k];

        t_179[k] = f_7 * isf_69[k]
                   + f_1 * ksd0_71[k]
                   - f_2 * ksd1_71[k]
                   + f_3 * pc_z[k] * ksf_119[k];

        t_180[k] = f_14 * isf_120[k]
                   + f_1 * ksd0_72[k]
                   - f_2 * ksd1_72[k]
                   + f_3 * pc_x[k] * ksf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, isf_70, isf_80, isf_82, \
                         isf_123, ksd0_75, ksd1_75, ksf_120, ksf_122, \
                         ksf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * isf_80[k]
                   + f_3 * pc_y[k] * ksf_120[k];

        t_182[k] = f_8 * isf_70[k]
                   + f_3 * pc_z[k] * ksf_120[k];

        t_183[k] = f_14 * isf_123[k]
                   + f_4 * ksd0_75[k]
                   - f_5 * ksd1_75[k]
                   + f_3 * pc_x[k] * ksf_123[k];

        t_184[k] = f_8 * isf_82[k]
                   + f_3 * pc_y[k] * ksf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, isf_125, isf_126, isf_127, isf_128, \
                         ksd0_77, ksd1_77, ksf_125, ksf_126, ksf_127, \
                         ksf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_14 * isf_125[k]
                   + f_4 * ksd0_77[k]
                   - f_5 * ksd1_77[k]
                   + f_3 * pc_x[k] * ksf_125[k];

        t_186[k] = f_14 * isf_126[k]
                   + f_3 * pc_x[k] * ksf_126[k];

        t_187[k] = f_14 * isf_127[k]
                   + f_3 * pc_x[k] * ksf_127[k];

        t_188[k] = f_14 * isf_128[k]
                   + f_3 * pc_x[k] * ksf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, isf_76, isf_86, isf_129, \
                         ksd0_75, ksd1_75, ksf_126, ksf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_14 * isf_129[k]
                   + f_3 * pc_x[k] * ksf_129[k];

        t_190[k] = f_8 * isf_86[k]
                   + f_1 * ksd0_75[k]
                   - f_2 * ksd1_75[k]
                   + f_3 * pc_y[k] * ksf_126[k];

        t_191[k] = f_8 * isf_76[k]
                   + f_3 * pc_z[k] * ksf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, isg0_135, isf_79, \
                         isf_88, isf_89, isg1_135, ksd0_77, ksd1_77, ksf_128, \
                         ksf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * isf_88[k]
                   + f_4 * ksd0_77[k]
                   - f_5 * ksd1_77[k]
                   + f_3 * pc_y[k] * ksf_128[k];

        t_193[k] = f_8 * isf_89[k]
                   + f_3 * pc_y[k] * ksf_129[k];

        t_194[k] = f_8 * isf_79[k]
                   + f_1 * ksd0_77[k]
                   - f_2 * ksd1_77[k]
                   + f_3 * pc_z[k] * ksf_129[k];

        t_195[k] = pa_y[k] * isg0_135[k]
                   - f_6 * pc_y[k] * isg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, isg0_138, isf_80, \
                         isf_90, isf_91, isf_92, isg1_138, ksf_130, \
                         ksf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * isf_90[k]
                   + f_3 * pc_y[k] * ksf_130[k];

        t_197[k] = f_14 * isf_80[k]
                   + f_3 * pc_z[k] * ksf_130[k];

        t_198[k] = pa_y[k] * isg0_138[k]
                   + f_8 * isf_91[k]
                   - f_6 * pc_y[k] * isg1_138[k];

        t_199[k] = f_7 * isf_92[k]
                   + f_3 * pc_y[k] * ksf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, isg0_140, isf_136, \
                         isf_137, isf_138, isg1_140, ksf_136, ksf_137, \
                         ksf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * isg0_140[k]
                   - f_6 * pc_y[k] * isg1_140[k];

        t_201[k] = f_14 * isf_136[k]
                   + f_3 * pc_x[k] * ksf_136[k];

        t_202[k] = f_14 * isf_137[k]
                   + f_3 * pc_x[k] * ksf_137[k];

        t_203[k] = f_14 * isf_138[k]
                   + f_3 * pc_x[k] * ksf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, isf_86, isf_96, isf_139, \
                         ksd0_81, ksd1_81, ksf_136, ksf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_14 * isf_139[k]
                   + f_3 * pc_x[k] * ksf_139[k];

        t_205[k] = f_7 * isf_96[k]
                   + f_1 * ksd0_81[k]
                   - f_2 * ksd1_81[k]
                   + f_3 * pc_y[k] * ksf_136[k];

        t_206[k] = f_14 * isf_86[k]
                   + f_3 * pc_z[k] * ksf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, isg0_149, isf_98, isf_99, isg1_149, \
                         ksd0_83, ksd1_83, ksf_138, ksf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * isf_98[k]
                   + f_4 * ksd0_83[k]
                   - f_5 * ksd1_83[k]
                   + f_3 * pc_y[k] * ksf_138[k];

        t_208[k] = f_7 * isf_99[k]
                   + f_3 * pc_y[k] * ksf_139[k];

        t_209[k] = pa_y[k] * isg0_149[k]
                   - f_6 * pc_y[k] * isg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, isf_90, isf_140, \
                         ksd0_84, ksd1_84, ksf_140, ksf_141, ksf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_14 * isf_140[k]
                   + f_1 * ksd0_84[k]
                   - f_2 * ksd1_84[k]
                   + f_3 * pc_x[k] * ksf_140[k];

        t_211[k] = f_3 * pc_y[k] * ksf_140[k];

        t_212[k] = f_13 * isf_90[k]
                   + f_3 * pc_z[k] * ksf_140[k];

        t_213[k] = f_4 * ksd0_84[k]
                   - f_5 * ksd1_84[k]
                   + f_3 * pc_y[k] * ksf_141[k];

        t_214[k] = f_3 * pc_y[k] * ksf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, isf_145, isf_146, isf_147, \
                         ksd0_89, ksd1_89, ksf_145, ksf_146, ksf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_14 * isf_145[k]
                   + f_4 * ksd0_89[k]
                   - f_5 * ksd1_89[k]
                   + f_3 * pc_x[k] * ksf_145[k];

        t_216[k] = f_14 * isf_146[k]
                   + f_3 * pc_x[k] * ksf_146[k];

        t_217[k] = f_14 * isf_147[k]
                   + f_3 * pc_x[k] * ksf_147[k];

        t_218[k] = f_3 * pc_y[k] * ksf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, isf_149, ksd0_87, ksd0_88, ksd1_87, \
                         ksd1_88, ksf_146, ksf_147, ksf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_14 * isf_149[k]
                   + f_3 * pc_x[k] * ksf_149[k];

        t_220[k] = f_1 * ksd0_87[k]
                   - f_2 * ksd1_87[k]
                   + f_3 * pc_y[k] * ksf_146[k];

        t_221[k] = f_10 * ksd0_88[k]
                   - f_11 * ksd1_88[k]
                   + f_3 * pc_y[k] * ksf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, isf_99, isf_150, \
                         ksd0_89, ksd0_90, ksd1_89, ksd1_90, ksf_148, ksf_149, \
                         ksf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * ksd0_89[k]
                   - f_5 * ksd1_89[k]
                   + f_3 * pc_y[k] * ksf_148[k];

        t_223[k] = f_3 * pc_y[k] * ksf_149[k];

        t_224[k] = f_13 * isf_99[k]
                   + f_1 * ksd0_89[k]
                   - f_2 * ksd1_89[k]
                   + f_3 * pc_z[k] * ksf_149[k];

        t_225[k] = f_8 * isf_150[k]
                   + f_1 * ksd0_90[k]
                   - f_2 * ksd1_90[k]
                   + f_3 * pc_x[k] * ksf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, isf_100, isf_153, \
                         ksd0_93, ksd1_93, ksf_150, ksf_151, ksf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_12 * isf_100[k]
                   + f_3 * pc_y[k] * ksf_150[k];

        t_227[k] = f_3 * pc_z[k] * ksf_150[k];

        t_228[k] = f_8 * isf_153[k]
                   + f_4 * ksd0_93[k]
                   - f_5 * ksd1_93[k]
                   + f_3 * pc_x[k] * ksf_153[k];

        t_229[k] = f_3 * pc_z[k] * ksf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, isf_156, isf_158, ksd0_90, \
                         ksd1_90, ksf_152, ksf_153, ksf_156, ksf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * ksd0_90[k]
                   - f_5 * ksd1_90[k]
                   + f_3 * pc_z[k] * ksf_152[k];

        t_231[k] = f_8 * isf_156[k]
                   + f_3 * pc_x[k] * ksf_156[k];

        t_232[k] = f_3 * pc_z[k] * ksf_153[k];

        t_233[k] = f_8 * isf_158[k]
                   + f_3 * pc_x[k] * ksf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, isf_106, \
                         isf_109, isf_159, ksd0_93, ksd1_93, ksf_156, ksf_157, \
                         ksf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_8 * isf_159[k]
                   + f_3 * pc_x[k] * ksf_159[k];

        t_235[k] = f_12 * isf_106[k]
                   + f_1 * ksd0_93[k]
                   - f_2 * ksd1_93[k]
                   + f_3 * pc_y[k] * ksf_156[k];

        t_236[k] = f_3 * pc_z[k] * ksf_156[k];

        t_237[k] = f_4 * ksd0_93[k]
                   - f_5 * ksd1_93[k]
                   + f_3 * pc_z[k] * ksf_157[k];

        t_238[k] = f_12 * isf_109[k]
                   + f_3 * pc_y[k] * ksf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, isg0_150, isf_100, \
                         isf_110, isg1_150, ksd0_95, ksd1_95, ksf_159, \
                         ksf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * ksd0_95[k]
                   - f_2 * ksd1_95[k]
                   + f_3 * pc_z[k] * ksf_159[k];

        t_240[k] = pa_z[k] * isg0_150[k]
                   - f_6 * pc_z[k] * isg1_150[k];

        t_241[k] = f_13 * isf_110[k]
                   + f_3 * pc_y[k] * ksf_160[k];

        t_242[k] = f_7 * isf_100[k]
                   + f_3 * pc_z[k] * ksf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, isg0_153, isf_112, \
                         isf_165, isg1_153, ksd0_101, ksd1_101, ksf_162, \
                         ksf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * isg0_153[k]
                   - f_6 * pc_z[k] * isg1_153[k];

        t_244[k] = f_13 * isf_112[k]
                   + f_3 * pc_y[k] * ksf_162[k];

        t_245[k] = f_8 * isf_165[k]
                   + f_4 * ksd0_101[k]
                   - f_5 * ksd1_101[k]
                   + f_3 * pc_x[k] * ksf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, isf_166, isf_167, isf_168, isf_169, \
                         ksf_166, ksf_167, ksf_168, ksf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_8 * isf_166[k]
                   + f_3 * pc_x[k] * ksf_166[k];

        t_247[k] = f_8 * isf_167[k]
                   + f_3 * pc_x[k] * ksf_167[k];

        t_248[k] = f_8 * isf_168[k]
                   + f_3 * pc_x[k] * ksf_168[k];

        t_249[k] = f_8 * isf_169[k]
                   + f_3 * pc_x[k] * ksf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, isg0_160, isf_106, isf_118, \
                         isg1_160, ksd0_101, ksd1_101, ksf_166, \
                         ksf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * isg0_160[k]
                   - f_6 * pc_z[k] * isg1_160[k];

        t_251[k] = f_7 * isf_106[k]
                   + f_3 * pc_z[k] * ksf_166[k];

        t_252[k] = f_13 * isf_118[k]
                   + f_4 * ksd0_101[k]
                   - f_5 * ksd1_101[k]
                   + f_3 * pc_y[k] * ksf_168[k];
    }
}

static auto
compute_prim_ksg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isg0,
                                                          const size_t isf, const size_t isg1,
                                                          const size_t ksd0, const size_t ksd1,
                                                          const size_t ksf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;

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
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isg0_210 = buffer.data(isg0 + 210);
    const auto *isg0_213 = buffer.data(isg0 + 213);
    const auto *isg0_215 = buffer.data(isg0 + 215);
    const auto *isg0_224 = buffer.data(isg0 + 224);
    const auto *isg0_225 = buffer.data(isg0 + 225);
    const auto *isg0_228 = buffer.data(isg0 + 228);
    const auto *isg0_315 = buffer.data(isg0 + 315);
    const auto *isg0_318 = buffer.data(isg0 + 318);
    const auto *isg0_325 = buffer.data(isg0 + 325);
    const auto *isg0_327 = buffer.data(isg0 + 327);
    const auto *isg0_329 = buffer.data(isg0 + 329);
    const auto *isg0_335 = buffer.data(isg0 + 335);
    const auto *isg0_340 = buffer.data(isg0 + 340);
    const auto *isg0_342 = buffer.data(isg0 + 342);
    const auto *isg0_344 = buffer.data(isg0 + 344);
    const auto *isg0_345 = buffer.data(isg0 + 345);
    const auto *isg0_348 = buffer.data(isg0 + 348);
    const auto *isg0_350 = buffer.data(isg0 + 350);
    const auto *isg0_355 = buffer.data(isg0 + 355);
    const auto *isg0_357 = buffer.data(isg0 + 357);
    const auto *isg0_359 = buffer.data(isg0 + 359);
    const auto *isg0_360 = buffer.data(isg0 + 360);
    const auto *isg0_363 = buffer.data(isg0 + 363);
    const auto *isg0_365 = buffer.data(isg0 + 365);
    const auto *isg0_370 = buffer.data(isg0 + 370);
    const auto *isg0_372 = buffer.data(isg0 + 372);

    const auto *isf_109 = buffer.data(isf + 109);
    const auto *isf_110 = buffer.data(isf + 110);
    const auto *isf_116 = buffer.data(isf + 116);
    const auto *isf_119 = buffer.data(isf + 119);
    const auto *isf_120 = buffer.data(isf + 120);
    const auto *isf_122 = buffer.data(isf + 122);
    const auto *isf_126 = buffer.data(isf + 126);
    const auto *isf_128 = buffer.data(isf + 128);
    const auto *isf_129 = buffer.data(isf + 129);
    const auto *isf_130 = buffer.data(isf + 130);
    const auto *isf_132 = buffer.data(isf + 132);
    const auto *isf_136 = buffer.data(isf + 136);
    const auto *isf_138 = buffer.data(isf + 138);
    const auto *isf_139 = buffer.data(isf + 139);
    const auto *isf_140 = buffer.data(isf + 140);
    const auto *isf_141 = buffer.data(isf + 141);
    const auto *isf_142 = buffer.data(isf + 142);
    const auto *isf_146 = buffer.data(isf + 146);
    const auto *isf_148 = buffer.data(isf + 148);
    const auto *isf_149 = buffer.data(isf + 149);
    const auto *isf_150 = buffer.data(isf + 150);
    const auto *isf_156 = buffer.data(isf + 156);
    const auto *isf_159 = buffer.data(isf + 159);
    const auto *isf_160 = buffer.data(isf + 160);
    const auto *isf_162 = buffer.data(isf + 162);
    const auto *isf_166 = buffer.data(isf + 166);
    const auto *isf_169 = buffer.data(isf + 169);
    const auto *isf_170 = buffer.data(isf + 170);
    const auto *isf_172 = buffer.data(isf + 172);
    const auto *isf_173 = buffer.data(isf + 173);
    const auto *isf_175 = buffer.data(isf + 175);
    const auto *isf_176 = buffer.data(isf + 176);
    const auto *isf_177 = buffer.data(isf + 177);
    const auto *isf_178 = buffer.data(isf + 178);
    const auto *isf_179 = buffer.data(isf + 179);
    const auto *isf_180 = buffer.data(isf + 180);
    const auto *isf_182 = buffer.data(isf + 182);
    const auto *isf_183 = buffer.data(isf + 183);
    const auto *isf_185 = buffer.data(isf + 185);
    const auto *isf_186 = buffer.data(isf + 186);
    const auto *isf_187 = buffer.data(isf + 187);
    const auto *isf_188 = buffer.data(isf + 188);
    const auto *isf_189 = buffer.data(isf + 189);
    const auto *isf_196 = buffer.data(isf + 196);
    const auto *isf_197 = buffer.data(isf + 197);
    const auto *isf_198 = buffer.data(isf + 198);
    const auto *isf_199 = buffer.data(isf + 199);
    const auto *isf_200 = buffer.data(isf + 200);
    const auto *isf_205 = buffer.data(isf + 205);
    const auto *isf_206 = buffer.data(isf + 206);
    const auto *isf_207 = buffer.data(isf + 207);
    const auto *isf_209 = buffer.data(isf + 209);
    const auto *isf_210 = buffer.data(isf + 210);
    const auto *isf_213 = buffer.data(isf + 213);
    const auto *isf_216 = buffer.data(isf + 216);
    const auto *isf_218 = buffer.data(isf + 218);
    const auto *isf_219 = buffer.data(isf + 219);
    const auto *isf_225 = buffer.data(isf + 225);
    const auto *isf_226 = buffer.data(isf + 226);
    const auto *isf_227 = buffer.data(isf + 227);
    const auto *isf_228 = buffer.data(isf + 228);
    const auto *isf_229 = buffer.data(isf + 229);
    const auto *isf_230 = buffer.data(isf + 230);
    const auto *isf_233 = buffer.data(isf + 233);
    const auto *isf_235 = buffer.data(isf + 235);
    const auto *isf_236 = buffer.data(isf + 236);
    const auto *isf_237 = buffer.data(isf + 237);
    const auto *isf_238 = buffer.data(isf + 238);
    const auto *isf_239 = buffer.data(isf + 239);
    const auto *isf_240 = buffer.data(isf + 240);
    const auto *isf_243 = buffer.data(isf + 243);
    const auto *isf_245 = buffer.data(isf + 245);
    const auto *isf_246 = buffer.data(isf + 246);
    const auto *isf_247 = buffer.data(isf + 247);
    const auto *isf_248 = buffer.data(isf + 248);
    const auto *isf_249 = buffer.data(isf + 249);

    const auto *isg1_210 = buffer.data(isg1 + 210);
    const auto *isg1_213 = buffer.data(isg1 + 213);
    const auto *isg1_215 = buffer.data(isg1 + 215);
    const auto *isg1_224 = buffer.data(isg1 + 224);
    const auto *isg1_225 = buffer.data(isg1 + 225);
    const auto *isg1_228 = buffer.data(isg1 + 228);
    const auto *isg1_315 = buffer.data(isg1 + 315);
    const auto *isg1_318 = buffer.data(isg1 + 318);
    const auto *isg1_325 = buffer.data(isg1 + 325);
    const auto *isg1_327 = buffer.data(isg1 + 327);
    const auto *isg1_329 = buffer.data(isg1 + 329);
    const auto *isg1_335 = buffer.data(isg1 + 335);
    const auto *isg1_340 = buffer.data(isg1 + 340);
    const auto *isg1_342 = buffer.data(isg1 + 342);
    const auto *isg1_344 = buffer.data(isg1 + 344);
    const auto *isg1_345 = buffer.data(isg1 + 345);
    const auto *isg1_348 = buffer.data(isg1 + 348);
    const auto *isg1_350 = buffer.data(isg1 + 350);
    const auto *isg1_355 = buffer.data(isg1 + 355);
    const auto *isg1_357 = buffer.data(isg1 + 357);
    const auto *isg1_359 = buffer.data(isg1 + 359);
    const auto *isg1_360 = buffer.data(isg1 + 360);
    const auto *isg1_363 = buffer.data(isg1 + 363);
    const auto *isg1_365 = buffer.data(isg1 + 365);
    const auto *isg1_370 = buffer.data(isg1 + 370);
    const auto *isg1_372 = buffer.data(isg1 + 372);

    const auto *ksd0_101 = buffer.data(ksd0 + 101);
    const auto *ksd0_102 = buffer.data(ksd0 + 102);
    const auto *ksd0_105 = buffer.data(ksd0 + 105);
    const auto *ksd0_107 = buffer.data(ksd0 + 107);
    const auto *ksd0_108 = buffer.data(ksd0 + 108);
    const auto *ksd0_111 = buffer.data(ksd0 + 111);
    const auto *ksd0_113 = buffer.data(ksd0 + 113);
    const auto *ksd0_117 = buffer.data(ksd0 + 117);
    const auto *ksd0_119 = buffer.data(ksd0 + 119);
    const auto *ksd0_120 = buffer.data(ksd0 + 120);
    const auto *ksd0_123 = buffer.data(ksd0 + 123);
    const auto *ksd0_124 = buffer.data(ksd0 + 124);
    const auto *ksd0_125 = buffer.data(ksd0 + 125);
    const auto *ksd0_126 = buffer.data(ksd0 + 126);

    const auto *ksd1_101 = buffer.data(ksd1 + 101);
    const auto *ksd1_102 = buffer.data(ksd1 + 102);
    const auto *ksd1_105 = buffer.data(ksd1 + 105);
    const auto *ksd1_107 = buffer.data(ksd1 + 107);
    const auto *ksd1_108 = buffer.data(ksd1 + 108);
    const auto *ksd1_111 = buffer.data(ksd1 + 111);
    const auto *ksd1_113 = buffer.data(ksd1 + 113);
    const auto *ksd1_117 = buffer.data(ksd1 + 117);
    const auto *ksd1_119 = buffer.data(ksd1 + 119);
    const auto *ksd1_120 = buffer.data(ksd1 + 120);
    const auto *ksd1_123 = buffer.data(ksd1 + 123);
    const auto *ksd1_124 = buffer.data(ksd1 + 124);
    const auto *ksd1_125 = buffer.data(ksd1 + 125);
    const auto *ksd1_126 = buffer.data(ksd1 + 126);

    const auto *ksf_169 = buffer.data(ksf + 169);
    const auto *ksf_170 = buffer.data(ksf + 170);
    const auto *ksf_172 = buffer.data(ksf + 172);
    const auto *ksf_173 = buffer.data(ksf + 173);
    const auto *ksf_175 = buffer.data(ksf + 175);
    const auto *ksf_176 = buffer.data(ksf + 176);
    const auto *ksf_177 = buffer.data(ksf + 177);
    const auto *ksf_178 = buffer.data(ksf + 178);
    const auto *ksf_179 = buffer.data(ksf + 179);
    const auto *ksf_180 = buffer.data(ksf + 180);
    const auto *ksf_182 = buffer.data(ksf + 182);
    const auto *ksf_183 = buffer.data(ksf + 183);
    const auto *ksf_185 = buffer.data(ksf + 185);
    const auto *ksf_186 = buffer.data(ksf + 186);
    const auto *ksf_187 = buffer.data(ksf + 187);
    const auto *ksf_188 = buffer.data(ksf + 188);
    const auto *ksf_189 = buffer.data(ksf + 189);
    const auto *ksf_190 = buffer.data(ksf + 190);
    const auto *ksf_192 = buffer.data(ksf + 192);
    const auto *ksf_196 = buffer.data(ksf + 196);
    const auto *ksf_197 = buffer.data(ksf + 197);
    const auto *ksf_198 = buffer.data(ksf + 198);
    const auto *ksf_199 = buffer.data(ksf + 199);
    const auto *ksf_200 = buffer.data(ksf + 200);
    const auto *ksf_201 = buffer.data(ksf + 201);
    const auto *ksf_202 = buffer.data(ksf + 202);
    const auto *ksf_205 = buffer.data(ksf + 205);
    const auto *ksf_206 = buffer.data(ksf + 206);
    const auto *ksf_207 = buffer.data(ksf + 207);
    const auto *ksf_208 = buffer.data(ksf + 208);
    const auto *ksf_209 = buffer.data(ksf + 209);
    const auto *ksf_210 = buffer.data(ksf + 210);
    const auto *ksf_211 = buffer.data(ksf + 211);
    const auto *ksf_212 = buffer.data(ksf + 212);
    const auto *ksf_213 = buffer.data(ksf + 213);
    const auto *ksf_216 = buffer.data(ksf + 216);
    const auto *ksf_218 = buffer.data(ksf + 218);
    const auto *ksf_219 = buffer.data(ksf + 219);
    const auto *ksf_220 = buffer.data(ksf + 220);
    const auto *ksf_222 = buffer.data(ksf + 222);
    const auto *ksf_226 = buffer.data(ksf + 226);
    const auto *ksf_227 = buffer.data(ksf + 227);
    const auto *ksf_228 = buffer.data(ksf + 228);
    const auto *ksf_229 = buffer.data(ksf + 229);
    const auto *ksf_230 = buffer.data(ksf + 230);
    const auto *ksf_232 = buffer.data(ksf + 232);
    const auto *ksf_236 = buffer.data(ksf + 236);
    const auto *ksf_237 = buffer.data(ksf + 237);
    const auto *ksf_238 = buffer.data(ksf + 238);
    const auto *ksf_239 = buffer.data(ksf + 239);
    const auto *ksf_240 = buffer.data(ksf + 240);
    const auto *ksf_242 = buffer.data(ksf + 242);
    const auto *ksf_246 = buffer.data(ksf + 246);
    const auto *ksf_247 = buffer.data(ksf + 247);
    const auto *ksf_248 = buffer.data(ksf + 248);
    const auto *ksf_249 = buffer.data(ksf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, isf_109, isf_119, isf_170, \
                         ksd0_101, ksd0_102, ksd1_101, ksd1_102, ksf_169, \
                         ksf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_13 * isf_119[k]
                   + f_3 * pc_y[k] * ksf_169[k];

        t_254[k] = f_7 * isf_109[k]
                   + f_1 * ksd0_101[k]
                   - f_2 * ksd1_101[k]
                   + f_3 * pc_z[k] * ksf_169[k];

        t_255[k] = f_8 * isf_170[k]
                   + f_1 * ksd0_102[k]
                   - f_2 * ksd1_102[k]
                   + f_3 * pc_x[k] * ksf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, isf_110, isf_120, \
                         isf_122, isf_173, ksd0_105, ksd1_105, ksf_170, ksf_172, \
                         ksf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * isf_120[k]
                   + f_3 * pc_y[k] * ksf_170[k];

        t_257[k] = f_8 * isf_110[k]
                   + f_3 * pc_z[k] * ksf_170[k];

        t_258[k] = f_8 * isf_173[k]
                   + f_4 * ksd0_105[k]
                   - f_5 * ksd1_105[k]
                   + f_3 * pc_x[k] * ksf_173[k];

        t_259[k] = f_14 * isf_122[k]
                   + f_3 * pc_y[k] * ksf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, isf_175, isf_176, isf_177, isf_178, \
                         ksd0_107, ksd1_107, ksf_175, ksf_176, ksf_177, \
                         ksf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_8 * isf_175[k]
                   + f_4 * ksd0_107[k]
                   - f_5 * ksd1_107[k]
                   + f_3 * pc_x[k] * ksf_175[k];

        t_261[k] = f_8 * isf_176[k]
                   + f_3 * pc_x[k] * ksf_176[k];

        t_262[k] = f_8 * isf_177[k]
                   + f_3 * pc_x[k] * ksf_177[k];

        t_263[k] = f_8 * isf_178[k]
                   + f_3 * pc_x[k] * ksf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, isf_116, isf_126, isf_179, \
                         ksd0_105, ksd1_105, ksf_176, ksf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * isf_179[k]
                   + f_3 * pc_x[k] * ksf_179[k];

        t_265[k] = f_14 * isf_126[k]
                   + f_1 * ksd0_105[k]
                   - f_2 * ksd1_105[k]
                   + f_3 * pc_y[k] * ksf_176[k];

        t_266[k] = f_8 * isf_116[k]
                   + f_3 * pc_z[k] * ksf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, isf_119, isf_128, isf_129, ksd0_107, \
                         ksd1_107, ksf_178, ksf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * isf_128[k]
                   + f_4 * ksd0_107[k]
                   - f_5 * ksd1_107[k]
                   + f_3 * pc_y[k] * ksf_178[k];

        t_268[k] = f_14 * isf_129[k]
                   + f_3 * pc_y[k] * ksf_179[k];

        t_269[k] = f_8 * isf_119[k]
                   + f_1 * ksd0_107[k]
                   - f_2 * ksd1_107[k]
                   + f_3 * pc_z[k] * ksf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, isf_120, isf_130, isf_180, \
                         ksd0_108, ksd1_108, ksf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_8 * isf_180[k]
                   + f_1 * ksd0_108[k]
                   - f_2 * ksd1_108[k]
                   + f_3 * pc_x[k] * ksf_180[k];

        t_271[k] = f_8 * isf_130[k]
                   + f_3 * pc_y[k] * ksf_180[k];

        t_272[k] = f_14 * isf_120[k]
                   + f_3 * pc_z[k] * ksf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, isf_132, isf_183, isf_185, ksd0_111, \
                         ksd0_113, ksd1_111, ksd1_113, ksf_182, ksf_183, \
                         ksf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_8 * isf_183[k]
                   + f_4 * ksd0_111[k]
                   - f_5 * ksd1_111[k]
                   + f_3 * pc_x[k] * ksf_183[k];

        t_274[k] = f_8 * isf_132[k]
                   + f_3 * pc_y[k] * ksf_182[k];

        t_275[k] = f_8 * isf_185[k]
                   + f_4 * ksd0_113[k]
                   - f_5 * ksd1_113[k]
                   + f_3 * pc_x[k] * ksf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, isf_186, isf_187, isf_188, isf_189, \
                         ksf_186, ksf_187, ksf_188, ksf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * isf_186[k]
                   + f_3 * pc_x[k] * ksf_186[k];

        t_277[k] = f_8 * isf_187[k]
                   + f_3 * pc_x[k] * ksf_187[k];

        t_278[k] = f_8 * isf_188[k]
                   + f_3 * pc_x[k] * ksf_188[k];

        t_279[k] = f_8 * isf_189[k]
                   + f_3 * pc_x[k] * ksf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, isf_126, isf_136, isf_138, ksd0_111, \
                         ksd0_113, ksd1_111, ksd1_113, ksf_186, \
                         ksf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * isf_136[k]
                   + f_1 * ksd0_111[k]
                   - f_2 * ksd1_111[k]
                   + f_3 * pc_y[k] * ksf_186[k];

        t_281[k] = f_14 * isf_126[k]
                   + f_3 * pc_z[k] * ksf_186[k];

        t_282[k] = f_8 * isf_138[k]
                   + f_4 * ksd0_113[k]
                   - f_5 * ksd1_113[k]
                   + f_3 * pc_y[k] * ksf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, isg0_210, isf_129, \
                         isf_139, isf_140, isg1_210, ksd0_113, ksd1_113, ksf_189, \
                         ksf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * isf_139[k]
                   + f_3 * pc_y[k] * ksf_189[k];

        t_284[k] = f_14 * isf_129[k]
                   + f_1 * ksd0_113[k]
                   - f_2 * ksd1_113[k]
                   + f_3 * pc_z[k] * ksf_189[k];

        t_285[k] = pa_y[k] * isg0_210[k]
                   - f_6 * pc_y[k] * isg1_210[k];

        t_286[k] = f_7 * isf_140[k]
                   + f_3 * pc_y[k] * ksf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, isg0_213, isg0_215, \
                         isf_130, isf_141, isf_142, isg1_213, isg1_215, ksf_190, \
                         ksf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_13 * isf_130[k]
                   + f_3 * pc_z[k] * ksf_190[k];

        t_288[k] = pa_y[k] * isg0_213[k]
                   + f_8 * isf_141[k]
                   - f_6 * pc_y[k] * isg1_213[k];

        t_289[k] = f_7 * isf_142[k]
                   + f_3 * pc_y[k] * ksf_192[k];

        t_290[k] = pa_y[k] * isg0_215[k]
                   - f_6 * pc_y[k] * isg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, isf_196, isf_197, isf_198, isf_199, \
                         ksf_196, ksf_197, ksf_198, ksf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_8 * isf_196[k]
                   + f_3 * pc_x[k] * ksf_196[k];

        t_292[k] = f_8 * isf_197[k]
                   + f_3 * pc_x[k] * ksf_197[k];

        t_293[k] = f_8 * isf_198[k]
                   + f_3 * pc_x[k] * ksf_198[k];

        t_294[k] = f_8 * isf_199[k]
                   + f_3 * pc_x[k] * ksf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, isf_136, isf_146, isf_148, ksd0_117, \
                         ksd0_119, ksd1_117, ksd1_119, ksf_196, \
                         ksf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * isf_146[k]
                   + f_1 * ksd0_117[k]
                   - f_2 * ksd1_117[k]
                   + f_3 * pc_y[k] * ksf_196[k];

        t_296[k] = f_13 * isf_136[k]
                   + f_3 * pc_z[k] * ksf_196[k];

        t_297[k] = f_7 * isf_148[k]
                   + f_4 * ksd0_119[k]
                   - f_5 * ksd1_119[k]
                   + f_3 * pc_y[k] * ksf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, isg0_224, isf_149, \
                         isf_200, isg1_224, ksd0_120, ksd1_120, ksf_199, \
                         ksf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * isf_149[k]
                   + f_3 * pc_y[k] * ksf_199[k];

        t_299[k] = pa_y[k] * isg0_224[k]
                   - f_6 * pc_y[k] * isg1_224[k];

        t_300[k] = f_8 * isf_200[k]
                   + f_1 * ksd0_120[k]
                   - f_2 * ksd1_120[k]
                   + f_3 * pc_x[k] * ksf_200[k];

        t_301[k] = f_3 * pc_y[k] * ksf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, isf_140, ksd0_120, ksd1_120, \
                         ksf_200, ksf_201, ksf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_12 * isf_140[k]
                   + f_3 * pc_z[k] * ksf_200[k];

        t_303[k] = f_4 * ksd0_120[k]
                   - f_5 * ksd1_120[k]
                   + f_3 * pc_y[k] * ksf_201[k];

        t_304[k] = f_3 * pc_y[k] * ksf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, isf_205, isf_206, isf_207, \
                         ksd0_125, ksd1_125, ksf_205, ksf_206, \
                         ksf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * isf_205[k]
                   + f_4 * ksd0_125[k]
                   - f_5 * ksd1_125[k]
                   + f_3 * pc_x[k] * ksf_205[k];

        t_306[k] = f_8 * isf_206[k]
                   + f_3 * pc_x[k] * ksf_206[k];

        t_307[k] = f_8 * isf_207[k]
                   + f_3 * pc_x[k] * ksf_207[k];

        t_308[k] = f_3 * pc_y[k] * ksf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, isf_209, ksd0_123, ksd0_124, \
                         ksd1_123, ksd1_124, ksf_206, ksf_207, \
                         ksf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_8 * isf_209[k]
                   + f_3 * pc_x[k] * ksf_209[k];

        t_310[k] = f_1 * ksd0_123[k]
                   - f_2 * ksd1_123[k]
                   + f_3 * pc_y[k] * ksf_206[k];

        t_311[k] = f_10 * ksd0_124[k]
                   - f_11 * ksd1_124[k]
                   + f_3 * pc_y[k] * ksf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pa_x, pc_x, pc_y, pc_z, isg0_315, \
                         isf_149, isf_210, isg1_315, ksd0_125, ksd1_125, ksf_208, \
                         ksf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * ksd0_125[k]
                   - f_5 * ksd1_125[k]
                   + f_3 * pc_y[k] * ksf_208[k];

        t_313[k] = f_3 * pc_y[k] * ksf_209[k];

        t_314[k] = f_12 * isf_149[k]
                   + f_1 * ksd0_125[k]
                   - f_2 * ksd1_125[k]
                   + f_3 * pc_z[k] * ksf_209[k];

        t_315[k] = pa_x[k] * isg0_315[k]
                   + f_13 * isf_210[k]
                   - f_6 * pc_x[k] * isg1_315[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_x, pc_x, pc_y, pc_z, isg0_318, \
                         isf_150, isf_213, isg1_318, ksf_210, ksf_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_9 * isf_150[k]
                   + f_3 * pc_y[k] * ksf_210[k];

        t_317[k] = f_3 * pc_z[k] * ksf_210[k];

        t_318[k] = pa_x[k] * isg0_318[k]
                   + f_8 * isf_213[k]
                   - f_6 * pc_x[k] * isg1_318[k];

        t_319[k] = f_3 * pc_z[k] * ksf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, isf_216, isf_218, ksd0_126, \
                         ksd1_126, ksf_212, ksf_213, ksf_216, ksf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * ksd0_126[k]
                   - f_5 * ksd1_126[k]
                   + f_3 * pc_z[k] * ksf_212[k];

        t_321[k] = f_7 * isf_216[k]
                   + f_3 * pc_x[k] * ksf_216[k];

        t_322[k] = f_3 * pc_z[k] * ksf_213[k];

        t_323[k] = f_7 * isf_218[k]
                   + f_3 * pc_x[k] * ksf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_x, pc_x, pc_z, isg0_325, isg0_327, \
                         isf_219, isg1_325, isg1_327, ksf_216, \
                         ksf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_7 * isf_219[k]
                   + f_3 * pc_x[k] * ksf_219[k];

        t_325[k] = pa_x[k] * isg0_325[k]
                   - f_6 * pc_x[k] * isg1_325[k];

        t_326[k] = f_3 * pc_z[k] * ksf_216[k];

        t_327[k] = pa_x[k] * isg0_327[k]
                   - f_6 * pc_x[k] * isg1_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_x, pa_z, pc_x, pc_y, pc_z, isg0_225, \
                         isg0_329, isf_159, isg1_225, isg1_329, \
                         ksf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_9 * isf_159[k]
                   + f_3 * pc_y[k] * ksf_219[k];

        t_329[k] = pa_x[k] * isg0_329[k]
                   - f_6 * pc_x[k] * isg1_329[k];

        t_330[k] = pa_z[k] * isg0_225[k]
                   - f_6 * pc_z[k] * isg1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pa_z, pc_y, pc_z, isg0_228, isf_150, \
                         isf_160, isf_162, isg1_228, ksf_220, ksf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_12 * isf_160[k]
                   + f_3 * pc_y[k] * ksf_220[k];

        t_332[k] = f_7 * isf_150[k]
                   + f_3 * pc_z[k] * ksf_220[k];

        t_333[k] = pa_z[k] * isg0_228[k]
                   - f_6 * pc_z[k] * isg1_228[k];

        t_334[k] = f_12 * isf_162[k]
                   + f_3 * pc_y[k] * ksf_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pa_x, pc_x, isg0_335, isf_225, isf_226, \
                         isf_227, isf_228, isg1_335, ksf_226, ksf_227, \
                         ksf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pa_x[k] * isg0_335[k]
                   + f_8 * isf_225[k]
                   - f_6 * pc_x[k] * isg1_335[k];

        t_336[k] = f_7 * isf_226[k]
                   + f_3 * pc_x[k] * ksf_226[k];

        t_337[k] = f_7 * isf_227[k]
                   + f_3 * pc_x[k] * ksf_227[k];

        t_338[k] = f_7 * isf_228[k]
                   + f_3 * pc_x[k] * ksf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pa_x, pc_x, pc_z, isg0_340, isg0_342, \
                         isf_156, isf_229, isg1_340, isg1_342, ksf_226, \
                         ksf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_7 * isf_229[k]
                   + f_3 * pc_x[k] * ksf_229[k];

        t_340[k] = pa_x[k] * isg0_340[k]
                   - f_6 * pc_x[k] * isg1_340[k];

        t_341[k] = f_7 * isf_156[k]
                   + f_3 * pc_z[k] * ksf_226[k];

        t_342[k] = pa_x[k] * isg0_342[k]
                   - f_6 * pc_x[k] * isg1_342[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pa_x, pc_x, pc_y, isg0_344, isg0_345, \
                         isf_169, isf_170, isf_230, isg1_344, isg1_345, ksf_229, \
                         ksf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_12 * isf_169[k]
                   + f_3 * pc_y[k] * ksf_229[k];

        t_344[k] = pa_x[k] * isg0_344[k]
                   - f_6 * pc_x[k] * isg1_344[k];

        t_345[k] = pa_x[k] * isg0_345[k]
                   + f_13 * isf_230[k]
                   - f_6 * pc_x[k] * isg1_345[k];

        t_346[k] = f_13 * isf_170[k]
                   + f_3 * pc_y[k] * ksf_230[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pa_x, pc_x, pc_y, pc_z, isg0_348, isf_160, \
                         isf_172, isf_233, isg1_348, ksf_230, ksf_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_8 * isf_160[k]
                   + f_3 * pc_z[k] * ksf_230[k];

        t_348[k] = pa_x[k] * isg0_348[k]
                   + f_8 * isf_233[k]
                   - f_6 * pc_x[k] * isg1_348[k];

        t_349[k] = f_13 * isf_172[k]
                   + f_3 * pc_y[k] * ksf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_x, pc_x, isg0_350, isf_235, isf_236, \
                         isf_237, isf_238, isg1_350, ksf_236, ksf_237, \
                         ksf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pa_x[k] * isg0_350[k]
                   + f_8 * isf_235[k]
                   - f_6 * pc_x[k] * isg1_350[k];

        t_351[k] = f_7 * isf_236[k]
                   + f_3 * pc_x[k] * ksf_236[k];

        t_352[k] = f_7 * isf_237[k]
                   + f_3 * pc_x[k] * ksf_237[k];

        t_353[k] = f_7 * isf_238[k]
                   + f_3 * pc_x[k] * ksf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_x, pc_x, pc_z, isg0_355, isg0_357, \
                         isf_166, isf_239, isg1_355, isg1_357, ksf_236, \
                         ksf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_7 * isf_239[k]
                   + f_3 * pc_x[k] * ksf_239[k];

        t_355[k] = pa_x[k] * isg0_355[k]
                   - f_6 * pc_x[k] * isg1_355[k];

        t_356[k] = f_8 * isf_166[k]
                   + f_3 * pc_z[k] * ksf_236[k];

        t_357[k] = pa_x[k] * isg0_357[k]
                   - f_6 * pc_x[k] * isg1_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pa_x, pc_x, pc_y, isg0_359, isg0_360, \
                         isf_179, isf_180, isf_240, isg1_359, isg1_360, ksf_239, \
                         ksf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_13 * isf_179[k]
                   + f_3 * pc_y[k] * ksf_239[k];

        t_359[k] = pa_x[k] * isg0_359[k]
                   - f_6 * pc_x[k] * isg1_359[k];

        t_360[k] = pa_x[k] * isg0_360[k]
                   + f_13 * isf_240[k]
                   - f_6 * pc_x[k] * isg1_360[k];

        t_361[k] = f_14 * isf_180[k]
                   + f_3 * pc_y[k] * ksf_240[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pa_x, pc_x, pc_y, pc_z, isg0_363, isf_170, \
                         isf_182, isf_243, isg1_363, ksf_240, ksf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_14 * isf_170[k]
                   + f_3 * pc_z[k] * ksf_240[k];

        t_363[k] = pa_x[k] * isg0_363[k]
                   + f_8 * isf_243[k]
                   - f_6 * pc_x[k] * isg1_363[k];

        t_364[k] = f_14 * isf_182[k]
                   + f_3 * pc_y[k] * ksf_242[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pa_x, pc_x, isg0_365, isf_245, isf_246, \
                         isf_247, isf_248, isg1_365, ksf_246, ksf_247, \
                         ksf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = pa_x[k] * isg0_365[k]
                   + f_8 * isf_245[k]
                   - f_6 * pc_x[k] * isg1_365[k];

        t_366[k] = f_7 * isf_246[k]
                   + f_3 * pc_x[k] * ksf_246[k];

        t_367[k] = f_7 * isf_247[k]
                   + f_3 * pc_x[k] * ksf_247[k];

        t_368[k] = f_7 * isf_248[k]
                   + f_3 * pc_x[k] * ksf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pa_x, pc_x, pc_z, isg0_370, isg0_372, \
                         isf_176, isf_249, isg1_370, isg1_372, ksf_246, \
                         ksf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_7 * isf_249[k]
                   + f_3 * pc_x[k] * ksf_249[k];

        t_370[k] = pa_x[k] * isg0_370[k]
                   - f_6 * pc_x[k] * isg1_370[k];

        t_371[k] = f_14 * isf_176[k]
                   + f_3 * pc_z[k] * ksf_246[k];

        t_372[k] = pa_x[k] * isg0_372[k]
                   - f_6 * pc_x[k] * isg1_372[k];
    }
}

static auto
compute_prim_ksg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isg0,
                                                          const size_t isf, const size_t isg1,
                                                          const size_t ksd0, const size_t ksd1,
                                                          const size_t ksf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 1.5 / q;

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
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isg0_300 = buffer.data(isg0 + 300);
    const auto *isg0_305 = buffer.data(isg0 + 305);
    const auto *isg0_315 = buffer.data(isg0 + 315);
    const auto *isg0_316 = buffer.data(isg0 + 316);
    const auto *isg0_318 = buffer.data(isg0 + 318);
    const auto *isg0_325 = buffer.data(isg0 + 325);
    const auto *isg0_327 = buffer.data(isg0 + 327);
    const auto *isg0_374 = buffer.data(isg0 + 374);
    const auto *isg0_375 = buffer.data(isg0 + 375);
    const auto *isg0_378 = buffer.data(isg0 + 378);
    const auto *isg0_380 = buffer.data(isg0 + 380);
    const auto *isg0_385 = buffer.data(isg0 + 385);
    const auto *isg0_387 = buffer.data(isg0 + 387);
    const auto *isg0_389 = buffer.data(isg0 + 389);
    const auto *isg0_393 = buffer.data(isg0 + 393);
    const auto *isg0_400 = buffer.data(isg0 + 400);
    const auto *isg0_402 = buffer.data(isg0 + 402);
    const auto *isg0_404 = buffer.data(isg0 + 404);
    const auto *isg0_405 = buffer.data(isg0 + 405);
    const auto *isg0_410 = buffer.data(isg0 + 410);
    const auto *isg0_415 = buffer.data(isg0 + 415);
    const auto *isg0_416 = buffer.data(isg0 + 416);
    const auto *isg0_417 = buffer.data(isg0 + 417);
    const auto *isg0_419 = buffer.data(isg0 + 419);

    const auto *isf_180 = buffer.data(isf + 180);
    const auto *isf_186 = buffer.data(isf + 186);
    const auto *isf_189 = buffer.data(isf + 189);
    const auto *isf_190 = buffer.data(isf + 190);
    const auto *isf_192 = buffer.data(isf + 192);
    const auto *isf_196 = buffer.data(isf + 196);
    const auto *isf_199 = buffer.data(isf + 199);
    const auto *isf_200 = buffer.data(isf + 200);
    const auto *isf_202 = buffer.data(isf + 202);
    const auto *isf_209 = buffer.data(isf + 209);
    const auto *isf_216 = buffer.data(isf + 216);
    const auto *isf_217 = buffer.data(isf + 217);
    const auto *isf_219 = buffer.data(isf + 219);
    const auto *isf_226 = buffer.data(isf + 226);
    const auto *isf_229 = buffer.data(isf + 229);
    const auto *isf_236 = buffer.data(isf + 236);
    const auto *isf_238 = buffer.data(isf + 238);
    const auto *isf_239 = buffer.data(isf + 239);
    const auto *isf_246 = buffer.data(isf + 246);
    const auto *isf_248 = buffer.data(isf + 248);
    const auto *isf_249 = buffer.data(isf + 249);
    const auto *isf_250 = buffer.data(isf + 250);
    const auto *isf_253 = buffer.data(isf + 253);
    const auto *isf_255 = buffer.data(isf + 255);
    const auto *isf_256 = buffer.data(isf + 256);
    const auto *isf_257 = buffer.data(isf + 257);
    const auto *isf_258 = buffer.data(isf + 258);
    const auto *isf_259 = buffer.data(isf + 259);
    const auto *isf_263 = buffer.data(isf + 263);
    const auto *isf_266 = buffer.data(isf + 266);
    const auto *isf_267 = buffer.data(isf + 267);
    const auto *isf_268 = buffer.data(isf + 268);
    const auto *isf_269 = buffer.data(isf + 269);
    const auto *isf_270 = buffer.data(isf + 270);
    const auto *isf_275 = buffer.data(isf + 275);
    const auto *isf_276 = buffer.data(isf + 276);
    const auto *isf_277 = buffer.data(isf + 277);
    const auto *isf_279 = buffer.data(isf + 279);

    const auto *isg1_300 = buffer.data(isg1 + 300);
    const auto *isg1_305 = buffer.data(isg1 + 305);
    const auto *isg1_315 = buffer.data(isg1 + 315);
    const auto *isg1_316 = buffer.data(isg1 + 316);
    const auto *isg1_318 = buffer.data(isg1 + 318);
    const auto *isg1_325 = buffer.data(isg1 + 325);
    const auto *isg1_327 = buffer.data(isg1 + 327);
    const auto *isg1_374 = buffer.data(isg1 + 374);
    const auto *isg1_375 = buffer.data(isg1 + 375);
    const auto *isg1_378 = buffer.data(isg1 + 378);
    const auto *isg1_380 = buffer.data(isg1 + 380);
    const auto *isg1_385 = buffer.data(isg1 + 385);
    const auto *isg1_387 = buffer.data(isg1 + 387);
    const auto *isg1_389 = buffer.data(isg1 + 389);
    const auto *isg1_393 = buffer.data(isg1 + 393);
    const auto *isg1_400 = buffer.data(isg1 + 400);
    const auto *isg1_402 = buffer.data(isg1 + 402);
    const auto *isg1_404 = buffer.data(isg1 + 404);
    const auto *isg1_405 = buffer.data(isg1 + 405);
    const auto *isg1_410 = buffer.data(isg1 + 410);
    const auto *isg1_415 = buffer.data(isg1 + 415);
    const auto *isg1_416 = buffer.data(isg1 + 416);
    const auto *isg1_417 = buffer.data(isg1 + 417);
    const auto *isg1_419 = buffer.data(isg1 + 419);

    const auto *ksd0_162 = buffer.data(ksd0 + 162);
    const auto *ksd0_168 = buffer.data(ksd0 + 168);
    const auto *ksd0_169 = buffer.data(ksd0 + 169);
    const auto *ksd0_171 = buffer.data(ksd0 + 171);
    const auto *ksd0_173 = buffer.data(ksd0 + 173);
    const auto *ksd0_176 = buffer.data(ksd0 + 176);
    const auto *ksd0_178 = buffer.data(ksd0 + 178);
    const auto *ksd0_179 = buffer.data(ksd0 + 179);
    const auto *ksd0_180 = buffer.data(ksd0 + 180);
    const auto *ksd0_181 = buffer.data(ksd0 + 181);
    const auto *ksd0_182 = buffer.data(ksd0 + 182);
    const auto *ksd0_183 = buffer.data(ksd0 + 183);
    const auto *ksd0_184 = buffer.data(ksd0 + 184);
    const auto *ksd0_185 = buffer.data(ksd0 + 185);
    const auto *ksd0_186 = buffer.data(ksd0 + 186);
    const auto *ksd0_187 = buffer.data(ksd0 + 187);
    const auto *ksd0_188 = buffer.data(ksd0 + 188);
    const auto *ksd0_189 = buffer.data(ksd0 + 189);
    const auto *ksd0_190 = buffer.data(ksd0 + 190);
    const auto *ksd0_191 = buffer.data(ksd0 + 191);
    const auto *ksd0_192 = buffer.data(ksd0 + 192);
    const auto *ksd0_193 = buffer.data(ksd0 + 193);
    const auto *ksd0_194 = buffer.data(ksd0 + 194);
    const auto *ksd0_195 = buffer.data(ksd0 + 195);
    const auto *ksd0_196 = buffer.data(ksd0 + 196);
    const auto *ksd0_197 = buffer.data(ksd0 + 197);
    const auto *ksd0_198 = buffer.data(ksd0 + 198);
    const auto *ksd0_199 = buffer.data(ksd0 + 199);
    const auto *ksd0_200 = buffer.data(ksd0 + 200);

    const auto *ksd1_162 = buffer.data(ksd1 + 162);
    const auto *ksd1_168 = buffer.data(ksd1 + 168);
    const auto *ksd1_169 = buffer.data(ksd1 + 169);
    const auto *ksd1_171 = buffer.data(ksd1 + 171);
    const auto *ksd1_173 = buffer.data(ksd1 + 173);
    const auto *ksd1_176 = buffer.data(ksd1 + 176);
    const auto *ksd1_178 = buffer.data(ksd1 + 178);
    const auto *ksd1_179 = buffer.data(ksd1 + 179);
    const auto *ksd1_180 = buffer.data(ksd1 + 180);
    const auto *ksd1_181 = buffer.data(ksd1 + 181);
    const auto *ksd1_182 = buffer.data(ksd1 + 182);
    const auto *ksd1_183 = buffer.data(ksd1 + 183);
    const auto *ksd1_184 = buffer.data(ksd1 + 184);
    const auto *ksd1_185 = buffer.data(ksd1 + 185);
    const auto *ksd1_186 = buffer.data(ksd1 + 186);
    const auto *ksd1_187 = buffer.data(ksd1 + 187);
    const auto *ksd1_188 = buffer.data(ksd1 + 188);
    const auto *ksd1_189 = buffer.data(ksd1 + 189);
    const auto *ksd1_190 = buffer.data(ksd1 + 190);
    const auto *ksd1_191 = buffer.data(ksd1 + 191);
    const auto *ksd1_192 = buffer.data(ksd1 + 192);
    const auto *ksd1_193 = buffer.data(ksd1 + 193);
    const auto *ksd1_194 = buffer.data(ksd1 + 194);
    const auto *ksd1_195 = buffer.data(ksd1 + 195);
    const auto *ksd1_196 = buffer.data(ksd1 + 196);
    const auto *ksd1_197 = buffer.data(ksd1 + 197);
    const auto *ksd1_198 = buffer.data(ksd1 + 198);
    const auto *ksd1_199 = buffer.data(ksd1 + 199);
    const auto *ksd1_200 = buffer.data(ksd1 + 200);

    const auto *ksf_249 = buffer.data(ksf + 249);
    const auto *ksf_250 = buffer.data(ksf + 250);
    const auto *ksf_252 = buffer.data(ksf + 252);
    const auto *ksf_256 = buffer.data(ksf + 256);
    const auto *ksf_257 = buffer.data(ksf + 257);
    const auto *ksf_258 = buffer.data(ksf + 258);
    const auto *ksf_259 = buffer.data(ksf + 259);
    const auto *ksf_260 = buffer.data(ksf + 260);
    const auto *ksf_262 = buffer.data(ksf + 262);
    const auto *ksf_266 = buffer.data(ksf + 266);
    const auto *ksf_267 = buffer.data(ksf + 267);
    const auto *ksf_268 = buffer.data(ksf + 268);
    const auto *ksf_269 = buffer.data(ksf + 269);
    const auto *ksf_270 = buffer.data(ksf + 270);
    const auto *ksf_271 = buffer.data(ksf + 271);
    const auto *ksf_272 = buffer.data(ksf + 272);
    const auto *ksf_275 = buffer.data(ksf + 275);
    const auto *ksf_276 = buffer.data(ksf + 276);
    const auto *ksf_277 = buffer.data(ksf + 277);
    const auto *ksf_279 = buffer.data(ksf + 279);
    const auto *ksf_280 = buffer.data(ksf + 280);
    const auto *ksf_281 = buffer.data(ksf + 281);
    const auto *ksf_283 = buffer.data(ksf + 283);
    const auto *ksf_285 = buffer.data(ksf + 285);
    const auto *ksf_286 = buffer.data(ksf + 286);
    const auto *ksf_287 = buffer.data(ksf + 287);
    const auto *ksf_288 = buffer.data(ksf + 288);
    const auto *ksf_289 = buffer.data(ksf + 289);
    const auto *ksf_292 = buffer.data(ksf + 292);
    const auto *ksf_294 = buffer.data(ksf + 294);
    const auto *ksf_295 = buffer.data(ksf + 295);
    const auto *ksf_296 = buffer.data(ksf + 296);
    const auto *ksf_297 = buffer.data(ksf + 297);
    const auto *ksf_298 = buffer.data(ksf + 298);
    const auto *ksf_299 = buffer.data(ksf + 299);
    const auto *ksf_300 = buffer.data(ksf + 300);
    const auto *ksf_301 = buffer.data(ksf + 301);
    const auto *ksf_302 = buffer.data(ksf + 302);
    const auto *ksf_303 = buffer.data(ksf + 303);
    const auto *ksf_304 = buffer.data(ksf + 304);
    const auto *ksf_305 = buffer.data(ksf + 305);
    const auto *ksf_306 = buffer.data(ksf + 306);
    const auto *ksf_307 = buffer.data(ksf + 307);
    const auto *ksf_308 = buffer.data(ksf + 308);
    const auto *ksf_309 = buffer.data(ksf + 309);
    const auto *ksf_310 = buffer.data(ksf + 310);
    const auto *ksf_311 = buffer.data(ksf + 311);
    const auto *ksf_312 = buffer.data(ksf + 312);
    const auto *ksf_313 = buffer.data(ksf + 313);
    const auto *ksf_314 = buffer.data(ksf + 314);
    const auto *ksf_315 = buffer.data(ksf + 315);
    const auto *ksf_316 = buffer.data(ksf + 316);
    const auto *ksf_317 = buffer.data(ksf + 317);
    const auto *ksf_318 = buffer.data(ksf + 318);
    const auto *ksf_319 = buffer.data(ksf + 319);
    const auto *ksf_320 = buffer.data(ksf + 320);
    const auto *ksf_321 = buffer.data(ksf + 321);
    const auto *ksf_322 = buffer.data(ksf + 322);
    const auto *ksf_323 = buffer.data(ksf + 323);
    const auto *ksf_324 = buffer.data(ksf + 324);
    const auto *ksf_325 = buffer.data(ksf + 325);
    const auto *ksf_326 = buffer.data(ksf + 326);
    const auto *ksf_327 = buffer.data(ksf + 327);
    const auto *ksf_328 = buffer.data(ksf + 328);
    const auto *ksf_329 = buffer.data(ksf + 329);
    const auto *ksf_330 = buffer.data(ksf + 330);
    const auto *ksf_331 = buffer.data(ksf + 331);
    const auto *ksf_332 = buffer.data(ksf + 332);

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pa_x, pc_x, pc_y, isg0_374, isg0_375, \
                         isf_189, isf_190, isf_250, isg1_374, isg1_375, ksf_249, \
                         ksf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * isf_189[k]
                   + f_3 * pc_y[k] * ksf_249[k];

        t_374[k] = pa_x[k] * isg0_374[k]
                   - f_6 * pc_x[k] * isg1_374[k];

        t_375[k] = pa_x[k] * isg0_375[k]
                   + f_13 * isf_250[k]
                   - f_6 * pc_x[k] * isg1_375[k];

        t_376[k] = f_8 * isf_190[k]
                   + f_3 * pc_y[k] * ksf_250[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pa_x, pc_x, pc_y, pc_z, isg0_378, isf_180, \
                         isf_192, isf_253, isg1_378, ksf_250, ksf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_13 * isf_180[k]
                   + f_3 * pc_z[k] * ksf_250[k];

        t_378[k] = pa_x[k] * isg0_378[k]
                   + f_8 * isf_253[k]
                   - f_6 * pc_x[k] * isg1_378[k];

        t_379[k] = f_8 * isf_192[k]
                   + f_3 * pc_y[k] * ksf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pa_x, pc_x, isg0_380, isf_255, isf_256, \
                         isf_257, isf_258, isg1_380, ksf_256, ksf_257, \
                         ksf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = pa_x[k] * isg0_380[k]
                   + f_8 * isf_255[k]
                   - f_6 * pc_x[k] * isg1_380[k];

        t_381[k] = f_7 * isf_256[k]
                   + f_3 * pc_x[k] * ksf_256[k];

        t_382[k] = f_7 * isf_257[k]
                   + f_3 * pc_x[k] * ksf_257[k];

        t_383[k] = f_7 * isf_258[k]
                   + f_3 * pc_x[k] * ksf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pa_x, pc_x, pc_z, isg0_385, isg0_387, \
                         isf_186, isf_259, isg1_385, isg1_387, ksf_256, \
                         ksf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_7 * isf_259[k]
                   + f_3 * pc_x[k] * ksf_259[k];

        t_385[k] = pa_x[k] * isg0_385[k]
                   - f_6 * pc_x[k] * isg1_385[k];

        t_386[k] = f_13 * isf_186[k]
                   + f_3 * pc_z[k] * ksf_256[k];

        t_387[k] = pa_x[k] * isg0_387[k]
                   - f_6 * pc_x[k] * isg1_387[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pa_x, pa_y, pc_x, pc_y, isg0_300, \
                         isg0_389, isf_199, isf_200, isg1_300, isg1_389, ksf_259, \
                         ksf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * isf_199[k]
                   + f_3 * pc_y[k] * ksf_259[k];

        t_389[k] = pa_x[k] * isg0_389[k]
                   - f_6 * pc_x[k] * isg1_389[k];

        t_390[k] = pa_y[k] * isg0_300[k]
                   - f_6 * pc_y[k] * isg1_300[k];

        t_391[k] = f_7 * isf_200[k]
                   + f_3 * pc_y[k] * ksf_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, pa_x, pc_x, pc_y, pc_z, isg0_393, isf_190, \
                         isf_202, isf_263, isg1_393, ksf_260, ksf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_12 * isf_190[k]
                   + f_3 * pc_z[k] * ksf_260[k];

        t_393[k] = pa_x[k] * isg0_393[k]
                   + f_8 * isf_263[k]
                   - f_6 * pc_x[k] * isg1_393[k];

        t_394[k] = f_7 * isf_202[k]
                   + f_3 * pc_y[k] * ksf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, isg0_305, isf_266, \
                         isf_267, isf_268, isg1_305, ksf_266, ksf_267, \
                         ksf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * isg0_305[k]
                   - f_6 * pc_y[k] * isg1_305[k];

        t_396[k] = f_7 * isf_266[k]
                   + f_3 * pc_x[k] * ksf_266[k];

        t_397[k] = f_7 * isf_267[k]
                   + f_3 * pc_x[k] * ksf_267[k];

        t_398[k] = f_7 * isf_268[k]
                   + f_3 * pc_x[k] * ksf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pa_x, pc_x, pc_z, isg0_400, isg0_402, \
                         isf_196, isf_269, isg1_400, isg1_402, ksf_266, \
                         ksf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_7 * isf_269[k]
                   + f_3 * pc_x[k] * ksf_269[k];

        t_400[k] = pa_x[k] * isg0_400[k]
                   - f_6 * pc_x[k] * isg1_400[k];

        t_401[k] = f_12 * isf_196[k]
                   + f_3 * pc_z[k] * ksf_266[k];

        t_402[k] = pa_x[k] * isg0_402[k]
                   - f_6 * pc_x[k] * isg1_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pa_x, pc_x, pc_y, isg0_404, isg0_405, \
                         isf_209, isf_270, isg1_404, isg1_405, ksf_269, \
                         ksf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * isf_209[k]
                   + f_3 * pc_y[k] * ksf_269[k];

        t_404[k] = pa_x[k] * isg0_404[k]
                   - f_6 * pc_x[k] * isg1_404[k];

        t_405[k] = pa_x[k] * isg0_405[k]
                   + f_13 * isf_270[k]
                   - f_6 * pc_x[k] * isg1_405[k];

        t_406[k] = f_3 * pc_y[k] * ksf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_y, pc_z, isf_200, ksd0_162, ksd1_162, \
                         ksf_270, ksf_271, ksf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_9 * isf_200[k]
                   + f_3 * pc_z[k] * ksf_270[k];

        t_408[k] = f_4 * ksd0_162[k]
                   - f_5 * ksd1_162[k]
                   + f_3 * pc_y[k] * ksf_271[k];

        t_409[k] = f_3 * pc_y[k] * ksf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pc_x, pc_y, isg0_410, isf_275, \
                         isf_276, isf_277, isg1_410, ksf_275, ksf_276, \
                         ksf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pa_x[k] * isg0_410[k]
                   + f_8 * isf_275[k]
                   - f_6 * pc_x[k] * isg1_410[k];

        t_411[k] = f_7 * isf_276[k]
                   + f_3 * pc_x[k] * ksf_276[k];

        t_412[k] = f_7 * isf_277[k]
                   + f_3 * pc_x[k] * ksf_277[k];

        t_413[k] = f_3 * pc_y[k] * ksf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, pa_x, pc_x, pc_y, isg0_415, \
                         isg0_416, isg0_417, isf_279, isg1_415, isg1_416, isg1_417, \
                         ksf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_7 * isf_279[k]
                   + f_3 * pc_x[k] * ksf_279[k];

        t_415[k] = pa_x[k] * isg0_415[k]
                   - f_6 * pc_x[k] * isg1_415[k];

        t_416[k] = pa_x[k] * isg0_416[k]
                   - f_6 * pc_x[k] * isg1_416[k];

        t_417[k] = pa_x[k] * isg0_417[k]
                   - f_6 * pc_x[k] * isg1_417[k];

        t_418[k] = f_3 * pc_y[k] * ksf_279[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pa_x, pc_x, pc_z, isg0_419, isg1_419, \
                         ksd0_168, ksd0_169, ksd1_168, ksd1_169, ksf_280, \
                         ksf_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = pa_x[k] * isg0_419[k]
                   - f_6 * pc_x[k] * isg1_419[k];

        t_420[k] = f_1 * ksd0_168[k]
                   - f_2 * ksd1_168[k]
                   + f_3 * pc_x[k] * ksf_280[k];

        t_421[k] = f_10 * ksd0_169[k]
                   - f_11 * ksd1_169[k]
                   + f_3 * pc_x[k] * ksf_281[k];

        t_422[k] = f_3 * pc_z[k] * ksf_280[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, t_427, pc_x, pc_z, ksd0_171, ksd0_173, \
                         ksd1_171, ksd1_173, ksf_281, ksf_283, ksf_285, ksf_286, \
                         ksf_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * ksd0_171[k]
                   - f_5 * ksd1_171[k]
                   + f_3 * pc_x[k] * ksf_283[k];

        t_424[k] = f_3 * pc_z[k] * ksf_281[k];

        t_425[k] = f_4 * ksd0_173[k]
                   - f_5 * ksd1_173[k]
                   + f_3 * pc_x[k] * ksf_285[k];

        t_426[k] = f_3 * pc_x[k] * ksf_286[k];

        t_427[k] = f_3 * pc_x[k] * ksf_287[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, t_432, pc_x, pc_y, pc_z, isf_216, \
                         ksd0_171, ksd1_171, ksf_286, ksf_287, ksf_288, \
                         ksf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_3 * pc_x[k] * ksf_288[k];

        t_429[k] = f_3 * pc_x[k] * ksf_289[k];

        t_430[k] = f_0 * isf_216[k]
                   + f_1 * ksd0_171[k]
                   - f_2 * ksd1_171[k]
                   + f_3 * pc_y[k] * ksf_286[k];

        t_431[k] = f_3 * pc_z[k] * ksf_286[k];

        t_432[k] = f_4 * ksd0_171[k]
                   - f_5 * ksd1_171[k]
                   + f_3 * pc_z[k] * ksf_287[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pa_z, pc_y, pc_z, isg0_315, isg0_316, \
                         isf_219, isg1_315, isg1_316, ksd0_173, ksd1_173, \
                         ksf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_0 * isf_219[k]
                   + f_3 * pc_y[k] * ksf_289[k];

        t_434[k] = f_1 * ksd0_173[k]
                   - f_2 * ksd1_173[k]
                   + f_3 * pc_z[k] * ksf_289[k];

        t_435[k] = pa_z[k] * isg0_315[k]
                   - f_6 * pc_z[k] * isg1_315[k];

        t_436[k] = pa_z[k] * isg0_316[k]
                   - f_6 * pc_z[k] * isg1_316[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_z, pc_x, pc_z, isg0_318, isg1_318, ksd0_176, \
                         ksd0_178, ksd1_176, ksd1_178, ksf_292, \
                         ksf_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_10 * ksd0_176[k]
                   - f_11 * ksd1_176[k]
                   + f_3 * pc_x[k] * ksf_292[k];

        t_438[k] = pa_z[k] * isg0_318[k]
                   - f_6 * pc_z[k] * isg1_318[k];

        t_439[k] = f_4 * ksd0_178[k]
                   - f_5 * ksd1_178[k]
                   + f_3 * pc_x[k] * ksf_294[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, ksd0_179, ksd1_179, ksf_295, \
                         ksf_296, ksf_297, ksf_298, ksf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_4 * ksd0_179[k]
                   - f_5 * ksd1_179[k]
                   + f_3 * pc_x[k] * ksf_295[k];

        t_441[k] = f_3 * pc_x[k] * ksf_296[k];

        t_442[k] = f_3 * pc_x[k] * ksf_297[k];

        t_443[k] = f_3 * pc_x[k] * ksf_298[k];

        t_444[k] = f_3 * pc_x[k] * ksf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, isg0_325, isg0_327, \
                         isf_216, isf_217, isf_229, isg1_325, isg1_327, ksf_296, \
                         ksf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * isg0_325[k]
                   - f_6 * pc_z[k] * isg1_325[k];

        t_446[k] = f_7 * isf_216[k]
                   + f_3 * pc_z[k] * ksf_296[k];

        t_447[k] = pa_z[k] * isg0_327[k]
                   + f_8 * isf_217[k]
                   - f_6 * pc_z[k] * isg1_327[k];

        t_448[k] = f_9 * isf_229[k]
                   + f_3 * pc_y[k] * ksf_299[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, pc_x, pc_z, isf_219, ksd0_179, ksd0_180, \
                         ksd0_181, ksd1_179, ksd1_180, ksd1_181, ksf_299, ksf_300, \
                         ksf_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_7 * isf_219[k]
                   + f_1 * ksd0_179[k]
                   - f_2 * ksd1_179[k]
                   + f_3 * pc_z[k] * ksf_299[k];

        t_450[k] = f_1 * ksd0_180[k]
                   - f_2 * ksd1_180[k]
                   + f_3 * pc_x[k] * ksf_300[k];

        t_451[k] = f_10 * ksd0_181[k]
                   - f_11 * ksd1_181[k]
                   + f_3 * pc_x[k] * ksf_301[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pc_x, ksd0_182, ksd0_183, ksd0_184, ksd1_182, \
                         ksd1_183, ksd1_184, ksf_302, ksf_303, \
                         ksf_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_10 * ksd0_182[k]
                   - f_11 * ksd1_182[k]
                   + f_3 * pc_x[k] * ksf_302[k];

        t_453[k] = f_4 * ksd0_183[k]
                   - f_5 * ksd1_183[k]
                   + f_3 * pc_x[k] * ksf_303[k];

        t_454[k] = f_4 * ksd0_184[k]
                   - f_5 * ksd1_184[k]
                   + f_3 * pc_x[k] * ksf_304[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, t_459, pc_x, ksd0_185, ksd1_185, ksf_305, \
                         ksf_306, ksf_307, ksf_308, ksf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_4 * ksd0_185[k]
                   - f_5 * ksd1_185[k]
                   + f_3 * pc_x[k] * ksf_305[k];

        t_456[k] = f_3 * pc_x[k] * ksf_306[k];

        t_457[k] = f_3 * pc_x[k] * ksf_307[k];

        t_458[k] = f_3 * pc_x[k] * ksf_308[k];

        t_459[k] = f_3 * pc_x[k] * ksf_309[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pc_y, pc_z, isf_226, isf_236, isf_238, ksd0_183, \
                         ksd0_185, ksd1_183, ksd1_185, ksf_306, \
                         ksf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_12 * isf_236[k]
                   + f_1 * ksd0_183[k]
                   - f_2 * ksd1_183[k]
                   + f_3 * pc_y[k] * ksf_306[k];

        t_461[k] = f_8 * isf_226[k]
                   + f_3 * pc_z[k] * ksf_306[k];

        t_462[k] = f_12 * isf_238[k]
                   + f_4 * ksd0_185[k]
                   - f_5 * ksd1_185[k]
                   + f_3 * pc_y[k] * ksf_308[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, pc_x, pc_y, pc_z, isf_229, isf_239, ksd0_185, \
                         ksd0_186, ksd1_185, ksd1_186, ksf_309, \
                         ksf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_12 * isf_239[k]
                   + f_3 * pc_y[k] * ksf_309[k];

        t_464[k] = f_8 * isf_229[k]
                   + f_1 * ksd0_185[k]
                   - f_2 * ksd1_185[k]
                   + f_3 * pc_z[k] * ksf_309[k];

        t_465[k] = f_1 * ksd0_186[k]
                   - f_2 * ksd1_186[k]
                   + f_3 * pc_x[k] * ksf_310[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, pc_x, ksd0_187, ksd0_188, ksd0_189, ksd1_187, \
                         ksd1_188, ksd1_189, ksf_311, ksf_312, \
                         ksf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_10 * ksd0_187[k]
                   - f_11 * ksd1_187[k]
                   + f_3 * pc_x[k] * ksf_311[k];

        t_467[k] = f_10 * ksd0_188[k]
                   - f_11 * ksd1_188[k]
                   + f_3 * pc_x[k] * ksf_312[k];

        t_468[k] = f_4 * ksd0_189[k]
                   - f_5 * ksd1_189[k]
                   + f_3 * pc_x[k] * ksf_313[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, pc_x, ksd0_190, ksd0_191, \
                         ksd1_190, ksd1_191, ksf_314, ksf_315, ksf_316, ksf_317, \
                         ksf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_4 * ksd0_190[k]
                   - f_5 * ksd1_190[k]
                   + f_3 * pc_x[k] * ksf_314[k];

        t_470[k] = f_4 * ksd0_191[k]
                   - f_5 * ksd1_191[k]
                   + f_3 * pc_x[k] * ksf_315[k];

        t_471[k] = f_3 * pc_x[k] * ksf_316[k];

        t_472[k] = f_3 * pc_x[k] * ksf_317[k];

        t_473[k] = f_3 * pc_x[k] * ksf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, isf_236, isf_246, ksd0_189, \
                         ksd1_189, ksf_316, ksf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_3 * pc_x[k] * ksf_319[k];

        t_475[k] = f_13 * isf_246[k]
                   + f_1 * ksd0_189[k]
                   - f_2 * ksd1_189[k]
                   + f_3 * pc_y[k] * ksf_316[k];

        t_476[k] = f_14 * isf_236[k]
                   + f_3 * pc_z[k] * ksf_316[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, isf_239, isf_248, isf_249, ksd0_191, \
                         ksd1_191, ksf_318, ksf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_13 * isf_248[k]
                   + f_4 * ksd0_191[k]
                   - f_5 * ksd1_191[k]
                   + f_3 * pc_y[k] * ksf_318[k];

        t_478[k] = f_13 * isf_249[k]
                   + f_3 * pc_y[k] * ksf_319[k];

        t_479[k] = f_14 * isf_239[k]
                   + f_1 * ksd0_191[k]
                   - f_2 * ksd1_191[k]
                   + f_3 * pc_z[k] * ksf_319[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pc_x, ksd0_192, ksd0_193, ksd0_194, ksd1_192, \
                         ksd1_193, ksd1_194, ksf_320, ksf_321, \
                         ksf_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_1 * ksd0_192[k]
                   - f_2 * ksd1_192[k]
                   + f_3 * pc_x[k] * ksf_320[k];

        t_481[k] = f_10 * ksd0_193[k]
                   - f_11 * ksd1_193[k]
                   + f_3 * pc_x[k] * ksf_321[k];

        t_482[k] = f_10 * ksd0_194[k]
                   - f_11 * ksd1_194[k]
                   + f_3 * pc_x[k] * ksf_322[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pc_x, ksd0_195, ksd0_196, ksd0_197, \
                         ksd1_195, ksd1_196, ksd1_197, ksf_323, ksf_324, ksf_325, \
                         ksf_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_4 * ksd0_195[k]
                   - f_5 * ksd1_195[k]
                   + f_3 * pc_x[k] * ksf_323[k];

        t_484[k] = f_4 * ksd0_196[k]
                   - f_5 * ksd1_196[k]
                   + f_3 * pc_x[k] * ksf_324[k];

        t_485[k] = f_4 * ksd0_197[k]
                   - f_5 * ksd1_197[k]
                   + f_3 * pc_x[k] * ksf_325[k];

        t_486[k] = f_3 * pc_x[k] * ksf_326[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, pc_x, pc_y, pc_z, isf_246, \
                         isf_256, ksd0_195, ksd1_195, ksf_326, ksf_327, ksf_328, \
                         ksf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_3 * pc_x[k] * ksf_327[k];

        t_488[k] = f_3 * pc_x[k] * ksf_328[k];

        t_489[k] = f_3 * pc_x[k] * ksf_329[k];

        t_490[k] = f_14 * isf_256[k]
                   + f_1 * ksd0_195[k]
                   - f_2 * ksd1_195[k]
                   + f_3 * pc_y[k] * ksf_326[k];

        t_491[k] = f_13 * isf_246[k]
                   + f_3 * pc_z[k] * ksf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, isf_249, isf_258, isf_259, ksd0_197, \
                         ksd1_197, ksf_328, ksf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * isf_258[k]
                   + f_4 * ksd0_197[k]
                   - f_5 * ksd1_197[k]
                   + f_3 * pc_y[k] * ksf_328[k];

        t_493[k] = f_14 * isf_259[k]
                   + f_3 * pc_y[k] * ksf_329[k];

        t_494[k] = f_13 * isf_249[k]
                   + f_1 * ksd0_197[k]
                   - f_2 * ksd1_197[k]
                   + f_3 * pc_z[k] * ksf_329[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_x, ksd0_198, ksd0_199, ksd0_200, ksd1_198, \
                         ksd1_199, ksd1_200, ksf_330, ksf_331, \
                         ksf_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_1 * ksd0_198[k]
                   - f_2 * ksd1_198[k]
                   + f_3 * pc_x[k] * ksf_330[k];

        t_496[k] = f_10 * ksd0_199[k]
                   - f_11 * ksd1_199[k]
                   + f_3 * pc_x[k] * ksf_331[k];

        t_497[k] = f_10 * ksd0_200[k]
                   - f_11 * ksd1_200[k]
                   + f_3 * pc_x[k] * ksf_332[k];
    }
}

static auto
compute_prim_ksg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isg0,
                                                          const size_t isf, const size_t isg1,
                                                          const size_t ksd0, const size_t ksd1,
                                                          const size_t ksf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 2.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isg0_405 = buffer.data(isg0 + 405);
    const auto *isg0_407 = buffer.data(isg0 + 407);
    const auto *isg0_410 = buffer.data(isg0 + 410);
    const auto *isg0_415 = buffer.data(isg0 + 415);
    const auto *isg0_417 = buffer.data(isg0 + 417);
    const auto *isg0_419 = buffer.data(isg0 + 419);

    const auto *isf_256 = buffer.data(isf + 256);
    const auto *isf_259 = buffer.data(isf + 259);
    const auto *isf_266 = buffer.data(isf + 266);
    const auto *isf_268 = buffer.data(isf + 268);
    const auto *isf_269 = buffer.data(isf + 269);
    const auto *isf_276 = buffer.data(isf + 276);
    const auto *isf_278 = buffer.data(isf + 278);
    const auto *isf_279 = buffer.data(isf + 279);

    const auto *isg1_405 = buffer.data(isg1 + 405);
    const auto *isg1_407 = buffer.data(isg1 + 407);
    const auto *isg1_410 = buffer.data(isg1 + 410);
    const auto *isg1_415 = buffer.data(isg1 + 415);
    const auto *isg1_417 = buffer.data(isg1 + 417);
    const auto *isg1_419 = buffer.data(isg1 + 419);

    const auto *ksd0_201 = buffer.data(ksd0 + 201);
    const auto *ksd0_202 = buffer.data(ksd0 + 202);
    const auto *ksd0_203 = buffer.data(ksd0 + 203);
    const auto *ksd0_205 = buffer.data(ksd0 + 205);
    const auto *ksd0_207 = buffer.data(ksd0 + 207);
    const auto *ksd0_208 = buffer.data(ksd0 + 208);
    const auto *ksd0_210 = buffer.data(ksd0 + 210);
    const auto *ksd0_212 = buffer.data(ksd0 + 212);
    const auto *ksd0_213 = buffer.data(ksd0 + 213);
    const auto *ksd0_214 = buffer.data(ksd0 + 214);
    const auto *ksd0_215 = buffer.data(ksd0 + 215);

    const auto *ksd1_201 = buffer.data(ksd1 + 201);
    const auto *ksd1_202 = buffer.data(ksd1 + 202);
    const auto *ksd1_203 = buffer.data(ksd1 + 203);
    const auto *ksd1_205 = buffer.data(ksd1 + 205);
    const auto *ksd1_207 = buffer.data(ksd1 + 207);
    const auto *ksd1_208 = buffer.data(ksd1 + 208);
    const auto *ksd1_210 = buffer.data(ksd1 + 210);
    const auto *ksd1_212 = buffer.data(ksd1 + 212);
    const auto *ksd1_213 = buffer.data(ksd1 + 213);
    const auto *ksd1_214 = buffer.data(ksd1 + 214);
    const auto *ksd1_215 = buffer.data(ksd1 + 215);

    const auto *ksf_333 = buffer.data(ksf + 333);
    const auto *ksf_334 = buffer.data(ksf + 334);
    const auto *ksf_335 = buffer.data(ksf + 335);
    const auto *ksf_336 = buffer.data(ksf + 336);
    const auto *ksf_337 = buffer.data(ksf + 337);
    const auto *ksf_338 = buffer.data(ksf + 338);
    const auto *ksf_339 = buffer.data(ksf + 339);
    const auto *ksf_341 = buffer.data(ksf + 341);
    const auto *ksf_343 = buffer.data(ksf + 343);
    const auto *ksf_344 = buffer.data(ksf + 344);
    const auto *ksf_346 = buffer.data(ksf + 346);
    const auto *ksf_347 = buffer.data(ksf + 347);
    const auto *ksf_348 = buffer.data(ksf + 348);
    const auto *ksf_349 = buffer.data(ksf + 349);
    const auto *ksf_350 = buffer.data(ksf + 350);
    const auto *ksf_352 = buffer.data(ksf + 352);
    const auto *ksf_353 = buffer.data(ksf + 353);
    const auto *ksf_355 = buffer.data(ksf + 355);
    const auto *ksf_356 = buffer.data(ksf + 356);
    const auto *ksf_357 = buffer.data(ksf + 357);
    const auto *ksf_358 = buffer.data(ksf + 358);
    const auto *ksf_359 = buffer.data(ksf + 359);

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pc_x, ksd0_201, ksd0_202, ksd0_203, \
                         ksd1_201, ksd1_202, ksd1_203, ksf_333, ksf_334, ksf_335, \
                         ksf_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_4 * ksd0_201[k]
                   - f_5 * ksd1_201[k]
                   + f_3 * pc_x[k] * ksf_333[k];

        t_499[k] = f_4 * ksd0_202[k]
                   - f_5 * ksd1_202[k]
                   + f_3 * pc_x[k] * ksf_334[k];

        t_500[k] = f_4 * ksd0_203[k]
                   - f_5 * ksd1_203[k]
                   + f_3 * pc_x[k] * ksf_335[k];

        t_501[k] = f_3 * pc_x[k] * ksf_336[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, t_506, pc_x, pc_y, pc_z, isf_256, \
                         isf_266, ksd0_201, ksd1_201, ksf_336, ksf_337, ksf_338, \
                         ksf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_3 * pc_x[k] * ksf_337[k];

        t_503[k] = f_3 * pc_x[k] * ksf_338[k];

        t_504[k] = f_3 * pc_x[k] * ksf_339[k];

        t_505[k] = f_8 * isf_266[k]
                   + f_1 * ksd0_201[k]
                   - f_2 * ksd1_201[k]
                   + f_3 * pc_y[k] * ksf_336[k];

        t_506[k] = f_12 * isf_256[k]
                   + f_3 * pc_z[k] * ksf_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pa_y, pc_y, pc_z, isg0_405, isf_259, \
                         isf_268, isf_269, isg1_405, ksd0_203, ksd1_203, ksf_338, \
                         ksf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * isf_268[k]
                   + f_4 * ksd0_203[k]
                   - f_5 * ksd1_203[k]
                   + f_3 * pc_y[k] * ksf_338[k];

        t_508[k] = f_8 * isf_269[k]
                   + f_3 * pc_y[k] * ksf_339[k];

        t_509[k] = f_12 * isf_259[k]
                   + f_1 * ksd0_203[k]
                   - f_2 * ksd1_203[k]
                   + f_3 * pc_z[k] * ksf_339[k];

        t_510[k] = pa_y[k] * isg0_405[k]
                   - f_6 * pc_y[k] * isg1_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pa_y, pc_x, pc_y, isg0_407, isg1_407, ksd0_205, \
                         ksd0_207, ksd1_205, ksd1_207, ksf_341, \
                         ksf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_10 * ksd0_205[k]
                   - f_11 * ksd1_205[k]
                   + f_3 * pc_x[k] * ksf_341[k];

        t_512[k] = pa_y[k] * isg0_407[k]
                   - f_6 * pc_y[k] * isg1_407[k];

        t_513[k] = f_4 * ksd0_207[k]
                   - f_5 * ksd1_207[k]
                   + f_3 * pc_x[k] * ksf_343[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, pa_y, pc_x, pc_y, isg0_410, \
                         isg1_410, ksd0_208, ksd1_208, ksf_344, ksf_346, ksf_347, \
                         ksf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_4 * ksd0_208[k]
                   - f_5 * ksd1_208[k]
                   + f_3 * pc_x[k] * ksf_344[k];

        t_515[k] = pa_y[k] * isg0_410[k]
                   - f_6 * pc_y[k] * isg1_410[k];

        t_516[k] = f_3 * pc_x[k] * ksf_346[k];

        t_517[k] = f_3 * pc_x[k] * ksf_347[k];

        t_518[k] = f_3 * pc_x[k] * ksf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pa_y, pc_x, pc_y, pc_z, isg0_415, isf_266, \
                         isf_276, isg1_415, ksf_346, ksf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_3 * pc_x[k] * ksf_349[k];

        t_520[k] = pa_y[k] * isg0_415[k]
                   + f_13 * isf_276[k]
                   - f_6 * pc_y[k] * isg1_415[k];

        t_521[k] = f_9 * isf_266[k]
                   + f_3 * pc_z[k] * ksf_346[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pa_y, pc_y, isg0_417, isg0_419, isf_278, \
                         isf_279, isg1_417, isg1_419, ksf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = pa_y[k] * isg0_417[k]
                   + f_8 * isf_278[k]
                   - f_6 * pc_y[k] * isg1_417[k];

        t_523[k] = f_7 * isf_279[k]
                   + f_3 * pc_y[k] * ksf_349[k];

        t_524[k] = pa_y[k] * isg0_419[k]
                   - f_6 * pc_y[k] * isg1_419[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, pc_y, ksd0_210, ksd0_212, \
                         ksd0_213, ksd1_210, ksd1_212, ksd1_213, ksf_350, ksf_352, \
                         ksf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_1 * ksd0_210[k]
                   - f_2 * ksd1_210[k]
                   + f_3 * pc_x[k] * ksf_350[k];

        t_526[k] = f_3 * pc_y[k] * ksf_350[k];

        t_527[k] = f_10 * ksd0_212[k]
                   - f_11 * ksd1_212[k]
                   + f_3 * pc_x[k] * ksf_352[k];

        t_528[k] = f_4 * ksd0_213[k]
                   - f_5 * ksd1_213[k]
                   + f_3 * pc_x[k] * ksf_353[k];

        t_529[k] = f_3 * pc_y[k] * ksf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, t_534, pc_x, ksd0_215, ksd1_215, ksf_355, \
                         ksf_356, ksf_357, ksf_358, ksf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_4 * ksd0_215[k]
                   - f_5 * ksd1_215[k]
                   + f_3 * pc_x[k] * ksf_355[k];

        t_531[k] = f_3 * pc_x[k] * ksf_356[k];

        t_532[k] = f_3 * pc_x[k] * ksf_357[k];

        t_533[k] = f_3 * pc_x[k] * ksf_358[k];

        t_534[k] = f_3 * pc_x[k] * ksf_359[k];
    }

#pragma omp simd aligned(t_535, t_536, t_537, t_538, pc_y, ksd0_213, ksd0_214, ksd0_215, \
                         ksd1_213, ksd1_214, ksd1_215, ksf_356, ksf_357, ksf_358, \
                         ksf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_535[k] = f_1 * ksd0_213[k]
                   - f_2 * ksd1_213[k]
                   + f_3 * pc_y[k] * ksf_356[k];

        t_536[k] = f_10 * ksd0_214[k]
                   - f_11 * ksd1_214[k]
                   + f_3 * pc_y[k] * ksf_357[k];

        t_537[k] = f_4 * ksd0_215[k]
                   - f_5 * ksd1_215[k]
                   + f_3 * pc_y[k] * ksf_358[k];

        t_538[k] = f_3 * pc_y[k] * ksf_359[k];
    }

#pragma omp simd aligned(t_539, pc_z, isf_279, ksd0_215, ksd1_215, \
                         ksf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_0 * isf_279[k]
                   + f_1 * ksd0_215[k]
                   - f_2 * ksd1_215[k]
                   + f_3 * pc_z[k] * ksf_359[k];
    }
}

auto
compute_prim_ksg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isg0, const size_t isf,
                                                   const size_t isg1, const size_t ksd0,
                                                   const size_t ksd1, const size_t ksf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isg0, isf,
                                                              isg1, ksd0, ksd1, ksf, ncols,
                                                              gamma, p, q);

    compute_prim_ksg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isg0, isf,
                                                              isg1, ksd0, ksd1, ksf, ncols,
                                                              gamma, p, q);

    compute_prim_ksg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, isg0, isf,
                                                              isg1, ksd0, ksd1, ksf, ncols,
                                                              gamma, p, q);

    compute_prim_ksg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, isg0, isf,
                                                              isg1, ksd0, ksd1, ksf, ncols,
                                                              gamma, p, q);

    compute_prim_ksg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, isg0, isf,
                                                              isg1, ksd0, ksd1, ksf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
