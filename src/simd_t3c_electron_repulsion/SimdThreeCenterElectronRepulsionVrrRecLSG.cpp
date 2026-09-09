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


#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
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

    const auto *ksg0_0 = buffer.data(ksg0 + 0);
    const auto *ksg0_3 = buffer.data(ksg0 + 3);
    const auto *ksg0_5 = buffer.data(ksg0 + 5);
    const auto *ksg0_10 = buffer.data(ksg0 + 10);
    const auto *ksg0_14 = buffer.data(ksg0 + 14);
    const auto *ksg0_18 = buffer.data(ksg0 + 18);
    const auto *ksg0_25 = buffer.data(ksg0 + 25);
    const auto *ksg0_30 = buffer.data(ksg0 + 30);
    const auto *ksg0_35 = buffer.data(ksg0 + 35);
    const auto *ksg0_44 = buffer.data(ksg0 + 44);
    const auto *ksg0_45 = buffer.data(ksg0 + 45);
    const auto *ksg0_48 = buffer.data(ksg0 + 48);
    const auto *ksg0_55 = buffer.data(ksg0 + 55);
    const auto *ksg0_75 = buffer.data(ksg0 + 75);
    const auto *ksg0_78 = buffer.data(ksg0 + 78);
    const auto *ksg0_80 = buffer.data(ksg0 + 80);

    const auto *ksf_0 = buffer.data(ksf + 0);
    const auto *ksf_1 = buffer.data(ksf + 1);
    const auto *ksf_2 = buffer.data(ksf + 2);
    const auto *ksf_6 = buffer.data(ksf + 6);
    const auto *ksf_9 = buffer.data(ksf + 9);
    const auto *ksf_10 = buffer.data(ksf + 10);
    const auto *ksf_16 = buffer.data(ksf + 16);
    const auto *ksf_18 = buffer.data(ksf + 18);
    const auto *ksf_19 = buffer.data(ksf + 19);
    const auto *ksf_20 = buffer.data(ksf + 20);
    const auto *ksf_22 = buffer.data(ksf + 22);
    const auto *ksf_26 = buffer.data(ksf + 26);
    const auto *ksf_27 = buffer.data(ksf + 27);
    const auto *ksf_28 = buffer.data(ksf + 28);
    const auto *ksf_29 = buffer.data(ksf + 29);
    const auto *ksf_30 = buffer.data(ksf + 30);
    const auto *ksf_33 = buffer.data(ksf + 33);
    const auto *ksf_36 = buffer.data(ksf + 36);
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
    const auto *ksf_59 = buffer.data(ksf + 59);
    const auto *ksf_60 = buffer.data(ksf + 60);
    const auto *ksf_63 = buffer.data(ksf + 63);
    const auto *ksf_66 = buffer.data(ksf + 66);
    const auto *ksf_68 = buffer.data(ksf + 68);
    const auto *ksf_69 = buffer.data(ksf + 69);
    const auto *ksf_75 = buffer.data(ksf + 75);
    const auto *ksf_76 = buffer.data(ksf + 76);
    const auto *ksf_77 = buffer.data(ksf + 77);
    const auto *ksf_78 = buffer.data(ksf + 78);
    const auto *ksf_79 = buffer.data(ksf + 79);
    const auto *ksf_86 = buffer.data(ksf + 86);
    const auto *ksf_87 = buffer.data(ksf + 87);
    const auto *ksf_88 = buffer.data(ksf + 88);
    const auto *ksf_89 = buffer.data(ksf + 89);

    const auto *ksg1_0 = buffer.data(ksg1 + 0);
    const auto *ksg1_3 = buffer.data(ksg1 + 3);
    const auto *ksg1_5 = buffer.data(ksg1 + 5);
    const auto *ksg1_10 = buffer.data(ksg1 + 10);
    const auto *ksg1_14 = buffer.data(ksg1 + 14);
    const auto *ksg1_18 = buffer.data(ksg1 + 18);
    const auto *ksg1_25 = buffer.data(ksg1 + 25);
    const auto *ksg1_30 = buffer.data(ksg1 + 30);
    const auto *ksg1_35 = buffer.data(ksg1 + 35);
    const auto *ksg1_44 = buffer.data(ksg1 + 44);
    const auto *ksg1_45 = buffer.data(ksg1 + 45);
    const auto *ksg1_48 = buffer.data(ksg1 + 48);
    const auto *ksg1_55 = buffer.data(ksg1 + 55);
    const auto *ksg1_75 = buffer.data(ksg1 + 75);
    const auto *ksg1_78 = buffer.data(ksg1 + 78);
    const auto *ksg1_80 = buffer.data(ksg1 + 80);

    const auto *lsd0_0 = buffer.data(lsd0 + 0);
    const auto *lsd0_3 = buffer.data(lsd0 + 3);
    const auto *lsd0_5 = buffer.data(lsd0 + 5);
    const auto *lsd0_9 = buffer.data(lsd0 + 9);
    const auto *lsd0_16 = buffer.data(lsd0 + 16);
    const auto *lsd0_17 = buffer.data(lsd0 + 17);
    const auto *lsd0_18 = buffer.data(lsd0 + 18);
    const auto *lsd0_21 = buffer.data(lsd0 + 21);
    const auto *lsd0_23 = buffer.data(lsd0 + 23);
    const auto *lsd0_29 = buffer.data(lsd0 + 29);
    const auto *lsd0_30 = buffer.data(lsd0 + 30);
    const auto *lsd0_33 = buffer.data(lsd0 + 33);
    const auto *lsd0_34 = buffer.data(lsd0 + 34);
    const auto *lsd0_35 = buffer.data(lsd0 + 35);
    const auto *lsd0_36 = buffer.data(lsd0 + 36);
    const auto *lsd0_39 = buffer.data(lsd0 + 39);
    const auto *lsd0_41 = buffer.data(lsd0 + 41);
    const auto *lsd0_47 = buffer.data(lsd0 + 47);

    const auto *lsd1_0 = buffer.data(lsd1 + 0);
    const auto *lsd1_3 = buffer.data(lsd1 + 3);
    const auto *lsd1_5 = buffer.data(lsd1 + 5);
    const auto *lsd1_9 = buffer.data(lsd1 + 9);
    const auto *lsd1_16 = buffer.data(lsd1 + 16);
    const auto *lsd1_17 = buffer.data(lsd1 + 17);
    const auto *lsd1_18 = buffer.data(lsd1 + 18);
    const auto *lsd1_21 = buffer.data(lsd1 + 21);
    const auto *lsd1_23 = buffer.data(lsd1 + 23);
    const auto *lsd1_29 = buffer.data(lsd1 + 29);
    const auto *lsd1_30 = buffer.data(lsd1 + 30);
    const auto *lsd1_33 = buffer.data(lsd1 + 33);
    const auto *lsd1_34 = buffer.data(lsd1 + 34);
    const auto *lsd1_35 = buffer.data(lsd1 + 35);
    const auto *lsd1_36 = buffer.data(lsd1 + 36);
    const auto *lsd1_39 = buffer.data(lsd1 + 39);
    const auto *lsd1_41 = buffer.data(lsd1 + 41);
    const auto *lsd1_47 = buffer.data(lsd1 + 47);

    const auto *lsf_0 = buffer.data(lsf + 0);
    const auto *lsf_1 = buffer.data(lsf + 1);
    const auto *lsf_2 = buffer.data(lsf + 2);
    const auto *lsf_3 = buffer.data(lsf + 3);
    const auto *lsf_5 = buffer.data(lsf + 5);
    const auto *lsf_6 = buffer.data(lsf + 6);
    const auto *lsf_8 = buffer.data(lsf + 8);
    const auto *lsf_9 = buffer.data(lsf + 9);
    const auto *lsf_10 = buffer.data(lsf + 10);
    const auto *lsf_11 = buffer.data(lsf + 11);
    const auto *lsf_13 = buffer.data(lsf + 13);
    const auto *lsf_16 = buffer.data(lsf + 16);
    const auto *lsf_17 = buffer.data(lsf + 17);
    const auto *lsf_18 = buffer.data(lsf + 18);
    const auto *lsf_19 = buffer.data(lsf + 19);
    const auto *lsf_20 = buffer.data(lsf + 20);
    const auto *lsf_22 = buffer.data(lsf + 22);
    const auto *lsf_25 = buffer.data(lsf + 25);
    const auto *lsf_26 = buffer.data(lsf + 26);
    const auto *lsf_27 = buffer.data(lsf + 27);
    const auto *lsf_28 = buffer.data(lsf + 28);
    const auto *lsf_29 = buffer.data(lsf + 29);
    const auto *lsf_30 = buffer.data(lsf + 30);
    const auto *lsf_31 = buffer.data(lsf + 31);
    const auto *lsf_32 = buffer.data(lsf + 32);
    const auto *lsf_33 = buffer.data(lsf + 33);
    const auto *lsf_36 = buffer.data(lsf + 36);
    const auto *lsf_37 = buffer.data(lsf + 37);
    const auto *lsf_38 = buffer.data(lsf + 38);
    const auto *lsf_39 = buffer.data(lsf + 39);
    const auto *lsf_40 = buffer.data(lsf + 40);
    const auto *lsf_42 = buffer.data(lsf + 42);
    const auto *lsf_46 = buffer.data(lsf + 46);
    const auto *lsf_47 = buffer.data(lsf + 47);
    const auto *lsf_48 = buffer.data(lsf + 48);
    const auto *lsf_49 = buffer.data(lsf + 49);
    const auto *lsf_50 = buffer.data(lsf + 50);
    const auto *lsf_51 = buffer.data(lsf + 51);
    const auto *lsf_52 = buffer.data(lsf + 52);
    const auto *lsf_55 = buffer.data(lsf + 55);
    const auto *lsf_56 = buffer.data(lsf + 56);
    const auto *lsf_57 = buffer.data(lsf + 57);
    const auto *lsf_58 = buffer.data(lsf + 58);
    const auto *lsf_59 = buffer.data(lsf + 59);
    const auto *lsf_60 = buffer.data(lsf + 60);
    const auto *lsf_61 = buffer.data(lsf + 61);
    const auto *lsf_62 = buffer.data(lsf + 62);
    const auto *lsf_63 = buffer.data(lsf + 63);
    const auto *lsf_66 = buffer.data(lsf + 66);
    const auto *lsf_67 = buffer.data(lsf + 67);
    const auto *lsf_68 = buffer.data(lsf + 68);
    const auto *lsf_69 = buffer.data(lsf + 69);
    const auto *lsf_70 = buffer.data(lsf + 70);
    const auto *lsf_72 = buffer.data(lsf + 72);
    const auto *lsf_75 = buffer.data(lsf + 75);
    const auto *lsf_76 = buffer.data(lsf + 76);
    const auto *lsf_77 = buffer.data(lsf + 77);
    const auto *lsf_78 = buffer.data(lsf + 78);
    const auto *lsf_79 = buffer.data(lsf + 79);
    const auto *lsf_80 = buffer.data(lsf + 80);
    const auto *lsf_82 = buffer.data(lsf + 82);
    const auto *lsf_86 = buffer.data(lsf + 86);
    const auto *lsf_87 = buffer.data(lsf + 87);
    const auto *lsf_88 = buffer.data(lsf + 88);
    const auto *lsf_89 = buffer.data(lsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksf_0, lsd0_0, \
                         lsd1_0, lsf_0, lsf_1, lsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksf_0[k]
                 + f_1 * lsd0_0[k]
                 - f_2 * lsd1_0[k]
                 + f_3 * pc_x[k] * lsf_0[k];

        t_1[k] = f_3 * pc_y[k] * lsf_0[k];

        t_2[k] = f_3 * pc_z[k] * lsf_0[k];

        t_3[k] = f_4 * lsd0_0[k]
                 - f_5 * lsd1_0[k]
                 + f_3 * pc_y[k] * lsf_1[k];

        t_4[k] = f_3 * pc_y[k] * lsf_2[k];

        t_5[k] = f_4 * lsd0_0[k]
                 - f_5 * lsd1_0[k]
                 + f_3 * pc_z[k] * lsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, ksf_6, ksf_9, lsd0_3, \
                         lsd1_3, lsf_3, lsf_5, lsf_6, lsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ksf_6[k]
                 + f_3 * pc_x[k] * lsf_6[k];

        t_7[k] = f_3 * pc_z[k] * lsf_3[k];

        t_8[k] = f_3 * pc_y[k] * lsf_5[k];

        t_9[k] = f_0 * ksf_9[k]
                 + f_3 * pc_x[k] * lsf_9[k];

        t_10[k] = f_1 * lsd0_3[k]
                  - f_2 * lsd1_3[k]
                  + f_3 * pc_y[k] * lsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, ksg0_0, ksg1_0, \
                         lsd0_5, lsd1_5, lsf_6, lsf_8, lsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * lsf_6[k];

        t_12[k] = f_4 * lsd0_5[k]
                  - f_5 * lsd1_5[k]
                  + f_3 * pc_y[k] * lsf_8[k];

        t_13[k] = f_3 * pc_y[k] * lsf_9[k];

        t_14[k] = f_1 * lsd0_5[k]
                  - f_2 * lsd1_5[k]
                  + f_3 * pc_z[k] * lsf_9[k];

        t_15[k] = pa_y[k] * ksg0_0[k]
                  - f_6 * pc_y[k] * ksg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, ksg0_3, ksg0_5, \
                         ksf_0, ksf_1, ksg1_3, ksg1_5, lsf_10, lsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * ksf_0[k]
                  + f_3 * pc_y[k] * lsf_10[k];

        t_17[k] = f_3 * pc_z[k] * lsf_10[k];

        t_18[k] = pa_y[k] * ksg0_3[k]
                  + f_8 * ksf_1[k]
                  - f_6 * pc_y[k] * ksg1_3[k];

        t_19[k] = f_3 * pc_z[k] * lsf_11[k];

        t_20[k] = pa_y[k] * ksg0_5[k]
                  - f_6 * pc_y[k] * ksg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, ksf_16, ksf_18, ksf_19, lsf_13, \
                         lsf_16, lsf_18, lsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * ksf_16[k]
                  + f_3 * pc_x[k] * lsf_16[k];

        t_22[k] = f_3 * pc_z[k] * lsf_13[k];

        t_23[k] = f_9 * ksf_18[k]
                  + f_3 * pc_x[k] * lsf_18[k];

        t_24[k] = f_9 * ksf_19[k]
                  + f_3 * pc_x[k] * lsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, ksf_6, ksf_9, lsd0_9, lsd1_9, \
                         lsf_16, lsf_17, lsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * ksf_6[k]
                  + f_1 * lsd0_9[k]
                  - f_2 * lsd1_9[k]
                  + f_3 * pc_y[k] * lsf_16[k];

        t_26[k] = f_3 * pc_z[k] * lsf_16[k];

        t_27[k] = f_4 * lsd0_9[k]
                  - f_5 * lsd1_9[k]
                  + f_3 * pc_z[k] * lsf_17[k];

        t_28[k] = f_7 * ksf_9[k]
                  + f_3 * pc_y[k] * lsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, ksg0_0, ksg0_14, \
                         ksf_0, ksg1_0, ksg1_14, lsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * ksg0_14[k]
                  - f_6 * pc_y[k] * ksg1_14[k];

        t_30[k] = pa_z[k] * ksg0_0[k]
                  - f_6 * pc_z[k] * ksg1_0[k];

        t_31[k] = f_3 * pc_y[k] * lsf_20[k];

        t_32[k] = f_7 * ksf_0[k]
                  + f_3 * pc_z[k] * lsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, ksg0_3, ksg0_5, \
                         ksf_2, ksf_26, ksg1_3, ksg1_5, lsf_22, \
                         lsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * ksg0_3[k]
                  - f_6 * pc_z[k] * ksg1_3[k];

        t_34[k] = f_3 * pc_y[k] * lsf_22[k];

        t_35[k] = pa_z[k] * ksg0_5[k]
                  + f_8 * ksf_2[k]
                  - f_6 * pc_z[k] * ksg1_5[k];

        t_36[k] = f_9 * ksf_26[k]
                  + f_3 * pc_x[k] * lsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, ksg0_10, ksf_27, \
                         ksf_29, ksg1_10, lsf_25, lsf_27, lsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * ksf_27[k]
                  + f_3 * pc_x[k] * lsf_27[k];

        t_38[k] = f_3 * pc_y[k] * lsf_25[k];

        t_39[k] = f_9 * ksf_29[k]
                  + f_3 * pc_x[k] * lsf_29[k];

        t_40[k] = pa_z[k] * ksg0_10[k]
                  - f_6 * pc_z[k] * ksg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, ksf_9, lsd0_16, lsd0_17, lsd1_16, \
                         lsd1_17, lsf_27, lsf_28, lsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * lsd0_16[k]
                  - f_11 * lsd1_16[k]
                  + f_3 * pc_y[k] * lsf_27[k];

        t_42[k] = f_4 * lsd0_17[k]
                  - f_5 * lsd1_17[k]
                  + f_3 * pc_y[k] * lsf_28[k];

        t_43[k] = f_3 * pc_y[k] * lsf_29[k];

        t_44[k] = f_7 * ksf_9[k]
                  + f_1 * lsd0_17[k]
                  - f_2 * lsd1_17[k]
                  + f_3 * pc_z[k] * lsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, ksf_10, ksf_30, ksf_33, \
                         lsd0_18, lsd0_21, lsd1_18, lsd1_21, lsf_30, \
                         lsf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * ksf_30[k]
                  + f_1 * lsd0_18[k]
                  - f_2 * lsd1_18[k]
                  + f_3 * pc_x[k] * lsf_30[k];

        t_46[k] = f_8 * ksf_10[k]
                  + f_3 * pc_y[k] * lsf_30[k];

        t_47[k] = f_3 * pc_z[k] * lsf_30[k];

        t_48[k] = f_12 * ksf_33[k]
                  + f_4 * lsd0_21[k]
                  - f_5 * lsd1_21[k]
                  + f_3 * pc_x[k] * lsf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, ksf_36, ksf_38, lsd0_18, \
                         lsd1_18, lsf_31, lsf_32, lsf_33, lsf_36, \
                         lsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * lsf_31[k];

        t_50[k] = f_4 * lsd0_18[k]
                  - f_5 * lsd1_18[k]
                  + f_3 * pc_z[k] * lsf_32[k];

        t_51[k] = f_12 * ksf_36[k]
                  + f_3 * pc_x[k] * lsf_36[k];

        t_52[k] = f_3 * pc_z[k] * lsf_33[k];

        t_53[k] = f_12 * ksf_38[k]
                  + f_3 * pc_x[k] * lsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, ksf_16, ksf_19, \
                         ksf_39, lsd0_21, lsd1_21, lsf_36, lsf_37, \
                         lsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * ksf_39[k]
                  + f_3 * pc_x[k] * lsf_39[k];

        t_55[k] = f_8 * ksf_16[k]
                  + f_1 * lsd0_21[k]
                  - f_2 * lsd1_21[k]
                  + f_3 * pc_y[k] * lsf_36[k];

        t_56[k] = f_3 * pc_z[k] * lsf_36[k];

        t_57[k] = f_4 * lsd0_21[k]
                  - f_5 * lsd1_21[k]
                  + f_3 * pc_z[k] * lsf_37[k];

        t_58[k] = f_8 * ksf_19[k]
                  + f_3 * pc_y[k] * lsf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, ksg0_30, ksf_10, ksf_20, \
                         ksg1_30, lsd0_23, lsd1_23, lsf_39, lsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * lsd0_23[k]
                  - f_2 * lsd1_23[k]
                  + f_3 * pc_z[k] * lsf_39[k];

        t_60[k] = pa_y[k] * ksg0_30[k]
                  - f_6 * pc_y[k] * ksg1_30[k];

        t_61[k] = f_7 * ksf_20[k]
                  + f_3 * pc_y[k] * lsf_40[k];

        t_62[k] = f_7 * ksf_10[k]
                  + f_3 * pc_z[k] * lsf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, ksg0_18, ksg0_35, ksf_22, \
                         ksg1_18, ksg1_35, lsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * ksg0_18[k]
                  - f_6 * pc_z[k] * ksg1_18[k];

        t_64[k] = f_7 * ksf_22[k]
                  + f_3 * pc_y[k] * lsf_42[k];

        t_65[k] = pa_y[k] * ksg0_35[k]
                  - f_6 * pc_y[k] * ksg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, ksf_46, ksf_47, ksf_48, ksf_49, lsf_46, \
                         lsf_47, lsf_48, lsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * ksf_46[k]
                  + f_3 * pc_x[k] * lsf_46[k];

        t_67[k] = f_12 * ksf_47[k]
                  + f_3 * pc_x[k] * lsf_47[k];

        t_68[k] = f_12 * ksf_48[k]
                  + f_3 * pc_x[k] * lsf_48[k];

        t_69[k] = f_12 * ksf_49[k]
                  + f_3 * pc_x[k] * lsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, ksg0_25, ksf_16, ksf_28, ksg1_25, \
                         lsd0_29, lsd1_29, lsf_46, lsf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * ksg0_25[k]
                  - f_6 * pc_z[k] * ksg1_25[k];

        t_71[k] = f_7 * ksf_16[k]
                  + f_3 * pc_z[k] * lsf_46[k];

        t_72[k] = f_7 * ksf_28[k]
                  + f_4 * lsd0_29[k]
                  - f_5 * lsd1_29[k]
                  + f_3 * pc_y[k] * lsf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, ksg0_44, ksf_29, ksf_50, \
                         ksg1_44, lsd0_30, lsd1_30, lsf_49, lsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * ksf_29[k]
                  + f_3 * pc_y[k] * lsf_49[k];

        t_74[k] = pa_y[k] * ksg0_44[k]
                  - f_6 * pc_y[k] * ksg1_44[k];

        t_75[k] = f_12 * ksf_50[k]
                  + f_1 * lsd0_30[k]
                  - f_2 * lsd1_30[k]
                  + f_3 * pc_x[k] * lsf_50[k];

        t_76[k] = f_3 * pc_y[k] * lsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, ksf_20, lsd0_30, lsd1_30, lsf_50, \
                         lsf_51, lsf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * ksf_20[k]
                  + f_3 * pc_z[k] * lsf_50[k];

        t_78[k] = f_4 * lsd0_30[k]
                  - f_5 * lsd1_30[k]
                  + f_3 * pc_y[k] * lsf_51[k];

        t_79[k] = f_3 * pc_y[k] * lsf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, ksf_55, ksf_56, ksf_57, lsd0_35, \
                         lsd1_35, lsf_55, lsf_56, lsf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * ksf_55[k]
                  + f_4 * lsd0_35[k]
                  - f_5 * lsd1_35[k]
                  + f_3 * pc_x[k] * lsf_55[k];

        t_81[k] = f_12 * ksf_56[k]
                  + f_3 * pc_x[k] * lsf_56[k];

        t_82[k] = f_12 * ksf_57[k]
                  + f_3 * pc_x[k] * lsf_57[k];

        t_83[k] = f_3 * pc_y[k] * lsf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, ksf_59, lsd0_33, lsd0_34, lsd1_33, \
                         lsd1_34, lsf_56, lsf_57, lsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * ksf_59[k]
                  + f_3 * pc_x[k] * lsf_59[k];

        t_85[k] = f_1 * lsd0_33[k]
                  - f_2 * lsd1_33[k]
                  + f_3 * pc_y[k] * lsf_56[k];

        t_86[k] = f_10 * lsd0_34[k]
                  - f_11 * lsd1_34[k]
                  + f_3 * pc_y[k] * lsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, ksf_29, ksf_60, lsd0_35, \
                         lsd0_36, lsd1_35, lsd1_36, lsf_58, lsf_59, \
                         lsf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * lsd0_35[k]
                  - f_5 * lsd1_35[k]
                  + f_3 * pc_y[k] * lsf_58[k];

        t_88[k] = f_3 * pc_y[k] * lsf_59[k];

        t_89[k] = f_8 * ksf_29[k]
                  + f_1 * lsd0_35[k]
                  - f_2 * lsd1_35[k]
                  + f_3 * pc_z[k] * lsf_59[k];

        t_90[k] = f_13 * ksf_60[k]
                  + f_1 * lsd0_36[k]
                  - f_2 * lsd1_36[k]
                  + f_3 * pc_x[k] * lsf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, ksf_30, ksf_63, lsd0_39, \
                         lsd1_39, lsf_60, lsf_61, lsf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * ksf_30[k]
                  + f_3 * pc_y[k] * lsf_60[k];

        t_92[k] = f_3 * pc_z[k] * lsf_60[k];

        t_93[k] = f_13 * ksf_63[k]
                  + f_4 * lsd0_39[k]
                  - f_5 * lsd1_39[k]
                  + f_3 * pc_x[k] * lsf_63[k];

        t_94[k] = f_3 * pc_z[k] * lsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, ksf_66, ksf_68, lsd0_36, lsd1_36, \
                         lsf_62, lsf_63, lsf_66, lsf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * lsd0_36[k]
                  - f_5 * lsd1_36[k]
                  + f_3 * pc_z[k] * lsf_62[k];

        t_96[k] = f_13 * ksf_66[k]
                  + f_3 * pc_x[k] * lsf_66[k];

        t_97[k] = f_3 * pc_z[k] * lsf_63[k];

        t_98[k] = f_13 * ksf_68[k]
                  + f_3 * pc_x[k] * lsf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, ksf_36, ksf_39, \
                         ksf_69, lsd0_39, lsd1_39, lsf_66, lsf_67, \
                         lsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * ksf_69[k]
                  + f_3 * pc_x[k] * lsf_69[k];

        t_100[k] = f_14 * ksf_36[k]
                   + f_1 * lsd0_39[k]
                   - f_2 * lsd1_39[k]
                   + f_3 * pc_y[k] * lsf_66[k];

        t_101[k] = f_3 * pc_z[k] * lsf_66[k];

        t_102[k] = f_4 * lsd0_39[k]
                   - f_5 * lsd1_39[k]
                   + f_3 * pc_z[k] * lsf_67[k];

        t_103[k] = f_14 * ksf_39[k]
                   + f_3 * pc_y[k] * lsf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, ksg0_45, ksf_30, \
                         ksf_40, ksg1_45, lsd0_41, lsd1_41, lsf_69, \
                         lsf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * lsd0_41[k]
                   - f_2 * lsd1_41[k]
                   + f_3 * pc_z[k] * lsf_69[k];

        t_105[k] = pa_z[k] * ksg0_45[k]
                   - f_6 * pc_z[k] * ksg1_45[k];

        t_106[k] = f_8 * ksf_40[k]
                   + f_3 * pc_y[k] * lsf_70[k];

        t_107[k] = f_7 * ksf_30[k]
                   + f_3 * pc_z[k] * lsf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, ksg0_48, ksf_42, ksf_75, \
                         ksg1_48, lsd0_47, lsd1_47, lsf_72, lsf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * ksg0_48[k]
                   - f_6 * pc_z[k] * ksg1_48[k];

        t_109[k] = f_8 * ksf_42[k]
                   + f_3 * pc_y[k] * lsf_72[k];

        t_110[k] = f_13 * ksf_75[k]
                   + f_4 * lsd0_47[k]
                   - f_5 * lsd1_47[k]
                   + f_3 * pc_x[k] * lsf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, ksf_76, ksf_77, ksf_78, ksf_79, \
                         lsf_76, lsf_77, lsf_78, lsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * ksf_76[k]
                   + f_3 * pc_x[k] * lsf_76[k];

        t_112[k] = f_13 * ksf_77[k]
                   + f_3 * pc_x[k] * lsf_77[k];

        t_113[k] = f_13 * ksf_78[k]
                   + f_3 * pc_x[k] * lsf_78[k];

        t_114[k] = f_13 * ksf_79[k]
                   + f_3 * pc_x[k] * lsf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, ksg0_55, ksf_36, ksf_48, \
                         ksg1_55, lsd0_47, lsd1_47, lsf_76, lsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * ksg0_55[k]
                   - f_6 * pc_z[k] * ksg1_55[k];

        t_116[k] = f_7 * ksf_36[k]
                   + f_3 * pc_z[k] * lsf_76[k];

        t_117[k] = f_8 * ksf_48[k]
                   + f_4 * lsd0_47[k]
                   - f_5 * lsd1_47[k]
                   + f_3 * pc_y[k] * lsf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, ksg0_75, ksf_39, \
                         ksf_49, ksf_50, ksg1_75, lsd0_47, lsd1_47, lsf_79, \
                         lsf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * ksf_49[k]
                   + f_3 * pc_y[k] * lsf_79[k];

        t_119[k] = f_7 * ksf_39[k]
                   + f_1 * lsd0_47[k]
                   - f_2 * lsd1_47[k]
                   + f_3 * pc_z[k] * lsf_79[k];

        t_120[k] = pa_y[k] * ksg0_75[k]
                   - f_6 * pc_y[k] * ksg1_75[k];

        t_121[k] = f_7 * ksf_50[k]
                   + f_3 * pc_y[k] * lsf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, ksg0_78, ksg0_80, \
                         ksf_40, ksf_51, ksf_52, ksg1_78, ksg1_80, lsf_80, \
                         lsf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * ksf_40[k]
                   + f_3 * pc_z[k] * lsf_80[k];

        t_123[k] = pa_y[k] * ksg0_78[k]
                   + f_8 * ksf_51[k]
                   - f_6 * pc_y[k] * ksg1_78[k];

        t_124[k] = f_7 * ksf_52[k]
                   + f_3 * pc_y[k] * lsf_82[k];

        t_125[k] = pa_y[k] * ksg0_80[k]
                   - f_6 * pc_y[k] * ksg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, ksf_86, ksf_87, ksf_88, ksf_89, \
                         lsf_86, lsf_87, lsf_88, lsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * ksf_86[k]
                   + f_3 * pc_x[k] * lsf_86[k];

        t_127[k] = f_13 * ksf_87[k]
                   + f_3 * pc_x[k] * lsf_87[k];

        t_128[k] = f_13 * ksf_88[k]
                   + f_3 * pc_x[k] * lsf_88[k];

        t_129[k] = f_13 * ksf_89[k]
                   + f_3 * pc_x[k] * lsf_89[k];
    }
}

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
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
    const auto f_13 = 2.5 / q;
    const auto f_14 = 1.5 / q;
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
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksg0_89 = buffer.data(ksg0 + 89);
    const auto *ksg0_90 = buffer.data(ksg0 + 90);
    const auto *ksg0_93 = buffer.data(ksg0 + 93);
    const auto *ksg0_100 = buffer.data(ksg0 + 100);
    const auto *ksg0_135 = buffer.data(ksg0 + 135);
    const auto *ksg0_138 = buffer.data(ksg0 + 138);
    const auto *ksg0_140 = buffer.data(ksg0 + 140);
    const auto *ksg0_149 = buffer.data(ksg0 + 149);
    const auto *ksg0_150 = buffer.data(ksg0 + 150);
    const auto *ksg0_153 = buffer.data(ksg0 + 153);
    const auto *ksg0_160 = buffer.data(ksg0 + 160);

    const auto *ksf_46 = buffer.data(ksf + 46);
    const auto *ksf_50 = buffer.data(ksf + 50);
    const auto *ksf_56 = buffer.data(ksf + 56);
    const auto *ksf_58 = buffer.data(ksf + 58);
    const auto *ksf_59 = buffer.data(ksf + 59);
    const auto *ksf_60 = buffer.data(ksf + 60);
    const auto *ksf_66 = buffer.data(ksf + 66);
    const auto *ksf_69 = buffer.data(ksf + 69);
    const auto *ksf_70 = buffer.data(ksf + 70);
    const auto *ksf_72 = buffer.data(ksf + 72);
    const auto *ksf_76 = buffer.data(ksf + 76);
    const auto *ksf_78 = buffer.data(ksf + 78);
    const auto *ksf_79 = buffer.data(ksf + 79);
    const auto *ksf_80 = buffer.data(ksf + 80);
    const auto *ksf_82 = buffer.data(ksf + 82);
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
    const auto *ksf_103 = buffer.data(ksf + 103);
    const auto *ksf_106 = buffer.data(ksf + 106);
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
    const auto *ksf_123 = buffer.data(ksf + 123);
    const auto *ksf_125 = buffer.data(ksf + 125);
    const auto *ksf_126 = buffer.data(ksf + 126);
    const auto *ksf_127 = buffer.data(ksf + 127);
    const auto *ksf_128 = buffer.data(ksf + 128);
    const auto *ksf_129 = buffer.data(ksf + 129);
    const auto *ksf_136 = buffer.data(ksf + 136);
    const auto *ksf_137 = buffer.data(ksf + 137);
    const auto *ksf_138 = buffer.data(ksf + 138);
    const auto *ksf_139 = buffer.data(ksf + 139);
    const auto *ksf_140 = buffer.data(ksf + 140);
    const auto *ksf_145 = buffer.data(ksf + 145);
    const auto *ksf_146 = buffer.data(ksf + 146);
    const auto *ksf_147 = buffer.data(ksf + 147);
    const auto *ksf_149 = buffer.data(ksf + 149);
    const auto *ksf_150 = buffer.data(ksf + 150);
    const auto *ksf_153 = buffer.data(ksf + 153);
    const auto *ksf_156 = buffer.data(ksf + 156);
    const auto *ksf_158 = buffer.data(ksf + 158);
    const auto *ksf_159 = buffer.data(ksf + 159);
    const auto *ksf_165 = buffer.data(ksf + 165);
    const auto *ksf_166 = buffer.data(ksf + 166);
    const auto *ksf_167 = buffer.data(ksf + 167);
    const auto *ksf_168 = buffer.data(ksf + 168);
    const auto *ksf_169 = buffer.data(ksf + 169);

    const auto *ksg1_89 = buffer.data(ksg1 + 89);
    const auto *ksg1_90 = buffer.data(ksg1 + 90);
    const auto *ksg1_93 = buffer.data(ksg1 + 93);
    const auto *ksg1_100 = buffer.data(ksg1 + 100);
    const auto *ksg1_135 = buffer.data(ksg1 + 135);
    const auto *ksg1_138 = buffer.data(ksg1 + 138);
    const auto *ksg1_140 = buffer.data(ksg1 + 140);
    const auto *ksg1_149 = buffer.data(ksg1 + 149);
    const auto *ksg1_150 = buffer.data(ksg1 + 150);
    const auto *ksg1_153 = buffer.data(ksg1 + 153);
    const auto *ksg1_160 = buffer.data(ksg1 + 160);

    const auto *lsd0_51 = buffer.data(lsd0 + 51);
    const auto *lsd0_53 = buffer.data(lsd0 + 53);
    const auto *lsd0_54 = buffer.data(lsd0 + 54);
    const auto *lsd0_57 = buffer.data(lsd0 + 57);
    const auto *lsd0_58 = buffer.data(lsd0 + 58);
    const auto *lsd0_59 = buffer.data(lsd0 + 59);
    const auto *lsd0_60 = buffer.data(lsd0 + 60);
    const auto *lsd0_63 = buffer.data(lsd0 + 63);
    const auto *lsd0_65 = buffer.data(lsd0 + 65);
    const auto *lsd0_71 = buffer.data(lsd0 + 71);
    const auto *lsd0_72 = buffer.data(lsd0 + 72);
    const auto *lsd0_75 = buffer.data(lsd0 + 75);
    const auto *lsd0_77 = buffer.data(lsd0 + 77);
    const auto *lsd0_81 = buffer.data(lsd0 + 81);
    const auto *lsd0_83 = buffer.data(lsd0 + 83);
    const auto *lsd0_84 = buffer.data(lsd0 + 84);
    const auto *lsd0_87 = buffer.data(lsd0 + 87);
    const auto *lsd0_88 = buffer.data(lsd0 + 88);
    const auto *lsd0_89 = buffer.data(lsd0 + 89);
    const auto *lsd0_90 = buffer.data(lsd0 + 90);
    const auto *lsd0_93 = buffer.data(lsd0 + 93);
    const auto *lsd0_95 = buffer.data(lsd0 + 95);
    const auto *lsd0_101 = buffer.data(lsd0 + 101);

    const auto *lsd1_51 = buffer.data(lsd1 + 51);
    const auto *lsd1_53 = buffer.data(lsd1 + 53);
    const auto *lsd1_54 = buffer.data(lsd1 + 54);
    const auto *lsd1_57 = buffer.data(lsd1 + 57);
    const auto *lsd1_58 = buffer.data(lsd1 + 58);
    const auto *lsd1_59 = buffer.data(lsd1 + 59);
    const auto *lsd1_60 = buffer.data(lsd1 + 60);
    const auto *lsd1_63 = buffer.data(lsd1 + 63);
    const auto *lsd1_65 = buffer.data(lsd1 + 65);
    const auto *lsd1_71 = buffer.data(lsd1 + 71);
    const auto *lsd1_72 = buffer.data(lsd1 + 72);
    const auto *lsd1_75 = buffer.data(lsd1 + 75);
    const auto *lsd1_77 = buffer.data(lsd1 + 77);
    const auto *lsd1_81 = buffer.data(lsd1 + 81);
    const auto *lsd1_83 = buffer.data(lsd1 + 83);
    const auto *lsd1_84 = buffer.data(lsd1 + 84);
    const auto *lsd1_87 = buffer.data(lsd1 + 87);
    const auto *lsd1_88 = buffer.data(lsd1 + 88);
    const auto *lsd1_89 = buffer.data(lsd1 + 89);
    const auto *lsd1_90 = buffer.data(lsd1 + 90);
    const auto *lsd1_93 = buffer.data(lsd1 + 93);
    const auto *lsd1_95 = buffer.data(lsd1 + 95);
    const auto *lsd1_101 = buffer.data(lsd1 + 101);

    const auto *lsf_86 = buffer.data(lsf + 86);
    const auto *lsf_88 = buffer.data(lsf + 88);
    const auto *lsf_89 = buffer.data(lsf + 89);
    const auto *lsf_90 = buffer.data(lsf + 90);
    const auto *lsf_91 = buffer.data(lsf + 91);
    const auto *lsf_92 = buffer.data(lsf + 92);
    const auto *lsf_95 = buffer.data(lsf + 95);
    const auto *lsf_96 = buffer.data(lsf + 96);
    const auto *lsf_97 = buffer.data(lsf + 97);
    const auto *lsf_98 = buffer.data(lsf + 98);
    const auto *lsf_99 = buffer.data(lsf + 99);
    const auto *lsf_100 = buffer.data(lsf + 100);
    const auto *lsf_101 = buffer.data(lsf + 101);
    const auto *lsf_102 = buffer.data(lsf + 102);
    const auto *lsf_103 = buffer.data(lsf + 103);
    const auto *lsf_106 = buffer.data(lsf + 106);
    const auto *lsf_107 = buffer.data(lsf + 107);
    const auto *lsf_108 = buffer.data(lsf + 108);
    const auto *lsf_109 = buffer.data(lsf + 109);
    const auto *lsf_110 = buffer.data(lsf + 110);
    const auto *lsf_112 = buffer.data(lsf + 112);
    const auto *lsf_115 = buffer.data(lsf + 115);
    const auto *lsf_116 = buffer.data(lsf + 116);
    const auto *lsf_117 = buffer.data(lsf + 117);
    const auto *lsf_118 = buffer.data(lsf + 118);
    const auto *lsf_119 = buffer.data(lsf + 119);
    const auto *lsf_120 = buffer.data(lsf + 120);
    const auto *lsf_122 = buffer.data(lsf + 122);
    const auto *lsf_123 = buffer.data(lsf + 123);
    const auto *lsf_125 = buffer.data(lsf + 125);
    const auto *lsf_126 = buffer.data(lsf + 126);
    const auto *lsf_127 = buffer.data(lsf + 127);
    const auto *lsf_128 = buffer.data(lsf + 128);
    const auto *lsf_129 = buffer.data(lsf + 129);
    const auto *lsf_130 = buffer.data(lsf + 130);
    const auto *lsf_132 = buffer.data(lsf + 132);
    const auto *lsf_136 = buffer.data(lsf + 136);
    const auto *lsf_137 = buffer.data(lsf + 137);
    const auto *lsf_138 = buffer.data(lsf + 138);
    const auto *lsf_139 = buffer.data(lsf + 139);
    const auto *lsf_140 = buffer.data(lsf + 140);
    const auto *lsf_141 = buffer.data(lsf + 141);
    const auto *lsf_142 = buffer.data(lsf + 142);
    const auto *lsf_145 = buffer.data(lsf + 145);
    const auto *lsf_146 = buffer.data(lsf + 146);
    const auto *lsf_147 = buffer.data(lsf + 147);
    const auto *lsf_148 = buffer.data(lsf + 148);
    const auto *lsf_149 = buffer.data(lsf + 149);
    const auto *lsf_150 = buffer.data(lsf + 150);
    const auto *lsf_151 = buffer.data(lsf + 151);
    const auto *lsf_152 = buffer.data(lsf + 152);
    const auto *lsf_153 = buffer.data(lsf + 153);
    const auto *lsf_156 = buffer.data(lsf + 156);
    const auto *lsf_157 = buffer.data(lsf + 157);
    const auto *lsf_158 = buffer.data(lsf + 158);
    const auto *lsf_159 = buffer.data(lsf + 159);
    const auto *lsf_160 = buffer.data(lsf + 160);
    const auto *lsf_162 = buffer.data(lsf + 162);
    const auto *lsf_165 = buffer.data(lsf + 165);
    const auto *lsf_166 = buffer.data(lsf + 166);
    const auto *lsf_167 = buffer.data(lsf + 167);
    const auto *lsf_168 = buffer.data(lsf + 168);
    const auto *lsf_169 = buffer.data(lsf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, ksf_46, ksf_56, ksf_58, lsd0_51, \
                         lsd0_53, lsd1_51, lsd1_53, lsf_86, lsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * ksf_56[k]
                   + f_1 * lsd0_51[k]
                   - f_2 * lsd1_51[k]
                   + f_3 * pc_y[k] * lsf_86[k];

        t_131[k] = f_8 * ksf_46[k]
                   + f_3 * pc_z[k] * lsf_86[k];

        t_132[k] = f_7 * ksf_58[k]
                   + f_4 * lsd0_53[k]
                   - f_5 * lsd1_53[k]
                   + f_3 * pc_y[k] * lsf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, ksg0_89, ksf_59, \
                         ksf_90, ksg1_89, lsd0_54, lsd1_54, lsf_89, \
                         lsf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * ksf_59[k]
                   + f_3 * pc_y[k] * lsf_89[k];

        t_134[k] = pa_y[k] * ksg0_89[k]
                   - f_6 * pc_y[k] * ksg1_89[k];

        t_135[k] = f_13 * ksf_90[k]
                   + f_1 * lsd0_54[k]
                   - f_2 * lsd1_54[k]
                   + f_3 * pc_x[k] * lsf_90[k];

        t_136[k] = f_3 * pc_y[k] * lsf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, ksf_50, lsd0_54, lsd1_54, lsf_90, \
                         lsf_91, lsf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * ksf_50[k]
                   + f_3 * pc_z[k] * lsf_90[k];

        t_138[k] = f_4 * lsd0_54[k]
                   - f_5 * lsd1_54[k]
                   + f_3 * pc_y[k] * lsf_91[k];

        t_139[k] = f_3 * pc_y[k] * lsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, ksf_95, ksf_96, ksf_97, \
                         lsd0_59, lsd1_59, lsf_95, lsf_96, lsf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * ksf_95[k]
                   + f_4 * lsd0_59[k]
                   - f_5 * lsd1_59[k]
                   + f_3 * pc_x[k] * lsf_95[k];

        t_141[k] = f_13 * ksf_96[k]
                   + f_3 * pc_x[k] * lsf_96[k];

        t_142[k] = f_13 * ksf_97[k]
                   + f_3 * pc_x[k] * lsf_97[k];

        t_143[k] = f_3 * pc_y[k] * lsf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, ksf_99, lsd0_57, lsd0_58, lsd1_57, \
                         lsd1_58, lsf_96, lsf_97, lsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * ksf_99[k]
                   + f_3 * pc_x[k] * lsf_99[k];

        t_145[k] = f_1 * lsd0_57[k]
                   - f_2 * lsd1_57[k]
                   + f_3 * pc_y[k] * lsf_96[k];

        t_146[k] = f_10 * lsd0_58[k]
                   - f_11 * lsd1_58[k]
                   + f_3 * pc_y[k] * lsf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, ksf_59, ksf_100, \
                         lsd0_59, lsd0_60, lsd1_59, lsd1_60, lsf_98, lsf_99, \
                         lsf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * lsd0_59[k]
                   - f_5 * lsd1_59[k]
                   + f_3 * pc_y[k] * lsf_98[k];

        t_148[k] = f_3 * pc_y[k] * lsf_99[k];

        t_149[k] = f_14 * ksf_59[k]
                   + f_1 * lsd0_59[k]
                   - f_2 * lsd1_59[k]
                   + f_3 * pc_z[k] * lsf_99[k];

        t_150[k] = f_15 * ksf_100[k]
                   + f_1 * lsd0_60[k]
                   - f_2 * lsd1_60[k]
                   + f_3 * pc_x[k] * lsf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, ksf_60, ksf_103, \
                         lsd0_63, lsd1_63, lsf_100, lsf_101, lsf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_15 * ksf_60[k]
                   + f_3 * pc_y[k] * lsf_100[k];

        t_152[k] = f_3 * pc_z[k] * lsf_100[k];

        t_153[k] = f_15 * ksf_103[k]
                   + f_4 * lsd0_63[k]
                   - f_5 * lsd1_63[k]
                   + f_3 * pc_x[k] * lsf_103[k];

        t_154[k] = f_3 * pc_z[k] * lsf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, ksf_106, ksf_108, lsd0_60, \
                         lsd1_60, lsf_102, lsf_103, lsf_106, lsf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * lsd0_60[k]
                   - f_5 * lsd1_60[k]
                   + f_3 * pc_z[k] * lsf_102[k];

        t_156[k] = f_15 * ksf_106[k]
                   + f_3 * pc_x[k] * lsf_106[k];

        t_157[k] = f_3 * pc_z[k] * lsf_103[k];

        t_158[k] = f_15 * ksf_108[k]
                   + f_3 * pc_x[k] * lsf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, ksf_66, ksf_69, \
                         ksf_109, lsd0_63, lsd1_63, lsf_106, lsf_107, \
                         lsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_15 * ksf_109[k]
                   + f_3 * pc_x[k] * lsf_109[k];

        t_160[k] = f_15 * ksf_66[k]
                   + f_1 * lsd0_63[k]
                   - f_2 * lsd1_63[k]
                   + f_3 * pc_y[k] * lsf_106[k];

        t_161[k] = f_3 * pc_z[k] * lsf_106[k];

        t_162[k] = f_4 * lsd0_63[k]
                   - f_5 * lsd1_63[k]
                   + f_3 * pc_z[k] * lsf_107[k];

        t_163[k] = f_15 * ksf_69[k]
                   + f_3 * pc_y[k] * lsf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, ksg0_90, ksf_60, \
                         ksf_70, ksg1_90, lsd0_65, lsd1_65, lsf_109, \
                         lsf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * lsd0_65[k]
                   - f_2 * lsd1_65[k]
                   + f_3 * pc_z[k] * lsf_109[k];

        t_165[k] = pa_z[k] * ksg0_90[k]
                   - f_6 * pc_z[k] * ksg1_90[k];

        t_166[k] = f_14 * ksf_70[k]
                   + f_3 * pc_y[k] * lsf_110[k];

        t_167[k] = f_7 * ksf_60[k]
                   + f_3 * pc_z[k] * lsf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, ksg0_93, ksf_72, \
                         ksf_115, ksg1_93, lsd0_71, lsd1_71, lsf_112, \
                         lsf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * ksg0_93[k]
                   - f_6 * pc_z[k] * ksg1_93[k];

        t_169[k] = f_14 * ksf_72[k]
                   + f_3 * pc_y[k] * lsf_112[k];

        t_170[k] = f_15 * ksf_115[k]
                   + f_4 * lsd0_71[k]
                   - f_5 * lsd1_71[k]
                   + f_3 * pc_x[k] * lsf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, ksf_116, ksf_117, ksf_118, ksf_119, \
                         lsf_116, lsf_117, lsf_118, lsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * ksf_116[k]
                   + f_3 * pc_x[k] * lsf_116[k];

        t_172[k] = f_15 * ksf_117[k]
                   + f_3 * pc_x[k] * lsf_117[k];

        t_173[k] = f_15 * ksf_118[k]
                   + f_3 * pc_x[k] * lsf_118[k];

        t_174[k] = f_15 * ksf_119[k]
                   + f_3 * pc_x[k] * lsf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, ksg0_100, ksf_66, ksf_78, \
                         ksg1_100, lsd0_71, lsd1_71, lsf_116, lsf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * ksg0_100[k]
                   - f_6 * pc_z[k] * ksg1_100[k];

        t_176[k] = f_7 * ksf_66[k]
                   + f_3 * pc_z[k] * lsf_116[k];

        t_177[k] = f_14 * ksf_78[k]
                   + f_4 * lsd0_71[k]
                   - f_5 * lsd1_71[k]
                   + f_3 * pc_y[k] * lsf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, ksf_69, ksf_79, ksf_120, \
                         lsd0_71, lsd0_72, lsd1_71, lsd1_72, lsf_119, \
                         lsf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * ksf_79[k]
                   + f_3 * pc_y[k] * lsf_119[k];

        t_179[k] = f_7 * ksf_69[k]
                   + f_1 * lsd0_71[k]
                   - f_2 * lsd1_71[k]
                   + f_3 * pc_z[k] * lsf_119[k];

        t_180[k] = f_15 * ksf_120[k]
                   + f_1 * lsd0_72[k]
                   - f_2 * lsd1_72[k]
                   + f_3 * pc_x[k] * lsf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, ksf_70, ksf_80, ksf_82, \
                         ksf_123, lsd0_75, lsd1_75, lsf_120, lsf_122, \
                         lsf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * ksf_80[k]
                   + f_3 * pc_y[k] * lsf_120[k];

        t_182[k] = f_8 * ksf_70[k]
                   + f_3 * pc_z[k] * lsf_120[k];

        t_183[k] = f_15 * ksf_123[k]
                   + f_4 * lsd0_75[k]
                   - f_5 * lsd1_75[k]
                   + f_3 * pc_x[k] * lsf_123[k];

        t_184[k] = f_8 * ksf_82[k]
                   + f_3 * pc_y[k] * lsf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, ksf_125, ksf_126, ksf_127, ksf_128, \
                         lsd0_77, lsd1_77, lsf_125, lsf_126, lsf_127, \
                         lsf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * ksf_125[k]
                   + f_4 * lsd0_77[k]
                   - f_5 * lsd1_77[k]
                   + f_3 * pc_x[k] * lsf_125[k];

        t_186[k] = f_15 * ksf_126[k]
                   + f_3 * pc_x[k] * lsf_126[k];

        t_187[k] = f_15 * ksf_127[k]
                   + f_3 * pc_x[k] * lsf_127[k];

        t_188[k] = f_15 * ksf_128[k]
                   + f_3 * pc_x[k] * lsf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, ksf_76, ksf_86, ksf_129, \
                         lsd0_75, lsd1_75, lsf_126, lsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * ksf_129[k]
                   + f_3 * pc_x[k] * lsf_129[k];

        t_190[k] = f_8 * ksf_86[k]
                   + f_1 * lsd0_75[k]
                   - f_2 * lsd1_75[k]
                   + f_3 * pc_y[k] * lsf_126[k];

        t_191[k] = f_8 * ksf_76[k]
                   + f_3 * pc_z[k] * lsf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, ksg0_135, ksf_79, \
                         ksf_88, ksf_89, ksg1_135, lsd0_77, lsd1_77, lsf_128, \
                         lsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * ksf_88[k]
                   + f_4 * lsd0_77[k]
                   - f_5 * lsd1_77[k]
                   + f_3 * pc_y[k] * lsf_128[k];

        t_193[k] = f_8 * ksf_89[k]
                   + f_3 * pc_y[k] * lsf_129[k];

        t_194[k] = f_8 * ksf_79[k]
                   + f_1 * lsd0_77[k]
                   - f_2 * lsd1_77[k]
                   + f_3 * pc_z[k] * lsf_129[k];

        t_195[k] = pa_y[k] * ksg0_135[k]
                   - f_6 * pc_y[k] * ksg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, ksg0_138, ksf_80, \
                         ksf_90, ksf_91, ksf_92, ksg1_138, lsf_130, \
                         lsf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * ksf_90[k]
                   + f_3 * pc_y[k] * lsf_130[k];

        t_197[k] = f_14 * ksf_80[k]
                   + f_3 * pc_z[k] * lsf_130[k];

        t_198[k] = pa_y[k] * ksg0_138[k]
                   + f_8 * ksf_91[k]
                   - f_6 * pc_y[k] * ksg1_138[k];

        t_199[k] = f_7 * ksf_92[k]
                   + f_3 * pc_y[k] * lsf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, ksg0_140, ksf_136, \
                         ksf_137, ksf_138, ksg1_140, lsf_136, lsf_137, \
                         lsf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * ksg0_140[k]
                   - f_6 * pc_y[k] * ksg1_140[k];

        t_201[k] = f_15 * ksf_136[k]
                   + f_3 * pc_x[k] * lsf_136[k];

        t_202[k] = f_15 * ksf_137[k]
                   + f_3 * pc_x[k] * lsf_137[k];

        t_203[k] = f_15 * ksf_138[k]
                   + f_3 * pc_x[k] * lsf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, ksf_86, ksf_96, ksf_139, \
                         lsd0_81, lsd1_81, lsf_136, lsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_15 * ksf_139[k]
                   + f_3 * pc_x[k] * lsf_139[k];

        t_205[k] = f_7 * ksf_96[k]
                   + f_1 * lsd0_81[k]
                   - f_2 * lsd1_81[k]
                   + f_3 * pc_y[k] * lsf_136[k];

        t_206[k] = f_14 * ksf_86[k]
                   + f_3 * pc_z[k] * lsf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, ksg0_149, ksf_98, ksf_99, ksg1_149, \
                         lsd0_83, lsd1_83, lsf_138, lsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * ksf_98[k]
                   + f_4 * lsd0_83[k]
                   - f_5 * lsd1_83[k]
                   + f_3 * pc_y[k] * lsf_138[k];

        t_208[k] = f_7 * ksf_99[k]
                   + f_3 * pc_y[k] * lsf_139[k];

        t_209[k] = pa_y[k] * ksg0_149[k]
                   - f_6 * pc_y[k] * ksg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, ksf_90, ksf_140, \
                         lsd0_84, lsd1_84, lsf_140, lsf_141, lsf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * ksf_140[k]
                   + f_1 * lsd0_84[k]
                   - f_2 * lsd1_84[k]
                   + f_3 * pc_x[k] * lsf_140[k];

        t_211[k] = f_3 * pc_y[k] * lsf_140[k];

        t_212[k] = f_15 * ksf_90[k]
                   + f_3 * pc_z[k] * lsf_140[k];

        t_213[k] = f_4 * lsd0_84[k]
                   - f_5 * lsd1_84[k]
                   + f_3 * pc_y[k] * lsf_141[k];

        t_214[k] = f_3 * pc_y[k] * lsf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, ksf_145, ksf_146, ksf_147, \
                         lsd0_89, lsd1_89, lsf_145, lsf_146, lsf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * ksf_145[k]
                   + f_4 * lsd0_89[k]
                   - f_5 * lsd1_89[k]
                   + f_3 * pc_x[k] * lsf_145[k];

        t_216[k] = f_15 * ksf_146[k]
                   + f_3 * pc_x[k] * lsf_146[k];

        t_217[k] = f_15 * ksf_147[k]
                   + f_3 * pc_x[k] * lsf_147[k];

        t_218[k] = f_3 * pc_y[k] * lsf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, ksf_149, lsd0_87, lsd0_88, lsd1_87, \
                         lsd1_88, lsf_146, lsf_147, lsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * ksf_149[k]
                   + f_3 * pc_x[k] * lsf_149[k];

        t_220[k] = f_1 * lsd0_87[k]
                   - f_2 * lsd1_87[k]
                   + f_3 * pc_y[k] * lsf_146[k];

        t_221[k] = f_10 * lsd0_88[k]
                   - f_11 * lsd1_88[k]
                   + f_3 * pc_y[k] * lsf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, ksf_99, ksf_150, \
                         lsd0_89, lsd0_90, lsd1_89, lsd1_90, lsf_148, lsf_149, \
                         lsf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * lsd0_89[k]
                   - f_5 * lsd1_89[k]
                   + f_3 * pc_y[k] * lsf_148[k];

        t_223[k] = f_3 * pc_y[k] * lsf_149[k];

        t_224[k] = f_15 * ksf_99[k]
                   + f_1 * lsd0_89[k]
                   - f_2 * lsd1_89[k]
                   + f_3 * pc_z[k] * lsf_149[k];

        t_225[k] = f_14 * ksf_150[k]
                   + f_1 * lsd0_90[k]
                   - f_2 * lsd1_90[k]
                   + f_3 * pc_x[k] * lsf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, ksf_100, ksf_153, \
                         lsd0_93, lsd1_93, lsf_150, lsf_151, lsf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_13 * ksf_100[k]
                   + f_3 * pc_y[k] * lsf_150[k];

        t_227[k] = f_3 * pc_z[k] * lsf_150[k];

        t_228[k] = f_14 * ksf_153[k]
                   + f_4 * lsd0_93[k]
                   - f_5 * lsd1_93[k]
                   + f_3 * pc_x[k] * lsf_153[k];

        t_229[k] = f_3 * pc_z[k] * lsf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, ksf_156, ksf_158, lsd0_90, \
                         lsd1_90, lsf_152, lsf_153, lsf_156, lsf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * lsd0_90[k]
                   - f_5 * lsd1_90[k]
                   + f_3 * pc_z[k] * lsf_152[k];

        t_231[k] = f_14 * ksf_156[k]
                   + f_3 * pc_x[k] * lsf_156[k];

        t_232[k] = f_3 * pc_z[k] * lsf_153[k];

        t_233[k] = f_14 * ksf_158[k]
                   + f_3 * pc_x[k] * lsf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, ksf_106, \
                         ksf_109, ksf_159, lsd0_93, lsd1_93, lsf_156, lsf_157, \
                         lsf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_14 * ksf_159[k]
                   + f_3 * pc_x[k] * lsf_159[k];

        t_235[k] = f_13 * ksf_106[k]
                   + f_1 * lsd0_93[k]
                   - f_2 * lsd1_93[k]
                   + f_3 * pc_y[k] * lsf_156[k];

        t_236[k] = f_3 * pc_z[k] * lsf_156[k];

        t_237[k] = f_4 * lsd0_93[k]
                   - f_5 * lsd1_93[k]
                   + f_3 * pc_z[k] * lsf_157[k];

        t_238[k] = f_13 * ksf_109[k]
                   + f_3 * pc_y[k] * lsf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, ksg0_150, ksf_100, \
                         ksf_110, ksg1_150, lsd0_95, lsd1_95, lsf_159, \
                         lsf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * lsd0_95[k]
                   - f_2 * lsd1_95[k]
                   + f_3 * pc_z[k] * lsf_159[k];

        t_240[k] = pa_z[k] * ksg0_150[k]
                   - f_6 * pc_z[k] * ksg1_150[k];

        t_241[k] = f_15 * ksf_110[k]
                   + f_3 * pc_y[k] * lsf_160[k];

        t_242[k] = f_7 * ksf_100[k]
                   + f_3 * pc_z[k] * lsf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, ksg0_153, ksf_112, \
                         ksf_165, ksg1_153, lsd0_101, lsd1_101, lsf_162, \
                         lsf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * ksg0_153[k]
                   - f_6 * pc_z[k] * ksg1_153[k];

        t_244[k] = f_15 * ksf_112[k]
                   + f_3 * pc_y[k] * lsf_162[k];

        t_245[k] = f_14 * ksf_165[k]
                   + f_4 * lsd0_101[k]
                   - f_5 * lsd1_101[k]
                   + f_3 * pc_x[k] * lsf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, ksf_166, ksf_167, ksf_168, ksf_169, \
                         lsf_166, lsf_167, lsf_168, lsf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_14 * ksf_166[k]
                   + f_3 * pc_x[k] * lsf_166[k];

        t_247[k] = f_14 * ksf_167[k]
                   + f_3 * pc_x[k] * lsf_167[k];

        t_248[k] = f_14 * ksf_168[k]
                   + f_3 * pc_x[k] * lsf_168[k];

        t_249[k] = f_14 * ksf_169[k]
                   + f_3 * pc_x[k] * lsf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, ksg0_160, ksf_106, ksf_118, \
                         ksg1_160, lsd0_101, lsd1_101, lsf_166, \
                         lsf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * ksg0_160[k]
                   - f_6 * pc_z[k] * ksg1_160[k];

        t_251[k] = f_7 * ksf_106[k]
                   + f_3 * pc_z[k] * lsf_166[k];

        t_252[k] = f_15 * ksf_118[k]
                   + f_4 * lsd0_101[k]
                   - f_5 * lsd1_101[k]
                   + f_3 * pc_y[k] * lsf_168[k];
    }
}

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
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
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksg0_210 = buffer.data(ksg0 + 210);
    const auto *ksg0_213 = buffer.data(ksg0 + 213);
    const auto *ksg0_215 = buffer.data(ksg0 + 215);
    const auto *ksg0_224 = buffer.data(ksg0 + 224);
    const auto *ksg0_225 = buffer.data(ksg0 + 225);
    const auto *ksg0_228 = buffer.data(ksg0 + 228);
    const auto *ksg0_235 = buffer.data(ksg0 + 235);

    const auto *ksf_109 = buffer.data(ksf + 109);
    const auto *ksf_110 = buffer.data(ksf + 110);
    const auto *ksf_116 = buffer.data(ksf + 116);
    const auto *ksf_119 = buffer.data(ksf + 119);
    const auto *ksf_120 = buffer.data(ksf + 120);
    const auto *ksf_122 = buffer.data(ksf + 122);
    const auto *ksf_126 = buffer.data(ksf + 126);
    const auto *ksf_128 = buffer.data(ksf + 128);
    const auto *ksf_129 = buffer.data(ksf + 129);
    const auto *ksf_130 = buffer.data(ksf + 130);
    const auto *ksf_132 = buffer.data(ksf + 132);
    const auto *ksf_136 = buffer.data(ksf + 136);
    const auto *ksf_138 = buffer.data(ksf + 138);
    const auto *ksf_139 = buffer.data(ksf + 139);
    const auto *ksf_140 = buffer.data(ksf + 140);
    const auto *ksf_141 = buffer.data(ksf + 141);
    const auto *ksf_142 = buffer.data(ksf + 142);
    const auto *ksf_146 = buffer.data(ksf + 146);
    const auto *ksf_148 = buffer.data(ksf + 148);
    const auto *ksf_149 = buffer.data(ksf + 149);
    const auto *ksf_150 = buffer.data(ksf + 150);
    const auto *ksf_156 = buffer.data(ksf + 156);
    const auto *ksf_159 = buffer.data(ksf + 159);
    const auto *ksf_160 = buffer.data(ksf + 160);
    const auto *ksf_162 = buffer.data(ksf + 162);
    const auto *ksf_166 = buffer.data(ksf + 166);
    const auto *ksf_168 = buffer.data(ksf + 168);
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
    const auto *ksf_196 = buffer.data(ksf + 196);
    const auto *ksf_197 = buffer.data(ksf + 197);
    const auto *ksf_198 = buffer.data(ksf + 198);
    const auto *ksf_199 = buffer.data(ksf + 199);
    const auto *ksf_200 = buffer.data(ksf + 200);
    const auto *ksf_205 = buffer.data(ksf + 205);
    const auto *ksf_206 = buffer.data(ksf + 206);
    const auto *ksf_207 = buffer.data(ksf + 207);
    const auto *ksf_209 = buffer.data(ksf + 209);
    const auto *ksf_210 = buffer.data(ksf + 210);
    const auto *ksf_213 = buffer.data(ksf + 213);
    const auto *ksf_216 = buffer.data(ksf + 216);
    const auto *ksf_218 = buffer.data(ksf + 218);
    const auto *ksf_219 = buffer.data(ksf + 219);
    const auto *ksf_225 = buffer.data(ksf + 225);
    const auto *ksf_226 = buffer.data(ksf + 226);
    const auto *ksf_227 = buffer.data(ksf + 227);
    const auto *ksf_228 = buffer.data(ksf + 228);
    const auto *ksf_229 = buffer.data(ksf + 229);
    const auto *ksf_230 = buffer.data(ksf + 230);
    const auto *ksf_233 = buffer.data(ksf + 233);
    const auto *ksf_235 = buffer.data(ksf + 235);
    const auto *ksf_236 = buffer.data(ksf + 236);
    const auto *ksf_237 = buffer.data(ksf + 237);
    const auto *ksf_238 = buffer.data(ksf + 238);
    const auto *ksf_239 = buffer.data(ksf + 239);
    const auto *ksf_240 = buffer.data(ksf + 240);
    const auto *ksf_243 = buffer.data(ksf + 243);
    const auto *ksf_245 = buffer.data(ksf + 245);
    const auto *ksf_246 = buffer.data(ksf + 246);
    const auto *ksf_247 = buffer.data(ksf + 247);
    const auto *ksf_248 = buffer.data(ksf + 248);
    const auto *ksf_249 = buffer.data(ksf + 249);

    const auto *ksg1_210 = buffer.data(ksg1 + 210);
    const auto *ksg1_213 = buffer.data(ksg1 + 213);
    const auto *ksg1_215 = buffer.data(ksg1 + 215);
    const auto *ksg1_224 = buffer.data(ksg1 + 224);
    const auto *ksg1_225 = buffer.data(ksg1 + 225);
    const auto *ksg1_228 = buffer.data(ksg1 + 228);
    const auto *ksg1_235 = buffer.data(ksg1 + 235);

    const auto *lsd0_101 = buffer.data(lsd0 + 101);
    const auto *lsd0_102 = buffer.data(lsd0 + 102);
    const auto *lsd0_105 = buffer.data(lsd0 + 105);
    const auto *lsd0_107 = buffer.data(lsd0 + 107);
    const auto *lsd0_108 = buffer.data(lsd0 + 108);
    const auto *lsd0_111 = buffer.data(lsd0 + 111);
    const auto *lsd0_113 = buffer.data(lsd0 + 113);
    const auto *lsd0_117 = buffer.data(lsd0 + 117);
    const auto *lsd0_119 = buffer.data(lsd0 + 119);
    const auto *lsd0_120 = buffer.data(lsd0 + 120);
    const auto *lsd0_123 = buffer.data(lsd0 + 123);
    const auto *lsd0_124 = buffer.data(lsd0 + 124);
    const auto *lsd0_125 = buffer.data(lsd0 + 125);
    const auto *lsd0_126 = buffer.data(lsd0 + 126);
    const auto *lsd0_129 = buffer.data(lsd0 + 129);
    const auto *lsd0_131 = buffer.data(lsd0 + 131);
    const auto *lsd0_137 = buffer.data(lsd0 + 137);
    const auto *lsd0_138 = buffer.data(lsd0 + 138);
    const auto *lsd0_141 = buffer.data(lsd0 + 141);
    const auto *lsd0_143 = buffer.data(lsd0 + 143);
    const auto *lsd0_144 = buffer.data(lsd0 + 144);
    const auto *lsd0_147 = buffer.data(lsd0 + 147);
    const auto *lsd0_149 = buffer.data(lsd0 + 149);

    const auto *lsd1_101 = buffer.data(lsd1 + 101);
    const auto *lsd1_102 = buffer.data(lsd1 + 102);
    const auto *lsd1_105 = buffer.data(lsd1 + 105);
    const auto *lsd1_107 = buffer.data(lsd1 + 107);
    const auto *lsd1_108 = buffer.data(lsd1 + 108);
    const auto *lsd1_111 = buffer.data(lsd1 + 111);
    const auto *lsd1_113 = buffer.data(lsd1 + 113);
    const auto *lsd1_117 = buffer.data(lsd1 + 117);
    const auto *lsd1_119 = buffer.data(lsd1 + 119);
    const auto *lsd1_120 = buffer.data(lsd1 + 120);
    const auto *lsd1_123 = buffer.data(lsd1 + 123);
    const auto *lsd1_124 = buffer.data(lsd1 + 124);
    const auto *lsd1_125 = buffer.data(lsd1 + 125);
    const auto *lsd1_126 = buffer.data(lsd1 + 126);
    const auto *lsd1_129 = buffer.data(lsd1 + 129);
    const auto *lsd1_131 = buffer.data(lsd1 + 131);
    const auto *lsd1_137 = buffer.data(lsd1 + 137);
    const auto *lsd1_138 = buffer.data(lsd1 + 138);
    const auto *lsd1_141 = buffer.data(lsd1 + 141);
    const auto *lsd1_143 = buffer.data(lsd1 + 143);
    const auto *lsd1_144 = buffer.data(lsd1 + 144);
    const auto *lsd1_147 = buffer.data(lsd1 + 147);
    const auto *lsd1_149 = buffer.data(lsd1 + 149);

    const auto *lsf_169 = buffer.data(lsf + 169);
    const auto *lsf_170 = buffer.data(lsf + 170);
    const auto *lsf_172 = buffer.data(lsf + 172);
    const auto *lsf_173 = buffer.data(lsf + 173);
    const auto *lsf_175 = buffer.data(lsf + 175);
    const auto *lsf_176 = buffer.data(lsf + 176);
    const auto *lsf_177 = buffer.data(lsf + 177);
    const auto *lsf_178 = buffer.data(lsf + 178);
    const auto *lsf_179 = buffer.data(lsf + 179);
    const auto *lsf_180 = buffer.data(lsf + 180);
    const auto *lsf_182 = buffer.data(lsf + 182);
    const auto *lsf_183 = buffer.data(lsf + 183);
    const auto *lsf_185 = buffer.data(lsf + 185);
    const auto *lsf_186 = buffer.data(lsf + 186);
    const auto *lsf_187 = buffer.data(lsf + 187);
    const auto *lsf_188 = buffer.data(lsf + 188);
    const auto *lsf_189 = buffer.data(lsf + 189);
    const auto *lsf_190 = buffer.data(lsf + 190);
    const auto *lsf_192 = buffer.data(lsf + 192);
    const auto *lsf_196 = buffer.data(lsf + 196);
    const auto *lsf_197 = buffer.data(lsf + 197);
    const auto *lsf_198 = buffer.data(lsf + 198);
    const auto *lsf_199 = buffer.data(lsf + 199);
    const auto *lsf_200 = buffer.data(lsf + 200);
    const auto *lsf_201 = buffer.data(lsf + 201);
    const auto *lsf_202 = buffer.data(lsf + 202);
    const auto *lsf_205 = buffer.data(lsf + 205);
    const auto *lsf_206 = buffer.data(lsf + 206);
    const auto *lsf_207 = buffer.data(lsf + 207);
    const auto *lsf_208 = buffer.data(lsf + 208);
    const auto *lsf_209 = buffer.data(lsf + 209);
    const auto *lsf_210 = buffer.data(lsf + 210);
    const auto *lsf_211 = buffer.data(lsf + 211);
    const auto *lsf_212 = buffer.data(lsf + 212);
    const auto *lsf_213 = buffer.data(lsf + 213);
    const auto *lsf_216 = buffer.data(lsf + 216);
    const auto *lsf_217 = buffer.data(lsf + 217);
    const auto *lsf_218 = buffer.data(lsf + 218);
    const auto *lsf_219 = buffer.data(lsf + 219);
    const auto *lsf_220 = buffer.data(lsf + 220);
    const auto *lsf_222 = buffer.data(lsf + 222);
    const auto *lsf_225 = buffer.data(lsf + 225);
    const auto *lsf_226 = buffer.data(lsf + 226);
    const auto *lsf_227 = buffer.data(lsf + 227);
    const auto *lsf_228 = buffer.data(lsf + 228);
    const auto *lsf_229 = buffer.data(lsf + 229);
    const auto *lsf_230 = buffer.data(lsf + 230);
    const auto *lsf_232 = buffer.data(lsf + 232);
    const auto *lsf_233 = buffer.data(lsf + 233);
    const auto *lsf_235 = buffer.data(lsf + 235);
    const auto *lsf_236 = buffer.data(lsf + 236);
    const auto *lsf_237 = buffer.data(lsf + 237);
    const auto *lsf_238 = buffer.data(lsf + 238);
    const auto *lsf_239 = buffer.data(lsf + 239);
    const auto *lsf_240 = buffer.data(lsf + 240);
    const auto *lsf_242 = buffer.data(lsf + 242);
    const auto *lsf_243 = buffer.data(lsf + 243);
    const auto *lsf_245 = buffer.data(lsf + 245);
    const auto *lsf_246 = buffer.data(lsf + 246);
    const auto *lsf_247 = buffer.data(lsf + 247);
    const auto *lsf_248 = buffer.data(lsf + 248);
    const auto *lsf_249 = buffer.data(lsf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, ksf_109, ksf_119, ksf_170, \
                         lsd0_101, lsd0_102, lsd1_101, lsd1_102, lsf_169, \
                         lsf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_15 * ksf_119[k]
                   + f_3 * pc_y[k] * lsf_169[k];

        t_254[k] = f_7 * ksf_109[k]
                   + f_1 * lsd0_101[k]
                   - f_2 * lsd1_101[k]
                   + f_3 * pc_z[k] * lsf_169[k];

        t_255[k] = f_14 * ksf_170[k]
                   + f_1 * lsd0_102[k]
                   - f_2 * lsd1_102[k]
                   + f_3 * pc_x[k] * lsf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, ksf_110, ksf_120, \
                         ksf_122, ksf_173, lsd0_105, lsd1_105, lsf_170, lsf_172, \
                         lsf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * ksf_120[k]
                   + f_3 * pc_y[k] * lsf_170[k];

        t_257[k] = f_8 * ksf_110[k]
                   + f_3 * pc_z[k] * lsf_170[k];

        t_258[k] = f_14 * ksf_173[k]
                   + f_4 * lsd0_105[k]
                   - f_5 * lsd1_105[k]
                   + f_3 * pc_x[k] * lsf_173[k];

        t_259[k] = f_14 * ksf_122[k]
                   + f_3 * pc_y[k] * lsf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, ksf_175, ksf_176, ksf_177, ksf_178, \
                         lsd0_107, lsd1_107, lsf_175, lsf_176, lsf_177, \
                         lsf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_14 * ksf_175[k]
                   + f_4 * lsd0_107[k]
                   - f_5 * lsd1_107[k]
                   + f_3 * pc_x[k] * lsf_175[k];

        t_261[k] = f_14 * ksf_176[k]
                   + f_3 * pc_x[k] * lsf_176[k];

        t_262[k] = f_14 * ksf_177[k]
                   + f_3 * pc_x[k] * lsf_177[k];

        t_263[k] = f_14 * ksf_178[k]
                   + f_3 * pc_x[k] * lsf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, ksf_116, ksf_126, ksf_179, \
                         lsd0_105, lsd1_105, lsf_176, lsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_14 * ksf_179[k]
                   + f_3 * pc_x[k] * lsf_179[k];

        t_265[k] = f_14 * ksf_126[k]
                   + f_1 * lsd0_105[k]
                   - f_2 * lsd1_105[k]
                   + f_3 * pc_y[k] * lsf_176[k];

        t_266[k] = f_8 * ksf_116[k]
                   + f_3 * pc_z[k] * lsf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, ksf_119, ksf_128, ksf_129, lsd0_107, \
                         lsd1_107, lsf_178, lsf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * ksf_128[k]
                   + f_4 * lsd0_107[k]
                   - f_5 * lsd1_107[k]
                   + f_3 * pc_y[k] * lsf_178[k];

        t_268[k] = f_14 * ksf_129[k]
                   + f_3 * pc_y[k] * lsf_179[k];

        t_269[k] = f_8 * ksf_119[k]
                   + f_1 * lsd0_107[k]
                   - f_2 * lsd1_107[k]
                   + f_3 * pc_z[k] * lsf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, ksf_120, ksf_130, ksf_180, \
                         lsd0_108, lsd1_108, lsf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_14 * ksf_180[k]
                   + f_1 * lsd0_108[k]
                   - f_2 * lsd1_108[k]
                   + f_3 * pc_x[k] * lsf_180[k];

        t_271[k] = f_8 * ksf_130[k]
                   + f_3 * pc_y[k] * lsf_180[k];

        t_272[k] = f_14 * ksf_120[k]
                   + f_3 * pc_z[k] * lsf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, ksf_132, ksf_183, ksf_185, lsd0_111, \
                         lsd0_113, lsd1_111, lsd1_113, lsf_182, lsf_183, \
                         lsf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_14 * ksf_183[k]
                   + f_4 * lsd0_111[k]
                   - f_5 * lsd1_111[k]
                   + f_3 * pc_x[k] * lsf_183[k];

        t_274[k] = f_8 * ksf_132[k]
                   + f_3 * pc_y[k] * lsf_182[k];

        t_275[k] = f_14 * ksf_185[k]
                   + f_4 * lsd0_113[k]
                   - f_5 * lsd1_113[k]
                   + f_3 * pc_x[k] * lsf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, ksf_186, ksf_187, ksf_188, ksf_189, \
                         lsf_186, lsf_187, lsf_188, lsf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_14 * ksf_186[k]
                   + f_3 * pc_x[k] * lsf_186[k];

        t_277[k] = f_14 * ksf_187[k]
                   + f_3 * pc_x[k] * lsf_187[k];

        t_278[k] = f_14 * ksf_188[k]
                   + f_3 * pc_x[k] * lsf_188[k];

        t_279[k] = f_14 * ksf_189[k]
                   + f_3 * pc_x[k] * lsf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, ksf_126, ksf_136, ksf_138, lsd0_111, \
                         lsd0_113, lsd1_111, lsd1_113, lsf_186, \
                         lsf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * ksf_136[k]
                   + f_1 * lsd0_111[k]
                   - f_2 * lsd1_111[k]
                   + f_3 * pc_y[k] * lsf_186[k];

        t_281[k] = f_14 * ksf_126[k]
                   + f_3 * pc_z[k] * lsf_186[k];

        t_282[k] = f_8 * ksf_138[k]
                   + f_4 * lsd0_113[k]
                   - f_5 * lsd1_113[k]
                   + f_3 * pc_y[k] * lsf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, ksg0_210, ksf_129, \
                         ksf_139, ksf_140, ksg1_210, lsd0_113, lsd1_113, lsf_189, \
                         lsf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * ksf_139[k]
                   + f_3 * pc_y[k] * lsf_189[k];

        t_284[k] = f_14 * ksf_129[k]
                   + f_1 * lsd0_113[k]
                   - f_2 * lsd1_113[k]
                   + f_3 * pc_z[k] * lsf_189[k];

        t_285[k] = pa_y[k] * ksg0_210[k]
                   - f_6 * pc_y[k] * ksg1_210[k];

        t_286[k] = f_7 * ksf_140[k]
                   + f_3 * pc_y[k] * lsf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, ksg0_213, ksg0_215, \
                         ksf_130, ksf_141, ksf_142, ksg1_213, ksg1_215, lsf_190, \
                         lsf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_15 * ksf_130[k]
                   + f_3 * pc_z[k] * lsf_190[k];

        t_288[k] = pa_y[k] * ksg0_213[k]
                   + f_8 * ksf_141[k]
                   - f_6 * pc_y[k] * ksg1_213[k];

        t_289[k] = f_7 * ksf_142[k]
                   + f_3 * pc_y[k] * lsf_192[k];

        t_290[k] = pa_y[k] * ksg0_215[k]
                   - f_6 * pc_y[k] * ksg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, ksf_196, ksf_197, ksf_198, ksf_199, \
                         lsf_196, lsf_197, lsf_198, lsf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_14 * ksf_196[k]
                   + f_3 * pc_x[k] * lsf_196[k];

        t_292[k] = f_14 * ksf_197[k]
                   + f_3 * pc_x[k] * lsf_197[k];

        t_293[k] = f_14 * ksf_198[k]
                   + f_3 * pc_x[k] * lsf_198[k];

        t_294[k] = f_14 * ksf_199[k]
                   + f_3 * pc_x[k] * lsf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, ksf_136, ksf_146, ksf_148, lsd0_117, \
                         lsd0_119, lsd1_117, lsd1_119, lsf_196, \
                         lsf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * ksf_146[k]
                   + f_1 * lsd0_117[k]
                   - f_2 * lsd1_117[k]
                   + f_3 * pc_y[k] * lsf_196[k];

        t_296[k] = f_15 * ksf_136[k]
                   + f_3 * pc_z[k] * lsf_196[k];

        t_297[k] = f_7 * ksf_148[k]
                   + f_4 * lsd0_119[k]
                   - f_5 * lsd1_119[k]
                   + f_3 * pc_y[k] * lsf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, ksg0_224, ksf_149, \
                         ksf_200, ksg1_224, lsd0_120, lsd1_120, lsf_199, \
                         lsf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * ksf_149[k]
                   + f_3 * pc_y[k] * lsf_199[k];

        t_299[k] = pa_y[k] * ksg0_224[k]
                   - f_6 * pc_y[k] * ksg1_224[k];

        t_300[k] = f_14 * ksf_200[k]
                   + f_1 * lsd0_120[k]
                   - f_2 * lsd1_120[k]
                   + f_3 * pc_x[k] * lsf_200[k];

        t_301[k] = f_3 * pc_y[k] * lsf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, ksf_140, lsd0_120, lsd1_120, \
                         lsf_200, lsf_201, lsf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_13 * ksf_140[k]
                   + f_3 * pc_z[k] * lsf_200[k];

        t_303[k] = f_4 * lsd0_120[k]
                   - f_5 * lsd1_120[k]
                   + f_3 * pc_y[k] * lsf_201[k];

        t_304[k] = f_3 * pc_y[k] * lsf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, ksf_205, ksf_206, ksf_207, \
                         lsd0_125, lsd1_125, lsf_205, lsf_206, \
                         lsf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_14 * ksf_205[k]
                   + f_4 * lsd0_125[k]
                   - f_5 * lsd1_125[k]
                   + f_3 * pc_x[k] * lsf_205[k];

        t_306[k] = f_14 * ksf_206[k]
                   + f_3 * pc_x[k] * lsf_206[k];

        t_307[k] = f_14 * ksf_207[k]
                   + f_3 * pc_x[k] * lsf_207[k];

        t_308[k] = f_3 * pc_y[k] * lsf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, ksf_209, lsd0_123, lsd0_124, \
                         lsd1_123, lsd1_124, lsf_206, lsf_207, \
                         lsf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_14 * ksf_209[k]
                   + f_3 * pc_x[k] * lsf_209[k];

        t_310[k] = f_1 * lsd0_123[k]
                   - f_2 * lsd1_123[k]
                   + f_3 * pc_y[k] * lsf_206[k];

        t_311[k] = f_10 * lsd0_124[k]
                   - f_11 * lsd1_124[k]
                   + f_3 * pc_y[k] * lsf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, ksf_149, ksf_210, \
                         lsd0_125, lsd0_126, lsd1_125, lsd1_126, lsf_208, lsf_209, \
                         lsf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * lsd0_125[k]
                   - f_5 * lsd1_125[k]
                   + f_3 * pc_y[k] * lsf_208[k];

        t_313[k] = f_3 * pc_y[k] * lsf_209[k];

        t_314[k] = f_13 * ksf_149[k]
                   + f_1 * lsd0_125[k]
                   - f_2 * lsd1_125[k]
                   + f_3 * pc_z[k] * lsf_209[k];

        t_315[k] = f_8 * ksf_210[k]
                   + f_1 * lsd0_126[k]
                   - f_2 * lsd1_126[k]
                   + f_3 * pc_x[k] * lsf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, ksf_150, ksf_213, \
                         lsd0_129, lsd1_129, lsf_210, lsf_211, \
                         lsf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_12 * ksf_150[k]
                   + f_3 * pc_y[k] * lsf_210[k];

        t_317[k] = f_3 * pc_z[k] * lsf_210[k];

        t_318[k] = f_8 * ksf_213[k]
                   + f_4 * lsd0_129[k]
                   - f_5 * lsd1_129[k]
                   + f_3 * pc_x[k] * lsf_213[k];

        t_319[k] = f_3 * pc_z[k] * lsf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, ksf_216, ksf_218, lsd0_126, \
                         lsd1_126, lsf_212, lsf_213, lsf_216, lsf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * lsd0_126[k]
                   - f_5 * lsd1_126[k]
                   + f_3 * pc_z[k] * lsf_212[k];

        t_321[k] = f_8 * ksf_216[k]
                   + f_3 * pc_x[k] * lsf_216[k];

        t_322[k] = f_3 * pc_z[k] * lsf_213[k];

        t_323[k] = f_8 * ksf_218[k]
                   + f_3 * pc_x[k] * lsf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, ksf_156, \
                         ksf_159, ksf_219, lsd0_129, lsd1_129, lsf_216, lsf_217, \
                         lsf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_8 * ksf_219[k]
                   + f_3 * pc_x[k] * lsf_219[k];

        t_325[k] = f_12 * ksf_156[k]
                   + f_1 * lsd0_129[k]
                   - f_2 * lsd1_129[k]
                   + f_3 * pc_y[k] * lsf_216[k];

        t_326[k] = f_3 * pc_z[k] * lsf_216[k];

        t_327[k] = f_4 * lsd0_129[k]
                   - f_5 * lsd1_129[k]
                   + f_3 * pc_z[k] * lsf_217[k];

        t_328[k] = f_12 * ksf_159[k]
                   + f_3 * pc_y[k] * lsf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_z, pc_y, pc_z, ksg0_225, ksf_150, \
                         ksf_160, ksg1_225, lsd0_131, lsd1_131, lsf_219, \
                         lsf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * lsd0_131[k]
                   - f_2 * lsd1_131[k]
                   + f_3 * pc_z[k] * lsf_219[k];

        t_330[k] = pa_z[k] * ksg0_225[k]
                   - f_6 * pc_z[k] * ksg1_225[k];

        t_331[k] = f_13 * ksf_160[k]
                   + f_3 * pc_y[k] * lsf_220[k];

        t_332[k] = f_7 * ksf_150[k]
                   + f_3 * pc_z[k] * lsf_220[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_x, pc_y, pc_z, ksg0_228, ksf_162, \
                         ksf_225, ksg1_228, lsd0_137, lsd1_137, lsf_222, \
                         lsf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * ksg0_228[k]
                   - f_6 * pc_z[k] * ksg1_228[k];

        t_334[k] = f_13 * ksf_162[k]
                   + f_3 * pc_y[k] * lsf_222[k];

        t_335[k] = f_8 * ksf_225[k]
                   + f_4 * lsd0_137[k]
                   - f_5 * lsd1_137[k]
                   + f_3 * pc_x[k] * lsf_225[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_x, ksf_226, ksf_227, ksf_228, ksf_229, \
                         lsf_226, lsf_227, lsf_228, lsf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_8 * ksf_226[k]
                   + f_3 * pc_x[k] * lsf_226[k];

        t_337[k] = f_8 * ksf_227[k]
                   + f_3 * pc_x[k] * lsf_227[k];

        t_338[k] = f_8 * ksf_228[k]
                   + f_3 * pc_x[k] * lsf_228[k];

        t_339[k] = f_8 * ksf_229[k]
                   + f_3 * pc_x[k] * lsf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pc_y, pc_z, ksg0_235, ksf_156, ksf_168, \
                         ksg1_235, lsd0_137, lsd1_137, lsf_226, \
                         lsf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * ksg0_235[k]
                   - f_6 * pc_z[k] * ksg1_235[k];

        t_341[k] = f_7 * ksf_156[k]
                   + f_3 * pc_z[k] * lsf_226[k];

        t_342[k] = f_13 * ksf_168[k]
                   + f_4 * lsd0_137[k]
                   - f_5 * lsd1_137[k]
                   + f_3 * pc_y[k] * lsf_228[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, ksf_159, ksf_169, ksf_230, \
                         lsd0_137, lsd0_138, lsd1_137, lsd1_138, lsf_229, \
                         lsf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_13 * ksf_169[k]
                   + f_3 * pc_y[k] * lsf_229[k];

        t_344[k] = f_7 * ksf_159[k]
                   + f_1 * lsd0_137[k]
                   - f_2 * lsd1_137[k]
                   + f_3 * pc_z[k] * lsf_229[k];

        t_345[k] = f_8 * ksf_230[k]
                   + f_1 * lsd0_138[k]
                   - f_2 * lsd1_138[k]
                   + f_3 * pc_x[k] * lsf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pc_x, pc_y, pc_z, ksf_160, ksf_170, \
                         ksf_172, ksf_233, lsd0_141, lsd1_141, lsf_230, lsf_232, \
                         lsf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_15 * ksf_170[k]
                   + f_3 * pc_y[k] * lsf_230[k];

        t_347[k] = f_8 * ksf_160[k]
                   + f_3 * pc_z[k] * lsf_230[k];

        t_348[k] = f_8 * ksf_233[k]
                   + f_4 * lsd0_141[k]
                   - f_5 * lsd1_141[k]
                   + f_3 * pc_x[k] * lsf_233[k];

        t_349[k] = f_15 * ksf_172[k]
                   + f_3 * pc_y[k] * lsf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, ksf_235, ksf_236, ksf_237, ksf_238, \
                         lsd0_143, lsd1_143, lsf_235, lsf_236, lsf_237, \
                         lsf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_8 * ksf_235[k]
                   + f_4 * lsd0_143[k]
                   - f_5 * lsd1_143[k]
                   + f_3 * pc_x[k] * lsf_235[k];

        t_351[k] = f_8 * ksf_236[k]
                   + f_3 * pc_x[k] * lsf_236[k];

        t_352[k] = f_8 * ksf_237[k]
                   + f_3 * pc_x[k] * lsf_237[k];

        t_353[k] = f_8 * ksf_238[k]
                   + f_3 * pc_x[k] * lsf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, ksf_166, ksf_176, ksf_239, \
                         lsd0_141, lsd1_141, lsf_236, lsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_8 * ksf_239[k]
                   + f_3 * pc_x[k] * lsf_239[k];

        t_355[k] = f_15 * ksf_176[k]
                   + f_1 * lsd0_141[k]
                   - f_2 * lsd1_141[k]
                   + f_3 * pc_y[k] * lsf_236[k];

        t_356[k] = f_8 * ksf_166[k]
                   + f_3 * pc_z[k] * lsf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, ksf_169, ksf_178, ksf_179, lsd0_143, \
                         lsd1_143, lsf_238, lsf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_15 * ksf_178[k]
                   + f_4 * lsd0_143[k]
                   - f_5 * lsd1_143[k]
                   + f_3 * pc_y[k] * lsf_238[k];

        t_358[k] = f_15 * ksf_179[k]
                   + f_3 * pc_y[k] * lsf_239[k];

        t_359[k] = f_8 * ksf_169[k]
                   + f_1 * lsd0_143[k]
                   - f_2 * lsd1_143[k]
                   + f_3 * pc_z[k] * lsf_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_x, pc_y, pc_z, ksf_170, ksf_180, ksf_240, \
                         lsd0_144, lsd1_144, lsf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_8 * ksf_240[k]
                   + f_1 * lsd0_144[k]
                   - f_2 * lsd1_144[k]
                   + f_3 * pc_x[k] * lsf_240[k];

        t_361[k] = f_14 * ksf_180[k]
                   + f_3 * pc_y[k] * lsf_240[k];

        t_362[k] = f_14 * ksf_170[k]
                   + f_3 * pc_z[k] * lsf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_x, pc_y, ksf_182, ksf_243, ksf_245, lsd0_147, \
                         lsd0_149, lsd1_147, lsd1_149, lsf_242, lsf_243, \
                         lsf_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_8 * ksf_243[k]
                   + f_4 * lsd0_147[k]
                   - f_5 * lsd1_147[k]
                   + f_3 * pc_x[k] * lsf_243[k];

        t_364[k] = f_14 * ksf_182[k]
                   + f_3 * pc_y[k] * lsf_242[k];

        t_365[k] = f_8 * ksf_245[k]
                   + f_4 * lsd0_149[k]
                   - f_5 * lsd1_149[k]
                   + f_3 * pc_x[k] * lsf_245[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, ksf_246, ksf_247, ksf_248, ksf_249, \
                         lsf_246, lsf_247, lsf_248, lsf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_8 * ksf_246[k]
                   + f_3 * pc_x[k] * lsf_246[k];

        t_367[k] = f_8 * ksf_247[k]
                   + f_3 * pc_x[k] * lsf_247[k];

        t_368[k] = f_8 * ksf_248[k]
                   + f_3 * pc_x[k] * lsf_248[k];

        t_369[k] = f_8 * ksf_249[k]
                   + f_3 * pc_x[k] * lsf_249[k];
    }
}

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
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
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 1.5 / q;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksg0_300 = buffer.data(ksg0 + 300);
    const auto *ksg0_303 = buffer.data(ksg0 + 303);
    const auto *ksg0_305 = buffer.data(ksg0 + 305);
    const auto *ksg0_314 = buffer.data(ksg0 + 314);
    const auto *ksg0_315 = buffer.data(ksg0 + 315);
    const auto *ksg0_318 = buffer.data(ksg0 + 318);
    const auto *ksg0_420 = buffer.data(ksg0 + 420);
    const auto *ksg0_423 = buffer.data(ksg0 + 423);
    const auto *ksg0_430 = buffer.data(ksg0 + 430);
    const auto *ksg0_432 = buffer.data(ksg0 + 432);
    const auto *ksg0_434 = buffer.data(ksg0 + 434);
    const auto *ksg0_440 = buffer.data(ksg0 + 440);
    const auto *ksg0_445 = buffer.data(ksg0 + 445);
    const auto *ksg0_447 = buffer.data(ksg0 + 447);
    const auto *ksg0_449 = buffer.data(ksg0 + 449);
    const auto *ksg0_450 = buffer.data(ksg0 + 450);
    const auto *ksg0_453 = buffer.data(ksg0 + 453);
    const auto *ksg0_455 = buffer.data(ksg0 + 455);
    const auto *ksg0_460 = buffer.data(ksg0 + 460);
    const auto *ksg0_462 = buffer.data(ksg0 + 462);
    const auto *ksg0_464 = buffer.data(ksg0 + 464);
    const auto *ksg0_465 = buffer.data(ksg0 + 465);
    const auto *ksg0_468 = buffer.data(ksg0 + 468);
    const auto *ksg0_470 = buffer.data(ksg0 + 470);
    const auto *ksg0_475 = buffer.data(ksg0 + 475);
    const auto *ksg0_477 = buffer.data(ksg0 + 477);
    const auto *ksg0_479 = buffer.data(ksg0 + 479);
    const auto *ksg0_480 = buffer.data(ksg0 + 480);
    const auto *ksg0_483 = buffer.data(ksg0 + 483);
    const auto *ksg0_485 = buffer.data(ksg0 + 485);
    const auto *ksg0_490 = buffer.data(ksg0 + 490);
    const auto *ksg0_492 = buffer.data(ksg0 + 492);

    const auto *ksf_176 = buffer.data(ksf + 176);
    const auto *ksf_179 = buffer.data(ksf + 179);
    const auto *ksf_180 = buffer.data(ksf + 180);
    const auto *ksf_186 = buffer.data(ksf + 186);
    const auto *ksf_188 = buffer.data(ksf + 188);
    const auto *ksf_189 = buffer.data(ksf + 189);
    const auto *ksf_190 = buffer.data(ksf + 190);
    const auto *ksf_192 = buffer.data(ksf + 192);
    const auto *ksf_196 = buffer.data(ksf + 196);
    const auto *ksf_198 = buffer.data(ksf + 198);
    const auto *ksf_199 = buffer.data(ksf + 199);
    const auto *ksf_200 = buffer.data(ksf + 200);
    const auto *ksf_201 = buffer.data(ksf + 201);
    const auto *ksf_202 = buffer.data(ksf + 202);
    const auto *ksf_206 = buffer.data(ksf + 206);
    const auto *ksf_208 = buffer.data(ksf + 208);
    const auto *ksf_209 = buffer.data(ksf + 209);
    const auto *ksf_210 = buffer.data(ksf + 210);
    const auto *ksf_216 = buffer.data(ksf + 216);
    const auto *ksf_219 = buffer.data(ksf + 219);
    const auto *ksf_220 = buffer.data(ksf + 220);
    const auto *ksf_222 = buffer.data(ksf + 222);
    const auto *ksf_226 = buffer.data(ksf + 226);
    const auto *ksf_229 = buffer.data(ksf + 229);
    const auto *ksf_230 = buffer.data(ksf + 230);
    const auto *ksf_232 = buffer.data(ksf + 232);
    const auto *ksf_236 = buffer.data(ksf + 236);
    const auto *ksf_239 = buffer.data(ksf + 239);
    const auto *ksf_240 = buffer.data(ksf + 240);
    const auto *ksf_242 = buffer.data(ksf + 242);
    const auto *ksf_246 = buffer.data(ksf + 246);
    const auto *ksf_249 = buffer.data(ksf + 249);
    const auto *ksf_250 = buffer.data(ksf + 250);
    const auto *ksf_252 = buffer.data(ksf + 252);
    const auto *ksf_253 = buffer.data(ksf + 253);
    const auto *ksf_255 = buffer.data(ksf + 255);
    const auto *ksf_256 = buffer.data(ksf + 256);
    const auto *ksf_257 = buffer.data(ksf + 257);
    const auto *ksf_258 = buffer.data(ksf + 258);
    const auto *ksf_259 = buffer.data(ksf + 259);
    const auto *ksf_266 = buffer.data(ksf + 266);
    const auto *ksf_267 = buffer.data(ksf + 267);
    const auto *ksf_268 = buffer.data(ksf + 268);
    const auto *ksf_269 = buffer.data(ksf + 269);
    const auto *ksf_270 = buffer.data(ksf + 270);
    const auto *ksf_275 = buffer.data(ksf + 275);
    const auto *ksf_276 = buffer.data(ksf + 276);
    const auto *ksf_277 = buffer.data(ksf + 277);
    const auto *ksf_279 = buffer.data(ksf + 279);
    const auto *ksf_280 = buffer.data(ksf + 280);
    const auto *ksf_283 = buffer.data(ksf + 283);
    const auto *ksf_286 = buffer.data(ksf + 286);
    const auto *ksf_288 = buffer.data(ksf + 288);
    const auto *ksf_289 = buffer.data(ksf + 289);
    const auto *ksf_295 = buffer.data(ksf + 295);
    const auto *ksf_296 = buffer.data(ksf + 296);
    const auto *ksf_297 = buffer.data(ksf + 297);
    const auto *ksf_298 = buffer.data(ksf + 298);
    const auto *ksf_299 = buffer.data(ksf + 299);
    const auto *ksf_300 = buffer.data(ksf + 300);
    const auto *ksf_303 = buffer.data(ksf + 303);
    const auto *ksf_305 = buffer.data(ksf + 305);
    const auto *ksf_306 = buffer.data(ksf + 306);
    const auto *ksf_307 = buffer.data(ksf + 307);
    const auto *ksf_308 = buffer.data(ksf + 308);
    const auto *ksf_309 = buffer.data(ksf + 309);
    const auto *ksf_310 = buffer.data(ksf + 310);
    const auto *ksf_313 = buffer.data(ksf + 313);
    const auto *ksf_315 = buffer.data(ksf + 315);
    const auto *ksf_316 = buffer.data(ksf + 316);
    const auto *ksf_317 = buffer.data(ksf + 317);
    const auto *ksf_318 = buffer.data(ksf + 318);
    const auto *ksf_319 = buffer.data(ksf + 319);
    const auto *ksf_320 = buffer.data(ksf + 320);
    const auto *ksf_323 = buffer.data(ksf + 323);
    const auto *ksf_325 = buffer.data(ksf + 325);
    const auto *ksf_326 = buffer.data(ksf + 326);
    const auto *ksf_327 = buffer.data(ksf + 327);
    const auto *ksf_328 = buffer.data(ksf + 328);
    const auto *ksf_329 = buffer.data(ksf + 329);

    const auto *ksg1_300 = buffer.data(ksg1 + 300);
    const auto *ksg1_303 = buffer.data(ksg1 + 303);
    const auto *ksg1_305 = buffer.data(ksg1 + 305);
    const auto *ksg1_314 = buffer.data(ksg1 + 314);
    const auto *ksg1_315 = buffer.data(ksg1 + 315);
    const auto *ksg1_318 = buffer.data(ksg1 + 318);
    const auto *ksg1_420 = buffer.data(ksg1 + 420);
    const auto *ksg1_423 = buffer.data(ksg1 + 423);
    const auto *ksg1_430 = buffer.data(ksg1 + 430);
    const auto *ksg1_432 = buffer.data(ksg1 + 432);
    const auto *ksg1_434 = buffer.data(ksg1 + 434);
    const auto *ksg1_440 = buffer.data(ksg1 + 440);
    const auto *ksg1_445 = buffer.data(ksg1 + 445);
    const auto *ksg1_447 = buffer.data(ksg1 + 447);
    const auto *ksg1_449 = buffer.data(ksg1 + 449);
    const auto *ksg1_450 = buffer.data(ksg1 + 450);
    const auto *ksg1_453 = buffer.data(ksg1 + 453);
    const auto *ksg1_455 = buffer.data(ksg1 + 455);
    const auto *ksg1_460 = buffer.data(ksg1 + 460);
    const auto *ksg1_462 = buffer.data(ksg1 + 462);
    const auto *ksg1_464 = buffer.data(ksg1 + 464);
    const auto *ksg1_465 = buffer.data(ksg1 + 465);
    const auto *ksg1_468 = buffer.data(ksg1 + 468);
    const auto *ksg1_470 = buffer.data(ksg1 + 470);
    const auto *ksg1_475 = buffer.data(ksg1 + 475);
    const auto *ksg1_477 = buffer.data(ksg1 + 477);
    const auto *ksg1_479 = buffer.data(ksg1 + 479);
    const auto *ksg1_480 = buffer.data(ksg1 + 480);
    const auto *ksg1_483 = buffer.data(ksg1 + 483);
    const auto *ksg1_485 = buffer.data(ksg1 + 485);
    const auto *ksg1_490 = buffer.data(ksg1 + 490);
    const auto *ksg1_492 = buffer.data(ksg1 + 492);

    const auto *lsd0_147 = buffer.data(lsd0 + 147);
    const auto *lsd0_149 = buffer.data(lsd0 + 149);
    const auto *lsd0_150 = buffer.data(lsd0 + 150);
    const auto *lsd0_153 = buffer.data(lsd0 + 153);
    const auto *lsd0_155 = buffer.data(lsd0 + 155);
    const auto *lsd0_159 = buffer.data(lsd0 + 159);
    const auto *lsd0_161 = buffer.data(lsd0 + 161);
    const auto *lsd0_162 = buffer.data(lsd0 + 162);
    const auto *lsd0_165 = buffer.data(lsd0 + 165);
    const auto *lsd0_166 = buffer.data(lsd0 + 166);
    const auto *lsd0_167 = buffer.data(lsd0 + 167);
    const auto *lsd0_168 = buffer.data(lsd0 + 168);

    const auto *lsd1_147 = buffer.data(lsd1 + 147);
    const auto *lsd1_149 = buffer.data(lsd1 + 149);
    const auto *lsd1_150 = buffer.data(lsd1 + 150);
    const auto *lsd1_153 = buffer.data(lsd1 + 153);
    const auto *lsd1_155 = buffer.data(lsd1 + 155);
    const auto *lsd1_159 = buffer.data(lsd1 + 159);
    const auto *lsd1_161 = buffer.data(lsd1 + 161);
    const auto *lsd1_162 = buffer.data(lsd1 + 162);
    const auto *lsd1_165 = buffer.data(lsd1 + 165);
    const auto *lsd1_166 = buffer.data(lsd1 + 166);
    const auto *lsd1_167 = buffer.data(lsd1 + 167);
    const auto *lsd1_168 = buffer.data(lsd1 + 168);

    const auto *lsf_246 = buffer.data(lsf + 246);
    const auto *lsf_248 = buffer.data(lsf + 248);
    const auto *lsf_249 = buffer.data(lsf + 249);
    const auto *lsf_250 = buffer.data(lsf + 250);
    const auto *lsf_252 = buffer.data(lsf + 252);
    const auto *lsf_253 = buffer.data(lsf + 253);
    const auto *lsf_255 = buffer.data(lsf + 255);
    const auto *lsf_256 = buffer.data(lsf + 256);
    const auto *lsf_257 = buffer.data(lsf + 257);
    const auto *lsf_258 = buffer.data(lsf + 258);
    const auto *lsf_259 = buffer.data(lsf + 259);
    const auto *lsf_260 = buffer.data(lsf + 260);
    const auto *lsf_262 = buffer.data(lsf + 262);
    const auto *lsf_266 = buffer.data(lsf + 266);
    const auto *lsf_267 = buffer.data(lsf + 267);
    const auto *lsf_268 = buffer.data(lsf + 268);
    const auto *lsf_269 = buffer.data(lsf + 269);
    const auto *lsf_270 = buffer.data(lsf + 270);
    const auto *lsf_271 = buffer.data(lsf + 271);
    const auto *lsf_272 = buffer.data(lsf + 272);
    const auto *lsf_275 = buffer.data(lsf + 275);
    const auto *lsf_276 = buffer.data(lsf + 276);
    const auto *lsf_277 = buffer.data(lsf + 277);
    const auto *lsf_278 = buffer.data(lsf + 278);
    const auto *lsf_279 = buffer.data(lsf + 279);
    const auto *lsf_280 = buffer.data(lsf + 280);
    const auto *lsf_281 = buffer.data(lsf + 281);
    const auto *lsf_282 = buffer.data(lsf + 282);
    const auto *lsf_283 = buffer.data(lsf + 283);
    const auto *lsf_286 = buffer.data(lsf + 286);
    const auto *lsf_288 = buffer.data(lsf + 288);
    const auto *lsf_289 = buffer.data(lsf + 289);
    const auto *lsf_290 = buffer.data(lsf + 290);
    const auto *lsf_292 = buffer.data(lsf + 292);
    const auto *lsf_296 = buffer.data(lsf + 296);
    const auto *lsf_297 = buffer.data(lsf + 297);
    const auto *lsf_298 = buffer.data(lsf + 298);
    const auto *lsf_299 = buffer.data(lsf + 299);
    const auto *lsf_300 = buffer.data(lsf + 300);
    const auto *lsf_302 = buffer.data(lsf + 302);
    const auto *lsf_306 = buffer.data(lsf + 306);
    const auto *lsf_307 = buffer.data(lsf + 307);
    const auto *lsf_308 = buffer.data(lsf + 308);
    const auto *lsf_309 = buffer.data(lsf + 309);
    const auto *lsf_310 = buffer.data(lsf + 310);
    const auto *lsf_312 = buffer.data(lsf + 312);
    const auto *lsf_316 = buffer.data(lsf + 316);
    const auto *lsf_317 = buffer.data(lsf + 317);
    const auto *lsf_318 = buffer.data(lsf + 318);
    const auto *lsf_319 = buffer.data(lsf + 319);
    const auto *lsf_320 = buffer.data(lsf + 320);
    const auto *lsf_322 = buffer.data(lsf + 322);
    const auto *lsf_326 = buffer.data(lsf + 326);
    const auto *lsf_327 = buffer.data(lsf + 327);
    const auto *lsf_328 = buffer.data(lsf + 328);
    const auto *lsf_329 = buffer.data(lsf + 329);

#pragma omp simd aligned(t_370, t_371, t_372, pc_y, pc_z, ksf_176, ksf_186, ksf_188, lsd0_147, \
                         lsd0_149, lsd1_147, lsd1_149, lsf_246, \
                         lsf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * ksf_186[k]
                   + f_1 * lsd0_147[k]
                   - f_2 * lsd1_147[k]
                   + f_3 * pc_y[k] * lsf_246[k];

        t_371[k] = f_14 * ksf_176[k]
                   + f_3 * pc_z[k] * lsf_246[k];

        t_372[k] = f_14 * ksf_188[k]
                   + f_4 * lsd0_149[k]
                   - f_5 * lsd1_149[k]
                   + f_3 * pc_y[k] * lsf_248[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, pc_z, ksf_179, ksf_189, ksf_250, \
                         lsd0_149, lsd0_150, lsd1_149, lsd1_150, lsf_249, \
                         lsf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * ksf_189[k]
                   + f_3 * pc_y[k] * lsf_249[k];

        t_374[k] = f_14 * ksf_179[k]
                   + f_1 * lsd0_149[k]
                   - f_2 * lsd1_149[k]
                   + f_3 * pc_z[k] * lsf_249[k];

        t_375[k] = f_8 * ksf_250[k]
                   + f_1 * lsd0_150[k]
                   - f_2 * lsd1_150[k]
                   + f_3 * pc_x[k] * lsf_250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, pc_y, pc_z, ksf_180, ksf_190, \
                         ksf_192, ksf_253, lsd0_153, lsd1_153, lsf_250, lsf_252, \
                         lsf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * ksf_190[k]
                   + f_3 * pc_y[k] * lsf_250[k];

        t_377[k] = f_15 * ksf_180[k]
                   + f_3 * pc_z[k] * lsf_250[k];

        t_378[k] = f_8 * ksf_253[k]
                   + f_4 * lsd0_153[k]
                   - f_5 * lsd1_153[k]
                   + f_3 * pc_x[k] * lsf_253[k];

        t_379[k] = f_8 * ksf_192[k]
                   + f_3 * pc_y[k] * lsf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, ksf_255, ksf_256, ksf_257, ksf_258, \
                         lsd0_155, lsd1_155, lsf_255, lsf_256, lsf_257, \
                         lsf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_8 * ksf_255[k]
                   + f_4 * lsd0_155[k]
                   - f_5 * lsd1_155[k]
                   + f_3 * pc_x[k] * lsf_255[k];

        t_381[k] = f_8 * ksf_256[k]
                   + f_3 * pc_x[k] * lsf_256[k];

        t_382[k] = f_8 * ksf_257[k]
                   + f_3 * pc_x[k] * lsf_257[k];

        t_383[k] = f_8 * ksf_258[k]
                   + f_3 * pc_x[k] * lsf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, ksf_186, ksf_196, ksf_259, \
                         lsd0_153, lsd1_153, lsf_256, lsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_8 * ksf_259[k]
                   + f_3 * pc_x[k] * lsf_259[k];

        t_385[k] = f_8 * ksf_196[k]
                   + f_1 * lsd0_153[k]
                   - f_2 * lsd1_153[k]
                   + f_3 * pc_y[k] * lsf_256[k];

        t_386[k] = f_15 * ksf_186[k]
                   + f_3 * pc_z[k] * lsf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, ksg0_300, ksf_189, \
                         ksf_198, ksf_199, ksg1_300, lsd0_155, lsd1_155, lsf_258, \
                         lsf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * ksf_198[k]
                   + f_4 * lsd0_155[k]
                   - f_5 * lsd1_155[k]
                   + f_3 * pc_y[k] * lsf_258[k];

        t_388[k] = f_8 * ksf_199[k]
                   + f_3 * pc_y[k] * lsf_259[k];

        t_389[k] = f_15 * ksf_189[k]
                   + f_1 * lsd0_155[k]
                   - f_2 * lsd1_155[k]
                   + f_3 * pc_z[k] * lsf_259[k];

        t_390[k] = pa_y[k] * ksg0_300[k]
                   - f_6 * pc_y[k] * ksg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_y, pc_z, ksg0_303, ksf_190, \
                         ksf_200, ksf_201, ksf_202, ksg1_303, lsf_260, \
                         lsf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * ksf_200[k]
                   + f_3 * pc_y[k] * lsf_260[k];

        t_392[k] = f_13 * ksf_190[k]
                   + f_3 * pc_z[k] * lsf_260[k];

        t_393[k] = pa_y[k] * ksg0_303[k]
                   + f_8 * ksf_201[k]
                   - f_6 * pc_y[k] * ksg1_303[k];

        t_394[k] = f_7 * ksf_202[k]
                   + f_3 * pc_y[k] * lsf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, ksg0_305, ksf_266, \
                         ksf_267, ksf_268, ksg1_305, lsf_266, lsf_267, \
                         lsf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * ksg0_305[k]
                   - f_6 * pc_y[k] * ksg1_305[k];

        t_396[k] = f_8 * ksf_266[k]
                   + f_3 * pc_x[k] * lsf_266[k];

        t_397[k] = f_8 * ksf_267[k]
                   + f_3 * pc_x[k] * lsf_267[k];

        t_398[k] = f_8 * ksf_268[k]
                   + f_3 * pc_x[k] * lsf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, ksf_196, ksf_206, ksf_269, \
                         lsd0_159, lsd1_159, lsf_266, lsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_8 * ksf_269[k]
                   + f_3 * pc_x[k] * lsf_269[k];

        t_400[k] = f_7 * ksf_206[k]
                   + f_1 * lsd0_159[k]
                   - f_2 * lsd1_159[k]
                   + f_3 * pc_y[k] * lsf_266[k];

        t_401[k] = f_13 * ksf_196[k]
                   + f_3 * pc_z[k] * lsf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, ksg0_314, ksf_208, ksf_209, \
                         ksg1_314, lsd0_161, lsd1_161, lsf_268, \
                         lsf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * ksf_208[k]
                   + f_4 * lsd0_161[k]
                   - f_5 * lsd1_161[k]
                   + f_3 * pc_y[k] * lsf_268[k];

        t_403[k] = f_7 * ksf_209[k]
                   + f_3 * pc_y[k] * lsf_269[k];

        t_404[k] = pa_y[k] * ksg0_314[k]
                   - f_6 * pc_y[k] * ksg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, pc_z, ksf_200, \
                         ksf_270, lsd0_162, lsd1_162, lsf_270, lsf_271, \
                         lsf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_8 * ksf_270[k]
                   + f_1 * lsd0_162[k]
                   - f_2 * lsd1_162[k]
                   + f_3 * pc_x[k] * lsf_270[k];

        t_406[k] = f_3 * pc_y[k] * lsf_270[k];

        t_407[k] = f_12 * ksf_200[k]
                   + f_3 * pc_z[k] * lsf_270[k];

        t_408[k] = f_4 * lsd0_162[k]
                   - f_5 * lsd1_162[k]
                   + f_3 * pc_y[k] * lsf_271[k];

        t_409[k] = f_3 * pc_y[k] * lsf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, ksf_275, ksf_276, ksf_277, \
                         lsd0_167, lsd1_167, lsf_275, lsf_276, \
                         lsf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_8 * ksf_275[k]
                   + f_4 * lsd0_167[k]
                   - f_5 * lsd1_167[k]
                   + f_3 * pc_x[k] * lsf_275[k];

        t_411[k] = f_8 * ksf_276[k]
                   + f_3 * pc_x[k] * lsf_276[k];

        t_412[k] = f_8 * ksf_277[k]
                   + f_3 * pc_x[k] * lsf_277[k];

        t_413[k] = f_3 * pc_y[k] * lsf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, ksf_279, lsd0_165, lsd0_166, \
                         lsd1_165, lsd1_166, lsf_276, lsf_277, \
                         lsf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_8 * ksf_279[k]
                   + f_3 * pc_x[k] * lsf_279[k];

        t_415[k] = f_1 * lsd0_165[k]
                   - f_2 * lsd1_165[k]
                   + f_3 * pc_y[k] * lsf_276[k];

        t_416[k] = f_10 * lsd0_166[k]
                   - f_11 * lsd1_166[k]
                   + f_3 * pc_y[k] * lsf_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pa_x, pc_x, pc_y, pc_z, ksg0_420, \
                         ksf_209, ksf_280, ksg1_420, lsd0_167, lsd1_167, lsf_278, \
                         lsf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * lsd0_167[k]
                   - f_5 * lsd1_167[k]
                   + f_3 * pc_y[k] * lsf_278[k];

        t_418[k] = f_3 * pc_y[k] * lsf_279[k];

        t_419[k] = f_12 * ksf_209[k]
                   + f_1 * lsd0_167[k]
                   - f_2 * lsd1_167[k]
                   + f_3 * pc_z[k] * lsf_279[k];

        t_420[k] = pa_x[k] * ksg0_420[k]
                   + f_15 * ksf_280[k]
                   - f_6 * pc_x[k] * ksg1_420[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pa_x, pc_x, pc_y, pc_z, ksg0_423, \
                         ksf_210, ksf_283, ksg1_423, lsf_280, lsf_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_9 * ksf_210[k]
                   + f_3 * pc_y[k] * lsf_280[k];

        t_422[k] = f_3 * pc_z[k] * lsf_280[k];

        t_423[k] = pa_x[k] * ksg0_423[k]
                   + f_8 * ksf_283[k]
                   - f_6 * pc_x[k] * ksg1_423[k];

        t_424[k] = f_3 * pc_z[k] * lsf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_z, ksf_286, ksf_288, lsd0_168, \
                         lsd1_168, lsf_282, lsf_283, lsf_286, lsf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * lsd0_168[k]
                   - f_5 * lsd1_168[k]
                   + f_3 * pc_z[k] * lsf_282[k];

        t_426[k] = f_7 * ksf_286[k]
                   + f_3 * pc_x[k] * lsf_286[k];

        t_427[k] = f_3 * pc_z[k] * lsf_283[k];

        t_428[k] = f_7 * ksf_288[k]
                   + f_3 * pc_x[k] * lsf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_x, pc_x, pc_z, ksg0_430, ksg0_432, \
                         ksf_289, ksg1_430, ksg1_432, lsf_286, \
                         lsf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_7 * ksf_289[k]
                   + f_3 * pc_x[k] * lsf_289[k];

        t_430[k] = pa_x[k] * ksg0_430[k]
                   - f_6 * pc_x[k] * ksg1_430[k];

        t_431[k] = f_3 * pc_z[k] * lsf_286[k];

        t_432[k] = pa_x[k] * ksg0_432[k]
                   - f_6 * pc_x[k] * ksg1_432[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, pa_x, pa_z, pc_x, pc_y, pc_z, ksg0_315, \
                         ksg0_434, ksf_219, ksg1_315, ksg1_434, \
                         lsf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_9 * ksf_219[k]
                   + f_3 * pc_y[k] * lsf_289[k];

        t_434[k] = pa_x[k] * ksg0_434[k]
                   - f_6 * pc_x[k] * ksg1_434[k];

        t_435[k] = pa_z[k] * ksg0_315[k]
                   - f_6 * pc_z[k] * ksg1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pa_z, pc_y, pc_z, ksg0_318, ksf_210, \
                         ksf_220, ksf_222, ksg1_318, lsf_290, lsf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_12 * ksf_220[k]
                   + f_3 * pc_y[k] * lsf_290[k];

        t_437[k] = f_7 * ksf_210[k]
                   + f_3 * pc_z[k] * lsf_290[k];

        t_438[k] = pa_z[k] * ksg0_318[k]
                   - f_6 * pc_z[k] * ksg1_318[k];

        t_439[k] = f_12 * ksf_222[k]
                   + f_3 * pc_y[k] * lsf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pa_x, pc_x, ksg0_440, ksf_295, ksf_296, \
                         ksf_297, ksf_298, ksg1_440, lsf_296, lsf_297, \
                         lsf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = pa_x[k] * ksg0_440[k]
                   + f_8 * ksf_295[k]
                   - f_6 * pc_x[k] * ksg1_440[k];

        t_441[k] = f_7 * ksf_296[k]
                   + f_3 * pc_x[k] * lsf_296[k];

        t_442[k] = f_7 * ksf_297[k]
                   + f_3 * pc_x[k] * lsf_297[k];

        t_443[k] = f_7 * ksf_298[k]
                   + f_3 * pc_x[k] * lsf_298[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pa_x, pc_x, pc_z, ksg0_445, ksg0_447, \
                         ksf_216, ksf_299, ksg1_445, ksg1_447, lsf_296, \
                         lsf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_7 * ksf_299[k]
                   + f_3 * pc_x[k] * lsf_299[k];

        t_445[k] = pa_x[k] * ksg0_445[k]
                   - f_6 * pc_x[k] * ksg1_445[k];

        t_446[k] = f_7 * ksf_216[k]
                   + f_3 * pc_z[k] * lsf_296[k];

        t_447[k] = pa_x[k] * ksg0_447[k]
                   - f_6 * pc_x[k] * ksg1_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pa_x, pc_x, pc_y, ksg0_449, ksg0_450, \
                         ksf_229, ksf_230, ksf_300, ksg1_449, ksg1_450, lsf_299, \
                         lsf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_12 * ksf_229[k]
                   + f_3 * pc_y[k] * lsf_299[k];

        t_449[k] = pa_x[k] * ksg0_449[k]
                   - f_6 * pc_x[k] * ksg1_449[k];

        t_450[k] = pa_x[k] * ksg0_450[k]
                   + f_15 * ksf_300[k]
                   - f_6 * pc_x[k] * ksg1_450[k];

        t_451[k] = f_13 * ksf_230[k]
                   + f_3 * pc_y[k] * lsf_300[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pa_x, pc_x, pc_y, pc_z, ksg0_453, ksf_220, \
                         ksf_232, ksf_303, ksg1_453, lsf_300, lsf_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_8 * ksf_220[k]
                   + f_3 * pc_z[k] * lsf_300[k];

        t_453[k] = pa_x[k] * ksg0_453[k]
                   + f_8 * ksf_303[k]
                   - f_6 * pc_x[k] * ksg1_453[k];

        t_454[k] = f_13 * ksf_232[k]
                   + f_3 * pc_y[k] * lsf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pa_x, pc_x, ksg0_455, ksf_305, ksf_306, \
                         ksf_307, ksf_308, ksg1_455, lsf_306, lsf_307, \
                         lsf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = pa_x[k] * ksg0_455[k]
                   + f_8 * ksf_305[k]
                   - f_6 * pc_x[k] * ksg1_455[k];

        t_456[k] = f_7 * ksf_306[k]
                   + f_3 * pc_x[k] * lsf_306[k];

        t_457[k] = f_7 * ksf_307[k]
                   + f_3 * pc_x[k] * lsf_307[k];

        t_458[k] = f_7 * ksf_308[k]
                   + f_3 * pc_x[k] * lsf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, t_462, pa_x, pc_x, pc_z, ksg0_460, ksg0_462, \
                         ksf_226, ksf_309, ksg1_460, ksg1_462, lsf_306, \
                         lsf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_7 * ksf_309[k]
                   + f_3 * pc_x[k] * lsf_309[k];

        t_460[k] = pa_x[k] * ksg0_460[k]
                   - f_6 * pc_x[k] * ksg1_460[k];

        t_461[k] = f_8 * ksf_226[k]
                   + f_3 * pc_z[k] * lsf_306[k];

        t_462[k] = pa_x[k] * ksg0_462[k]
                   - f_6 * pc_x[k] * ksg1_462[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pa_x, pc_x, pc_y, ksg0_464, ksg0_465, \
                         ksf_239, ksf_240, ksf_310, ksg1_464, ksg1_465, lsf_309, \
                         lsf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_13 * ksf_239[k]
                   + f_3 * pc_y[k] * lsf_309[k];

        t_464[k] = pa_x[k] * ksg0_464[k]
                   - f_6 * pc_x[k] * ksg1_464[k];

        t_465[k] = pa_x[k] * ksg0_465[k]
                   + f_15 * ksf_310[k]
                   - f_6 * pc_x[k] * ksg1_465[k];

        t_466[k] = f_15 * ksf_240[k]
                   + f_3 * pc_y[k] * lsf_310[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pa_x, pc_x, pc_y, pc_z, ksg0_468, ksf_230, \
                         ksf_242, ksf_313, ksg1_468, lsf_310, lsf_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_14 * ksf_230[k]
                   + f_3 * pc_z[k] * lsf_310[k];

        t_468[k] = pa_x[k] * ksg0_468[k]
                   + f_8 * ksf_313[k]
                   - f_6 * pc_x[k] * ksg1_468[k];

        t_469[k] = f_15 * ksf_242[k]
                   + f_3 * pc_y[k] * lsf_312[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pa_x, pc_x, ksg0_470, ksf_315, ksf_316, \
                         ksf_317, ksf_318, ksg1_470, lsf_316, lsf_317, \
                         lsf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = pa_x[k] * ksg0_470[k]
                   + f_8 * ksf_315[k]
                   - f_6 * pc_x[k] * ksg1_470[k];

        t_471[k] = f_7 * ksf_316[k]
                   + f_3 * pc_x[k] * lsf_316[k];

        t_472[k] = f_7 * ksf_317[k]
                   + f_3 * pc_x[k] * lsf_317[k];

        t_473[k] = f_7 * ksf_318[k]
                   + f_3 * pc_x[k] * lsf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pa_x, pc_x, pc_z, ksg0_475, ksg0_477, \
                         ksf_236, ksf_319, ksg1_475, ksg1_477, lsf_316, \
                         lsf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_7 * ksf_319[k]
                   + f_3 * pc_x[k] * lsf_319[k];

        t_475[k] = pa_x[k] * ksg0_475[k]
                   - f_6 * pc_x[k] * ksg1_475[k];

        t_476[k] = f_14 * ksf_236[k]
                   + f_3 * pc_z[k] * lsf_316[k];

        t_477[k] = pa_x[k] * ksg0_477[k]
                   - f_6 * pc_x[k] * ksg1_477[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pa_x, pc_x, pc_y, ksg0_479, ksg0_480, \
                         ksf_249, ksf_250, ksf_320, ksg1_479, ksg1_480, lsf_319, \
                         lsf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_15 * ksf_249[k]
                   + f_3 * pc_y[k] * lsf_319[k];

        t_479[k] = pa_x[k] * ksg0_479[k]
                   - f_6 * pc_x[k] * ksg1_479[k];

        t_480[k] = pa_x[k] * ksg0_480[k]
                   + f_15 * ksf_320[k]
                   - f_6 * pc_x[k] * ksg1_480[k];

        t_481[k] = f_14 * ksf_250[k]
                   + f_3 * pc_y[k] * lsf_320[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pa_x, pc_x, pc_y, pc_z, ksg0_483, ksf_240, \
                         ksf_252, ksf_323, ksg1_483, lsf_320, lsf_322 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_15 * ksf_240[k]
                   + f_3 * pc_z[k] * lsf_320[k];

        t_483[k] = pa_x[k] * ksg0_483[k]
                   + f_8 * ksf_323[k]
                   - f_6 * pc_x[k] * ksg1_483[k];

        t_484[k] = f_14 * ksf_252[k]
                   + f_3 * pc_y[k] * lsf_322[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pa_x, pc_x, ksg0_485, ksf_325, ksf_326, \
                         ksf_327, ksf_328, ksg1_485, lsf_326, lsf_327, \
                         lsf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_x[k] * ksg0_485[k]
                   + f_8 * ksf_325[k]
                   - f_6 * pc_x[k] * ksg1_485[k];

        t_486[k] = f_7 * ksf_326[k]
                   + f_3 * pc_x[k] * lsf_326[k];

        t_487[k] = f_7 * ksf_327[k]
                   + f_3 * pc_x[k] * lsf_327[k];

        t_488[k] = f_7 * ksf_328[k]
                   + f_3 * pc_x[k] * lsf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, pa_x, pc_x, pc_z, ksg0_490, ksg0_492, \
                         ksf_246, ksf_329, ksg1_490, ksg1_492, lsf_326, \
                         lsf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_7 * ksf_329[k]
                   + f_3 * pc_x[k] * lsf_329[k];

        t_490[k] = pa_x[k] * ksg0_490[k]
                   - f_6 * pc_x[k] * ksg1_490[k];

        t_491[k] = f_15 * ksf_246[k]
                   + f_3 * pc_z[k] * lsf_326[k];

        t_492[k] = pa_x[k] * ksg0_492[k]
                   - f_6 * pc_x[k] * ksg1_492[k];
    }
}

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.0 / q;

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
    auto *t_614 = buffer.data(target + 614);
    auto *t_615 = buffer.data(target + 615);
    auto *t_616 = buffer.data(target + 616);
    auto *t_617 = buffer.data(target + 617);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksg0_405 = buffer.data(ksg0 + 405);
    const auto *ksg0_410 = buffer.data(ksg0 + 410);
    const auto *ksg0_420 = buffer.data(ksg0 + 420);
    const auto *ksg0_421 = buffer.data(ksg0 + 421);
    const auto *ksg0_423 = buffer.data(ksg0 + 423);
    const auto *ksg0_430 = buffer.data(ksg0 + 430);
    const auto *ksg0_432 = buffer.data(ksg0 + 432);
    const auto *ksg0_494 = buffer.data(ksg0 + 494);
    const auto *ksg0_495 = buffer.data(ksg0 + 495);
    const auto *ksg0_498 = buffer.data(ksg0 + 498);
    const auto *ksg0_500 = buffer.data(ksg0 + 500);
    const auto *ksg0_505 = buffer.data(ksg0 + 505);
    const auto *ksg0_507 = buffer.data(ksg0 + 507);
    const auto *ksg0_509 = buffer.data(ksg0 + 509);
    const auto *ksg0_513 = buffer.data(ksg0 + 513);
    const auto *ksg0_520 = buffer.data(ksg0 + 520);
    const auto *ksg0_522 = buffer.data(ksg0 + 522);
    const auto *ksg0_524 = buffer.data(ksg0 + 524);
    const auto *ksg0_525 = buffer.data(ksg0 + 525);
    const auto *ksg0_530 = buffer.data(ksg0 + 530);
    const auto *ksg0_535 = buffer.data(ksg0 + 535);
    const auto *ksg0_536 = buffer.data(ksg0 + 536);
    const auto *ksg0_537 = buffer.data(ksg0 + 537);
    const auto *ksg0_539 = buffer.data(ksg0 + 539);

    const auto *ksf_250 = buffer.data(ksf + 250);
    const auto *ksf_256 = buffer.data(ksf + 256);
    const auto *ksf_259 = buffer.data(ksf + 259);
    const auto *ksf_260 = buffer.data(ksf + 260);
    const auto *ksf_262 = buffer.data(ksf + 262);
    const auto *ksf_266 = buffer.data(ksf + 266);
    const auto *ksf_269 = buffer.data(ksf + 269);
    const auto *ksf_270 = buffer.data(ksf + 270);
    const auto *ksf_272 = buffer.data(ksf + 272);
    const auto *ksf_279 = buffer.data(ksf + 279);
    const auto *ksf_286 = buffer.data(ksf + 286);
    const auto *ksf_287 = buffer.data(ksf + 287);
    const auto *ksf_289 = buffer.data(ksf + 289);
    const auto *ksf_296 = buffer.data(ksf + 296);
    const auto *ksf_299 = buffer.data(ksf + 299);
    const auto *ksf_306 = buffer.data(ksf + 306);
    const auto *ksf_308 = buffer.data(ksf + 308);
    const auto *ksf_309 = buffer.data(ksf + 309);
    const auto *ksf_316 = buffer.data(ksf + 316);
    const auto *ksf_318 = buffer.data(ksf + 318);
    const auto *ksf_319 = buffer.data(ksf + 319);
    const auto *ksf_326 = buffer.data(ksf + 326);
    const auto *ksf_328 = buffer.data(ksf + 328);
    const auto *ksf_329 = buffer.data(ksf + 329);
    const auto *ksf_330 = buffer.data(ksf + 330);
    const auto *ksf_333 = buffer.data(ksf + 333);
    const auto *ksf_335 = buffer.data(ksf + 335);
    const auto *ksf_336 = buffer.data(ksf + 336);
    const auto *ksf_337 = buffer.data(ksf + 337);
    const auto *ksf_338 = buffer.data(ksf + 338);
    const auto *ksf_339 = buffer.data(ksf + 339);
    const auto *ksf_343 = buffer.data(ksf + 343);
    const auto *ksf_346 = buffer.data(ksf + 346);
    const auto *ksf_347 = buffer.data(ksf + 347);
    const auto *ksf_348 = buffer.data(ksf + 348);
    const auto *ksf_349 = buffer.data(ksf + 349);
    const auto *ksf_350 = buffer.data(ksf + 350);
    const auto *ksf_355 = buffer.data(ksf + 355);
    const auto *ksf_356 = buffer.data(ksf + 356);
    const auto *ksf_357 = buffer.data(ksf + 357);
    const auto *ksf_359 = buffer.data(ksf + 359);

    const auto *ksg1_405 = buffer.data(ksg1 + 405);
    const auto *ksg1_410 = buffer.data(ksg1 + 410);
    const auto *ksg1_420 = buffer.data(ksg1 + 420);
    const auto *ksg1_421 = buffer.data(ksg1 + 421);
    const auto *ksg1_423 = buffer.data(ksg1 + 423);
    const auto *ksg1_430 = buffer.data(ksg1 + 430);
    const auto *ksg1_432 = buffer.data(ksg1 + 432);
    const auto *ksg1_494 = buffer.data(ksg1 + 494);
    const auto *ksg1_495 = buffer.data(ksg1 + 495);
    const auto *ksg1_498 = buffer.data(ksg1 + 498);
    const auto *ksg1_500 = buffer.data(ksg1 + 500);
    const auto *ksg1_505 = buffer.data(ksg1 + 505);
    const auto *ksg1_507 = buffer.data(ksg1 + 507);
    const auto *ksg1_509 = buffer.data(ksg1 + 509);
    const auto *ksg1_513 = buffer.data(ksg1 + 513);
    const auto *ksg1_520 = buffer.data(ksg1 + 520);
    const auto *ksg1_522 = buffer.data(ksg1 + 522);
    const auto *ksg1_524 = buffer.data(ksg1 + 524);
    const auto *ksg1_525 = buffer.data(ksg1 + 525);
    const auto *ksg1_530 = buffer.data(ksg1 + 530);
    const auto *ksg1_535 = buffer.data(ksg1 + 535);
    const auto *ksg1_536 = buffer.data(ksg1 + 536);
    const auto *ksg1_537 = buffer.data(ksg1 + 537);
    const auto *ksg1_539 = buffer.data(ksg1 + 539);

    const auto *lsd0_210 = buffer.data(lsd0 + 210);
    const auto *lsd0_216 = buffer.data(lsd0 + 216);
    const auto *lsd0_217 = buffer.data(lsd0 + 217);
    const auto *lsd0_219 = buffer.data(lsd0 + 219);
    const auto *lsd0_221 = buffer.data(lsd0 + 221);
    const auto *lsd0_224 = buffer.data(lsd0 + 224);
    const auto *lsd0_226 = buffer.data(lsd0 + 226);
    const auto *lsd0_227 = buffer.data(lsd0 + 227);
    const auto *lsd0_228 = buffer.data(lsd0 + 228);
    const auto *lsd0_229 = buffer.data(lsd0 + 229);
    const auto *lsd0_230 = buffer.data(lsd0 + 230);
    const auto *lsd0_231 = buffer.data(lsd0 + 231);
    const auto *lsd0_232 = buffer.data(lsd0 + 232);
    const auto *lsd0_233 = buffer.data(lsd0 + 233);
    const auto *lsd0_234 = buffer.data(lsd0 + 234);
    const auto *lsd0_235 = buffer.data(lsd0 + 235);
    const auto *lsd0_236 = buffer.data(lsd0 + 236);
    const auto *lsd0_237 = buffer.data(lsd0 + 237);
    const auto *lsd0_238 = buffer.data(lsd0 + 238);
    const auto *lsd0_239 = buffer.data(lsd0 + 239);
    const auto *lsd0_240 = buffer.data(lsd0 + 240);
    const auto *lsd0_241 = buffer.data(lsd0 + 241);
    const auto *lsd0_242 = buffer.data(lsd0 + 242);
    const auto *lsd0_243 = buffer.data(lsd0 + 243);
    const auto *lsd0_244 = buffer.data(lsd0 + 244);
    const auto *lsd0_245 = buffer.data(lsd0 + 245);
    const auto *lsd0_246 = buffer.data(lsd0 + 246);
    const auto *lsd0_247 = buffer.data(lsd0 + 247);
    const auto *lsd0_248 = buffer.data(lsd0 + 248);

    const auto *lsd1_210 = buffer.data(lsd1 + 210);
    const auto *lsd1_216 = buffer.data(lsd1 + 216);
    const auto *lsd1_217 = buffer.data(lsd1 + 217);
    const auto *lsd1_219 = buffer.data(lsd1 + 219);
    const auto *lsd1_221 = buffer.data(lsd1 + 221);
    const auto *lsd1_224 = buffer.data(lsd1 + 224);
    const auto *lsd1_226 = buffer.data(lsd1 + 226);
    const auto *lsd1_227 = buffer.data(lsd1 + 227);
    const auto *lsd1_228 = buffer.data(lsd1 + 228);
    const auto *lsd1_229 = buffer.data(lsd1 + 229);
    const auto *lsd1_230 = buffer.data(lsd1 + 230);
    const auto *lsd1_231 = buffer.data(lsd1 + 231);
    const auto *lsd1_232 = buffer.data(lsd1 + 232);
    const auto *lsd1_233 = buffer.data(lsd1 + 233);
    const auto *lsd1_234 = buffer.data(lsd1 + 234);
    const auto *lsd1_235 = buffer.data(lsd1 + 235);
    const auto *lsd1_236 = buffer.data(lsd1 + 236);
    const auto *lsd1_237 = buffer.data(lsd1 + 237);
    const auto *lsd1_238 = buffer.data(lsd1 + 238);
    const auto *lsd1_239 = buffer.data(lsd1 + 239);
    const auto *lsd1_240 = buffer.data(lsd1 + 240);
    const auto *lsd1_241 = buffer.data(lsd1 + 241);
    const auto *lsd1_242 = buffer.data(lsd1 + 242);
    const auto *lsd1_243 = buffer.data(lsd1 + 243);
    const auto *lsd1_244 = buffer.data(lsd1 + 244);
    const auto *lsd1_245 = buffer.data(lsd1 + 245);
    const auto *lsd1_246 = buffer.data(lsd1 + 246);
    const auto *lsd1_247 = buffer.data(lsd1 + 247);
    const auto *lsd1_248 = buffer.data(lsd1 + 248);

    const auto *lsf_329 = buffer.data(lsf + 329);
    const auto *lsf_330 = buffer.data(lsf + 330);
    const auto *lsf_332 = buffer.data(lsf + 332);
    const auto *lsf_336 = buffer.data(lsf + 336);
    const auto *lsf_337 = buffer.data(lsf + 337);
    const auto *lsf_338 = buffer.data(lsf + 338);
    const auto *lsf_339 = buffer.data(lsf + 339);
    const auto *lsf_340 = buffer.data(lsf + 340);
    const auto *lsf_342 = buffer.data(lsf + 342);
    const auto *lsf_346 = buffer.data(lsf + 346);
    const auto *lsf_347 = buffer.data(lsf + 347);
    const auto *lsf_348 = buffer.data(lsf + 348);
    const auto *lsf_349 = buffer.data(lsf + 349);
    const auto *lsf_350 = buffer.data(lsf + 350);
    const auto *lsf_351 = buffer.data(lsf + 351);
    const auto *lsf_352 = buffer.data(lsf + 352);
    const auto *lsf_355 = buffer.data(lsf + 355);
    const auto *lsf_356 = buffer.data(lsf + 356);
    const auto *lsf_357 = buffer.data(lsf + 357);
    const auto *lsf_359 = buffer.data(lsf + 359);
    const auto *lsf_360 = buffer.data(lsf + 360);
    const auto *lsf_361 = buffer.data(lsf + 361);
    const auto *lsf_363 = buffer.data(lsf + 363);
    const auto *lsf_365 = buffer.data(lsf + 365);
    const auto *lsf_366 = buffer.data(lsf + 366);
    const auto *lsf_367 = buffer.data(lsf + 367);
    const auto *lsf_368 = buffer.data(lsf + 368);
    const auto *lsf_369 = buffer.data(lsf + 369);
    const auto *lsf_372 = buffer.data(lsf + 372);
    const auto *lsf_374 = buffer.data(lsf + 374);
    const auto *lsf_375 = buffer.data(lsf + 375);
    const auto *lsf_376 = buffer.data(lsf + 376);
    const auto *lsf_377 = buffer.data(lsf + 377);
    const auto *lsf_378 = buffer.data(lsf + 378);
    const auto *lsf_379 = buffer.data(lsf + 379);
    const auto *lsf_380 = buffer.data(lsf + 380);
    const auto *lsf_381 = buffer.data(lsf + 381);
    const auto *lsf_382 = buffer.data(lsf + 382);
    const auto *lsf_383 = buffer.data(lsf + 383);
    const auto *lsf_384 = buffer.data(lsf + 384);
    const auto *lsf_385 = buffer.data(lsf + 385);
    const auto *lsf_386 = buffer.data(lsf + 386);
    const auto *lsf_387 = buffer.data(lsf + 387);
    const auto *lsf_388 = buffer.data(lsf + 388);
    const auto *lsf_389 = buffer.data(lsf + 389);
    const auto *lsf_390 = buffer.data(lsf + 390);
    const auto *lsf_391 = buffer.data(lsf + 391);
    const auto *lsf_392 = buffer.data(lsf + 392);
    const auto *lsf_393 = buffer.data(lsf + 393);
    const auto *lsf_394 = buffer.data(lsf + 394);
    const auto *lsf_395 = buffer.data(lsf + 395);
    const auto *lsf_396 = buffer.data(lsf + 396);
    const auto *lsf_397 = buffer.data(lsf + 397);
    const auto *lsf_398 = buffer.data(lsf + 398);
    const auto *lsf_399 = buffer.data(lsf + 399);
    const auto *lsf_400 = buffer.data(lsf + 400);
    const auto *lsf_401 = buffer.data(lsf + 401);
    const auto *lsf_402 = buffer.data(lsf + 402);
    const auto *lsf_403 = buffer.data(lsf + 403);
    const auto *lsf_404 = buffer.data(lsf + 404);
    const auto *lsf_405 = buffer.data(lsf + 405);
    const auto *lsf_406 = buffer.data(lsf + 406);
    const auto *lsf_407 = buffer.data(lsf + 407);
    const auto *lsf_408 = buffer.data(lsf + 408);
    const auto *lsf_409 = buffer.data(lsf + 409);
    const auto *lsf_410 = buffer.data(lsf + 410);
    const auto *lsf_411 = buffer.data(lsf + 411);
    const auto *lsf_412 = buffer.data(lsf + 412);

#pragma omp simd aligned(t_493, t_494, t_495, t_496, pa_x, pc_x, pc_y, ksg0_494, ksg0_495, \
                         ksf_259, ksf_260, ksf_330, ksg1_494, ksg1_495, lsf_329, \
                         lsf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_493[k] = f_14 * ksf_259[k]
                   + f_3 * pc_y[k] * lsf_329[k];

        t_494[k] = pa_x[k] * ksg0_494[k]
                   - f_6 * pc_x[k] * ksg1_494[k];

        t_495[k] = pa_x[k] * ksg0_495[k]
                   + f_15 * ksf_330[k]
                   - f_6 * pc_x[k] * ksg1_495[k];

        t_496[k] = f_8 * ksf_260[k]
                   + f_3 * pc_y[k] * lsf_330[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, pa_x, pc_x, pc_y, pc_z, ksg0_498, ksf_250, \
                         ksf_262, ksf_333, ksg1_498, lsf_330, lsf_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = f_13 * ksf_250[k]
                   + f_3 * pc_z[k] * lsf_330[k];

        t_498[k] = pa_x[k] * ksg0_498[k]
                   + f_8 * ksf_333[k]
                   - f_6 * pc_x[k] * ksg1_498[k];

        t_499[k] = f_8 * ksf_262[k]
                   + f_3 * pc_y[k] * lsf_332[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pa_x, pc_x, ksg0_500, ksf_335, ksf_336, \
                         ksf_337, ksf_338, ksg1_500, lsf_336, lsf_337, \
                         lsf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = pa_x[k] * ksg0_500[k]
                   + f_8 * ksf_335[k]
                   - f_6 * pc_x[k] * ksg1_500[k];

        t_501[k] = f_7 * ksf_336[k]
                   + f_3 * pc_x[k] * lsf_336[k];

        t_502[k] = f_7 * ksf_337[k]
                   + f_3 * pc_x[k] * lsf_337[k];

        t_503[k] = f_7 * ksf_338[k]
                   + f_3 * pc_x[k] * lsf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pa_x, pc_x, pc_z, ksg0_505, ksg0_507, \
                         ksf_256, ksf_339, ksg1_505, ksg1_507, lsf_336, \
                         lsf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_7 * ksf_339[k]
                   + f_3 * pc_x[k] * lsf_339[k];

        t_505[k] = pa_x[k] * ksg0_505[k]
                   - f_6 * pc_x[k] * ksg1_505[k];

        t_506[k] = f_13 * ksf_256[k]
                   + f_3 * pc_z[k] * lsf_336[k];

        t_507[k] = pa_x[k] * ksg0_507[k]
                   - f_6 * pc_x[k] * ksg1_507[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_x, pa_y, pc_x, pc_y, ksg0_405, \
                         ksg0_509, ksf_269, ksf_270, ksg1_405, ksg1_509, lsf_339, \
                         lsf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * ksf_269[k]
                   + f_3 * pc_y[k] * lsf_339[k];

        t_509[k] = pa_x[k] * ksg0_509[k]
                   - f_6 * pc_x[k] * ksg1_509[k];

        t_510[k] = pa_y[k] * ksg0_405[k]
                   - f_6 * pc_y[k] * ksg1_405[k];

        t_511[k] = f_7 * ksf_270[k]
                   + f_3 * pc_y[k] * lsf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_x, pc_x, pc_y, pc_z, ksg0_513, ksf_260, \
                         ksf_272, ksf_343, ksg1_513, lsf_340, lsf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * ksf_260[k]
                   + f_3 * pc_z[k] * lsf_340[k];

        t_513[k] = pa_x[k] * ksg0_513[k]
                   + f_8 * ksf_343[k]
                   - f_6 * pc_x[k] * ksg1_513[k];

        t_514[k] = f_7 * ksf_272[k]
                   + f_3 * pc_y[k] * lsf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pa_y, pc_x, pc_y, ksg0_410, ksf_346, \
                         ksf_347, ksf_348, ksg1_410, lsf_346, lsf_347, \
                         lsf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pa_y[k] * ksg0_410[k]
                   - f_6 * pc_y[k] * ksg1_410[k];

        t_516[k] = f_7 * ksf_346[k]
                   + f_3 * pc_x[k] * lsf_346[k];

        t_517[k] = f_7 * ksf_347[k]
                   + f_3 * pc_x[k] * lsf_347[k];

        t_518[k] = f_7 * ksf_348[k]
                   + f_3 * pc_x[k] * lsf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pa_x, pc_x, pc_z, ksg0_520, ksg0_522, \
                         ksf_266, ksf_349, ksg1_520, ksg1_522, lsf_346, \
                         lsf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_7 * ksf_349[k]
                   + f_3 * pc_x[k] * lsf_349[k];

        t_520[k] = pa_x[k] * ksg0_520[k]
                   - f_6 * pc_x[k] * ksg1_520[k];

        t_521[k] = f_12 * ksf_266[k]
                   + f_3 * pc_z[k] * lsf_346[k];

        t_522[k] = pa_x[k] * ksg0_522[k]
                   - f_6 * pc_x[k] * ksg1_522[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_x, pc_x, pc_y, ksg0_524, ksg0_525, \
                         ksf_279, ksf_350, ksg1_524, ksg1_525, lsf_349, \
                         lsf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * ksf_279[k]
                   + f_3 * pc_y[k] * lsf_349[k];

        t_524[k] = pa_x[k] * ksg0_524[k]
                   - f_6 * pc_x[k] * ksg1_524[k];

        t_525[k] = pa_x[k] * ksg0_525[k]
                   + f_15 * ksf_350[k]
                   - f_6 * pc_x[k] * ksg1_525[k];

        t_526[k] = f_3 * pc_y[k] * lsf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, pc_z, ksf_270, lsd0_210, lsd1_210, \
                         lsf_350, lsf_351, lsf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_9 * ksf_270[k]
                   + f_3 * pc_z[k] * lsf_350[k];

        t_528[k] = f_4 * lsd0_210[k]
                   - f_5 * lsd1_210[k]
                   + f_3 * pc_y[k] * lsf_351[k];

        t_529[k] = f_3 * pc_y[k] * lsf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_x, pc_x, pc_y, ksg0_530, ksf_355, \
                         ksf_356, ksf_357, ksg1_530, lsf_355, lsf_356, \
                         lsf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = pa_x[k] * ksg0_530[k]
                   + f_8 * ksf_355[k]
                   - f_6 * pc_x[k] * ksg1_530[k];

        t_531[k] = f_7 * ksf_356[k]
                   + f_3 * pc_x[k] * lsf_356[k];

        t_532[k] = f_7 * ksf_357[k]
                   + f_3 * pc_x[k] * lsf_357[k];

        t_533[k] = f_3 * pc_y[k] * lsf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, t_538, pa_x, pc_x, pc_y, ksg0_535, \
                         ksg0_536, ksg0_537, ksf_359, ksg1_535, ksg1_536, ksg1_537, \
                         lsf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_7 * ksf_359[k]
                   + f_3 * pc_x[k] * lsf_359[k];

        t_535[k] = pa_x[k] * ksg0_535[k]
                   - f_6 * pc_x[k] * ksg1_535[k];

        t_536[k] = pa_x[k] * ksg0_536[k]
                   - f_6 * pc_x[k] * ksg1_536[k];

        t_537[k] = pa_x[k] * ksg0_537[k]
                   - f_6 * pc_x[k] * ksg1_537[k];

        t_538[k] = f_3 * pc_y[k] * lsf_359[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, t_542, pa_x, pc_x, pc_z, ksg0_539, ksg1_539, \
                         lsd0_216, lsd0_217, lsd1_216, lsd1_217, lsf_360, \
                         lsf_361 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = pa_x[k] * ksg0_539[k]
                   - f_6 * pc_x[k] * ksg1_539[k];

        t_540[k] = f_1 * lsd0_216[k]
                   - f_2 * lsd1_216[k]
                   + f_3 * pc_x[k] * lsf_360[k];

        t_541[k] = f_10 * lsd0_217[k]
                   - f_11 * lsd1_217[k]
                   + f_3 * pc_x[k] * lsf_361[k];

        t_542[k] = f_3 * pc_z[k] * lsf_360[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, t_546, t_547, pc_x, pc_z, lsd0_219, lsd0_221, \
                         lsd1_219, lsd1_221, lsf_361, lsf_363, lsf_365, lsf_366, \
                         lsf_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_4 * lsd0_219[k]
                   - f_5 * lsd1_219[k]
                   + f_3 * pc_x[k] * lsf_363[k];

        t_544[k] = f_3 * pc_z[k] * lsf_361[k];

        t_545[k] = f_4 * lsd0_221[k]
                   - f_5 * lsd1_221[k]
                   + f_3 * pc_x[k] * lsf_365[k];

        t_546[k] = f_3 * pc_x[k] * lsf_366[k];

        t_547[k] = f_3 * pc_x[k] * lsf_367[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, pc_y, pc_z, ksf_286, \
                         lsd0_219, lsd1_219, lsf_366, lsf_367, lsf_368, \
                         lsf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_3 * pc_x[k] * lsf_368[k];

        t_549[k] = f_3 * pc_x[k] * lsf_369[k];

        t_550[k] = f_0 * ksf_286[k]
                   + f_1 * lsd0_219[k]
                   - f_2 * lsd1_219[k]
                   + f_3 * pc_y[k] * lsf_366[k];

        t_551[k] = f_3 * pc_z[k] * lsf_366[k];

        t_552[k] = f_4 * lsd0_219[k]
                   - f_5 * lsd1_219[k]
                   + f_3 * pc_z[k] * lsf_367[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, t_556, pa_z, pc_y, pc_z, ksg0_420, ksg0_421, \
                         ksf_289, ksg1_420, ksg1_421, lsd0_221, lsd1_221, \
                         lsf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_0 * ksf_289[k]
                   + f_3 * pc_y[k] * lsf_369[k];

        t_554[k] = f_1 * lsd0_221[k]
                   - f_2 * lsd1_221[k]
                   + f_3 * pc_z[k] * lsf_369[k];

        t_555[k] = pa_z[k] * ksg0_420[k]
                   - f_6 * pc_z[k] * ksg1_420[k];

        t_556[k] = pa_z[k] * ksg0_421[k]
                   - f_6 * pc_z[k] * ksg1_421[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pa_z, pc_x, pc_z, ksg0_423, ksg1_423, lsd0_224, \
                         lsd0_226, lsd1_224, lsd1_226, lsf_372, \
                         lsf_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_10 * lsd0_224[k]
                   - f_11 * lsd1_224[k]
                   + f_3 * pc_x[k] * lsf_372[k];

        t_558[k] = pa_z[k] * ksg0_423[k]
                   - f_6 * pc_z[k] * ksg1_423[k];

        t_559[k] = f_4 * lsd0_226[k]
                   - f_5 * lsd1_226[k]
                   + f_3 * pc_x[k] * lsf_374[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, pc_x, lsd0_227, lsd1_227, lsf_375, \
                         lsf_376, lsf_377, lsf_378, lsf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_4 * lsd0_227[k]
                   - f_5 * lsd1_227[k]
                   + f_3 * pc_x[k] * lsf_375[k];

        t_561[k] = f_3 * pc_x[k] * lsf_376[k];

        t_562[k] = f_3 * pc_x[k] * lsf_377[k];

        t_563[k] = f_3 * pc_x[k] * lsf_378[k];

        t_564[k] = f_3 * pc_x[k] * lsf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pa_z, pc_y, pc_z, ksg0_430, ksg0_432, \
                         ksf_286, ksf_287, ksf_299, ksg1_430, ksg1_432, lsf_376, \
                         lsf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pa_z[k] * ksg0_430[k]
                   - f_6 * pc_z[k] * ksg1_430[k];

        t_566[k] = f_7 * ksf_286[k]
                   + f_3 * pc_z[k] * lsf_376[k];

        t_567[k] = pa_z[k] * ksg0_432[k]
                   + f_8 * ksf_287[k]
                   - f_6 * pc_z[k] * ksg1_432[k];

        t_568[k] = f_9 * ksf_299[k]
                   + f_3 * pc_y[k] * lsf_379[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_z, ksf_289, lsd0_227, lsd0_228, \
                         lsd0_229, lsd1_227, lsd1_228, lsd1_229, lsf_379, lsf_380, \
                         lsf_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_7 * ksf_289[k]
                   + f_1 * lsd0_227[k]
                   - f_2 * lsd1_227[k]
                   + f_3 * pc_z[k] * lsf_379[k];

        t_570[k] = f_1 * lsd0_228[k]
                   - f_2 * lsd1_228[k]
                   + f_3 * pc_x[k] * lsf_380[k];

        t_571[k] = f_10 * lsd0_229[k]
                   - f_11 * lsd1_229[k]
                   + f_3 * pc_x[k] * lsf_381[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, lsd0_230, lsd0_231, lsd0_232, lsd1_230, \
                         lsd1_231, lsd1_232, lsf_382, lsf_383, \
                         lsf_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_10 * lsd0_230[k]
                   - f_11 * lsd1_230[k]
                   + f_3 * pc_x[k] * lsf_382[k];

        t_573[k] = f_4 * lsd0_231[k]
                   - f_5 * lsd1_231[k]
                   + f_3 * pc_x[k] * lsf_383[k];

        t_574[k] = f_4 * lsd0_232[k]
                   - f_5 * lsd1_232[k]
                   + f_3 * pc_x[k] * lsf_384[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, t_579, pc_x, lsd0_233, lsd1_233, lsf_385, \
                         lsf_386, lsf_387, lsf_388, lsf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_4 * lsd0_233[k]
                   - f_5 * lsd1_233[k]
                   + f_3 * pc_x[k] * lsf_385[k];

        t_576[k] = f_3 * pc_x[k] * lsf_386[k];

        t_577[k] = f_3 * pc_x[k] * lsf_387[k];

        t_578[k] = f_3 * pc_x[k] * lsf_388[k];

        t_579[k] = f_3 * pc_x[k] * lsf_389[k];
    }

#pragma omp simd aligned(t_580, t_581, t_582, pc_y, pc_z, ksf_296, ksf_306, ksf_308, lsd0_231, \
                         lsd0_233, lsd1_231, lsd1_233, lsf_386, \
                         lsf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_580[k] = f_12 * ksf_306[k]
                   + f_1 * lsd0_231[k]
                   - f_2 * lsd1_231[k]
                   + f_3 * pc_y[k] * lsf_386[k];

        t_581[k] = f_8 * ksf_296[k]
                   + f_3 * pc_z[k] * lsf_386[k];

        t_582[k] = f_12 * ksf_308[k]
                   + f_4 * lsd0_233[k]
                   - f_5 * lsd1_233[k]
                   + f_3 * pc_y[k] * lsf_388[k];
    }

#pragma omp simd aligned(t_583, t_584, t_585, pc_x, pc_y, pc_z, ksf_299, ksf_309, lsd0_233, \
                         lsd0_234, lsd1_233, lsd1_234, lsf_389, \
                         lsf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_583[k] = f_12 * ksf_309[k]
                   + f_3 * pc_y[k] * lsf_389[k];

        t_584[k] = f_8 * ksf_299[k]
                   + f_1 * lsd0_233[k]
                   - f_2 * lsd1_233[k]
                   + f_3 * pc_z[k] * lsf_389[k];

        t_585[k] = f_1 * lsd0_234[k]
                   - f_2 * lsd1_234[k]
                   + f_3 * pc_x[k] * lsf_390[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, pc_x, lsd0_235, lsd0_236, lsd0_237, lsd1_235, \
                         lsd1_236, lsd1_237, lsf_391, lsf_392, \
                         lsf_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_10 * lsd0_235[k]
                   - f_11 * lsd1_235[k]
                   + f_3 * pc_x[k] * lsf_391[k];

        t_587[k] = f_10 * lsd0_236[k]
                   - f_11 * lsd1_236[k]
                   + f_3 * pc_x[k] * lsf_392[k];

        t_588[k] = f_4 * lsd0_237[k]
                   - f_5 * lsd1_237[k]
                   + f_3 * pc_x[k] * lsf_393[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, t_593, pc_x, lsd0_238, lsd0_239, \
                         lsd1_238, lsd1_239, lsf_394, lsf_395, lsf_396, lsf_397, \
                         lsf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_4 * lsd0_238[k]
                   - f_5 * lsd1_238[k]
                   + f_3 * pc_x[k] * lsf_394[k];

        t_590[k] = f_4 * lsd0_239[k]
                   - f_5 * lsd1_239[k]
                   + f_3 * pc_x[k] * lsf_395[k];

        t_591[k] = f_3 * pc_x[k] * lsf_396[k];

        t_592[k] = f_3 * pc_x[k] * lsf_397[k];

        t_593[k] = f_3 * pc_x[k] * lsf_398[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pc_x, pc_y, pc_z, ksf_306, ksf_316, lsd0_237, \
                         lsd1_237, lsf_396, lsf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_3 * pc_x[k] * lsf_399[k];

        t_595[k] = f_13 * ksf_316[k]
                   + f_1 * lsd0_237[k]
                   - f_2 * lsd1_237[k]
                   + f_3 * pc_y[k] * lsf_396[k];

        t_596[k] = f_14 * ksf_306[k]
                   + f_3 * pc_z[k] * lsf_396[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pc_y, pc_z, ksf_309, ksf_318, ksf_319, lsd0_239, \
                         lsd1_239, lsf_398, lsf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = f_13 * ksf_318[k]
                   + f_4 * lsd0_239[k]
                   - f_5 * lsd1_239[k]
                   + f_3 * pc_y[k] * lsf_398[k];

        t_598[k] = f_13 * ksf_319[k]
                   + f_3 * pc_y[k] * lsf_399[k];

        t_599[k] = f_14 * ksf_309[k]
                   + f_1 * lsd0_239[k]
                   - f_2 * lsd1_239[k]
                   + f_3 * pc_z[k] * lsf_399[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pc_x, lsd0_240, lsd0_241, lsd0_242, lsd1_240, \
                         lsd1_241, lsd1_242, lsf_400, lsf_401, \
                         lsf_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_1 * lsd0_240[k]
                   - f_2 * lsd1_240[k]
                   + f_3 * pc_x[k] * lsf_400[k];

        t_601[k] = f_10 * lsd0_241[k]
                   - f_11 * lsd1_241[k]
                   + f_3 * pc_x[k] * lsf_401[k];

        t_602[k] = f_10 * lsd0_242[k]
                   - f_11 * lsd1_242[k]
                   + f_3 * pc_x[k] * lsf_402[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pc_x, lsd0_243, lsd0_244, lsd0_245, \
                         lsd1_243, lsd1_244, lsd1_245, lsf_403, lsf_404, lsf_405, \
                         lsf_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_4 * lsd0_243[k]
                   - f_5 * lsd1_243[k]
                   + f_3 * pc_x[k] * lsf_403[k];

        t_604[k] = f_4 * lsd0_244[k]
                   - f_5 * lsd1_244[k]
                   + f_3 * pc_x[k] * lsf_404[k];

        t_605[k] = f_4 * lsd0_245[k]
                   - f_5 * lsd1_245[k]
                   + f_3 * pc_x[k] * lsf_405[k];

        t_606[k] = f_3 * pc_x[k] * lsf_406[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, t_610, t_611, pc_x, pc_y, pc_z, ksf_316, \
                         ksf_326, lsd0_243, lsd1_243, lsf_406, lsf_407, lsf_408, \
                         lsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_3 * pc_x[k] * lsf_407[k];

        t_608[k] = f_3 * pc_x[k] * lsf_408[k];

        t_609[k] = f_3 * pc_x[k] * lsf_409[k];

        t_610[k] = f_15 * ksf_326[k]
                   + f_1 * lsd0_243[k]
                   - f_2 * lsd1_243[k]
                   + f_3 * pc_y[k] * lsf_406[k];

        t_611[k] = f_15 * ksf_316[k]
                   + f_3 * pc_z[k] * lsf_406[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pc_y, pc_z, ksf_319, ksf_328, ksf_329, lsd0_245, \
                         lsd1_245, lsf_408, lsf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_15 * ksf_328[k]
                   + f_4 * lsd0_245[k]
                   - f_5 * lsd1_245[k]
                   + f_3 * pc_y[k] * lsf_408[k];

        t_613[k] = f_15 * ksf_329[k]
                   + f_3 * pc_y[k] * lsf_409[k];

        t_614[k] = f_15 * ksf_319[k]
                   + f_1 * lsd0_245[k]
                   - f_2 * lsd1_245[k]
                   + f_3 * pc_z[k] * lsf_409[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, lsd0_246, lsd0_247, lsd0_248, lsd1_246, \
                         lsd1_247, lsd1_248, lsf_410, lsf_411, \
                         lsf_412 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_1 * lsd0_246[k]
                   - f_2 * lsd1_246[k]
                   + f_3 * pc_x[k] * lsf_410[k];

        t_616[k] = f_10 * lsd0_247[k]
                   - f_11 * lsd1_247[k]
                   + f_3 * pc_x[k] * lsf_411[k];

        t_617[k] = f_10 * lsd0_248[k]
                   - f_11 * lsd1_248[k]
                   + f_3 * pc_x[k] * lsf_412[k];
    }
}

static auto
compute_prim_lsg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksg0,
                                                          const size_t ksf, const size_t ksg1,
                                                          const size_t lsd0, const size_t lsd1,
                                                          const size_t lsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.5 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 3.0 / q;
    const auto f_13 = 2.5 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksg0_525 = buffer.data(ksg0 + 525);
    const auto *ksg0_527 = buffer.data(ksg0 + 527);
    const auto *ksg0_530 = buffer.data(ksg0 + 530);
    const auto *ksg0_535 = buffer.data(ksg0 + 535);
    const auto *ksg0_537 = buffer.data(ksg0 + 537);
    const auto *ksg0_539 = buffer.data(ksg0 + 539);

    const auto *ksf_326 = buffer.data(ksf + 326);
    const auto *ksf_329 = buffer.data(ksf + 329);
    const auto *ksf_336 = buffer.data(ksf + 336);
    const auto *ksf_338 = buffer.data(ksf + 338);
    const auto *ksf_339 = buffer.data(ksf + 339);
    const auto *ksf_346 = buffer.data(ksf + 346);
    const auto *ksf_348 = buffer.data(ksf + 348);
    const auto *ksf_349 = buffer.data(ksf + 349);
    const auto *ksf_356 = buffer.data(ksf + 356);
    const auto *ksf_358 = buffer.data(ksf + 358);
    const auto *ksf_359 = buffer.data(ksf + 359);

    const auto *ksg1_525 = buffer.data(ksg1 + 525);
    const auto *ksg1_527 = buffer.data(ksg1 + 527);
    const auto *ksg1_530 = buffer.data(ksg1 + 530);
    const auto *ksg1_535 = buffer.data(ksg1 + 535);
    const auto *ksg1_537 = buffer.data(ksg1 + 537);
    const auto *ksg1_539 = buffer.data(ksg1 + 539);

    const auto *lsd0_249 = buffer.data(lsd0 + 249);
    const auto *lsd0_250 = buffer.data(lsd0 + 250);
    const auto *lsd0_251 = buffer.data(lsd0 + 251);
    const auto *lsd0_252 = buffer.data(lsd0 + 252);
    const auto *lsd0_253 = buffer.data(lsd0 + 253);
    const auto *lsd0_254 = buffer.data(lsd0 + 254);
    const auto *lsd0_255 = buffer.data(lsd0 + 255);
    const auto *lsd0_256 = buffer.data(lsd0 + 256);
    const auto *lsd0_257 = buffer.data(lsd0 + 257);
    const auto *lsd0_259 = buffer.data(lsd0 + 259);
    const auto *lsd0_261 = buffer.data(lsd0 + 261);
    const auto *lsd0_262 = buffer.data(lsd0 + 262);
    const auto *lsd0_264 = buffer.data(lsd0 + 264);
    const auto *lsd0_266 = buffer.data(lsd0 + 266);
    const auto *lsd0_267 = buffer.data(lsd0 + 267);
    const auto *lsd0_268 = buffer.data(lsd0 + 268);
    const auto *lsd0_269 = buffer.data(lsd0 + 269);

    const auto *lsd1_249 = buffer.data(lsd1 + 249);
    const auto *lsd1_250 = buffer.data(lsd1 + 250);
    const auto *lsd1_251 = buffer.data(lsd1 + 251);
    const auto *lsd1_252 = buffer.data(lsd1 + 252);
    const auto *lsd1_253 = buffer.data(lsd1 + 253);
    const auto *lsd1_254 = buffer.data(lsd1 + 254);
    const auto *lsd1_255 = buffer.data(lsd1 + 255);
    const auto *lsd1_256 = buffer.data(lsd1 + 256);
    const auto *lsd1_257 = buffer.data(lsd1 + 257);
    const auto *lsd1_259 = buffer.data(lsd1 + 259);
    const auto *lsd1_261 = buffer.data(lsd1 + 261);
    const auto *lsd1_262 = buffer.data(lsd1 + 262);
    const auto *lsd1_264 = buffer.data(lsd1 + 264);
    const auto *lsd1_266 = buffer.data(lsd1 + 266);
    const auto *lsd1_267 = buffer.data(lsd1 + 267);
    const auto *lsd1_268 = buffer.data(lsd1 + 268);
    const auto *lsd1_269 = buffer.data(lsd1 + 269);

    const auto *lsf_413 = buffer.data(lsf + 413);
    const auto *lsf_414 = buffer.data(lsf + 414);
    const auto *lsf_415 = buffer.data(lsf + 415);
    const auto *lsf_416 = buffer.data(lsf + 416);
    const auto *lsf_417 = buffer.data(lsf + 417);
    const auto *lsf_418 = buffer.data(lsf + 418);
    const auto *lsf_419 = buffer.data(lsf + 419);
    const auto *lsf_420 = buffer.data(lsf + 420);
    const auto *lsf_421 = buffer.data(lsf + 421);
    const auto *lsf_422 = buffer.data(lsf + 422);
    const auto *lsf_423 = buffer.data(lsf + 423);
    const auto *lsf_424 = buffer.data(lsf + 424);
    const auto *lsf_425 = buffer.data(lsf + 425);
    const auto *lsf_426 = buffer.data(lsf + 426);
    const auto *lsf_427 = buffer.data(lsf + 427);
    const auto *lsf_428 = buffer.data(lsf + 428);
    const auto *lsf_429 = buffer.data(lsf + 429);
    const auto *lsf_431 = buffer.data(lsf + 431);
    const auto *lsf_433 = buffer.data(lsf + 433);
    const auto *lsf_434 = buffer.data(lsf + 434);
    const auto *lsf_436 = buffer.data(lsf + 436);
    const auto *lsf_437 = buffer.data(lsf + 437);
    const auto *lsf_438 = buffer.data(lsf + 438);
    const auto *lsf_439 = buffer.data(lsf + 439);
    const auto *lsf_440 = buffer.data(lsf + 440);
    const auto *lsf_442 = buffer.data(lsf + 442);
    const auto *lsf_443 = buffer.data(lsf + 443);
    const auto *lsf_445 = buffer.data(lsf + 445);
    const auto *lsf_446 = buffer.data(lsf + 446);
    const auto *lsf_447 = buffer.data(lsf + 447);
    const auto *lsf_448 = buffer.data(lsf + 448);
    const auto *lsf_449 = buffer.data(lsf + 449);

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pc_x, lsd0_249, lsd0_250, lsd0_251, \
                         lsd1_249, lsd1_250, lsd1_251, lsf_413, lsf_414, lsf_415, \
                         lsf_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_4 * lsd0_249[k]
                   - f_5 * lsd1_249[k]
                   + f_3 * pc_x[k] * lsf_413[k];

        t_619[k] = f_4 * lsd0_250[k]
                   - f_5 * lsd1_250[k]
                   + f_3 * pc_x[k] * lsf_414[k];

        t_620[k] = f_4 * lsd0_251[k]
                   - f_5 * lsd1_251[k]
                   + f_3 * pc_x[k] * lsf_415[k];

        t_621[k] = f_3 * pc_x[k] * lsf_416[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, t_626, pc_x, pc_y, pc_z, ksf_326, \
                         ksf_336, lsd0_249, lsd1_249, lsf_416, lsf_417, lsf_418, \
                         lsf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_3 * pc_x[k] * lsf_417[k];

        t_623[k] = f_3 * pc_x[k] * lsf_418[k];

        t_624[k] = f_3 * pc_x[k] * lsf_419[k];

        t_625[k] = f_14 * ksf_336[k]
                   + f_1 * lsd0_249[k]
                   - f_2 * lsd1_249[k]
                   + f_3 * pc_y[k] * lsf_416[k];

        t_626[k] = f_13 * ksf_326[k]
                   + f_3 * pc_z[k] * lsf_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pc_y, pc_z, ksf_329, ksf_338, ksf_339, lsd0_251, \
                         lsd1_251, lsf_418, lsf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_14 * ksf_338[k]
                   + f_4 * lsd0_251[k]
                   - f_5 * lsd1_251[k]
                   + f_3 * pc_y[k] * lsf_418[k];

        t_628[k] = f_14 * ksf_339[k]
                   + f_3 * pc_y[k] * lsf_419[k];

        t_629[k] = f_13 * ksf_329[k]
                   + f_1 * lsd0_251[k]
                   - f_2 * lsd1_251[k]
                   + f_3 * pc_z[k] * lsf_419[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pc_x, lsd0_252, lsd0_253, lsd0_254, lsd1_252, \
                         lsd1_253, lsd1_254, lsf_420, lsf_421, \
                         lsf_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_1 * lsd0_252[k]
                   - f_2 * lsd1_252[k]
                   + f_3 * pc_x[k] * lsf_420[k];

        t_631[k] = f_10 * lsd0_253[k]
                   - f_11 * lsd1_253[k]
                   + f_3 * pc_x[k] * lsf_421[k];

        t_632[k] = f_10 * lsd0_254[k]
                   - f_11 * lsd1_254[k]
                   + f_3 * pc_x[k] * lsf_422[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pc_x, lsd0_255, lsd0_256, lsd0_257, \
                         lsd1_255, lsd1_256, lsd1_257, lsf_423, lsf_424, lsf_425, \
                         lsf_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_4 * lsd0_255[k]
                   - f_5 * lsd1_255[k]
                   + f_3 * pc_x[k] * lsf_423[k];

        t_634[k] = f_4 * lsd0_256[k]
                   - f_5 * lsd1_256[k]
                   + f_3 * pc_x[k] * lsf_424[k];

        t_635[k] = f_4 * lsd0_257[k]
                   - f_5 * lsd1_257[k]
                   + f_3 * pc_x[k] * lsf_425[k];

        t_636[k] = f_3 * pc_x[k] * lsf_426[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, t_641, pc_x, pc_y, pc_z, ksf_336, \
                         ksf_346, lsd0_255, lsd1_255, lsf_426, lsf_427, lsf_428, \
                         lsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_3 * pc_x[k] * lsf_427[k];

        t_638[k] = f_3 * pc_x[k] * lsf_428[k];

        t_639[k] = f_3 * pc_x[k] * lsf_429[k];

        t_640[k] = f_8 * ksf_346[k]
                   + f_1 * lsd0_255[k]
                   - f_2 * lsd1_255[k]
                   + f_3 * pc_y[k] * lsf_426[k];

        t_641[k] = f_12 * ksf_336[k]
                   + f_3 * pc_z[k] * lsf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pc_y, pc_z, ksg0_525, ksf_339, \
                         ksf_348, ksf_349, ksg1_525, lsd0_257, lsd1_257, lsf_428, \
                         lsf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * ksf_348[k]
                   + f_4 * lsd0_257[k]
                   - f_5 * lsd1_257[k]
                   + f_3 * pc_y[k] * lsf_428[k];

        t_643[k] = f_8 * ksf_349[k]
                   + f_3 * pc_y[k] * lsf_429[k];

        t_644[k] = f_12 * ksf_339[k]
                   + f_1 * lsd0_257[k]
                   - f_2 * lsd1_257[k]
                   + f_3 * pc_z[k] * lsf_429[k];

        t_645[k] = pa_y[k] * ksg0_525[k]
                   - f_6 * pc_y[k] * ksg1_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, pa_y, pc_x, pc_y, ksg0_527, ksg1_527, lsd0_259, \
                         lsd0_261, lsd1_259, lsd1_261, lsf_431, \
                         lsf_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_10 * lsd0_259[k]
                   - f_11 * lsd1_259[k]
                   + f_3 * pc_x[k] * lsf_431[k];

        t_647[k] = pa_y[k] * ksg0_527[k]
                   - f_6 * pc_y[k] * ksg1_527[k];

        t_648[k] = f_4 * lsd0_261[k]
                   - f_5 * lsd1_261[k]
                   + f_3 * pc_x[k] * lsf_433[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, pa_y, pc_x, pc_y, ksg0_530, \
                         ksg1_530, lsd0_262, lsd1_262, lsf_434, lsf_436, lsf_437, \
                         lsf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_4 * lsd0_262[k]
                   - f_5 * lsd1_262[k]
                   + f_3 * pc_x[k] * lsf_434[k];

        t_650[k] = pa_y[k] * ksg0_530[k]
                   - f_6 * pc_y[k] * ksg1_530[k];

        t_651[k] = f_3 * pc_x[k] * lsf_436[k];

        t_652[k] = f_3 * pc_x[k] * lsf_437[k];

        t_653[k] = f_3 * pc_x[k] * lsf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pa_y, pc_x, pc_y, pc_z, ksg0_535, ksf_346, \
                         ksf_356, ksg1_535, lsf_436, lsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_3 * pc_x[k] * lsf_439[k];

        t_655[k] = pa_y[k] * ksg0_535[k]
                   + f_15 * ksf_356[k]
                   - f_6 * pc_y[k] * ksg1_535[k];

        t_656[k] = f_9 * ksf_346[k]
                   + f_3 * pc_z[k] * lsf_436[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_y, pc_y, ksg0_537, ksg0_539, ksf_358, \
                         ksf_359, ksg1_537, ksg1_539, lsf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pa_y[k] * ksg0_537[k]
                   + f_8 * ksf_358[k]
                   - f_6 * pc_y[k] * ksg1_537[k];

        t_658[k] = f_7 * ksf_359[k]
                   + f_3 * pc_y[k] * lsf_439[k];

        t_659[k] = pa_y[k] * ksg0_539[k]
                   - f_6 * pc_y[k] * ksg1_539[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, lsd0_264, lsd0_266, \
                         lsd0_267, lsd1_264, lsd1_266, lsd1_267, lsf_440, lsf_442, \
                         lsf_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_1 * lsd0_264[k]
                   - f_2 * lsd1_264[k]
                   + f_3 * pc_x[k] * lsf_440[k];

        t_661[k] = f_3 * pc_y[k] * lsf_440[k];

        t_662[k] = f_10 * lsd0_266[k]
                   - f_11 * lsd1_266[k]
                   + f_3 * pc_x[k] * lsf_442[k];

        t_663[k] = f_4 * lsd0_267[k]
                   - f_5 * lsd1_267[k]
                   + f_3 * pc_x[k] * lsf_443[k];

        t_664[k] = f_3 * pc_y[k] * lsf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, pc_x, lsd0_269, lsd1_269, lsf_445, \
                         lsf_446, lsf_447, lsf_448, lsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_4 * lsd0_269[k]
                   - f_5 * lsd1_269[k]
                   + f_3 * pc_x[k] * lsf_445[k];

        t_666[k] = f_3 * pc_x[k] * lsf_446[k];

        t_667[k] = f_3 * pc_x[k] * lsf_447[k];

        t_668[k] = f_3 * pc_x[k] * lsf_448[k];

        t_669[k] = f_3 * pc_x[k] * lsf_449[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, lsd0_267, lsd0_268, lsd0_269, \
                         lsd1_267, lsd1_268, lsd1_269, lsf_446, lsf_447, lsf_448, \
                         lsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_1 * lsd0_267[k]
                   - f_2 * lsd1_267[k]
                   + f_3 * pc_y[k] * lsf_446[k];

        t_671[k] = f_10 * lsd0_268[k]
                   - f_11 * lsd1_268[k]
                   + f_3 * pc_y[k] * lsf_447[k];

        t_672[k] = f_4 * lsd0_269[k]
                   - f_5 * lsd1_269[k]
                   + f_3 * pc_y[k] * lsf_448[k];

        t_673[k] = f_3 * pc_y[k] * lsf_449[k];
    }

#pragma omp simd aligned(t_674, pc_z, ksf_359, lsd0_269, lsd1_269, \
                         lsf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * ksf_359[k]
                   + f_1 * lsd0_269[k]
                   - f_2 * lsd1_269[k]
                   + f_3 * pc_z[k] * lsf_449[k];
    }
}

auto
compute_prim_lsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksg0, const size_t ksf,
                                                   const size_t ksg1, const size_t lsd0,
                                                   const size_t lsd1, const size_t lsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);

    compute_prim_lsg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);

    compute_prim_lsg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);

    compute_prim_lsg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);

    compute_prim_lsg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);

    compute_prim_lsg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ksg0, ksf,
                                                              ksg1, lsd0, lsd1, lsf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
