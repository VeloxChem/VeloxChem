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


#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksd0,
                                                          const size_t ksp, const size_t ksd1,
                                                          const size_t lss0, const size_t lss1,
                                                          const size_t lsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksd0_0 = buffer.data(ksd0 + 0);
    const auto *ksd0_3 = buffer.data(ksd0 + 3);
    const auto *ksd0_5 = buffer.data(ksd0 + 5);
    const auto *ksd0_9 = buffer.data(ksd0 + 9);
    const auto *ksd0_12 = buffer.data(ksd0 + 12);
    const auto *ksd0_17 = buffer.data(ksd0 + 17);
    const auto *ksd0_18 = buffer.data(ksd0 + 18);
    const auto *ksd0_21 = buffer.data(ksd0 + 21);
    const auto *ksd0_30 = buffer.data(ksd0 + 30);
    const auto *ksd0_35 = buffer.data(ksd0 + 35);
    const auto *ksd0_36 = buffer.data(ksd0 + 36);
    const auto *ksd0_39 = buffer.data(ksd0 + 39);
    const auto *ksd0_54 = buffer.data(ksd0 + 54);
    const auto *ksd0_59 = buffer.data(ksd0 + 59);
    const auto *ksd0_60 = buffer.data(ksd0 + 60);
    const auto *ksd0_63 = buffer.data(ksd0 + 63);
    const auto *ksd0_84 = buffer.data(ksd0 + 84);
    const auto *ksd0_89 = buffer.data(ksd0 + 89);

    const auto *ksp_0 = buffer.data(ksp + 0);
    const auto *ksp_1 = buffer.data(ksp + 1);
    const auto *ksp_2 = buffer.data(ksp + 2);
    const auto *ksp_4 = buffer.data(ksp + 4);
    const auto *ksp_8 = buffer.data(ksp + 8);
    const auto *ksp_9 = buffer.data(ksp + 9);
    const auto *ksp_10 = buffer.data(ksp + 10);
    const auto *ksp_11 = buffer.data(ksp + 11);
    const auto *ksp_13 = buffer.data(ksp + 13);
    const auto *ksp_14 = buffer.data(ksp + 14);
    const auto *ksp_15 = buffer.data(ksp + 15);
    const auto *ksp_16 = buffer.data(ksp + 16);
    const auto *ksp_17 = buffer.data(ksp + 17);
    const auto *ksp_18 = buffer.data(ksp + 18);
    const auto *ksp_19 = buffer.data(ksp + 19);
    const auto *ksp_20 = buffer.data(ksp + 20);
    const auto *ksp_22 = buffer.data(ksp + 22);
    const auto *ksp_23 = buffer.data(ksp + 23);
    const auto *ksp_25 = buffer.data(ksp + 25);
    const auto *ksp_26 = buffer.data(ksp + 26);
    const auto *ksp_27 = buffer.data(ksp + 27);
    const auto *ksp_28 = buffer.data(ksp + 28);
    const auto *ksp_29 = buffer.data(ksp + 29);
    const auto *ksp_30 = buffer.data(ksp + 30);
    const auto *ksp_31 = buffer.data(ksp + 31);
    const auto *ksp_32 = buffer.data(ksp + 32);
    const auto *ksp_34 = buffer.data(ksp + 34);
    const auto *ksp_35 = buffer.data(ksp + 35);
    const auto *ksp_36 = buffer.data(ksp + 36);
    const auto *ksp_37 = buffer.data(ksp + 37);
    const auto *ksp_38 = buffer.data(ksp + 38);
    const auto *ksp_40 = buffer.data(ksp + 40);
    const auto *ksp_41 = buffer.data(ksp + 41);
    const auto *ksp_42 = buffer.data(ksp + 42);
    const auto *ksp_43 = buffer.data(ksp + 43);
    const auto *ksp_44 = buffer.data(ksp + 44);
    const auto *ksp_45 = buffer.data(ksp + 45);
    const auto *ksp_46 = buffer.data(ksp + 46);
    const auto *ksp_49 = buffer.data(ksp + 49);
    const auto *ksp_50 = buffer.data(ksp + 50);
    const auto *ksp_51 = buffer.data(ksp + 51);
    const auto *ksp_52 = buffer.data(ksp + 52);
    const auto *ksp_53 = buffer.data(ksp + 53);
    const auto *ksp_54 = buffer.data(ksp + 54);
    const auto *ksp_55 = buffer.data(ksp + 55);
    const auto *ksp_56 = buffer.data(ksp + 56);
    const auto *ksp_58 = buffer.data(ksp + 58);
    const auto *ksp_59 = buffer.data(ksp + 59);
    const auto *ksp_60 = buffer.data(ksp + 60);
    const auto *ksp_62 = buffer.data(ksp + 62);

    const auto *ksd1_0 = buffer.data(ksd1 + 0);
    const auto *ksd1_3 = buffer.data(ksd1 + 3);
    const auto *ksd1_5 = buffer.data(ksd1 + 5);
    const auto *ksd1_9 = buffer.data(ksd1 + 9);
    const auto *ksd1_12 = buffer.data(ksd1 + 12);
    const auto *ksd1_17 = buffer.data(ksd1 + 17);
    const auto *ksd1_18 = buffer.data(ksd1 + 18);
    const auto *ksd1_21 = buffer.data(ksd1 + 21);
    const auto *ksd1_30 = buffer.data(ksd1 + 30);
    const auto *ksd1_35 = buffer.data(ksd1 + 35);
    const auto *ksd1_36 = buffer.data(ksd1 + 36);
    const auto *ksd1_39 = buffer.data(ksd1 + 39);
    const auto *ksd1_54 = buffer.data(ksd1 + 54);
    const auto *ksd1_59 = buffer.data(ksd1 + 59);
    const auto *ksd1_60 = buffer.data(ksd1 + 60);
    const auto *ksd1_63 = buffer.data(ksd1 + 63);
    const auto *ksd1_84 = buffer.data(ksd1 + 84);
    const auto *ksd1_89 = buffer.data(ksd1 + 89);

    const auto *lss0_0 = buffer.data(lss0 + 0);
    const auto *lss0_1 = buffer.data(lss0 + 1);
    const auto *lss0_2 = buffer.data(lss0 + 2);
    const auto *lss0_3 = buffer.data(lss0 + 3);
    const auto *lss0_5 = buffer.data(lss0 + 5);
    const auto *lss0_6 = buffer.data(lss0 + 6);
    const auto *lss0_7 = buffer.data(lss0 + 7);
    const auto *lss0_8 = buffer.data(lss0 + 8);
    const auto *lss0_9 = buffer.data(lss0 + 9);
    const auto *lss0_10 = buffer.data(lss0 + 10);
    const auto *lss0_11 = buffer.data(lss0 + 11);
    const auto *lss0_12 = buffer.data(lss0 + 12);
    const auto *lss0_13 = buffer.data(lss0 + 13);
    const auto *lss0_14 = buffer.data(lss0 + 14);
    const auto *lss0_15 = buffer.data(lss0 + 15);
    const auto *lss0_16 = buffer.data(lss0 + 16);
    const auto *lss0_17 = buffer.data(lss0 + 17);
    const auto *lss0_18 = buffer.data(lss0 + 18);
    const auto *lss0_19 = buffer.data(lss0 + 19);
    const auto *lss0_20 = buffer.data(lss0 + 20);

    const auto *lss1_0 = buffer.data(lss1 + 0);
    const auto *lss1_1 = buffer.data(lss1 + 1);
    const auto *lss1_2 = buffer.data(lss1 + 2);
    const auto *lss1_3 = buffer.data(lss1 + 3);
    const auto *lss1_5 = buffer.data(lss1 + 5);
    const auto *lss1_6 = buffer.data(lss1 + 6);
    const auto *lss1_7 = buffer.data(lss1 + 7);
    const auto *lss1_8 = buffer.data(lss1 + 8);
    const auto *lss1_9 = buffer.data(lss1 + 9);
    const auto *lss1_10 = buffer.data(lss1 + 10);
    const auto *lss1_11 = buffer.data(lss1 + 11);
    const auto *lss1_12 = buffer.data(lss1 + 12);
    const auto *lss1_13 = buffer.data(lss1 + 13);
    const auto *lss1_14 = buffer.data(lss1 + 14);
    const auto *lss1_15 = buffer.data(lss1 + 15);
    const auto *lss1_16 = buffer.data(lss1 + 16);
    const auto *lss1_17 = buffer.data(lss1 + 17);
    const auto *lss1_18 = buffer.data(lss1 + 18);
    const auto *lss1_19 = buffer.data(lss1 + 19);
    const auto *lss1_20 = buffer.data(lss1 + 20);

    const auto *lsp_0 = buffer.data(lsp + 0);
    const auto *lsp_1 = buffer.data(lsp + 1);
    const auto *lsp_2 = buffer.data(lsp + 2);
    const auto *lsp_3 = buffer.data(lsp + 3);
    const auto *lsp_4 = buffer.data(lsp + 4);
    const auto *lsp_6 = buffer.data(lsp + 6);
    const auto *lsp_8 = buffer.data(lsp + 8);
    const auto *lsp_9 = buffer.data(lsp + 9);
    const auto *lsp_10 = buffer.data(lsp + 10);
    const auto *lsp_11 = buffer.data(lsp + 11);
    const auto *lsp_13 = buffer.data(lsp + 13);
    const auto *lsp_14 = buffer.data(lsp + 14);
    const auto *lsp_15 = buffer.data(lsp + 15);
    const auto *lsp_16 = buffer.data(lsp + 16);
    const auto *lsp_17 = buffer.data(lsp + 17);
    const auto *lsp_18 = buffer.data(lsp + 18);
    const auto *lsp_19 = buffer.data(lsp + 19);
    const auto *lsp_20 = buffer.data(lsp + 20);
    const auto *lsp_22 = buffer.data(lsp + 22);
    const auto *lsp_23 = buffer.data(lsp + 23);
    const auto *lsp_25 = buffer.data(lsp + 25);
    const auto *lsp_26 = buffer.data(lsp + 26);
    const auto *lsp_27 = buffer.data(lsp + 27);
    const auto *lsp_28 = buffer.data(lsp + 28);
    const auto *lsp_29 = buffer.data(lsp + 29);
    const auto *lsp_30 = buffer.data(lsp + 30);
    const auto *lsp_31 = buffer.data(lsp + 31);
    const auto *lsp_32 = buffer.data(lsp + 32);
    const auto *lsp_34 = buffer.data(lsp + 34);
    const auto *lsp_35 = buffer.data(lsp + 35);
    const auto *lsp_36 = buffer.data(lsp + 36);
    const auto *lsp_37 = buffer.data(lsp + 37);
    const auto *lsp_38 = buffer.data(lsp + 38);
    const auto *lsp_40 = buffer.data(lsp + 40);
    const auto *lsp_41 = buffer.data(lsp + 41);
    const auto *lsp_42 = buffer.data(lsp + 42);
    const auto *lsp_43 = buffer.data(lsp + 43);
    const auto *lsp_44 = buffer.data(lsp + 44);
    const auto *lsp_45 = buffer.data(lsp + 45);
    const auto *lsp_46 = buffer.data(lsp + 46);
    const auto *lsp_47 = buffer.data(lsp + 47);
    const auto *lsp_49 = buffer.data(lsp + 49);
    const auto *lsp_50 = buffer.data(lsp + 50);
    const auto *lsp_51 = buffer.data(lsp + 51);
    const auto *lsp_52 = buffer.data(lsp + 52);
    const auto *lsp_53 = buffer.data(lsp + 53);
    const auto *lsp_54 = buffer.data(lsp + 54);
    const auto *lsp_55 = buffer.data(lsp + 55);
    const auto *lsp_56 = buffer.data(lsp + 56);
    const auto *lsp_58 = buffer.data(lsp + 58);
    const auto *lsp_59 = buffer.data(lsp + 59);
    const auto *lsp_60 = buffer.data(lsp + 60);
    const auto *lsp_61 = buffer.data(lsp + 61);
    const auto *lsp_62 = buffer.data(lsp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksp_0, lss0_0, \
                         lss1_0, lsp_0, lsp_1, lsp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksp_0[k]
                 + f_1 * lss0_0[k]
                 - f_2 * lss1_0[k]
                 + f_3 * pc_x[k] * lsp_0[k];

        t_1[k] = f_3 * pc_y[k] * lsp_0[k];

        t_2[k] = f_3 * pc_z[k] * lsp_0[k];

        t_3[k] = f_1 * lss0_0[k]
                 - f_2 * lss1_0[k]
                 + f_3 * pc_y[k] * lsp_1[k];

        t_4[k] = f_3 * pc_y[k] * lsp_2[k];

        t_5[k] = f_1 * lss0_0[k]
                 - f_2 * lss1_0[k]
                 + f_3 * pc_z[k] * lsp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, ksd0_0, ksp_1, ksp_4, \
                         ksd1_0, lss0_1, lss1_1, lsp_3, lsp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * ksd0_0[k]
                 - f_4 * pc_y[k] * ksd1_0[k];

        t_7[k] = f_5 * ksp_4[k]
                 + f_3 * pc_x[k] * lsp_4[k];

        t_8[k] = f_3 * pc_z[k] * lsp_3[k];

        t_9[k] = f_6 * ksp_1[k]
                 + f_1 * lss0_1[k]
                 - f_2 * lss1_1[k]
                 + f_3 * pc_y[k] * lsp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, ksd0_0, ksd0_5, \
                         ksd1_0, ksd1_5, lsp_4, lsp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * lsp_4[k];

        t_11[k] = pa_y[k] * ksd0_5[k]
                  - f_4 * pc_y[k] * ksd1_5[k];

        t_12[k] = pa_z[k] * ksd0_0[k]
                  - f_4 * pc_z[k] * ksd1_0[k];

        t_13[k] = f_3 * pc_y[k] * lsp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, ksd0_3, ksp_2, ksp_8, \
                         ksd1_3, lss0_2, lss1_2, lsp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * ksp_8[k]
                  + f_3 * pc_x[k] * lsp_8[k];

        t_15[k] = pa_z[k] * ksd0_3[k]
                  - f_4 * pc_z[k] * ksd1_3[k];

        t_16[k] = f_3 * pc_y[k] * lsp_8[k];

        t_17[k] = f_6 * ksp_2[k]
                  + f_1 * lss0_2[k]
                  - f_2 * lss1_2[k]
                  + f_3 * pc_z[k] * lsp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, ksp_4, ksp_9, ksp_10, \
                         lss0_3, lss1_3, lsp_9, lsp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * ksp_9[k]
                  + f_1 * lss0_3[k]
                  - f_2 * lss1_3[k]
                  + f_3 * pc_x[k] * lsp_9[k];

        t_19[k] = f_7 * ksp_10[k]
                  + f_3 * pc_x[k] * lsp_10[k];

        t_20[k] = f_3 * pc_z[k] * lsp_9[k];

        t_21[k] = f_8 * ksp_4[k]
                  + f_1 * lss0_3[k]
                  - f_2 * lss1_3[k]
                  + f_3 * pc_y[k] * lsp_10[k];

        t_22[k] = f_3 * pc_z[k] * lsp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, ksd0_12, ksp_13, ksd1_12, \
                         lss0_3, lss1_3, lsp_11, lsp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * lss0_3[k]
                  - f_2 * lss1_3[k]
                  + f_3 * pc_z[k] * lsp_11[k];

        t_24[k] = pa_y[k] * ksd0_12[k]
                  - f_4 * pc_y[k] * ksd1_12[k];

        t_25[k] = f_7 * ksp_13[k]
                  + f_3 * pc_x[k] * lsp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, ksd0_9, \
                         ksd0_17, ksp_8, ksp_14, ksd1_9, ksd1_17, \
                         lsp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * ksp_14[k]
                  + f_3 * pc_x[k] * lsp_14[k];

        t_27[k] = pa_z[k] * ksd0_9[k]
                  - f_4 * pc_z[k] * ksd1_9[k];

        t_28[k] = f_6 * ksp_8[k]
                  + f_3 * pc_y[k] * lsp_14[k];

        t_29[k] = pa_y[k] * ksd0_17[k]
                  - f_4 * pc_y[k] * ksd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, ksp_15, ksp_17, lss0_5, \
                         lss1_5, lsp_15, lsp_16, lsp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * ksp_15[k]
                  + f_1 * lss0_5[k]
                  - f_2 * lss1_5[k]
                  + f_3 * pc_x[k] * lsp_15[k];

        t_31[k] = f_3 * pc_y[k] * lsp_15[k];

        t_32[k] = f_7 * ksp_17[k]
                  + f_3 * pc_x[k] * lsp_17[k];

        t_33[k] = f_1 * lss0_5[k]
                  - f_2 * lss1_5[k]
                  + f_3 * pc_y[k] * lsp_16[k];

        t_34[k] = f_3 * pc_y[k] * lsp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, ksp_8, ksp_18, ksp_19, lss0_5, \
                         lss0_6, lss1_5, lss1_6, lsp_17, lsp_18, \
                         lsp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * ksp_8[k]
                  + f_1 * lss0_5[k]
                  - f_2 * lss1_5[k]
                  + f_3 * pc_z[k] * lsp_17[k];

        t_36[k] = f_9 * ksp_18[k]
                  + f_1 * lss0_6[k]
                  - f_2 * lss1_6[k]
                  + f_3 * pc_x[k] * lsp_18[k];

        t_37[k] = f_9 * ksp_19[k]
                  + f_3 * pc_x[k] * lsp_19[k];

        t_38[k] = f_3 * pc_z[k] * lsp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, ksd0_18, ksp_10, ksd1_18, \
                         lss0_6, lss1_6, lsp_19, lsp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * ksp_10[k]
                  + f_1 * lss0_6[k]
                  - f_2 * lss1_6[k]
                  + f_3 * pc_y[k] * lsp_19[k];

        t_40[k] = f_3 * pc_z[k] * lsp_19[k];

        t_41[k] = f_1 * lss0_6[k]
                  - f_2 * lss1_6[k]
                  + f_3 * pc_z[k] * lsp_20[k];

        t_42[k] = pa_z[k] * ksd0_18[k]
                  - f_4 * pc_z[k] * ksd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, ksd0_21, ksp_14, \
                         ksp_22, ksp_23, ksd1_21, lsp_22, lsp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * ksp_22[k]
                  + f_3 * pc_x[k] * lsp_22[k];

        t_44[k] = f_9 * ksp_23[k]
                  + f_3 * pc_x[k] * lsp_23[k];

        t_45[k] = pa_z[k] * ksd0_21[k]
                  - f_4 * pc_z[k] * ksd1_21[k];

        t_46[k] = f_8 * ksp_14[k]
                  + f_3 * pc_y[k] * lsp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, ksd0_30, ksp_11, ksp_25, \
                         ksd1_30, lss0_7, lss1_7, lsp_23, lsp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * ksp_11[k]
                  + f_1 * lss0_7[k]
                  - f_2 * lss1_7[k]
                  + f_3 * pc_z[k] * lsp_23[k];

        t_48[k] = pa_y[k] * ksd0_30[k]
                  - f_4 * pc_y[k] * ksd1_30[k];

        t_49[k] = f_9 * ksp_25[k]
                  + f_3 * pc_x[k] * lsp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, ksd0_35, ksp_16, ksp_17, \
                         ksp_26, ksd1_35, lss0_8, lss1_8, lsp_25, \
                         lsp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * ksp_26[k]
                  + f_3 * pc_x[k] * lsp_26[k];

        t_51[k] = f_6 * ksp_16[k]
                  + f_1 * lss0_8[k]
                  - f_2 * lss1_8[k]
                  + f_3 * pc_y[k] * lsp_25[k];

        t_52[k] = f_6 * ksp_17[k]
                  + f_3 * pc_y[k] * lsp_26[k];

        t_53[k] = pa_y[k] * ksd0_35[k]
                  - f_4 * pc_y[k] * ksd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, ksp_27, ksp_29, lss0_9, \
                         lss1_9, lsp_27, lsp_28, lsp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * ksp_27[k]
                  + f_1 * lss0_9[k]
                  - f_2 * lss1_9[k]
                  + f_3 * pc_x[k] * lsp_27[k];

        t_55[k] = f_3 * pc_y[k] * lsp_27[k];

        t_56[k] = f_9 * ksp_29[k]
                  + f_3 * pc_x[k] * lsp_29[k];

        t_57[k] = f_1 * lss0_9[k]
                  - f_2 * lss1_9[k]
                  + f_3 * pc_y[k] * lsp_28[k];

        t_58[k] = f_3 * pc_y[k] * lsp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, ksp_17, ksp_30, ksp_31, lss0_9, \
                         lss0_10, lss1_9, lss1_10, lsp_29, lsp_30, \
                         lsp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * ksp_17[k]
                  + f_1 * lss0_9[k]
                  - f_2 * lss1_9[k]
                  + f_3 * pc_z[k] * lsp_29[k];

        t_60[k] = f_11 * ksp_30[k]
                  + f_1 * lss0_10[k]
                  - f_2 * lss1_10[k]
                  + f_3 * pc_x[k] * lsp_30[k];

        t_61[k] = f_11 * ksp_31[k]
                  + f_3 * pc_x[k] * lsp_31[k];

        t_62[k] = f_3 * pc_z[k] * lsp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, ksd0_36, ksp_19, ksd1_36, \
                         lss0_10, lss1_10, lsp_31, lsp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_11 * ksp_19[k]
                  + f_1 * lss0_10[k]
                  - f_2 * lss1_10[k]
                  + f_3 * pc_y[k] * lsp_31[k];

        t_64[k] = f_3 * pc_z[k] * lsp_31[k];

        t_65[k] = f_1 * lss0_10[k]
                  - f_2 * lss1_10[k]
                  + f_3 * pc_z[k] * lsp_32[k];

        t_66[k] = pa_z[k] * ksd0_36[k]
                  - f_4 * pc_z[k] * ksd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, ksd0_39, ksp_23, \
                         ksp_34, ksp_35, ksd1_39, lsp_34, lsp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * ksp_34[k]
                  + f_3 * pc_x[k] * lsp_34[k];

        t_68[k] = f_11 * ksp_35[k]
                  + f_3 * pc_x[k] * lsp_35[k];

        t_69[k] = pa_z[k] * ksd0_39[k]
                  - f_4 * pc_z[k] * ksd1_39[k];

        t_70[k] = f_10 * ksp_23[k]
                  + f_3 * pc_y[k] * lsp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, ksp_20, ksp_36, ksp_37, lss0_11, \
                         lss0_12, lss1_11, lss1_12, lsp_35, lsp_36, \
                         lsp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * ksp_20[k]
                  + f_1 * lss0_11[k]
                  - f_2 * lss1_11[k]
                  + f_3 * pc_z[k] * lsp_35[k];

        t_72[k] = f_11 * ksp_36[k]
                  + f_1 * lss0_12[k]
                  - f_2 * lss1_12[k]
                  + f_3 * pc_x[k] * lsp_36[k];

        t_73[k] = f_11 * ksp_37[k]
                  + f_3 * pc_x[k] * lsp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, ksp_23, ksp_25, ksp_26, \
                         ksp_38, lss0_12, lss1_12, lsp_37, lsp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * ksp_38[k]
                  + f_3 * pc_x[k] * lsp_38[k];

        t_75[k] = f_8 * ksp_25[k]
                  + f_1 * lss0_12[k]
                  - f_2 * lss1_12[k]
                  + f_3 * pc_y[k] * lsp_37[k];

        t_76[k] = f_8 * ksp_26[k]
                  + f_3 * pc_y[k] * lsp_38[k];

        t_77[k] = f_8 * ksp_23[k]
                  + f_1 * lss0_12[k]
                  - f_2 * lss1_12[k]
                  + f_3 * pc_z[k] * lsp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, ksd0_54, ksp_28, ksp_40, \
                         ksp_41, ksd1_54, lss0_13, lss1_13, lsp_40, \
                         lsp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * ksd0_54[k]
                  - f_4 * pc_y[k] * ksd1_54[k];

        t_79[k] = f_11 * ksp_40[k]
                  + f_3 * pc_x[k] * lsp_40[k];

        t_80[k] = f_11 * ksp_41[k]
                  + f_3 * pc_x[k] * lsp_41[k];

        t_81[k] = f_6 * ksp_28[k]
                  + f_1 * lss0_13[k]
                  - f_2 * lss1_13[k]
                  + f_3 * pc_y[k] * lsp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, ksd0_59, ksp_29, ksp_42, \
                         ksd1_59, lss0_14, lss1_14, lsp_41, lsp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * ksp_29[k]
                  + f_3 * pc_y[k] * lsp_41[k];

        t_83[k] = pa_y[k] * ksd0_59[k]
                  - f_4 * pc_y[k] * ksd1_59[k];

        t_84[k] = f_11 * ksp_42[k]
                  + f_1 * lss0_14[k]
                  - f_2 * lss1_14[k]
                  + f_3 * pc_x[k] * lsp_42[k];

        t_85[k] = f_3 * pc_y[k] * lsp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, ksp_29, ksp_44, lss0_14, \
                         lss1_14, lsp_43, lsp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * ksp_44[k]
                  + f_3 * pc_x[k] * lsp_44[k];

        t_87[k] = f_1 * lss0_14[k]
                  - f_2 * lss1_14[k]
                  + f_3 * pc_y[k] * lsp_43[k];

        t_88[k] = f_3 * pc_y[k] * lsp_44[k];

        t_89[k] = f_11 * ksp_29[k]
                  + f_1 * lss0_14[k]
                  - f_2 * lss1_14[k]
                  + f_3 * pc_z[k] * lsp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, ksp_31, ksp_45, \
                         ksp_46, lss0_15, lss1_15, lsp_45, lsp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_10 * ksp_45[k]
                  + f_1 * lss0_15[k]
                  - f_2 * lss1_15[k]
                  + f_3 * pc_x[k] * lsp_45[k];

        t_91[k] = f_10 * ksp_46[k]
                  + f_3 * pc_x[k] * lsp_46[k];

        t_92[k] = f_3 * pc_z[k] * lsp_45[k];

        t_93[k] = f_9 * ksp_31[k]
                  + f_1 * lss0_15[k]
                  - f_2 * lss1_15[k]
                  + f_3 * pc_y[k] * lsp_46[k];

        t_94[k] = f_3 * pc_z[k] * lsp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, ksd0_60, ksp_49, ksp_50, \
                         ksd1_60, lss0_15, lss1_15, lsp_47, lsp_49, \
                         lsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * lss0_15[k]
                  - f_2 * lss1_15[k]
                  + f_3 * pc_z[k] * lsp_47[k];

        t_96[k] = pa_z[k] * ksd0_60[k]
                  - f_4 * pc_z[k] * ksd1_60[k];

        t_97[k] = f_10 * ksp_49[k]
                  + f_3 * pc_x[k] * lsp_49[k];

        t_98[k] = f_10 * ksp_50[k]
                  + f_3 * pc_x[k] * lsp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, ksd0_63, ksp_32, ksp_35, \
                         ksd1_63, lss0_16, lss1_16, lsp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * ksd0_63[k]
                  - f_4 * pc_z[k] * ksd1_63[k];

        t_100[k] = f_11 * ksp_35[k]
                   + f_3 * pc_y[k] * lsp_50[k];

        t_101[k] = f_6 * ksp_32[k]
                   + f_1 * lss0_16[k]
                   - f_2 * lss1_16[k]
                   + f_3 * pc_z[k] * lsp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, ksp_37, ksp_51, ksp_52, \
                         ksp_53, lss0_17, lss1_17, lsp_51, lsp_52, \
                         lsp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_10 * ksp_51[k]
                   + f_1 * lss0_17[k]
                   - f_2 * lss1_17[k]
                   + f_3 * pc_x[k] * lsp_51[k];

        t_103[k] = f_10 * ksp_52[k]
                   + f_3 * pc_x[k] * lsp_52[k];

        t_104[k] = f_10 * ksp_53[k]
                   + f_3 * pc_x[k] * lsp_53[k];

        t_105[k] = f_10 * ksp_37[k]
                   + f_1 * lss0_17[k]
                   - f_2 * lss1_17[k]
                   + f_3 * pc_y[k] * lsp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, ksp_35, ksp_38, ksp_54, \
                         lss0_17, lss0_18, lss1_17, lss1_18, lsp_53, \
                         lsp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * ksp_38[k]
                   + f_3 * pc_y[k] * lsp_53[k];

        t_107[k] = f_8 * ksp_35[k]
                   + f_1 * lss0_17[k]
                   - f_2 * lss1_17[k]
                   + f_3 * pc_z[k] * lsp_53[k];

        t_108[k] = f_10 * ksp_54[k]
                   + f_1 * lss0_18[k]
                   - f_2 * lss1_18[k]
                   + f_3 * pc_x[k] * lsp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, ksp_40, ksp_41, ksp_55, \
                         ksp_56, lss0_18, lss1_18, lsp_55, lsp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_10 * ksp_55[k]
                   + f_3 * pc_x[k] * lsp_55[k];

        t_110[k] = f_10 * ksp_56[k]
                   + f_3 * pc_x[k] * lsp_56[k];

        t_111[k] = f_8 * ksp_40[k]
                   + f_1 * lss0_18[k]
                   - f_2 * lss1_18[k]
                   + f_3 * pc_y[k] * lsp_55[k];

        t_112[k] = f_8 * ksp_41[k]
                   + f_3 * pc_y[k] * lsp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, ksd0_84, ksp_38, ksp_58, \
                         ksd1_84, lss0_18, lss1_18, lsp_56, lsp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * ksp_38[k]
                   + f_1 * lss0_18[k]
                   - f_2 * lss1_18[k]
                   + f_3 * pc_z[k] * lsp_56[k];

        t_114[k] = pa_y[k] * ksd0_84[k]
                   - f_4 * pc_y[k] * ksd1_84[k];

        t_115[k] = f_10 * ksp_58[k]
                   + f_3 * pc_x[k] * lsp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, ksd0_89, ksp_43, \
                         ksp_44, ksp_59, ksd1_89, lss0_19, lss1_19, lsp_58, \
                         lsp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_10 * ksp_59[k]
                   + f_3 * pc_x[k] * lsp_59[k];

        t_117[k] = f_6 * ksp_43[k]
                   + f_1 * lss0_19[k]
                   - f_2 * lss1_19[k]
                   + f_3 * pc_y[k] * lsp_58[k];

        t_118[k] = f_6 * ksp_44[k]
                   + f_3 * pc_y[k] * lsp_59[k];

        t_119[k] = pa_y[k] * ksd0_89[k]
                   - f_4 * pc_y[k] * ksd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, ksp_60, ksp_62, \
                         lss0_20, lss1_20, lsp_60, lsp_61, lsp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_10 * ksp_60[k]
                   + f_1 * lss0_20[k]
                   - f_2 * lss1_20[k]
                   + f_3 * pc_x[k] * lsp_60[k];

        t_121[k] = f_3 * pc_y[k] * lsp_60[k];

        t_122[k] = f_10 * ksp_62[k]
                   + f_3 * pc_x[k] * lsp_62[k];

        t_123[k] = f_1 * lss0_20[k]
                   - f_2 * lss1_20[k]
                   + f_3 * pc_y[k] * lsp_61[k];

        t_124[k] = f_3 * pc_y[k] * lsp_62[k];
    }
}

static auto
compute_prim_lsd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksd0,
                                                          const size_t ksp, const size_t ksd1,
                                                          const size_t lss0, const size_t lss1,
                                                          const size_t lsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 3.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksd0_90 = buffer.data(ksd0 + 90);
    const auto *ksd0_93 = buffer.data(ksd0 + 93);
    const auto *ksd0_120 = buffer.data(ksd0 + 120);
    const auto *ksd0_125 = buffer.data(ksd0 + 125);
    const auto *ksd0_126 = buffer.data(ksd0 + 126);
    const auto *ksd0_162 = buffer.data(ksd0 + 162);
    const auto *ksd0_168 = buffer.data(ksd0 + 168);
    const auto *ksd0_171 = buffer.data(ksd0 + 171);
    const auto *ksd0_173 = buffer.data(ksd0 + 173);
    const auto *ksd0_177 = buffer.data(ksd0 + 177);
    const auto *ksd0_179 = buffer.data(ksd0 + 179);
    const auto *ksd0_180 = buffer.data(ksd0 + 180);
    const auto *ksd0_183 = buffer.data(ksd0 + 183);
    const auto *ksd0_185 = buffer.data(ksd0 + 185);
    const auto *ksd0_186 = buffer.data(ksd0 + 186);
    const auto *ksd0_189 = buffer.data(ksd0 + 189);
    const auto *ksd0_191 = buffer.data(ksd0 + 191);
    const auto *ksd0_192 = buffer.data(ksd0 + 192);
    const auto *ksd0_195 = buffer.data(ksd0 + 195);
    const auto *ksd0_197 = buffer.data(ksd0 + 197);
    const auto *ksd0_198 = buffer.data(ksd0 + 198);
    const auto *ksd0_201 = buffer.data(ksd0 + 201);
    const auto *ksd0_203 = buffer.data(ksd0 + 203);
    const auto *ksd0_207 = buffer.data(ksd0 + 207);
    const auto *ksd0_209 = buffer.data(ksd0 + 209);
    const auto *ksd0_210 = buffer.data(ksd0 + 210);
    const auto *ksd0_213 = buffer.data(ksd0 + 213);
    const auto *ksd0_215 = buffer.data(ksd0 + 215);

    const auto *ksp_44 = buffer.data(ksp + 44);
    const auto *ksp_46 = buffer.data(ksp + 46);
    const auto *ksp_47 = buffer.data(ksp + 47);
    const auto *ksp_50 = buffer.data(ksp + 50);
    const auto *ksp_52 = buffer.data(ksp + 52);
    const auto *ksp_53 = buffer.data(ksp + 53);
    const auto *ksp_55 = buffer.data(ksp + 55);
    const auto *ksp_56 = buffer.data(ksp + 56);
    const auto *ksp_58 = buffer.data(ksp + 58);
    const auto *ksp_59 = buffer.data(ksp + 59);
    const auto *ksp_61 = buffer.data(ksp + 61);
    const auto *ksp_62 = buffer.data(ksp + 62);
    const auto *ksp_63 = buffer.data(ksp + 63);
    const auto *ksp_64 = buffer.data(ksp + 64);
    const auto *ksp_67 = buffer.data(ksp + 67);
    const auto *ksp_68 = buffer.data(ksp + 68);
    const auto *ksp_69 = buffer.data(ksp + 69);
    const auto *ksp_70 = buffer.data(ksp + 70);
    const auto *ksp_71 = buffer.data(ksp + 71);
    const auto *ksp_72 = buffer.data(ksp + 72);
    const auto *ksp_73 = buffer.data(ksp + 73);
    const auto *ksp_74 = buffer.data(ksp + 74);
    const auto *ksp_75 = buffer.data(ksp + 75);
    const auto *ksp_76 = buffer.data(ksp + 76);
    const auto *ksp_77 = buffer.data(ksp + 77);
    const auto *ksp_79 = buffer.data(ksp + 79);
    const auto *ksp_80 = buffer.data(ksp + 80);
    const auto *ksp_81 = buffer.data(ksp + 81);
    const auto *ksp_83 = buffer.data(ksp + 83);
    const auto *ksp_84 = buffer.data(ksp + 84);
    const auto *ksp_85 = buffer.data(ksp + 85);
    const auto *ksp_86 = buffer.data(ksp + 86);
    const auto *ksp_88 = buffer.data(ksp + 88);
    const auto *ksp_89 = buffer.data(ksp + 89);
    const auto *ksp_90 = buffer.data(ksp + 90);
    const auto *ksp_91 = buffer.data(ksp + 91);
    const auto *ksp_92 = buffer.data(ksp + 92);
    const auto *ksp_93 = buffer.data(ksp + 93);
    const auto *ksp_94 = buffer.data(ksp + 94);
    const auto *ksp_95 = buffer.data(ksp + 95);
    const auto *ksp_96 = buffer.data(ksp + 96);
    const auto *ksp_97 = buffer.data(ksp + 97);
    const auto *ksp_98 = buffer.data(ksp + 98);
    const auto *ksp_99 = buffer.data(ksp + 99);
    const auto *ksp_100 = buffer.data(ksp + 100);
    const auto *ksp_101 = buffer.data(ksp + 101);
    const auto *ksp_103 = buffer.data(ksp + 103);
    const auto *ksp_104 = buffer.data(ksp + 104);
    const auto *ksp_105 = buffer.data(ksp + 105);
    const auto *ksp_107 = buffer.data(ksp + 107);

    const auto *ksd1_90 = buffer.data(ksd1 + 90);
    const auto *ksd1_93 = buffer.data(ksd1 + 93);
    const auto *ksd1_120 = buffer.data(ksd1 + 120);
    const auto *ksd1_125 = buffer.data(ksd1 + 125);
    const auto *ksd1_126 = buffer.data(ksd1 + 126);
    const auto *ksd1_162 = buffer.data(ksd1 + 162);
    const auto *ksd1_168 = buffer.data(ksd1 + 168);
    const auto *ksd1_171 = buffer.data(ksd1 + 171);
    const auto *ksd1_173 = buffer.data(ksd1 + 173);
    const auto *ksd1_177 = buffer.data(ksd1 + 177);
    const auto *ksd1_179 = buffer.data(ksd1 + 179);
    const auto *ksd1_180 = buffer.data(ksd1 + 180);
    const auto *ksd1_183 = buffer.data(ksd1 + 183);
    const auto *ksd1_185 = buffer.data(ksd1 + 185);
    const auto *ksd1_186 = buffer.data(ksd1 + 186);
    const auto *ksd1_189 = buffer.data(ksd1 + 189);
    const auto *ksd1_191 = buffer.data(ksd1 + 191);
    const auto *ksd1_192 = buffer.data(ksd1 + 192);
    const auto *ksd1_195 = buffer.data(ksd1 + 195);
    const auto *ksd1_197 = buffer.data(ksd1 + 197);
    const auto *ksd1_198 = buffer.data(ksd1 + 198);
    const auto *ksd1_201 = buffer.data(ksd1 + 201);
    const auto *ksd1_203 = buffer.data(ksd1 + 203);
    const auto *ksd1_207 = buffer.data(ksd1 + 207);
    const auto *ksd1_209 = buffer.data(ksd1 + 209);
    const auto *ksd1_210 = buffer.data(ksd1 + 210);
    const auto *ksd1_213 = buffer.data(ksd1 + 213);
    const auto *ksd1_215 = buffer.data(ksd1 + 215);

    const auto *lss0_20 = buffer.data(lss0 + 20);
    const auto *lss0_21 = buffer.data(lss0 + 21);
    const auto *lss0_22 = buffer.data(lss0 + 22);
    const auto *lss0_23 = buffer.data(lss0 + 23);
    const auto *lss0_24 = buffer.data(lss0 + 24);
    const auto *lss0_25 = buffer.data(lss0 + 25);
    const auto *lss0_26 = buffer.data(lss0 + 26);
    const auto *lss0_27 = buffer.data(lss0 + 27);
    const auto *lss0_36 = buffer.data(lss0 + 36);
    const auto *lss0_37 = buffer.data(lss0 + 37);
    const auto *lss0_38 = buffer.data(lss0 + 38);
    const auto *lss0_39 = buffer.data(lss0 + 39);
    const auto *lss0_40 = buffer.data(lss0 + 40);
    const auto *lss0_41 = buffer.data(lss0 + 41);

    const auto *lss1_20 = buffer.data(lss1 + 20);
    const auto *lss1_21 = buffer.data(lss1 + 21);
    const auto *lss1_22 = buffer.data(lss1 + 22);
    const auto *lss1_23 = buffer.data(lss1 + 23);
    const auto *lss1_24 = buffer.data(lss1 + 24);
    const auto *lss1_25 = buffer.data(lss1 + 25);
    const auto *lss1_26 = buffer.data(lss1 + 26);
    const auto *lss1_27 = buffer.data(lss1 + 27);
    const auto *lss1_36 = buffer.data(lss1 + 36);
    const auto *lss1_37 = buffer.data(lss1 + 37);
    const auto *lss1_38 = buffer.data(lss1 + 38);
    const auto *lss1_39 = buffer.data(lss1 + 39);
    const auto *lss1_40 = buffer.data(lss1 + 40);
    const auto *lss1_41 = buffer.data(lss1 + 41);

    const auto *lsp_62 = buffer.data(lsp + 62);
    const auto *lsp_63 = buffer.data(lsp + 63);
    const auto *lsp_64 = buffer.data(lsp + 64);
    const auto *lsp_65 = buffer.data(lsp + 65);
    const auto *lsp_67 = buffer.data(lsp + 67);
    const auto *lsp_68 = buffer.data(lsp + 68);
    const auto *lsp_69 = buffer.data(lsp + 69);
    const auto *lsp_70 = buffer.data(lsp + 70);
    const auto *lsp_71 = buffer.data(lsp + 71);
    const auto *lsp_72 = buffer.data(lsp + 72);
    const auto *lsp_73 = buffer.data(lsp + 73);
    const auto *lsp_74 = buffer.data(lsp + 74);
    const auto *lsp_75 = buffer.data(lsp + 75);
    const auto *lsp_76 = buffer.data(lsp + 76);
    const auto *lsp_77 = buffer.data(lsp + 77);
    const auto *lsp_79 = buffer.data(lsp + 79);
    const auto *lsp_80 = buffer.data(lsp + 80);
    const auto *lsp_81 = buffer.data(lsp + 81);
    const auto *lsp_82 = buffer.data(lsp + 82);
    const auto *lsp_83 = buffer.data(lsp + 83);
    const auto *lsp_84 = buffer.data(lsp + 84);
    const auto *lsp_85 = buffer.data(lsp + 85);
    const auto *lsp_88 = buffer.data(lsp + 88);
    const auto *lsp_89 = buffer.data(lsp + 89);
    const auto *lsp_91 = buffer.data(lsp + 91);
    const auto *lsp_92 = buffer.data(lsp + 92);
    const auto *lsp_94 = buffer.data(lsp + 94);
    const auto *lsp_95 = buffer.data(lsp + 95);
    const auto *lsp_97 = buffer.data(lsp + 97);
    const auto *lsp_98 = buffer.data(lsp + 98);
    const auto *lsp_100 = buffer.data(lsp + 100);
    const auto *lsp_101 = buffer.data(lsp + 101);
    const auto *lsp_103 = buffer.data(lsp + 103);
    const auto *lsp_104 = buffer.data(lsp + 104);
    const auto *lsp_105 = buffer.data(lsp + 105);
    const auto *lsp_107 = buffer.data(lsp + 107);
    const auto *lsp_108 = buffer.data(lsp + 108);
    const auto *lsp_109 = buffer.data(lsp + 109);
    const auto *lsp_110 = buffer.data(lsp + 110);
    const auto *lsp_112 = buffer.data(lsp + 112);
    const auto *lsp_113 = buffer.data(lsp + 113);
    const auto *lsp_114 = buffer.data(lsp + 114);
    const auto *lsp_115 = buffer.data(lsp + 115);
    const auto *lsp_116 = buffer.data(lsp + 116);
    const auto *lsp_117 = buffer.data(lsp + 117);
    const auto *lsp_118 = buffer.data(lsp + 118);
    const auto *lsp_119 = buffer.data(lsp + 119);
    const auto *lsp_120 = buffer.data(lsp + 120);
    const auto *lsp_121 = buffer.data(lsp + 121);
    const auto *lsp_122 = buffer.data(lsp + 122);
    const auto *lsp_123 = buffer.data(lsp + 123);
    const auto *lsp_124 = buffer.data(lsp + 124);
    const auto *lsp_125 = buffer.data(lsp + 125);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_z, ksp_44, ksp_63, ksp_64, \
                         lss0_20, lss0_21, lss1_20, lss1_21, lsp_62, lsp_63, \
                         lsp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_9 * ksp_44[k]
                   + f_1 * lss0_20[k]
                   - f_2 * lss1_20[k]
                   + f_3 * pc_z[k] * lsp_62[k];

        t_126[k] = f_8 * ksp_63[k]
                   + f_1 * lss0_21[k]
                   - f_2 * lss1_21[k]
                   + f_3 * pc_x[k] * lsp_63[k];

        t_127[k] = f_8 * ksp_64[k]
                   + f_3 * pc_x[k] * lsp_64[k];

        t_128[k] = f_3 * pc_z[k] * lsp_63[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pc_y, pc_z, ksd0_90, ksp_46, \
                         ksd1_90, lss0_21, lss1_21, lsp_64, lsp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_7 * ksp_46[k]
                   + f_1 * lss0_21[k]
                   - f_2 * lss1_21[k]
                   + f_3 * pc_y[k] * lsp_64[k];

        t_130[k] = f_3 * pc_z[k] * lsp_64[k];

        t_131[k] = f_1 * lss0_21[k]
                   - f_2 * lss1_21[k]
                   + f_3 * pc_z[k] * lsp_65[k];

        t_132[k] = pa_z[k] * ksd0_90[k]
                   - f_4 * pc_z[k] * ksd1_90[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, ksd0_93, ksp_50, \
                         ksp_67, ksp_68, ksd1_93, lsp_67, lsp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_8 * ksp_67[k]
                   + f_3 * pc_x[k] * lsp_67[k];

        t_134[k] = f_8 * ksp_68[k]
                   + f_3 * pc_x[k] * lsp_68[k];

        t_135[k] = pa_z[k] * ksd0_93[k]
                   - f_4 * pc_z[k] * ksd1_93[k];

        t_136[k] = f_9 * ksp_50[k]
                   + f_3 * pc_y[k] * lsp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_x, pc_z, ksp_47, ksp_69, ksp_70, lss0_22, \
                         lss0_23, lss1_22, lss1_23, lsp_68, lsp_69, \
                         lsp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * ksp_47[k]
                   + f_1 * lss0_22[k]
                   - f_2 * lss1_22[k]
                   + f_3 * pc_z[k] * lsp_68[k];

        t_138[k] = f_8 * ksp_69[k]
                   + f_1 * lss0_23[k]
                   - f_2 * lss1_23[k]
                   + f_3 * pc_x[k] * lsp_69[k];

        t_139[k] = f_8 * ksp_70[k]
                   + f_3 * pc_x[k] * lsp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, ksp_50, ksp_52, ksp_53, \
                         ksp_71, lss0_23, lss1_23, lsp_70, lsp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_8 * ksp_71[k]
                   + f_3 * pc_x[k] * lsp_71[k];

        t_141[k] = f_11 * ksp_52[k]
                   + f_1 * lss0_23[k]
                   - f_2 * lss1_23[k]
                   + f_3 * pc_y[k] * lsp_70[k];

        t_142[k] = f_11 * ksp_53[k]
                   + f_3 * pc_y[k] * lsp_71[k];

        t_143[k] = f_8 * ksp_50[k]
                   + f_1 * lss0_23[k]
                   - f_2 * lss1_23[k]
                   + f_3 * pc_z[k] * lsp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pc_x, pc_y, ksp_55, ksp_72, ksp_73, \
                         ksp_74, lss0_24, lss1_24, lsp_72, lsp_73, \
                         lsp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_8 * ksp_72[k]
                   + f_1 * lss0_24[k]
                   - f_2 * lss1_24[k]
                   + f_3 * pc_x[k] * lsp_72[k];

        t_145[k] = f_8 * ksp_73[k]
                   + f_3 * pc_x[k] * lsp_73[k];

        t_146[k] = f_8 * ksp_74[k]
                   + f_3 * pc_x[k] * lsp_74[k];

        t_147[k] = f_10 * ksp_55[k]
                   + f_1 * lss0_24[k]
                   - f_2 * lss1_24[k]
                   + f_3 * pc_y[k] * lsp_73[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, ksp_53, ksp_56, ksp_75, \
                         lss0_24, lss0_25, lss1_24, lss1_25, lsp_74, \
                         lsp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * ksp_56[k]
                   + f_3 * pc_y[k] * lsp_74[k];

        t_149[k] = f_10 * ksp_53[k]
                   + f_1 * lss0_24[k]
                   - f_2 * lss1_24[k]
                   + f_3 * pc_z[k] * lsp_74[k];

        t_150[k] = f_8 * ksp_75[k]
                   + f_1 * lss0_25[k]
                   - f_2 * lss1_25[k]
                   + f_3 * pc_x[k] * lsp_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, ksp_58, ksp_59, ksp_76, \
                         ksp_77, lss0_25, lss1_25, lsp_76, lsp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_8 * ksp_76[k]
                   + f_3 * pc_x[k] * lsp_76[k];

        t_152[k] = f_8 * ksp_77[k]
                   + f_3 * pc_x[k] * lsp_77[k];

        t_153[k] = f_8 * ksp_58[k]
                   + f_1 * lss0_25[k]
                   - f_2 * lss1_25[k]
                   + f_3 * pc_y[k] * lsp_76[k];

        t_154[k] = f_8 * ksp_59[k]
                   + f_3 * pc_y[k] * lsp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, ksd0_120, ksp_56, \
                         ksp_79, ksd1_120, lss0_25, lss1_25, lsp_77, \
                         lsp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_11 * ksp_56[k]
                   + f_1 * lss0_25[k]
                   - f_2 * lss1_25[k]
                   + f_3 * pc_z[k] * lsp_77[k];

        t_156[k] = pa_y[k] * ksd0_120[k]
                   - f_4 * pc_y[k] * ksd1_120[k];

        t_157[k] = f_8 * ksp_79[k]
                   + f_3 * pc_x[k] * lsp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, ksd0_125, ksp_61, \
                         ksp_62, ksp_80, ksd1_125, lss0_26, lss1_26, lsp_79, \
                         lsp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_8 * ksp_80[k]
                   + f_3 * pc_x[k] * lsp_80[k];

        t_159[k] = f_6 * ksp_61[k]
                   + f_1 * lss0_26[k]
                   - f_2 * lss1_26[k]
                   + f_3 * pc_y[k] * lsp_79[k];

        t_160[k] = f_6 * ksp_62[k]
                   + f_3 * pc_y[k] * lsp_80[k];

        t_161[k] = pa_y[k] * ksd0_125[k]
                   - f_4 * pc_y[k] * ksd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, ksp_81, ksp_83, \
                         lss0_27, lss1_27, lsp_81, lsp_82, lsp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_8 * ksp_81[k]
                   + f_1 * lss0_27[k]
                   - f_2 * lss1_27[k]
                   + f_3 * pc_x[k] * lsp_81[k];

        t_163[k] = f_3 * pc_y[k] * lsp_81[k];

        t_164[k] = f_8 * ksp_83[k]
                   + f_3 * pc_x[k] * lsp_83[k];

        t_165[k] = f_1 * lss0_27[k]
                   - f_2 * lss1_27[k]
                   + f_3 * pc_y[k] * lsp_82[k];

        t_166[k] = f_3 * pc_y[k] * lsp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_x, pc_x, pc_z, ksd0_168, ksp_62, ksp_84, \
                         ksp_85, ksd1_168, lss0_27, lss1_27, lsp_83, \
                         lsp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_7 * ksp_62[k]
                   + f_1 * lss0_27[k]
                   - f_2 * lss1_27[k]
                   + f_3 * pc_z[k] * lsp_83[k];

        t_168[k] = pa_x[k] * ksd0_168[k]
                   + f_8 * ksp_84[k]
                   - f_4 * pc_x[k] * ksd1_168[k];

        t_169[k] = f_6 * ksp_85[k]
                   + f_3 * pc_x[k] * lsp_85[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pa_x, pc_x, pc_z, ksd0_171, ksd0_173, \
                         ksd1_171, ksd1_173, lsp_84, lsp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_3 * pc_z[k] * lsp_84[k];

        t_171[k] = pa_x[k] * ksd0_171[k]
                   - f_4 * pc_x[k] * ksd1_171[k];

        t_172[k] = f_3 * pc_z[k] * lsp_85[k];

        t_173[k] = pa_x[k] * ksd0_173[k]
                   - f_4 * pc_x[k] * ksd1_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pa_x, pa_z, pc_x, pc_z, ksd0_126, \
                         ksd0_177, ksp_88, ksp_89, ksd1_126, ksd1_177, lsp_88, \
                         lsp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_z[k] * ksd0_126[k]
                   - f_4 * pc_z[k] * ksd1_126[k];

        t_175[k] = f_6 * ksp_88[k]
                   + f_3 * pc_x[k] * lsp_88[k];

        t_176[k] = f_6 * ksp_89[k]
                   + f_3 * pc_x[k] * lsp_89[k];

        t_177[k] = pa_x[k] * ksd0_177[k]
                   - f_4 * pc_x[k] * ksd1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pa_x, pc_x, pc_y, ksd0_179, ksd0_180, \
                         ksp_68, ksp_90, ksp_91, ksd1_179, ksd1_180, lsp_89, \
                         lsp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_7 * ksp_68[k]
                   + f_3 * pc_y[k] * lsp_89[k];

        t_179[k] = pa_x[k] * ksd0_179[k]
                   - f_4 * pc_x[k] * ksd1_179[k];

        t_180[k] = pa_x[k] * ksd0_180[k]
                   + f_8 * ksp_90[k]
                   - f_4 * pc_x[k] * ksd1_180[k];

        t_181[k] = f_6 * ksp_91[k]
                   + f_3 * pc_x[k] * lsp_91[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pa_x, pc_x, pc_y, ksd0_183, ksd0_185, \
                         ksp_71, ksp_92, ksd1_183, ksd1_185, lsp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_6 * ksp_92[k]
                   + f_3 * pc_x[k] * lsp_92[k];

        t_183[k] = pa_x[k] * ksd0_183[k]
                   - f_4 * pc_x[k] * ksd1_183[k];

        t_184[k] = f_9 * ksp_71[k]
                   + f_3 * pc_y[k] * lsp_92[k];

        t_185[k] = pa_x[k] * ksd0_185[k]
                   - f_4 * pc_x[k] * ksd1_185[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pa_x, pc_x, ksd0_186, ksd0_189, ksp_93, \
                         ksp_94, ksp_95, ksd1_186, ksd1_189, lsp_94, \
                         lsp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_x[k] * ksd0_186[k]
                   + f_8 * ksp_93[k]
                   - f_4 * pc_x[k] * ksd1_186[k];

        t_187[k] = f_6 * ksp_94[k]
                   + f_3 * pc_x[k] * lsp_94[k];

        t_188[k] = f_6 * ksp_95[k]
                   + f_3 * pc_x[k] * lsp_95[k];

        t_189[k] = pa_x[k] * ksd0_189[k]
                   - f_4 * pc_x[k] * ksd1_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_x, pc_x, pc_y, ksd0_191, ksd0_192, \
                         ksp_74, ksp_96, ksp_97, ksd1_191, ksd1_192, lsp_95, \
                         lsp_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_11 * ksp_74[k]
                   + f_3 * pc_y[k] * lsp_95[k];

        t_191[k] = pa_x[k] * ksd0_191[k]
                   - f_4 * pc_x[k] * ksd1_191[k];

        t_192[k] = pa_x[k] * ksd0_192[k]
                   + f_8 * ksp_96[k]
                   - f_4 * pc_x[k] * ksd1_192[k];

        t_193[k] = f_6 * ksp_97[k]
                   + f_3 * pc_x[k] * lsp_97[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_x, pc_x, pc_y, ksd0_195, ksd0_197, \
                         ksp_77, ksp_98, ksd1_195, ksd1_197, lsp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_6 * ksp_98[k]
                   + f_3 * pc_x[k] * lsp_98[k];

        t_195[k] = pa_x[k] * ksd0_195[k]
                   - f_4 * pc_x[k] * ksd1_195[k];

        t_196[k] = f_10 * ksp_77[k]
                   + f_3 * pc_y[k] * lsp_98[k];

        t_197[k] = pa_x[k] * ksd0_197[k]
                   - f_4 * pc_x[k] * ksd1_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_x, pc_x, ksd0_198, ksd0_201, ksp_99, \
                         ksp_100, ksp_101, ksd1_198, ksd1_201, lsp_100, \
                         lsp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pa_x[k] * ksd0_198[k]
                   + f_8 * ksp_99[k]
                   - f_4 * pc_x[k] * ksd1_198[k];

        t_199[k] = f_6 * ksp_100[k]
                   + f_3 * pc_x[k] * lsp_100[k];

        t_200[k] = f_6 * ksp_101[k]
                   + f_3 * pc_x[k] * lsp_101[k];

        t_201[k] = pa_x[k] * ksd0_201[k]
                   - f_4 * pc_x[k] * ksd1_201[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pa_x, pa_y, pc_x, pc_y, ksd0_162, \
                         ksd0_203, ksp_80, ksp_103, ksd1_162, ksd1_203, lsp_101, \
                         lsp_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_8 * ksp_80[k]
                   + f_3 * pc_y[k] * lsp_101[k];

        t_203[k] = pa_x[k] * ksd0_203[k]
                   - f_4 * pc_x[k] * ksd1_203[k];

        t_204[k] = pa_y[k] * ksd0_162[k]
                   - f_4 * pc_y[k] * ksd1_162[k];

        t_205[k] = f_6 * ksp_103[k]
                   + f_3 * pc_x[k] * lsp_103[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pc_x, pc_y, ksd0_207, ksd0_209, \
                         ksp_83, ksp_104, ksd1_207, ksd1_209, lsp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_6 * ksp_104[k]
                   + f_3 * pc_x[k] * lsp_104[k];

        t_207[k] = pa_x[k] * ksd0_207[k]
                   - f_4 * pc_x[k] * ksd1_207[k];

        t_208[k] = f_6 * ksp_83[k]
                   + f_3 * pc_y[k] * lsp_104[k];

        t_209[k] = pa_x[k] * ksd0_209[k]
                   - f_4 * pc_x[k] * ksd1_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pa_x, pc_x, pc_y, ksd0_210, \
                         ksd0_213, ksp_105, ksp_107, ksd1_210, ksd1_213, lsp_105, \
                         lsp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = pa_x[k] * ksd0_210[k]
                   + f_8 * ksp_105[k]
                   - f_4 * pc_x[k] * ksd1_210[k];

        t_211[k] = f_3 * pc_y[k] * lsp_105[k];

        t_212[k] = f_6 * ksp_107[k]
                   + f_3 * pc_x[k] * lsp_107[k];

        t_213[k] = pa_x[k] * ksd0_213[k]
                   - f_4 * pc_x[k] * ksd1_213[k];

        t_214[k] = f_3 * pc_y[k] * lsp_107[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pa_x, pc_x, pc_y, ksd0_215, \
                         ksp_85, ksd1_215, lss0_36, lss1_36, lsp_108, lsp_109, \
                         lsp_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pa_x[k] * ksd0_215[k]
                   - f_4 * pc_x[k] * ksd1_215[k];

        t_216[k] = f_1 * lss0_36[k]
                   - f_2 * lss1_36[k]
                   + f_3 * pc_x[k] * lsp_108[k];

        t_217[k] = f_3 * pc_x[k] * lsp_109[k];

        t_218[k] = f_3 * pc_x[k] * lsp_110[k];

        t_219[k] = f_0 * ksp_85[k]
                   + f_1 * lss0_36[k]
                   - f_2 * lss1_36[k]
                   + f_3 * pc_y[k] * lsp_109[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pa_z, pc_x, pc_z, ksd0_168, \
                         ksd1_168, lss0_36, lss1_36, lsp_109, lsp_110, lsp_112, \
                         lsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_3 * pc_z[k] * lsp_109[k];

        t_221[k] = f_1 * lss0_36[k]
                   - f_2 * lss1_36[k]
                   + f_3 * pc_z[k] * lsp_110[k];

        t_222[k] = pa_z[k] * ksd0_168[k]
                   - f_4 * pc_z[k] * ksd1_168[k];

        t_223[k] = f_3 * pc_x[k] * lsp_112[k];

        t_224[k] = f_3 * pc_x[k] * lsp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_z, pc_y, pc_z, ksd0_171, ksp_86, ksp_89, \
                         ksd1_171, lss0_37, lss1_37, lsp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_z[k] * ksd0_171[k]
                   - f_4 * pc_z[k] * ksd1_171[k];

        t_226[k] = f_5 * ksp_89[k]
                   + f_3 * pc_y[k] * lsp_113[k];

        t_227[k] = f_6 * ksp_86[k]
                   + f_1 * lss0_37[k]
                   - f_2 * lss1_37[k]
                   + f_3 * pc_z[k] * lsp_113[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, t_232, pc_x, pc_y, ksp_91, ksp_92, \
                         lss0_38, lss1_38, lsp_114, lsp_115, lsp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_1 * lss0_38[k]
                   - f_2 * lss1_38[k]
                   + f_3 * pc_x[k] * lsp_114[k];

        t_229[k] = f_3 * pc_x[k] * lsp_115[k];

        t_230[k] = f_3 * pc_x[k] * lsp_116[k];

        t_231[k] = f_7 * ksp_91[k]
                   + f_1 * lss0_38[k]
                   - f_2 * lss1_38[k]
                   + f_3 * pc_y[k] * lsp_115[k];

        t_232[k] = f_7 * ksp_92[k]
                   + f_3 * pc_y[k] * lsp_116[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pc_x, pc_z, ksp_89, lss0_38, lss0_39, \
                         lss1_38, lss1_39, lsp_116, lsp_117, lsp_118, \
                         lsp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_8 * ksp_89[k]
                   + f_1 * lss0_38[k]
                   - f_2 * lss1_38[k]
                   + f_3 * pc_z[k] * lsp_116[k];

        t_234[k] = f_1 * lss0_39[k]
                   - f_2 * lss1_39[k]
                   + f_3 * pc_x[k] * lsp_117[k];

        t_235[k] = f_3 * pc_x[k] * lsp_118[k];

        t_236[k] = f_3 * pc_x[k] * lsp_119[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pc_y, pc_z, ksp_92, ksp_94, ksp_95, lss0_39, \
                         lss1_39, lsp_118, lsp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_9 * ksp_94[k]
                   + f_1 * lss0_39[k]
                   - f_2 * lss1_39[k]
                   + f_3 * pc_y[k] * lsp_118[k];

        t_238[k] = f_9 * ksp_95[k]
                   + f_3 * pc_y[k] * lsp_119[k];

        t_239[k] = f_10 * ksp_92[k]
                   + f_1 * lss0_39[k]
                   - f_2 * lss1_39[k]
                   + f_3 * pc_z[k] * lsp_119[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, ksp_97, ksp_98, \
                         lss0_40, lss1_40, lsp_120, lsp_121, lsp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_1 * lss0_40[k]
                   - f_2 * lss1_40[k]
                   + f_3 * pc_x[k] * lsp_120[k];

        t_241[k] = f_3 * pc_x[k] * lsp_121[k];

        t_242[k] = f_3 * pc_x[k] * lsp_122[k];

        t_243[k] = f_11 * ksp_97[k]
                   + f_1 * lss0_40[k]
                   - f_2 * lss1_40[k]
                   + f_3 * pc_y[k] * lsp_121[k];

        t_244[k] = f_11 * ksp_98[k]
                   + f_3 * pc_y[k] * lsp_122[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_z, ksp_95, lss0_40, lss0_41, \
                         lss1_40, lss1_41, lsp_122, lsp_123, lsp_124, \
                         lsp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * ksp_95[k]
                   + f_1 * lss0_40[k]
                   - f_2 * lss1_40[k]
                   + f_3 * pc_z[k] * lsp_122[k];

        t_246[k] = f_1 * lss0_41[k]
                   - f_2 * lss1_41[k]
                   + f_3 * pc_x[k] * lsp_123[k];

        t_247[k] = f_3 * pc_x[k] * lsp_124[k];

        t_248[k] = f_3 * pc_x[k] * lsp_125[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_y, pc_z, ksp_98, ksp_100, ksp_101, lss0_41, \
                         lss1_41, lsp_124, lsp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_10 * ksp_100[k]
                   + f_1 * lss0_41[k]
                   - f_2 * lss1_41[k]
                   + f_3 * pc_y[k] * lsp_124[k];

        t_250[k] = f_10 * ksp_101[k]
                   + f_3 * pc_y[k] * lsp_125[k];

        t_251[k] = f_9 * ksp_98[k]
                   + f_1 * lss0_41[k]
                   - f_2 * lss1_41[k]
                   + f_3 * pc_z[k] * lsp_125[k];
    }
}

static auto
compute_prim_lsd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksd0,
                                                          const size_t ksp, const size_t ksd1,
                                                          const size_t lss0, const size_t lss1,
                                                          const size_t lsp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.0 / q;
    const auto f_8 = 1.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksd0_210 = buffer.data(ksd0 + 210);
    const auto *ksd0_213 = buffer.data(ksd0 + 213);
    const auto *ksd0_215 = buffer.data(ksd0 + 215);

    const auto *ksp_101 = buffer.data(ksp + 101);
    const auto *ksp_103 = buffer.data(ksp + 103);
    const auto *ksp_104 = buffer.data(ksp + 104);
    const auto *ksp_106 = buffer.data(ksp + 106);
    const auto *ksp_107 = buffer.data(ksp + 107);

    const auto *ksd1_210 = buffer.data(ksd1 + 210);
    const auto *ksd1_213 = buffer.data(ksd1 + 213);
    const auto *ksd1_215 = buffer.data(ksd1 + 215);

    const auto *lss0_42 = buffer.data(lss0 + 42);
    const auto *lss0_44 = buffer.data(lss0 + 44);

    const auto *lss1_42 = buffer.data(lss1 + 42);
    const auto *lss1_44 = buffer.data(lss1 + 44);

    const auto *lsp_126 = buffer.data(lsp + 126);
    const auto *lsp_127 = buffer.data(lsp + 127);
    const auto *lsp_128 = buffer.data(lsp + 128);
    const auto *lsp_130 = buffer.data(lsp + 130);
    const auto *lsp_131 = buffer.data(lsp + 131);
    const auto *lsp_132 = buffer.data(lsp + 132);
    const auto *lsp_133 = buffer.data(lsp + 133);
    const auto *lsp_134 = buffer.data(lsp + 134);

#pragma omp simd aligned(t_252, t_253, t_254, t_255, t_256, pc_x, pc_y, ksp_103, ksp_104, \
                         lss0_42, lss1_42, lsp_126, lsp_127, lsp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_1 * lss0_42[k]
                   - f_2 * lss1_42[k]
                   + f_3 * pc_x[k] * lsp_126[k];

        t_253[k] = f_3 * pc_x[k] * lsp_127[k];

        t_254[k] = f_3 * pc_x[k] * lsp_128[k];

        t_255[k] = f_8 * ksp_103[k]
                   + f_1 * lss0_42[k]
                   - f_2 * lss1_42[k]
                   + f_3 * pc_y[k] * lsp_127[k];

        t_256[k] = f_8 * ksp_104[k]
                   + f_3 * pc_y[k] * lsp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_y, pc_x, pc_y, pc_z, ksd0_210, \
                         ksp_101, ksd1_210, lss0_42, lss1_42, lsp_128, lsp_130, \
                         lsp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_7 * ksp_101[k]
                   + f_1 * lss0_42[k]
                   - f_2 * lss1_42[k]
                   + f_3 * pc_z[k] * lsp_128[k];

        t_258[k] = pa_y[k] * ksd0_210[k]
                   - f_4 * pc_y[k] * ksd1_210[k];

        t_259[k] = f_3 * pc_x[k] * lsp_130[k];

        t_260[k] = f_3 * pc_x[k] * lsp_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pa_y, pc_y, ksd0_213, ksd0_215, ksp_106, \
                         ksp_107, ksd1_213, ksd1_215, lsp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_y[k] * ksd0_213[k]
                   + f_8 * ksp_106[k]
                   - f_4 * pc_y[k] * ksd1_213[k];

        t_262[k] = f_6 * ksp_107[k]
                   + f_3 * pc_y[k] * lsp_131[k];

        t_263[k] = pa_y[k] * ksd0_215[k]
                   - f_4 * pc_y[k] * ksd1_215[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, t_269, pc_x, pc_y, pc_z, ksp_107, \
                         lss0_44, lss1_44, lsp_132, lsp_133, lsp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_1 * lss0_44[k]
                   - f_2 * lss1_44[k]
                   + f_3 * pc_x[k] * lsp_132[k];

        t_265[k] = f_3 * pc_x[k] * lsp_133[k];

        t_266[k] = f_3 * pc_x[k] * lsp_134[k];

        t_267[k] = f_1 * lss0_44[k]
                   - f_2 * lss1_44[k]
                   + f_3 * pc_y[k] * lsp_133[k];

        t_268[k] = f_3 * pc_y[k] * lsp_134[k];

        t_269[k] = f_0 * ksp_107[k]
                   + f_1 * lss0_44[k]
                   - f_2 * lss1_44[k]
                   + f_3 * pc_z[k] * lsp_134[k];
    }
}

auto
compute_prim_lsd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksd0, const size_t ksp,
                                                   const size_t ksd1, const size_t lss0,
                                                   const size_t lss1, const size_t lsp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksd0, ksp,
                                                              ksd1, lss0, lss1, lsp, ncols,
                                                              gamma, p, q);

    compute_prim_lsd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksd0, ksp,
                                                              ksd1, lss0, lss1, lsp, ncols,
                                                              gamma, p, q);

    compute_prim_lsd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksd0, ksp,
                                                              ksd1, lss0, lss1, lsp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
