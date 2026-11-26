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


#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_msd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsd0,
                                                          const size_t lsp, const size_t lsd1,
                                                          const size_t mss0, const size_t mss1,
                                                          const size_t msp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsd0_0 = buffer.data(lsd0 + 0);
    const auto *lsd0_3 = buffer.data(lsd0 + 3);
    const auto *lsd0_5 = buffer.data(lsd0 + 5);
    const auto *lsd0_9 = buffer.data(lsd0 + 9);
    const auto *lsd0_12 = buffer.data(lsd0 + 12);
    const auto *lsd0_17 = buffer.data(lsd0 + 17);
    const auto *lsd0_18 = buffer.data(lsd0 + 18);
    const auto *lsd0_21 = buffer.data(lsd0 + 21);
    const auto *lsd0_30 = buffer.data(lsd0 + 30);
    const auto *lsd0_35 = buffer.data(lsd0 + 35);
    const auto *lsd0_36 = buffer.data(lsd0 + 36);
    const auto *lsd0_39 = buffer.data(lsd0 + 39);
    const auto *lsd0_54 = buffer.data(lsd0 + 54);
    const auto *lsd0_59 = buffer.data(lsd0 + 59);
    const auto *lsd0_60 = buffer.data(lsd0 + 60);
    const auto *lsd0_63 = buffer.data(lsd0 + 63);
    const auto *lsd0_84 = buffer.data(lsd0 + 84);
    const auto *lsd0_89 = buffer.data(lsd0 + 89);

    const auto *lsp_0 = buffer.data(lsp + 0);
    const auto *lsp_1 = buffer.data(lsp + 1);
    const auto *lsp_2 = buffer.data(lsp + 2);
    const auto *lsp_4 = buffer.data(lsp + 4);
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
    const auto *lsp_62 = buffer.data(lsp + 62);

    const auto *lsd1_0 = buffer.data(lsd1 + 0);
    const auto *lsd1_3 = buffer.data(lsd1 + 3);
    const auto *lsd1_5 = buffer.data(lsd1 + 5);
    const auto *lsd1_9 = buffer.data(lsd1 + 9);
    const auto *lsd1_12 = buffer.data(lsd1 + 12);
    const auto *lsd1_17 = buffer.data(lsd1 + 17);
    const auto *lsd1_18 = buffer.data(lsd1 + 18);
    const auto *lsd1_21 = buffer.data(lsd1 + 21);
    const auto *lsd1_30 = buffer.data(lsd1 + 30);
    const auto *lsd1_35 = buffer.data(lsd1 + 35);
    const auto *lsd1_36 = buffer.data(lsd1 + 36);
    const auto *lsd1_39 = buffer.data(lsd1 + 39);
    const auto *lsd1_54 = buffer.data(lsd1 + 54);
    const auto *lsd1_59 = buffer.data(lsd1 + 59);
    const auto *lsd1_60 = buffer.data(lsd1 + 60);
    const auto *lsd1_63 = buffer.data(lsd1 + 63);
    const auto *lsd1_84 = buffer.data(lsd1 + 84);
    const auto *lsd1_89 = buffer.data(lsd1 + 89);

    const auto *mss0_0 = buffer.data(mss0 + 0);
    const auto *mss0_1 = buffer.data(mss0 + 1);
    const auto *mss0_2 = buffer.data(mss0 + 2);
    const auto *mss0_3 = buffer.data(mss0 + 3);
    const auto *mss0_5 = buffer.data(mss0 + 5);
    const auto *mss0_6 = buffer.data(mss0 + 6);
    const auto *mss0_7 = buffer.data(mss0 + 7);
    const auto *mss0_8 = buffer.data(mss0 + 8);
    const auto *mss0_9 = buffer.data(mss0 + 9);
    const auto *mss0_10 = buffer.data(mss0 + 10);
    const auto *mss0_11 = buffer.data(mss0 + 11);
    const auto *mss0_12 = buffer.data(mss0 + 12);
    const auto *mss0_13 = buffer.data(mss0 + 13);
    const auto *mss0_14 = buffer.data(mss0 + 14);
    const auto *mss0_15 = buffer.data(mss0 + 15);
    const auto *mss0_16 = buffer.data(mss0 + 16);
    const auto *mss0_17 = buffer.data(mss0 + 17);
    const auto *mss0_18 = buffer.data(mss0 + 18);
    const auto *mss0_19 = buffer.data(mss0 + 19);
    const auto *mss0_20 = buffer.data(mss0 + 20);

    const auto *mss1_0 = buffer.data(mss1 + 0);
    const auto *mss1_1 = buffer.data(mss1 + 1);
    const auto *mss1_2 = buffer.data(mss1 + 2);
    const auto *mss1_3 = buffer.data(mss1 + 3);
    const auto *mss1_5 = buffer.data(mss1 + 5);
    const auto *mss1_6 = buffer.data(mss1 + 6);
    const auto *mss1_7 = buffer.data(mss1 + 7);
    const auto *mss1_8 = buffer.data(mss1 + 8);
    const auto *mss1_9 = buffer.data(mss1 + 9);
    const auto *mss1_10 = buffer.data(mss1 + 10);
    const auto *mss1_11 = buffer.data(mss1 + 11);
    const auto *mss1_12 = buffer.data(mss1 + 12);
    const auto *mss1_13 = buffer.data(mss1 + 13);
    const auto *mss1_14 = buffer.data(mss1 + 14);
    const auto *mss1_15 = buffer.data(mss1 + 15);
    const auto *mss1_16 = buffer.data(mss1 + 16);
    const auto *mss1_17 = buffer.data(mss1 + 17);
    const auto *mss1_18 = buffer.data(mss1 + 18);
    const auto *mss1_19 = buffer.data(mss1 + 19);
    const auto *mss1_20 = buffer.data(mss1 + 20);

    const auto *msp_0 = buffer.data(msp + 0);
    const auto *msp_1 = buffer.data(msp + 1);
    const auto *msp_2 = buffer.data(msp + 2);
    const auto *msp_3 = buffer.data(msp + 3);
    const auto *msp_4 = buffer.data(msp + 4);
    const auto *msp_6 = buffer.data(msp + 6);
    const auto *msp_8 = buffer.data(msp + 8);
    const auto *msp_9 = buffer.data(msp + 9);
    const auto *msp_10 = buffer.data(msp + 10);
    const auto *msp_11 = buffer.data(msp + 11);
    const auto *msp_13 = buffer.data(msp + 13);
    const auto *msp_14 = buffer.data(msp + 14);
    const auto *msp_15 = buffer.data(msp + 15);
    const auto *msp_16 = buffer.data(msp + 16);
    const auto *msp_17 = buffer.data(msp + 17);
    const auto *msp_18 = buffer.data(msp + 18);
    const auto *msp_19 = buffer.data(msp + 19);
    const auto *msp_20 = buffer.data(msp + 20);
    const auto *msp_22 = buffer.data(msp + 22);
    const auto *msp_23 = buffer.data(msp + 23);
    const auto *msp_25 = buffer.data(msp + 25);
    const auto *msp_26 = buffer.data(msp + 26);
    const auto *msp_27 = buffer.data(msp + 27);
    const auto *msp_28 = buffer.data(msp + 28);
    const auto *msp_29 = buffer.data(msp + 29);
    const auto *msp_30 = buffer.data(msp + 30);
    const auto *msp_31 = buffer.data(msp + 31);
    const auto *msp_32 = buffer.data(msp + 32);
    const auto *msp_34 = buffer.data(msp + 34);
    const auto *msp_35 = buffer.data(msp + 35);
    const auto *msp_36 = buffer.data(msp + 36);
    const auto *msp_37 = buffer.data(msp + 37);
    const auto *msp_38 = buffer.data(msp + 38);
    const auto *msp_40 = buffer.data(msp + 40);
    const auto *msp_41 = buffer.data(msp + 41);
    const auto *msp_42 = buffer.data(msp + 42);
    const auto *msp_43 = buffer.data(msp + 43);
    const auto *msp_44 = buffer.data(msp + 44);
    const auto *msp_45 = buffer.data(msp + 45);
    const auto *msp_46 = buffer.data(msp + 46);
    const auto *msp_47 = buffer.data(msp + 47);
    const auto *msp_49 = buffer.data(msp + 49);
    const auto *msp_50 = buffer.data(msp + 50);
    const auto *msp_51 = buffer.data(msp + 51);
    const auto *msp_52 = buffer.data(msp + 52);
    const auto *msp_53 = buffer.data(msp + 53);
    const auto *msp_54 = buffer.data(msp + 54);
    const auto *msp_55 = buffer.data(msp + 55);
    const auto *msp_56 = buffer.data(msp + 56);
    const auto *msp_58 = buffer.data(msp + 58);
    const auto *msp_59 = buffer.data(msp + 59);
    const auto *msp_60 = buffer.data(msp + 60);
    const auto *msp_61 = buffer.data(msp + 61);
    const auto *msp_62 = buffer.data(msp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, lsp_0, mss0_0, \
                         mss1_0, msp_0, msp_1, msp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lsp_0[k]
                 + f_1 * mss0_0[k]
                 - f_2 * mss1_0[k]
                 + f_3 * pc_x[k] * msp_0[k];

        t_1[k] = f_3 * pc_y[k] * msp_0[k];

        t_2[k] = f_3 * pc_z[k] * msp_0[k];

        t_3[k] = f_1 * mss0_0[k]
                 - f_2 * mss1_0[k]
                 + f_3 * pc_y[k] * msp_1[k];

        t_4[k] = f_3 * pc_y[k] * msp_2[k];

        t_5[k] = f_1 * mss0_0[k]
                 - f_2 * mss1_0[k]
                 + f_3 * pc_z[k] * msp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, lsd0_0, lsp_1, lsp_4, \
                         lsd1_0, mss0_1, mss1_1, msp_3, msp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * lsd0_0[k]
                 - f_4 * pc_y[k] * lsd1_0[k];

        t_7[k] = f_5 * lsp_4[k]
                 + f_3 * pc_x[k] * msp_4[k];

        t_8[k] = f_3 * pc_z[k] * msp_3[k];

        t_9[k] = f_6 * lsp_1[k]
                 + f_1 * mss0_1[k]
                 - f_2 * mss1_1[k]
                 + f_3 * pc_y[k] * msp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, lsd0_0, lsd0_5, \
                         lsd1_0, lsd1_5, msp_4, msp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * msp_4[k];

        t_11[k] = pa_y[k] * lsd0_5[k]
                  - f_4 * pc_y[k] * lsd1_5[k];

        t_12[k] = pa_z[k] * lsd0_0[k]
                  - f_4 * pc_z[k] * lsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * msp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, lsd0_3, lsp_2, lsp_8, \
                         lsd1_3, mss0_2, mss1_2, msp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * lsp_8[k]
                  + f_3 * pc_x[k] * msp_8[k];

        t_15[k] = pa_z[k] * lsd0_3[k]
                  - f_4 * pc_z[k] * lsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * msp_8[k];

        t_17[k] = f_6 * lsp_2[k]
                  + f_1 * mss0_2[k]
                  - f_2 * mss1_2[k]
                  + f_3 * pc_z[k] * msp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, lsp_4, lsp_9, lsp_10, \
                         mss0_3, mss1_3, msp_9, msp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * lsp_9[k]
                  + f_1 * mss0_3[k]
                  - f_2 * mss1_3[k]
                  + f_3 * pc_x[k] * msp_9[k];

        t_19[k] = f_7 * lsp_10[k]
                  + f_3 * pc_x[k] * msp_10[k];

        t_20[k] = f_3 * pc_z[k] * msp_9[k];

        t_21[k] = f_8 * lsp_4[k]
                  + f_1 * mss0_3[k]
                  - f_2 * mss1_3[k]
                  + f_3 * pc_y[k] * msp_10[k];

        t_22[k] = f_3 * pc_z[k] * msp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, lsd0_12, lsp_13, lsd1_12, \
                         mss0_3, mss1_3, msp_11, msp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * mss0_3[k]
                  - f_2 * mss1_3[k]
                  + f_3 * pc_z[k] * msp_11[k];

        t_24[k] = pa_y[k] * lsd0_12[k]
                  - f_4 * pc_y[k] * lsd1_12[k];

        t_25[k] = f_7 * lsp_13[k]
                  + f_3 * pc_x[k] * msp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, lsd0_9, \
                         lsd0_17, lsp_8, lsp_14, lsd1_9, lsd1_17, \
                         msp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * lsp_14[k]
                  + f_3 * pc_x[k] * msp_14[k];

        t_27[k] = pa_z[k] * lsd0_9[k]
                  - f_4 * pc_z[k] * lsd1_9[k];

        t_28[k] = f_6 * lsp_8[k]
                  + f_3 * pc_y[k] * msp_14[k];

        t_29[k] = pa_y[k] * lsd0_17[k]
                  - f_4 * pc_y[k] * lsd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, lsp_15, lsp_17, mss0_5, \
                         mss1_5, msp_15, msp_16, msp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * lsp_15[k]
                  + f_1 * mss0_5[k]
                  - f_2 * mss1_5[k]
                  + f_3 * pc_x[k] * msp_15[k];

        t_31[k] = f_3 * pc_y[k] * msp_15[k];

        t_32[k] = f_7 * lsp_17[k]
                  + f_3 * pc_x[k] * msp_17[k];

        t_33[k] = f_1 * mss0_5[k]
                  - f_2 * mss1_5[k]
                  + f_3 * pc_y[k] * msp_16[k];

        t_34[k] = f_3 * pc_y[k] * msp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, lsp_8, lsp_18, lsp_19, mss0_5, \
                         mss0_6, mss1_5, mss1_6, msp_17, msp_18, \
                         msp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * lsp_8[k]
                  + f_1 * mss0_5[k]
                  - f_2 * mss1_5[k]
                  + f_3 * pc_z[k] * msp_17[k];

        t_36[k] = f_9 * lsp_18[k]
                  + f_1 * mss0_6[k]
                  - f_2 * mss1_6[k]
                  + f_3 * pc_x[k] * msp_18[k];

        t_37[k] = f_9 * lsp_19[k]
                  + f_3 * pc_x[k] * msp_19[k];

        t_38[k] = f_3 * pc_z[k] * msp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, lsd0_18, lsp_10, lsd1_18, \
                         mss0_6, mss1_6, msp_19, msp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_10 * lsp_10[k]
                  + f_1 * mss0_6[k]
                  - f_2 * mss1_6[k]
                  + f_3 * pc_y[k] * msp_19[k];

        t_40[k] = f_3 * pc_z[k] * msp_19[k];

        t_41[k] = f_1 * mss0_6[k]
                  - f_2 * mss1_6[k]
                  + f_3 * pc_z[k] * msp_20[k];

        t_42[k] = pa_z[k] * lsd0_18[k]
                  - f_4 * pc_z[k] * lsd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, lsd0_21, lsp_14, \
                         lsp_22, lsp_23, lsd1_21, msp_22, msp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * lsp_22[k]
                  + f_3 * pc_x[k] * msp_22[k];

        t_44[k] = f_9 * lsp_23[k]
                  + f_3 * pc_x[k] * msp_23[k];

        t_45[k] = pa_z[k] * lsd0_21[k]
                  - f_4 * pc_z[k] * lsd1_21[k];

        t_46[k] = f_8 * lsp_14[k]
                  + f_3 * pc_y[k] * msp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, lsd0_30, lsp_11, lsp_25, \
                         lsd1_30, mss0_7, mss1_7, msp_23, msp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * lsp_11[k]
                  + f_1 * mss0_7[k]
                  - f_2 * mss1_7[k]
                  + f_3 * pc_z[k] * msp_23[k];

        t_48[k] = pa_y[k] * lsd0_30[k]
                  - f_4 * pc_y[k] * lsd1_30[k];

        t_49[k] = f_9 * lsp_25[k]
                  + f_3 * pc_x[k] * msp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, lsd0_35, lsp_16, lsp_17, \
                         lsp_26, lsd1_35, mss0_8, mss1_8, msp_25, \
                         msp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * lsp_26[k]
                  + f_3 * pc_x[k] * msp_26[k];

        t_51[k] = f_6 * lsp_16[k]
                  + f_1 * mss0_8[k]
                  - f_2 * mss1_8[k]
                  + f_3 * pc_y[k] * msp_25[k];

        t_52[k] = f_6 * lsp_17[k]
                  + f_3 * pc_y[k] * msp_26[k];

        t_53[k] = pa_y[k] * lsd0_35[k]
                  - f_4 * pc_y[k] * lsd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, lsp_27, lsp_29, mss0_9, \
                         mss1_9, msp_27, msp_28, msp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * lsp_27[k]
                  + f_1 * mss0_9[k]
                  - f_2 * mss1_9[k]
                  + f_3 * pc_x[k] * msp_27[k];

        t_55[k] = f_3 * pc_y[k] * msp_27[k];

        t_56[k] = f_9 * lsp_29[k]
                  + f_3 * pc_x[k] * msp_29[k];

        t_57[k] = f_1 * mss0_9[k]
                  - f_2 * mss1_9[k]
                  + f_3 * pc_y[k] * msp_28[k];

        t_58[k] = f_3 * pc_y[k] * msp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, lsp_17, lsp_30, lsp_31, mss0_9, \
                         mss0_10, mss1_9, mss1_10, msp_29, msp_30, \
                         msp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_10 * lsp_17[k]
                  + f_1 * mss0_9[k]
                  - f_2 * mss1_9[k]
                  + f_3 * pc_z[k] * msp_29[k];

        t_60[k] = f_11 * lsp_30[k]
                  + f_1 * mss0_10[k]
                  - f_2 * mss1_10[k]
                  + f_3 * pc_x[k] * msp_30[k];

        t_61[k] = f_11 * lsp_31[k]
                  + f_3 * pc_x[k] * msp_31[k];

        t_62[k] = f_3 * pc_z[k] * msp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, lsd0_36, lsp_19, lsd1_36, \
                         mss0_10, mss1_10, msp_31, msp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_12 * lsp_19[k]
                  + f_1 * mss0_10[k]
                  - f_2 * mss1_10[k]
                  + f_3 * pc_y[k] * msp_31[k];

        t_64[k] = f_3 * pc_z[k] * msp_31[k];

        t_65[k] = f_1 * mss0_10[k]
                  - f_2 * mss1_10[k]
                  + f_3 * pc_z[k] * msp_32[k];

        t_66[k] = pa_z[k] * lsd0_36[k]
                  - f_4 * pc_z[k] * lsd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, lsd0_39, lsp_23, \
                         lsp_34, lsp_35, lsd1_39, msp_34, msp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_11 * lsp_34[k]
                  + f_3 * pc_x[k] * msp_34[k];

        t_68[k] = f_11 * lsp_35[k]
                  + f_3 * pc_x[k] * msp_35[k];

        t_69[k] = pa_z[k] * lsd0_39[k]
                  - f_4 * pc_z[k] * lsd1_39[k];

        t_70[k] = f_10 * lsp_23[k]
                  + f_3 * pc_y[k] * msp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, lsp_20, lsp_36, lsp_37, mss0_11, \
                         mss0_12, mss1_11, mss1_12, msp_35, msp_36, \
                         msp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * lsp_20[k]
                  + f_1 * mss0_11[k]
                  - f_2 * mss1_11[k]
                  + f_3 * pc_z[k] * msp_35[k];

        t_72[k] = f_11 * lsp_36[k]
                  + f_1 * mss0_12[k]
                  - f_2 * mss1_12[k]
                  + f_3 * pc_x[k] * msp_36[k];

        t_73[k] = f_11 * lsp_37[k]
                  + f_3 * pc_x[k] * msp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, lsp_23, lsp_25, lsp_26, \
                         lsp_38, mss0_12, mss1_12, msp_37, msp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * lsp_38[k]
                  + f_3 * pc_x[k] * msp_38[k];

        t_75[k] = f_8 * lsp_25[k]
                  + f_1 * mss0_12[k]
                  - f_2 * mss1_12[k]
                  + f_3 * pc_y[k] * msp_37[k];

        t_76[k] = f_8 * lsp_26[k]
                  + f_3 * pc_y[k] * msp_38[k];

        t_77[k] = f_8 * lsp_23[k]
                  + f_1 * mss0_12[k]
                  - f_2 * mss1_12[k]
                  + f_3 * pc_z[k] * msp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, lsd0_54, lsp_28, lsp_40, \
                         lsp_41, lsd1_54, mss0_13, mss1_13, msp_40, \
                         msp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * lsd0_54[k]
                  - f_4 * pc_y[k] * lsd1_54[k];

        t_79[k] = f_11 * lsp_40[k]
                  + f_3 * pc_x[k] * msp_40[k];

        t_80[k] = f_11 * lsp_41[k]
                  + f_3 * pc_x[k] * msp_41[k];

        t_81[k] = f_6 * lsp_28[k]
                  + f_1 * mss0_13[k]
                  - f_2 * mss1_13[k]
                  + f_3 * pc_y[k] * msp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, lsd0_59, lsp_29, lsp_42, \
                         lsd1_59, mss0_14, mss1_14, msp_41, msp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * lsp_29[k]
                  + f_3 * pc_y[k] * msp_41[k];

        t_83[k] = pa_y[k] * lsd0_59[k]
                  - f_4 * pc_y[k] * lsd1_59[k];

        t_84[k] = f_11 * lsp_42[k]
                  + f_1 * mss0_14[k]
                  - f_2 * mss1_14[k]
                  + f_3 * pc_x[k] * msp_42[k];

        t_85[k] = f_3 * pc_y[k] * msp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, lsp_29, lsp_44, mss0_14, \
                         mss1_14, msp_43, msp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_11 * lsp_44[k]
                  + f_3 * pc_x[k] * msp_44[k];

        t_87[k] = f_1 * mss0_14[k]
                  - f_2 * mss1_14[k]
                  + f_3 * pc_y[k] * msp_43[k];

        t_88[k] = f_3 * pc_y[k] * msp_44[k];

        t_89[k] = f_12 * lsp_29[k]
                  + f_1 * mss0_14[k]
                  - f_2 * mss1_14[k]
                  + f_3 * pc_z[k] * msp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, lsp_31, lsp_45, \
                         lsp_46, mss0_15, mss1_15, msp_45, msp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_12 * lsp_45[k]
                  + f_1 * mss0_15[k]
                  - f_2 * mss1_15[k]
                  + f_3 * pc_x[k] * msp_45[k];

        t_91[k] = f_12 * lsp_46[k]
                  + f_3 * pc_x[k] * msp_46[k];

        t_92[k] = f_3 * pc_z[k] * msp_45[k];

        t_93[k] = f_11 * lsp_31[k]
                  + f_1 * mss0_15[k]
                  - f_2 * mss1_15[k]
                  + f_3 * pc_y[k] * msp_46[k];

        t_94[k] = f_3 * pc_z[k] * msp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_z, pc_x, pc_z, lsd0_60, lsp_49, lsp_50, \
                         lsd1_60, mss0_15, mss1_15, msp_47, msp_49, \
                         msp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_1 * mss0_15[k]
                  - f_2 * mss1_15[k]
                  + f_3 * pc_z[k] * msp_47[k];

        t_96[k] = pa_z[k] * lsd0_60[k]
                  - f_4 * pc_z[k] * lsd1_60[k];

        t_97[k] = f_12 * lsp_49[k]
                  + f_3 * pc_x[k] * msp_49[k];

        t_98[k] = f_12 * lsp_50[k]
                  + f_3 * pc_x[k] * msp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_z, pc_y, pc_z, lsd0_63, lsp_32, lsp_35, \
                         lsd1_63, mss0_16, mss1_16, msp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_z[k] * lsd0_63[k]
                  - f_4 * pc_z[k] * lsd1_63[k];

        t_100[k] = f_12 * lsp_35[k]
                   + f_3 * pc_y[k] * msp_50[k];

        t_101[k] = f_6 * lsp_32[k]
                   + f_1 * mss0_16[k]
                   - f_2 * mss1_16[k]
                   + f_3 * pc_z[k] * msp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_y, lsp_37, lsp_51, lsp_52, \
                         lsp_53, mss0_17, mss1_17, msp_51, msp_52, \
                         msp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * lsp_51[k]
                   + f_1 * mss0_17[k]
                   - f_2 * mss1_17[k]
                   + f_3 * pc_x[k] * msp_51[k];

        t_103[k] = f_12 * lsp_52[k]
                   + f_3 * pc_x[k] * msp_52[k];

        t_104[k] = f_12 * lsp_53[k]
                   + f_3 * pc_x[k] * msp_53[k];

        t_105[k] = f_10 * lsp_37[k]
                   + f_1 * mss0_17[k]
                   - f_2 * mss1_17[k]
                   + f_3 * pc_y[k] * msp_52[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pc_x, pc_y, pc_z, lsp_35, lsp_38, lsp_54, \
                         mss0_17, mss0_18, mss1_17, mss1_18, msp_53, \
                         msp_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_10 * lsp_38[k]
                   + f_3 * pc_y[k] * msp_53[k];

        t_107[k] = f_8 * lsp_35[k]
                   + f_1 * mss0_17[k]
                   - f_2 * mss1_17[k]
                   + f_3 * pc_z[k] * msp_53[k];

        t_108[k] = f_12 * lsp_54[k]
                   + f_1 * mss0_18[k]
                   - f_2 * mss1_18[k]
                   + f_3 * pc_x[k] * msp_54[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pc_x, pc_y, lsp_40, lsp_41, lsp_55, \
                         lsp_56, mss0_18, mss1_18, msp_55, msp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_12 * lsp_55[k]
                   + f_3 * pc_x[k] * msp_55[k];

        t_110[k] = f_12 * lsp_56[k]
                   + f_3 * pc_x[k] * msp_56[k];

        t_111[k] = f_8 * lsp_40[k]
                   + f_1 * mss0_18[k]
                   - f_2 * mss1_18[k]
                   + f_3 * pc_y[k] * msp_55[k];

        t_112[k] = f_8 * lsp_41[k]
                   + f_3 * pc_y[k] * msp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, lsd0_84, lsp_38, lsp_58, \
                         lsd1_84, mss0_18, mss1_18, msp_56, msp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_10 * lsp_38[k]
                   + f_1 * mss0_18[k]
                   - f_2 * mss1_18[k]
                   + f_3 * pc_z[k] * msp_56[k];

        t_114[k] = pa_y[k] * lsd0_84[k]
                   - f_4 * pc_y[k] * lsd1_84[k];

        t_115[k] = f_12 * lsp_58[k]
                   + f_3 * pc_x[k] * msp_58[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_x, pc_y, lsd0_89, lsp_43, \
                         lsp_44, lsp_59, lsd1_89, mss0_19, mss1_19, msp_58, \
                         msp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_12 * lsp_59[k]
                   + f_3 * pc_x[k] * msp_59[k];

        t_117[k] = f_6 * lsp_43[k]
                   + f_1 * mss0_19[k]
                   - f_2 * mss1_19[k]
                   + f_3 * pc_y[k] * msp_58[k];

        t_118[k] = f_6 * lsp_44[k]
                   + f_3 * pc_y[k] * msp_59[k];

        t_119[k] = pa_y[k] * lsd0_89[k]
                   - f_4 * pc_y[k] * lsd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, lsp_60, lsp_62, \
                         mss0_20, mss1_20, msp_60, msp_61, msp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_12 * lsp_60[k]
                   + f_1 * mss0_20[k]
                   - f_2 * mss1_20[k]
                   + f_3 * pc_x[k] * msp_60[k];

        t_121[k] = f_3 * pc_y[k] * msp_60[k];

        t_122[k] = f_12 * lsp_62[k]
                   + f_3 * pc_x[k] * msp_62[k];

        t_123[k] = f_1 * mss0_20[k]
                   - f_2 * mss1_20[k]
                   + f_3 * pc_y[k] * msp_61[k];

        t_124[k] = f_3 * pc_y[k] * msp_62[k];
    }
}

static auto
compute_prim_msd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsd0,
                                                          const size_t lsp, const size_t lsd1,
                                                          const size_t mss0, const size_t mss1,
                                                          const size_t msp, const size_t ncols,
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
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsd0_90 = buffer.data(lsd0 + 90);
    const auto *lsd0_93 = buffer.data(lsd0 + 93);
    const auto *lsd0_120 = buffer.data(lsd0 + 120);
    const auto *lsd0_125 = buffer.data(lsd0 + 125);
    const auto *lsd0_126 = buffer.data(lsd0 + 126);
    const auto *lsd0_129 = buffer.data(lsd0 + 129);
    const auto *lsd0_162 = buffer.data(lsd0 + 162);
    const auto *lsd0_167 = buffer.data(lsd0 + 167);
    const auto *lsd0_168 = buffer.data(lsd0 + 168);
    const auto *lsd0_216 = buffer.data(lsd0 + 216);
    const auto *lsd0_219 = buffer.data(lsd0 + 219);
    const auto *lsd0_221 = buffer.data(lsd0 + 221);
    const auto *lsd0_225 = buffer.data(lsd0 + 225);
    const auto *lsd0_227 = buffer.data(lsd0 + 227);
    const auto *lsd0_228 = buffer.data(lsd0 + 228);
    const auto *lsd0_231 = buffer.data(lsd0 + 231);
    const auto *lsd0_233 = buffer.data(lsd0 + 233);
    const auto *lsd0_234 = buffer.data(lsd0 + 234);
    const auto *lsd0_237 = buffer.data(lsd0 + 237);
    const auto *lsd0_239 = buffer.data(lsd0 + 239);
    const auto *lsd0_240 = buffer.data(lsd0 + 240);
    const auto *lsd0_243 = buffer.data(lsd0 + 243);

    const auto *lsp_44 = buffer.data(lsp + 44);
    const auto *lsp_46 = buffer.data(lsp + 46);
    const auto *lsp_47 = buffer.data(lsp + 47);
    const auto *lsp_50 = buffer.data(lsp + 50);
    const auto *lsp_52 = buffer.data(lsp + 52);
    const auto *lsp_53 = buffer.data(lsp + 53);
    const auto *lsp_55 = buffer.data(lsp + 55);
    const auto *lsp_56 = buffer.data(lsp + 56);
    const auto *lsp_58 = buffer.data(lsp + 58);
    const auto *lsp_59 = buffer.data(lsp + 59);
    const auto *lsp_61 = buffer.data(lsp + 61);
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
    const auto *lsp_90 = buffer.data(lsp + 90);
    const auto *lsp_91 = buffer.data(lsp + 91);
    const auto *lsp_92 = buffer.data(lsp + 92);
    const auto *lsp_93 = buffer.data(lsp + 93);
    const auto *lsp_94 = buffer.data(lsp + 94);
    const auto *lsp_95 = buffer.data(lsp + 95);
    const auto *lsp_96 = buffer.data(lsp + 96);
    const auto *lsp_97 = buffer.data(lsp + 97);
    const auto *lsp_98 = buffer.data(lsp + 98);
    const auto *lsp_99 = buffer.data(lsp + 99);
    const auto *lsp_100 = buffer.data(lsp + 100);
    const auto *lsp_101 = buffer.data(lsp + 101);
    const auto *lsp_103 = buffer.data(lsp + 103);
    const auto *lsp_104 = buffer.data(lsp + 104);
    const auto *lsp_105 = buffer.data(lsp + 105);
    const auto *lsp_107 = buffer.data(lsp + 107);
    const auto *lsp_108 = buffer.data(lsp + 108);
    const auto *lsp_109 = buffer.data(lsp + 109);
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

    const auto *lsd1_90 = buffer.data(lsd1 + 90);
    const auto *lsd1_93 = buffer.data(lsd1 + 93);
    const auto *lsd1_120 = buffer.data(lsd1 + 120);
    const auto *lsd1_125 = buffer.data(lsd1 + 125);
    const auto *lsd1_126 = buffer.data(lsd1 + 126);
    const auto *lsd1_129 = buffer.data(lsd1 + 129);
    const auto *lsd1_162 = buffer.data(lsd1 + 162);
    const auto *lsd1_167 = buffer.data(lsd1 + 167);
    const auto *lsd1_168 = buffer.data(lsd1 + 168);
    const auto *lsd1_216 = buffer.data(lsd1 + 216);
    const auto *lsd1_219 = buffer.data(lsd1 + 219);
    const auto *lsd1_221 = buffer.data(lsd1 + 221);
    const auto *lsd1_225 = buffer.data(lsd1 + 225);
    const auto *lsd1_227 = buffer.data(lsd1 + 227);
    const auto *lsd1_228 = buffer.data(lsd1 + 228);
    const auto *lsd1_231 = buffer.data(lsd1 + 231);
    const auto *lsd1_233 = buffer.data(lsd1 + 233);
    const auto *lsd1_234 = buffer.data(lsd1 + 234);
    const auto *lsd1_237 = buffer.data(lsd1 + 237);
    const auto *lsd1_239 = buffer.data(lsd1 + 239);
    const auto *lsd1_240 = buffer.data(lsd1 + 240);
    const auto *lsd1_243 = buffer.data(lsd1 + 243);

    const auto *mss0_20 = buffer.data(mss0 + 20);
    const auto *mss0_21 = buffer.data(mss0 + 21);
    const auto *mss0_22 = buffer.data(mss0 + 22);
    const auto *mss0_23 = buffer.data(mss0 + 23);
    const auto *mss0_24 = buffer.data(mss0 + 24);
    const auto *mss0_25 = buffer.data(mss0 + 25);
    const auto *mss0_26 = buffer.data(mss0 + 26);
    const auto *mss0_27 = buffer.data(mss0 + 27);
    const auto *mss0_28 = buffer.data(mss0 + 28);
    const auto *mss0_29 = buffer.data(mss0 + 29);
    const auto *mss0_30 = buffer.data(mss0 + 30);
    const auto *mss0_31 = buffer.data(mss0 + 31);
    const auto *mss0_32 = buffer.data(mss0 + 32);
    const auto *mss0_33 = buffer.data(mss0 + 33);
    const auto *mss0_34 = buffer.data(mss0 + 34);
    const auto *mss0_35 = buffer.data(mss0 + 35);

    const auto *mss1_20 = buffer.data(mss1 + 20);
    const auto *mss1_21 = buffer.data(mss1 + 21);
    const auto *mss1_22 = buffer.data(mss1 + 22);
    const auto *mss1_23 = buffer.data(mss1 + 23);
    const auto *mss1_24 = buffer.data(mss1 + 24);
    const auto *mss1_25 = buffer.data(mss1 + 25);
    const auto *mss1_26 = buffer.data(mss1 + 26);
    const auto *mss1_27 = buffer.data(mss1 + 27);
    const auto *mss1_28 = buffer.data(mss1 + 28);
    const auto *mss1_29 = buffer.data(mss1 + 29);
    const auto *mss1_30 = buffer.data(mss1 + 30);
    const auto *mss1_31 = buffer.data(mss1 + 31);
    const auto *mss1_32 = buffer.data(mss1 + 32);
    const auto *mss1_33 = buffer.data(mss1 + 33);
    const auto *mss1_34 = buffer.data(mss1 + 34);
    const auto *mss1_35 = buffer.data(mss1 + 35);

    const auto *msp_62 = buffer.data(msp + 62);
    const auto *msp_63 = buffer.data(msp + 63);
    const auto *msp_64 = buffer.data(msp + 64);
    const auto *msp_65 = buffer.data(msp + 65);
    const auto *msp_67 = buffer.data(msp + 67);
    const auto *msp_68 = buffer.data(msp + 68);
    const auto *msp_69 = buffer.data(msp + 69);
    const auto *msp_70 = buffer.data(msp + 70);
    const auto *msp_71 = buffer.data(msp + 71);
    const auto *msp_72 = buffer.data(msp + 72);
    const auto *msp_73 = buffer.data(msp + 73);
    const auto *msp_74 = buffer.data(msp + 74);
    const auto *msp_75 = buffer.data(msp + 75);
    const auto *msp_76 = buffer.data(msp + 76);
    const auto *msp_77 = buffer.data(msp + 77);
    const auto *msp_79 = buffer.data(msp + 79);
    const auto *msp_80 = buffer.data(msp + 80);
    const auto *msp_81 = buffer.data(msp + 81);
    const auto *msp_82 = buffer.data(msp + 82);
    const auto *msp_83 = buffer.data(msp + 83);
    const auto *msp_84 = buffer.data(msp + 84);
    const auto *msp_85 = buffer.data(msp + 85);
    const auto *msp_86 = buffer.data(msp + 86);
    const auto *msp_88 = buffer.data(msp + 88);
    const auto *msp_89 = buffer.data(msp + 89);
    const auto *msp_90 = buffer.data(msp + 90);
    const auto *msp_91 = buffer.data(msp + 91);
    const auto *msp_92 = buffer.data(msp + 92);
    const auto *msp_93 = buffer.data(msp + 93);
    const auto *msp_94 = buffer.data(msp + 94);
    const auto *msp_95 = buffer.data(msp + 95);
    const auto *msp_96 = buffer.data(msp + 96);
    const auto *msp_97 = buffer.data(msp + 97);
    const auto *msp_98 = buffer.data(msp + 98);
    const auto *msp_99 = buffer.data(msp + 99);
    const auto *msp_100 = buffer.data(msp + 100);
    const auto *msp_101 = buffer.data(msp + 101);
    const auto *msp_103 = buffer.data(msp + 103);
    const auto *msp_104 = buffer.data(msp + 104);
    const auto *msp_105 = buffer.data(msp + 105);
    const auto *msp_106 = buffer.data(msp + 106);
    const auto *msp_107 = buffer.data(msp + 107);
    const auto *msp_108 = buffer.data(msp + 108);
    const auto *msp_109 = buffer.data(msp + 109);
    const auto *msp_112 = buffer.data(msp + 112);
    const auto *msp_113 = buffer.data(msp + 113);
    const auto *msp_115 = buffer.data(msp + 115);
    const auto *msp_116 = buffer.data(msp + 116);
    const auto *msp_118 = buffer.data(msp + 118);
    const auto *msp_119 = buffer.data(msp + 119);
    const auto *msp_121 = buffer.data(msp + 121);
    const auto *msp_122 = buffer.data(msp + 122);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_z, lsp_44, lsp_63, lsp_64, \
                         mss0_20, mss0_21, mss1_20, mss1_21, msp_62, msp_63, \
                         msp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_11 * lsp_44[k]
                   + f_1 * mss0_20[k]
                   - f_2 * mss1_20[k]
                   + f_3 * pc_z[k] * msp_62[k];

        t_126[k] = f_10 * lsp_63[k]
                   + f_1 * mss0_21[k]
                   - f_2 * mss1_21[k]
                   + f_3 * pc_x[k] * msp_63[k];

        t_127[k] = f_10 * lsp_64[k]
                   + f_3 * pc_x[k] * msp_64[k];

        t_128[k] = f_3 * pc_z[k] * msp_63[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_z, pc_y, pc_z, lsd0_90, lsp_46, \
                         lsd1_90, mss0_21, mss1_21, msp_64, msp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_9 * lsp_46[k]
                   + f_1 * mss0_21[k]
                   - f_2 * mss1_21[k]
                   + f_3 * pc_y[k] * msp_64[k];

        t_130[k] = f_3 * pc_z[k] * msp_64[k];

        t_131[k] = f_1 * mss0_21[k]
                   - f_2 * mss1_21[k]
                   + f_3 * pc_z[k] * msp_65[k];

        t_132[k] = pa_z[k] * lsd0_90[k]
                   - f_4 * pc_z[k] * lsd1_90[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, lsd0_93, lsp_50, \
                         lsp_67, lsp_68, lsd1_93, msp_67, msp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_10 * lsp_67[k]
                   + f_3 * pc_x[k] * msp_67[k];

        t_134[k] = f_10 * lsp_68[k]
                   + f_3 * pc_x[k] * msp_68[k];

        t_135[k] = pa_z[k] * lsd0_93[k]
                   - f_4 * pc_z[k] * lsd1_93[k];

        t_136[k] = f_11 * lsp_50[k]
                   + f_3 * pc_y[k] * msp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_x, pc_z, lsp_47, lsp_69, lsp_70, mss0_22, \
                         mss0_23, mss1_22, mss1_23, msp_68, msp_69, \
                         msp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * lsp_47[k]
                   + f_1 * mss0_22[k]
                   - f_2 * mss1_22[k]
                   + f_3 * pc_z[k] * msp_68[k];

        t_138[k] = f_10 * lsp_69[k]
                   + f_1 * mss0_23[k]
                   - f_2 * mss1_23[k]
                   + f_3 * pc_x[k] * msp_69[k];

        t_139[k] = f_10 * lsp_70[k]
                   + f_3 * pc_x[k] * msp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, pc_z, lsp_50, lsp_52, lsp_53, \
                         lsp_71, mss0_23, mss1_23, msp_70, msp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_10 * lsp_71[k]
                   + f_3 * pc_x[k] * msp_71[k];

        t_141[k] = f_12 * lsp_52[k]
                   + f_1 * mss0_23[k]
                   - f_2 * mss1_23[k]
                   + f_3 * pc_y[k] * msp_70[k];

        t_142[k] = f_12 * lsp_53[k]
                   + f_3 * pc_y[k] * msp_71[k];

        t_143[k] = f_8 * lsp_50[k]
                   + f_1 * mss0_23[k]
                   - f_2 * mss1_23[k]
                   + f_3 * pc_z[k] * msp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pc_x, pc_y, lsp_55, lsp_72, lsp_73, \
                         lsp_74, mss0_24, mss1_24, msp_72, msp_73, \
                         msp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_10 * lsp_72[k]
                   + f_1 * mss0_24[k]
                   - f_2 * mss1_24[k]
                   + f_3 * pc_x[k] * msp_72[k];

        t_145[k] = f_10 * lsp_73[k]
                   + f_3 * pc_x[k] * msp_73[k];

        t_146[k] = f_10 * lsp_74[k]
                   + f_3 * pc_x[k] * msp_74[k];

        t_147[k] = f_10 * lsp_55[k]
                   + f_1 * mss0_24[k]
                   - f_2 * mss1_24[k]
                   + f_3 * pc_y[k] * msp_73[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pc_x, pc_y, pc_z, lsp_53, lsp_56, lsp_75, \
                         mss0_24, mss0_25, mss1_24, mss1_25, msp_74, \
                         msp_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * lsp_56[k]
                   + f_3 * pc_y[k] * msp_74[k];

        t_149[k] = f_10 * lsp_53[k]
                   + f_1 * mss0_24[k]
                   - f_2 * mss1_24[k]
                   + f_3 * pc_z[k] * msp_74[k];

        t_150[k] = f_10 * lsp_75[k]
                   + f_1 * mss0_25[k]
                   - f_2 * mss1_25[k]
                   + f_3 * pc_x[k] * msp_75[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, lsp_58, lsp_59, lsp_76, \
                         lsp_77, mss0_25, mss1_25, msp_76, msp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_10 * lsp_76[k]
                   + f_3 * pc_x[k] * msp_76[k];

        t_152[k] = f_10 * lsp_77[k]
                   + f_3 * pc_x[k] * msp_77[k];

        t_153[k] = f_8 * lsp_58[k]
                   + f_1 * mss0_25[k]
                   - f_2 * mss1_25[k]
                   + f_3 * pc_y[k] * msp_76[k];

        t_154[k] = f_8 * lsp_59[k]
                   + f_3 * pc_y[k] * msp_77[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, pa_y, pc_x, pc_y, pc_z, lsd0_120, lsp_56, \
                         lsp_79, lsd1_120, mss0_25, mss1_25, msp_77, \
                         msp_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * lsp_56[k]
                   + f_1 * mss0_25[k]
                   - f_2 * mss1_25[k]
                   + f_3 * pc_z[k] * msp_77[k];

        t_156[k] = pa_y[k] * lsd0_120[k]
                   - f_4 * pc_y[k] * lsd1_120[k];

        t_157[k] = f_10 * lsp_79[k]
                   + f_3 * pc_x[k] * msp_79[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, lsd0_125, lsp_61, \
                         lsp_62, lsp_80, lsd1_125, mss0_26, mss1_26, msp_79, \
                         msp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_10 * lsp_80[k]
                   + f_3 * pc_x[k] * msp_80[k];

        t_159[k] = f_6 * lsp_61[k]
                   + f_1 * mss0_26[k]
                   - f_2 * mss1_26[k]
                   + f_3 * pc_y[k] * msp_79[k];

        t_160[k] = f_6 * lsp_62[k]
                   + f_3 * pc_y[k] * msp_80[k];

        t_161[k] = pa_y[k] * lsd0_125[k]
                   - f_4 * pc_y[k] * lsd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pc_x, pc_y, lsp_81, lsp_83, \
                         mss0_27, mss1_27, msp_81, msp_82, msp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_10 * lsp_81[k]
                   + f_1 * mss0_27[k]
                   - f_2 * mss1_27[k]
                   + f_3 * pc_x[k] * msp_81[k];

        t_163[k] = f_3 * pc_y[k] * msp_81[k];

        t_164[k] = f_10 * lsp_83[k]
                   + f_3 * pc_x[k] * msp_83[k];

        t_165[k] = f_1 * mss0_27[k]
                   - f_2 * mss1_27[k]
                   + f_3 * pc_y[k] * msp_82[k];

        t_166[k] = f_3 * pc_y[k] * msp_83[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pc_x, pc_z, lsp_62, lsp_84, lsp_85, \
                         mss0_27, mss0_28, mss1_27, mss1_28, msp_83, msp_84, \
                         msp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_9 * lsp_62[k]
                   + f_1 * mss0_27[k]
                   - f_2 * mss1_27[k]
                   + f_3 * pc_z[k] * msp_83[k];

        t_168[k] = f_8 * lsp_84[k]
                   + f_1 * mss0_28[k]
                   - f_2 * mss1_28[k]
                   + f_3 * pc_x[k] * msp_84[k];

        t_169[k] = f_8 * lsp_85[k]
                   + f_3 * pc_x[k] * msp_85[k];

        t_170[k] = f_3 * pc_z[k] * msp_84[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pa_z, pc_y, pc_z, lsd0_126, lsp_64, \
                         lsd1_126, mss0_28, mss1_28, msp_85, msp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_7 * lsp_64[k]
                   + f_1 * mss0_28[k]
                   - f_2 * mss1_28[k]
                   + f_3 * pc_y[k] * msp_85[k];

        t_172[k] = f_3 * pc_z[k] * msp_85[k];

        t_173[k] = f_1 * mss0_28[k]
                   - f_2 * mss1_28[k]
                   + f_3 * pc_z[k] * msp_86[k];

        t_174[k] = pa_z[k] * lsd0_126[k]
                   - f_4 * pc_z[k] * lsd1_126[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pc_x, pc_y, pc_z, lsd0_129, lsp_68, \
                         lsp_88, lsp_89, lsd1_129, msp_88, msp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_8 * lsp_88[k]
                   + f_3 * pc_x[k] * msp_88[k];

        t_176[k] = f_8 * lsp_89[k]
                   + f_3 * pc_x[k] * msp_89[k];

        t_177[k] = pa_z[k] * lsd0_129[k]
                   - f_4 * pc_z[k] * lsd1_129[k];

        t_178[k] = f_9 * lsp_68[k]
                   + f_3 * pc_y[k] * msp_89[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, lsp_65, lsp_90, lsp_91, mss0_29, \
                         mss0_30, mss1_29, mss1_30, msp_89, msp_90, \
                         msp_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_6 * lsp_65[k]
                   + f_1 * mss0_29[k]
                   - f_2 * mss1_29[k]
                   + f_3 * pc_z[k] * msp_89[k];

        t_180[k] = f_8 * lsp_90[k]
                   + f_1 * mss0_30[k]
                   - f_2 * mss1_30[k]
                   + f_3 * pc_x[k] * msp_90[k];

        t_181[k] = f_8 * lsp_91[k]
                   + f_3 * pc_x[k] * msp_91[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_y, pc_z, lsp_68, lsp_70, lsp_71, \
                         lsp_92, mss0_30, mss1_30, msp_91, msp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_8 * lsp_92[k]
                   + f_3 * pc_x[k] * msp_92[k];

        t_183[k] = f_11 * lsp_70[k]
                   + f_1 * mss0_30[k]
                   - f_2 * mss1_30[k]
                   + f_3 * pc_y[k] * msp_91[k];

        t_184[k] = f_11 * lsp_71[k]
                   + f_3 * pc_y[k] * msp_92[k];

        t_185[k] = f_8 * lsp_68[k]
                   + f_1 * mss0_30[k]
                   - f_2 * mss1_30[k]
                   + f_3 * pc_z[k] * msp_92[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, pc_y, lsp_73, lsp_93, lsp_94, \
                         lsp_95, mss0_31, mss1_31, msp_93, msp_94, \
                         msp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_8 * lsp_93[k]
                   + f_1 * mss0_31[k]
                   - f_2 * mss1_31[k]
                   + f_3 * pc_x[k] * msp_93[k];

        t_187[k] = f_8 * lsp_94[k]
                   + f_3 * pc_x[k] * msp_94[k];

        t_188[k] = f_8 * lsp_95[k]
                   + f_3 * pc_x[k] * msp_95[k];

        t_189[k] = f_12 * lsp_73[k]
                   + f_1 * mss0_31[k]
                   - f_2 * mss1_31[k]
                   + f_3 * pc_y[k] * msp_94[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_x, pc_y, pc_z, lsp_71, lsp_74, lsp_96, \
                         mss0_31, mss0_32, mss1_31, mss1_32, msp_95, \
                         msp_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_12 * lsp_74[k]
                   + f_3 * pc_y[k] * msp_95[k];

        t_191[k] = f_10 * lsp_71[k]
                   + f_1 * mss0_31[k]
                   - f_2 * mss1_31[k]
                   + f_3 * pc_z[k] * msp_95[k];

        t_192[k] = f_8 * lsp_96[k]
                   + f_1 * mss0_32[k]
                   - f_2 * mss1_32[k]
                   + f_3 * pc_x[k] * msp_96[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, lsp_76, lsp_77, lsp_97, \
                         lsp_98, mss0_32, mss1_32, msp_97, msp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * lsp_97[k]
                   + f_3 * pc_x[k] * msp_97[k];

        t_194[k] = f_8 * lsp_98[k]
                   + f_3 * pc_x[k] * msp_98[k];

        t_195[k] = f_10 * lsp_76[k]
                   + f_1 * mss0_32[k]
                   - f_2 * mss1_32[k]
                   + f_3 * pc_y[k] * msp_97[k];

        t_196[k] = f_10 * lsp_77[k]
                   + f_3 * pc_y[k] * msp_98[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pc_x, pc_z, lsp_74, lsp_99, lsp_100, mss0_32, \
                         mss0_33, mss1_32, mss1_33, msp_98, msp_99, \
                         msp_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * lsp_74[k]
                   + f_1 * mss0_32[k]
                   - f_2 * mss1_32[k]
                   + f_3 * pc_z[k] * msp_98[k];

        t_198[k] = f_8 * lsp_99[k]
                   + f_1 * mss0_33[k]
                   - f_2 * mss1_33[k]
                   + f_3 * pc_x[k] * msp_99[k];

        t_199[k] = f_8 * lsp_100[k]
                   + f_3 * pc_x[k] * msp_100[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, pc_y, pc_z, lsp_77, lsp_79, lsp_80, \
                         lsp_101, mss0_33, mss1_33, msp_100, msp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_8 * lsp_101[k]
                   + f_3 * pc_x[k] * msp_101[k];

        t_201[k] = f_8 * lsp_79[k]
                   + f_1 * mss0_33[k]
                   - f_2 * mss1_33[k]
                   + f_3 * pc_y[k] * msp_100[k];

        t_202[k] = f_8 * lsp_80[k]
                   + f_3 * pc_y[k] * msp_101[k];

        t_203[k] = f_11 * lsp_77[k]
                   + f_1 * mss0_33[k]
                   - f_2 * mss1_33[k]
                   + f_3 * pc_z[k] * msp_101[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_y, pc_x, pc_y, lsd0_162, lsp_82, \
                         lsp_103, lsp_104, lsd1_162, mss0_34, mss1_34, msp_103, \
                         msp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * lsd0_162[k]
                   - f_4 * pc_y[k] * lsd1_162[k];

        t_205[k] = f_8 * lsp_103[k]
                   + f_3 * pc_x[k] * msp_103[k];

        t_206[k] = f_8 * lsp_104[k]
                   + f_3 * pc_x[k] * msp_104[k];

        t_207[k] = f_6 * lsp_82[k]
                   + f_1 * mss0_34[k]
                   - f_2 * mss1_34[k]
                   + f_3 * pc_y[k] * msp_103[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_y, pc_x, pc_y, lsd0_167, lsp_83, \
                         lsp_105, lsd1_167, mss0_35, mss1_35, msp_104, \
                         msp_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_6 * lsp_83[k]
                   + f_3 * pc_y[k] * msp_104[k];

        t_209[k] = pa_y[k] * lsd0_167[k]
                   - f_4 * pc_y[k] * lsd1_167[k];

        t_210[k] = f_8 * lsp_105[k]
                   + f_1 * mss0_35[k]
                   - f_2 * mss1_35[k]
                   + f_3 * pc_x[k] * msp_105[k];

        t_211[k] = f_3 * pc_y[k] * msp_105[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, lsp_83, lsp_107, \
                         mss0_35, mss1_35, msp_106, msp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_8 * lsp_107[k]
                   + f_3 * pc_x[k] * msp_107[k];

        t_213[k] = f_1 * mss0_35[k]
                   - f_2 * mss1_35[k]
                   + f_3 * pc_y[k] * msp_106[k];

        t_214[k] = f_3 * pc_y[k] * msp_107[k];

        t_215[k] = f_7 * lsp_83[k]
                   + f_1 * mss0_35[k]
                   - f_2 * mss1_35[k]
                   + f_3 * pc_z[k] * msp_107[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, pa_x, pc_x, pc_z, lsd0_216, \
                         lsd0_219, lsp_108, lsp_109, lsd1_216, lsd1_219, msp_108, \
                         msp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pa_x[k] * lsd0_216[k]
                   + f_8 * lsp_108[k]
                   - f_4 * pc_x[k] * lsd1_216[k];

        t_217[k] = f_6 * lsp_109[k]
                   + f_3 * pc_x[k] * msp_109[k];

        t_218[k] = f_3 * pc_z[k] * msp_108[k];

        t_219[k] = pa_x[k] * lsd0_219[k]
                   - f_4 * pc_x[k] * lsd1_219[k];

        t_220[k] = f_3 * pc_z[k] * msp_109[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_x, pa_z, pc_x, pc_z, lsd0_168, \
                         lsd0_221, lsp_112, lsp_113, lsd1_168, lsd1_221, msp_112, \
                         msp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_x[k] * lsd0_221[k]
                   - f_4 * pc_x[k] * lsd1_221[k];

        t_222[k] = pa_z[k] * lsd0_168[k]
                   - f_4 * pc_z[k] * lsd1_168[k];

        t_223[k] = f_6 * lsp_112[k]
                   + f_3 * pc_x[k] * msp_112[k];

        t_224[k] = f_6 * lsp_113[k]
                   + f_3 * pc_x[k] * msp_113[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_x, pc_x, pc_y, lsd0_225, lsd0_227, \
                         lsd0_228, lsp_89, lsp_114, lsd1_225, lsd1_227, lsd1_228, \
                         msp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = pa_x[k] * lsd0_225[k]
                   - f_4 * pc_x[k] * lsd1_225[k];

        t_226[k] = f_7 * lsp_89[k]
                   + f_3 * pc_y[k] * msp_113[k];

        t_227[k] = pa_x[k] * lsd0_227[k]
                   - f_4 * pc_x[k] * lsd1_227[k];

        t_228[k] = pa_x[k] * lsd0_228[k]
                   + f_8 * lsp_114[k]
                   - f_4 * pc_x[k] * lsd1_228[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_x, pc_x, pc_y, lsd0_231, lsp_92, \
                         lsp_115, lsp_116, lsd1_231, msp_115, msp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_6 * lsp_115[k]
                   + f_3 * pc_x[k] * msp_115[k];

        t_230[k] = f_6 * lsp_116[k]
                   + f_3 * pc_x[k] * msp_116[k];

        t_231[k] = pa_x[k] * lsd0_231[k]
                   - f_4 * pc_x[k] * lsd1_231[k];

        t_232[k] = f_9 * lsp_92[k]
                   + f_3 * pc_y[k] * msp_116[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_x, pc_x, lsd0_233, lsd0_234, lsp_117, \
                         lsp_118, lsp_119, lsd1_233, lsd1_234, msp_118, \
                         msp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_x[k] * lsd0_233[k]
                   - f_4 * pc_x[k] * lsd1_233[k];

        t_234[k] = pa_x[k] * lsd0_234[k]
                   + f_8 * lsp_117[k]
                   - f_4 * pc_x[k] * lsd1_234[k];

        t_235[k] = f_6 * lsp_118[k]
                   + f_3 * pc_x[k] * msp_118[k];

        t_236[k] = f_6 * lsp_119[k]
                   + f_3 * pc_x[k] * msp_119[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pa_x, pc_x, pc_y, lsd0_237, lsd0_239, \
                         lsd0_240, lsp_95, lsp_120, lsd1_237, lsd1_239, lsd1_240, \
                         msp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pa_x[k] * lsd0_237[k]
                   - f_4 * pc_x[k] * lsd1_237[k];

        t_238[k] = f_11 * lsp_95[k]
                   + f_3 * pc_y[k] * msp_119[k];

        t_239[k] = pa_x[k] * lsd0_239[k]
                   - f_4 * pc_x[k] * lsd1_239[k];

        t_240[k] = pa_x[k] * lsd0_240[k]
                   + f_8 * lsp_120[k]
                   - f_4 * pc_x[k] * lsd1_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pa_x, pc_x, pc_y, lsd0_243, lsp_98, \
                         lsp_121, lsp_122, lsd1_243, msp_121, msp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_6 * lsp_121[k]
                   + f_3 * pc_x[k] * msp_121[k];

        t_242[k] = f_6 * lsp_122[k]
                   + f_3 * pc_x[k] * msp_122[k];

        t_243[k] = pa_x[k] * lsd0_243[k]
                   - f_4 * pc_x[k] * lsd1_243[k];

        t_244[k] = f_12 * lsp_98[k]
                   + f_3 * pc_y[k] * msp_122[k];
    }
}

static auto
compute_prim_msd_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t lsd0,
                                                          const size_t lsp, const size_t lsd1,
                                                          const size_t mss0, const size_t mss1,
                                                          const size_t msp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 4.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 3.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *lsd0_210 = buffer.data(lsd0 + 210);
    const auto *lsd0_216 = buffer.data(lsd0 + 216);
    const auto *lsd0_219 = buffer.data(lsd0 + 219);
    const auto *lsd0_245 = buffer.data(lsd0 + 245);
    const auto *lsd0_246 = buffer.data(lsd0 + 246);
    const auto *lsd0_249 = buffer.data(lsd0 + 249);
    const auto *lsd0_251 = buffer.data(lsd0 + 251);
    const auto *lsd0_252 = buffer.data(lsd0 + 252);
    const auto *lsd0_255 = buffer.data(lsd0 + 255);
    const auto *lsd0_257 = buffer.data(lsd0 + 257);
    const auto *lsd0_261 = buffer.data(lsd0 + 261);
    const auto *lsd0_263 = buffer.data(lsd0 + 263);
    const auto *lsd0_264 = buffer.data(lsd0 + 264);
    const auto *lsd0_267 = buffer.data(lsd0 + 267);
    const auto *lsd0_269 = buffer.data(lsd0 + 269);

    const auto *lsp_101 = buffer.data(lsp + 101);
    const auto *lsp_104 = buffer.data(lsp + 104);
    const auto *lsp_107 = buffer.data(lsp + 107);
    const auto *lsp_109 = buffer.data(lsp + 109);
    const auto *lsp_110 = buffer.data(lsp + 110);
    const auto *lsp_113 = buffer.data(lsp + 113);
    const auto *lsp_115 = buffer.data(lsp + 115);
    const auto *lsp_116 = buffer.data(lsp + 116);
    const auto *lsp_118 = buffer.data(lsp + 118);
    const auto *lsp_119 = buffer.data(lsp + 119);
    const auto *lsp_121 = buffer.data(lsp + 121);
    const auto *lsp_122 = buffer.data(lsp + 122);
    const auto *lsp_123 = buffer.data(lsp + 123);
    const auto *lsp_124 = buffer.data(lsp + 124);
    const auto *lsp_125 = buffer.data(lsp + 125);
    const auto *lsp_126 = buffer.data(lsp + 126);
    const auto *lsp_127 = buffer.data(lsp + 127);
    const auto *lsp_128 = buffer.data(lsp + 128);
    const auto *lsp_130 = buffer.data(lsp + 130);
    const auto *lsp_131 = buffer.data(lsp + 131);
    const auto *lsp_132 = buffer.data(lsp + 132);
    const auto *lsp_133 = buffer.data(lsp + 133);
    const auto *lsp_134 = buffer.data(lsp + 134);

    const auto *lsd1_210 = buffer.data(lsd1 + 210);
    const auto *lsd1_216 = buffer.data(lsd1 + 216);
    const auto *lsd1_219 = buffer.data(lsd1 + 219);
    const auto *lsd1_245 = buffer.data(lsd1 + 245);
    const auto *lsd1_246 = buffer.data(lsd1 + 246);
    const auto *lsd1_249 = buffer.data(lsd1 + 249);
    const auto *lsd1_251 = buffer.data(lsd1 + 251);
    const auto *lsd1_252 = buffer.data(lsd1 + 252);
    const auto *lsd1_255 = buffer.data(lsd1 + 255);
    const auto *lsd1_257 = buffer.data(lsd1 + 257);
    const auto *lsd1_261 = buffer.data(lsd1 + 261);
    const auto *lsd1_263 = buffer.data(lsd1 + 263);
    const auto *lsd1_264 = buffer.data(lsd1 + 264);
    const auto *lsd1_267 = buffer.data(lsd1 + 267);
    const auto *lsd1_269 = buffer.data(lsd1 + 269);

    const auto *mss0_45 = buffer.data(mss0 + 45);
    const auto *mss0_46 = buffer.data(mss0 + 46);
    const auto *mss0_47 = buffer.data(mss0 + 47);
    const auto *mss0_48 = buffer.data(mss0 + 48);
    const auto *mss0_49 = buffer.data(mss0 + 49);
    const auto *mss0_50 = buffer.data(mss0 + 50);
    const auto *mss0_51 = buffer.data(mss0 + 51);
    const auto *mss0_52 = buffer.data(mss0 + 52);
    const auto *mss0_54 = buffer.data(mss0 + 54);

    const auto *mss1_45 = buffer.data(mss1 + 45);
    const auto *mss1_46 = buffer.data(mss1 + 46);
    const auto *mss1_47 = buffer.data(mss1 + 47);
    const auto *mss1_48 = buffer.data(mss1 + 48);
    const auto *mss1_49 = buffer.data(mss1 + 49);
    const auto *mss1_50 = buffer.data(mss1 + 50);
    const auto *mss1_51 = buffer.data(mss1 + 51);
    const auto *mss1_52 = buffer.data(mss1 + 52);
    const auto *mss1_54 = buffer.data(mss1 + 54);

    const auto *msp_124 = buffer.data(msp + 124);
    const auto *msp_125 = buffer.data(msp + 125);
    const auto *msp_127 = buffer.data(msp + 127);
    const auto *msp_128 = buffer.data(msp + 128);
    const auto *msp_130 = buffer.data(msp + 130);
    const auto *msp_131 = buffer.data(msp + 131);
    const auto *msp_132 = buffer.data(msp + 132);
    const auto *msp_134 = buffer.data(msp + 134);
    const auto *msp_135 = buffer.data(msp + 135);
    const auto *msp_136 = buffer.data(msp + 136);
    const auto *msp_137 = buffer.data(msp + 137);
    const auto *msp_139 = buffer.data(msp + 139);
    const auto *msp_140 = buffer.data(msp + 140);
    const auto *msp_141 = buffer.data(msp + 141);
    const auto *msp_142 = buffer.data(msp + 142);
    const auto *msp_143 = buffer.data(msp + 143);
    const auto *msp_144 = buffer.data(msp + 144);
    const auto *msp_145 = buffer.data(msp + 145);
    const auto *msp_146 = buffer.data(msp + 146);
    const auto *msp_147 = buffer.data(msp + 147);
    const auto *msp_148 = buffer.data(msp + 148);
    const auto *msp_149 = buffer.data(msp + 149);
    const auto *msp_150 = buffer.data(msp + 150);
    const auto *msp_151 = buffer.data(msp + 151);
    const auto *msp_152 = buffer.data(msp + 152);
    const auto *msp_153 = buffer.data(msp + 153);
    const auto *msp_154 = buffer.data(msp + 154);
    const auto *msp_155 = buffer.data(msp + 155);
    const auto *msp_156 = buffer.data(msp + 156);
    const auto *msp_157 = buffer.data(msp + 157);
    const auto *msp_158 = buffer.data(msp + 158);
    const auto *msp_160 = buffer.data(msp + 160);
    const auto *msp_161 = buffer.data(msp + 161);
    const auto *msp_162 = buffer.data(msp + 162);
    const auto *msp_163 = buffer.data(msp + 163);
    const auto *msp_164 = buffer.data(msp + 164);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pa_x, pc_x, lsd0_245, lsd0_246, lsp_123, \
                         lsp_124, lsp_125, lsd1_245, lsd1_246, msp_124, \
                         msp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = pa_x[k] * lsd0_245[k]
                   - f_4 * pc_x[k] * lsd1_245[k];

        t_246[k] = pa_x[k] * lsd0_246[k]
                   + f_8 * lsp_123[k]
                   - f_4 * pc_x[k] * lsd1_246[k];

        t_247[k] = f_6 * lsp_124[k]
                   + f_3 * pc_x[k] * msp_124[k];

        t_248[k] = f_6 * lsp_125[k]
                   + f_3 * pc_x[k] * msp_125[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, pa_x, pc_x, pc_y, lsd0_249, lsd0_251, \
                         lsd0_252, lsp_101, lsp_126, lsd1_249, lsd1_251, lsd1_252, \
                         msp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_x[k] * lsd0_249[k]
                   - f_4 * pc_x[k] * lsd1_249[k];

        t_250[k] = f_10 * lsp_101[k]
                   + f_3 * pc_y[k] * msp_125[k];

        t_251[k] = pa_x[k] * lsd0_251[k]
                   - f_4 * pc_x[k] * lsd1_251[k];

        t_252[k] = pa_x[k] * lsd0_252[k]
                   + f_8 * lsp_126[k]
                   - f_4 * pc_x[k] * lsd1_252[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_x, pc_x, pc_y, lsd0_255, lsp_104, \
                         lsp_127, lsp_128, lsd1_255, msp_127, msp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_6 * lsp_127[k]
                   + f_3 * pc_x[k] * msp_127[k];

        t_254[k] = f_6 * lsp_128[k]
                   + f_3 * pc_x[k] * msp_128[k];

        t_255[k] = pa_x[k] * lsd0_255[k]
                   - f_4 * pc_x[k] * lsd1_255[k];

        t_256[k] = f_8 * lsp_104[k]
                   + f_3 * pc_y[k] * msp_128[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pa_x, pa_y, pc_x, pc_y, lsd0_210, \
                         lsd0_257, lsp_130, lsp_131, lsd1_210, lsd1_257, msp_130, \
                         msp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pa_x[k] * lsd0_257[k]
                   - f_4 * pc_x[k] * lsd1_257[k];

        t_258[k] = pa_y[k] * lsd0_210[k]
                   - f_4 * pc_y[k] * lsd1_210[k];

        t_259[k] = f_6 * lsp_130[k]
                   + f_3 * pc_x[k] * msp_130[k];

        t_260[k] = f_6 * lsp_131[k]
                   + f_3 * pc_x[k] * msp_131[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pa_x, pc_x, pc_y, lsd0_261, lsd0_263, \
                         lsd0_264, lsp_107, lsp_132, lsd1_261, lsd1_263, lsd1_264, \
                         msp_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pa_x[k] * lsd0_261[k]
                   - f_4 * pc_x[k] * lsd1_261[k];

        t_262[k] = f_6 * lsp_107[k]
                   + f_3 * pc_y[k] * msp_131[k];

        t_263[k] = pa_x[k] * lsd0_263[k]
                   - f_4 * pc_x[k] * lsd1_263[k];

        t_264[k] = pa_x[k] * lsd0_264[k]
                   + f_8 * lsp_132[k]
                   - f_4 * pc_x[k] * lsd1_264[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pa_x, pc_x, pc_y, lsd0_267, \
                         lsd0_269, lsp_134, lsd1_267, lsd1_269, msp_132, \
                         msp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_3 * pc_y[k] * msp_132[k];

        t_266[k] = f_6 * lsp_134[k]
                   + f_3 * pc_x[k] * msp_134[k];

        t_267[k] = pa_x[k] * lsd0_267[k]
                   - f_4 * pc_x[k] * lsd1_267[k];

        t_268[k] = f_3 * pc_y[k] * msp_134[k];

        t_269[k] = pa_x[k] * lsd0_269[k]
                   - f_4 * pc_x[k] * lsd1_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, lsp_109, \
                         mss0_45, mss1_45, msp_135, msp_136, msp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * mss0_45[k]
                   - f_2 * mss1_45[k]
                   + f_3 * pc_x[k] * msp_135[k];

        t_271[k] = f_3 * pc_x[k] * msp_136[k];

        t_272[k] = f_3 * pc_x[k] * msp_137[k];

        t_273[k] = f_0 * lsp_109[k]
                   + f_1 * mss0_45[k]
                   - f_2 * mss1_45[k]
                   + f_3 * pc_y[k] * msp_136[k];

        t_274[k] = f_3 * pc_z[k] * msp_136[k];

        t_275[k] = f_1 * mss0_45[k]
                   - f_2 * mss1_45[k]
                   + f_3 * pc_z[k] * msp_137[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, pa_z, pc_x, pc_y, pc_z, lsd0_216, \
                         lsd0_219, lsp_113, lsd1_216, lsd1_219, msp_139, \
                         msp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_z[k] * lsd0_216[k]
                   - f_4 * pc_z[k] * lsd1_216[k];

        t_277[k] = f_3 * pc_x[k] * msp_139[k];

        t_278[k] = f_3 * pc_x[k] * msp_140[k];

        t_279[k] = pa_z[k] * lsd0_219[k]
                   - f_4 * pc_z[k] * lsd1_219[k];

        t_280[k] = f_5 * lsp_113[k]
                   + f_3 * pc_y[k] * msp_140[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pc_x, pc_z, lsp_110, mss0_46, mss0_47, \
                         mss1_46, mss1_47, msp_140, msp_141, msp_142, \
                         msp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_6 * lsp_110[k]
                   + f_1 * mss0_46[k]
                   - f_2 * mss1_46[k]
                   + f_3 * pc_z[k] * msp_140[k];

        t_282[k] = f_1 * mss0_47[k]
                   - f_2 * mss1_47[k]
                   + f_3 * pc_x[k] * msp_141[k];

        t_283[k] = f_3 * pc_x[k] * msp_142[k];

        t_284[k] = f_3 * pc_x[k] * msp_143[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pc_y, pc_z, lsp_113, lsp_115, lsp_116, mss0_47, \
                         mss1_47, msp_142, msp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_7 * lsp_115[k]
                   + f_1 * mss0_47[k]
                   - f_2 * mss1_47[k]
                   + f_3 * pc_y[k] * msp_142[k];

        t_286[k] = f_7 * lsp_116[k]
                   + f_3 * pc_y[k] * msp_143[k];

        t_287[k] = f_8 * lsp_113[k]
                   + f_1 * mss0_47[k]
                   - f_2 * mss1_47[k]
                   + f_3 * pc_z[k] * msp_143[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pc_x, pc_y, lsp_118, lsp_119, \
                         mss0_48, mss1_48, msp_144, msp_145, msp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_1 * mss0_48[k]
                   - f_2 * mss1_48[k]
                   + f_3 * pc_x[k] * msp_144[k];

        t_289[k] = f_3 * pc_x[k] * msp_145[k];

        t_290[k] = f_3 * pc_x[k] * msp_146[k];

        t_291[k] = f_9 * lsp_118[k]
                   + f_1 * mss0_48[k]
                   - f_2 * mss1_48[k]
                   + f_3 * pc_y[k] * msp_145[k];

        t_292[k] = f_9 * lsp_119[k]
                   + f_3 * pc_y[k] * msp_146[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pc_x, pc_z, lsp_116, mss0_48, mss0_49, \
                         mss1_48, mss1_49, msp_146, msp_147, msp_148, \
                         msp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_10 * lsp_116[k]
                   + f_1 * mss0_48[k]
                   - f_2 * mss1_48[k]
                   + f_3 * pc_z[k] * msp_146[k];

        t_294[k] = f_1 * mss0_49[k]
                   - f_2 * mss1_49[k]
                   + f_3 * pc_x[k] * msp_147[k];

        t_295[k] = f_3 * pc_x[k] * msp_148[k];

        t_296[k] = f_3 * pc_x[k] * msp_149[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, lsp_119, lsp_121, lsp_122, mss0_49, \
                         mss1_49, msp_148, msp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_11 * lsp_121[k]
                   + f_1 * mss0_49[k]
                   - f_2 * mss1_49[k]
                   + f_3 * pc_y[k] * msp_148[k];

        t_298[k] = f_11 * lsp_122[k]
                   + f_3 * pc_y[k] * msp_149[k];

        t_299[k] = f_12 * lsp_119[k]
                   + f_1 * mss0_49[k]
                   - f_2 * mss1_49[k]
                   + f_3 * pc_z[k] * msp_149[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, pc_y, lsp_124, lsp_125, \
                         mss0_50, mss1_50, msp_150, msp_151, msp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * mss0_50[k]
                   - f_2 * mss1_50[k]
                   + f_3 * pc_x[k] * msp_150[k];

        t_301[k] = f_3 * pc_x[k] * msp_151[k];

        t_302[k] = f_3 * pc_x[k] * msp_152[k];

        t_303[k] = f_12 * lsp_124[k]
                   + f_1 * mss0_50[k]
                   - f_2 * mss1_50[k]
                   + f_3 * pc_y[k] * msp_151[k];

        t_304[k] = f_12 * lsp_125[k]
                   + f_3 * pc_y[k] * msp_152[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_z, lsp_122, mss0_50, mss0_51, \
                         mss1_50, mss1_51, msp_152, msp_153, msp_154, \
                         msp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_11 * lsp_122[k]
                   + f_1 * mss0_50[k]
                   - f_2 * mss1_50[k]
                   + f_3 * pc_z[k] * msp_152[k];

        t_306[k] = f_1 * mss0_51[k]
                   - f_2 * mss1_51[k]
                   + f_3 * pc_x[k] * msp_153[k];

        t_307[k] = f_3 * pc_x[k] * msp_154[k];

        t_308[k] = f_3 * pc_x[k] * msp_155[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_y, pc_z, lsp_125, lsp_127, lsp_128, mss0_51, \
                         mss1_51, msp_154, msp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_10 * lsp_127[k]
                   + f_1 * mss0_51[k]
                   - f_2 * mss1_51[k]
                   + f_3 * pc_y[k] * msp_154[k];

        t_310[k] = f_10 * lsp_128[k]
                   + f_3 * pc_y[k] * msp_155[k];

        t_311[k] = f_9 * lsp_125[k]
                   + f_1 * mss0_51[k]
                   - f_2 * mss1_51[k]
                   + f_3 * pc_z[k] * msp_155[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, t_316, pc_x, pc_y, lsp_130, lsp_131, \
                         mss0_52, mss1_52, msp_156, msp_157, msp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_1 * mss0_52[k]
                   - f_2 * mss1_52[k]
                   + f_3 * pc_x[k] * msp_156[k];

        t_313[k] = f_3 * pc_x[k] * msp_157[k];

        t_314[k] = f_3 * pc_x[k] * msp_158[k];

        t_315[k] = f_8 * lsp_130[k]
                   + f_1 * mss0_52[k]
                   - f_2 * mss1_52[k]
                   + f_3 * pc_y[k] * msp_157[k];

        t_316[k] = f_8 * lsp_131[k]
                   + f_3 * pc_y[k] * msp_158[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pa_y, pc_x, pc_y, pc_z, lsd0_264, \
                         lsp_128, lsd1_264, mss0_52, mss1_52, msp_158, msp_160, \
                         msp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_7 * lsp_128[k]
                   + f_1 * mss0_52[k]
                   - f_2 * mss1_52[k]
                   + f_3 * pc_z[k] * msp_158[k];

        t_318[k] = pa_y[k] * lsd0_264[k]
                   - f_4 * pc_y[k] * lsd1_264[k];

        t_319[k] = f_3 * pc_x[k] * msp_160[k];

        t_320[k] = f_3 * pc_x[k] * msp_161[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_y, pc_y, lsd0_267, lsd0_269, lsp_133, \
                         lsp_134, lsd1_267, lsd1_269, msp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = pa_y[k] * lsd0_267[k]
                   + f_8 * lsp_133[k]
                   - f_4 * pc_y[k] * lsd1_267[k];

        t_322[k] = f_6 * lsp_134[k]
                   + f_3 * pc_y[k] * msp_161[k];

        t_323[k] = pa_y[k] * lsd0_269[k]
                   - f_4 * pc_y[k] * lsd1_269[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pc_x, pc_y, pc_z, lsp_134, \
                         mss0_54, mss1_54, msp_162, msp_163, msp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * mss0_54[k]
                   - f_2 * mss1_54[k]
                   + f_3 * pc_x[k] * msp_162[k];

        t_325[k] = f_3 * pc_x[k] * msp_163[k];

        t_326[k] = f_3 * pc_x[k] * msp_164[k];

        t_327[k] = f_1 * mss0_54[k]
                   - f_2 * mss1_54[k]
                   + f_3 * pc_y[k] * msp_163[k];

        t_328[k] = f_3 * pc_y[k] * msp_164[k];

        t_329[k] = f_0 * lsp_134[k]
                   + f_1 * mss0_54[k]
                   - f_2 * mss1_54[k]
                   + f_3 * pc_z[k] * msp_164[k];
    }
}

auto
compute_prim_msd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t lsd0, const size_t lsp,
                                                   const size_t lsd1, const size_t mss0,
                                                   const size_t mss1, const size_t msp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_msd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, lsd0, lsp,
                                                              lsd1, mss0, mss1, msp, ncols,
                                                              gamma, p, q);

    compute_prim_msd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, lsd0, lsp,
                                                              lsd1, mss0, mss1, msp, ncols,
                                                              gamma, p, q);

    compute_prim_msd_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, lsd0, lsp,
                                                              lsd1, mss0, mss1, msp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
