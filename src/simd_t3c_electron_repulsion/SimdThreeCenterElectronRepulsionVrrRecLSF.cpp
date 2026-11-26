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


#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksf0,
                                                          const size_t ksd, const size_t ksf1,
                                                          const size_t lsp0, const size_t lsp1,
                                                          const size_t lsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksf0_0 = buffer.data(ksf0 + 0);
    const auto *ksf0_6 = buffer.data(ksf0 + 6);
    const auto *ksf0_9 = buffer.data(ksf0 + 9);
    const auto *ksf0_16 = buffer.data(ksf0 + 16);
    const auto *ksf0_20 = buffer.data(ksf0 + 20);
    const auto *ksf0_29 = buffer.data(ksf0 + 29);
    const auto *ksf0_30 = buffer.data(ksf0 + 30);
    const auto *ksf0_36 = buffer.data(ksf0 + 36);
    const auto *ksf0_50 = buffer.data(ksf0 + 50);
    const auto *ksf0_59 = buffer.data(ksf0 + 59);
    const auto *ksf0_60 = buffer.data(ksf0 + 60);
    const auto *ksf0_66 = buffer.data(ksf0 + 66);
    const auto *ksf0_90 = buffer.data(ksf0 + 90);

    const auto *ksd_0 = buffer.data(ksd + 0);
    const auto *ksd_3 = buffer.data(ksd + 3);
    const auto *ksd_5 = buffer.data(ksd + 5);
    const auto *ksd_6 = buffer.data(ksd + 6);
    const auto *ksd_9 = buffer.data(ksd + 9);
    const auto *ksd_11 = buffer.data(ksd + 11);
    const auto *ksd_12 = buffer.data(ksd + 12);
    const auto *ksd_15 = buffer.data(ksd + 15);
    const auto *ksd_17 = buffer.data(ksd + 17);
    const auto *ksd_18 = buffer.data(ksd + 18);
    const auto *ksd_21 = buffer.data(ksd + 21);
    const auto *ksd_23 = buffer.data(ksd + 23);
    const auto *ksd_24 = buffer.data(ksd + 24);
    const auto *ksd_27 = buffer.data(ksd + 27);
    const auto *ksd_28 = buffer.data(ksd + 28);
    const auto *ksd_29 = buffer.data(ksd + 29);
    const auto *ksd_30 = buffer.data(ksd + 30);
    const auto *ksd_33 = buffer.data(ksd + 33);
    const auto *ksd_35 = buffer.data(ksd + 35);
    const auto *ksd_36 = buffer.data(ksd + 36);
    const auto *ksd_39 = buffer.data(ksd + 39);
    const auto *ksd_41 = buffer.data(ksd + 41);
    const auto *ksd_42 = buffer.data(ksd + 42);
    const auto *ksd_45 = buffer.data(ksd + 45);
    const auto *ksd_46 = buffer.data(ksd + 46);
    const auto *ksd_47 = buffer.data(ksd + 47);
    const auto *ksd_48 = buffer.data(ksd + 48);
    const auto *ksd_51 = buffer.data(ksd + 51);
    const auto *ksd_52 = buffer.data(ksd + 52);
    const auto *ksd_53 = buffer.data(ksd + 53);
    const auto *ksd_54 = buffer.data(ksd + 54);
    const auto *ksd_57 = buffer.data(ksd + 57);
    const auto *ksd_59 = buffer.data(ksd + 59);
    const auto *ksd_60 = buffer.data(ksd + 60);
    const auto *ksd_63 = buffer.data(ksd + 63);
    const auto *ksd_65 = buffer.data(ksd + 65);
    const auto *ksd_69 = buffer.data(ksd + 69);
    const auto *ksd_70 = buffer.data(ksd + 70);
    const auto *ksd_71 = buffer.data(ksd + 71);
    const auto *ksd_72 = buffer.data(ksd + 72);
    const auto *ksd_75 = buffer.data(ksd + 75);
    const auto *ksd_76 = buffer.data(ksd + 76);
    const auto *ksd_77 = buffer.data(ksd + 77);

    const auto *ksf1_0 = buffer.data(ksf1 + 0);
    const auto *ksf1_6 = buffer.data(ksf1 + 6);
    const auto *ksf1_9 = buffer.data(ksf1 + 9);
    const auto *ksf1_16 = buffer.data(ksf1 + 16);
    const auto *ksf1_20 = buffer.data(ksf1 + 20);
    const auto *ksf1_29 = buffer.data(ksf1 + 29);
    const auto *ksf1_30 = buffer.data(ksf1 + 30);
    const auto *ksf1_36 = buffer.data(ksf1 + 36);
    const auto *ksf1_50 = buffer.data(ksf1 + 50);
    const auto *ksf1_59 = buffer.data(ksf1 + 59);
    const auto *ksf1_60 = buffer.data(ksf1 + 60);
    const auto *ksf1_66 = buffer.data(ksf1 + 66);
    const auto *ksf1_90 = buffer.data(ksf1 + 90);

    const auto *lsp0_0 = buffer.data(lsp0 + 0);
    const auto *lsp0_1 = buffer.data(lsp0 + 1);
    const auto *lsp0_2 = buffer.data(lsp0 + 2);
    const auto *lsp0_4 = buffer.data(lsp0 + 4);
    const auto *lsp0_8 = buffer.data(lsp0 + 8);
    const auto *lsp0_9 = buffer.data(lsp0 + 9);
    const auto *lsp0_10 = buffer.data(lsp0 + 10);
    const auto *lsp0_11 = buffer.data(lsp0 + 11);
    const auto *lsp0_15 = buffer.data(lsp0 + 15);
    const auto *lsp0_16 = buffer.data(lsp0 + 16);
    const auto *lsp0_17 = buffer.data(lsp0 + 17);
    const auto *lsp0_18 = buffer.data(lsp0 + 18);
    const auto *lsp0_19 = buffer.data(lsp0 + 19);
    const auto *lsp0_20 = buffer.data(lsp0 + 20);
    const auto *lsp0_23 = buffer.data(lsp0 + 23);
    const auto *lsp0_25 = buffer.data(lsp0 + 25);
    const auto *lsp0_27 = buffer.data(lsp0 + 27);
    const auto *lsp0_28 = buffer.data(lsp0 + 28);
    const auto *lsp0_29 = buffer.data(lsp0 + 29);
    const auto *lsp0_30 = buffer.data(lsp0 + 30);
    const auto *lsp0_31 = buffer.data(lsp0 + 31);
    const auto *lsp0_32 = buffer.data(lsp0 + 32);
    const auto *lsp0_35 = buffer.data(lsp0 + 35);
    const auto *lsp0_36 = buffer.data(lsp0 + 36);
    const auto *lsp0_37 = buffer.data(lsp0 + 37);
    const auto *lsp0_38 = buffer.data(lsp0 + 38);

    const auto *lsp1_0 = buffer.data(lsp1 + 0);
    const auto *lsp1_1 = buffer.data(lsp1 + 1);
    const auto *lsp1_2 = buffer.data(lsp1 + 2);
    const auto *lsp1_4 = buffer.data(lsp1 + 4);
    const auto *lsp1_8 = buffer.data(lsp1 + 8);
    const auto *lsp1_9 = buffer.data(lsp1 + 9);
    const auto *lsp1_10 = buffer.data(lsp1 + 10);
    const auto *lsp1_11 = buffer.data(lsp1 + 11);
    const auto *lsp1_15 = buffer.data(lsp1 + 15);
    const auto *lsp1_16 = buffer.data(lsp1 + 16);
    const auto *lsp1_17 = buffer.data(lsp1 + 17);
    const auto *lsp1_18 = buffer.data(lsp1 + 18);
    const auto *lsp1_19 = buffer.data(lsp1 + 19);
    const auto *lsp1_20 = buffer.data(lsp1 + 20);
    const auto *lsp1_23 = buffer.data(lsp1 + 23);
    const auto *lsp1_25 = buffer.data(lsp1 + 25);
    const auto *lsp1_27 = buffer.data(lsp1 + 27);
    const auto *lsp1_28 = buffer.data(lsp1 + 28);
    const auto *lsp1_29 = buffer.data(lsp1 + 29);
    const auto *lsp1_30 = buffer.data(lsp1 + 30);
    const auto *lsp1_31 = buffer.data(lsp1 + 31);
    const auto *lsp1_32 = buffer.data(lsp1 + 32);
    const auto *lsp1_35 = buffer.data(lsp1 + 35);
    const auto *lsp1_36 = buffer.data(lsp1 + 36);
    const auto *lsp1_37 = buffer.data(lsp1 + 37);
    const auto *lsp1_38 = buffer.data(lsp1 + 38);

    const auto *lsd_0 = buffer.data(lsd + 0);
    const auto *lsd_2 = buffer.data(lsd + 2);
    const auto *lsd_3 = buffer.data(lsd + 3);
    const auto *lsd_5 = buffer.data(lsd + 5);
    const auto *lsd_6 = buffer.data(lsd + 6);
    const auto *lsd_7 = buffer.data(lsd + 7);
    const auto *lsd_9 = buffer.data(lsd + 9);
    const auto *lsd_11 = buffer.data(lsd + 11);
    const auto *lsd_12 = buffer.data(lsd + 12);
    const auto *lsd_14 = buffer.data(lsd + 14);
    const auto *lsd_15 = buffer.data(lsd + 15);
    const auto *lsd_16 = buffer.data(lsd + 16);
    const auto *lsd_17 = buffer.data(lsd + 17);
    const auto *lsd_18 = buffer.data(lsd + 18);
    const auto *lsd_19 = buffer.data(lsd + 19);
    const auto *lsd_21 = buffer.data(lsd + 21);
    const auto *lsd_23 = buffer.data(lsd + 23);
    const auto *lsd_24 = buffer.data(lsd + 24);
    const auto *lsd_27 = buffer.data(lsd + 27);
    const auto *lsd_28 = buffer.data(lsd + 28);
    const auto *lsd_29 = buffer.data(lsd + 29);
    const auto *lsd_30 = buffer.data(lsd + 30);
    const auto *lsd_32 = buffer.data(lsd + 32);
    const auto *lsd_33 = buffer.data(lsd + 33);
    const auto *lsd_34 = buffer.data(lsd + 34);
    const auto *lsd_35 = buffer.data(lsd + 35);
    const auto *lsd_36 = buffer.data(lsd + 36);
    const auto *lsd_37 = buffer.data(lsd + 37);
    const auto *lsd_39 = buffer.data(lsd + 39);
    const auto *lsd_41 = buffer.data(lsd + 41);
    const auto *lsd_42 = buffer.data(lsd + 42);
    const auto *lsd_45 = buffer.data(lsd + 45);
    const auto *lsd_46 = buffer.data(lsd + 46);
    const auto *lsd_47 = buffer.data(lsd + 47);
    const auto *lsd_48 = buffer.data(lsd + 48);
    const auto *lsd_51 = buffer.data(lsd + 51);
    const auto *lsd_52 = buffer.data(lsd + 52);
    const auto *lsd_53 = buffer.data(lsd + 53);
    const auto *lsd_54 = buffer.data(lsd + 54);
    const auto *lsd_56 = buffer.data(lsd + 56);
    const auto *lsd_57 = buffer.data(lsd + 57);
    const auto *lsd_58 = buffer.data(lsd + 58);
    const auto *lsd_59 = buffer.data(lsd + 59);
    const auto *lsd_60 = buffer.data(lsd + 60);
    const auto *lsd_61 = buffer.data(lsd + 61);
    const auto *lsd_63 = buffer.data(lsd + 63);
    const auto *lsd_65 = buffer.data(lsd + 65);
    const auto *lsd_66 = buffer.data(lsd + 66);
    const auto *lsd_69 = buffer.data(lsd + 69);
    const auto *lsd_70 = buffer.data(lsd + 70);
    const auto *lsd_71 = buffer.data(lsd + 71);
    const auto *lsd_72 = buffer.data(lsd + 72);
    const auto *lsd_75 = buffer.data(lsd + 75);
    const auto *lsd_76 = buffer.data(lsd + 76);
    const auto *lsd_77 = buffer.data(lsd + 77);
    const auto *lsd_78 = buffer.data(lsd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ksd_0, ksd_3, lsp0_0, \
                         lsp1_0, lsd_0, lsd_2, lsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksd_0[k]
                 + f_1 * lsp0_0[k]
                 - f_2 * lsp1_0[k]
                 + f_3 * pc_x[k] * lsd_0[k];

        t_1[k] = f_3 * pc_y[k] * lsd_0[k];

        t_2[k] = f_3 * pc_z[k] * lsd_0[k];

        t_3[k] = f_0 * ksd_3[k]
                 + f_3 * pc_x[k] * lsd_3[k];

        t_4[k] = f_3 * pc_y[k] * lsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, ksd_5, lsp0_1, lsp0_2, \
                         lsp1_1, lsp1_2, lsd_3, lsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * ksd_5[k]
                 + f_3 * pc_x[k] * lsd_5[k];

        t_6[k] = f_1 * lsp0_1[k]
                 - f_2 * lsp1_1[k]
                 + f_3 * pc_y[k] * lsd_3[k];

        t_7[k] = f_3 * pc_z[k] * lsd_3[k];

        t_8[k] = f_3 * pc_y[k] * lsd_5[k];

        t_9[k] = f_1 * lsp0_2[k]
                 - f_2 * lsp1_2[k]
                 + f_3 * pc_z[k] * lsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, ksf0_0, ksd_0, \
                         ksd_9, ksf1_0, lsd_6, lsd_7, lsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * ksf0_0[k]
                  - f_4 * pc_y[k] * ksf1_0[k];

        t_11[k] = f_5 * ksd_0[k]
                  + f_3 * pc_y[k] * lsd_6[k];

        t_12[k] = f_3 * pc_z[k] * lsd_6[k];

        t_13[k] = f_6 * ksd_9[k]
                  + f_3 * pc_x[k] * lsd_9[k];

        t_14[k] = f_3 * pc_z[k] * lsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, ksd_3, ksd_5, ksd_11, \
                         lsp0_4, lsp1_4, lsd_9, lsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * ksd_11[k]
                  + f_3 * pc_x[k] * lsd_11[k];

        t_16[k] = f_5 * ksd_3[k]
                  + f_1 * lsp0_4[k]
                  - f_2 * lsp1_4[k]
                  + f_3 * pc_y[k] * lsd_9[k];

        t_17[k] = f_3 * pc_z[k] * lsd_9[k];

        t_18[k] = f_5 * ksd_5[k]
                  + f_3 * pc_y[k] * lsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, ksf0_0, ksf0_9, \
                         ksd_0, ksf1_0, ksf1_9, lsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * ksf0_9[k]
                  - f_4 * pc_y[k] * ksf1_9[k];

        t_20[k] = pa_z[k] * ksf0_0[k]
                  - f_4 * pc_z[k] * ksf1_0[k];

        t_21[k] = f_3 * pc_y[k] * lsd_12[k];

        t_22[k] = f_5 * ksd_0[k]
                  + f_3 * pc_z[k] * lsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, ksf0_6, ksd_15, \
                         ksd_17, ksf1_6, lsd_14, lsd_15, lsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * ksd_15[k]
                  + f_3 * pc_x[k] * lsd_15[k];

        t_24[k] = f_3 * pc_y[k] * lsd_14[k];

        t_25[k] = f_6 * ksd_17[k]
                  + f_3 * pc_x[k] * lsd_17[k];

        t_26[k] = pa_z[k] * ksf0_6[k]
                  - f_4 * pc_z[k] * ksf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, ksd_5, ksd_18, lsp0_8, \
                         lsp0_9, lsp1_8, lsp1_9, lsd_16, lsd_17, \
                         lsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * lsp0_8[k]
                  - f_8 * lsp1_8[k]
                  + f_3 * pc_y[k] * lsd_16[k];

        t_28[k] = f_3 * pc_y[k] * lsd_17[k];

        t_29[k] = f_5 * ksd_5[k]
                  + f_1 * lsp0_8[k]
                  - f_2 * lsp1_8[k]
                  + f_3 * pc_z[k] * lsd_17[k];

        t_30[k] = f_9 * ksd_18[k]
                  + f_1 * lsp0_9[k]
                  - f_2 * lsp1_9[k]
                  + f_3 * pc_x[k] * lsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, ksd_6, ksd_21, \
                         ksd_23, lsd_18, lsd_19, lsd_21, lsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * ksd_6[k]
                  + f_3 * pc_y[k] * lsd_18[k];

        t_32[k] = f_3 * pc_z[k] * lsd_18[k];

        t_33[k] = f_9 * ksd_21[k]
                  + f_3 * pc_x[k] * lsd_21[k];

        t_34[k] = f_3 * pc_z[k] * lsd_19[k];

        t_35[k] = f_9 * ksd_23[k]
                  + f_3 * pc_x[k] * lsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, ksd_9, ksd_11, lsp0_10, lsp0_11, \
                         lsp1_10, lsp1_11, lsd_21, lsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * ksd_9[k]
                  + f_1 * lsp0_10[k]
                  - f_2 * lsp1_10[k]
                  + f_3 * pc_y[k] * lsd_21[k];

        t_37[k] = f_3 * pc_z[k] * lsd_21[k];

        t_38[k] = f_10 * ksd_11[k]
                  + f_3 * pc_y[k] * lsd_23[k];

        t_39[k] = f_1 * lsp0_11[k]
                  - f_2 * lsp1_11[k]
                  + f_3 * pc_z[k] * lsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, ksf0_20, ksd_6, \
                         ksd_12, ksd_27, ksf1_20, lsd_24, lsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * ksf0_20[k]
                  - f_4 * pc_y[k] * ksf1_20[k];

        t_41[k] = f_5 * ksd_12[k]
                  + f_3 * pc_y[k] * lsd_24[k];

        t_42[k] = f_5 * ksd_6[k]
                  + f_3 * pc_z[k] * lsd_24[k];

        t_43[k] = f_9 * ksd_27[k]
                  + f_3 * pc_x[k] * lsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, ksf0_16, ksd_9, ksd_28, \
                         ksd_29, ksf1_16, lsd_27, lsd_28, lsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * ksd_28[k]
                  + f_3 * pc_x[k] * lsd_28[k];

        t_45[k] = f_9 * ksd_29[k]
                  + f_3 * pc_x[k] * lsd_29[k];

        t_46[k] = pa_z[k] * ksf0_16[k]
                  - f_4 * pc_z[k] * ksf1_16[k];

        t_47[k] = f_5 * ksd_9[k]
                  + f_3 * pc_z[k] * lsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, ksf0_29, ksd_17, ksd_30, \
                         ksf1_29, lsp0_15, lsp1_15, lsd_29, lsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * ksd_17[k]
                  + f_3 * pc_y[k] * lsd_29[k];

        t_49[k] = pa_y[k] * ksf0_29[k]
                  - f_4 * pc_y[k] * ksf1_29[k];

        t_50[k] = f_9 * ksd_30[k]
                  + f_1 * lsp0_15[k]
                  - f_2 * lsp1_15[k]
                  + f_3 * pc_x[k] * lsd_30[k];

        t_51[k] = f_3 * pc_y[k] * lsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, ksd_12, ksd_33, ksd_35, \
                         lsd_30, lsd_32, lsd_33, lsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * ksd_12[k]
                  + f_3 * pc_z[k] * lsd_30[k];

        t_53[k] = f_9 * ksd_33[k]
                  + f_3 * pc_x[k] * lsd_33[k];

        t_54[k] = f_3 * pc_y[k] * lsd_32[k];

        t_55[k] = f_9 * ksd_35[k]
                  + f_3 * pc_x[k] * lsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, ksd_17, lsp0_16, lsp0_17, \
                         lsp1_16, lsp1_17, lsd_33, lsd_34, lsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * lsp0_16[k]
                  - f_2 * lsp1_16[k]
                  + f_3 * pc_y[k] * lsd_33[k];

        t_57[k] = f_7 * lsp0_17[k]
                  - f_8 * lsp1_17[k]
                  + f_3 * pc_y[k] * lsd_34[k];

        t_58[k] = f_3 * pc_y[k] * lsd_35[k];

        t_59[k] = f_10 * ksd_17[k]
                  + f_1 * lsp0_17[k]
                  - f_2 * lsp1_17[k]
                  + f_3 * pc_z[k] * lsd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, ksd_18, ksd_36, \
                         ksd_39, lsp0_18, lsp1_18, lsd_36, lsd_37, \
                         lsd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * ksd_36[k]
                  + f_1 * lsp0_18[k]
                  - f_2 * lsp1_18[k]
                  + f_3 * pc_x[k] * lsd_36[k];

        t_61[k] = f_12 * ksd_18[k]
                  + f_3 * pc_y[k] * lsd_36[k];

        t_62[k] = f_3 * pc_z[k] * lsd_36[k];

        t_63[k] = f_11 * ksd_39[k]
                  + f_3 * pc_x[k] * lsd_39[k];

        t_64[k] = f_3 * pc_z[k] * lsd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, ksd_21, ksd_23, ksd_41, \
                         lsp0_19, lsp1_19, lsd_39, lsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * ksd_41[k]
                  + f_3 * pc_x[k] * lsd_41[k];

        t_66[k] = f_12 * ksd_21[k]
                  + f_1 * lsp0_19[k]
                  - f_2 * lsp1_19[k]
                  + f_3 * pc_y[k] * lsd_39[k];

        t_67[k] = f_3 * pc_z[k] * lsd_39[k];

        t_68[k] = f_12 * ksd_23[k]
                  + f_3 * pc_y[k] * lsd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, ksf0_30, ksd_18, ksd_24, \
                         ksf1_30, lsp0_20, lsp1_20, lsd_41, lsd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * lsp0_20[k]
                  - f_2 * lsp1_20[k]
                  + f_3 * pc_z[k] * lsd_41[k];

        t_70[k] = pa_z[k] * ksf0_30[k]
                  - f_4 * pc_z[k] * ksf1_30[k];

        t_71[k] = f_10 * ksd_24[k]
                  + f_3 * pc_y[k] * lsd_42[k];

        t_72[k] = f_5 * ksd_18[k]
                  + f_3 * pc_z[k] * lsd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, ksf0_36, ksd_45, ksd_46, \
                         ksd_47, ksf1_36, lsd_45, lsd_46, lsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * ksd_45[k]
                  + f_3 * pc_x[k] * lsd_45[k];

        t_74[k] = f_11 * ksd_46[k]
                  + f_3 * pc_x[k] * lsd_46[k];

        t_75[k] = f_11 * ksd_47[k]
                  + f_3 * pc_x[k] * lsd_47[k];

        t_76[k] = pa_z[k] * ksf0_36[k]
                  - f_4 * pc_z[k] * ksf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, ksf0_50, ksd_21, ksd_23, \
                         ksd_29, ksf1_50, lsp0_23, lsp1_23, lsd_45, \
                         lsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * ksd_21[k]
                  + f_3 * pc_z[k] * lsd_45[k];

        t_78[k] = f_10 * ksd_29[k]
                  + f_3 * pc_y[k] * lsd_47[k];

        t_79[k] = f_5 * ksd_23[k]
                  + f_1 * lsp0_23[k]
                  - f_2 * lsp1_23[k]
                  + f_3 * pc_z[k] * lsd_47[k];

        t_80[k] = pa_y[k] * ksf0_50[k]
                  - f_4 * pc_y[k] * ksf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, ksd_24, ksd_30, ksd_51, \
                         ksd_52, lsd_48, lsd_51, lsd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * ksd_30[k]
                  + f_3 * pc_y[k] * lsd_48[k];

        t_82[k] = f_10 * ksd_24[k]
                  + f_3 * pc_z[k] * lsd_48[k];

        t_83[k] = f_11 * ksd_51[k]
                  + f_3 * pc_x[k] * lsd_51[k];

        t_84[k] = f_11 * ksd_52[k]
                  + f_3 * pc_x[k] * lsd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, ksd_27, ksd_33, ksd_35, \
                         ksd_53, lsp0_25, lsp1_25, lsd_51, lsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * ksd_53[k]
                  + f_3 * pc_x[k] * lsd_53[k];

        t_86[k] = f_5 * ksd_33[k]
                  + f_1 * lsp0_25[k]
                  - f_2 * lsp1_25[k]
                  + f_3 * pc_y[k] * lsd_51[k];

        t_87[k] = f_10 * ksd_27[k]
                  + f_3 * pc_z[k] * lsd_51[k];

        t_88[k] = f_5 * ksd_35[k]
                  + f_3 * pc_y[k] * lsd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, ksf0_59, ksd_30, \
                         ksd_54, ksf1_59, lsp0_27, lsp1_27, lsd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * ksf0_59[k]
                  - f_4 * pc_y[k] * ksf1_59[k];

        t_90[k] = f_11 * ksd_54[k]
                  + f_1 * lsp0_27[k]
                  - f_2 * lsp1_27[k]
                  + f_3 * pc_x[k] * lsd_54[k];

        t_91[k] = f_3 * pc_y[k] * lsd_54[k];

        t_92[k] = f_12 * ksd_30[k]
                  + f_3 * pc_z[k] * lsd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, ksd_57, ksd_59, lsp0_28, lsp1_28, \
                         lsd_56, lsd_57, lsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * ksd_57[k]
                  + f_3 * pc_x[k] * lsd_57[k];

        t_94[k] = f_3 * pc_y[k] * lsd_56[k];

        t_95[k] = f_11 * ksd_59[k]
                  + f_3 * pc_x[k] * lsd_59[k];

        t_96[k] = f_1 * lsp0_28[k]
                  - f_2 * lsp1_28[k]
                  + f_3 * pc_y[k] * lsd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, ksd_35, ksd_60, lsp0_29, \
                         lsp0_30, lsp1_29, lsp1_30, lsd_58, lsd_59, \
                         lsd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * lsp0_29[k]
                  - f_8 * lsp1_29[k]
                  + f_3 * pc_y[k] * lsd_58[k];

        t_98[k] = f_3 * pc_y[k] * lsd_59[k];

        t_99[k] = f_12 * ksd_35[k]
                  + f_1 * lsp0_29[k]
                  - f_2 * lsp1_29[k]
                  + f_3 * pc_z[k] * lsd_59[k];

        t_100[k] = f_13 * ksd_60[k]
                   + f_1 * lsp0_30[k]
                   - f_2 * lsp1_30[k]
                   + f_3 * pc_x[k] * lsd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, ksd_36, ksd_63, \
                         ksd_65, lsd_60, lsd_61, lsd_63, lsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_13 * ksd_36[k]
                   + f_3 * pc_y[k] * lsd_60[k];

        t_102[k] = f_3 * pc_z[k] * lsd_60[k];

        t_103[k] = f_13 * ksd_63[k]
                   + f_3 * pc_x[k] * lsd_63[k];

        t_104[k] = f_3 * pc_z[k] * lsd_61[k];

        t_105[k] = f_13 * ksd_65[k]
                   + f_3 * pc_x[k] * lsd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, ksd_39, ksd_41, lsp0_31, \
                         lsp0_32, lsp1_31, lsp1_32, lsd_63, lsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_13 * ksd_39[k]
                   + f_1 * lsp0_31[k]
                   - f_2 * lsp1_31[k]
                   + f_3 * pc_y[k] * lsd_63[k];

        t_107[k] = f_3 * pc_z[k] * lsd_63[k];

        t_108[k] = f_13 * ksd_41[k]
                   + f_3 * pc_y[k] * lsd_65[k];

        t_109[k] = f_1 * lsp0_32[k]
                   - f_2 * lsp1_32[k]
                   + f_3 * pc_z[k] * lsd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, ksf0_60, ksd_36, \
                         ksd_42, ksd_69, ksf1_60, lsd_66, lsd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * ksf0_60[k]
                   - f_4 * pc_z[k] * ksf1_60[k];

        t_111[k] = f_12 * ksd_42[k]
                   + f_3 * pc_y[k] * lsd_66[k];

        t_112[k] = f_5 * ksd_36[k]
                   + f_3 * pc_z[k] * lsd_66[k];

        t_113[k] = f_13 * ksd_69[k]
                   + f_3 * pc_x[k] * lsd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, ksf0_66, ksd_39, \
                         ksd_70, ksd_71, ksf1_66, lsd_69, lsd_70, \
                         lsd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * ksd_70[k]
                   + f_3 * pc_x[k] * lsd_70[k];

        t_115[k] = f_13 * ksd_71[k]
                   + f_3 * pc_x[k] * lsd_71[k];

        t_116[k] = pa_z[k] * ksf0_66[k]
                   - f_4 * pc_z[k] * ksf1_66[k];

        t_117[k] = f_5 * ksd_39[k]
                   + f_3 * pc_z[k] * lsd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, ksd_41, ksd_47, ksd_72, \
                         lsp0_35, lsp0_36, lsp1_35, lsp1_36, lsd_71, \
                         lsd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * ksd_47[k]
                   + f_3 * pc_y[k] * lsd_71[k];

        t_119[k] = f_5 * ksd_41[k]
                   + f_1 * lsp0_35[k]
                   - f_2 * lsp1_35[k]
                   + f_3 * pc_z[k] * lsd_71[k];

        t_120[k] = f_13 * ksd_72[k]
                   + f_1 * lsp0_36[k]
                   - f_2 * lsp1_36[k]
                   + f_3 * pc_x[k] * lsd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, ksd_42, ksd_48, ksd_75, \
                         ksd_76, lsd_72, lsd_75, lsd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * ksd_48[k]
                   + f_3 * pc_y[k] * lsd_72[k];

        t_122[k] = f_10 * ksd_42[k]
                   + f_3 * pc_z[k] * lsd_72[k];

        t_123[k] = f_13 * ksd_75[k]
                   + f_3 * pc_x[k] * lsd_75[k];

        t_124[k] = f_13 * ksd_76[k]
                   + f_3 * pc_x[k] * lsd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, ksd_45, ksd_51, ksd_53, \
                         ksd_77, lsp0_37, lsp1_37, lsd_75, lsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * ksd_77[k]
                   + f_3 * pc_x[k] * lsd_77[k];

        t_126[k] = f_10 * ksd_51[k]
                   + f_1 * lsp0_37[k]
                   - f_2 * lsp1_37[k]
                   + f_3 * pc_y[k] * lsd_75[k];

        t_127[k] = f_10 * ksd_45[k]
                   + f_3 * pc_z[k] * lsd_75[k];

        t_128[k] = f_10 * ksd_53[k]
                   + f_3 * pc_y[k] * lsd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, ksf0_90, ksd_47, \
                         ksd_48, ksd_54, ksf1_90, lsp0_38, lsp1_38, lsd_77, \
                         lsd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * ksd_47[k]
                   + f_1 * lsp0_38[k]
                   - f_2 * lsp1_38[k]
                   + f_3 * pc_z[k] * lsd_77[k];

        t_130[k] = pa_y[k] * ksf0_90[k]
                   - f_4 * pc_y[k] * ksf1_90[k];

        t_131[k] = f_5 * ksd_54[k]
                   + f_3 * pc_y[k] * lsd_78[k];

        t_132[k] = f_12 * ksd_48[k]
                   + f_3 * pc_z[k] * lsd_78[k];
    }
}

static auto
compute_prim_lsf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksf0,
                                                          const size_t ksd, const size_t ksf1,
                                                          const size_t lsp0, const size_t lsp1,
                                                          const size_t lsd, const size_t ncols,
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
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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
    auto *t_255 = buffer.data(target + 255);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksf0_99 = buffer.data(ksf0 + 99);
    const auto *ksf0_100 = buffer.data(ksf0 + 100);
    const auto *ksf0_106 = buffer.data(ksf0 + 106);
    const auto *ksf0_140 = buffer.data(ksf0 + 140);
    const auto *ksf0_149 = buffer.data(ksf0 + 149);
    const auto *ksf0_150 = buffer.data(ksf0 + 150);
    const auto *ksf0_156 = buffer.data(ksf0 + 156);

    const auto *ksd_51 = buffer.data(ksd + 51);
    const auto *ksd_54 = buffer.data(ksd + 54);
    const auto *ksd_57 = buffer.data(ksd + 57);
    const auto *ksd_59 = buffer.data(ksd + 59);
    const auto *ksd_60 = buffer.data(ksd + 60);
    const auto *ksd_63 = buffer.data(ksd + 63);
    const auto *ksd_65 = buffer.data(ksd + 65);
    const auto *ksd_66 = buffer.data(ksd + 66);
    const auto *ksd_69 = buffer.data(ksd + 69);
    const auto *ksd_71 = buffer.data(ksd + 71);
    const auto *ksd_72 = buffer.data(ksd + 72);
    const auto *ksd_75 = buffer.data(ksd + 75);
    const auto *ksd_77 = buffer.data(ksd + 77);
    const auto *ksd_78 = buffer.data(ksd + 78);
    const auto *ksd_81 = buffer.data(ksd + 81);
    const auto *ksd_82 = buffer.data(ksd + 82);
    const auto *ksd_83 = buffer.data(ksd + 83);
    const auto *ksd_84 = buffer.data(ksd + 84);
    const auto *ksd_87 = buffer.data(ksd + 87);
    const auto *ksd_89 = buffer.data(ksd + 89);
    const auto *ksd_90 = buffer.data(ksd + 90);
    const auto *ksd_93 = buffer.data(ksd + 93);
    const auto *ksd_95 = buffer.data(ksd + 95);
    const auto *ksd_96 = buffer.data(ksd + 96);
    const auto *ksd_99 = buffer.data(ksd + 99);
    const auto *ksd_100 = buffer.data(ksd + 100);
    const auto *ksd_101 = buffer.data(ksd + 101);
    const auto *ksd_102 = buffer.data(ksd + 102);
    const auto *ksd_105 = buffer.data(ksd + 105);
    const auto *ksd_106 = buffer.data(ksd + 106);
    const auto *ksd_107 = buffer.data(ksd + 107);
    const auto *ksd_108 = buffer.data(ksd + 108);
    const auto *ksd_111 = buffer.data(ksd + 111);
    const auto *ksd_112 = buffer.data(ksd + 112);
    const auto *ksd_113 = buffer.data(ksd + 113);
    const auto *ksd_114 = buffer.data(ksd + 114);
    const auto *ksd_117 = buffer.data(ksd + 117);
    const auto *ksd_118 = buffer.data(ksd + 118);
    const auto *ksd_119 = buffer.data(ksd + 119);
    const auto *ksd_120 = buffer.data(ksd + 120);
    const auto *ksd_123 = buffer.data(ksd + 123);
    const auto *ksd_125 = buffer.data(ksd + 125);
    const auto *ksd_126 = buffer.data(ksd + 126);
    const auto *ksd_129 = buffer.data(ksd + 129);
    const auto *ksd_131 = buffer.data(ksd + 131);
    const auto *ksd_135 = buffer.data(ksd + 135);
    const auto *ksd_136 = buffer.data(ksd + 136);
    const auto *ksd_137 = buffer.data(ksd + 137);
    const auto *ksd_138 = buffer.data(ksd + 138);
    const auto *ksd_141 = buffer.data(ksd + 141);
    const auto *ksd_142 = buffer.data(ksd + 142);
    const auto *ksd_143 = buffer.data(ksd + 143);
    const auto *ksd_144 = buffer.data(ksd + 144);
    const auto *ksd_147 = buffer.data(ksd + 147);
    const auto *ksd_148 = buffer.data(ksd + 148);
    const auto *ksd_149 = buffer.data(ksd + 149);
    const auto *ksd_150 = buffer.data(ksd + 150);
    const auto *ksd_153 = buffer.data(ksd + 153);
    const auto *ksd_154 = buffer.data(ksd + 154);
    const auto *ksd_155 = buffer.data(ksd + 155);

    const auto *ksf1_99 = buffer.data(ksf1 + 99);
    const auto *ksf1_100 = buffer.data(ksf1 + 100);
    const auto *ksf1_106 = buffer.data(ksf1 + 106);
    const auto *ksf1_140 = buffer.data(ksf1 + 140);
    const auto *ksf1_149 = buffer.data(ksf1 + 149);
    const auto *ksf1_150 = buffer.data(ksf1 + 150);
    const auto *ksf1_156 = buffer.data(ksf1 + 156);

    const auto *lsp0_40 = buffer.data(lsp0 + 40);
    const auto *lsp0_42 = buffer.data(lsp0 + 42);
    const auto *lsp0_43 = buffer.data(lsp0 + 43);
    const auto *lsp0_44 = buffer.data(lsp0 + 44);
    const auto *lsp0_45 = buffer.data(lsp0 + 45);
    const auto *lsp0_46 = buffer.data(lsp0 + 46);
    const auto *lsp0_47 = buffer.data(lsp0 + 47);
    const auto *lsp0_50 = buffer.data(lsp0 + 50);
    const auto *lsp0_51 = buffer.data(lsp0 + 51);
    const auto *lsp0_52 = buffer.data(lsp0 + 52);
    const auto *lsp0_53 = buffer.data(lsp0 + 53);
    const auto *lsp0_54 = buffer.data(lsp0 + 54);
    const auto *lsp0_55 = buffer.data(lsp0 + 55);
    const auto *lsp0_56 = buffer.data(lsp0 + 56);
    const auto *lsp0_58 = buffer.data(lsp0 + 58);
    const auto *lsp0_60 = buffer.data(lsp0 + 60);
    const auto *lsp0_61 = buffer.data(lsp0 + 61);
    const auto *lsp0_62 = buffer.data(lsp0 + 62);
    const auto *lsp0_63 = buffer.data(lsp0 + 63);
    const auto *lsp0_64 = buffer.data(lsp0 + 64);
    const auto *lsp0_65 = buffer.data(lsp0 + 65);
    const auto *lsp0_68 = buffer.data(lsp0 + 68);
    const auto *lsp0_69 = buffer.data(lsp0 + 69);
    const auto *lsp0_70 = buffer.data(lsp0 + 70);
    const auto *lsp0_71 = buffer.data(lsp0 + 71);
    const auto *lsp0_72 = buffer.data(lsp0 + 72);
    const auto *lsp0_73 = buffer.data(lsp0 + 73);
    const auto *lsp0_74 = buffer.data(lsp0 + 74);
    const auto *lsp0_75 = buffer.data(lsp0 + 75);

    const auto *lsp1_40 = buffer.data(lsp1 + 40);
    const auto *lsp1_42 = buffer.data(lsp1 + 42);
    const auto *lsp1_43 = buffer.data(lsp1 + 43);
    const auto *lsp1_44 = buffer.data(lsp1 + 44);
    const auto *lsp1_45 = buffer.data(lsp1 + 45);
    const auto *lsp1_46 = buffer.data(lsp1 + 46);
    const auto *lsp1_47 = buffer.data(lsp1 + 47);
    const auto *lsp1_50 = buffer.data(lsp1 + 50);
    const auto *lsp1_51 = buffer.data(lsp1 + 51);
    const auto *lsp1_52 = buffer.data(lsp1 + 52);
    const auto *lsp1_53 = buffer.data(lsp1 + 53);
    const auto *lsp1_54 = buffer.data(lsp1 + 54);
    const auto *lsp1_55 = buffer.data(lsp1 + 55);
    const auto *lsp1_56 = buffer.data(lsp1 + 56);
    const auto *lsp1_58 = buffer.data(lsp1 + 58);
    const auto *lsp1_60 = buffer.data(lsp1 + 60);
    const auto *lsp1_61 = buffer.data(lsp1 + 61);
    const auto *lsp1_62 = buffer.data(lsp1 + 62);
    const auto *lsp1_63 = buffer.data(lsp1 + 63);
    const auto *lsp1_64 = buffer.data(lsp1 + 64);
    const auto *lsp1_65 = buffer.data(lsp1 + 65);
    const auto *lsp1_68 = buffer.data(lsp1 + 68);
    const auto *lsp1_69 = buffer.data(lsp1 + 69);
    const auto *lsp1_70 = buffer.data(lsp1 + 70);
    const auto *lsp1_71 = buffer.data(lsp1 + 71);
    const auto *lsp1_72 = buffer.data(lsp1 + 72);
    const auto *lsp1_73 = buffer.data(lsp1 + 73);
    const auto *lsp1_74 = buffer.data(lsp1 + 74);
    const auto *lsp1_75 = buffer.data(lsp1 + 75);

    const auto *lsd_81 = buffer.data(lsd + 81);
    const auto *lsd_82 = buffer.data(lsd + 82);
    const auto *lsd_83 = buffer.data(lsd + 83);
    const auto *lsd_84 = buffer.data(lsd + 84);
    const auto *lsd_86 = buffer.data(lsd + 86);
    const auto *lsd_87 = buffer.data(lsd + 87);
    const auto *lsd_88 = buffer.data(lsd + 88);
    const auto *lsd_89 = buffer.data(lsd + 89);
    const auto *lsd_90 = buffer.data(lsd + 90);
    const auto *lsd_91 = buffer.data(lsd + 91);
    const auto *lsd_93 = buffer.data(lsd + 93);
    const auto *lsd_95 = buffer.data(lsd + 95);
    const auto *lsd_96 = buffer.data(lsd + 96);
    const auto *lsd_99 = buffer.data(lsd + 99);
    const auto *lsd_100 = buffer.data(lsd + 100);
    const auto *lsd_101 = buffer.data(lsd + 101);
    const auto *lsd_102 = buffer.data(lsd + 102);
    const auto *lsd_105 = buffer.data(lsd + 105);
    const auto *lsd_106 = buffer.data(lsd + 106);
    const auto *lsd_107 = buffer.data(lsd + 107);
    const auto *lsd_108 = buffer.data(lsd + 108);
    const auto *lsd_111 = buffer.data(lsd + 111);
    const auto *lsd_112 = buffer.data(lsd + 112);
    const auto *lsd_113 = buffer.data(lsd + 113);
    const auto *lsd_114 = buffer.data(lsd + 114);
    const auto *lsd_117 = buffer.data(lsd + 117);
    const auto *lsd_118 = buffer.data(lsd + 118);
    const auto *lsd_119 = buffer.data(lsd + 119);
    const auto *lsd_120 = buffer.data(lsd + 120);
    const auto *lsd_122 = buffer.data(lsd + 122);
    const auto *lsd_123 = buffer.data(lsd + 123);
    const auto *lsd_124 = buffer.data(lsd + 124);
    const auto *lsd_125 = buffer.data(lsd + 125);
    const auto *lsd_126 = buffer.data(lsd + 126);
    const auto *lsd_127 = buffer.data(lsd + 127);
    const auto *lsd_129 = buffer.data(lsd + 129);
    const auto *lsd_131 = buffer.data(lsd + 131);
    const auto *lsd_132 = buffer.data(lsd + 132);
    const auto *lsd_135 = buffer.data(lsd + 135);
    const auto *lsd_136 = buffer.data(lsd + 136);
    const auto *lsd_137 = buffer.data(lsd + 137);
    const auto *lsd_138 = buffer.data(lsd + 138);
    const auto *lsd_141 = buffer.data(lsd + 141);
    const auto *lsd_142 = buffer.data(lsd + 142);
    const auto *lsd_143 = buffer.data(lsd + 143);
    const auto *lsd_144 = buffer.data(lsd + 144);
    const auto *lsd_147 = buffer.data(lsd + 147);
    const auto *lsd_148 = buffer.data(lsd + 148);
    const auto *lsd_149 = buffer.data(lsd + 149);
    const auto *lsd_150 = buffer.data(lsd + 150);
    const auto *lsd_153 = buffer.data(lsd + 153);
    const auto *lsd_154 = buffer.data(lsd + 154);
    const auto *lsd_155 = buffer.data(lsd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, ksd_57, ksd_81, ksd_82, \
                         ksd_83, lsp0_40, lsp1_40, lsd_81, lsd_82, \
                         lsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_13 * ksd_81[k]
                   + f_3 * pc_x[k] * lsd_81[k];

        t_134[k] = f_13 * ksd_82[k]
                   + f_3 * pc_x[k] * lsd_82[k];

        t_135[k] = f_13 * ksd_83[k]
                   + f_3 * pc_x[k] * lsd_83[k];

        t_136[k] = f_5 * ksd_57[k]
                   + f_1 * lsp0_40[k]
                   - f_2 * lsp1_40[k]
                   + f_3 * pc_y[k] * lsd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, ksf0_99, ksd_51, ksd_59, \
                         ksf1_99, lsd_81, lsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * ksd_51[k]
                   + f_3 * pc_z[k] * lsd_81[k];

        t_138[k] = f_5 * ksd_59[k]
                   + f_3 * pc_y[k] * lsd_83[k];

        t_139[k] = pa_y[k] * ksf0_99[k]
                   - f_4 * pc_y[k] * ksf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, ksd_54, ksd_84, \
                         ksd_87, lsp0_42, lsp1_42, lsd_84, lsd_86, \
                         lsd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * ksd_84[k]
                   + f_1 * lsp0_42[k]
                   - f_2 * lsp1_42[k]
                   + f_3 * pc_x[k] * lsd_84[k];

        t_141[k] = f_3 * pc_y[k] * lsd_84[k];

        t_142[k] = f_13 * ksd_54[k]
                   + f_3 * pc_z[k] * lsd_84[k];

        t_143[k] = f_13 * ksd_87[k]
                   + f_3 * pc_x[k] * lsd_87[k];

        t_144[k] = f_3 * pc_y[k] * lsd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, ksd_89, lsp0_43, lsp0_44, \
                         lsp1_43, lsp1_44, lsd_87, lsd_88, lsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * ksd_89[k]
                   + f_3 * pc_x[k] * lsd_89[k];

        t_146[k] = f_1 * lsp0_43[k]
                   - f_2 * lsp1_43[k]
                   + f_3 * pc_y[k] * lsd_87[k];

        t_147[k] = f_7 * lsp0_44[k]
                   - f_8 * lsp1_44[k]
                   + f_3 * pc_y[k] * lsd_88[k];

        t_148[k] = f_3 * pc_y[k] * lsd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, ksd_59, ksd_60, ksd_90, \
                         lsp0_44, lsp0_45, lsp1_44, lsp1_45, lsd_89, \
                         lsd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_13 * ksd_59[k]
                   + f_1 * lsp0_44[k]
                   - f_2 * lsp1_44[k]
                   + f_3 * pc_z[k] * lsd_89[k];

        t_150[k] = f_12 * ksd_90[k]
                   + f_1 * lsp0_45[k]
                   - f_2 * lsp1_45[k]
                   + f_3 * pc_x[k] * lsd_90[k];

        t_151[k] = f_11 * ksd_60[k]
                   + f_3 * pc_y[k] * lsd_90[k];

        t_152[k] = f_3 * pc_z[k] * lsd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, ksd_63, ksd_93, \
                         ksd_95, lsp0_46, lsp1_46, lsd_91, lsd_93, \
                         lsd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_12 * ksd_93[k]
                   + f_3 * pc_x[k] * lsd_93[k];

        t_154[k] = f_3 * pc_z[k] * lsd_91[k];

        t_155[k] = f_12 * ksd_95[k]
                   + f_3 * pc_x[k] * lsd_95[k];

        t_156[k] = f_11 * ksd_63[k]
                   + f_1 * lsp0_46[k]
                   - f_2 * lsp1_46[k]
                   + f_3 * pc_y[k] * lsd_93[k];

        t_157[k] = f_3 * pc_z[k] * lsd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, ksf0_100, ksd_65, \
                         ksd_66, ksf1_100, lsp0_47, lsp1_47, lsd_95, \
                         lsd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_11 * ksd_65[k]
                   + f_3 * pc_y[k] * lsd_95[k];

        t_159[k] = f_1 * lsp0_47[k]
                   - f_2 * lsp1_47[k]
                   + f_3 * pc_z[k] * lsd_95[k];

        t_160[k] = pa_z[k] * ksf0_100[k]
                   - f_4 * pc_z[k] * ksf1_100[k];

        t_161[k] = f_13 * ksd_66[k]
                   + f_3 * pc_y[k] * lsd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, ksd_60, ksd_99, ksd_100, \
                         ksd_101, lsd_96, lsd_99, lsd_100, lsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * ksd_60[k]
                   + f_3 * pc_z[k] * lsd_96[k];

        t_163[k] = f_12 * ksd_99[k]
                   + f_3 * pc_x[k] * lsd_99[k];

        t_164[k] = f_12 * ksd_100[k]
                   + f_3 * pc_x[k] * lsd_100[k];

        t_165[k] = f_12 * ksd_101[k]
                   + f_3 * pc_x[k] * lsd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, ksf0_106, ksd_63, \
                         ksd_65, ksd_71, ksf1_106, lsp0_50, lsp1_50, lsd_99, \
                         lsd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * ksf0_106[k]
                   - f_4 * pc_z[k] * ksf1_106[k];

        t_167[k] = f_5 * ksd_63[k]
                   + f_3 * pc_z[k] * lsd_99[k];

        t_168[k] = f_13 * ksd_71[k]
                   + f_3 * pc_y[k] * lsd_101[k];

        t_169[k] = f_5 * ksd_65[k]
                   + f_1 * lsp0_50[k]
                   - f_2 * lsp1_50[k]
                   + f_3 * pc_z[k] * lsd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, ksd_66, ksd_72, \
                         ksd_102, ksd_105, lsp0_51, lsp1_51, lsd_102, \
                         lsd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_12 * ksd_102[k]
                   + f_1 * lsp0_51[k]
                   - f_2 * lsp1_51[k]
                   + f_3 * pc_x[k] * lsd_102[k];

        t_171[k] = f_12 * ksd_72[k]
                   + f_3 * pc_y[k] * lsd_102[k];

        t_172[k] = f_10 * ksd_66[k]
                   + f_3 * pc_z[k] * lsd_102[k];

        t_173[k] = f_12 * ksd_105[k]
                   + f_3 * pc_x[k] * lsd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, ksd_69, ksd_75, \
                         ksd_106, ksd_107, lsp0_52, lsp1_52, lsd_105, lsd_106, \
                         lsd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * ksd_106[k]
                   + f_3 * pc_x[k] * lsd_106[k];

        t_175[k] = f_12 * ksd_107[k]
                   + f_3 * pc_x[k] * lsd_107[k];

        t_176[k] = f_12 * ksd_75[k]
                   + f_1 * lsp0_52[k]
                   - f_2 * lsp1_52[k]
                   + f_3 * pc_y[k] * lsd_105[k];

        t_177[k] = f_10 * ksd_69[k]
                   + f_3 * pc_z[k] * lsd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, ksd_71, ksd_77, ksd_108, \
                         lsp0_53, lsp0_54, lsp1_53, lsp1_54, lsd_107, \
                         lsd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * ksd_77[k]
                   + f_3 * pc_y[k] * lsd_107[k];

        t_179[k] = f_10 * ksd_71[k]
                   + f_1 * lsp0_53[k]
                   - f_2 * lsp1_53[k]
                   + f_3 * pc_z[k] * lsd_107[k];

        t_180[k] = f_12 * ksd_108[k]
                   + f_1 * lsp0_54[k]
                   - f_2 * lsp1_54[k]
                   + f_3 * pc_x[k] * lsd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, ksd_72, ksd_78, \
                         ksd_111, ksd_112, lsd_108, lsd_111, lsd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * ksd_78[k]
                   + f_3 * pc_y[k] * lsd_108[k];

        t_182[k] = f_12 * ksd_72[k]
                   + f_3 * pc_z[k] * lsd_108[k];

        t_183[k] = f_12 * ksd_111[k]
                   + f_3 * pc_x[k] * lsd_111[k];

        t_184[k] = f_12 * ksd_112[k]
                   + f_3 * pc_x[k] * lsd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, ksd_75, ksd_81, ksd_83, \
                         ksd_113, lsp0_55, lsp1_55, lsd_111, lsd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_12 * ksd_113[k]
                   + f_3 * pc_x[k] * lsd_113[k];

        t_186[k] = f_10 * ksd_81[k]
                   + f_1 * lsp0_55[k]
                   - f_2 * lsp1_55[k]
                   + f_3 * pc_y[k] * lsd_111[k];

        t_187[k] = f_12 * ksd_75[k]
                   + f_3 * pc_z[k] * lsd_111[k];

        t_188[k] = f_10 * ksd_83[k]
                   + f_3 * pc_y[k] * lsd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, ksf0_140, ksd_77, \
                         ksd_78, ksd_84, ksf1_140, lsp0_56, lsp1_56, lsd_113, \
                         lsd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * ksd_77[k]
                   + f_1 * lsp0_56[k]
                   - f_2 * lsp1_56[k]
                   + f_3 * pc_z[k] * lsd_113[k];

        t_190[k] = pa_y[k] * ksf0_140[k]
                   - f_4 * pc_y[k] * ksf1_140[k];

        t_191[k] = f_5 * ksd_84[k]
                   + f_3 * pc_y[k] * lsd_114[k];

        t_192[k] = f_13 * ksd_78[k]
                   + f_3 * pc_z[k] * lsd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, ksd_87, ksd_117, ksd_118, \
                         ksd_119, lsp0_58, lsp1_58, lsd_117, lsd_118, \
                         lsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_12 * ksd_117[k]
                   + f_3 * pc_x[k] * lsd_117[k];

        t_194[k] = f_12 * ksd_118[k]
                   + f_3 * pc_x[k] * lsd_118[k];

        t_195[k] = f_12 * ksd_119[k]
                   + f_3 * pc_x[k] * lsd_119[k];

        t_196[k] = f_5 * ksd_87[k]
                   + f_1 * lsp0_58[k]
                   - f_2 * lsp1_58[k]
                   + f_3 * pc_y[k] * lsd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, ksf0_149, ksd_81, ksd_89, \
                         ksf1_149, lsd_117, lsd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_13 * ksd_81[k]
                   + f_3 * pc_z[k] * lsd_117[k];

        t_198[k] = f_5 * ksd_89[k]
                   + f_3 * pc_y[k] * lsd_119[k];

        t_199[k] = pa_y[k] * ksf0_149[k]
                   - f_4 * pc_y[k] * ksf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, ksd_84, ksd_120, \
                         ksd_123, lsp0_60, lsp1_60, lsd_120, lsd_122, \
                         lsd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_12 * ksd_120[k]
                   + f_1 * lsp0_60[k]
                   - f_2 * lsp1_60[k]
                   + f_3 * pc_x[k] * lsd_120[k];

        t_201[k] = f_3 * pc_y[k] * lsd_120[k];

        t_202[k] = f_11 * ksd_84[k]
                   + f_3 * pc_z[k] * lsd_120[k];

        t_203[k] = f_12 * ksd_123[k]
                   + f_3 * pc_x[k] * lsd_123[k];

        t_204[k] = f_3 * pc_y[k] * lsd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, ksd_125, lsp0_61, lsp0_62, \
                         lsp1_61, lsp1_62, lsd_123, lsd_124, lsd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_12 * ksd_125[k]
                   + f_3 * pc_x[k] * lsd_125[k];

        t_206[k] = f_1 * lsp0_61[k]
                   - f_2 * lsp1_61[k]
                   + f_3 * pc_y[k] * lsd_123[k];

        t_207[k] = f_7 * lsp0_62[k]
                   - f_8 * lsp1_62[k]
                   + f_3 * pc_y[k] * lsd_124[k];

        t_208[k] = f_3 * pc_y[k] * lsd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, ksd_89, ksd_90, \
                         ksd_126, lsp0_62, lsp0_63, lsp1_62, lsp1_63, lsd_125, \
                         lsd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_11 * ksd_89[k]
                   + f_1 * lsp0_62[k]
                   - f_2 * lsp1_62[k]
                   + f_3 * pc_z[k] * lsd_125[k];

        t_210[k] = f_10 * ksd_126[k]
                   + f_1 * lsp0_63[k]
                   - f_2 * lsp1_63[k]
                   + f_3 * pc_x[k] * lsd_126[k];

        t_211[k] = f_9 * ksd_90[k]
                   + f_3 * pc_y[k] * lsd_126[k];

        t_212[k] = f_3 * pc_z[k] * lsd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pc_x, pc_y, pc_z, ksd_93, ksd_129, \
                         ksd_131, lsp0_64, lsp1_64, lsd_127, lsd_129, \
                         lsd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_10 * ksd_129[k]
                   + f_3 * pc_x[k] * lsd_129[k];

        t_214[k] = f_3 * pc_z[k] * lsd_127[k];

        t_215[k] = f_10 * ksd_131[k]
                   + f_3 * pc_x[k] * lsd_131[k];

        t_216[k] = f_9 * ksd_93[k]
                   + f_1 * lsp0_64[k]
                   - f_2 * lsp1_64[k]
                   + f_3 * pc_y[k] * lsd_129[k];

        t_217[k] = f_3 * pc_z[k] * lsd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pc_y, pc_z, ksf0_150, ksd_95, \
                         ksd_96, ksf1_150, lsp0_65, lsp1_65, lsd_131, \
                         lsd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_9 * ksd_95[k]
                   + f_3 * pc_y[k] * lsd_131[k];

        t_219[k] = f_1 * lsp0_65[k]
                   - f_2 * lsp1_65[k]
                   + f_3 * pc_z[k] * lsd_131[k];

        t_220[k] = pa_z[k] * ksf0_150[k]
                   - f_4 * pc_z[k] * ksf1_150[k];

        t_221[k] = f_11 * ksd_96[k]
                   + f_3 * pc_y[k] * lsd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, ksd_90, ksd_135, ksd_136, \
                         ksd_137, lsd_132, lsd_135, lsd_136, lsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_5 * ksd_90[k]
                   + f_3 * pc_z[k] * lsd_132[k];

        t_223[k] = f_10 * ksd_135[k]
                   + f_3 * pc_x[k] * lsd_135[k];

        t_224[k] = f_10 * ksd_136[k]
                   + f_3 * pc_x[k] * lsd_136[k];

        t_225[k] = f_10 * ksd_137[k]
                   + f_3 * pc_x[k] * lsd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pc_y, pc_z, ksf0_156, ksd_93, \
                         ksd_95, ksd_101, ksf1_156, lsp0_68, lsp1_68, lsd_135, \
                         lsd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * ksf0_156[k]
                   - f_4 * pc_z[k] * ksf1_156[k];

        t_227[k] = f_5 * ksd_93[k]
                   + f_3 * pc_z[k] * lsd_135[k];

        t_228[k] = f_11 * ksd_101[k]
                   + f_3 * pc_y[k] * lsd_137[k];

        t_229[k] = f_5 * ksd_95[k]
                   + f_1 * lsp0_68[k]
                   - f_2 * lsp1_68[k]
                   + f_3 * pc_z[k] * lsd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, ksd_96, ksd_102, \
                         ksd_138, ksd_141, lsp0_69, lsp1_69, lsd_138, \
                         lsd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_10 * ksd_138[k]
                   + f_1 * lsp0_69[k]
                   - f_2 * lsp1_69[k]
                   + f_3 * pc_x[k] * lsd_138[k];

        t_231[k] = f_13 * ksd_102[k]
                   + f_3 * pc_y[k] * lsd_138[k];

        t_232[k] = f_10 * ksd_96[k]
                   + f_3 * pc_z[k] * lsd_138[k];

        t_233[k] = f_10 * ksd_141[k]
                   + f_3 * pc_x[k] * lsd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, ksd_99, ksd_105, \
                         ksd_142, ksd_143, lsp0_70, lsp1_70, lsd_141, lsd_142, \
                         lsd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_10 * ksd_142[k]
                   + f_3 * pc_x[k] * lsd_142[k];

        t_235[k] = f_10 * ksd_143[k]
                   + f_3 * pc_x[k] * lsd_143[k];

        t_236[k] = f_13 * ksd_105[k]
                   + f_1 * lsp0_70[k]
                   - f_2 * lsp1_70[k]
                   + f_3 * pc_y[k] * lsd_141[k];

        t_237[k] = f_10 * ksd_99[k]
                   + f_3 * pc_z[k] * lsd_141[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pc_x, pc_y, pc_z, ksd_101, ksd_107, ksd_144, \
                         lsp0_71, lsp0_72, lsp1_71, lsp1_72, lsd_143, \
                         lsd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_13 * ksd_107[k]
                   + f_3 * pc_y[k] * lsd_143[k];

        t_239[k] = f_10 * ksd_101[k]
                   + f_1 * lsp0_71[k]
                   - f_2 * lsp1_71[k]
                   + f_3 * pc_z[k] * lsd_143[k];

        t_240[k] = f_10 * ksd_144[k]
                   + f_1 * lsp0_72[k]
                   - f_2 * lsp1_72[k]
                   + f_3 * pc_x[k] * lsd_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, ksd_102, ksd_108, \
                         ksd_147, ksd_148, lsd_144, lsd_147, lsd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * ksd_108[k]
                   + f_3 * pc_y[k] * lsd_144[k];

        t_242[k] = f_12 * ksd_102[k]
                   + f_3 * pc_z[k] * lsd_144[k];

        t_243[k] = f_10 * ksd_147[k]
                   + f_3 * pc_x[k] * lsd_147[k];

        t_244[k] = f_10 * ksd_148[k]
                   + f_3 * pc_x[k] * lsd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, ksd_105, ksd_111, \
                         ksd_113, ksd_149, lsp0_73, lsp1_73, lsd_147, \
                         lsd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_10 * ksd_149[k]
                   + f_3 * pc_x[k] * lsd_149[k];

        t_246[k] = f_12 * ksd_111[k]
                   + f_1 * lsp0_73[k]
                   - f_2 * lsp1_73[k]
                   + f_3 * pc_y[k] * lsd_147[k];

        t_247[k] = f_12 * ksd_105[k]
                   + f_3 * pc_z[k] * lsd_147[k];

        t_248[k] = f_12 * ksd_113[k]
                   + f_3 * pc_y[k] * lsd_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, ksd_107, ksd_114, ksd_150, \
                         lsp0_74, lsp0_75, lsp1_74, lsp1_75, lsd_149, \
                         lsd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * ksd_107[k]
                   + f_1 * lsp0_74[k]
                   - f_2 * lsp1_74[k]
                   + f_3 * pc_z[k] * lsd_149[k];

        t_250[k] = f_10 * ksd_150[k]
                   + f_1 * lsp0_75[k]
                   - f_2 * lsp1_75[k]
                   + f_3 * pc_x[k] * lsd_150[k];

        t_251[k] = f_10 * ksd_114[k]
                   + f_3 * pc_y[k] * lsd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, ksd_108, ksd_153, ksd_154, \
                         ksd_155, lsd_150, lsd_153, lsd_154, lsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_13 * ksd_108[k]
                   + f_3 * pc_z[k] * lsd_150[k];

        t_253[k] = f_10 * ksd_153[k]
                   + f_3 * pc_x[k] * lsd_153[k];

        t_254[k] = f_10 * ksd_154[k]
                   + f_3 * pc_x[k] * lsd_154[k];

        t_255[k] = f_10 * ksd_155[k]
                   + f_3 * pc_x[k] * lsd_155[k];
    }
}

static auto
compute_prim_lsf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksf0,
                                                          const size_t ksd, const size_t ksf1,
                                                          const size_t lsp0, const size_t lsp1,
                                                          const size_t lsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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
    auto *t_380 = buffer.data(target + 380);
    auto *t_381 = buffer.data(target + 381);
    auto *t_382 = buffer.data(target + 382);
    auto *t_383 = buffer.data(target + 383);
    auto *t_384 = buffer.data(target + 384);
    auto *t_385 = buffer.data(target + 385);
    auto *t_386 = buffer.data(target + 386);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksf0_200 = buffer.data(ksf0 + 200);
    const auto *ksf0_209 = buffer.data(ksf0 + 209);
    const auto *ksf0_210 = buffer.data(ksf0 + 210);
    const auto *ksf0_270 = buffer.data(ksf0 + 270);
    const auto *ksf0_280 = buffer.data(ksf0 + 280);
    const auto *ksf0_281 = buffer.data(ksf0 + 281);
    const auto *ksf0_286 = buffer.data(ksf0 + 286);
    const auto *ksf0_289 = buffer.data(ksf0 + 289);
    const auto *ksf0_296 = buffer.data(ksf0 + 296);
    const auto *ksf0_299 = buffer.data(ksf0 + 299);
    const auto *ksf0_300 = buffer.data(ksf0 + 300);
    const auto *ksf0_306 = buffer.data(ksf0 + 306);
    const auto *ksf0_309 = buffer.data(ksf0 + 309);
    const auto *ksf0_310 = buffer.data(ksf0 + 310);
    const auto *ksf0_316 = buffer.data(ksf0 + 316);
    const auto *ksf0_319 = buffer.data(ksf0 + 319);
    const auto *ksf0_320 = buffer.data(ksf0 + 320);
    const auto *ksf0_326 = buffer.data(ksf0 + 326);
    const auto *ksf0_329 = buffer.data(ksf0 + 329);
    const auto *ksf0_330 = buffer.data(ksf0 + 330);
    const auto *ksf0_336 = buffer.data(ksf0 + 336);
    const auto *ksf0_339 = buffer.data(ksf0 + 339);
    const auto *ksf0_346 = buffer.data(ksf0 + 346);
    const auto *ksf0_349 = buffer.data(ksf0 + 349);
    const auto *ksf0_350 = buffer.data(ksf0 + 350);
    const auto *ksf0_356 = buffer.data(ksf0 + 356);
    const auto *ksf0_357 = buffer.data(ksf0 + 357);
    const auto *ksf0_359 = buffer.data(ksf0 + 359);

    const auto *ksd_111 = buffer.data(ksd + 111);
    const auto *ksd_113 = buffer.data(ksd + 113);
    const auto *ksd_114 = buffer.data(ksd + 114);
    const auto *ksd_117 = buffer.data(ksd + 117);
    const auto *ksd_119 = buffer.data(ksd + 119);
    const auto *ksd_120 = buffer.data(ksd + 120);
    const auto *ksd_123 = buffer.data(ksd + 123);
    const auto *ksd_125 = buffer.data(ksd + 125);
    const auto *ksd_126 = buffer.data(ksd + 126);
    const auto *ksd_129 = buffer.data(ksd + 129);
    const auto *ksd_131 = buffer.data(ksd + 131);
    const auto *ksd_132 = buffer.data(ksd + 132);
    const auto *ksd_135 = buffer.data(ksd + 135);
    const auto *ksd_137 = buffer.data(ksd + 137);
    const auto *ksd_138 = buffer.data(ksd + 138);
    const auto *ksd_141 = buffer.data(ksd + 141);
    const auto *ksd_143 = buffer.data(ksd + 143);
    const auto *ksd_144 = buffer.data(ksd + 144);
    const auto *ksd_147 = buffer.data(ksd + 147);
    const auto *ksd_149 = buffer.data(ksd + 149);
    const auto *ksd_150 = buffer.data(ksd + 150);
    const auto *ksd_153 = buffer.data(ksd + 153);
    const auto *ksd_155 = buffer.data(ksd + 155);
    const auto *ksd_156 = buffer.data(ksd + 156);
    const auto *ksd_159 = buffer.data(ksd + 159);
    const auto *ksd_160 = buffer.data(ksd + 160);
    const auto *ksd_161 = buffer.data(ksd + 161);
    const auto *ksd_162 = buffer.data(ksd + 162);
    const auto *ksd_165 = buffer.data(ksd + 165);
    const auto *ksd_167 = buffer.data(ksd + 167);
    const auto *ksd_168 = buffer.data(ksd + 168);
    const auto *ksd_171 = buffer.data(ksd + 171);
    const auto *ksd_173 = buffer.data(ksd + 173);
    const auto *ksd_177 = buffer.data(ksd + 177);
    const auto *ksd_178 = buffer.data(ksd + 178);
    const auto *ksd_179 = buffer.data(ksd + 179);
    const auto *ksd_180 = buffer.data(ksd + 180);
    const auto *ksd_183 = buffer.data(ksd + 183);
    const auto *ksd_184 = buffer.data(ksd + 184);
    const auto *ksd_185 = buffer.data(ksd + 185);
    const auto *ksd_186 = buffer.data(ksd + 186);
    const auto *ksd_189 = buffer.data(ksd + 189);
    const auto *ksd_190 = buffer.data(ksd + 190);
    const auto *ksd_191 = buffer.data(ksd + 191);
    const auto *ksd_192 = buffer.data(ksd + 192);
    const auto *ksd_195 = buffer.data(ksd + 195);
    const auto *ksd_196 = buffer.data(ksd + 196);
    const auto *ksd_197 = buffer.data(ksd + 197);
    const auto *ksd_198 = buffer.data(ksd + 198);
    const auto *ksd_201 = buffer.data(ksd + 201);
    const auto *ksd_202 = buffer.data(ksd + 202);
    const auto *ksd_203 = buffer.data(ksd + 203);
    const auto *ksd_207 = buffer.data(ksd + 207);
    const auto *ksd_208 = buffer.data(ksd + 208);
    const auto *ksd_209 = buffer.data(ksd + 209);
    const auto *ksd_210 = buffer.data(ksd + 210);
    const auto *ksd_213 = buffer.data(ksd + 213);
    const auto *ksd_215 = buffer.data(ksd + 215);

    const auto *ksf1_200 = buffer.data(ksf1 + 200);
    const auto *ksf1_209 = buffer.data(ksf1 + 209);
    const auto *ksf1_210 = buffer.data(ksf1 + 210);
    const auto *ksf1_270 = buffer.data(ksf1 + 270);
    const auto *ksf1_280 = buffer.data(ksf1 + 280);
    const auto *ksf1_281 = buffer.data(ksf1 + 281);
    const auto *ksf1_286 = buffer.data(ksf1 + 286);
    const auto *ksf1_289 = buffer.data(ksf1 + 289);
    const auto *ksf1_296 = buffer.data(ksf1 + 296);
    const auto *ksf1_299 = buffer.data(ksf1 + 299);
    const auto *ksf1_300 = buffer.data(ksf1 + 300);
    const auto *ksf1_306 = buffer.data(ksf1 + 306);
    const auto *ksf1_309 = buffer.data(ksf1 + 309);
    const auto *ksf1_310 = buffer.data(ksf1 + 310);
    const auto *ksf1_316 = buffer.data(ksf1 + 316);
    const auto *ksf1_319 = buffer.data(ksf1 + 319);
    const auto *ksf1_320 = buffer.data(ksf1 + 320);
    const auto *ksf1_326 = buffer.data(ksf1 + 326);
    const auto *ksf1_329 = buffer.data(ksf1 + 329);
    const auto *ksf1_330 = buffer.data(ksf1 + 330);
    const auto *ksf1_336 = buffer.data(ksf1 + 336);
    const auto *ksf1_339 = buffer.data(ksf1 + 339);
    const auto *ksf1_346 = buffer.data(ksf1 + 346);
    const auto *ksf1_349 = buffer.data(ksf1 + 349);
    const auto *ksf1_350 = buffer.data(ksf1 + 350);
    const auto *ksf1_356 = buffer.data(ksf1 + 356);
    const auto *ksf1_357 = buffer.data(ksf1 + 357);
    const auto *ksf1_359 = buffer.data(ksf1 + 359);

    const auto *lsp0_76 = buffer.data(lsp0 + 76);
    const auto *lsp0_77 = buffer.data(lsp0 + 77);
    const auto *lsp0_79 = buffer.data(lsp0 + 79);
    const auto *lsp0_81 = buffer.data(lsp0 + 81);
    const auto *lsp0_82 = buffer.data(lsp0 + 82);
    const auto *lsp0_83 = buffer.data(lsp0 + 83);
    const auto *lsp0_108 = buffer.data(lsp0 + 108);
    const auto *lsp0_109 = buffer.data(lsp0 + 109);
    const auto *lsp0_110 = buffer.data(lsp0 + 110);
    const auto *lsp0_113 = buffer.data(lsp0 + 113);
    const auto *lsp0_114 = buffer.data(lsp0 + 114);
    const auto *lsp0_115 = buffer.data(lsp0 + 115);
    const auto *lsp0_116 = buffer.data(lsp0 + 116);

    const auto *lsp1_76 = buffer.data(lsp1 + 76);
    const auto *lsp1_77 = buffer.data(lsp1 + 77);
    const auto *lsp1_79 = buffer.data(lsp1 + 79);
    const auto *lsp1_81 = buffer.data(lsp1 + 81);
    const auto *lsp1_82 = buffer.data(lsp1 + 82);
    const auto *lsp1_83 = buffer.data(lsp1 + 83);
    const auto *lsp1_108 = buffer.data(lsp1 + 108);
    const auto *lsp1_109 = buffer.data(lsp1 + 109);
    const auto *lsp1_110 = buffer.data(lsp1 + 110);
    const auto *lsp1_113 = buffer.data(lsp1 + 113);
    const auto *lsp1_114 = buffer.data(lsp1 + 114);
    const auto *lsp1_115 = buffer.data(lsp1 + 115);
    const auto *lsp1_116 = buffer.data(lsp1 + 116);

    const auto *lsd_153 = buffer.data(lsd + 153);
    const auto *lsd_155 = buffer.data(lsd + 155);
    const auto *lsd_156 = buffer.data(lsd + 156);
    const auto *lsd_159 = buffer.data(lsd + 159);
    const auto *lsd_160 = buffer.data(lsd + 160);
    const auto *lsd_161 = buffer.data(lsd + 161);
    const auto *lsd_162 = buffer.data(lsd + 162);
    const auto *lsd_164 = buffer.data(lsd + 164);
    const auto *lsd_165 = buffer.data(lsd + 165);
    const auto *lsd_166 = buffer.data(lsd + 166);
    const auto *lsd_167 = buffer.data(lsd + 167);
    const auto *lsd_168 = buffer.data(lsd + 168);
    const auto *lsd_169 = buffer.data(lsd + 169);
    const auto *lsd_171 = buffer.data(lsd + 171);
    const auto *lsd_173 = buffer.data(lsd + 173);
    const auto *lsd_174 = buffer.data(lsd + 174);
    const auto *lsd_177 = buffer.data(lsd + 177);
    const auto *lsd_178 = buffer.data(lsd + 178);
    const auto *lsd_179 = buffer.data(lsd + 179);
    const auto *lsd_180 = buffer.data(lsd + 180);
    const auto *lsd_183 = buffer.data(lsd + 183);
    const auto *lsd_184 = buffer.data(lsd + 184);
    const auto *lsd_185 = buffer.data(lsd + 185);
    const auto *lsd_186 = buffer.data(lsd + 186);
    const auto *lsd_189 = buffer.data(lsd + 189);
    const auto *lsd_190 = buffer.data(lsd + 190);
    const auto *lsd_191 = buffer.data(lsd + 191);
    const auto *lsd_192 = buffer.data(lsd + 192);
    const auto *lsd_195 = buffer.data(lsd + 195);
    const auto *lsd_196 = buffer.data(lsd + 196);
    const auto *lsd_197 = buffer.data(lsd + 197);
    const auto *lsd_198 = buffer.data(lsd + 198);
    const auto *lsd_201 = buffer.data(lsd + 201);
    const auto *lsd_202 = buffer.data(lsd + 202);
    const auto *lsd_203 = buffer.data(lsd + 203);
    const auto *lsd_204 = buffer.data(lsd + 204);
    const auto *lsd_207 = buffer.data(lsd + 207);
    const auto *lsd_208 = buffer.data(lsd + 208);
    const auto *lsd_209 = buffer.data(lsd + 209);
    const auto *lsd_210 = buffer.data(lsd + 210);
    const auto *lsd_212 = buffer.data(lsd + 212);
    const auto *lsd_213 = buffer.data(lsd + 213);
    const auto *lsd_215 = buffer.data(lsd + 215);
    const auto *lsd_216 = buffer.data(lsd + 216);
    const auto *lsd_217 = buffer.data(lsd + 217);
    const auto *lsd_219 = buffer.data(lsd + 219);
    const auto *lsd_220 = buffer.data(lsd + 220);
    const auto *lsd_221 = buffer.data(lsd + 221);
    const auto *lsd_224 = buffer.data(lsd + 224);
    const auto *lsd_225 = buffer.data(lsd + 225);
    const auto *lsd_226 = buffer.data(lsd + 226);
    const auto *lsd_227 = buffer.data(lsd + 227);
    const auto *lsd_228 = buffer.data(lsd + 228);
    const auto *lsd_229 = buffer.data(lsd + 229);
    const auto *lsd_230 = buffer.data(lsd + 230);
    const auto *lsd_231 = buffer.data(lsd + 231);
    const auto *lsd_232 = buffer.data(lsd + 232);
    const auto *lsd_233 = buffer.data(lsd + 233);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_y, pc_z, ksd_111, ksd_113, ksd_117, \
                         ksd_119, lsp0_76, lsp0_77, lsp1_76, lsp1_77, lsd_153, \
                         lsd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * ksd_117[k]
                   + f_1 * lsp0_76[k]
                   - f_2 * lsp1_76[k]
                   + f_3 * pc_y[k] * lsd_153[k];

        t_257[k] = f_13 * ksd_111[k]
                   + f_3 * pc_z[k] * lsd_153[k];

        t_258[k] = f_10 * ksd_119[k]
                   + f_3 * pc_y[k] * lsd_155[k];

        t_259[k] = f_13 * ksd_113[k]
                   + f_1 * lsp0_77[k]
                   - f_2 * lsp1_77[k]
                   + f_3 * pc_z[k] * lsd_155[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, ksf0_200, \
                         ksd_114, ksd_120, ksd_159, ksf1_200, lsd_156, \
                         lsd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * ksf0_200[k]
                   - f_4 * pc_y[k] * ksf1_200[k];

        t_261[k] = f_5 * ksd_120[k]
                   + f_3 * pc_y[k] * lsd_156[k];

        t_262[k] = f_11 * ksd_114[k]
                   + f_3 * pc_z[k] * lsd_156[k];

        t_263[k] = f_10 * ksd_159[k]
                   + f_3 * pc_x[k] * lsd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, pc_z, ksd_117, ksd_123, \
                         ksd_160, ksd_161, lsp0_79, lsp1_79, lsd_159, lsd_160, \
                         lsd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_10 * ksd_160[k]
                   + f_3 * pc_x[k] * lsd_160[k];

        t_265[k] = f_10 * ksd_161[k]
                   + f_3 * pc_x[k] * lsd_161[k];

        t_266[k] = f_5 * ksd_123[k]
                   + f_1 * lsp0_79[k]
                   - f_2 * lsp1_79[k]
                   + f_3 * pc_y[k] * lsd_159[k];

        t_267[k] = f_11 * ksd_117[k]
                   + f_3 * pc_z[k] * lsd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pc_x, pc_y, ksf0_209, ksd_125, \
                         ksd_162, ksf1_209, lsp0_81, lsp1_81, lsd_161, \
                         lsd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * ksd_125[k]
                   + f_3 * pc_y[k] * lsd_161[k];

        t_269[k] = pa_y[k] * ksf0_209[k]
                   - f_4 * pc_y[k] * ksf1_209[k];

        t_270[k] = f_10 * ksd_162[k]
                   + f_1 * lsp0_81[k]
                   - f_2 * lsp1_81[k]
                   + f_3 * pc_x[k] * lsd_162[k];

        t_271[k] = f_3 * pc_y[k] * lsd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, ksd_120, ksd_165, \
                         ksd_167, lsd_162, lsd_164, lsd_165, lsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_9 * ksd_120[k]
                   + f_3 * pc_z[k] * lsd_162[k];

        t_273[k] = f_10 * ksd_165[k]
                   + f_3 * pc_x[k] * lsd_165[k];

        t_274[k] = f_3 * pc_y[k] * lsd_164[k];

        t_275[k] = f_10 * ksd_167[k]
                   + f_3 * pc_x[k] * lsd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, ksd_125, lsp0_82, lsp0_83, \
                         lsp1_82, lsp1_83, lsd_165, lsd_166, lsd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * lsp0_82[k]
                   - f_2 * lsp1_82[k]
                   + f_3 * pc_y[k] * lsd_165[k];

        t_277[k] = f_7 * lsp0_83[k]
                   - f_8 * lsp1_83[k]
                   + f_3 * pc_y[k] * lsd_166[k];

        t_278[k] = f_3 * pc_y[k] * lsd_167[k];

        t_279[k] = f_9 * ksd_125[k]
                   + f_1 * lsp0_83[k]
                   - f_2 * lsp1_83[k]
                   + f_3 * pc_z[k] * lsd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_x, pc_x, pc_y, pc_z, ksf0_280, \
                         ksd_126, ksd_168, ksd_171, ksf1_280, lsd_168, \
                         lsd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_x[k] * ksf0_280[k]
                   + f_12 * ksd_168[k]
                   - f_4 * pc_x[k] * ksf1_280[k];

        t_281[k] = f_6 * ksd_126[k]
                   + f_3 * pc_y[k] * lsd_168[k];

        t_282[k] = f_3 * pc_z[k] * lsd_168[k];

        t_283[k] = f_5 * ksd_171[k]
                   + f_3 * pc_x[k] * lsd_171[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, t_288, pa_x, pc_x, pc_y, pc_z, ksf0_286, \
                         ksd_131, ksd_173, ksf1_286, lsd_169, lsd_171, \
                         lsd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * lsd_169[k];

        t_285[k] = f_5 * ksd_173[k]
                   + f_3 * pc_x[k] * lsd_173[k];

        t_286[k] = pa_x[k] * ksf0_286[k]
                   - f_4 * pc_x[k] * ksf1_286[k];

        t_287[k] = f_3 * pc_z[k] * lsd_171[k];

        t_288[k] = f_6 * ksd_131[k]
                   + f_3 * pc_y[k] * lsd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_x, pa_z, pc_x, pc_y, pc_z, ksf0_210, \
                         ksf0_289, ksd_126, ksd_132, ksf1_210, ksf1_289, \
                         lsd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = pa_x[k] * ksf0_289[k]
                   - f_4 * pc_x[k] * ksf1_289[k];

        t_290[k] = pa_z[k] * ksf0_210[k]
                   - f_4 * pc_z[k] * ksf1_210[k];

        t_291[k] = f_9 * ksd_132[k]
                   + f_3 * pc_y[k] * lsd_174[k];

        t_292[k] = f_5 * ksd_126[k]
                   + f_3 * pc_z[k] * lsd_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_x, pc_x, ksf0_296, ksd_177, ksd_178, \
                         ksd_179, ksf1_296, lsd_177, lsd_178, lsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_5 * ksd_177[k]
                   + f_3 * pc_x[k] * lsd_177[k];

        t_294[k] = f_5 * ksd_178[k]
                   + f_3 * pc_x[k] * lsd_178[k];

        t_295[k] = f_5 * ksd_179[k]
                   + f_3 * pc_x[k] * lsd_179[k];

        t_296[k] = pa_x[k] * ksf0_296[k]
                   - f_4 * pc_x[k] * ksf1_296[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_x, pc_x, pc_y, pc_z, ksf0_299, ksd_129, \
                         ksd_137, ksf1_299, lsd_177, lsd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * ksd_129[k]
                   + f_3 * pc_z[k] * lsd_177[k];

        t_298[k] = f_9 * ksd_137[k]
                   + f_3 * pc_y[k] * lsd_179[k];

        t_299[k] = pa_x[k] * ksf0_299[k]
                   - f_4 * pc_x[k] * ksf1_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pc_x, pc_y, pc_z, ksf0_300, \
                         ksd_132, ksd_138, ksd_180, ksd_183, ksf1_300, lsd_180, \
                         lsd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = pa_x[k] * ksf0_300[k]
                   + f_12 * ksd_180[k]
                   - f_4 * pc_x[k] * ksf1_300[k];

        t_301[k] = f_11 * ksd_138[k]
                   + f_3 * pc_y[k] * lsd_180[k];

        t_302[k] = f_10 * ksd_132[k]
                   + f_3 * pc_z[k] * lsd_180[k];

        t_303[k] = f_5 * ksd_183[k]
                   + f_3 * pc_x[k] * lsd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_x, pc_x, pc_z, ksf0_306, ksd_135, \
                         ksd_184, ksd_185, ksf1_306, lsd_183, lsd_184, \
                         lsd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_5 * ksd_184[k]
                   + f_3 * pc_x[k] * lsd_184[k];

        t_305[k] = f_5 * ksd_185[k]
                   + f_3 * pc_x[k] * lsd_185[k];

        t_306[k] = pa_x[k] * ksf0_306[k]
                   - f_4 * pc_x[k] * ksf1_306[k];

        t_307[k] = f_10 * ksd_135[k]
                   + f_3 * pc_z[k] * lsd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_x, pc_x, pc_y, ksf0_309, ksf0_310, \
                         ksd_143, ksd_144, ksd_186, ksf1_309, ksf1_310, lsd_185, \
                         lsd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_11 * ksd_143[k]
                   + f_3 * pc_y[k] * lsd_185[k];

        t_309[k] = pa_x[k] * ksf0_309[k]
                   - f_4 * pc_x[k] * ksf1_309[k];

        t_310[k] = pa_x[k] * ksf0_310[k]
                   + f_12 * ksd_186[k]
                   - f_4 * pc_x[k] * ksf1_310[k];

        t_311[k] = f_13 * ksd_144[k]
                   + f_3 * pc_y[k] * lsd_186[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_z, ksd_138, ksd_189, ksd_190, \
                         ksd_191, lsd_186, lsd_189, lsd_190, lsd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_12 * ksd_138[k]
                   + f_3 * pc_z[k] * lsd_186[k];

        t_313[k] = f_5 * ksd_189[k]
                   + f_3 * pc_x[k] * lsd_189[k];

        t_314[k] = f_5 * ksd_190[k]
                   + f_3 * pc_x[k] * lsd_190[k];

        t_315[k] = f_5 * ksd_191[k]
                   + f_3 * pc_x[k] * lsd_191[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_x, pc_x, pc_y, pc_z, ksf0_316, \
                         ksf0_319, ksd_141, ksd_149, ksf1_316, ksf1_319, lsd_189, \
                         lsd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = pa_x[k] * ksf0_316[k]
                   - f_4 * pc_x[k] * ksf1_316[k];

        t_317[k] = f_12 * ksd_141[k]
                   + f_3 * pc_z[k] * lsd_189[k];

        t_318[k] = f_13 * ksd_149[k]
                   + f_3 * pc_y[k] * lsd_191[k];

        t_319[k] = pa_x[k] * ksf0_319[k]
                   - f_4 * pc_x[k] * ksf1_319[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pa_x, pc_x, pc_y, pc_z, ksf0_320, \
                         ksd_144, ksd_150, ksd_192, ksd_195, ksf1_320, lsd_192, \
                         lsd_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pa_x[k] * ksf0_320[k]
                   + f_12 * ksd_192[k]
                   - f_4 * pc_x[k] * ksf1_320[k];

        t_321[k] = f_12 * ksd_150[k]
                   + f_3 * pc_y[k] * lsd_192[k];

        t_322[k] = f_13 * ksd_144[k]
                   + f_3 * pc_z[k] * lsd_192[k];

        t_323[k] = f_5 * ksd_195[k]
                   + f_3 * pc_x[k] * lsd_195[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pa_x, pc_x, pc_z, ksf0_326, ksd_147, \
                         ksd_196, ksd_197, ksf1_326, lsd_195, lsd_196, \
                         lsd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_5 * ksd_196[k]
                   + f_3 * pc_x[k] * lsd_196[k];

        t_325[k] = f_5 * ksd_197[k]
                   + f_3 * pc_x[k] * lsd_197[k];

        t_326[k] = pa_x[k] * ksf0_326[k]
                   - f_4 * pc_x[k] * ksf1_326[k];

        t_327[k] = f_13 * ksd_147[k]
                   + f_3 * pc_z[k] * lsd_195[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_x, pc_x, pc_y, ksf0_329, ksf0_330, \
                         ksd_155, ksd_156, ksd_198, ksf1_329, ksf1_330, lsd_197, \
                         lsd_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_12 * ksd_155[k]
                   + f_3 * pc_y[k] * lsd_197[k];

        t_329[k] = pa_x[k] * ksf0_329[k]
                   - f_4 * pc_x[k] * ksf1_329[k];

        t_330[k] = pa_x[k] * ksf0_330[k]
                   + f_12 * ksd_198[k]
                   - f_4 * pc_x[k] * ksf1_330[k];

        t_331[k] = f_10 * ksd_156[k]
                   + f_3 * pc_y[k] * lsd_198[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pc_x, pc_z, ksd_150, ksd_201, ksd_202, \
                         ksd_203, lsd_198, lsd_201, lsd_202, lsd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_11 * ksd_150[k]
                   + f_3 * pc_z[k] * lsd_198[k];

        t_333[k] = f_5 * ksd_201[k]
                   + f_3 * pc_x[k] * lsd_201[k];

        t_334[k] = f_5 * ksd_202[k]
                   + f_3 * pc_x[k] * lsd_202[k];

        t_335[k] = f_5 * ksd_203[k]
                   + f_3 * pc_x[k] * lsd_203[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_x, pc_x, pc_y, pc_z, ksf0_336, \
                         ksf0_339, ksd_153, ksd_161, ksf1_336, ksf1_339, lsd_201, \
                         lsd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_x[k] * ksf0_336[k]
                   - f_4 * pc_x[k] * ksf1_336[k];

        t_337[k] = f_11 * ksd_153[k]
                   + f_3 * pc_z[k] * lsd_201[k];

        t_338[k] = f_10 * ksd_161[k]
                   + f_3 * pc_y[k] * lsd_203[k];

        t_339[k] = pa_x[k] * ksf0_339[k]
                   - f_4 * pc_x[k] * ksf1_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_y, pc_x, pc_y, pc_z, ksf0_270, \
                         ksd_156, ksd_162, ksd_207, ksf1_270, lsd_204, \
                         lsd_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_y[k] * ksf0_270[k]
                   - f_4 * pc_y[k] * ksf1_270[k];

        t_341[k] = f_5 * ksd_162[k]
                   + f_3 * pc_y[k] * lsd_204[k];

        t_342[k] = f_9 * ksd_156[k]
                   + f_3 * pc_z[k] * lsd_204[k];

        t_343[k] = f_5 * ksd_207[k]
                   + f_3 * pc_x[k] * lsd_207[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_x, pc_x, pc_z, ksf0_346, ksd_159, \
                         ksd_208, ksd_209, ksf1_346, lsd_207, lsd_208, \
                         lsd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_5 * ksd_208[k]
                   + f_3 * pc_x[k] * lsd_208[k];

        t_345[k] = f_5 * ksd_209[k]
                   + f_3 * pc_x[k] * lsd_209[k];

        t_346[k] = pa_x[k] * ksf0_346[k]
                   - f_4 * pc_x[k] * ksf1_346[k];

        t_347[k] = f_9 * ksd_159[k]
                   + f_3 * pc_z[k] * lsd_207[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, t_351, pa_x, pc_x, pc_y, ksf0_349, ksf0_350, \
                         ksd_167, ksd_210, ksf1_349, ksf1_350, lsd_209, \
                         lsd_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_5 * ksd_167[k]
                   + f_3 * pc_y[k] * lsd_209[k];

        t_349[k] = pa_x[k] * ksf0_349[k]
                   - f_4 * pc_x[k] * ksf1_349[k];

        t_350[k] = pa_x[k] * ksf0_350[k]
                   + f_12 * ksd_210[k]
                   - f_4 * pc_x[k] * ksf1_350[k];

        t_351[k] = f_3 * pc_y[k] * lsd_210[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pc_x, pc_y, pc_z, ksd_162, ksd_213, \
                         ksd_215, lsd_210, lsd_212, lsd_213, lsd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_6 * ksd_162[k]
                   + f_3 * pc_z[k] * lsd_210[k];

        t_353[k] = f_5 * ksd_213[k]
                   + f_3 * pc_x[k] * lsd_213[k];

        t_354[k] = f_3 * pc_y[k] * lsd_212[k];

        t_355[k] = f_5 * ksd_215[k]
                   + f_3 * pc_x[k] * lsd_215[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pa_x, pc_x, pc_y, ksf0_356, ksf0_357, \
                         ksf0_359, ksf1_356, ksf1_357, ksf1_359, \
                         lsd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pa_x[k] * ksf0_356[k]
                   - f_4 * pc_x[k] * ksf1_356[k];

        t_357[k] = pa_x[k] * ksf0_357[k]
                   - f_4 * pc_x[k] * ksf1_357[k];

        t_358[k] = f_3 * pc_y[k] * lsd_215[k];

        t_359[k] = pa_x[k] * ksf0_359[k]
                   - f_4 * pc_x[k] * ksf1_359[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, t_363, t_364, pc_x, pc_z, lsp0_108, lsp0_109, \
                         lsp1_108, lsp1_109, lsd_216, lsd_217, lsd_219, \
                         lsd_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * lsp0_108[k]
                   - f_2 * lsp1_108[k]
                   + f_3 * pc_x[k] * lsd_216[k];

        t_361[k] = f_7 * lsp0_109[k]
                   - f_8 * lsp1_109[k]
                   + f_3 * pc_x[k] * lsd_217[k];

        t_362[k] = f_3 * pc_z[k] * lsd_216[k];

        t_363[k] = f_3 * pc_x[k] * lsd_219[k];

        t_364[k] = f_3 * pc_x[k] * lsd_220[k];
    }

#pragma omp simd aligned(t_365, t_366, t_367, t_368, t_369, pc_x, pc_y, pc_z, ksd_171, \
                         ksd_173, lsp0_109, lsp0_110, lsp1_109, lsp1_110, lsd_219, \
                         lsd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = f_3 * pc_x[k] * lsd_221[k];

        t_366[k] = f_0 * ksd_171[k]
                   + f_1 * lsp0_109[k]
                   - f_2 * lsp1_109[k]
                   + f_3 * pc_y[k] * lsd_219[k];

        t_367[k] = f_3 * pc_z[k] * lsd_219[k];

        t_368[k] = f_0 * ksd_173[k]
                   + f_3 * pc_y[k] * lsd_221[k];

        t_369[k] = f_1 * lsp0_110[k]
                   - f_2 * lsp1_110[k]
                   + f_3 * pc_z[k] * lsd_221[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_z, pc_x, pc_z, ksf0_280, ksf0_281, \
                         ksf1_280, ksf1_281, lsp0_113, lsp1_113, lsd_224, \
                         lsd_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_z[k] * ksf0_280[k]
                   - f_4 * pc_z[k] * ksf1_280[k];

        t_371[k] = pa_z[k] * ksf0_281[k]
                   - f_4 * pc_z[k] * ksf1_281[k];

        t_372[k] = f_7 * lsp0_113[k]
                   - f_8 * lsp1_113[k]
                   + f_3 * pc_x[k] * lsd_224[k];

        t_373[k] = f_3 * pc_x[k] * lsd_225[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, t_377, t_378, pa_z, pc_x, pc_y, pc_z, ksf0_286, \
                         ksd_171, ksd_179, ksf1_286, lsd_225, lsd_226, \
                         lsd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_3 * pc_x[k] * lsd_226[k];

        t_375[k] = f_3 * pc_x[k] * lsd_227[k];

        t_376[k] = pa_z[k] * ksf0_286[k]
                   - f_4 * pc_z[k] * ksf1_286[k];

        t_377[k] = f_5 * ksd_171[k]
                   + f_3 * pc_z[k] * lsd_225[k];

        t_378[k] = f_6 * ksd_179[k]
                   + f_3 * pc_y[k] * lsd_227[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pc_x, pc_z, ksd_173, lsp0_113, lsp0_114, \
                         lsp0_115, lsp1_113, lsp1_114, lsp1_115, lsd_227, lsd_228, \
                         lsd_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = f_5 * ksd_173[k]
                   + f_1 * lsp0_113[k]
                   - f_2 * lsp1_113[k]
                   + f_3 * pc_z[k] * lsd_227[k];

        t_380[k] = f_1 * lsp0_114[k]
                   - f_2 * lsp1_114[k]
                   + f_3 * pc_x[k] * lsd_228[k];

        t_381[k] = f_7 * lsp0_115[k]
                   - f_8 * lsp1_115[k]
                   + f_3 * pc_x[k] * lsd_229[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, t_386, pc_x, pc_y, ksd_183, lsp0_115, \
                         lsp0_116, lsp1_115, lsp1_116, lsd_230, lsd_231, lsd_232, \
                         lsd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_7 * lsp0_116[k]
                   - f_8 * lsp1_116[k]
                   + f_3 * pc_x[k] * lsd_230[k];

        t_383[k] = f_3 * pc_x[k] * lsd_231[k];

        t_384[k] = f_3 * pc_x[k] * lsd_232[k];

        t_385[k] = f_3 * pc_x[k] * lsd_233[k];

        t_386[k] = f_9 * ksd_183[k]
                   + f_1 * lsp0_115[k]
                   - f_2 * lsp1_115[k]
                   + f_3 * pc_y[k] * lsd_231[k];
    }
}

static auto
compute_prim_lsf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksf0,
                                                          const size_t ksd, const size_t ksf1,
                                                          const size_t lsp0, const size_t lsp1,
                                                          const size_t lsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 3.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 3.0 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 2.5 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksf0_350 = buffer.data(ksf0 + 350);
    const auto *ksf0_352 = buffer.data(ksf0 + 352);
    const auto *ksf0_356 = buffer.data(ksf0 + 356);
    const auto *ksf0_359 = buffer.data(ksf0 + 359);

    const auto *ksd_177 = buffer.data(ksd + 177);
    const auto *ksd_179 = buffer.data(ksd + 179);
    const auto *ksd_183 = buffer.data(ksd + 183);
    const auto *ksd_185 = buffer.data(ksd + 185);
    const auto *ksd_189 = buffer.data(ksd + 189);
    const auto *ksd_191 = buffer.data(ksd + 191);
    const auto *ksd_195 = buffer.data(ksd + 195);
    const auto *ksd_197 = buffer.data(ksd + 197);
    const auto *ksd_201 = buffer.data(ksd + 201);
    const auto *ksd_203 = buffer.data(ksd + 203);
    const auto *ksd_207 = buffer.data(ksd + 207);
    const auto *ksd_209 = buffer.data(ksd + 209);
    const auto *ksd_213 = buffer.data(ksd + 213);
    const auto *ksd_215 = buffer.data(ksd + 215);

    const auto *ksf1_350 = buffer.data(ksf1 + 350);
    const auto *ksf1_352 = buffer.data(ksf1 + 352);
    const auto *ksf1_356 = buffer.data(ksf1 + 356);
    const auto *ksf1_359 = buffer.data(ksf1 + 359);

    const auto *lsp0_116 = buffer.data(lsp0 + 116);
    const auto *lsp0_117 = buffer.data(lsp0 + 117);
    const auto *lsp0_118 = buffer.data(lsp0 + 118);
    const auto *lsp0_119 = buffer.data(lsp0 + 119);
    const auto *lsp0_120 = buffer.data(lsp0 + 120);
    const auto *lsp0_121 = buffer.data(lsp0 + 121);
    const auto *lsp0_122 = buffer.data(lsp0 + 122);
    const auto *lsp0_123 = buffer.data(lsp0 + 123);
    const auto *lsp0_124 = buffer.data(lsp0 + 124);
    const auto *lsp0_125 = buffer.data(lsp0 + 125);
    const auto *lsp0_126 = buffer.data(lsp0 + 126);
    const auto *lsp0_127 = buffer.data(lsp0 + 127);
    const auto *lsp0_128 = buffer.data(lsp0 + 128);
    const auto *lsp0_130 = buffer.data(lsp0 + 130);
    const auto *lsp0_132 = buffer.data(lsp0 + 132);
    const auto *lsp0_133 = buffer.data(lsp0 + 133);
    const auto *lsp0_134 = buffer.data(lsp0 + 134);

    const auto *lsp1_116 = buffer.data(lsp1 + 116);
    const auto *lsp1_117 = buffer.data(lsp1 + 117);
    const auto *lsp1_118 = buffer.data(lsp1 + 118);
    const auto *lsp1_119 = buffer.data(lsp1 + 119);
    const auto *lsp1_120 = buffer.data(lsp1 + 120);
    const auto *lsp1_121 = buffer.data(lsp1 + 121);
    const auto *lsp1_122 = buffer.data(lsp1 + 122);
    const auto *lsp1_123 = buffer.data(lsp1 + 123);
    const auto *lsp1_124 = buffer.data(lsp1 + 124);
    const auto *lsp1_125 = buffer.data(lsp1 + 125);
    const auto *lsp1_126 = buffer.data(lsp1 + 126);
    const auto *lsp1_127 = buffer.data(lsp1 + 127);
    const auto *lsp1_128 = buffer.data(lsp1 + 128);
    const auto *lsp1_130 = buffer.data(lsp1 + 130);
    const auto *lsp1_132 = buffer.data(lsp1 + 132);
    const auto *lsp1_133 = buffer.data(lsp1 + 133);
    const auto *lsp1_134 = buffer.data(lsp1 + 134);

    const auto *lsd_231 = buffer.data(lsd + 231);
    const auto *lsd_233 = buffer.data(lsd + 233);
    const auto *lsd_234 = buffer.data(lsd + 234);
    const auto *lsd_235 = buffer.data(lsd + 235);
    const auto *lsd_236 = buffer.data(lsd + 236);
    const auto *lsd_237 = buffer.data(lsd + 237);
    const auto *lsd_238 = buffer.data(lsd + 238);
    const auto *lsd_239 = buffer.data(lsd + 239);
    const auto *lsd_240 = buffer.data(lsd + 240);
    const auto *lsd_241 = buffer.data(lsd + 241);
    const auto *lsd_242 = buffer.data(lsd + 242);
    const auto *lsd_243 = buffer.data(lsd + 243);
    const auto *lsd_244 = buffer.data(lsd + 244);
    const auto *lsd_245 = buffer.data(lsd + 245);
    const auto *lsd_246 = buffer.data(lsd + 246);
    const auto *lsd_247 = buffer.data(lsd + 247);
    const auto *lsd_248 = buffer.data(lsd + 248);
    const auto *lsd_249 = buffer.data(lsd + 249);
    const auto *lsd_250 = buffer.data(lsd + 250);
    const auto *lsd_251 = buffer.data(lsd + 251);
    const auto *lsd_252 = buffer.data(lsd + 252);
    const auto *lsd_253 = buffer.data(lsd + 253);
    const auto *lsd_254 = buffer.data(lsd + 254);
    const auto *lsd_255 = buffer.data(lsd + 255);
    const auto *lsd_256 = buffer.data(lsd + 256);
    const auto *lsd_257 = buffer.data(lsd + 257);
    const auto *lsd_259 = buffer.data(lsd + 259);
    const auto *lsd_261 = buffer.data(lsd + 261);
    const auto *lsd_262 = buffer.data(lsd + 262);
    const auto *lsd_263 = buffer.data(lsd + 263);
    const auto *lsd_264 = buffer.data(lsd + 264);
    const auto *lsd_266 = buffer.data(lsd + 266);
    const auto *lsd_267 = buffer.data(lsd + 267);
    const auto *lsd_268 = buffer.data(lsd + 268);
    const auto *lsd_269 = buffer.data(lsd + 269);

#pragma omp simd aligned(t_387, t_388, t_389, pc_y, pc_z, ksd_177, ksd_179, ksd_185, lsp0_116, \
                         lsp1_116, lsd_231, lsd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_10 * ksd_177[k]
                   + f_3 * pc_z[k] * lsd_231[k];

        t_388[k] = f_9 * ksd_185[k]
                   + f_3 * pc_y[k] * lsd_233[k];

        t_389[k] = f_10 * ksd_179[k]
                   + f_1 * lsp0_116[k]
                   - f_2 * lsp1_116[k]
                   + f_3 * pc_z[k] * lsd_233[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pc_x, lsp0_117, lsp0_118, lsp0_119, \
                         lsp1_117, lsp1_118, lsp1_119, lsd_234, lsd_235, lsd_236, \
                         lsd_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_1 * lsp0_117[k]
                   - f_2 * lsp1_117[k]
                   + f_3 * pc_x[k] * lsd_234[k];

        t_391[k] = f_7 * lsp0_118[k]
                   - f_8 * lsp1_118[k]
                   + f_3 * pc_x[k] * lsd_235[k];

        t_392[k] = f_7 * lsp0_119[k]
                   - f_8 * lsp1_119[k]
                   + f_3 * pc_x[k] * lsd_236[k];

        t_393[k] = f_3 * pc_x[k] * lsd_237[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, ksd_183, \
                         ksd_189, ksd_191, lsp0_118, lsp1_118, lsd_237, lsd_238, \
                         lsd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_3 * pc_x[k] * lsd_238[k];

        t_395[k] = f_3 * pc_x[k] * lsd_239[k];

        t_396[k] = f_11 * ksd_189[k]
                   + f_1 * lsp0_118[k]
                   - f_2 * lsp1_118[k]
                   + f_3 * pc_y[k] * lsd_237[k];

        t_397[k] = f_12 * ksd_183[k]
                   + f_3 * pc_z[k] * lsd_237[k];

        t_398[k] = f_11 * ksd_191[k]
                   + f_3 * pc_y[k] * lsd_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_z, ksd_185, lsp0_119, lsp0_120, \
                         lsp0_121, lsp1_119, lsp1_120, lsp1_121, lsd_239, lsd_240, \
                         lsd_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_12 * ksd_185[k]
                   + f_1 * lsp0_119[k]
                   - f_2 * lsp1_119[k]
                   + f_3 * pc_z[k] * lsd_239[k];

        t_400[k] = f_1 * lsp0_120[k]
                   - f_2 * lsp1_120[k]
                   + f_3 * pc_x[k] * lsd_240[k];

        t_401[k] = f_7 * lsp0_121[k]
                   - f_8 * lsp1_121[k]
                   + f_3 * pc_x[k] * lsd_241[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, t_406, pc_x, pc_y, ksd_195, lsp0_121, \
                         lsp0_122, lsp1_121, lsp1_122, lsd_242, lsd_243, lsd_244, \
                         lsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * lsp0_122[k]
                   - f_8 * lsp1_122[k]
                   + f_3 * pc_x[k] * lsd_242[k];

        t_403[k] = f_3 * pc_x[k] * lsd_243[k];

        t_404[k] = f_3 * pc_x[k] * lsd_244[k];

        t_405[k] = f_3 * pc_x[k] * lsd_245[k];

        t_406[k] = f_13 * ksd_195[k]
                   + f_1 * lsp0_121[k]
                   - f_2 * lsp1_121[k]
                   + f_3 * pc_y[k] * lsd_243[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pc_y, pc_z, ksd_189, ksd_191, ksd_197, lsp0_122, \
                         lsp1_122, lsd_243, lsd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_13 * ksd_189[k]
                   + f_3 * pc_z[k] * lsd_243[k];

        t_408[k] = f_13 * ksd_197[k]
                   + f_3 * pc_y[k] * lsd_245[k];

        t_409[k] = f_13 * ksd_191[k]
                   + f_1 * lsp0_122[k]
                   - f_2 * lsp1_122[k]
                   + f_3 * pc_z[k] * lsd_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, lsp0_123, lsp0_124, lsp0_125, \
                         lsp1_123, lsp1_124, lsp1_125, lsd_246, lsd_247, lsd_248, \
                         lsd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_1 * lsp0_123[k]
                   - f_2 * lsp1_123[k]
                   + f_3 * pc_x[k] * lsd_246[k];

        t_411[k] = f_7 * lsp0_124[k]
                   - f_8 * lsp1_124[k]
                   + f_3 * pc_x[k] * lsd_247[k];

        t_412[k] = f_7 * lsp0_125[k]
                   - f_8 * lsp1_125[k]
                   + f_3 * pc_x[k] * lsd_248[k];

        t_413[k] = f_3 * pc_x[k] * lsd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, t_418, pc_x, pc_y, pc_z, ksd_195, \
                         ksd_201, ksd_203, lsp0_124, lsp1_124, lsd_249, lsd_250, \
                         lsd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_3 * pc_x[k] * lsd_250[k];

        t_415[k] = f_3 * pc_x[k] * lsd_251[k];

        t_416[k] = f_12 * ksd_201[k]
                   + f_1 * lsp0_124[k]
                   - f_2 * lsp1_124[k]
                   + f_3 * pc_y[k] * lsd_249[k];

        t_417[k] = f_11 * ksd_195[k]
                   + f_3 * pc_z[k] * lsd_249[k];

        t_418[k] = f_12 * ksd_203[k]
                   + f_3 * pc_y[k] * lsd_251[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, pc_x, pc_z, ksd_197, lsp0_125, lsp0_126, \
                         lsp0_127, lsp1_125, lsp1_126, lsp1_127, lsd_251, lsd_252, \
                         lsd_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_11 * ksd_197[k]
                   + f_1 * lsp0_125[k]
                   - f_2 * lsp1_125[k]
                   + f_3 * pc_z[k] * lsd_251[k];

        t_420[k] = f_1 * lsp0_126[k]
                   - f_2 * lsp1_126[k]
                   + f_3 * pc_x[k] * lsd_252[k];

        t_421[k] = f_7 * lsp0_127[k]
                   - f_8 * lsp1_127[k]
                   + f_3 * pc_x[k] * lsd_253[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, t_426, pc_x, pc_y, ksd_207, lsp0_127, \
                         lsp0_128, lsp1_127, lsp1_128, lsd_254, lsd_255, lsd_256, \
                         lsd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_7 * lsp0_128[k]
                   - f_8 * lsp1_128[k]
                   + f_3 * pc_x[k] * lsd_254[k];

        t_423[k] = f_3 * pc_x[k] * lsd_255[k];

        t_424[k] = f_3 * pc_x[k] * lsd_256[k];

        t_425[k] = f_3 * pc_x[k] * lsd_257[k];

        t_426[k] = f_10 * ksd_207[k]
                   + f_1 * lsp0_127[k]
                   - f_2 * lsp1_127[k]
                   + f_3 * pc_y[k] * lsd_255[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, pa_y, pc_y, pc_z, ksf0_350, ksd_201, \
                         ksd_203, ksd_209, ksf1_350, lsp0_128, lsp1_128, lsd_255, \
                         lsd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_9 * ksd_201[k]
                   + f_3 * pc_z[k] * lsd_255[k];

        t_428[k] = f_10 * ksd_209[k]
                   + f_3 * pc_y[k] * lsd_257[k];

        t_429[k] = f_9 * ksd_203[k]
                   + f_1 * lsp0_128[k]
                   - f_2 * lsp1_128[k]
                   + f_3 * pc_z[k] * lsd_257[k];

        t_430[k] = pa_y[k] * ksf0_350[k]
                   - f_4 * pc_y[k] * ksf1_350[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, t_434, t_435, pa_y, pc_x, pc_y, ksf0_352, \
                         ksf1_352, lsp0_130, lsp1_130, lsd_259, lsd_261, lsd_262, \
                         lsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_7 * lsp0_130[k]
                   - f_8 * lsp1_130[k]
                   + f_3 * pc_x[k] * lsd_259[k];

        t_432[k] = pa_y[k] * ksf0_352[k]
                   - f_4 * pc_y[k] * ksf1_352[k];

        t_433[k] = f_3 * pc_x[k] * lsd_261[k];

        t_434[k] = f_3 * pc_x[k] * lsd_262[k];

        t_435[k] = f_3 * pc_x[k] * lsd_263[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pa_y, pc_y, pc_z, ksf0_356, ksf0_359, \
                         ksd_207, ksd_213, ksd_215, ksf1_356, ksf1_359, lsd_261, \
                         lsd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = pa_y[k] * ksf0_356[k]
                   + f_12 * ksd_213[k]
                   - f_4 * pc_y[k] * ksf1_356[k];

        t_437[k] = f_6 * ksd_207[k]
                   + f_3 * pc_z[k] * lsd_261[k];

        t_438[k] = f_5 * ksd_215[k]
                   + f_3 * pc_y[k] * lsd_263[k];

        t_439[k] = pa_y[k] * ksf0_359[k]
                   - f_4 * pc_y[k] * ksf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pc_y, lsp0_132, lsp0_134, \
                         lsp1_132, lsp1_134, lsd_264, lsd_266, lsd_267, \
                         lsd_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_1 * lsp0_132[k]
                   - f_2 * lsp1_132[k]
                   + f_3 * pc_x[k] * lsd_264[k];

        t_441[k] = f_3 * pc_y[k] * lsd_264[k];

        t_442[k] = f_7 * lsp0_134[k]
                   - f_8 * lsp1_134[k]
                   + f_3 * pc_x[k] * lsd_266[k];

        t_443[k] = f_3 * pc_x[k] * lsd_267[k];

        t_444[k] = f_3 * pc_x[k] * lsd_268[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, t_449, pc_x, pc_y, pc_z, ksd_215, \
                         lsp0_133, lsp0_134, lsp1_133, lsp1_134, lsd_267, lsd_268, \
                         lsd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_3 * pc_x[k] * lsd_269[k];

        t_446[k] = f_1 * lsp0_133[k]
                   - f_2 * lsp1_133[k]
                   + f_3 * pc_y[k] * lsd_267[k];

        t_447[k] = f_7 * lsp0_134[k]
                   - f_8 * lsp1_134[k]
                   + f_3 * pc_y[k] * lsd_268[k];

        t_448[k] = f_3 * pc_y[k] * lsd_269[k];

        t_449[k] = f_0 * ksd_215[k]
                   + f_1 * lsp0_134[k]
                   - f_2 * lsp1_134[k]
                   + f_3 * pc_z[k] * lsd_269[k];
    }
}

auto
compute_prim_lsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksf0, const size_t ksd,
                                                   const size_t ksf1, const size_t lsp0,
                                                   const size_t lsp1, const size_t lsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksf0, ksd,
                                                              ksf1, lsp0, lsp1, lsd, ncols,
                                                              gamma, p, q);

    compute_prim_lsf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksf0, ksd,
                                                              ksf1, lsp0, lsp1, lsd, ncols,
                                                              gamma, p, q);

    compute_prim_lsf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksf0, ksd,
                                                              ksf1, lsp0, lsp1, lsd, ncols,
                                                              gamma, p, q);

    compute_prim_lsf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksf0, ksd,
                                                              ksf1, lsp0, lsp1, lsd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
