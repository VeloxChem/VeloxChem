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


#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_0 = buffer.data(ksi0 + 0);
    const auto *ksi0_3 = buffer.data(ksi0 + 3);
    const auto *ksi0_5 = buffer.data(ksi0 + 5);
    const auto *ksi0_6 = buffer.data(ksi0 + 6);
    const auto *ksi0_9 = buffer.data(ksi0 + 9);
    const auto *ksi0_10 = buffer.data(ksi0 + 10);
    const auto *ksi0_14 = buffer.data(ksi0 + 14);
    const auto *ksi0_21 = buffer.data(ksi0 + 21);
    const auto *ksi0_27 = buffer.data(ksi0 + 27);
    const auto *ksi0_31 = buffer.data(ksi0 + 31);
    const auto *ksi0_34 = buffer.data(ksi0 + 34);
    const auto *ksi0_38 = buffer.data(ksi0 + 38);
    const auto *ksi0_56 = buffer.data(ksi0 + 56);
    const auto *ksi0_61 = buffer.data(ksi0 + 61);
    const auto *ksi0_65 = buffer.data(ksi0 + 65);
    const auto *ksi0_68 = buffer.data(ksi0 + 68);
    const auto *ksi0_70 = buffer.data(ksi0 + 70);

    const auto *ksh_0 = buffer.data(ksh + 0);
    const auto *ksh_1 = buffer.data(ksh + 1);
    const auto *ksh_2 = buffer.data(ksh + 2);
    const auto *ksh_3 = buffer.data(ksh + 3);
    const auto *ksh_5 = buffer.data(ksh + 5);
    const auto *ksh_6 = buffer.data(ksh + 6);
    const auto *ksh_9 = buffer.data(ksh + 9);
    const auto *ksh_15 = buffer.data(ksh + 15);
    const auto *ksh_17 = buffer.data(ksh + 17);
    const auto *ksh_18 = buffer.data(ksh + 18);
    const auto *ksh_20 = buffer.data(ksh + 20);
    const auto *ksh_21 = buffer.data(ksh + 21);
    const auto *ksh_24 = buffer.data(ksh + 24);
    const auto *ksh_26 = buffer.data(ksh + 26);
    const auto *ksh_27 = buffer.data(ksh + 27);
    const auto *ksh_30 = buffer.data(ksh + 30);
    const auto *ksh_36 = buffer.data(ksh + 36);
    const auto *ksh_38 = buffer.data(ksh + 38);
    const auto *ksh_39 = buffer.data(ksh + 39);
    const auto *ksh_40 = buffer.data(ksh + 40);
    const auto *ksh_41 = buffer.data(ksh + 41);
    const auto *ksh_42 = buffer.data(ksh + 42);
    const auto *ksh_44 = buffer.data(ksh + 44);
    const auto *ksh_47 = buffer.data(ksh + 47);
    const auto *ksh_50 = buffer.data(ksh + 50);
    const auto *ksh_51 = buffer.data(ksh + 51);
    const auto *ksh_57 = buffer.data(ksh + 57);
    const auto *ksh_58 = buffer.data(ksh + 58);
    const auto *ksh_59 = buffer.data(ksh + 59);
    const auto *ksh_60 = buffer.data(ksh + 60);
    const auto *ksh_62 = buffer.data(ksh + 62);
    const auto *ksh_63 = buffer.data(ksh + 63);
    const auto *ksh_66 = buffer.data(ksh + 66);
    const auto *ksh_69 = buffer.data(ksh + 69);
    const auto *ksh_73 = buffer.data(ksh + 73);
    const auto *ksh_78 = buffer.data(ksh + 78);
    const auto *ksh_80 = buffer.data(ksh + 80);
    const auto *ksh_81 = buffer.data(ksh + 81);
    const auto *ksh_82 = buffer.data(ksh + 82);
    const auto *ksh_83 = buffer.data(ksh + 83);
    const auto *ksh_99 = buffer.data(ksh + 99);

    const auto *ksi1_0 = buffer.data(ksi1 + 0);
    const auto *ksi1_3 = buffer.data(ksi1 + 3);
    const auto *ksi1_5 = buffer.data(ksi1 + 5);
    const auto *ksi1_6 = buffer.data(ksi1 + 6);
    const auto *ksi1_9 = buffer.data(ksi1 + 9);
    const auto *ksi1_10 = buffer.data(ksi1 + 10);
    const auto *ksi1_14 = buffer.data(ksi1 + 14);
    const auto *ksi1_21 = buffer.data(ksi1 + 21);
    const auto *ksi1_27 = buffer.data(ksi1 + 27);
    const auto *ksi1_31 = buffer.data(ksi1 + 31);
    const auto *ksi1_34 = buffer.data(ksi1 + 34);
    const auto *ksi1_38 = buffer.data(ksi1 + 38);
    const auto *ksi1_56 = buffer.data(ksi1 + 56);
    const auto *ksi1_61 = buffer.data(ksi1 + 61);
    const auto *ksi1_65 = buffer.data(ksi1 + 65);
    const auto *ksi1_68 = buffer.data(ksi1 + 68);
    const auto *ksi1_70 = buffer.data(ksi1 + 70);

    const auto *lsg0_0 = buffer.data(lsg0 + 0);
    const auto *lsg0_1 = buffer.data(lsg0 + 1);
    const auto *lsg0_2 = buffer.data(lsg0 + 2);
    const auto *lsg0_3 = buffer.data(lsg0 + 3);
    const auto *lsg0_5 = buffer.data(lsg0 + 5);
    const auto *lsg0_10 = buffer.data(lsg0 + 10);
    const auto *lsg0_12 = buffer.data(lsg0 + 12);
    const auto *lsg0_13 = buffer.data(lsg0 + 13);
    const auto *lsg0_14 = buffer.data(lsg0 + 14);
    const auto *lsg0_18 = buffer.data(lsg0 + 18);
    const auto *lsg0_25 = buffer.data(lsg0 + 25);
    const auto *lsg0_26 = buffer.data(lsg0 + 26);
    const auto *lsg0_27 = buffer.data(lsg0 + 27);
    const auto *lsg0_32 = buffer.data(lsg0 + 32);
    const auto *lsg0_34 = buffer.data(lsg0 + 34);
    const auto *lsg0_35 = buffer.data(lsg0 + 35);
    const auto *lsg0_41 = buffer.data(lsg0 + 41);
    const auto *lsg0_42 = buffer.data(lsg0 + 42);
    const auto *lsg0_43 = buffer.data(lsg0 + 43);
    const auto *lsg0_44 = buffer.data(lsg0 + 44);
    const auto *lsg0_45 = buffer.data(lsg0 + 45);
    const auto *lsg0_47 = buffer.data(lsg0 + 47);
    const auto *lsg0_48 = buffer.data(lsg0 + 48);
    const auto *lsg0_50 = buffer.data(lsg0 + 50);
    const auto *lsg0_51 = buffer.data(lsg0 + 51);
    const auto *lsg0_55 = buffer.data(lsg0 + 55);
    const auto *lsg0_56 = buffer.data(lsg0 + 56);
    const auto *lsg0_57 = buffer.data(lsg0 + 57);
    const auto *lsg0_59 = buffer.data(lsg0 + 59);

    const auto *lsg1_0 = buffer.data(lsg1 + 0);
    const auto *lsg1_1 = buffer.data(lsg1 + 1);
    const auto *lsg1_2 = buffer.data(lsg1 + 2);
    const auto *lsg1_3 = buffer.data(lsg1 + 3);
    const auto *lsg1_5 = buffer.data(lsg1 + 5);
    const auto *lsg1_10 = buffer.data(lsg1 + 10);
    const auto *lsg1_12 = buffer.data(lsg1 + 12);
    const auto *lsg1_13 = buffer.data(lsg1 + 13);
    const auto *lsg1_14 = buffer.data(lsg1 + 14);
    const auto *lsg1_18 = buffer.data(lsg1 + 18);
    const auto *lsg1_25 = buffer.data(lsg1 + 25);
    const auto *lsg1_26 = buffer.data(lsg1 + 26);
    const auto *lsg1_27 = buffer.data(lsg1 + 27);
    const auto *lsg1_32 = buffer.data(lsg1 + 32);
    const auto *lsg1_34 = buffer.data(lsg1 + 34);
    const auto *lsg1_35 = buffer.data(lsg1 + 35);
    const auto *lsg1_41 = buffer.data(lsg1 + 41);
    const auto *lsg1_42 = buffer.data(lsg1 + 42);
    const auto *lsg1_43 = buffer.data(lsg1 + 43);
    const auto *lsg1_44 = buffer.data(lsg1 + 44);
    const auto *lsg1_45 = buffer.data(lsg1 + 45);
    const auto *lsg1_47 = buffer.data(lsg1 + 47);
    const auto *lsg1_48 = buffer.data(lsg1 + 48);
    const auto *lsg1_50 = buffer.data(lsg1 + 50);
    const auto *lsg1_51 = buffer.data(lsg1 + 51);
    const auto *lsg1_55 = buffer.data(lsg1 + 55);
    const auto *lsg1_56 = buffer.data(lsg1 + 56);
    const auto *lsg1_57 = buffer.data(lsg1 + 57);
    const auto *lsg1_59 = buffer.data(lsg1 + 59);

    const auto *lsh_0 = buffer.data(lsh + 0);
    const auto *lsh_1 = buffer.data(lsh + 1);
    const auto *lsh_2 = buffer.data(lsh + 2);
    const auto *lsh_3 = buffer.data(lsh + 3);
    const auto *lsh_5 = buffer.data(lsh + 5);
    const auto *lsh_6 = buffer.data(lsh + 6);
    const auto *lsh_8 = buffer.data(lsh + 8);
    const auto *lsh_9 = buffer.data(lsh + 9);
    const auto *lsh_10 = buffer.data(lsh + 10);
    const auto *lsh_14 = buffer.data(lsh + 14);
    const auto *lsh_15 = buffer.data(lsh + 15);
    const auto *lsh_17 = buffer.data(lsh + 17);
    const auto *lsh_18 = buffer.data(lsh + 18);
    const auto *lsh_19 = buffer.data(lsh + 19);
    const auto *lsh_20 = buffer.data(lsh + 20);
    const auto *lsh_21 = buffer.data(lsh + 21);
    const auto *lsh_22 = buffer.data(lsh + 22);
    const auto *lsh_24 = buffer.data(lsh + 24);
    const auto *lsh_26 = buffer.data(lsh + 26);
    const auto *lsh_27 = buffer.data(lsh + 27);
    const auto *lsh_28 = buffer.data(lsh + 28);
    const auto *lsh_30 = buffer.data(lsh + 30);
    const auto *lsh_31 = buffer.data(lsh + 31);
    const auto *lsh_36 = buffer.data(lsh + 36);
    const auto *lsh_37 = buffer.data(lsh + 37);
    const auto *lsh_38 = buffer.data(lsh + 38);
    const auto *lsh_39 = buffer.data(lsh + 39);
    const auto *lsh_40 = buffer.data(lsh + 40);
    const auto *lsh_41 = buffer.data(lsh + 41);
    const auto *lsh_42 = buffer.data(lsh + 42);
    const auto *lsh_44 = buffer.data(lsh + 44);
    const auto *lsh_46 = buffer.data(lsh + 46);
    const auto *lsh_47 = buffer.data(lsh + 47);
    const auto *lsh_49 = buffer.data(lsh + 49);
    const auto *lsh_50 = buffer.data(lsh + 50);
    const auto *lsh_51 = buffer.data(lsh + 51);
    const auto *lsh_56 = buffer.data(lsh + 56);
    const auto *lsh_57 = buffer.data(lsh + 57);
    const auto *lsh_58 = buffer.data(lsh + 58);
    const auto *lsh_59 = buffer.data(lsh + 59);
    const auto *lsh_60 = buffer.data(lsh + 60);
    const auto *lsh_61 = buffer.data(lsh + 61);
    const auto *lsh_62 = buffer.data(lsh + 62);
    const auto *lsh_63 = buffer.data(lsh + 63);
    const auto *lsh_64 = buffer.data(lsh + 64);
    const auto *lsh_65 = buffer.data(lsh + 65);
    const auto *lsh_66 = buffer.data(lsh + 66);
    const auto *lsh_68 = buffer.data(lsh + 68);
    const auto *lsh_69 = buffer.data(lsh + 69);
    const auto *lsh_70 = buffer.data(lsh + 70);
    const auto *lsh_72 = buffer.data(lsh + 72);
    const auto *lsh_73 = buffer.data(lsh + 73);
    const auto *lsh_78 = buffer.data(lsh + 78);
    const auto *lsh_79 = buffer.data(lsh + 79);
    const auto *lsh_80 = buffer.data(lsh + 80);
    const auto *lsh_81 = buffer.data(lsh + 81);
    const auto *lsh_82 = buffer.data(lsh + 82);
    const auto *lsh_83 = buffer.data(lsh + 83);
    const auto *lsh_84 = buffer.data(lsh + 84);
    const auto *lsh_86 = buffer.data(lsh + 86);
    const auto *lsh_87 = buffer.data(lsh + 87);
    const auto *lsh_89 = buffer.data(lsh + 89);
    const auto *lsh_90 = buffer.data(lsh + 90);
    const auto *lsh_93 = buffer.data(lsh + 93);
    const auto *lsh_99 = buffer.data(lsh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ksh_0, lsg0_0, \
                         lsg1_0, lsh_0, lsh_1, lsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ksh_0[k]
                 + f_1 * lsg0_0[k]
                 - f_2 * lsg1_0[k]
                 + f_3 * pc_x[k] * lsh_0[k];

        t_1[k] = f_3 * pc_y[k] * lsh_0[k];

        t_2[k] = f_3 * pc_z[k] * lsh_0[k];

        t_3[k] = f_4 * lsg0_0[k]
                 - f_5 * lsg1_0[k]
                 + f_3 * pc_y[k] * lsh_1[k];

        t_4[k] = f_3 * pc_y[k] * lsh_2[k];

        t_5[k] = f_4 * lsg0_0[k]
                 - f_5 * lsg1_0[k]
                 + f_3 * pc_z[k] * lsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, lsg0_1, lsg0_2, lsg0_3, lsg1_1, \
                         lsg1_2, lsg1_3, lsh_3, lsh_5, lsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * lsg0_1[k]
                 - f_7 * lsg1_1[k]
                 + f_3 * pc_y[k] * lsh_3[k];

        t_7[k] = f_3 * pc_z[k] * lsh_3[k];

        t_8[k] = f_3 * pc_y[k] * lsh_5[k];

        t_9[k] = f_6 * lsg0_2[k]
                 - f_7 * lsg1_2[k]
                 + f_3 * pc_z[k] * lsh_5[k];

        t_10[k] = f_8 * lsg0_3[k]
                  - f_9 * lsg1_3[k]
                  + f_3 * pc_y[k] * lsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, ksh_15, lsg0_5, \
                         lsg1_5, lsh_6, lsh_8, lsh_9, lsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * lsh_6[k];

        t_12[k] = f_4 * lsg0_5[k]
                  - f_5 * lsg1_5[k]
                  + f_3 * pc_y[k] * lsh_8[k];

        t_13[k] = f_3 * pc_y[k] * lsh_9[k];

        t_14[k] = f_8 * lsg0_5[k]
                  - f_9 * lsg1_5[k]
                  + f_3 * pc_z[k] * lsh_9[k];

        t_15[k] = f_0 * ksh_15[k]
                  + f_3 * pc_x[k] * lsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, ksh_17, ksh_18, \
                         ksh_20, lsh_10, lsh_14, lsh_17, lsh_18, \
                         lsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * lsh_10[k];

        t_17[k] = f_0 * ksh_17[k]
                  + f_3 * pc_x[k] * lsh_17[k];

        t_18[k] = f_0 * ksh_18[k]
                  + f_3 * pc_x[k] * lsh_18[k];

        t_19[k] = f_3 * pc_y[k] * lsh_14[k];

        t_20[k] = f_0 * ksh_20[k]
                  + f_3 * pc_x[k] * lsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, lsg0_10, lsg0_12, lsg0_13, \
                         lsg1_10, lsg1_12, lsg1_13, lsh_15, lsh_17, \
                         lsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * lsg0_10[k]
                  - f_2 * lsg1_10[k]
                  + f_3 * pc_y[k] * lsh_15[k];

        t_22[k] = f_3 * pc_z[k] * lsh_15[k];

        t_23[k] = f_8 * lsg0_12[k]
                  - f_9 * lsg1_12[k]
                  + f_3 * pc_y[k] * lsh_17[k];

        t_24[k] = f_6 * lsg0_13[k]
                  - f_7 * lsg1_13[k]
                  + f_3 * pc_y[k] * lsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, ksi0_0, ksh_0, \
                         ksi1_0, lsg0_14, lsg1_14, lsh_19, lsh_20, \
                         lsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * lsg0_14[k]
                  - f_5 * lsg1_14[k]
                  + f_3 * pc_y[k] * lsh_19[k];

        t_26[k] = f_3 * pc_y[k] * lsh_20[k];

        t_27[k] = f_1 * lsg0_14[k]
                  - f_2 * lsg1_14[k]
                  + f_3 * pc_z[k] * lsh_20[k];

        t_28[k] = pa_y[k] * ksi0_0[k]
                  - f_10 * pc_y[k] * ksi1_0[k];

        t_29[k] = f_11 * ksh_0[k]
                  + f_3 * pc_y[k] * lsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, ksi0_3, ksi0_5, ksh_1, \
                         ksi1_3, ksi1_5, lsh_21, lsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * lsh_21[k];

        t_31[k] = pa_y[k] * ksi0_3[k]
                  + f_12 * ksh_1[k]
                  - f_10 * pc_y[k] * ksi1_3[k];

        t_32[k] = f_3 * pc_z[k] * lsh_22[k];

        t_33[k] = pa_y[k] * ksi0_5[k]
                  - f_10 * pc_y[k] * ksi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, ksi0_6, ksi0_9, ksh_3, \
                         ksh_5, ksi1_6, ksi1_9, lsh_24, lsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * ksi0_6[k]
                  + f_13 * ksh_3[k]
                  - f_10 * pc_y[k] * ksi1_6[k];

        t_35[k] = f_3 * pc_z[k] * lsh_24[k];

        t_36[k] = f_11 * ksh_5[k]
                  + f_3 * pc_y[k] * lsh_26[k];

        t_37[k] = pa_y[k] * ksi0_9[k]
                  - f_10 * pc_y[k] * ksi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, ksi0_10, ksh_6, ksh_9, \
                         ksi1_10, lsg0_18, lsg1_18, lsh_27, lsh_28, \
                         lsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * ksi0_10[k]
                  + f_14 * ksh_6[k]
                  - f_10 * pc_y[k] * ksi1_10[k];

        t_39[k] = f_3 * pc_z[k] * lsh_27[k];

        t_40[k] = f_4 * lsg0_18[k]
                  - f_5 * lsg1_18[k]
                  + f_3 * pc_z[k] * lsh_28[k];

        t_41[k] = f_11 * ksh_9[k]
                  + f_3 * pc_y[k] * lsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, ksi0_14, ksh_36, \
                         ksh_38, ksi1_14, lsh_31, lsh_36, lsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * ksi0_14[k]
                  - f_10 * pc_y[k] * ksi1_14[k];

        t_43[k] = f_15 * ksh_36[k]
                  + f_3 * pc_x[k] * lsh_36[k];

        t_44[k] = f_3 * pc_z[k] * lsh_31[k];

        t_45[k] = f_15 * ksh_38[k]
                  + f_3 * pc_x[k] * lsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, ksh_15, ksh_39, ksh_40, ksh_41, \
                         lsg0_25, lsg1_25, lsh_36, lsh_39, lsh_40, \
                         lsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * ksh_39[k]
                  + f_3 * pc_x[k] * lsh_39[k];

        t_47[k] = f_15 * ksh_40[k]
                  + f_3 * pc_x[k] * lsh_40[k];

        t_48[k] = f_15 * ksh_41[k]
                  + f_3 * pc_x[k] * lsh_41[k];

        t_49[k] = f_11 * ksh_15[k]
                  + f_1 * lsg0_25[k]
                  - f_2 * lsg1_25[k]
                  + f_3 * pc_y[k] * lsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, lsg0_25, lsg0_26, lsg0_27, lsg1_25, \
                         lsg1_26, lsg1_27, lsh_36, lsh_37, lsh_38, \
                         lsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * lsh_36[k];

        t_51[k] = f_4 * lsg0_25[k]
                  - f_5 * lsg1_25[k]
                  + f_3 * pc_z[k] * lsh_37[k];

        t_52[k] = f_6 * lsg0_26[k]
                  - f_7 * lsg1_26[k]
                  + f_3 * pc_z[k] * lsh_38[k];

        t_53[k] = f_8 * lsg0_27[k]
                  - f_9 * lsg1_27[k]
                  + f_3 * pc_z[k] * lsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, ksi0_0, ksi0_27, \
                         ksh_20, ksi1_0, ksi1_27, lsh_41, lsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * ksh_20[k]
                  + f_3 * pc_y[k] * lsh_41[k];

        t_55[k] = pa_y[k] * ksi0_27[k]
                  - f_10 * pc_y[k] * ksi1_27[k];

        t_56[k] = pa_z[k] * ksi0_0[k]
                  - f_10 * pc_z[k] * ksi1_0[k];

        t_57[k] = f_3 * pc_y[k] * lsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, ksi0_3, ksi0_5, ksh_0, \
                         ksh_2, ksi1_3, ksi1_5, lsh_42, lsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * ksh_0[k]
                  + f_3 * pc_z[k] * lsh_42[k];

        t_59[k] = pa_z[k] * ksi0_3[k]
                  - f_10 * pc_z[k] * ksi1_3[k];

        t_60[k] = f_3 * pc_y[k] * lsh_44[k];

        t_61[k] = pa_z[k] * ksi0_5[k]
                  + f_12 * ksh_2[k]
                  - f_10 * pc_z[k] * ksi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, ksi0_6, ksi0_9, ksh_5, \
                         ksi1_6, ksi1_9, lsg0_32, lsg1_32, lsh_46, \
                         lsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * ksi0_6[k]
                  - f_10 * pc_z[k] * ksi1_6[k];

        t_63[k] = f_4 * lsg0_32[k]
                  - f_5 * lsg1_32[k]
                  + f_3 * pc_y[k] * lsh_46[k];

        t_64[k] = f_3 * pc_y[k] * lsh_47[k];

        t_65[k] = pa_z[k] * ksi0_9[k]
                  + f_13 * ksh_5[k]
                  - f_10 * pc_z[k] * ksi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, ksi0_10, ksi1_10, lsg0_34, \
                         lsg0_35, lsg1_34, lsg1_35, lsh_49, lsh_50, \
                         lsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * ksi0_10[k]
                  - f_10 * pc_z[k] * ksi1_10[k];

        t_67[k] = f_6 * lsg0_34[k]
                  - f_7 * lsg1_34[k]
                  + f_3 * pc_y[k] * lsh_49[k];

        t_68[k] = f_4 * lsg0_35[k]
                  - f_5 * lsg1_35[k]
                  + f_3 * pc_y[k] * lsh_50[k];

        t_69[k] = f_3 * pc_y[k] * lsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, ksi0_14, ksh_9, ksh_57, \
                         ksh_58, ksh_59, ksi1_14, lsh_57, lsh_58, \
                         lsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * ksi0_14[k]
                  + f_14 * ksh_9[k]
                  - f_10 * pc_z[k] * ksi1_14[k];

        t_71[k] = f_15 * ksh_57[k]
                  + f_3 * pc_x[k] * lsh_57[k];

        t_72[k] = f_15 * ksh_58[k]
                  + f_3 * pc_x[k] * lsh_58[k];

        t_73[k] = f_15 * ksh_59[k]
                  + f_3 * pc_x[k] * lsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, ksi0_21, ksh_60, \
                         ksh_62, ksi1_21, lsh_56, lsh_60, lsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * ksh_60[k]
                  + f_3 * pc_x[k] * lsh_60[k];

        t_75[k] = f_3 * pc_y[k] * lsh_56[k];

        t_76[k] = f_15 * ksh_62[k]
                  + f_3 * pc_x[k] * lsh_62[k];

        t_77[k] = pa_z[k] * ksi0_21[k]
                  - f_10 * pc_z[k] * ksi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, lsg0_41, lsg0_42, lsg0_43, lsg1_41, lsg1_42, \
                         lsg1_43, lsh_58, lsh_59, lsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * lsg0_41[k]
                  - f_17 * lsg1_41[k]
                  + f_3 * pc_y[k] * lsh_58[k];

        t_79[k] = f_8 * lsg0_42[k]
                  - f_9 * lsg1_42[k]
                  + f_3 * pc_y[k] * lsh_59[k];

        t_80[k] = f_6 * lsg0_43[k]
                  - f_7 * lsg1_43[k]
                  + f_3 * pc_y[k] * lsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, ksh_20, ksh_63, lsg0_44, \
                         lsg0_45, lsg1_44, lsg1_45, lsh_61, lsh_62, \
                         lsh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * lsg0_44[k]
                  - f_5 * lsg1_44[k]
                  + f_3 * pc_y[k] * lsh_61[k];

        t_82[k] = f_3 * pc_y[k] * lsh_62[k];

        t_83[k] = f_11 * ksh_20[k]
                  + f_1 * lsg0_44[k]
                  - f_2 * lsg1_44[k]
                  + f_3 * pc_z[k] * lsh_62[k];

        t_84[k] = f_18 * ksh_63[k]
                  + f_1 * lsg0_45[k]
                  - f_2 * lsg1_45[k]
                  + f_3 * pc_x[k] * lsh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, ksh_21, ksh_66, lsg0_48, \
                         lsg1_48, lsh_63, lsh_64, lsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * ksh_21[k]
                  + f_3 * pc_y[k] * lsh_63[k];

        t_86[k] = f_3 * pc_z[k] * lsh_63[k];

        t_87[k] = f_18 * ksh_66[k]
                  + f_8 * lsg0_48[k]
                  - f_9 * lsg1_48[k]
                  + f_3 * pc_x[k] * lsh_66[k];

        t_88[k] = f_3 * pc_z[k] * lsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, ksh_69, lsg0_45, lsg0_51, lsg1_45, \
                         lsg1_51, lsh_65, lsh_66, lsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * lsg0_45[k]
                  - f_5 * lsg1_45[k]
                  + f_3 * pc_z[k] * lsh_65[k];

        t_90[k] = f_18 * ksh_69[k]
                  + f_6 * lsg0_51[k]
                  - f_7 * lsg1_51[k]
                  + f_3 * pc_x[k] * lsh_69[k];

        t_91[k] = f_3 * pc_z[k] * lsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, ksh_26, ksh_73, lsg0_47, \
                         lsg0_55, lsg1_47, lsg1_55, lsh_68, lsh_69, \
                         lsh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * ksh_26[k]
                  + f_3 * pc_y[k] * lsh_68[k];

        t_93[k] = f_6 * lsg0_47[k]
                  - f_7 * lsg1_47[k]
                  + f_3 * pc_z[k] * lsh_68[k];

        t_94[k] = f_18 * ksh_73[k]
                  + f_4 * lsg0_55[k]
                  - f_5 * lsg1_55[k]
                  + f_3 * pc_x[k] * lsh_73[k];

        t_95[k] = f_3 * pc_z[k] * lsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, ksh_30, ksh_78, lsg0_48, \
                         lsg0_50, lsg1_48, lsg1_50, lsh_70, lsh_72, \
                         lsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * lsg0_48[k]
                  - f_5 * lsg1_48[k]
                  + f_3 * pc_z[k] * lsh_70[k];

        t_97[k] = f_12 * ksh_30[k]
                  + f_3 * pc_y[k] * lsh_72[k];

        t_98[k] = f_8 * lsg0_50[k]
                  - f_9 * lsg1_50[k]
                  + f_3 * pc_z[k] * lsh_72[k];

        t_99[k] = f_18 * ksh_78[k]
                  + f_3 * pc_x[k] * lsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, ksh_80, ksh_81, \
                         ksh_82, ksh_83, lsh_73, lsh_80, lsh_81, lsh_82, \
                         lsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * lsh_73[k];

        t_101[k] = f_18 * ksh_80[k]
                   + f_3 * pc_x[k] * lsh_80[k];

        t_102[k] = f_18 * ksh_81[k]
                   + f_3 * pc_x[k] * lsh_81[k];

        t_103[k] = f_18 * ksh_82[k]
                   + f_3 * pc_x[k] * lsh_82[k];

        t_104[k] = f_18 * ksh_83[k]
                   + f_3 * pc_x[k] * lsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, ksh_36, lsg0_55, lsg0_56, \
                         lsg1_55, lsg1_56, lsh_78, lsh_79, lsh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * ksh_36[k]
                   + f_1 * lsg0_55[k]
                   - f_2 * lsg1_55[k]
                   + f_3 * pc_y[k] * lsh_78[k];

        t_106[k] = f_3 * pc_z[k] * lsh_78[k];

        t_107[k] = f_4 * lsg0_55[k]
                   - f_5 * lsg1_55[k]
                   + f_3 * pc_z[k] * lsh_79[k];

        t_108[k] = f_6 * lsg0_56[k]
                   - f_7 * lsg1_56[k]
                   + f_3 * pc_z[k] * lsh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, ksi0_56, ksh_41, \
                         ksi1_56, lsg0_57, lsg0_59, lsg1_57, lsg1_59, lsh_81, \
                         lsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * lsg0_57[k]
                   - f_9 * lsg1_57[k]
                   + f_3 * pc_z[k] * lsh_81[k];

        t_110[k] = f_12 * ksh_41[k]
                   + f_3 * pc_y[k] * lsh_83[k];

        t_111[k] = f_1 * lsg0_59[k]
                   - f_2 * lsg1_59[k]
                   + f_3 * pc_z[k] * lsh_83[k];

        t_112[k] = pa_y[k] * ksi0_56[k]
                   - f_10 * pc_y[k] * ksi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, ksi0_31, ksh_21, \
                         ksh_42, ksh_44, ksi1_31, lsh_84, lsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * ksh_42[k]
                   + f_3 * pc_y[k] * lsh_84[k];

        t_114[k] = f_11 * ksh_21[k]
                   + f_3 * pc_z[k] * lsh_84[k];

        t_115[k] = pa_z[k] * ksi0_31[k]
                   - f_10 * pc_z[k] * ksi1_31[k];

        t_116[k] = f_11 * ksh_44[k]
                   + f_3 * pc_y[k] * lsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, ksi0_34, ksi0_61, \
                         ksh_24, ksh_47, ksi1_34, ksi1_61, lsh_87, \
                         lsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * ksi0_61[k]
                   - f_10 * pc_y[k] * ksi1_61[k];

        t_118[k] = pa_z[k] * ksi0_34[k]
                   - f_10 * pc_z[k] * ksi1_34[k];

        t_119[k] = f_11 * ksh_24[k]
                   + f_3 * pc_z[k] * lsh_87[k];

        t_120[k] = f_11 * ksh_47[k]
                   + f_3 * pc_y[k] * lsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, ksi0_38, ksi0_65, \
                         ksh_27, ksi1_38, ksi1_65, lsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * ksi0_65[k]
                   - f_10 * pc_y[k] * ksi1_65[k];

        t_122[k] = pa_z[k] * ksi0_38[k]
                   - f_10 * pc_z[k] * ksi1_38[k];

        t_123[k] = f_11 * ksh_27[k]
                   + f_3 * pc_z[k] * lsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, ksi0_68, ksi0_70, \
                         ksh_50, ksh_51, ksh_99, ksi1_68, ksi1_70, lsh_93, \
                         lsh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * ksi0_68[k]
                   + f_12 * ksh_50[k]
                   - f_10 * pc_y[k] * ksi1_68[k];

        t_125[k] = f_11 * ksh_51[k]
                   + f_3 * pc_y[k] * lsh_93[k];

        t_126[k] = pa_y[k] * ksi0_70[k]
                   - f_10 * pc_y[k] * ksi1_70[k];

        t_127[k] = f_18 * ksh_99[k]
                   + f_3 * pc_x[k] * lsh_99[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_49 = buffer.data(ksi0 + 49);
    const auto *ksi0_83 = buffer.data(ksi0 + 83);
    const auto *ksi0_84 = buffer.data(ksi0 + 84);
    const auto *ksi0_87 = buffer.data(ksi0 + 87);
    const auto *ksi0_90 = buffer.data(ksi0 + 90);
    const auto *ksi0_94 = buffer.data(ksi0 + 94);
    const auto *ksi0_96 = buffer.data(ksi0 + 96);
    const auto *ksi0_105 = buffer.data(ksi0 + 105);
    const auto *ksi0_140 = buffer.data(ksi0 + 140);
    const auto *ksi0_143 = buffer.data(ksi0 + 143);
    const auto *ksi0_145 = buffer.data(ksi0 + 145);
    const auto *ksi0_146 = buffer.data(ksi0 + 146);
    const auto *ksi0_149 = buffer.data(ksi0 + 149);
    const auto *ksi0_150 = buffer.data(ksi0 + 150);
    const auto *ksi0_152 = buffer.data(ksi0 + 152);
    const auto *ksi0_154 = buffer.data(ksi0 + 154);

    const auto *ksh_36 = buffer.data(ksh + 36);
    const auto *ksh_42 = buffer.data(ksh + 42);
    const auto *ksh_59 = buffer.data(ksh + 59);
    const auto *ksh_60 = buffer.data(ksh + 60);
    const auto *ksh_61 = buffer.data(ksh + 61);
    const auto *ksh_62 = buffer.data(ksh + 62);
    const auto *ksh_63 = buffer.data(ksh + 63);
    const auto *ksh_66 = buffer.data(ksh + 66);
    const auto *ksh_68 = buffer.data(ksh + 68);
    const auto *ksh_69 = buffer.data(ksh + 69);
    const auto *ksh_70 = buffer.data(ksh + 70);
    const auto *ksh_72 = buffer.data(ksh + 72);
    const auto *ksh_78 = buffer.data(ksh + 78);
    const auto *ksh_83 = buffer.data(ksh + 83);
    const auto *ksh_84 = buffer.data(ksh + 84);
    const auto *ksh_86 = buffer.data(ksh + 86);
    const auto *ksh_87 = buffer.data(ksh + 87);
    const auto *ksh_89 = buffer.data(ksh + 89);
    const auto *ksh_90 = buffer.data(ksh + 90);
    const auto *ksh_93 = buffer.data(ksh + 93);
    const auto *ksh_99 = buffer.data(ksh + 99);
    const auto *ksh_100 = buffer.data(ksh + 100);
    const auto *ksh_101 = buffer.data(ksh + 101);
    const auto *ksh_102 = buffer.data(ksh + 102);
    const auto *ksh_103 = buffer.data(ksh + 103);
    const auto *ksh_104 = buffer.data(ksh + 104);
    const auto *ksh_105 = buffer.data(ksh + 105);
    const auto *ksh_106 = buffer.data(ksh + 106);
    const auto *ksh_107 = buffer.data(ksh + 107);
    const auto *ksh_108 = buffer.data(ksh + 108);
    const auto *ksh_110 = buffer.data(ksh + 110);
    const auto *ksh_111 = buffer.data(ksh + 111);
    const auto *ksh_113 = buffer.data(ksh + 113);
    const auto *ksh_114 = buffer.data(ksh + 114);
    const auto *ksh_119 = buffer.data(ksh + 119);
    const auto *ksh_120 = buffer.data(ksh + 120);
    const auto *ksh_121 = buffer.data(ksh + 121);
    const auto *ksh_122 = buffer.data(ksh + 122);
    const auto *ksh_123 = buffer.data(ksh + 123);
    const auto *ksh_125 = buffer.data(ksh + 125);
    const auto *ksh_126 = buffer.data(ksh + 126);
    const auto *ksh_129 = buffer.data(ksh + 129);
    const auto *ksh_132 = buffer.data(ksh + 132);
    const auto *ksh_136 = buffer.data(ksh + 136);
    const auto *ksh_141 = buffer.data(ksh + 141);
    const auto *ksh_143 = buffer.data(ksh + 143);
    const auto *ksh_144 = buffer.data(ksh + 144);
    const auto *ksh_145 = buffer.data(ksh + 145);
    const auto *ksh_146 = buffer.data(ksh + 146);
    const auto *ksh_152 = buffer.data(ksh + 152);
    const auto *ksh_156 = buffer.data(ksh + 156);
    const auto *ksh_161 = buffer.data(ksh + 161);
    const auto *ksh_162 = buffer.data(ksh + 162);
    const auto *ksh_163 = buffer.data(ksh + 163);
    const auto *ksh_164 = buffer.data(ksh + 164);
    const auto *ksh_165 = buffer.data(ksh + 165);
    const auto *ksh_166 = buffer.data(ksh + 166);
    const auto *ksh_167 = buffer.data(ksh + 167);
    const auto *ksh_183 = buffer.data(ksh + 183);
    const auto *ksh_184 = buffer.data(ksh + 184);
    const auto *ksh_185 = buffer.data(ksh + 185);
    const auto *ksh_186 = buffer.data(ksh + 186);
    const auto *ksh_187 = buffer.data(ksh + 187);
    const auto *ksh_188 = buffer.data(ksh + 188);

    const auto *ksi1_49 = buffer.data(ksi1 + 49);
    const auto *ksi1_83 = buffer.data(ksi1 + 83);
    const auto *ksi1_84 = buffer.data(ksi1 + 84);
    const auto *ksi1_87 = buffer.data(ksi1 + 87);
    const auto *ksi1_90 = buffer.data(ksi1 + 90);
    const auto *ksi1_94 = buffer.data(ksi1 + 94);
    const auto *ksi1_96 = buffer.data(ksi1 + 96);
    const auto *ksi1_105 = buffer.data(ksi1 + 105);
    const auto *ksi1_140 = buffer.data(ksi1 + 140);
    const auto *ksi1_143 = buffer.data(ksi1 + 143);
    const auto *ksi1_145 = buffer.data(ksi1 + 145);
    const auto *ksi1_146 = buffer.data(ksi1 + 146);
    const auto *ksi1_149 = buffer.data(ksi1 + 149);
    const auto *ksi1_150 = buffer.data(ksi1 + 150);
    const auto *ksi1_152 = buffer.data(ksi1 + 152);
    const auto *ksi1_154 = buffer.data(ksi1 + 154);

    const auto *lsg0_72 = buffer.data(lsg0 + 72);
    const auto *lsg0_73 = buffer.data(lsg0 + 73);
    const auto *lsg0_74 = buffer.data(lsg0 + 74);
    const auto *lsg0_75 = buffer.data(lsg0 + 75);
    const auto *lsg0_76 = buffer.data(lsg0 + 76);
    const auto *lsg0_77 = buffer.data(lsg0 + 77);
    const auto *lsg0_78 = buffer.data(lsg0 + 78);
    const auto *lsg0_79 = buffer.data(lsg0 + 79);
    const auto *lsg0_80 = buffer.data(lsg0 + 80);
    const auto *lsg0_84 = buffer.data(lsg0 + 84);
    const auto *lsg0_85 = buffer.data(lsg0 + 85);
    const auto *lsg0_86 = buffer.data(lsg0 + 86);
    const auto *lsg0_87 = buffer.data(lsg0 + 87);
    const auto *lsg0_88 = buffer.data(lsg0 + 88);
    const auto *lsg0_89 = buffer.data(lsg0 + 89);
    const auto *lsg0_90 = buffer.data(lsg0 + 90);
    const auto *lsg0_92 = buffer.data(lsg0 + 92);
    const auto *lsg0_93 = buffer.data(lsg0 + 93);
    const auto *lsg0_95 = buffer.data(lsg0 + 95);
    const auto *lsg0_96 = buffer.data(lsg0 + 96);
    const auto *lsg0_100 = buffer.data(lsg0 + 100);
    const auto *lsg0_101 = buffer.data(lsg0 + 101);
    const auto *lsg0_102 = buffer.data(lsg0 + 102);
    const auto *lsg0_104 = buffer.data(lsg0 + 104);
    const auto *lsg0_110 = buffer.data(lsg0 + 110);
    const auto *lsg0_114 = buffer.data(lsg0 + 114);
    const auto *lsg0_117 = buffer.data(lsg0 + 117);
    const auto *lsg0_118 = buffer.data(lsg0 + 118);
    const auto *lsg0_119 = buffer.data(lsg0 + 119);
    const auto *lsg0_130 = buffer.data(lsg0 + 130);
    const auto *lsg0_132 = buffer.data(lsg0 + 132);

    const auto *lsg1_72 = buffer.data(lsg1 + 72);
    const auto *lsg1_73 = buffer.data(lsg1 + 73);
    const auto *lsg1_74 = buffer.data(lsg1 + 74);
    const auto *lsg1_75 = buffer.data(lsg1 + 75);
    const auto *lsg1_76 = buffer.data(lsg1 + 76);
    const auto *lsg1_77 = buffer.data(lsg1 + 77);
    const auto *lsg1_78 = buffer.data(lsg1 + 78);
    const auto *lsg1_79 = buffer.data(lsg1 + 79);
    const auto *lsg1_80 = buffer.data(lsg1 + 80);
    const auto *lsg1_84 = buffer.data(lsg1 + 84);
    const auto *lsg1_85 = buffer.data(lsg1 + 85);
    const auto *lsg1_86 = buffer.data(lsg1 + 86);
    const auto *lsg1_87 = buffer.data(lsg1 + 87);
    const auto *lsg1_88 = buffer.data(lsg1 + 88);
    const auto *lsg1_89 = buffer.data(lsg1 + 89);
    const auto *lsg1_90 = buffer.data(lsg1 + 90);
    const auto *lsg1_92 = buffer.data(lsg1 + 92);
    const auto *lsg1_93 = buffer.data(lsg1 + 93);
    const auto *lsg1_95 = buffer.data(lsg1 + 95);
    const auto *lsg1_96 = buffer.data(lsg1 + 96);
    const auto *lsg1_100 = buffer.data(lsg1 + 100);
    const auto *lsg1_101 = buffer.data(lsg1 + 101);
    const auto *lsg1_102 = buffer.data(lsg1 + 102);
    const auto *lsg1_104 = buffer.data(lsg1 + 104);
    const auto *lsg1_110 = buffer.data(lsg1 + 110);
    const auto *lsg1_114 = buffer.data(lsg1 + 114);
    const auto *lsg1_117 = buffer.data(lsg1 + 117);
    const auto *lsg1_118 = buffer.data(lsg1 + 118);
    const auto *lsg1_119 = buffer.data(lsg1 + 119);
    const auto *lsg1_130 = buffer.data(lsg1 + 130);
    const auto *lsg1_132 = buffer.data(lsg1 + 132);

    const auto *lsh_99 = buffer.data(lsh + 99);
    const auto *lsh_100 = buffer.data(lsh + 100);
    const auto *lsh_101 = buffer.data(lsh + 101);
    const auto *lsh_102 = buffer.data(lsh + 102);
    const auto *lsh_103 = buffer.data(lsh + 103);
    const auto *lsh_104 = buffer.data(lsh + 104);
    const auto *lsh_105 = buffer.data(lsh + 105);
    const auto *lsh_106 = buffer.data(lsh + 106);
    const auto *lsh_107 = buffer.data(lsh + 107);
    const auto *lsh_108 = buffer.data(lsh + 108);
    const auto *lsh_109 = buffer.data(lsh + 109);
    const auto *lsh_110 = buffer.data(lsh + 110);
    const auto *lsh_111 = buffer.data(lsh + 111);
    const auto *lsh_112 = buffer.data(lsh + 112);
    const auto *lsh_113 = buffer.data(lsh + 113);
    const auto *lsh_114 = buffer.data(lsh + 114);
    const auto *lsh_119 = buffer.data(lsh + 119);
    const auto *lsh_120 = buffer.data(lsh + 120);
    const auto *lsh_121 = buffer.data(lsh + 121);
    const auto *lsh_122 = buffer.data(lsh + 122);
    const auto *lsh_123 = buffer.data(lsh + 123);
    const auto *lsh_124 = buffer.data(lsh + 124);
    const auto *lsh_125 = buffer.data(lsh + 125);
    const auto *lsh_126 = buffer.data(lsh + 126);
    const auto *lsh_127 = buffer.data(lsh + 127);
    const auto *lsh_128 = buffer.data(lsh + 128);
    const auto *lsh_129 = buffer.data(lsh + 129);
    const auto *lsh_131 = buffer.data(lsh + 131);
    const auto *lsh_132 = buffer.data(lsh + 132);
    const auto *lsh_133 = buffer.data(lsh + 133);
    const auto *lsh_135 = buffer.data(lsh + 135);
    const auto *lsh_136 = buffer.data(lsh + 136);
    const auto *lsh_141 = buffer.data(lsh + 141);
    const auto *lsh_142 = buffer.data(lsh + 142);
    const auto *lsh_143 = buffer.data(lsh + 143);
    const auto *lsh_144 = buffer.data(lsh + 144);
    const auto *lsh_145 = buffer.data(lsh + 145);
    const auto *lsh_146 = buffer.data(lsh + 146);
    const auto *lsh_147 = buffer.data(lsh + 147);
    const auto *lsh_149 = buffer.data(lsh + 149);
    const auto *lsh_150 = buffer.data(lsh + 150);
    const auto *lsh_152 = buffer.data(lsh + 152);
    const auto *lsh_153 = buffer.data(lsh + 153);
    const auto *lsh_156 = buffer.data(lsh + 156);
    const auto *lsh_161 = buffer.data(lsh + 161);
    const auto *lsh_162 = buffer.data(lsh + 162);
    const auto *lsh_163 = buffer.data(lsh + 163);
    const auto *lsh_164 = buffer.data(lsh + 164);
    const auto *lsh_165 = buffer.data(lsh + 165);
    const auto *lsh_166 = buffer.data(lsh + 166);
    const auto *lsh_167 = buffer.data(lsh + 167);
    const auto *lsh_168 = buffer.data(lsh + 168);
    const auto *lsh_170 = buffer.data(lsh + 170);
    const auto *lsh_171 = buffer.data(lsh + 171);
    const auto *lsh_173 = buffer.data(lsh + 173);
    const auto *lsh_174 = buffer.data(lsh + 174);
    const auto *lsh_177 = buffer.data(lsh + 177);
    const auto *lsh_183 = buffer.data(lsh + 183);
    const auto *lsh_184 = buffer.data(lsh + 184);
    const auto *lsh_185 = buffer.data(lsh + 185);
    const auto *lsh_186 = buffer.data(lsh + 186);
    const auto *lsh_187 = buffer.data(lsh + 187);
    const auto *lsh_188 = buffer.data(lsh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, ksh_100, ksh_101, ksh_102, \
                         ksh_103, ksh_104, lsh_100, lsh_101, lsh_102, lsh_103, \
                         lsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * ksh_100[k]
                   + f_3 * pc_x[k] * lsh_100[k];

        t_129[k] = f_18 * ksh_101[k]
                   + f_3 * pc_x[k] * lsh_101[k];

        t_130[k] = f_18 * ksh_102[k]
                   + f_3 * pc_x[k] * lsh_102[k];

        t_131[k] = f_18 * ksh_103[k]
                   + f_3 * pc_x[k] * lsh_103[k];

        t_132[k] = f_18 * ksh_104[k]
                   + f_3 * pc_x[k] * lsh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, ksi0_49, ksh_36, ksh_59, \
                         ksi1_49, lsg0_72, lsg1_72, lsh_99, lsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * ksi0_49[k]
                   - f_10 * pc_z[k] * ksi1_49[k];

        t_134[k] = f_11 * ksh_36[k]
                   + f_3 * pc_z[k] * lsh_99[k];

        t_135[k] = f_11 * ksh_59[k]
                   + f_8 * lsg0_72[k]
                   - f_9 * lsg1_72[k]
                   + f_3 * pc_y[k] * lsh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, ksh_60, ksh_61, ksh_62, lsg0_73, lsg0_74, \
                         lsg1_73, lsg1_74, lsh_102, lsh_103, lsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * ksh_60[k]
                   + f_6 * lsg0_73[k]
                   - f_7 * lsg1_73[k]
                   + f_3 * pc_y[k] * lsh_102[k];

        t_137[k] = f_11 * ksh_61[k]
                   + f_4 * lsg0_74[k]
                   - f_5 * lsg1_74[k]
                   + f_3 * pc_y[k] * lsh_103[k];

        t_138[k] = f_11 * ksh_62[k]
                   + f_3 * pc_y[k] * lsh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, ksi0_83, ksh_42, \
                         ksh_105, ksi1_83, lsg0_75, lsg1_75, lsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * ksi0_83[k]
                   - f_10 * pc_y[k] * ksi1_83[k];

        t_140[k] = f_18 * ksh_105[k]
                   + f_1 * lsg0_75[k]
                   - f_2 * lsg1_75[k]
                   + f_3 * pc_x[k] * lsh_105[k];

        t_141[k] = f_3 * pc_y[k] * lsh_105[k];

        t_142[k] = f_12 * ksh_42[k]
                   + f_3 * pc_z[k] * lsh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, ksh_110, lsg0_75, lsg0_80, lsg1_75, \
                         lsg1_80, lsh_106, lsh_107, lsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * lsg0_75[k]
                   - f_5 * lsg1_75[k]
                   + f_3 * pc_y[k] * lsh_106[k];

        t_144[k] = f_3 * pc_y[k] * lsh_107[k];

        t_145[k] = f_18 * ksh_110[k]
                   + f_8 * lsg0_80[k]
                   - f_9 * lsg1_80[k]
                   + f_3 * pc_x[k] * lsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, lsg0_76, lsg0_77, lsg1_76, lsg1_77, \
                         lsh_108, lsh_109, lsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * lsg0_76[k]
                   - f_7 * lsg1_76[k]
                   + f_3 * pc_y[k] * lsh_108[k];

        t_147[k] = f_4 * lsg0_77[k]
                   - f_5 * lsg1_77[k]
                   + f_3 * pc_y[k] * lsh_109[k];

        t_148[k] = f_3 * pc_y[k] * lsh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, ksh_114, lsg0_78, lsg0_79, lsg0_84, \
                         lsg1_78, lsg1_79, lsg1_84, lsh_111, lsh_112, \
                         lsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * ksh_114[k]
                   + f_6 * lsg0_84[k]
                   - f_7 * lsg1_84[k]
                   + f_3 * pc_x[k] * lsh_114[k];

        t_150[k] = f_8 * lsg0_78[k]
                   - f_9 * lsg1_78[k]
                   + f_3 * pc_y[k] * lsh_111[k];

        t_151[k] = f_6 * lsg0_79[k]
                   - f_7 * lsg1_79[k]
                   + f_3 * pc_y[k] * lsh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, ksh_119, ksh_120, lsg0_80, \
                         lsg0_89, lsg1_80, lsg1_89, lsh_113, lsh_114, lsh_119, \
                         lsh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * lsg0_80[k]
                   - f_5 * lsg1_80[k]
                   + f_3 * pc_y[k] * lsh_113[k];

        t_153[k] = f_3 * pc_y[k] * lsh_114[k];

        t_154[k] = f_18 * ksh_119[k]
                   + f_4 * lsg0_89[k]
                   - f_5 * lsg1_89[k]
                   + f_3 * pc_x[k] * lsh_119[k];

        t_155[k] = f_18 * ksh_120[k]
                   + f_3 * pc_x[k] * lsh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, ksh_121, ksh_122, \
                         ksh_123, ksh_125, lsh_119, lsh_121, lsh_122, lsh_123, \
                         lsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * ksh_121[k]
                   + f_3 * pc_x[k] * lsh_121[k];

        t_157[k] = f_18 * ksh_122[k]
                   + f_3 * pc_x[k] * lsh_122[k];

        t_158[k] = f_18 * ksh_123[k]
                   + f_3 * pc_x[k] * lsh_123[k];

        t_159[k] = f_3 * pc_y[k] * lsh_119[k];

        t_160[k] = f_18 * ksh_125[k]
                   + f_3 * pc_x[k] * lsh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, lsg0_85, lsg0_86, lsg0_87, lsg1_85, \
                         lsg1_86, lsg1_87, lsh_120, lsh_121, lsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * lsg0_85[k]
                   - f_2 * lsg1_85[k]
                   + f_3 * pc_y[k] * lsh_120[k];

        t_162[k] = f_16 * lsg0_86[k]
                   - f_17 * lsg1_86[k]
                   + f_3 * pc_y[k] * lsh_121[k];

        t_163[k] = f_8 * lsg0_87[k]
                   - f_9 * lsg1_87[k]
                   + f_3 * pc_y[k] * lsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, ksh_62, lsg0_88, lsg0_89, \
                         lsg1_88, lsg1_89, lsh_123, lsh_124, lsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * lsg0_88[k]
                   - f_7 * lsg1_88[k]
                   + f_3 * pc_y[k] * lsh_123[k];

        t_165[k] = f_4 * lsg0_89[k]
                   - f_5 * lsg1_89[k]
                   + f_3 * pc_y[k] * lsh_124[k];

        t_166[k] = f_3 * pc_y[k] * lsh_125[k];

        t_167[k] = f_12 * ksh_62[k]
                   + f_1 * lsg0_89[k]
                   - f_2 * lsg1_89[k]
                   + f_3 * pc_z[k] * lsh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, ksh_63, ksh_126, \
                         ksh_129, lsg0_90, lsg0_93, lsg1_90, lsg1_93, lsh_126, \
                         lsh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_19 * ksh_126[k]
                   + f_1 * lsg0_90[k]
                   - f_2 * lsg1_90[k]
                   + f_3 * pc_x[k] * lsh_126[k];

        t_169[k] = f_13 * ksh_63[k]
                   + f_3 * pc_y[k] * lsh_126[k];

        t_170[k] = f_3 * pc_z[k] * lsh_126[k];

        t_171[k] = f_19 * ksh_129[k]
                   + f_8 * lsg0_93[k]
                   - f_9 * lsg1_93[k]
                   + f_3 * pc_x[k] * lsh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, ksh_132, lsg0_90, lsg0_96, \
                         lsg1_90, lsg1_96, lsh_127, lsh_128, lsh_129, \
                         lsh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * lsh_127[k];

        t_173[k] = f_4 * lsg0_90[k]
                   - f_5 * lsg1_90[k]
                   + f_3 * pc_z[k] * lsh_128[k];

        t_174[k] = f_19 * ksh_132[k]
                   + f_6 * lsg0_96[k]
                   - f_7 * lsg1_96[k]
                   + f_3 * pc_x[k] * lsh_132[k];

        t_175[k] = f_3 * pc_z[k] * lsh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, ksh_68, ksh_136, \
                         lsg0_92, lsg0_100, lsg1_92, lsg1_100, lsh_131, lsh_132, \
                         lsh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * ksh_68[k]
                   + f_3 * pc_y[k] * lsh_131[k];

        t_177[k] = f_6 * lsg0_92[k]
                   - f_7 * lsg1_92[k]
                   + f_3 * pc_z[k] * lsh_131[k];

        t_178[k] = f_19 * ksh_136[k]
                   + f_4 * lsg0_100[k]
                   - f_5 * lsg1_100[k]
                   + f_3 * pc_x[k] * lsh_136[k];

        t_179[k] = f_3 * pc_z[k] * lsh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, ksh_72, ksh_141, \
                         lsg0_93, lsg0_95, lsg1_93, lsg1_95, lsh_133, lsh_135, \
                         lsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * lsg0_93[k]
                   - f_5 * lsg1_93[k]
                   + f_3 * pc_z[k] * lsh_133[k];

        t_181[k] = f_13 * ksh_72[k]
                   + f_3 * pc_y[k] * lsh_135[k];

        t_182[k] = f_8 * lsg0_95[k]
                   - f_9 * lsg1_95[k]
                   + f_3 * pc_z[k] * lsh_135[k];

        t_183[k] = f_19 * ksh_141[k]
                   + f_3 * pc_x[k] * lsh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, ksh_143, ksh_144, \
                         ksh_145, ksh_146, lsh_136, lsh_143, lsh_144, lsh_145, \
                         lsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * lsh_136[k];

        t_185[k] = f_19 * ksh_143[k]
                   + f_3 * pc_x[k] * lsh_143[k];

        t_186[k] = f_19 * ksh_144[k]
                   + f_3 * pc_x[k] * lsh_144[k];

        t_187[k] = f_19 * ksh_145[k]
                   + f_3 * pc_x[k] * lsh_145[k];

        t_188[k] = f_19 * ksh_146[k]
                   + f_3 * pc_x[k] * lsh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, ksh_78, lsg0_100, lsg0_101, \
                         lsg1_100, lsg1_101, lsh_141, lsh_142, \
                         lsh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * ksh_78[k]
                   + f_1 * lsg0_100[k]
                   - f_2 * lsg1_100[k]
                   + f_3 * pc_y[k] * lsh_141[k];

        t_190[k] = f_3 * pc_z[k] * lsh_141[k];

        t_191[k] = f_4 * lsg0_100[k]
                   - f_5 * lsg1_100[k]
                   + f_3 * pc_z[k] * lsh_142[k];

        t_192[k] = f_6 * lsg0_101[k]
                   - f_7 * lsg1_101[k]
                   + f_3 * pc_z[k] * lsh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, ksi0_84, ksh_83, \
                         ksi1_84, lsg0_102, lsg0_104, lsg1_102, lsg1_104, lsh_144, \
                         lsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * lsg0_102[k]
                   - f_9 * lsg1_102[k]
                   + f_3 * pc_z[k] * lsh_144[k];

        t_194[k] = f_13 * ksh_83[k]
                   + f_3 * pc_y[k] * lsh_146[k];

        t_195[k] = f_1 * lsg0_104[k]
                   - f_2 * lsg1_104[k]
                   + f_3 * pc_z[k] * lsh_146[k];

        t_196[k] = pa_z[k] * ksi0_84[k]
                   - f_10 * pc_z[k] * ksi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, ksi0_87, ksh_63, \
                         ksh_84, ksh_86, ksi1_87, lsh_147, lsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * ksh_84[k]
                   + f_3 * pc_y[k] * lsh_147[k];

        t_198[k] = f_11 * ksh_63[k]
                   + f_3 * pc_z[k] * lsh_147[k];

        t_199[k] = pa_z[k] * ksi0_87[k]
                   - f_10 * pc_z[k] * ksi1_87[k];

        t_200[k] = f_12 * ksh_86[k]
                   + f_3 * pc_y[k] * lsh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, ksi0_90, ksh_66, ksh_152, \
                         ksi1_90, lsg0_110, lsg1_110, lsh_150, \
                         lsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_19 * ksh_152[k]
                   + f_8 * lsg0_110[k]
                   - f_9 * lsg1_110[k]
                   + f_3 * pc_x[k] * lsh_152[k];

        t_202[k] = pa_z[k] * ksi0_90[k]
                   - f_10 * pc_z[k] * ksi1_90[k];

        t_203[k] = f_11 * ksh_66[k]
                   + f_3 * pc_z[k] * lsh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, ksi0_94, ksh_89, \
                         ksh_156, ksi1_94, lsg0_114, lsg1_114, lsh_152, \
                         lsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * ksh_89[k]
                   + f_3 * pc_y[k] * lsh_152[k];

        t_205[k] = f_19 * ksh_156[k]
                   + f_6 * lsg0_114[k]
                   - f_7 * lsg1_114[k]
                   + f_3 * pc_x[k] * lsh_156[k];

        t_206[k] = pa_z[k] * ksi0_94[k]
                   - f_10 * pc_z[k] * ksi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, ksi0_96, ksh_69, ksh_70, \
                         ksh_93, ksi1_96, lsh_153, lsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * ksh_69[k]
                   + f_3 * pc_z[k] * lsh_153[k];

        t_208[k] = pa_z[k] * ksi0_96[k]
                   + f_12 * ksh_70[k]
                   - f_10 * pc_z[k] * ksi1_96[k];

        t_209[k] = f_12 * ksh_93[k]
                   + f_3 * pc_y[k] * lsh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, ksh_161, ksh_162, ksh_163, ksh_164, \
                         lsg0_119, lsg1_119, lsh_161, lsh_162, lsh_163, \
                         lsh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_19 * ksh_161[k]
                   + f_4 * lsg0_119[k]
                   - f_5 * lsg1_119[k]
                   + f_3 * pc_x[k] * lsh_161[k];

        t_211[k] = f_19 * ksh_162[k]
                   + f_3 * pc_x[k] * lsh_162[k];

        t_212[k] = f_19 * ksh_163[k]
                   + f_3 * pc_x[k] * lsh_163[k];

        t_213[k] = f_19 * ksh_164[k]
                   + f_3 * pc_x[k] * lsh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, ksi0_105, ksh_165, \
                         ksh_166, ksh_167, ksi1_105, lsh_165, lsh_166, \
                         lsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_19 * ksh_165[k]
                   + f_3 * pc_x[k] * lsh_165[k];

        t_215[k] = f_19 * ksh_166[k]
                   + f_3 * pc_x[k] * lsh_166[k];

        t_216[k] = f_19 * ksh_167[k]
                   + f_3 * pc_x[k] * lsh_167[k];

        t_217[k] = pa_z[k] * ksi0_105[k]
                   - f_10 * pc_z[k] * ksi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, ksh_78, ksh_101, ksh_102, lsg0_117, \
                         lsg0_118, lsg1_117, lsg1_118, lsh_162, lsh_164, \
                         lsh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * ksh_78[k]
                   + f_3 * pc_z[k] * lsh_162[k];

        t_219[k] = f_12 * ksh_101[k]
                   + f_8 * lsg0_117[k]
                   - f_9 * lsg1_117[k]
                   + f_3 * pc_y[k] * lsh_164[k];

        t_220[k] = f_12 * ksh_102[k]
                   + f_6 * lsg0_118[k]
                   - f_7 * lsg1_118[k]
                   + f_3 * pc_y[k] * lsh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, ksi0_140, ksh_83, \
                         ksh_103, ksh_104, ksi1_140, lsg0_119, lsg1_119, lsh_166, \
                         lsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * ksh_103[k]
                   + f_4 * lsg0_119[k]
                   - f_5 * lsg1_119[k]
                   + f_3 * pc_y[k] * lsh_166[k];

        t_222[k] = f_12 * ksh_104[k]
                   + f_3 * pc_y[k] * lsh_167[k];

        t_223[k] = f_11 * ksh_83[k]
                   + f_1 * lsg0_119[k]
                   - f_2 * lsg1_119[k]
                   + f_3 * pc_z[k] * lsh_167[k];

        t_224[k] = pa_y[k] * ksi0_140[k]
                   - f_10 * pc_y[k] * ksi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, ksi0_143, ksh_84, \
                         ksh_105, ksh_106, ksh_107, ksi1_143, lsh_168, \
                         lsh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * ksh_105[k]
                   + f_3 * pc_y[k] * lsh_168[k];

        t_226[k] = f_12 * ksh_84[k]
                   + f_3 * pc_z[k] * lsh_168[k];

        t_227[k] = pa_y[k] * ksi0_143[k]
                   + f_12 * ksh_106[k]
                   - f_10 * pc_y[k] * ksi1_143[k];

        t_228[k] = f_11 * ksh_107[k]
                   + f_3 * pc_y[k] * lsh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, ksi0_145, ksi0_146, \
                         ksh_87, ksh_108, ksh_110, ksi1_145, ksi1_146, lsh_171, \
                         lsh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * ksi0_145[k]
                   - f_10 * pc_y[k] * ksi1_145[k];

        t_230[k] = pa_y[k] * ksi0_146[k]
                   + f_13 * ksh_108[k]
                   - f_10 * pc_y[k] * ksi1_146[k];

        t_231[k] = f_12 * ksh_87[k]
                   + f_3 * pc_z[k] * lsh_171[k];

        t_232[k] = f_11 * ksh_110[k]
                   + f_3 * pc_y[k] * lsh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, ksi0_149, ksi0_150, ksh_90, \
                         ksh_111, ksi1_149, ksi1_150, lsh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * ksi0_149[k]
                   - f_10 * pc_y[k] * ksi1_149[k];

        t_234[k] = pa_y[k] * ksi0_150[k]
                   + f_14 * ksh_111[k]
                   - f_10 * pc_y[k] * ksi1_150[k];

        t_235[k] = f_12 * ksh_90[k]
                   + f_3 * pc_z[k] * lsh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, ksi0_152, ksi0_154, \
                         ksh_113, ksh_114, ksh_183, ksi1_152, ksi1_154, lsh_177, \
                         lsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * ksi0_152[k]
                   + f_12 * ksh_113[k]
                   - f_10 * pc_y[k] * ksi1_152[k];

        t_237[k] = f_11 * ksh_114[k]
                   + f_3 * pc_y[k] * lsh_177[k];

        t_238[k] = pa_y[k] * ksi0_154[k]
                   - f_10 * pc_y[k] * ksi1_154[k];

        t_239[k] = f_19 * ksh_183[k]
                   + f_3 * pc_x[k] * lsh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, ksh_184, ksh_185, ksh_186, \
                         ksh_187, ksh_188, lsh_184, lsh_185, lsh_186, lsh_187, \
                         lsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_19 * ksh_184[k]
                   + f_3 * pc_x[k] * lsh_184[k];

        t_241[k] = f_19 * ksh_185[k]
                   + f_3 * pc_x[k] * lsh_185[k];

        t_242[k] = f_19 * ksh_186[k]
                   + f_3 * pc_x[k] * lsh_186[k];

        t_243[k] = f_19 * ksh_187[k]
                   + f_3 * pc_x[k] * lsh_187[k];

        t_244[k] = f_19 * ksh_188[k]
                   + f_3 * pc_x[k] * lsh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, ksh_99, ksh_120, ksh_122, lsg0_130, \
                         lsg0_132, lsg1_130, lsg1_132, lsh_183, \
                         lsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * ksh_120[k]
                   + f_1 * lsg0_130[k]
                   - f_2 * lsg1_130[k]
                   + f_3 * pc_y[k] * lsh_183[k];

        t_246[k] = f_12 * ksh_99[k]
                   + f_3 * pc_z[k] * lsh_183[k];

        t_247[k] = f_11 * ksh_122[k]
                   + f_8 * lsg0_132[k]
                   - f_9 * lsg1_132[k]
                   + f_3 * pc_y[k] * lsh_185[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_167 = buffer.data(ksi0 + 167);
    const auto *ksi0_168 = buffer.data(ksi0 + 168);
    const auto *ksi0_171 = buffer.data(ksi0 + 171);
    const auto *ksi0_174 = buffer.data(ksi0 + 174);
    const auto *ksi0_178 = buffer.data(ksi0 + 178);
    const auto *ksi0_180 = buffer.data(ksi0 + 180);
    const auto *ksi0_189 = buffer.data(ksi0 + 189);

    const auto *ksh_105 = buffer.data(ksh + 105);
    const auto *ksh_123 = buffer.data(ksh + 123);
    const auto *ksh_124 = buffer.data(ksh + 124);
    const auto *ksh_125 = buffer.data(ksh + 125);
    const auto *ksh_126 = buffer.data(ksh + 126);
    const auto *ksh_129 = buffer.data(ksh + 129);
    const auto *ksh_131 = buffer.data(ksh + 131);
    const auto *ksh_132 = buffer.data(ksh + 132);
    const auto *ksh_133 = buffer.data(ksh + 133);
    const auto *ksh_135 = buffer.data(ksh + 135);
    const auto *ksh_141 = buffer.data(ksh + 141);
    const auto *ksh_146 = buffer.data(ksh + 146);
    const auto *ksh_147 = buffer.data(ksh + 147);
    const auto *ksh_149 = buffer.data(ksh + 149);
    const auto *ksh_150 = buffer.data(ksh + 150);
    const auto *ksh_152 = buffer.data(ksh + 152);
    const auto *ksh_153 = buffer.data(ksh + 153);
    const auto *ksh_156 = buffer.data(ksh + 156);
    const auto *ksh_162 = buffer.data(ksh + 162);
    const auto *ksh_164 = buffer.data(ksh + 164);
    const auto *ksh_165 = buffer.data(ksh + 165);
    const auto *ksh_166 = buffer.data(ksh + 166);
    const auto *ksh_167 = buffer.data(ksh + 167);
    const auto *ksh_168 = buffer.data(ksh + 168);
    const auto *ksh_170 = buffer.data(ksh + 170);
    const auto *ksh_173 = buffer.data(ksh + 173);
    const auto *ksh_177 = buffer.data(ksh + 177);
    const auto *ksh_183 = buffer.data(ksh + 183);
    const auto *ksh_185 = buffer.data(ksh + 185);
    const auto *ksh_186 = buffer.data(ksh + 186);
    const auto *ksh_187 = buffer.data(ksh + 187);
    const auto *ksh_189 = buffer.data(ksh + 189);
    const auto *ksh_194 = buffer.data(ksh + 194);
    const auto *ksh_198 = buffer.data(ksh + 198);
    const auto *ksh_203 = buffer.data(ksh + 203);
    const auto *ksh_204 = buffer.data(ksh + 204);
    const auto *ksh_205 = buffer.data(ksh + 205);
    const auto *ksh_206 = buffer.data(ksh + 206);
    const auto *ksh_207 = buffer.data(ksh + 207);
    const auto *ksh_209 = buffer.data(ksh + 209);
    const auto *ksh_210 = buffer.data(ksh + 210);
    const auto *ksh_213 = buffer.data(ksh + 213);
    const auto *ksh_216 = buffer.data(ksh + 216);
    const auto *ksh_220 = buffer.data(ksh + 220);
    const auto *ksh_225 = buffer.data(ksh + 225);
    const auto *ksh_227 = buffer.data(ksh + 227);
    const auto *ksh_228 = buffer.data(ksh + 228);
    const auto *ksh_229 = buffer.data(ksh + 229);
    const auto *ksh_230 = buffer.data(ksh + 230);
    const auto *ksh_236 = buffer.data(ksh + 236);
    const auto *ksh_240 = buffer.data(ksh + 240);
    const auto *ksh_245 = buffer.data(ksh + 245);
    const auto *ksh_246 = buffer.data(ksh + 246);
    const auto *ksh_247 = buffer.data(ksh + 247);
    const auto *ksh_248 = buffer.data(ksh + 248);
    const auto *ksh_249 = buffer.data(ksh + 249);
    const auto *ksh_250 = buffer.data(ksh + 250);
    const auto *ksh_251 = buffer.data(ksh + 251);
    const auto *ksh_252 = buffer.data(ksh + 252);
    const auto *ksh_255 = buffer.data(ksh + 255);
    const auto *ksh_257 = buffer.data(ksh + 257);
    const auto *ksh_258 = buffer.data(ksh + 258);
    const auto *ksh_261 = buffer.data(ksh + 261);
    const auto *ksh_262 = buffer.data(ksh + 262);
    const auto *ksh_264 = buffer.data(ksh + 264);
    const auto *ksh_266 = buffer.data(ksh + 266);
    const auto *ksh_267 = buffer.data(ksh + 267);
    const auto *ksh_268 = buffer.data(ksh + 268);
    const auto *ksh_269 = buffer.data(ksh + 269);
    const auto *ksh_270 = buffer.data(ksh + 270);
    const auto *ksh_271 = buffer.data(ksh + 271);
    const auto *ksh_272 = buffer.data(ksh + 272);

    const auto *ksi1_167 = buffer.data(ksi1 + 167);
    const auto *ksi1_168 = buffer.data(ksi1 + 168);
    const auto *ksi1_171 = buffer.data(ksi1 + 171);
    const auto *ksi1_174 = buffer.data(ksi1 + 174);
    const auto *ksi1_178 = buffer.data(ksi1 + 178);
    const auto *ksi1_180 = buffer.data(ksi1 + 180);
    const auto *ksi1_189 = buffer.data(ksi1 + 189);

    const auto *lsg0_133 = buffer.data(lsg0 + 133);
    const auto *lsg0_134 = buffer.data(lsg0 + 134);
    const auto *lsg0_135 = buffer.data(lsg0 + 135);
    const auto *lsg0_136 = buffer.data(lsg0 + 136);
    const auto *lsg0_137 = buffer.data(lsg0 + 137);
    const auto *lsg0_138 = buffer.data(lsg0 + 138);
    const auto *lsg0_139 = buffer.data(lsg0 + 139);
    const auto *lsg0_140 = buffer.data(lsg0 + 140);
    const auto *lsg0_144 = buffer.data(lsg0 + 144);
    const auto *lsg0_145 = buffer.data(lsg0 + 145);
    const auto *lsg0_146 = buffer.data(lsg0 + 146);
    const auto *lsg0_147 = buffer.data(lsg0 + 147);
    const auto *lsg0_148 = buffer.data(lsg0 + 148);
    const auto *lsg0_149 = buffer.data(lsg0 + 149);
    const auto *lsg0_150 = buffer.data(lsg0 + 150);
    const auto *lsg0_152 = buffer.data(lsg0 + 152);
    const auto *lsg0_153 = buffer.data(lsg0 + 153);
    const auto *lsg0_155 = buffer.data(lsg0 + 155);
    const auto *lsg0_156 = buffer.data(lsg0 + 156);
    const auto *lsg0_160 = buffer.data(lsg0 + 160);
    const auto *lsg0_161 = buffer.data(lsg0 + 161);
    const auto *lsg0_162 = buffer.data(lsg0 + 162);
    const auto *lsg0_164 = buffer.data(lsg0 + 164);
    const auto *lsg0_170 = buffer.data(lsg0 + 170);
    const auto *lsg0_174 = buffer.data(lsg0 + 174);
    const auto *lsg0_177 = buffer.data(lsg0 + 177);
    const auto *lsg0_178 = buffer.data(lsg0 + 178);
    const auto *lsg0_179 = buffer.data(lsg0 + 179);
    const auto *lsg0_180 = buffer.data(lsg0 + 180);
    const auto *lsg0_183 = buffer.data(lsg0 + 183);
    const auto *lsg0_185 = buffer.data(lsg0 + 185);
    const auto *lsg0_186 = buffer.data(lsg0 + 186);
    const auto *lsg0_189 = buffer.data(lsg0 + 189);
    const auto *lsg0_190 = buffer.data(lsg0 + 190);
    const auto *lsg0_192 = buffer.data(lsg0 + 192);
    const auto *lsg0_193 = buffer.data(lsg0 + 193);
    const auto *lsg0_194 = buffer.data(lsg0 + 194);

    const auto *lsg1_133 = buffer.data(lsg1 + 133);
    const auto *lsg1_134 = buffer.data(lsg1 + 134);
    const auto *lsg1_135 = buffer.data(lsg1 + 135);
    const auto *lsg1_136 = buffer.data(lsg1 + 136);
    const auto *lsg1_137 = buffer.data(lsg1 + 137);
    const auto *lsg1_138 = buffer.data(lsg1 + 138);
    const auto *lsg1_139 = buffer.data(lsg1 + 139);
    const auto *lsg1_140 = buffer.data(lsg1 + 140);
    const auto *lsg1_144 = buffer.data(lsg1 + 144);
    const auto *lsg1_145 = buffer.data(lsg1 + 145);
    const auto *lsg1_146 = buffer.data(lsg1 + 146);
    const auto *lsg1_147 = buffer.data(lsg1 + 147);
    const auto *lsg1_148 = buffer.data(lsg1 + 148);
    const auto *lsg1_149 = buffer.data(lsg1 + 149);
    const auto *lsg1_150 = buffer.data(lsg1 + 150);
    const auto *lsg1_152 = buffer.data(lsg1 + 152);
    const auto *lsg1_153 = buffer.data(lsg1 + 153);
    const auto *lsg1_155 = buffer.data(lsg1 + 155);
    const auto *lsg1_156 = buffer.data(lsg1 + 156);
    const auto *lsg1_160 = buffer.data(lsg1 + 160);
    const auto *lsg1_161 = buffer.data(lsg1 + 161);
    const auto *lsg1_162 = buffer.data(lsg1 + 162);
    const auto *lsg1_164 = buffer.data(lsg1 + 164);
    const auto *lsg1_170 = buffer.data(lsg1 + 170);
    const auto *lsg1_174 = buffer.data(lsg1 + 174);
    const auto *lsg1_177 = buffer.data(lsg1 + 177);
    const auto *lsg1_178 = buffer.data(lsg1 + 178);
    const auto *lsg1_179 = buffer.data(lsg1 + 179);
    const auto *lsg1_180 = buffer.data(lsg1 + 180);
    const auto *lsg1_183 = buffer.data(lsg1 + 183);
    const auto *lsg1_185 = buffer.data(lsg1 + 185);
    const auto *lsg1_186 = buffer.data(lsg1 + 186);
    const auto *lsg1_189 = buffer.data(lsg1 + 189);
    const auto *lsg1_190 = buffer.data(lsg1 + 190);
    const auto *lsg1_192 = buffer.data(lsg1 + 192);
    const auto *lsg1_193 = buffer.data(lsg1 + 193);
    const auto *lsg1_194 = buffer.data(lsg1 + 194);

    const auto *lsh_186 = buffer.data(lsh + 186);
    const auto *lsh_187 = buffer.data(lsh + 187);
    const auto *lsh_188 = buffer.data(lsh + 188);
    const auto *lsh_189 = buffer.data(lsh + 189);
    const auto *lsh_190 = buffer.data(lsh + 190);
    const auto *lsh_191 = buffer.data(lsh + 191);
    const auto *lsh_192 = buffer.data(lsh + 192);
    const auto *lsh_193 = buffer.data(lsh + 193);
    const auto *lsh_194 = buffer.data(lsh + 194);
    const auto *lsh_195 = buffer.data(lsh + 195);
    const auto *lsh_196 = buffer.data(lsh + 196);
    const auto *lsh_197 = buffer.data(lsh + 197);
    const auto *lsh_198 = buffer.data(lsh + 198);
    const auto *lsh_203 = buffer.data(lsh + 203);
    const auto *lsh_204 = buffer.data(lsh + 204);
    const auto *lsh_205 = buffer.data(lsh + 205);
    const auto *lsh_206 = buffer.data(lsh + 206);
    const auto *lsh_207 = buffer.data(lsh + 207);
    const auto *lsh_208 = buffer.data(lsh + 208);
    const auto *lsh_209 = buffer.data(lsh + 209);
    const auto *lsh_210 = buffer.data(lsh + 210);
    const auto *lsh_211 = buffer.data(lsh + 211);
    const auto *lsh_212 = buffer.data(lsh + 212);
    const auto *lsh_213 = buffer.data(lsh + 213);
    const auto *lsh_215 = buffer.data(lsh + 215);
    const auto *lsh_216 = buffer.data(lsh + 216);
    const auto *lsh_217 = buffer.data(lsh + 217);
    const auto *lsh_219 = buffer.data(lsh + 219);
    const auto *lsh_220 = buffer.data(lsh + 220);
    const auto *lsh_225 = buffer.data(lsh + 225);
    const auto *lsh_226 = buffer.data(lsh + 226);
    const auto *lsh_227 = buffer.data(lsh + 227);
    const auto *lsh_228 = buffer.data(lsh + 228);
    const auto *lsh_229 = buffer.data(lsh + 229);
    const auto *lsh_230 = buffer.data(lsh + 230);
    const auto *lsh_231 = buffer.data(lsh + 231);
    const auto *lsh_233 = buffer.data(lsh + 233);
    const auto *lsh_234 = buffer.data(lsh + 234);
    const auto *lsh_236 = buffer.data(lsh + 236);
    const auto *lsh_237 = buffer.data(lsh + 237);
    const auto *lsh_240 = buffer.data(lsh + 240);
    const auto *lsh_245 = buffer.data(lsh + 245);
    const auto *lsh_246 = buffer.data(lsh + 246);
    const auto *lsh_247 = buffer.data(lsh + 247);
    const auto *lsh_248 = buffer.data(lsh + 248);
    const auto *lsh_249 = buffer.data(lsh + 249);
    const auto *lsh_250 = buffer.data(lsh + 250);
    const auto *lsh_251 = buffer.data(lsh + 251);
    const auto *lsh_252 = buffer.data(lsh + 252);
    const auto *lsh_254 = buffer.data(lsh + 254);
    const auto *lsh_255 = buffer.data(lsh + 255);
    const auto *lsh_257 = buffer.data(lsh + 257);
    const auto *lsh_258 = buffer.data(lsh + 258);
    const auto *lsh_261 = buffer.data(lsh + 261);
    const auto *lsh_262 = buffer.data(lsh + 262);
    const auto *lsh_264 = buffer.data(lsh + 264);
    const auto *lsh_266 = buffer.data(lsh + 266);
    const auto *lsh_267 = buffer.data(lsh + 267);
    const auto *lsh_268 = buffer.data(lsh + 268);
    const auto *lsh_269 = buffer.data(lsh + 269);
    const auto *lsh_270 = buffer.data(lsh + 270);
    const auto *lsh_271 = buffer.data(lsh + 271);
    const auto *lsh_272 = buffer.data(lsh + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, ksh_123, ksh_124, ksh_125, lsg0_133, \
                         lsg0_134, lsg1_133, lsg1_134, lsh_186, lsh_187, \
                         lsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * ksh_123[k]
                   + f_6 * lsg0_133[k]
                   - f_7 * lsg1_133[k]
                   + f_3 * pc_y[k] * lsh_186[k];

        t_249[k] = f_11 * ksh_124[k]
                   + f_4 * lsg0_134[k]
                   - f_5 * lsg1_134[k]
                   + f_3 * pc_y[k] * lsh_187[k];

        t_250[k] = f_11 * ksh_125[k]
                   + f_3 * pc_y[k] * lsh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, ksi0_167, \
                         ksh_105, ksh_189, ksi1_167, lsg0_135, lsg1_135, \
                         lsh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * ksi0_167[k]
                   - f_10 * pc_y[k] * ksi1_167[k];

        t_252[k] = f_19 * ksh_189[k]
                   + f_1 * lsg0_135[k]
                   - f_2 * lsg1_135[k]
                   + f_3 * pc_x[k] * lsh_189[k];

        t_253[k] = f_3 * pc_y[k] * lsh_189[k];

        t_254[k] = f_13 * ksh_105[k]
                   + f_3 * pc_z[k] * lsh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, ksh_194, lsg0_135, lsg0_140, \
                         lsg1_135, lsg1_140, lsh_190, lsh_191, \
                         lsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * lsg0_135[k]
                   - f_5 * lsg1_135[k]
                   + f_3 * pc_y[k] * lsh_190[k];

        t_256[k] = f_3 * pc_y[k] * lsh_191[k];

        t_257[k] = f_19 * ksh_194[k]
                   + f_8 * lsg0_140[k]
                   - f_9 * lsg1_140[k]
                   + f_3 * pc_x[k] * lsh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, lsg0_136, lsg0_137, lsg1_136, lsg1_137, \
                         lsh_192, lsh_193, lsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * lsg0_136[k]
                   - f_7 * lsg1_136[k]
                   + f_3 * pc_y[k] * lsh_192[k];

        t_259[k] = f_4 * lsg0_137[k]
                   - f_5 * lsg1_137[k]
                   + f_3 * pc_y[k] * lsh_193[k];

        t_260[k] = f_3 * pc_y[k] * lsh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, ksh_198, lsg0_138, lsg0_139, \
                         lsg0_144, lsg1_138, lsg1_139, lsg1_144, lsh_195, lsh_196, \
                         lsh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_19 * ksh_198[k]
                   + f_6 * lsg0_144[k]
                   - f_7 * lsg1_144[k]
                   + f_3 * pc_x[k] * lsh_198[k];

        t_262[k] = f_8 * lsg0_138[k]
                   - f_9 * lsg1_138[k]
                   + f_3 * pc_y[k] * lsh_195[k];

        t_263[k] = f_6 * lsg0_139[k]
                   - f_7 * lsg1_139[k]
                   + f_3 * pc_y[k] * lsh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, ksh_203, ksh_204, lsg0_140, \
                         lsg0_149, lsg1_140, lsg1_149, lsh_197, lsh_198, lsh_203, \
                         lsh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * lsg0_140[k]
                   - f_5 * lsg1_140[k]
                   + f_3 * pc_y[k] * lsh_197[k];

        t_265[k] = f_3 * pc_y[k] * lsh_198[k];

        t_266[k] = f_19 * ksh_203[k]
                   + f_4 * lsg0_149[k]
                   - f_5 * lsg1_149[k]
                   + f_3 * pc_x[k] * lsh_203[k];

        t_267[k] = f_19 * ksh_204[k]
                   + f_3 * pc_x[k] * lsh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, ksh_205, ksh_206, \
                         ksh_207, ksh_209, lsh_203, lsh_205, lsh_206, lsh_207, \
                         lsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_19 * ksh_205[k]
                   + f_3 * pc_x[k] * lsh_205[k];

        t_269[k] = f_19 * ksh_206[k]
                   + f_3 * pc_x[k] * lsh_206[k];

        t_270[k] = f_19 * ksh_207[k]
                   + f_3 * pc_x[k] * lsh_207[k];

        t_271[k] = f_3 * pc_y[k] * lsh_203[k];

        t_272[k] = f_19 * ksh_209[k]
                   + f_3 * pc_x[k] * lsh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, lsg0_145, lsg0_146, lsg0_147, lsg1_145, \
                         lsg1_146, lsg1_147, lsh_204, lsh_205, \
                         lsh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * lsg0_145[k]
                   - f_2 * lsg1_145[k]
                   + f_3 * pc_y[k] * lsh_204[k];

        t_274[k] = f_16 * lsg0_146[k]
                   - f_17 * lsg1_146[k]
                   + f_3 * pc_y[k] * lsh_205[k];

        t_275[k] = f_8 * lsg0_147[k]
                   - f_9 * lsg1_147[k]
                   + f_3 * pc_y[k] * lsh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, ksh_125, lsg0_148, lsg0_149, \
                         lsg1_148, lsg1_149, lsh_207, lsh_208, \
                         lsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * lsg0_148[k]
                   - f_7 * lsg1_148[k]
                   + f_3 * pc_y[k] * lsh_207[k];

        t_277[k] = f_4 * lsg0_149[k]
                   - f_5 * lsg1_149[k]
                   + f_3 * pc_y[k] * lsh_208[k];

        t_278[k] = f_3 * pc_y[k] * lsh_209[k];

        t_279[k] = f_13 * ksh_125[k]
                   + f_1 * lsg0_149[k]
                   - f_2 * lsg1_149[k]
                   + f_3 * pc_z[k] * lsh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, ksh_126, ksh_210, \
                         ksh_213, lsg0_150, lsg0_153, lsg1_150, lsg1_153, lsh_210, \
                         lsh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * ksh_210[k]
                   + f_1 * lsg0_150[k]
                   - f_2 * lsg1_150[k]
                   + f_3 * pc_x[k] * lsh_210[k];

        t_281[k] = f_14 * ksh_126[k]
                   + f_3 * pc_y[k] * lsh_210[k];

        t_282[k] = f_3 * pc_z[k] * lsh_210[k];

        t_283[k] = f_14 * ksh_213[k]
                   + f_8 * lsg0_153[k]
                   - f_9 * lsg1_153[k]
                   + f_3 * pc_x[k] * lsh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, ksh_216, lsg0_150, lsg0_156, \
                         lsg1_150, lsg1_156, lsh_211, lsh_212, lsh_213, \
                         lsh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * lsh_211[k];

        t_285[k] = f_4 * lsg0_150[k]
                   - f_5 * lsg1_150[k]
                   + f_3 * pc_z[k] * lsh_212[k];

        t_286[k] = f_14 * ksh_216[k]
                   + f_6 * lsg0_156[k]
                   - f_7 * lsg1_156[k]
                   + f_3 * pc_x[k] * lsh_216[k];

        t_287[k] = f_3 * pc_z[k] * lsh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, ksh_131, ksh_220, \
                         lsg0_152, lsg0_160, lsg1_152, lsg1_160, lsh_215, lsh_216, \
                         lsh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * ksh_131[k]
                   + f_3 * pc_y[k] * lsh_215[k];

        t_289[k] = f_6 * lsg0_152[k]
                   - f_7 * lsg1_152[k]
                   + f_3 * pc_z[k] * lsh_215[k];

        t_290[k] = f_14 * ksh_220[k]
                   + f_4 * lsg0_160[k]
                   - f_5 * lsg1_160[k]
                   + f_3 * pc_x[k] * lsh_220[k];

        t_291[k] = f_3 * pc_z[k] * lsh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, ksh_135, ksh_225, \
                         lsg0_153, lsg0_155, lsg1_153, lsg1_155, lsh_217, lsh_219, \
                         lsh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * lsg0_153[k]
                   - f_5 * lsg1_153[k]
                   + f_3 * pc_z[k] * lsh_217[k];

        t_293[k] = f_14 * ksh_135[k]
                   + f_3 * pc_y[k] * lsh_219[k];

        t_294[k] = f_8 * lsg0_155[k]
                   - f_9 * lsg1_155[k]
                   + f_3 * pc_z[k] * lsh_219[k];

        t_295[k] = f_14 * ksh_225[k]
                   + f_3 * pc_x[k] * lsh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, ksh_227, ksh_228, \
                         ksh_229, ksh_230, lsh_220, lsh_227, lsh_228, lsh_229, \
                         lsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * lsh_220[k];

        t_297[k] = f_14 * ksh_227[k]
                   + f_3 * pc_x[k] * lsh_227[k];

        t_298[k] = f_14 * ksh_228[k]
                   + f_3 * pc_x[k] * lsh_228[k];

        t_299[k] = f_14 * ksh_229[k]
                   + f_3 * pc_x[k] * lsh_229[k];

        t_300[k] = f_14 * ksh_230[k]
                   + f_3 * pc_x[k] * lsh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, ksh_141, lsg0_160, lsg0_161, \
                         lsg1_160, lsg1_161, lsh_225, lsh_226, \
                         lsh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * ksh_141[k]
                   + f_1 * lsg0_160[k]
                   - f_2 * lsg1_160[k]
                   + f_3 * pc_y[k] * lsh_225[k];

        t_302[k] = f_3 * pc_z[k] * lsh_225[k];

        t_303[k] = f_4 * lsg0_160[k]
                   - f_5 * lsg1_160[k]
                   + f_3 * pc_z[k] * lsh_226[k];

        t_304[k] = f_6 * lsg0_161[k]
                   - f_7 * lsg1_161[k]
                   + f_3 * pc_z[k] * lsh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, ksi0_168, ksh_146, \
                         ksi1_168, lsg0_162, lsg0_164, lsg1_162, lsg1_164, lsh_228, \
                         lsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * lsg0_162[k]
                   - f_9 * lsg1_162[k]
                   + f_3 * pc_z[k] * lsh_228[k];

        t_306[k] = f_14 * ksh_146[k]
                   + f_3 * pc_y[k] * lsh_230[k];

        t_307[k] = f_1 * lsg0_164[k]
                   - f_2 * lsg1_164[k]
                   + f_3 * pc_z[k] * lsh_230[k];

        t_308[k] = pa_z[k] * ksi0_168[k]
                   - f_10 * pc_z[k] * ksi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, ksi0_171, ksh_126, \
                         ksh_147, ksh_149, ksi1_171, lsh_231, lsh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * ksh_147[k]
                   + f_3 * pc_y[k] * lsh_231[k];

        t_310[k] = f_11 * ksh_126[k]
                   + f_3 * pc_z[k] * lsh_231[k];

        t_311[k] = pa_z[k] * ksi0_171[k]
                   - f_10 * pc_z[k] * ksi1_171[k];

        t_312[k] = f_13 * ksh_149[k]
                   + f_3 * pc_y[k] * lsh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, ksi0_174, ksh_129, ksh_236, \
                         ksi1_174, lsg0_170, lsg1_170, lsh_234, \
                         lsh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_14 * ksh_236[k]
                   + f_8 * lsg0_170[k]
                   - f_9 * lsg1_170[k]
                   + f_3 * pc_x[k] * lsh_236[k];

        t_314[k] = pa_z[k] * ksi0_174[k]
                   - f_10 * pc_z[k] * ksi1_174[k];

        t_315[k] = f_11 * ksh_129[k]
                   + f_3 * pc_z[k] * lsh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, ksi0_178, ksh_152, \
                         ksh_240, ksi1_178, lsg0_174, lsg1_174, lsh_236, \
                         lsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * ksh_152[k]
                   + f_3 * pc_y[k] * lsh_236[k];

        t_317[k] = f_14 * ksh_240[k]
                   + f_6 * lsg0_174[k]
                   - f_7 * lsg1_174[k]
                   + f_3 * pc_x[k] * lsh_240[k];

        t_318[k] = pa_z[k] * ksi0_178[k]
                   - f_10 * pc_z[k] * ksi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, ksi0_180, ksh_132, ksh_133, \
                         ksh_156, ksi1_180, lsh_237, lsh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * ksh_132[k]
                   + f_3 * pc_z[k] * lsh_237[k];

        t_320[k] = pa_z[k] * ksi0_180[k]
                   + f_12 * ksh_133[k]
                   - f_10 * pc_z[k] * ksi1_180[k];

        t_321[k] = f_13 * ksh_156[k]
                   + f_3 * pc_y[k] * lsh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, ksh_245, ksh_246, ksh_247, ksh_248, \
                         lsg0_179, lsg1_179, lsh_245, lsh_246, lsh_247, \
                         lsh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * ksh_245[k]
                   + f_4 * lsg0_179[k]
                   - f_5 * lsg1_179[k]
                   + f_3 * pc_x[k] * lsh_245[k];

        t_323[k] = f_14 * ksh_246[k]
                   + f_3 * pc_x[k] * lsh_246[k];

        t_324[k] = f_14 * ksh_247[k]
                   + f_3 * pc_x[k] * lsh_247[k];

        t_325[k] = f_14 * ksh_248[k]
                   + f_3 * pc_x[k] * lsh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, ksi0_189, ksh_249, \
                         ksh_250, ksh_251, ksi1_189, lsh_249, lsh_250, \
                         lsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_14 * ksh_249[k]
                   + f_3 * pc_x[k] * lsh_249[k];

        t_327[k] = f_14 * ksh_250[k]
                   + f_3 * pc_x[k] * lsh_250[k];

        t_328[k] = f_14 * ksh_251[k]
                   + f_3 * pc_x[k] * lsh_251[k];

        t_329[k] = pa_z[k] * ksi0_189[k]
                   - f_10 * pc_z[k] * ksi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, ksh_141, ksh_164, ksh_165, lsg0_177, \
                         lsg0_178, lsg1_177, lsg1_178, lsh_246, lsh_248, \
                         lsh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * ksh_141[k]
                   + f_3 * pc_z[k] * lsh_246[k];

        t_331[k] = f_13 * ksh_164[k]
                   + f_8 * lsg0_177[k]
                   - f_9 * lsg1_177[k]
                   + f_3 * pc_y[k] * lsh_248[k];

        t_332[k] = f_13 * ksh_165[k]
                   + f_6 * lsg0_178[k]
                   - f_7 * lsg1_178[k]
                   + f_3 * pc_y[k] * lsh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, ksh_146, ksh_166, ksh_167, lsg0_179, \
                         lsg1_179, lsh_250, lsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * ksh_166[k]
                   + f_4 * lsg0_179[k]
                   - f_5 * lsg1_179[k]
                   + f_3 * pc_y[k] * lsh_250[k];

        t_334[k] = f_13 * ksh_167[k]
                   + f_3 * pc_y[k] * lsh_251[k];

        t_335[k] = f_11 * ksh_146[k]
                   + f_1 * lsg0_179[k]
                   - f_2 * lsg1_179[k]
                   + f_3 * pc_z[k] * lsh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, ksh_147, ksh_168, ksh_252, \
                         lsg0_180, lsg1_180, lsh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_14 * ksh_252[k]
                   + f_1 * lsg0_180[k]
                   - f_2 * lsg1_180[k]
                   + f_3 * pc_x[k] * lsh_252[k];

        t_337[k] = f_12 * ksh_168[k]
                   + f_3 * pc_y[k] * lsh_252[k];

        t_338[k] = f_12 * ksh_147[k]
                   + f_3 * pc_z[k] * lsh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, ksh_170, ksh_255, ksh_257, lsg0_183, \
                         lsg0_185, lsg1_183, lsg1_185, lsh_254, lsh_255, \
                         lsh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_14 * ksh_255[k]
                   + f_8 * lsg0_183[k]
                   - f_9 * lsg1_183[k]
                   + f_3 * pc_x[k] * lsh_255[k];

        t_340[k] = f_12 * ksh_170[k]
                   + f_3 * pc_y[k] * lsh_254[k];

        t_341[k] = f_14 * ksh_257[k]
                   + f_8 * lsg0_185[k]
                   - f_9 * lsg1_185[k]
                   + f_3 * pc_x[k] * lsh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, ksh_150, ksh_173, ksh_258, \
                         lsg0_186, lsg1_186, lsh_255, lsh_257, \
                         lsh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_14 * ksh_258[k]
                   + f_6 * lsg0_186[k]
                   - f_7 * lsg1_186[k]
                   + f_3 * pc_x[k] * lsh_258[k];

        t_343[k] = f_12 * ksh_150[k]
                   + f_3 * pc_z[k] * lsh_255[k];

        t_344[k] = f_12 * ksh_173[k]
                   + f_3 * pc_y[k] * lsh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, ksh_153, ksh_261, ksh_262, lsg0_189, \
                         lsg0_190, lsg1_189, lsg1_190, lsh_258, lsh_261, \
                         lsh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_14 * ksh_261[k]
                   + f_6 * lsg0_189[k]
                   - f_7 * lsg1_189[k]
                   + f_3 * pc_x[k] * lsh_261[k];

        t_346[k] = f_14 * ksh_262[k]
                   + f_4 * lsg0_190[k]
                   - f_5 * lsg1_190[k]
                   + f_3 * pc_x[k] * lsh_262[k];

        t_347[k] = f_12 * ksh_153[k]
                   + f_3 * pc_z[k] * lsh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, ksh_177, ksh_264, ksh_266, lsg0_192, \
                         lsg0_194, lsg1_192, lsg1_194, lsh_261, lsh_264, \
                         lsh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_14 * ksh_264[k]
                   + f_4 * lsg0_192[k]
                   - f_5 * lsg1_192[k]
                   + f_3 * pc_x[k] * lsh_264[k];

        t_349[k] = f_12 * ksh_177[k]
                   + f_3 * pc_y[k] * lsh_261[k];

        t_350[k] = f_14 * ksh_266[k]
                   + f_4 * lsg0_194[k]
                   - f_5 * lsg1_194[k]
                   + f_3 * pc_x[k] * lsh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, ksh_267, ksh_268, ksh_269, \
                         ksh_270, ksh_271, lsh_267, lsh_268, lsh_269, lsh_270, \
                         lsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_14 * ksh_267[k]
                   + f_3 * pc_x[k] * lsh_267[k];

        t_352[k] = f_14 * ksh_268[k]
                   + f_3 * pc_x[k] * lsh_268[k];

        t_353[k] = f_14 * ksh_269[k]
                   + f_3 * pc_x[k] * lsh_269[k];

        t_354[k] = f_14 * ksh_270[k]
                   + f_3 * pc_x[k] * lsh_270[k];

        t_355[k] = f_14 * ksh_271[k]
                   + f_3 * pc_x[k] * lsh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, ksh_162, ksh_183, ksh_272, \
                         lsg0_190, lsg1_190, lsh_267, lsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_14 * ksh_272[k]
                   + f_3 * pc_x[k] * lsh_272[k];

        t_357[k] = f_12 * ksh_183[k]
                   + f_1 * lsg0_190[k]
                   - f_2 * lsg1_190[k]
                   + f_3 * pc_y[k] * lsh_267[k];

        t_358[k] = f_12 * ksh_162[k]
                   + f_3 * pc_z[k] * lsh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, ksh_185, ksh_186, ksh_187, lsg0_192, \
                         lsg0_193, lsg0_194, lsg1_192, lsg1_193, lsg1_194, lsh_269, lsh_270, \
                         lsh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * ksh_185[k]
                   + f_8 * lsg0_192[k]
                   - f_9 * lsg1_192[k]
                   + f_3 * pc_y[k] * lsh_269[k];

        t_360[k] = f_12 * ksh_186[k]
                   + f_6 * lsg0_193[k]
                   - f_7 * lsg1_193[k]
                   + f_3 * pc_y[k] * lsh_270[k];

        t_361[k] = f_12 * ksh_187[k]
                   + f_4 * lsg0_194[k]
                   - f_5 * lsg1_194[k]
                   + f_3 * pc_y[k] * lsh_271[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_252 = buffer.data(ksi0 + 252);
    const auto *ksi0_255 = buffer.data(ksi0 + 255);
    const auto *ksi0_257 = buffer.data(ksi0 + 257);
    const auto *ksi0_258 = buffer.data(ksi0 + 258);
    const auto *ksi0_261 = buffer.data(ksi0 + 261);
    const auto *ksi0_262 = buffer.data(ksi0 + 262);
    const auto *ksi0_264 = buffer.data(ksi0 + 264);
    const auto *ksi0_266 = buffer.data(ksi0 + 266);
    const auto *ksi0_279 = buffer.data(ksi0 + 279);
    const auto *ksi0_280 = buffer.data(ksi0 + 280);
    const auto *ksi0_283 = buffer.data(ksi0 + 283);
    const auto *ksi0_286 = buffer.data(ksi0 + 286);
    const auto *ksi0_290 = buffer.data(ksi0 + 290);
    const auto *ksi0_292 = buffer.data(ksi0 + 292);
    const auto *ksi0_301 = buffer.data(ksi0 + 301);

    const auto *ksh_167 = buffer.data(ksh + 167);
    const auto *ksh_168 = buffer.data(ksh + 168);
    const auto *ksh_171 = buffer.data(ksh + 171);
    const auto *ksh_174 = buffer.data(ksh + 174);
    const auto *ksh_183 = buffer.data(ksh + 183);
    const auto *ksh_188 = buffer.data(ksh + 188);
    const auto *ksh_189 = buffer.data(ksh + 189);
    const auto *ksh_190 = buffer.data(ksh + 190);
    const auto *ksh_191 = buffer.data(ksh + 191);
    const auto *ksh_192 = buffer.data(ksh + 192);
    const auto *ksh_194 = buffer.data(ksh + 194);
    const auto *ksh_195 = buffer.data(ksh + 195);
    const auto *ksh_197 = buffer.data(ksh + 197);
    const auto *ksh_198 = buffer.data(ksh + 198);
    const auto *ksh_204 = buffer.data(ksh + 204);
    const auto *ksh_206 = buffer.data(ksh + 206);
    const auto *ksh_207 = buffer.data(ksh + 207);
    const auto *ksh_208 = buffer.data(ksh + 208);
    const auto *ksh_209 = buffer.data(ksh + 209);
    const auto *ksh_210 = buffer.data(ksh + 210);
    const auto *ksh_213 = buffer.data(ksh + 213);
    const auto *ksh_215 = buffer.data(ksh + 215);
    const auto *ksh_216 = buffer.data(ksh + 216);
    const auto *ksh_217 = buffer.data(ksh + 217);
    const auto *ksh_219 = buffer.data(ksh + 219);
    const auto *ksh_225 = buffer.data(ksh + 225);
    const auto *ksh_230 = buffer.data(ksh + 230);
    const auto *ksh_231 = buffer.data(ksh + 231);
    const auto *ksh_233 = buffer.data(ksh + 233);
    const auto *ksh_236 = buffer.data(ksh + 236);
    const auto *ksh_240 = buffer.data(ksh + 240);
    const auto *ksh_248 = buffer.data(ksh + 248);
    const auto *ksh_249 = buffer.data(ksh + 249);
    const auto *ksh_250 = buffer.data(ksh + 250);
    const auto *ksh_251 = buffer.data(ksh + 251);
    const auto *ksh_252 = buffer.data(ksh + 252);
    const auto *ksh_288 = buffer.data(ksh + 288);
    const auto *ksh_289 = buffer.data(ksh + 289);
    const auto *ksh_290 = buffer.data(ksh + 290);
    const auto *ksh_291 = buffer.data(ksh + 291);
    const auto *ksh_292 = buffer.data(ksh + 292);
    const auto *ksh_293 = buffer.data(ksh + 293);
    const auto *ksh_294 = buffer.data(ksh + 294);
    const auto *ksh_299 = buffer.data(ksh + 299);
    const auto *ksh_303 = buffer.data(ksh + 303);
    const auto *ksh_308 = buffer.data(ksh + 308);
    const auto *ksh_309 = buffer.data(ksh + 309);
    const auto *ksh_310 = buffer.data(ksh + 310);
    const auto *ksh_311 = buffer.data(ksh + 311);
    const auto *ksh_312 = buffer.data(ksh + 312);
    const auto *ksh_314 = buffer.data(ksh + 314);
    const auto *ksh_315 = buffer.data(ksh + 315);
    const auto *ksh_318 = buffer.data(ksh + 318);
    const auto *ksh_321 = buffer.data(ksh + 321);
    const auto *ksh_325 = buffer.data(ksh + 325);
    const auto *ksh_330 = buffer.data(ksh + 330);
    const auto *ksh_332 = buffer.data(ksh + 332);
    const auto *ksh_333 = buffer.data(ksh + 333);
    const auto *ksh_334 = buffer.data(ksh + 334);
    const auto *ksh_335 = buffer.data(ksh + 335);
    const auto *ksh_341 = buffer.data(ksh + 341);
    const auto *ksh_345 = buffer.data(ksh + 345);
    const auto *ksh_350 = buffer.data(ksh + 350);
    const auto *ksh_351 = buffer.data(ksh + 351);
    const auto *ksh_352 = buffer.data(ksh + 352);
    const auto *ksh_353 = buffer.data(ksh + 353);
    const auto *ksh_354 = buffer.data(ksh + 354);
    const auto *ksh_355 = buffer.data(ksh + 355);
    const auto *ksh_356 = buffer.data(ksh + 356);
    const auto *ksh_357 = buffer.data(ksh + 357);

    const auto *ksi1_252 = buffer.data(ksi1 + 252);
    const auto *ksi1_255 = buffer.data(ksi1 + 255);
    const auto *ksi1_257 = buffer.data(ksi1 + 257);
    const auto *ksi1_258 = buffer.data(ksi1 + 258);
    const auto *ksi1_261 = buffer.data(ksi1 + 261);
    const auto *ksi1_262 = buffer.data(ksi1 + 262);
    const auto *ksi1_264 = buffer.data(ksi1 + 264);
    const auto *ksi1_266 = buffer.data(ksi1 + 266);
    const auto *ksi1_279 = buffer.data(ksi1 + 279);
    const auto *ksi1_280 = buffer.data(ksi1 + 280);
    const auto *ksi1_283 = buffer.data(ksi1 + 283);
    const auto *ksi1_286 = buffer.data(ksi1 + 286);
    const auto *ksi1_290 = buffer.data(ksi1 + 290);
    const auto *ksi1_292 = buffer.data(ksi1 + 292);
    const auto *ksi1_301 = buffer.data(ksi1 + 301);

    const auto *lsg0_194 = buffer.data(lsg0 + 194);
    const auto *lsg0_205 = buffer.data(lsg0 + 205);
    const auto *lsg0_207 = buffer.data(lsg0 + 207);
    const auto *lsg0_208 = buffer.data(lsg0 + 208);
    const auto *lsg0_209 = buffer.data(lsg0 + 209);
    const auto *lsg0_210 = buffer.data(lsg0 + 210);
    const auto *lsg0_211 = buffer.data(lsg0 + 211);
    const auto *lsg0_212 = buffer.data(lsg0 + 212);
    const auto *lsg0_213 = buffer.data(lsg0 + 213);
    const auto *lsg0_214 = buffer.data(lsg0 + 214);
    const auto *lsg0_215 = buffer.data(lsg0 + 215);
    const auto *lsg0_219 = buffer.data(lsg0 + 219);
    const auto *lsg0_220 = buffer.data(lsg0 + 220);
    const auto *lsg0_221 = buffer.data(lsg0 + 221);
    const auto *lsg0_222 = buffer.data(lsg0 + 222);
    const auto *lsg0_223 = buffer.data(lsg0 + 223);
    const auto *lsg0_224 = buffer.data(lsg0 + 224);
    const auto *lsg0_225 = buffer.data(lsg0 + 225);
    const auto *lsg0_227 = buffer.data(lsg0 + 227);
    const auto *lsg0_228 = buffer.data(lsg0 + 228);
    const auto *lsg0_230 = buffer.data(lsg0 + 230);
    const auto *lsg0_231 = buffer.data(lsg0 + 231);
    const auto *lsg0_235 = buffer.data(lsg0 + 235);
    const auto *lsg0_236 = buffer.data(lsg0 + 236);
    const auto *lsg0_237 = buffer.data(lsg0 + 237);
    const auto *lsg0_239 = buffer.data(lsg0 + 239);
    const auto *lsg0_245 = buffer.data(lsg0 + 245);
    const auto *lsg0_249 = buffer.data(lsg0 + 249);
    const auto *lsg0_252 = buffer.data(lsg0 + 252);
    const auto *lsg0_253 = buffer.data(lsg0 + 253);
    const auto *lsg0_254 = buffer.data(lsg0 + 254);
    const auto *lsg0_255 = buffer.data(lsg0 + 255);

    const auto *lsg1_194 = buffer.data(lsg1 + 194);
    const auto *lsg1_205 = buffer.data(lsg1 + 205);
    const auto *lsg1_207 = buffer.data(lsg1 + 207);
    const auto *lsg1_208 = buffer.data(lsg1 + 208);
    const auto *lsg1_209 = buffer.data(lsg1 + 209);
    const auto *lsg1_210 = buffer.data(lsg1 + 210);
    const auto *lsg1_211 = buffer.data(lsg1 + 211);
    const auto *lsg1_212 = buffer.data(lsg1 + 212);
    const auto *lsg1_213 = buffer.data(lsg1 + 213);
    const auto *lsg1_214 = buffer.data(lsg1 + 214);
    const auto *lsg1_215 = buffer.data(lsg1 + 215);
    const auto *lsg1_219 = buffer.data(lsg1 + 219);
    const auto *lsg1_220 = buffer.data(lsg1 + 220);
    const auto *lsg1_221 = buffer.data(lsg1 + 221);
    const auto *lsg1_222 = buffer.data(lsg1 + 222);
    const auto *lsg1_223 = buffer.data(lsg1 + 223);
    const auto *lsg1_224 = buffer.data(lsg1 + 224);
    const auto *lsg1_225 = buffer.data(lsg1 + 225);
    const auto *lsg1_227 = buffer.data(lsg1 + 227);
    const auto *lsg1_228 = buffer.data(lsg1 + 228);
    const auto *lsg1_230 = buffer.data(lsg1 + 230);
    const auto *lsg1_231 = buffer.data(lsg1 + 231);
    const auto *lsg1_235 = buffer.data(lsg1 + 235);
    const auto *lsg1_236 = buffer.data(lsg1 + 236);
    const auto *lsg1_237 = buffer.data(lsg1 + 237);
    const auto *lsg1_239 = buffer.data(lsg1 + 239);
    const auto *lsg1_245 = buffer.data(lsg1 + 245);
    const auto *lsg1_249 = buffer.data(lsg1 + 249);
    const auto *lsg1_252 = buffer.data(lsg1 + 252);
    const auto *lsg1_253 = buffer.data(lsg1 + 253);
    const auto *lsg1_254 = buffer.data(lsg1 + 254);
    const auto *lsg1_255 = buffer.data(lsg1 + 255);

    const auto *lsh_272 = buffer.data(lsh + 272);
    const auto *lsh_273 = buffer.data(lsh + 273);
    const auto *lsh_275 = buffer.data(lsh + 275);
    const auto *lsh_276 = buffer.data(lsh + 276);
    const auto *lsh_278 = buffer.data(lsh + 278);
    const auto *lsh_279 = buffer.data(lsh + 279);
    const auto *lsh_282 = buffer.data(lsh + 282);
    const auto *lsh_288 = buffer.data(lsh + 288);
    const auto *lsh_289 = buffer.data(lsh + 289);
    const auto *lsh_290 = buffer.data(lsh + 290);
    const auto *lsh_291 = buffer.data(lsh + 291);
    const auto *lsh_292 = buffer.data(lsh + 292);
    const auto *lsh_293 = buffer.data(lsh + 293);
    const auto *lsh_294 = buffer.data(lsh + 294);
    const auto *lsh_295 = buffer.data(lsh + 295);
    const auto *lsh_296 = buffer.data(lsh + 296);
    const auto *lsh_297 = buffer.data(lsh + 297);
    const auto *lsh_298 = buffer.data(lsh + 298);
    const auto *lsh_299 = buffer.data(lsh + 299);
    const auto *lsh_300 = buffer.data(lsh + 300);
    const auto *lsh_301 = buffer.data(lsh + 301);
    const auto *lsh_302 = buffer.data(lsh + 302);
    const auto *lsh_303 = buffer.data(lsh + 303);
    const auto *lsh_308 = buffer.data(lsh + 308);
    const auto *lsh_309 = buffer.data(lsh + 309);
    const auto *lsh_310 = buffer.data(lsh + 310);
    const auto *lsh_311 = buffer.data(lsh + 311);
    const auto *lsh_312 = buffer.data(lsh + 312);
    const auto *lsh_313 = buffer.data(lsh + 313);
    const auto *lsh_314 = buffer.data(lsh + 314);
    const auto *lsh_315 = buffer.data(lsh + 315);
    const auto *lsh_316 = buffer.data(lsh + 316);
    const auto *lsh_317 = buffer.data(lsh + 317);
    const auto *lsh_318 = buffer.data(lsh + 318);
    const auto *lsh_320 = buffer.data(lsh + 320);
    const auto *lsh_321 = buffer.data(lsh + 321);
    const auto *lsh_322 = buffer.data(lsh + 322);
    const auto *lsh_324 = buffer.data(lsh + 324);
    const auto *lsh_325 = buffer.data(lsh + 325);
    const auto *lsh_330 = buffer.data(lsh + 330);
    const auto *lsh_331 = buffer.data(lsh + 331);
    const auto *lsh_332 = buffer.data(lsh + 332);
    const auto *lsh_333 = buffer.data(lsh + 333);
    const auto *lsh_334 = buffer.data(lsh + 334);
    const auto *lsh_335 = buffer.data(lsh + 335);
    const auto *lsh_336 = buffer.data(lsh + 336);
    const auto *lsh_338 = buffer.data(lsh + 338);
    const auto *lsh_339 = buffer.data(lsh + 339);
    const auto *lsh_341 = buffer.data(lsh + 341);
    const auto *lsh_342 = buffer.data(lsh + 342);
    const auto *lsh_345 = buffer.data(lsh + 345);
    const auto *lsh_350 = buffer.data(lsh + 350);
    const auto *lsh_351 = buffer.data(lsh + 351);
    const auto *lsh_352 = buffer.data(lsh + 352);
    const auto *lsh_353 = buffer.data(lsh + 353);
    const auto *lsh_354 = buffer.data(lsh + 354);
    const auto *lsh_355 = buffer.data(lsh + 355);
    const auto *lsh_356 = buffer.data(lsh + 356);
    const auto *lsh_357 = buffer.data(lsh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, ksi0_252, ksh_167, \
                         ksh_188, ksh_189, ksi1_252, lsg0_194, lsg1_194, lsh_272, \
                         lsh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * ksh_188[k]
                   + f_3 * pc_y[k] * lsh_272[k];

        t_363[k] = f_12 * ksh_167[k]
                   + f_1 * lsg0_194[k]
                   - f_2 * lsg1_194[k]
                   + f_3 * pc_z[k] * lsh_272[k];

        t_364[k] = pa_y[k] * ksi0_252[k]
                   - f_10 * pc_y[k] * ksi1_252[k];

        t_365[k] = f_11 * ksh_189[k]
                   + f_3 * pc_y[k] * lsh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, ksi0_255, ksi0_257, \
                         ksh_168, ksh_190, ksh_191, ksi1_255, ksi1_257, lsh_273, \
                         lsh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * ksh_168[k]
                   + f_3 * pc_z[k] * lsh_273[k];

        t_367[k] = pa_y[k] * ksi0_255[k]
                   + f_12 * ksh_190[k]
                   - f_10 * pc_y[k] * ksi1_255[k];

        t_368[k] = f_11 * ksh_191[k]
                   + f_3 * pc_y[k] * lsh_275[k];

        t_369[k] = pa_y[k] * ksi0_257[k]
                   - f_10 * pc_y[k] * ksi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, ksi0_258, ksi0_261, \
                         ksh_171, ksh_192, ksh_194, ksi1_258, ksi1_261, lsh_276, \
                         lsh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * ksi0_258[k]
                   + f_13 * ksh_192[k]
                   - f_10 * pc_y[k] * ksi1_258[k];

        t_371[k] = f_13 * ksh_171[k]
                   + f_3 * pc_z[k] * lsh_276[k];

        t_372[k] = f_11 * ksh_194[k]
                   + f_3 * pc_y[k] * lsh_278[k];

        t_373[k] = pa_y[k] * ksi0_261[k]
                   - f_10 * pc_y[k] * ksi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, ksi0_262, ksi0_264, ksh_174, \
                         ksh_195, ksh_197, ksi1_262, ksi1_264, \
                         lsh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * ksi0_262[k]
                   + f_14 * ksh_195[k]
                   - f_10 * pc_y[k] * ksi1_262[k];

        t_375[k] = f_13 * ksh_174[k]
                   + f_3 * pc_z[k] * lsh_279[k];

        t_376[k] = pa_y[k] * ksi0_264[k]
                   + f_12 * ksh_197[k]
                   - f_10 * pc_y[k] * ksi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, ksi0_266, ksh_198, \
                         ksh_288, ksh_289, ksi1_266, lsh_282, lsh_288, \
                         lsh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * ksh_198[k]
                   + f_3 * pc_y[k] * lsh_282[k];

        t_378[k] = pa_y[k] * ksi0_266[k]
                   - f_10 * pc_y[k] * ksi1_266[k];

        t_379[k] = f_14 * ksh_288[k]
                   + f_3 * pc_x[k] * lsh_288[k];

        t_380[k] = f_14 * ksh_289[k]
                   + f_3 * pc_x[k] * lsh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, ksh_290, ksh_291, ksh_292, ksh_293, \
                         lsh_290, lsh_291, lsh_292, lsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_14 * ksh_290[k]
                   + f_3 * pc_x[k] * lsh_290[k];

        t_382[k] = f_14 * ksh_291[k]
                   + f_3 * pc_x[k] * lsh_291[k];

        t_383[k] = f_14 * ksh_292[k]
                   + f_3 * pc_x[k] * lsh_292[k];

        t_384[k] = f_14 * ksh_293[k]
                   + f_3 * pc_x[k] * lsh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, ksh_183, ksh_204, ksh_206, lsg0_205, \
                         lsg0_207, lsg1_205, lsg1_207, lsh_288, \
                         lsh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * ksh_204[k]
                   + f_1 * lsg0_205[k]
                   - f_2 * lsg1_205[k]
                   + f_3 * pc_y[k] * lsh_288[k];

        t_386[k] = f_13 * ksh_183[k]
                   + f_3 * pc_z[k] * lsh_288[k];

        t_387[k] = f_11 * ksh_206[k]
                   + f_8 * lsg0_207[k]
                   - f_9 * lsg1_207[k]
                   + f_3 * pc_y[k] * lsh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, ksh_207, ksh_208, ksh_209, lsg0_208, \
                         lsg0_209, lsg1_208, lsg1_209, lsh_291, lsh_292, \
                         lsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * ksh_207[k]
                   + f_6 * lsg0_208[k]
                   - f_7 * lsg1_208[k]
                   + f_3 * pc_y[k] * lsh_291[k];

        t_389[k] = f_11 * ksh_208[k]
                   + f_4 * lsg0_209[k]
                   - f_5 * lsg1_209[k]
                   + f_3 * pc_y[k] * lsh_292[k];

        t_390[k] = f_11 * ksh_209[k]
                   + f_3 * pc_y[k] * lsh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, ksi0_279, \
                         ksh_189, ksh_294, ksi1_279, lsg0_210, lsg1_210, \
                         lsh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * ksi0_279[k]
                   - f_10 * pc_y[k] * ksi1_279[k];

        t_392[k] = f_14 * ksh_294[k]
                   + f_1 * lsg0_210[k]
                   - f_2 * lsg1_210[k]
                   + f_3 * pc_x[k] * lsh_294[k];

        t_393[k] = f_3 * pc_y[k] * lsh_294[k];

        t_394[k] = f_14 * ksh_189[k]
                   + f_3 * pc_z[k] * lsh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, ksh_299, lsg0_210, lsg0_215, \
                         lsg1_210, lsg1_215, lsh_295, lsh_296, \
                         lsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * lsg0_210[k]
                   - f_5 * lsg1_210[k]
                   + f_3 * pc_y[k] * lsh_295[k];

        t_396[k] = f_3 * pc_y[k] * lsh_296[k];

        t_397[k] = f_14 * ksh_299[k]
                   + f_8 * lsg0_215[k]
                   - f_9 * lsg1_215[k]
                   + f_3 * pc_x[k] * lsh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, lsg0_211, lsg0_212, lsg1_211, lsg1_212, \
                         lsh_297, lsh_298, lsh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * lsg0_211[k]
                   - f_7 * lsg1_211[k]
                   + f_3 * pc_y[k] * lsh_297[k];

        t_399[k] = f_4 * lsg0_212[k]
                   - f_5 * lsg1_212[k]
                   + f_3 * pc_y[k] * lsh_298[k];

        t_400[k] = f_3 * pc_y[k] * lsh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, ksh_303, lsg0_213, lsg0_214, \
                         lsg0_219, lsg1_213, lsg1_214, lsg1_219, lsh_300, lsh_301, \
                         lsh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_14 * ksh_303[k]
                   + f_6 * lsg0_219[k]
                   - f_7 * lsg1_219[k]
                   + f_3 * pc_x[k] * lsh_303[k];

        t_402[k] = f_8 * lsg0_213[k]
                   - f_9 * lsg1_213[k]
                   + f_3 * pc_y[k] * lsh_300[k];

        t_403[k] = f_6 * lsg0_214[k]
                   - f_7 * lsg1_214[k]
                   + f_3 * pc_y[k] * lsh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, ksh_308, ksh_309, lsg0_215, \
                         lsg0_224, lsg1_215, lsg1_224, lsh_302, lsh_303, lsh_308, \
                         lsh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * lsg0_215[k]
                   - f_5 * lsg1_215[k]
                   + f_3 * pc_y[k] * lsh_302[k];

        t_405[k] = f_3 * pc_y[k] * lsh_303[k];

        t_406[k] = f_14 * ksh_308[k]
                   + f_4 * lsg0_224[k]
                   - f_5 * lsg1_224[k]
                   + f_3 * pc_x[k] * lsh_308[k];

        t_407[k] = f_14 * ksh_309[k]
                   + f_3 * pc_x[k] * lsh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, ksh_310, ksh_311, \
                         ksh_312, ksh_314, lsh_308, lsh_310, lsh_311, lsh_312, \
                         lsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_14 * ksh_310[k]
                   + f_3 * pc_x[k] * lsh_310[k];

        t_409[k] = f_14 * ksh_311[k]
                   + f_3 * pc_x[k] * lsh_311[k];

        t_410[k] = f_14 * ksh_312[k]
                   + f_3 * pc_x[k] * lsh_312[k];

        t_411[k] = f_3 * pc_y[k] * lsh_308[k];

        t_412[k] = f_14 * ksh_314[k]
                   + f_3 * pc_x[k] * lsh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, lsg0_220, lsg0_221, lsg0_222, lsg1_220, \
                         lsg1_221, lsg1_222, lsh_309, lsh_310, \
                         lsh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * lsg0_220[k]
                   - f_2 * lsg1_220[k]
                   + f_3 * pc_y[k] * lsh_309[k];

        t_414[k] = f_16 * lsg0_221[k]
                   - f_17 * lsg1_221[k]
                   + f_3 * pc_y[k] * lsh_310[k];

        t_415[k] = f_8 * lsg0_222[k]
                   - f_9 * lsg1_222[k]
                   + f_3 * pc_y[k] * lsh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, ksh_209, lsg0_223, lsg0_224, \
                         lsg1_223, lsg1_224, lsh_312, lsh_313, \
                         lsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * lsg0_223[k]
                   - f_7 * lsg1_223[k]
                   + f_3 * pc_y[k] * lsh_312[k];

        t_417[k] = f_4 * lsg0_224[k]
                   - f_5 * lsg1_224[k]
                   + f_3 * pc_y[k] * lsh_313[k];

        t_418[k] = f_3 * pc_y[k] * lsh_314[k];

        t_419[k] = f_14 * ksh_209[k]
                   + f_1 * lsg0_224[k]
                   - f_2 * lsg1_224[k]
                   + f_3 * pc_z[k] * lsh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, ksh_210, ksh_315, \
                         ksh_318, lsg0_225, lsg0_228, lsg1_225, lsg1_228, lsh_315, \
                         lsh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_13 * ksh_315[k]
                   + f_1 * lsg0_225[k]
                   - f_2 * lsg1_225[k]
                   + f_3 * pc_x[k] * lsh_315[k];

        t_421[k] = f_19 * ksh_210[k]
                   + f_3 * pc_y[k] * lsh_315[k];

        t_422[k] = f_3 * pc_z[k] * lsh_315[k];

        t_423[k] = f_13 * ksh_318[k]
                   + f_8 * lsg0_228[k]
                   - f_9 * lsg1_228[k]
                   + f_3 * pc_x[k] * lsh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, ksh_321, lsg0_225, lsg0_231, \
                         lsg1_225, lsg1_231, lsh_316, lsh_317, lsh_318, \
                         lsh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * lsh_316[k];

        t_425[k] = f_4 * lsg0_225[k]
                   - f_5 * lsg1_225[k]
                   + f_3 * pc_z[k] * lsh_317[k];

        t_426[k] = f_13 * ksh_321[k]
                   + f_6 * lsg0_231[k]
                   - f_7 * lsg1_231[k]
                   + f_3 * pc_x[k] * lsh_321[k];

        t_427[k] = f_3 * pc_z[k] * lsh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, ksh_215, ksh_325, \
                         lsg0_227, lsg0_235, lsg1_227, lsg1_235, lsh_320, lsh_321, \
                         lsh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_19 * ksh_215[k]
                   + f_3 * pc_y[k] * lsh_320[k];

        t_429[k] = f_6 * lsg0_227[k]
                   - f_7 * lsg1_227[k]
                   + f_3 * pc_z[k] * lsh_320[k];

        t_430[k] = f_13 * ksh_325[k]
                   + f_4 * lsg0_235[k]
                   - f_5 * lsg1_235[k]
                   + f_3 * pc_x[k] * lsh_325[k];

        t_431[k] = f_3 * pc_z[k] * lsh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, ksh_219, ksh_330, \
                         lsg0_228, lsg0_230, lsg1_228, lsg1_230, lsh_322, lsh_324, \
                         lsh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * lsg0_228[k]
                   - f_5 * lsg1_228[k]
                   + f_3 * pc_z[k] * lsh_322[k];

        t_433[k] = f_19 * ksh_219[k]
                   + f_3 * pc_y[k] * lsh_324[k];

        t_434[k] = f_8 * lsg0_230[k]
                   - f_9 * lsg1_230[k]
                   + f_3 * pc_z[k] * lsh_324[k];

        t_435[k] = f_13 * ksh_330[k]
                   + f_3 * pc_x[k] * lsh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, ksh_332, ksh_333, \
                         ksh_334, ksh_335, lsh_325, lsh_332, lsh_333, lsh_334, \
                         lsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * lsh_325[k];

        t_437[k] = f_13 * ksh_332[k]
                   + f_3 * pc_x[k] * lsh_332[k];

        t_438[k] = f_13 * ksh_333[k]
                   + f_3 * pc_x[k] * lsh_333[k];

        t_439[k] = f_13 * ksh_334[k]
                   + f_3 * pc_x[k] * lsh_334[k];

        t_440[k] = f_13 * ksh_335[k]
                   + f_3 * pc_x[k] * lsh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, ksh_225, lsg0_235, lsg0_236, \
                         lsg1_235, lsg1_236, lsh_330, lsh_331, \
                         lsh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_19 * ksh_225[k]
                   + f_1 * lsg0_235[k]
                   - f_2 * lsg1_235[k]
                   + f_3 * pc_y[k] * lsh_330[k];

        t_442[k] = f_3 * pc_z[k] * lsh_330[k];

        t_443[k] = f_4 * lsg0_235[k]
                   - f_5 * lsg1_235[k]
                   + f_3 * pc_z[k] * lsh_331[k];

        t_444[k] = f_6 * lsg0_236[k]
                   - f_7 * lsg1_236[k]
                   + f_3 * pc_z[k] * lsh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, ksi0_280, ksh_230, \
                         ksi1_280, lsg0_237, lsg0_239, lsg1_237, lsg1_239, lsh_333, \
                         lsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * lsg0_237[k]
                   - f_9 * lsg1_237[k]
                   + f_3 * pc_z[k] * lsh_333[k];

        t_446[k] = f_19 * ksh_230[k]
                   + f_3 * pc_y[k] * lsh_335[k];

        t_447[k] = f_1 * lsg0_239[k]
                   - f_2 * lsg1_239[k]
                   + f_3 * pc_z[k] * lsh_335[k];

        t_448[k] = pa_z[k] * ksi0_280[k]
                   - f_10 * pc_z[k] * ksi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, ksi0_283, ksh_210, \
                         ksh_231, ksh_233, ksi1_283, lsh_336, lsh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * ksh_231[k]
                   + f_3 * pc_y[k] * lsh_336[k];

        t_450[k] = f_11 * ksh_210[k]
                   + f_3 * pc_z[k] * lsh_336[k];

        t_451[k] = pa_z[k] * ksi0_283[k]
                   - f_10 * pc_z[k] * ksi1_283[k];

        t_452[k] = f_14 * ksh_233[k]
                   + f_3 * pc_y[k] * lsh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, ksi0_286, ksh_213, ksh_341, \
                         ksi1_286, lsg0_245, lsg1_245, lsh_339, \
                         lsh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_13 * ksh_341[k]
                   + f_8 * lsg0_245[k]
                   - f_9 * lsg1_245[k]
                   + f_3 * pc_x[k] * lsh_341[k];

        t_454[k] = pa_z[k] * ksi0_286[k]
                   - f_10 * pc_z[k] * ksi1_286[k];

        t_455[k] = f_11 * ksh_213[k]
                   + f_3 * pc_z[k] * lsh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, ksi0_290, ksh_236, \
                         ksh_345, ksi1_290, lsg0_249, lsg1_249, lsh_341, \
                         lsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * ksh_236[k]
                   + f_3 * pc_y[k] * lsh_341[k];

        t_457[k] = f_13 * ksh_345[k]
                   + f_6 * lsg0_249[k]
                   - f_7 * lsg1_249[k]
                   + f_3 * pc_x[k] * lsh_345[k];

        t_458[k] = pa_z[k] * ksi0_290[k]
                   - f_10 * pc_z[k] * ksi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, ksi0_292, ksh_216, ksh_217, \
                         ksh_240, ksi1_292, lsh_342, lsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * ksh_216[k]
                   + f_3 * pc_z[k] * lsh_342[k];

        t_460[k] = pa_z[k] * ksi0_292[k]
                   + f_12 * ksh_217[k]
                   - f_10 * pc_z[k] * ksi1_292[k];

        t_461[k] = f_14 * ksh_240[k]
                   + f_3 * pc_y[k] * lsh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, ksh_350, ksh_351, ksh_352, ksh_353, \
                         lsg0_254, lsg1_254, lsh_350, lsh_351, lsh_352, \
                         lsh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_13 * ksh_350[k]
                   + f_4 * lsg0_254[k]
                   - f_5 * lsg1_254[k]
                   + f_3 * pc_x[k] * lsh_350[k];

        t_463[k] = f_13 * ksh_351[k]
                   + f_3 * pc_x[k] * lsh_351[k];

        t_464[k] = f_13 * ksh_352[k]
                   + f_3 * pc_x[k] * lsh_352[k];

        t_465[k] = f_13 * ksh_353[k]
                   + f_3 * pc_x[k] * lsh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, ksi0_301, ksh_354, \
                         ksh_355, ksh_356, ksi1_301, lsh_354, lsh_355, \
                         lsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_13 * ksh_354[k]
                   + f_3 * pc_x[k] * lsh_354[k];

        t_467[k] = f_13 * ksh_355[k]
                   + f_3 * pc_x[k] * lsh_355[k];

        t_468[k] = f_13 * ksh_356[k]
                   + f_3 * pc_x[k] * lsh_356[k];

        t_469[k] = pa_z[k] * ksi0_301[k]
                   - f_10 * pc_z[k] * ksi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, ksh_225, ksh_248, ksh_249, lsg0_252, \
                         lsg0_253, lsg1_252, lsg1_253, lsh_351, lsh_353, \
                         lsh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * ksh_225[k]
                   + f_3 * pc_z[k] * lsh_351[k];

        t_471[k] = f_14 * ksh_248[k]
                   + f_8 * lsg0_252[k]
                   - f_9 * lsg1_252[k]
                   + f_3 * pc_y[k] * lsh_353[k];

        t_472[k] = f_14 * ksh_249[k]
                   + f_6 * lsg0_253[k]
                   - f_7 * lsg1_253[k]
                   + f_3 * pc_y[k] * lsh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, ksh_230, ksh_250, ksh_251, lsg0_254, \
                         lsg1_254, lsh_355, lsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * ksh_250[k]
                   + f_4 * lsg0_254[k]
                   - f_5 * lsg1_254[k]
                   + f_3 * pc_y[k] * lsh_355[k];

        t_474[k] = f_14 * ksh_251[k]
                   + f_3 * pc_y[k] * lsh_356[k];

        t_475[k] = f_11 * ksh_230[k]
                   + f_1 * lsg0_254[k]
                   - f_2 * lsg1_254[k]
                   + f_3 * pc_z[k] * lsh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, ksh_231, ksh_252, ksh_357, \
                         lsg0_255, lsg1_255, lsh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_13 * ksh_357[k]
                   + f_1 * lsg0_255[k]
                   - f_2 * lsg1_255[k]
                   + f_3 * pc_x[k] * lsh_357[k];

        t_477[k] = f_13 * ksh_252[k]
                   + f_3 * pc_y[k] * lsh_357[k];

        t_478[k] = f_12 * ksh_231[k]
                   + f_3 * pc_z[k] * lsh_357[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_392 = buffer.data(ksi0 + 392);
    const auto *ksi0_395 = buffer.data(ksi0 + 395);
    const auto *ksi0_397 = buffer.data(ksi0 + 397);
    const auto *ksi0_398 = buffer.data(ksi0 + 398);
    const auto *ksi0_401 = buffer.data(ksi0 + 401);
    const auto *ksi0_402 = buffer.data(ksi0 + 402);
    const auto *ksi0_404 = buffer.data(ksi0 + 404);
    const auto *ksi0_406 = buffer.data(ksi0 + 406);
    const auto *ksi0_419 = buffer.data(ksi0 + 419);

    const auto *ksh_234 = buffer.data(ksh + 234);
    const auto *ksh_237 = buffer.data(ksh + 237);
    const auto *ksh_246 = buffer.data(ksh + 246);
    const auto *ksh_251 = buffer.data(ksh + 251);
    const auto *ksh_252 = buffer.data(ksh + 252);
    const auto *ksh_254 = buffer.data(ksh + 254);
    const auto *ksh_255 = buffer.data(ksh + 255);
    const auto *ksh_257 = buffer.data(ksh + 257);
    const auto *ksh_258 = buffer.data(ksh + 258);
    const auto *ksh_261 = buffer.data(ksh + 261);
    const auto *ksh_267 = buffer.data(ksh + 267);
    const auto *ksh_269 = buffer.data(ksh + 269);
    const auto *ksh_270 = buffer.data(ksh + 270);
    const auto *ksh_271 = buffer.data(ksh + 271);
    const auto *ksh_272 = buffer.data(ksh + 272);
    const auto *ksh_273 = buffer.data(ksh + 273);
    const auto *ksh_275 = buffer.data(ksh + 275);
    const auto *ksh_276 = buffer.data(ksh + 276);
    const auto *ksh_278 = buffer.data(ksh + 278);
    const auto *ksh_279 = buffer.data(ksh + 279);
    const auto *ksh_282 = buffer.data(ksh + 282);
    const auto *ksh_288 = buffer.data(ksh + 288);
    const auto *ksh_290 = buffer.data(ksh + 290);
    const auto *ksh_291 = buffer.data(ksh + 291);
    const auto *ksh_292 = buffer.data(ksh + 292);
    const auto *ksh_293 = buffer.data(ksh + 293);
    const auto *ksh_294 = buffer.data(ksh + 294);
    const auto *ksh_295 = buffer.data(ksh + 295);
    const auto *ksh_296 = buffer.data(ksh + 296);
    const auto *ksh_297 = buffer.data(ksh + 297);
    const auto *ksh_299 = buffer.data(ksh + 299);
    const auto *ksh_300 = buffer.data(ksh + 300);
    const auto *ksh_302 = buffer.data(ksh + 302);
    const auto *ksh_303 = buffer.data(ksh + 303);
    const auto *ksh_309 = buffer.data(ksh + 309);
    const auto *ksh_311 = buffer.data(ksh + 311);
    const auto *ksh_312 = buffer.data(ksh + 312);
    const auto *ksh_313 = buffer.data(ksh + 313);
    const auto *ksh_314 = buffer.data(ksh + 314);
    const auto *ksh_315 = buffer.data(ksh + 315);
    const auto *ksh_360 = buffer.data(ksh + 360);
    const auto *ksh_362 = buffer.data(ksh + 362);
    const auto *ksh_363 = buffer.data(ksh + 363);
    const auto *ksh_366 = buffer.data(ksh + 366);
    const auto *ksh_367 = buffer.data(ksh + 367);
    const auto *ksh_369 = buffer.data(ksh + 369);
    const auto *ksh_371 = buffer.data(ksh + 371);
    const auto *ksh_372 = buffer.data(ksh + 372);
    const auto *ksh_373 = buffer.data(ksh + 373);
    const auto *ksh_374 = buffer.data(ksh + 374);
    const auto *ksh_375 = buffer.data(ksh + 375);
    const auto *ksh_376 = buffer.data(ksh + 376);
    const auto *ksh_377 = buffer.data(ksh + 377);
    const auto *ksh_378 = buffer.data(ksh + 378);
    const auto *ksh_381 = buffer.data(ksh + 381);
    const auto *ksh_383 = buffer.data(ksh + 383);
    const auto *ksh_384 = buffer.data(ksh + 384);
    const auto *ksh_387 = buffer.data(ksh + 387);
    const auto *ksh_388 = buffer.data(ksh + 388);
    const auto *ksh_390 = buffer.data(ksh + 390);
    const auto *ksh_392 = buffer.data(ksh + 392);
    const auto *ksh_393 = buffer.data(ksh + 393);
    const auto *ksh_394 = buffer.data(ksh + 394);
    const auto *ksh_395 = buffer.data(ksh + 395);
    const auto *ksh_396 = buffer.data(ksh + 396);
    const auto *ksh_397 = buffer.data(ksh + 397);
    const auto *ksh_398 = buffer.data(ksh + 398);
    const auto *ksh_414 = buffer.data(ksh + 414);
    const auto *ksh_415 = buffer.data(ksh + 415);
    const auto *ksh_416 = buffer.data(ksh + 416);
    const auto *ksh_417 = buffer.data(ksh + 417);
    const auto *ksh_418 = buffer.data(ksh + 418);
    const auto *ksh_419 = buffer.data(ksh + 419);
    const auto *ksh_420 = buffer.data(ksh + 420);
    const auto *ksh_425 = buffer.data(ksh + 425);
    const auto *ksh_429 = buffer.data(ksh + 429);
    const auto *ksh_434 = buffer.data(ksh + 434);
    const auto *ksh_435 = buffer.data(ksh + 435);
    const auto *ksh_436 = buffer.data(ksh + 436);
    const auto *ksh_437 = buffer.data(ksh + 437);
    const auto *ksh_438 = buffer.data(ksh + 438);
    const auto *ksh_440 = buffer.data(ksh + 440);
    const auto *ksh_441 = buffer.data(ksh + 441);
    const auto *ksh_444 = buffer.data(ksh + 444);

    const auto *ksi1_392 = buffer.data(ksi1 + 392);
    const auto *ksi1_395 = buffer.data(ksi1 + 395);
    const auto *ksi1_397 = buffer.data(ksi1 + 397);
    const auto *ksi1_398 = buffer.data(ksi1 + 398);
    const auto *ksi1_401 = buffer.data(ksi1 + 401);
    const auto *ksi1_402 = buffer.data(ksi1 + 402);
    const auto *ksi1_404 = buffer.data(ksi1 + 404);
    const auto *ksi1_406 = buffer.data(ksi1 + 406);
    const auto *ksi1_419 = buffer.data(ksi1 + 419);

    const auto *lsg0_258 = buffer.data(lsg0 + 258);
    const auto *lsg0_260 = buffer.data(lsg0 + 260);
    const auto *lsg0_261 = buffer.data(lsg0 + 261);
    const auto *lsg0_264 = buffer.data(lsg0 + 264);
    const auto *lsg0_265 = buffer.data(lsg0 + 265);
    const auto *lsg0_267 = buffer.data(lsg0 + 267);
    const auto *lsg0_268 = buffer.data(lsg0 + 268);
    const auto *lsg0_269 = buffer.data(lsg0 + 269);
    const auto *lsg0_270 = buffer.data(lsg0 + 270);
    const auto *lsg0_273 = buffer.data(lsg0 + 273);
    const auto *lsg0_275 = buffer.data(lsg0 + 275);
    const auto *lsg0_276 = buffer.data(lsg0 + 276);
    const auto *lsg0_279 = buffer.data(lsg0 + 279);
    const auto *lsg0_280 = buffer.data(lsg0 + 280);
    const auto *lsg0_282 = buffer.data(lsg0 + 282);
    const auto *lsg0_283 = buffer.data(lsg0 + 283);
    const auto *lsg0_284 = buffer.data(lsg0 + 284);
    const auto *lsg0_295 = buffer.data(lsg0 + 295);
    const auto *lsg0_297 = buffer.data(lsg0 + 297);
    const auto *lsg0_298 = buffer.data(lsg0 + 298);
    const auto *lsg0_299 = buffer.data(lsg0 + 299);
    const auto *lsg0_300 = buffer.data(lsg0 + 300);
    const auto *lsg0_301 = buffer.data(lsg0 + 301);
    const auto *lsg0_302 = buffer.data(lsg0 + 302);
    const auto *lsg0_303 = buffer.data(lsg0 + 303);
    const auto *lsg0_304 = buffer.data(lsg0 + 304);
    const auto *lsg0_305 = buffer.data(lsg0 + 305);
    const auto *lsg0_309 = buffer.data(lsg0 + 309);
    const auto *lsg0_310 = buffer.data(lsg0 + 310);
    const auto *lsg0_311 = buffer.data(lsg0 + 311);
    const auto *lsg0_312 = buffer.data(lsg0 + 312);
    const auto *lsg0_313 = buffer.data(lsg0 + 313);
    const auto *lsg0_314 = buffer.data(lsg0 + 314);
    const auto *lsg0_315 = buffer.data(lsg0 + 315);
    const auto *lsg0_318 = buffer.data(lsg0 + 318);

    const auto *lsg1_258 = buffer.data(lsg1 + 258);
    const auto *lsg1_260 = buffer.data(lsg1 + 260);
    const auto *lsg1_261 = buffer.data(lsg1 + 261);
    const auto *lsg1_264 = buffer.data(lsg1 + 264);
    const auto *lsg1_265 = buffer.data(lsg1 + 265);
    const auto *lsg1_267 = buffer.data(lsg1 + 267);
    const auto *lsg1_268 = buffer.data(lsg1 + 268);
    const auto *lsg1_269 = buffer.data(lsg1 + 269);
    const auto *lsg1_270 = buffer.data(lsg1 + 270);
    const auto *lsg1_273 = buffer.data(lsg1 + 273);
    const auto *lsg1_275 = buffer.data(lsg1 + 275);
    const auto *lsg1_276 = buffer.data(lsg1 + 276);
    const auto *lsg1_279 = buffer.data(lsg1 + 279);
    const auto *lsg1_280 = buffer.data(lsg1 + 280);
    const auto *lsg1_282 = buffer.data(lsg1 + 282);
    const auto *lsg1_283 = buffer.data(lsg1 + 283);
    const auto *lsg1_284 = buffer.data(lsg1 + 284);
    const auto *lsg1_295 = buffer.data(lsg1 + 295);
    const auto *lsg1_297 = buffer.data(lsg1 + 297);
    const auto *lsg1_298 = buffer.data(lsg1 + 298);
    const auto *lsg1_299 = buffer.data(lsg1 + 299);
    const auto *lsg1_300 = buffer.data(lsg1 + 300);
    const auto *lsg1_301 = buffer.data(lsg1 + 301);
    const auto *lsg1_302 = buffer.data(lsg1 + 302);
    const auto *lsg1_303 = buffer.data(lsg1 + 303);
    const auto *lsg1_304 = buffer.data(lsg1 + 304);
    const auto *lsg1_305 = buffer.data(lsg1 + 305);
    const auto *lsg1_309 = buffer.data(lsg1 + 309);
    const auto *lsg1_310 = buffer.data(lsg1 + 310);
    const auto *lsg1_311 = buffer.data(lsg1 + 311);
    const auto *lsg1_312 = buffer.data(lsg1 + 312);
    const auto *lsg1_313 = buffer.data(lsg1 + 313);
    const auto *lsg1_314 = buffer.data(lsg1 + 314);
    const auto *lsg1_315 = buffer.data(lsg1 + 315);
    const auto *lsg1_318 = buffer.data(lsg1 + 318);

    const auto *lsh_359 = buffer.data(lsh + 359);
    const auto *lsh_360 = buffer.data(lsh + 360);
    const auto *lsh_362 = buffer.data(lsh + 362);
    const auto *lsh_363 = buffer.data(lsh + 363);
    const auto *lsh_366 = buffer.data(lsh + 366);
    const auto *lsh_367 = buffer.data(lsh + 367);
    const auto *lsh_369 = buffer.data(lsh + 369);
    const auto *lsh_371 = buffer.data(lsh + 371);
    const auto *lsh_372 = buffer.data(lsh + 372);
    const auto *lsh_373 = buffer.data(lsh + 373);
    const auto *lsh_374 = buffer.data(lsh + 374);
    const auto *lsh_375 = buffer.data(lsh + 375);
    const auto *lsh_376 = buffer.data(lsh + 376);
    const auto *lsh_377 = buffer.data(lsh + 377);
    const auto *lsh_378 = buffer.data(lsh + 378);
    const auto *lsh_380 = buffer.data(lsh + 380);
    const auto *lsh_381 = buffer.data(lsh + 381);
    const auto *lsh_383 = buffer.data(lsh + 383);
    const auto *lsh_384 = buffer.data(lsh + 384);
    const auto *lsh_387 = buffer.data(lsh + 387);
    const auto *lsh_388 = buffer.data(lsh + 388);
    const auto *lsh_390 = buffer.data(lsh + 390);
    const auto *lsh_392 = buffer.data(lsh + 392);
    const auto *lsh_393 = buffer.data(lsh + 393);
    const auto *lsh_394 = buffer.data(lsh + 394);
    const auto *lsh_395 = buffer.data(lsh + 395);
    const auto *lsh_396 = buffer.data(lsh + 396);
    const auto *lsh_397 = buffer.data(lsh + 397);
    const auto *lsh_398 = buffer.data(lsh + 398);
    const auto *lsh_399 = buffer.data(lsh + 399);
    const auto *lsh_401 = buffer.data(lsh + 401);
    const auto *lsh_402 = buffer.data(lsh + 402);
    const auto *lsh_404 = buffer.data(lsh + 404);
    const auto *lsh_405 = buffer.data(lsh + 405);
    const auto *lsh_408 = buffer.data(lsh + 408);
    const auto *lsh_414 = buffer.data(lsh + 414);
    const auto *lsh_415 = buffer.data(lsh + 415);
    const auto *lsh_416 = buffer.data(lsh + 416);
    const auto *lsh_417 = buffer.data(lsh + 417);
    const auto *lsh_418 = buffer.data(lsh + 418);
    const auto *lsh_419 = buffer.data(lsh + 419);
    const auto *lsh_420 = buffer.data(lsh + 420);
    const auto *lsh_421 = buffer.data(lsh + 421);
    const auto *lsh_422 = buffer.data(lsh + 422);
    const auto *lsh_423 = buffer.data(lsh + 423);
    const auto *lsh_424 = buffer.data(lsh + 424);
    const auto *lsh_425 = buffer.data(lsh + 425);
    const auto *lsh_426 = buffer.data(lsh + 426);
    const auto *lsh_427 = buffer.data(lsh + 427);
    const auto *lsh_428 = buffer.data(lsh + 428);
    const auto *lsh_429 = buffer.data(lsh + 429);
    const auto *lsh_434 = buffer.data(lsh + 434);
    const auto *lsh_435 = buffer.data(lsh + 435);
    const auto *lsh_436 = buffer.data(lsh + 436);
    const auto *lsh_437 = buffer.data(lsh + 437);
    const auto *lsh_438 = buffer.data(lsh + 438);
    const auto *lsh_439 = buffer.data(lsh + 439);
    const auto *lsh_440 = buffer.data(lsh + 440);
    const auto *lsh_441 = buffer.data(lsh + 441);
    const auto *lsh_444 = buffer.data(lsh + 444);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, ksh_254, ksh_360, ksh_362, lsg0_258, \
                         lsg0_260, lsg1_258, lsg1_260, lsh_359, lsh_360, \
                         lsh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_13 * ksh_360[k]
                   + f_8 * lsg0_258[k]
                   - f_9 * lsg1_258[k]
                   + f_3 * pc_x[k] * lsh_360[k];

        t_480[k] = f_13 * ksh_254[k]
                   + f_3 * pc_y[k] * lsh_359[k];

        t_481[k] = f_13 * ksh_362[k]
                   + f_8 * lsg0_260[k]
                   - f_9 * lsg1_260[k]
                   + f_3 * pc_x[k] * lsh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, ksh_234, ksh_257, ksh_363, \
                         lsg0_261, lsg1_261, lsh_360, lsh_362, \
                         lsh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_13 * ksh_363[k]
                   + f_6 * lsg0_261[k]
                   - f_7 * lsg1_261[k]
                   + f_3 * pc_x[k] * lsh_363[k];

        t_483[k] = f_12 * ksh_234[k]
                   + f_3 * pc_z[k] * lsh_360[k];

        t_484[k] = f_13 * ksh_257[k]
                   + f_3 * pc_y[k] * lsh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, ksh_237, ksh_366, ksh_367, lsg0_264, \
                         lsg0_265, lsg1_264, lsg1_265, lsh_363, lsh_366, \
                         lsh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_13 * ksh_366[k]
                   + f_6 * lsg0_264[k]
                   - f_7 * lsg1_264[k]
                   + f_3 * pc_x[k] * lsh_366[k];

        t_486[k] = f_13 * ksh_367[k]
                   + f_4 * lsg0_265[k]
                   - f_5 * lsg1_265[k]
                   + f_3 * pc_x[k] * lsh_367[k];

        t_487[k] = f_12 * ksh_237[k]
                   + f_3 * pc_z[k] * lsh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, ksh_261, ksh_369, ksh_371, lsg0_267, \
                         lsg0_269, lsg1_267, lsg1_269, lsh_366, lsh_369, \
                         lsh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_13 * ksh_369[k]
                   + f_4 * lsg0_267[k]
                   - f_5 * lsg1_267[k]
                   + f_3 * pc_x[k] * lsh_369[k];

        t_489[k] = f_13 * ksh_261[k]
                   + f_3 * pc_y[k] * lsh_366[k];

        t_490[k] = f_13 * ksh_371[k]
                   + f_4 * lsg0_269[k]
                   - f_5 * lsg1_269[k]
                   + f_3 * pc_x[k] * lsh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, ksh_372, ksh_373, ksh_374, \
                         ksh_375, ksh_376, lsh_372, lsh_373, lsh_374, lsh_375, \
                         lsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_13 * ksh_372[k]
                   + f_3 * pc_x[k] * lsh_372[k];

        t_492[k] = f_13 * ksh_373[k]
                   + f_3 * pc_x[k] * lsh_373[k];

        t_493[k] = f_13 * ksh_374[k]
                   + f_3 * pc_x[k] * lsh_374[k];

        t_494[k] = f_13 * ksh_375[k]
                   + f_3 * pc_x[k] * lsh_375[k];

        t_495[k] = f_13 * ksh_376[k]
                   + f_3 * pc_x[k] * lsh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, ksh_246, ksh_267, ksh_377, \
                         lsg0_265, lsg1_265, lsh_372, lsh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_13 * ksh_377[k]
                   + f_3 * pc_x[k] * lsh_377[k];

        t_497[k] = f_13 * ksh_267[k]
                   + f_1 * lsg0_265[k]
                   - f_2 * lsg1_265[k]
                   + f_3 * pc_y[k] * lsh_372[k];

        t_498[k] = f_12 * ksh_246[k]
                   + f_3 * pc_z[k] * lsh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, ksh_269, ksh_270, ksh_271, lsg0_267, \
                         lsg0_268, lsg0_269, lsg1_267, lsg1_268, lsg1_269, lsh_374, lsh_375, \
                         lsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * ksh_269[k]
                   + f_8 * lsg0_267[k]
                   - f_9 * lsg1_267[k]
                   + f_3 * pc_y[k] * lsh_374[k];

        t_500[k] = f_13 * ksh_270[k]
                   + f_6 * lsg0_268[k]
                   - f_7 * lsg1_268[k]
                   + f_3 * pc_y[k] * lsh_375[k];

        t_501[k] = f_13 * ksh_271[k]
                   + f_4 * lsg0_269[k]
                   - f_5 * lsg1_269[k]
                   + f_3 * pc_y[k] * lsh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, ksh_251, ksh_272, ksh_378, \
                         lsg0_269, lsg0_270, lsg1_269, lsg1_270, lsh_377, \
                         lsh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * ksh_272[k]
                   + f_3 * pc_y[k] * lsh_377[k];

        t_503[k] = f_12 * ksh_251[k]
                   + f_1 * lsg0_269[k]
                   - f_2 * lsg1_269[k]
                   + f_3 * pc_z[k] * lsh_377[k];

        t_504[k] = f_13 * ksh_378[k]
                   + f_1 * lsg0_270[k]
                   - f_2 * lsg1_270[k]
                   + f_3 * pc_x[k] * lsh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, ksh_252, ksh_273, \
                         ksh_275, ksh_381, lsg0_273, lsg1_273, lsh_378, lsh_380, \
                         lsh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * ksh_273[k]
                   + f_3 * pc_y[k] * lsh_378[k];

        t_506[k] = f_13 * ksh_252[k]
                   + f_3 * pc_z[k] * lsh_378[k];

        t_507[k] = f_13 * ksh_381[k]
                   + f_8 * lsg0_273[k]
                   - f_9 * lsg1_273[k]
                   + f_3 * pc_x[k] * lsh_381[k];

        t_508[k] = f_12 * ksh_275[k]
                   + f_3 * pc_y[k] * lsh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, ksh_255, ksh_383, ksh_384, lsg0_275, \
                         lsg0_276, lsg1_275, lsg1_276, lsh_381, lsh_383, \
                         lsh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_13 * ksh_383[k]
                   + f_8 * lsg0_275[k]
                   - f_9 * lsg1_275[k]
                   + f_3 * pc_x[k] * lsh_383[k];

        t_510[k] = f_13 * ksh_384[k]
                   + f_6 * lsg0_276[k]
                   - f_7 * lsg1_276[k]
                   + f_3 * pc_x[k] * lsh_384[k];

        t_511[k] = f_13 * ksh_255[k]
                   + f_3 * pc_z[k] * lsh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, ksh_278, ksh_387, ksh_388, lsg0_279, \
                         lsg0_280, lsg1_279, lsg1_280, lsh_383, lsh_387, \
                         lsh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * ksh_278[k]
                   + f_3 * pc_y[k] * lsh_383[k];

        t_513[k] = f_13 * ksh_387[k]
                   + f_6 * lsg0_279[k]
                   - f_7 * lsg1_279[k]
                   + f_3 * pc_x[k] * lsh_387[k];

        t_514[k] = f_13 * ksh_388[k]
                   + f_4 * lsg0_280[k]
                   - f_5 * lsg1_280[k]
                   + f_3 * pc_x[k] * lsh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, ksh_258, ksh_282, ksh_390, \
                         lsg0_282, lsg1_282, lsh_384, lsh_387, \
                         lsh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * ksh_258[k]
                   + f_3 * pc_z[k] * lsh_384[k];

        t_516[k] = f_13 * ksh_390[k]
                   + f_4 * lsg0_282[k]
                   - f_5 * lsg1_282[k]
                   + f_3 * pc_x[k] * lsh_390[k];

        t_517[k] = f_12 * ksh_282[k]
                   + f_3 * pc_y[k] * lsh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, ksh_392, ksh_393, ksh_394, ksh_395, \
                         lsg0_284, lsg1_284, lsh_392, lsh_393, lsh_394, \
                         lsh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_13 * ksh_392[k]
                   + f_4 * lsg0_284[k]
                   - f_5 * lsg1_284[k]
                   + f_3 * pc_x[k] * lsh_392[k];

        t_519[k] = f_13 * ksh_393[k]
                   + f_3 * pc_x[k] * lsh_393[k];

        t_520[k] = f_13 * ksh_394[k]
                   + f_3 * pc_x[k] * lsh_394[k];

        t_521[k] = f_13 * ksh_395[k]
                   + f_3 * pc_x[k] * lsh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, ksh_288, ksh_396, ksh_397, \
                         ksh_398, lsg0_280, lsg1_280, lsh_393, lsh_396, lsh_397, \
                         lsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_13 * ksh_396[k]
                   + f_3 * pc_x[k] * lsh_396[k];

        t_523[k] = f_13 * ksh_397[k]
                   + f_3 * pc_x[k] * lsh_397[k];

        t_524[k] = f_13 * ksh_398[k]
                   + f_3 * pc_x[k] * lsh_398[k];

        t_525[k] = f_12 * ksh_288[k]
                   + f_1 * lsg0_280[k]
                   - f_2 * lsg1_280[k]
                   + f_3 * pc_y[k] * lsh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, ksh_267, ksh_290, ksh_291, lsg0_282, \
                         lsg0_283, lsg1_282, lsg1_283, lsh_393, lsh_395, \
                         lsh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * ksh_267[k]
                   + f_3 * pc_z[k] * lsh_393[k];

        t_527[k] = f_12 * ksh_290[k]
                   + f_8 * lsg0_282[k]
                   - f_9 * lsg1_282[k]
                   + f_3 * pc_y[k] * lsh_395[k];

        t_528[k] = f_12 * ksh_291[k]
                   + f_6 * lsg0_283[k]
                   - f_7 * lsg1_283[k]
                   + f_3 * pc_y[k] * lsh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, ksi0_392, ksh_272, \
                         ksh_292, ksh_293, ksi1_392, lsg0_284, lsg1_284, lsh_397, \
                         lsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * ksh_292[k]
                   + f_4 * lsg0_284[k]
                   - f_5 * lsg1_284[k]
                   + f_3 * pc_y[k] * lsh_397[k];

        t_530[k] = f_12 * ksh_293[k]
                   + f_3 * pc_y[k] * lsh_398[k];

        t_531[k] = f_13 * ksh_272[k]
                   + f_1 * lsg0_284[k]
                   - f_2 * lsg1_284[k]
                   + f_3 * pc_z[k] * lsh_398[k];

        t_532[k] = pa_y[k] * ksi0_392[k]
                   - f_10 * pc_y[k] * ksi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, ksi0_395, ksh_273, \
                         ksh_294, ksh_295, ksh_296, ksi1_395, lsh_399, \
                         lsh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * ksh_294[k]
                   + f_3 * pc_y[k] * lsh_399[k];

        t_534[k] = f_14 * ksh_273[k]
                   + f_3 * pc_z[k] * lsh_399[k];

        t_535[k] = pa_y[k] * ksi0_395[k]
                   + f_12 * ksh_295[k]
                   - f_10 * pc_y[k] * ksi1_395[k];

        t_536[k] = f_11 * ksh_296[k]
                   + f_3 * pc_y[k] * lsh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, ksi0_397, ksi0_398, \
                         ksh_276, ksh_297, ksh_299, ksi1_397, ksi1_398, lsh_402, \
                         lsh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * ksi0_397[k]
                   - f_10 * pc_y[k] * ksi1_397[k];

        t_538[k] = pa_y[k] * ksi0_398[k]
                   + f_13 * ksh_297[k]
                   - f_10 * pc_y[k] * ksi1_398[k];

        t_539[k] = f_14 * ksh_276[k]
                   + f_3 * pc_z[k] * lsh_402[k];

        t_540[k] = f_11 * ksh_299[k]
                   + f_3 * pc_y[k] * lsh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, ksi0_401, ksi0_402, ksh_279, \
                         ksh_300, ksi1_401, ksi1_402, lsh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * ksi0_401[k]
                   - f_10 * pc_y[k] * ksi1_401[k];

        t_542[k] = pa_y[k] * ksi0_402[k]
                   + f_14 * ksh_300[k]
                   - f_10 * pc_y[k] * ksi1_402[k];

        t_543[k] = f_14 * ksh_279[k]
                   + f_3 * pc_z[k] * lsh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, ksi0_404, ksi0_406, \
                         ksh_302, ksh_303, ksh_414, ksi1_404, ksi1_406, lsh_408, \
                         lsh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * ksi0_404[k]
                   + f_12 * ksh_302[k]
                   - f_10 * pc_y[k] * ksi1_404[k];

        t_545[k] = f_11 * ksh_303[k]
                   + f_3 * pc_y[k] * lsh_408[k];

        t_546[k] = pa_y[k] * ksi0_406[k]
                   - f_10 * pc_y[k] * ksi1_406[k];

        t_547[k] = f_13 * ksh_414[k]
                   + f_3 * pc_x[k] * lsh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, ksh_415, ksh_416, ksh_417, \
                         ksh_418, ksh_419, lsh_415, lsh_416, lsh_417, lsh_418, \
                         lsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_13 * ksh_415[k]
                   + f_3 * pc_x[k] * lsh_415[k];

        t_549[k] = f_13 * ksh_416[k]
                   + f_3 * pc_x[k] * lsh_416[k];

        t_550[k] = f_13 * ksh_417[k]
                   + f_3 * pc_x[k] * lsh_417[k];

        t_551[k] = f_13 * ksh_418[k]
                   + f_3 * pc_x[k] * lsh_418[k];

        t_552[k] = f_13 * ksh_419[k]
                   + f_3 * pc_x[k] * lsh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, ksh_288, ksh_309, ksh_311, lsg0_295, \
                         lsg0_297, lsg1_295, lsg1_297, lsh_414, \
                         lsh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * ksh_309[k]
                   + f_1 * lsg0_295[k]
                   - f_2 * lsg1_295[k]
                   + f_3 * pc_y[k] * lsh_414[k];

        t_554[k] = f_14 * ksh_288[k]
                   + f_3 * pc_z[k] * lsh_414[k];

        t_555[k] = f_11 * ksh_311[k]
                   + f_8 * lsg0_297[k]
                   - f_9 * lsg1_297[k]
                   + f_3 * pc_y[k] * lsh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, ksh_312, ksh_313, ksh_314, lsg0_298, \
                         lsg0_299, lsg1_298, lsg1_299, lsh_417, lsh_418, \
                         lsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * ksh_312[k]
                   + f_6 * lsg0_298[k]
                   - f_7 * lsg1_298[k]
                   + f_3 * pc_y[k] * lsh_417[k];

        t_557[k] = f_11 * ksh_313[k]
                   + f_4 * lsg0_299[k]
                   - f_5 * lsg1_299[k]
                   + f_3 * pc_y[k] * lsh_418[k];

        t_558[k] = f_11 * ksh_314[k]
                   + f_3 * pc_y[k] * lsh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, ksi0_419, \
                         ksh_294, ksh_420, ksi1_419, lsg0_300, lsg1_300, \
                         lsh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * ksi0_419[k]
                   - f_10 * pc_y[k] * ksi1_419[k];

        t_560[k] = f_13 * ksh_420[k]
                   + f_1 * lsg0_300[k]
                   - f_2 * lsg1_300[k]
                   + f_3 * pc_x[k] * lsh_420[k];

        t_561[k] = f_3 * pc_y[k] * lsh_420[k];

        t_562[k] = f_19 * ksh_294[k]
                   + f_3 * pc_z[k] * lsh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, ksh_425, lsg0_300, lsg0_305, \
                         lsg1_300, lsg1_305, lsh_421, lsh_422, \
                         lsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * lsg0_300[k]
                   - f_5 * lsg1_300[k]
                   + f_3 * pc_y[k] * lsh_421[k];

        t_564[k] = f_3 * pc_y[k] * lsh_422[k];

        t_565[k] = f_13 * ksh_425[k]
                   + f_8 * lsg0_305[k]
                   - f_9 * lsg1_305[k]
                   + f_3 * pc_x[k] * lsh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, lsg0_301, lsg0_302, lsg1_301, lsg1_302, \
                         lsh_423, lsh_424, lsh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * lsg0_301[k]
                   - f_7 * lsg1_301[k]
                   + f_3 * pc_y[k] * lsh_423[k];

        t_567[k] = f_4 * lsg0_302[k]
                   - f_5 * lsg1_302[k]
                   + f_3 * pc_y[k] * lsh_424[k];

        t_568[k] = f_3 * pc_y[k] * lsh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, ksh_429, lsg0_303, lsg0_304, \
                         lsg0_309, lsg1_303, lsg1_304, lsg1_309, lsh_426, lsh_427, \
                         lsh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_13 * ksh_429[k]
                   + f_6 * lsg0_309[k]
                   - f_7 * lsg1_309[k]
                   + f_3 * pc_x[k] * lsh_429[k];

        t_570[k] = f_8 * lsg0_303[k]
                   - f_9 * lsg1_303[k]
                   + f_3 * pc_y[k] * lsh_426[k];

        t_571[k] = f_6 * lsg0_304[k]
                   - f_7 * lsg1_304[k]
                   + f_3 * pc_y[k] * lsh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, ksh_434, ksh_435, lsg0_305, \
                         lsg0_314, lsg1_305, lsg1_314, lsh_428, lsh_429, lsh_434, \
                         lsh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * lsg0_305[k]
                   - f_5 * lsg1_305[k]
                   + f_3 * pc_y[k] * lsh_428[k];

        t_573[k] = f_3 * pc_y[k] * lsh_429[k];

        t_574[k] = f_13 * ksh_434[k]
                   + f_4 * lsg0_314[k]
                   - f_5 * lsg1_314[k]
                   + f_3 * pc_x[k] * lsh_434[k];

        t_575[k] = f_13 * ksh_435[k]
                   + f_3 * pc_x[k] * lsh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, ksh_436, ksh_437, \
                         ksh_438, ksh_440, lsh_434, lsh_436, lsh_437, lsh_438, \
                         lsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_13 * ksh_436[k]
                   + f_3 * pc_x[k] * lsh_436[k];

        t_577[k] = f_13 * ksh_437[k]
                   + f_3 * pc_x[k] * lsh_437[k];

        t_578[k] = f_13 * ksh_438[k]
                   + f_3 * pc_x[k] * lsh_438[k];

        t_579[k] = f_3 * pc_y[k] * lsh_434[k];

        t_580[k] = f_13 * ksh_440[k]
                   + f_3 * pc_x[k] * lsh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, lsg0_310, lsg0_311, lsg0_312, lsg1_310, \
                         lsg1_311, lsg1_312, lsh_435, lsh_436, \
                         lsh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * lsg0_310[k]
                   - f_2 * lsg1_310[k]
                   + f_3 * pc_y[k] * lsh_435[k];

        t_582[k] = f_16 * lsg0_311[k]
                   - f_17 * lsg1_311[k]
                   + f_3 * pc_y[k] * lsh_436[k];

        t_583[k] = f_8 * lsg0_312[k]
                   - f_9 * lsg1_312[k]
                   + f_3 * pc_y[k] * lsh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, ksh_314, lsg0_313, lsg0_314, \
                         lsg1_313, lsg1_314, lsh_438, lsh_439, \
                         lsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * lsg0_313[k]
                   - f_7 * lsg1_313[k]
                   + f_3 * pc_y[k] * lsh_438[k];

        t_585[k] = f_4 * lsg0_314[k]
                   - f_5 * lsg1_314[k]
                   + f_3 * pc_y[k] * lsh_439[k];

        t_586[k] = f_3 * pc_y[k] * lsh_440[k];

        t_587[k] = f_19 * ksh_314[k]
                   + f_1 * lsg0_314[k]
                   - f_2 * lsg1_314[k]
                   + f_3 * pc_z[k] * lsh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pc_x, pc_y, pc_z, ksh_315, ksh_441, \
                         ksh_444, lsg0_315, lsg0_318, lsg1_315, lsg1_318, lsh_441, \
                         lsh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_12 * ksh_441[k]
                   + f_1 * lsg0_315[k]
                   - f_2 * lsg1_315[k]
                   + f_3 * pc_x[k] * lsh_441[k];

        t_589[k] = f_18 * ksh_315[k]
                   + f_3 * pc_y[k] * lsh_441[k];

        t_590[k] = f_3 * pc_z[k] * lsh_441[k];

        t_591[k] = f_12 * ksh_444[k]
                   + f_8 * lsg0_318[k]
                   - f_9 * lsg1_318[k]
                   + f_3 * pc_x[k] * lsh_444[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_420 = buffer.data(ksi0 + 420);
    const auto *ksi0_423 = buffer.data(ksi0 + 423);
    const auto *ksi0_426 = buffer.data(ksi0 + 426);
    const auto *ksi0_430 = buffer.data(ksi0 + 430);
    const auto *ksi0_432 = buffer.data(ksi0 + 432);
    const auto *ksi0_441 = buffer.data(ksi0 + 441);

    const auto *ksh_315 = buffer.data(ksh + 315);
    const auto *ksh_318 = buffer.data(ksh + 318);
    const auto *ksh_320 = buffer.data(ksh + 320);
    const auto *ksh_321 = buffer.data(ksh + 321);
    const auto *ksh_322 = buffer.data(ksh + 322);
    const auto *ksh_324 = buffer.data(ksh + 324);
    const auto *ksh_330 = buffer.data(ksh + 330);
    const auto *ksh_335 = buffer.data(ksh + 335);
    const auto *ksh_336 = buffer.data(ksh + 336);
    const auto *ksh_338 = buffer.data(ksh + 338);
    const auto *ksh_339 = buffer.data(ksh + 339);
    const auto *ksh_341 = buffer.data(ksh + 341);
    const auto *ksh_342 = buffer.data(ksh + 342);
    const auto *ksh_345 = buffer.data(ksh + 345);
    const auto *ksh_351 = buffer.data(ksh + 351);
    const auto *ksh_353 = buffer.data(ksh + 353);
    const auto *ksh_354 = buffer.data(ksh + 354);
    const auto *ksh_355 = buffer.data(ksh + 355);
    const auto *ksh_356 = buffer.data(ksh + 356);
    const auto *ksh_357 = buffer.data(ksh + 357);
    const auto *ksh_359 = buffer.data(ksh + 359);
    const auto *ksh_360 = buffer.data(ksh + 360);
    const auto *ksh_362 = buffer.data(ksh + 362);
    const auto *ksh_363 = buffer.data(ksh + 363);
    const auto *ksh_366 = buffer.data(ksh + 366);
    const auto *ksh_372 = buffer.data(ksh + 372);
    const auto *ksh_374 = buffer.data(ksh + 374);
    const auto *ksh_375 = buffer.data(ksh + 375);
    const auto *ksh_376 = buffer.data(ksh + 376);
    const auto *ksh_377 = buffer.data(ksh + 377);
    const auto *ksh_378 = buffer.data(ksh + 378);
    const auto *ksh_380 = buffer.data(ksh + 380);
    const auto *ksh_383 = buffer.data(ksh + 383);
    const auto *ksh_387 = buffer.data(ksh + 387);
    const auto *ksh_393 = buffer.data(ksh + 393);
    const auto *ksh_395 = buffer.data(ksh + 395);
    const auto *ksh_396 = buffer.data(ksh + 396);
    const auto *ksh_397 = buffer.data(ksh + 397);
    const auto *ksh_398 = buffer.data(ksh + 398);
    const auto *ksh_399 = buffer.data(ksh + 399);
    const auto *ksh_447 = buffer.data(ksh + 447);
    const auto *ksh_451 = buffer.data(ksh + 451);
    const auto *ksh_456 = buffer.data(ksh + 456);
    const auto *ksh_458 = buffer.data(ksh + 458);
    const auto *ksh_459 = buffer.data(ksh + 459);
    const auto *ksh_460 = buffer.data(ksh + 460);
    const auto *ksh_461 = buffer.data(ksh + 461);
    const auto *ksh_467 = buffer.data(ksh + 467);
    const auto *ksh_471 = buffer.data(ksh + 471);
    const auto *ksh_476 = buffer.data(ksh + 476);
    const auto *ksh_477 = buffer.data(ksh + 477);
    const auto *ksh_478 = buffer.data(ksh + 478);
    const auto *ksh_479 = buffer.data(ksh + 479);
    const auto *ksh_480 = buffer.data(ksh + 480);
    const auto *ksh_481 = buffer.data(ksh + 481);
    const auto *ksh_482 = buffer.data(ksh + 482);
    const auto *ksh_483 = buffer.data(ksh + 483);
    const auto *ksh_486 = buffer.data(ksh + 486);
    const auto *ksh_488 = buffer.data(ksh + 488);
    const auto *ksh_489 = buffer.data(ksh + 489);
    const auto *ksh_492 = buffer.data(ksh + 492);
    const auto *ksh_493 = buffer.data(ksh + 493);
    const auto *ksh_495 = buffer.data(ksh + 495);
    const auto *ksh_497 = buffer.data(ksh + 497);
    const auto *ksh_498 = buffer.data(ksh + 498);
    const auto *ksh_499 = buffer.data(ksh + 499);
    const auto *ksh_500 = buffer.data(ksh + 500);
    const auto *ksh_501 = buffer.data(ksh + 501);
    const auto *ksh_502 = buffer.data(ksh + 502);
    const auto *ksh_503 = buffer.data(ksh + 503);
    const auto *ksh_504 = buffer.data(ksh + 504);
    const auto *ksh_507 = buffer.data(ksh + 507);
    const auto *ksh_509 = buffer.data(ksh + 509);
    const auto *ksh_510 = buffer.data(ksh + 510);
    const auto *ksh_513 = buffer.data(ksh + 513);
    const auto *ksh_514 = buffer.data(ksh + 514);
    const auto *ksh_516 = buffer.data(ksh + 516);
    const auto *ksh_518 = buffer.data(ksh + 518);
    const auto *ksh_519 = buffer.data(ksh + 519);
    const auto *ksh_520 = buffer.data(ksh + 520);
    const auto *ksh_521 = buffer.data(ksh + 521);
    const auto *ksh_522 = buffer.data(ksh + 522);
    const auto *ksh_523 = buffer.data(ksh + 523);
    const auto *ksh_524 = buffer.data(ksh + 524);
    const auto *ksh_525 = buffer.data(ksh + 525);

    const auto *ksi1_420 = buffer.data(ksi1 + 420);
    const auto *ksi1_423 = buffer.data(ksi1 + 423);
    const auto *ksi1_426 = buffer.data(ksi1 + 426);
    const auto *ksi1_430 = buffer.data(ksi1 + 430);
    const auto *ksi1_432 = buffer.data(ksi1 + 432);
    const auto *ksi1_441 = buffer.data(ksi1 + 441);

    const auto *lsg0_315 = buffer.data(lsg0 + 315);
    const auto *lsg0_317 = buffer.data(lsg0 + 317);
    const auto *lsg0_318 = buffer.data(lsg0 + 318);
    const auto *lsg0_320 = buffer.data(lsg0 + 320);
    const auto *lsg0_321 = buffer.data(lsg0 + 321);
    const auto *lsg0_325 = buffer.data(lsg0 + 325);
    const auto *lsg0_326 = buffer.data(lsg0 + 326);
    const auto *lsg0_327 = buffer.data(lsg0 + 327);
    const auto *lsg0_329 = buffer.data(lsg0 + 329);
    const auto *lsg0_335 = buffer.data(lsg0 + 335);
    const auto *lsg0_339 = buffer.data(lsg0 + 339);
    const auto *lsg0_342 = buffer.data(lsg0 + 342);
    const auto *lsg0_343 = buffer.data(lsg0 + 343);
    const auto *lsg0_344 = buffer.data(lsg0 + 344);
    const auto *lsg0_345 = buffer.data(lsg0 + 345);
    const auto *lsg0_348 = buffer.data(lsg0 + 348);
    const auto *lsg0_350 = buffer.data(lsg0 + 350);
    const auto *lsg0_351 = buffer.data(lsg0 + 351);
    const auto *lsg0_354 = buffer.data(lsg0 + 354);
    const auto *lsg0_355 = buffer.data(lsg0 + 355);
    const auto *lsg0_357 = buffer.data(lsg0 + 357);
    const auto *lsg0_358 = buffer.data(lsg0 + 358);
    const auto *lsg0_359 = buffer.data(lsg0 + 359);
    const auto *lsg0_360 = buffer.data(lsg0 + 360);
    const auto *lsg0_363 = buffer.data(lsg0 + 363);
    const auto *lsg0_365 = buffer.data(lsg0 + 365);
    const auto *lsg0_366 = buffer.data(lsg0 + 366);
    const auto *lsg0_369 = buffer.data(lsg0 + 369);
    const auto *lsg0_370 = buffer.data(lsg0 + 370);
    const auto *lsg0_372 = buffer.data(lsg0 + 372);
    const auto *lsg0_373 = buffer.data(lsg0 + 373);
    const auto *lsg0_374 = buffer.data(lsg0 + 374);
    const auto *lsg0_375 = buffer.data(lsg0 + 375);

    const auto *lsg1_315 = buffer.data(lsg1 + 315);
    const auto *lsg1_317 = buffer.data(lsg1 + 317);
    const auto *lsg1_318 = buffer.data(lsg1 + 318);
    const auto *lsg1_320 = buffer.data(lsg1 + 320);
    const auto *lsg1_321 = buffer.data(lsg1 + 321);
    const auto *lsg1_325 = buffer.data(lsg1 + 325);
    const auto *lsg1_326 = buffer.data(lsg1 + 326);
    const auto *lsg1_327 = buffer.data(lsg1 + 327);
    const auto *lsg1_329 = buffer.data(lsg1 + 329);
    const auto *lsg1_335 = buffer.data(lsg1 + 335);
    const auto *lsg1_339 = buffer.data(lsg1 + 339);
    const auto *lsg1_342 = buffer.data(lsg1 + 342);
    const auto *lsg1_343 = buffer.data(lsg1 + 343);
    const auto *lsg1_344 = buffer.data(lsg1 + 344);
    const auto *lsg1_345 = buffer.data(lsg1 + 345);
    const auto *lsg1_348 = buffer.data(lsg1 + 348);
    const auto *lsg1_350 = buffer.data(lsg1 + 350);
    const auto *lsg1_351 = buffer.data(lsg1 + 351);
    const auto *lsg1_354 = buffer.data(lsg1 + 354);
    const auto *lsg1_355 = buffer.data(lsg1 + 355);
    const auto *lsg1_357 = buffer.data(lsg1 + 357);
    const auto *lsg1_358 = buffer.data(lsg1 + 358);
    const auto *lsg1_359 = buffer.data(lsg1 + 359);
    const auto *lsg1_360 = buffer.data(lsg1 + 360);
    const auto *lsg1_363 = buffer.data(lsg1 + 363);
    const auto *lsg1_365 = buffer.data(lsg1 + 365);
    const auto *lsg1_366 = buffer.data(lsg1 + 366);
    const auto *lsg1_369 = buffer.data(lsg1 + 369);
    const auto *lsg1_370 = buffer.data(lsg1 + 370);
    const auto *lsg1_372 = buffer.data(lsg1 + 372);
    const auto *lsg1_373 = buffer.data(lsg1 + 373);
    const auto *lsg1_374 = buffer.data(lsg1 + 374);
    const auto *lsg1_375 = buffer.data(lsg1 + 375);

    const auto *lsh_442 = buffer.data(lsh + 442);
    const auto *lsh_443 = buffer.data(lsh + 443);
    const auto *lsh_444 = buffer.data(lsh + 444);
    const auto *lsh_446 = buffer.data(lsh + 446);
    const auto *lsh_447 = buffer.data(lsh + 447);
    const auto *lsh_448 = buffer.data(lsh + 448);
    const auto *lsh_450 = buffer.data(lsh + 450);
    const auto *lsh_451 = buffer.data(lsh + 451);
    const auto *lsh_456 = buffer.data(lsh + 456);
    const auto *lsh_457 = buffer.data(lsh + 457);
    const auto *lsh_458 = buffer.data(lsh + 458);
    const auto *lsh_459 = buffer.data(lsh + 459);
    const auto *lsh_460 = buffer.data(lsh + 460);
    const auto *lsh_461 = buffer.data(lsh + 461);
    const auto *lsh_462 = buffer.data(lsh + 462);
    const auto *lsh_464 = buffer.data(lsh + 464);
    const auto *lsh_465 = buffer.data(lsh + 465);
    const auto *lsh_467 = buffer.data(lsh + 467);
    const auto *lsh_468 = buffer.data(lsh + 468);
    const auto *lsh_471 = buffer.data(lsh + 471);
    const auto *lsh_476 = buffer.data(lsh + 476);
    const auto *lsh_477 = buffer.data(lsh + 477);
    const auto *lsh_478 = buffer.data(lsh + 478);
    const auto *lsh_479 = buffer.data(lsh + 479);
    const auto *lsh_480 = buffer.data(lsh + 480);
    const auto *lsh_481 = buffer.data(lsh + 481);
    const auto *lsh_482 = buffer.data(lsh + 482);
    const auto *lsh_483 = buffer.data(lsh + 483);
    const auto *lsh_485 = buffer.data(lsh + 485);
    const auto *lsh_486 = buffer.data(lsh + 486);
    const auto *lsh_488 = buffer.data(lsh + 488);
    const auto *lsh_489 = buffer.data(lsh + 489);
    const auto *lsh_492 = buffer.data(lsh + 492);
    const auto *lsh_493 = buffer.data(lsh + 493);
    const auto *lsh_495 = buffer.data(lsh + 495);
    const auto *lsh_497 = buffer.data(lsh + 497);
    const auto *lsh_498 = buffer.data(lsh + 498);
    const auto *lsh_499 = buffer.data(lsh + 499);
    const auto *lsh_500 = buffer.data(lsh + 500);
    const auto *lsh_501 = buffer.data(lsh + 501);
    const auto *lsh_502 = buffer.data(lsh + 502);
    const auto *lsh_503 = buffer.data(lsh + 503);
    const auto *lsh_504 = buffer.data(lsh + 504);
    const auto *lsh_506 = buffer.data(lsh + 506);
    const auto *lsh_507 = buffer.data(lsh + 507);
    const auto *lsh_509 = buffer.data(lsh + 509);
    const auto *lsh_510 = buffer.data(lsh + 510);
    const auto *lsh_513 = buffer.data(lsh + 513);
    const auto *lsh_514 = buffer.data(lsh + 514);
    const auto *lsh_516 = buffer.data(lsh + 516);
    const auto *lsh_518 = buffer.data(lsh + 518);
    const auto *lsh_519 = buffer.data(lsh + 519);
    const auto *lsh_520 = buffer.data(lsh + 520);
    const auto *lsh_521 = buffer.data(lsh + 521);
    const auto *lsh_522 = buffer.data(lsh + 522);
    const auto *lsh_523 = buffer.data(lsh + 523);
    const auto *lsh_524 = buffer.data(lsh + 524);
    const auto *lsh_525 = buffer.data(lsh + 525);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pc_x, pc_z, ksh_447, lsg0_315, lsg0_321, \
                         lsg1_315, lsg1_321, lsh_442, lsh_443, lsh_444, \
                         lsh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * lsh_442[k];

        t_593[k] = f_4 * lsg0_315[k]
                   - f_5 * lsg1_315[k]
                   + f_3 * pc_z[k] * lsh_443[k];

        t_594[k] = f_12 * ksh_447[k]
                   + f_6 * lsg0_321[k]
                   - f_7 * lsg1_321[k]
                   + f_3 * pc_x[k] * lsh_447[k];

        t_595[k] = f_3 * pc_z[k] * lsh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pc_x, pc_y, pc_z, ksh_320, ksh_451, \
                         lsg0_317, lsg0_325, lsg1_317, lsg1_325, lsh_446, lsh_447, \
                         lsh_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_18 * ksh_320[k]
                   + f_3 * pc_y[k] * lsh_446[k];

        t_597[k] = f_6 * lsg0_317[k]
                   - f_7 * lsg1_317[k]
                   + f_3 * pc_z[k] * lsh_446[k];

        t_598[k] = f_12 * ksh_451[k]
                   + f_4 * lsg0_325[k]
                   - f_5 * lsg1_325[k]
                   + f_3 * pc_x[k] * lsh_451[k];

        t_599[k] = f_3 * pc_z[k] * lsh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, ksh_324, ksh_456, \
                         lsg0_318, lsg0_320, lsg1_318, lsg1_320, lsh_448, lsh_450, \
                         lsh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * lsg0_318[k]
                   - f_5 * lsg1_318[k]
                   + f_3 * pc_z[k] * lsh_448[k];

        t_601[k] = f_18 * ksh_324[k]
                   + f_3 * pc_y[k] * lsh_450[k];

        t_602[k] = f_8 * lsg0_320[k]
                   - f_9 * lsg1_320[k]
                   + f_3 * pc_z[k] * lsh_450[k];

        t_603[k] = f_12 * ksh_456[k]
                   + f_3 * pc_x[k] * lsh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, ksh_458, ksh_459, \
                         ksh_460, ksh_461, lsh_451, lsh_458, lsh_459, lsh_460, \
                         lsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * lsh_451[k];

        t_605[k] = f_12 * ksh_458[k]
                   + f_3 * pc_x[k] * lsh_458[k];

        t_606[k] = f_12 * ksh_459[k]
                   + f_3 * pc_x[k] * lsh_459[k];

        t_607[k] = f_12 * ksh_460[k]
                   + f_3 * pc_x[k] * lsh_460[k];

        t_608[k] = f_12 * ksh_461[k]
                   + f_3 * pc_x[k] * lsh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pc_y, pc_z, ksh_330, lsg0_325, lsg0_326, \
                         lsg1_325, lsg1_326, lsh_456, lsh_457, \
                         lsh_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_18 * ksh_330[k]
                   + f_1 * lsg0_325[k]
                   - f_2 * lsg1_325[k]
                   + f_3 * pc_y[k] * lsh_456[k];

        t_610[k] = f_3 * pc_z[k] * lsh_456[k];

        t_611[k] = f_4 * lsg0_325[k]
                   - f_5 * lsg1_325[k]
                   + f_3 * pc_z[k] * lsh_457[k];

        t_612[k] = f_6 * lsg0_326[k]
                   - f_7 * lsg1_326[k]
                   + f_3 * pc_z[k] * lsh_458[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, t_616, pa_z, pc_y, pc_z, ksi0_420, ksh_335, \
                         ksi1_420, lsg0_327, lsg0_329, lsg1_327, lsg1_329, lsh_459, \
                         lsh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_8 * lsg0_327[k]
                   - f_9 * lsg1_327[k]
                   + f_3 * pc_z[k] * lsh_459[k];

        t_614[k] = f_18 * ksh_335[k]
                   + f_3 * pc_y[k] * lsh_461[k];

        t_615[k] = f_1 * lsg0_329[k]
                   - f_2 * lsg1_329[k]
                   + f_3 * pc_z[k] * lsh_461[k];

        t_616[k] = pa_z[k] * ksi0_420[k]
                   - f_10 * pc_z[k] * ksi1_420[k];
    }

#pragma omp simd aligned(t_617, t_618, t_619, t_620, pa_z, pc_y, pc_z, ksi0_423, ksh_315, \
                         ksh_336, ksh_338, ksi1_423, lsh_462, lsh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_617[k] = f_19 * ksh_336[k]
                   + f_3 * pc_y[k] * lsh_462[k];

        t_618[k] = f_11 * ksh_315[k]
                   + f_3 * pc_z[k] * lsh_462[k];

        t_619[k] = pa_z[k] * ksi0_423[k]
                   - f_10 * pc_z[k] * ksi1_423[k];

        t_620[k] = f_19 * ksh_338[k]
                   + f_3 * pc_y[k] * lsh_464[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, pa_z, pc_x, pc_z, ksi0_426, ksh_318, ksh_467, \
                         ksi1_426, lsg0_335, lsg1_335, lsh_465, \
                         lsh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_12 * ksh_467[k]
                   + f_8 * lsg0_335[k]
                   - f_9 * lsg1_335[k]
                   + f_3 * pc_x[k] * lsh_467[k];

        t_622[k] = pa_z[k] * ksi0_426[k]
                   - f_10 * pc_z[k] * ksi1_426[k];

        t_623[k] = f_11 * ksh_318[k]
                   + f_3 * pc_z[k] * lsh_465[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_z, pc_x, pc_y, pc_z, ksi0_430, ksh_341, \
                         ksh_471, ksi1_430, lsg0_339, lsg1_339, lsh_467, \
                         lsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_19 * ksh_341[k]
                   + f_3 * pc_y[k] * lsh_467[k];

        t_625[k] = f_12 * ksh_471[k]
                   + f_6 * lsg0_339[k]
                   - f_7 * lsg1_339[k]
                   + f_3 * pc_x[k] * lsh_471[k];

        t_626[k] = pa_z[k] * ksi0_430[k]
                   - f_10 * pc_z[k] * ksi1_430[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_z, pc_y, pc_z, ksi0_432, ksh_321, ksh_322, \
                         ksh_345, ksi1_432, lsh_468, lsh_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = f_11 * ksh_321[k]
                   + f_3 * pc_z[k] * lsh_468[k];

        t_628[k] = pa_z[k] * ksi0_432[k]
                   + f_12 * ksh_322[k]
                   - f_10 * pc_z[k] * ksi1_432[k];

        t_629[k] = f_19 * ksh_345[k]
                   + f_3 * pc_y[k] * lsh_471[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pc_x, ksh_476, ksh_477, ksh_478, ksh_479, \
                         lsg0_344, lsg1_344, lsh_476, lsh_477, lsh_478, \
                         lsh_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = f_12 * ksh_476[k]
                   + f_4 * lsg0_344[k]
                   - f_5 * lsg1_344[k]
                   + f_3 * pc_x[k] * lsh_476[k];

        t_631[k] = f_12 * ksh_477[k]
                   + f_3 * pc_x[k] * lsh_477[k];

        t_632[k] = f_12 * ksh_478[k]
                   + f_3 * pc_x[k] * lsh_478[k];

        t_633[k] = f_12 * ksh_479[k]
                   + f_3 * pc_x[k] * lsh_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_z, pc_x, pc_z, ksi0_441, ksh_480, \
                         ksh_481, ksh_482, ksi1_441, lsh_480, lsh_481, \
                         lsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_12 * ksh_480[k]
                   + f_3 * pc_x[k] * lsh_480[k];

        t_635[k] = f_12 * ksh_481[k]
                   + f_3 * pc_x[k] * lsh_481[k];

        t_636[k] = f_12 * ksh_482[k]
                   + f_3 * pc_x[k] * lsh_482[k];

        t_637[k] = pa_z[k] * ksi0_441[k]
                   - f_10 * pc_z[k] * ksi1_441[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, pc_y, pc_z, ksh_330, ksh_353, ksh_354, lsg0_342, \
                         lsg0_343, lsg1_342, lsg1_343, lsh_477, lsh_479, \
                         lsh_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_11 * ksh_330[k]
                   + f_3 * pc_z[k] * lsh_477[k];

        t_639[k] = f_19 * ksh_353[k]
                   + f_8 * lsg0_342[k]
                   - f_9 * lsg1_342[k]
                   + f_3 * pc_y[k] * lsh_479[k];

        t_640[k] = f_19 * ksh_354[k]
                   + f_6 * lsg0_343[k]
                   - f_7 * lsg1_343[k]
                   + f_3 * pc_y[k] * lsh_480[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, pc_y, pc_z, ksh_335, ksh_355, ksh_356, lsg0_344, \
                         lsg1_344, lsh_481, lsh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = f_19 * ksh_355[k]
                   + f_4 * lsg0_344[k]
                   - f_5 * lsg1_344[k]
                   + f_3 * pc_y[k] * lsh_481[k];

        t_642[k] = f_19 * ksh_356[k]
                   + f_3 * pc_y[k] * lsh_482[k];

        t_643[k] = f_11 * ksh_335[k]
                   + f_1 * lsg0_344[k]
                   - f_2 * lsg1_344[k]
                   + f_3 * pc_z[k] * lsh_482[k];
    }

#pragma omp simd aligned(t_644, t_645, t_646, pc_x, pc_y, pc_z, ksh_336, ksh_357, ksh_483, \
                         lsg0_345, lsg1_345, lsh_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_644[k] = f_12 * ksh_483[k]
                   + f_1 * lsg0_345[k]
                   - f_2 * lsg1_345[k]
                   + f_3 * pc_x[k] * lsh_483[k];

        t_645[k] = f_14 * ksh_357[k]
                   + f_3 * pc_y[k] * lsh_483[k];

        t_646[k] = f_12 * ksh_336[k]
                   + f_3 * pc_z[k] * lsh_483[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, pc_x, pc_y, ksh_359, ksh_486, ksh_488, lsg0_348, \
                         lsg0_350, lsg1_348, lsg1_350, lsh_485, lsh_486, \
                         lsh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = f_12 * ksh_486[k]
                   + f_8 * lsg0_348[k]
                   - f_9 * lsg1_348[k]
                   + f_3 * pc_x[k] * lsh_486[k];

        t_648[k] = f_14 * ksh_359[k]
                   + f_3 * pc_y[k] * lsh_485[k];

        t_649[k] = f_12 * ksh_488[k]
                   + f_8 * lsg0_350[k]
                   - f_9 * lsg1_350[k]
                   + f_3 * pc_x[k] * lsh_488[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, pc_x, pc_y, pc_z, ksh_339, ksh_362, ksh_489, \
                         lsg0_351, lsg1_351, lsh_486, lsh_488, \
                         lsh_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = f_12 * ksh_489[k]
                   + f_6 * lsg0_351[k]
                   - f_7 * lsg1_351[k]
                   + f_3 * pc_x[k] * lsh_489[k];

        t_651[k] = f_12 * ksh_339[k]
                   + f_3 * pc_z[k] * lsh_486[k];

        t_652[k] = f_14 * ksh_362[k]
                   + f_3 * pc_y[k] * lsh_488[k];
    }

#pragma omp simd aligned(t_653, t_654, t_655, pc_x, pc_z, ksh_342, ksh_492, ksh_493, lsg0_354, \
                         lsg0_355, lsg1_354, lsg1_355, lsh_489, lsh_492, \
                         lsh_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_653[k] = f_12 * ksh_492[k]
                   + f_6 * lsg0_354[k]
                   - f_7 * lsg1_354[k]
                   + f_3 * pc_x[k] * lsh_492[k];

        t_654[k] = f_12 * ksh_493[k]
                   + f_4 * lsg0_355[k]
                   - f_5 * lsg1_355[k]
                   + f_3 * pc_x[k] * lsh_493[k];

        t_655[k] = f_12 * ksh_342[k]
                   + f_3 * pc_z[k] * lsh_489[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, pc_x, pc_y, ksh_366, ksh_495, ksh_497, lsg0_357, \
                         lsg0_359, lsg1_357, lsg1_359, lsh_492, lsh_495, \
                         lsh_497 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_12 * ksh_495[k]
                   + f_4 * lsg0_357[k]
                   - f_5 * lsg1_357[k]
                   + f_3 * pc_x[k] * lsh_495[k];

        t_657[k] = f_14 * ksh_366[k]
                   + f_3 * pc_y[k] * lsh_492[k];

        t_658[k] = f_12 * ksh_497[k]
                   + f_4 * lsg0_359[k]
                   - f_5 * lsg1_359[k]
                   + f_3 * pc_x[k] * lsh_497[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, t_663, pc_x, ksh_498, ksh_499, ksh_500, \
                         ksh_501, ksh_502, lsh_498, lsh_499, lsh_500, lsh_501, \
                         lsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = f_12 * ksh_498[k]
                   + f_3 * pc_x[k] * lsh_498[k];

        t_660[k] = f_12 * ksh_499[k]
                   + f_3 * pc_x[k] * lsh_499[k];

        t_661[k] = f_12 * ksh_500[k]
                   + f_3 * pc_x[k] * lsh_500[k];

        t_662[k] = f_12 * ksh_501[k]
                   + f_3 * pc_x[k] * lsh_501[k];

        t_663[k] = f_12 * ksh_502[k]
                   + f_3 * pc_x[k] * lsh_502[k];
    }

#pragma omp simd aligned(t_664, t_665, t_666, pc_x, pc_y, pc_z, ksh_351, ksh_372, ksh_503, \
                         lsg0_355, lsg1_355, lsh_498, lsh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_664[k] = f_12 * ksh_503[k]
                   + f_3 * pc_x[k] * lsh_503[k];

        t_665[k] = f_14 * ksh_372[k]
                   + f_1 * lsg0_355[k]
                   - f_2 * lsg1_355[k]
                   + f_3 * pc_y[k] * lsh_498[k];

        t_666[k] = f_12 * ksh_351[k]
                   + f_3 * pc_z[k] * lsh_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, ksh_374, ksh_375, ksh_376, lsg0_357, \
                         lsg0_358, lsg0_359, lsg1_357, lsg1_358, lsg1_359, lsh_500, lsh_501, \
                         lsh_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * ksh_374[k]
                   + f_8 * lsg0_357[k]
                   - f_9 * lsg1_357[k]
                   + f_3 * pc_y[k] * lsh_500[k];

        t_668[k] = f_14 * ksh_375[k]
                   + f_6 * lsg0_358[k]
                   - f_7 * lsg1_358[k]
                   + f_3 * pc_y[k] * lsh_501[k];

        t_669[k] = f_14 * ksh_376[k]
                   + f_4 * lsg0_359[k]
                   - f_5 * lsg1_359[k]
                   + f_3 * pc_y[k] * lsh_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, ksh_356, ksh_377, ksh_504, \
                         lsg0_359, lsg0_360, lsg1_359, lsg1_360, lsh_503, \
                         lsh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * ksh_377[k]
                   + f_3 * pc_y[k] * lsh_503[k];

        t_671[k] = f_12 * ksh_356[k]
                   + f_1 * lsg0_359[k]
                   - f_2 * lsg1_359[k]
                   + f_3 * pc_z[k] * lsh_503[k];

        t_672[k] = f_12 * ksh_504[k]
                   + f_1 * lsg0_360[k]
                   - f_2 * lsg1_360[k]
                   + f_3 * pc_x[k] * lsh_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, t_676, pc_x, pc_y, pc_z, ksh_357, ksh_378, \
                         ksh_380, ksh_507, lsg0_363, lsg1_363, lsh_504, lsh_506, \
                         lsh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_13 * ksh_378[k]
                   + f_3 * pc_y[k] * lsh_504[k];

        t_674[k] = f_13 * ksh_357[k]
                   + f_3 * pc_z[k] * lsh_504[k];

        t_675[k] = f_12 * ksh_507[k]
                   + f_8 * lsg0_363[k]
                   - f_9 * lsg1_363[k]
                   + f_3 * pc_x[k] * lsh_507[k];

        t_676[k] = f_13 * ksh_380[k]
                   + f_3 * pc_y[k] * lsh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_x, pc_z, ksh_360, ksh_509, ksh_510, lsg0_365, \
                         lsg0_366, lsg1_365, lsg1_366, lsh_507, lsh_509, \
                         lsh_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_12 * ksh_509[k]
                   + f_8 * lsg0_365[k]
                   - f_9 * lsg1_365[k]
                   + f_3 * pc_x[k] * lsh_509[k];

        t_678[k] = f_12 * ksh_510[k]
                   + f_6 * lsg0_366[k]
                   - f_7 * lsg1_366[k]
                   + f_3 * pc_x[k] * lsh_510[k];

        t_679[k] = f_13 * ksh_360[k]
                   + f_3 * pc_z[k] * lsh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pc_x, pc_y, ksh_383, ksh_513, ksh_514, lsg0_369, \
                         lsg0_370, lsg1_369, lsg1_370, lsh_509, lsh_513, \
                         lsh_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * ksh_383[k]
                   + f_3 * pc_y[k] * lsh_509[k];

        t_681[k] = f_12 * ksh_513[k]
                   + f_6 * lsg0_369[k]
                   - f_7 * lsg1_369[k]
                   + f_3 * pc_x[k] * lsh_513[k];

        t_682[k] = f_12 * ksh_514[k]
                   + f_4 * lsg0_370[k]
                   - f_5 * lsg1_370[k]
                   + f_3 * pc_x[k] * lsh_514[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pc_x, pc_y, pc_z, ksh_363, ksh_387, ksh_516, \
                         lsg0_372, lsg1_372, lsh_510, lsh_513, \
                         lsh_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * ksh_363[k]
                   + f_3 * pc_z[k] * lsh_510[k];

        t_684[k] = f_12 * ksh_516[k]
                   + f_4 * lsg0_372[k]
                   - f_5 * lsg1_372[k]
                   + f_3 * pc_x[k] * lsh_516[k];

        t_685[k] = f_13 * ksh_387[k]
                   + f_3 * pc_y[k] * lsh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pc_x, ksh_518, ksh_519, ksh_520, ksh_521, \
                         lsg0_374, lsg1_374, lsh_518, lsh_519, lsh_520, \
                         lsh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_12 * ksh_518[k]
                   + f_4 * lsg0_374[k]
                   - f_5 * lsg1_374[k]
                   + f_3 * pc_x[k] * lsh_518[k];

        t_687[k] = f_12 * ksh_519[k]
                   + f_3 * pc_x[k] * lsh_519[k];

        t_688[k] = f_12 * ksh_520[k]
                   + f_3 * pc_x[k] * lsh_520[k];

        t_689[k] = f_12 * ksh_521[k]
                   + f_3 * pc_x[k] * lsh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pc_x, pc_y, ksh_393, ksh_522, ksh_523, \
                         ksh_524, lsg0_370, lsg1_370, lsh_519, lsh_522, lsh_523, \
                         lsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_12 * ksh_522[k]
                   + f_3 * pc_x[k] * lsh_522[k];

        t_691[k] = f_12 * ksh_523[k]
                   + f_3 * pc_x[k] * lsh_523[k];

        t_692[k] = f_12 * ksh_524[k]
                   + f_3 * pc_x[k] * lsh_524[k];

        t_693[k] = f_13 * ksh_393[k]
                   + f_1 * lsg0_370[k]
                   - f_2 * lsg1_370[k]
                   + f_3 * pc_y[k] * lsh_519[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pc_y, pc_z, ksh_372, ksh_395, ksh_396, lsg0_372, \
                         lsg0_373, lsg1_372, lsg1_373, lsh_519, lsh_521, \
                         lsh_522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * ksh_372[k]
                   + f_3 * pc_z[k] * lsh_519[k];

        t_695[k] = f_13 * ksh_395[k]
                   + f_8 * lsg0_372[k]
                   - f_9 * lsg1_372[k]
                   + f_3 * pc_y[k] * lsh_521[k];

        t_696[k] = f_13 * ksh_396[k]
                   + f_6 * lsg0_373[k]
                   - f_7 * lsg1_373[k]
                   + f_3 * pc_y[k] * lsh_522[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, ksh_377, ksh_397, ksh_398, lsg0_374, \
                         lsg1_374, lsh_523, lsh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_13 * ksh_397[k]
                   + f_4 * lsg0_374[k]
                   - f_5 * lsg1_374[k]
                   + f_3 * pc_y[k] * lsh_523[k];

        t_698[k] = f_13 * ksh_398[k]
                   + f_3 * pc_y[k] * lsh_524[k];

        t_699[k] = f_13 * ksh_377[k]
                   + f_1 * lsg0_374[k]
                   - f_2 * lsg1_374[k]
                   + f_3 * pc_z[k] * lsh_524[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, pc_y, pc_z, ksh_378, ksh_399, ksh_525, \
                         lsg0_375, lsg1_375, lsh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_12 * ksh_525[k]
                   + f_1 * lsg0_375[k]
                   - f_2 * lsg1_375[k]
                   + f_3 * pc_x[k] * lsh_525[k];

        t_701[k] = f_12 * ksh_399[k]
                   + f_3 * pc_y[k] * lsh_525[k];

        t_702[k] = f_14 * ksh_378[k]
                   + f_3 * pc_z[k] * lsh_525[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_560 = buffer.data(ksi0 + 560);
    const auto *ksi0_563 = buffer.data(ksi0 + 563);
    const auto *ksi0_565 = buffer.data(ksi0 + 565);
    const auto *ksi0_566 = buffer.data(ksi0 + 566);
    const auto *ksi0_569 = buffer.data(ksi0 + 569);
    const auto *ksi0_570 = buffer.data(ksi0 + 570);
    const auto *ksi0_572 = buffer.data(ksi0 + 572);
    const auto *ksi0_574 = buffer.data(ksi0 + 574);
    const auto *ksi0_587 = buffer.data(ksi0 + 587);
    const auto *ksi0_588 = buffer.data(ksi0 + 588);
    const auto *ksi0_591 = buffer.data(ksi0 + 591);
    const auto *ksi0_594 = buffer.data(ksi0 + 594);
    const auto *ksi0_784 = buffer.data(ksi0 + 784);
    const auto *ksi0_787 = buffer.data(ksi0 + 787);
    const auto *ksi0_790 = buffer.data(ksi0 + 790);
    const auto *ksi0_794 = buffer.data(ksi0 + 794);
    const auto *ksi0_805 = buffer.data(ksi0 + 805);
    const auto *ksi0_807 = buffer.data(ksi0 + 807);
    const auto *ksi0_808 = buffer.data(ksi0 + 808);
    const auto *ksi0_809 = buffer.data(ksi0 + 809);
    const auto *ksi0_811 = buffer.data(ksi0 + 811);
    const auto *ksi0_817 = buffer.data(ksi0 + 817);
    const auto *ksi0_821 = buffer.data(ksi0 + 821);

    const auto *ksh_381 = buffer.data(ksh + 381);
    const auto *ksh_384 = buffer.data(ksh + 384);
    const auto *ksh_393 = buffer.data(ksh + 393);
    const auto *ksh_398 = buffer.data(ksh + 398);
    const auto *ksh_399 = buffer.data(ksh + 399);
    const auto *ksh_401 = buffer.data(ksh + 401);
    const auto *ksh_402 = buffer.data(ksh + 402);
    const auto *ksh_404 = buffer.data(ksh + 404);
    const auto *ksh_405 = buffer.data(ksh + 405);
    const auto *ksh_408 = buffer.data(ksh + 408);
    const auto *ksh_414 = buffer.data(ksh + 414);
    const auto *ksh_416 = buffer.data(ksh + 416);
    const auto *ksh_417 = buffer.data(ksh + 417);
    const auto *ksh_418 = buffer.data(ksh + 418);
    const auto *ksh_419 = buffer.data(ksh + 419);
    const auto *ksh_420 = buffer.data(ksh + 420);
    const auto *ksh_421 = buffer.data(ksh + 421);
    const auto *ksh_422 = buffer.data(ksh + 422);
    const auto *ksh_423 = buffer.data(ksh + 423);
    const auto *ksh_425 = buffer.data(ksh + 425);
    const auto *ksh_426 = buffer.data(ksh + 426);
    const auto *ksh_428 = buffer.data(ksh + 428);
    const auto *ksh_429 = buffer.data(ksh + 429);
    const auto *ksh_435 = buffer.data(ksh + 435);
    const auto *ksh_437 = buffer.data(ksh + 437);
    const auto *ksh_438 = buffer.data(ksh + 438);
    const auto *ksh_439 = buffer.data(ksh + 439);
    const auto *ksh_440 = buffer.data(ksh + 440);
    const auto *ksh_441 = buffer.data(ksh + 441);
    const auto *ksh_444 = buffer.data(ksh + 444);
    const auto *ksh_446 = buffer.data(ksh + 446);
    const auto *ksh_450 = buffer.data(ksh + 450);
    const auto *ksh_461 = buffer.data(ksh + 461);
    const auto *ksh_462 = buffer.data(ksh + 462);
    const auto *ksh_464 = buffer.data(ksh + 464);
    const auto *ksh_467 = buffer.data(ksh + 467);
    const auto *ksh_528 = buffer.data(ksh + 528);
    const auto *ksh_530 = buffer.data(ksh + 530);
    const auto *ksh_531 = buffer.data(ksh + 531);
    const auto *ksh_534 = buffer.data(ksh + 534);
    const auto *ksh_535 = buffer.data(ksh + 535);
    const auto *ksh_537 = buffer.data(ksh + 537);
    const auto *ksh_539 = buffer.data(ksh + 539);
    const auto *ksh_540 = buffer.data(ksh + 540);
    const auto *ksh_541 = buffer.data(ksh + 541);
    const auto *ksh_542 = buffer.data(ksh + 542);
    const auto *ksh_543 = buffer.data(ksh + 543);
    const auto *ksh_544 = buffer.data(ksh + 544);
    const auto *ksh_545 = buffer.data(ksh + 545);
    const auto *ksh_561 = buffer.data(ksh + 561);
    const auto *ksh_562 = buffer.data(ksh + 562);
    const auto *ksh_563 = buffer.data(ksh + 563);
    const auto *ksh_564 = buffer.data(ksh + 564);
    const auto *ksh_565 = buffer.data(ksh + 565);
    const auto *ksh_566 = buffer.data(ksh + 566);
    const auto *ksh_567 = buffer.data(ksh + 567);
    const auto *ksh_572 = buffer.data(ksh + 572);
    const auto *ksh_576 = buffer.data(ksh + 576);
    const auto *ksh_581 = buffer.data(ksh + 581);
    const auto *ksh_582 = buffer.data(ksh + 582);
    const auto *ksh_583 = buffer.data(ksh + 583);
    const auto *ksh_584 = buffer.data(ksh + 584);
    const auto *ksh_585 = buffer.data(ksh + 585);
    const auto *ksh_587 = buffer.data(ksh + 587);
    const auto *ksh_588 = buffer.data(ksh + 588);
    const auto *ksh_591 = buffer.data(ksh + 591);
    const auto *ksh_594 = buffer.data(ksh + 594);
    const auto *ksh_598 = buffer.data(ksh + 598);
    const auto *ksh_603 = buffer.data(ksh + 603);
    const auto *ksh_605 = buffer.data(ksh + 605);
    const auto *ksh_606 = buffer.data(ksh + 606);
    const auto *ksh_607 = buffer.data(ksh + 607);
    const auto *ksh_608 = buffer.data(ksh + 608);
    const auto *ksh_614 = buffer.data(ksh + 614);
    const auto *ksh_618 = buffer.data(ksh + 618);

    const auto *ksi1_560 = buffer.data(ksi1 + 560);
    const auto *ksi1_563 = buffer.data(ksi1 + 563);
    const auto *ksi1_565 = buffer.data(ksi1 + 565);
    const auto *ksi1_566 = buffer.data(ksi1 + 566);
    const auto *ksi1_569 = buffer.data(ksi1 + 569);
    const auto *ksi1_570 = buffer.data(ksi1 + 570);
    const auto *ksi1_572 = buffer.data(ksi1 + 572);
    const auto *ksi1_574 = buffer.data(ksi1 + 574);
    const auto *ksi1_587 = buffer.data(ksi1 + 587);
    const auto *ksi1_588 = buffer.data(ksi1 + 588);
    const auto *ksi1_591 = buffer.data(ksi1 + 591);
    const auto *ksi1_594 = buffer.data(ksi1 + 594);
    const auto *ksi1_784 = buffer.data(ksi1 + 784);
    const auto *ksi1_787 = buffer.data(ksi1 + 787);
    const auto *ksi1_790 = buffer.data(ksi1 + 790);
    const auto *ksi1_794 = buffer.data(ksi1 + 794);
    const auto *ksi1_805 = buffer.data(ksi1 + 805);
    const auto *ksi1_807 = buffer.data(ksi1 + 807);
    const auto *ksi1_808 = buffer.data(ksi1 + 808);
    const auto *ksi1_809 = buffer.data(ksi1 + 809);
    const auto *ksi1_811 = buffer.data(ksi1 + 811);
    const auto *ksi1_817 = buffer.data(ksi1 + 817);
    const auto *ksi1_821 = buffer.data(ksi1 + 821);

    const auto *lsg0_378 = buffer.data(lsg0 + 378);
    const auto *lsg0_380 = buffer.data(lsg0 + 380);
    const auto *lsg0_381 = buffer.data(lsg0 + 381);
    const auto *lsg0_384 = buffer.data(lsg0 + 384);
    const auto *lsg0_385 = buffer.data(lsg0 + 385);
    const auto *lsg0_387 = buffer.data(lsg0 + 387);
    const auto *lsg0_388 = buffer.data(lsg0 + 388);
    const auto *lsg0_389 = buffer.data(lsg0 + 389);
    const auto *lsg0_400 = buffer.data(lsg0 + 400);
    const auto *lsg0_402 = buffer.data(lsg0 + 402);
    const auto *lsg0_403 = buffer.data(lsg0 + 403);
    const auto *lsg0_404 = buffer.data(lsg0 + 404);
    const auto *lsg0_405 = buffer.data(lsg0 + 405);
    const auto *lsg0_406 = buffer.data(lsg0 + 406);
    const auto *lsg0_407 = buffer.data(lsg0 + 407);
    const auto *lsg0_408 = buffer.data(lsg0 + 408);
    const auto *lsg0_409 = buffer.data(lsg0 + 409);
    const auto *lsg0_410 = buffer.data(lsg0 + 410);
    const auto *lsg0_414 = buffer.data(lsg0 + 414);
    const auto *lsg0_415 = buffer.data(lsg0 + 415);
    const auto *lsg0_416 = buffer.data(lsg0 + 416);
    const auto *lsg0_417 = buffer.data(lsg0 + 417);
    const auto *lsg0_418 = buffer.data(lsg0 + 418);
    const auto *lsg0_419 = buffer.data(lsg0 + 419);
    const auto *lsg0_420 = buffer.data(lsg0 + 420);
    const auto *lsg0_422 = buffer.data(lsg0 + 422);
    const auto *lsg0_423 = buffer.data(lsg0 + 423);
    const auto *lsg0_425 = buffer.data(lsg0 + 425);

    const auto *lsg1_378 = buffer.data(lsg1 + 378);
    const auto *lsg1_380 = buffer.data(lsg1 + 380);
    const auto *lsg1_381 = buffer.data(lsg1 + 381);
    const auto *lsg1_384 = buffer.data(lsg1 + 384);
    const auto *lsg1_385 = buffer.data(lsg1 + 385);
    const auto *lsg1_387 = buffer.data(lsg1 + 387);
    const auto *lsg1_388 = buffer.data(lsg1 + 388);
    const auto *lsg1_389 = buffer.data(lsg1 + 389);
    const auto *lsg1_400 = buffer.data(lsg1 + 400);
    const auto *lsg1_402 = buffer.data(lsg1 + 402);
    const auto *lsg1_403 = buffer.data(lsg1 + 403);
    const auto *lsg1_404 = buffer.data(lsg1 + 404);
    const auto *lsg1_405 = buffer.data(lsg1 + 405);
    const auto *lsg1_406 = buffer.data(lsg1 + 406);
    const auto *lsg1_407 = buffer.data(lsg1 + 407);
    const auto *lsg1_408 = buffer.data(lsg1 + 408);
    const auto *lsg1_409 = buffer.data(lsg1 + 409);
    const auto *lsg1_410 = buffer.data(lsg1 + 410);
    const auto *lsg1_414 = buffer.data(lsg1 + 414);
    const auto *lsg1_415 = buffer.data(lsg1 + 415);
    const auto *lsg1_416 = buffer.data(lsg1 + 416);
    const auto *lsg1_417 = buffer.data(lsg1 + 417);
    const auto *lsg1_418 = buffer.data(lsg1 + 418);
    const auto *lsg1_419 = buffer.data(lsg1 + 419);
    const auto *lsg1_420 = buffer.data(lsg1 + 420);
    const auto *lsg1_422 = buffer.data(lsg1 + 422);
    const auto *lsg1_423 = buffer.data(lsg1 + 423);
    const auto *lsg1_425 = buffer.data(lsg1 + 425);

    const auto *lsh_527 = buffer.data(lsh + 527);
    const auto *lsh_528 = buffer.data(lsh + 528);
    const auto *lsh_530 = buffer.data(lsh + 530);
    const auto *lsh_531 = buffer.data(lsh + 531);
    const auto *lsh_534 = buffer.data(lsh + 534);
    const auto *lsh_535 = buffer.data(lsh + 535);
    const auto *lsh_537 = buffer.data(lsh + 537);
    const auto *lsh_539 = buffer.data(lsh + 539);
    const auto *lsh_540 = buffer.data(lsh + 540);
    const auto *lsh_541 = buffer.data(lsh + 541);
    const auto *lsh_542 = buffer.data(lsh + 542);
    const auto *lsh_543 = buffer.data(lsh + 543);
    const auto *lsh_544 = buffer.data(lsh + 544);
    const auto *lsh_545 = buffer.data(lsh + 545);
    const auto *lsh_546 = buffer.data(lsh + 546);
    const auto *lsh_548 = buffer.data(lsh + 548);
    const auto *lsh_549 = buffer.data(lsh + 549);
    const auto *lsh_551 = buffer.data(lsh + 551);
    const auto *lsh_552 = buffer.data(lsh + 552);
    const auto *lsh_555 = buffer.data(lsh + 555);
    const auto *lsh_561 = buffer.data(lsh + 561);
    const auto *lsh_562 = buffer.data(lsh + 562);
    const auto *lsh_563 = buffer.data(lsh + 563);
    const auto *lsh_564 = buffer.data(lsh + 564);
    const auto *lsh_565 = buffer.data(lsh + 565);
    const auto *lsh_566 = buffer.data(lsh + 566);
    const auto *lsh_567 = buffer.data(lsh + 567);
    const auto *lsh_568 = buffer.data(lsh + 568);
    const auto *lsh_569 = buffer.data(lsh + 569);
    const auto *lsh_570 = buffer.data(lsh + 570);
    const auto *lsh_571 = buffer.data(lsh + 571);
    const auto *lsh_572 = buffer.data(lsh + 572);
    const auto *lsh_573 = buffer.data(lsh + 573);
    const auto *lsh_574 = buffer.data(lsh + 574);
    const auto *lsh_575 = buffer.data(lsh + 575);
    const auto *lsh_576 = buffer.data(lsh + 576);
    const auto *lsh_581 = buffer.data(lsh + 581);
    const auto *lsh_582 = buffer.data(lsh + 582);
    const auto *lsh_583 = buffer.data(lsh + 583);
    const auto *lsh_584 = buffer.data(lsh + 584);
    const auto *lsh_585 = buffer.data(lsh + 585);
    const auto *lsh_586 = buffer.data(lsh + 586);
    const auto *lsh_587 = buffer.data(lsh + 587);
    const auto *lsh_588 = buffer.data(lsh + 588);
    const auto *lsh_589 = buffer.data(lsh + 589);
    const auto *lsh_590 = buffer.data(lsh + 590);
    const auto *lsh_591 = buffer.data(lsh + 591);
    const auto *lsh_593 = buffer.data(lsh + 593);
    const auto *lsh_594 = buffer.data(lsh + 594);
    const auto *lsh_595 = buffer.data(lsh + 595);
    const auto *lsh_597 = buffer.data(lsh + 597);
    const auto *lsh_598 = buffer.data(lsh + 598);
    const auto *lsh_603 = buffer.data(lsh + 603);
    const auto *lsh_605 = buffer.data(lsh + 605);
    const auto *lsh_606 = buffer.data(lsh + 606);
    const auto *lsh_607 = buffer.data(lsh + 607);
    const auto *lsh_608 = buffer.data(lsh + 608);
    const auto *lsh_609 = buffer.data(lsh + 609);
    const auto *lsh_611 = buffer.data(lsh + 611);
    const auto *lsh_612 = buffer.data(lsh + 612);
    const auto *lsh_614 = buffer.data(lsh + 614);

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, ksh_401, ksh_528, ksh_530, lsg0_378, \
                         lsg0_380, lsg1_378, lsg1_380, lsh_527, lsh_528, \
                         lsh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_12 * ksh_528[k]
                   + f_8 * lsg0_378[k]
                   - f_9 * lsg1_378[k]
                   + f_3 * pc_x[k] * lsh_528[k];

        t_704[k] = f_12 * ksh_401[k]
                   + f_3 * pc_y[k] * lsh_527[k];

        t_705[k] = f_12 * ksh_530[k]
                   + f_8 * lsg0_380[k]
                   - f_9 * lsg1_380[k]
                   + f_3 * pc_x[k] * lsh_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, pc_y, pc_z, ksh_381, ksh_404, ksh_531, \
                         lsg0_381, lsg1_381, lsh_528, lsh_530, \
                         lsh_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_12 * ksh_531[k]
                   + f_6 * lsg0_381[k]
                   - f_7 * lsg1_381[k]
                   + f_3 * pc_x[k] * lsh_531[k];

        t_707[k] = f_14 * ksh_381[k]
                   + f_3 * pc_z[k] * lsh_528[k];

        t_708[k] = f_12 * ksh_404[k]
                   + f_3 * pc_y[k] * lsh_530[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, ksh_384, ksh_534, ksh_535, lsg0_384, \
                         lsg0_385, lsg1_384, lsg1_385, lsh_531, lsh_534, \
                         lsh_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_12 * ksh_534[k]
                   + f_6 * lsg0_384[k]
                   - f_7 * lsg1_384[k]
                   + f_3 * pc_x[k] * lsh_534[k];

        t_710[k] = f_12 * ksh_535[k]
                   + f_4 * lsg0_385[k]
                   - f_5 * lsg1_385[k]
                   + f_3 * pc_x[k] * lsh_535[k];

        t_711[k] = f_14 * ksh_384[k]
                   + f_3 * pc_z[k] * lsh_531[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pc_x, pc_y, ksh_408, ksh_537, ksh_539, lsg0_387, \
                         lsg0_389, lsg1_387, lsg1_389, lsh_534, lsh_537, \
                         lsh_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_12 * ksh_537[k]
                   + f_4 * lsg0_387[k]
                   - f_5 * lsg1_387[k]
                   + f_3 * pc_x[k] * lsh_537[k];

        t_713[k] = f_12 * ksh_408[k]
                   + f_3 * pc_y[k] * lsh_534[k];

        t_714[k] = f_12 * ksh_539[k]
                   + f_4 * lsg0_389[k]
                   - f_5 * lsg1_389[k]
                   + f_3 * pc_x[k] * lsh_539[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, t_718, t_719, pc_x, ksh_540, ksh_541, ksh_542, \
                         ksh_543, ksh_544, lsh_540, lsh_541, lsh_542, lsh_543, \
                         lsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = f_12 * ksh_540[k]
                   + f_3 * pc_x[k] * lsh_540[k];

        t_716[k] = f_12 * ksh_541[k]
                   + f_3 * pc_x[k] * lsh_541[k];

        t_717[k] = f_12 * ksh_542[k]
                   + f_3 * pc_x[k] * lsh_542[k];

        t_718[k] = f_12 * ksh_543[k]
                   + f_3 * pc_x[k] * lsh_543[k];

        t_719[k] = f_12 * ksh_544[k]
                   + f_3 * pc_x[k] * lsh_544[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, ksh_393, ksh_414, ksh_545, \
                         lsg0_385, lsg1_385, lsh_540, lsh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_12 * ksh_545[k]
                   + f_3 * pc_x[k] * lsh_545[k];

        t_721[k] = f_12 * ksh_414[k]
                   + f_1 * lsg0_385[k]
                   - f_2 * lsg1_385[k]
                   + f_3 * pc_y[k] * lsh_540[k];

        t_722[k] = f_14 * ksh_393[k]
                   + f_3 * pc_z[k] * lsh_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, ksh_416, ksh_417, ksh_418, lsg0_387, \
                         lsg0_388, lsg0_389, lsg1_387, lsg1_388, lsg1_389, lsh_542, lsh_543, \
                         lsh_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * ksh_416[k]
                   + f_8 * lsg0_387[k]
                   - f_9 * lsg1_387[k]
                   + f_3 * pc_y[k] * lsh_542[k];

        t_724[k] = f_12 * ksh_417[k]
                   + f_6 * lsg0_388[k]
                   - f_7 * lsg1_388[k]
                   + f_3 * pc_y[k] * lsh_543[k];

        t_725[k] = f_12 * ksh_418[k]
                   + f_4 * lsg0_389[k]
                   - f_5 * lsg1_389[k]
                   + f_3 * pc_y[k] * lsh_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_y, pc_y, pc_z, ksi0_560, ksh_398, \
                         ksh_419, ksh_420, ksi1_560, lsg0_389, lsg1_389, lsh_545, \
                         lsh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * ksh_419[k]
                   + f_3 * pc_y[k] * lsh_545[k];

        t_727[k] = f_14 * ksh_398[k]
                   + f_1 * lsg0_389[k]
                   - f_2 * lsg1_389[k]
                   + f_3 * pc_z[k] * lsh_545[k];

        t_728[k] = pa_y[k] * ksi0_560[k]
                   - f_10 * pc_y[k] * ksi1_560[k];

        t_729[k] = f_11 * ksh_420[k]
                   + f_3 * pc_y[k] * lsh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, t_733, pa_y, pc_y, pc_z, ksi0_563, ksi0_565, \
                         ksh_399, ksh_421, ksh_422, ksi1_563, ksi1_565, lsh_546, \
                         lsh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_19 * ksh_399[k]
                   + f_3 * pc_z[k] * lsh_546[k];

        t_731[k] = pa_y[k] * ksi0_563[k]
                   + f_12 * ksh_421[k]
                   - f_10 * pc_y[k] * ksi1_563[k];

        t_732[k] = f_11 * ksh_422[k]
                   + f_3 * pc_y[k] * lsh_548[k];

        t_733[k] = pa_y[k] * ksi0_565[k]
                   - f_10 * pc_y[k] * ksi1_565[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, t_737, pa_y, pc_y, pc_z, ksi0_566, ksi0_569, \
                         ksh_402, ksh_423, ksh_425, ksi1_566, ksi1_569, lsh_549, \
                         lsh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pa_y[k] * ksi0_566[k]
                   + f_13 * ksh_423[k]
                   - f_10 * pc_y[k] * ksi1_566[k];

        t_735[k] = f_19 * ksh_402[k]
                   + f_3 * pc_z[k] * lsh_549[k];

        t_736[k] = f_11 * ksh_425[k]
                   + f_3 * pc_y[k] * lsh_551[k];

        t_737[k] = pa_y[k] * ksi0_569[k]
                   - f_10 * pc_y[k] * ksi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pa_y, pc_y, pc_z, ksi0_570, ksi0_572, ksh_405, \
                         ksh_426, ksh_428, ksi1_570, ksi1_572, \
                         lsh_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = pa_y[k] * ksi0_570[k]
                   + f_14 * ksh_426[k]
                   - f_10 * pc_y[k] * ksi1_570[k];

        t_739[k] = f_19 * ksh_405[k]
                   + f_3 * pc_z[k] * lsh_552[k];

        t_740[k] = pa_y[k] * ksi0_572[k]
                   + f_12 * ksh_428[k]
                   - f_10 * pc_y[k] * ksi1_572[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, pa_y, pc_x, pc_y, ksi0_574, ksh_429, \
                         ksh_561, ksh_562, ksi1_574, lsh_555, lsh_561, \
                         lsh_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_11 * ksh_429[k]
                   + f_3 * pc_y[k] * lsh_555[k];

        t_742[k] = pa_y[k] * ksi0_574[k]
                   - f_10 * pc_y[k] * ksi1_574[k];

        t_743[k] = f_12 * ksh_561[k]
                   + f_3 * pc_x[k] * lsh_561[k];

        t_744[k] = f_12 * ksh_562[k]
                   + f_3 * pc_x[k] * lsh_562[k];
    }

#pragma omp simd aligned(t_745, t_746, t_747, t_748, pc_x, ksh_563, ksh_564, ksh_565, ksh_566, \
                         lsh_563, lsh_564, lsh_565, lsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_745[k] = f_12 * ksh_563[k]
                   + f_3 * pc_x[k] * lsh_563[k];

        t_746[k] = f_12 * ksh_564[k]
                   + f_3 * pc_x[k] * lsh_564[k];

        t_747[k] = f_12 * ksh_565[k]
                   + f_3 * pc_x[k] * lsh_565[k];

        t_748[k] = f_12 * ksh_566[k]
                   + f_3 * pc_x[k] * lsh_566[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_y, pc_z, ksh_414, ksh_435, ksh_437, lsg0_400, \
                         lsg0_402, lsg1_400, lsg1_402, lsh_561, \
                         lsh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * ksh_435[k]
                   + f_1 * lsg0_400[k]
                   - f_2 * lsg1_400[k]
                   + f_3 * pc_y[k] * lsh_561[k];

        t_750[k] = f_19 * ksh_414[k]
                   + f_3 * pc_z[k] * lsh_561[k];

        t_751[k] = f_11 * ksh_437[k]
                   + f_8 * lsg0_402[k]
                   - f_9 * lsg1_402[k]
                   + f_3 * pc_y[k] * lsh_563[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, pc_y, ksh_438, ksh_439, ksh_440, lsg0_403, \
                         lsg0_404, lsg1_403, lsg1_404, lsh_564, lsh_565, \
                         lsh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_11 * ksh_438[k]
                   + f_6 * lsg0_403[k]
                   - f_7 * lsg1_403[k]
                   + f_3 * pc_y[k] * lsh_564[k];

        t_753[k] = f_11 * ksh_439[k]
                   + f_4 * lsg0_404[k]
                   - f_5 * lsg1_404[k]
                   + f_3 * pc_y[k] * lsh_565[k];

        t_754[k] = f_11 * ksh_440[k]
                   + f_3 * pc_y[k] * lsh_566[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pa_y, pc_x, pc_y, pc_z, ksi0_587, \
                         ksh_420, ksh_567, ksi1_587, lsg0_405, lsg1_405, \
                         lsh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = pa_y[k] * ksi0_587[k]
                   - f_10 * pc_y[k] * ksi1_587[k];

        t_756[k] = f_12 * ksh_567[k]
                   + f_1 * lsg0_405[k]
                   - f_2 * lsg1_405[k]
                   + f_3 * pc_x[k] * lsh_567[k];

        t_757[k] = f_3 * pc_y[k] * lsh_567[k];

        t_758[k] = f_18 * ksh_420[k]
                   + f_3 * pc_z[k] * lsh_567[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pc_x, pc_y, ksh_572, lsg0_405, lsg0_410, \
                         lsg1_405, lsg1_410, lsh_568, lsh_569, \
                         lsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = f_4 * lsg0_405[k]
                   - f_5 * lsg1_405[k]
                   + f_3 * pc_y[k] * lsh_568[k];

        t_760[k] = f_3 * pc_y[k] * lsh_569[k];

        t_761[k] = f_12 * ksh_572[k]
                   + f_8 * lsg0_410[k]
                   - f_9 * lsg1_410[k]
                   + f_3 * pc_x[k] * lsh_572[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pc_y, lsg0_406, lsg0_407, lsg1_406, lsg1_407, \
                         lsh_570, lsh_571, lsh_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = f_6 * lsg0_406[k]
                   - f_7 * lsg1_406[k]
                   + f_3 * pc_y[k] * lsh_570[k];

        t_763[k] = f_4 * lsg0_407[k]
                   - f_5 * lsg1_407[k]
                   + f_3 * pc_y[k] * lsh_571[k];

        t_764[k] = f_3 * pc_y[k] * lsh_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, pc_y, ksh_576, lsg0_408, lsg0_409, \
                         lsg0_414, lsg1_408, lsg1_409, lsg1_414, lsh_573, lsh_574, \
                         lsh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_12 * ksh_576[k]
                   + f_6 * lsg0_414[k]
                   - f_7 * lsg1_414[k]
                   + f_3 * pc_x[k] * lsh_576[k];

        t_766[k] = f_8 * lsg0_408[k]
                   - f_9 * lsg1_408[k]
                   + f_3 * pc_y[k] * lsh_573[k];

        t_767[k] = f_6 * lsg0_409[k]
                   - f_7 * lsg1_409[k]
                   + f_3 * pc_y[k] * lsh_574[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pc_x, pc_y, ksh_581, ksh_582, lsg0_410, \
                         lsg0_419, lsg1_410, lsg1_419, lsh_575, lsh_576, lsh_581, \
                         lsh_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * lsg0_410[k]
                   - f_5 * lsg1_410[k]
                   + f_3 * pc_y[k] * lsh_575[k];

        t_769[k] = f_3 * pc_y[k] * lsh_576[k];

        t_770[k] = f_12 * ksh_581[k]
                   + f_4 * lsg0_419[k]
                   - f_5 * lsg1_419[k]
                   + f_3 * pc_x[k] * lsh_581[k];

        t_771[k] = f_12 * ksh_582[k]
                   + f_3 * pc_x[k] * lsh_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, pc_y, ksh_583, ksh_584, \
                         ksh_585, ksh_587, lsh_581, lsh_583, lsh_584, lsh_585, \
                         lsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_12 * ksh_583[k]
                   + f_3 * pc_x[k] * lsh_583[k];

        t_773[k] = f_12 * ksh_584[k]
                   + f_3 * pc_x[k] * lsh_584[k];

        t_774[k] = f_12 * ksh_585[k]
                   + f_3 * pc_x[k] * lsh_585[k];

        t_775[k] = f_3 * pc_y[k] * lsh_581[k];

        t_776[k] = f_12 * ksh_587[k]
                   + f_3 * pc_x[k] * lsh_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, lsg0_415, lsg0_416, lsg0_417, lsg1_415, \
                         lsg1_416, lsg1_417, lsh_582, lsh_583, \
                         lsh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_1 * lsg0_415[k]
                   - f_2 * lsg1_415[k]
                   + f_3 * pc_y[k] * lsh_582[k];

        t_778[k] = f_16 * lsg0_416[k]
                   - f_17 * lsg1_416[k]
                   + f_3 * pc_y[k] * lsh_583[k];

        t_779[k] = f_8 * lsg0_417[k]
                   - f_9 * lsg1_417[k]
                   + f_3 * pc_y[k] * lsh_584[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, t_783, pc_y, pc_z, ksh_440, lsg0_418, lsg0_419, \
                         lsg1_418, lsg1_419, lsh_585, lsh_586, \
                         lsh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_6 * lsg0_418[k]
                   - f_7 * lsg1_418[k]
                   + f_3 * pc_y[k] * lsh_585[k];

        t_781[k] = f_4 * lsg0_419[k]
                   - f_5 * lsg1_419[k]
                   + f_3 * pc_y[k] * lsh_586[k];

        t_782[k] = f_3 * pc_y[k] * lsh_587[k];

        t_783[k] = f_18 * ksh_440[k]
                   + f_1 * lsg0_419[k]
                   - f_2 * lsg1_419[k]
                   + f_3 * pc_z[k] * lsh_587[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pa_x, pc_x, pc_y, pc_z, ksi0_784, \
                         ksi0_787, ksh_441, ksh_588, ksh_591, ksi1_784, ksi1_787, \
                         lsh_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = pa_x[k] * ksi0_784[k]
                   + f_18 * ksh_588[k]
                   - f_10 * pc_x[k] * ksi1_784[k];

        t_785[k] = f_15 * ksh_441[k]
                   + f_3 * pc_y[k] * lsh_588[k];

        t_786[k] = f_3 * pc_z[k] * lsh_588[k];

        t_787[k] = pa_x[k] * ksi0_787[k]
                   + f_14 * ksh_591[k]
                   - f_10 * pc_x[k] * ksi1_787[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, pa_x, pc_x, pc_z, ksi0_790, ksh_594, \
                         ksi1_790, lsg0_420, lsg1_420, lsh_589, lsh_590, \
                         lsh_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_3 * pc_z[k] * lsh_589[k];

        t_789[k] = f_4 * lsg0_420[k]
                   - f_5 * lsg1_420[k]
                   + f_3 * pc_z[k] * lsh_590[k];

        t_790[k] = pa_x[k] * ksi0_790[k]
                   + f_13 * ksh_594[k]
                   - f_10 * pc_x[k] * ksi1_790[k];

        t_791[k] = f_3 * pc_z[k] * lsh_591[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, t_795, pa_x, pc_x, pc_y, pc_z, ksi0_794, \
                         ksh_446, ksh_598, ksi1_794, lsg0_422, lsg1_422, lsh_593, \
                         lsh_594 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_15 * ksh_446[k]
                   + f_3 * pc_y[k] * lsh_593[k];

        t_793[k] = f_6 * lsg0_422[k]
                   - f_7 * lsg1_422[k]
                   + f_3 * pc_z[k] * lsh_593[k];

        t_794[k] = pa_x[k] * ksi0_794[k]
                   + f_12 * ksh_598[k]
                   - f_10 * pc_x[k] * ksi1_794[k];

        t_795[k] = f_3 * pc_z[k] * lsh_594[k];
    }

#pragma omp simd aligned(t_796, t_797, t_798, t_799, pc_x, pc_y, pc_z, ksh_450, ksh_603, \
                         lsg0_423, lsg0_425, lsg1_423, lsg1_425, lsh_595, lsh_597, \
                         lsh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_796[k] = f_4 * lsg0_423[k]
                   - f_5 * lsg1_423[k]
                   + f_3 * pc_z[k] * lsh_595[k];

        t_797[k] = f_15 * ksh_450[k]
                   + f_3 * pc_y[k] * lsh_597[k];

        t_798[k] = f_8 * lsg0_425[k]
                   - f_9 * lsg1_425[k]
                   + f_3 * pc_z[k] * lsh_597[k];

        t_799[k] = f_11 * ksh_603[k]
                   + f_3 * pc_x[k] * lsh_603[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, t_804, pc_x, pc_z, ksh_605, ksh_606, \
                         ksh_607, ksh_608, lsh_598, lsh_605, lsh_606, lsh_607, \
                         lsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_3 * pc_z[k] * lsh_598[k];

        t_801[k] = f_11 * ksh_605[k]
                   + f_3 * pc_x[k] * lsh_605[k];

        t_802[k] = f_11 * ksh_606[k]
                   + f_3 * pc_x[k] * lsh_606[k];

        t_803[k] = f_11 * ksh_607[k]
                   + f_3 * pc_x[k] * lsh_607[k];

        t_804[k] = f_11 * ksh_608[k]
                   + f_3 * pc_x[k] * lsh_608[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, t_808, pa_x, pc_x, pc_z, ksi0_805, ksi0_807, \
                         ksi0_808, ksi1_805, ksi1_807, ksi1_808, \
                         lsh_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = pa_x[k] * ksi0_805[k]
                   - f_10 * pc_x[k] * ksi1_805[k];

        t_806[k] = f_3 * pc_z[k] * lsh_603[k];

        t_807[k] = pa_x[k] * ksi0_807[k]
                   - f_10 * pc_x[k] * ksi1_807[k];

        t_808[k] = pa_x[k] * ksi0_808[k]
                   - f_10 * pc_x[k] * ksi1_808[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pa_x, pc_x, pc_y, ksi0_809, ksi0_811, ksh_461, \
                         ksi1_809, ksi1_811, lsh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = pa_x[k] * ksi0_809[k]
                   - f_10 * pc_x[k] * ksi1_809[k];

        t_810[k] = f_15 * ksh_461[k]
                   + f_3 * pc_y[k] * lsh_608[k];

        t_811[k] = pa_x[k] * ksi0_811[k]
                   - f_10 * pc_x[k] * ksi1_811[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pa_z, pc_y, pc_z, ksi0_588, ksi0_591, \
                         ksh_441, ksh_462, ksi1_588, ksi1_591, \
                         lsh_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pa_z[k] * ksi0_588[k]
                   - f_10 * pc_z[k] * ksi1_588[k];

        t_813[k] = f_18 * ksh_462[k]
                   + f_3 * pc_y[k] * lsh_609[k];

        t_814[k] = f_11 * ksh_441[k]
                   + f_3 * pc_z[k] * lsh_609[k];

        t_815[k] = pa_z[k] * ksi0_591[k]
                   - f_10 * pc_z[k] * ksi1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pa_x, pa_z, pc_x, pc_y, pc_z, ksi0_594, \
                         ksi0_817, ksh_464, ksh_614, ksi1_594, ksi1_817, \
                         lsh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_18 * ksh_464[k]
                   + f_3 * pc_y[k] * lsh_611[k];

        t_817[k] = pa_x[k] * ksi0_817[k]
                   + f_14 * ksh_614[k]
                   - f_10 * pc_x[k] * ksi1_817[k];

        t_818[k] = pa_z[k] * ksi0_594[k]
                   - f_10 * pc_z[k] * ksi1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pa_x, pc_x, pc_y, pc_z, ksi0_821, ksh_444, \
                         ksh_467, ksh_618, ksi1_821, lsh_612, lsh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_11 * ksh_444[k]
                   + f_3 * pc_z[k] * lsh_612[k];

        t_820[k] = f_18 * ksh_467[k]
                   + f_3 * pc_y[k] * lsh_614[k];

        t_821[k] = pa_x[k] * ksi0_821[k]
                   + f_13 * ksh_618[k]
                   - f_10 * pc_x[k] * ksi1_821[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_598 = buffer.data(ksi0 + 598);
    const auto *ksi0_824 = buffer.data(ksi0 + 824);
    const auto *ksi0_826 = buffer.data(ksi0 + 826);
    const auto *ksi0_833 = buffer.data(ksi0 + 833);
    const auto *ksi0_835 = buffer.data(ksi0 + 835);
    const auto *ksi0_836 = buffer.data(ksi0 + 836);
    const auto *ksi0_837 = buffer.data(ksi0 + 837);
    const auto *ksi0_839 = buffer.data(ksi0 + 839);
    const auto *ksi0_840 = buffer.data(ksi0 + 840);
    const auto *ksi0_843 = buffer.data(ksi0 + 843);
    const auto *ksi0_845 = buffer.data(ksi0 + 845);
    const auto *ksi0_846 = buffer.data(ksi0 + 846);
    const auto *ksi0_849 = buffer.data(ksi0 + 849);
    const auto *ksi0_850 = buffer.data(ksi0 + 850);
    const auto *ksi0_852 = buffer.data(ksi0 + 852);
    const auto *ksi0_854 = buffer.data(ksi0 + 854);
    const auto *ksi0_861 = buffer.data(ksi0 + 861);
    const auto *ksi0_863 = buffer.data(ksi0 + 863);
    const auto *ksi0_864 = buffer.data(ksi0 + 864);
    const auto *ksi0_865 = buffer.data(ksi0 + 865);
    const auto *ksi0_867 = buffer.data(ksi0 + 867);
    const auto *ksi0_868 = buffer.data(ksi0 + 868);
    const auto *ksi0_871 = buffer.data(ksi0 + 871);
    const auto *ksi0_873 = buffer.data(ksi0 + 873);
    const auto *ksi0_874 = buffer.data(ksi0 + 874);
    const auto *ksi0_877 = buffer.data(ksi0 + 877);
    const auto *ksi0_878 = buffer.data(ksi0 + 878);
    const auto *ksi0_880 = buffer.data(ksi0 + 880);
    const auto *ksi0_882 = buffer.data(ksi0 + 882);
    const auto *ksi0_889 = buffer.data(ksi0 + 889);
    const auto *ksi0_891 = buffer.data(ksi0 + 891);
    const auto *ksi0_892 = buffer.data(ksi0 + 892);
    const auto *ksi0_893 = buffer.data(ksi0 + 893);
    const auto *ksi0_895 = buffer.data(ksi0 + 895);
    const auto *ksi0_896 = buffer.data(ksi0 + 896);
    const auto *ksi0_899 = buffer.data(ksi0 + 899);
    const auto *ksi0_901 = buffer.data(ksi0 + 901);
    const auto *ksi0_902 = buffer.data(ksi0 + 902);
    const auto *ksi0_905 = buffer.data(ksi0 + 905);
    const auto *ksi0_906 = buffer.data(ksi0 + 906);
    const auto *ksi0_908 = buffer.data(ksi0 + 908);
    const auto *ksi0_910 = buffer.data(ksi0 + 910);
    const auto *ksi0_917 = buffer.data(ksi0 + 917);
    const auto *ksi0_919 = buffer.data(ksi0 + 919);
    const auto *ksi0_920 = buffer.data(ksi0 + 920);
    const auto *ksi0_921 = buffer.data(ksi0 + 921);
    const auto *ksi0_923 = buffer.data(ksi0 + 923);
    const auto *ksi0_924 = buffer.data(ksi0 + 924);
    const auto *ksi0_927 = buffer.data(ksi0 + 927);
    const auto *ksi0_929 = buffer.data(ksi0 + 929);
    const auto *ksi0_930 = buffer.data(ksi0 + 930);
    const auto *ksi0_933 = buffer.data(ksi0 + 933);
    const auto *ksi0_934 = buffer.data(ksi0 + 934);
    const auto *ksi0_936 = buffer.data(ksi0 + 936);
    const auto *ksi0_938 = buffer.data(ksi0 + 938);

    const auto *ksh_447 = buffer.data(ksh + 447);
    const auto *ksh_456 = buffer.data(ksh + 456);
    const auto *ksh_462 = buffer.data(ksh + 462);
    const auto *ksh_465 = buffer.data(ksh + 465);
    const auto *ksh_468 = buffer.data(ksh + 468);
    const auto *ksh_471 = buffer.data(ksh + 471);
    const auto *ksh_477 = buffer.data(ksh + 477);
    const auto *ksh_482 = buffer.data(ksh + 482);
    const auto *ksh_483 = buffer.data(ksh + 483);
    const auto *ksh_485 = buffer.data(ksh + 485);
    const auto *ksh_486 = buffer.data(ksh + 486);
    const auto *ksh_488 = buffer.data(ksh + 488);
    const auto *ksh_489 = buffer.data(ksh + 489);
    const auto *ksh_492 = buffer.data(ksh + 492);
    const auto *ksh_498 = buffer.data(ksh + 498);
    const auto *ksh_503 = buffer.data(ksh + 503);
    const auto *ksh_504 = buffer.data(ksh + 504);
    const auto *ksh_506 = buffer.data(ksh + 506);
    const auto *ksh_507 = buffer.data(ksh + 507);
    const auto *ksh_509 = buffer.data(ksh + 509);
    const auto *ksh_510 = buffer.data(ksh + 510);
    const auto *ksh_513 = buffer.data(ksh + 513);
    const auto *ksh_519 = buffer.data(ksh + 519);
    const auto *ksh_524 = buffer.data(ksh + 524);
    const auto *ksh_525 = buffer.data(ksh + 525);
    const auto *ksh_527 = buffer.data(ksh + 527);
    const auto *ksh_528 = buffer.data(ksh + 528);
    const auto *ksh_530 = buffer.data(ksh + 530);
    const auto *ksh_531 = buffer.data(ksh + 531);
    const auto *ksh_534 = buffer.data(ksh + 534);
    const auto *ksh_545 = buffer.data(ksh + 545);
    const auto *ksh_546 = buffer.data(ksh + 546);
    const auto *ksh_548 = buffer.data(ksh + 548);
    const auto *ksh_551 = buffer.data(ksh + 551);
    const auto *ksh_555 = buffer.data(ksh + 555);
    const auto *ksh_621 = buffer.data(ksh + 621);
    const auto *ksh_623 = buffer.data(ksh + 623);
    const auto *ksh_624 = buffer.data(ksh + 624);
    const auto *ksh_625 = buffer.data(ksh + 625);
    const auto *ksh_626 = buffer.data(ksh + 626);
    const auto *ksh_627 = buffer.data(ksh + 627);
    const auto *ksh_628 = buffer.data(ksh + 628);
    const auto *ksh_629 = buffer.data(ksh + 629);
    const auto *ksh_630 = buffer.data(ksh + 630);
    const auto *ksh_633 = buffer.data(ksh + 633);
    const auto *ksh_635 = buffer.data(ksh + 635);
    const auto *ksh_636 = buffer.data(ksh + 636);
    const auto *ksh_639 = buffer.data(ksh + 639);
    const auto *ksh_640 = buffer.data(ksh + 640);
    const auto *ksh_642 = buffer.data(ksh + 642);
    const auto *ksh_644 = buffer.data(ksh + 644);
    const auto *ksh_645 = buffer.data(ksh + 645);
    const auto *ksh_646 = buffer.data(ksh + 646);
    const auto *ksh_647 = buffer.data(ksh + 647);
    const auto *ksh_648 = buffer.data(ksh + 648);
    const auto *ksh_649 = buffer.data(ksh + 649);
    const auto *ksh_650 = buffer.data(ksh + 650);
    const auto *ksh_651 = buffer.data(ksh + 651);
    const auto *ksh_654 = buffer.data(ksh + 654);
    const auto *ksh_656 = buffer.data(ksh + 656);
    const auto *ksh_657 = buffer.data(ksh + 657);
    const auto *ksh_660 = buffer.data(ksh + 660);
    const auto *ksh_661 = buffer.data(ksh + 661);
    const auto *ksh_663 = buffer.data(ksh + 663);
    const auto *ksh_665 = buffer.data(ksh + 665);
    const auto *ksh_666 = buffer.data(ksh + 666);
    const auto *ksh_667 = buffer.data(ksh + 667);
    const auto *ksh_668 = buffer.data(ksh + 668);
    const auto *ksh_669 = buffer.data(ksh + 669);
    const auto *ksh_670 = buffer.data(ksh + 670);
    const auto *ksh_671 = buffer.data(ksh + 671);
    const auto *ksh_672 = buffer.data(ksh + 672);
    const auto *ksh_675 = buffer.data(ksh + 675);
    const auto *ksh_677 = buffer.data(ksh + 677);
    const auto *ksh_678 = buffer.data(ksh + 678);
    const auto *ksh_681 = buffer.data(ksh + 681);
    const auto *ksh_682 = buffer.data(ksh + 682);
    const auto *ksh_684 = buffer.data(ksh + 684);
    const auto *ksh_686 = buffer.data(ksh + 686);
    const auto *ksh_687 = buffer.data(ksh + 687);
    const auto *ksh_688 = buffer.data(ksh + 688);
    const auto *ksh_689 = buffer.data(ksh + 689);
    const auto *ksh_690 = buffer.data(ksh + 690);
    const auto *ksh_691 = buffer.data(ksh + 691);
    const auto *ksh_692 = buffer.data(ksh + 692);
    const auto *ksh_693 = buffer.data(ksh + 693);
    const auto *ksh_696 = buffer.data(ksh + 696);
    const auto *ksh_698 = buffer.data(ksh + 698);
    const auto *ksh_699 = buffer.data(ksh + 699);
    const auto *ksh_702 = buffer.data(ksh + 702);
    const auto *ksh_703 = buffer.data(ksh + 703);
    const auto *ksh_705 = buffer.data(ksh + 705);
    const auto *ksh_707 = buffer.data(ksh + 707);
    const auto *ksh_708 = buffer.data(ksh + 708);
    const auto *ksh_709 = buffer.data(ksh + 709);
    const auto *ksh_710 = buffer.data(ksh + 710);

    const auto *ksi1_598 = buffer.data(ksi1 + 598);
    const auto *ksi1_824 = buffer.data(ksi1 + 824);
    const auto *ksi1_826 = buffer.data(ksi1 + 826);
    const auto *ksi1_833 = buffer.data(ksi1 + 833);
    const auto *ksi1_835 = buffer.data(ksi1 + 835);
    const auto *ksi1_836 = buffer.data(ksi1 + 836);
    const auto *ksi1_837 = buffer.data(ksi1 + 837);
    const auto *ksi1_839 = buffer.data(ksi1 + 839);
    const auto *ksi1_840 = buffer.data(ksi1 + 840);
    const auto *ksi1_843 = buffer.data(ksi1 + 843);
    const auto *ksi1_845 = buffer.data(ksi1 + 845);
    const auto *ksi1_846 = buffer.data(ksi1 + 846);
    const auto *ksi1_849 = buffer.data(ksi1 + 849);
    const auto *ksi1_850 = buffer.data(ksi1 + 850);
    const auto *ksi1_852 = buffer.data(ksi1 + 852);
    const auto *ksi1_854 = buffer.data(ksi1 + 854);
    const auto *ksi1_861 = buffer.data(ksi1 + 861);
    const auto *ksi1_863 = buffer.data(ksi1 + 863);
    const auto *ksi1_864 = buffer.data(ksi1 + 864);
    const auto *ksi1_865 = buffer.data(ksi1 + 865);
    const auto *ksi1_867 = buffer.data(ksi1 + 867);
    const auto *ksi1_868 = buffer.data(ksi1 + 868);
    const auto *ksi1_871 = buffer.data(ksi1 + 871);
    const auto *ksi1_873 = buffer.data(ksi1 + 873);
    const auto *ksi1_874 = buffer.data(ksi1 + 874);
    const auto *ksi1_877 = buffer.data(ksi1 + 877);
    const auto *ksi1_878 = buffer.data(ksi1 + 878);
    const auto *ksi1_880 = buffer.data(ksi1 + 880);
    const auto *ksi1_882 = buffer.data(ksi1 + 882);
    const auto *ksi1_889 = buffer.data(ksi1 + 889);
    const auto *ksi1_891 = buffer.data(ksi1 + 891);
    const auto *ksi1_892 = buffer.data(ksi1 + 892);
    const auto *ksi1_893 = buffer.data(ksi1 + 893);
    const auto *ksi1_895 = buffer.data(ksi1 + 895);
    const auto *ksi1_896 = buffer.data(ksi1 + 896);
    const auto *ksi1_899 = buffer.data(ksi1 + 899);
    const auto *ksi1_901 = buffer.data(ksi1 + 901);
    const auto *ksi1_902 = buffer.data(ksi1 + 902);
    const auto *ksi1_905 = buffer.data(ksi1 + 905);
    const auto *ksi1_906 = buffer.data(ksi1 + 906);
    const auto *ksi1_908 = buffer.data(ksi1 + 908);
    const auto *ksi1_910 = buffer.data(ksi1 + 910);
    const auto *ksi1_917 = buffer.data(ksi1 + 917);
    const auto *ksi1_919 = buffer.data(ksi1 + 919);
    const auto *ksi1_920 = buffer.data(ksi1 + 920);
    const auto *ksi1_921 = buffer.data(ksi1 + 921);
    const auto *ksi1_923 = buffer.data(ksi1 + 923);
    const auto *ksi1_924 = buffer.data(ksi1 + 924);
    const auto *ksi1_927 = buffer.data(ksi1 + 927);
    const auto *ksi1_929 = buffer.data(ksi1 + 929);
    const auto *ksi1_930 = buffer.data(ksi1 + 930);
    const auto *ksi1_933 = buffer.data(ksi1 + 933);
    const auto *ksi1_934 = buffer.data(ksi1 + 934);
    const auto *ksi1_936 = buffer.data(ksi1 + 936);
    const auto *ksi1_938 = buffer.data(ksi1 + 938);

    const auto *lsh_615 = buffer.data(lsh + 615);
    const auto *lsh_618 = buffer.data(lsh + 618);
    const auto *lsh_624 = buffer.data(lsh + 624);
    const auto *lsh_625 = buffer.data(lsh + 625);
    const auto *lsh_626 = buffer.data(lsh + 626);
    const auto *lsh_627 = buffer.data(lsh + 627);
    const auto *lsh_628 = buffer.data(lsh + 628);
    const auto *lsh_629 = buffer.data(lsh + 629);
    const auto *lsh_630 = buffer.data(lsh + 630);
    const auto *lsh_632 = buffer.data(lsh + 632);
    const auto *lsh_633 = buffer.data(lsh + 633);
    const auto *lsh_635 = buffer.data(lsh + 635);
    const auto *lsh_636 = buffer.data(lsh + 636);
    const auto *lsh_639 = buffer.data(lsh + 639);
    const auto *lsh_645 = buffer.data(lsh + 645);
    const auto *lsh_646 = buffer.data(lsh + 646);
    const auto *lsh_647 = buffer.data(lsh + 647);
    const auto *lsh_648 = buffer.data(lsh + 648);
    const auto *lsh_649 = buffer.data(lsh + 649);
    const auto *lsh_650 = buffer.data(lsh + 650);
    const auto *lsh_651 = buffer.data(lsh + 651);
    const auto *lsh_653 = buffer.data(lsh + 653);
    const auto *lsh_654 = buffer.data(lsh + 654);
    const auto *lsh_656 = buffer.data(lsh + 656);
    const auto *lsh_657 = buffer.data(lsh + 657);
    const auto *lsh_660 = buffer.data(lsh + 660);
    const auto *lsh_666 = buffer.data(lsh + 666);
    const auto *lsh_667 = buffer.data(lsh + 667);
    const auto *lsh_668 = buffer.data(lsh + 668);
    const auto *lsh_669 = buffer.data(lsh + 669);
    const auto *lsh_670 = buffer.data(lsh + 670);
    const auto *lsh_671 = buffer.data(lsh + 671);
    const auto *lsh_672 = buffer.data(lsh + 672);
    const auto *lsh_674 = buffer.data(lsh + 674);
    const auto *lsh_675 = buffer.data(lsh + 675);
    const auto *lsh_677 = buffer.data(lsh + 677);
    const auto *lsh_678 = buffer.data(lsh + 678);
    const auto *lsh_681 = buffer.data(lsh + 681);
    const auto *lsh_687 = buffer.data(lsh + 687);
    const auto *lsh_688 = buffer.data(lsh + 688);
    const auto *lsh_689 = buffer.data(lsh + 689);
    const auto *lsh_690 = buffer.data(lsh + 690);
    const auto *lsh_691 = buffer.data(lsh + 691);
    const auto *lsh_692 = buffer.data(lsh + 692);
    const auto *lsh_693 = buffer.data(lsh + 693);
    const auto *lsh_695 = buffer.data(lsh + 695);
    const auto *lsh_696 = buffer.data(lsh + 696);
    const auto *lsh_698 = buffer.data(lsh + 698);
    const auto *lsh_699 = buffer.data(lsh + 699);
    const auto *lsh_702 = buffer.data(lsh + 702);
    const auto *lsh_708 = buffer.data(lsh + 708);
    const auto *lsh_709 = buffer.data(lsh + 709);
    const auto *lsh_710 = buffer.data(lsh + 710);

#pragma omp simd aligned(t_822, t_823, t_824, pa_x, pa_z, pc_x, pc_z, ksi0_598, ksi0_824, \
                         ksh_447, ksh_621, ksi1_598, ksi1_824, \
                         lsh_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_z[k] * ksi0_598[k]
                   - f_10 * pc_z[k] * ksi1_598[k];

        t_823[k] = f_11 * ksh_447[k]
                   + f_3 * pc_z[k] * lsh_615[k];

        t_824[k] = pa_x[k] * ksi0_824[k]
                   + f_12 * ksh_621[k]
                   - f_10 * pc_x[k] * ksi1_824[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, pa_x, pc_x, pc_y, ksi0_826, ksh_471, \
                         ksh_623, ksh_624, ksh_625, ksi1_826, lsh_618, lsh_624, \
                         lsh_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_18 * ksh_471[k]
                   + f_3 * pc_y[k] * lsh_618[k];

        t_826[k] = pa_x[k] * ksi0_826[k]
                   + f_12 * ksh_623[k]
                   - f_10 * pc_x[k] * ksi1_826[k];

        t_827[k] = f_11 * ksh_624[k]
                   + f_3 * pc_x[k] * lsh_624[k];

        t_828[k] = f_11 * ksh_625[k]
                   + f_3 * pc_x[k] * lsh_625[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, pc_x, ksh_626, ksh_627, ksh_628, ksh_629, \
                         lsh_626, lsh_627, lsh_628, lsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_11 * ksh_626[k]
                   + f_3 * pc_x[k] * lsh_626[k];

        t_830[k] = f_11 * ksh_627[k]
                   + f_3 * pc_x[k] * lsh_627[k];

        t_831[k] = f_11 * ksh_628[k]
                   + f_3 * pc_x[k] * lsh_628[k];

        t_832[k] = f_11 * ksh_629[k]
                   + f_3 * pc_x[k] * lsh_629[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, t_836, pa_x, pc_x, pc_z, ksi0_833, ksi0_835, \
                         ksi0_836, ksh_456, ksi1_833, ksi1_835, ksi1_836, \
                         lsh_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = pa_x[k] * ksi0_833[k]
                   - f_10 * pc_x[k] * ksi1_833[k];

        t_834[k] = f_11 * ksh_456[k]
                   + f_3 * pc_z[k] * lsh_624[k];

        t_835[k] = pa_x[k] * ksi0_835[k]
                   - f_10 * pc_x[k] * ksi1_835[k];

        t_836[k] = pa_x[k] * ksi0_836[k]
                   - f_10 * pc_x[k] * ksi1_836[k];
    }

#pragma omp simd aligned(t_837, t_838, t_839, t_840, pa_x, pc_x, pc_y, ksi0_837, ksi0_839, \
                         ksi0_840, ksh_482, ksh_630, ksi1_837, ksi1_839, ksi1_840, \
                         lsh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_837[k] = pa_x[k] * ksi0_837[k]
                   - f_10 * pc_x[k] * ksi1_837[k];

        t_838[k] = f_18 * ksh_482[k]
                   + f_3 * pc_y[k] * lsh_629[k];

        t_839[k] = pa_x[k] * ksi0_839[k]
                   - f_10 * pc_x[k] * ksi1_839[k];

        t_840[k] = pa_x[k] * ksi0_840[k]
                   + f_18 * ksh_630[k]
                   - f_10 * pc_x[k] * ksi1_840[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pa_x, pc_x, pc_y, pc_z, ksi0_843, \
                         ksh_462, ksh_483, ksh_485, ksh_633, ksi1_843, lsh_630, \
                         lsh_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_19 * ksh_483[k]
                   + f_3 * pc_y[k] * lsh_630[k];

        t_842[k] = f_12 * ksh_462[k]
                   + f_3 * pc_z[k] * lsh_630[k];

        t_843[k] = pa_x[k] * ksi0_843[k]
                   + f_14 * ksh_633[k]
                   - f_10 * pc_x[k] * ksi1_843[k];

        t_844[k] = f_19 * ksh_485[k]
                   + f_3 * pc_y[k] * lsh_632[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pa_x, pc_x, pc_z, ksi0_845, ksi0_846, ksh_465, \
                         ksh_635, ksh_636, ksi1_845, ksi1_846, \
                         lsh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pa_x[k] * ksi0_845[k]
                   + f_14 * ksh_635[k]
                   - f_10 * pc_x[k] * ksi1_845[k];

        t_846[k] = pa_x[k] * ksi0_846[k]
                   + f_13 * ksh_636[k]
                   - f_10 * pc_x[k] * ksi1_846[k];

        t_847[k] = f_12 * ksh_465[k]
                   + f_3 * pc_z[k] * lsh_633[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pa_x, pc_x, pc_y, ksi0_849, ksi0_850, ksh_488, \
                         ksh_639, ksh_640, ksi1_849, ksi1_850, \
                         lsh_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_19 * ksh_488[k]
                   + f_3 * pc_y[k] * lsh_635[k];

        t_849[k] = pa_x[k] * ksi0_849[k]
                   + f_13 * ksh_639[k]
                   - f_10 * pc_x[k] * ksi1_849[k];

        t_850[k] = pa_x[k] * ksi0_850[k]
                   + f_12 * ksh_640[k]
                   - f_10 * pc_x[k] * ksi1_850[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pa_x, pc_x, pc_y, pc_z, ksi0_852, ksh_468, \
                         ksh_492, ksh_642, ksi1_852, lsh_636, lsh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_12 * ksh_468[k]
                   + f_3 * pc_z[k] * lsh_636[k];

        t_852[k] = pa_x[k] * ksi0_852[k]
                   + f_12 * ksh_642[k]
                   - f_10 * pc_x[k] * ksi1_852[k];

        t_853[k] = f_19 * ksh_492[k]
                   + f_3 * pc_y[k] * lsh_639[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pa_x, pc_x, ksi0_854, ksh_644, ksh_645, \
                         ksh_646, ksh_647, ksi1_854, lsh_645, lsh_646, \
                         lsh_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = pa_x[k] * ksi0_854[k]
                   + f_12 * ksh_644[k]
                   - f_10 * pc_x[k] * ksi1_854[k];

        t_855[k] = f_11 * ksh_645[k]
                   + f_3 * pc_x[k] * lsh_645[k];

        t_856[k] = f_11 * ksh_646[k]
                   + f_3 * pc_x[k] * lsh_646[k];

        t_857[k] = f_11 * ksh_647[k]
                   + f_3 * pc_x[k] * lsh_647[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pa_x, pc_x, ksi0_861, ksh_648, ksh_649, \
                         ksh_650, ksi1_861, lsh_648, lsh_649, lsh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_11 * ksh_648[k]
                   + f_3 * pc_x[k] * lsh_648[k];

        t_859[k] = f_11 * ksh_649[k]
                   + f_3 * pc_x[k] * lsh_649[k];

        t_860[k] = f_11 * ksh_650[k]
                   + f_3 * pc_x[k] * lsh_650[k];

        t_861[k] = pa_x[k] * ksi0_861[k]
                   - f_10 * pc_x[k] * ksi1_861[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pa_x, pc_x, pc_z, ksi0_863, ksi0_864, \
                         ksi0_865, ksh_477, ksi1_863, ksi1_864, ksi1_865, \
                         lsh_645 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_12 * ksh_477[k]
                   + f_3 * pc_z[k] * lsh_645[k];

        t_863[k] = pa_x[k] * ksi0_863[k]
                   - f_10 * pc_x[k] * ksi1_863[k];

        t_864[k] = pa_x[k] * ksi0_864[k]
                   - f_10 * pc_x[k] * ksi1_864[k];

        t_865[k] = pa_x[k] * ksi0_865[k]
                   - f_10 * pc_x[k] * ksi1_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pa_x, pc_x, pc_y, ksi0_867, ksi0_868, \
                         ksh_503, ksh_504, ksh_651, ksi1_867, ksi1_868, lsh_650, \
                         lsh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_19 * ksh_503[k]
                   + f_3 * pc_y[k] * lsh_650[k];

        t_867[k] = pa_x[k] * ksi0_867[k]
                   - f_10 * pc_x[k] * ksi1_867[k];

        t_868[k] = pa_x[k] * ksi0_868[k]
                   + f_18 * ksh_651[k]
                   - f_10 * pc_x[k] * ksi1_868[k];

        t_869[k] = f_14 * ksh_504[k]
                   + f_3 * pc_y[k] * lsh_651[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pa_x, pc_x, pc_y, pc_z, ksi0_871, ksh_483, \
                         ksh_506, ksh_654, ksi1_871, lsh_651, lsh_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = f_13 * ksh_483[k]
                   + f_3 * pc_z[k] * lsh_651[k];

        t_871[k] = pa_x[k] * ksi0_871[k]
                   + f_14 * ksh_654[k]
                   - f_10 * pc_x[k] * ksi1_871[k];

        t_872[k] = f_14 * ksh_506[k]
                   + f_3 * pc_y[k] * lsh_653[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, pa_x, pc_x, pc_z, ksi0_873, ksi0_874, ksh_486, \
                         ksh_656, ksh_657, ksi1_873, ksi1_874, \
                         lsh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pa_x[k] * ksi0_873[k]
                   + f_14 * ksh_656[k]
                   - f_10 * pc_x[k] * ksi1_873[k];

        t_874[k] = pa_x[k] * ksi0_874[k]
                   + f_13 * ksh_657[k]
                   - f_10 * pc_x[k] * ksi1_874[k];

        t_875[k] = f_13 * ksh_486[k]
                   + f_3 * pc_z[k] * lsh_654[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_x, pc_x, pc_y, ksi0_877, ksi0_878, ksh_509, \
                         ksh_660, ksh_661, ksi1_877, ksi1_878, \
                         lsh_656 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = f_14 * ksh_509[k]
                   + f_3 * pc_y[k] * lsh_656[k];

        t_877[k] = pa_x[k] * ksi0_877[k]
                   + f_13 * ksh_660[k]
                   - f_10 * pc_x[k] * ksi1_877[k];

        t_878[k] = pa_x[k] * ksi0_878[k]
                   + f_12 * ksh_661[k]
                   - f_10 * pc_x[k] * ksi1_878[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pa_x, pc_x, pc_y, pc_z, ksi0_880, ksh_489, \
                         ksh_513, ksh_663, ksi1_880, lsh_657, lsh_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_13 * ksh_489[k]
                   + f_3 * pc_z[k] * lsh_657[k];

        t_880[k] = pa_x[k] * ksi0_880[k]
                   + f_12 * ksh_663[k]
                   - f_10 * pc_x[k] * ksi1_880[k];

        t_881[k] = f_14 * ksh_513[k]
                   + f_3 * pc_y[k] * lsh_660[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, t_885, pa_x, pc_x, ksi0_882, ksh_665, ksh_666, \
                         ksh_667, ksh_668, ksi1_882, lsh_666, lsh_667, \
                         lsh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pa_x[k] * ksi0_882[k]
                   + f_12 * ksh_665[k]
                   - f_10 * pc_x[k] * ksi1_882[k];

        t_883[k] = f_11 * ksh_666[k]
                   + f_3 * pc_x[k] * lsh_666[k];

        t_884[k] = f_11 * ksh_667[k]
                   + f_3 * pc_x[k] * lsh_667[k];

        t_885[k] = f_11 * ksh_668[k]
                   + f_3 * pc_x[k] * lsh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, pa_x, pc_x, ksi0_889, ksh_669, ksh_670, \
                         ksh_671, ksi1_889, lsh_669, lsh_670, lsh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_11 * ksh_669[k]
                   + f_3 * pc_x[k] * lsh_669[k];

        t_887[k] = f_11 * ksh_670[k]
                   + f_3 * pc_x[k] * lsh_670[k];

        t_888[k] = f_11 * ksh_671[k]
                   + f_3 * pc_x[k] * lsh_671[k];

        t_889[k] = pa_x[k] * ksi0_889[k]
                   - f_10 * pc_x[k] * ksi1_889[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pa_x, pc_x, pc_z, ksi0_891, ksi0_892, \
                         ksi0_893, ksh_498, ksi1_891, ksi1_892, ksi1_893, \
                         lsh_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_13 * ksh_498[k]
                   + f_3 * pc_z[k] * lsh_666[k];

        t_891[k] = pa_x[k] * ksi0_891[k]
                   - f_10 * pc_x[k] * ksi1_891[k];

        t_892[k] = pa_x[k] * ksi0_892[k]
                   - f_10 * pc_x[k] * ksi1_892[k];

        t_893[k] = pa_x[k] * ksi0_893[k]
                   - f_10 * pc_x[k] * ksi1_893[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pa_x, pc_x, pc_y, ksi0_895, ksi0_896, \
                         ksh_524, ksh_525, ksh_672, ksi1_895, ksi1_896, lsh_671, \
                         lsh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_14 * ksh_524[k]
                   + f_3 * pc_y[k] * lsh_671[k];

        t_895[k] = pa_x[k] * ksi0_895[k]
                   - f_10 * pc_x[k] * ksi1_895[k];

        t_896[k] = pa_x[k] * ksi0_896[k]
                   + f_18 * ksh_672[k]
                   - f_10 * pc_x[k] * ksi1_896[k];

        t_897[k] = f_13 * ksh_525[k]
                   + f_3 * pc_y[k] * lsh_672[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, pa_x, pc_x, pc_y, pc_z, ksi0_899, ksh_504, \
                         ksh_527, ksh_675, ksi1_899, lsh_672, lsh_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_14 * ksh_504[k]
                   + f_3 * pc_z[k] * lsh_672[k];

        t_899[k] = pa_x[k] * ksi0_899[k]
                   + f_14 * ksh_675[k]
                   - f_10 * pc_x[k] * ksi1_899[k];

        t_900[k] = f_13 * ksh_527[k]
                   + f_3 * pc_y[k] * lsh_674[k];
    }

#pragma omp simd aligned(t_901, t_902, t_903, pa_x, pc_x, pc_z, ksi0_901, ksi0_902, ksh_507, \
                         ksh_677, ksh_678, ksi1_901, ksi1_902, \
                         lsh_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_901[k] = pa_x[k] * ksi0_901[k]
                   + f_14 * ksh_677[k]
                   - f_10 * pc_x[k] * ksi1_901[k];

        t_902[k] = pa_x[k] * ksi0_902[k]
                   + f_13 * ksh_678[k]
                   - f_10 * pc_x[k] * ksi1_902[k];

        t_903[k] = f_14 * ksh_507[k]
                   + f_3 * pc_z[k] * lsh_675[k];
    }

#pragma omp simd aligned(t_904, t_905, t_906, pa_x, pc_x, pc_y, ksi0_905, ksi0_906, ksh_530, \
                         ksh_681, ksh_682, ksi1_905, ksi1_906, \
                         lsh_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_904[k] = f_13 * ksh_530[k]
                   + f_3 * pc_y[k] * lsh_677[k];

        t_905[k] = pa_x[k] * ksi0_905[k]
                   + f_13 * ksh_681[k]
                   - f_10 * pc_x[k] * ksi1_905[k];

        t_906[k] = pa_x[k] * ksi0_906[k]
                   + f_12 * ksh_682[k]
                   - f_10 * pc_x[k] * ksi1_906[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, pa_x, pc_x, pc_y, pc_z, ksi0_908, ksh_510, \
                         ksh_534, ksh_684, ksi1_908, lsh_678, lsh_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_14 * ksh_510[k]
                   + f_3 * pc_z[k] * lsh_678[k];

        t_908[k] = pa_x[k] * ksi0_908[k]
                   + f_12 * ksh_684[k]
                   - f_10 * pc_x[k] * ksi1_908[k];

        t_909[k] = f_13 * ksh_534[k]
                   + f_3 * pc_y[k] * lsh_681[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pa_x, pc_x, ksi0_910, ksh_686, ksh_687, \
                         ksh_688, ksh_689, ksi1_910, lsh_687, lsh_688, \
                         lsh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = pa_x[k] * ksi0_910[k]
                   + f_12 * ksh_686[k]
                   - f_10 * pc_x[k] * ksi1_910[k];

        t_911[k] = f_11 * ksh_687[k]
                   + f_3 * pc_x[k] * lsh_687[k];

        t_912[k] = f_11 * ksh_688[k]
                   + f_3 * pc_x[k] * lsh_688[k];

        t_913[k] = f_11 * ksh_689[k]
                   + f_3 * pc_x[k] * lsh_689[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, pa_x, pc_x, ksi0_917, ksh_690, ksh_691, \
                         ksh_692, ksi1_917, lsh_690, lsh_691, lsh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_11 * ksh_690[k]
                   + f_3 * pc_x[k] * lsh_690[k];

        t_915[k] = f_11 * ksh_691[k]
                   + f_3 * pc_x[k] * lsh_691[k];

        t_916[k] = f_11 * ksh_692[k]
                   + f_3 * pc_x[k] * lsh_692[k];

        t_917[k] = pa_x[k] * ksi0_917[k]
                   - f_10 * pc_x[k] * ksi1_917[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, t_921, pa_x, pc_x, pc_z, ksi0_919, ksi0_920, \
                         ksi0_921, ksh_519, ksi1_919, ksi1_920, ksi1_921, \
                         lsh_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_14 * ksh_519[k]
                   + f_3 * pc_z[k] * lsh_687[k];

        t_919[k] = pa_x[k] * ksi0_919[k]
                   - f_10 * pc_x[k] * ksi1_919[k];

        t_920[k] = pa_x[k] * ksi0_920[k]
                   - f_10 * pc_x[k] * ksi1_920[k];

        t_921[k] = pa_x[k] * ksi0_921[k]
                   - f_10 * pc_x[k] * ksi1_921[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, pa_x, pc_x, pc_y, ksi0_923, ksi0_924, \
                         ksh_545, ksh_546, ksh_693, ksi1_923, ksi1_924, lsh_692, \
                         lsh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * ksh_545[k]
                   + f_3 * pc_y[k] * lsh_692[k];

        t_923[k] = pa_x[k] * ksi0_923[k]
                   - f_10 * pc_x[k] * ksi1_923[k];

        t_924[k] = pa_x[k] * ksi0_924[k]
                   + f_18 * ksh_693[k]
                   - f_10 * pc_x[k] * ksi1_924[k];

        t_925[k] = f_12 * ksh_546[k]
                   + f_3 * pc_y[k] * lsh_693[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, pa_x, pc_x, pc_y, pc_z, ksi0_927, ksh_525, \
                         ksh_548, ksh_696, ksi1_927, lsh_693, lsh_695 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_19 * ksh_525[k]
                   + f_3 * pc_z[k] * lsh_693[k];

        t_927[k] = pa_x[k] * ksi0_927[k]
                   + f_14 * ksh_696[k]
                   - f_10 * pc_x[k] * ksi1_927[k];

        t_928[k] = f_12 * ksh_548[k]
                   + f_3 * pc_y[k] * lsh_695[k];
    }

#pragma omp simd aligned(t_929, t_930, t_931, pa_x, pc_x, pc_z, ksi0_929, ksi0_930, ksh_528, \
                         ksh_698, ksh_699, ksi1_929, ksi1_930, \
                         lsh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_929[k] = pa_x[k] * ksi0_929[k]
                   + f_14 * ksh_698[k]
                   - f_10 * pc_x[k] * ksi1_929[k];

        t_930[k] = pa_x[k] * ksi0_930[k]
                   + f_13 * ksh_699[k]
                   - f_10 * pc_x[k] * ksi1_930[k];

        t_931[k] = f_19 * ksh_528[k]
                   + f_3 * pc_z[k] * lsh_696[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pa_x, pc_x, pc_y, ksi0_933, ksi0_934, ksh_551, \
                         ksh_702, ksh_703, ksi1_933, ksi1_934, \
                         lsh_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_12 * ksh_551[k]
                   + f_3 * pc_y[k] * lsh_698[k];

        t_933[k] = pa_x[k] * ksi0_933[k]
                   + f_13 * ksh_702[k]
                   - f_10 * pc_x[k] * ksi1_933[k];

        t_934[k] = pa_x[k] * ksi0_934[k]
                   + f_12 * ksh_703[k]
                   - f_10 * pc_x[k] * ksi1_934[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, pa_x, pc_x, pc_y, pc_z, ksi0_936, ksh_531, \
                         ksh_555, ksh_705, ksi1_936, lsh_699, lsh_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_19 * ksh_531[k]
                   + f_3 * pc_z[k] * lsh_699[k];

        t_936[k] = pa_x[k] * ksi0_936[k]
                   + f_12 * ksh_705[k]
                   - f_10 * pc_x[k] * ksi1_936[k];

        t_937[k] = f_12 * ksh_555[k]
                   + f_3 * pc_y[k] * lsh_702[k];
    }

#pragma omp simd aligned(t_938, t_939, t_940, t_941, pa_x, pc_x, ksi0_938, ksh_707, ksh_708, \
                         ksh_709, ksh_710, ksi1_938, lsh_708, lsh_709, \
                         lsh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_938[k] = pa_x[k] * ksi0_938[k]
                   + f_12 * ksh_707[k]
                   - f_10 * pc_x[k] * ksi1_938[k];

        t_939[k] = f_11 * ksh_708[k]
                   + f_3 * pc_x[k] * lsh_708[k];

        t_940[k] = f_11 * ksh_709[k]
                   + f_3 * pc_x[k] * lsh_709[k];

        t_941[k] = f_11 * ksh_710[k]
                   + f_3 * pc_x[k] * lsh_710[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t ksi0,
                                                          const size_t ksh, const size_t ksi1,
                                                          const size_t lsg0, const size_t lsg1,
                                                          const size_t lsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_942 = buffer.data(target + 942);
    auto *t_943 = buffer.data(target + 943);
    auto *t_944 = buffer.data(target + 944);
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);
    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_756 = buffer.data(ksi0 + 756);
    const auto *ksi0_761 = buffer.data(ksi0 + 761);
    const auto *ksi0_765 = buffer.data(ksi0 + 765);
    const auto *ksi0_770 = buffer.data(ksi0 + 770);
    const auto *ksi0_784 = buffer.data(ksi0 + 784);
    const auto *ksi0_785 = buffer.data(ksi0 + 785);
    const auto *ksi0_787 = buffer.data(ksi0 + 787);
    const auto *ksi0_790 = buffer.data(ksi0 + 790);
    const auto *ksi0_794 = buffer.data(ksi0 + 794);
    const auto *ksi0_805 = buffer.data(ksi0 + 805);
    const auto *ksi0_807 = buffer.data(ksi0 + 807);
    const auto *ksi0_808 = buffer.data(ksi0 + 808);
    const auto *ksi0_809 = buffer.data(ksi0 + 809);
    const auto *ksi0_945 = buffer.data(ksi0 + 945);
    const auto *ksi0_947 = buffer.data(ksi0 + 947);
    const auto *ksi0_948 = buffer.data(ksi0 + 948);
    const auto *ksi0_949 = buffer.data(ksi0 + 949);
    const auto *ksi0_951 = buffer.data(ksi0 + 951);
    const auto *ksi0_955 = buffer.data(ksi0 + 955);
    const auto *ksi0_958 = buffer.data(ksi0 + 958);
    const auto *ksi0_962 = buffer.data(ksi0 + 962);
    const auto *ksi0_964 = buffer.data(ksi0 + 964);
    const auto *ksi0_973 = buffer.data(ksi0 + 973);
    const auto *ksi0_975 = buffer.data(ksi0 + 975);
    const auto *ksi0_976 = buffer.data(ksi0 + 976);
    const auto *ksi0_977 = buffer.data(ksi0 + 977);
    const auto *ksi0_979 = buffer.data(ksi0 + 979);
    const auto *ksi0_980 = buffer.data(ksi0 + 980);
    const auto *ksi0_985 = buffer.data(ksi0 + 985);
    const auto *ksi0_989 = buffer.data(ksi0 + 989);
    const auto *ksi0_994 = buffer.data(ksi0 + 994);
    const auto *ksi0_1001 = buffer.data(ksi0 + 1001);
    const auto *ksi0_1002 = buffer.data(ksi0 + 1002);
    const auto *ksi0_1003 = buffer.data(ksi0 + 1003);
    const auto *ksi0_1004 = buffer.data(ksi0 + 1004);
    const auto *ksi0_1005 = buffer.data(ksi0 + 1005);
    const auto *ksi0_1007 = buffer.data(ksi0 + 1007);

    const auto *ksh_540 = buffer.data(ksh + 540);
    const auto *ksh_546 = buffer.data(ksh + 546);
    const auto *ksh_549 = buffer.data(ksh + 549);
    const auto *ksh_552 = buffer.data(ksh + 552);
    const auto *ksh_561 = buffer.data(ksh + 561);
    const auto *ksh_566 = buffer.data(ksh + 566);
    const auto *ksh_567 = buffer.data(ksh + 567);
    const auto *ksh_569 = buffer.data(ksh + 569);
    const auto *ksh_572 = buffer.data(ksh + 572);
    const auto *ksh_576 = buffer.data(ksh + 576);
    const auto *ksh_587 = buffer.data(ksh + 587);
    const auto *ksh_603 = buffer.data(ksh + 603);
    const auto *ksh_604 = buffer.data(ksh + 604);
    const auto *ksh_605 = buffer.data(ksh + 605);
    const auto *ksh_606 = buffer.data(ksh + 606);
    const auto *ksh_608 = buffer.data(ksh + 608);
    const auto *ksh_629 = buffer.data(ksh + 629);
    const auto *ksh_711 = buffer.data(ksh + 711);
    const auto *ksh_712 = buffer.data(ksh + 712);
    const auto *ksh_713 = buffer.data(ksh + 713);
    const auto *ksh_717 = buffer.data(ksh + 717);
    const auto *ksh_720 = buffer.data(ksh + 720);
    const auto *ksh_724 = buffer.data(ksh + 724);
    const auto *ksh_726 = buffer.data(ksh + 726);
    const auto *ksh_729 = buffer.data(ksh + 729);
    const auto *ksh_730 = buffer.data(ksh + 730);
    const auto *ksh_731 = buffer.data(ksh + 731);
    const auto *ksh_732 = buffer.data(ksh + 732);
    const auto *ksh_733 = buffer.data(ksh + 733);
    const auto *ksh_734 = buffer.data(ksh + 734);
    const auto *ksh_735 = buffer.data(ksh + 735);
    const auto *ksh_740 = buffer.data(ksh + 740);
    const auto *ksh_744 = buffer.data(ksh + 744);
    const auto *ksh_749 = buffer.data(ksh + 749);
    const auto *ksh_750 = buffer.data(ksh + 750);
    const auto *ksh_751 = buffer.data(ksh + 751);
    const auto *ksh_752 = buffer.data(ksh + 752);
    const auto *ksh_753 = buffer.data(ksh + 753);
    const auto *ksh_755 = buffer.data(ksh + 755);

    const auto *ksi1_756 = buffer.data(ksi1 + 756);
    const auto *ksi1_761 = buffer.data(ksi1 + 761);
    const auto *ksi1_765 = buffer.data(ksi1 + 765);
    const auto *ksi1_770 = buffer.data(ksi1 + 770);
    const auto *ksi1_784 = buffer.data(ksi1 + 784);
    const auto *ksi1_785 = buffer.data(ksi1 + 785);
    const auto *ksi1_787 = buffer.data(ksi1 + 787);
    const auto *ksi1_790 = buffer.data(ksi1 + 790);
    const auto *ksi1_794 = buffer.data(ksi1 + 794);
    const auto *ksi1_805 = buffer.data(ksi1 + 805);
    const auto *ksi1_807 = buffer.data(ksi1 + 807);
    const auto *ksi1_808 = buffer.data(ksi1 + 808);
    const auto *ksi1_809 = buffer.data(ksi1 + 809);
    const auto *ksi1_945 = buffer.data(ksi1 + 945);
    const auto *ksi1_947 = buffer.data(ksi1 + 947);
    const auto *ksi1_948 = buffer.data(ksi1 + 948);
    const auto *ksi1_949 = buffer.data(ksi1 + 949);
    const auto *ksi1_951 = buffer.data(ksi1 + 951);
    const auto *ksi1_955 = buffer.data(ksi1 + 955);
    const auto *ksi1_958 = buffer.data(ksi1 + 958);
    const auto *ksi1_962 = buffer.data(ksi1 + 962);
    const auto *ksi1_964 = buffer.data(ksi1 + 964);
    const auto *ksi1_973 = buffer.data(ksi1 + 973);
    const auto *ksi1_975 = buffer.data(ksi1 + 975);
    const auto *ksi1_976 = buffer.data(ksi1 + 976);
    const auto *ksi1_977 = buffer.data(ksi1 + 977);
    const auto *ksi1_979 = buffer.data(ksi1 + 979);
    const auto *ksi1_980 = buffer.data(ksi1 + 980);
    const auto *ksi1_985 = buffer.data(ksi1 + 985);
    const auto *ksi1_989 = buffer.data(ksi1 + 989);
    const auto *ksi1_994 = buffer.data(ksi1 + 994);
    const auto *ksi1_1001 = buffer.data(ksi1 + 1001);
    const auto *ksi1_1002 = buffer.data(ksi1 + 1002);
    const auto *ksi1_1003 = buffer.data(ksi1 + 1003);
    const auto *ksi1_1004 = buffer.data(ksi1 + 1004);
    const auto *ksi1_1005 = buffer.data(ksi1 + 1005);
    const auto *ksi1_1007 = buffer.data(ksi1 + 1007);

    const auto *lsg0_525 = buffer.data(lsg0 + 525);
    const auto *lsg0_526 = buffer.data(lsg0 + 526);
    const auto *lsg0_527 = buffer.data(lsg0 + 527);
    const auto *lsg0_528 = buffer.data(lsg0 + 528);
    const auto *lsg0_529 = buffer.data(lsg0 + 529);
    const auto *lsg0_530 = buffer.data(lsg0 + 530);
    const auto *lsg0_540 = buffer.data(lsg0 + 540);
    const auto *lsg0_541 = buffer.data(lsg0 + 541);
    const auto *lsg0_543 = buffer.data(lsg0 + 543);
    const auto *lsg0_545 = buffer.data(lsg0 + 545);
    const auto *lsg0_546 = buffer.data(lsg0 + 546);
    const auto *lsg0_548 = buffer.data(lsg0 + 548);
    const auto *lsg0_549 = buffer.data(lsg0 + 549);
    const auto *lsg0_550 = buffer.data(lsg0 + 550);
    const auto *lsg0_551 = buffer.data(lsg0 + 551);
    const auto *lsg0_552 = buffer.data(lsg0 + 552);
    const auto *lsg0_553 = buffer.data(lsg0 + 553);
    const auto *lsg0_554 = buffer.data(lsg0 + 554);
    const auto *lsg0_557 = buffer.data(lsg0 + 557);
    const auto *lsg0_559 = buffer.data(lsg0 + 559);
    const auto *lsg0_560 = buffer.data(lsg0 + 560);
    const auto *lsg0_562 = buffer.data(lsg0 + 562);
    const auto *lsg0_563 = buffer.data(lsg0 + 563);
    const auto *lsg0_564 = buffer.data(lsg0 + 564);
    const auto *lsg0_566 = buffer.data(lsg0 + 566);
    const auto *lsg0_567 = buffer.data(lsg0 + 567);
    const auto *lsg0_568 = buffer.data(lsg0 + 568);
    const auto *lsg0_569 = buffer.data(lsg0 + 569);
    const auto *lsg0_570 = buffer.data(lsg0 + 570);
    const auto *lsg0_571 = buffer.data(lsg0 + 571);
    const auto *lsg0_572 = buffer.data(lsg0 + 572);
    const auto *lsg0_573 = buffer.data(lsg0 + 573);

    const auto *lsg1_525 = buffer.data(lsg1 + 525);
    const auto *lsg1_526 = buffer.data(lsg1 + 526);
    const auto *lsg1_527 = buffer.data(lsg1 + 527);
    const auto *lsg1_528 = buffer.data(lsg1 + 528);
    const auto *lsg1_529 = buffer.data(lsg1 + 529);
    const auto *lsg1_530 = buffer.data(lsg1 + 530);
    const auto *lsg1_540 = buffer.data(lsg1 + 540);
    const auto *lsg1_541 = buffer.data(lsg1 + 541);
    const auto *lsg1_543 = buffer.data(lsg1 + 543);
    const auto *lsg1_545 = buffer.data(lsg1 + 545);
    const auto *lsg1_546 = buffer.data(lsg1 + 546);
    const auto *lsg1_548 = buffer.data(lsg1 + 548);
    const auto *lsg1_549 = buffer.data(lsg1 + 549);
    const auto *lsg1_550 = buffer.data(lsg1 + 550);
    const auto *lsg1_551 = buffer.data(lsg1 + 551);
    const auto *lsg1_552 = buffer.data(lsg1 + 552);
    const auto *lsg1_553 = buffer.data(lsg1 + 553);
    const auto *lsg1_554 = buffer.data(lsg1 + 554);
    const auto *lsg1_557 = buffer.data(lsg1 + 557);
    const auto *lsg1_559 = buffer.data(lsg1 + 559);
    const auto *lsg1_560 = buffer.data(lsg1 + 560);
    const auto *lsg1_562 = buffer.data(lsg1 + 562);
    const auto *lsg1_563 = buffer.data(lsg1 + 563);
    const auto *lsg1_564 = buffer.data(lsg1 + 564);
    const auto *lsg1_566 = buffer.data(lsg1 + 566);
    const auto *lsg1_567 = buffer.data(lsg1 + 567);
    const auto *lsg1_568 = buffer.data(lsg1 + 568);
    const auto *lsg1_569 = buffer.data(lsg1 + 569);
    const auto *lsg1_570 = buffer.data(lsg1 + 570);
    const auto *lsg1_571 = buffer.data(lsg1 + 571);
    const auto *lsg1_572 = buffer.data(lsg1 + 572);
    const auto *lsg1_573 = buffer.data(lsg1 + 573);

    const auto *lsh_708 = buffer.data(lsh + 708);
    const auto *lsh_711 = buffer.data(lsh + 711);
    const auto *lsh_712 = buffer.data(lsh + 712);
    const auto *lsh_713 = buffer.data(lsh + 713);
    const auto *lsh_714 = buffer.data(lsh + 714);
    const auto *lsh_716 = buffer.data(lsh + 716);
    const auto *lsh_717 = buffer.data(lsh + 717);
    const auto *lsh_719 = buffer.data(lsh + 719);
    const auto *lsh_720 = buffer.data(lsh + 720);
    const auto *lsh_723 = buffer.data(lsh + 723);
    const auto *lsh_729 = buffer.data(lsh + 729);
    const auto *lsh_730 = buffer.data(lsh + 730);
    const auto *lsh_731 = buffer.data(lsh + 731);
    const auto *lsh_732 = buffer.data(lsh + 732);
    const auto *lsh_733 = buffer.data(lsh + 733);
    const auto *lsh_734 = buffer.data(lsh + 734);
    const auto *lsh_735 = buffer.data(lsh + 735);
    const auto *lsh_736 = buffer.data(lsh + 736);
    const auto *lsh_737 = buffer.data(lsh + 737);
    const auto *lsh_738 = buffer.data(lsh + 738);
    const auto *lsh_739 = buffer.data(lsh + 739);
    const auto *lsh_740 = buffer.data(lsh + 740);
    const auto *lsh_741 = buffer.data(lsh + 741);
    const auto *lsh_742 = buffer.data(lsh + 742);
    const auto *lsh_743 = buffer.data(lsh + 743);
    const auto *lsh_744 = buffer.data(lsh + 744);
    const auto *lsh_749 = buffer.data(lsh + 749);
    const auto *lsh_750 = buffer.data(lsh + 750);
    const auto *lsh_751 = buffer.data(lsh + 751);
    const auto *lsh_752 = buffer.data(lsh + 752);
    const auto *lsh_753 = buffer.data(lsh + 753);
    const auto *lsh_755 = buffer.data(lsh + 755);
    const auto *lsh_756 = buffer.data(lsh + 756);
    const auto *lsh_757 = buffer.data(lsh + 757);
    const auto *lsh_759 = buffer.data(lsh + 759);
    const auto *lsh_761 = buffer.data(lsh + 761);
    const auto *lsh_762 = buffer.data(lsh + 762);
    const auto *lsh_764 = buffer.data(lsh + 764);
    const auto *lsh_765 = buffer.data(lsh + 765);
    const auto *lsh_766 = buffer.data(lsh + 766);
    const auto *lsh_768 = buffer.data(lsh + 768);
    const auto *lsh_769 = buffer.data(lsh + 769);
    const auto *lsh_770 = buffer.data(lsh + 770);
    const auto *lsh_771 = buffer.data(lsh + 771);
    const auto *lsh_772 = buffer.data(lsh + 772);
    const auto *lsh_773 = buffer.data(lsh + 773);
    const auto *lsh_774 = buffer.data(lsh + 774);
    const auto *lsh_775 = buffer.data(lsh + 775);
    const auto *lsh_776 = buffer.data(lsh + 776);
    const auto *lsh_779 = buffer.data(lsh + 779);
    const auto *lsh_781 = buffer.data(lsh + 781);
    const auto *lsh_782 = buffer.data(lsh + 782);
    const auto *lsh_784 = buffer.data(lsh + 784);
    const auto *lsh_785 = buffer.data(lsh + 785);
    const auto *lsh_786 = buffer.data(lsh + 786);
    const auto *lsh_788 = buffer.data(lsh + 788);
    const auto *lsh_789 = buffer.data(lsh + 789);
    const auto *lsh_790 = buffer.data(lsh + 790);
    const auto *lsh_791 = buffer.data(lsh + 791);
    const auto *lsh_792 = buffer.data(lsh + 792);
    const auto *lsh_793 = buffer.data(lsh + 793);
    const auto *lsh_794 = buffer.data(lsh + 794);
    const auto *lsh_795 = buffer.data(lsh + 795);
    const auto *lsh_796 = buffer.data(lsh + 796);
    const auto *lsh_797 = buffer.data(lsh + 797);
    const auto *lsh_798 = buffer.data(lsh + 798);
    const auto *lsh_799 = buffer.data(lsh + 799);
    const auto *lsh_800 = buffer.data(lsh + 800);
    const auto *lsh_801 = buffer.data(lsh + 801);

#pragma omp simd aligned(t_942, t_943, t_944, t_945, pa_x, pc_x, ksi0_945, ksh_711, ksh_712, \
                         ksh_713, ksi1_945, lsh_711, lsh_712, lsh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_11 * ksh_711[k]
                   + f_3 * pc_x[k] * lsh_711[k];

        t_943[k] = f_11 * ksh_712[k]
                   + f_3 * pc_x[k] * lsh_712[k];

        t_944[k] = f_11 * ksh_713[k]
                   + f_3 * pc_x[k] * lsh_713[k];

        t_945[k] = pa_x[k] * ksi0_945[k]
                   - f_10 * pc_x[k] * ksi1_945[k];
    }

#pragma omp simd aligned(t_946, t_947, t_948, t_949, pa_x, pc_x, pc_z, ksi0_947, ksi0_948, \
                         ksi0_949, ksh_540, ksi1_947, ksi1_948, ksi1_949, \
                         lsh_708 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_946[k] = f_19 * ksh_540[k]
                   + f_3 * pc_z[k] * lsh_708[k];

        t_947[k] = pa_x[k] * ksi0_947[k]
                   - f_10 * pc_x[k] * ksi1_947[k];

        t_948[k] = pa_x[k] * ksi0_948[k]
                   - f_10 * pc_x[k] * ksi1_948[k];

        t_949[k] = pa_x[k] * ksi0_949[k]
                   - f_10 * pc_x[k] * ksi1_949[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pa_x, pa_y, pc_x, pc_y, ksi0_756, \
                         ksi0_951, ksh_566, ksh_567, ksi1_756, ksi1_951, lsh_713, \
                         lsh_714 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_12 * ksh_566[k]
                   + f_3 * pc_y[k] * lsh_713[k];

        t_951[k] = pa_x[k] * ksi0_951[k]
                   - f_10 * pc_x[k] * ksi1_951[k];

        t_952[k] = pa_y[k] * ksi0_756[k]
                   - f_10 * pc_y[k] * ksi1_756[k];

        t_953[k] = f_11 * ksh_567[k]
                   + f_3 * pc_y[k] * lsh_714[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pa_x, pc_x, pc_y, pc_z, ksi0_955, ksh_546, \
                         ksh_569, ksh_717, ksi1_955, lsh_714, lsh_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_18 * ksh_546[k]
                   + f_3 * pc_z[k] * lsh_714[k];

        t_955[k] = pa_x[k] * ksi0_955[k]
                   + f_14 * ksh_717[k]
                   - f_10 * pc_x[k] * ksi1_955[k];

        t_956[k] = f_11 * ksh_569[k]
                   + f_3 * pc_y[k] * lsh_716[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pa_x, pa_y, pc_x, pc_y, pc_z, ksi0_761, \
                         ksi0_958, ksh_549, ksh_720, ksi1_761, ksi1_958, \
                         lsh_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = pa_y[k] * ksi0_761[k]
                   - f_10 * pc_y[k] * ksi1_761[k];

        t_958[k] = pa_x[k] * ksi0_958[k]
                   + f_13 * ksh_720[k]
                   - f_10 * pc_x[k] * ksi1_958[k];

        t_959[k] = f_18 * ksh_549[k]
                   + f_3 * pc_z[k] * lsh_717[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pa_x, pa_y, pc_x, pc_y, ksi0_765, ksi0_962, \
                         ksh_572, ksh_724, ksi1_765, ksi1_962, \
                         lsh_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_11 * ksh_572[k]
                   + f_3 * pc_y[k] * lsh_719[k];

        t_961[k] = pa_y[k] * ksi0_765[k]
                   - f_10 * pc_y[k] * ksi1_765[k];

        t_962[k] = pa_x[k] * ksi0_962[k]
                   + f_12 * ksh_724[k]
                   - f_10 * pc_x[k] * ksi1_962[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pa_x, pc_x, pc_y, pc_z, ksi0_964, ksh_552, \
                         ksh_576, ksh_726, ksi1_964, lsh_720, lsh_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_18 * ksh_552[k]
                   + f_3 * pc_z[k] * lsh_720[k];

        t_964[k] = pa_x[k] * ksi0_964[k]
                   + f_12 * ksh_726[k]
                   - f_10 * pc_x[k] * ksi1_964[k];

        t_965[k] = f_11 * ksh_576[k]
                   + f_3 * pc_y[k] * lsh_723[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, t_969, pa_y, pc_x, pc_y, ksi0_770, ksh_729, \
                         ksh_730, ksh_731, ksi1_770, lsh_729, lsh_730, \
                         lsh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = pa_y[k] * ksi0_770[k]
                   - f_10 * pc_y[k] * ksi1_770[k];

        t_967[k] = f_11 * ksh_729[k]
                   + f_3 * pc_x[k] * lsh_729[k];

        t_968[k] = f_11 * ksh_730[k]
                   + f_3 * pc_x[k] * lsh_730[k];

        t_969[k] = f_11 * ksh_731[k]
                   + f_3 * pc_x[k] * lsh_731[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, pa_x, pc_x, ksi0_973, ksh_732, ksh_733, \
                         ksh_734, ksi1_973, lsh_732, lsh_733, lsh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_11 * ksh_732[k]
                   + f_3 * pc_x[k] * lsh_732[k];

        t_971[k] = f_11 * ksh_733[k]
                   + f_3 * pc_x[k] * lsh_733[k];

        t_972[k] = f_11 * ksh_734[k]
                   + f_3 * pc_x[k] * lsh_734[k];

        t_973[k] = pa_x[k] * ksi0_973[k]
                   - f_10 * pc_x[k] * ksi1_973[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, pa_x, pc_x, pc_z, ksi0_975, ksi0_976, \
                         ksi0_977, ksh_561, ksi1_975, ksi1_976, ksi1_977, \
                         lsh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_18 * ksh_561[k]
                   + f_3 * pc_z[k] * lsh_729[k];

        t_975[k] = pa_x[k] * ksi0_975[k]
                   - f_10 * pc_x[k] * ksi1_975[k];

        t_976[k] = pa_x[k] * ksi0_976[k]
                   - f_10 * pc_x[k] * ksi1_976[k];

        t_977[k] = pa_x[k] * ksi0_977[k]
                   - f_10 * pc_x[k] * ksi1_977[k];
    }

#pragma omp simd aligned(t_978, t_979, t_980, t_981, pa_x, pc_x, pc_y, ksi0_979, ksi0_980, \
                         ksh_587, ksh_735, ksi1_979, ksi1_980, lsh_734, \
                         lsh_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_978[k] = f_11 * ksh_587[k]
                   + f_3 * pc_y[k] * lsh_734[k];

        t_979[k] = pa_x[k] * ksi0_979[k]
                   - f_10 * pc_x[k] * ksi1_979[k];

        t_980[k] = pa_x[k] * ksi0_980[k]
                   + f_18 * ksh_735[k]
                   - f_10 * pc_x[k] * ksi1_980[k];

        t_981[k] = f_3 * pc_y[k] * lsh_735[k];
    }

#pragma omp simd aligned(t_982, t_983, t_984, pc_y, pc_z, ksh_567, lsg0_525, lsg1_525, \
                         lsh_735, lsh_736, lsh_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_982[k] = f_15 * ksh_567[k]
                   + f_3 * pc_z[k] * lsh_735[k];

        t_983[k] = f_4 * lsg0_525[k]
                   - f_5 * lsg1_525[k]
                   + f_3 * pc_y[k] * lsh_736[k];

        t_984[k] = f_3 * pc_y[k] * lsh_737[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, pa_x, pc_x, pc_y, ksi0_985, ksh_740, ksi1_985, \
                         lsg0_526, lsg0_527, lsg1_526, lsg1_527, lsh_738, \
                         lsh_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = pa_x[k] * ksi0_985[k]
                   + f_14 * ksh_740[k]
                   - f_10 * pc_x[k] * ksi1_985[k];

        t_986[k] = f_6 * lsg0_526[k]
                   - f_7 * lsg1_526[k]
                   + f_3 * pc_y[k] * lsh_738[k];

        t_987[k] = f_4 * lsg0_527[k]
                   - f_5 * lsg1_527[k]
                   + f_3 * pc_y[k] * lsh_739[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pa_x, pc_x, pc_y, ksi0_989, ksh_744, ksi1_989, \
                         lsg0_528, lsg1_528, lsh_740, lsh_741 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_3 * pc_y[k] * lsh_740[k];

        t_989[k] = pa_x[k] * ksi0_989[k]
                   + f_13 * ksh_744[k]
                   - f_10 * pc_x[k] * ksi1_989[k];

        t_990[k] = f_8 * lsg0_528[k]
                   - f_9 * lsg1_528[k]
                   + f_3 * pc_y[k] * lsh_741[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, pc_y, lsg0_529, lsg0_530, lsg1_529, lsg1_530, \
                         lsh_742, lsh_743, lsh_744 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_6 * lsg0_529[k]
                   - f_7 * lsg1_529[k]
                   + f_3 * pc_y[k] * lsh_742[k];

        t_992[k] = f_4 * lsg0_530[k]
                   - f_5 * lsg1_530[k]
                   + f_3 * pc_y[k] * lsh_743[k];

        t_993[k] = f_3 * pc_y[k] * lsh_744[k];
    }

#pragma omp simd aligned(t_994, t_995, t_996, t_997, pa_x, pc_x, ksi0_994, ksh_749, ksh_750, \
                         ksh_751, ksh_752, ksi1_994, lsh_750, lsh_751, \
                         lsh_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_994[k] = pa_x[k] * ksi0_994[k]
                   + f_12 * ksh_749[k]
                   - f_10 * pc_x[k] * ksi1_994[k];

        t_995[k] = f_11 * ksh_750[k]
                   + f_3 * pc_x[k] * lsh_750[k];

        t_996[k] = f_11 * ksh_751[k]
                   + f_3 * pc_x[k] * lsh_751[k];

        t_997[k] = f_11 * ksh_752[k]
                   + f_3 * pc_x[k] * lsh_752[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, pa_x, pc_x, pc_y, ksi0_1001, ksh_753, \
                         ksh_755, ksi1_1001, lsh_749, lsh_753, \
                         lsh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_11 * ksh_753[k]
                   + f_3 * pc_x[k] * lsh_753[k];

        t_999[k] = f_3 * pc_y[k] * lsh_749[k];

        t_1000[k] = f_11 * ksh_755[k]
                    + f_3 * pc_x[k] * lsh_755[k];

        t_1001[k] = pa_x[k] * ksi0_1001[k]
                    - f_10 * pc_x[k] * ksi1_1001[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, t_1005, pa_x, pc_x, ksi0_1002, ksi0_1003, \
                         ksi0_1004, ksi0_1005, ksi1_1002, ksi1_1003, ksi1_1004, \
                         ksi1_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = pa_x[k] * ksi0_1002[k]
                    - f_10 * pc_x[k] * ksi1_1002[k];

        t_1003[k] = pa_x[k] * ksi0_1003[k]
                    - f_10 * pc_x[k] * ksi1_1003[k];

        t_1004[k] = pa_x[k] * ksi0_1004[k]
                    - f_10 * pc_x[k] * ksi1_1004[k];

        t_1005[k] = pa_x[k] * ksi0_1005[k]
                    - f_10 * pc_x[k] * ksi1_1005[k];
    }

#pragma omp simd aligned(t_1006, t_1007, t_1008, t_1009, pa_x, pc_x, pc_y, ksi0_1007, \
                         ksi1_1007, lsg0_540, lsg0_541, lsg1_540, lsg1_541, lsh_755, lsh_756, \
                         lsh_757 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1006[k] = f_3 * pc_y[k] * lsh_755[k];

        t_1007[k] = pa_x[k] * ksi0_1007[k]
                    - f_10 * pc_x[k] * ksi1_1007[k];

        t_1008[k] = f_1 * lsg0_540[k]
                    - f_2 * lsg1_540[k]
                    + f_3 * pc_x[k] * lsh_756[k];

        t_1009[k] = f_16 * lsg0_541[k]
                    - f_17 * lsg1_541[k]
                    + f_3 * pc_x[k] * lsh_757[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, pc_x, pc_z, lsg0_543, lsg0_545, \
                         lsg1_543, lsg1_545, lsh_756, lsh_757, lsh_759, \
                         lsh_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_3 * pc_z[k] * lsh_756[k];

        t_1011[k] = f_8 * lsg0_543[k]
                    - f_9 * lsg1_543[k]
                    + f_3 * pc_x[k] * lsh_759[k];

        t_1012[k] = f_3 * pc_z[k] * lsh_757[k];

        t_1013[k] = f_8 * lsg0_545[k]
                    - f_9 * lsg1_545[k]
                    + f_3 * pc_x[k] * lsh_761[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, t_1017, pc_x, pc_z, lsg0_546, lsg0_548, \
                         lsg0_549, lsg1_546, lsg1_548, lsg1_549, lsh_759, lsh_762, lsh_764, \
                         lsh_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_6 * lsg0_546[k]
                    - f_7 * lsg1_546[k]
                    + f_3 * pc_x[k] * lsh_762[k];

        t_1015[k] = f_3 * pc_z[k] * lsh_759[k];

        t_1016[k] = f_6 * lsg0_548[k]
                    - f_7 * lsg1_548[k]
                    + f_3 * pc_x[k] * lsh_764[k];

        t_1017[k] = f_6 * lsg0_549[k]
                    - f_7 * lsg1_549[k]
                    + f_3 * pc_x[k] * lsh_765[k];
    }

#pragma omp simd aligned(t_1018, t_1019, t_1020, t_1021, pc_x, pc_z, lsg0_550, lsg0_552, \
                         lsg0_553, lsg1_550, lsg1_552, lsg1_553, lsh_762, lsh_766, lsh_768, \
                         lsh_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1018[k] = f_4 * lsg0_550[k]
                    - f_5 * lsg1_550[k]
                    + f_3 * pc_x[k] * lsh_766[k];

        t_1019[k] = f_3 * pc_z[k] * lsh_762[k];

        t_1020[k] = f_4 * lsg0_552[k]
                    - f_5 * lsg1_552[k]
                    + f_3 * pc_x[k] * lsh_768[k];

        t_1021[k] = f_4 * lsg0_553[k]
                    - f_5 * lsg1_553[k]
                    + f_3 * pc_x[k] * lsh_769[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, t_1025, t_1026, t_1027, pc_x, lsg0_554, \
                         lsg1_554, lsh_770, lsh_771, lsh_772, lsh_773, lsh_774, \
                         lsh_775 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_4 * lsg0_554[k]
                    - f_5 * lsg1_554[k]
                    + f_3 * pc_x[k] * lsh_770[k];

        t_1023[k] = f_3 * pc_x[k] * lsh_771[k];

        t_1024[k] = f_3 * pc_x[k] * lsh_772[k];

        t_1025[k] = f_3 * pc_x[k] * lsh_773[k];

        t_1026[k] = f_3 * pc_x[k] * lsh_774[k];

        t_1027[k] = f_3 * pc_x[k] * lsh_775[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, t_1031, pc_x, pc_y, pc_z, ksh_603, lsg0_550, \
                         lsg1_550, lsh_771, lsh_772, lsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_3 * pc_x[k] * lsh_776[k];

        t_1029[k] = f_0 * ksh_603[k]
                    + f_1 * lsg0_550[k]
                    - f_2 * lsg1_550[k]
                    + f_3 * pc_y[k] * lsh_771[k];

        t_1030[k] = f_3 * pc_z[k] * lsh_771[k];

        t_1031[k] = f_4 * lsg0_550[k]
                    - f_5 * lsg1_550[k]
                    + f_3 * pc_z[k] * lsh_772[k];
    }

#pragma omp simd aligned(t_1032, t_1033, t_1034, t_1035, pc_y, pc_z, ksh_608, lsg0_551, \
                         lsg0_552, lsg0_554, lsg1_551, lsg1_552, lsg1_554, lsh_773, lsh_774, \
                         lsh_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1032[k] = f_6 * lsg0_551[k]
                    - f_7 * lsg1_551[k]
                    + f_3 * pc_z[k] * lsh_773[k];

        t_1033[k] = f_8 * lsg0_552[k]
                    - f_9 * lsg1_552[k]
                    + f_3 * pc_z[k] * lsh_774[k];

        t_1034[k] = f_0 * ksh_608[k]
                    + f_3 * pc_y[k] * lsh_776[k];

        t_1035[k] = f_1 * lsg0_554[k]
                    - f_2 * lsg1_554[k]
                    + f_3 * pc_z[k] * lsh_776[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, t_1039, pa_z, pc_x, pc_z, ksi0_784, ksi0_785, \
                         ksi0_787, ksi1_784, ksi1_785, ksi1_787, lsg0_557, lsg1_557, \
                         lsh_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = pa_z[k] * ksi0_784[k]
                    - f_10 * pc_z[k] * ksi1_784[k];

        t_1037[k] = pa_z[k] * ksi0_785[k]
                    - f_10 * pc_z[k] * ksi1_785[k];

        t_1038[k] = f_16 * lsg0_557[k]
                    - f_17 * lsg1_557[k]
                    + f_3 * pc_x[k] * lsh_779[k];

        t_1039[k] = pa_z[k] * ksi0_787[k]
                    - f_10 * pc_z[k] * ksi1_787[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pa_z, pc_x, pc_z, ksi0_790, ksi1_790, \
                         lsg0_559, lsg0_560, lsg1_559, lsg1_560, lsh_781, \
                         lsh_782 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_8 * lsg0_559[k]
                    - f_9 * lsg1_559[k]
                    + f_3 * pc_x[k] * lsh_781[k];

        t_1041[k] = f_8 * lsg0_560[k]
                    - f_9 * lsg1_560[k]
                    + f_3 * pc_x[k] * lsh_782[k];

        t_1042[k] = pa_z[k] * ksi0_790[k]
                    - f_10 * pc_z[k] * ksi1_790[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, lsg0_562, lsg0_563, lsg0_564, lsg1_562, \
                         lsg1_563, lsg1_564, lsh_784, lsh_785, \
                         lsh_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_6 * lsg0_562[k]
                    - f_7 * lsg1_562[k]
                    + f_3 * pc_x[k] * lsh_784[k];

        t_1044[k] = f_6 * lsg0_563[k]
                    - f_7 * lsg1_563[k]
                    + f_3 * pc_x[k] * lsh_785[k];

        t_1045[k] = f_6 * lsg0_564[k]
                    - f_7 * lsg1_564[k]
                    + f_3 * pc_x[k] * lsh_786[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pa_z, pc_x, pc_z, ksi0_794, ksi1_794, \
                         lsg0_566, lsg0_567, lsg1_566, lsg1_567, lsh_788, \
                         lsh_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = pa_z[k] * ksi0_794[k]
                    - f_10 * pc_z[k] * ksi1_794[k];

        t_1047[k] = f_4 * lsg0_566[k]
                    - f_5 * lsg1_566[k]
                    + f_3 * pc_x[k] * lsh_788[k];

        t_1048[k] = f_4 * lsg0_567[k]
                    - f_5 * lsg1_567[k]
                    + f_3 * pc_x[k] * lsh_789[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, t_1052, t_1053, pc_x, lsg0_568, lsg0_569, \
                         lsg1_568, lsg1_569, lsh_790, lsh_791, lsh_792, lsh_793, \
                         lsh_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_4 * lsg0_568[k]
                    - f_5 * lsg1_568[k]
                    + f_3 * pc_x[k] * lsh_790[k];

        t_1050[k] = f_4 * lsg0_569[k]
                    - f_5 * lsg1_569[k]
                    + f_3 * pc_x[k] * lsh_791[k];

        t_1051[k] = f_3 * pc_x[k] * lsh_792[k];

        t_1052[k] = f_3 * pc_x[k] * lsh_793[k];

        t_1053[k] = f_3 * pc_x[k] * lsh_794[k];
    }

#pragma omp simd aligned(t_1054, t_1055, t_1056, t_1057, t_1058, pa_z, pc_x, pc_z, ksi0_805, \
                         ksh_603, ksi1_805, lsh_792, lsh_795, lsh_796, \
                         lsh_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1054[k] = f_3 * pc_x[k] * lsh_795[k];

        t_1055[k] = f_3 * pc_x[k] * lsh_796[k];

        t_1056[k] = f_3 * pc_x[k] * lsh_797[k];

        t_1057[k] = pa_z[k] * ksi0_805[k]
                    - f_10 * pc_z[k] * ksi1_805[k];

        t_1058[k] = f_11 * ksh_603[k]
                    + f_3 * pc_z[k] * lsh_792[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pa_z, pc_z, ksi0_807, ksi0_808, ksi0_809, \
                         ksh_604, ksh_605, ksh_606, ksi1_807, ksi1_808, \
                         ksi1_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = pa_z[k] * ksi0_807[k]
                    + f_12 * ksh_604[k]
                    - f_10 * pc_z[k] * ksi1_807[k];

        t_1060[k] = pa_z[k] * ksi0_808[k]
                    + f_13 * ksh_605[k]
                    - f_10 * pc_z[k] * ksi1_808[k];

        t_1061[k] = pa_z[k] * ksi0_809[k]
                    + f_14 * ksh_606[k]
                    - f_10 * pc_z[k] * ksi1_809[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, pc_x, pc_y, pc_z, ksh_608, ksh_629, lsg0_569, \
                         lsg0_570, lsg1_569, lsg1_570, lsh_797, \
                         lsh_798 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_15 * ksh_629[k]
                    + f_3 * pc_y[k] * lsh_797[k];

        t_1063[k] = f_11 * ksh_608[k]
                    + f_1 * lsg0_569[k]
                    - f_2 * lsg1_569[k]
                    + f_3 * pc_z[k] * lsh_797[k];

        t_1064[k] = f_1 * lsg0_570[k]
                    - f_2 * lsg1_570[k]
                    + f_3 * pc_x[k] * lsh_798[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, pc_x, lsg0_571, lsg0_572, lsg0_573, lsg1_571, \
                         lsg1_572, lsg1_573, lsh_799, lsh_800, \
                         lsh_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_16 * lsg0_571[k]
                    - f_17 * lsg1_571[k]
                    + f_3 * pc_x[k] * lsh_799[k];

        t_1066[k] = f_16 * lsg0_572[k]
                    - f_17 * lsg1_572[k]
                    + f_3 * pc_x[k] * lsh_800[k];

        t_1067[k] = f_8 * lsg0_573[k]
                    - f_9 * lsg1_573[k]
                    + f_3 * pc_x[k] * lsh_801[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t ksh, const size_t lsg0,
                                                          const size_t lsg1, const size_t lsh,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);
    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);
    auto *t_1170 = buffer.data(target + 1170);
    auto *t_1171 = buffer.data(target + 1171);
    auto *t_1172 = buffer.data(target + 1172);
    auto *t_1173 = buffer.data(target + 1173);
    auto *t_1174 = buffer.data(target + 1174);
    auto *t_1175 = buffer.data(target + 1175);
    auto *t_1176 = buffer.data(target + 1176);
    auto *t_1177 = buffer.data(target + 1177);
    auto *t_1178 = buffer.data(target + 1178);
    auto *t_1179 = buffer.data(target + 1179);
    auto *t_1180 = buffer.data(target + 1180);
    auto *t_1181 = buffer.data(target + 1181);
    auto *t_1182 = buffer.data(target + 1182);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksh_624 = buffer.data(ksh + 624);
    const auto *ksh_629 = buffer.data(ksh + 629);
    const auto *ksh_645 = buffer.data(ksh + 645);
    const auto *ksh_647 = buffer.data(ksh + 647);
    const auto *ksh_648 = buffer.data(ksh + 648);
    const auto *ksh_649 = buffer.data(ksh + 649);
    const auto *ksh_650 = buffer.data(ksh + 650);
    const auto *ksh_666 = buffer.data(ksh + 666);
    const auto *ksh_668 = buffer.data(ksh + 668);
    const auto *ksh_669 = buffer.data(ksh + 669);
    const auto *ksh_670 = buffer.data(ksh + 670);
    const auto *ksh_671 = buffer.data(ksh + 671);
    const auto *ksh_687 = buffer.data(ksh + 687);
    const auto *ksh_689 = buffer.data(ksh + 689);
    const auto *ksh_690 = buffer.data(ksh + 690);
    const auto *ksh_691 = buffer.data(ksh + 691);
    const auto *ksh_692 = buffer.data(ksh + 692);
    const auto *ksh_708 = buffer.data(ksh + 708);
    const auto *ksh_710 = buffer.data(ksh + 710);
    const auto *ksh_711 = buffer.data(ksh + 711);
    const auto *ksh_712 = buffer.data(ksh + 712);
    const auto *ksh_713 = buffer.data(ksh + 713);

    const auto *lsg0_574 = buffer.data(lsg0 + 574);
    const auto *lsg0_575 = buffer.data(lsg0 + 575);
    const auto *lsg0_576 = buffer.data(lsg0 + 576);
    const auto *lsg0_577 = buffer.data(lsg0 + 577);
    const auto *lsg0_578 = buffer.data(lsg0 + 578);
    const auto *lsg0_579 = buffer.data(lsg0 + 579);
    const auto *lsg0_580 = buffer.data(lsg0 + 580);
    const auto *lsg0_581 = buffer.data(lsg0 + 581);
    const auto *lsg0_582 = buffer.data(lsg0 + 582);
    const auto *lsg0_583 = buffer.data(lsg0 + 583);
    const auto *lsg0_584 = buffer.data(lsg0 + 584);
    const auto *lsg0_585 = buffer.data(lsg0 + 585);
    const auto *lsg0_586 = buffer.data(lsg0 + 586);
    const auto *lsg0_587 = buffer.data(lsg0 + 587);
    const auto *lsg0_588 = buffer.data(lsg0 + 588);
    const auto *lsg0_589 = buffer.data(lsg0 + 589);
    const auto *lsg0_590 = buffer.data(lsg0 + 590);
    const auto *lsg0_591 = buffer.data(lsg0 + 591);
    const auto *lsg0_592 = buffer.data(lsg0 + 592);
    const auto *lsg0_593 = buffer.data(lsg0 + 593);
    const auto *lsg0_594 = buffer.data(lsg0 + 594);
    const auto *lsg0_595 = buffer.data(lsg0 + 595);
    const auto *lsg0_596 = buffer.data(lsg0 + 596);
    const auto *lsg0_597 = buffer.data(lsg0 + 597);
    const auto *lsg0_598 = buffer.data(lsg0 + 598);
    const auto *lsg0_599 = buffer.data(lsg0 + 599);
    const auto *lsg0_600 = buffer.data(lsg0 + 600);
    const auto *lsg0_601 = buffer.data(lsg0 + 601);
    const auto *lsg0_602 = buffer.data(lsg0 + 602);
    const auto *lsg0_603 = buffer.data(lsg0 + 603);
    const auto *lsg0_604 = buffer.data(lsg0 + 604);
    const auto *lsg0_605 = buffer.data(lsg0 + 605);
    const auto *lsg0_606 = buffer.data(lsg0 + 606);
    const auto *lsg0_607 = buffer.data(lsg0 + 607);
    const auto *lsg0_608 = buffer.data(lsg0 + 608);
    const auto *lsg0_609 = buffer.data(lsg0 + 609);
    const auto *lsg0_610 = buffer.data(lsg0 + 610);
    const auto *lsg0_611 = buffer.data(lsg0 + 611);
    const auto *lsg0_612 = buffer.data(lsg0 + 612);
    const auto *lsg0_613 = buffer.data(lsg0 + 613);
    const auto *lsg0_614 = buffer.data(lsg0 + 614);
    const auto *lsg0_615 = buffer.data(lsg0 + 615);
    const auto *lsg0_616 = buffer.data(lsg0 + 616);
    const auto *lsg0_617 = buffer.data(lsg0 + 617);
    const auto *lsg0_618 = buffer.data(lsg0 + 618);
    const auto *lsg0_619 = buffer.data(lsg0 + 619);
    const auto *lsg0_620 = buffer.data(lsg0 + 620);
    const auto *lsg0_621 = buffer.data(lsg0 + 621);
    const auto *lsg0_622 = buffer.data(lsg0 + 622);
    const auto *lsg0_623 = buffer.data(lsg0 + 623);
    const auto *lsg0_624 = buffer.data(lsg0 + 624);
    const auto *lsg0_625 = buffer.data(lsg0 + 625);
    const auto *lsg0_626 = buffer.data(lsg0 + 626);
    const auto *lsg0_627 = buffer.data(lsg0 + 627);
    const auto *lsg0_628 = buffer.data(lsg0 + 628);
    const auto *lsg0_629 = buffer.data(lsg0 + 629);
    const auto *lsg0_630 = buffer.data(lsg0 + 630);
    const auto *lsg0_631 = buffer.data(lsg0 + 631);
    const auto *lsg0_632 = buffer.data(lsg0 + 632);
    const auto *lsg0_633 = buffer.data(lsg0 + 633);
    const auto *lsg0_634 = buffer.data(lsg0 + 634);
    const auto *lsg0_635 = buffer.data(lsg0 + 635);
    const auto *lsg0_636 = buffer.data(lsg0 + 636);

    const auto *lsg1_574 = buffer.data(lsg1 + 574);
    const auto *lsg1_575 = buffer.data(lsg1 + 575);
    const auto *lsg1_576 = buffer.data(lsg1 + 576);
    const auto *lsg1_577 = buffer.data(lsg1 + 577);
    const auto *lsg1_578 = buffer.data(lsg1 + 578);
    const auto *lsg1_579 = buffer.data(lsg1 + 579);
    const auto *lsg1_580 = buffer.data(lsg1 + 580);
    const auto *lsg1_581 = buffer.data(lsg1 + 581);
    const auto *lsg1_582 = buffer.data(lsg1 + 582);
    const auto *lsg1_583 = buffer.data(lsg1 + 583);
    const auto *lsg1_584 = buffer.data(lsg1 + 584);
    const auto *lsg1_585 = buffer.data(lsg1 + 585);
    const auto *lsg1_586 = buffer.data(lsg1 + 586);
    const auto *lsg1_587 = buffer.data(lsg1 + 587);
    const auto *lsg1_588 = buffer.data(lsg1 + 588);
    const auto *lsg1_589 = buffer.data(lsg1 + 589);
    const auto *lsg1_590 = buffer.data(lsg1 + 590);
    const auto *lsg1_591 = buffer.data(lsg1 + 591);
    const auto *lsg1_592 = buffer.data(lsg1 + 592);
    const auto *lsg1_593 = buffer.data(lsg1 + 593);
    const auto *lsg1_594 = buffer.data(lsg1 + 594);
    const auto *lsg1_595 = buffer.data(lsg1 + 595);
    const auto *lsg1_596 = buffer.data(lsg1 + 596);
    const auto *lsg1_597 = buffer.data(lsg1 + 597);
    const auto *lsg1_598 = buffer.data(lsg1 + 598);
    const auto *lsg1_599 = buffer.data(lsg1 + 599);
    const auto *lsg1_600 = buffer.data(lsg1 + 600);
    const auto *lsg1_601 = buffer.data(lsg1 + 601);
    const auto *lsg1_602 = buffer.data(lsg1 + 602);
    const auto *lsg1_603 = buffer.data(lsg1 + 603);
    const auto *lsg1_604 = buffer.data(lsg1 + 604);
    const auto *lsg1_605 = buffer.data(lsg1 + 605);
    const auto *lsg1_606 = buffer.data(lsg1 + 606);
    const auto *lsg1_607 = buffer.data(lsg1 + 607);
    const auto *lsg1_608 = buffer.data(lsg1 + 608);
    const auto *lsg1_609 = buffer.data(lsg1 + 609);
    const auto *lsg1_610 = buffer.data(lsg1 + 610);
    const auto *lsg1_611 = buffer.data(lsg1 + 611);
    const auto *lsg1_612 = buffer.data(lsg1 + 612);
    const auto *lsg1_613 = buffer.data(lsg1 + 613);
    const auto *lsg1_614 = buffer.data(lsg1 + 614);
    const auto *lsg1_615 = buffer.data(lsg1 + 615);
    const auto *lsg1_616 = buffer.data(lsg1 + 616);
    const auto *lsg1_617 = buffer.data(lsg1 + 617);
    const auto *lsg1_618 = buffer.data(lsg1 + 618);
    const auto *lsg1_619 = buffer.data(lsg1 + 619);
    const auto *lsg1_620 = buffer.data(lsg1 + 620);
    const auto *lsg1_621 = buffer.data(lsg1 + 621);
    const auto *lsg1_622 = buffer.data(lsg1 + 622);
    const auto *lsg1_623 = buffer.data(lsg1 + 623);
    const auto *lsg1_624 = buffer.data(lsg1 + 624);
    const auto *lsg1_625 = buffer.data(lsg1 + 625);
    const auto *lsg1_626 = buffer.data(lsg1 + 626);
    const auto *lsg1_627 = buffer.data(lsg1 + 627);
    const auto *lsg1_628 = buffer.data(lsg1 + 628);
    const auto *lsg1_629 = buffer.data(lsg1 + 629);
    const auto *lsg1_630 = buffer.data(lsg1 + 630);
    const auto *lsg1_631 = buffer.data(lsg1 + 631);
    const auto *lsg1_632 = buffer.data(lsg1 + 632);
    const auto *lsg1_633 = buffer.data(lsg1 + 633);
    const auto *lsg1_634 = buffer.data(lsg1 + 634);
    const auto *lsg1_635 = buffer.data(lsg1 + 635);
    const auto *lsg1_636 = buffer.data(lsg1 + 636);

    const auto *lsh_802 = buffer.data(lsh + 802);
    const auto *lsh_803 = buffer.data(lsh + 803);
    const auto *lsh_804 = buffer.data(lsh + 804);
    const auto *lsh_805 = buffer.data(lsh + 805);
    const auto *lsh_806 = buffer.data(lsh + 806);
    const auto *lsh_807 = buffer.data(lsh + 807);
    const auto *lsh_808 = buffer.data(lsh + 808);
    const auto *lsh_809 = buffer.data(lsh + 809);
    const auto *lsh_810 = buffer.data(lsh + 810);
    const auto *lsh_811 = buffer.data(lsh + 811);
    const auto *lsh_812 = buffer.data(lsh + 812);
    const auto *lsh_813 = buffer.data(lsh + 813);
    const auto *lsh_814 = buffer.data(lsh + 814);
    const auto *lsh_815 = buffer.data(lsh + 815);
    const auto *lsh_816 = buffer.data(lsh + 816);
    const auto *lsh_817 = buffer.data(lsh + 817);
    const auto *lsh_818 = buffer.data(lsh + 818);
    const auto *lsh_819 = buffer.data(lsh + 819);
    const auto *lsh_820 = buffer.data(lsh + 820);
    const auto *lsh_821 = buffer.data(lsh + 821);
    const auto *lsh_822 = buffer.data(lsh + 822);
    const auto *lsh_823 = buffer.data(lsh + 823);
    const auto *lsh_824 = buffer.data(lsh + 824);
    const auto *lsh_825 = buffer.data(lsh + 825);
    const auto *lsh_826 = buffer.data(lsh + 826);
    const auto *lsh_827 = buffer.data(lsh + 827);
    const auto *lsh_828 = buffer.data(lsh + 828);
    const auto *lsh_829 = buffer.data(lsh + 829);
    const auto *lsh_830 = buffer.data(lsh + 830);
    const auto *lsh_831 = buffer.data(lsh + 831);
    const auto *lsh_832 = buffer.data(lsh + 832);
    const auto *lsh_833 = buffer.data(lsh + 833);
    const auto *lsh_834 = buffer.data(lsh + 834);
    const auto *lsh_835 = buffer.data(lsh + 835);
    const auto *lsh_836 = buffer.data(lsh + 836);
    const auto *lsh_837 = buffer.data(lsh + 837);
    const auto *lsh_838 = buffer.data(lsh + 838);
    const auto *lsh_839 = buffer.data(lsh + 839);
    const auto *lsh_840 = buffer.data(lsh + 840);
    const auto *lsh_841 = buffer.data(lsh + 841);
    const auto *lsh_842 = buffer.data(lsh + 842);
    const auto *lsh_843 = buffer.data(lsh + 843);
    const auto *lsh_844 = buffer.data(lsh + 844);
    const auto *lsh_845 = buffer.data(lsh + 845);
    const auto *lsh_846 = buffer.data(lsh + 846);
    const auto *lsh_847 = buffer.data(lsh + 847);
    const auto *lsh_848 = buffer.data(lsh + 848);
    const auto *lsh_849 = buffer.data(lsh + 849);
    const auto *lsh_850 = buffer.data(lsh + 850);
    const auto *lsh_851 = buffer.data(lsh + 851);
    const auto *lsh_852 = buffer.data(lsh + 852);
    const auto *lsh_853 = buffer.data(lsh + 853);
    const auto *lsh_854 = buffer.data(lsh + 854);
    const auto *lsh_855 = buffer.data(lsh + 855);
    const auto *lsh_856 = buffer.data(lsh + 856);
    const auto *lsh_857 = buffer.data(lsh + 857);
    const auto *lsh_858 = buffer.data(lsh + 858);
    const auto *lsh_859 = buffer.data(lsh + 859);
    const auto *lsh_860 = buffer.data(lsh + 860);
    const auto *lsh_861 = buffer.data(lsh + 861);
    const auto *lsh_862 = buffer.data(lsh + 862);
    const auto *lsh_863 = buffer.data(lsh + 863);
    const auto *lsh_864 = buffer.data(lsh + 864);
    const auto *lsh_865 = buffer.data(lsh + 865);
    const auto *lsh_866 = buffer.data(lsh + 866);
    const auto *lsh_867 = buffer.data(lsh + 867);
    const auto *lsh_868 = buffer.data(lsh + 868);
    const auto *lsh_869 = buffer.data(lsh + 869);
    const auto *lsh_870 = buffer.data(lsh + 870);
    const auto *lsh_871 = buffer.data(lsh + 871);
    const auto *lsh_872 = buffer.data(lsh + 872);
    const auto *lsh_873 = buffer.data(lsh + 873);
    const auto *lsh_874 = buffer.data(lsh + 874);
    const auto *lsh_875 = buffer.data(lsh + 875);
    const auto *lsh_876 = buffer.data(lsh + 876);
    const auto *lsh_877 = buffer.data(lsh + 877);
    const auto *lsh_878 = buffer.data(lsh + 878);
    const auto *lsh_879 = buffer.data(lsh + 879);
    const auto *lsh_880 = buffer.data(lsh + 880);
    const auto *lsh_881 = buffer.data(lsh + 881);
    const auto *lsh_882 = buffer.data(lsh + 882);
    const auto *lsh_883 = buffer.data(lsh + 883);
    const auto *lsh_884 = buffer.data(lsh + 884);
    const auto *lsh_885 = buffer.data(lsh + 885);
    const auto *lsh_886 = buffer.data(lsh + 886);
    const auto *lsh_887 = buffer.data(lsh + 887);
    const auto *lsh_888 = buffer.data(lsh + 888);

#pragma omp simd aligned(t_1068, t_1069, t_1070, pc_x, lsg0_574, lsg0_575, lsg0_576, lsg1_574, \
                         lsg1_575, lsg1_576, lsh_802, lsh_803, \
                         lsh_804 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_8 * lsg0_574[k]
                    - f_9 * lsg1_574[k]
                    + f_3 * pc_x[k] * lsh_802[k];

        t_1069[k] = f_8 * lsg0_575[k]
                    - f_9 * lsg1_575[k]
                    + f_3 * pc_x[k] * lsh_803[k];

        t_1070[k] = f_6 * lsg0_576[k]
                    - f_7 * lsg1_576[k]
                    + f_3 * pc_x[k] * lsh_804[k];
    }

#pragma omp simd aligned(t_1071, t_1072, t_1073, pc_x, lsg0_577, lsg0_578, lsg0_579, lsg1_577, \
                         lsg1_578, lsg1_579, lsh_805, lsh_806, \
                         lsh_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1071[k] = f_6 * lsg0_577[k]
                    - f_7 * lsg1_577[k]
                    + f_3 * pc_x[k] * lsh_805[k];

        t_1072[k] = f_6 * lsg0_578[k]
                    - f_7 * lsg1_578[k]
                    + f_3 * pc_x[k] * lsh_806[k];

        t_1073[k] = f_6 * lsg0_579[k]
                    - f_7 * lsg1_579[k]
                    + f_3 * pc_x[k] * lsh_807[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, pc_x, lsg0_580, lsg0_581, lsg0_582, lsg1_580, \
                         lsg1_581, lsg1_582, lsh_808, lsh_809, \
                         lsh_810 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = f_4 * lsg0_580[k]
                    - f_5 * lsg1_580[k]
                    + f_3 * pc_x[k] * lsh_808[k];

        t_1075[k] = f_4 * lsg0_581[k]
                    - f_5 * lsg1_581[k]
                    + f_3 * pc_x[k] * lsh_809[k];

        t_1076[k] = f_4 * lsg0_582[k]
                    - f_5 * lsg1_582[k]
                    + f_3 * pc_x[k] * lsh_810[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, t_1080, t_1081, pc_x, lsg0_583, lsg0_584, \
                         lsg1_583, lsg1_584, lsh_811, lsh_812, lsh_813, lsh_814, \
                         lsh_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_4 * lsg0_583[k]
                    - f_5 * lsg1_583[k]
                    + f_3 * pc_x[k] * lsh_811[k];

        t_1078[k] = f_4 * lsg0_584[k]
                    - f_5 * lsg1_584[k]
                    + f_3 * pc_x[k] * lsh_812[k];

        t_1079[k] = f_3 * pc_x[k] * lsh_813[k];

        t_1080[k] = f_3 * pc_x[k] * lsh_814[k];

        t_1081[k] = f_3 * pc_x[k] * lsh_815[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, t_1085, t_1086, pc_x, pc_y, pc_z, ksh_624, \
                         ksh_645, lsg0_580, lsg1_580, lsh_813, lsh_816, lsh_817, \
                         lsh_818 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_3 * pc_x[k] * lsh_816[k];

        t_1083[k] = f_3 * pc_x[k] * lsh_817[k];

        t_1084[k] = f_3 * pc_x[k] * lsh_818[k];

        t_1085[k] = f_18 * ksh_645[k]
                    + f_1 * lsg0_580[k]
                    - f_2 * lsg1_580[k]
                    + f_3 * pc_y[k] * lsh_813[k];

        t_1086[k] = f_12 * ksh_624[k]
                    + f_3 * pc_z[k] * lsh_813[k];
    }

#pragma omp simd aligned(t_1087, t_1088, t_1089, pc_y, ksh_647, ksh_648, ksh_649, lsg0_582, \
                         lsg0_583, lsg0_584, lsg1_582, lsg1_583, lsg1_584, lsh_815, lsh_816, \
                         lsh_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_18 * ksh_647[k]
                    + f_8 * lsg0_582[k]
                    - f_9 * lsg1_582[k]
                    + f_3 * pc_y[k] * lsh_815[k];

        t_1088[k] = f_18 * ksh_648[k]
                    + f_6 * lsg0_583[k]
                    - f_7 * lsg1_583[k]
                    + f_3 * pc_y[k] * lsh_816[k];

        t_1089[k] = f_18 * ksh_649[k]
                    + f_4 * lsg0_584[k]
                    - f_5 * lsg1_584[k]
                    + f_3 * pc_y[k] * lsh_817[k];
    }

#pragma omp simd aligned(t_1090, t_1091, t_1092, pc_x, pc_y, pc_z, ksh_629, ksh_650, lsg0_584, \
                         lsg0_585, lsg1_584, lsg1_585, lsh_818, \
                         lsh_819 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1090[k] = f_18 * ksh_650[k]
                    + f_3 * pc_y[k] * lsh_818[k];

        t_1091[k] = f_12 * ksh_629[k]
                    + f_1 * lsg0_584[k]
                    - f_2 * lsg1_584[k]
                    + f_3 * pc_z[k] * lsh_818[k];

        t_1092[k] = f_1 * lsg0_585[k]
                    - f_2 * lsg1_585[k]
                    + f_3 * pc_x[k] * lsh_819[k];
    }

#pragma omp simd aligned(t_1093, t_1094, t_1095, pc_x, lsg0_586, lsg0_587, lsg0_588, lsg1_586, \
                         lsg1_587, lsg1_588, lsh_820, lsh_821, \
                         lsh_822 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1093[k] = f_16 * lsg0_586[k]
                    - f_17 * lsg1_586[k]
                    + f_3 * pc_x[k] * lsh_820[k];

        t_1094[k] = f_16 * lsg0_587[k]
                    - f_17 * lsg1_587[k]
                    + f_3 * pc_x[k] * lsh_821[k];

        t_1095[k] = f_8 * lsg0_588[k]
                    - f_9 * lsg1_588[k]
                    + f_3 * pc_x[k] * lsh_822[k];
    }

#pragma omp simd aligned(t_1096, t_1097, t_1098, pc_x, lsg0_589, lsg0_590, lsg0_591, lsg1_589, \
                         lsg1_590, lsg1_591, lsh_823, lsh_824, \
                         lsh_825 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1096[k] = f_8 * lsg0_589[k]
                    - f_9 * lsg1_589[k]
                    + f_3 * pc_x[k] * lsh_823[k];

        t_1097[k] = f_8 * lsg0_590[k]
                    - f_9 * lsg1_590[k]
                    + f_3 * pc_x[k] * lsh_824[k];

        t_1098[k] = f_6 * lsg0_591[k]
                    - f_7 * lsg1_591[k]
                    + f_3 * pc_x[k] * lsh_825[k];
    }

#pragma omp simd aligned(t_1099, t_1100, t_1101, pc_x, lsg0_592, lsg0_593, lsg0_594, lsg1_592, \
                         lsg1_593, lsg1_594, lsh_826, lsh_827, \
                         lsh_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1099[k] = f_6 * lsg0_592[k]
                    - f_7 * lsg1_592[k]
                    + f_3 * pc_x[k] * lsh_826[k];

        t_1100[k] = f_6 * lsg0_593[k]
                    - f_7 * lsg1_593[k]
                    + f_3 * pc_x[k] * lsh_827[k];

        t_1101[k] = f_6 * lsg0_594[k]
                    - f_7 * lsg1_594[k]
                    + f_3 * pc_x[k] * lsh_828[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, pc_x, lsg0_595, lsg0_596, lsg0_597, lsg1_595, \
                         lsg1_596, lsg1_597, lsh_829, lsh_830, \
                         lsh_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = f_4 * lsg0_595[k]
                    - f_5 * lsg1_595[k]
                    + f_3 * pc_x[k] * lsh_829[k];

        t_1103[k] = f_4 * lsg0_596[k]
                    - f_5 * lsg1_596[k]
                    + f_3 * pc_x[k] * lsh_830[k];

        t_1104[k] = f_4 * lsg0_597[k]
                    - f_5 * lsg1_597[k]
                    + f_3 * pc_x[k] * lsh_831[k];
    }

#pragma omp simd aligned(t_1105, t_1106, t_1107, t_1108, t_1109, pc_x, lsg0_598, lsg0_599, \
                         lsg1_598, lsg1_599, lsh_832, lsh_833, lsh_834, lsh_835, \
                         lsh_836 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1105[k] = f_4 * lsg0_598[k]
                    - f_5 * lsg1_598[k]
                    + f_3 * pc_x[k] * lsh_832[k];

        t_1106[k] = f_4 * lsg0_599[k]
                    - f_5 * lsg1_599[k]
                    + f_3 * pc_x[k] * lsh_833[k];

        t_1107[k] = f_3 * pc_x[k] * lsh_834[k];

        t_1108[k] = f_3 * pc_x[k] * lsh_835[k];

        t_1109[k] = f_3 * pc_x[k] * lsh_836[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, pc_y, pc_z, ksh_645, \
                         ksh_666, lsg0_595, lsg1_595, lsh_834, lsh_837, lsh_838, \
                         lsh_839 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_3 * pc_x[k] * lsh_837[k];

        t_1111[k] = f_3 * pc_x[k] * lsh_838[k];

        t_1112[k] = f_3 * pc_x[k] * lsh_839[k];

        t_1113[k] = f_19 * ksh_666[k]
                    + f_1 * lsg0_595[k]
                    - f_2 * lsg1_595[k]
                    + f_3 * pc_y[k] * lsh_834[k];

        t_1114[k] = f_13 * ksh_645[k]
                    + f_3 * pc_z[k] * lsh_834[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_y, ksh_668, ksh_669, ksh_670, lsg0_597, \
                         lsg0_598, lsg0_599, lsg1_597, lsg1_598, lsg1_599, lsh_836, lsh_837, \
                         lsh_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_19 * ksh_668[k]
                    + f_8 * lsg0_597[k]
                    - f_9 * lsg1_597[k]
                    + f_3 * pc_y[k] * lsh_836[k];

        t_1116[k] = f_19 * ksh_669[k]
                    + f_6 * lsg0_598[k]
                    - f_7 * lsg1_598[k]
                    + f_3 * pc_y[k] * lsh_837[k];

        t_1117[k] = f_19 * ksh_670[k]
                    + f_4 * lsg0_599[k]
                    - f_5 * lsg1_599[k]
                    + f_3 * pc_y[k] * lsh_838[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_x, pc_y, pc_z, ksh_650, ksh_671, lsg0_599, \
                         lsg0_600, lsg1_599, lsg1_600, lsh_839, \
                         lsh_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_19 * ksh_671[k]
                    + f_3 * pc_y[k] * lsh_839[k];

        t_1119[k] = f_13 * ksh_650[k]
                    + f_1 * lsg0_599[k]
                    - f_2 * lsg1_599[k]
                    + f_3 * pc_z[k] * lsh_839[k];

        t_1120[k] = f_1 * lsg0_600[k]
                    - f_2 * lsg1_600[k]
                    + f_3 * pc_x[k] * lsh_840[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_x, lsg0_601, lsg0_602, lsg0_603, lsg1_601, \
                         lsg1_602, lsg1_603, lsh_841, lsh_842, \
                         lsh_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_16 * lsg0_601[k]
                    - f_17 * lsg1_601[k]
                    + f_3 * pc_x[k] * lsh_841[k];

        t_1122[k] = f_16 * lsg0_602[k]
                    - f_17 * lsg1_602[k]
                    + f_3 * pc_x[k] * lsh_842[k];

        t_1123[k] = f_8 * lsg0_603[k]
                    - f_9 * lsg1_603[k]
                    + f_3 * pc_x[k] * lsh_843[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, lsg0_604, lsg0_605, lsg0_606, lsg1_604, \
                         lsg1_605, lsg1_606, lsh_844, lsh_845, \
                         lsh_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_8 * lsg0_604[k]
                    - f_9 * lsg1_604[k]
                    + f_3 * pc_x[k] * lsh_844[k];

        t_1125[k] = f_8 * lsg0_605[k]
                    - f_9 * lsg1_605[k]
                    + f_3 * pc_x[k] * lsh_845[k];

        t_1126[k] = f_6 * lsg0_606[k]
                    - f_7 * lsg1_606[k]
                    + f_3 * pc_x[k] * lsh_846[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, lsg0_607, lsg0_608, lsg0_609, lsg1_607, \
                         lsg1_608, lsg1_609, lsh_847, lsh_848, \
                         lsh_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_6 * lsg0_607[k]
                    - f_7 * lsg1_607[k]
                    + f_3 * pc_x[k] * lsh_847[k];

        t_1128[k] = f_6 * lsg0_608[k]
                    - f_7 * lsg1_608[k]
                    + f_3 * pc_x[k] * lsh_848[k];

        t_1129[k] = f_6 * lsg0_609[k]
                    - f_7 * lsg1_609[k]
                    + f_3 * pc_x[k] * lsh_849[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, lsg0_610, lsg0_611, lsg0_612, lsg1_610, \
                         lsg1_611, lsg1_612, lsh_850, lsh_851, \
                         lsh_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_4 * lsg0_610[k]
                    - f_5 * lsg1_610[k]
                    + f_3 * pc_x[k] * lsh_850[k];

        t_1131[k] = f_4 * lsg0_611[k]
                    - f_5 * lsg1_611[k]
                    + f_3 * pc_x[k] * lsh_851[k];

        t_1132[k] = f_4 * lsg0_612[k]
                    - f_5 * lsg1_612[k]
                    + f_3 * pc_x[k] * lsh_852[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, t_1136, t_1137, pc_x, lsg0_613, lsg0_614, \
                         lsg1_613, lsg1_614, lsh_853, lsh_854, lsh_855, lsh_856, \
                         lsh_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_4 * lsg0_613[k]
                    - f_5 * lsg1_613[k]
                    + f_3 * pc_x[k] * lsh_853[k];

        t_1134[k] = f_4 * lsg0_614[k]
                    - f_5 * lsg1_614[k]
                    + f_3 * pc_x[k] * lsh_854[k];

        t_1135[k] = f_3 * pc_x[k] * lsh_855[k];

        t_1136[k] = f_3 * pc_x[k] * lsh_856[k];

        t_1137[k] = f_3 * pc_x[k] * lsh_857[k];
    }

#pragma omp simd aligned(t_1138, t_1139, t_1140, t_1141, t_1142, pc_x, pc_y, pc_z, ksh_666, \
                         ksh_687, lsg0_610, lsg1_610, lsh_855, lsh_858, lsh_859, \
                         lsh_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1138[k] = f_3 * pc_x[k] * lsh_858[k];

        t_1139[k] = f_3 * pc_x[k] * lsh_859[k];

        t_1140[k] = f_3 * pc_x[k] * lsh_860[k];

        t_1141[k] = f_14 * ksh_687[k]
                    + f_1 * lsg0_610[k]
                    - f_2 * lsg1_610[k]
                    + f_3 * pc_y[k] * lsh_855[k];

        t_1142[k] = f_14 * ksh_666[k]
                    + f_3 * pc_z[k] * lsh_855[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_y, ksh_689, ksh_690, ksh_691, lsg0_612, \
                         lsg0_613, lsg0_614, lsg1_612, lsg1_613, lsg1_614, lsh_857, lsh_858, \
                         lsh_859 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_14 * ksh_689[k]
                    + f_8 * lsg0_612[k]
                    - f_9 * lsg1_612[k]
                    + f_3 * pc_y[k] * lsh_857[k];

        t_1144[k] = f_14 * ksh_690[k]
                    + f_6 * lsg0_613[k]
                    - f_7 * lsg1_613[k]
                    + f_3 * pc_y[k] * lsh_858[k];

        t_1145[k] = f_14 * ksh_691[k]
                    + f_4 * lsg0_614[k]
                    - f_5 * lsg1_614[k]
                    + f_3 * pc_y[k] * lsh_859[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_y, pc_z, ksh_671, ksh_692, lsg0_614, \
                         lsg0_615, lsg1_614, lsg1_615, lsh_860, \
                         lsh_861 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_14 * ksh_692[k]
                    + f_3 * pc_y[k] * lsh_860[k];

        t_1147[k] = f_14 * ksh_671[k]
                    + f_1 * lsg0_614[k]
                    - f_2 * lsg1_614[k]
                    + f_3 * pc_z[k] * lsh_860[k];

        t_1148[k] = f_1 * lsg0_615[k]
                    - f_2 * lsg1_615[k]
                    + f_3 * pc_x[k] * lsh_861[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, pc_x, lsg0_616, lsg0_617, lsg0_618, lsg1_616, \
                         lsg1_617, lsg1_618, lsh_862, lsh_863, \
                         lsh_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_16 * lsg0_616[k]
                    - f_17 * lsg1_616[k]
                    + f_3 * pc_x[k] * lsh_862[k];

        t_1150[k] = f_16 * lsg0_617[k]
                    - f_17 * lsg1_617[k]
                    + f_3 * pc_x[k] * lsh_863[k];

        t_1151[k] = f_8 * lsg0_618[k]
                    - f_9 * lsg1_618[k]
                    + f_3 * pc_x[k] * lsh_864[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, pc_x, lsg0_619, lsg0_620, lsg0_621, lsg1_619, \
                         lsg1_620, lsg1_621, lsh_865, lsh_866, \
                         lsh_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = f_8 * lsg0_619[k]
                    - f_9 * lsg1_619[k]
                    + f_3 * pc_x[k] * lsh_865[k];

        t_1153[k] = f_8 * lsg0_620[k]
                    - f_9 * lsg1_620[k]
                    + f_3 * pc_x[k] * lsh_866[k];

        t_1154[k] = f_6 * lsg0_621[k]
                    - f_7 * lsg1_621[k]
                    + f_3 * pc_x[k] * lsh_867[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, pc_x, lsg0_622, lsg0_623, lsg0_624, lsg1_622, \
                         lsg1_623, lsg1_624, lsh_868, lsh_869, \
                         lsh_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_6 * lsg0_622[k]
                    - f_7 * lsg1_622[k]
                    + f_3 * pc_x[k] * lsh_868[k];

        t_1156[k] = f_6 * lsg0_623[k]
                    - f_7 * lsg1_623[k]
                    + f_3 * pc_x[k] * lsh_869[k];

        t_1157[k] = f_6 * lsg0_624[k]
                    - f_7 * lsg1_624[k]
                    + f_3 * pc_x[k] * lsh_870[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, pc_x, lsg0_625, lsg0_626, lsg0_627, lsg1_625, \
                         lsg1_626, lsg1_627, lsh_871, lsh_872, \
                         lsh_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = f_4 * lsg0_625[k]
                    - f_5 * lsg1_625[k]
                    + f_3 * pc_x[k] * lsh_871[k];

        t_1159[k] = f_4 * lsg0_626[k]
                    - f_5 * lsg1_626[k]
                    + f_3 * pc_x[k] * lsh_872[k];

        t_1160[k] = f_4 * lsg0_627[k]
                    - f_5 * lsg1_627[k]
                    + f_3 * pc_x[k] * lsh_873[k];
    }

#pragma omp simd aligned(t_1161, t_1162, t_1163, t_1164, t_1165, pc_x, lsg0_628, lsg0_629, \
                         lsg1_628, lsg1_629, lsh_874, lsh_875, lsh_876, lsh_877, \
                         lsh_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1161[k] = f_4 * lsg0_628[k]
                    - f_5 * lsg1_628[k]
                    + f_3 * pc_x[k] * lsh_874[k];

        t_1162[k] = f_4 * lsg0_629[k]
                    - f_5 * lsg1_629[k]
                    + f_3 * pc_x[k] * lsh_875[k];

        t_1163[k] = f_3 * pc_x[k] * lsh_876[k];

        t_1164[k] = f_3 * pc_x[k] * lsh_877[k];

        t_1165[k] = f_3 * pc_x[k] * lsh_878[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, t_1169, t_1170, pc_x, pc_y, pc_z, ksh_687, \
                         ksh_708, lsg0_625, lsg1_625, lsh_876, lsh_879, lsh_880, \
                         lsh_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_3 * pc_x[k] * lsh_879[k];

        t_1167[k] = f_3 * pc_x[k] * lsh_880[k];

        t_1168[k] = f_3 * pc_x[k] * lsh_881[k];

        t_1169[k] = f_13 * ksh_708[k]
                    + f_1 * lsg0_625[k]
                    - f_2 * lsg1_625[k]
                    + f_3 * pc_y[k] * lsh_876[k];

        t_1170[k] = f_19 * ksh_687[k]
                    + f_3 * pc_z[k] * lsh_876[k];
    }

#pragma omp simd aligned(t_1171, t_1172, t_1173, pc_y, ksh_710, ksh_711, ksh_712, lsg0_627, \
                         lsg0_628, lsg0_629, lsg1_627, lsg1_628, lsg1_629, lsh_878, lsh_879, \
                         lsh_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1171[k] = f_13 * ksh_710[k]
                    + f_8 * lsg0_627[k]
                    - f_9 * lsg1_627[k]
                    + f_3 * pc_y[k] * lsh_878[k];

        t_1172[k] = f_13 * ksh_711[k]
                    + f_6 * lsg0_628[k]
                    - f_7 * lsg1_628[k]
                    + f_3 * pc_y[k] * lsh_879[k];

        t_1173[k] = f_13 * ksh_712[k]
                    + f_4 * lsg0_629[k]
                    - f_5 * lsg1_629[k]
                    + f_3 * pc_y[k] * lsh_880[k];
    }

#pragma omp simd aligned(t_1174, t_1175, t_1176, pc_x, pc_y, pc_z, ksh_692, ksh_713, lsg0_629, \
                         lsg0_630, lsg1_629, lsg1_630, lsh_881, \
                         lsh_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1174[k] = f_13 * ksh_713[k]
                    + f_3 * pc_y[k] * lsh_881[k];

        t_1175[k] = f_19 * ksh_692[k]
                    + f_1 * lsg0_629[k]
                    - f_2 * lsg1_629[k]
                    + f_3 * pc_z[k] * lsh_881[k];

        t_1176[k] = f_1 * lsg0_630[k]
                    - f_2 * lsg1_630[k]
                    + f_3 * pc_x[k] * lsh_882[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, pc_x, lsg0_631, lsg0_632, lsg0_633, lsg1_631, \
                         lsg1_632, lsg1_633, lsh_883, lsh_884, \
                         lsh_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_16 * lsg0_631[k]
                    - f_17 * lsg1_631[k]
                    + f_3 * pc_x[k] * lsh_883[k];

        t_1178[k] = f_16 * lsg0_632[k]
                    - f_17 * lsg1_632[k]
                    + f_3 * pc_x[k] * lsh_884[k];

        t_1179[k] = f_8 * lsg0_633[k]
                    - f_9 * lsg1_633[k]
                    + f_3 * pc_x[k] * lsh_885[k];
    }

#pragma omp simd aligned(t_1180, t_1181, t_1182, pc_x, lsg0_634, lsg0_635, lsg0_636, lsg1_634, \
                         lsg1_635, lsg1_636, lsh_886, lsh_887, \
                         lsh_888 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1180[k] = f_8 * lsg0_634[k]
                    - f_9 * lsg1_634[k]
                    + f_3 * pc_x[k] * lsh_886[k];

        t_1181[k] = f_8 * lsg0_635[k]
                    - f_9 * lsg1_635[k]
                    + f_3 * pc_x[k] * lsh_887[k];

        t_1182[k] = f_6 * lsg0_636[k]
                    - f_7 * lsg1_636[k]
                    + f_3 * pc_x[k] * lsh_888[k];
    }
}

static auto
compute_prim_lsi_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pa,
                                                           const size_t pc, const size_t ksi0,
                                                           const size_t ksh, const size_t ksi1,
                                                           const size_t lsg0, const size_t lsg1,
                                                           const size_t lsh, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_1 = 2.5 / gamma;
    const auto f_2 = 2.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = gamma / q;
    const auto f_11 = 0.5 / q;
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 3.0 / q;

    auto *t_1183 = buffer.data(target + 1183);
    auto *t_1184 = buffer.data(target + 1184);
    auto *t_1185 = buffer.data(target + 1185);
    auto *t_1186 = buffer.data(target + 1186);
    auto *t_1187 = buffer.data(target + 1187);
    auto *t_1188 = buffer.data(target + 1188);
    auto *t_1189 = buffer.data(target + 1189);
    auto *t_1190 = buffer.data(target + 1190);
    auto *t_1191 = buffer.data(target + 1191);
    auto *t_1192 = buffer.data(target + 1192);
    auto *t_1193 = buffer.data(target + 1193);
    auto *t_1194 = buffer.data(target + 1194);
    auto *t_1195 = buffer.data(target + 1195);
    auto *t_1196 = buffer.data(target + 1196);
    auto *t_1197 = buffer.data(target + 1197);
    auto *t_1198 = buffer.data(target + 1198);
    auto *t_1199 = buffer.data(target + 1199);
    auto *t_1200 = buffer.data(target + 1200);
    auto *t_1201 = buffer.data(target + 1201);
    auto *t_1202 = buffer.data(target + 1202);
    auto *t_1203 = buffer.data(target + 1203);
    auto *t_1204 = buffer.data(target + 1204);
    auto *t_1205 = buffer.data(target + 1205);
    auto *t_1206 = buffer.data(target + 1206);
    auto *t_1207 = buffer.data(target + 1207);
    auto *t_1208 = buffer.data(target + 1208);
    auto *t_1209 = buffer.data(target + 1209);
    auto *t_1210 = buffer.data(target + 1210);
    auto *t_1211 = buffer.data(target + 1211);
    auto *t_1212 = buffer.data(target + 1212);
    auto *t_1213 = buffer.data(target + 1213);
    auto *t_1214 = buffer.data(target + 1214);
    auto *t_1215 = buffer.data(target + 1215);
    auto *t_1216 = buffer.data(target + 1216);
    auto *t_1217 = buffer.data(target + 1217);
    auto *t_1218 = buffer.data(target + 1218);
    auto *t_1219 = buffer.data(target + 1219);
    auto *t_1220 = buffer.data(target + 1220);
    auto *t_1221 = buffer.data(target + 1221);
    auto *t_1222 = buffer.data(target + 1222);
    auto *t_1223 = buffer.data(target + 1223);
    auto *t_1224 = buffer.data(target + 1224);
    auto *t_1225 = buffer.data(target + 1225);
    auto *t_1226 = buffer.data(target + 1226);
    auto *t_1227 = buffer.data(target + 1227);
    auto *t_1228 = buffer.data(target + 1228);
    auto *t_1229 = buffer.data(target + 1229);
    auto *t_1230 = buffer.data(target + 1230);
    auto *t_1231 = buffer.data(target + 1231);
    auto *t_1232 = buffer.data(target + 1232);
    auto *t_1233 = buffer.data(target + 1233);
    auto *t_1234 = buffer.data(target + 1234);
    auto *t_1235 = buffer.data(target + 1235);
    auto *t_1236 = buffer.data(target + 1236);
    auto *t_1237 = buffer.data(target + 1237);
    auto *t_1238 = buffer.data(target + 1238);
    auto *t_1239 = buffer.data(target + 1239);
    auto *t_1240 = buffer.data(target + 1240);
    auto *t_1241 = buffer.data(target + 1241);
    auto *t_1242 = buffer.data(target + 1242);
    auto *t_1243 = buffer.data(target + 1243);
    auto *t_1244 = buffer.data(target + 1244);
    auto *t_1245 = buffer.data(target + 1245);
    auto *t_1246 = buffer.data(target + 1246);
    auto *t_1247 = buffer.data(target + 1247);
    auto *t_1248 = buffer.data(target + 1248);
    auto *t_1249 = buffer.data(target + 1249);
    auto *t_1250 = buffer.data(target + 1250);
    auto *t_1251 = buffer.data(target + 1251);
    auto *t_1252 = buffer.data(target + 1252);
    auto *t_1253 = buffer.data(target + 1253);
    auto *t_1254 = buffer.data(target + 1254);
    auto *t_1255 = buffer.data(target + 1255);
    auto *t_1256 = buffer.data(target + 1256);
    auto *t_1257 = buffer.data(target + 1257);
    auto *t_1258 = buffer.data(target + 1258);
    auto *t_1259 = buffer.data(target + 1259);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ksi0_980 = buffer.data(ksi0 + 980);
    const auto *ksi0_982 = buffer.data(ksi0 + 982);
    const auto *ksi0_985 = buffer.data(ksi0 + 985);
    const auto *ksi0_989 = buffer.data(ksi0 + 989);
    const auto *ksi0_994 = buffer.data(ksi0 + 994);
    const auto *ksi0_1001 = buffer.data(ksi0 + 1001);
    const auto *ksi0_1003 = buffer.data(ksi0 + 1003);
    const auto *ksi0_1004 = buffer.data(ksi0 + 1004);
    const auto *ksi0_1005 = buffer.data(ksi0 + 1005);
    const auto *ksi0_1007 = buffer.data(ksi0 + 1007);

    const auto *ksh_708 = buffer.data(ksh + 708);
    const auto *ksh_713 = buffer.data(ksh + 713);
    const auto *ksh_729 = buffer.data(ksh + 729);
    const auto *ksh_731 = buffer.data(ksh + 731);
    const auto *ksh_732 = buffer.data(ksh + 732);
    const auto *ksh_733 = buffer.data(ksh + 733);
    const auto *ksh_734 = buffer.data(ksh + 734);
    const auto *ksh_750 = buffer.data(ksh + 750);
    const auto *ksh_752 = buffer.data(ksh + 752);
    const auto *ksh_753 = buffer.data(ksh + 753);
    const auto *ksh_754 = buffer.data(ksh + 754);
    const auto *ksh_755 = buffer.data(ksh + 755);

    const auto *ksi1_980 = buffer.data(ksi1 + 980);
    const auto *ksi1_982 = buffer.data(ksi1 + 982);
    const auto *ksi1_985 = buffer.data(ksi1 + 985);
    const auto *ksi1_989 = buffer.data(ksi1 + 989);
    const auto *ksi1_994 = buffer.data(ksi1 + 994);
    const auto *ksi1_1001 = buffer.data(ksi1 + 1001);
    const auto *ksi1_1003 = buffer.data(ksi1 + 1003);
    const auto *ksi1_1004 = buffer.data(ksi1 + 1004);
    const auto *ksi1_1005 = buffer.data(ksi1 + 1005);
    const auto *ksi1_1007 = buffer.data(ksi1 + 1007);

    const auto *lsg0_637 = buffer.data(lsg0 + 637);
    const auto *lsg0_638 = buffer.data(lsg0 + 638);
    const auto *lsg0_639 = buffer.data(lsg0 + 639);
    const auto *lsg0_640 = buffer.data(lsg0 + 640);
    const auto *lsg0_641 = buffer.data(lsg0 + 641);
    const auto *lsg0_642 = buffer.data(lsg0 + 642);
    const auto *lsg0_643 = buffer.data(lsg0 + 643);
    const auto *lsg0_644 = buffer.data(lsg0 + 644);
    const auto *lsg0_646 = buffer.data(lsg0 + 646);
    const auto *lsg0_648 = buffer.data(lsg0 + 648);
    const auto *lsg0_649 = buffer.data(lsg0 + 649);
    const auto *lsg0_651 = buffer.data(lsg0 + 651);
    const auto *lsg0_652 = buffer.data(lsg0 + 652);
    const auto *lsg0_653 = buffer.data(lsg0 + 653);
    const auto *lsg0_655 = buffer.data(lsg0 + 655);
    const auto *lsg0_656 = buffer.data(lsg0 + 656);
    const auto *lsg0_657 = buffer.data(lsg0 + 657);
    const auto *lsg0_658 = buffer.data(lsg0 + 658);
    const auto *lsg0_660 = buffer.data(lsg0 + 660);
    const auto *lsg0_662 = buffer.data(lsg0 + 662);
    const auto *lsg0_663 = buffer.data(lsg0 + 663);
    const auto *lsg0_665 = buffer.data(lsg0 + 665);
    const auto *lsg0_666 = buffer.data(lsg0 + 666);
    const auto *lsg0_667 = buffer.data(lsg0 + 667);
    const auto *lsg0_669 = buffer.data(lsg0 + 669);
    const auto *lsg0_670 = buffer.data(lsg0 + 670);
    const auto *lsg0_671 = buffer.data(lsg0 + 671);
    const auto *lsg0_672 = buffer.data(lsg0 + 672);
    const auto *lsg0_673 = buffer.data(lsg0 + 673);
    const auto *lsg0_674 = buffer.data(lsg0 + 674);

    const auto *lsg1_637 = buffer.data(lsg1 + 637);
    const auto *lsg1_638 = buffer.data(lsg1 + 638);
    const auto *lsg1_639 = buffer.data(lsg1 + 639);
    const auto *lsg1_640 = buffer.data(lsg1 + 640);
    const auto *lsg1_641 = buffer.data(lsg1 + 641);
    const auto *lsg1_642 = buffer.data(lsg1 + 642);
    const auto *lsg1_643 = buffer.data(lsg1 + 643);
    const auto *lsg1_644 = buffer.data(lsg1 + 644);
    const auto *lsg1_646 = buffer.data(lsg1 + 646);
    const auto *lsg1_648 = buffer.data(lsg1 + 648);
    const auto *lsg1_649 = buffer.data(lsg1 + 649);
    const auto *lsg1_651 = buffer.data(lsg1 + 651);
    const auto *lsg1_652 = buffer.data(lsg1 + 652);
    const auto *lsg1_653 = buffer.data(lsg1 + 653);
    const auto *lsg1_655 = buffer.data(lsg1 + 655);
    const auto *lsg1_656 = buffer.data(lsg1 + 656);
    const auto *lsg1_657 = buffer.data(lsg1 + 657);
    const auto *lsg1_658 = buffer.data(lsg1 + 658);
    const auto *lsg1_660 = buffer.data(lsg1 + 660);
    const auto *lsg1_662 = buffer.data(lsg1 + 662);
    const auto *lsg1_663 = buffer.data(lsg1 + 663);
    const auto *lsg1_665 = buffer.data(lsg1 + 665);
    const auto *lsg1_666 = buffer.data(lsg1 + 666);
    const auto *lsg1_667 = buffer.data(lsg1 + 667);
    const auto *lsg1_669 = buffer.data(lsg1 + 669);
    const auto *lsg1_670 = buffer.data(lsg1 + 670);
    const auto *lsg1_671 = buffer.data(lsg1 + 671);
    const auto *lsg1_672 = buffer.data(lsg1 + 672);
    const auto *lsg1_673 = buffer.data(lsg1 + 673);
    const auto *lsg1_674 = buffer.data(lsg1 + 674);

    const auto *lsh_889 = buffer.data(lsh + 889);
    const auto *lsh_890 = buffer.data(lsh + 890);
    const auto *lsh_891 = buffer.data(lsh + 891);
    const auto *lsh_892 = buffer.data(lsh + 892);
    const auto *lsh_893 = buffer.data(lsh + 893);
    const auto *lsh_894 = buffer.data(lsh + 894);
    const auto *lsh_895 = buffer.data(lsh + 895);
    const auto *lsh_896 = buffer.data(lsh + 896);
    const auto *lsh_897 = buffer.data(lsh + 897);
    const auto *lsh_898 = buffer.data(lsh + 898);
    const auto *lsh_899 = buffer.data(lsh + 899);
    const auto *lsh_900 = buffer.data(lsh + 900);
    const auto *lsh_901 = buffer.data(lsh + 901);
    const auto *lsh_902 = buffer.data(lsh + 902);
    const auto *lsh_904 = buffer.data(lsh + 904);
    const auto *lsh_906 = buffer.data(lsh + 906);
    const auto *lsh_907 = buffer.data(lsh + 907);
    const auto *lsh_909 = buffer.data(lsh + 909);
    const auto *lsh_910 = buffer.data(lsh + 910);
    const auto *lsh_911 = buffer.data(lsh + 911);
    const auto *lsh_913 = buffer.data(lsh + 913);
    const auto *lsh_914 = buffer.data(lsh + 914);
    const auto *lsh_915 = buffer.data(lsh + 915);
    const auto *lsh_916 = buffer.data(lsh + 916);
    const auto *lsh_918 = buffer.data(lsh + 918);
    const auto *lsh_919 = buffer.data(lsh + 919);
    const auto *lsh_920 = buffer.data(lsh + 920);
    const auto *lsh_921 = buffer.data(lsh + 921);
    const auto *lsh_922 = buffer.data(lsh + 922);
    const auto *lsh_923 = buffer.data(lsh + 923);
    const auto *lsh_924 = buffer.data(lsh + 924);
    const auto *lsh_926 = buffer.data(lsh + 926);
    const auto *lsh_927 = buffer.data(lsh + 927);
    const auto *lsh_929 = buffer.data(lsh + 929);
    const auto *lsh_930 = buffer.data(lsh + 930);
    const auto *lsh_931 = buffer.data(lsh + 931);
    const auto *lsh_933 = buffer.data(lsh + 933);
    const auto *lsh_934 = buffer.data(lsh + 934);
    const auto *lsh_935 = buffer.data(lsh + 935);
    const auto *lsh_936 = buffer.data(lsh + 936);
    const auto *lsh_938 = buffer.data(lsh + 938);
    const auto *lsh_939 = buffer.data(lsh + 939);
    const auto *lsh_940 = buffer.data(lsh + 940);
    const auto *lsh_941 = buffer.data(lsh + 941);
    const auto *lsh_942 = buffer.data(lsh + 942);
    const auto *lsh_943 = buffer.data(lsh + 943);
    const auto *lsh_944 = buffer.data(lsh + 944);

#pragma omp simd aligned(t_1183, t_1184, t_1185, pc_x, lsg0_637, lsg0_638, lsg0_639, lsg1_637, \
                         lsg1_638, lsg1_639, lsh_889, lsh_890, \
                         lsh_891 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1183[k] = f_6 * lsg0_637[k]
                    - f_7 * lsg1_637[k]
                    + f_3 * pc_x[k] * lsh_889[k];

        t_1184[k] = f_6 * lsg0_638[k]
                    - f_7 * lsg1_638[k]
                    + f_3 * pc_x[k] * lsh_890[k];

        t_1185[k] = f_6 * lsg0_639[k]
                    - f_7 * lsg1_639[k]
                    + f_3 * pc_x[k] * lsh_891[k];
    }

#pragma omp simd aligned(t_1186, t_1187, t_1188, pc_x, lsg0_640, lsg0_641, lsg0_642, lsg1_640, \
                         lsg1_641, lsg1_642, lsh_892, lsh_893, \
                         lsh_894 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1186[k] = f_4 * lsg0_640[k]
                    - f_5 * lsg1_640[k]
                    + f_3 * pc_x[k] * lsh_892[k];

        t_1187[k] = f_4 * lsg0_641[k]
                    - f_5 * lsg1_641[k]
                    + f_3 * pc_x[k] * lsh_893[k];

        t_1188[k] = f_4 * lsg0_642[k]
                    - f_5 * lsg1_642[k]
                    + f_3 * pc_x[k] * lsh_894[k];
    }

#pragma omp simd aligned(t_1189, t_1190, t_1191, t_1192, t_1193, pc_x, lsg0_643, lsg0_644, \
                         lsg1_643, lsg1_644, lsh_895, lsh_896, lsh_897, lsh_898, \
                         lsh_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1189[k] = f_4 * lsg0_643[k]
                    - f_5 * lsg1_643[k]
                    + f_3 * pc_x[k] * lsh_895[k];

        t_1190[k] = f_4 * lsg0_644[k]
                    - f_5 * lsg1_644[k]
                    + f_3 * pc_x[k] * lsh_896[k];

        t_1191[k] = f_3 * pc_x[k] * lsh_897[k];

        t_1192[k] = f_3 * pc_x[k] * lsh_898[k];

        t_1193[k] = f_3 * pc_x[k] * lsh_899[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, t_1197, t_1198, pc_x, pc_y, pc_z, ksh_708, \
                         ksh_729, lsg0_640, lsg1_640, lsh_897, lsh_900, lsh_901, \
                         lsh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = f_3 * pc_x[k] * lsh_900[k];

        t_1195[k] = f_3 * pc_x[k] * lsh_901[k];

        t_1196[k] = f_3 * pc_x[k] * lsh_902[k];

        t_1197[k] = f_12 * ksh_729[k]
                    + f_1 * lsg0_640[k]
                    - f_2 * lsg1_640[k]
                    + f_3 * pc_y[k] * lsh_897[k];

        t_1198[k] = f_18 * ksh_708[k]
                    + f_3 * pc_z[k] * lsh_897[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, pc_y, ksh_731, ksh_732, ksh_733, lsg0_642, \
                         lsg0_643, lsg0_644, lsg1_642, lsg1_643, lsg1_644, lsh_899, lsh_900, \
                         lsh_901 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_12 * ksh_731[k]
                    + f_8 * lsg0_642[k]
                    - f_9 * lsg1_642[k]
                    + f_3 * pc_y[k] * lsh_899[k];

        t_1200[k] = f_12 * ksh_732[k]
                    + f_6 * lsg0_643[k]
                    - f_7 * lsg1_643[k]
                    + f_3 * pc_y[k] * lsh_900[k];

        t_1201[k] = f_12 * ksh_733[k]
                    + f_4 * lsg0_644[k]
                    - f_5 * lsg1_644[k]
                    + f_3 * pc_y[k] * lsh_901[k];
    }

#pragma omp simd aligned(t_1202, t_1203, t_1204, pa_y, pc_y, pc_z, ksi0_980, ksh_713, ksh_734, \
                         ksi1_980, lsg0_644, lsg1_644, lsh_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1202[k] = f_12 * ksh_734[k]
                    + f_3 * pc_y[k] * lsh_902[k];

        t_1203[k] = f_18 * ksh_713[k]
                    + f_1 * lsg0_644[k]
                    - f_2 * lsg1_644[k]
                    + f_3 * pc_z[k] * lsh_902[k];

        t_1204[k] = pa_y[k] * ksi0_980[k]
                    - f_10 * pc_y[k] * ksi1_980[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, pa_y, pc_x, pc_y, ksi0_982, ksi1_982, \
                         lsg0_646, lsg0_648, lsg1_646, lsg1_648, lsh_904, \
                         lsh_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_16 * lsg0_646[k]
                    - f_17 * lsg1_646[k]
                    + f_3 * pc_x[k] * lsh_904[k];

        t_1206[k] = pa_y[k] * ksi0_982[k]
                    - f_10 * pc_y[k] * ksi1_982[k];

        t_1207[k] = f_8 * lsg0_648[k]
                    - f_9 * lsg1_648[k]
                    + f_3 * pc_x[k] * lsh_906[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, pa_y, pc_x, pc_y, ksi0_985, ksi1_985, \
                         lsg0_649, lsg0_651, lsg1_649, lsg1_651, lsh_907, \
                         lsh_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_8 * lsg0_649[k]
                    - f_9 * lsg1_649[k]
                    + f_3 * pc_x[k] * lsh_907[k];

        t_1209[k] = pa_y[k] * ksi0_985[k]
                    - f_10 * pc_y[k] * ksi1_985[k];

        t_1210[k] = f_6 * lsg0_651[k]
                    - f_7 * lsg1_651[k]
                    + f_3 * pc_x[k] * lsh_909[k];
    }

#pragma omp simd aligned(t_1211, t_1212, t_1213, pa_y, pc_x, pc_y, ksi0_989, ksi1_989, \
                         lsg0_652, lsg0_653, lsg1_652, lsg1_653, lsh_910, \
                         lsh_911 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_6 * lsg0_652[k]
                    - f_7 * lsg1_652[k]
                    + f_3 * pc_x[k] * lsh_910[k];

        t_1212[k] = f_6 * lsg0_653[k]
                    - f_7 * lsg1_653[k]
                    + f_3 * pc_x[k] * lsh_911[k];

        t_1213[k] = pa_y[k] * ksi0_989[k]
                    - f_10 * pc_y[k] * ksi1_989[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, pc_x, lsg0_655, lsg0_656, lsg0_657, lsg1_655, \
                         lsg1_656, lsg1_657, lsh_913, lsh_914, \
                         lsh_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = f_4 * lsg0_655[k]
                    - f_5 * lsg1_655[k]
                    + f_3 * pc_x[k] * lsh_913[k];

        t_1215[k] = f_4 * lsg0_656[k]
                    - f_5 * lsg1_656[k]
                    + f_3 * pc_x[k] * lsh_914[k];

        t_1216[k] = f_4 * lsg0_657[k]
                    - f_5 * lsg1_657[k]
                    + f_3 * pc_x[k] * lsh_915[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, t_1221, pa_y, pc_x, pc_y, ksi0_994, \
                         ksi1_994, lsg0_658, lsg1_658, lsh_916, lsh_918, lsh_919, \
                         lsh_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_4 * lsg0_658[k]
                    - f_5 * lsg1_658[k]
                    + f_3 * pc_x[k] * lsh_916[k];

        t_1218[k] = pa_y[k] * ksi0_994[k]
                    - f_10 * pc_y[k] * ksi1_994[k];

        t_1219[k] = f_3 * pc_x[k] * lsh_918[k];

        t_1220[k] = f_3 * pc_x[k] * lsh_919[k];

        t_1221[k] = f_3 * pc_x[k] * lsh_920[k];
    }

#pragma omp simd aligned(t_1222, t_1223, t_1224, t_1225, pa_y, pc_x, pc_y, ksi0_1001, ksh_750, \
                         ksi1_1001, lsh_921, lsh_922, lsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1222[k] = f_3 * pc_x[k] * lsh_921[k];

        t_1223[k] = f_3 * pc_x[k] * lsh_922[k];

        t_1224[k] = f_3 * pc_x[k] * lsh_923[k];

        t_1225[k] = pa_y[k] * ksi0_1001[k]
                    + f_18 * ksh_750[k]
                    - f_10 * pc_y[k] * ksi1_1001[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, pa_y, pc_y, pc_z, ksi0_1003, ksi0_1004, \
                         ksh_729, ksh_752, ksh_753, ksi1_1003, ksi1_1004, \
                         lsh_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_15 * ksh_729[k]
                    + f_3 * pc_z[k] * lsh_918[k];

        t_1227[k] = pa_y[k] * ksi0_1003[k]
                    + f_14 * ksh_752[k]
                    - f_10 * pc_y[k] * ksi1_1003[k];

        t_1228[k] = pa_y[k] * ksi0_1004[k]
                    + f_13 * ksh_753[k]
                    - f_10 * pc_y[k] * ksi1_1004[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, pa_y, pc_y, ksi0_1005, ksi0_1007, ksh_754, \
                         ksh_755, ksi1_1005, ksi1_1007, lsh_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = pa_y[k] * ksi0_1005[k]
                    + f_12 * ksh_754[k]
                    - f_10 * pc_y[k] * ksi1_1005[k];

        t_1230[k] = f_11 * ksh_755[k]
                    + f_3 * pc_y[k] * lsh_923[k];

        t_1231[k] = pa_y[k] * ksi0_1007[k]
                    - f_10 * pc_y[k] * ksi1_1007[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, t_1236, pc_x, pc_y, lsg0_660, \
                         lsg0_662, lsg0_663, lsg1_660, lsg1_662, lsg1_663, lsh_924, lsh_926, \
                         lsh_927 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_1 * lsg0_660[k]
                    - f_2 * lsg1_660[k]
                    + f_3 * pc_x[k] * lsh_924[k];

        t_1233[k] = f_3 * pc_y[k] * lsh_924[k];

        t_1234[k] = f_16 * lsg0_662[k]
                    - f_17 * lsg1_662[k]
                    + f_3 * pc_x[k] * lsh_926[k];

        t_1235[k] = f_8 * lsg0_663[k]
                    - f_9 * lsg1_663[k]
                    + f_3 * pc_x[k] * lsh_927[k];

        t_1236[k] = f_3 * pc_y[k] * lsh_926[k];
    }

#pragma omp simd aligned(t_1237, t_1238, t_1239, t_1240, pc_x, pc_y, lsg0_665, lsg0_666, \
                         lsg0_667, lsg1_665, lsg1_666, lsg1_667, lsh_929, lsh_930, \
                         lsh_931 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1237[k] = f_8 * lsg0_665[k]
                    - f_9 * lsg1_665[k]
                    + f_3 * pc_x[k] * lsh_929[k];

        t_1238[k] = f_6 * lsg0_666[k]
                    - f_7 * lsg1_666[k]
                    + f_3 * pc_x[k] * lsh_930[k];

        t_1239[k] = f_6 * lsg0_667[k]
                    - f_7 * lsg1_667[k]
                    + f_3 * pc_x[k] * lsh_931[k];

        t_1240[k] = f_3 * pc_y[k] * lsh_929[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, pc_x, lsg0_669, lsg0_670, lsg0_671, lsg1_669, \
                         lsg1_670, lsg1_671, lsh_933, lsh_934, \
                         lsh_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_6 * lsg0_669[k]
                    - f_7 * lsg1_669[k]
                    + f_3 * pc_x[k] * lsh_933[k];

        t_1242[k] = f_4 * lsg0_670[k]
                    - f_5 * lsg1_670[k]
                    + f_3 * pc_x[k] * lsh_934[k];

        t_1243[k] = f_4 * lsg0_671[k]
                    - f_5 * lsg1_671[k]
                    + f_3 * pc_x[k] * lsh_935[k];
    }

#pragma omp simd aligned(t_1244, t_1245, t_1246, t_1247, t_1248, pc_x, pc_y, lsg0_672, \
                         lsg0_674, lsg1_672, lsg1_674, lsh_933, lsh_936, lsh_938, lsh_939, \
                         lsh_940 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1244[k] = f_4 * lsg0_672[k]
                    - f_5 * lsg1_672[k]
                    + f_3 * pc_x[k] * lsh_936[k];

        t_1245[k] = f_3 * pc_y[k] * lsh_933[k];

        t_1246[k] = f_4 * lsg0_674[k]
                    - f_5 * lsg1_674[k]
                    + f_3 * pc_x[k] * lsh_938[k];

        t_1247[k] = f_3 * pc_x[k] * lsh_939[k];

        t_1248[k] = f_3 * pc_x[k] * lsh_940[k];
    }

#pragma omp simd aligned(t_1249, t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, lsg0_670, \
                         lsg1_670, lsh_939, lsh_941, lsh_942, lsh_943, \
                         lsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1249[k] = f_3 * pc_x[k] * lsh_941[k];

        t_1250[k] = f_3 * pc_x[k] * lsh_942[k];

        t_1251[k] = f_3 * pc_x[k] * lsh_943[k];

        t_1252[k] = f_3 * pc_x[k] * lsh_944[k];

        t_1253[k] = f_1 * lsg0_670[k]
                    - f_2 * lsg1_670[k]
                    + f_3 * pc_y[k] * lsh_939[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, lsg0_671, lsg0_672, lsg0_673, lsg1_671, \
                         lsg1_672, lsg1_673, lsh_940, lsh_941, \
                         lsh_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_16 * lsg0_671[k]
                    - f_17 * lsg1_671[k]
                    + f_3 * pc_y[k] * lsh_940[k];

        t_1255[k] = f_8 * lsg0_672[k]
                    - f_9 * lsg1_672[k]
                    + f_3 * pc_y[k] * lsh_941[k];

        t_1256[k] = f_6 * lsg0_673[k]
                    - f_7 * lsg1_673[k]
                    + f_3 * pc_y[k] * lsh_942[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, ksh_755, lsg0_674, lsg1_674, \
                         lsh_943, lsh_944 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_4 * lsg0_674[k]
                    - f_5 * lsg1_674[k]
                    + f_3 * pc_y[k] * lsh_943[k];

        t_1258[k] = f_3 * pc_y[k] * lsh_944[k];

        t_1259[k] = f_0 * ksh_755[k]
                    + f_1 * lsg0_674[k]
                    - f_2 * lsg1_674[k]
                    + f_3 * pc_z[k] * lsh_944[k];
    }
}

auto
compute_prim_lsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t ksi0, const size_t ksh,
                                                   const size_t ksi1, const size_t lsg0,
                                                   const size_t lsg1, const size_t lsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_lsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsh, ncols, gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, ksi0, ksh,
                                                              ksi1, lsg0, lsg1, lsh, ncols,
                                                              gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece9(buffer, target, pc, ksh, lsg0,
                                                              lsg1, lsh, ncols, gamma, p, q);

    compute_prim_lsi_three_center_electron_repulsion_0_piece10(buffer, target, pa, pc, ksi0,
                                                               ksh, ksi1, lsg0, lsg1, lsh,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
