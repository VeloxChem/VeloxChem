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


#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *hsi0_0 = buffer.data(hsi0 + 0);
    const auto *hsi0_3 = buffer.data(hsi0 + 3);
    const auto *hsi0_5 = buffer.data(hsi0 + 5);
    const auto *hsi0_6 = buffer.data(hsi0 + 6);
    const auto *hsi0_9 = buffer.data(hsi0 + 9);
    const auto *hsi0_10 = buffer.data(hsi0 + 10);
    const auto *hsi0_14 = buffer.data(hsi0 + 14);
    const auto *hsi0_21 = buffer.data(hsi0 + 21);
    const auto *hsi0_27 = buffer.data(hsi0 + 27);
    const auto *hsi0_31 = buffer.data(hsi0 + 31);
    const auto *hsi0_34 = buffer.data(hsi0 + 34);
    const auto *hsi0_38 = buffer.data(hsi0 + 38);
    const auto *hsi0_56 = buffer.data(hsi0 + 56);
    const auto *hsi0_61 = buffer.data(hsi0 + 61);
    const auto *hsi0_65 = buffer.data(hsi0 + 65);
    const auto *hsi0_68 = buffer.data(hsi0 + 68);
    const auto *hsi0_70 = buffer.data(hsi0 + 70);

    const auto *hsh_0 = buffer.data(hsh + 0);
    const auto *hsh_1 = buffer.data(hsh + 1);
    const auto *hsh_2 = buffer.data(hsh + 2);
    const auto *hsh_3 = buffer.data(hsh + 3);
    const auto *hsh_5 = buffer.data(hsh + 5);
    const auto *hsh_6 = buffer.data(hsh + 6);
    const auto *hsh_9 = buffer.data(hsh + 9);
    const auto *hsh_15 = buffer.data(hsh + 15);
    const auto *hsh_17 = buffer.data(hsh + 17);
    const auto *hsh_18 = buffer.data(hsh + 18);
    const auto *hsh_20 = buffer.data(hsh + 20);
    const auto *hsh_21 = buffer.data(hsh + 21);
    const auto *hsh_24 = buffer.data(hsh + 24);
    const auto *hsh_26 = buffer.data(hsh + 26);
    const auto *hsh_27 = buffer.data(hsh + 27);
    const auto *hsh_30 = buffer.data(hsh + 30);
    const auto *hsh_36 = buffer.data(hsh + 36);
    const auto *hsh_38 = buffer.data(hsh + 38);
    const auto *hsh_39 = buffer.data(hsh + 39);
    const auto *hsh_40 = buffer.data(hsh + 40);
    const auto *hsh_41 = buffer.data(hsh + 41);
    const auto *hsh_42 = buffer.data(hsh + 42);
    const auto *hsh_44 = buffer.data(hsh + 44);
    const auto *hsh_47 = buffer.data(hsh + 47);
    const auto *hsh_50 = buffer.data(hsh + 50);
    const auto *hsh_51 = buffer.data(hsh + 51);
    const auto *hsh_57 = buffer.data(hsh + 57);
    const auto *hsh_58 = buffer.data(hsh + 58);
    const auto *hsh_59 = buffer.data(hsh + 59);
    const auto *hsh_60 = buffer.data(hsh + 60);
    const auto *hsh_62 = buffer.data(hsh + 62);
    const auto *hsh_63 = buffer.data(hsh + 63);
    const auto *hsh_66 = buffer.data(hsh + 66);
    const auto *hsh_69 = buffer.data(hsh + 69);
    const auto *hsh_73 = buffer.data(hsh + 73);
    const auto *hsh_78 = buffer.data(hsh + 78);
    const auto *hsh_80 = buffer.data(hsh + 80);
    const auto *hsh_81 = buffer.data(hsh + 81);
    const auto *hsh_82 = buffer.data(hsh + 82);
    const auto *hsh_83 = buffer.data(hsh + 83);
    const auto *hsh_99 = buffer.data(hsh + 99);

    const auto *hsi1_0 = buffer.data(hsi1 + 0);
    const auto *hsi1_3 = buffer.data(hsi1 + 3);
    const auto *hsi1_5 = buffer.data(hsi1 + 5);
    const auto *hsi1_6 = buffer.data(hsi1 + 6);
    const auto *hsi1_9 = buffer.data(hsi1 + 9);
    const auto *hsi1_10 = buffer.data(hsi1 + 10);
    const auto *hsi1_14 = buffer.data(hsi1 + 14);
    const auto *hsi1_21 = buffer.data(hsi1 + 21);
    const auto *hsi1_27 = buffer.data(hsi1 + 27);
    const auto *hsi1_31 = buffer.data(hsi1 + 31);
    const auto *hsi1_34 = buffer.data(hsi1 + 34);
    const auto *hsi1_38 = buffer.data(hsi1 + 38);
    const auto *hsi1_56 = buffer.data(hsi1 + 56);
    const auto *hsi1_61 = buffer.data(hsi1 + 61);
    const auto *hsi1_65 = buffer.data(hsi1 + 65);
    const auto *hsi1_68 = buffer.data(hsi1 + 68);
    const auto *hsi1_70 = buffer.data(hsi1 + 70);

    const auto *isg0_0 = buffer.data(isg0 + 0);
    const auto *isg0_1 = buffer.data(isg0 + 1);
    const auto *isg0_2 = buffer.data(isg0 + 2);
    const auto *isg0_3 = buffer.data(isg0 + 3);
    const auto *isg0_5 = buffer.data(isg0 + 5);
    const auto *isg0_10 = buffer.data(isg0 + 10);
    const auto *isg0_12 = buffer.data(isg0 + 12);
    const auto *isg0_13 = buffer.data(isg0 + 13);
    const auto *isg0_14 = buffer.data(isg0 + 14);
    const auto *isg0_18 = buffer.data(isg0 + 18);
    const auto *isg0_25 = buffer.data(isg0 + 25);
    const auto *isg0_26 = buffer.data(isg0 + 26);
    const auto *isg0_27 = buffer.data(isg0 + 27);
    const auto *isg0_32 = buffer.data(isg0 + 32);
    const auto *isg0_34 = buffer.data(isg0 + 34);
    const auto *isg0_35 = buffer.data(isg0 + 35);
    const auto *isg0_41 = buffer.data(isg0 + 41);
    const auto *isg0_42 = buffer.data(isg0 + 42);
    const auto *isg0_43 = buffer.data(isg0 + 43);
    const auto *isg0_44 = buffer.data(isg0 + 44);
    const auto *isg0_45 = buffer.data(isg0 + 45);
    const auto *isg0_47 = buffer.data(isg0 + 47);
    const auto *isg0_48 = buffer.data(isg0 + 48);
    const auto *isg0_50 = buffer.data(isg0 + 50);
    const auto *isg0_51 = buffer.data(isg0 + 51);
    const auto *isg0_55 = buffer.data(isg0 + 55);
    const auto *isg0_56 = buffer.data(isg0 + 56);
    const auto *isg0_57 = buffer.data(isg0 + 57);
    const auto *isg0_59 = buffer.data(isg0 + 59);

    const auto *isg1_0 = buffer.data(isg1 + 0);
    const auto *isg1_1 = buffer.data(isg1 + 1);
    const auto *isg1_2 = buffer.data(isg1 + 2);
    const auto *isg1_3 = buffer.data(isg1 + 3);
    const auto *isg1_5 = buffer.data(isg1 + 5);
    const auto *isg1_10 = buffer.data(isg1 + 10);
    const auto *isg1_12 = buffer.data(isg1 + 12);
    const auto *isg1_13 = buffer.data(isg1 + 13);
    const auto *isg1_14 = buffer.data(isg1 + 14);
    const auto *isg1_18 = buffer.data(isg1 + 18);
    const auto *isg1_25 = buffer.data(isg1 + 25);
    const auto *isg1_26 = buffer.data(isg1 + 26);
    const auto *isg1_27 = buffer.data(isg1 + 27);
    const auto *isg1_32 = buffer.data(isg1 + 32);
    const auto *isg1_34 = buffer.data(isg1 + 34);
    const auto *isg1_35 = buffer.data(isg1 + 35);
    const auto *isg1_41 = buffer.data(isg1 + 41);
    const auto *isg1_42 = buffer.data(isg1 + 42);
    const auto *isg1_43 = buffer.data(isg1 + 43);
    const auto *isg1_44 = buffer.data(isg1 + 44);
    const auto *isg1_45 = buffer.data(isg1 + 45);
    const auto *isg1_47 = buffer.data(isg1 + 47);
    const auto *isg1_48 = buffer.data(isg1 + 48);
    const auto *isg1_50 = buffer.data(isg1 + 50);
    const auto *isg1_51 = buffer.data(isg1 + 51);
    const auto *isg1_55 = buffer.data(isg1 + 55);
    const auto *isg1_56 = buffer.data(isg1 + 56);
    const auto *isg1_57 = buffer.data(isg1 + 57);
    const auto *isg1_59 = buffer.data(isg1 + 59);

    const auto *ish_0 = buffer.data(ish + 0);
    const auto *ish_1 = buffer.data(ish + 1);
    const auto *ish_2 = buffer.data(ish + 2);
    const auto *ish_3 = buffer.data(ish + 3);
    const auto *ish_5 = buffer.data(ish + 5);
    const auto *ish_6 = buffer.data(ish + 6);
    const auto *ish_8 = buffer.data(ish + 8);
    const auto *ish_9 = buffer.data(ish + 9);
    const auto *ish_10 = buffer.data(ish + 10);
    const auto *ish_14 = buffer.data(ish + 14);
    const auto *ish_15 = buffer.data(ish + 15);
    const auto *ish_17 = buffer.data(ish + 17);
    const auto *ish_18 = buffer.data(ish + 18);
    const auto *ish_19 = buffer.data(ish + 19);
    const auto *ish_20 = buffer.data(ish + 20);
    const auto *ish_21 = buffer.data(ish + 21);
    const auto *ish_22 = buffer.data(ish + 22);
    const auto *ish_24 = buffer.data(ish + 24);
    const auto *ish_26 = buffer.data(ish + 26);
    const auto *ish_27 = buffer.data(ish + 27);
    const auto *ish_28 = buffer.data(ish + 28);
    const auto *ish_30 = buffer.data(ish + 30);
    const auto *ish_31 = buffer.data(ish + 31);
    const auto *ish_36 = buffer.data(ish + 36);
    const auto *ish_37 = buffer.data(ish + 37);
    const auto *ish_38 = buffer.data(ish + 38);
    const auto *ish_39 = buffer.data(ish + 39);
    const auto *ish_40 = buffer.data(ish + 40);
    const auto *ish_41 = buffer.data(ish + 41);
    const auto *ish_42 = buffer.data(ish + 42);
    const auto *ish_44 = buffer.data(ish + 44);
    const auto *ish_46 = buffer.data(ish + 46);
    const auto *ish_47 = buffer.data(ish + 47);
    const auto *ish_49 = buffer.data(ish + 49);
    const auto *ish_50 = buffer.data(ish + 50);
    const auto *ish_51 = buffer.data(ish + 51);
    const auto *ish_56 = buffer.data(ish + 56);
    const auto *ish_57 = buffer.data(ish + 57);
    const auto *ish_58 = buffer.data(ish + 58);
    const auto *ish_59 = buffer.data(ish + 59);
    const auto *ish_60 = buffer.data(ish + 60);
    const auto *ish_61 = buffer.data(ish + 61);
    const auto *ish_62 = buffer.data(ish + 62);
    const auto *ish_63 = buffer.data(ish + 63);
    const auto *ish_64 = buffer.data(ish + 64);
    const auto *ish_65 = buffer.data(ish + 65);
    const auto *ish_66 = buffer.data(ish + 66);
    const auto *ish_68 = buffer.data(ish + 68);
    const auto *ish_69 = buffer.data(ish + 69);
    const auto *ish_70 = buffer.data(ish + 70);
    const auto *ish_72 = buffer.data(ish + 72);
    const auto *ish_73 = buffer.data(ish + 73);
    const auto *ish_78 = buffer.data(ish + 78);
    const auto *ish_79 = buffer.data(ish + 79);
    const auto *ish_80 = buffer.data(ish + 80);
    const auto *ish_81 = buffer.data(ish + 81);
    const auto *ish_82 = buffer.data(ish + 82);
    const auto *ish_83 = buffer.data(ish + 83);
    const auto *ish_84 = buffer.data(ish + 84);
    const auto *ish_86 = buffer.data(ish + 86);
    const auto *ish_87 = buffer.data(ish + 87);
    const auto *ish_89 = buffer.data(ish + 89);
    const auto *ish_90 = buffer.data(ish + 90);
    const auto *ish_93 = buffer.data(ish + 93);
    const auto *ish_99 = buffer.data(ish + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsh_0, isg0_0, \
                         isg1_0, ish_0, ish_1, ish_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsh_0[k]
                 + f_1 * isg0_0[k]
                 - f_2 * isg1_0[k]
                 + f_3 * pc_x[k] * ish_0[k];

        t_1[k] = f_3 * pc_y[k] * ish_0[k];

        t_2[k] = f_3 * pc_z[k] * ish_0[k];

        t_3[k] = f_4 * isg0_0[k]
                 - f_5 * isg1_0[k]
                 + f_3 * pc_y[k] * ish_1[k];

        t_4[k] = f_3 * pc_y[k] * ish_2[k];

        t_5[k] = f_4 * isg0_0[k]
                 - f_5 * isg1_0[k]
                 + f_3 * pc_z[k] * ish_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, isg0_1, isg0_2, isg0_3, isg1_1, \
                         isg1_2, isg1_3, ish_3, ish_5, ish_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * isg0_1[k]
                 - f_7 * isg1_1[k]
                 + f_3 * pc_y[k] * ish_3[k];

        t_7[k] = f_3 * pc_z[k] * ish_3[k];

        t_8[k] = f_3 * pc_y[k] * ish_5[k];

        t_9[k] = f_6 * isg0_2[k]
                 - f_7 * isg1_2[k]
                 + f_3 * pc_z[k] * ish_5[k];

        t_10[k] = f_8 * isg0_3[k]
                  - f_9 * isg1_3[k]
                  + f_3 * pc_y[k] * ish_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, hsh_15, isg0_5, \
                         isg1_5, ish_6, ish_8, ish_9, ish_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ish_6[k];

        t_12[k] = f_4 * isg0_5[k]
                  - f_5 * isg1_5[k]
                  + f_3 * pc_y[k] * ish_8[k];

        t_13[k] = f_3 * pc_y[k] * ish_9[k];

        t_14[k] = f_8 * isg0_5[k]
                  - f_9 * isg1_5[k]
                  + f_3 * pc_z[k] * ish_9[k];

        t_15[k] = f_0 * hsh_15[k]
                  + f_3 * pc_x[k] * ish_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, hsh_17, hsh_18, \
                         hsh_20, ish_10, ish_14, ish_17, ish_18, \
                         ish_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * ish_10[k];

        t_17[k] = f_0 * hsh_17[k]
                  + f_3 * pc_x[k] * ish_17[k];

        t_18[k] = f_0 * hsh_18[k]
                  + f_3 * pc_x[k] * ish_18[k];

        t_19[k] = f_3 * pc_y[k] * ish_14[k];

        t_20[k] = f_0 * hsh_20[k]
                  + f_3 * pc_x[k] * ish_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, isg0_10, isg0_12, isg0_13, \
                         isg1_10, isg1_12, isg1_13, ish_15, ish_17, \
                         ish_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * isg0_10[k]
                  - f_2 * isg1_10[k]
                  + f_3 * pc_y[k] * ish_15[k];

        t_22[k] = f_3 * pc_z[k] * ish_15[k];

        t_23[k] = f_8 * isg0_12[k]
                  - f_9 * isg1_12[k]
                  + f_3 * pc_y[k] * ish_17[k];

        t_24[k] = f_6 * isg0_13[k]
                  - f_7 * isg1_13[k]
                  + f_3 * pc_y[k] * ish_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, hsi0_0, hsh_0, \
                         hsi1_0, isg0_14, isg1_14, ish_19, ish_20, \
                         ish_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * isg0_14[k]
                  - f_5 * isg1_14[k]
                  + f_3 * pc_y[k] * ish_19[k];

        t_26[k] = f_3 * pc_y[k] * ish_20[k];

        t_27[k] = f_1 * isg0_14[k]
                  - f_2 * isg1_14[k]
                  + f_3 * pc_z[k] * ish_20[k];

        t_28[k] = pa_y[k] * hsi0_0[k]
                  - f_10 * pc_y[k] * hsi1_0[k];

        t_29[k] = f_11 * hsh_0[k]
                  + f_3 * pc_y[k] * ish_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, hsi0_3, hsi0_5, hsh_1, \
                         hsi1_3, hsi1_5, ish_21, ish_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * ish_21[k];

        t_31[k] = pa_y[k] * hsi0_3[k]
                  + f_12 * hsh_1[k]
                  - f_10 * pc_y[k] * hsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * ish_22[k];

        t_33[k] = pa_y[k] * hsi0_5[k]
                  - f_10 * pc_y[k] * hsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, hsi0_6, hsi0_9, hsh_3, \
                         hsh_5, hsi1_6, hsi1_9, ish_24, ish_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * hsi0_6[k]
                  + f_13 * hsh_3[k]
                  - f_10 * pc_y[k] * hsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * ish_24[k];

        t_36[k] = f_11 * hsh_5[k]
                  + f_3 * pc_y[k] * ish_26[k];

        t_37[k] = pa_y[k] * hsi0_9[k]
                  - f_10 * pc_y[k] * hsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, hsi0_10, hsh_6, hsh_9, \
                         hsi1_10, isg0_18, isg1_18, ish_27, ish_28, \
                         ish_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * hsi0_10[k]
                  + f_14 * hsh_6[k]
                  - f_10 * pc_y[k] * hsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * ish_27[k];

        t_40[k] = f_4 * isg0_18[k]
                  - f_5 * isg1_18[k]
                  + f_3 * pc_z[k] * ish_28[k];

        t_41[k] = f_11 * hsh_9[k]
                  + f_3 * pc_y[k] * ish_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, hsi0_14, hsh_36, \
                         hsh_38, hsi1_14, ish_31, ish_36, ish_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * hsi0_14[k]
                  - f_10 * pc_y[k] * hsi1_14[k];

        t_43[k] = f_15 * hsh_36[k]
                  + f_3 * pc_x[k] * ish_36[k];

        t_44[k] = f_3 * pc_z[k] * ish_31[k];

        t_45[k] = f_15 * hsh_38[k]
                  + f_3 * pc_x[k] * ish_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, hsh_15, hsh_39, hsh_40, hsh_41, \
                         isg0_25, isg1_25, ish_36, ish_39, ish_40, \
                         ish_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * hsh_39[k]
                  + f_3 * pc_x[k] * ish_39[k];

        t_47[k] = f_15 * hsh_40[k]
                  + f_3 * pc_x[k] * ish_40[k];

        t_48[k] = f_15 * hsh_41[k]
                  + f_3 * pc_x[k] * ish_41[k];

        t_49[k] = f_11 * hsh_15[k]
                  + f_1 * isg0_25[k]
                  - f_2 * isg1_25[k]
                  + f_3 * pc_y[k] * ish_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, isg0_25, isg0_26, isg0_27, isg1_25, \
                         isg1_26, isg1_27, ish_36, ish_37, ish_38, \
                         ish_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * ish_36[k];

        t_51[k] = f_4 * isg0_25[k]
                  - f_5 * isg1_25[k]
                  + f_3 * pc_z[k] * ish_37[k];

        t_52[k] = f_6 * isg0_26[k]
                  - f_7 * isg1_26[k]
                  + f_3 * pc_z[k] * ish_38[k];

        t_53[k] = f_8 * isg0_27[k]
                  - f_9 * isg1_27[k]
                  + f_3 * pc_z[k] * ish_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, hsi0_0, hsi0_27, \
                         hsh_20, hsi1_0, hsi1_27, ish_41, ish_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * hsh_20[k]
                  + f_3 * pc_y[k] * ish_41[k];

        t_55[k] = pa_y[k] * hsi0_27[k]
                  - f_10 * pc_y[k] * hsi1_27[k];

        t_56[k] = pa_z[k] * hsi0_0[k]
                  - f_10 * pc_z[k] * hsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * ish_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, hsi0_3, hsi0_5, hsh_0, \
                         hsh_2, hsi1_3, hsi1_5, ish_42, ish_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * hsh_0[k]
                  + f_3 * pc_z[k] * ish_42[k];

        t_59[k] = pa_z[k] * hsi0_3[k]
                  - f_10 * pc_z[k] * hsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * ish_44[k];

        t_61[k] = pa_z[k] * hsi0_5[k]
                  + f_12 * hsh_2[k]
                  - f_10 * pc_z[k] * hsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, hsi0_6, hsi0_9, hsh_5, \
                         hsi1_6, hsi1_9, isg0_32, isg1_32, ish_46, \
                         ish_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * hsi0_6[k]
                  - f_10 * pc_z[k] * hsi1_6[k];

        t_63[k] = f_4 * isg0_32[k]
                  - f_5 * isg1_32[k]
                  + f_3 * pc_y[k] * ish_46[k];

        t_64[k] = f_3 * pc_y[k] * ish_47[k];

        t_65[k] = pa_z[k] * hsi0_9[k]
                  + f_13 * hsh_5[k]
                  - f_10 * pc_z[k] * hsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, hsi0_10, hsi1_10, isg0_34, \
                         isg0_35, isg1_34, isg1_35, ish_49, ish_50, \
                         ish_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * hsi0_10[k]
                  - f_10 * pc_z[k] * hsi1_10[k];

        t_67[k] = f_6 * isg0_34[k]
                  - f_7 * isg1_34[k]
                  + f_3 * pc_y[k] * ish_49[k];

        t_68[k] = f_4 * isg0_35[k]
                  - f_5 * isg1_35[k]
                  + f_3 * pc_y[k] * ish_50[k];

        t_69[k] = f_3 * pc_y[k] * ish_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, hsi0_14, hsh_9, hsh_57, \
                         hsh_58, hsh_59, hsi1_14, ish_57, ish_58, \
                         ish_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * hsi0_14[k]
                  + f_14 * hsh_9[k]
                  - f_10 * pc_z[k] * hsi1_14[k];

        t_71[k] = f_15 * hsh_57[k]
                  + f_3 * pc_x[k] * ish_57[k];

        t_72[k] = f_15 * hsh_58[k]
                  + f_3 * pc_x[k] * ish_58[k];

        t_73[k] = f_15 * hsh_59[k]
                  + f_3 * pc_x[k] * ish_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, hsi0_21, hsh_60, \
                         hsh_62, hsi1_21, ish_56, ish_60, ish_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * hsh_60[k]
                  + f_3 * pc_x[k] * ish_60[k];

        t_75[k] = f_3 * pc_y[k] * ish_56[k];

        t_76[k] = f_15 * hsh_62[k]
                  + f_3 * pc_x[k] * ish_62[k];

        t_77[k] = pa_z[k] * hsi0_21[k]
                  - f_10 * pc_z[k] * hsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, isg0_41, isg0_42, isg0_43, isg1_41, isg1_42, \
                         isg1_43, ish_58, ish_59, ish_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * isg0_41[k]
                  - f_17 * isg1_41[k]
                  + f_3 * pc_y[k] * ish_58[k];

        t_79[k] = f_8 * isg0_42[k]
                  - f_9 * isg1_42[k]
                  + f_3 * pc_y[k] * ish_59[k];

        t_80[k] = f_6 * isg0_43[k]
                  - f_7 * isg1_43[k]
                  + f_3 * pc_y[k] * ish_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, hsh_20, hsh_63, isg0_44, \
                         isg0_45, isg1_44, isg1_45, ish_61, ish_62, \
                         ish_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * isg0_44[k]
                  - f_5 * isg1_44[k]
                  + f_3 * pc_y[k] * ish_61[k];

        t_82[k] = f_3 * pc_y[k] * ish_62[k];

        t_83[k] = f_11 * hsh_20[k]
                  + f_1 * isg0_44[k]
                  - f_2 * isg1_44[k]
                  + f_3 * pc_z[k] * ish_62[k];

        t_84[k] = f_14 * hsh_63[k]
                  + f_1 * isg0_45[k]
                  - f_2 * isg1_45[k]
                  + f_3 * pc_x[k] * ish_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, hsh_21, hsh_66, isg0_48, \
                         isg1_48, ish_63, ish_64, ish_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * hsh_21[k]
                  + f_3 * pc_y[k] * ish_63[k];

        t_86[k] = f_3 * pc_z[k] * ish_63[k];

        t_87[k] = f_14 * hsh_66[k]
                  + f_8 * isg0_48[k]
                  - f_9 * isg1_48[k]
                  + f_3 * pc_x[k] * ish_66[k];

        t_88[k] = f_3 * pc_z[k] * ish_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, hsh_69, isg0_45, isg0_51, isg1_45, \
                         isg1_51, ish_65, ish_66, ish_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * isg0_45[k]
                  - f_5 * isg1_45[k]
                  + f_3 * pc_z[k] * ish_65[k];

        t_90[k] = f_14 * hsh_69[k]
                  + f_6 * isg0_51[k]
                  - f_7 * isg1_51[k]
                  + f_3 * pc_x[k] * ish_69[k];

        t_91[k] = f_3 * pc_z[k] * ish_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, hsh_26, hsh_73, isg0_47, \
                         isg0_55, isg1_47, isg1_55, ish_68, ish_69, \
                         ish_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * hsh_26[k]
                  + f_3 * pc_y[k] * ish_68[k];

        t_93[k] = f_6 * isg0_47[k]
                  - f_7 * isg1_47[k]
                  + f_3 * pc_z[k] * ish_68[k];

        t_94[k] = f_14 * hsh_73[k]
                  + f_4 * isg0_55[k]
                  - f_5 * isg1_55[k]
                  + f_3 * pc_x[k] * ish_73[k];

        t_95[k] = f_3 * pc_z[k] * ish_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, hsh_30, hsh_78, isg0_48, \
                         isg0_50, isg1_48, isg1_50, ish_70, ish_72, \
                         ish_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * isg0_48[k]
                  - f_5 * isg1_48[k]
                  + f_3 * pc_z[k] * ish_70[k];

        t_97[k] = f_12 * hsh_30[k]
                  + f_3 * pc_y[k] * ish_72[k];

        t_98[k] = f_8 * isg0_50[k]
                  - f_9 * isg1_50[k]
                  + f_3 * pc_z[k] * ish_72[k];

        t_99[k] = f_14 * hsh_78[k]
                  + f_3 * pc_x[k] * ish_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, hsh_80, hsh_81, \
                         hsh_82, hsh_83, ish_73, ish_80, ish_81, ish_82, \
                         ish_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * ish_73[k];

        t_101[k] = f_14 * hsh_80[k]
                   + f_3 * pc_x[k] * ish_80[k];

        t_102[k] = f_14 * hsh_81[k]
                   + f_3 * pc_x[k] * ish_81[k];

        t_103[k] = f_14 * hsh_82[k]
                   + f_3 * pc_x[k] * ish_82[k];

        t_104[k] = f_14 * hsh_83[k]
                   + f_3 * pc_x[k] * ish_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, hsh_36, isg0_55, isg0_56, \
                         isg1_55, isg1_56, ish_78, ish_79, ish_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * hsh_36[k]
                   + f_1 * isg0_55[k]
                   - f_2 * isg1_55[k]
                   + f_3 * pc_y[k] * ish_78[k];

        t_106[k] = f_3 * pc_z[k] * ish_78[k];

        t_107[k] = f_4 * isg0_55[k]
                   - f_5 * isg1_55[k]
                   + f_3 * pc_z[k] * ish_79[k];

        t_108[k] = f_6 * isg0_56[k]
                   - f_7 * isg1_56[k]
                   + f_3 * pc_z[k] * ish_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, hsi0_56, hsh_41, \
                         hsi1_56, isg0_57, isg0_59, isg1_57, isg1_59, ish_81, \
                         ish_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * isg0_57[k]
                   - f_9 * isg1_57[k]
                   + f_3 * pc_z[k] * ish_81[k];

        t_110[k] = f_12 * hsh_41[k]
                   + f_3 * pc_y[k] * ish_83[k];

        t_111[k] = f_1 * isg0_59[k]
                   - f_2 * isg1_59[k]
                   + f_3 * pc_z[k] * ish_83[k];

        t_112[k] = pa_y[k] * hsi0_56[k]
                   - f_10 * pc_y[k] * hsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, hsi0_31, hsh_21, \
                         hsh_42, hsh_44, hsi1_31, ish_84, ish_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * hsh_42[k]
                   + f_3 * pc_y[k] * ish_84[k];

        t_114[k] = f_11 * hsh_21[k]
                   + f_3 * pc_z[k] * ish_84[k];

        t_115[k] = pa_z[k] * hsi0_31[k]
                   - f_10 * pc_z[k] * hsi1_31[k];

        t_116[k] = f_11 * hsh_44[k]
                   + f_3 * pc_y[k] * ish_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, hsi0_34, hsi0_61, \
                         hsh_24, hsh_47, hsi1_34, hsi1_61, ish_87, \
                         ish_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * hsi0_61[k]
                   - f_10 * pc_y[k] * hsi1_61[k];

        t_118[k] = pa_z[k] * hsi0_34[k]
                   - f_10 * pc_z[k] * hsi1_34[k];

        t_119[k] = f_11 * hsh_24[k]
                   + f_3 * pc_z[k] * ish_87[k];

        t_120[k] = f_11 * hsh_47[k]
                   + f_3 * pc_y[k] * ish_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, hsi0_38, hsi0_65, \
                         hsh_27, hsi1_38, hsi1_65, ish_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * hsi0_65[k]
                   - f_10 * pc_y[k] * hsi1_65[k];

        t_122[k] = pa_z[k] * hsi0_38[k]
                   - f_10 * pc_z[k] * hsi1_38[k];

        t_123[k] = f_11 * hsh_27[k]
                   + f_3 * pc_z[k] * ish_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, hsi0_68, hsi0_70, \
                         hsh_50, hsh_51, hsh_99, hsi1_68, hsi1_70, ish_93, \
                         ish_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * hsi0_68[k]
                   + f_12 * hsh_50[k]
                   - f_10 * pc_y[k] * hsi1_68[k];

        t_125[k] = f_11 * hsh_51[k]
                   + f_3 * pc_y[k] * ish_93[k];

        t_126[k] = pa_y[k] * hsi0_70[k]
                   - f_10 * pc_y[k] * hsi1_70[k];

        t_127[k] = f_14 * hsh_99[k]
                   + f_3 * pc_x[k] * ish_99[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
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

    const auto *hsi0_49 = buffer.data(hsi0 + 49);
    const auto *hsi0_83 = buffer.data(hsi0 + 83);
    const auto *hsi0_84 = buffer.data(hsi0 + 84);
    const auto *hsi0_87 = buffer.data(hsi0 + 87);
    const auto *hsi0_90 = buffer.data(hsi0 + 90);
    const auto *hsi0_94 = buffer.data(hsi0 + 94);
    const auto *hsi0_96 = buffer.data(hsi0 + 96);
    const auto *hsi0_105 = buffer.data(hsi0 + 105);
    const auto *hsi0_140 = buffer.data(hsi0 + 140);
    const auto *hsi0_143 = buffer.data(hsi0 + 143);
    const auto *hsi0_145 = buffer.data(hsi0 + 145);
    const auto *hsi0_146 = buffer.data(hsi0 + 146);
    const auto *hsi0_149 = buffer.data(hsi0 + 149);
    const auto *hsi0_150 = buffer.data(hsi0 + 150);
    const auto *hsi0_152 = buffer.data(hsi0 + 152);
    const auto *hsi0_154 = buffer.data(hsi0 + 154);

    const auto *hsh_36 = buffer.data(hsh + 36);
    const auto *hsh_42 = buffer.data(hsh + 42);
    const auto *hsh_59 = buffer.data(hsh + 59);
    const auto *hsh_60 = buffer.data(hsh + 60);
    const auto *hsh_61 = buffer.data(hsh + 61);
    const auto *hsh_62 = buffer.data(hsh + 62);
    const auto *hsh_63 = buffer.data(hsh + 63);
    const auto *hsh_66 = buffer.data(hsh + 66);
    const auto *hsh_68 = buffer.data(hsh + 68);
    const auto *hsh_69 = buffer.data(hsh + 69);
    const auto *hsh_70 = buffer.data(hsh + 70);
    const auto *hsh_72 = buffer.data(hsh + 72);
    const auto *hsh_78 = buffer.data(hsh + 78);
    const auto *hsh_83 = buffer.data(hsh + 83);
    const auto *hsh_84 = buffer.data(hsh + 84);
    const auto *hsh_86 = buffer.data(hsh + 86);
    const auto *hsh_87 = buffer.data(hsh + 87);
    const auto *hsh_89 = buffer.data(hsh + 89);
    const auto *hsh_90 = buffer.data(hsh + 90);
    const auto *hsh_93 = buffer.data(hsh + 93);
    const auto *hsh_99 = buffer.data(hsh + 99);
    const auto *hsh_100 = buffer.data(hsh + 100);
    const auto *hsh_101 = buffer.data(hsh + 101);
    const auto *hsh_102 = buffer.data(hsh + 102);
    const auto *hsh_103 = buffer.data(hsh + 103);
    const auto *hsh_104 = buffer.data(hsh + 104);
    const auto *hsh_105 = buffer.data(hsh + 105);
    const auto *hsh_106 = buffer.data(hsh + 106);
    const auto *hsh_107 = buffer.data(hsh + 107);
    const auto *hsh_108 = buffer.data(hsh + 108);
    const auto *hsh_110 = buffer.data(hsh + 110);
    const auto *hsh_111 = buffer.data(hsh + 111);
    const auto *hsh_113 = buffer.data(hsh + 113);
    const auto *hsh_114 = buffer.data(hsh + 114);
    const auto *hsh_119 = buffer.data(hsh + 119);
    const auto *hsh_120 = buffer.data(hsh + 120);
    const auto *hsh_121 = buffer.data(hsh + 121);
    const auto *hsh_122 = buffer.data(hsh + 122);
    const auto *hsh_123 = buffer.data(hsh + 123);
    const auto *hsh_125 = buffer.data(hsh + 125);
    const auto *hsh_126 = buffer.data(hsh + 126);
    const auto *hsh_129 = buffer.data(hsh + 129);
    const auto *hsh_132 = buffer.data(hsh + 132);
    const auto *hsh_136 = buffer.data(hsh + 136);
    const auto *hsh_141 = buffer.data(hsh + 141);
    const auto *hsh_143 = buffer.data(hsh + 143);
    const auto *hsh_144 = buffer.data(hsh + 144);
    const auto *hsh_145 = buffer.data(hsh + 145);
    const auto *hsh_146 = buffer.data(hsh + 146);
    const auto *hsh_152 = buffer.data(hsh + 152);
    const auto *hsh_156 = buffer.data(hsh + 156);
    const auto *hsh_161 = buffer.data(hsh + 161);
    const auto *hsh_162 = buffer.data(hsh + 162);
    const auto *hsh_163 = buffer.data(hsh + 163);
    const auto *hsh_164 = buffer.data(hsh + 164);
    const auto *hsh_165 = buffer.data(hsh + 165);
    const auto *hsh_166 = buffer.data(hsh + 166);
    const auto *hsh_167 = buffer.data(hsh + 167);
    const auto *hsh_183 = buffer.data(hsh + 183);
    const auto *hsh_184 = buffer.data(hsh + 184);
    const auto *hsh_185 = buffer.data(hsh + 185);
    const auto *hsh_186 = buffer.data(hsh + 186);
    const auto *hsh_187 = buffer.data(hsh + 187);
    const auto *hsh_188 = buffer.data(hsh + 188);

    const auto *hsi1_49 = buffer.data(hsi1 + 49);
    const auto *hsi1_83 = buffer.data(hsi1 + 83);
    const auto *hsi1_84 = buffer.data(hsi1 + 84);
    const auto *hsi1_87 = buffer.data(hsi1 + 87);
    const auto *hsi1_90 = buffer.data(hsi1 + 90);
    const auto *hsi1_94 = buffer.data(hsi1 + 94);
    const auto *hsi1_96 = buffer.data(hsi1 + 96);
    const auto *hsi1_105 = buffer.data(hsi1 + 105);
    const auto *hsi1_140 = buffer.data(hsi1 + 140);
    const auto *hsi1_143 = buffer.data(hsi1 + 143);
    const auto *hsi1_145 = buffer.data(hsi1 + 145);
    const auto *hsi1_146 = buffer.data(hsi1 + 146);
    const auto *hsi1_149 = buffer.data(hsi1 + 149);
    const auto *hsi1_150 = buffer.data(hsi1 + 150);
    const auto *hsi1_152 = buffer.data(hsi1 + 152);
    const auto *hsi1_154 = buffer.data(hsi1 + 154);

    const auto *isg0_72 = buffer.data(isg0 + 72);
    const auto *isg0_73 = buffer.data(isg0 + 73);
    const auto *isg0_74 = buffer.data(isg0 + 74);
    const auto *isg0_75 = buffer.data(isg0 + 75);
    const auto *isg0_76 = buffer.data(isg0 + 76);
    const auto *isg0_77 = buffer.data(isg0 + 77);
    const auto *isg0_78 = buffer.data(isg0 + 78);
    const auto *isg0_79 = buffer.data(isg0 + 79);
    const auto *isg0_80 = buffer.data(isg0 + 80);
    const auto *isg0_84 = buffer.data(isg0 + 84);
    const auto *isg0_85 = buffer.data(isg0 + 85);
    const auto *isg0_86 = buffer.data(isg0 + 86);
    const auto *isg0_87 = buffer.data(isg0 + 87);
    const auto *isg0_88 = buffer.data(isg0 + 88);
    const auto *isg0_89 = buffer.data(isg0 + 89);
    const auto *isg0_90 = buffer.data(isg0 + 90);
    const auto *isg0_92 = buffer.data(isg0 + 92);
    const auto *isg0_93 = buffer.data(isg0 + 93);
    const auto *isg0_95 = buffer.data(isg0 + 95);
    const auto *isg0_96 = buffer.data(isg0 + 96);
    const auto *isg0_100 = buffer.data(isg0 + 100);
    const auto *isg0_101 = buffer.data(isg0 + 101);
    const auto *isg0_102 = buffer.data(isg0 + 102);
    const auto *isg0_104 = buffer.data(isg0 + 104);
    const auto *isg0_110 = buffer.data(isg0 + 110);
    const auto *isg0_114 = buffer.data(isg0 + 114);
    const auto *isg0_117 = buffer.data(isg0 + 117);
    const auto *isg0_118 = buffer.data(isg0 + 118);
    const auto *isg0_119 = buffer.data(isg0 + 119);
    const auto *isg0_130 = buffer.data(isg0 + 130);
    const auto *isg0_132 = buffer.data(isg0 + 132);

    const auto *isg1_72 = buffer.data(isg1 + 72);
    const auto *isg1_73 = buffer.data(isg1 + 73);
    const auto *isg1_74 = buffer.data(isg1 + 74);
    const auto *isg1_75 = buffer.data(isg1 + 75);
    const auto *isg1_76 = buffer.data(isg1 + 76);
    const auto *isg1_77 = buffer.data(isg1 + 77);
    const auto *isg1_78 = buffer.data(isg1 + 78);
    const auto *isg1_79 = buffer.data(isg1 + 79);
    const auto *isg1_80 = buffer.data(isg1 + 80);
    const auto *isg1_84 = buffer.data(isg1 + 84);
    const auto *isg1_85 = buffer.data(isg1 + 85);
    const auto *isg1_86 = buffer.data(isg1 + 86);
    const auto *isg1_87 = buffer.data(isg1 + 87);
    const auto *isg1_88 = buffer.data(isg1 + 88);
    const auto *isg1_89 = buffer.data(isg1 + 89);
    const auto *isg1_90 = buffer.data(isg1 + 90);
    const auto *isg1_92 = buffer.data(isg1 + 92);
    const auto *isg1_93 = buffer.data(isg1 + 93);
    const auto *isg1_95 = buffer.data(isg1 + 95);
    const auto *isg1_96 = buffer.data(isg1 + 96);
    const auto *isg1_100 = buffer.data(isg1 + 100);
    const auto *isg1_101 = buffer.data(isg1 + 101);
    const auto *isg1_102 = buffer.data(isg1 + 102);
    const auto *isg1_104 = buffer.data(isg1 + 104);
    const auto *isg1_110 = buffer.data(isg1 + 110);
    const auto *isg1_114 = buffer.data(isg1 + 114);
    const auto *isg1_117 = buffer.data(isg1 + 117);
    const auto *isg1_118 = buffer.data(isg1 + 118);
    const auto *isg1_119 = buffer.data(isg1 + 119);
    const auto *isg1_130 = buffer.data(isg1 + 130);
    const auto *isg1_132 = buffer.data(isg1 + 132);

    const auto *ish_99 = buffer.data(ish + 99);
    const auto *ish_100 = buffer.data(ish + 100);
    const auto *ish_101 = buffer.data(ish + 101);
    const auto *ish_102 = buffer.data(ish + 102);
    const auto *ish_103 = buffer.data(ish + 103);
    const auto *ish_104 = buffer.data(ish + 104);
    const auto *ish_105 = buffer.data(ish + 105);
    const auto *ish_106 = buffer.data(ish + 106);
    const auto *ish_107 = buffer.data(ish + 107);
    const auto *ish_108 = buffer.data(ish + 108);
    const auto *ish_109 = buffer.data(ish + 109);
    const auto *ish_110 = buffer.data(ish + 110);
    const auto *ish_111 = buffer.data(ish + 111);
    const auto *ish_112 = buffer.data(ish + 112);
    const auto *ish_113 = buffer.data(ish + 113);
    const auto *ish_114 = buffer.data(ish + 114);
    const auto *ish_119 = buffer.data(ish + 119);
    const auto *ish_120 = buffer.data(ish + 120);
    const auto *ish_121 = buffer.data(ish + 121);
    const auto *ish_122 = buffer.data(ish + 122);
    const auto *ish_123 = buffer.data(ish + 123);
    const auto *ish_124 = buffer.data(ish + 124);
    const auto *ish_125 = buffer.data(ish + 125);
    const auto *ish_126 = buffer.data(ish + 126);
    const auto *ish_127 = buffer.data(ish + 127);
    const auto *ish_128 = buffer.data(ish + 128);
    const auto *ish_129 = buffer.data(ish + 129);
    const auto *ish_131 = buffer.data(ish + 131);
    const auto *ish_132 = buffer.data(ish + 132);
    const auto *ish_133 = buffer.data(ish + 133);
    const auto *ish_135 = buffer.data(ish + 135);
    const auto *ish_136 = buffer.data(ish + 136);
    const auto *ish_141 = buffer.data(ish + 141);
    const auto *ish_142 = buffer.data(ish + 142);
    const auto *ish_143 = buffer.data(ish + 143);
    const auto *ish_144 = buffer.data(ish + 144);
    const auto *ish_145 = buffer.data(ish + 145);
    const auto *ish_146 = buffer.data(ish + 146);
    const auto *ish_147 = buffer.data(ish + 147);
    const auto *ish_149 = buffer.data(ish + 149);
    const auto *ish_150 = buffer.data(ish + 150);
    const auto *ish_152 = buffer.data(ish + 152);
    const auto *ish_153 = buffer.data(ish + 153);
    const auto *ish_156 = buffer.data(ish + 156);
    const auto *ish_161 = buffer.data(ish + 161);
    const auto *ish_162 = buffer.data(ish + 162);
    const auto *ish_163 = buffer.data(ish + 163);
    const auto *ish_164 = buffer.data(ish + 164);
    const auto *ish_165 = buffer.data(ish + 165);
    const auto *ish_166 = buffer.data(ish + 166);
    const auto *ish_167 = buffer.data(ish + 167);
    const auto *ish_168 = buffer.data(ish + 168);
    const auto *ish_170 = buffer.data(ish + 170);
    const auto *ish_171 = buffer.data(ish + 171);
    const auto *ish_173 = buffer.data(ish + 173);
    const auto *ish_174 = buffer.data(ish + 174);
    const auto *ish_177 = buffer.data(ish + 177);
    const auto *ish_183 = buffer.data(ish + 183);
    const auto *ish_184 = buffer.data(ish + 184);
    const auto *ish_185 = buffer.data(ish + 185);
    const auto *ish_186 = buffer.data(ish + 186);
    const auto *ish_187 = buffer.data(ish + 187);
    const auto *ish_188 = buffer.data(ish + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, hsh_100, hsh_101, hsh_102, \
                         hsh_103, hsh_104, ish_100, ish_101, ish_102, ish_103, \
                         ish_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_14 * hsh_100[k]
                   + f_3 * pc_x[k] * ish_100[k];

        t_129[k] = f_14 * hsh_101[k]
                   + f_3 * pc_x[k] * ish_101[k];

        t_130[k] = f_14 * hsh_102[k]
                   + f_3 * pc_x[k] * ish_102[k];

        t_131[k] = f_14 * hsh_103[k]
                   + f_3 * pc_x[k] * ish_103[k];

        t_132[k] = f_14 * hsh_104[k]
                   + f_3 * pc_x[k] * ish_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, hsi0_49, hsh_36, hsh_59, \
                         hsi1_49, isg0_72, isg1_72, ish_99, ish_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * hsi0_49[k]
                   - f_10 * pc_z[k] * hsi1_49[k];

        t_134[k] = f_11 * hsh_36[k]
                   + f_3 * pc_z[k] * ish_99[k];

        t_135[k] = f_11 * hsh_59[k]
                   + f_8 * isg0_72[k]
                   - f_9 * isg1_72[k]
                   + f_3 * pc_y[k] * ish_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, hsh_60, hsh_61, hsh_62, isg0_73, isg0_74, \
                         isg1_73, isg1_74, ish_102, ish_103, ish_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * hsh_60[k]
                   + f_6 * isg0_73[k]
                   - f_7 * isg1_73[k]
                   + f_3 * pc_y[k] * ish_102[k];

        t_137[k] = f_11 * hsh_61[k]
                   + f_4 * isg0_74[k]
                   - f_5 * isg1_74[k]
                   + f_3 * pc_y[k] * ish_103[k];

        t_138[k] = f_11 * hsh_62[k]
                   + f_3 * pc_y[k] * ish_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, hsi0_83, hsh_42, \
                         hsh_105, hsi1_83, isg0_75, isg1_75, ish_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * hsi0_83[k]
                   - f_10 * pc_y[k] * hsi1_83[k];

        t_140[k] = f_14 * hsh_105[k]
                   + f_1 * isg0_75[k]
                   - f_2 * isg1_75[k]
                   + f_3 * pc_x[k] * ish_105[k];

        t_141[k] = f_3 * pc_y[k] * ish_105[k];

        t_142[k] = f_12 * hsh_42[k]
                   + f_3 * pc_z[k] * ish_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, hsh_110, isg0_75, isg0_80, isg1_75, \
                         isg1_80, ish_106, ish_107, ish_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * isg0_75[k]
                   - f_5 * isg1_75[k]
                   + f_3 * pc_y[k] * ish_106[k];

        t_144[k] = f_3 * pc_y[k] * ish_107[k];

        t_145[k] = f_14 * hsh_110[k]
                   + f_8 * isg0_80[k]
                   - f_9 * isg1_80[k]
                   + f_3 * pc_x[k] * ish_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, isg0_76, isg0_77, isg1_76, isg1_77, \
                         ish_108, ish_109, ish_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * isg0_76[k]
                   - f_7 * isg1_76[k]
                   + f_3 * pc_y[k] * ish_108[k];

        t_147[k] = f_4 * isg0_77[k]
                   - f_5 * isg1_77[k]
                   + f_3 * pc_y[k] * ish_109[k];

        t_148[k] = f_3 * pc_y[k] * ish_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, hsh_114, isg0_78, isg0_79, isg0_84, \
                         isg1_78, isg1_79, isg1_84, ish_111, ish_112, \
                         ish_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * hsh_114[k]
                   + f_6 * isg0_84[k]
                   - f_7 * isg1_84[k]
                   + f_3 * pc_x[k] * ish_114[k];

        t_150[k] = f_8 * isg0_78[k]
                   - f_9 * isg1_78[k]
                   + f_3 * pc_y[k] * ish_111[k];

        t_151[k] = f_6 * isg0_79[k]
                   - f_7 * isg1_79[k]
                   + f_3 * pc_y[k] * ish_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, hsh_119, hsh_120, isg0_80, \
                         isg0_89, isg1_80, isg1_89, ish_113, ish_114, ish_119, \
                         ish_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * isg0_80[k]
                   - f_5 * isg1_80[k]
                   + f_3 * pc_y[k] * ish_113[k];

        t_153[k] = f_3 * pc_y[k] * ish_114[k];

        t_154[k] = f_14 * hsh_119[k]
                   + f_4 * isg0_89[k]
                   - f_5 * isg1_89[k]
                   + f_3 * pc_x[k] * ish_119[k];

        t_155[k] = f_14 * hsh_120[k]
                   + f_3 * pc_x[k] * ish_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, hsh_121, hsh_122, \
                         hsh_123, hsh_125, ish_119, ish_121, ish_122, ish_123, \
                         ish_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_14 * hsh_121[k]
                   + f_3 * pc_x[k] * ish_121[k];

        t_157[k] = f_14 * hsh_122[k]
                   + f_3 * pc_x[k] * ish_122[k];

        t_158[k] = f_14 * hsh_123[k]
                   + f_3 * pc_x[k] * ish_123[k];

        t_159[k] = f_3 * pc_y[k] * ish_119[k];

        t_160[k] = f_14 * hsh_125[k]
                   + f_3 * pc_x[k] * ish_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, isg0_85, isg0_86, isg0_87, isg1_85, \
                         isg1_86, isg1_87, ish_120, ish_121, ish_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * isg0_85[k]
                   - f_2 * isg1_85[k]
                   + f_3 * pc_y[k] * ish_120[k];

        t_162[k] = f_16 * isg0_86[k]
                   - f_17 * isg1_86[k]
                   + f_3 * pc_y[k] * ish_121[k];

        t_163[k] = f_8 * isg0_87[k]
                   - f_9 * isg1_87[k]
                   + f_3 * pc_y[k] * ish_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, hsh_62, isg0_88, isg0_89, \
                         isg1_88, isg1_89, ish_123, ish_124, ish_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * isg0_88[k]
                   - f_7 * isg1_88[k]
                   + f_3 * pc_y[k] * ish_123[k];

        t_165[k] = f_4 * isg0_89[k]
                   - f_5 * isg1_89[k]
                   + f_3 * pc_y[k] * ish_124[k];

        t_166[k] = f_3 * pc_y[k] * ish_125[k];

        t_167[k] = f_12 * hsh_62[k]
                   + f_1 * isg0_89[k]
                   - f_2 * isg1_89[k]
                   + f_3 * pc_z[k] * ish_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, hsh_63, hsh_126, \
                         hsh_129, isg0_90, isg0_93, isg1_90, isg1_93, ish_126, \
                         ish_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_13 * hsh_126[k]
                   + f_1 * isg0_90[k]
                   - f_2 * isg1_90[k]
                   + f_3 * pc_x[k] * ish_126[k];

        t_169[k] = f_13 * hsh_63[k]
                   + f_3 * pc_y[k] * ish_126[k];

        t_170[k] = f_3 * pc_z[k] * ish_126[k];

        t_171[k] = f_13 * hsh_129[k]
                   + f_8 * isg0_93[k]
                   - f_9 * isg1_93[k]
                   + f_3 * pc_x[k] * ish_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, hsh_132, isg0_90, isg0_96, \
                         isg1_90, isg1_96, ish_127, ish_128, ish_129, \
                         ish_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * ish_127[k];

        t_173[k] = f_4 * isg0_90[k]
                   - f_5 * isg1_90[k]
                   + f_3 * pc_z[k] * ish_128[k];

        t_174[k] = f_13 * hsh_132[k]
                   + f_6 * isg0_96[k]
                   - f_7 * isg1_96[k]
                   + f_3 * pc_x[k] * ish_132[k];

        t_175[k] = f_3 * pc_z[k] * ish_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, hsh_68, hsh_136, \
                         isg0_92, isg0_100, isg1_92, isg1_100, ish_131, ish_132, \
                         ish_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * hsh_68[k]
                   + f_3 * pc_y[k] * ish_131[k];

        t_177[k] = f_6 * isg0_92[k]
                   - f_7 * isg1_92[k]
                   + f_3 * pc_z[k] * ish_131[k];

        t_178[k] = f_13 * hsh_136[k]
                   + f_4 * isg0_100[k]
                   - f_5 * isg1_100[k]
                   + f_3 * pc_x[k] * ish_136[k];

        t_179[k] = f_3 * pc_z[k] * ish_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, hsh_72, hsh_141, \
                         isg0_93, isg0_95, isg1_93, isg1_95, ish_133, ish_135, \
                         ish_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * isg0_93[k]
                   - f_5 * isg1_93[k]
                   + f_3 * pc_z[k] * ish_133[k];

        t_181[k] = f_13 * hsh_72[k]
                   + f_3 * pc_y[k] * ish_135[k];

        t_182[k] = f_8 * isg0_95[k]
                   - f_9 * isg1_95[k]
                   + f_3 * pc_z[k] * ish_135[k];

        t_183[k] = f_13 * hsh_141[k]
                   + f_3 * pc_x[k] * ish_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, hsh_143, hsh_144, \
                         hsh_145, hsh_146, ish_136, ish_143, ish_144, ish_145, \
                         ish_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * ish_136[k];

        t_185[k] = f_13 * hsh_143[k]
                   + f_3 * pc_x[k] * ish_143[k];

        t_186[k] = f_13 * hsh_144[k]
                   + f_3 * pc_x[k] * ish_144[k];

        t_187[k] = f_13 * hsh_145[k]
                   + f_3 * pc_x[k] * ish_145[k];

        t_188[k] = f_13 * hsh_146[k]
                   + f_3 * pc_x[k] * ish_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, hsh_78, isg0_100, isg0_101, \
                         isg1_100, isg1_101, ish_141, ish_142, \
                         ish_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * hsh_78[k]
                   + f_1 * isg0_100[k]
                   - f_2 * isg1_100[k]
                   + f_3 * pc_y[k] * ish_141[k];

        t_190[k] = f_3 * pc_z[k] * ish_141[k];

        t_191[k] = f_4 * isg0_100[k]
                   - f_5 * isg1_100[k]
                   + f_3 * pc_z[k] * ish_142[k];

        t_192[k] = f_6 * isg0_101[k]
                   - f_7 * isg1_101[k]
                   + f_3 * pc_z[k] * ish_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, hsi0_84, hsh_83, \
                         hsi1_84, isg0_102, isg0_104, isg1_102, isg1_104, ish_144, \
                         ish_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * isg0_102[k]
                   - f_9 * isg1_102[k]
                   + f_3 * pc_z[k] * ish_144[k];

        t_194[k] = f_13 * hsh_83[k]
                   + f_3 * pc_y[k] * ish_146[k];

        t_195[k] = f_1 * isg0_104[k]
                   - f_2 * isg1_104[k]
                   + f_3 * pc_z[k] * ish_146[k];

        t_196[k] = pa_z[k] * hsi0_84[k]
                   - f_10 * pc_z[k] * hsi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, hsi0_87, hsh_63, \
                         hsh_84, hsh_86, hsi1_87, ish_147, ish_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * hsh_84[k]
                   + f_3 * pc_y[k] * ish_147[k];

        t_198[k] = f_11 * hsh_63[k]
                   + f_3 * pc_z[k] * ish_147[k];

        t_199[k] = pa_z[k] * hsi0_87[k]
                   - f_10 * pc_z[k] * hsi1_87[k];

        t_200[k] = f_12 * hsh_86[k]
                   + f_3 * pc_y[k] * ish_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, hsi0_90, hsh_66, hsh_152, \
                         hsi1_90, isg0_110, isg1_110, ish_150, \
                         ish_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_13 * hsh_152[k]
                   + f_8 * isg0_110[k]
                   - f_9 * isg1_110[k]
                   + f_3 * pc_x[k] * ish_152[k];

        t_202[k] = pa_z[k] * hsi0_90[k]
                   - f_10 * pc_z[k] * hsi1_90[k];

        t_203[k] = f_11 * hsh_66[k]
                   + f_3 * pc_z[k] * ish_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, hsi0_94, hsh_89, \
                         hsh_156, hsi1_94, isg0_114, isg1_114, ish_152, \
                         ish_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * hsh_89[k]
                   + f_3 * pc_y[k] * ish_152[k];

        t_205[k] = f_13 * hsh_156[k]
                   + f_6 * isg0_114[k]
                   - f_7 * isg1_114[k]
                   + f_3 * pc_x[k] * ish_156[k];

        t_206[k] = pa_z[k] * hsi0_94[k]
                   - f_10 * pc_z[k] * hsi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, hsi0_96, hsh_69, hsh_70, \
                         hsh_93, hsi1_96, ish_153, ish_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * hsh_69[k]
                   + f_3 * pc_z[k] * ish_153[k];

        t_208[k] = pa_z[k] * hsi0_96[k]
                   + f_12 * hsh_70[k]
                   - f_10 * pc_z[k] * hsi1_96[k];

        t_209[k] = f_12 * hsh_93[k]
                   + f_3 * pc_y[k] * ish_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, hsh_161, hsh_162, hsh_163, hsh_164, \
                         isg0_119, isg1_119, ish_161, ish_162, ish_163, \
                         ish_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_13 * hsh_161[k]
                   + f_4 * isg0_119[k]
                   - f_5 * isg1_119[k]
                   + f_3 * pc_x[k] * ish_161[k];

        t_211[k] = f_13 * hsh_162[k]
                   + f_3 * pc_x[k] * ish_162[k];

        t_212[k] = f_13 * hsh_163[k]
                   + f_3 * pc_x[k] * ish_163[k];

        t_213[k] = f_13 * hsh_164[k]
                   + f_3 * pc_x[k] * ish_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, hsi0_105, hsh_165, \
                         hsh_166, hsh_167, hsi1_105, ish_165, ish_166, \
                         ish_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_13 * hsh_165[k]
                   + f_3 * pc_x[k] * ish_165[k];

        t_215[k] = f_13 * hsh_166[k]
                   + f_3 * pc_x[k] * ish_166[k];

        t_216[k] = f_13 * hsh_167[k]
                   + f_3 * pc_x[k] * ish_167[k];

        t_217[k] = pa_z[k] * hsi0_105[k]
                   - f_10 * pc_z[k] * hsi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, hsh_78, hsh_101, hsh_102, isg0_117, \
                         isg0_118, isg1_117, isg1_118, ish_162, ish_164, \
                         ish_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * hsh_78[k]
                   + f_3 * pc_z[k] * ish_162[k];

        t_219[k] = f_12 * hsh_101[k]
                   + f_8 * isg0_117[k]
                   - f_9 * isg1_117[k]
                   + f_3 * pc_y[k] * ish_164[k];

        t_220[k] = f_12 * hsh_102[k]
                   + f_6 * isg0_118[k]
                   - f_7 * isg1_118[k]
                   + f_3 * pc_y[k] * ish_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, hsi0_140, hsh_83, \
                         hsh_103, hsh_104, hsi1_140, isg0_119, isg1_119, ish_166, \
                         ish_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * hsh_103[k]
                   + f_4 * isg0_119[k]
                   - f_5 * isg1_119[k]
                   + f_3 * pc_y[k] * ish_166[k];

        t_222[k] = f_12 * hsh_104[k]
                   + f_3 * pc_y[k] * ish_167[k];

        t_223[k] = f_11 * hsh_83[k]
                   + f_1 * isg0_119[k]
                   - f_2 * isg1_119[k]
                   + f_3 * pc_z[k] * ish_167[k];

        t_224[k] = pa_y[k] * hsi0_140[k]
                   - f_10 * pc_y[k] * hsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, hsi0_143, hsh_84, \
                         hsh_105, hsh_106, hsh_107, hsi1_143, ish_168, \
                         ish_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * hsh_105[k]
                   + f_3 * pc_y[k] * ish_168[k];

        t_226[k] = f_12 * hsh_84[k]
                   + f_3 * pc_z[k] * ish_168[k];

        t_227[k] = pa_y[k] * hsi0_143[k]
                   + f_12 * hsh_106[k]
                   - f_10 * pc_y[k] * hsi1_143[k];

        t_228[k] = f_11 * hsh_107[k]
                   + f_3 * pc_y[k] * ish_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, hsi0_145, hsi0_146, \
                         hsh_87, hsh_108, hsh_110, hsi1_145, hsi1_146, ish_171, \
                         ish_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * hsi0_145[k]
                   - f_10 * pc_y[k] * hsi1_145[k];

        t_230[k] = pa_y[k] * hsi0_146[k]
                   + f_13 * hsh_108[k]
                   - f_10 * pc_y[k] * hsi1_146[k];

        t_231[k] = f_12 * hsh_87[k]
                   + f_3 * pc_z[k] * ish_171[k];

        t_232[k] = f_11 * hsh_110[k]
                   + f_3 * pc_y[k] * ish_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, hsi0_149, hsi0_150, hsh_90, \
                         hsh_111, hsi1_149, hsi1_150, ish_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * hsi0_149[k]
                   - f_10 * pc_y[k] * hsi1_149[k];

        t_234[k] = pa_y[k] * hsi0_150[k]
                   + f_14 * hsh_111[k]
                   - f_10 * pc_y[k] * hsi1_150[k];

        t_235[k] = f_12 * hsh_90[k]
                   + f_3 * pc_z[k] * ish_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, hsi0_152, hsi0_154, \
                         hsh_113, hsh_114, hsh_183, hsi1_152, hsi1_154, ish_177, \
                         ish_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * hsi0_152[k]
                   + f_12 * hsh_113[k]
                   - f_10 * pc_y[k] * hsi1_152[k];

        t_237[k] = f_11 * hsh_114[k]
                   + f_3 * pc_y[k] * ish_177[k];

        t_238[k] = pa_y[k] * hsi0_154[k]
                   - f_10 * pc_y[k] * hsi1_154[k];

        t_239[k] = f_13 * hsh_183[k]
                   + f_3 * pc_x[k] * ish_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, hsh_184, hsh_185, hsh_186, \
                         hsh_187, hsh_188, ish_184, ish_185, ish_186, ish_187, \
                         ish_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_13 * hsh_184[k]
                   + f_3 * pc_x[k] * ish_184[k];

        t_241[k] = f_13 * hsh_185[k]
                   + f_3 * pc_x[k] * ish_185[k];

        t_242[k] = f_13 * hsh_186[k]
                   + f_3 * pc_x[k] * ish_186[k];

        t_243[k] = f_13 * hsh_187[k]
                   + f_3 * pc_x[k] * ish_187[k];

        t_244[k] = f_13 * hsh_188[k]
                   + f_3 * pc_x[k] * ish_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, hsh_99, hsh_120, hsh_122, isg0_130, \
                         isg0_132, isg1_130, isg1_132, ish_183, \
                         ish_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * hsh_120[k]
                   + f_1 * isg0_130[k]
                   - f_2 * isg1_130[k]
                   + f_3 * pc_y[k] * ish_183[k];

        t_246[k] = f_12 * hsh_99[k]
                   + f_3 * pc_z[k] * ish_183[k];

        t_247[k] = f_11 * hsh_122[k]
                   + f_8 * isg0_132[k]
                   - f_9 * isg1_132[k]
                   + f_3 * pc_y[k] * ish_185[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
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

    const auto *hsi0_167 = buffer.data(hsi0 + 167);
    const auto *hsi0_168 = buffer.data(hsi0 + 168);
    const auto *hsi0_171 = buffer.data(hsi0 + 171);
    const auto *hsi0_174 = buffer.data(hsi0 + 174);
    const auto *hsi0_178 = buffer.data(hsi0 + 178);
    const auto *hsi0_180 = buffer.data(hsi0 + 180);
    const auto *hsi0_189 = buffer.data(hsi0 + 189);

    const auto *hsh_105 = buffer.data(hsh + 105);
    const auto *hsh_123 = buffer.data(hsh + 123);
    const auto *hsh_124 = buffer.data(hsh + 124);
    const auto *hsh_125 = buffer.data(hsh + 125);
    const auto *hsh_126 = buffer.data(hsh + 126);
    const auto *hsh_129 = buffer.data(hsh + 129);
    const auto *hsh_131 = buffer.data(hsh + 131);
    const auto *hsh_132 = buffer.data(hsh + 132);
    const auto *hsh_133 = buffer.data(hsh + 133);
    const auto *hsh_135 = buffer.data(hsh + 135);
    const auto *hsh_141 = buffer.data(hsh + 141);
    const auto *hsh_146 = buffer.data(hsh + 146);
    const auto *hsh_147 = buffer.data(hsh + 147);
    const auto *hsh_149 = buffer.data(hsh + 149);
    const auto *hsh_150 = buffer.data(hsh + 150);
    const auto *hsh_152 = buffer.data(hsh + 152);
    const auto *hsh_153 = buffer.data(hsh + 153);
    const auto *hsh_156 = buffer.data(hsh + 156);
    const auto *hsh_162 = buffer.data(hsh + 162);
    const auto *hsh_164 = buffer.data(hsh + 164);
    const auto *hsh_165 = buffer.data(hsh + 165);
    const auto *hsh_166 = buffer.data(hsh + 166);
    const auto *hsh_167 = buffer.data(hsh + 167);
    const auto *hsh_168 = buffer.data(hsh + 168);
    const auto *hsh_170 = buffer.data(hsh + 170);
    const auto *hsh_173 = buffer.data(hsh + 173);
    const auto *hsh_177 = buffer.data(hsh + 177);
    const auto *hsh_183 = buffer.data(hsh + 183);
    const auto *hsh_185 = buffer.data(hsh + 185);
    const auto *hsh_186 = buffer.data(hsh + 186);
    const auto *hsh_187 = buffer.data(hsh + 187);
    const auto *hsh_189 = buffer.data(hsh + 189);
    const auto *hsh_194 = buffer.data(hsh + 194);
    const auto *hsh_198 = buffer.data(hsh + 198);
    const auto *hsh_203 = buffer.data(hsh + 203);
    const auto *hsh_204 = buffer.data(hsh + 204);
    const auto *hsh_205 = buffer.data(hsh + 205);
    const auto *hsh_206 = buffer.data(hsh + 206);
    const auto *hsh_207 = buffer.data(hsh + 207);
    const auto *hsh_209 = buffer.data(hsh + 209);
    const auto *hsh_210 = buffer.data(hsh + 210);
    const auto *hsh_213 = buffer.data(hsh + 213);
    const auto *hsh_216 = buffer.data(hsh + 216);
    const auto *hsh_220 = buffer.data(hsh + 220);
    const auto *hsh_225 = buffer.data(hsh + 225);
    const auto *hsh_227 = buffer.data(hsh + 227);
    const auto *hsh_228 = buffer.data(hsh + 228);
    const auto *hsh_229 = buffer.data(hsh + 229);
    const auto *hsh_230 = buffer.data(hsh + 230);
    const auto *hsh_236 = buffer.data(hsh + 236);
    const auto *hsh_240 = buffer.data(hsh + 240);
    const auto *hsh_245 = buffer.data(hsh + 245);
    const auto *hsh_246 = buffer.data(hsh + 246);
    const auto *hsh_247 = buffer.data(hsh + 247);
    const auto *hsh_248 = buffer.data(hsh + 248);
    const auto *hsh_249 = buffer.data(hsh + 249);
    const auto *hsh_250 = buffer.data(hsh + 250);
    const auto *hsh_251 = buffer.data(hsh + 251);
    const auto *hsh_252 = buffer.data(hsh + 252);
    const auto *hsh_255 = buffer.data(hsh + 255);
    const auto *hsh_257 = buffer.data(hsh + 257);
    const auto *hsh_258 = buffer.data(hsh + 258);
    const auto *hsh_261 = buffer.data(hsh + 261);
    const auto *hsh_262 = buffer.data(hsh + 262);
    const auto *hsh_264 = buffer.data(hsh + 264);
    const auto *hsh_266 = buffer.data(hsh + 266);
    const auto *hsh_267 = buffer.data(hsh + 267);
    const auto *hsh_268 = buffer.data(hsh + 268);
    const auto *hsh_269 = buffer.data(hsh + 269);
    const auto *hsh_270 = buffer.data(hsh + 270);
    const auto *hsh_271 = buffer.data(hsh + 271);
    const auto *hsh_272 = buffer.data(hsh + 272);

    const auto *hsi1_167 = buffer.data(hsi1 + 167);
    const auto *hsi1_168 = buffer.data(hsi1 + 168);
    const auto *hsi1_171 = buffer.data(hsi1 + 171);
    const auto *hsi1_174 = buffer.data(hsi1 + 174);
    const auto *hsi1_178 = buffer.data(hsi1 + 178);
    const auto *hsi1_180 = buffer.data(hsi1 + 180);
    const auto *hsi1_189 = buffer.data(hsi1 + 189);

    const auto *isg0_133 = buffer.data(isg0 + 133);
    const auto *isg0_134 = buffer.data(isg0 + 134);
    const auto *isg0_135 = buffer.data(isg0 + 135);
    const auto *isg0_136 = buffer.data(isg0 + 136);
    const auto *isg0_137 = buffer.data(isg0 + 137);
    const auto *isg0_138 = buffer.data(isg0 + 138);
    const auto *isg0_139 = buffer.data(isg0 + 139);
    const auto *isg0_140 = buffer.data(isg0 + 140);
    const auto *isg0_144 = buffer.data(isg0 + 144);
    const auto *isg0_145 = buffer.data(isg0 + 145);
    const auto *isg0_146 = buffer.data(isg0 + 146);
    const auto *isg0_147 = buffer.data(isg0 + 147);
    const auto *isg0_148 = buffer.data(isg0 + 148);
    const auto *isg0_149 = buffer.data(isg0 + 149);
    const auto *isg0_150 = buffer.data(isg0 + 150);
    const auto *isg0_152 = buffer.data(isg0 + 152);
    const auto *isg0_153 = buffer.data(isg0 + 153);
    const auto *isg0_155 = buffer.data(isg0 + 155);
    const auto *isg0_156 = buffer.data(isg0 + 156);
    const auto *isg0_160 = buffer.data(isg0 + 160);
    const auto *isg0_161 = buffer.data(isg0 + 161);
    const auto *isg0_162 = buffer.data(isg0 + 162);
    const auto *isg0_164 = buffer.data(isg0 + 164);
    const auto *isg0_170 = buffer.data(isg0 + 170);
    const auto *isg0_174 = buffer.data(isg0 + 174);
    const auto *isg0_177 = buffer.data(isg0 + 177);
    const auto *isg0_178 = buffer.data(isg0 + 178);
    const auto *isg0_179 = buffer.data(isg0 + 179);
    const auto *isg0_180 = buffer.data(isg0 + 180);
    const auto *isg0_183 = buffer.data(isg0 + 183);
    const auto *isg0_185 = buffer.data(isg0 + 185);
    const auto *isg0_186 = buffer.data(isg0 + 186);
    const auto *isg0_189 = buffer.data(isg0 + 189);
    const auto *isg0_190 = buffer.data(isg0 + 190);
    const auto *isg0_192 = buffer.data(isg0 + 192);
    const auto *isg0_193 = buffer.data(isg0 + 193);
    const auto *isg0_194 = buffer.data(isg0 + 194);

    const auto *isg1_133 = buffer.data(isg1 + 133);
    const auto *isg1_134 = buffer.data(isg1 + 134);
    const auto *isg1_135 = buffer.data(isg1 + 135);
    const auto *isg1_136 = buffer.data(isg1 + 136);
    const auto *isg1_137 = buffer.data(isg1 + 137);
    const auto *isg1_138 = buffer.data(isg1 + 138);
    const auto *isg1_139 = buffer.data(isg1 + 139);
    const auto *isg1_140 = buffer.data(isg1 + 140);
    const auto *isg1_144 = buffer.data(isg1 + 144);
    const auto *isg1_145 = buffer.data(isg1 + 145);
    const auto *isg1_146 = buffer.data(isg1 + 146);
    const auto *isg1_147 = buffer.data(isg1 + 147);
    const auto *isg1_148 = buffer.data(isg1 + 148);
    const auto *isg1_149 = buffer.data(isg1 + 149);
    const auto *isg1_150 = buffer.data(isg1 + 150);
    const auto *isg1_152 = buffer.data(isg1 + 152);
    const auto *isg1_153 = buffer.data(isg1 + 153);
    const auto *isg1_155 = buffer.data(isg1 + 155);
    const auto *isg1_156 = buffer.data(isg1 + 156);
    const auto *isg1_160 = buffer.data(isg1 + 160);
    const auto *isg1_161 = buffer.data(isg1 + 161);
    const auto *isg1_162 = buffer.data(isg1 + 162);
    const auto *isg1_164 = buffer.data(isg1 + 164);
    const auto *isg1_170 = buffer.data(isg1 + 170);
    const auto *isg1_174 = buffer.data(isg1 + 174);
    const auto *isg1_177 = buffer.data(isg1 + 177);
    const auto *isg1_178 = buffer.data(isg1 + 178);
    const auto *isg1_179 = buffer.data(isg1 + 179);
    const auto *isg1_180 = buffer.data(isg1 + 180);
    const auto *isg1_183 = buffer.data(isg1 + 183);
    const auto *isg1_185 = buffer.data(isg1 + 185);
    const auto *isg1_186 = buffer.data(isg1 + 186);
    const auto *isg1_189 = buffer.data(isg1 + 189);
    const auto *isg1_190 = buffer.data(isg1 + 190);
    const auto *isg1_192 = buffer.data(isg1 + 192);
    const auto *isg1_193 = buffer.data(isg1 + 193);
    const auto *isg1_194 = buffer.data(isg1 + 194);

    const auto *ish_186 = buffer.data(ish + 186);
    const auto *ish_187 = buffer.data(ish + 187);
    const auto *ish_188 = buffer.data(ish + 188);
    const auto *ish_189 = buffer.data(ish + 189);
    const auto *ish_190 = buffer.data(ish + 190);
    const auto *ish_191 = buffer.data(ish + 191);
    const auto *ish_192 = buffer.data(ish + 192);
    const auto *ish_193 = buffer.data(ish + 193);
    const auto *ish_194 = buffer.data(ish + 194);
    const auto *ish_195 = buffer.data(ish + 195);
    const auto *ish_196 = buffer.data(ish + 196);
    const auto *ish_197 = buffer.data(ish + 197);
    const auto *ish_198 = buffer.data(ish + 198);
    const auto *ish_203 = buffer.data(ish + 203);
    const auto *ish_204 = buffer.data(ish + 204);
    const auto *ish_205 = buffer.data(ish + 205);
    const auto *ish_206 = buffer.data(ish + 206);
    const auto *ish_207 = buffer.data(ish + 207);
    const auto *ish_208 = buffer.data(ish + 208);
    const auto *ish_209 = buffer.data(ish + 209);
    const auto *ish_210 = buffer.data(ish + 210);
    const auto *ish_211 = buffer.data(ish + 211);
    const auto *ish_212 = buffer.data(ish + 212);
    const auto *ish_213 = buffer.data(ish + 213);
    const auto *ish_215 = buffer.data(ish + 215);
    const auto *ish_216 = buffer.data(ish + 216);
    const auto *ish_217 = buffer.data(ish + 217);
    const auto *ish_219 = buffer.data(ish + 219);
    const auto *ish_220 = buffer.data(ish + 220);
    const auto *ish_225 = buffer.data(ish + 225);
    const auto *ish_226 = buffer.data(ish + 226);
    const auto *ish_227 = buffer.data(ish + 227);
    const auto *ish_228 = buffer.data(ish + 228);
    const auto *ish_229 = buffer.data(ish + 229);
    const auto *ish_230 = buffer.data(ish + 230);
    const auto *ish_231 = buffer.data(ish + 231);
    const auto *ish_233 = buffer.data(ish + 233);
    const auto *ish_234 = buffer.data(ish + 234);
    const auto *ish_236 = buffer.data(ish + 236);
    const auto *ish_237 = buffer.data(ish + 237);
    const auto *ish_240 = buffer.data(ish + 240);
    const auto *ish_245 = buffer.data(ish + 245);
    const auto *ish_246 = buffer.data(ish + 246);
    const auto *ish_247 = buffer.data(ish + 247);
    const auto *ish_248 = buffer.data(ish + 248);
    const auto *ish_249 = buffer.data(ish + 249);
    const auto *ish_250 = buffer.data(ish + 250);
    const auto *ish_251 = buffer.data(ish + 251);
    const auto *ish_252 = buffer.data(ish + 252);
    const auto *ish_254 = buffer.data(ish + 254);
    const auto *ish_255 = buffer.data(ish + 255);
    const auto *ish_257 = buffer.data(ish + 257);
    const auto *ish_258 = buffer.data(ish + 258);
    const auto *ish_261 = buffer.data(ish + 261);
    const auto *ish_262 = buffer.data(ish + 262);
    const auto *ish_264 = buffer.data(ish + 264);
    const auto *ish_266 = buffer.data(ish + 266);
    const auto *ish_267 = buffer.data(ish + 267);
    const auto *ish_268 = buffer.data(ish + 268);
    const auto *ish_269 = buffer.data(ish + 269);
    const auto *ish_270 = buffer.data(ish + 270);
    const auto *ish_271 = buffer.data(ish + 271);
    const auto *ish_272 = buffer.data(ish + 272);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, hsh_123, hsh_124, hsh_125, isg0_133, \
                         isg0_134, isg1_133, isg1_134, ish_186, ish_187, \
                         ish_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * hsh_123[k]
                   + f_6 * isg0_133[k]
                   - f_7 * isg1_133[k]
                   + f_3 * pc_y[k] * ish_186[k];

        t_249[k] = f_11 * hsh_124[k]
                   + f_4 * isg0_134[k]
                   - f_5 * isg1_134[k]
                   + f_3 * pc_y[k] * ish_187[k];

        t_250[k] = f_11 * hsh_125[k]
                   + f_3 * pc_y[k] * ish_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, hsi0_167, \
                         hsh_105, hsh_189, hsi1_167, isg0_135, isg1_135, \
                         ish_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * hsi0_167[k]
                   - f_10 * pc_y[k] * hsi1_167[k];

        t_252[k] = f_13 * hsh_189[k]
                   + f_1 * isg0_135[k]
                   - f_2 * isg1_135[k]
                   + f_3 * pc_x[k] * ish_189[k];

        t_253[k] = f_3 * pc_y[k] * ish_189[k];

        t_254[k] = f_13 * hsh_105[k]
                   + f_3 * pc_z[k] * ish_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, hsh_194, isg0_135, isg0_140, \
                         isg1_135, isg1_140, ish_190, ish_191, \
                         ish_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * isg0_135[k]
                   - f_5 * isg1_135[k]
                   + f_3 * pc_y[k] * ish_190[k];

        t_256[k] = f_3 * pc_y[k] * ish_191[k];

        t_257[k] = f_13 * hsh_194[k]
                   + f_8 * isg0_140[k]
                   - f_9 * isg1_140[k]
                   + f_3 * pc_x[k] * ish_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, isg0_136, isg0_137, isg1_136, isg1_137, \
                         ish_192, ish_193, ish_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * isg0_136[k]
                   - f_7 * isg1_136[k]
                   + f_3 * pc_y[k] * ish_192[k];

        t_259[k] = f_4 * isg0_137[k]
                   - f_5 * isg1_137[k]
                   + f_3 * pc_y[k] * ish_193[k];

        t_260[k] = f_3 * pc_y[k] * ish_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, hsh_198, isg0_138, isg0_139, \
                         isg0_144, isg1_138, isg1_139, isg1_144, ish_195, ish_196, \
                         ish_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_13 * hsh_198[k]
                   + f_6 * isg0_144[k]
                   - f_7 * isg1_144[k]
                   + f_3 * pc_x[k] * ish_198[k];

        t_262[k] = f_8 * isg0_138[k]
                   - f_9 * isg1_138[k]
                   + f_3 * pc_y[k] * ish_195[k];

        t_263[k] = f_6 * isg0_139[k]
                   - f_7 * isg1_139[k]
                   + f_3 * pc_y[k] * ish_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, hsh_203, hsh_204, isg0_140, \
                         isg0_149, isg1_140, isg1_149, ish_197, ish_198, ish_203, \
                         ish_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * isg0_140[k]
                   - f_5 * isg1_140[k]
                   + f_3 * pc_y[k] * ish_197[k];

        t_265[k] = f_3 * pc_y[k] * ish_198[k];

        t_266[k] = f_13 * hsh_203[k]
                   + f_4 * isg0_149[k]
                   - f_5 * isg1_149[k]
                   + f_3 * pc_x[k] * ish_203[k];

        t_267[k] = f_13 * hsh_204[k]
                   + f_3 * pc_x[k] * ish_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, hsh_205, hsh_206, \
                         hsh_207, hsh_209, ish_203, ish_205, ish_206, ish_207, \
                         ish_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_13 * hsh_205[k]
                   + f_3 * pc_x[k] * ish_205[k];

        t_269[k] = f_13 * hsh_206[k]
                   + f_3 * pc_x[k] * ish_206[k];

        t_270[k] = f_13 * hsh_207[k]
                   + f_3 * pc_x[k] * ish_207[k];

        t_271[k] = f_3 * pc_y[k] * ish_203[k];

        t_272[k] = f_13 * hsh_209[k]
                   + f_3 * pc_x[k] * ish_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, isg0_145, isg0_146, isg0_147, isg1_145, \
                         isg1_146, isg1_147, ish_204, ish_205, \
                         ish_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * isg0_145[k]
                   - f_2 * isg1_145[k]
                   + f_3 * pc_y[k] * ish_204[k];

        t_274[k] = f_16 * isg0_146[k]
                   - f_17 * isg1_146[k]
                   + f_3 * pc_y[k] * ish_205[k];

        t_275[k] = f_8 * isg0_147[k]
                   - f_9 * isg1_147[k]
                   + f_3 * pc_y[k] * ish_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, hsh_125, isg0_148, isg0_149, \
                         isg1_148, isg1_149, ish_207, ish_208, \
                         ish_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * isg0_148[k]
                   - f_7 * isg1_148[k]
                   + f_3 * pc_y[k] * ish_207[k];

        t_277[k] = f_4 * isg0_149[k]
                   - f_5 * isg1_149[k]
                   + f_3 * pc_y[k] * ish_208[k];

        t_278[k] = f_3 * pc_y[k] * ish_209[k];

        t_279[k] = f_13 * hsh_125[k]
                   + f_1 * isg0_149[k]
                   - f_2 * isg1_149[k]
                   + f_3 * pc_z[k] * ish_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, hsh_126, hsh_210, \
                         hsh_213, isg0_150, isg0_153, isg1_150, isg1_153, ish_210, \
                         ish_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_12 * hsh_210[k]
                   + f_1 * isg0_150[k]
                   - f_2 * isg1_150[k]
                   + f_3 * pc_x[k] * ish_210[k];

        t_281[k] = f_14 * hsh_126[k]
                   + f_3 * pc_y[k] * ish_210[k];

        t_282[k] = f_3 * pc_z[k] * ish_210[k];

        t_283[k] = f_12 * hsh_213[k]
                   + f_8 * isg0_153[k]
                   - f_9 * isg1_153[k]
                   + f_3 * pc_x[k] * ish_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, hsh_216, isg0_150, isg0_156, \
                         isg1_150, isg1_156, ish_211, ish_212, ish_213, \
                         ish_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * ish_211[k];

        t_285[k] = f_4 * isg0_150[k]
                   - f_5 * isg1_150[k]
                   + f_3 * pc_z[k] * ish_212[k];

        t_286[k] = f_12 * hsh_216[k]
                   + f_6 * isg0_156[k]
                   - f_7 * isg1_156[k]
                   + f_3 * pc_x[k] * ish_216[k];

        t_287[k] = f_3 * pc_z[k] * ish_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, hsh_131, hsh_220, \
                         isg0_152, isg0_160, isg1_152, isg1_160, ish_215, ish_216, \
                         ish_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * hsh_131[k]
                   + f_3 * pc_y[k] * ish_215[k];

        t_289[k] = f_6 * isg0_152[k]
                   - f_7 * isg1_152[k]
                   + f_3 * pc_z[k] * ish_215[k];

        t_290[k] = f_12 * hsh_220[k]
                   + f_4 * isg0_160[k]
                   - f_5 * isg1_160[k]
                   + f_3 * pc_x[k] * ish_220[k];

        t_291[k] = f_3 * pc_z[k] * ish_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, hsh_135, hsh_225, \
                         isg0_153, isg0_155, isg1_153, isg1_155, ish_217, ish_219, \
                         ish_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * isg0_153[k]
                   - f_5 * isg1_153[k]
                   + f_3 * pc_z[k] * ish_217[k];

        t_293[k] = f_14 * hsh_135[k]
                   + f_3 * pc_y[k] * ish_219[k];

        t_294[k] = f_8 * isg0_155[k]
                   - f_9 * isg1_155[k]
                   + f_3 * pc_z[k] * ish_219[k];

        t_295[k] = f_12 * hsh_225[k]
                   + f_3 * pc_x[k] * ish_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, hsh_227, hsh_228, \
                         hsh_229, hsh_230, ish_220, ish_227, ish_228, ish_229, \
                         ish_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * ish_220[k];

        t_297[k] = f_12 * hsh_227[k]
                   + f_3 * pc_x[k] * ish_227[k];

        t_298[k] = f_12 * hsh_228[k]
                   + f_3 * pc_x[k] * ish_228[k];

        t_299[k] = f_12 * hsh_229[k]
                   + f_3 * pc_x[k] * ish_229[k];

        t_300[k] = f_12 * hsh_230[k]
                   + f_3 * pc_x[k] * ish_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, hsh_141, isg0_160, isg0_161, \
                         isg1_160, isg1_161, ish_225, ish_226, \
                         ish_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * hsh_141[k]
                   + f_1 * isg0_160[k]
                   - f_2 * isg1_160[k]
                   + f_3 * pc_y[k] * ish_225[k];

        t_302[k] = f_3 * pc_z[k] * ish_225[k];

        t_303[k] = f_4 * isg0_160[k]
                   - f_5 * isg1_160[k]
                   + f_3 * pc_z[k] * ish_226[k];

        t_304[k] = f_6 * isg0_161[k]
                   - f_7 * isg1_161[k]
                   + f_3 * pc_z[k] * ish_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, hsi0_168, hsh_146, \
                         hsi1_168, isg0_162, isg0_164, isg1_162, isg1_164, ish_228, \
                         ish_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * isg0_162[k]
                   - f_9 * isg1_162[k]
                   + f_3 * pc_z[k] * ish_228[k];

        t_306[k] = f_14 * hsh_146[k]
                   + f_3 * pc_y[k] * ish_230[k];

        t_307[k] = f_1 * isg0_164[k]
                   - f_2 * isg1_164[k]
                   + f_3 * pc_z[k] * ish_230[k];

        t_308[k] = pa_z[k] * hsi0_168[k]
                   - f_10 * pc_z[k] * hsi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, hsi0_171, hsh_126, \
                         hsh_147, hsh_149, hsi1_171, ish_231, ish_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * hsh_147[k]
                   + f_3 * pc_y[k] * ish_231[k];

        t_310[k] = f_11 * hsh_126[k]
                   + f_3 * pc_z[k] * ish_231[k];

        t_311[k] = pa_z[k] * hsi0_171[k]
                   - f_10 * pc_z[k] * hsi1_171[k];

        t_312[k] = f_13 * hsh_149[k]
                   + f_3 * pc_y[k] * ish_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, hsi0_174, hsh_129, hsh_236, \
                         hsi1_174, isg0_170, isg1_170, ish_234, \
                         ish_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_12 * hsh_236[k]
                   + f_8 * isg0_170[k]
                   - f_9 * isg1_170[k]
                   + f_3 * pc_x[k] * ish_236[k];

        t_314[k] = pa_z[k] * hsi0_174[k]
                   - f_10 * pc_z[k] * hsi1_174[k];

        t_315[k] = f_11 * hsh_129[k]
                   + f_3 * pc_z[k] * ish_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, hsi0_178, hsh_152, \
                         hsh_240, hsi1_178, isg0_174, isg1_174, ish_236, \
                         ish_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * hsh_152[k]
                   + f_3 * pc_y[k] * ish_236[k];

        t_317[k] = f_12 * hsh_240[k]
                   + f_6 * isg0_174[k]
                   - f_7 * isg1_174[k]
                   + f_3 * pc_x[k] * ish_240[k];

        t_318[k] = pa_z[k] * hsi0_178[k]
                   - f_10 * pc_z[k] * hsi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, hsi0_180, hsh_132, hsh_133, \
                         hsh_156, hsi1_180, ish_237, ish_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * hsh_132[k]
                   + f_3 * pc_z[k] * ish_237[k];

        t_320[k] = pa_z[k] * hsi0_180[k]
                   + f_12 * hsh_133[k]
                   - f_10 * pc_z[k] * hsi1_180[k];

        t_321[k] = f_13 * hsh_156[k]
                   + f_3 * pc_y[k] * ish_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, hsh_245, hsh_246, hsh_247, hsh_248, \
                         isg0_179, isg1_179, ish_245, ish_246, ish_247, \
                         ish_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_12 * hsh_245[k]
                   + f_4 * isg0_179[k]
                   - f_5 * isg1_179[k]
                   + f_3 * pc_x[k] * ish_245[k];

        t_323[k] = f_12 * hsh_246[k]
                   + f_3 * pc_x[k] * ish_246[k];

        t_324[k] = f_12 * hsh_247[k]
                   + f_3 * pc_x[k] * ish_247[k];

        t_325[k] = f_12 * hsh_248[k]
                   + f_3 * pc_x[k] * ish_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, hsi0_189, hsh_249, \
                         hsh_250, hsh_251, hsi1_189, ish_249, ish_250, \
                         ish_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * hsh_249[k]
                   + f_3 * pc_x[k] * ish_249[k];

        t_327[k] = f_12 * hsh_250[k]
                   + f_3 * pc_x[k] * ish_250[k];

        t_328[k] = f_12 * hsh_251[k]
                   + f_3 * pc_x[k] * ish_251[k];

        t_329[k] = pa_z[k] * hsi0_189[k]
                   - f_10 * pc_z[k] * hsi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, hsh_141, hsh_164, hsh_165, isg0_177, \
                         isg0_178, isg1_177, isg1_178, ish_246, ish_248, \
                         ish_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * hsh_141[k]
                   + f_3 * pc_z[k] * ish_246[k];

        t_331[k] = f_13 * hsh_164[k]
                   + f_8 * isg0_177[k]
                   - f_9 * isg1_177[k]
                   + f_3 * pc_y[k] * ish_248[k];

        t_332[k] = f_13 * hsh_165[k]
                   + f_6 * isg0_178[k]
                   - f_7 * isg1_178[k]
                   + f_3 * pc_y[k] * ish_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, hsh_146, hsh_166, hsh_167, isg0_179, \
                         isg1_179, ish_250, ish_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * hsh_166[k]
                   + f_4 * isg0_179[k]
                   - f_5 * isg1_179[k]
                   + f_3 * pc_y[k] * ish_250[k];

        t_334[k] = f_13 * hsh_167[k]
                   + f_3 * pc_y[k] * ish_251[k];

        t_335[k] = f_11 * hsh_146[k]
                   + f_1 * isg0_179[k]
                   - f_2 * isg1_179[k]
                   + f_3 * pc_z[k] * ish_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, hsh_147, hsh_168, hsh_252, \
                         isg0_180, isg1_180, ish_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_12 * hsh_252[k]
                   + f_1 * isg0_180[k]
                   - f_2 * isg1_180[k]
                   + f_3 * pc_x[k] * ish_252[k];

        t_337[k] = f_12 * hsh_168[k]
                   + f_3 * pc_y[k] * ish_252[k];

        t_338[k] = f_12 * hsh_147[k]
                   + f_3 * pc_z[k] * ish_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, hsh_170, hsh_255, hsh_257, isg0_183, \
                         isg0_185, isg1_183, isg1_185, ish_254, ish_255, \
                         ish_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_12 * hsh_255[k]
                   + f_8 * isg0_183[k]
                   - f_9 * isg1_183[k]
                   + f_3 * pc_x[k] * ish_255[k];

        t_340[k] = f_12 * hsh_170[k]
                   + f_3 * pc_y[k] * ish_254[k];

        t_341[k] = f_12 * hsh_257[k]
                   + f_8 * isg0_185[k]
                   - f_9 * isg1_185[k]
                   + f_3 * pc_x[k] * ish_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, hsh_150, hsh_173, hsh_258, \
                         isg0_186, isg1_186, ish_255, ish_257, \
                         ish_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_12 * hsh_258[k]
                   + f_6 * isg0_186[k]
                   - f_7 * isg1_186[k]
                   + f_3 * pc_x[k] * ish_258[k];

        t_343[k] = f_12 * hsh_150[k]
                   + f_3 * pc_z[k] * ish_255[k];

        t_344[k] = f_12 * hsh_173[k]
                   + f_3 * pc_y[k] * ish_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, hsh_153, hsh_261, hsh_262, isg0_189, \
                         isg0_190, isg1_189, isg1_190, ish_258, ish_261, \
                         ish_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_12 * hsh_261[k]
                   + f_6 * isg0_189[k]
                   - f_7 * isg1_189[k]
                   + f_3 * pc_x[k] * ish_261[k];

        t_346[k] = f_12 * hsh_262[k]
                   + f_4 * isg0_190[k]
                   - f_5 * isg1_190[k]
                   + f_3 * pc_x[k] * ish_262[k];

        t_347[k] = f_12 * hsh_153[k]
                   + f_3 * pc_z[k] * ish_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, hsh_177, hsh_264, hsh_266, isg0_192, \
                         isg0_194, isg1_192, isg1_194, ish_261, ish_264, \
                         ish_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_12 * hsh_264[k]
                   + f_4 * isg0_192[k]
                   - f_5 * isg1_192[k]
                   + f_3 * pc_x[k] * ish_264[k];

        t_349[k] = f_12 * hsh_177[k]
                   + f_3 * pc_y[k] * ish_261[k];

        t_350[k] = f_12 * hsh_266[k]
                   + f_4 * isg0_194[k]
                   - f_5 * isg1_194[k]
                   + f_3 * pc_x[k] * ish_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, hsh_267, hsh_268, hsh_269, \
                         hsh_270, hsh_271, ish_267, ish_268, ish_269, ish_270, \
                         ish_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_12 * hsh_267[k]
                   + f_3 * pc_x[k] * ish_267[k];

        t_352[k] = f_12 * hsh_268[k]
                   + f_3 * pc_x[k] * ish_268[k];

        t_353[k] = f_12 * hsh_269[k]
                   + f_3 * pc_x[k] * ish_269[k];

        t_354[k] = f_12 * hsh_270[k]
                   + f_3 * pc_x[k] * ish_270[k];

        t_355[k] = f_12 * hsh_271[k]
                   + f_3 * pc_x[k] * ish_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, hsh_162, hsh_183, hsh_272, \
                         isg0_190, isg1_190, ish_267, ish_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_12 * hsh_272[k]
                   + f_3 * pc_x[k] * ish_272[k];

        t_357[k] = f_12 * hsh_183[k]
                   + f_1 * isg0_190[k]
                   - f_2 * isg1_190[k]
                   + f_3 * pc_y[k] * ish_267[k];

        t_358[k] = f_12 * hsh_162[k]
                   + f_3 * pc_z[k] * ish_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, hsh_185, hsh_186, hsh_187, isg0_192, \
                         isg0_193, isg0_194, isg1_192, isg1_193, isg1_194, ish_269, ish_270, \
                         ish_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * hsh_185[k]
                   + f_8 * isg0_192[k]
                   - f_9 * isg1_192[k]
                   + f_3 * pc_y[k] * ish_269[k];

        t_360[k] = f_12 * hsh_186[k]
                   + f_6 * isg0_193[k]
                   - f_7 * isg1_193[k]
                   + f_3 * pc_y[k] * ish_270[k];

        t_361[k] = f_12 * hsh_187[k]
                   + f_4 * isg0_194[k]
                   - f_5 * isg1_194[k]
                   + f_3 * pc_y[k] * ish_271[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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
    auto *t_479 = buffer.data(target + 479);
    auto *t_480 = buffer.data(target + 480);
    auto *t_481 = buffer.data(target + 481);
    auto *t_482 = buffer.data(target + 482);
    auto *t_483 = buffer.data(target + 483);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsi0_252 = buffer.data(hsi0 + 252);
    const auto *hsi0_255 = buffer.data(hsi0 + 255);
    const auto *hsi0_257 = buffer.data(hsi0 + 257);
    const auto *hsi0_258 = buffer.data(hsi0 + 258);
    const auto *hsi0_261 = buffer.data(hsi0 + 261);
    const auto *hsi0_262 = buffer.data(hsi0 + 262);
    const auto *hsi0_264 = buffer.data(hsi0 + 264);
    const auto *hsi0_266 = buffer.data(hsi0 + 266);
    const auto *hsi0_279 = buffer.data(hsi0 + 279);
    const auto *hsi0_280 = buffer.data(hsi0 + 280);
    const auto *hsi0_283 = buffer.data(hsi0 + 283);
    const auto *hsi0_286 = buffer.data(hsi0 + 286);
    const auto *hsi0_290 = buffer.data(hsi0 + 290);
    const auto *hsi0_420 = buffer.data(hsi0 + 420);
    const auto *hsi0_423 = buffer.data(hsi0 + 423);
    const auto *hsi0_426 = buffer.data(hsi0 + 426);
    const auto *hsi0_430 = buffer.data(hsi0 + 430);
    const auto *hsi0_441 = buffer.data(hsi0 + 441);
    const auto *hsi0_443 = buffer.data(hsi0 + 443);
    const auto *hsi0_444 = buffer.data(hsi0 + 444);
    const auto *hsi0_445 = buffer.data(hsi0 + 445);
    const auto *hsi0_447 = buffer.data(hsi0 + 447);
    const auto *hsi0_453 = buffer.data(hsi0 + 453);
    const auto *hsi0_457 = buffer.data(hsi0 + 457);
    const auto *hsi0_460 = buffer.data(hsi0 + 460);
    const auto *hsi0_462 = buffer.data(hsi0 + 462);
    const auto *hsi0_469 = buffer.data(hsi0 + 469);
    const auto *hsi0_471 = buffer.data(hsi0 + 471);
    const auto *hsi0_472 = buffer.data(hsi0 + 472);
    const auto *hsi0_473 = buffer.data(hsi0 + 473);
    const auto *hsi0_475 = buffer.data(hsi0 + 475);
    const auto *hsi0_476 = buffer.data(hsi0 + 476);
    const auto *hsi0_479 = buffer.data(hsi0 + 479);
    const auto *hsi0_481 = buffer.data(hsi0 + 481);
    const auto *hsi0_482 = buffer.data(hsi0 + 482);

    const auto *hsh_167 = buffer.data(hsh + 167);
    const auto *hsh_168 = buffer.data(hsh + 168);
    const auto *hsh_171 = buffer.data(hsh + 171);
    const auto *hsh_174 = buffer.data(hsh + 174);
    const auto *hsh_183 = buffer.data(hsh + 183);
    const auto *hsh_188 = buffer.data(hsh + 188);
    const auto *hsh_189 = buffer.data(hsh + 189);
    const auto *hsh_190 = buffer.data(hsh + 190);
    const auto *hsh_191 = buffer.data(hsh + 191);
    const auto *hsh_192 = buffer.data(hsh + 192);
    const auto *hsh_194 = buffer.data(hsh + 194);
    const auto *hsh_195 = buffer.data(hsh + 195);
    const auto *hsh_197 = buffer.data(hsh + 197);
    const auto *hsh_198 = buffer.data(hsh + 198);
    const auto *hsh_204 = buffer.data(hsh + 204);
    const auto *hsh_206 = buffer.data(hsh + 206);
    const auto *hsh_207 = buffer.data(hsh + 207);
    const auto *hsh_208 = buffer.data(hsh + 208);
    const auto *hsh_209 = buffer.data(hsh + 209);
    const auto *hsh_210 = buffer.data(hsh + 210);
    const auto *hsh_213 = buffer.data(hsh + 213);
    const auto *hsh_215 = buffer.data(hsh + 215);
    const auto *hsh_216 = buffer.data(hsh + 216);
    const auto *hsh_219 = buffer.data(hsh + 219);
    const auto *hsh_225 = buffer.data(hsh + 225);
    const auto *hsh_230 = buffer.data(hsh + 230);
    const auto *hsh_231 = buffer.data(hsh + 231);
    const auto *hsh_233 = buffer.data(hsh + 233);
    const auto *hsh_234 = buffer.data(hsh + 234);
    const auto *hsh_236 = buffer.data(hsh + 236);
    const auto *hsh_240 = buffer.data(hsh + 240);
    const auto *hsh_251 = buffer.data(hsh + 251);
    const auto *hsh_252 = buffer.data(hsh + 252);
    const auto *hsh_254 = buffer.data(hsh + 254);
    const auto *hsh_288 = buffer.data(hsh + 288);
    const auto *hsh_289 = buffer.data(hsh + 289);
    const auto *hsh_290 = buffer.data(hsh + 290);
    const auto *hsh_291 = buffer.data(hsh + 291);
    const auto *hsh_292 = buffer.data(hsh + 292);
    const auto *hsh_293 = buffer.data(hsh + 293);
    const auto *hsh_294 = buffer.data(hsh + 294);
    const auto *hsh_299 = buffer.data(hsh + 299);
    const auto *hsh_303 = buffer.data(hsh + 303);
    const auto *hsh_308 = buffer.data(hsh + 308);
    const auto *hsh_309 = buffer.data(hsh + 309);
    const auto *hsh_310 = buffer.data(hsh + 310);
    const auto *hsh_311 = buffer.data(hsh + 311);
    const auto *hsh_312 = buffer.data(hsh + 312);
    const auto *hsh_314 = buffer.data(hsh + 314);
    const auto *hsh_315 = buffer.data(hsh + 315);
    const auto *hsh_318 = buffer.data(hsh + 318);
    const auto *hsh_321 = buffer.data(hsh + 321);
    const auto *hsh_325 = buffer.data(hsh + 325);
    const auto *hsh_330 = buffer.data(hsh + 330);
    const auto *hsh_332 = buffer.data(hsh + 332);
    const auto *hsh_333 = buffer.data(hsh + 333);
    const auto *hsh_334 = buffer.data(hsh + 334);
    const auto *hsh_335 = buffer.data(hsh + 335);
    const auto *hsh_341 = buffer.data(hsh + 341);
    const auto *hsh_345 = buffer.data(hsh + 345);
    const auto *hsh_348 = buffer.data(hsh + 348);
    const auto *hsh_350 = buffer.data(hsh + 350);
    const auto *hsh_351 = buffer.data(hsh + 351);
    const auto *hsh_352 = buffer.data(hsh + 352);
    const auto *hsh_353 = buffer.data(hsh + 353);
    const auto *hsh_354 = buffer.data(hsh + 354);
    const auto *hsh_355 = buffer.data(hsh + 355);
    const auto *hsh_356 = buffer.data(hsh + 356);
    const auto *hsh_357 = buffer.data(hsh + 357);
    const auto *hsh_360 = buffer.data(hsh + 360);
    const auto *hsh_362 = buffer.data(hsh + 362);
    const auto *hsh_363 = buffer.data(hsh + 363);

    const auto *hsi1_252 = buffer.data(hsi1 + 252);
    const auto *hsi1_255 = buffer.data(hsi1 + 255);
    const auto *hsi1_257 = buffer.data(hsi1 + 257);
    const auto *hsi1_258 = buffer.data(hsi1 + 258);
    const auto *hsi1_261 = buffer.data(hsi1 + 261);
    const auto *hsi1_262 = buffer.data(hsi1 + 262);
    const auto *hsi1_264 = buffer.data(hsi1 + 264);
    const auto *hsi1_266 = buffer.data(hsi1 + 266);
    const auto *hsi1_279 = buffer.data(hsi1 + 279);
    const auto *hsi1_280 = buffer.data(hsi1 + 280);
    const auto *hsi1_283 = buffer.data(hsi1 + 283);
    const auto *hsi1_286 = buffer.data(hsi1 + 286);
    const auto *hsi1_290 = buffer.data(hsi1 + 290);
    const auto *hsi1_420 = buffer.data(hsi1 + 420);
    const auto *hsi1_423 = buffer.data(hsi1 + 423);
    const auto *hsi1_426 = buffer.data(hsi1 + 426);
    const auto *hsi1_430 = buffer.data(hsi1 + 430);
    const auto *hsi1_441 = buffer.data(hsi1 + 441);
    const auto *hsi1_443 = buffer.data(hsi1 + 443);
    const auto *hsi1_444 = buffer.data(hsi1 + 444);
    const auto *hsi1_445 = buffer.data(hsi1 + 445);
    const auto *hsi1_447 = buffer.data(hsi1 + 447);
    const auto *hsi1_453 = buffer.data(hsi1 + 453);
    const auto *hsi1_457 = buffer.data(hsi1 + 457);
    const auto *hsi1_460 = buffer.data(hsi1 + 460);
    const auto *hsi1_462 = buffer.data(hsi1 + 462);
    const auto *hsi1_469 = buffer.data(hsi1 + 469);
    const auto *hsi1_471 = buffer.data(hsi1 + 471);
    const auto *hsi1_472 = buffer.data(hsi1 + 472);
    const auto *hsi1_473 = buffer.data(hsi1 + 473);
    const auto *hsi1_475 = buffer.data(hsi1 + 475);
    const auto *hsi1_476 = buffer.data(hsi1 + 476);
    const auto *hsi1_479 = buffer.data(hsi1 + 479);
    const auto *hsi1_481 = buffer.data(hsi1 + 481);
    const auto *hsi1_482 = buffer.data(hsi1 + 482);

    const auto *isg0_194 = buffer.data(isg0 + 194);
    const auto *isg0_205 = buffer.data(isg0 + 205);
    const auto *isg0_207 = buffer.data(isg0 + 207);
    const auto *isg0_208 = buffer.data(isg0 + 208);
    const auto *isg0_209 = buffer.data(isg0 + 209);
    const auto *isg0_210 = buffer.data(isg0 + 210);
    const auto *isg0_211 = buffer.data(isg0 + 211);
    const auto *isg0_212 = buffer.data(isg0 + 212);
    const auto *isg0_213 = buffer.data(isg0 + 213);
    const auto *isg0_214 = buffer.data(isg0 + 214);
    const auto *isg0_215 = buffer.data(isg0 + 215);
    const auto *isg0_219 = buffer.data(isg0 + 219);
    const auto *isg0_220 = buffer.data(isg0 + 220);
    const auto *isg0_221 = buffer.data(isg0 + 221);
    const auto *isg0_222 = buffer.data(isg0 + 222);
    const auto *isg0_223 = buffer.data(isg0 + 223);
    const auto *isg0_224 = buffer.data(isg0 + 224);
    const auto *isg0_225 = buffer.data(isg0 + 225);
    const auto *isg0_227 = buffer.data(isg0 + 227);
    const auto *isg0_228 = buffer.data(isg0 + 228);
    const auto *isg0_230 = buffer.data(isg0 + 230);

    const auto *isg1_194 = buffer.data(isg1 + 194);
    const auto *isg1_205 = buffer.data(isg1 + 205);
    const auto *isg1_207 = buffer.data(isg1 + 207);
    const auto *isg1_208 = buffer.data(isg1 + 208);
    const auto *isg1_209 = buffer.data(isg1 + 209);
    const auto *isg1_210 = buffer.data(isg1 + 210);
    const auto *isg1_211 = buffer.data(isg1 + 211);
    const auto *isg1_212 = buffer.data(isg1 + 212);
    const auto *isg1_213 = buffer.data(isg1 + 213);
    const auto *isg1_214 = buffer.data(isg1 + 214);
    const auto *isg1_215 = buffer.data(isg1 + 215);
    const auto *isg1_219 = buffer.data(isg1 + 219);
    const auto *isg1_220 = buffer.data(isg1 + 220);
    const auto *isg1_221 = buffer.data(isg1 + 221);
    const auto *isg1_222 = buffer.data(isg1 + 222);
    const auto *isg1_223 = buffer.data(isg1 + 223);
    const auto *isg1_224 = buffer.data(isg1 + 224);
    const auto *isg1_225 = buffer.data(isg1 + 225);
    const auto *isg1_227 = buffer.data(isg1 + 227);
    const auto *isg1_228 = buffer.data(isg1 + 228);
    const auto *isg1_230 = buffer.data(isg1 + 230);

    const auto *ish_272 = buffer.data(ish + 272);
    const auto *ish_273 = buffer.data(ish + 273);
    const auto *ish_275 = buffer.data(ish + 275);
    const auto *ish_276 = buffer.data(ish + 276);
    const auto *ish_278 = buffer.data(ish + 278);
    const auto *ish_279 = buffer.data(ish + 279);
    const auto *ish_282 = buffer.data(ish + 282);
    const auto *ish_288 = buffer.data(ish + 288);
    const auto *ish_289 = buffer.data(ish + 289);
    const auto *ish_290 = buffer.data(ish + 290);
    const auto *ish_291 = buffer.data(ish + 291);
    const auto *ish_292 = buffer.data(ish + 292);
    const auto *ish_293 = buffer.data(ish + 293);
    const auto *ish_294 = buffer.data(ish + 294);
    const auto *ish_295 = buffer.data(ish + 295);
    const auto *ish_296 = buffer.data(ish + 296);
    const auto *ish_297 = buffer.data(ish + 297);
    const auto *ish_298 = buffer.data(ish + 298);
    const auto *ish_299 = buffer.data(ish + 299);
    const auto *ish_300 = buffer.data(ish + 300);
    const auto *ish_301 = buffer.data(ish + 301);
    const auto *ish_302 = buffer.data(ish + 302);
    const auto *ish_303 = buffer.data(ish + 303);
    const auto *ish_308 = buffer.data(ish + 308);
    const auto *ish_309 = buffer.data(ish + 309);
    const auto *ish_310 = buffer.data(ish + 310);
    const auto *ish_311 = buffer.data(ish + 311);
    const auto *ish_312 = buffer.data(ish + 312);
    const auto *ish_313 = buffer.data(ish + 313);
    const auto *ish_314 = buffer.data(ish + 314);
    const auto *ish_315 = buffer.data(ish + 315);
    const auto *ish_316 = buffer.data(ish + 316);
    const auto *ish_317 = buffer.data(ish + 317);
    const auto *ish_318 = buffer.data(ish + 318);
    const auto *ish_320 = buffer.data(ish + 320);
    const auto *ish_321 = buffer.data(ish + 321);
    const auto *ish_322 = buffer.data(ish + 322);
    const auto *ish_324 = buffer.data(ish + 324);
    const auto *ish_325 = buffer.data(ish + 325);
    const auto *ish_330 = buffer.data(ish + 330);
    const auto *ish_332 = buffer.data(ish + 332);
    const auto *ish_333 = buffer.data(ish + 333);
    const auto *ish_334 = buffer.data(ish + 334);
    const auto *ish_335 = buffer.data(ish + 335);
    const auto *ish_336 = buffer.data(ish + 336);
    const auto *ish_338 = buffer.data(ish + 338);
    const auto *ish_339 = buffer.data(ish + 339);
    const auto *ish_341 = buffer.data(ish + 341);
    const auto *ish_342 = buffer.data(ish + 342);
    const auto *ish_345 = buffer.data(ish + 345);
    const auto *ish_351 = buffer.data(ish + 351);
    const auto *ish_352 = buffer.data(ish + 352);
    const auto *ish_353 = buffer.data(ish + 353);
    const auto *ish_354 = buffer.data(ish + 354);
    const auto *ish_355 = buffer.data(ish + 355);
    const auto *ish_356 = buffer.data(ish + 356);
    const auto *ish_357 = buffer.data(ish + 357);
    const auto *ish_359 = buffer.data(ish + 359);
    const auto *ish_360 = buffer.data(ish + 360);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, hsi0_252, hsh_167, \
                         hsh_188, hsh_189, hsi1_252, isg0_194, isg1_194, ish_272, \
                         ish_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * hsh_188[k]
                   + f_3 * pc_y[k] * ish_272[k];

        t_363[k] = f_12 * hsh_167[k]
                   + f_1 * isg0_194[k]
                   - f_2 * isg1_194[k]
                   + f_3 * pc_z[k] * ish_272[k];

        t_364[k] = pa_y[k] * hsi0_252[k]
                   - f_10 * pc_y[k] * hsi1_252[k];

        t_365[k] = f_11 * hsh_189[k]
                   + f_3 * pc_y[k] * ish_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, hsi0_255, hsi0_257, \
                         hsh_168, hsh_190, hsh_191, hsi1_255, hsi1_257, ish_273, \
                         ish_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * hsh_168[k]
                   + f_3 * pc_z[k] * ish_273[k];

        t_367[k] = pa_y[k] * hsi0_255[k]
                   + f_12 * hsh_190[k]
                   - f_10 * pc_y[k] * hsi1_255[k];

        t_368[k] = f_11 * hsh_191[k]
                   + f_3 * pc_y[k] * ish_275[k];

        t_369[k] = pa_y[k] * hsi0_257[k]
                   - f_10 * pc_y[k] * hsi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, hsi0_258, hsi0_261, \
                         hsh_171, hsh_192, hsh_194, hsi1_258, hsi1_261, ish_276, \
                         ish_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * hsi0_258[k]
                   + f_13 * hsh_192[k]
                   - f_10 * pc_y[k] * hsi1_258[k];

        t_371[k] = f_13 * hsh_171[k]
                   + f_3 * pc_z[k] * ish_276[k];

        t_372[k] = f_11 * hsh_194[k]
                   + f_3 * pc_y[k] * ish_278[k];

        t_373[k] = pa_y[k] * hsi0_261[k]
                   - f_10 * pc_y[k] * hsi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, hsi0_262, hsi0_264, hsh_174, \
                         hsh_195, hsh_197, hsi1_262, hsi1_264, \
                         ish_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * hsi0_262[k]
                   + f_14 * hsh_195[k]
                   - f_10 * pc_y[k] * hsi1_262[k];

        t_375[k] = f_13 * hsh_174[k]
                   + f_3 * pc_z[k] * ish_279[k];

        t_376[k] = pa_y[k] * hsi0_264[k]
                   + f_12 * hsh_197[k]
                   - f_10 * pc_y[k] * hsi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, hsi0_266, hsh_198, \
                         hsh_288, hsh_289, hsi1_266, ish_282, ish_288, \
                         ish_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * hsh_198[k]
                   + f_3 * pc_y[k] * ish_282[k];

        t_378[k] = pa_y[k] * hsi0_266[k]
                   - f_10 * pc_y[k] * hsi1_266[k];

        t_379[k] = f_12 * hsh_288[k]
                   + f_3 * pc_x[k] * ish_288[k];

        t_380[k] = f_12 * hsh_289[k]
                   + f_3 * pc_x[k] * ish_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, hsh_290, hsh_291, hsh_292, hsh_293, \
                         ish_290, ish_291, ish_292, ish_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_12 * hsh_290[k]
                   + f_3 * pc_x[k] * ish_290[k];

        t_382[k] = f_12 * hsh_291[k]
                   + f_3 * pc_x[k] * ish_291[k];

        t_383[k] = f_12 * hsh_292[k]
                   + f_3 * pc_x[k] * ish_292[k];

        t_384[k] = f_12 * hsh_293[k]
                   + f_3 * pc_x[k] * ish_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, hsh_183, hsh_204, hsh_206, isg0_205, \
                         isg0_207, isg1_205, isg1_207, ish_288, \
                         ish_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * hsh_204[k]
                   + f_1 * isg0_205[k]
                   - f_2 * isg1_205[k]
                   + f_3 * pc_y[k] * ish_288[k];

        t_386[k] = f_13 * hsh_183[k]
                   + f_3 * pc_z[k] * ish_288[k];

        t_387[k] = f_11 * hsh_206[k]
                   + f_8 * isg0_207[k]
                   - f_9 * isg1_207[k]
                   + f_3 * pc_y[k] * ish_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, hsh_207, hsh_208, hsh_209, isg0_208, \
                         isg0_209, isg1_208, isg1_209, ish_291, ish_292, \
                         ish_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * hsh_207[k]
                   + f_6 * isg0_208[k]
                   - f_7 * isg1_208[k]
                   + f_3 * pc_y[k] * ish_291[k];

        t_389[k] = f_11 * hsh_208[k]
                   + f_4 * isg0_209[k]
                   - f_5 * isg1_209[k]
                   + f_3 * pc_y[k] * ish_292[k];

        t_390[k] = f_11 * hsh_209[k]
                   + f_3 * pc_y[k] * ish_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, hsi0_279, \
                         hsh_189, hsh_294, hsi1_279, isg0_210, isg1_210, \
                         ish_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * hsi0_279[k]
                   - f_10 * pc_y[k] * hsi1_279[k];

        t_392[k] = f_12 * hsh_294[k]
                   + f_1 * isg0_210[k]
                   - f_2 * isg1_210[k]
                   + f_3 * pc_x[k] * ish_294[k];

        t_393[k] = f_3 * pc_y[k] * ish_294[k];

        t_394[k] = f_14 * hsh_189[k]
                   + f_3 * pc_z[k] * ish_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, hsh_299, isg0_210, isg0_215, \
                         isg1_210, isg1_215, ish_295, ish_296, \
                         ish_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * isg0_210[k]
                   - f_5 * isg1_210[k]
                   + f_3 * pc_y[k] * ish_295[k];

        t_396[k] = f_3 * pc_y[k] * ish_296[k];

        t_397[k] = f_12 * hsh_299[k]
                   + f_8 * isg0_215[k]
                   - f_9 * isg1_215[k]
                   + f_3 * pc_x[k] * ish_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, isg0_211, isg0_212, isg1_211, isg1_212, \
                         ish_297, ish_298, ish_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * isg0_211[k]
                   - f_7 * isg1_211[k]
                   + f_3 * pc_y[k] * ish_297[k];

        t_399[k] = f_4 * isg0_212[k]
                   - f_5 * isg1_212[k]
                   + f_3 * pc_y[k] * ish_298[k];

        t_400[k] = f_3 * pc_y[k] * ish_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, hsh_303, isg0_213, isg0_214, \
                         isg0_219, isg1_213, isg1_214, isg1_219, ish_300, ish_301, \
                         ish_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_12 * hsh_303[k]
                   + f_6 * isg0_219[k]
                   - f_7 * isg1_219[k]
                   + f_3 * pc_x[k] * ish_303[k];

        t_402[k] = f_8 * isg0_213[k]
                   - f_9 * isg1_213[k]
                   + f_3 * pc_y[k] * ish_300[k];

        t_403[k] = f_6 * isg0_214[k]
                   - f_7 * isg1_214[k]
                   + f_3 * pc_y[k] * ish_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, hsh_308, hsh_309, isg0_215, \
                         isg0_224, isg1_215, isg1_224, ish_302, ish_303, ish_308, \
                         ish_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * isg0_215[k]
                   - f_5 * isg1_215[k]
                   + f_3 * pc_y[k] * ish_302[k];

        t_405[k] = f_3 * pc_y[k] * ish_303[k];

        t_406[k] = f_12 * hsh_308[k]
                   + f_4 * isg0_224[k]
                   - f_5 * isg1_224[k]
                   + f_3 * pc_x[k] * ish_308[k];

        t_407[k] = f_12 * hsh_309[k]
                   + f_3 * pc_x[k] * ish_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, hsh_310, hsh_311, \
                         hsh_312, hsh_314, ish_308, ish_310, ish_311, ish_312, \
                         ish_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_12 * hsh_310[k]
                   + f_3 * pc_x[k] * ish_310[k];

        t_409[k] = f_12 * hsh_311[k]
                   + f_3 * pc_x[k] * ish_311[k];

        t_410[k] = f_12 * hsh_312[k]
                   + f_3 * pc_x[k] * ish_312[k];

        t_411[k] = f_3 * pc_y[k] * ish_308[k];

        t_412[k] = f_12 * hsh_314[k]
                   + f_3 * pc_x[k] * ish_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, isg0_220, isg0_221, isg0_222, isg1_220, \
                         isg1_221, isg1_222, ish_309, ish_310, \
                         ish_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * isg0_220[k]
                   - f_2 * isg1_220[k]
                   + f_3 * pc_y[k] * ish_309[k];

        t_414[k] = f_16 * isg0_221[k]
                   - f_17 * isg1_221[k]
                   + f_3 * pc_y[k] * ish_310[k];

        t_415[k] = f_8 * isg0_222[k]
                   - f_9 * isg1_222[k]
                   + f_3 * pc_y[k] * ish_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, hsh_209, isg0_223, isg0_224, \
                         isg1_223, isg1_224, ish_312, ish_313, \
                         ish_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * isg0_223[k]
                   - f_7 * isg1_223[k]
                   + f_3 * pc_y[k] * ish_312[k];

        t_417[k] = f_4 * isg0_224[k]
                   - f_5 * isg1_224[k]
                   + f_3 * pc_y[k] * ish_313[k];

        t_418[k] = f_3 * pc_y[k] * ish_314[k];

        t_419[k] = f_14 * hsh_209[k]
                   + f_1 * isg0_224[k]
                   - f_2 * isg1_224[k]
                   + f_3 * pc_z[k] * ish_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pa_x, pc_x, pc_y, pc_z, hsi0_420, \
                         hsi0_423, hsh_210, hsh_315, hsh_318, hsi1_420, hsi1_423, \
                         ish_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pa_x[k] * hsi0_420[k]
                   + f_0 * hsh_315[k]
                   - f_10 * pc_x[k] * hsi1_420[k];

        t_421[k] = f_15 * hsh_210[k]
                   + f_3 * pc_y[k] * ish_315[k];

        t_422[k] = f_3 * pc_z[k] * ish_315[k];

        t_423[k] = pa_x[k] * hsi0_423[k]
                   + f_14 * hsh_318[k]
                   - f_10 * pc_x[k] * hsi1_423[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pa_x, pc_x, pc_z, hsi0_426, hsh_321, \
                         hsi1_426, isg0_225, isg1_225, ish_316, ish_317, \
                         ish_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * ish_316[k];

        t_425[k] = f_4 * isg0_225[k]
                   - f_5 * isg1_225[k]
                   + f_3 * pc_z[k] * ish_317[k];

        t_426[k] = pa_x[k] * hsi0_426[k]
                   + f_13 * hsh_321[k]
                   - f_10 * pc_x[k] * hsi1_426[k];

        t_427[k] = f_3 * pc_z[k] * ish_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pa_x, pc_x, pc_y, pc_z, hsi0_430, \
                         hsh_215, hsh_325, hsi1_430, isg0_227, isg1_227, ish_320, \
                         ish_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_15 * hsh_215[k]
                   + f_3 * pc_y[k] * ish_320[k];

        t_429[k] = f_6 * isg0_227[k]
                   - f_7 * isg1_227[k]
                   + f_3 * pc_z[k] * ish_320[k];

        t_430[k] = pa_x[k] * hsi0_430[k]
                   + f_12 * hsh_325[k]
                   - f_10 * pc_x[k] * hsi1_430[k];

        t_431[k] = f_3 * pc_z[k] * ish_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, hsh_219, hsh_330, \
                         isg0_228, isg0_230, isg1_228, isg1_230, ish_322, ish_324, \
                         ish_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * isg0_228[k]
                   - f_5 * isg1_228[k]
                   + f_3 * pc_z[k] * ish_322[k];

        t_433[k] = f_15 * hsh_219[k]
                   + f_3 * pc_y[k] * ish_324[k];

        t_434[k] = f_8 * isg0_230[k]
                   - f_9 * isg1_230[k]
                   + f_3 * pc_z[k] * ish_324[k];

        t_435[k] = f_11 * hsh_330[k]
                   + f_3 * pc_x[k] * ish_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, hsh_332, hsh_333, \
                         hsh_334, hsh_335, ish_325, ish_332, ish_333, ish_334, \
                         ish_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * ish_325[k];

        t_437[k] = f_11 * hsh_332[k]
                   + f_3 * pc_x[k] * ish_332[k];

        t_438[k] = f_11 * hsh_333[k]
                   + f_3 * pc_x[k] * ish_333[k];

        t_439[k] = f_11 * hsh_334[k]
                   + f_3 * pc_x[k] * ish_334[k];

        t_440[k] = f_11 * hsh_335[k]
                   + f_3 * pc_x[k] * ish_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_x, pc_x, pc_z, hsi0_441, hsi0_443, \
                         hsi0_444, hsi1_441, hsi1_443, hsi1_444, \
                         ish_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_x[k] * hsi0_441[k]
                   - f_10 * pc_x[k] * hsi1_441[k];

        t_442[k] = f_3 * pc_z[k] * ish_330[k];

        t_443[k] = pa_x[k] * hsi0_443[k]
                   - f_10 * pc_x[k] * hsi1_443[k];

        t_444[k] = pa_x[k] * hsi0_444[k]
                   - f_10 * pc_x[k] * hsi1_444[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_x, pc_x, pc_y, hsi0_445, hsi0_447, hsh_230, \
                         hsi1_445, hsi1_447, ish_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_x[k] * hsi0_445[k]
                   - f_10 * pc_x[k] * hsi1_445[k];

        t_446[k] = f_15 * hsh_230[k]
                   + f_3 * pc_y[k] * ish_335[k];

        t_447[k] = pa_x[k] * hsi0_447[k]
                   - f_10 * pc_x[k] * hsi1_447[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pa_z, pc_y, pc_z, hsi0_280, hsi0_283, \
                         hsh_210, hsh_231, hsi1_280, hsi1_283, \
                         ish_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pa_z[k] * hsi0_280[k]
                   - f_10 * pc_z[k] * hsi1_280[k];

        t_449[k] = f_14 * hsh_231[k]
                   + f_3 * pc_y[k] * ish_336[k];

        t_450[k] = f_11 * hsh_210[k]
                   + f_3 * pc_z[k] * ish_336[k];

        t_451[k] = pa_z[k] * hsi0_283[k]
                   - f_10 * pc_z[k] * hsi1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pa_x, pa_z, pc_x, pc_y, pc_z, hsi0_286, \
                         hsi0_453, hsh_233, hsh_341, hsi1_286, hsi1_453, \
                         ish_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_14 * hsh_233[k]
                   + f_3 * pc_y[k] * ish_338[k];

        t_453[k] = pa_x[k] * hsi0_453[k]
                   + f_14 * hsh_341[k]
                   - f_10 * pc_x[k] * hsi1_453[k];

        t_454[k] = pa_z[k] * hsi0_286[k]
                   - f_10 * pc_z[k] * hsi1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pa_x, pc_x, pc_y, pc_z, hsi0_457, hsh_213, \
                         hsh_236, hsh_345, hsi1_457, ish_339, ish_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_11 * hsh_213[k]
                   + f_3 * pc_z[k] * ish_339[k];

        t_456[k] = f_14 * hsh_236[k]
                   + f_3 * pc_y[k] * ish_341[k];

        t_457[k] = pa_x[k] * hsi0_457[k]
                   + f_13 * hsh_345[k]
                   - f_10 * pc_x[k] * hsi1_457[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pa_x, pa_z, pc_x, pc_z, hsi0_290, hsi0_460, \
                         hsh_216, hsh_348, hsi1_290, hsi1_460, \
                         ish_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pa_z[k] * hsi0_290[k]
                   - f_10 * pc_z[k] * hsi1_290[k];

        t_459[k] = f_11 * hsh_216[k]
                   + f_3 * pc_z[k] * ish_342[k];

        t_460[k] = pa_x[k] * hsi0_460[k]
                   + f_12 * hsh_348[k]
                   - f_10 * pc_x[k] * hsi1_460[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pa_x, pc_x, pc_y, hsi0_462, hsh_240, \
                         hsh_350, hsh_351, hsh_352, hsi1_462, ish_345, ish_351, \
                         ish_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_14 * hsh_240[k]
                   + f_3 * pc_y[k] * ish_345[k];

        t_462[k] = pa_x[k] * hsi0_462[k]
                   + f_12 * hsh_350[k]
                   - f_10 * pc_x[k] * hsi1_462[k];

        t_463[k] = f_11 * hsh_351[k]
                   + f_3 * pc_x[k] * ish_351[k];

        t_464[k] = f_11 * hsh_352[k]
                   + f_3 * pc_x[k] * ish_352[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pc_x, hsh_353, hsh_354, hsh_355, hsh_356, \
                         ish_353, ish_354, ish_355, ish_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_11 * hsh_353[k]
                   + f_3 * pc_x[k] * ish_353[k];

        t_466[k] = f_11 * hsh_354[k]
                   + f_3 * pc_x[k] * ish_354[k];

        t_467[k] = f_11 * hsh_355[k]
                   + f_3 * pc_x[k] * ish_355[k];

        t_468[k] = f_11 * hsh_356[k]
                   + f_3 * pc_x[k] * ish_356[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, pa_x, pc_x, pc_z, hsi0_469, hsi0_471, \
                         hsi0_472, hsh_225, hsi1_469, hsi1_471, hsi1_472, \
                         ish_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = pa_x[k] * hsi0_469[k]
                   - f_10 * pc_x[k] * hsi1_469[k];

        t_470[k] = f_11 * hsh_225[k]
                   + f_3 * pc_z[k] * ish_351[k];

        t_471[k] = pa_x[k] * hsi0_471[k]
                   - f_10 * pc_x[k] * hsi1_471[k];

        t_472[k] = pa_x[k] * hsi0_472[k]
                   - f_10 * pc_x[k] * hsi1_472[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, pa_x, pc_x, pc_y, hsi0_473, hsi0_475, \
                         hsi0_476, hsh_251, hsh_357, hsi1_473, hsi1_475, hsi1_476, \
                         ish_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = pa_x[k] * hsi0_473[k]
                   - f_10 * pc_x[k] * hsi1_473[k];

        t_474[k] = f_14 * hsh_251[k]
                   + f_3 * pc_y[k] * ish_356[k];

        t_475[k] = pa_x[k] * hsi0_475[k]
                   - f_10 * pc_x[k] * hsi1_475[k];

        t_476[k] = pa_x[k] * hsi0_476[k]
                   + f_0 * hsh_357[k]
                   - f_10 * pc_x[k] * hsi1_476[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pa_x, pc_x, pc_y, pc_z, hsi0_479, \
                         hsh_231, hsh_252, hsh_254, hsh_360, hsi1_479, ish_357, \
                         ish_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_13 * hsh_252[k]
                   + f_3 * pc_y[k] * ish_357[k];

        t_478[k] = f_12 * hsh_231[k]
                   + f_3 * pc_z[k] * ish_357[k];

        t_479[k] = pa_x[k] * hsi0_479[k]
                   + f_14 * hsh_360[k]
                   - f_10 * pc_x[k] * hsi1_479[k];

        t_480[k] = f_13 * hsh_254[k]
                   + f_3 * pc_y[k] * ish_359[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, pa_x, pc_x, pc_z, hsi0_481, hsi0_482, hsh_234, \
                         hsh_362, hsh_363, hsi1_481, hsi1_482, \
                         ish_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = pa_x[k] * hsi0_481[k]
                   + f_14 * hsh_362[k]
                   - f_10 * pc_x[k] * hsi1_481[k];

        t_482[k] = pa_x[k] * hsi0_482[k]
                   + f_13 * hsh_363[k]
                   - f_10 * pc_x[k] * hsi1_482[k];

        t_483[k] = f_12 * hsh_234[k]
                   + f_3 * pc_z[k] * ish_360[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsi0_392 = buffer.data(hsi0 + 392);
    const auto *hsi0_397 = buffer.data(hsi0 + 397);
    const auto *hsi0_401 = buffer.data(hsi0 + 401);
    const auto *hsi0_406 = buffer.data(hsi0 + 406);
    const auto *hsi0_485 = buffer.data(hsi0 + 485);
    const auto *hsi0_486 = buffer.data(hsi0 + 486);
    const auto *hsi0_488 = buffer.data(hsi0 + 488);
    const auto *hsi0_490 = buffer.data(hsi0 + 490);
    const auto *hsi0_497 = buffer.data(hsi0 + 497);
    const auto *hsi0_499 = buffer.data(hsi0 + 499);
    const auto *hsi0_500 = buffer.data(hsi0 + 500);
    const auto *hsi0_501 = buffer.data(hsi0 + 501);
    const auto *hsi0_503 = buffer.data(hsi0 + 503);
    const auto *hsi0_504 = buffer.data(hsi0 + 504);
    const auto *hsi0_507 = buffer.data(hsi0 + 507);
    const auto *hsi0_509 = buffer.data(hsi0 + 509);
    const auto *hsi0_510 = buffer.data(hsi0 + 510);
    const auto *hsi0_513 = buffer.data(hsi0 + 513);
    const auto *hsi0_514 = buffer.data(hsi0 + 514);
    const auto *hsi0_516 = buffer.data(hsi0 + 516);
    const auto *hsi0_518 = buffer.data(hsi0 + 518);
    const auto *hsi0_525 = buffer.data(hsi0 + 525);
    const auto *hsi0_527 = buffer.data(hsi0 + 527);
    const auto *hsi0_528 = buffer.data(hsi0 + 528);
    const auto *hsi0_529 = buffer.data(hsi0 + 529);
    const auto *hsi0_531 = buffer.data(hsi0 + 531);
    const auto *hsi0_535 = buffer.data(hsi0 + 535);
    const auto *hsi0_538 = buffer.data(hsi0 + 538);
    const auto *hsi0_542 = buffer.data(hsi0 + 542);
    const auto *hsi0_544 = buffer.data(hsi0 + 544);
    const auto *hsi0_553 = buffer.data(hsi0 + 553);
    const auto *hsi0_555 = buffer.data(hsi0 + 555);
    const auto *hsi0_556 = buffer.data(hsi0 + 556);
    const auto *hsi0_557 = buffer.data(hsi0 + 557);
    const auto *hsi0_559 = buffer.data(hsi0 + 559);
    const auto *hsi0_560 = buffer.data(hsi0 + 560);
    const auto *hsi0_565 = buffer.data(hsi0 + 565);
    const auto *hsi0_569 = buffer.data(hsi0 + 569);
    const auto *hsi0_574 = buffer.data(hsi0 + 574);
    const auto *hsi0_581 = buffer.data(hsi0 + 581);
    const auto *hsi0_582 = buffer.data(hsi0 + 582);
    const auto *hsi0_583 = buffer.data(hsi0 + 583);
    const auto *hsi0_584 = buffer.data(hsi0 + 584);
    const auto *hsi0_585 = buffer.data(hsi0 + 585);
    const auto *hsi0_587 = buffer.data(hsi0 + 587);

    const auto *hsh_237 = buffer.data(hsh + 237);
    const auto *hsh_246 = buffer.data(hsh + 246);
    const auto *hsh_252 = buffer.data(hsh + 252);
    const auto *hsh_255 = buffer.data(hsh + 255);
    const auto *hsh_257 = buffer.data(hsh + 257);
    const auto *hsh_258 = buffer.data(hsh + 258);
    const auto *hsh_261 = buffer.data(hsh + 261);
    const auto *hsh_267 = buffer.data(hsh + 267);
    const auto *hsh_272 = buffer.data(hsh + 272);
    const auto *hsh_273 = buffer.data(hsh + 273);
    const auto *hsh_275 = buffer.data(hsh + 275);
    const auto *hsh_276 = buffer.data(hsh + 276);
    const auto *hsh_278 = buffer.data(hsh + 278);
    const auto *hsh_279 = buffer.data(hsh + 279);
    const auto *hsh_282 = buffer.data(hsh + 282);
    const auto *hsh_288 = buffer.data(hsh + 288);
    const auto *hsh_293 = buffer.data(hsh + 293);
    const auto *hsh_294 = buffer.data(hsh + 294);
    const auto *hsh_296 = buffer.data(hsh + 296);
    const auto *hsh_299 = buffer.data(hsh + 299);
    const auto *hsh_303 = buffer.data(hsh + 303);
    const auto *hsh_314 = buffer.data(hsh + 314);
    const auto *hsh_366 = buffer.data(hsh + 366);
    const auto *hsh_367 = buffer.data(hsh + 367);
    const auto *hsh_369 = buffer.data(hsh + 369);
    const auto *hsh_371 = buffer.data(hsh + 371);
    const auto *hsh_372 = buffer.data(hsh + 372);
    const auto *hsh_373 = buffer.data(hsh + 373);
    const auto *hsh_374 = buffer.data(hsh + 374);
    const auto *hsh_375 = buffer.data(hsh + 375);
    const auto *hsh_376 = buffer.data(hsh + 376);
    const auto *hsh_377 = buffer.data(hsh + 377);
    const auto *hsh_378 = buffer.data(hsh + 378);
    const auto *hsh_381 = buffer.data(hsh + 381);
    const auto *hsh_383 = buffer.data(hsh + 383);
    const auto *hsh_384 = buffer.data(hsh + 384);
    const auto *hsh_387 = buffer.data(hsh + 387);
    const auto *hsh_388 = buffer.data(hsh + 388);
    const auto *hsh_390 = buffer.data(hsh + 390);
    const auto *hsh_392 = buffer.data(hsh + 392);
    const auto *hsh_393 = buffer.data(hsh + 393);
    const auto *hsh_394 = buffer.data(hsh + 394);
    const auto *hsh_395 = buffer.data(hsh + 395);
    const auto *hsh_396 = buffer.data(hsh + 396);
    const auto *hsh_397 = buffer.data(hsh + 397);
    const auto *hsh_398 = buffer.data(hsh + 398);
    const auto *hsh_402 = buffer.data(hsh + 402);
    const auto *hsh_405 = buffer.data(hsh + 405);
    const auto *hsh_409 = buffer.data(hsh + 409);
    const auto *hsh_411 = buffer.data(hsh + 411);
    const auto *hsh_414 = buffer.data(hsh + 414);
    const auto *hsh_415 = buffer.data(hsh + 415);
    const auto *hsh_416 = buffer.data(hsh + 416);
    const auto *hsh_417 = buffer.data(hsh + 417);
    const auto *hsh_418 = buffer.data(hsh + 418);
    const auto *hsh_419 = buffer.data(hsh + 419);
    const auto *hsh_420 = buffer.data(hsh + 420);
    const auto *hsh_425 = buffer.data(hsh + 425);
    const auto *hsh_429 = buffer.data(hsh + 429);
    const auto *hsh_434 = buffer.data(hsh + 434);
    const auto *hsh_435 = buffer.data(hsh + 435);
    const auto *hsh_436 = buffer.data(hsh + 436);
    const auto *hsh_437 = buffer.data(hsh + 437);
    const auto *hsh_438 = buffer.data(hsh + 438);
    const auto *hsh_440 = buffer.data(hsh + 440);

    const auto *hsi1_392 = buffer.data(hsi1 + 392);
    const auto *hsi1_397 = buffer.data(hsi1 + 397);
    const auto *hsi1_401 = buffer.data(hsi1 + 401);
    const auto *hsi1_406 = buffer.data(hsi1 + 406);
    const auto *hsi1_485 = buffer.data(hsi1 + 485);
    const auto *hsi1_486 = buffer.data(hsi1 + 486);
    const auto *hsi1_488 = buffer.data(hsi1 + 488);
    const auto *hsi1_490 = buffer.data(hsi1 + 490);
    const auto *hsi1_497 = buffer.data(hsi1 + 497);
    const auto *hsi1_499 = buffer.data(hsi1 + 499);
    const auto *hsi1_500 = buffer.data(hsi1 + 500);
    const auto *hsi1_501 = buffer.data(hsi1 + 501);
    const auto *hsi1_503 = buffer.data(hsi1 + 503);
    const auto *hsi1_504 = buffer.data(hsi1 + 504);
    const auto *hsi1_507 = buffer.data(hsi1 + 507);
    const auto *hsi1_509 = buffer.data(hsi1 + 509);
    const auto *hsi1_510 = buffer.data(hsi1 + 510);
    const auto *hsi1_513 = buffer.data(hsi1 + 513);
    const auto *hsi1_514 = buffer.data(hsi1 + 514);
    const auto *hsi1_516 = buffer.data(hsi1 + 516);
    const auto *hsi1_518 = buffer.data(hsi1 + 518);
    const auto *hsi1_525 = buffer.data(hsi1 + 525);
    const auto *hsi1_527 = buffer.data(hsi1 + 527);
    const auto *hsi1_528 = buffer.data(hsi1 + 528);
    const auto *hsi1_529 = buffer.data(hsi1 + 529);
    const auto *hsi1_531 = buffer.data(hsi1 + 531);
    const auto *hsi1_535 = buffer.data(hsi1 + 535);
    const auto *hsi1_538 = buffer.data(hsi1 + 538);
    const auto *hsi1_542 = buffer.data(hsi1 + 542);
    const auto *hsi1_544 = buffer.data(hsi1 + 544);
    const auto *hsi1_553 = buffer.data(hsi1 + 553);
    const auto *hsi1_555 = buffer.data(hsi1 + 555);
    const auto *hsi1_556 = buffer.data(hsi1 + 556);
    const auto *hsi1_557 = buffer.data(hsi1 + 557);
    const auto *hsi1_559 = buffer.data(hsi1 + 559);
    const auto *hsi1_560 = buffer.data(hsi1 + 560);
    const auto *hsi1_565 = buffer.data(hsi1 + 565);
    const auto *hsi1_569 = buffer.data(hsi1 + 569);
    const auto *hsi1_574 = buffer.data(hsi1 + 574);
    const auto *hsi1_581 = buffer.data(hsi1 + 581);
    const auto *hsi1_582 = buffer.data(hsi1 + 582);
    const auto *hsi1_583 = buffer.data(hsi1 + 583);
    const auto *hsi1_584 = buffer.data(hsi1 + 584);
    const auto *hsi1_585 = buffer.data(hsi1 + 585);
    const auto *hsi1_587 = buffer.data(hsi1 + 587);

    const auto *isg0_300 = buffer.data(isg0 + 300);
    const auto *isg0_301 = buffer.data(isg0 + 301);
    const auto *isg0_302 = buffer.data(isg0 + 302);
    const auto *isg0_303 = buffer.data(isg0 + 303);
    const auto *isg0_304 = buffer.data(isg0 + 304);
    const auto *isg0_305 = buffer.data(isg0 + 305);
    const auto *isg0_315 = buffer.data(isg0 + 315);
    const auto *isg0_316 = buffer.data(isg0 + 316);
    const auto *isg0_318 = buffer.data(isg0 + 318);
    const auto *isg0_320 = buffer.data(isg0 + 320);
    const auto *isg0_321 = buffer.data(isg0 + 321);
    const auto *isg0_323 = buffer.data(isg0 + 323);
    const auto *isg0_324 = buffer.data(isg0 + 324);
    const auto *isg0_325 = buffer.data(isg0 + 325);
    const auto *isg0_327 = buffer.data(isg0 + 327);
    const auto *isg0_328 = buffer.data(isg0 + 328);
    const auto *isg0_329 = buffer.data(isg0 + 329);

    const auto *isg1_300 = buffer.data(isg1 + 300);
    const auto *isg1_301 = buffer.data(isg1 + 301);
    const auto *isg1_302 = buffer.data(isg1 + 302);
    const auto *isg1_303 = buffer.data(isg1 + 303);
    const auto *isg1_304 = buffer.data(isg1 + 304);
    const auto *isg1_305 = buffer.data(isg1 + 305);
    const auto *isg1_315 = buffer.data(isg1 + 315);
    const auto *isg1_316 = buffer.data(isg1 + 316);
    const auto *isg1_318 = buffer.data(isg1 + 318);
    const auto *isg1_320 = buffer.data(isg1 + 320);
    const auto *isg1_321 = buffer.data(isg1 + 321);
    const auto *isg1_323 = buffer.data(isg1 + 323);
    const auto *isg1_324 = buffer.data(isg1 + 324);
    const auto *isg1_325 = buffer.data(isg1 + 325);
    const auto *isg1_327 = buffer.data(isg1 + 327);
    const auto *isg1_328 = buffer.data(isg1 + 328);
    const auto *isg1_329 = buffer.data(isg1 + 329);

    const auto *ish_362 = buffer.data(ish + 362);
    const auto *ish_363 = buffer.data(ish + 363);
    const auto *ish_366 = buffer.data(ish + 366);
    const auto *ish_372 = buffer.data(ish + 372);
    const auto *ish_373 = buffer.data(ish + 373);
    const auto *ish_374 = buffer.data(ish + 374);
    const auto *ish_375 = buffer.data(ish + 375);
    const auto *ish_376 = buffer.data(ish + 376);
    const auto *ish_377 = buffer.data(ish + 377);
    const auto *ish_378 = buffer.data(ish + 378);
    const auto *ish_380 = buffer.data(ish + 380);
    const auto *ish_381 = buffer.data(ish + 381);
    const auto *ish_383 = buffer.data(ish + 383);
    const auto *ish_384 = buffer.data(ish + 384);
    const auto *ish_387 = buffer.data(ish + 387);
    const auto *ish_393 = buffer.data(ish + 393);
    const auto *ish_394 = buffer.data(ish + 394);
    const auto *ish_395 = buffer.data(ish + 395);
    const auto *ish_396 = buffer.data(ish + 396);
    const auto *ish_397 = buffer.data(ish + 397);
    const auto *ish_398 = buffer.data(ish + 398);
    const auto *ish_399 = buffer.data(ish + 399);
    const auto *ish_401 = buffer.data(ish + 401);
    const auto *ish_402 = buffer.data(ish + 402);
    const auto *ish_404 = buffer.data(ish + 404);
    const auto *ish_405 = buffer.data(ish + 405);
    const auto *ish_408 = buffer.data(ish + 408);
    const auto *ish_414 = buffer.data(ish + 414);
    const auto *ish_415 = buffer.data(ish + 415);
    const auto *ish_416 = buffer.data(ish + 416);
    const auto *ish_417 = buffer.data(ish + 417);
    const auto *ish_418 = buffer.data(ish + 418);
    const auto *ish_419 = buffer.data(ish + 419);
    const auto *ish_420 = buffer.data(ish + 420);
    const auto *ish_421 = buffer.data(ish + 421);
    const auto *ish_422 = buffer.data(ish + 422);
    const auto *ish_423 = buffer.data(ish + 423);
    const auto *ish_424 = buffer.data(ish + 424);
    const auto *ish_425 = buffer.data(ish + 425);
    const auto *ish_426 = buffer.data(ish + 426);
    const auto *ish_427 = buffer.data(ish + 427);
    const auto *ish_428 = buffer.data(ish + 428);
    const auto *ish_429 = buffer.data(ish + 429);
    const auto *ish_434 = buffer.data(ish + 434);
    const auto *ish_435 = buffer.data(ish + 435);
    const auto *ish_436 = buffer.data(ish + 436);
    const auto *ish_437 = buffer.data(ish + 437);
    const auto *ish_438 = buffer.data(ish + 438);
    const auto *ish_440 = buffer.data(ish + 440);
    const auto *ish_441 = buffer.data(ish + 441);
    const auto *ish_442 = buffer.data(ish + 442);
    const auto *ish_444 = buffer.data(ish + 444);
    const auto *ish_446 = buffer.data(ish + 446);
    const auto *ish_447 = buffer.data(ish + 447);
    const auto *ish_449 = buffer.data(ish + 449);
    const auto *ish_450 = buffer.data(ish + 450);
    const auto *ish_451 = buffer.data(ish + 451);
    const auto *ish_453 = buffer.data(ish + 453);
    const auto *ish_454 = buffer.data(ish + 454);
    const auto *ish_455 = buffer.data(ish + 455);
    const auto *ish_456 = buffer.data(ish + 456);
    const auto *ish_457 = buffer.data(ish + 457);
    const auto *ish_458 = buffer.data(ish + 458);
    const auto *ish_459 = buffer.data(ish + 459);
    const auto *ish_460 = buffer.data(ish + 460);

#pragma omp simd aligned(t_484, t_485, t_486, pa_x, pc_x, pc_y, hsi0_485, hsi0_486, hsh_257, \
                         hsh_366, hsh_367, hsi1_485, hsi1_486, \
                         ish_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_13 * hsh_257[k]
                   + f_3 * pc_y[k] * ish_362[k];

        t_485[k] = pa_x[k] * hsi0_485[k]
                   + f_13 * hsh_366[k]
                   - f_10 * pc_x[k] * hsi1_485[k];

        t_486[k] = pa_x[k] * hsi0_486[k]
                   + f_12 * hsh_367[k]
                   - f_10 * pc_x[k] * hsi1_486[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, pa_x, pc_x, pc_y, pc_z, hsi0_488, hsh_237, \
                         hsh_261, hsh_369, hsi1_488, ish_363, ish_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_12 * hsh_237[k]
                   + f_3 * pc_z[k] * ish_363[k];

        t_488[k] = pa_x[k] * hsi0_488[k]
                   + f_12 * hsh_369[k]
                   - f_10 * pc_x[k] * hsi1_488[k];

        t_489[k] = f_13 * hsh_261[k]
                   + f_3 * pc_y[k] * ish_366[k];
    }

#pragma omp simd aligned(t_490, t_491, t_492, t_493, pa_x, pc_x, hsi0_490, hsh_371, hsh_372, \
                         hsh_373, hsh_374, hsi1_490, ish_372, ish_373, \
                         ish_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = pa_x[k] * hsi0_490[k]
                   + f_12 * hsh_371[k]
                   - f_10 * pc_x[k] * hsi1_490[k];

        t_491[k] = f_11 * hsh_372[k]
                   + f_3 * pc_x[k] * ish_372[k];

        t_492[k] = f_11 * hsh_373[k]
                   + f_3 * pc_x[k] * ish_373[k];

        t_493[k] = f_11 * hsh_374[k]
                   + f_3 * pc_x[k] * ish_374[k];
    }

#pragma omp simd aligned(t_494, t_495, t_496, t_497, pa_x, pc_x, hsi0_497, hsh_375, hsh_376, \
                         hsh_377, hsi1_497, ish_375, ish_376, ish_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_11 * hsh_375[k]
                   + f_3 * pc_x[k] * ish_375[k];

        t_495[k] = f_11 * hsh_376[k]
                   + f_3 * pc_x[k] * ish_376[k];

        t_496[k] = f_11 * hsh_377[k]
                   + f_3 * pc_x[k] * ish_377[k];

        t_497[k] = pa_x[k] * hsi0_497[k]
                   - f_10 * pc_x[k] * hsi1_497[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, t_501, pa_x, pc_x, pc_z, hsi0_499, hsi0_500, \
                         hsi0_501, hsh_246, hsi1_499, hsi1_500, hsi1_501, \
                         ish_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_12 * hsh_246[k]
                   + f_3 * pc_z[k] * ish_372[k];

        t_499[k] = pa_x[k] * hsi0_499[k]
                   - f_10 * pc_x[k] * hsi1_499[k];

        t_500[k] = pa_x[k] * hsi0_500[k]
                   - f_10 * pc_x[k] * hsi1_500[k];

        t_501[k] = pa_x[k] * hsi0_501[k]
                   - f_10 * pc_x[k] * hsi1_501[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, t_505, pa_x, pc_x, pc_y, hsi0_503, hsi0_504, \
                         hsh_272, hsh_273, hsh_378, hsi1_503, hsi1_504, ish_377, \
                         ish_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * hsh_272[k]
                   + f_3 * pc_y[k] * ish_377[k];

        t_503[k] = pa_x[k] * hsi0_503[k]
                   - f_10 * pc_x[k] * hsi1_503[k];

        t_504[k] = pa_x[k] * hsi0_504[k]
                   + f_0 * hsh_378[k]
                   - f_10 * pc_x[k] * hsi1_504[k];

        t_505[k] = f_12 * hsh_273[k]
                   + f_3 * pc_y[k] * ish_378[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pa_x, pc_x, pc_y, pc_z, hsi0_507, hsh_252, \
                         hsh_275, hsh_381, hsi1_507, ish_378, ish_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_13 * hsh_252[k]
                   + f_3 * pc_z[k] * ish_378[k];

        t_507[k] = pa_x[k] * hsi0_507[k]
                   + f_14 * hsh_381[k]
                   - f_10 * pc_x[k] * hsi1_507[k];

        t_508[k] = f_12 * hsh_275[k]
                   + f_3 * pc_y[k] * ish_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pa_x, pc_x, pc_z, hsi0_509, hsi0_510, hsh_255, \
                         hsh_383, hsh_384, hsi1_509, hsi1_510, \
                         ish_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pa_x[k] * hsi0_509[k]
                   + f_14 * hsh_383[k]
                   - f_10 * pc_x[k] * hsi1_509[k];

        t_510[k] = pa_x[k] * hsi0_510[k]
                   + f_13 * hsh_384[k]
                   - f_10 * pc_x[k] * hsi1_510[k];

        t_511[k] = f_13 * hsh_255[k]
                   + f_3 * pc_z[k] * ish_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pa_x, pc_x, pc_y, hsi0_513, hsi0_514, hsh_278, \
                         hsh_387, hsh_388, hsi1_513, hsi1_514, \
                         ish_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * hsh_278[k]
                   + f_3 * pc_y[k] * ish_383[k];

        t_513[k] = pa_x[k] * hsi0_513[k]
                   + f_13 * hsh_387[k]
                   - f_10 * pc_x[k] * hsi1_513[k];

        t_514[k] = pa_x[k] * hsi0_514[k]
                   + f_12 * hsh_388[k]
                   - f_10 * pc_x[k] * hsi1_514[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pa_x, pc_x, pc_y, pc_z, hsi0_516, hsh_258, \
                         hsh_282, hsh_390, hsi1_516, ish_384, ish_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * hsh_258[k]
                   + f_3 * pc_z[k] * ish_384[k];

        t_516[k] = pa_x[k] * hsi0_516[k]
                   + f_12 * hsh_390[k]
                   - f_10 * pc_x[k] * hsi1_516[k];

        t_517[k] = f_12 * hsh_282[k]
                   + f_3 * pc_y[k] * ish_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pa_x, pc_x, hsi0_518, hsh_392, hsh_393, \
                         hsh_394, hsh_395, hsi1_518, ish_393, ish_394, \
                         ish_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pa_x[k] * hsi0_518[k]
                   + f_12 * hsh_392[k]
                   - f_10 * pc_x[k] * hsi1_518[k];

        t_519[k] = f_11 * hsh_393[k]
                   + f_3 * pc_x[k] * ish_393[k];

        t_520[k] = f_11 * hsh_394[k]
                   + f_3 * pc_x[k] * ish_394[k];

        t_521[k] = f_11 * hsh_395[k]
                   + f_3 * pc_x[k] * ish_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pa_x, pc_x, hsi0_525, hsh_396, hsh_397, \
                         hsh_398, hsi1_525, ish_396, ish_397, ish_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_11 * hsh_396[k]
                   + f_3 * pc_x[k] * ish_396[k];

        t_523[k] = f_11 * hsh_397[k]
                   + f_3 * pc_x[k] * ish_397[k];

        t_524[k] = f_11 * hsh_398[k]
                   + f_3 * pc_x[k] * ish_398[k];

        t_525[k] = pa_x[k] * hsi0_525[k]
                   - f_10 * pc_x[k] * hsi1_525[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pa_x, pc_x, pc_z, hsi0_527, hsi0_528, \
                         hsi0_529, hsh_267, hsi1_527, hsi1_528, hsi1_529, \
                         ish_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * hsh_267[k]
                   + f_3 * pc_z[k] * ish_393[k];

        t_527[k] = pa_x[k] * hsi0_527[k]
                   - f_10 * pc_x[k] * hsi1_527[k];

        t_528[k] = pa_x[k] * hsi0_528[k]
                   - f_10 * pc_x[k] * hsi1_528[k];

        t_529[k] = pa_x[k] * hsi0_529[k]
                   - f_10 * pc_x[k] * hsi1_529[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_x, pa_y, pc_x, pc_y, hsi0_392, \
                         hsi0_531, hsh_293, hsh_294, hsi1_392, hsi1_531, ish_398, \
                         ish_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_12 * hsh_293[k]
                   + f_3 * pc_y[k] * ish_398[k];

        t_531[k] = pa_x[k] * hsi0_531[k]
                   - f_10 * pc_x[k] * hsi1_531[k];

        t_532[k] = pa_y[k] * hsi0_392[k]
                   - f_10 * pc_y[k] * hsi1_392[k];

        t_533[k] = f_11 * hsh_294[k]
                   + f_3 * pc_y[k] * ish_399[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pa_x, pc_x, pc_y, pc_z, hsi0_535, hsh_273, \
                         hsh_296, hsh_402, hsi1_535, ish_399, ish_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_14 * hsh_273[k]
                   + f_3 * pc_z[k] * ish_399[k];

        t_535[k] = pa_x[k] * hsi0_535[k]
                   + f_14 * hsh_402[k]
                   - f_10 * pc_x[k] * hsi1_535[k];

        t_536[k] = f_11 * hsh_296[k]
                   + f_3 * pc_y[k] * ish_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_x, pa_y, pc_x, pc_y, pc_z, hsi0_397, \
                         hsi0_538, hsh_276, hsh_405, hsi1_397, hsi1_538, \
                         ish_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * hsi0_397[k]
                   - f_10 * pc_y[k] * hsi1_397[k];

        t_538[k] = pa_x[k] * hsi0_538[k]
                   + f_13 * hsh_405[k]
                   - f_10 * pc_x[k] * hsi1_538[k];

        t_539[k] = f_14 * hsh_276[k]
                   + f_3 * pc_z[k] * ish_402[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pa_x, pa_y, pc_x, pc_y, hsi0_401, hsi0_542, \
                         hsh_299, hsh_409, hsi1_401, hsi1_542, \
                         ish_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = f_11 * hsh_299[k]
                   + f_3 * pc_y[k] * ish_404[k];

        t_541[k] = pa_y[k] * hsi0_401[k]
                   - f_10 * pc_y[k] * hsi1_401[k];

        t_542[k] = pa_x[k] * hsi0_542[k]
                   + f_12 * hsh_409[k]
                   - f_10 * pc_x[k] * hsi1_542[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pa_x, pc_x, pc_y, pc_z, hsi0_544, hsh_279, \
                         hsh_303, hsh_411, hsi1_544, ish_405, ish_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = f_14 * hsh_279[k]
                   + f_3 * pc_z[k] * ish_405[k];

        t_544[k] = pa_x[k] * hsi0_544[k]
                   + f_12 * hsh_411[k]
                   - f_10 * pc_x[k] * hsi1_544[k];

        t_545[k] = f_11 * hsh_303[k]
                   + f_3 * pc_y[k] * ish_408[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pa_y, pc_x, pc_y, hsi0_406, hsh_414, \
                         hsh_415, hsh_416, hsi1_406, ish_414, ish_415, \
                         ish_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * hsi0_406[k]
                   - f_10 * pc_y[k] * hsi1_406[k];

        t_547[k] = f_11 * hsh_414[k]
                   + f_3 * pc_x[k] * ish_414[k];

        t_548[k] = f_11 * hsh_415[k]
                   + f_3 * pc_x[k] * ish_415[k];

        t_549[k] = f_11 * hsh_416[k]
                   + f_3 * pc_x[k] * ish_416[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pa_x, pc_x, hsi0_553, hsh_417, hsh_418, \
                         hsh_419, hsi1_553, ish_417, ish_418, ish_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_11 * hsh_417[k]
                   + f_3 * pc_x[k] * ish_417[k];

        t_551[k] = f_11 * hsh_418[k]
                   + f_3 * pc_x[k] * ish_418[k];

        t_552[k] = f_11 * hsh_419[k]
                   + f_3 * pc_x[k] * ish_419[k];

        t_553[k] = pa_x[k] * hsi0_553[k]
                   - f_10 * pc_x[k] * hsi1_553[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, pa_x, pc_x, pc_z, hsi0_555, hsi0_556, \
                         hsi0_557, hsh_288, hsi1_555, hsi1_556, hsi1_557, \
                         ish_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_14 * hsh_288[k]
                   + f_3 * pc_z[k] * ish_414[k];

        t_555[k] = pa_x[k] * hsi0_555[k]
                   - f_10 * pc_x[k] * hsi1_555[k];

        t_556[k] = pa_x[k] * hsi0_556[k]
                   - f_10 * pc_x[k] * hsi1_556[k];

        t_557[k] = pa_x[k] * hsi0_557[k]
                   - f_10 * pc_x[k] * hsi1_557[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pa_x, pc_x, pc_y, hsi0_559, hsi0_560, \
                         hsh_314, hsh_420, hsi1_559, hsi1_560, ish_419, \
                         ish_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_11 * hsh_314[k]
                   + f_3 * pc_y[k] * ish_419[k];

        t_559[k] = pa_x[k] * hsi0_559[k]
                   - f_10 * pc_x[k] * hsi1_559[k];

        t_560[k] = pa_x[k] * hsi0_560[k]
                   + f_0 * hsh_420[k]
                   - f_10 * pc_x[k] * hsi1_560[k];

        t_561[k] = f_3 * pc_y[k] * ish_420[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, pc_y, pc_z, hsh_294, isg0_300, isg1_300, \
                         ish_420, ish_421, ish_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_15 * hsh_294[k]
                   + f_3 * pc_z[k] * ish_420[k];

        t_563[k] = f_4 * isg0_300[k]
                   - f_5 * isg1_300[k]
                   + f_3 * pc_y[k] * ish_421[k];

        t_564[k] = f_3 * pc_y[k] * ish_422[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pa_x, pc_x, pc_y, hsi0_565, hsh_425, hsi1_565, \
                         isg0_301, isg0_302, isg1_301, isg1_302, ish_423, \
                         ish_424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pa_x[k] * hsi0_565[k]
                   + f_14 * hsh_425[k]
                   - f_10 * pc_x[k] * hsi1_565[k];

        t_566[k] = f_6 * isg0_301[k]
                   - f_7 * isg1_301[k]
                   + f_3 * pc_y[k] * ish_423[k];

        t_567[k] = f_4 * isg0_302[k]
                   - f_5 * isg1_302[k]
                   + f_3 * pc_y[k] * ish_424[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pa_x, pc_x, pc_y, hsi0_569, hsh_429, hsi1_569, \
                         isg0_303, isg1_303, ish_425, ish_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_3 * pc_y[k] * ish_425[k];

        t_569[k] = pa_x[k] * hsi0_569[k]
                   + f_13 * hsh_429[k]
                   - f_10 * pc_x[k] * hsi1_569[k];

        t_570[k] = f_8 * isg0_303[k]
                   - f_9 * isg1_303[k]
                   + f_3 * pc_y[k] * ish_426[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, pc_y, isg0_304, isg0_305, isg1_304, isg1_305, \
                         ish_427, ish_428, ish_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_6 * isg0_304[k]
                   - f_7 * isg1_304[k]
                   + f_3 * pc_y[k] * ish_427[k];

        t_572[k] = f_4 * isg0_305[k]
                   - f_5 * isg1_305[k]
                   + f_3 * pc_y[k] * ish_428[k];

        t_573[k] = f_3 * pc_y[k] * ish_429[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_x, pc_x, hsi0_574, hsh_434, hsh_435, \
                         hsh_436, hsh_437, hsi1_574, ish_435, ish_436, \
                         ish_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = pa_x[k] * hsi0_574[k]
                   + f_12 * hsh_434[k]
                   - f_10 * pc_x[k] * hsi1_574[k];

        t_575[k] = f_11 * hsh_435[k]
                   + f_3 * pc_x[k] * ish_435[k];

        t_576[k] = f_11 * hsh_436[k]
                   + f_3 * pc_x[k] * ish_436[k];

        t_577[k] = f_11 * hsh_437[k]
                   + f_3 * pc_x[k] * ish_437[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pa_x, pc_x, pc_y, hsi0_581, hsh_438, \
                         hsh_440, hsi1_581, ish_434, ish_438, ish_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_11 * hsh_438[k]
                   + f_3 * pc_x[k] * ish_438[k];

        t_579[k] = f_3 * pc_y[k] * ish_434[k];

        t_580[k] = f_11 * hsh_440[k]
                   + f_3 * pc_x[k] * ish_440[k];

        t_581[k] = pa_x[k] * hsi0_581[k]
                   - f_10 * pc_x[k] * hsi1_581[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pa_x, pc_x, hsi0_582, hsi0_583, hsi0_584, \
                         hsi0_585, hsi1_582, hsi1_583, hsi1_584, \
                         hsi1_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = pa_x[k] * hsi0_582[k]
                   - f_10 * pc_x[k] * hsi1_582[k];

        t_583[k] = pa_x[k] * hsi0_583[k]
                   - f_10 * pc_x[k] * hsi1_583[k];

        t_584[k] = pa_x[k] * hsi0_584[k]
                   - f_10 * pc_x[k] * hsi1_584[k];

        t_585[k] = pa_x[k] * hsi0_585[k]
                   - f_10 * pc_x[k] * hsi1_585[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pa_x, pc_x, pc_y, hsi0_587, hsi1_587, \
                         isg0_315, isg0_316, isg1_315, isg1_316, ish_440, ish_441, \
                         ish_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = f_3 * pc_y[k] * ish_440[k];

        t_587[k] = pa_x[k] * hsi0_587[k]
                   - f_10 * pc_x[k] * hsi1_587[k];

        t_588[k] = f_1 * isg0_315[k]
                   - f_2 * isg1_315[k]
                   + f_3 * pc_x[k] * ish_441[k];

        t_589[k] = f_16 * isg0_316[k]
                   - f_17 * isg1_316[k]
                   + f_3 * pc_x[k] * ish_442[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pc_x, pc_z, isg0_318, isg0_320, isg1_318, \
                         isg1_320, ish_441, ish_442, ish_444, ish_446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = f_3 * pc_z[k] * ish_441[k];

        t_591[k] = f_8 * isg0_318[k]
                   - f_9 * isg1_318[k]
                   + f_3 * pc_x[k] * ish_444[k];

        t_592[k] = f_3 * pc_z[k] * ish_442[k];

        t_593[k] = f_8 * isg0_320[k]
                   - f_9 * isg1_320[k]
                   + f_3 * pc_x[k] * ish_446[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pc_x, pc_z, isg0_321, isg0_323, isg0_324, \
                         isg1_321, isg1_323, isg1_324, ish_444, ish_447, ish_449, \
                         ish_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_6 * isg0_321[k]
                   - f_7 * isg1_321[k]
                   + f_3 * pc_x[k] * ish_447[k];

        t_595[k] = f_3 * pc_z[k] * ish_444[k];

        t_596[k] = f_6 * isg0_323[k]
                   - f_7 * isg1_323[k]
                   + f_3 * pc_x[k] * ish_449[k];

        t_597[k] = f_6 * isg0_324[k]
                   - f_7 * isg1_324[k]
                   + f_3 * pc_x[k] * ish_450[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pc_x, pc_z, isg0_325, isg0_327, isg0_328, \
                         isg1_325, isg1_327, isg1_328, ish_447, ish_451, ish_453, \
                         ish_454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_4 * isg0_325[k]
                   - f_5 * isg1_325[k]
                   + f_3 * pc_x[k] * ish_451[k];

        t_599[k] = f_3 * pc_z[k] * ish_447[k];

        t_600[k] = f_4 * isg0_327[k]
                   - f_5 * isg1_327[k]
                   + f_3 * pc_x[k] * ish_453[k];

        t_601[k] = f_4 * isg0_328[k]
                   - f_5 * isg1_328[k]
                   + f_3 * pc_x[k] * ish_454[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, t_606, t_607, pc_x, isg0_329, isg1_329, \
                         ish_455, ish_456, ish_457, ish_458, ish_459, \
                         ish_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_4 * isg0_329[k]
                   - f_5 * isg1_329[k]
                   + f_3 * pc_x[k] * ish_455[k];

        t_603[k] = f_3 * pc_x[k] * ish_456[k];

        t_604[k] = f_3 * pc_x[k] * ish_457[k];

        t_605[k] = f_3 * pc_x[k] * ish_458[k];

        t_606[k] = f_3 * pc_x[k] * ish_459[k];

        t_607[k] = f_3 * pc_x[k] * ish_460[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsi0_420 = buffer.data(hsi0 + 420);
    const auto *hsi0_421 = buffer.data(hsi0 + 421);
    const auto *hsi0_423 = buffer.data(hsi0 + 423);
    const auto *hsi0_426 = buffer.data(hsi0 + 426);
    const auto *hsi0_430 = buffer.data(hsi0 + 430);
    const auto *hsi0_441 = buffer.data(hsi0 + 441);
    const auto *hsi0_443 = buffer.data(hsi0 + 443);
    const auto *hsi0_444 = buffer.data(hsi0 + 444);
    const auto *hsi0_445 = buffer.data(hsi0 + 445);

    const auto *hsh_330 = buffer.data(hsh + 330);
    const auto *hsh_331 = buffer.data(hsh + 331);
    const auto *hsh_332 = buffer.data(hsh + 332);
    const auto *hsh_333 = buffer.data(hsh + 333);
    const auto *hsh_335 = buffer.data(hsh + 335);
    const auto *hsh_351 = buffer.data(hsh + 351);
    const auto *hsh_356 = buffer.data(hsh + 356);
    const auto *hsh_372 = buffer.data(hsh + 372);
    const auto *hsh_374 = buffer.data(hsh + 374);
    const auto *hsh_375 = buffer.data(hsh + 375);
    const auto *hsh_376 = buffer.data(hsh + 376);
    const auto *hsh_377 = buffer.data(hsh + 377);
    const auto *hsh_393 = buffer.data(hsh + 393);
    const auto *hsh_395 = buffer.data(hsh + 395);
    const auto *hsh_396 = buffer.data(hsh + 396);
    const auto *hsh_397 = buffer.data(hsh + 397);
    const auto *hsh_398 = buffer.data(hsh + 398);
    const auto *hsh_414 = buffer.data(hsh + 414);
    const auto *hsh_416 = buffer.data(hsh + 416);
    const auto *hsh_417 = buffer.data(hsh + 417);
    const auto *hsh_418 = buffer.data(hsh + 418);

    const auto *hsi1_420 = buffer.data(hsi1 + 420);
    const auto *hsi1_421 = buffer.data(hsi1 + 421);
    const auto *hsi1_423 = buffer.data(hsi1 + 423);
    const auto *hsi1_426 = buffer.data(hsi1 + 426);
    const auto *hsi1_430 = buffer.data(hsi1 + 430);
    const auto *hsi1_441 = buffer.data(hsi1 + 441);
    const auto *hsi1_443 = buffer.data(hsi1 + 443);
    const auto *hsi1_444 = buffer.data(hsi1 + 444);
    const auto *hsi1_445 = buffer.data(hsi1 + 445);

    const auto *isg0_325 = buffer.data(isg0 + 325);
    const auto *isg0_326 = buffer.data(isg0 + 326);
    const auto *isg0_327 = buffer.data(isg0 + 327);
    const auto *isg0_329 = buffer.data(isg0 + 329);
    const auto *isg0_332 = buffer.data(isg0 + 332);
    const auto *isg0_334 = buffer.data(isg0 + 334);
    const auto *isg0_335 = buffer.data(isg0 + 335);
    const auto *isg0_337 = buffer.data(isg0 + 337);
    const auto *isg0_338 = buffer.data(isg0 + 338);
    const auto *isg0_339 = buffer.data(isg0 + 339);
    const auto *isg0_341 = buffer.data(isg0 + 341);
    const auto *isg0_342 = buffer.data(isg0 + 342);
    const auto *isg0_343 = buffer.data(isg0 + 343);
    const auto *isg0_344 = buffer.data(isg0 + 344);
    const auto *isg0_345 = buffer.data(isg0 + 345);
    const auto *isg0_346 = buffer.data(isg0 + 346);
    const auto *isg0_347 = buffer.data(isg0 + 347);
    const auto *isg0_348 = buffer.data(isg0 + 348);
    const auto *isg0_349 = buffer.data(isg0 + 349);
    const auto *isg0_350 = buffer.data(isg0 + 350);
    const auto *isg0_351 = buffer.data(isg0 + 351);
    const auto *isg0_352 = buffer.data(isg0 + 352);
    const auto *isg0_353 = buffer.data(isg0 + 353);
    const auto *isg0_354 = buffer.data(isg0 + 354);
    const auto *isg0_355 = buffer.data(isg0 + 355);
    const auto *isg0_356 = buffer.data(isg0 + 356);
    const auto *isg0_357 = buffer.data(isg0 + 357);
    const auto *isg0_358 = buffer.data(isg0 + 358);
    const auto *isg0_359 = buffer.data(isg0 + 359);
    const auto *isg0_360 = buffer.data(isg0 + 360);
    const auto *isg0_361 = buffer.data(isg0 + 361);
    const auto *isg0_362 = buffer.data(isg0 + 362);
    const auto *isg0_363 = buffer.data(isg0 + 363);
    const auto *isg0_364 = buffer.data(isg0 + 364);
    const auto *isg0_365 = buffer.data(isg0 + 365);
    const auto *isg0_366 = buffer.data(isg0 + 366);
    const auto *isg0_367 = buffer.data(isg0 + 367);
    const auto *isg0_368 = buffer.data(isg0 + 368);
    const auto *isg0_369 = buffer.data(isg0 + 369);
    const auto *isg0_370 = buffer.data(isg0 + 370);
    const auto *isg0_371 = buffer.data(isg0 + 371);
    const auto *isg0_372 = buffer.data(isg0 + 372);
    const auto *isg0_373 = buffer.data(isg0 + 373);
    const auto *isg0_374 = buffer.data(isg0 + 374);
    const auto *isg0_375 = buffer.data(isg0 + 375);
    const auto *isg0_376 = buffer.data(isg0 + 376);
    const auto *isg0_377 = buffer.data(isg0 + 377);
    const auto *isg0_378 = buffer.data(isg0 + 378);
    const auto *isg0_379 = buffer.data(isg0 + 379);
    const auto *isg0_380 = buffer.data(isg0 + 380);
    const auto *isg0_381 = buffer.data(isg0 + 381);
    const auto *isg0_382 = buffer.data(isg0 + 382);
    const auto *isg0_383 = buffer.data(isg0 + 383);
    const auto *isg0_384 = buffer.data(isg0 + 384);
    const auto *isg0_385 = buffer.data(isg0 + 385);
    const auto *isg0_386 = buffer.data(isg0 + 386);
    const auto *isg0_387 = buffer.data(isg0 + 387);
    const auto *isg0_388 = buffer.data(isg0 + 388);
    const auto *isg0_389 = buffer.data(isg0 + 389);

    const auto *isg1_325 = buffer.data(isg1 + 325);
    const auto *isg1_326 = buffer.data(isg1 + 326);
    const auto *isg1_327 = buffer.data(isg1 + 327);
    const auto *isg1_329 = buffer.data(isg1 + 329);
    const auto *isg1_332 = buffer.data(isg1 + 332);
    const auto *isg1_334 = buffer.data(isg1 + 334);
    const auto *isg1_335 = buffer.data(isg1 + 335);
    const auto *isg1_337 = buffer.data(isg1 + 337);
    const auto *isg1_338 = buffer.data(isg1 + 338);
    const auto *isg1_339 = buffer.data(isg1 + 339);
    const auto *isg1_341 = buffer.data(isg1 + 341);
    const auto *isg1_342 = buffer.data(isg1 + 342);
    const auto *isg1_343 = buffer.data(isg1 + 343);
    const auto *isg1_344 = buffer.data(isg1 + 344);
    const auto *isg1_345 = buffer.data(isg1 + 345);
    const auto *isg1_346 = buffer.data(isg1 + 346);
    const auto *isg1_347 = buffer.data(isg1 + 347);
    const auto *isg1_348 = buffer.data(isg1 + 348);
    const auto *isg1_349 = buffer.data(isg1 + 349);
    const auto *isg1_350 = buffer.data(isg1 + 350);
    const auto *isg1_351 = buffer.data(isg1 + 351);
    const auto *isg1_352 = buffer.data(isg1 + 352);
    const auto *isg1_353 = buffer.data(isg1 + 353);
    const auto *isg1_354 = buffer.data(isg1 + 354);
    const auto *isg1_355 = buffer.data(isg1 + 355);
    const auto *isg1_356 = buffer.data(isg1 + 356);
    const auto *isg1_357 = buffer.data(isg1 + 357);
    const auto *isg1_358 = buffer.data(isg1 + 358);
    const auto *isg1_359 = buffer.data(isg1 + 359);
    const auto *isg1_360 = buffer.data(isg1 + 360);
    const auto *isg1_361 = buffer.data(isg1 + 361);
    const auto *isg1_362 = buffer.data(isg1 + 362);
    const auto *isg1_363 = buffer.data(isg1 + 363);
    const auto *isg1_364 = buffer.data(isg1 + 364);
    const auto *isg1_365 = buffer.data(isg1 + 365);
    const auto *isg1_366 = buffer.data(isg1 + 366);
    const auto *isg1_367 = buffer.data(isg1 + 367);
    const auto *isg1_368 = buffer.data(isg1 + 368);
    const auto *isg1_369 = buffer.data(isg1 + 369);
    const auto *isg1_370 = buffer.data(isg1 + 370);
    const auto *isg1_371 = buffer.data(isg1 + 371);
    const auto *isg1_372 = buffer.data(isg1 + 372);
    const auto *isg1_373 = buffer.data(isg1 + 373);
    const auto *isg1_374 = buffer.data(isg1 + 374);
    const auto *isg1_375 = buffer.data(isg1 + 375);
    const auto *isg1_376 = buffer.data(isg1 + 376);
    const auto *isg1_377 = buffer.data(isg1 + 377);
    const auto *isg1_378 = buffer.data(isg1 + 378);
    const auto *isg1_379 = buffer.data(isg1 + 379);
    const auto *isg1_380 = buffer.data(isg1 + 380);
    const auto *isg1_381 = buffer.data(isg1 + 381);
    const auto *isg1_382 = buffer.data(isg1 + 382);
    const auto *isg1_383 = buffer.data(isg1 + 383);
    const auto *isg1_384 = buffer.data(isg1 + 384);
    const auto *isg1_385 = buffer.data(isg1 + 385);
    const auto *isg1_386 = buffer.data(isg1 + 386);
    const auto *isg1_387 = buffer.data(isg1 + 387);
    const auto *isg1_388 = buffer.data(isg1 + 388);
    const auto *isg1_389 = buffer.data(isg1 + 389);

    const auto *ish_456 = buffer.data(ish + 456);
    const auto *ish_457 = buffer.data(ish + 457);
    const auto *ish_458 = buffer.data(ish + 458);
    const auto *ish_459 = buffer.data(ish + 459);
    const auto *ish_461 = buffer.data(ish + 461);
    const auto *ish_464 = buffer.data(ish + 464);
    const auto *ish_466 = buffer.data(ish + 466);
    const auto *ish_467 = buffer.data(ish + 467);
    const auto *ish_469 = buffer.data(ish + 469);
    const auto *ish_470 = buffer.data(ish + 470);
    const auto *ish_471 = buffer.data(ish + 471);
    const auto *ish_473 = buffer.data(ish + 473);
    const auto *ish_474 = buffer.data(ish + 474);
    const auto *ish_475 = buffer.data(ish + 475);
    const auto *ish_476 = buffer.data(ish + 476);
    const auto *ish_477 = buffer.data(ish + 477);
    const auto *ish_478 = buffer.data(ish + 478);
    const auto *ish_479 = buffer.data(ish + 479);
    const auto *ish_480 = buffer.data(ish + 480);
    const auto *ish_481 = buffer.data(ish + 481);
    const auto *ish_482 = buffer.data(ish + 482);
    const auto *ish_483 = buffer.data(ish + 483);
    const auto *ish_484 = buffer.data(ish + 484);
    const auto *ish_485 = buffer.data(ish + 485);
    const auto *ish_486 = buffer.data(ish + 486);
    const auto *ish_487 = buffer.data(ish + 487);
    const auto *ish_488 = buffer.data(ish + 488);
    const auto *ish_489 = buffer.data(ish + 489);
    const auto *ish_490 = buffer.data(ish + 490);
    const auto *ish_491 = buffer.data(ish + 491);
    const auto *ish_492 = buffer.data(ish + 492);
    const auto *ish_493 = buffer.data(ish + 493);
    const auto *ish_494 = buffer.data(ish + 494);
    const auto *ish_495 = buffer.data(ish + 495);
    const auto *ish_496 = buffer.data(ish + 496);
    const auto *ish_497 = buffer.data(ish + 497);
    const auto *ish_498 = buffer.data(ish + 498);
    const auto *ish_499 = buffer.data(ish + 499);
    const auto *ish_500 = buffer.data(ish + 500);
    const auto *ish_501 = buffer.data(ish + 501);
    const auto *ish_502 = buffer.data(ish + 502);
    const auto *ish_503 = buffer.data(ish + 503);
    const auto *ish_504 = buffer.data(ish + 504);
    const auto *ish_505 = buffer.data(ish + 505);
    const auto *ish_506 = buffer.data(ish + 506);
    const auto *ish_507 = buffer.data(ish + 507);
    const auto *ish_508 = buffer.data(ish + 508);
    const auto *ish_509 = buffer.data(ish + 509);
    const auto *ish_510 = buffer.data(ish + 510);
    const auto *ish_511 = buffer.data(ish + 511);
    const auto *ish_512 = buffer.data(ish + 512);
    const auto *ish_513 = buffer.data(ish + 513);
    const auto *ish_514 = buffer.data(ish + 514);
    const auto *ish_515 = buffer.data(ish + 515);
    const auto *ish_516 = buffer.data(ish + 516);
    const auto *ish_517 = buffer.data(ish + 517);
    const auto *ish_518 = buffer.data(ish + 518);
    const auto *ish_519 = buffer.data(ish + 519);
    const auto *ish_520 = buffer.data(ish + 520);
    const auto *ish_521 = buffer.data(ish + 521);
    const auto *ish_522 = buffer.data(ish + 522);
    const auto *ish_523 = buffer.data(ish + 523);
    const auto *ish_524 = buffer.data(ish + 524);
    const auto *ish_525 = buffer.data(ish + 525);
    const auto *ish_526 = buffer.data(ish + 526);
    const auto *ish_527 = buffer.data(ish + 527);
    const auto *ish_528 = buffer.data(ish + 528);
    const auto *ish_529 = buffer.data(ish + 529);
    const auto *ish_530 = buffer.data(ish + 530);
    const auto *ish_531 = buffer.data(ish + 531);
    const auto *ish_532 = buffer.data(ish + 532);
    const auto *ish_533 = buffer.data(ish + 533);
    const auto *ish_534 = buffer.data(ish + 534);
    const auto *ish_535 = buffer.data(ish + 535);
    const auto *ish_536 = buffer.data(ish + 536);
    const auto *ish_537 = buffer.data(ish + 537);
    const auto *ish_538 = buffer.data(ish + 538);
    const auto *ish_539 = buffer.data(ish + 539);
    const auto *ish_540 = buffer.data(ish + 540);
    const auto *ish_541 = buffer.data(ish + 541);
    const auto *ish_542 = buffer.data(ish + 542);
    const auto *ish_543 = buffer.data(ish + 543);
    const auto *ish_544 = buffer.data(ish + 544);
    const auto *ish_545 = buffer.data(ish + 545);

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pc_x, pc_y, pc_z, hsh_330, isg0_325, \
                         isg1_325, ish_456, ish_457, ish_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_3 * pc_x[k] * ish_461[k];

        t_609[k] = f_0 * hsh_330[k]
                   + f_1 * isg0_325[k]
                   - f_2 * isg1_325[k]
                   + f_3 * pc_y[k] * ish_456[k];

        t_610[k] = f_3 * pc_z[k] * ish_456[k];

        t_611[k] = f_4 * isg0_325[k]
                   - f_5 * isg1_325[k]
                   + f_3 * pc_z[k] * ish_457[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pc_y, pc_z, hsh_335, isg0_326, isg0_327, \
                         isg0_329, isg1_326, isg1_327, isg1_329, ish_458, ish_459, \
                         ish_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_6 * isg0_326[k]
                   - f_7 * isg1_326[k]
                   + f_3 * pc_z[k] * ish_458[k];

        t_613[k] = f_8 * isg0_327[k]
                   - f_9 * isg1_327[k]
                   + f_3 * pc_z[k] * ish_459[k];

        t_614[k] = f_0 * hsh_335[k]
                   + f_3 * pc_y[k] * ish_461[k];

        t_615[k] = f_1 * isg0_329[k]
                   - f_2 * isg1_329[k]
                   + f_3 * pc_z[k] * ish_461[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pa_z, pc_x, pc_z, hsi0_420, hsi0_421, \
                         hsi0_423, hsi1_420, hsi1_421, hsi1_423, isg0_332, isg1_332, \
                         ish_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pa_z[k] * hsi0_420[k]
                   - f_10 * pc_z[k] * hsi1_420[k];

        t_617[k] = pa_z[k] * hsi0_421[k]
                   - f_10 * pc_z[k] * hsi1_421[k];

        t_618[k] = f_16 * isg0_332[k]
                   - f_17 * isg1_332[k]
                   + f_3 * pc_x[k] * ish_464[k];

        t_619[k] = pa_z[k] * hsi0_423[k]
                   - f_10 * pc_z[k] * hsi1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_z, pc_x, pc_z, hsi0_426, hsi1_426, isg0_334, \
                         isg0_335, isg1_334, isg1_335, ish_466, \
                         ish_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_8 * isg0_334[k]
                   - f_9 * isg1_334[k]
                   + f_3 * pc_x[k] * ish_466[k];

        t_621[k] = f_8 * isg0_335[k]
                   - f_9 * isg1_335[k]
                   + f_3 * pc_x[k] * ish_467[k];

        t_622[k] = pa_z[k] * hsi0_426[k]
                   - f_10 * pc_z[k] * hsi1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, isg0_337, isg0_338, isg0_339, isg1_337, \
                         isg1_338, isg1_339, ish_469, ish_470, \
                         ish_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_6 * isg0_337[k]
                   - f_7 * isg1_337[k]
                   + f_3 * pc_x[k] * ish_469[k];

        t_624[k] = f_6 * isg0_338[k]
                   - f_7 * isg1_338[k]
                   + f_3 * pc_x[k] * ish_470[k];

        t_625[k] = f_6 * isg0_339[k]
                   - f_7 * isg1_339[k]
                   + f_3 * pc_x[k] * ish_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_z, pc_x, pc_z, hsi0_430, hsi1_430, isg0_341, \
                         isg0_342, isg1_341, isg1_342, ish_473, \
                         ish_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_z[k] * hsi0_430[k]
                   - f_10 * pc_z[k] * hsi1_430[k];

        t_627[k] = f_4 * isg0_341[k]
                   - f_5 * isg1_341[k]
                   + f_3 * pc_x[k] * ish_473[k];

        t_628[k] = f_4 * isg0_342[k]
                   - f_5 * isg1_342[k]
                   + f_3 * pc_x[k] * ish_474[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, pc_x, isg0_343, isg0_344, \
                         isg1_343, isg1_344, ish_475, ish_476, ish_477, ish_478, \
                         ish_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_4 * isg0_343[k]
                   - f_5 * isg1_343[k]
                   + f_3 * pc_x[k] * ish_475[k];

        t_630[k] = f_4 * isg0_344[k]
                   - f_5 * isg1_344[k]
                   + f_3 * pc_x[k] * ish_476[k];

        t_631[k] = f_3 * pc_x[k] * ish_477[k];

        t_632[k] = f_3 * pc_x[k] * ish_478[k];

        t_633[k] = f_3 * pc_x[k] * ish_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pa_z, pc_x, pc_z, hsi0_441, \
                         hsh_330, hsi1_441, ish_477, ish_480, ish_481, \
                         ish_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_3 * pc_x[k] * ish_480[k];

        t_635[k] = f_3 * pc_x[k] * ish_481[k];

        t_636[k] = f_3 * pc_x[k] * ish_482[k];

        t_637[k] = pa_z[k] * hsi0_441[k]
                   - f_10 * pc_z[k] * hsi1_441[k];

        t_638[k] = f_11 * hsh_330[k]
                   + f_3 * pc_z[k] * ish_477[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_z, pc_z, hsi0_443, hsi0_444, hsi0_445, \
                         hsh_331, hsh_332, hsh_333, hsi1_443, hsi1_444, \
                         hsi1_445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pa_z[k] * hsi0_443[k]
                   + f_12 * hsh_331[k]
                   - f_10 * pc_z[k] * hsi1_443[k];

        t_640[k] = pa_z[k] * hsi0_444[k]
                   + f_13 * hsh_332[k]
                   - f_10 * pc_z[k] * hsi1_444[k];

        t_641[k] = pa_z[k] * hsi0_445[k]
                   + f_14 * hsh_333[k]
                   - f_10 * pc_z[k] * hsi1_445[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, pc_z, hsh_335, hsh_356, isg0_344, \
                         isg0_345, isg1_344, isg1_345, ish_482, \
                         ish_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_15 * hsh_356[k]
                   + f_3 * pc_y[k] * ish_482[k];

        t_643[k] = f_11 * hsh_335[k]
                   + f_1 * isg0_344[k]
                   - f_2 * isg1_344[k]
                   + f_3 * pc_z[k] * ish_482[k];

        t_644[k] = f_1 * isg0_345[k]
                   - f_2 * isg1_345[k]
                   + f_3 * pc_x[k] * ish_483[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, isg0_346, isg0_347, isg0_348, isg1_346, \
                         isg1_347, isg1_348, ish_484, ish_485, \
                         ish_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_16 * isg0_346[k]
                   - f_17 * isg1_346[k]
                   + f_3 * pc_x[k] * ish_484[k];

        t_646[k] = f_16 * isg0_347[k]
                   - f_17 * isg1_347[k]
                   + f_3 * pc_x[k] * ish_485[k];

        t_647[k] = f_8 * isg0_348[k]
                   - f_9 * isg1_348[k]
                   + f_3 * pc_x[k] * ish_486[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, isg0_349, isg0_350, isg0_351, isg1_349, \
                         isg1_350, isg1_351, ish_487, ish_488, \
                         ish_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_8 * isg0_349[k]
                   - f_9 * isg1_349[k]
                   + f_3 * pc_x[k] * ish_487[k];

        t_649[k] = f_8 * isg0_350[k]
                   - f_9 * isg1_350[k]
                   + f_3 * pc_x[k] * ish_488[k];

        t_650[k] = f_6 * isg0_351[k]
                   - f_7 * isg1_351[k]
                   + f_3 * pc_x[k] * ish_489[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, isg0_352, isg0_353, isg0_354, isg1_352, \
                         isg1_353, isg1_354, ish_490, ish_491, \
                         ish_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_6 * isg0_352[k]
                   - f_7 * isg1_352[k]
                   + f_3 * pc_x[k] * ish_490[k];

        t_652[k] = f_6 * isg0_353[k]
                   - f_7 * isg1_353[k]
                   + f_3 * pc_x[k] * ish_491[k];

        t_653[k] = f_6 * isg0_354[k]
                   - f_7 * isg1_354[k]
                   + f_3 * pc_x[k] * ish_492[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, isg0_355, isg0_356, isg0_357, isg1_355, \
                         isg1_356, isg1_357, ish_493, ish_494, \
                         ish_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_4 * isg0_355[k]
                   - f_5 * isg1_355[k]
                   + f_3 * pc_x[k] * ish_493[k];

        t_655[k] = f_4 * isg0_356[k]
                   - f_5 * isg1_356[k]
                   + f_3 * pc_x[k] * ish_494[k];

        t_656[k] = f_4 * isg0_357[k]
                   - f_5 * isg1_357[k]
                   + f_3 * pc_x[k] * ish_495[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, pc_x, isg0_358, isg0_359, \
                         isg1_358, isg1_359, ish_496, ish_497, ish_498, ish_499, \
                         ish_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_4 * isg0_358[k]
                   - f_5 * isg1_358[k]
                   + f_3 * pc_x[k] * ish_496[k];

        t_658[k] = f_4 * isg0_359[k]
                   - f_5 * isg1_359[k]
                   + f_3 * pc_x[k] * ish_497[k];

        t_659[k] = f_3 * pc_x[k] * ish_498[k];

        t_660[k] = f_3 * pc_x[k] * ish_499[k];

        t_661[k] = f_3 * pc_x[k] * ish_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, t_666, pc_x, pc_y, pc_z, hsh_351, \
                         hsh_372, isg0_355, isg1_355, ish_498, ish_501, ish_502, \
                         ish_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_3 * pc_x[k] * ish_501[k];

        t_663[k] = f_3 * pc_x[k] * ish_502[k];

        t_664[k] = f_3 * pc_x[k] * ish_503[k];

        t_665[k] = f_14 * hsh_372[k]
                   + f_1 * isg0_355[k]
                   - f_2 * isg1_355[k]
                   + f_3 * pc_y[k] * ish_498[k];

        t_666[k] = f_12 * hsh_351[k]
                   + f_3 * pc_z[k] * ish_498[k];
    }

#pragma omp simd aligned(t_667, t_668, t_669, pc_y, hsh_374, hsh_375, hsh_376, isg0_357, \
                         isg0_358, isg0_359, isg1_357, isg1_358, isg1_359, ish_500, ish_501, \
                         ish_502 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_667[k] = f_14 * hsh_374[k]
                   + f_8 * isg0_357[k]
                   - f_9 * isg1_357[k]
                   + f_3 * pc_y[k] * ish_500[k];

        t_668[k] = f_14 * hsh_375[k]
                   + f_6 * isg0_358[k]
                   - f_7 * isg1_358[k]
                   + f_3 * pc_y[k] * ish_501[k];

        t_669[k] = f_14 * hsh_376[k]
                   + f_4 * isg0_359[k]
                   - f_5 * isg1_359[k]
                   + f_3 * pc_y[k] * ish_502[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, pc_x, pc_y, pc_z, hsh_356, hsh_377, isg0_359, \
                         isg0_360, isg1_359, isg1_360, ish_503, \
                         ish_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * hsh_377[k]
                   + f_3 * pc_y[k] * ish_503[k];

        t_671[k] = f_12 * hsh_356[k]
                   + f_1 * isg0_359[k]
                   - f_2 * isg1_359[k]
                   + f_3 * pc_z[k] * ish_503[k];

        t_672[k] = f_1 * isg0_360[k]
                   - f_2 * isg1_360[k]
                   + f_3 * pc_x[k] * ish_504[k];
    }

#pragma omp simd aligned(t_673, t_674, t_675, pc_x, isg0_361, isg0_362, isg0_363, isg1_361, \
                         isg1_362, isg1_363, ish_505, ish_506, \
                         ish_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_673[k] = f_16 * isg0_361[k]
                   - f_17 * isg1_361[k]
                   + f_3 * pc_x[k] * ish_505[k];

        t_674[k] = f_16 * isg0_362[k]
                   - f_17 * isg1_362[k]
                   + f_3 * pc_x[k] * ish_506[k];

        t_675[k] = f_8 * isg0_363[k]
                   - f_9 * isg1_363[k]
                   + f_3 * pc_x[k] * ish_507[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pc_x, isg0_364, isg0_365, isg0_366, isg1_364, \
                         isg1_365, isg1_366, ish_508, ish_509, \
                         ish_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_8 * isg0_364[k]
                   - f_9 * isg1_364[k]
                   + f_3 * pc_x[k] * ish_508[k];

        t_677[k] = f_8 * isg0_365[k]
                   - f_9 * isg1_365[k]
                   + f_3 * pc_x[k] * ish_509[k];

        t_678[k] = f_6 * isg0_366[k]
                   - f_7 * isg1_366[k]
                   + f_3 * pc_x[k] * ish_510[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pc_x, isg0_367, isg0_368, isg0_369, isg1_367, \
                         isg1_368, isg1_369, ish_511, ish_512, \
                         ish_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = f_6 * isg0_367[k]
                   - f_7 * isg1_367[k]
                   + f_3 * pc_x[k] * ish_511[k];

        t_680[k] = f_6 * isg0_368[k]
                   - f_7 * isg1_368[k]
                   + f_3 * pc_x[k] * ish_512[k];

        t_681[k] = f_6 * isg0_369[k]
                   - f_7 * isg1_369[k]
                   + f_3 * pc_x[k] * ish_513[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, pc_x, isg0_370, isg0_371, isg0_372, isg1_370, \
                         isg1_371, isg1_372, ish_514, ish_515, \
                         ish_516 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = f_4 * isg0_370[k]
                   - f_5 * isg1_370[k]
                   + f_3 * pc_x[k] * ish_514[k];

        t_683[k] = f_4 * isg0_371[k]
                   - f_5 * isg1_371[k]
                   + f_3 * pc_x[k] * ish_515[k];

        t_684[k] = f_4 * isg0_372[k]
                   - f_5 * isg1_372[k]
                   + f_3 * pc_x[k] * ish_516[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, t_689, pc_x, isg0_373, isg0_374, \
                         isg1_373, isg1_374, ish_517, ish_518, ish_519, ish_520, \
                         ish_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = f_4 * isg0_373[k]
                   - f_5 * isg1_373[k]
                   + f_3 * pc_x[k] * ish_517[k];

        t_686[k] = f_4 * isg0_374[k]
                   - f_5 * isg1_374[k]
                   + f_3 * pc_x[k] * ish_518[k];

        t_687[k] = f_3 * pc_x[k] * ish_519[k];

        t_688[k] = f_3 * pc_x[k] * ish_520[k];

        t_689[k] = f_3 * pc_x[k] * ish_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, t_694, pc_x, pc_y, pc_z, hsh_372, \
                         hsh_393, isg0_370, isg1_370, ish_519, ish_522, ish_523, \
                         ish_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_3 * pc_x[k] * ish_522[k];

        t_691[k] = f_3 * pc_x[k] * ish_523[k];

        t_692[k] = f_3 * pc_x[k] * ish_524[k];

        t_693[k] = f_13 * hsh_393[k]
                   + f_1 * isg0_370[k]
                   - f_2 * isg1_370[k]
                   + f_3 * pc_y[k] * ish_519[k];

        t_694[k] = f_13 * hsh_372[k]
                   + f_3 * pc_z[k] * ish_519[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pc_y, hsh_395, hsh_396, hsh_397, isg0_372, \
                         isg0_373, isg0_374, isg1_372, isg1_373, isg1_374, ish_521, ish_522, \
                         ish_523 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = f_13 * hsh_395[k]
                   + f_8 * isg0_372[k]
                   - f_9 * isg1_372[k]
                   + f_3 * pc_y[k] * ish_521[k];

        t_696[k] = f_13 * hsh_396[k]
                   + f_6 * isg0_373[k]
                   - f_7 * isg1_373[k]
                   + f_3 * pc_y[k] * ish_522[k];

        t_697[k] = f_13 * hsh_397[k]
                   + f_4 * isg0_374[k]
                   - f_5 * isg1_374[k]
                   + f_3 * pc_y[k] * ish_523[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pc_x, pc_y, pc_z, hsh_377, hsh_398, isg0_374, \
                         isg0_375, isg1_374, isg1_375, ish_524, \
                         ish_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_13 * hsh_398[k]
                   + f_3 * pc_y[k] * ish_524[k];

        t_699[k] = f_13 * hsh_377[k]
                   + f_1 * isg0_374[k]
                   - f_2 * isg1_374[k]
                   + f_3 * pc_z[k] * ish_524[k];

        t_700[k] = f_1 * isg0_375[k]
                   - f_2 * isg1_375[k]
                   + f_3 * pc_x[k] * ish_525[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, pc_x, isg0_376, isg0_377, isg0_378, isg1_376, \
                         isg1_377, isg1_378, ish_526, ish_527, \
                         ish_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_16 * isg0_376[k]
                   - f_17 * isg1_376[k]
                   + f_3 * pc_x[k] * ish_526[k];

        t_702[k] = f_16 * isg0_377[k]
                   - f_17 * isg1_377[k]
                   + f_3 * pc_x[k] * ish_527[k];

        t_703[k] = f_8 * isg0_378[k]
                   - f_9 * isg1_378[k]
                   + f_3 * pc_x[k] * ish_528[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, pc_x, isg0_379, isg0_380, isg0_381, isg1_379, \
                         isg1_380, isg1_381, ish_529, ish_530, \
                         ish_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_8 * isg0_379[k]
                   - f_9 * isg1_379[k]
                   + f_3 * pc_x[k] * ish_529[k];

        t_705[k] = f_8 * isg0_380[k]
                   - f_9 * isg1_380[k]
                   + f_3 * pc_x[k] * ish_530[k];

        t_706[k] = f_6 * isg0_381[k]
                   - f_7 * isg1_381[k]
                   + f_3 * pc_x[k] * ish_531[k];
    }

#pragma omp simd aligned(t_707, t_708, t_709, pc_x, isg0_382, isg0_383, isg0_384, isg1_382, \
                         isg1_383, isg1_384, ish_532, ish_533, \
                         ish_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_707[k] = f_6 * isg0_382[k]
                   - f_7 * isg1_382[k]
                   + f_3 * pc_x[k] * ish_532[k];

        t_708[k] = f_6 * isg0_383[k]
                   - f_7 * isg1_383[k]
                   + f_3 * pc_x[k] * ish_533[k];

        t_709[k] = f_6 * isg0_384[k]
                   - f_7 * isg1_384[k]
                   + f_3 * pc_x[k] * ish_534[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, pc_x, isg0_385, isg0_386, isg0_387, isg1_385, \
                         isg1_386, isg1_387, ish_535, ish_536, \
                         ish_537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_4 * isg0_385[k]
                   - f_5 * isg1_385[k]
                   + f_3 * pc_x[k] * ish_535[k];

        t_711[k] = f_4 * isg0_386[k]
                   - f_5 * isg1_386[k]
                   + f_3 * pc_x[k] * ish_536[k];

        t_712[k] = f_4 * isg0_387[k]
                   - f_5 * isg1_387[k]
                   + f_3 * pc_x[k] * ish_537[k];
    }

#pragma omp simd aligned(t_713, t_714, t_715, t_716, t_717, pc_x, isg0_388, isg0_389, \
                         isg1_388, isg1_389, ish_538, ish_539, ish_540, ish_541, \
                         ish_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_713[k] = f_4 * isg0_388[k]
                   - f_5 * isg1_388[k]
                   + f_3 * pc_x[k] * ish_538[k];

        t_714[k] = f_4 * isg0_389[k]
                   - f_5 * isg1_389[k]
                   + f_3 * pc_x[k] * ish_539[k];

        t_715[k] = f_3 * pc_x[k] * ish_540[k];

        t_716[k] = f_3 * pc_x[k] * ish_541[k];

        t_717[k] = f_3 * pc_x[k] * ish_542[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, t_722, pc_x, pc_y, pc_z, hsh_393, \
                         hsh_414, isg0_385, isg1_385, ish_540, ish_543, ish_544, \
                         ish_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_3 * pc_x[k] * ish_543[k];

        t_719[k] = f_3 * pc_x[k] * ish_544[k];

        t_720[k] = f_3 * pc_x[k] * ish_545[k];

        t_721[k] = f_12 * hsh_414[k]
                   + f_1 * isg0_385[k]
                   - f_2 * isg1_385[k]
                   + f_3 * pc_y[k] * ish_540[k];

        t_722[k] = f_14 * hsh_393[k]
                   + f_3 * pc_z[k] * ish_540[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_y, hsh_416, hsh_417, hsh_418, isg0_387, \
                         isg0_388, isg0_389, isg1_387, isg1_388, isg1_389, ish_542, ish_543, \
                         ish_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_12 * hsh_416[k]
                   + f_8 * isg0_387[k]
                   - f_9 * isg1_387[k]
                   + f_3 * pc_y[k] * ish_542[k];

        t_724[k] = f_12 * hsh_417[k]
                   + f_6 * isg0_388[k]
                   - f_7 * isg1_388[k]
                   + f_3 * pc_y[k] * ish_543[k];

        t_725[k] = f_12 * hsh_418[k]
                   + f_4 * isg0_389[k]
                   - f_5 * isg1_389[k]
                   + f_3 * pc_y[k] * ish_544[k];
    }
}

static auto
compute_prim_isi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsi0,
                                                          const size_t hsh, const size_t hsi1,
                                                          const size_t isg0, const size_t isg1,
                                                          const size_t ish, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_15 = 2.5 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsi0_560 = buffer.data(hsi0 + 560);
    const auto *hsi0_562 = buffer.data(hsi0 + 562);
    const auto *hsi0_565 = buffer.data(hsi0 + 565);
    const auto *hsi0_569 = buffer.data(hsi0 + 569);
    const auto *hsi0_574 = buffer.data(hsi0 + 574);
    const auto *hsi0_581 = buffer.data(hsi0 + 581);
    const auto *hsi0_583 = buffer.data(hsi0 + 583);
    const auto *hsi0_584 = buffer.data(hsi0 + 584);
    const auto *hsi0_585 = buffer.data(hsi0 + 585);
    const auto *hsi0_587 = buffer.data(hsi0 + 587);

    const auto *hsh_398 = buffer.data(hsh + 398);
    const auto *hsh_414 = buffer.data(hsh + 414);
    const auto *hsh_419 = buffer.data(hsh + 419);
    const auto *hsh_435 = buffer.data(hsh + 435);
    const auto *hsh_437 = buffer.data(hsh + 437);
    const auto *hsh_438 = buffer.data(hsh + 438);
    const auto *hsh_439 = buffer.data(hsh + 439);
    const auto *hsh_440 = buffer.data(hsh + 440);

    const auto *hsi1_560 = buffer.data(hsi1 + 560);
    const auto *hsi1_562 = buffer.data(hsi1 + 562);
    const auto *hsi1_565 = buffer.data(hsi1 + 565);
    const auto *hsi1_569 = buffer.data(hsi1 + 569);
    const auto *hsi1_574 = buffer.data(hsi1 + 574);
    const auto *hsi1_581 = buffer.data(hsi1 + 581);
    const auto *hsi1_583 = buffer.data(hsi1 + 583);
    const auto *hsi1_584 = buffer.data(hsi1 + 584);
    const auto *hsi1_585 = buffer.data(hsi1 + 585);
    const auto *hsi1_587 = buffer.data(hsi1 + 587);

    const auto *isg0_389 = buffer.data(isg0 + 389);
    const auto *isg0_391 = buffer.data(isg0 + 391);
    const auto *isg0_393 = buffer.data(isg0 + 393);
    const auto *isg0_394 = buffer.data(isg0 + 394);
    const auto *isg0_396 = buffer.data(isg0 + 396);
    const auto *isg0_397 = buffer.data(isg0 + 397);
    const auto *isg0_398 = buffer.data(isg0 + 398);
    const auto *isg0_400 = buffer.data(isg0 + 400);
    const auto *isg0_401 = buffer.data(isg0 + 401);
    const auto *isg0_402 = buffer.data(isg0 + 402);
    const auto *isg0_403 = buffer.data(isg0 + 403);
    const auto *isg0_405 = buffer.data(isg0 + 405);
    const auto *isg0_407 = buffer.data(isg0 + 407);
    const auto *isg0_408 = buffer.data(isg0 + 408);
    const auto *isg0_410 = buffer.data(isg0 + 410);
    const auto *isg0_411 = buffer.data(isg0 + 411);
    const auto *isg0_412 = buffer.data(isg0 + 412);
    const auto *isg0_414 = buffer.data(isg0 + 414);
    const auto *isg0_415 = buffer.data(isg0 + 415);
    const auto *isg0_416 = buffer.data(isg0 + 416);
    const auto *isg0_417 = buffer.data(isg0 + 417);
    const auto *isg0_418 = buffer.data(isg0 + 418);
    const auto *isg0_419 = buffer.data(isg0 + 419);

    const auto *isg1_389 = buffer.data(isg1 + 389);
    const auto *isg1_391 = buffer.data(isg1 + 391);
    const auto *isg1_393 = buffer.data(isg1 + 393);
    const auto *isg1_394 = buffer.data(isg1 + 394);
    const auto *isg1_396 = buffer.data(isg1 + 396);
    const auto *isg1_397 = buffer.data(isg1 + 397);
    const auto *isg1_398 = buffer.data(isg1 + 398);
    const auto *isg1_400 = buffer.data(isg1 + 400);
    const auto *isg1_401 = buffer.data(isg1 + 401);
    const auto *isg1_402 = buffer.data(isg1 + 402);
    const auto *isg1_403 = buffer.data(isg1 + 403);
    const auto *isg1_405 = buffer.data(isg1 + 405);
    const auto *isg1_407 = buffer.data(isg1 + 407);
    const auto *isg1_408 = buffer.data(isg1 + 408);
    const auto *isg1_410 = buffer.data(isg1 + 410);
    const auto *isg1_411 = buffer.data(isg1 + 411);
    const auto *isg1_412 = buffer.data(isg1 + 412);
    const auto *isg1_414 = buffer.data(isg1 + 414);
    const auto *isg1_415 = buffer.data(isg1 + 415);
    const auto *isg1_416 = buffer.data(isg1 + 416);
    const auto *isg1_417 = buffer.data(isg1 + 417);
    const auto *isg1_418 = buffer.data(isg1 + 418);
    const auto *isg1_419 = buffer.data(isg1 + 419);

    const auto *ish_545 = buffer.data(ish + 545);
    const auto *ish_547 = buffer.data(ish + 547);
    const auto *ish_549 = buffer.data(ish + 549);
    const auto *ish_550 = buffer.data(ish + 550);
    const auto *ish_552 = buffer.data(ish + 552);
    const auto *ish_553 = buffer.data(ish + 553);
    const auto *ish_554 = buffer.data(ish + 554);
    const auto *ish_556 = buffer.data(ish + 556);
    const auto *ish_557 = buffer.data(ish + 557);
    const auto *ish_558 = buffer.data(ish + 558);
    const auto *ish_559 = buffer.data(ish + 559);
    const auto *ish_561 = buffer.data(ish + 561);
    const auto *ish_562 = buffer.data(ish + 562);
    const auto *ish_563 = buffer.data(ish + 563);
    const auto *ish_564 = buffer.data(ish + 564);
    const auto *ish_565 = buffer.data(ish + 565);
    const auto *ish_566 = buffer.data(ish + 566);
    const auto *ish_567 = buffer.data(ish + 567);
    const auto *ish_569 = buffer.data(ish + 569);
    const auto *ish_570 = buffer.data(ish + 570);
    const auto *ish_572 = buffer.data(ish + 572);
    const auto *ish_573 = buffer.data(ish + 573);
    const auto *ish_574 = buffer.data(ish + 574);
    const auto *ish_576 = buffer.data(ish + 576);
    const auto *ish_577 = buffer.data(ish + 577);
    const auto *ish_578 = buffer.data(ish + 578);
    const auto *ish_579 = buffer.data(ish + 579);
    const auto *ish_581 = buffer.data(ish + 581);
    const auto *ish_582 = buffer.data(ish + 582);
    const auto *ish_583 = buffer.data(ish + 583);
    const auto *ish_584 = buffer.data(ish + 584);
    const auto *ish_585 = buffer.data(ish + 585);
    const auto *ish_586 = buffer.data(ish + 586);
    const auto *ish_587 = buffer.data(ish + 587);

#pragma omp simd aligned(t_726, t_727, t_728, pa_y, pc_y, pc_z, hsi0_560, hsh_398, hsh_419, \
                         hsi1_560, isg0_389, isg1_389, ish_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * hsh_419[k]
                   + f_3 * pc_y[k] * ish_545[k];

        t_727[k] = f_14 * hsh_398[k]
                   + f_1 * isg0_389[k]
                   - f_2 * isg1_389[k]
                   + f_3 * pc_z[k] * ish_545[k];

        t_728[k] = pa_y[k] * hsi0_560[k]
                   - f_10 * pc_y[k] * hsi1_560[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pa_y, pc_x, pc_y, hsi0_562, hsi1_562, isg0_391, \
                         isg0_393, isg1_391, isg1_393, ish_547, \
                         ish_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_16 * isg0_391[k]
                   - f_17 * isg1_391[k]
                   + f_3 * pc_x[k] * ish_547[k];

        t_730[k] = pa_y[k] * hsi0_562[k]
                   - f_10 * pc_y[k] * hsi1_562[k];

        t_731[k] = f_8 * isg0_393[k]
                   - f_9 * isg1_393[k]
                   + f_3 * pc_x[k] * ish_549[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_x, pc_y, hsi0_565, hsi1_565, isg0_394, \
                         isg0_396, isg1_394, isg1_396, ish_550, \
                         ish_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_8 * isg0_394[k]
                   - f_9 * isg1_394[k]
                   + f_3 * pc_x[k] * ish_550[k];

        t_733[k] = pa_y[k] * hsi0_565[k]
                   - f_10 * pc_y[k] * hsi1_565[k];

        t_734[k] = f_6 * isg0_396[k]
                   - f_7 * isg1_396[k]
                   + f_3 * pc_x[k] * ish_552[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pa_y, pc_x, pc_y, hsi0_569, hsi1_569, isg0_397, \
                         isg0_398, isg1_397, isg1_398, ish_553, \
                         ish_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_6 * isg0_397[k]
                   - f_7 * isg1_397[k]
                   + f_3 * pc_x[k] * ish_553[k];

        t_736[k] = f_6 * isg0_398[k]
                   - f_7 * isg1_398[k]
                   + f_3 * pc_x[k] * ish_554[k];

        t_737[k] = pa_y[k] * hsi0_569[k]
                   - f_10 * pc_y[k] * hsi1_569[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pc_x, isg0_400, isg0_401, isg0_402, isg1_400, \
                         isg1_401, isg1_402, ish_556, ish_557, \
                         ish_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_4 * isg0_400[k]
                   - f_5 * isg1_400[k]
                   + f_3 * pc_x[k] * ish_556[k];

        t_739[k] = f_4 * isg0_401[k]
                   - f_5 * isg1_401[k]
                   + f_3 * pc_x[k] * ish_557[k];

        t_740[k] = f_4 * isg0_402[k]
                   - f_5 * isg1_402[k]
                   + f_3 * pc_x[k] * ish_558[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, t_745, pa_y, pc_x, pc_y, hsi0_574, \
                         hsi1_574, isg0_403, isg1_403, ish_559, ish_561, ish_562, \
                         ish_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_4 * isg0_403[k]
                   - f_5 * isg1_403[k]
                   + f_3 * pc_x[k] * ish_559[k];

        t_742[k] = pa_y[k] * hsi0_574[k]
                   - f_10 * pc_y[k] * hsi1_574[k];

        t_743[k] = f_3 * pc_x[k] * ish_561[k];

        t_744[k] = f_3 * pc_x[k] * ish_562[k];

        t_745[k] = f_3 * pc_x[k] * ish_563[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pa_y, pc_x, pc_y, hsi0_581, hsh_435, \
                         hsi1_581, ish_564, ish_565, ish_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_3 * pc_x[k] * ish_564[k];

        t_747[k] = f_3 * pc_x[k] * ish_565[k];

        t_748[k] = f_3 * pc_x[k] * ish_566[k];

        t_749[k] = pa_y[k] * hsi0_581[k]
                   + f_0 * hsh_435[k]
                   - f_10 * pc_y[k] * hsi1_581[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pa_y, pc_y, pc_z, hsi0_583, hsi0_584, hsh_414, \
                         hsh_437, hsh_438, hsi1_583, hsi1_584, \
                         ish_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_15 * hsh_414[k]
                   + f_3 * pc_z[k] * ish_561[k];

        t_751[k] = pa_y[k] * hsi0_583[k]
                   + f_14 * hsh_437[k]
                   - f_10 * pc_y[k] * hsi1_583[k];

        t_752[k] = pa_y[k] * hsi0_584[k]
                   + f_13 * hsh_438[k]
                   - f_10 * pc_y[k] * hsi1_584[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pa_y, pc_y, hsi0_585, hsi0_587, hsh_439, \
                         hsh_440, hsi1_585, hsi1_587, ish_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = pa_y[k] * hsi0_585[k]
                   + f_12 * hsh_439[k]
                   - f_10 * pc_y[k] * hsi1_585[k];

        t_754[k] = f_11 * hsh_440[k]
                   + f_3 * pc_y[k] * ish_566[k];

        t_755[k] = pa_y[k] * hsi0_587[k]
                   - f_10 * pc_y[k] * hsi1_587[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, t_760, pc_x, pc_y, isg0_405, isg0_407, \
                         isg0_408, isg1_405, isg1_407, isg1_408, ish_567, ish_569, \
                         ish_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_1 * isg0_405[k]
                   - f_2 * isg1_405[k]
                   + f_3 * pc_x[k] * ish_567[k];

        t_757[k] = f_3 * pc_y[k] * ish_567[k];

        t_758[k] = f_16 * isg0_407[k]
                   - f_17 * isg1_407[k]
                   + f_3 * pc_x[k] * ish_569[k];

        t_759[k] = f_8 * isg0_408[k]
                   - f_9 * isg1_408[k]
                   + f_3 * pc_x[k] * ish_570[k];

        t_760[k] = f_3 * pc_y[k] * ish_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, pc_x, pc_y, isg0_410, isg0_411, isg0_412, \
                         isg1_410, isg1_411, isg1_412, ish_572, ish_573, \
                         ish_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_8 * isg0_410[k]
                   - f_9 * isg1_410[k]
                   + f_3 * pc_x[k] * ish_572[k];

        t_762[k] = f_6 * isg0_411[k]
                   - f_7 * isg1_411[k]
                   + f_3 * pc_x[k] * ish_573[k];

        t_763[k] = f_6 * isg0_412[k]
                   - f_7 * isg1_412[k]
                   + f_3 * pc_x[k] * ish_574[k];

        t_764[k] = f_3 * pc_y[k] * ish_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pc_x, isg0_414, isg0_415, isg0_416, isg1_414, \
                         isg1_415, isg1_416, ish_576, ish_577, \
                         ish_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = f_6 * isg0_414[k]
                   - f_7 * isg1_414[k]
                   + f_3 * pc_x[k] * ish_576[k];

        t_766[k] = f_4 * isg0_415[k]
                   - f_5 * isg1_415[k]
                   + f_3 * pc_x[k] * ish_577[k];

        t_767[k] = f_4 * isg0_416[k]
                   - f_5 * isg1_416[k]
                   + f_3 * pc_x[k] * ish_578[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, t_772, pc_x, pc_y, isg0_417, isg0_419, \
                         isg1_417, isg1_419, ish_576, ish_579, ish_581, ish_582, \
                         ish_583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = f_4 * isg0_417[k]
                   - f_5 * isg1_417[k]
                   + f_3 * pc_x[k] * ish_579[k];

        t_769[k] = f_3 * pc_y[k] * ish_576[k];

        t_770[k] = f_4 * isg0_419[k]
                   - f_5 * isg1_419[k]
                   + f_3 * pc_x[k] * ish_581[k];

        t_771[k] = f_3 * pc_x[k] * ish_582[k];

        t_772[k] = f_3 * pc_x[k] * ish_583[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, t_776, t_777, pc_x, pc_y, isg0_415, isg1_415, \
                         ish_582, ish_584, ish_585, ish_586, ish_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_3 * pc_x[k] * ish_584[k];

        t_774[k] = f_3 * pc_x[k] * ish_585[k];

        t_775[k] = f_3 * pc_x[k] * ish_586[k];

        t_776[k] = f_3 * pc_x[k] * ish_587[k];

        t_777[k] = f_1 * isg0_415[k]
                   - f_2 * isg1_415[k]
                   + f_3 * pc_y[k] * ish_582[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pc_y, isg0_416, isg0_417, isg0_418, isg1_416, \
                         isg1_417, isg1_418, ish_583, ish_584, \
                         ish_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_16 * isg0_416[k]
                   - f_17 * isg1_416[k]
                   + f_3 * pc_y[k] * ish_583[k];

        t_779[k] = f_8 * isg0_417[k]
                   - f_9 * isg1_417[k]
                   + f_3 * pc_y[k] * ish_584[k];

        t_780[k] = f_6 * isg0_418[k]
                   - f_7 * isg1_418[k]
                   + f_3 * pc_y[k] * ish_585[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pc_y, pc_z, hsh_440, isg0_419, isg1_419, \
                         ish_586, ish_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_4 * isg0_419[k]
                   - f_5 * isg1_419[k]
                   + f_3 * pc_y[k] * ish_586[k];

        t_782[k] = f_3 * pc_y[k] * ish_587[k];

        t_783[k] = f_0 * hsh_440[k]
                   + f_1 * isg0_419[k]
                   - f_2 * isg1_419[k]
                   + f_3 * pc_z[k] * ish_587[k];
    }
}

auto
compute_prim_isi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsi0, const size_t hsh,
                                                   const size_t hsi1, const size_t isg0,
                                                   const size_t isg1, const size_t ish,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);

    compute_prim_isi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, hsi0, hsh,
                                                              hsi1, isg0, isg1, ish, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
