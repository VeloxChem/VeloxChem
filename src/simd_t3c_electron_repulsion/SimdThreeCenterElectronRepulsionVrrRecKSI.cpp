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


#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 2.5 / q;

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

    const auto *isi0_0 = buffer.data(isi0 + 0);
    const auto *isi0_3 = buffer.data(isi0 + 3);
    const auto *isi0_5 = buffer.data(isi0 + 5);
    const auto *isi0_6 = buffer.data(isi0 + 6);
    const auto *isi0_9 = buffer.data(isi0 + 9);
    const auto *isi0_10 = buffer.data(isi0 + 10);
    const auto *isi0_14 = buffer.data(isi0 + 14);
    const auto *isi0_21 = buffer.data(isi0 + 21);
    const auto *isi0_27 = buffer.data(isi0 + 27);
    const auto *isi0_31 = buffer.data(isi0 + 31);
    const auto *isi0_34 = buffer.data(isi0 + 34);
    const auto *isi0_38 = buffer.data(isi0 + 38);
    const auto *isi0_56 = buffer.data(isi0 + 56);
    const auto *isi0_61 = buffer.data(isi0 + 61);
    const auto *isi0_65 = buffer.data(isi0 + 65);
    const auto *isi0_68 = buffer.data(isi0 + 68);
    const auto *isi0_70 = buffer.data(isi0 + 70);

    const auto *ish_0 = buffer.data(ish + 0);
    const auto *ish_1 = buffer.data(ish + 1);
    const auto *ish_2 = buffer.data(ish + 2);
    const auto *ish_3 = buffer.data(ish + 3);
    const auto *ish_5 = buffer.data(ish + 5);
    const auto *ish_6 = buffer.data(ish + 6);
    const auto *ish_9 = buffer.data(ish + 9);
    const auto *ish_15 = buffer.data(ish + 15);
    const auto *ish_17 = buffer.data(ish + 17);
    const auto *ish_18 = buffer.data(ish + 18);
    const auto *ish_20 = buffer.data(ish + 20);
    const auto *ish_21 = buffer.data(ish + 21);
    const auto *ish_24 = buffer.data(ish + 24);
    const auto *ish_26 = buffer.data(ish + 26);
    const auto *ish_27 = buffer.data(ish + 27);
    const auto *ish_30 = buffer.data(ish + 30);
    const auto *ish_36 = buffer.data(ish + 36);
    const auto *ish_38 = buffer.data(ish + 38);
    const auto *ish_39 = buffer.data(ish + 39);
    const auto *ish_40 = buffer.data(ish + 40);
    const auto *ish_41 = buffer.data(ish + 41);
    const auto *ish_42 = buffer.data(ish + 42);
    const auto *ish_44 = buffer.data(ish + 44);
    const auto *ish_47 = buffer.data(ish + 47);
    const auto *ish_50 = buffer.data(ish + 50);
    const auto *ish_51 = buffer.data(ish + 51);
    const auto *ish_57 = buffer.data(ish + 57);
    const auto *ish_58 = buffer.data(ish + 58);
    const auto *ish_59 = buffer.data(ish + 59);
    const auto *ish_60 = buffer.data(ish + 60);
    const auto *ish_62 = buffer.data(ish + 62);
    const auto *ish_63 = buffer.data(ish + 63);
    const auto *ish_66 = buffer.data(ish + 66);
    const auto *ish_69 = buffer.data(ish + 69);
    const auto *ish_73 = buffer.data(ish + 73);
    const auto *ish_78 = buffer.data(ish + 78);
    const auto *ish_80 = buffer.data(ish + 80);
    const auto *ish_81 = buffer.data(ish + 81);
    const auto *ish_82 = buffer.data(ish + 82);
    const auto *ish_83 = buffer.data(ish + 83);
    const auto *ish_99 = buffer.data(ish + 99);

    const auto *isi1_0 = buffer.data(isi1 + 0);
    const auto *isi1_3 = buffer.data(isi1 + 3);
    const auto *isi1_5 = buffer.data(isi1 + 5);
    const auto *isi1_6 = buffer.data(isi1 + 6);
    const auto *isi1_9 = buffer.data(isi1 + 9);
    const auto *isi1_10 = buffer.data(isi1 + 10);
    const auto *isi1_14 = buffer.data(isi1 + 14);
    const auto *isi1_21 = buffer.data(isi1 + 21);
    const auto *isi1_27 = buffer.data(isi1 + 27);
    const auto *isi1_31 = buffer.data(isi1 + 31);
    const auto *isi1_34 = buffer.data(isi1 + 34);
    const auto *isi1_38 = buffer.data(isi1 + 38);
    const auto *isi1_56 = buffer.data(isi1 + 56);
    const auto *isi1_61 = buffer.data(isi1 + 61);
    const auto *isi1_65 = buffer.data(isi1 + 65);
    const auto *isi1_68 = buffer.data(isi1 + 68);
    const auto *isi1_70 = buffer.data(isi1 + 70);

    const auto *ksg0_0 = buffer.data(ksg0 + 0);
    const auto *ksg0_1 = buffer.data(ksg0 + 1);
    const auto *ksg0_2 = buffer.data(ksg0 + 2);
    const auto *ksg0_3 = buffer.data(ksg0 + 3);
    const auto *ksg0_5 = buffer.data(ksg0 + 5);
    const auto *ksg0_10 = buffer.data(ksg0 + 10);
    const auto *ksg0_12 = buffer.data(ksg0 + 12);
    const auto *ksg0_13 = buffer.data(ksg0 + 13);
    const auto *ksg0_14 = buffer.data(ksg0 + 14);
    const auto *ksg0_18 = buffer.data(ksg0 + 18);
    const auto *ksg0_25 = buffer.data(ksg0 + 25);
    const auto *ksg0_26 = buffer.data(ksg0 + 26);
    const auto *ksg0_27 = buffer.data(ksg0 + 27);
    const auto *ksg0_32 = buffer.data(ksg0 + 32);
    const auto *ksg0_34 = buffer.data(ksg0 + 34);
    const auto *ksg0_35 = buffer.data(ksg0 + 35);
    const auto *ksg0_41 = buffer.data(ksg0 + 41);
    const auto *ksg0_42 = buffer.data(ksg0 + 42);
    const auto *ksg0_43 = buffer.data(ksg0 + 43);
    const auto *ksg0_44 = buffer.data(ksg0 + 44);
    const auto *ksg0_45 = buffer.data(ksg0 + 45);
    const auto *ksg0_47 = buffer.data(ksg0 + 47);
    const auto *ksg0_48 = buffer.data(ksg0 + 48);
    const auto *ksg0_50 = buffer.data(ksg0 + 50);
    const auto *ksg0_51 = buffer.data(ksg0 + 51);
    const auto *ksg0_55 = buffer.data(ksg0 + 55);
    const auto *ksg0_56 = buffer.data(ksg0 + 56);
    const auto *ksg0_57 = buffer.data(ksg0 + 57);
    const auto *ksg0_59 = buffer.data(ksg0 + 59);

    const auto *ksg1_0 = buffer.data(ksg1 + 0);
    const auto *ksg1_1 = buffer.data(ksg1 + 1);
    const auto *ksg1_2 = buffer.data(ksg1 + 2);
    const auto *ksg1_3 = buffer.data(ksg1 + 3);
    const auto *ksg1_5 = buffer.data(ksg1 + 5);
    const auto *ksg1_10 = buffer.data(ksg1 + 10);
    const auto *ksg1_12 = buffer.data(ksg1 + 12);
    const auto *ksg1_13 = buffer.data(ksg1 + 13);
    const auto *ksg1_14 = buffer.data(ksg1 + 14);
    const auto *ksg1_18 = buffer.data(ksg1 + 18);
    const auto *ksg1_25 = buffer.data(ksg1 + 25);
    const auto *ksg1_26 = buffer.data(ksg1 + 26);
    const auto *ksg1_27 = buffer.data(ksg1 + 27);
    const auto *ksg1_32 = buffer.data(ksg1 + 32);
    const auto *ksg1_34 = buffer.data(ksg1 + 34);
    const auto *ksg1_35 = buffer.data(ksg1 + 35);
    const auto *ksg1_41 = buffer.data(ksg1 + 41);
    const auto *ksg1_42 = buffer.data(ksg1 + 42);
    const auto *ksg1_43 = buffer.data(ksg1 + 43);
    const auto *ksg1_44 = buffer.data(ksg1 + 44);
    const auto *ksg1_45 = buffer.data(ksg1 + 45);
    const auto *ksg1_47 = buffer.data(ksg1 + 47);
    const auto *ksg1_48 = buffer.data(ksg1 + 48);
    const auto *ksg1_50 = buffer.data(ksg1 + 50);
    const auto *ksg1_51 = buffer.data(ksg1 + 51);
    const auto *ksg1_55 = buffer.data(ksg1 + 55);
    const auto *ksg1_56 = buffer.data(ksg1 + 56);
    const auto *ksg1_57 = buffer.data(ksg1 + 57);
    const auto *ksg1_59 = buffer.data(ksg1 + 59);

    const auto *ksh_0 = buffer.data(ksh + 0);
    const auto *ksh_1 = buffer.data(ksh + 1);
    const auto *ksh_2 = buffer.data(ksh + 2);
    const auto *ksh_3 = buffer.data(ksh + 3);
    const auto *ksh_5 = buffer.data(ksh + 5);
    const auto *ksh_6 = buffer.data(ksh + 6);
    const auto *ksh_8 = buffer.data(ksh + 8);
    const auto *ksh_9 = buffer.data(ksh + 9);
    const auto *ksh_10 = buffer.data(ksh + 10);
    const auto *ksh_14 = buffer.data(ksh + 14);
    const auto *ksh_15 = buffer.data(ksh + 15);
    const auto *ksh_17 = buffer.data(ksh + 17);
    const auto *ksh_18 = buffer.data(ksh + 18);
    const auto *ksh_19 = buffer.data(ksh + 19);
    const auto *ksh_20 = buffer.data(ksh + 20);
    const auto *ksh_21 = buffer.data(ksh + 21);
    const auto *ksh_22 = buffer.data(ksh + 22);
    const auto *ksh_24 = buffer.data(ksh + 24);
    const auto *ksh_26 = buffer.data(ksh + 26);
    const auto *ksh_27 = buffer.data(ksh + 27);
    const auto *ksh_28 = buffer.data(ksh + 28);
    const auto *ksh_30 = buffer.data(ksh + 30);
    const auto *ksh_31 = buffer.data(ksh + 31);
    const auto *ksh_36 = buffer.data(ksh + 36);
    const auto *ksh_37 = buffer.data(ksh + 37);
    const auto *ksh_38 = buffer.data(ksh + 38);
    const auto *ksh_39 = buffer.data(ksh + 39);
    const auto *ksh_40 = buffer.data(ksh + 40);
    const auto *ksh_41 = buffer.data(ksh + 41);
    const auto *ksh_42 = buffer.data(ksh + 42);
    const auto *ksh_44 = buffer.data(ksh + 44);
    const auto *ksh_46 = buffer.data(ksh + 46);
    const auto *ksh_47 = buffer.data(ksh + 47);
    const auto *ksh_49 = buffer.data(ksh + 49);
    const auto *ksh_50 = buffer.data(ksh + 50);
    const auto *ksh_51 = buffer.data(ksh + 51);
    const auto *ksh_56 = buffer.data(ksh + 56);
    const auto *ksh_57 = buffer.data(ksh + 57);
    const auto *ksh_58 = buffer.data(ksh + 58);
    const auto *ksh_59 = buffer.data(ksh + 59);
    const auto *ksh_60 = buffer.data(ksh + 60);
    const auto *ksh_61 = buffer.data(ksh + 61);
    const auto *ksh_62 = buffer.data(ksh + 62);
    const auto *ksh_63 = buffer.data(ksh + 63);
    const auto *ksh_64 = buffer.data(ksh + 64);
    const auto *ksh_65 = buffer.data(ksh + 65);
    const auto *ksh_66 = buffer.data(ksh + 66);
    const auto *ksh_68 = buffer.data(ksh + 68);
    const auto *ksh_69 = buffer.data(ksh + 69);
    const auto *ksh_70 = buffer.data(ksh + 70);
    const auto *ksh_72 = buffer.data(ksh + 72);
    const auto *ksh_73 = buffer.data(ksh + 73);
    const auto *ksh_78 = buffer.data(ksh + 78);
    const auto *ksh_79 = buffer.data(ksh + 79);
    const auto *ksh_80 = buffer.data(ksh + 80);
    const auto *ksh_81 = buffer.data(ksh + 81);
    const auto *ksh_82 = buffer.data(ksh + 82);
    const auto *ksh_83 = buffer.data(ksh + 83);
    const auto *ksh_84 = buffer.data(ksh + 84);
    const auto *ksh_86 = buffer.data(ksh + 86);
    const auto *ksh_87 = buffer.data(ksh + 87);
    const auto *ksh_89 = buffer.data(ksh + 89);
    const auto *ksh_90 = buffer.data(ksh + 90);
    const auto *ksh_93 = buffer.data(ksh + 93);
    const auto *ksh_99 = buffer.data(ksh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, ish_0, ksg0_0, \
                         ksg1_0, ksh_0, ksh_1, ksh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ish_0[k]
                 + f_1 * ksg0_0[k]
                 - f_2 * ksg1_0[k]
                 + f_3 * pc_x[k] * ksh_0[k];

        t_1[k] = f_3 * pc_y[k] * ksh_0[k];

        t_2[k] = f_3 * pc_z[k] * ksh_0[k];

        t_3[k] = f_4 * ksg0_0[k]
                 - f_5 * ksg1_0[k]
                 + f_3 * pc_y[k] * ksh_1[k];

        t_4[k] = f_3 * pc_y[k] * ksh_2[k];

        t_5[k] = f_4 * ksg0_0[k]
                 - f_5 * ksg1_0[k]
                 + f_3 * pc_z[k] * ksh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, ksg0_1, ksg0_2, ksg0_3, ksg1_1, \
                         ksg1_2, ksg1_3, ksh_3, ksh_5, ksh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * ksg0_1[k]
                 - f_7 * ksg1_1[k]
                 + f_3 * pc_y[k] * ksh_3[k];

        t_7[k] = f_3 * pc_z[k] * ksh_3[k];

        t_8[k] = f_3 * pc_y[k] * ksh_5[k];

        t_9[k] = f_6 * ksg0_2[k]
                 - f_7 * ksg1_2[k]
                 + f_3 * pc_z[k] * ksh_5[k];

        t_10[k] = f_8 * ksg0_3[k]
                  - f_9 * ksg1_3[k]
                  + f_3 * pc_y[k] * ksh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, ish_15, ksg0_5, \
                         ksg1_5, ksh_6, ksh_8, ksh_9, ksh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * ksh_6[k];

        t_12[k] = f_4 * ksg0_5[k]
                  - f_5 * ksg1_5[k]
                  + f_3 * pc_y[k] * ksh_8[k];

        t_13[k] = f_3 * pc_y[k] * ksh_9[k];

        t_14[k] = f_8 * ksg0_5[k]
                  - f_9 * ksg1_5[k]
                  + f_3 * pc_z[k] * ksh_9[k];

        t_15[k] = f_0 * ish_15[k]
                  + f_3 * pc_x[k] * ksh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, ish_17, ish_18, \
                         ish_20, ksh_10, ksh_14, ksh_17, ksh_18, \
                         ksh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * ksh_10[k];

        t_17[k] = f_0 * ish_17[k]
                  + f_3 * pc_x[k] * ksh_17[k];

        t_18[k] = f_0 * ish_18[k]
                  + f_3 * pc_x[k] * ksh_18[k];

        t_19[k] = f_3 * pc_y[k] * ksh_14[k];

        t_20[k] = f_0 * ish_20[k]
                  + f_3 * pc_x[k] * ksh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, ksg0_10, ksg0_12, ksg0_13, \
                         ksg1_10, ksg1_12, ksg1_13, ksh_15, ksh_17, \
                         ksh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * ksg0_10[k]
                  - f_2 * ksg1_10[k]
                  + f_3 * pc_y[k] * ksh_15[k];

        t_22[k] = f_3 * pc_z[k] * ksh_15[k];

        t_23[k] = f_8 * ksg0_12[k]
                  - f_9 * ksg1_12[k]
                  + f_3 * pc_y[k] * ksh_17[k];

        t_24[k] = f_6 * ksg0_13[k]
                  - f_7 * ksg1_13[k]
                  + f_3 * pc_y[k] * ksh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, isi0_0, ish_0, \
                         isi1_0, ksg0_14, ksg1_14, ksh_19, ksh_20, \
                         ksh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * ksg0_14[k]
                  - f_5 * ksg1_14[k]
                  + f_3 * pc_y[k] * ksh_19[k];

        t_26[k] = f_3 * pc_y[k] * ksh_20[k];

        t_27[k] = f_1 * ksg0_14[k]
                  - f_2 * ksg1_14[k]
                  + f_3 * pc_z[k] * ksh_20[k];

        t_28[k] = pa_y[k] * isi0_0[k]
                  - f_10 * pc_y[k] * isi1_0[k];

        t_29[k] = f_11 * ish_0[k]
                  + f_3 * pc_y[k] * ksh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, isi0_3, isi0_5, ish_1, \
                         isi1_3, isi1_5, ksh_21, ksh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * ksh_21[k];

        t_31[k] = pa_y[k] * isi0_3[k]
                  + f_12 * ish_1[k]
                  - f_10 * pc_y[k] * isi1_3[k];

        t_32[k] = f_3 * pc_z[k] * ksh_22[k];

        t_33[k] = pa_y[k] * isi0_5[k]
                  - f_10 * pc_y[k] * isi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, isi0_6, isi0_9, ish_3, \
                         ish_5, isi1_6, isi1_9, ksh_24, ksh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * isi0_6[k]
                  + f_13 * ish_3[k]
                  - f_10 * pc_y[k] * isi1_6[k];

        t_35[k] = f_3 * pc_z[k] * ksh_24[k];

        t_36[k] = f_11 * ish_5[k]
                  + f_3 * pc_y[k] * ksh_26[k];

        t_37[k] = pa_y[k] * isi0_9[k]
                  - f_10 * pc_y[k] * isi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, isi0_10, ish_6, ish_9, \
                         isi1_10, ksg0_18, ksg1_18, ksh_27, ksh_28, \
                         ksh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * isi0_10[k]
                  + f_14 * ish_6[k]
                  - f_10 * pc_y[k] * isi1_10[k];

        t_39[k] = f_3 * pc_z[k] * ksh_27[k];

        t_40[k] = f_4 * ksg0_18[k]
                  - f_5 * ksg1_18[k]
                  + f_3 * pc_z[k] * ksh_28[k];

        t_41[k] = f_11 * ish_9[k]
                  + f_3 * pc_y[k] * ksh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, isi0_14, ish_36, \
                         ish_38, isi1_14, ksh_31, ksh_36, ksh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * isi0_14[k]
                  - f_10 * pc_y[k] * isi1_14[k];

        t_43[k] = f_15 * ish_36[k]
                  + f_3 * pc_x[k] * ksh_36[k];

        t_44[k] = f_3 * pc_z[k] * ksh_31[k];

        t_45[k] = f_15 * ish_38[k]
                  + f_3 * pc_x[k] * ksh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, ish_15, ish_39, ish_40, ish_41, \
                         ksg0_25, ksg1_25, ksh_36, ksh_39, ksh_40, \
                         ksh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_15 * ish_39[k]
                  + f_3 * pc_x[k] * ksh_39[k];

        t_47[k] = f_15 * ish_40[k]
                  + f_3 * pc_x[k] * ksh_40[k];

        t_48[k] = f_15 * ish_41[k]
                  + f_3 * pc_x[k] * ksh_41[k];

        t_49[k] = f_11 * ish_15[k]
                  + f_1 * ksg0_25[k]
                  - f_2 * ksg1_25[k]
                  + f_3 * pc_y[k] * ksh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, ksg0_25, ksg0_26, ksg0_27, ksg1_25, \
                         ksg1_26, ksg1_27, ksh_36, ksh_37, ksh_38, \
                         ksh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * ksh_36[k];

        t_51[k] = f_4 * ksg0_25[k]
                  - f_5 * ksg1_25[k]
                  + f_3 * pc_z[k] * ksh_37[k];

        t_52[k] = f_6 * ksg0_26[k]
                  - f_7 * ksg1_26[k]
                  + f_3 * pc_z[k] * ksh_38[k];

        t_53[k] = f_8 * ksg0_27[k]
                  - f_9 * ksg1_27[k]
                  + f_3 * pc_z[k] * ksh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, isi0_0, isi0_27, \
                         ish_20, isi1_0, isi1_27, ksh_41, ksh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * ish_20[k]
                  + f_3 * pc_y[k] * ksh_41[k];

        t_55[k] = pa_y[k] * isi0_27[k]
                  - f_10 * pc_y[k] * isi1_27[k];

        t_56[k] = pa_z[k] * isi0_0[k]
                  - f_10 * pc_z[k] * isi1_0[k];

        t_57[k] = f_3 * pc_y[k] * ksh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, isi0_3, isi0_5, ish_0, \
                         ish_2, isi1_3, isi1_5, ksh_42, ksh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * ish_0[k]
                  + f_3 * pc_z[k] * ksh_42[k];

        t_59[k] = pa_z[k] * isi0_3[k]
                  - f_10 * pc_z[k] * isi1_3[k];

        t_60[k] = f_3 * pc_y[k] * ksh_44[k];

        t_61[k] = pa_z[k] * isi0_5[k]
                  + f_12 * ish_2[k]
                  - f_10 * pc_z[k] * isi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, isi0_6, isi0_9, ish_5, \
                         isi1_6, isi1_9, ksg0_32, ksg1_32, ksh_46, \
                         ksh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * isi0_6[k]
                  - f_10 * pc_z[k] * isi1_6[k];

        t_63[k] = f_4 * ksg0_32[k]
                  - f_5 * ksg1_32[k]
                  + f_3 * pc_y[k] * ksh_46[k];

        t_64[k] = f_3 * pc_y[k] * ksh_47[k];

        t_65[k] = pa_z[k] * isi0_9[k]
                  + f_13 * ish_5[k]
                  - f_10 * pc_z[k] * isi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, isi0_10, isi1_10, ksg0_34, \
                         ksg0_35, ksg1_34, ksg1_35, ksh_49, ksh_50, \
                         ksh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * isi0_10[k]
                  - f_10 * pc_z[k] * isi1_10[k];

        t_67[k] = f_6 * ksg0_34[k]
                  - f_7 * ksg1_34[k]
                  + f_3 * pc_y[k] * ksh_49[k];

        t_68[k] = f_4 * ksg0_35[k]
                  - f_5 * ksg1_35[k]
                  + f_3 * pc_y[k] * ksh_50[k];

        t_69[k] = f_3 * pc_y[k] * ksh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, isi0_14, ish_9, ish_57, \
                         ish_58, ish_59, isi1_14, ksh_57, ksh_58, \
                         ksh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * isi0_14[k]
                  + f_14 * ish_9[k]
                  - f_10 * pc_z[k] * isi1_14[k];

        t_71[k] = f_15 * ish_57[k]
                  + f_3 * pc_x[k] * ksh_57[k];

        t_72[k] = f_15 * ish_58[k]
                  + f_3 * pc_x[k] * ksh_58[k];

        t_73[k] = f_15 * ish_59[k]
                  + f_3 * pc_x[k] * ksh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, isi0_21, ish_60, \
                         ish_62, isi1_21, ksh_56, ksh_60, ksh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_15 * ish_60[k]
                  + f_3 * pc_x[k] * ksh_60[k];

        t_75[k] = f_3 * pc_y[k] * ksh_56[k];

        t_76[k] = f_15 * ish_62[k]
                  + f_3 * pc_x[k] * ksh_62[k];

        t_77[k] = pa_z[k] * isi0_21[k]
                  - f_10 * pc_z[k] * isi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, ksg0_41, ksg0_42, ksg0_43, ksg1_41, ksg1_42, \
                         ksg1_43, ksh_58, ksh_59, ksh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_16 * ksg0_41[k]
                  - f_17 * ksg1_41[k]
                  + f_3 * pc_y[k] * ksh_58[k];

        t_79[k] = f_8 * ksg0_42[k]
                  - f_9 * ksg1_42[k]
                  + f_3 * pc_y[k] * ksh_59[k];

        t_80[k] = f_6 * ksg0_43[k]
                  - f_7 * ksg1_43[k]
                  + f_3 * pc_y[k] * ksh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, ish_20, ish_63, ksg0_44, \
                         ksg0_45, ksg1_44, ksg1_45, ksh_61, ksh_62, \
                         ksh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * ksg0_44[k]
                  - f_5 * ksg1_44[k]
                  + f_3 * pc_y[k] * ksh_61[k];

        t_82[k] = f_3 * pc_y[k] * ksh_62[k];

        t_83[k] = f_11 * ish_20[k]
                  + f_1 * ksg0_44[k]
                  - f_2 * ksg1_44[k]
                  + f_3 * pc_z[k] * ksh_62[k];

        t_84[k] = f_18 * ish_63[k]
                  + f_1 * ksg0_45[k]
                  - f_2 * ksg1_45[k]
                  + f_3 * pc_x[k] * ksh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, ish_21, ish_66, ksg0_48, \
                         ksg1_48, ksh_63, ksh_64, ksh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * ish_21[k]
                  + f_3 * pc_y[k] * ksh_63[k];

        t_86[k] = f_3 * pc_z[k] * ksh_63[k];

        t_87[k] = f_18 * ish_66[k]
                  + f_8 * ksg0_48[k]
                  - f_9 * ksg1_48[k]
                  + f_3 * pc_x[k] * ksh_66[k];

        t_88[k] = f_3 * pc_z[k] * ksh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, ish_69, ksg0_45, ksg0_51, ksg1_45, \
                         ksg1_51, ksh_65, ksh_66, ksh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * ksg0_45[k]
                  - f_5 * ksg1_45[k]
                  + f_3 * pc_z[k] * ksh_65[k];

        t_90[k] = f_18 * ish_69[k]
                  + f_6 * ksg0_51[k]
                  - f_7 * ksg1_51[k]
                  + f_3 * pc_x[k] * ksh_69[k];

        t_91[k] = f_3 * pc_z[k] * ksh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, ish_26, ish_73, ksg0_47, \
                         ksg0_55, ksg1_47, ksg1_55, ksh_68, ksh_69, \
                         ksh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * ish_26[k]
                  + f_3 * pc_y[k] * ksh_68[k];

        t_93[k] = f_6 * ksg0_47[k]
                  - f_7 * ksg1_47[k]
                  + f_3 * pc_z[k] * ksh_68[k];

        t_94[k] = f_18 * ish_73[k]
                  + f_4 * ksg0_55[k]
                  - f_5 * ksg1_55[k]
                  + f_3 * pc_x[k] * ksh_73[k];

        t_95[k] = f_3 * pc_z[k] * ksh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, ish_30, ish_78, ksg0_48, \
                         ksg0_50, ksg1_48, ksg1_50, ksh_70, ksh_72, \
                         ksh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * ksg0_48[k]
                  - f_5 * ksg1_48[k]
                  + f_3 * pc_z[k] * ksh_70[k];

        t_97[k] = f_12 * ish_30[k]
                  + f_3 * pc_y[k] * ksh_72[k];

        t_98[k] = f_8 * ksg0_50[k]
                  - f_9 * ksg1_50[k]
                  + f_3 * pc_z[k] * ksh_72[k];

        t_99[k] = f_18 * ish_78[k]
                  + f_3 * pc_x[k] * ksh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, ish_80, ish_81, \
                         ish_82, ish_83, ksh_73, ksh_80, ksh_81, ksh_82, \
                         ksh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * ksh_73[k];

        t_101[k] = f_18 * ish_80[k]
                   + f_3 * pc_x[k] * ksh_80[k];

        t_102[k] = f_18 * ish_81[k]
                   + f_3 * pc_x[k] * ksh_81[k];

        t_103[k] = f_18 * ish_82[k]
                   + f_3 * pc_x[k] * ksh_82[k];

        t_104[k] = f_18 * ish_83[k]
                   + f_3 * pc_x[k] * ksh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, ish_36, ksg0_55, ksg0_56, \
                         ksg1_55, ksg1_56, ksh_78, ksh_79, ksh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * ish_36[k]
                   + f_1 * ksg0_55[k]
                   - f_2 * ksg1_55[k]
                   + f_3 * pc_y[k] * ksh_78[k];

        t_106[k] = f_3 * pc_z[k] * ksh_78[k];

        t_107[k] = f_4 * ksg0_55[k]
                   - f_5 * ksg1_55[k]
                   + f_3 * pc_z[k] * ksh_79[k];

        t_108[k] = f_6 * ksg0_56[k]
                   - f_7 * ksg1_56[k]
                   + f_3 * pc_z[k] * ksh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, isi0_56, ish_41, \
                         isi1_56, ksg0_57, ksg0_59, ksg1_57, ksg1_59, ksh_81, \
                         ksh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * ksg0_57[k]
                   - f_9 * ksg1_57[k]
                   + f_3 * pc_z[k] * ksh_81[k];

        t_110[k] = f_12 * ish_41[k]
                   + f_3 * pc_y[k] * ksh_83[k];

        t_111[k] = f_1 * ksg0_59[k]
                   - f_2 * ksg1_59[k]
                   + f_3 * pc_z[k] * ksh_83[k];

        t_112[k] = pa_y[k] * isi0_56[k]
                   - f_10 * pc_y[k] * isi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, isi0_31, ish_21, \
                         ish_42, ish_44, isi1_31, ksh_84, ksh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * ish_42[k]
                   + f_3 * pc_y[k] * ksh_84[k];

        t_114[k] = f_11 * ish_21[k]
                   + f_3 * pc_z[k] * ksh_84[k];

        t_115[k] = pa_z[k] * isi0_31[k]
                   - f_10 * pc_z[k] * isi1_31[k];

        t_116[k] = f_11 * ish_44[k]
                   + f_3 * pc_y[k] * ksh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, isi0_34, isi0_61, \
                         ish_24, ish_47, isi1_34, isi1_61, ksh_87, \
                         ksh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * isi0_61[k]
                   - f_10 * pc_y[k] * isi1_61[k];

        t_118[k] = pa_z[k] * isi0_34[k]
                   - f_10 * pc_z[k] * isi1_34[k];

        t_119[k] = f_11 * ish_24[k]
                   + f_3 * pc_z[k] * ksh_87[k];

        t_120[k] = f_11 * ish_47[k]
                   + f_3 * pc_y[k] * ksh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, isi0_38, isi0_65, \
                         ish_27, isi1_38, isi1_65, ksh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * isi0_65[k]
                   - f_10 * pc_y[k] * isi1_65[k];

        t_122[k] = pa_z[k] * isi0_38[k]
                   - f_10 * pc_z[k] * isi1_38[k];

        t_123[k] = f_11 * ish_27[k]
                   + f_3 * pc_z[k] * ksh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, isi0_68, isi0_70, \
                         ish_50, ish_51, ish_99, isi1_68, isi1_70, ksh_93, \
                         ksh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * isi0_68[k]
                   + f_12 * ish_50[k]
                   - f_10 * pc_y[k] * isi1_68[k];

        t_125[k] = f_11 * ish_51[k]
                   + f_3 * pc_y[k] * ksh_93[k];

        t_126[k] = pa_y[k] * isi0_70[k]
                   - f_10 * pc_y[k] * isi1_70[k];

        t_127[k] = f_18 * ish_99[k]
                   + f_3 * pc_x[k] * ksh_99[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
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
    const auto f_18 = 2.5 / q;

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

    const auto *isi0_49 = buffer.data(isi0 + 49);
    const auto *isi0_83 = buffer.data(isi0 + 83);
    const auto *isi0_84 = buffer.data(isi0 + 84);
    const auto *isi0_87 = buffer.data(isi0 + 87);
    const auto *isi0_90 = buffer.data(isi0 + 90);
    const auto *isi0_94 = buffer.data(isi0 + 94);
    const auto *isi0_96 = buffer.data(isi0 + 96);
    const auto *isi0_105 = buffer.data(isi0 + 105);
    const auto *isi0_140 = buffer.data(isi0 + 140);
    const auto *isi0_143 = buffer.data(isi0 + 143);
    const auto *isi0_145 = buffer.data(isi0 + 145);
    const auto *isi0_146 = buffer.data(isi0 + 146);
    const auto *isi0_149 = buffer.data(isi0 + 149);
    const auto *isi0_150 = buffer.data(isi0 + 150);
    const auto *isi0_152 = buffer.data(isi0 + 152);
    const auto *isi0_154 = buffer.data(isi0 + 154);

    const auto *ish_36 = buffer.data(ish + 36);
    const auto *ish_42 = buffer.data(ish + 42);
    const auto *ish_59 = buffer.data(ish + 59);
    const auto *ish_60 = buffer.data(ish + 60);
    const auto *ish_61 = buffer.data(ish + 61);
    const auto *ish_62 = buffer.data(ish + 62);
    const auto *ish_63 = buffer.data(ish + 63);
    const auto *ish_66 = buffer.data(ish + 66);
    const auto *ish_68 = buffer.data(ish + 68);
    const auto *ish_69 = buffer.data(ish + 69);
    const auto *ish_70 = buffer.data(ish + 70);
    const auto *ish_72 = buffer.data(ish + 72);
    const auto *ish_78 = buffer.data(ish + 78);
    const auto *ish_83 = buffer.data(ish + 83);
    const auto *ish_84 = buffer.data(ish + 84);
    const auto *ish_86 = buffer.data(ish + 86);
    const auto *ish_87 = buffer.data(ish + 87);
    const auto *ish_89 = buffer.data(ish + 89);
    const auto *ish_90 = buffer.data(ish + 90);
    const auto *ish_93 = buffer.data(ish + 93);
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
    const auto *ish_110 = buffer.data(ish + 110);
    const auto *ish_111 = buffer.data(ish + 111);
    const auto *ish_113 = buffer.data(ish + 113);
    const auto *ish_114 = buffer.data(ish + 114);
    const auto *ish_119 = buffer.data(ish + 119);
    const auto *ish_120 = buffer.data(ish + 120);
    const auto *ish_121 = buffer.data(ish + 121);
    const auto *ish_122 = buffer.data(ish + 122);
    const auto *ish_123 = buffer.data(ish + 123);
    const auto *ish_125 = buffer.data(ish + 125);
    const auto *ish_126 = buffer.data(ish + 126);
    const auto *ish_129 = buffer.data(ish + 129);
    const auto *ish_132 = buffer.data(ish + 132);
    const auto *ish_136 = buffer.data(ish + 136);
    const auto *ish_141 = buffer.data(ish + 141);
    const auto *ish_143 = buffer.data(ish + 143);
    const auto *ish_144 = buffer.data(ish + 144);
    const auto *ish_145 = buffer.data(ish + 145);
    const auto *ish_146 = buffer.data(ish + 146);
    const auto *ish_152 = buffer.data(ish + 152);
    const auto *ish_156 = buffer.data(ish + 156);
    const auto *ish_161 = buffer.data(ish + 161);
    const auto *ish_162 = buffer.data(ish + 162);
    const auto *ish_163 = buffer.data(ish + 163);
    const auto *ish_164 = buffer.data(ish + 164);
    const auto *ish_165 = buffer.data(ish + 165);
    const auto *ish_166 = buffer.data(ish + 166);
    const auto *ish_167 = buffer.data(ish + 167);
    const auto *ish_183 = buffer.data(ish + 183);
    const auto *ish_184 = buffer.data(ish + 184);
    const auto *ish_185 = buffer.data(ish + 185);
    const auto *ish_186 = buffer.data(ish + 186);
    const auto *ish_187 = buffer.data(ish + 187);
    const auto *ish_188 = buffer.data(ish + 188);

    const auto *isi1_49 = buffer.data(isi1 + 49);
    const auto *isi1_83 = buffer.data(isi1 + 83);
    const auto *isi1_84 = buffer.data(isi1 + 84);
    const auto *isi1_87 = buffer.data(isi1 + 87);
    const auto *isi1_90 = buffer.data(isi1 + 90);
    const auto *isi1_94 = buffer.data(isi1 + 94);
    const auto *isi1_96 = buffer.data(isi1 + 96);
    const auto *isi1_105 = buffer.data(isi1 + 105);
    const auto *isi1_140 = buffer.data(isi1 + 140);
    const auto *isi1_143 = buffer.data(isi1 + 143);
    const auto *isi1_145 = buffer.data(isi1 + 145);
    const auto *isi1_146 = buffer.data(isi1 + 146);
    const auto *isi1_149 = buffer.data(isi1 + 149);
    const auto *isi1_150 = buffer.data(isi1 + 150);
    const auto *isi1_152 = buffer.data(isi1 + 152);
    const auto *isi1_154 = buffer.data(isi1 + 154);

    const auto *ksg0_72 = buffer.data(ksg0 + 72);
    const auto *ksg0_73 = buffer.data(ksg0 + 73);
    const auto *ksg0_74 = buffer.data(ksg0 + 74);
    const auto *ksg0_75 = buffer.data(ksg0 + 75);
    const auto *ksg0_76 = buffer.data(ksg0 + 76);
    const auto *ksg0_77 = buffer.data(ksg0 + 77);
    const auto *ksg0_78 = buffer.data(ksg0 + 78);
    const auto *ksg0_79 = buffer.data(ksg0 + 79);
    const auto *ksg0_80 = buffer.data(ksg0 + 80);
    const auto *ksg0_84 = buffer.data(ksg0 + 84);
    const auto *ksg0_85 = buffer.data(ksg0 + 85);
    const auto *ksg0_86 = buffer.data(ksg0 + 86);
    const auto *ksg0_87 = buffer.data(ksg0 + 87);
    const auto *ksg0_88 = buffer.data(ksg0 + 88);
    const auto *ksg0_89 = buffer.data(ksg0 + 89);
    const auto *ksg0_90 = buffer.data(ksg0 + 90);
    const auto *ksg0_92 = buffer.data(ksg0 + 92);
    const auto *ksg0_93 = buffer.data(ksg0 + 93);
    const auto *ksg0_95 = buffer.data(ksg0 + 95);
    const auto *ksg0_96 = buffer.data(ksg0 + 96);
    const auto *ksg0_100 = buffer.data(ksg0 + 100);
    const auto *ksg0_101 = buffer.data(ksg0 + 101);
    const auto *ksg0_102 = buffer.data(ksg0 + 102);
    const auto *ksg0_104 = buffer.data(ksg0 + 104);
    const auto *ksg0_110 = buffer.data(ksg0 + 110);
    const auto *ksg0_114 = buffer.data(ksg0 + 114);
    const auto *ksg0_117 = buffer.data(ksg0 + 117);
    const auto *ksg0_118 = buffer.data(ksg0 + 118);
    const auto *ksg0_119 = buffer.data(ksg0 + 119);
    const auto *ksg0_130 = buffer.data(ksg0 + 130);
    const auto *ksg0_132 = buffer.data(ksg0 + 132);

    const auto *ksg1_72 = buffer.data(ksg1 + 72);
    const auto *ksg1_73 = buffer.data(ksg1 + 73);
    const auto *ksg1_74 = buffer.data(ksg1 + 74);
    const auto *ksg1_75 = buffer.data(ksg1 + 75);
    const auto *ksg1_76 = buffer.data(ksg1 + 76);
    const auto *ksg1_77 = buffer.data(ksg1 + 77);
    const auto *ksg1_78 = buffer.data(ksg1 + 78);
    const auto *ksg1_79 = buffer.data(ksg1 + 79);
    const auto *ksg1_80 = buffer.data(ksg1 + 80);
    const auto *ksg1_84 = buffer.data(ksg1 + 84);
    const auto *ksg1_85 = buffer.data(ksg1 + 85);
    const auto *ksg1_86 = buffer.data(ksg1 + 86);
    const auto *ksg1_87 = buffer.data(ksg1 + 87);
    const auto *ksg1_88 = buffer.data(ksg1 + 88);
    const auto *ksg1_89 = buffer.data(ksg1 + 89);
    const auto *ksg1_90 = buffer.data(ksg1 + 90);
    const auto *ksg1_92 = buffer.data(ksg1 + 92);
    const auto *ksg1_93 = buffer.data(ksg1 + 93);
    const auto *ksg1_95 = buffer.data(ksg1 + 95);
    const auto *ksg1_96 = buffer.data(ksg1 + 96);
    const auto *ksg1_100 = buffer.data(ksg1 + 100);
    const auto *ksg1_101 = buffer.data(ksg1 + 101);
    const auto *ksg1_102 = buffer.data(ksg1 + 102);
    const auto *ksg1_104 = buffer.data(ksg1 + 104);
    const auto *ksg1_110 = buffer.data(ksg1 + 110);
    const auto *ksg1_114 = buffer.data(ksg1 + 114);
    const auto *ksg1_117 = buffer.data(ksg1 + 117);
    const auto *ksg1_118 = buffer.data(ksg1 + 118);
    const auto *ksg1_119 = buffer.data(ksg1 + 119);
    const auto *ksg1_130 = buffer.data(ksg1 + 130);
    const auto *ksg1_132 = buffer.data(ksg1 + 132);

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
    const auto *ksh_109 = buffer.data(ksh + 109);
    const auto *ksh_110 = buffer.data(ksh + 110);
    const auto *ksh_111 = buffer.data(ksh + 111);
    const auto *ksh_112 = buffer.data(ksh + 112);
    const auto *ksh_113 = buffer.data(ksh + 113);
    const auto *ksh_114 = buffer.data(ksh + 114);
    const auto *ksh_119 = buffer.data(ksh + 119);
    const auto *ksh_120 = buffer.data(ksh + 120);
    const auto *ksh_121 = buffer.data(ksh + 121);
    const auto *ksh_122 = buffer.data(ksh + 122);
    const auto *ksh_123 = buffer.data(ksh + 123);
    const auto *ksh_124 = buffer.data(ksh + 124);
    const auto *ksh_125 = buffer.data(ksh + 125);
    const auto *ksh_126 = buffer.data(ksh + 126);
    const auto *ksh_127 = buffer.data(ksh + 127);
    const auto *ksh_128 = buffer.data(ksh + 128);
    const auto *ksh_129 = buffer.data(ksh + 129);
    const auto *ksh_131 = buffer.data(ksh + 131);
    const auto *ksh_132 = buffer.data(ksh + 132);
    const auto *ksh_133 = buffer.data(ksh + 133);
    const auto *ksh_135 = buffer.data(ksh + 135);
    const auto *ksh_136 = buffer.data(ksh + 136);
    const auto *ksh_141 = buffer.data(ksh + 141);
    const auto *ksh_142 = buffer.data(ksh + 142);
    const auto *ksh_143 = buffer.data(ksh + 143);
    const auto *ksh_144 = buffer.data(ksh + 144);
    const auto *ksh_145 = buffer.data(ksh + 145);
    const auto *ksh_146 = buffer.data(ksh + 146);
    const auto *ksh_147 = buffer.data(ksh + 147);
    const auto *ksh_149 = buffer.data(ksh + 149);
    const auto *ksh_150 = buffer.data(ksh + 150);
    const auto *ksh_152 = buffer.data(ksh + 152);
    const auto *ksh_153 = buffer.data(ksh + 153);
    const auto *ksh_156 = buffer.data(ksh + 156);
    const auto *ksh_161 = buffer.data(ksh + 161);
    const auto *ksh_162 = buffer.data(ksh + 162);
    const auto *ksh_163 = buffer.data(ksh + 163);
    const auto *ksh_164 = buffer.data(ksh + 164);
    const auto *ksh_165 = buffer.data(ksh + 165);
    const auto *ksh_166 = buffer.data(ksh + 166);
    const auto *ksh_167 = buffer.data(ksh + 167);
    const auto *ksh_168 = buffer.data(ksh + 168);
    const auto *ksh_170 = buffer.data(ksh + 170);
    const auto *ksh_171 = buffer.data(ksh + 171);
    const auto *ksh_173 = buffer.data(ksh + 173);
    const auto *ksh_174 = buffer.data(ksh + 174);
    const auto *ksh_177 = buffer.data(ksh + 177);
    const auto *ksh_183 = buffer.data(ksh + 183);
    const auto *ksh_184 = buffer.data(ksh + 184);
    const auto *ksh_185 = buffer.data(ksh + 185);
    const auto *ksh_186 = buffer.data(ksh + 186);
    const auto *ksh_187 = buffer.data(ksh + 187);
    const auto *ksh_188 = buffer.data(ksh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, ish_100, ish_101, ish_102, \
                         ish_103, ish_104, ksh_100, ksh_101, ksh_102, ksh_103, \
                         ksh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_18 * ish_100[k]
                   + f_3 * pc_x[k] * ksh_100[k];

        t_129[k] = f_18 * ish_101[k]
                   + f_3 * pc_x[k] * ksh_101[k];

        t_130[k] = f_18 * ish_102[k]
                   + f_3 * pc_x[k] * ksh_102[k];

        t_131[k] = f_18 * ish_103[k]
                   + f_3 * pc_x[k] * ksh_103[k];

        t_132[k] = f_18 * ish_104[k]
                   + f_3 * pc_x[k] * ksh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, isi0_49, ish_36, ish_59, \
                         isi1_49, ksg0_72, ksg1_72, ksh_99, ksh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * isi0_49[k]
                   - f_10 * pc_z[k] * isi1_49[k];

        t_134[k] = f_11 * ish_36[k]
                   + f_3 * pc_z[k] * ksh_99[k];

        t_135[k] = f_11 * ish_59[k]
                   + f_8 * ksg0_72[k]
                   - f_9 * ksg1_72[k]
                   + f_3 * pc_y[k] * ksh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, ish_60, ish_61, ish_62, ksg0_73, ksg0_74, \
                         ksg1_73, ksg1_74, ksh_102, ksh_103, ksh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * ish_60[k]
                   + f_6 * ksg0_73[k]
                   - f_7 * ksg1_73[k]
                   + f_3 * pc_y[k] * ksh_102[k];

        t_137[k] = f_11 * ish_61[k]
                   + f_4 * ksg0_74[k]
                   - f_5 * ksg1_74[k]
                   + f_3 * pc_y[k] * ksh_103[k];

        t_138[k] = f_11 * ish_62[k]
                   + f_3 * pc_y[k] * ksh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, isi0_83, ish_42, \
                         ish_105, isi1_83, ksg0_75, ksg1_75, ksh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * isi0_83[k]
                   - f_10 * pc_y[k] * isi1_83[k];

        t_140[k] = f_18 * ish_105[k]
                   + f_1 * ksg0_75[k]
                   - f_2 * ksg1_75[k]
                   + f_3 * pc_x[k] * ksh_105[k];

        t_141[k] = f_3 * pc_y[k] * ksh_105[k];

        t_142[k] = f_12 * ish_42[k]
                   + f_3 * pc_z[k] * ksh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, ish_110, ksg0_75, ksg0_80, ksg1_75, \
                         ksg1_80, ksh_106, ksh_107, ksh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * ksg0_75[k]
                   - f_5 * ksg1_75[k]
                   + f_3 * pc_y[k] * ksh_106[k];

        t_144[k] = f_3 * pc_y[k] * ksh_107[k];

        t_145[k] = f_18 * ish_110[k]
                   + f_8 * ksg0_80[k]
                   - f_9 * ksg1_80[k]
                   + f_3 * pc_x[k] * ksh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, ksg0_76, ksg0_77, ksg1_76, ksg1_77, \
                         ksh_108, ksh_109, ksh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * ksg0_76[k]
                   - f_7 * ksg1_76[k]
                   + f_3 * pc_y[k] * ksh_108[k];

        t_147[k] = f_4 * ksg0_77[k]
                   - f_5 * ksg1_77[k]
                   + f_3 * pc_y[k] * ksh_109[k];

        t_148[k] = f_3 * pc_y[k] * ksh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, ish_114, ksg0_78, ksg0_79, ksg0_84, \
                         ksg1_78, ksg1_79, ksg1_84, ksh_111, ksh_112, \
                         ksh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_18 * ish_114[k]
                   + f_6 * ksg0_84[k]
                   - f_7 * ksg1_84[k]
                   + f_3 * pc_x[k] * ksh_114[k];

        t_150[k] = f_8 * ksg0_78[k]
                   - f_9 * ksg1_78[k]
                   + f_3 * pc_y[k] * ksh_111[k];

        t_151[k] = f_6 * ksg0_79[k]
                   - f_7 * ksg1_79[k]
                   + f_3 * pc_y[k] * ksh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, ish_119, ish_120, ksg0_80, \
                         ksg0_89, ksg1_80, ksg1_89, ksh_113, ksh_114, ksh_119, \
                         ksh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * ksg0_80[k]
                   - f_5 * ksg1_80[k]
                   + f_3 * pc_y[k] * ksh_113[k];

        t_153[k] = f_3 * pc_y[k] * ksh_114[k];

        t_154[k] = f_18 * ish_119[k]
                   + f_4 * ksg0_89[k]
                   - f_5 * ksg1_89[k]
                   + f_3 * pc_x[k] * ksh_119[k];

        t_155[k] = f_18 * ish_120[k]
                   + f_3 * pc_x[k] * ksh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, ish_121, ish_122, \
                         ish_123, ish_125, ksh_119, ksh_121, ksh_122, ksh_123, \
                         ksh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * ish_121[k]
                   + f_3 * pc_x[k] * ksh_121[k];

        t_157[k] = f_18 * ish_122[k]
                   + f_3 * pc_x[k] * ksh_122[k];

        t_158[k] = f_18 * ish_123[k]
                   + f_3 * pc_x[k] * ksh_123[k];

        t_159[k] = f_3 * pc_y[k] * ksh_119[k];

        t_160[k] = f_18 * ish_125[k]
                   + f_3 * pc_x[k] * ksh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, ksg0_85, ksg0_86, ksg0_87, ksg1_85, \
                         ksg1_86, ksg1_87, ksh_120, ksh_121, ksh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * ksg0_85[k]
                   - f_2 * ksg1_85[k]
                   + f_3 * pc_y[k] * ksh_120[k];

        t_162[k] = f_16 * ksg0_86[k]
                   - f_17 * ksg1_86[k]
                   + f_3 * pc_y[k] * ksh_121[k];

        t_163[k] = f_8 * ksg0_87[k]
                   - f_9 * ksg1_87[k]
                   + f_3 * pc_y[k] * ksh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, ish_62, ksg0_88, ksg0_89, \
                         ksg1_88, ksg1_89, ksh_123, ksh_124, ksh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * ksg0_88[k]
                   - f_7 * ksg1_88[k]
                   + f_3 * pc_y[k] * ksh_123[k];

        t_165[k] = f_4 * ksg0_89[k]
                   - f_5 * ksg1_89[k]
                   + f_3 * pc_y[k] * ksh_124[k];

        t_166[k] = f_3 * pc_y[k] * ksh_125[k];

        t_167[k] = f_12 * ish_62[k]
                   + f_1 * ksg0_89[k]
                   - f_2 * ksg1_89[k]
                   + f_3 * pc_z[k] * ksh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, ish_63, ish_126, \
                         ish_129, ksg0_90, ksg0_93, ksg1_90, ksg1_93, ksh_126, \
                         ksh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_14 * ish_126[k]
                   + f_1 * ksg0_90[k]
                   - f_2 * ksg1_90[k]
                   + f_3 * pc_x[k] * ksh_126[k];

        t_169[k] = f_13 * ish_63[k]
                   + f_3 * pc_y[k] * ksh_126[k];

        t_170[k] = f_3 * pc_z[k] * ksh_126[k];

        t_171[k] = f_14 * ish_129[k]
                   + f_8 * ksg0_93[k]
                   - f_9 * ksg1_93[k]
                   + f_3 * pc_x[k] * ksh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, ish_132, ksg0_90, ksg0_96, \
                         ksg1_90, ksg1_96, ksh_127, ksh_128, ksh_129, \
                         ksh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * ksh_127[k];

        t_173[k] = f_4 * ksg0_90[k]
                   - f_5 * ksg1_90[k]
                   + f_3 * pc_z[k] * ksh_128[k];

        t_174[k] = f_14 * ish_132[k]
                   + f_6 * ksg0_96[k]
                   - f_7 * ksg1_96[k]
                   + f_3 * pc_x[k] * ksh_132[k];

        t_175[k] = f_3 * pc_z[k] * ksh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, ish_68, ish_136, \
                         ksg0_92, ksg0_100, ksg1_92, ksg1_100, ksh_131, ksh_132, \
                         ksh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * ish_68[k]
                   + f_3 * pc_y[k] * ksh_131[k];

        t_177[k] = f_6 * ksg0_92[k]
                   - f_7 * ksg1_92[k]
                   + f_3 * pc_z[k] * ksh_131[k];

        t_178[k] = f_14 * ish_136[k]
                   + f_4 * ksg0_100[k]
                   - f_5 * ksg1_100[k]
                   + f_3 * pc_x[k] * ksh_136[k];

        t_179[k] = f_3 * pc_z[k] * ksh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, ish_72, ish_141, \
                         ksg0_93, ksg0_95, ksg1_93, ksg1_95, ksh_133, ksh_135, \
                         ksh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * ksg0_93[k]
                   - f_5 * ksg1_93[k]
                   + f_3 * pc_z[k] * ksh_133[k];

        t_181[k] = f_13 * ish_72[k]
                   + f_3 * pc_y[k] * ksh_135[k];

        t_182[k] = f_8 * ksg0_95[k]
                   - f_9 * ksg1_95[k]
                   + f_3 * pc_z[k] * ksh_135[k];

        t_183[k] = f_14 * ish_141[k]
                   + f_3 * pc_x[k] * ksh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, ish_143, ish_144, \
                         ish_145, ish_146, ksh_136, ksh_143, ksh_144, ksh_145, \
                         ksh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * ksh_136[k];

        t_185[k] = f_14 * ish_143[k]
                   + f_3 * pc_x[k] * ksh_143[k];

        t_186[k] = f_14 * ish_144[k]
                   + f_3 * pc_x[k] * ksh_144[k];

        t_187[k] = f_14 * ish_145[k]
                   + f_3 * pc_x[k] * ksh_145[k];

        t_188[k] = f_14 * ish_146[k]
                   + f_3 * pc_x[k] * ksh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, ish_78, ksg0_100, ksg0_101, \
                         ksg1_100, ksg1_101, ksh_141, ksh_142, \
                         ksh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * ish_78[k]
                   + f_1 * ksg0_100[k]
                   - f_2 * ksg1_100[k]
                   + f_3 * pc_y[k] * ksh_141[k];

        t_190[k] = f_3 * pc_z[k] * ksh_141[k];

        t_191[k] = f_4 * ksg0_100[k]
                   - f_5 * ksg1_100[k]
                   + f_3 * pc_z[k] * ksh_142[k];

        t_192[k] = f_6 * ksg0_101[k]
                   - f_7 * ksg1_101[k]
                   + f_3 * pc_z[k] * ksh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, isi0_84, ish_83, \
                         isi1_84, ksg0_102, ksg0_104, ksg1_102, ksg1_104, ksh_144, \
                         ksh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * ksg0_102[k]
                   - f_9 * ksg1_102[k]
                   + f_3 * pc_z[k] * ksh_144[k];

        t_194[k] = f_13 * ish_83[k]
                   + f_3 * pc_y[k] * ksh_146[k];

        t_195[k] = f_1 * ksg0_104[k]
                   - f_2 * ksg1_104[k]
                   + f_3 * pc_z[k] * ksh_146[k];

        t_196[k] = pa_z[k] * isi0_84[k]
                   - f_10 * pc_z[k] * isi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, isi0_87, ish_63, \
                         ish_84, ish_86, isi1_87, ksh_147, ksh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * ish_84[k]
                   + f_3 * pc_y[k] * ksh_147[k];

        t_198[k] = f_11 * ish_63[k]
                   + f_3 * pc_z[k] * ksh_147[k];

        t_199[k] = pa_z[k] * isi0_87[k]
                   - f_10 * pc_z[k] * isi1_87[k];

        t_200[k] = f_12 * ish_86[k]
                   + f_3 * pc_y[k] * ksh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, isi0_90, ish_66, ish_152, \
                         isi1_90, ksg0_110, ksg1_110, ksh_150, \
                         ksh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_14 * ish_152[k]
                   + f_8 * ksg0_110[k]
                   - f_9 * ksg1_110[k]
                   + f_3 * pc_x[k] * ksh_152[k];

        t_202[k] = pa_z[k] * isi0_90[k]
                   - f_10 * pc_z[k] * isi1_90[k];

        t_203[k] = f_11 * ish_66[k]
                   + f_3 * pc_z[k] * ksh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, isi0_94, ish_89, \
                         ish_156, isi1_94, ksg0_114, ksg1_114, ksh_152, \
                         ksh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * ish_89[k]
                   + f_3 * pc_y[k] * ksh_152[k];

        t_205[k] = f_14 * ish_156[k]
                   + f_6 * ksg0_114[k]
                   - f_7 * ksg1_114[k]
                   + f_3 * pc_x[k] * ksh_156[k];

        t_206[k] = pa_z[k] * isi0_94[k]
                   - f_10 * pc_z[k] * isi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, isi0_96, ish_69, ish_70, \
                         ish_93, isi1_96, ksh_153, ksh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * ish_69[k]
                   + f_3 * pc_z[k] * ksh_153[k];

        t_208[k] = pa_z[k] * isi0_96[k]
                   + f_12 * ish_70[k]
                   - f_10 * pc_z[k] * isi1_96[k];

        t_209[k] = f_12 * ish_93[k]
                   + f_3 * pc_y[k] * ksh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, ish_161, ish_162, ish_163, ish_164, \
                         ksg0_119, ksg1_119, ksh_161, ksh_162, ksh_163, \
                         ksh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_14 * ish_161[k]
                   + f_4 * ksg0_119[k]
                   - f_5 * ksg1_119[k]
                   + f_3 * pc_x[k] * ksh_161[k];

        t_211[k] = f_14 * ish_162[k]
                   + f_3 * pc_x[k] * ksh_162[k];

        t_212[k] = f_14 * ish_163[k]
                   + f_3 * pc_x[k] * ksh_163[k];

        t_213[k] = f_14 * ish_164[k]
                   + f_3 * pc_x[k] * ksh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, isi0_105, ish_165, \
                         ish_166, ish_167, isi1_105, ksh_165, ksh_166, \
                         ksh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_14 * ish_165[k]
                   + f_3 * pc_x[k] * ksh_165[k];

        t_215[k] = f_14 * ish_166[k]
                   + f_3 * pc_x[k] * ksh_166[k];

        t_216[k] = f_14 * ish_167[k]
                   + f_3 * pc_x[k] * ksh_167[k];

        t_217[k] = pa_z[k] * isi0_105[k]
                   - f_10 * pc_z[k] * isi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, ish_78, ish_101, ish_102, ksg0_117, \
                         ksg0_118, ksg1_117, ksg1_118, ksh_162, ksh_164, \
                         ksh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * ish_78[k]
                   + f_3 * pc_z[k] * ksh_162[k];

        t_219[k] = f_12 * ish_101[k]
                   + f_8 * ksg0_117[k]
                   - f_9 * ksg1_117[k]
                   + f_3 * pc_y[k] * ksh_164[k];

        t_220[k] = f_12 * ish_102[k]
                   + f_6 * ksg0_118[k]
                   - f_7 * ksg1_118[k]
                   + f_3 * pc_y[k] * ksh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, isi0_140, ish_83, \
                         ish_103, ish_104, isi1_140, ksg0_119, ksg1_119, ksh_166, \
                         ksh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * ish_103[k]
                   + f_4 * ksg0_119[k]
                   - f_5 * ksg1_119[k]
                   + f_3 * pc_y[k] * ksh_166[k];

        t_222[k] = f_12 * ish_104[k]
                   + f_3 * pc_y[k] * ksh_167[k];

        t_223[k] = f_11 * ish_83[k]
                   + f_1 * ksg0_119[k]
                   - f_2 * ksg1_119[k]
                   + f_3 * pc_z[k] * ksh_167[k];

        t_224[k] = pa_y[k] * isi0_140[k]
                   - f_10 * pc_y[k] * isi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, isi0_143, ish_84, \
                         ish_105, ish_106, ish_107, isi1_143, ksh_168, \
                         ksh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * ish_105[k]
                   + f_3 * pc_y[k] * ksh_168[k];

        t_226[k] = f_12 * ish_84[k]
                   + f_3 * pc_z[k] * ksh_168[k];

        t_227[k] = pa_y[k] * isi0_143[k]
                   + f_12 * ish_106[k]
                   - f_10 * pc_y[k] * isi1_143[k];

        t_228[k] = f_11 * ish_107[k]
                   + f_3 * pc_y[k] * ksh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, isi0_145, isi0_146, \
                         ish_87, ish_108, ish_110, isi1_145, isi1_146, ksh_171, \
                         ksh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * isi0_145[k]
                   - f_10 * pc_y[k] * isi1_145[k];

        t_230[k] = pa_y[k] * isi0_146[k]
                   + f_13 * ish_108[k]
                   - f_10 * pc_y[k] * isi1_146[k];

        t_231[k] = f_12 * ish_87[k]
                   + f_3 * pc_z[k] * ksh_171[k];

        t_232[k] = f_11 * ish_110[k]
                   + f_3 * pc_y[k] * ksh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, isi0_149, isi0_150, ish_90, \
                         ish_111, isi1_149, isi1_150, ksh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * isi0_149[k]
                   - f_10 * pc_y[k] * isi1_149[k];

        t_234[k] = pa_y[k] * isi0_150[k]
                   + f_14 * ish_111[k]
                   - f_10 * pc_y[k] * isi1_150[k];

        t_235[k] = f_12 * ish_90[k]
                   + f_3 * pc_z[k] * ksh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, isi0_152, isi0_154, \
                         ish_113, ish_114, ish_183, isi1_152, isi1_154, ksh_177, \
                         ksh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * isi0_152[k]
                   + f_12 * ish_113[k]
                   - f_10 * pc_y[k] * isi1_152[k];

        t_237[k] = f_11 * ish_114[k]
                   + f_3 * pc_y[k] * ksh_177[k];

        t_238[k] = pa_y[k] * isi0_154[k]
                   - f_10 * pc_y[k] * isi1_154[k];

        t_239[k] = f_14 * ish_183[k]
                   + f_3 * pc_x[k] * ksh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, ish_184, ish_185, ish_186, \
                         ish_187, ish_188, ksh_184, ksh_185, ksh_186, ksh_187, \
                         ksh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_14 * ish_184[k]
                   + f_3 * pc_x[k] * ksh_184[k];

        t_241[k] = f_14 * ish_185[k]
                   + f_3 * pc_x[k] * ksh_185[k];

        t_242[k] = f_14 * ish_186[k]
                   + f_3 * pc_x[k] * ksh_186[k];

        t_243[k] = f_14 * ish_187[k]
                   + f_3 * pc_x[k] * ksh_187[k];

        t_244[k] = f_14 * ish_188[k]
                   + f_3 * pc_x[k] * ksh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, ish_99, ish_120, ish_122, ksg0_130, \
                         ksg0_132, ksg1_130, ksg1_132, ksh_183, \
                         ksh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * ish_120[k]
                   + f_1 * ksg0_130[k]
                   - f_2 * ksg1_130[k]
                   + f_3 * pc_y[k] * ksh_183[k];

        t_246[k] = f_12 * ish_99[k]
                   + f_3 * pc_z[k] * ksh_183[k];

        t_247[k] = f_11 * ish_122[k]
                   + f_8 * ksg0_132[k]
                   - f_9 * ksg1_132[k]
                   + f_3 * pc_y[k] * ksh_185[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
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

    const auto *isi0_167 = buffer.data(isi0 + 167);
    const auto *isi0_168 = buffer.data(isi0 + 168);
    const auto *isi0_171 = buffer.data(isi0 + 171);
    const auto *isi0_174 = buffer.data(isi0 + 174);
    const auto *isi0_178 = buffer.data(isi0 + 178);
    const auto *isi0_180 = buffer.data(isi0 + 180);
    const auto *isi0_189 = buffer.data(isi0 + 189);

    const auto *ish_105 = buffer.data(ish + 105);
    const auto *ish_123 = buffer.data(ish + 123);
    const auto *ish_124 = buffer.data(ish + 124);
    const auto *ish_125 = buffer.data(ish + 125);
    const auto *ish_126 = buffer.data(ish + 126);
    const auto *ish_129 = buffer.data(ish + 129);
    const auto *ish_131 = buffer.data(ish + 131);
    const auto *ish_132 = buffer.data(ish + 132);
    const auto *ish_133 = buffer.data(ish + 133);
    const auto *ish_135 = buffer.data(ish + 135);
    const auto *ish_141 = buffer.data(ish + 141);
    const auto *ish_146 = buffer.data(ish + 146);
    const auto *ish_147 = buffer.data(ish + 147);
    const auto *ish_149 = buffer.data(ish + 149);
    const auto *ish_150 = buffer.data(ish + 150);
    const auto *ish_152 = buffer.data(ish + 152);
    const auto *ish_153 = buffer.data(ish + 153);
    const auto *ish_156 = buffer.data(ish + 156);
    const auto *ish_162 = buffer.data(ish + 162);
    const auto *ish_164 = buffer.data(ish + 164);
    const auto *ish_165 = buffer.data(ish + 165);
    const auto *ish_166 = buffer.data(ish + 166);
    const auto *ish_167 = buffer.data(ish + 167);
    const auto *ish_168 = buffer.data(ish + 168);
    const auto *ish_170 = buffer.data(ish + 170);
    const auto *ish_173 = buffer.data(ish + 173);
    const auto *ish_177 = buffer.data(ish + 177);
    const auto *ish_183 = buffer.data(ish + 183);
    const auto *ish_185 = buffer.data(ish + 185);
    const auto *ish_186 = buffer.data(ish + 186);
    const auto *ish_187 = buffer.data(ish + 187);
    const auto *ish_189 = buffer.data(ish + 189);
    const auto *ish_194 = buffer.data(ish + 194);
    const auto *ish_198 = buffer.data(ish + 198);
    const auto *ish_203 = buffer.data(ish + 203);
    const auto *ish_204 = buffer.data(ish + 204);
    const auto *ish_205 = buffer.data(ish + 205);
    const auto *ish_206 = buffer.data(ish + 206);
    const auto *ish_207 = buffer.data(ish + 207);
    const auto *ish_209 = buffer.data(ish + 209);
    const auto *ish_210 = buffer.data(ish + 210);
    const auto *ish_213 = buffer.data(ish + 213);
    const auto *ish_216 = buffer.data(ish + 216);
    const auto *ish_220 = buffer.data(ish + 220);
    const auto *ish_225 = buffer.data(ish + 225);
    const auto *ish_227 = buffer.data(ish + 227);
    const auto *ish_228 = buffer.data(ish + 228);
    const auto *ish_229 = buffer.data(ish + 229);
    const auto *ish_230 = buffer.data(ish + 230);
    const auto *ish_236 = buffer.data(ish + 236);
    const auto *ish_240 = buffer.data(ish + 240);
    const auto *ish_245 = buffer.data(ish + 245);
    const auto *ish_246 = buffer.data(ish + 246);
    const auto *ish_247 = buffer.data(ish + 247);
    const auto *ish_248 = buffer.data(ish + 248);
    const auto *ish_249 = buffer.data(ish + 249);
    const auto *ish_250 = buffer.data(ish + 250);
    const auto *ish_251 = buffer.data(ish + 251);
    const auto *ish_252 = buffer.data(ish + 252);
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

    const auto *isi1_167 = buffer.data(isi1 + 167);
    const auto *isi1_168 = buffer.data(isi1 + 168);
    const auto *isi1_171 = buffer.data(isi1 + 171);
    const auto *isi1_174 = buffer.data(isi1 + 174);
    const auto *isi1_178 = buffer.data(isi1 + 178);
    const auto *isi1_180 = buffer.data(isi1 + 180);
    const auto *isi1_189 = buffer.data(isi1 + 189);

    const auto *ksg0_133 = buffer.data(ksg0 + 133);
    const auto *ksg0_134 = buffer.data(ksg0 + 134);
    const auto *ksg0_135 = buffer.data(ksg0 + 135);
    const auto *ksg0_136 = buffer.data(ksg0 + 136);
    const auto *ksg0_137 = buffer.data(ksg0 + 137);
    const auto *ksg0_138 = buffer.data(ksg0 + 138);
    const auto *ksg0_139 = buffer.data(ksg0 + 139);
    const auto *ksg0_140 = buffer.data(ksg0 + 140);
    const auto *ksg0_144 = buffer.data(ksg0 + 144);
    const auto *ksg0_145 = buffer.data(ksg0 + 145);
    const auto *ksg0_146 = buffer.data(ksg0 + 146);
    const auto *ksg0_147 = buffer.data(ksg0 + 147);
    const auto *ksg0_148 = buffer.data(ksg0 + 148);
    const auto *ksg0_149 = buffer.data(ksg0 + 149);
    const auto *ksg0_150 = buffer.data(ksg0 + 150);
    const auto *ksg0_152 = buffer.data(ksg0 + 152);
    const auto *ksg0_153 = buffer.data(ksg0 + 153);
    const auto *ksg0_155 = buffer.data(ksg0 + 155);
    const auto *ksg0_156 = buffer.data(ksg0 + 156);
    const auto *ksg0_160 = buffer.data(ksg0 + 160);
    const auto *ksg0_161 = buffer.data(ksg0 + 161);
    const auto *ksg0_162 = buffer.data(ksg0 + 162);
    const auto *ksg0_164 = buffer.data(ksg0 + 164);
    const auto *ksg0_170 = buffer.data(ksg0 + 170);
    const auto *ksg0_174 = buffer.data(ksg0 + 174);
    const auto *ksg0_177 = buffer.data(ksg0 + 177);
    const auto *ksg0_178 = buffer.data(ksg0 + 178);
    const auto *ksg0_179 = buffer.data(ksg0 + 179);
    const auto *ksg0_180 = buffer.data(ksg0 + 180);
    const auto *ksg0_183 = buffer.data(ksg0 + 183);
    const auto *ksg0_185 = buffer.data(ksg0 + 185);
    const auto *ksg0_186 = buffer.data(ksg0 + 186);
    const auto *ksg0_189 = buffer.data(ksg0 + 189);
    const auto *ksg0_190 = buffer.data(ksg0 + 190);
    const auto *ksg0_192 = buffer.data(ksg0 + 192);
    const auto *ksg0_193 = buffer.data(ksg0 + 193);
    const auto *ksg0_194 = buffer.data(ksg0 + 194);

    const auto *ksg1_133 = buffer.data(ksg1 + 133);
    const auto *ksg1_134 = buffer.data(ksg1 + 134);
    const auto *ksg1_135 = buffer.data(ksg1 + 135);
    const auto *ksg1_136 = buffer.data(ksg1 + 136);
    const auto *ksg1_137 = buffer.data(ksg1 + 137);
    const auto *ksg1_138 = buffer.data(ksg1 + 138);
    const auto *ksg1_139 = buffer.data(ksg1 + 139);
    const auto *ksg1_140 = buffer.data(ksg1 + 140);
    const auto *ksg1_144 = buffer.data(ksg1 + 144);
    const auto *ksg1_145 = buffer.data(ksg1 + 145);
    const auto *ksg1_146 = buffer.data(ksg1 + 146);
    const auto *ksg1_147 = buffer.data(ksg1 + 147);
    const auto *ksg1_148 = buffer.data(ksg1 + 148);
    const auto *ksg1_149 = buffer.data(ksg1 + 149);
    const auto *ksg1_150 = buffer.data(ksg1 + 150);
    const auto *ksg1_152 = buffer.data(ksg1 + 152);
    const auto *ksg1_153 = buffer.data(ksg1 + 153);
    const auto *ksg1_155 = buffer.data(ksg1 + 155);
    const auto *ksg1_156 = buffer.data(ksg1 + 156);
    const auto *ksg1_160 = buffer.data(ksg1 + 160);
    const auto *ksg1_161 = buffer.data(ksg1 + 161);
    const auto *ksg1_162 = buffer.data(ksg1 + 162);
    const auto *ksg1_164 = buffer.data(ksg1 + 164);
    const auto *ksg1_170 = buffer.data(ksg1 + 170);
    const auto *ksg1_174 = buffer.data(ksg1 + 174);
    const auto *ksg1_177 = buffer.data(ksg1 + 177);
    const auto *ksg1_178 = buffer.data(ksg1 + 178);
    const auto *ksg1_179 = buffer.data(ksg1 + 179);
    const auto *ksg1_180 = buffer.data(ksg1 + 180);
    const auto *ksg1_183 = buffer.data(ksg1 + 183);
    const auto *ksg1_185 = buffer.data(ksg1 + 185);
    const auto *ksg1_186 = buffer.data(ksg1 + 186);
    const auto *ksg1_189 = buffer.data(ksg1 + 189);
    const auto *ksg1_190 = buffer.data(ksg1 + 190);
    const auto *ksg1_192 = buffer.data(ksg1 + 192);
    const auto *ksg1_193 = buffer.data(ksg1 + 193);
    const auto *ksg1_194 = buffer.data(ksg1 + 194);

    const auto *ksh_186 = buffer.data(ksh + 186);
    const auto *ksh_187 = buffer.data(ksh + 187);
    const auto *ksh_188 = buffer.data(ksh + 188);
    const auto *ksh_189 = buffer.data(ksh + 189);
    const auto *ksh_190 = buffer.data(ksh + 190);
    const auto *ksh_191 = buffer.data(ksh + 191);
    const auto *ksh_192 = buffer.data(ksh + 192);
    const auto *ksh_193 = buffer.data(ksh + 193);
    const auto *ksh_194 = buffer.data(ksh + 194);
    const auto *ksh_195 = buffer.data(ksh + 195);
    const auto *ksh_196 = buffer.data(ksh + 196);
    const auto *ksh_197 = buffer.data(ksh + 197);
    const auto *ksh_198 = buffer.data(ksh + 198);
    const auto *ksh_203 = buffer.data(ksh + 203);
    const auto *ksh_204 = buffer.data(ksh + 204);
    const auto *ksh_205 = buffer.data(ksh + 205);
    const auto *ksh_206 = buffer.data(ksh + 206);
    const auto *ksh_207 = buffer.data(ksh + 207);
    const auto *ksh_208 = buffer.data(ksh + 208);
    const auto *ksh_209 = buffer.data(ksh + 209);
    const auto *ksh_210 = buffer.data(ksh + 210);
    const auto *ksh_211 = buffer.data(ksh + 211);
    const auto *ksh_212 = buffer.data(ksh + 212);
    const auto *ksh_213 = buffer.data(ksh + 213);
    const auto *ksh_215 = buffer.data(ksh + 215);
    const auto *ksh_216 = buffer.data(ksh + 216);
    const auto *ksh_217 = buffer.data(ksh + 217);
    const auto *ksh_219 = buffer.data(ksh + 219);
    const auto *ksh_220 = buffer.data(ksh + 220);
    const auto *ksh_225 = buffer.data(ksh + 225);
    const auto *ksh_226 = buffer.data(ksh + 226);
    const auto *ksh_227 = buffer.data(ksh + 227);
    const auto *ksh_228 = buffer.data(ksh + 228);
    const auto *ksh_229 = buffer.data(ksh + 229);
    const auto *ksh_230 = buffer.data(ksh + 230);
    const auto *ksh_231 = buffer.data(ksh + 231);
    const auto *ksh_233 = buffer.data(ksh + 233);
    const auto *ksh_234 = buffer.data(ksh + 234);
    const auto *ksh_236 = buffer.data(ksh + 236);
    const auto *ksh_237 = buffer.data(ksh + 237);
    const auto *ksh_240 = buffer.data(ksh + 240);
    const auto *ksh_245 = buffer.data(ksh + 245);
    const auto *ksh_246 = buffer.data(ksh + 246);
    const auto *ksh_247 = buffer.data(ksh + 247);
    const auto *ksh_248 = buffer.data(ksh + 248);
    const auto *ksh_249 = buffer.data(ksh + 249);
    const auto *ksh_250 = buffer.data(ksh + 250);
    const auto *ksh_251 = buffer.data(ksh + 251);
    const auto *ksh_252 = buffer.data(ksh + 252);
    const auto *ksh_254 = buffer.data(ksh + 254);
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

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, ish_123, ish_124, ish_125, ksg0_133, \
                         ksg0_134, ksg1_133, ksg1_134, ksh_186, ksh_187, \
                         ksh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * ish_123[k]
                   + f_6 * ksg0_133[k]
                   - f_7 * ksg1_133[k]
                   + f_3 * pc_y[k] * ksh_186[k];

        t_249[k] = f_11 * ish_124[k]
                   + f_4 * ksg0_134[k]
                   - f_5 * ksg1_134[k]
                   + f_3 * pc_y[k] * ksh_187[k];

        t_250[k] = f_11 * ish_125[k]
                   + f_3 * pc_y[k] * ksh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, isi0_167, \
                         ish_105, ish_189, isi1_167, ksg0_135, ksg1_135, \
                         ksh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * isi0_167[k]
                   - f_10 * pc_y[k] * isi1_167[k];

        t_252[k] = f_14 * ish_189[k]
                   + f_1 * ksg0_135[k]
                   - f_2 * ksg1_135[k]
                   + f_3 * pc_x[k] * ksh_189[k];

        t_253[k] = f_3 * pc_y[k] * ksh_189[k];

        t_254[k] = f_13 * ish_105[k]
                   + f_3 * pc_z[k] * ksh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, ish_194, ksg0_135, ksg0_140, \
                         ksg1_135, ksg1_140, ksh_190, ksh_191, \
                         ksh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * ksg0_135[k]
                   - f_5 * ksg1_135[k]
                   + f_3 * pc_y[k] * ksh_190[k];

        t_256[k] = f_3 * pc_y[k] * ksh_191[k];

        t_257[k] = f_14 * ish_194[k]
                   + f_8 * ksg0_140[k]
                   - f_9 * ksg1_140[k]
                   + f_3 * pc_x[k] * ksh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, ksg0_136, ksg0_137, ksg1_136, ksg1_137, \
                         ksh_192, ksh_193, ksh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * ksg0_136[k]
                   - f_7 * ksg1_136[k]
                   + f_3 * pc_y[k] * ksh_192[k];

        t_259[k] = f_4 * ksg0_137[k]
                   - f_5 * ksg1_137[k]
                   + f_3 * pc_y[k] * ksh_193[k];

        t_260[k] = f_3 * pc_y[k] * ksh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, ish_198, ksg0_138, ksg0_139, \
                         ksg0_144, ksg1_138, ksg1_139, ksg1_144, ksh_195, ksh_196, \
                         ksh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_14 * ish_198[k]
                   + f_6 * ksg0_144[k]
                   - f_7 * ksg1_144[k]
                   + f_3 * pc_x[k] * ksh_198[k];

        t_262[k] = f_8 * ksg0_138[k]
                   - f_9 * ksg1_138[k]
                   + f_3 * pc_y[k] * ksh_195[k];

        t_263[k] = f_6 * ksg0_139[k]
                   - f_7 * ksg1_139[k]
                   + f_3 * pc_y[k] * ksh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, ish_203, ish_204, ksg0_140, \
                         ksg0_149, ksg1_140, ksg1_149, ksh_197, ksh_198, ksh_203, \
                         ksh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * ksg0_140[k]
                   - f_5 * ksg1_140[k]
                   + f_3 * pc_y[k] * ksh_197[k];

        t_265[k] = f_3 * pc_y[k] * ksh_198[k];

        t_266[k] = f_14 * ish_203[k]
                   + f_4 * ksg0_149[k]
                   - f_5 * ksg1_149[k]
                   + f_3 * pc_x[k] * ksh_203[k];

        t_267[k] = f_14 * ish_204[k]
                   + f_3 * pc_x[k] * ksh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, ish_205, ish_206, \
                         ish_207, ish_209, ksh_203, ksh_205, ksh_206, ksh_207, \
                         ksh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_14 * ish_205[k]
                   + f_3 * pc_x[k] * ksh_205[k];

        t_269[k] = f_14 * ish_206[k]
                   + f_3 * pc_x[k] * ksh_206[k];

        t_270[k] = f_14 * ish_207[k]
                   + f_3 * pc_x[k] * ksh_207[k];

        t_271[k] = f_3 * pc_y[k] * ksh_203[k];

        t_272[k] = f_14 * ish_209[k]
                   + f_3 * pc_x[k] * ksh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, ksg0_145, ksg0_146, ksg0_147, ksg1_145, \
                         ksg1_146, ksg1_147, ksh_204, ksh_205, \
                         ksh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * ksg0_145[k]
                   - f_2 * ksg1_145[k]
                   + f_3 * pc_y[k] * ksh_204[k];

        t_274[k] = f_16 * ksg0_146[k]
                   - f_17 * ksg1_146[k]
                   + f_3 * pc_y[k] * ksh_205[k];

        t_275[k] = f_8 * ksg0_147[k]
                   - f_9 * ksg1_147[k]
                   + f_3 * pc_y[k] * ksh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, ish_125, ksg0_148, ksg0_149, \
                         ksg1_148, ksg1_149, ksh_207, ksh_208, \
                         ksh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * ksg0_148[k]
                   - f_7 * ksg1_148[k]
                   + f_3 * pc_y[k] * ksh_207[k];

        t_277[k] = f_4 * ksg0_149[k]
                   - f_5 * ksg1_149[k]
                   + f_3 * pc_y[k] * ksh_208[k];

        t_278[k] = f_3 * pc_y[k] * ksh_209[k];

        t_279[k] = f_13 * ish_125[k]
                   + f_1 * ksg0_149[k]
                   - f_2 * ksg1_149[k]
                   + f_3 * pc_z[k] * ksh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pc_x, pc_y, pc_z, ish_126, ish_210, \
                         ish_213, ksg0_150, ksg0_153, ksg1_150, ksg1_153, ksh_210, \
                         ksh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_13 * ish_210[k]
                   + f_1 * ksg0_150[k]
                   - f_2 * ksg1_150[k]
                   + f_3 * pc_x[k] * ksh_210[k];

        t_281[k] = f_14 * ish_126[k]
                   + f_3 * pc_y[k] * ksh_210[k];

        t_282[k] = f_3 * pc_z[k] * ksh_210[k];

        t_283[k] = f_13 * ish_213[k]
                   + f_8 * ksg0_153[k]
                   - f_9 * ksg1_153[k]
                   + f_3 * pc_x[k] * ksh_213[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, pc_z, ish_216, ksg0_150, ksg0_156, \
                         ksg1_150, ksg1_156, ksh_211, ksh_212, ksh_213, \
                         ksh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * ksh_211[k];

        t_285[k] = f_4 * ksg0_150[k]
                   - f_5 * ksg1_150[k]
                   + f_3 * pc_z[k] * ksh_212[k];

        t_286[k] = f_13 * ish_216[k]
                   + f_6 * ksg0_156[k]
                   - f_7 * ksg1_156[k]
                   + f_3 * pc_x[k] * ksh_216[k];

        t_287[k] = f_3 * pc_z[k] * ksh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pc_x, pc_y, pc_z, ish_131, ish_220, \
                         ksg0_152, ksg0_160, ksg1_152, ksg1_160, ksh_215, ksh_216, \
                         ksh_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * ish_131[k]
                   + f_3 * pc_y[k] * ksh_215[k];

        t_289[k] = f_6 * ksg0_152[k]
                   - f_7 * ksg1_152[k]
                   + f_3 * pc_z[k] * ksh_215[k];

        t_290[k] = f_13 * ish_220[k]
                   + f_4 * ksg0_160[k]
                   - f_5 * ksg1_160[k]
                   + f_3 * pc_x[k] * ksh_220[k];

        t_291[k] = f_3 * pc_z[k] * ksh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, ish_135, ish_225, \
                         ksg0_153, ksg0_155, ksg1_153, ksg1_155, ksh_217, ksh_219, \
                         ksh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * ksg0_153[k]
                   - f_5 * ksg1_153[k]
                   + f_3 * pc_z[k] * ksh_217[k];

        t_293[k] = f_14 * ish_135[k]
                   + f_3 * pc_y[k] * ksh_219[k];

        t_294[k] = f_8 * ksg0_155[k]
                   - f_9 * ksg1_155[k]
                   + f_3 * pc_z[k] * ksh_219[k];

        t_295[k] = f_13 * ish_225[k]
                   + f_3 * pc_x[k] * ksh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, ish_227, ish_228, \
                         ish_229, ish_230, ksh_220, ksh_227, ksh_228, ksh_229, \
                         ksh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * ksh_220[k];

        t_297[k] = f_13 * ish_227[k]
                   + f_3 * pc_x[k] * ksh_227[k];

        t_298[k] = f_13 * ish_228[k]
                   + f_3 * pc_x[k] * ksh_228[k];

        t_299[k] = f_13 * ish_229[k]
                   + f_3 * pc_x[k] * ksh_229[k];

        t_300[k] = f_13 * ish_230[k]
                   + f_3 * pc_x[k] * ksh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pc_y, pc_z, ish_141, ksg0_160, ksg0_161, \
                         ksg1_160, ksg1_161, ksh_225, ksh_226, \
                         ksh_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_14 * ish_141[k]
                   + f_1 * ksg0_160[k]
                   - f_2 * ksg1_160[k]
                   + f_3 * pc_y[k] * ksh_225[k];

        t_302[k] = f_3 * pc_z[k] * ksh_225[k];

        t_303[k] = f_4 * ksg0_160[k]
                   - f_5 * ksg1_160[k]
                   + f_3 * pc_z[k] * ksh_226[k];

        t_304[k] = f_6 * ksg0_161[k]
                   - f_7 * ksg1_161[k]
                   + f_3 * pc_z[k] * ksh_227[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_z, pc_y, pc_z, isi0_168, ish_146, \
                         isi1_168, ksg0_162, ksg0_164, ksg1_162, ksg1_164, ksh_228, \
                         ksh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_8 * ksg0_162[k]
                   - f_9 * ksg1_162[k]
                   + f_3 * pc_z[k] * ksh_228[k];

        t_306[k] = f_14 * ish_146[k]
                   + f_3 * pc_y[k] * ksh_230[k];

        t_307[k] = f_1 * ksg0_164[k]
                   - f_2 * ksg1_164[k]
                   + f_3 * pc_z[k] * ksh_230[k];

        t_308[k] = pa_z[k] * isi0_168[k]
                   - f_10 * pc_z[k] * isi1_168[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_z, pc_y, pc_z, isi0_171, ish_126, \
                         ish_147, ish_149, isi1_171, ksh_231, ksh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_13 * ish_147[k]
                   + f_3 * pc_y[k] * ksh_231[k];

        t_310[k] = f_11 * ish_126[k]
                   + f_3 * pc_z[k] * ksh_231[k];

        t_311[k] = pa_z[k] * isi0_171[k]
                   - f_10 * pc_z[k] * isi1_171[k];

        t_312[k] = f_13 * ish_149[k]
                   + f_3 * pc_y[k] * ksh_233[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_z, pc_x, pc_z, isi0_174, ish_129, ish_236, \
                         isi1_174, ksg0_170, ksg1_170, ksh_234, \
                         ksh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_13 * ish_236[k]
                   + f_8 * ksg0_170[k]
                   - f_9 * ksg1_170[k]
                   + f_3 * pc_x[k] * ksh_236[k];

        t_314[k] = pa_z[k] * isi0_174[k]
                   - f_10 * pc_z[k] * isi1_174[k];

        t_315[k] = f_11 * ish_129[k]
                   + f_3 * pc_z[k] * ksh_234[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pa_z, pc_x, pc_y, pc_z, isi0_178, ish_152, \
                         ish_240, isi1_178, ksg0_174, ksg1_174, ksh_236, \
                         ksh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_13 * ish_152[k]
                   + f_3 * pc_y[k] * ksh_236[k];

        t_317[k] = f_13 * ish_240[k]
                   + f_6 * ksg0_174[k]
                   - f_7 * ksg1_174[k]
                   + f_3 * pc_x[k] * ksh_240[k];

        t_318[k] = pa_z[k] * isi0_178[k]
                   - f_10 * pc_z[k] * isi1_178[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pa_z, pc_y, pc_z, isi0_180, ish_132, ish_133, \
                         ish_156, isi1_180, ksh_237, ksh_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_11 * ish_132[k]
                   + f_3 * pc_z[k] * ksh_237[k];

        t_320[k] = pa_z[k] * isi0_180[k]
                   + f_12 * ish_133[k]
                   - f_10 * pc_z[k] * isi1_180[k];

        t_321[k] = f_13 * ish_156[k]
                   + f_3 * pc_y[k] * ksh_240[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, ish_245, ish_246, ish_247, ish_248, \
                         ksg0_179, ksg1_179, ksh_245, ksh_246, ksh_247, \
                         ksh_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_13 * ish_245[k]
                   + f_4 * ksg0_179[k]
                   - f_5 * ksg1_179[k]
                   + f_3 * pc_x[k] * ksh_245[k];

        t_323[k] = f_13 * ish_246[k]
                   + f_3 * pc_x[k] * ksh_246[k];

        t_324[k] = f_13 * ish_247[k]
                   + f_3 * pc_x[k] * ksh_247[k];

        t_325[k] = f_13 * ish_248[k]
                   + f_3 * pc_x[k] * ksh_248[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pa_z, pc_x, pc_z, isi0_189, ish_249, \
                         ish_250, ish_251, isi1_189, ksh_249, ksh_250, \
                         ksh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_13 * ish_249[k]
                   + f_3 * pc_x[k] * ksh_249[k];

        t_327[k] = f_13 * ish_250[k]
                   + f_3 * pc_x[k] * ksh_250[k];

        t_328[k] = f_13 * ish_251[k]
                   + f_3 * pc_x[k] * ksh_251[k];

        t_329[k] = pa_z[k] * isi0_189[k]
                   - f_10 * pc_z[k] * isi1_189[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, pc_z, ish_141, ish_164, ish_165, ksg0_177, \
                         ksg0_178, ksg1_177, ksg1_178, ksh_246, ksh_248, \
                         ksh_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_11 * ish_141[k]
                   + f_3 * pc_z[k] * ksh_246[k];

        t_331[k] = f_13 * ish_164[k]
                   + f_8 * ksg0_177[k]
                   - f_9 * ksg1_177[k]
                   + f_3 * pc_y[k] * ksh_248[k];

        t_332[k] = f_13 * ish_165[k]
                   + f_6 * ksg0_178[k]
                   - f_7 * ksg1_178[k]
                   + f_3 * pc_y[k] * ksh_249[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pc_y, pc_z, ish_146, ish_166, ish_167, ksg0_179, \
                         ksg1_179, ksh_250, ksh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_13 * ish_166[k]
                   + f_4 * ksg0_179[k]
                   - f_5 * ksg1_179[k]
                   + f_3 * pc_y[k] * ksh_250[k];

        t_334[k] = f_13 * ish_167[k]
                   + f_3 * pc_y[k] * ksh_251[k];

        t_335[k] = f_11 * ish_146[k]
                   + f_1 * ksg0_179[k]
                   - f_2 * ksg1_179[k]
                   + f_3 * pc_z[k] * ksh_251[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, pc_x, pc_y, pc_z, ish_147, ish_168, ish_252, \
                         ksg0_180, ksg1_180, ksh_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_13 * ish_252[k]
                   + f_1 * ksg0_180[k]
                   - f_2 * ksg1_180[k]
                   + f_3 * pc_x[k] * ksh_252[k];

        t_337[k] = f_12 * ish_168[k]
                   + f_3 * pc_y[k] * ksh_252[k];

        t_338[k] = f_12 * ish_147[k]
                   + f_3 * pc_z[k] * ksh_252[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, ish_170, ish_255, ish_257, ksg0_183, \
                         ksg0_185, ksg1_183, ksg1_185, ksh_254, ksh_255, \
                         ksh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_13 * ish_255[k]
                   + f_8 * ksg0_183[k]
                   - f_9 * ksg1_183[k]
                   + f_3 * pc_x[k] * ksh_255[k];

        t_340[k] = f_12 * ish_170[k]
                   + f_3 * pc_y[k] * ksh_254[k];

        t_341[k] = f_13 * ish_257[k]
                   + f_8 * ksg0_185[k]
                   - f_9 * ksg1_185[k]
                   + f_3 * pc_x[k] * ksh_257[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pc_x, pc_y, pc_z, ish_150, ish_173, ish_258, \
                         ksg0_186, ksg1_186, ksh_255, ksh_257, \
                         ksh_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_13 * ish_258[k]
                   + f_6 * ksg0_186[k]
                   - f_7 * ksg1_186[k]
                   + f_3 * pc_x[k] * ksh_258[k];

        t_343[k] = f_12 * ish_150[k]
                   + f_3 * pc_z[k] * ksh_255[k];

        t_344[k] = f_12 * ish_173[k]
                   + f_3 * pc_y[k] * ksh_257[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_z, ish_153, ish_261, ish_262, ksg0_189, \
                         ksg0_190, ksg1_189, ksg1_190, ksh_258, ksh_261, \
                         ksh_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_13 * ish_261[k]
                   + f_6 * ksg0_189[k]
                   - f_7 * ksg1_189[k]
                   + f_3 * pc_x[k] * ksh_261[k];

        t_346[k] = f_13 * ish_262[k]
                   + f_4 * ksg0_190[k]
                   - f_5 * ksg1_190[k]
                   + f_3 * pc_x[k] * ksh_262[k];

        t_347[k] = f_12 * ish_153[k]
                   + f_3 * pc_z[k] * ksh_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pc_x, pc_y, ish_177, ish_264, ish_266, ksg0_192, \
                         ksg0_194, ksg1_192, ksg1_194, ksh_261, ksh_264, \
                         ksh_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_13 * ish_264[k]
                   + f_4 * ksg0_192[k]
                   - f_5 * ksg1_192[k]
                   + f_3 * pc_x[k] * ksh_264[k];

        t_349[k] = f_12 * ish_177[k]
                   + f_3 * pc_y[k] * ksh_261[k];

        t_350[k] = f_13 * ish_266[k]
                   + f_4 * ksg0_194[k]
                   - f_5 * ksg1_194[k]
                   + f_3 * pc_x[k] * ksh_266[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, t_355, pc_x, ish_267, ish_268, ish_269, \
                         ish_270, ish_271, ksh_267, ksh_268, ksh_269, ksh_270, \
                         ksh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_13 * ish_267[k]
                   + f_3 * pc_x[k] * ksh_267[k];

        t_352[k] = f_13 * ish_268[k]
                   + f_3 * pc_x[k] * ksh_268[k];

        t_353[k] = f_13 * ish_269[k]
                   + f_3 * pc_x[k] * ksh_269[k];

        t_354[k] = f_13 * ish_270[k]
                   + f_3 * pc_x[k] * ksh_270[k];

        t_355[k] = f_13 * ish_271[k]
                   + f_3 * pc_x[k] * ksh_271[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_x, pc_y, pc_z, ish_162, ish_183, ish_272, \
                         ksg0_190, ksg1_190, ksh_267, ksh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_13 * ish_272[k]
                   + f_3 * pc_x[k] * ksh_272[k];

        t_357[k] = f_12 * ish_183[k]
                   + f_1 * ksg0_190[k]
                   - f_2 * ksg1_190[k]
                   + f_3 * pc_y[k] * ksh_267[k];

        t_358[k] = f_12 * ish_162[k]
                   + f_3 * pc_z[k] * ksh_267[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, ish_185, ish_186, ish_187, ksg0_192, \
                         ksg0_193, ksg0_194, ksg1_192, ksg1_193, ksg1_194, ksh_269, ksh_270, \
                         ksh_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_12 * ish_185[k]
                   + f_8 * ksg0_192[k]
                   - f_9 * ksg1_192[k]
                   + f_3 * pc_y[k] * ksh_269[k];

        t_360[k] = f_12 * ish_186[k]
                   + f_6 * ksg0_193[k]
                   - f_7 * ksg1_193[k]
                   + f_3 * pc_y[k] * ksh_270[k];

        t_361[k] = f_12 * ish_187[k]
                   + f_4 * ksg0_194[k]
                   - f_5 * ksg1_194[k]
                   + f_3 * pc_y[k] * ksh_271[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
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
    const auto f_18 = 2.5 / q;

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

    const auto *isi0_252 = buffer.data(isi0 + 252);
    const auto *isi0_255 = buffer.data(isi0 + 255);
    const auto *isi0_257 = buffer.data(isi0 + 257);
    const auto *isi0_258 = buffer.data(isi0 + 258);
    const auto *isi0_261 = buffer.data(isi0 + 261);
    const auto *isi0_262 = buffer.data(isi0 + 262);
    const auto *isi0_264 = buffer.data(isi0 + 264);
    const auto *isi0_266 = buffer.data(isi0 + 266);
    const auto *isi0_279 = buffer.data(isi0 + 279);
    const auto *isi0_280 = buffer.data(isi0 + 280);
    const auto *isi0_283 = buffer.data(isi0 + 283);
    const auto *isi0_286 = buffer.data(isi0 + 286);
    const auto *isi0_290 = buffer.data(isi0 + 290);
    const auto *isi0_292 = buffer.data(isi0 + 292);
    const auto *isi0_301 = buffer.data(isi0 + 301);

    const auto *ish_167 = buffer.data(ish + 167);
    const auto *ish_168 = buffer.data(ish + 168);
    const auto *ish_171 = buffer.data(ish + 171);
    const auto *ish_174 = buffer.data(ish + 174);
    const auto *ish_183 = buffer.data(ish + 183);
    const auto *ish_188 = buffer.data(ish + 188);
    const auto *ish_189 = buffer.data(ish + 189);
    const auto *ish_190 = buffer.data(ish + 190);
    const auto *ish_191 = buffer.data(ish + 191);
    const auto *ish_192 = buffer.data(ish + 192);
    const auto *ish_194 = buffer.data(ish + 194);
    const auto *ish_195 = buffer.data(ish + 195);
    const auto *ish_197 = buffer.data(ish + 197);
    const auto *ish_198 = buffer.data(ish + 198);
    const auto *ish_204 = buffer.data(ish + 204);
    const auto *ish_206 = buffer.data(ish + 206);
    const auto *ish_207 = buffer.data(ish + 207);
    const auto *ish_208 = buffer.data(ish + 208);
    const auto *ish_209 = buffer.data(ish + 209);
    const auto *ish_210 = buffer.data(ish + 210);
    const auto *ish_213 = buffer.data(ish + 213);
    const auto *ish_215 = buffer.data(ish + 215);
    const auto *ish_216 = buffer.data(ish + 216);
    const auto *ish_217 = buffer.data(ish + 217);
    const auto *ish_219 = buffer.data(ish + 219);
    const auto *ish_225 = buffer.data(ish + 225);
    const auto *ish_230 = buffer.data(ish + 230);
    const auto *ish_231 = buffer.data(ish + 231);
    const auto *ish_233 = buffer.data(ish + 233);
    const auto *ish_236 = buffer.data(ish + 236);
    const auto *ish_240 = buffer.data(ish + 240);
    const auto *ish_248 = buffer.data(ish + 248);
    const auto *ish_249 = buffer.data(ish + 249);
    const auto *ish_250 = buffer.data(ish + 250);
    const auto *ish_251 = buffer.data(ish + 251);
    const auto *ish_252 = buffer.data(ish + 252);
    const auto *ish_288 = buffer.data(ish + 288);
    const auto *ish_289 = buffer.data(ish + 289);
    const auto *ish_290 = buffer.data(ish + 290);
    const auto *ish_291 = buffer.data(ish + 291);
    const auto *ish_292 = buffer.data(ish + 292);
    const auto *ish_293 = buffer.data(ish + 293);
    const auto *ish_294 = buffer.data(ish + 294);
    const auto *ish_299 = buffer.data(ish + 299);
    const auto *ish_303 = buffer.data(ish + 303);
    const auto *ish_308 = buffer.data(ish + 308);
    const auto *ish_309 = buffer.data(ish + 309);
    const auto *ish_310 = buffer.data(ish + 310);
    const auto *ish_311 = buffer.data(ish + 311);
    const auto *ish_312 = buffer.data(ish + 312);
    const auto *ish_314 = buffer.data(ish + 314);
    const auto *ish_315 = buffer.data(ish + 315);
    const auto *ish_318 = buffer.data(ish + 318);
    const auto *ish_321 = buffer.data(ish + 321);
    const auto *ish_325 = buffer.data(ish + 325);
    const auto *ish_330 = buffer.data(ish + 330);
    const auto *ish_332 = buffer.data(ish + 332);
    const auto *ish_333 = buffer.data(ish + 333);
    const auto *ish_334 = buffer.data(ish + 334);
    const auto *ish_335 = buffer.data(ish + 335);
    const auto *ish_341 = buffer.data(ish + 341);
    const auto *ish_345 = buffer.data(ish + 345);
    const auto *ish_350 = buffer.data(ish + 350);
    const auto *ish_351 = buffer.data(ish + 351);
    const auto *ish_352 = buffer.data(ish + 352);
    const auto *ish_353 = buffer.data(ish + 353);
    const auto *ish_354 = buffer.data(ish + 354);
    const auto *ish_355 = buffer.data(ish + 355);
    const auto *ish_356 = buffer.data(ish + 356);
    const auto *ish_357 = buffer.data(ish + 357);

    const auto *isi1_252 = buffer.data(isi1 + 252);
    const auto *isi1_255 = buffer.data(isi1 + 255);
    const auto *isi1_257 = buffer.data(isi1 + 257);
    const auto *isi1_258 = buffer.data(isi1 + 258);
    const auto *isi1_261 = buffer.data(isi1 + 261);
    const auto *isi1_262 = buffer.data(isi1 + 262);
    const auto *isi1_264 = buffer.data(isi1 + 264);
    const auto *isi1_266 = buffer.data(isi1 + 266);
    const auto *isi1_279 = buffer.data(isi1 + 279);
    const auto *isi1_280 = buffer.data(isi1 + 280);
    const auto *isi1_283 = buffer.data(isi1 + 283);
    const auto *isi1_286 = buffer.data(isi1 + 286);
    const auto *isi1_290 = buffer.data(isi1 + 290);
    const auto *isi1_292 = buffer.data(isi1 + 292);
    const auto *isi1_301 = buffer.data(isi1 + 301);

    const auto *ksg0_194 = buffer.data(ksg0 + 194);
    const auto *ksg0_205 = buffer.data(ksg0 + 205);
    const auto *ksg0_207 = buffer.data(ksg0 + 207);
    const auto *ksg0_208 = buffer.data(ksg0 + 208);
    const auto *ksg0_209 = buffer.data(ksg0 + 209);
    const auto *ksg0_210 = buffer.data(ksg0 + 210);
    const auto *ksg0_211 = buffer.data(ksg0 + 211);
    const auto *ksg0_212 = buffer.data(ksg0 + 212);
    const auto *ksg0_213 = buffer.data(ksg0 + 213);
    const auto *ksg0_214 = buffer.data(ksg0 + 214);
    const auto *ksg0_215 = buffer.data(ksg0 + 215);
    const auto *ksg0_219 = buffer.data(ksg0 + 219);
    const auto *ksg0_220 = buffer.data(ksg0 + 220);
    const auto *ksg0_221 = buffer.data(ksg0 + 221);
    const auto *ksg0_222 = buffer.data(ksg0 + 222);
    const auto *ksg0_223 = buffer.data(ksg0 + 223);
    const auto *ksg0_224 = buffer.data(ksg0 + 224);
    const auto *ksg0_225 = buffer.data(ksg0 + 225);
    const auto *ksg0_227 = buffer.data(ksg0 + 227);
    const auto *ksg0_228 = buffer.data(ksg0 + 228);
    const auto *ksg0_230 = buffer.data(ksg0 + 230);
    const auto *ksg0_231 = buffer.data(ksg0 + 231);
    const auto *ksg0_235 = buffer.data(ksg0 + 235);
    const auto *ksg0_236 = buffer.data(ksg0 + 236);
    const auto *ksg0_237 = buffer.data(ksg0 + 237);
    const auto *ksg0_239 = buffer.data(ksg0 + 239);
    const auto *ksg0_245 = buffer.data(ksg0 + 245);
    const auto *ksg0_249 = buffer.data(ksg0 + 249);
    const auto *ksg0_252 = buffer.data(ksg0 + 252);
    const auto *ksg0_253 = buffer.data(ksg0 + 253);
    const auto *ksg0_254 = buffer.data(ksg0 + 254);
    const auto *ksg0_255 = buffer.data(ksg0 + 255);

    const auto *ksg1_194 = buffer.data(ksg1 + 194);
    const auto *ksg1_205 = buffer.data(ksg1 + 205);
    const auto *ksg1_207 = buffer.data(ksg1 + 207);
    const auto *ksg1_208 = buffer.data(ksg1 + 208);
    const auto *ksg1_209 = buffer.data(ksg1 + 209);
    const auto *ksg1_210 = buffer.data(ksg1 + 210);
    const auto *ksg1_211 = buffer.data(ksg1 + 211);
    const auto *ksg1_212 = buffer.data(ksg1 + 212);
    const auto *ksg1_213 = buffer.data(ksg1 + 213);
    const auto *ksg1_214 = buffer.data(ksg1 + 214);
    const auto *ksg1_215 = buffer.data(ksg1 + 215);
    const auto *ksg1_219 = buffer.data(ksg1 + 219);
    const auto *ksg1_220 = buffer.data(ksg1 + 220);
    const auto *ksg1_221 = buffer.data(ksg1 + 221);
    const auto *ksg1_222 = buffer.data(ksg1 + 222);
    const auto *ksg1_223 = buffer.data(ksg1 + 223);
    const auto *ksg1_224 = buffer.data(ksg1 + 224);
    const auto *ksg1_225 = buffer.data(ksg1 + 225);
    const auto *ksg1_227 = buffer.data(ksg1 + 227);
    const auto *ksg1_228 = buffer.data(ksg1 + 228);
    const auto *ksg1_230 = buffer.data(ksg1 + 230);
    const auto *ksg1_231 = buffer.data(ksg1 + 231);
    const auto *ksg1_235 = buffer.data(ksg1 + 235);
    const auto *ksg1_236 = buffer.data(ksg1 + 236);
    const auto *ksg1_237 = buffer.data(ksg1 + 237);
    const auto *ksg1_239 = buffer.data(ksg1 + 239);
    const auto *ksg1_245 = buffer.data(ksg1 + 245);
    const auto *ksg1_249 = buffer.data(ksg1 + 249);
    const auto *ksg1_252 = buffer.data(ksg1 + 252);
    const auto *ksg1_253 = buffer.data(ksg1 + 253);
    const auto *ksg1_254 = buffer.data(ksg1 + 254);
    const auto *ksg1_255 = buffer.data(ksg1 + 255);

    const auto *ksh_272 = buffer.data(ksh + 272);
    const auto *ksh_273 = buffer.data(ksh + 273);
    const auto *ksh_275 = buffer.data(ksh + 275);
    const auto *ksh_276 = buffer.data(ksh + 276);
    const auto *ksh_278 = buffer.data(ksh + 278);
    const auto *ksh_279 = buffer.data(ksh + 279);
    const auto *ksh_282 = buffer.data(ksh + 282);
    const auto *ksh_288 = buffer.data(ksh + 288);
    const auto *ksh_289 = buffer.data(ksh + 289);
    const auto *ksh_290 = buffer.data(ksh + 290);
    const auto *ksh_291 = buffer.data(ksh + 291);
    const auto *ksh_292 = buffer.data(ksh + 292);
    const auto *ksh_293 = buffer.data(ksh + 293);
    const auto *ksh_294 = buffer.data(ksh + 294);
    const auto *ksh_295 = buffer.data(ksh + 295);
    const auto *ksh_296 = buffer.data(ksh + 296);
    const auto *ksh_297 = buffer.data(ksh + 297);
    const auto *ksh_298 = buffer.data(ksh + 298);
    const auto *ksh_299 = buffer.data(ksh + 299);
    const auto *ksh_300 = buffer.data(ksh + 300);
    const auto *ksh_301 = buffer.data(ksh + 301);
    const auto *ksh_302 = buffer.data(ksh + 302);
    const auto *ksh_303 = buffer.data(ksh + 303);
    const auto *ksh_308 = buffer.data(ksh + 308);
    const auto *ksh_309 = buffer.data(ksh + 309);
    const auto *ksh_310 = buffer.data(ksh + 310);
    const auto *ksh_311 = buffer.data(ksh + 311);
    const auto *ksh_312 = buffer.data(ksh + 312);
    const auto *ksh_313 = buffer.data(ksh + 313);
    const auto *ksh_314 = buffer.data(ksh + 314);
    const auto *ksh_315 = buffer.data(ksh + 315);
    const auto *ksh_316 = buffer.data(ksh + 316);
    const auto *ksh_317 = buffer.data(ksh + 317);
    const auto *ksh_318 = buffer.data(ksh + 318);
    const auto *ksh_320 = buffer.data(ksh + 320);
    const auto *ksh_321 = buffer.data(ksh + 321);
    const auto *ksh_322 = buffer.data(ksh + 322);
    const auto *ksh_324 = buffer.data(ksh + 324);
    const auto *ksh_325 = buffer.data(ksh + 325);
    const auto *ksh_330 = buffer.data(ksh + 330);
    const auto *ksh_331 = buffer.data(ksh + 331);
    const auto *ksh_332 = buffer.data(ksh + 332);
    const auto *ksh_333 = buffer.data(ksh + 333);
    const auto *ksh_334 = buffer.data(ksh + 334);
    const auto *ksh_335 = buffer.data(ksh + 335);
    const auto *ksh_336 = buffer.data(ksh + 336);
    const auto *ksh_338 = buffer.data(ksh + 338);
    const auto *ksh_339 = buffer.data(ksh + 339);
    const auto *ksh_341 = buffer.data(ksh + 341);
    const auto *ksh_342 = buffer.data(ksh + 342);
    const auto *ksh_345 = buffer.data(ksh + 345);
    const auto *ksh_350 = buffer.data(ksh + 350);
    const auto *ksh_351 = buffer.data(ksh + 351);
    const auto *ksh_352 = buffer.data(ksh + 352);
    const auto *ksh_353 = buffer.data(ksh + 353);
    const auto *ksh_354 = buffer.data(ksh + 354);
    const auto *ksh_355 = buffer.data(ksh + 355);
    const auto *ksh_356 = buffer.data(ksh + 356);
    const auto *ksh_357 = buffer.data(ksh + 357);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_y, pc_y, pc_z, isi0_252, ish_167, \
                         ish_188, ish_189, isi1_252, ksg0_194, ksg1_194, ksh_272, \
                         ksh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * ish_188[k]
                   + f_3 * pc_y[k] * ksh_272[k];

        t_363[k] = f_12 * ish_167[k]
                   + f_1 * ksg0_194[k]
                   - f_2 * ksg1_194[k]
                   + f_3 * pc_z[k] * ksh_272[k];

        t_364[k] = pa_y[k] * isi0_252[k]
                   - f_10 * pc_y[k] * isi1_252[k];

        t_365[k] = f_11 * ish_189[k]
                   + f_3 * pc_y[k] * ksh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pa_y, pc_y, pc_z, isi0_255, isi0_257, \
                         ish_168, ish_190, ish_191, isi1_255, isi1_257, ksh_273, \
                         ksh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * ish_168[k]
                   + f_3 * pc_z[k] * ksh_273[k];

        t_367[k] = pa_y[k] * isi0_255[k]
                   + f_12 * ish_190[k]
                   - f_10 * pc_y[k] * isi1_255[k];

        t_368[k] = f_11 * ish_191[k]
                   + f_3 * pc_y[k] * ksh_275[k];

        t_369[k] = pa_y[k] * isi0_257[k]
                   - f_10 * pc_y[k] * isi1_257[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pa_y, pc_y, pc_z, isi0_258, isi0_261, \
                         ish_171, ish_192, ish_194, isi1_258, isi1_261, ksh_276, \
                         ksh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_y[k] * isi0_258[k]
                   + f_13 * ish_192[k]
                   - f_10 * pc_y[k] * isi1_258[k];

        t_371[k] = f_13 * ish_171[k]
                   + f_3 * pc_z[k] * ksh_276[k];

        t_372[k] = f_11 * ish_194[k]
                   + f_3 * pc_y[k] * ksh_278[k];

        t_373[k] = pa_y[k] * isi0_261[k]
                   - f_10 * pc_y[k] * isi1_261[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pa_y, pc_y, pc_z, isi0_262, isi0_264, ish_174, \
                         ish_195, ish_197, isi1_262, isi1_264, \
                         ksh_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pa_y[k] * isi0_262[k]
                   + f_14 * ish_195[k]
                   - f_10 * pc_y[k] * isi1_262[k];

        t_375[k] = f_13 * ish_174[k]
                   + f_3 * pc_z[k] * ksh_279[k];

        t_376[k] = pa_y[k] * isi0_264[k]
                   + f_12 * ish_197[k]
                   - f_10 * pc_y[k] * isi1_264[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, t_380, pa_y, pc_x, pc_y, isi0_266, ish_198, \
                         ish_288, ish_289, isi1_266, ksh_282, ksh_288, \
                         ksh_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * ish_198[k]
                   + f_3 * pc_y[k] * ksh_282[k];

        t_378[k] = pa_y[k] * isi0_266[k]
                   - f_10 * pc_y[k] * isi1_266[k];

        t_379[k] = f_13 * ish_288[k]
                   + f_3 * pc_x[k] * ksh_288[k];

        t_380[k] = f_13 * ish_289[k]
                   + f_3 * pc_x[k] * ksh_289[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, ish_290, ish_291, ish_292, ish_293, \
                         ksh_290, ksh_291, ksh_292, ksh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_13 * ish_290[k]
                   + f_3 * pc_x[k] * ksh_290[k];

        t_382[k] = f_13 * ish_291[k]
                   + f_3 * pc_x[k] * ksh_291[k];

        t_383[k] = f_13 * ish_292[k]
                   + f_3 * pc_x[k] * ksh_292[k];

        t_384[k] = f_13 * ish_293[k]
                   + f_3 * pc_x[k] * ksh_293[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pc_y, pc_z, ish_183, ish_204, ish_206, ksg0_205, \
                         ksg0_207, ksg1_205, ksg1_207, ksh_288, \
                         ksh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * ish_204[k]
                   + f_1 * ksg0_205[k]
                   - f_2 * ksg1_205[k]
                   + f_3 * pc_y[k] * ksh_288[k];

        t_386[k] = f_13 * ish_183[k]
                   + f_3 * pc_z[k] * ksh_288[k];

        t_387[k] = f_11 * ish_206[k]
                   + f_8 * ksg0_207[k]
                   - f_9 * ksg1_207[k]
                   + f_3 * pc_y[k] * ksh_290[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_y, ish_207, ish_208, ish_209, ksg0_208, \
                         ksg0_209, ksg1_208, ksg1_209, ksh_291, ksh_292, \
                         ksh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_11 * ish_207[k]
                   + f_6 * ksg0_208[k]
                   - f_7 * ksg1_208[k]
                   + f_3 * pc_y[k] * ksh_291[k];

        t_389[k] = f_11 * ish_208[k]
                   + f_4 * ksg0_209[k]
                   - f_5 * ksg1_209[k]
                   + f_3 * pc_y[k] * ksh_292[k];

        t_390[k] = f_11 * ish_209[k]
                   + f_3 * pc_y[k] * ksh_293[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_x, pc_y, pc_z, isi0_279, \
                         ish_189, ish_294, isi1_279, ksg0_210, ksg1_210, \
                         ksh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = pa_y[k] * isi0_279[k]
                   - f_10 * pc_y[k] * isi1_279[k];

        t_392[k] = f_13 * ish_294[k]
                   + f_1 * ksg0_210[k]
                   - f_2 * ksg1_210[k]
                   + f_3 * pc_x[k] * ksh_294[k];

        t_393[k] = f_3 * pc_y[k] * ksh_294[k];

        t_394[k] = f_14 * ish_189[k]
                   + f_3 * pc_z[k] * ksh_294[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, pc_x, pc_y, ish_299, ksg0_210, ksg0_215, \
                         ksg1_210, ksg1_215, ksh_295, ksh_296, \
                         ksh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_4 * ksg0_210[k]
                   - f_5 * ksg1_210[k]
                   + f_3 * pc_y[k] * ksh_295[k];

        t_396[k] = f_3 * pc_y[k] * ksh_296[k];

        t_397[k] = f_13 * ish_299[k]
                   + f_8 * ksg0_215[k]
                   - f_9 * ksg1_215[k]
                   + f_3 * pc_x[k] * ksh_299[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, ksg0_211, ksg0_212, ksg1_211, ksg1_212, \
                         ksh_297, ksh_298, ksh_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_6 * ksg0_211[k]
                   - f_7 * ksg1_211[k]
                   + f_3 * pc_y[k] * ksh_297[k];

        t_399[k] = f_4 * ksg0_212[k]
                   - f_5 * ksg1_212[k]
                   + f_3 * pc_y[k] * ksh_298[k];

        t_400[k] = f_3 * pc_y[k] * ksh_299[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, ish_303, ksg0_213, ksg0_214, \
                         ksg0_219, ksg1_213, ksg1_214, ksg1_219, ksh_300, ksh_301, \
                         ksh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_13 * ish_303[k]
                   + f_6 * ksg0_219[k]
                   - f_7 * ksg1_219[k]
                   + f_3 * pc_x[k] * ksh_303[k];

        t_402[k] = f_8 * ksg0_213[k]
                   - f_9 * ksg1_213[k]
                   + f_3 * pc_y[k] * ksh_300[k];

        t_403[k] = f_6 * ksg0_214[k]
                   - f_7 * ksg1_214[k]
                   + f_3 * pc_y[k] * ksh_301[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pc_x, pc_y, ish_308, ish_309, ksg0_215, \
                         ksg0_224, ksg1_215, ksg1_224, ksh_302, ksh_303, ksh_308, \
                         ksh_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = f_4 * ksg0_215[k]
                   - f_5 * ksg1_215[k]
                   + f_3 * pc_y[k] * ksh_302[k];

        t_405[k] = f_3 * pc_y[k] * ksh_303[k];

        t_406[k] = f_13 * ish_308[k]
                   + f_4 * ksg0_224[k]
                   - f_5 * ksg1_224[k]
                   + f_3 * pc_x[k] * ksh_308[k];

        t_407[k] = f_13 * ish_309[k]
                   + f_3 * pc_x[k] * ksh_309[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, t_412, pc_x, pc_y, ish_310, ish_311, \
                         ish_312, ish_314, ksh_308, ksh_310, ksh_311, ksh_312, \
                         ksh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_13 * ish_310[k]
                   + f_3 * pc_x[k] * ksh_310[k];

        t_409[k] = f_13 * ish_311[k]
                   + f_3 * pc_x[k] * ksh_311[k];

        t_410[k] = f_13 * ish_312[k]
                   + f_3 * pc_x[k] * ksh_312[k];

        t_411[k] = f_3 * pc_y[k] * ksh_308[k];

        t_412[k] = f_13 * ish_314[k]
                   + f_3 * pc_x[k] * ksh_314[k];
    }

#pragma omp simd aligned(t_413, t_414, t_415, pc_y, ksg0_220, ksg0_221, ksg0_222, ksg1_220, \
                         ksg1_221, ksg1_222, ksh_309, ksh_310, \
                         ksh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_413[k] = f_1 * ksg0_220[k]
                   - f_2 * ksg1_220[k]
                   + f_3 * pc_y[k] * ksh_309[k];

        t_414[k] = f_16 * ksg0_221[k]
                   - f_17 * ksg1_221[k]
                   + f_3 * pc_y[k] * ksh_310[k];

        t_415[k] = f_8 * ksg0_222[k]
                   - f_9 * ksg1_222[k]
                   + f_3 * pc_y[k] * ksh_311[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pc_y, pc_z, ish_209, ksg0_223, ksg0_224, \
                         ksg1_223, ksg1_224, ksh_312, ksh_313, \
                         ksh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = f_6 * ksg0_223[k]
                   - f_7 * ksg1_223[k]
                   + f_3 * pc_y[k] * ksh_312[k];

        t_417[k] = f_4 * ksg0_224[k]
                   - f_5 * ksg1_224[k]
                   + f_3 * pc_y[k] * ksh_313[k];

        t_418[k] = f_3 * pc_y[k] * ksh_314[k];

        t_419[k] = f_14 * ish_209[k]
                   + f_1 * ksg0_224[k]
                   - f_2 * ksg1_224[k]
                   + f_3 * pc_z[k] * ksh_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pc_x, pc_y, pc_z, ish_210, ish_315, \
                         ish_318, ksg0_225, ksg0_228, ksg1_225, ksg1_228, ksh_315, \
                         ksh_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_12 * ish_315[k]
                   + f_1 * ksg0_225[k]
                   - f_2 * ksg1_225[k]
                   + f_3 * pc_x[k] * ksh_315[k];

        t_421[k] = f_18 * ish_210[k]
                   + f_3 * pc_y[k] * ksh_315[k];

        t_422[k] = f_3 * pc_z[k] * ksh_315[k];

        t_423[k] = f_12 * ish_318[k]
                   + f_8 * ksg0_228[k]
                   - f_9 * ksg1_228[k]
                   + f_3 * pc_x[k] * ksh_318[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, pc_x, pc_z, ish_321, ksg0_225, ksg0_231, \
                         ksg1_225, ksg1_231, ksh_316, ksh_317, ksh_318, \
                         ksh_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_3 * pc_z[k] * ksh_316[k];

        t_425[k] = f_4 * ksg0_225[k]
                   - f_5 * ksg1_225[k]
                   + f_3 * pc_z[k] * ksh_317[k];

        t_426[k] = f_12 * ish_321[k]
                   + f_6 * ksg0_231[k]
                   - f_7 * ksg1_231[k]
                   + f_3 * pc_x[k] * ksh_321[k];

        t_427[k] = f_3 * pc_z[k] * ksh_318[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, ish_215, ish_325, \
                         ksg0_227, ksg0_235, ksg1_227, ksg1_235, ksh_320, ksh_321, \
                         ksh_325 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_18 * ish_215[k]
                   + f_3 * pc_y[k] * ksh_320[k];

        t_429[k] = f_6 * ksg0_227[k]
                   - f_7 * ksg1_227[k]
                   + f_3 * pc_z[k] * ksh_320[k];

        t_430[k] = f_12 * ish_325[k]
                   + f_4 * ksg0_235[k]
                   - f_5 * ksg1_235[k]
                   + f_3 * pc_x[k] * ksh_325[k];

        t_431[k] = f_3 * pc_z[k] * ksh_321[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, pc_y, pc_z, ish_219, ish_330, \
                         ksg0_228, ksg0_230, ksg1_228, ksg1_230, ksh_322, ksh_324, \
                         ksh_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_4 * ksg0_228[k]
                   - f_5 * ksg1_228[k]
                   + f_3 * pc_z[k] * ksh_322[k];

        t_433[k] = f_18 * ish_219[k]
                   + f_3 * pc_y[k] * ksh_324[k];

        t_434[k] = f_8 * ksg0_230[k]
                   - f_9 * ksg1_230[k]
                   + f_3 * pc_z[k] * ksh_324[k];

        t_435[k] = f_12 * ish_330[k]
                   + f_3 * pc_x[k] * ksh_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, pc_z, ish_332, ish_333, \
                         ish_334, ish_335, ksh_325, ksh_332, ksh_333, ksh_334, \
                         ksh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_3 * pc_z[k] * ksh_325[k];

        t_437[k] = f_12 * ish_332[k]
                   + f_3 * pc_x[k] * ksh_332[k];

        t_438[k] = f_12 * ish_333[k]
                   + f_3 * pc_x[k] * ksh_333[k];

        t_439[k] = f_12 * ish_334[k]
                   + f_3 * pc_x[k] * ksh_334[k];

        t_440[k] = f_12 * ish_335[k]
                   + f_3 * pc_x[k] * ksh_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, ish_225, ksg0_235, ksg0_236, \
                         ksg1_235, ksg1_236, ksh_330, ksh_331, \
                         ksh_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_18 * ish_225[k]
                   + f_1 * ksg0_235[k]
                   - f_2 * ksg1_235[k]
                   + f_3 * pc_y[k] * ksh_330[k];

        t_442[k] = f_3 * pc_z[k] * ksh_330[k];

        t_443[k] = f_4 * ksg0_235[k]
                   - f_5 * ksg1_235[k]
                   + f_3 * pc_z[k] * ksh_331[k];

        t_444[k] = f_6 * ksg0_236[k]
                   - f_7 * ksg1_236[k]
                   + f_3 * pc_z[k] * ksh_332[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pa_z, pc_y, pc_z, isi0_280, ish_230, \
                         isi1_280, ksg0_237, ksg0_239, ksg1_237, ksg1_239, ksh_333, \
                         ksh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * ksg0_237[k]
                   - f_9 * ksg1_237[k]
                   + f_3 * pc_z[k] * ksh_333[k];

        t_446[k] = f_18 * ish_230[k]
                   + f_3 * pc_y[k] * ksh_335[k];

        t_447[k] = f_1 * ksg0_239[k]
                   - f_2 * ksg1_239[k]
                   + f_3 * pc_z[k] * ksh_335[k];

        t_448[k] = pa_z[k] * isi0_280[k]
                   - f_10 * pc_z[k] * isi1_280[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pa_z, pc_y, pc_z, isi0_283, ish_210, \
                         ish_231, ish_233, isi1_283, ksh_336, ksh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_14 * ish_231[k]
                   + f_3 * pc_y[k] * ksh_336[k];

        t_450[k] = f_11 * ish_210[k]
                   + f_3 * pc_z[k] * ksh_336[k];

        t_451[k] = pa_z[k] * isi0_283[k]
                   - f_10 * pc_z[k] * isi1_283[k];

        t_452[k] = f_14 * ish_233[k]
                   + f_3 * pc_y[k] * ksh_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pa_z, pc_x, pc_z, isi0_286, ish_213, ish_341, \
                         isi1_286, ksg0_245, ksg1_245, ksh_339, \
                         ksh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_12 * ish_341[k]
                   + f_8 * ksg0_245[k]
                   - f_9 * ksg1_245[k]
                   + f_3 * pc_x[k] * ksh_341[k];

        t_454[k] = pa_z[k] * isi0_286[k]
                   - f_10 * pc_z[k] * isi1_286[k];

        t_455[k] = f_11 * ish_213[k]
                   + f_3 * pc_z[k] * ksh_339[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pa_z, pc_x, pc_y, pc_z, isi0_290, ish_236, \
                         ish_345, isi1_290, ksg0_249, ksg1_249, ksh_341, \
                         ksh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_14 * ish_236[k]
                   + f_3 * pc_y[k] * ksh_341[k];

        t_457[k] = f_12 * ish_345[k]
                   + f_6 * ksg0_249[k]
                   - f_7 * ksg1_249[k]
                   + f_3 * pc_x[k] * ksh_345[k];

        t_458[k] = pa_z[k] * isi0_290[k]
                   - f_10 * pc_z[k] * isi1_290[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pa_z, pc_y, pc_z, isi0_292, ish_216, ish_217, \
                         ish_240, isi1_292, ksh_342, ksh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_11 * ish_216[k]
                   + f_3 * pc_z[k] * ksh_342[k];

        t_460[k] = pa_z[k] * isi0_292[k]
                   + f_12 * ish_217[k]
                   - f_10 * pc_z[k] * isi1_292[k];

        t_461[k] = f_14 * ish_240[k]
                   + f_3 * pc_y[k] * ksh_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, ish_350, ish_351, ish_352, ish_353, \
                         ksg0_254, ksg1_254, ksh_350, ksh_351, ksh_352, \
                         ksh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_12 * ish_350[k]
                   + f_4 * ksg0_254[k]
                   - f_5 * ksg1_254[k]
                   + f_3 * pc_x[k] * ksh_350[k];

        t_463[k] = f_12 * ish_351[k]
                   + f_3 * pc_x[k] * ksh_351[k];

        t_464[k] = f_12 * ish_352[k]
                   + f_3 * pc_x[k] * ksh_352[k];

        t_465[k] = f_12 * ish_353[k]
                   + f_3 * pc_x[k] * ksh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_x, pc_z, isi0_301, ish_354, \
                         ish_355, ish_356, isi1_301, ksh_354, ksh_355, \
                         ksh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_12 * ish_354[k]
                   + f_3 * pc_x[k] * ksh_354[k];

        t_467[k] = f_12 * ish_355[k]
                   + f_3 * pc_x[k] * ksh_355[k];

        t_468[k] = f_12 * ish_356[k]
                   + f_3 * pc_x[k] * ksh_356[k];

        t_469[k] = pa_z[k] * isi0_301[k]
                   - f_10 * pc_z[k] * isi1_301[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pc_y, pc_z, ish_225, ish_248, ish_249, ksg0_252, \
                         ksg0_253, ksg1_252, ksg1_253, ksh_351, ksh_353, \
                         ksh_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_11 * ish_225[k]
                   + f_3 * pc_z[k] * ksh_351[k];

        t_471[k] = f_14 * ish_248[k]
                   + f_8 * ksg0_252[k]
                   - f_9 * ksg1_252[k]
                   + f_3 * pc_y[k] * ksh_353[k];

        t_472[k] = f_14 * ish_249[k]
                   + f_6 * ksg0_253[k]
                   - f_7 * ksg1_253[k]
                   + f_3 * pc_y[k] * ksh_354[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, pc_y, pc_z, ish_230, ish_250, ish_251, ksg0_254, \
                         ksg1_254, ksh_355, ksh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_14 * ish_250[k]
                   + f_4 * ksg0_254[k]
                   - f_5 * ksg1_254[k]
                   + f_3 * pc_y[k] * ksh_355[k];

        t_474[k] = f_14 * ish_251[k]
                   + f_3 * pc_y[k] * ksh_356[k];

        t_475[k] = f_11 * ish_230[k]
                   + f_1 * ksg0_254[k]
                   - f_2 * ksg1_254[k]
                   + f_3 * pc_z[k] * ksh_356[k];
    }

#pragma omp simd aligned(t_476, t_477, t_478, pc_x, pc_y, pc_z, ish_231, ish_252, ish_357, \
                         ksg0_255, ksg1_255, ksh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_476[k] = f_12 * ish_357[k]
                   + f_1 * ksg0_255[k]
                   - f_2 * ksg1_255[k]
                   + f_3 * pc_x[k] * ksh_357[k];

        t_477[k] = f_13 * ish_252[k]
                   + f_3 * pc_y[k] * ksh_357[k];

        t_478[k] = f_12 * ish_231[k]
                   + f_3 * pc_z[k] * ksh_357[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isi0_392 = buffer.data(isi0 + 392);
    const auto *isi0_395 = buffer.data(isi0 + 395);
    const auto *isi0_397 = buffer.data(isi0 + 397);
    const auto *isi0_398 = buffer.data(isi0 + 398);
    const auto *isi0_401 = buffer.data(isi0 + 401);
    const auto *isi0_402 = buffer.data(isi0 + 402);
    const auto *isi0_404 = buffer.data(isi0 + 404);
    const auto *isi0_406 = buffer.data(isi0 + 406);
    const auto *isi0_419 = buffer.data(isi0 + 419);
    const auto *isi0_588 = buffer.data(isi0 + 588);
    const auto *isi0_591 = buffer.data(isi0 + 591);

    const auto *ish_234 = buffer.data(ish + 234);
    const auto *ish_237 = buffer.data(ish + 237);
    const auto *ish_246 = buffer.data(ish + 246);
    const auto *ish_251 = buffer.data(ish + 251);
    const auto *ish_252 = buffer.data(ish + 252);
    const auto *ish_254 = buffer.data(ish + 254);
    const auto *ish_255 = buffer.data(ish + 255);
    const auto *ish_257 = buffer.data(ish + 257);
    const auto *ish_258 = buffer.data(ish + 258);
    const auto *ish_261 = buffer.data(ish + 261);
    const auto *ish_267 = buffer.data(ish + 267);
    const auto *ish_269 = buffer.data(ish + 269);
    const auto *ish_270 = buffer.data(ish + 270);
    const auto *ish_271 = buffer.data(ish + 271);
    const auto *ish_272 = buffer.data(ish + 272);
    const auto *ish_273 = buffer.data(ish + 273);
    const auto *ish_275 = buffer.data(ish + 275);
    const auto *ish_276 = buffer.data(ish + 276);
    const auto *ish_278 = buffer.data(ish + 278);
    const auto *ish_279 = buffer.data(ish + 279);
    const auto *ish_282 = buffer.data(ish + 282);
    const auto *ish_288 = buffer.data(ish + 288);
    const auto *ish_290 = buffer.data(ish + 290);
    const auto *ish_291 = buffer.data(ish + 291);
    const auto *ish_292 = buffer.data(ish + 292);
    const auto *ish_293 = buffer.data(ish + 293);
    const auto *ish_294 = buffer.data(ish + 294);
    const auto *ish_295 = buffer.data(ish + 295);
    const auto *ish_296 = buffer.data(ish + 296);
    const auto *ish_297 = buffer.data(ish + 297);
    const auto *ish_299 = buffer.data(ish + 299);
    const auto *ish_300 = buffer.data(ish + 300);
    const auto *ish_302 = buffer.data(ish + 302);
    const auto *ish_303 = buffer.data(ish + 303);
    const auto *ish_309 = buffer.data(ish + 309);
    const auto *ish_311 = buffer.data(ish + 311);
    const auto *ish_312 = buffer.data(ish + 312);
    const auto *ish_313 = buffer.data(ish + 313);
    const auto *ish_314 = buffer.data(ish + 314);
    const auto *ish_315 = buffer.data(ish + 315);
    const auto *ish_360 = buffer.data(ish + 360);
    const auto *ish_362 = buffer.data(ish + 362);
    const auto *ish_363 = buffer.data(ish + 363);
    const auto *ish_366 = buffer.data(ish + 366);
    const auto *ish_367 = buffer.data(ish + 367);
    const auto *ish_369 = buffer.data(ish + 369);
    const auto *ish_371 = buffer.data(ish + 371);
    const auto *ish_372 = buffer.data(ish + 372);
    const auto *ish_373 = buffer.data(ish + 373);
    const auto *ish_374 = buffer.data(ish + 374);
    const auto *ish_375 = buffer.data(ish + 375);
    const auto *ish_376 = buffer.data(ish + 376);
    const auto *ish_377 = buffer.data(ish + 377);
    const auto *ish_378 = buffer.data(ish + 378);
    const auto *ish_381 = buffer.data(ish + 381);
    const auto *ish_383 = buffer.data(ish + 383);
    const auto *ish_384 = buffer.data(ish + 384);
    const auto *ish_387 = buffer.data(ish + 387);
    const auto *ish_388 = buffer.data(ish + 388);
    const auto *ish_390 = buffer.data(ish + 390);
    const auto *ish_392 = buffer.data(ish + 392);
    const auto *ish_393 = buffer.data(ish + 393);
    const auto *ish_394 = buffer.data(ish + 394);
    const auto *ish_395 = buffer.data(ish + 395);
    const auto *ish_396 = buffer.data(ish + 396);
    const auto *ish_397 = buffer.data(ish + 397);
    const auto *ish_398 = buffer.data(ish + 398);
    const auto *ish_414 = buffer.data(ish + 414);
    const auto *ish_415 = buffer.data(ish + 415);
    const auto *ish_416 = buffer.data(ish + 416);
    const auto *ish_417 = buffer.data(ish + 417);
    const auto *ish_418 = buffer.data(ish + 418);
    const auto *ish_419 = buffer.data(ish + 419);
    const auto *ish_420 = buffer.data(ish + 420);
    const auto *ish_425 = buffer.data(ish + 425);
    const auto *ish_429 = buffer.data(ish + 429);
    const auto *ish_434 = buffer.data(ish + 434);
    const auto *ish_435 = buffer.data(ish + 435);
    const auto *ish_436 = buffer.data(ish + 436);
    const auto *ish_437 = buffer.data(ish + 437);
    const auto *ish_438 = buffer.data(ish + 438);
    const auto *ish_440 = buffer.data(ish + 440);
    const auto *ish_441 = buffer.data(ish + 441);
    const auto *ish_444 = buffer.data(ish + 444);

    const auto *isi1_392 = buffer.data(isi1 + 392);
    const auto *isi1_395 = buffer.data(isi1 + 395);
    const auto *isi1_397 = buffer.data(isi1 + 397);
    const auto *isi1_398 = buffer.data(isi1 + 398);
    const auto *isi1_401 = buffer.data(isi1 + 401);
    const auto *isi1_402 = buffer.data(isi1 + 402);
    const auto *isi1_404 = buffer.data(isi1 + 404);
    const auto *isi1_406 = buffer.data(isi1 + 406);
    const auto *isi1_419 = buffer.data(isi1 + 419);
    const auto *isi1_588 = buffer.data(isi1 + 588);
    const auto *isi1_591 = buffer.data(isi1 + 591);

    const auto *ksg0_258 = buffer.data(ksg0 + 258);
    const auto *ksg0_260 = buffer.data(ksg0 + 260);
    const auto *ksg0_261 = buffer.data(ksg0 + 261);
    const auto *ksg0_264 = buffer.data(ksg0 + 264);
    const auto *ksg0_265 = buffer.data(ksg0 + 265);
    const auto *ksg0_267 = buffer.data(ksg0 + 267);
    const auto *ksg0_268 = buffer.data(ksg0 + 268);
    const auto *ksg0_269 = buffer.data(ksg0 + 269);
    const auto *ksg0_270 = buffer.data(ksg0 + 270);
    const auto *ksg0_273 = buffer.data(ksg0 + 273);
    const auto *ksg0_275 = buffer.data(ksg0 + 275);
    const auto *ksg0_276 = buffer.data(ksg0 + 276);
    const auto *ksg0_279 = buffer.data(ksg0 + 279);
    const auto *ksg0_280 = buffer.data(ksg0 + 280);
    const auto *ksg0_282 = buffer.data(ksg0 + 282);
    const auto *ksg0_283 = buffer.data(ksg0 + 283);
    const auto *ksg0_284 = buffer.data(ksg0 + 284);
    const auto *ksg0_295 = buffer.data(ksg0 + 295);
    const auto *ksg0_297 = buffer.data(ksg0 + 297);
    const auto *ksg0_298 = buffer.data(ksg0 + 298);
    const auto *ksg0_299 = buffer.data(ksg0 + 299);
    const auto *ksg0_300 = buffer.data(ksg0 + 300);
    const auto *ksg0_301 = buffer.data(ksg0 + 301);
    const auto *ksg0_302 = buffer.data(ksg0 + 302);
    const auto *ksg0_303 = buffer.data(ksg0 + 303);
    const auto *ksg0_304 = buffer.data(ksg0 + 304);
    const auto *ksg0_305 = buffer.data(ksg0 + 305);
    const auto *ksg0_309 = buffer.data(ksg0 + 309);
    const auto *ksg0_310 = buffer.data(ksg0 + 310);
    const auto *ksg0_311 = buffer.data(ksg0 + 311);
    const auto *ksg0_312 = buffer.data(ksg0 + 312);
    const auto *ksg0_313 = buffer.data(ksg0 + 313);
    const auto *ksg0_314 = buffer.data(ksg0 + 314);

    const auto *ksg1_258 = buffer.data(ksg1 + 258);
    const auto *ksg1_260 = buffer.data(ksg1 + 260);
    const auto *ksg1_261 = buffer.data(ksg1 + 261);
    const auto *ksg1_264 = buffer.data(ksg1 + 264);
    const auto *ksg1_265 = buffer.data(ksg1 + 265);
    const auto *ksg1_267 = buffer.data(ksg1 + 267);
    const auto *ksg1_268 = buffer.data(ksg1 + 268);
    const auto *ksg1_269 = buffer.data(ksg1 + 269);
    const auto *ksg1_270 = buffer.data(ksg1 + 270);
    const auto *ksg1_273 = buffer.data(ksg1 + 273);
    const auto *ksg1_275 = buffer.data(ksg1 + 275);
    const auto *ksg1_276 = buffer.data(ksg1 + 276);
    const auto *ksg1_279 = buffer.data(ksg1 + 279);
    const auto *ksg1_280 = buffer.data(ksg1 + 280);
    const auto *ksg1_282 = buffer.data(ksg1 + 282);
    const auto *ksg1_283 = buffer.data(ksg1 + 283);
    const auto *ksg1_284 = buffer.data(ksg1 + 284);
    const auto *ksg1_295 = buffer.data(ksg1 + 295);
    const auto *ksg1_297 = buffer.data(ksg1 + 297);
    const auto *ksg1_298 = buffer.data(ksg1 + 298);
    const auto *ksg1_299 = buffer.data(ksg1 + 299);
    const auto *ksg1_300 = buffer.data(ksg1 + 300);
    const auto *ksg1_301 = buffer.data(ksg1 + 301);
    const auto *ksg1_302 = buffer.data(ksg1 + 302);
    const auto *ksg1_303 = buffer.data(ksg1 + 303);
    const auto *ksg1_304 = buffer.data(ksg1 + 304);
    const auto *ksg1_305 = buffer.data(ksg1 + 305);
    const auto *ksg1_309 = buffer.data(ksg1 + 309);
    const auto *ksg1_310 = buffer.data(ksg1 + 310);
    const auto *ksg1_311 = buffer.data(ksg1 + 311);
    const auto *ksg1_312 = buffer.data(ksg1 + 312);
    const auto *ksg1_313 = buffer.data(ksg1 + 313);
    const auto *ksg1_314 = buffer.data(ksg1 + 314);

    const auto *ksh_359 = buffer.data(ksh + 359);
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
    const auto *ksh_380 = buffer.data(ksh + 380);
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
    const auto *ksh_399 = buffer.data(ksh + 399);
    const auto *ksh_401 = buffer.data(ksh + 401);
    const auto *ksh_402 = buffer.data(ksh + 402);
    const auto *ksh_404 = buffer.data(ksh + 404);
    const auto *ksh_405 = buffer.data(ksh + 405);
    const auto *ksh_408 = buffer.data(ksh + 408);
    const auto *ksh_414 = buffer.data(ksh + 414);
    const auto *ksh_415 = buffer.data(ksh + 415);
    const auto *ksh_416 = buffer.data(ksh + 416);
    const auto *ksh_417 = buffer.data(ksh + 417);
    const auto *ksh_418 = buffer.data(ksh + 418);
    const auto *ksh_419 = buffer.data(ksh + 419);
    const auto *ksh_420 = buffer.data(ksh + 420);
    const auto *ksh_421 = buffer.data(ksh + 421);
    const auto *ksh_422 = buffer.data(ksh + 422);
    const auto *ksh_423 = buffer.data(ksh + 423);
    const auto *ksh_424 = buffer.data(ksh + 424);
    const auto *ksh_425 = buffer.data(ksh + 425);
    const auto *ksh_426 = buffer.data(ksh + 426);
    const auto *ksh_427 = buffer.data(ksh + 427);
    const auto *ksh_428 = buffer.data(ksh + 428);
    const auto *ksh_429 = buffer.data(ksh + 429);
    const auto *ksh_434 = buffer.data(ksh + 434);
    const auto *ksh_435 = buffer.data(ksh + 435);
    const auto *ksh_436 = buffer.data(ksh + 436);
    const auto *ksh_437 = buffer.data(ksh + 437);
    const auto *ksh_438 = buffer.data(ksh + 438);
    const auto *ksh_439 = buffer.data(ksh + 439);
    const auto *ksh_440 = buffer.data(ksh + 440);
    const auto *ksh_441 = buffer.data(ksh + 441);

#pragma omp simd aligned(t_479, t_480, t_481, pc_x, pc_y, ish_254, ish_360, ish_362, ksg0_258, \
                         ksg0_260, ksg1_258, ksg1_260, ksh_359, ksh_360, \
                         ksh_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_12 * ish_360[k]
                   + f_8 * ksg0_258[k]
                   - f_9 * ksg1_258[k]
                   + f_3 * pc_x[k] * ksh_360[k];

        t_480[k] = f_13 * ish_254[k]
                   + f_3 * pc_y[k] * ksh_359[k];

        t_481[k] = f_12 * ish_362[k]
                   + f_8 * ksg0_260[k]
                   - f_9 * ksg1_260[k]
                   + f_3 * pc_x[k] * ksh_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, pc_x, pc_y, pc_z, ish_234, ish_257, ish_363, \
                         ksg0_261, ksg1_261, ksh_360, ksh_362, \
                         ksh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_12 * ish_363[k]
                   + f_6 * ksg0_261[k]
                   - f_7 * ksg1_261[k]
                   + f_3 * pc_x[k] * ksh_363[k];

        t_483[k] = f_12 * ish_234[k]
                   + f_3 * pc_z[k] * ksh_360[k];

        t_484[k] = f_13 * ish_257[k]
                   + f_3 * pc_y[k] * ksh_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_z, ish_237, ish_366, ish_367, ksg0_264, \
                         ksg0_265, ksg1_264, ksg1_265, ksh_363, ksh_366, \
                         ksh_367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_12 * ish_366[k]
                   + f_6 * ksg0_264[k]
                   - f_7 * ksg1_264[k]
                   + f_3 * pc_x[k] * ksh_366[k];

        t_486[k] = f_12 * ish_367[k]
                   + f_4 * ksg0_265[k]
                   - f_5 * ksg1_265[k]
                   + f_3 * pc_x[k] * ksh_367[k];

        t_487[k] = f_12 * ish_237[k]
                   + f_3 * pc_z[k] * ksh_363[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pc_x, pc_y, ish_261, ish_369, ish_371, ksg0_267, \
                         ksg0_269, ksg1_267, ksg1_269, ksh_366, ksh_369, \
                         ksh_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_12 * ish_369[k]
                   + f_4 * ksg0_267[k]
                   - f_5 * ksg1_267[k]
                   + f_3 * pc_x[k] * ksh_369[k];

        t_489[k] = f_13 * ish_261[k]
                   + f_3 * pc_y[k] * ksh_366[k];

        t_490[k] = f_12 * ish_371[k]
                   + f_4 * ksg0_269[k]
                   - f_5 * ksg1_269[k]
                   + f_3 * pc_x[k] * ksh_371[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, t_495, pc_x, ish_372, ish_373, ish_374, \
                         ish_375, ish_376, ksh_372, ksh_373, ksh_374, ksh_375, \
                         ksh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_12 * ish_372[k]
                   + f_3 * pc_x[k] * ksh_372[k];

        t_492[k] = f_12 * ish_373[k]
                   + f_3 * pc_x[k] * ksh_373[k];

        t_493[k] = f_12 * ish_374[k]
                   + f_3 * pc_x[k] * ksh_374[k];

        t_494[k] = f_12 * ish_375[k]
                   + f_3 * pc_x[k] * ksh_375[k];

        t_495[k] = f_12 * ish_376[k]
                   + f_3 * pc_x[k] * ksh_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, pc_x, pc_y, pc_z, ish_246, ish_267, ish_377, \
                         ksg0_265, ksg1_265, ksh_372, ksh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_12 * ish_377[k]
                   + f_3 * pc_x[k] * ksh_377[k];

        t_497[k] = f_13 * ish_267[k]
                   + f_1 * ksg0_265[k]
                   - f_2 * ksg1_265[k]
                   + f_3 * pc_y[k] * ksh_372[k];

        t_498[k] = f_12 * ish_246[k]
                   + f_3 * pc_z[k] * ksh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, ish_269, ish_270, ish_271, ksg0_267, \
                         ksg0_268, ksg0_269, ksg1_267, ksg1_268, ksg1_269, ksh_374, ksh_375, \
                         ksh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * ish_269[k]
                   + f_8 * ksg0_267[k]
                   - f_9 * ksg1_267[k]
                   + f_3 * pc_y[k] * ksh_374[k];

        t_500[k] = f_13 * ish_270[k]
                   + f_6 * ksg0_268[k]
                   - f_7 * ksg1_268[k]
                   + f_3 * pc_y[k] * ksh_375[k];

        t_501[k] = f_13 * ish_271[k]
                   + f_4 * ksg0_269[k]
                   - f_5 * ksg1_269[k]
                   + f_3 * pc_y[k] * ksh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, ish_251, ish_272, ish_378, \
                         ksg0_269, ksg0_270, ksg1_269, ksg1_270, ksh_377, \
                         ksh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * ish_272[k]
                   + f_3 * pc_y[k] * ksh_377[k];

        t_503[k] = f_12 * ish_251[k]
                   + f_1 * ksg0_269[k]
                   - f_2 * ksg1_269[k]
                   + f_3 * pc_z[k] * ksh_377[k];

        t_504[k] = f_12 * ish_378[k]
                   + f_1 * ksg0_270[k]
                   - f_2 * ksg1_270[k]
                   + f_3 * pc_x[k] * ksh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pc_x, pc_y, pc_z, ish_252, ish_273, \
                         ish_275, ish_381, ksg0_273, ksg1_273, ksh_378, ksh_380, \
                         ksh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_12 * ish_273[k]
                   + f_3 * pc_y[k] * ksh_378[k];

        t_506[k] = f_13 * ish_252[k]
                   + f_3 * pc_z[k] * ksh_378[k];

        t_507[k] = f_12 * ish_381[k]
                   + f_8 * ksg0_273[k]
                   - f_9 * ksg1_273[k]
                   + f_3 * pc_x[k] * ksh_381[k];

        t_508[k] = f_12 * ish_275[k]
                   + f_3 * pc_y[k] * ksh_380[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pc_x, pc_z, ish_255, ish_383, ish_384, ksg0_275, \
                         ksg0_276, ksg1_275, ksg1_276, ksh_381, ksh_383, \
                         ksh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_12 * ish_383[k]
                   + f_8 * ksg0_275[k]
                   - f_9 * ksg1_275[k]
                   + f_3 * pc_x[k] * ksh_383[k];

        t_510[k] = f_12 * ish_384[k]
                   + f_6 * ksg0_276[k]
                   - f_7 * ksg1_276[k]
                   + f_3 * pc_x[k] * ksh_384[k];

        t_511[k] = f_13 * ish_255[k]
                   + f_3 * pc_z[k] * ksh_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pc_x, pc_y, ish_278, ish_387, ish_388, ksg0_279, \
                         ksg0_280, ksg1_279, ksg1_280, ksh_383, ksh_387, \
                         ksh_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_12 * ish_278[k]
                   + f_3 * pc_y[k] * ksh_383[k];

        t_513[k] = f_12 * ish_387[k]
                   + f_6 * ksg0_279[k]
                   - f_7 * ksg1_279[k]
                   + f_3 * pc_x[k] * ksh_387[k];

        t_514[k] = f_12 * ish_388[k]
                   + f_4 * ksg0_280[k]
                   - f_5 * ksg1_280[k]
                   + f_3 * pc_x[k] * ksh_388[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pc_x, pc_y, pc_z, ish_258, ish_282, ish_390, \
                         ksg0_282, ksg1_282, ksh_384, ksh_387, \
                         ksh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_13 * ish_258[k]
                   + f_3 * pc_z[k] * ksh_384[k];

        t_516[k] = f_12 * ish_390[k]
                   + f_4 * ksg0_282[k]
                   - f_5 * ksg1_282[k]
                   + f_3 * pc_x[k] * ksh_390[k];

        t_517[k] = f_12 * ish_282[k]
                   + f_3 * pc_y[k] * ksh_387[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pc_x, ish_392, ish_393, ish_394, ish_395, \
                         ksg0_284, ksg1_284, ksh_392, ksh_393, ksh_394, \
                         ksh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = f_12 * ish_392[k]
                   + f_4 * ksg0_284[k]
                   - f_5 * ksg1_284[k]
                   + f_3 * pc_x[k] * ksh_392[k];

        t_519[k] = f_12 * ish_393[k]
                   + f_3 * pc_x[k] * ksh_393[k];

        t_520[k] = f_12 * ish_394[k]
                   + f_3 * pc_x[k] * ksh_394[k];

        t_521[k] = f_12 * ish_395[k]
                   + f_3 * pc_x[k] * ksh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_y, ish_288, ish_396, ish_397, \
                         ish_398, ksg0_280, ksg1_280, ksh_393, ksh_396, ksh_397, \
                         ksh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_12 * ish_396[k]
                   + f_3 * pc_x[k] * ksh_396[k];

        t_523[k] = f_12 * ish_397[k]
                   + f_3 * pc_x[k] * ksh_397[k];

        t_524[k] = f_12 * ish_398[k]
                   + f_3 * pc_x[k] * ksh_398[k];

        t_525[k] = f_12 * ish_288[k]
                   + f_1 * ksg0_280[k]
                   - f_2 * ksg1_280[k]
                   + f_3 * pc_y[k] * ksh_393[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_y, pc_z, ish_267, ish_290, ish_291, ksg0_282, \
                         ksg0_283, ksg1_282, ksg1_283, ksh_393, ksh_395, \
                         ksh_396 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_13 * ish_267[k]
                   + f_3 * pc_z[k] * ksh_393[k];

        t_527[k] = f_12 * ish_290[k]
                   + f_8 * ksg0_282[k]
                   - f_9 * ksg1_282[k]
                   + f_3 * pc_y[k] * ksh_395[k];

        t_528[k] = f_12 * ish_291[k]
                   + f_6 * ksg0_283[k]
                   - f_7 * ksg1_283[k]
                   + f_3 * pc_y[k] * ksh_396[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, pa_y, pc_y, pc_z, isi0_392, ish_272, \
                         ish_292, ish_293, isi1_392, ksg0_284, ksg1_284, ksh_397, \
                         ksh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_12 * ish_292[k]
                   + f_4 * ksg0_284[k]
                   - f_5 * ksg1_284[k]
                   + f_3 * pc_y[k] * ksh_397[k];

        t_530[k] = f_12 * ish_293[k]
                   + f_3 * pc_y[k] * ksh_398[k];

        t_531[k] = f_13 * ish_272[k]
                   + f_1 * ksg0_284[k]
                   - f_2 * ksg1_284[k]
                   + f_3 * pc_z[k] * ksh_398[k];

        t_532[k] = pa_y[k] * isi0_392[k]
                   - f_10 * pc_y[k] * isi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pa_y, pc_y, pc_z, isi0_395, ish_273, \
                         ish_294, ish_295, ish_296, isi1_395, ksh_399, \
                         ksh_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_11 * ish_294[k]
                   + f_3 * pc_y[k] * ksh_399[k];

        t_534[k] = f_14 * ish_273[k]
                   + f_3 * pc_z[k] * ksh_399[k];

        t_535[k] = pa_y[k] * isi0_395[k]
                   + f_12 * ish_295[k]
                   - f_10 * pc_y[k] * isi1_395[k];

        t_536[k] = f_11 * ish_296[k]
                   + f_3 * pc_y[k] * ksh_401[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pa_y, pc_y, pc_z, isi0_397, isi0_398, \
                         ish_276, ish_297, ish_299, isi1_397, isi1_398, ksh_402, \
                         ksh_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = pa_y[k] * isi0_397[k]
                   - f_10 * pc_y[k] * isi1_397[k];

        t_538[k] = pa_y[k] * isi0_398[k]
                   + f_13 * ish_297[k]
                   - f_10 * pc_y[k] * isi1_398[k];

        t_539[k] = f_14 * ish_276[k]
                   + f_3 * pc_z[k] * ksh_402[k];

        t_540[k] = f_11 * ish_299[k]
                   + f_3 * pc_y[k] * ksh_404[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, pa_y, pc_y, pc_z, isi0_401, isi0_402, ish_279, \
                         ish_300, isi1_401, isi1_402, ksh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = pa_y[k] * isi0_401[k]
                   - f_10 * pc_y[k] * isi1_401[k];

        t_542[k] = pa_y[k] * isi0_402[k]
                   + f_14 * ish_300[k]
                   - f_10 * pc_y[k] * isi1_402[k];

        t_543[k] = f_14 * ish_279[k]
                   + f_3 * pc_z[k] * ksh_405[k];
    }

#pragma omp simd aligned(t_544, t_545, t_546, t_547, pa_y, pc_x, pc_y, isi0_404, isi0_406, \
                         ish_302, ish_303, ish_414, isi1_404, isi1_406, ksh_408, \
                         ksh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_544[k] = pa_y[k] * isi0_404[k]
                   + f_12 * ish_302[k]
                   - f_10 * pc_y[k] * isi1_404[k];

        t_545[k] = f_11 * ish_303[k]
                   + f_3 * pc_y[k] * ksh_408[k];

        t_546[k] = pa_y[k] * isi0_406[k]
                   - f_10 * pc_y[k] * isi1_406[k];

        t_547[k] = f_12 * ish_414[k]
                   + f_3 * pc_x[k] * ksh_414[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, t_551, t_552, pc_x, ish_415, ish_416, ish_417, \
                         ish_418, ish_419, ksh_415, ksh_416, ksh_417, ksh_418, \
                         ksh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_12 * ish_415[k]
                   + f_3 * pc_x[k] * ksh_415[k];

        t_549[k] = f_12 * ish_416[k]
                   + f_3 * pc_x[k] * ksh_416[k];

        t_550[k] = f_12 * ish_417[k]
                   + f_3 * pc_x[k] * ksh_417[k];

        t_551[k] = f_12 * ish_418[k]
                   + f_3 * pc_x[k] * ksh_418[k];

        t_552[k] = f_12 * ish_419[k]
                   + f_3 * pc_x[k] * ksh_419[k];
    }

#pragma omp simd aligned(t_553, t_554, t_555, pc_y, pc_z, ish_288, ish_309, ish_311, ksg0_295, \
                         ksg0_297, ksg1_295, ksg1_297, ksh_414, \
                         ksh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_553[k] = f_11 * ish_309[k]
                   + f_1 * ksg0_295[k]
                   - f_2 * ksg1_295[k]
                   + f_3 * pc_y[k] * ksh_414[k];

        t_554[k] = f_14 * ish_288[k]
                   + f_3 * pc_z[k] * ksh_414[k];

        t_555[k] = f_11 * ish_311[k]
                   + f_8 * ksg0_297[k]
                   - f_9 * ksg1_297[k]
                   + f_3 * pc_y[k] * ksh_416[k];
    }

#pragma omp simd aligned(t_556, t_557, t_558, pc_y, ish_312, ish_313, ish_314, ksg0_298, \
                         ksg0_299, ksg1_298, ksg1_299, ksh_417, ksh_418, \
                         ksh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_556[k] = f_11 * ish_312[k]
                   + f_6 * ksg0_298[k]
                   - f_7 * ksg1_298[k]
                   + f_3 * pc_y[k] * ksh_417[k];

        t_557[k] = f_11 * ish_313[k]
                   + f_4 * ksg0_299[k]
                   - f_5 * ksg1_299[k]
                   + f_3 * pc_y[k] * ksh_418[k];

        t_558[k] = f_11 * ish_314[k]
                   + f_3 * pc_y[k] * ksh_419[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_y, pc_x, pc_y, pc_z, isi0_419, \
                         ish_294, ish_420, isi1_419, ksg0_300, ksg1_300, \
                         ksh_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_y[k] * isi0_419[k]
                   - f_10 * pc_y[k] * isi1_419[k];

        t_560[k] = f_12 * ish_420[k]
                   + f_1 * ksg0_300[k]
                   - f_2 * ksg1_300[k]
                   + f_3 * pc_x[k] * ksh_420[k];

        t_561[k] = f_3 * pc_y[k] * ksh_420[k];

        t_562[k] = f_18 * ish_294[k]
                   + f_3 * pc_z[k] * ksh_420[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, pc_y, ish_425, ksg0_300, ksg0_305, \
                         ksg1_300, ksg1_305, ksh_421, ksh_422, \
                         ksh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_4 * ksg0_300[k]
                   - f_5 * ksg1_300[k]
                   + f_3 * pc_y[k] * ksh_421[k];

        t_564[k] = f_3 * pc_y[k] * ksh_422[k];

        t_565[k] = f_12 * ish_425[k]
                   + f_8 * ksg0_305[k]
                   - f_9 * ksg1_305[k]
                   + f_3 * pc_x[k] * ksh_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pc_y, ksg0_301, ksg0_302, ksg1_301, ksg1_302, \
                         ksh_423, ksh_424, ksh_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_6 * ksg0_301[k]
                   - f_7 * ksg1_301[k]
                   + f_3 * pc_y[k] * ksh_423[k];

        t_567[k] = f_4 * ksg0_302[k]
                   - f_5 * ksg1_302[k]
                   + f_3 * pc_y[k] * ksh_424[k];

        t_568[k] = f_3 * pc_y[k] * ksh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, pc_y, ish_429, ksg0_303, ksg0_304, \
                         ksg0_309, ksg1_303, ksg1_304, ksg1_309, ksh_426, ksh_427, \
                         ksh_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_12 * ish_429[k]
                   + f_6 * ksg0_309[k]
                   - f_7 * ksg1_309[k]
                   + f_3 * pc_x[k] * ksh_429[k];

        t_570[k] = f_8 * ksg0_303[k]
                   - f_9 * ksg1_303[k]
                   + f_3 * pc_y[k] * ksh_426[k];

        t_571[k] = f_6 * ksg0_304[k]
                   - f_7 * ksg1_304[k]
                   + f_3 * pc_y[k] * ksh_427[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, pc_y, ish_434, ish_435, ksg0_305, \
                         ksg0_314, ksg1_305, ksg1_314, ksh_428, ksh_429, ksh_434, \
                         ksh_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * ksg0_305[k]
                   - f_5 * ksg1_305[k]
                   + f_3 * pc_y[k] * ksh_428[k];

        t_573[k] = f_3 * pc_y[k] * ksh_429[k];

        t_574[k] = f_12 * ish_434[k]
                   + f_4 * ksg0_314[k]
                   - f_5 * ksg1_314[k]
                   + f_3 * pc_x[k] * ksh_434[k];

        t_575[k] = f_12 * ish_435[k]
                   + f_3 * pc_x[k] * ksh_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, pc_x, pc_y, ish_436, ish_437, \
                         ish_438, ish_440, ksh_434, ksh_436, ksh_437, ksh_438, \
                         ksh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_12 * ish_436[k]
                   + f_3 * pc_x[k] * ksh_436[k];

        t_577[k] = f_12 * ish_437[k]
                   + f_3 * pc_x[k] * ksh_437[k];

        t_578[k] = f_12 * ish_438[k]
                   + f_3 * pc_x[k] * ksh_438[k];

        t_579[k] = f_3 * pc_y[k] * ksh_434[k];

        t_580[k] = f_12 * ish_440[k]
                   + f_3 * pc_x[k] * ksh_440[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, ksg0_310, ksg0_311, ksg0_312, ksg1_310, \
                         ksg1_311, ksg1_312, ksh_435, ksh_436, \
                         ksh_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_1 * ksg0_310[k]
                   - f_2 * ksg1_310[k]
                   + f_3 * pc_y[k] * ksh_435[k];

        t_582[k] = f_16 * ksg0_311[k]
                   - f_17 * ksg1_311[k]
                   + f_3 * pc_y[k] * ksh_436[k];

        t_583[k] = f_8 * ksg0_312[k]
                   - f_9 * ksg1_312[k]
                   + f_3 * pc_y[k] * ksh_437[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pc_y, pc_z, ish_314, ksg0_313, ksg0_314, \
                         ksg1_313, ksg1_314, ksh_438, ksh_439, \
                         ksh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_6 * ksg0_313[k]
                   - f_7 * ksg1_313[k]
                   + f_3 * pc_y[k] * ksh_438[k];

        t_585[k] = f_4 * ksg0_314[k]
                   - f_5 * ksg1_314[k]
                   + f_3 * pc_y[k] * ksh_439[k];

        t_586[k] = f_3 * pc_y[k] * ksh_440[k];

        t_587[k] = f_18 * ish_314[k]
                   + f_1 * ksg0_314[k]
                   - f_2 * ksg1_314[k]
                   + f_3 * pc_z[k] * ksh_440[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pa_x, pc_x, pc_y, pc_z, isi0_588, \
                         isi0_591, ish_315, ish_441, ish_444, isi1_588, isi1_591, \
                         ksh_441 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_x[k] * isi0_588[k]
                   + f_15 * ish_441[k]
                   - f_10 * pc_x[k] * isi1_588[k];

        t_589[k] = f_15 * ish_315[k]
                   + f_3 * pc_y[k] * ksh_441[k];

        t_590[k] = f_3 * pc_z[k] * ksh_441[k];

        t_591[k] = pa_x[k] * isi0_591[k]
                   + f_14 * ish_444[k]
                   - f_10 * pc_x[k] * isi1_591[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

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
    const auto f_15 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isi0_420 = buffer.data(isi0 + 420);
    const auto *isi0_423 = buffer.data(isi0 + 423);
    const auto *isi0_426 = buffer.data(isi0 + 426);
    const auto *isi0_430 = buffer.data(isi0 + 430);
    const auto *isi0_594 = buffer.data(isi0 + 594);
    const auto *isi0_598 = buffer.data(isi0 + 598);
    const auto *isi0_609 = buffer.data(isi0 + 609);
    const auto *isi0_611 = buffer.data(isi0 + 611);
    const auto *isi0_612 = buffer.data(isi0 + 612);
    const auto *isi0_613 = buffer.data(isi0 + 613);
    const auto *isi0_615 = buffer.data(isi0 + 615);
    const auto *isi0_621 = buffer.data(isi0 + 621);
    const auto *isi0_625 = buffer.data(isi0 + 625);
    const auto *isi0_628 = buffer.data(isi0 + 628);
    const auto *isi0_630 = buffer.data(isi0 + 630);
    const auto *isi0_637 = buffer.data(isi0 + 637);
    const auto *isi0_639 = buffer.data(isi0 + 639);
    const auto *isi0_640 = buffer.data(isi0 + 640);
    const auto *isi0_641 = buffer.data(isi0 + 641);
    const auto *isi0_643 = buffer.data(isi0 + 643);
    const auto *isi0_644 = buffer.data(isi0 + 644);
    const auto *isi0_647 = buffer.data(isi0 + 647);
    const auto *isi0_649 = buffer.data(isi0 + 649);
    const auto *isi0_650 = buffer.data(isi0 + 650);
    const auto *isi0_653 = buffer.data(isi0 + 653);
    const auto *isi0_654 = buffer.data(isi0 + 654);
    const auto *isi0_656 = buffer.data(isi0 + 656);
    const auto *isi0_658 = buffer.data(isi0 + 658);
    const auto *isi0_665 = buffer.data(isi0 + 665);
    const auto *isi0_667 = buffer.data(isi0 + 667);
    const auto *isi0_668 = buffer.data(isi0 + 668);
    const auto *isi0_669 = buffer.data(isi0 + 669);
    const auto *isi0_671 = buffer.data(isi0 + 671);
    const auto *isi0_672 = buffer.data(isi0 + 672);
    const auto *isi0_675 = buffer.data(isi0 + 675);
    const auto *isi0_677 = buffer.data(isi0 + 677);
    const auto *isi0_678 = buffer.data(isi0 + 678);
    const auto *isi0_681 = buffer.data(isi0 + 681);
    const auto *isi0_682 = buffer.data(isi0 + 682);
    const auto *isi0_684 = buffer.data(isi0 + 684);
    const auto *isi0_686 = buffer.data(isi0 + 686);
    const auto *isi0_693 = buffer.data(isi0 + 693);
    const auto *isi0_695 = buffer.data(isi0 + 695);
    const auto *isi0_696 = buffer.data(isi0 + 696);
    const auto *isi0_697 = buffer.data(isi0 + 697);
    const auto *isi0_699 = buffer.data(isi0 + 699);
    const auto *isi0_700 = buffer.data(isi0 + 700);
    const auto *isi0_703 = buffer.data(isi0 + 703);
    const auto *isi0_705 = buffer.data(isi0 + 705);
    const auto *isi0_706 = buffer.data(isi0 + 706);
    const auto *isi0_709 = buffer.data(isi0 + 709);
    const auto *isi0_710 = buffer.data(isi0 + 710);
    const auto *isi0_712 = buffer.data(isi0 + 712);

    const auto *ish_315 = buffer.data(ish + 315);
    const auto *ish_318 = buffer.data(ish + 318);
    const auto *ish_320 = buffer.data(ish + 320);
    const auto *ish_321 = buffer.data(ish + 321);
    const auto *ish_324 = buffer.data(ish + 324);
    const auto *ish_330 = buffer.data(ish + 330);
    const auto *ish_335 = buffer.data(ish + 335);
    const auto *ish_336 = buffer.data(ish + 336);
    const auto *ish_338 = buffer.data(ish + 338);
    const auto *ish_339 = buffer.data(ish + 339);
    const auto *ish_341 = buffer.data(ish + 341);
    const auto *ish_342 = buffer.data(ish + 342);
    const auto *ish_345 = buffer.data(ish + 345);
    const auto *ish_351 = buffer.data(ish + 351);
    const auto *ish_356 = buffer.data(ish + 356);
    const auto *ish_357 = buffer.data(ish + 357);
    const auto *ish_359 = buffer.data(ish + 359);
    const auto *ish_360 = buffer.data(ish + 360);
    const auto *ish_362 = buffer.data(ish + 362);
    const auto *ish_363 = buffer.data(ish + 363);
    const auto *ish_366 = buffer.data(ish + 366);
    const auto *ish_372 = buffer.data(ish + 372);
    const auto *ish_377 = buffer.data(ish + 377);
    const auto *ish_378 = buffer.data(ish + 378);
    const auto *ish_380 = buffer.data(ish + 380);
    const auto *ish_381 = buffer.data(ish + 381);
    const auto *ish_383 = buffer.data(ish + 383);
    const auto *ish_384 = buffer.data(ish + 384);
    const auto *ish_387 = buffer.data(ish + 387);
    const auto *ish_398 = buffer.data(ish + 398);
    const auto *ish_399 = buffer.data(ish + 399);
    const auto *ish_401 = buffer.data(ish + 401);
    const auto *ish_404 = buffer.data(ish + 404);
    const auto *ish_408 = buffer.data(ish + 408);
    const auto *ish_447 = buffer.data(ish + 447);
    const auto *ish_451 = buffer.data(ish + 451);
    const auto *ish_456 = buffer.data(ish + 456);
    const auto *ish_458 = buffer.data(ish + 458);
    const auto *ish_459 = buffer.data(ish + 459);
    const auto *ish_460 = buffer.data(ish + 460);
    const auto *ish_461 = buffer.data(ish + 461);
    const auto *ish_467 = buffer.data(ish + 467);
    const auto *ish_471 = buffer.data(ish + 471);
    const auto *ish_474 = buffer.data(ish + 474);
    const auto *ish_476 = buffer.data(ish + 476);
    const auto *ish_477 = buffer.data(ish + 477);
    const auto *ish_478 = buffer.data(ish + 478);
    const auto *ish_479 = buffer.data(ish + 479);
    const auto *ish_480 = buffer.data(ish + 480);
    const auto *ish_481 = buffer.data(ish + 481);
    const auto *ish_482 = buffer.data(ish + 482);
    const auto *ish_483 = buffer.data(ish + 483);
    const auto *ish_486 = buffer.data(ish + 486);
    const auto *ish_488 = buffer.data(ish + 488);
    const auto *ish_489 = buffer.data(ish + 489);
    const auto *ish_492 = buffer.data(ish + 492);
    const auto *ish_493 = buffer.data(ish + 493);
    const auto *ish_495 = buffer.data(ish + 495);
    const auto *ish_497 = buffer.data(ish + 497);
    const auto *ish_498 = buffer.data(ish + 498);
    const auto *ish_499 = buffer.data(ish + 499);
    const auto *ish_500 = buffer.data(ish + 500);
    const auto *ish_501 = buffer.data(ish + 501);
    const auto *ish_502 = buffer.data(ish + 502);
    const auto *ish_503 = buffer.data(ish + 503);
    const auto *ish_504 = buffer.data(ish + 504);
    const auto *ish_507 = buffer.data(ish + 507);
    const auto *ish_509 = buffer.data(ish + 509);
    const auto *ish_510 = buffer.data(ish + 510);
    const auto *ish_513 = buffer.data(ish + 513);
    const auto *ish_514 = buffer.data(ish + 514);
    const auto *ish_516 = buffer.data(ish + 516);
    const auto *ish_518 = buffer.data(ish + 518);
    const auto *ish_519 = buffer.data(ish + 519);
    const auto *ish_520 = buffer.data(ish + 520);
    const auto *ish_521 = buffer.data(ish + 521);
    const auto *ish_522 = buffer.data(ish + 522);
    const auto *ish_523 = buffer.data(ish + 523);
    const auto *ish_524 = buffer.data(ish + 524);
    const auto *ish_525 = buffer.data(ish + 525);
    const auto *ish_528 = buffer.data(ish + 528);
    const auto *ish_530 = buffer.data(ish + 530);
    const auto *ish_531 = buffer.data(ish + 531);
    const auto *ish_534 = buffer.data(ish + 534);
    const auto *ish_535 = buffer.data(ish + 535);
    const auto *ish_537 = buffer.data(ish + 537);

    const auto *isi1_420 = buffer.data(isi1 + 420);
    const auto *isi1_423 = buffer.data(isi1 + 423);
    const auto *isi1_426 = buffer.data(isi1 + 426);
    const auto *isi1_430 = buffer.data(isi1 + 430);
    const auto *isi1_594 = buffer.data(isi1 + 594);
    const auto *isi1_598 = buffer.data(isi1 + 598);
    const auto *isi1_609 = buffer.data(isi1 + 609);
    const auto *isi1_611 = buffer.data(isi1 + 611);
    const auto *isi1_612 = buffer.data(isi1 + 612);
    const auto *isi1_613 = buffer.data(isi1 + 613);
    const auto *isi1_615 = buffer.data(isi1 + 615);
    const auto *isi1_621 = buffer.data(isi1 + 621);
    const auto *isi1_625 = buffer.data(isi1 + 625);
    const auto *isi1_628 = buffer.data(isi1 + 628);
    const auto *isi1_630 = buffer.data(isi1 + 630);
    const auto *isi1_637 = buffer.data(isi1 + 637);
    const auto *isi1_639 = buffer.data(isi1 + 639);
    const auto *isi1_640 = buffer.data(isi1 + 640);
    const auto *isi1_641 = buffer.data(isi1 + 641);
    const auto *isi1_643 = buffer.data(isi1 + 643);
    const auto *isi1_644 = buffer.data(isi1 + 644);
    const auto *isi1_647 = buffer.data(isi1 + 647);
    const auto *isi1_649 = buffer.data(isi1 + 649);
    const auto *isi1_650 = buffer.data(isi1 + 650);
    const auto *isi1_653 = buffer.data(isi1 + 653);
    const auto *isi1_654 = buffer.data(isi1 + 654);
    const auto *isi1_656 = buffer.data(isi1 + 656);
    const auto *isi1_658 = buffer.data(isi1 + 658);
    const auto *isi1_665 = buffer.data(isi1 + 665);
    const auto *isi1_667 = buffer.data(isi1 + 667);
    const auto *isi1_668 = buffer.data(isi1 + 668);
    const auto *isi1_669 = buffer.data(isi1 + 669);
    const auto *isi1_671 = buffer.data(isi1 + 671);
    const auto *isi1_672 = buffer.data(isi1 + 672);
    const auto *isi1_675 = buffer.data(isi1 + 675);
    const auto *isi1_677 = buffer.data(isi1 + 677);
    const auto *isi1_678 = buffer.data(isi1 + 678);
    const auto *isi1_681 = buffer.data(isi1 + 681);
    const auto *isi1_682 = buffer.data(isi1 + 682);
    const auto *isi1_684 = buffer.data(isi1 + 684);
    const auto *isi1_686 = buffer.data(isi1 + 686);
    const auto *isi1_693 = buffer.data(isi1 + 693);
    const auto *isi1_695 = buffer.data(isi1 + 695);
    const auto *isi1_696 = buffer.data(isi1 + 696);
    const auto *isi1_697 = buffer.data(isi1 + 697);
    const auto *isi1_699 = buffer.data(isi1 + 699);
    const auto *isi1_700 = buffer.data(isi1 + 700);
    const auto *isi1_703 = buffer.data(isi1 + 703);
    const auto *isi1_705 = buffer.data(isi1 + 705);
    const auto *isi1_706 = buffer.data(isi1 + 706);
    const auto *isi1_709 = buffer.data(isi1 + 709);
    const auto *isi1_710 = buffer.data(isi1 + 710);
    const auto *isi1_712 = buffer.data(isi1 + 712);

    const auto *ksg0_315 = buffer.data(ksg0 + 315);
    const auto *ksg0_317 = buffer.data(ksg0 + 317);
    const auto *ksg0_318 = buffer.data(ksg0 + 318);
    const auto *ksg0_320 = buffer.data(ksg0 + 320);

    const auto *ksg1_315 = buffer.data(ksg1 + 315);
    const auto *ksg1_317 = buffer.data(ksg1 + 317);
    const auto *ksg1_318 = buffer.data(ksg1 + 318);
    const auto *ksg1_320 = buffer.data(ksg1 + 320);

    const auto *ksh_442 = buffer.data(ksh + 442);
    const auto *ksh_443 = buffer.data(ksh + 443);
    const auto *ksh_444 = buffer.data(ksh + 444);
    const auto *ksh_446 = buffer.data(ksh + 446);
    const auto *ksh_447 = buffer.data(ksh + 447);
    const auto *ksh_448 = buffer.data(ksh + 448);
    const auto *ksh_450 = buffer.data(ksh + 450);
    const auto *ksh_451 = buffer.data(ksh + 451);
    const auto *ksh_456 = buffer.data(ksh + 456);
    const auto *ksh_458 = buffer.data(ksh + 458);
    const auto *ksh_459 = buffer.data(ksh + 459);
    const auto *ksh_460 = buffer.data(ksh + 460);
    const auto *ksh_461 = buffer.data(ksh + 461);
    const auto *ksh_462 = buffer.data(ksh + 462);
    const auto *ksh_464 = buffer.data(ksh + 464);
    const auto *ksh_465 = buffer.data(ksh + 465);
    const auto *ksh_467 = buffer.data(ksh + 467);
    const auto *ksh_468 = buffer.data(ksh + 468);
    const auto *ksh_471 = buffer.data(ksh + 471);
    const auto *ksh_477 = buffer.data(ksh + 477);
    const auto *ksh_478 = buffer.data(ksh + 478);
    const auto *ksh_479 = buffer.data(ksh + 479);
    const auto *ksh_480 = buffer.data(ksh + 480);
    const auto *ksh_481 = buffer.data(ksh + 481);
    const auto *ksh_482 = buffer.data(ksh + 482);
    const auto *ksh_483 = buffer.data(ksh + 483);
    const auto *ksh_485 = buffer.data(ksh + 485);
    const auto *ksh_486 = buffer.data(ksh + 486);
    const auto *ksh_488 = buffer.data(ksh + 488);
    const auto *ksh_489 = buffer.data(ksh + 489);
    const auto *ksh_492 = buffer.data(ksh + 492);
    const auto *ksh_498 = buffer.data(ksh + 498);
    const auto *ksh_499 = buffer.data(ksh + 499);
    const auto *ksh_500 = buffer.data(ksh + 500);
    const auto *ksh_501 = buffer.data(ksh + 501);
    const auto *ksh_502 = buffer.data(ksh + 502);
    const auto *ksh_503 = buffer.data(ksh + 503);
    const auto *ksh_504 = buffer.data(ksh + 504);
    const auto *ksh_506 = buffer.data(ksh + 506);
    const auto *ksh_507 = buffer.data(ksh + 507);
    const auto *ksh_509 = buffer.data(ksh + 509);
    const auto *ksh_510 = buffer.data(ksh + 510);
    const auto *ksh_513 = buffer.data(ksh + 513);
    const auto *ksh_519 = buffer.data(ksh + 519);
    const auto *ksh_520 = buffer.data(ksh + 520);
    const auto *ksh_521 = buffer.data(ksh + 521);
    const auto *ksh_522 = buffer.data(ksh + 522);
    const auto *ksh_523 = buffer.data(ksh + 523);
    const auto *ksh_524 = buffer.data(ksh + 524);
    const auto *ksh_525 = buffer.data(ksh + 525);
    const auto *ksh_527 = buffer.data(ksh + 527);
    const auto *ksh_528 = buffer.data(ksh + 528);
    const auto *ksh_530 = buffer.data(ksh + 530);
    const auto *ksh_531 = buffer.data(ksh + 531);
    const auto *ksh_534 = buffer.data(ksh + 534);

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pa_x, pc_x, pc_z, isi0_594, ish_447, \
                         isi1_594, ksg0_315, ksg1_315, ksh_442, ksh_443, \
                         ksh_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_3 * pc_z[k] * ksh_442[k];

        t_593[k] = f_4 * ksg0_315[k]
                   - f_5 * ksg1_315[k]
                   + f_3 * pc_z[k] * ksh_443[k];

        t_594[k] = pa_x[k] * isi0_594[k]
                   + f_13 * ish_447[k]
                   - f_10 * pc_x[k] * isi1_594[k];

        t_595[k] = f_3 * pc_z[k] * ksh_444[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pa_x, pc_x, pc_y, pc_z, isi0_598, \
                         ish_320, ish_451, isi1_598, ksg0_317, ksg1_317, ksh_446, \
                         ksh_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_15 * ish_320[k]
                   + f_3 * pc_y[k] * ksh_446[k];

        t_597[k] = f_6 * ksg0_317[k]
                   - f_7 * ksg1_317[k]
                   + f_3 * pc_z[k] * ksh_446[k];

        t_598[k] = pa_x[k] * isi0_598[k]
                   + f_12 * ish_451[k]
                   - f_10 * pc_x[k] * isi1_598[k];

        t_599[k] = f_3 * pc_z[k] * ksh_447[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pc_x, pc_y, pc_z, ish_324, ish_456, \
                         ksg0_318, ksg0_320, ksg1_318, ksg1_320, ksh_448, ksh_450, \
                         ksh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_4 * ksg0_318[k]
                   - f_5 * ksg1_318[k]
                   + f_3 * pc_z[k] * ksh_448[k];

        t_601[k] = f_15 * ish_324[k]
                   + f_3 * pc_y[k] * ksh_450[k];

        t_602[k] = f_8 * ksg0_320[k]
                   - f_9 * ksg1_320[k]
                   + f_3 * pc_z[k] * ksh_450[k];

        t_603[k] = f_11 * ish_456[k]
                   + f_3 * pc_x[k] * ksh_456[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, t_608, pc_x, pc_z, ish_458, ish_459, \
                         ish_460, ish_461, ksh_451, ksh_458, ksh_459, ksh_460, \
                         ksh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = f_3 * pc_z[k] * ksh_451[k];

        t_605[k] = f_11 * ish_458[k]
                   + f_3 * pc_x[k] * ksh_458[k];

        t_606[k] = f_11 * ish_459[k]
                   + f_3 * pc_x[k] * ksh_459[k];

        t_607[k] = f_11 * ish_460[k]
                   + f_3 * pc_x[k] * ksh_460[k];

        t_608[k] = f_11 * ish_461[k]
                   + f_3 * pc_x[k] * ksh_461[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, t_612, pa_x, pc_x, pc_z, isi0_609, isi0_611, \
                         isi0_612, isi1_609, isi1_611, isi1_612, \
                         ksh_456 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = pa_x[k] * isi0_609[k]
                   - f_10 * pc_x[k] * isi1_609[k];

        t_610[k] = f_3 * pc_z[k] * ksh_456[k];

        t_611[k] = pa_x[k] * isi0_611[k]
                   - f_10 * pc_x[k] * isi1_611[k];

        t_612[k] = pa_x[k] * isi0_612[k]
                   - f_10 * pc_x[k] * isi1_612[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pa_x, pc_x, pc_y, isi0_613, isi0_615, ish_335, \
                         isi1_613, isi1_615, ksh_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = pa_x[k] * isi0_613[k]
                   - f_10 * pc_x[k] * isi1_613[k];

        t_614[k] = f_15 * ish_335[k]
                   + f_3 * pc_y[k] * ksh_461[k];

        t_615[k] = pa_x[k] * isi0_615[k]
                   - f_10 * pc_x[k] * isi1_615[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pa_z, pc_y, pc_z, isi0_420, isi0_423, \
                         ish_315, ish_336, isi1_420, isi1_423, \
                         ksh_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pa_z[k] * isi0_420[k]
                   - f_10 * pc_z[k] * isi1_420[k];

        t_617[k] = f_18 * ish_336[k]
                   + f_3 * pc_y[k] * ksh_462[k];

        t_618[k] = f_11 * ish_315[k]
                   + f_3 * pc_z[k] * ksh_462[k];

        t_619[k] = pa_z[k] * isi0_423[k]
                   - f_10 * pc_z[k] * isi1_423[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_x, pa_z, pc_x, pc_y, pc_z, isi0_426, \
                         isi0_621, ish_338, ish_467, isi1_426, isi1_621, \
                         ksh_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_18 * ish_338[k]
                   + f_3 * pc_y[k] * ksh_464[k];

        t_621[k] = pa_x[k] * isi0_621[k]
                   + f_14 * ish_467[k]
                   - f_10 * pc_x[k] * isi1_621[k];

        t_622[k] = pa_z[k] * isi0_426[k]
                   - f_10 * pc_z[k] * isi1_426[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pa_x, pc_x, pc_y, pc_z, isi0_625, ish_318, \
                         ish_341, ish_471, isi1_625, ksh_465, ksh_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_11 * ish_318[k]
                   + f_3 * pc_z[k] * ksh_465[k];

        t_624[k] = f_18 * ish_341[k]
                   + f_3 * pc_y[k] * ksh_467[k];

        t_625[k] = pa_x[k] * isi0_625[k]
                   + f_13 * ish_471[k]
                   - f_10 * pc_x[k] * isi1_625[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_x, pa_z, pc_x, pc_z, isi0_430, isi0_628, \
                         ish_321, ish_474, isi1_430, isi1_628, \
                         ksh_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_z[k] * isi0_430[k]
                   - f_10 * pc_z[k] * isi1_430[k];

        t_627[k] = f_11 * ish_321[k]
                   + f_3 * pc_z[k] * ksh_468[k];

        t_628[k] = pa_x[k] * isi0_628[k]
                   + f_12 * ish_474[k]
                   - f_10 * pc_x[k] * isi1_628[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pa_x, pc_x, pc_y, isi0_630, ish_345, \
                         ish_476, ish_477, ish_478, isi1_630, ksh_471, ksh_477, \
                         ksh_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_18 * ish_345[k]
                   + f_3 * pc_y[k] * ksh_471[k];

        t_630[k] = pa_x[k] * isi0_630[k]
                   + f_12 * ish_476[k]
                   - f_10 * pc_x[k] * isi1_630[k];

        t_631[k] = f_11 * ish_477[k]
                   + f_3 * pc_x[k] * ksh_477[k];

        t_632[k] = f_11 * ish_478[k]
                   + f_3 * pc_x[k] * ksh_478[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pc_x, ish_479, ish_480, ish_481, ish_482, \
                         ksh_479, ksh_480, ksh_481, ksh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_11 * ish_479[k]
                   + f_3 * pc_x[k] * ksh_479[k];

        t_634[k] = f_11 * ish_480[k]
                   + f_3 * pc_x[k] * ksh_480[k];

        t_635[k] = f_11 * ish_481[k]
                   + f_3 * pc_x[k] * ksh_481[k];

        t_636[k] = f_11 * ish_482[k]
                   + f_3 * pc_x[k] * ksh_482[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pa_x, pc_x, pc_z, isi0_637, isi0_639, \
                         isi0_640, ish_330, isi1_637, isi1_639, isi1_640, \
                         ksh_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = pa_x[k] * isi0_637[k]
                   - f_10 * pc_x[k] * isi1_637[k];

        t_638[k] = f_11 * ish_330[k]
                   + f_3 * pc_z[k] * ksh_477[k];

        t_639[k] = pa_x[k] * isi0_639[k]
                   - f_10 * pc_x[k] * isi1_639[k];

        t_640[k] = pa_x[k] * isi0_640[k]
                   - f_10 * pc_x[k] * isi1_640[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pa_x, pc_x, pc_y, isi0_641, isi0_643, \
                         isi0_644, ish_356, ish_483, isi1_641, isi1_643, isi1_644, \
                         ksh_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = pa_x[k] * isi0_641[k]
                   - f_10 * pc_x[k] * isi1_641[k];

        t_642[k] = f_18 * ish_356[k]
                   + f_3 * pc_y[k] * ksh_482[k];

        t_643[k] = pa_x[k] * isi0_643[k]
                   - f_10 * pc_x[k] * isi1_643[k];

        t_644[k] = pa_x[k] * isi0_644[k]
                   + f_15 * ish_483[k]
                   - f_10 * pc_x[k] * isi1_644[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pa_x, pc_x, pc_y, pc_z, isi0_647, \
                         ish_336, ish_357, ish_359, ish_486, isi1_647, ksh_483, \
                         ksh_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_14 * ish_357[k]
                   + f_3 * pc_y[k] * ksh_483[k];

        t_646[k] = f_12 * ish_336[k]
                   + f_3 * pc_z[k] * ksh_483[k];

        t_647[k] = pa_x[k] * isi0_647[k]
                   + f_14 * ish_486[k]
                   - f_10 * pc_x[k] * isi1_647[k];

        t_648[k] = f_14 * ish_359[k]
                   + f_3 * pc_y[k] * ksh_485[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, pa_x, pc_x, pc_z, isi0_649, isi0_650, ish_339, \
                         ish_488, ish_489, isi1_649, isi1_650, \
                         ksh_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = pa_x[k] * isi0_649[k]
                   + f_14 * ish_488[k]
                   - f_10 * pc_x[k] * isi1_649[k];

        t_650[k] = pa_x[k] * isi0_650[k]
                   + f_13 * ish_489[k]
                   - f_10 * pc_x[k] * isi1_650[k];

        t_651[k] = f_12 * ish_339[k]
                   + f_3 * pc_z[k] * ksh_486[k];
    }

#pragma omp simd aligned(t_652, t_653, t_654, pa_x, pc_x, pc_y, isi0_653, isi0_654, ish_362, \
                         ish_492, ish_493, isi1_653, isi1_654, \
                         ksh_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_652[k] = f_14 * ish_362[k]
                   + f_3 * pc_y[k] * ksh_488[k];

        t_653[k] = pa_x[k] * isi0_653[k]
                   + f_13 * ish_492[k]
                   - f_10 * pc_x[k] * isi1_653[k];

        t_654[k] = pa_x[k] * isi0_654[k]
                   + f_12 * ish_493[k]
                   - f_10 * pc_x[k] * isi1_654[k];
    }

#pragma omp simd aligned(t_655, t_656, t_657, pa_x, pc_x, pc_y, pc_z, isi0_656, ish_342, \
                         ish_366, ish_495, isi1_656, ksh_489, ksh_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_655[k] = f_12 * ish_342[k]
                   + f_3 * pc_z[k] * ksh_489[k];

        t_656[k] = pa_x[k] * isi0_656[k]
                   + f_12 * ish_495[k]
                   - f_10 * pc_x[k] * isi1_656[k];

        t_657[k] = f_14 * ish_366[k]
                   + f_3 * pc_y[k] * ksh_492[k];
    }

#pragma omp simd aligned(t_658, t_659, t_660, t_661, pa_x, pc_x, isi0_658, ish_497, ish_498, \
                         ish_499, ish_500, isi1_658, ksh_498, ksh_499, \
                         ksh_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_658[k] = pa_x[k] * isi0_658[k]
                   + f_12 * ish_497[k]
                   - f_10 * pc_x[k] * isi1_658[k];

        t_659[k] = f_11 * ish_498[k]
                   + f_3 * pc_x[k] * ksh_498[k];

        t_660[k] = f_11 * ish_499[k]
                   + f_3 * pc_x[k] * ksh_499[k];

        t_661[k] = f_11 * ish_500[k]
                   + f_3 * pc_x[k] * ksh_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pa_x, pc_x, isi0_665, ish_501, ish_502, \
                         ish_503, isi1_665, ksh_501, ksh_502, ksh_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_11 * ish_501[k]
                   + f_3 * pc_x[k] * ksh_501[k];

        t_663[k] = f_11 * ish_502[k]
                   + f_3 * pc_x[k] * ksh_502[k];

        t_664[k] = f_11 * ish_503[k]
                   + f_3 * pc_x[k] * ksh_503[k];

        t_665[k] = pa_x[k] * isi0_665[k]
                   - f_10 * pc_x[k] * isi1_665[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pa_x, pc_x, pc_z, isi0_667, isi0_668, \
                         isi0_669, ish_351, isi1_667, isi1_668, isi1_669, \
                         ksh_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_12 * ish_351[k]
                   + f_3 * pc_z[k] * ksh_498[k];

        t_667[k] = pa_x[k] * isi0_667[k]
                   - f_10 * pc_x[k] * isi1_667[k];

        t_668[k] = pa_x[k] * isi0_668[k]
                   - f_10 * pc_x[k] * isi1_668[k];

        t_669[k] = pa_x[k] * isi0_669[k]
                   - f_10 * pc_x[k] * isi1_669[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pa_x, pc_x, pc_y, isi0_671, isi0_672, \
                         ish_377, ish_378, ish_504, isi1_671, isi1_672, ksh_503, \
                         ksh_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_14 * ish_377[k]
                   + f_3 * pc_y[k] * ksh_503[k];

        t_671[k] = pa_x[k] * isi0_671[k]
                   - f_10 * pc_x[k] * isi1_671[k];

        t_672[k] = pa_x[k] * isi0_672[k]
                   + f_15 * ish_504[k]
                   - f_10 * pc_x[k] * isi1_672[k];

        t_673[k] = f_13 * ish_378[k]
                   + f_3 * pc_y[k] * ksh_504[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, pa_x, pc_x, pc_y, pc_z, isi0_675, ish_357, \
                         ish_380, ish_507, isi1_675, ksh_504, ksh_506 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_13 * ish_357[k]
                   + f_3 * pc_z[k] * ksh_504[k];

        t_675[k] = pa_x[k] * isi0_675[k]
                   + f_14 * ish_507[k]
                   - f_10 * pc_x[k] * isi1_675[k];

        t_676[k] = f_13 * ish_380[k]
                   + f_3 * pc_y[k] * ksh_506[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pa_x, pc_x, pc_z, isi0_677, isi0_678, ish_360, \
                         ish_509, ish_510, isi1_677, isi1_678, \
                         ksh_507 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = pa_x[k] * isi0_677[k]
                   + f_14 * ish_509[k]
                   - f_10 * pc_x[k] * isi1_677[k];

        t_678[k] = pa_x[k] * isi0_678[k]
                   + f_13 * ish_510[k]
                   - f_10 * pc_x[k] * isi1_678[k];

        t_679[k] = f_13 * ish_360[k]
                   + f_3 * pc_z[k] * ksh_507[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pa_x, pc_x, pc_y, isi0_681, isi0_682, ish_383, \
                         ish_513, ish_514, isi1_681, isi1_682, \
                         ksh_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_13 * ish_383[k]
                   + f_3 * pc_y[k] * ksh_509[k];

        t_681[k] = pa_x[k] * isi0_681[k]
                   + f_13 * ish_513[k]
                   - f_10 * pc_x[k] * isi1_681[k];

        t_682[k] = pa_x[k] * isi0_682[k]
                   + f_12 * ish_514[k]
                   - f_10 * pc_x[k] * isi1_682[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pa_x, pc_x, pc_y, pc_z, isi0_684, ish_363, \
                         ish_387, ish_516, isi1_684, ksh_510, ksh_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_13 * ish_363[k]
                   + f_3 * pc_z[k] * ksh_510[k];

        t_684[k] = pa_x[k] * isi0_684[k]
                   + f_12 * ish_516[k]
                   - f_10 * pc_x[k] * isi1_684[k];

        t_685[k] = f_13 * ish_387[k]
                   + f_3 * pc_y[k] * ksh_513[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, t_689, pa_x, pc_x, isi0_686, ish_518, ish_519, \
                         ish_520, ish_521, isi1_686, ksh_519, ksh_520, \
                         ksh_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = pa_x[k] * isi0_686[k]
                   + f_12 * ish_518[k]
                   - f_10 * pc_x[k] * isi1_686[k];

        t_687[k] = f_11 * ish_519[k]
                   + f_3 * pc_x[k] * ksh_519[k];

        t_688[k] = f_11 * ish_520[k]
                   + f_3 * pc_x[k] * ksh_520[k];

        t_689[k] = f_11 * ish_521[k]
                   + f_3 * pc_x[k] * ksh_521[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, t_693, pa_x, pc_x, isi0_693, ish_522, ish_523, \
                         ish_524, isi1_693, ksh_522, ksh_523, ksh_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_11 * ish_522[k]
                   + f_3 * pc_x[k] * ksh_522[k];

        t_691[k] = f_11 * ish_523[k]
                   + f_3 * pc_x[k] * ksh_523[k];

        t_692[k] = f_11 * ish_524[k]
                   + f_3 * pc_x[k] * ksh_524[k];

        t_693[k] = pa_x[k] * isi0_693[k]
                   - f_10 * pc_x[k] * isi1_693[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, t_697, pa_x, pc_x, pc_z, isi0_695, isi0_696, \
                         isi0_697, ish_372, isi1_695, isi1_696, isi1_697, \
                         ksh_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = f_13 * ish_372[k]
                   + f_3 * pc_z[k] * ksh_519[k];

        t_695[k] = pa_x[k] * isi0_695[k]
                   - f_10 * pc_x[k] * isi1_695[k];

        t_696[k] = pa_x[k] * isi0_696[k]
                   - f_10 * pc_x[k] * isi1_696[k];

        t_697[k] = pa_x[k] * isi0_697[k]
                   - f_10 * pc_x[k] * isi1_697[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, t_701, pa_x, pc_x, pc_y, isi0_699, isi0_700, \
                         ish_398, ish_399, ish_525, isi1_699, isi1_700, ksh_524, \
                         ksh_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = f_13 * ish_398[k]
                   + f_3 * pc_y[k] * ksh_524[k];

        t_699[k] = pa_x[k] * isi0_699[k]
                   - f_10 * pc_x[k] * isi1_699[k];

        t_700[k] = pa_x[k] * isi0_700[k]
                   + f_15 * ish_525[k]
                   - f_10 * pc_x[k] * isi1_700[k];

        t_701[k] = f_12 * ish_399[k]
                   + f_3 * pc_y[k] * ksh_525[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, pa_x, pc_x, pc_y, pc_z, isi0_703, ish_378, \
                         ish_401, ish_528, isi1_703, ksh_525, ksh_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_14 * ish_378[k]
                   + f_3 * pc_z[k] * ksh_525[k];

        t_703[k] = pa_x[k] * isi0_703[k]
                   + f_14 * ish_528[k]
                   - f_10 * pc_x[k] * isi1_703[k];

        t_704[k] = f_12 * ish_401[k]
                   + f_3 * pc_y[k] * ksh_527[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, pa_x, pc_x, pc_z, isi0_705, isi0_706, ish_381, \
                         ish_530, ish_531, isi1_705, isi1_706, \
                         ksh_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = pa_x[k] * isi0_705[k]
                   + f_14 * ish_530[k]
                   - f_10 * pc_x[k] * isi1_705[k];

        t_706[k] = pa_x[k] * isi0_706[k]
                   + f_13 * ish_531[k]
                   - f_10 * pc_x[k] * isi1_706[k];

        t_707[k] = f_14 * ish_381[k]
                   + f_3 * pc_z[k] * ksh_528[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, pa_x, pc_x, pc_y, isi0_709, isi0_710, ish_404, \
                         ish_534, ish_535, isi1_709, isi1_710, \
                         ksh_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = f_12 * ish_404[k]
                   + f_3 * pc_y[k] * ksh_530[k];

        t_709[k] = pa_x[k] * isi0_709[k]
                   + f_13 * ish_534[k]
                   - f_10 * pc_x[k] * isi1_709[k];

        t_710[k] = pa_x[k] * isi0_710[k]
                   + f_12 * ish_535[k]
                   - f_10 * pc_x[k] * isi1_710[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pa_x, pc_x, pc_y, pc_z, isi0_712, ish_384, \
                         ish_408, ish_537, isi1_712, ksh_531, ksh_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_14 * ish_384[k]
                   + f_3 * pc_z[k] * ksh_531[k];

        t_712[k] = pa_x[k] * isi0_712[k]
                   + f_12 * ish_537[k]
                   - f_10 * pc_x[k] * isi1_712[k];

        t_713[k] = f_12 * ish_408[k]
                   + f_3 * pc_y[k] * ksh_534[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isi0_560 = buffer.data(isi0 + 560);
    const auto *isi0_565 = buffer.data(isi0 + 565);
    const auto *isi0_569 = buffer.data(isi0 + 569);
    const auto *isi0_574 = buffer.data(isi0 + 574);
    const auto *isi0_588 = buffer.data(isi0 + 588);
    const auto *isi0_589 = buffer.data(isi0 + 589);
    const auto *isi0_591 = buffer.data(isi0 + 591);
    const auto *isi0_594 = buffer.data(isi0 + 594);
    const auto *isi0_598 = buffer.data(isi0 + 598);
    const auto *isi0_609 = buffer.data(isi0 + 609);
    const auto *isi0_611 = buffer.data(isi0 + 611);
    const auto *isi0_612 = buffer.data(isi0 + 612);
    const auto *isi0_613 = buffer.data(isi0 + 613);
    const auto *isi0_714 = buffer.data(isi0 + 714);
    const auto *isi0_721 = buffer.data(isi0 + 721);
    const auto *isi0_723 = buffer.data(isi0 + 723);
    const auto *isi0_724 = buffer.data(isi0 + 724);
    const auto *isi0_725 = buffer.data(isi0 + 725);
    const auto *isi0_727 = buffer.data(isi0 + 727);
    const auto *isi0_731 = buffer.data(isi0 + 731);
    const auto *isi0_734 = buffer.data(isi0 + 734);
    const auto *isi0_738 = buffer.data(isi0 + 738);
    const auto *isi0_740 = buffer.data(isi0 + 740);
    const auto *isi0_749 = buffer.data(isi0 + 749);
    const auto *isi0_751 = buffer.data(isi0 + 751);
    const auto *isi0_752 = buffer.data(isi0 + 752);
    const auto *isi0_753 = buffer.data(isi0 + 753);
    const auto *isi0_755 = buffer.data(isi0 + 755);
    const auto *isi0_756 = buffer.data(isi0 + 756);
    const auto *isi0_761 = buffer.data(isi0 + 761);
    const auto *isi0_765 = buffer.data(isi0 + 765);
    const auto *isi0_770 = buffer.data(isi0 + 770);
    const auto *isi0_777 = buffer.data(isi0 + 777);
    const auto *isi0_778 = buffer.data(isi0 + 778);
    const auto *isi0_779 = buffer.data(isi0 + 779);
    const auto *isi0_780 = buffer.data(isi0 + 780);
    const auto *isi0_781 = buffer.data(isi0 + 781);
    const auto *isi0_783 = buffer.data(isi0 + 783);

    const auto *ish_393 = buffer.data(ish + 393);
    const auto *ish_399 = buffer.data(ish + 399);
    const auto *ish_402 = buffer.data(ish + 402);
    const auto *ish_405 = buffer.data(ish + 405);
    const auto *ish_414 = buffer.data(ish + 414);
    const auto *ish_419 = buffer.data(ish + 419);
    const auto *ish_420 = buffer.data(ish + 420);
    const auto *ish_422 = buffer.data(ish + 422);
    const auto *ish_425 = buffer.data(ish + 425);
    const auto *ish_429 = buffer.data(ish + 429);
    const auto *ish_440 = buffer.data(ish + 440);
    const auto *ish_456 = buffer.data(ish + 456);
    const auto *ish_457 = buffer.data(ish + 457);
    const auto *ish_458 = buffer.data(ish + 458);
    const auto *ish_459 = buffer.data(ish + 459);
    const auto *ish_461 = buffer.data(ish + 461);
    const auto *ish_482 = buffer.data(ish + 482);
    const auto *ish_539 = buffer.data(ish + 539);
    const auto *ish_540 = buffer.data(ish + 540);
    const auto *ish_541 = buffer.data(ish + 541);
    const auto *ish_542 = buffer.data(ish + 542);
    const auto *ish_543 = buffer.data(ish + 543);
    const auto *ish_544 = buffer.data(ish + 544);
    const auto *ish_545 = buffer.data(ish + 545);
    const auto *ish_549 = buffer.data(ish + 549);
    const auto *ish_552 = buffer.data(ish + 552);
    const auto *ish_556 = buffer.data(ish + 556);
    const auto *ish_558 = buffer.data(ish + 558);
    const auto *ish_561 = buffer.data(ish + 561);
    const auto *ish_562 = buffer.data(ish + 562);
    const auto *ish_563 = buffer.data(ish + 563);
    const auto *ish_564 = buffer.data(ish + 564);
    const auto *ish_565 = buffer.data(ish + 565);
    const auto *ish_566 = buffer.data(ish + 566);
    const auto *ish_567 = buffer.data(ish + 567);
    const auto *ish_572 = buffer.data(ish + 572);
    const auto *ish_576 = buffer.data(ish + 576);
    const auto *ish_581 = buffer.data(ish + 581);
    const auto *ish_582 = buffer.data(ish + 582);
    const auto *ish_583 = buffer.data(ish + 583);
    const auto *ish_584 = buffer.data(ish + 584);
    const auto *ish_585 = buffer.data(ish + 585);
    const auto *ish_587 = buffer.data(ish + 587);

    const auto *isi1_560 = buffer.data(isi1 + 560);
    const auto *isi1_565 = buffer.data(isi1 + 565);
    const auto *isi1_569 = buffer.data(isi1 + 569);
    const auto *isi1_574 = buffer.data(isi1 + 574);
    const auto *isi1_588 = buffer.data(isi1 + 588);
    const auto *isi1_589 = buffer.data(isi1 + 589);
    const auto *isi1_591 = buffer.data(isi1 + 591);
    const auto *isi1_594 = buffer.data(isi1 + 594);
    const auto *isi1_598 = buffer.data(isi1 + 598);
    const auto *isi1_609 = buffer.data(isi1 + 609);
    const auto *isi1_611 = buffer.data(isi1 + 611);
    const auto *isi1_612 = buffer.data(isi1 + 612);
    const auto *isi1_613 = buffer.data(isi1 + 613);
    const auto *isi1_714 = buffer.data(isi1 + 714);
    const auto *isi1_721 = buffer.data(isi1 + 721);
    const auto *isi1_723 = buffer.data(isi1 + 723);
    const auto *isi1_724 = buffer.data(isi1 + 724);
    const auto *isi1_725 = buffer.data(isi1 + 725);
    const auto *isi1_727 = buffer.data(isi1 + 727);
    const auto *isi1_731 = buffer.data(isi1 + 731);
    const auto *isi1_734 = buffer.data(isi1 + 734);
    const auto *isi1_738 = buffer.data(isi1 + 738);
    const auto *isi1_740 = buffer.data(isi1 + 740);
    const auto *isi1_749 = buffer.data(isi1 + 749);
    const auto *isi1_751 = buffer.data(isi1 + 751);
    const auto *isi1_752 = buffer.data(isi1 + 752);
    const auto *isi1_753 = buffer.data(isi1 + 753);
    const auto *isi1_755 = buffer.data(isi1 + 755);
    const auto *isi1_756 = buffer.data(isi1 + 756);
    const auto *isi1_761 = buffer.data(isi1 + 761);
    const auto *isi1_765 = buffer.data(isi1 + 765);
    const auto *isi1_770 = buffer.data(isi1 + 770);
    const auto *isi1_777 = buffer.data(isi1 + 777);
    const auto *isi1_778 = buffer.data(isi1 + 778);
    const auto *isi1_779 = buffer.data(isi1 + 779);
    const auto *isi1_780 = buffer.data(isi1 + 780);
    const auto *isi1_781 = buffer.data(isi1 + 781);
    const auto *isi1_783 = buffer.data(isi1 + 783);

    const auto *ksg0_405 = buffer.data(ksg0 + 405);
    const auto *ksg0_406 = buffer.data(ksg0 + 406);
    const auto *ksg0_407 = buffer.data(ksg0 + 407);
    const auto *ksg0_408 = buffer.data(ksg0 + 408);
    const auto *ksg0_409 = buffer.data(ksg0 + 409);
    const auto *ksg0_410 = buffer.data(ksg0 + 410);
    const auto *ksg0_420 = buffer.data(ksg0 + 420);
    const auto *ksg0_421 = buffer.data(ksg0 + 421);
    const auto *ksg0_423 = buffer.data(ksg0 + 423);
    const auto *ksg0_425 = buffer.data(ksg0 + 425);
    const auto *ksg0_426 = buffer.data(ksg0 + 426);
    const auto *ksg0_428 = buffer.data(ksg0 + 428);
    const auto *ksg0_429 = buffer.data(ksg0 + 429);
    const auto *ksg0_430 = buffer.data(ksg0 + 430);
    const auto *ksg0_431 = buffer.data(ksg0 + 431);
    const auto *ksg0_432 = buffer.data(ksg0 + 432);
    const auto *ksg0_433 = buffer.data(ksg0 + 433);
    const auto *ksg0_434 = buffer.data(ksg0 + 434);
    const auto *ksg0_437 = buffer.data(ksg0 + 437);
    const auto *ksg0_439 = buffer.data(ksg0 + 439);
    const auto *ksg0_440 = buffer.data(ksg0 + 440);
    const auto *ksg0_442 = buffer.data(ksg0 + 442);
    const auto *ksg0_443 = buffer.data(ksg0 + 443);
    const auto *ksg0_444 = buffer.data(ksg0 + 444);
    const auto *ksg0_446 = buffer.data(ksg0 + 446);
    const auto *ksg0_447 = buffer.data(ksg0 + 447);
    const auto *ksg0_448 = buffer.data(ksg0 + 448);
    const auto *ksg0_449 = buffer.data(ksg0 + 449);
    const auto *ksg0_450 = buffer.data(ksg0 + 450);

    const auto *ksg1_405 = buffer.data(ksg1 + 405);
    const auto *ksg1_406 = buffer.data(ksg1 + 406);
    const auto *ksg1_407 = buffer.data(ksg1 + 407);
    const auto *ksg1_408 = buffer.data(ksg1 + 408);
    const auto *ksg1_409 = buffer.data(ksg1 + 409);
    const auto *ksg1_410 = buffer.data(ksg1 + 410);
    const auto *ksg1_420 = buffer.data(ksg1 + 420);
    const auto *ksg1_421 = buffer.data(ksg1 + 421);
    const auto *ksg1_423 = buffer.data(ksg1 + 423);
    const auto *ksg1_425 = buffer.data(ksg1 + 425);
    const auto *ksg1_426 = buffer.data(ksg1 + 426);
    const auto *ksg1_428 = buffer.data(ksg1 + 428);
    const auto *ksg1_429 = buffer.data(ksg1 + 429);
    const auto *ksg1_430 = buffer.data(ksg1 + 430);
    const auto *ksg1_431 = buffer.data(ksg1 + 431);
    const auto *ksg1_432 = buffer.data(ksg1 + 432);
    const auto *ksg1_433 = buffer.data(ksg1 + 433);
    const auto *ksg1_434 = buffer.data(ksg1 + 434);
    const auto *ksg1_437 = buffer.data(ksg1 + 437);
    const auto *ksg1_439 = buffer.data(ksg1 + 439);
    const auto *ksg1_440 = buffer.data(ksg1 + 440);
    const auto *ksg1_442 = buffer.data(ksg1 + 442);
    const auto *ksg1_443 = buffer.data(ksg1 + 443);
    const auto *ksg1_444 = buffer.data(ksg1 + 444);
    const auto *ksg1_446 = buffer.data(ksg1 + 446);
    const auto *ksg1_447 = buffer.data(ksg1 + 447);
    const auto *ksg1_448 = buffer.data(ksg1 + 448);
    const auto *ksg1_449 = buffer.data(ksg1 + 449);
    const auto *ksg1_450 = buffer.data(ksg1 + 450);

    const auto *ksh_540 = buffer.data(ksh + 540);
    const auto *ksh_541 = buffer.data(ksh + 541);
    const auto *ksh_542 = buffer.data(ksh + 542);
    const auto *ksh_543 = buffer.data(ksh + 543);
    const auto *ksh_544 = buffer.data(ksh + 544);
    const auto *ksh_545 = buffer.data(ksh + 545);
    const auto *ksh_546 = buffer.data(ksh + 546);
    const auto *ksh_548 = buffer.data(ksh + 548);
    const auto *ksh_549 = buffer.data(ksh + 549);
    const auto *ksh_551 = buffer.data(ksh + 551);
    const auto *ksh_552 = buffer.data(ksh + 552);
    const auto *ksh_555 = buffer.data(ksh + 555);
    const auto *ksh_561 = buffer.data(ksh + 561);
    const auto *ksh_562 = buffer.data(ksh + 562);
    const auto *ksh_563 = buffer.data(ksh + 563);
    const auto *ksh_564 = buffer.data(ksh + 564);
    const auto *ksh_565 = buffer.data(ksh + 565);
    const auto *ksh_566 = buffer.data(ksh + 566);
    const auto *ksh_567 = buffer.data(ksh + 567);
    const auto *ksh_568 = buffer.data(ksh + 568);
    const auto *ksh_569 = buffer.data(ksh + 569);
    const auto *ksh_570 = buffer.data(ksh + 570);
    const auto *ksh_571 = buffer.data(ksh + 571);
    const auto *ksh_572 = buffer.data(ksh + 572);
    const auto *ksh_573 = buffer.data(ksh + 573);
    const auto *ksh_574 = buffer.data(ksh + 574);
    const auto *ksh_575 = buffer.data(ksh + 575);
    const auto *ksh_576 = buffer.data(ksh + 576);
    const auto *ksh_581 = buffer.data(ksh + 581);
    const auto *ksh_582 = buffer.data(ksh + 582);
    const auto *ksh_583 = buffer.data(ksh + 583);
    const auto *ksh_584 = buffer.data(ksh + 584);
    const auto *ksh_585 = buffer.data(ksh + 585);
    const auto *ksh_587 = buffer.data(ksh + 587);
    const auto *ksh_588 = buffer.data(ksh + 588);
    const auto *ksh_589 = buffer.data(ksh + 589);
    const auto *ksh_591 = buffer.data(ksh + 591);
    const auto *ksh_593 = buffer.data(ksh + 593);
    const auto *ksh_594 = buffer.data(ksh + 594);
    const auto *ksh_596 = buffer.data(ksh + 596);
    const auto *ksh_597 = buffer.data(ksh + 597);
    const auto *ksh_598 = buffer.data(ksh + 598);
    const auto *ksh_600 = buffer.data(ksh + 600);
    const auto *ksh_601 = buffer.data(ksh + 601);
    const auto *ksh_602 = buffer.data(ksh + 602);
    const auto *ksh_603 = buffer.data(ksh + 603);
    const auto *ksh_604 = buffer.data(ksh + 604);
    const auto *ksh_605 = buffer.data(ksh + 605);
    const auto *ksh_606 = buffer.data(ksh + 606);
    const auto *ksh_607 = buffer.data(ksh + 607);
    const auto *ksh_608 = buffer.data(ksh + 608);
    const auto *ksh_611 = buffer.data(ksh + 611);
    const auto *ksh_613 = buffer.data(ksh + 613);
    const auto *ksh_614 = buffer.data(ksh + 614);
    const auto *ksh_616 = buffer.data(ksh + 616);
    const auto *ksh_617 = buffer.data(ksh + 617);
    const auto *ksh_618 = buffer.data(ksh + 618);
    const auto *ksh_620 = buffer.data(ksh + 620);
    const auto *ksh_621 = buffer.data(ksh + 621);
    const auto *ksh_622 = buffer.data(ksh + 622);
    const auto *ksh_623 = buffer.data(ksh + 623);
    const auto *ksh_624 = buffer.data(ksh + 624);
    const auto *ksh_625 = buffer.data(ksh + 625);
    const auto *ksh_626 = buffer.data(ksh + 626);
    const auto *ksh_627 = buffer.data(ksh + 627);
    const auto *ksh_628 = buffer.data(ksh + 628);
    const auto *ksh_629 = buffer.data(ksh + 629);
    const auto *ksh_630 = buffer.data(ksh + 630);

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pa_x, pc_x, isi0_714, ish_539, ish_540, \
                         ish_541, ish_542, isi1_714, ksh_540, ksh_541, \
                         ksh_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pa_x[k] * isi0_714[k]
                   + f_12 * ish_539[k]
                   - f_10 * pc_x[k] * isi1_714[k];

        t_715[k] = f_11 * ish_540[k]
                   + f_3 * pc_x[k] * ksh_540[k];

        t_716[k] = f_11 * ish_541[k]
                   + f_3 * pc_x[k] * ksh_541[k];

        t_717[k] = f_11 * ish_542[k]
                   + f_3 * pc_x[k] * ksh_542[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, t_721, pa_x, pc_x, isi0_721, ish_543, ish_544, \
                         ish_545, isi1_721, ksh_543, ksh_544, ksh_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_11 * ish_543[k]
                   + f_3 * pc_x[k] * ksh_543[k];

        t_719[k] = f_11 * ish_544[k]
                   + f_3 * pc_x[k] * ksh_544[k];

        t_720[k] = f_11 * ish_545[k]
                   + f_3 * pc_x[k] * ksh_545[k];

        t_721[k] = pa_x[k] * isi0_721[k]
                   - f_10 * pc_x[k] * isi1_721[k];
    }

#pragma omp simd aligned(t_722, t_723, t_724, t_725, pa_x, pc_x, pc_z, isi0_723, isi0_724, \
                         isi0_725, ish_393, isi1_723, isi1_724, isi1_725, \
                         ksh_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_722[k] = f_14 * ish_393[k]
                   + f_3 * pc_z[k] * ksh_540[k];

        t_723[k] = pa_x[k] * isi0_723[k]
                   - f_10 * pc_x[k] * isi1_723[k];

        t_724[k] = pa_x[k] * isi0_724[k]
                   - f_10 * pc_x[k] * isi1_724[k];

        t_725[k] = pa_x[k] * isi0_725[k]
                   - f_10 * pc_x[k] * isi1_725[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pa_x, pa_y, pc_x, pc_y, isi0_560, \
                         isi0_727, ish_419, ish_420, isi1_560, isi1_727, ksh_545, \
                         ksh_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_12 * ish_419[k]
                   + f_3 * pc_y[k] * ksh_545[k];

        t_727[k] = pa_x[k] * isi0_727[k]
                   - f_10 * pc_x[k] * isi1_727[k];

        t_728[k] = pa_y[k] * isi0_560[k]
                   - f_10 * pc_y[k] * isi1_560[k];

        t_729[k] = f_11 * ish_420[k]
                   + f_3 * pc_y[k] * ksh_546[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pa_x, pc_x, pc_y, pc_z, isi0_731, ish_399, \
                         ish_422, ish_549, isi1_731, ksh_546, ksh_548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_18 * ish_399[k]
                   + f_3 * pc_z[k] * ksh_546[k];

        t_731[k] = pa_x[k] * isi0_731[k]
                   + f_14 * ish_549[k]
                   - f_10 * pc_x[k] * isi1_731[k];

        t_732[k] = f_11 * ish_422[k]
                   + f_3 * pc_y[k] * ksh_548[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pa_x, pa_y, pc_x, pc_y, pc_z, isi0_565, \
                         isi0_734, ish_402, ish_552, isi1_565, isi1_734, \
                         ksh_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = pa_y[k] * isi0_565[k]
                   - f_10 * pc_y[k] * isi1_565[k];

        t_734[k] = pa_x[k] * isi0_734[k]
                   + f_13 * ish_552[k]
                   - f_10 * pc_x[k] * isi1_734[k];

        t_735[k] = f_18 * ish_402[k]
                   + f_3 * pc_z[k] * ksh_549[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pa_x, pa_y, pc_x, pc_y, isi0_569, isi0_738, \
                         ish_425, ish_556, isi1_569, isi1_738, \
                         ksh_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_11 * ish_425[k]
                   + f_3 * pc_y[k] * ksh_551[k];

        t_737[k] = pa_y[k] * isi0_569[k]
                   - f_10 * pc_y[k] * isi1_569[k];

        t_738[k] = pa_x[k] * isi0_738[k]
                   + f_12 * ish_556[k]
                   - f_10 * pc_x[k] * isi1_738[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pa_x, pc_x, pc_y, pc_z, isi0_740, ish_405, \
                         ish_429, ish_558, isi1_740, ksh_552, ksh_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_18 * ish_405[k]
                   + f_3 * pc_z[k] * ksh_552[k];

        t_740[k] = pa_x[k] * isi0_740[k]
                   + f_12 * ish_558[k]
                   - f_10 * pc_x[k] * isi1_740[k];

        t_741[k] = f_11 * ish_429[k]
                   + f_3 * pc_y[k] * ksh_555[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pa_y, pc_x, pc_y, isi0_574, ish_561, \
                         ish_562, ish_563, isi1_574, ksh_561, ksh_562, \
                         ksh_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = pa_y[k] * isi0_574[k]
                   - f_10 * pc_y[k] * isi1_574[k];

        t_743[k] = f_11 * ish_561[k]
                   + f_3 * pc_x[k] * ksh_561[k];

        t_744[k] = f_11 * ish_562[k]
                   + f_3 * pc_x[k] * ksh_562[k];

        t_745[k] = f_11 * ish_563[k]
                   + f_3 * pc_x[k] * ksh_563[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pa_x, pc_x, isi0_749, ish_564, ish_565, \
                         ish_566, isi1_749, ksh_564, ksh_565, ksh_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_11 * ish_564[k]
                   + f_3 * pc_x[k] * ksh_564[k];

        t_747[k] = f_11 * ish_565[k]
                   + f_3 * pc_x[k] * ksh_565[k];

        t_748[k] = f_11 * ish_566[k]
                   + f_3 * pc_x[k] * ksh_566[k];

        t_749[k] = pa_x[k] * isi0_749[k]
                   - f_10 * pc_x[k] * isi1_749[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, pa_x, pc_x, pc_z, isi0_751, isi0_752, \
                         isi0_753, ish_414, isi1_751, isi1_752, isi1_753, \
                         ksh_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_18 * ish_414[k]
                   + f_3 * pc_z[k] * ksh_561[k];

        t_751[k] = pa_x[k] * isi0_751[k]
                   - f_10 * pc_x[k] * isi1_751[k];

        t_752[k] = pa_x[k] * isi0_752[k]
                   - f_10 * pc_x[k] * isi1_752[k];

        t_753[k] = pa_x[k] * isi0_753[k]
                   - f_10 * pc_x[k] * isi1_753[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, t_757, pa_x, pc_x, pc_y, isi0_755, isi0_756, \
                         ish_440, ish_567, isi1_755, isi1_756, ksh_566, \
                         ksh_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_11 * ish_440[k]
                   + f_3 * pc_y[k] * ksh_566[k];

        t_755[k] = pa_x[k] * isi0_755[k]
                   - f_10 * pc_x[k] * isi1_755[k];

        t_756[k] = pa_x[k] * isi0_756[k]
                   + f_15 * ish_567[k]
                   - f_10 * pc_x[k] * isi1_756[k];

        t_757[k] = f_3 * pc_y[k] * ksh_567[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, pc_z, ish_420, ksg0_405, ksg1_405, \
                         ksh_567, ksh_568, ksh_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_15 * ish_420[k]
                   + f_3 * pc_z[k] * ksh_567[k];

        t_759[k] = f_4 * ksg0_405[k]
                   - f_5 * ksg1_405[k]
                   + f_3 * pc_y[k] * ksh_568[k];

        t_760[k] = f_3 * pc_y[k] * ksh_569[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pa_x, pc_x, pc_y, isi0_761, ish_572, isi1_761, \
                         ksg0_406, ksg0_407, ksg1_406, ksg1_407, ksh_570, \
                         ksh_571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = pa_x[k] * isi0_761[k]
                   + f_14 * ish_572[k]
                   - f_10 * pc_x[k] * isi1_761[k];

        t_762[k] = f_6 * ksg0_406[k]
                   - f_7 * ksg1_406[k]
                   + f_3 * pc_y[k] * ksh_570[k];

        t_763[k] = f_4 * ksg0_407[k]
                   - f_5 * ksg1_407[k]
                   + f_3 * pc_y[k] * ksh_571[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pa_x, pc_x, pc_y, isi0_765, ish_576, isi1_765, \
                         ksg0_408, ksg1_408, ksh_572, ksh_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_3 * pc_y[k] * ksh_572[k];

        t_765[k] = pa_x[k] * isi0_765[k]
                   + f_13 * ish_576[k]
                   - f_10 * pc_x[k] * isi1_765[k];

        t_766[k] = f_8 * ksg0_408[k]
                   - f_9 * ksg1_408[k]
                   + f_3 * pc_y[k] * ksh_573[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_y, ksg0_409, ksg0_410, ksg1_409, ksg1_410, \
                         ksh_574, ksh_575, ksh_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_6 * ksg0_409[k]
                   - f_7 * ksg1_409[k]
                   + f_3 * pc_y[k] * ksh_574[k];

        t_768[k] = f_4 * ksg0_410[k]
                   - f_5 * ksg1_410[k]
                   + f_3 * pc_y[k] * ksh_575[k];

        t_769[k] = f_3 * pc_y[k] * ksh_576[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pa_x, pc_x, isi0_770, ish_581, ish_582, \
                         ish_583, ish_584, isi1_770, ksh_582, ksh_583, \
                         ksh_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pa_x[k] * isi0_770[k]
                   + f_12 * ish_581[k]
                   - f_10 * pc_x[k] * isi1_770[k];

        t_771[k] = f_11 * ish_582[k]
                   + f_3 * pc_x[k] * ksh_582[k];

        t_772[k] = f_11 * ish_583[k]
                   + f_3 * pc_x[k] * ksh_583[k];

        t_773[k] = f_11 * ish_584[k]
                   + f_3 * pc_x[k] * ksh_584[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, t_777, pa_x, pc_x, pc_y, isi0_777, ish_585, \
                         ish_587, isi1_777, ksh_581, ksh_585, ksh_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_11 * ish_585[k]
                   + f_3 * pc_x[k] * ksh_585[k];

        t_775[k] = f_3 * pc_y[k] * ksh_581[k];

        t_776[k] = f_11 * ish_587[k]
                   + f_3 * pc_x[k] * ksh_587[k];

        t_777[k] = pa_x[k] * isi0_777[k]
                   - f_10 * pc_x[k] * isi1_777[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, t_781, pa_x, pc_x, isi0_778, isi0_779, isi0_780, \
                         isi0_781, isi1_778, isi1_779, isi1_780, \
                         isi1_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = pa_x[k] * isi0_778[k]
                   - f_10 * pc_x[k] * isi1_778[k];

        t_779[k] = pa_x[k] * isi0_779[k]
                   - f_10 * pc_x[k] * isi1_779[k];

        t_780[k] = pa_x[k] * isi0_780[k]
                   - f_10 * pc_x[k] * isi1_780[k];

        t_781[k] = pa_x[k] * isi0_781[k]
                   - f_10 * pc_x[k] * isi1_781[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, t_785, pa_x, pc_x, pc_y, isi0_783, isi1_783, \
                         ksg0_420, ksg0_421, ksg1_420, ksg1_421, ksh_587, ksh_588, \
                         ksh_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_3 * pc_y[k] * ksh_587[k];

        t_783[k] = pa_x[k] * isi0_783[k]
                   - f_10 * pc_x[k] * isi1_783[k];

        t_784[k] = f_1 * ksg0_420[k]
                   - f_2 * ksg1_420[k]
                   + f_3 * pc_x[k] * ksh_588[k];

        t_785[k] = f_16 * ksg0_421[k]
                   - f_17 * ksg1_421[k]
                   + f_3 * pc_x[k] * ksh_589[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, pc_z, ksg0_423, ksg0_425, ksg1_423, \
                         ksg1_425, ksh_588, ksh_589, ksh_591, ksh_593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_3 * pc_z[k] * ksh_588[k];

        t_787[k] = f_8 * ksg0_423[k]
                   - f_9 * ksg1_423[k]
                   + f_3 * pc_x[k] * ksh_591[k];

        t_788[k] = f_3 * pc_z[k] * ksh_589[k];

        t_789[k] = f_8 * ksg0_425[k]
                   - f_9 * ksg1_425[k]
                   + f_3 * pc_x[k] * ksh_593[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, t_793, pc_x, pc_z, ksg0_426, ksg0_428, ksg0_429, \
                         ksg1_426, ksg1_428, ksg1_429, ksh_591, ksh_594, ksh_596, \
                         ksh_597 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_6 * ksg0_426[k]
                   - f_7 * ksg1_426[k]
                   + f_3 * pc_x[k] * ksh_594[k];

        t_791[k] = f_3 * pc_z[k] * ksh_591[k];

        t_792[k] = f_6 * ksg0_428[k]
                   - f_7 * ksg1_428[k]
                   + f_3 * pc_x[k] * ksh_596[k];

        t_793[k] = f_6 * ksg0_429[k]
                   - f_7 * ksg1_429[k]
                   + f_3 * pc_x[k] * ksh_597[k];
    }

#pragma omp simd aligned(t_794, t_795, t_796, t_797, pc_x, pc_z, ksg0_430, ksg0_432, ksg0_433, \
                         ksg1_430, ksg1_432, ksg1_433, ksh_594, ksh_598, ksh_600, \
                         ksh_601 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_794[k] = f_4 * ksg0_430[k]
                   - f_5 * ksg1_430[k]
                   + f_3 * pc_x[k] * ksh_598[k];

        t_795[k] = f_3 * pc_z[k] * ksh_594[k];

        t_796[k] = f_4 * ksg0_432[k]
                   - f_5 * ksg1_432[k]
                   + f_3 * pc_x[k] * ksh_600[k];

        t_797[k] = f_4 * ksg0_433[k]
                   - f_5 * ksg1_433[k]
                   + f_3 * pc_x[k] * ksh_601[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, t_803, pc_x, ksg0_434, ksg1_434, \
                         ksh_602, ksh_603, ksh_604, ksh_605, ksh_606, \
                         ksh_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = f_4 * ksg0_434[k]
                   - f_5 * ksg1_434[k]
                   + f_3 * pc_x[k] * ksh_602[k];

        t_799[k] = f_3 * pc_x[k] * ksh_603[k];

        t_800[k] = f_3 * pc_x[k] * ksh_604[k];

        t_801[k] = f_3 * pc_x[k] * ksh_605[k];

        t_802[k] = f_3 * pc_x[k] * ksh_606[k];

        t_803[k] = f_3 * pc_x[k] * ksh_607[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pc_x, pc_y, pc_z, ish_456, ksg0_430, \
                         ksg1_430, ksh_603, ksh_604, ksh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_3 * pc_x[k] * ksh_608[k];

        t_805[k] = f_0 * ish_456[k]
                   + f_1 * ksg0_430[k]
                   - f_2 * ksg1_430[k]
                   + f_3 * pc_y[k] * ksh_603[k];

        t_806[k] = f_3 * pc_z[k] * ksh_603[k];

        t_807[k] = f_4 * ksg0_430[k]
                   - f_5 * ksg1_430[k]
                   + f_3 * pc_z[k] * ksh_604[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pc_y, pc_z, ish_461, ksg0_431, ksg0_432, \
                         ksg0_434, ksg1_431, ksg1_432, ksg1_434, ksh_605, ksh_606, \
                         ksh_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_6 * ksg0_431[k]
                   - f_7 * ksg1_431[k]
                   + f_3 * pc_z[k] * ksh_605[k];

        t_809[k] = f_8 * ksg0_432[k]
                   - f_9 * ksg1_432[k]
                   + f_3 * pc_z[k] * ksh_606[k];

        t_810[k] = f_0 * ish_461[k]
                   + f_3 * pc_y[k] * ksh_608[k];

        t_811[k] = f_1 * ksg0_434[k]
                   - f_2 * ksg1_434[k]
                   + f_3 * pc_z[k] * ksh_608[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, t_815, pa_z, pc_x, pc_z, isi0_588, isi0_589, \
                         isi0_591, isi1_588, isi1_589, isi1_591, ksg0_437, ksg1_437, \
                         ksh_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = pa_z[k] * isi0_588[k]
                   - f_10 * pc_z[k] * isi1_588[k];

        t_813[k] = pa_z[k] * isi0_589[k]
                   - f_10 * pc_z[k] * isi1_589[k];

        t_814[k] = f_16 * ksg0_437[k]
                   - f_17 * ksg1_437[k]
                   + f_3 * pc_x[k] * ksh_611[k];

        t_815[k] = pa_z[k] * isi0_591[k]
                   - f_10 * pc_z[k] * isi1_591[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pa_z, pc_x, pc_z, isi0_594, isi1_594, ksg0_439, \
                         ksg0_440, ksg1_439, ksg1_440, ksh_613, \
                         ksh_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_8 * ksg0_439[k]
                   - f_9 * ksg1_439[k]
                   + f_3 * pc_x[k] * ksh_613[k];

        t_817[k] = f_8 * ksg0_440[k]
                   - f_9 * ksg1_440[k]
                   + f_3 * pc_x[k] * ksh_614[k];

        t_818[k] = pa_z[k] * isi0_594[k]
                   - f_10 * pc_z[k] * isi1_594[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, ksg0_442, ksg0_443, ksg0_444, ksg1_442, \
                         ksg1_443, ksg1_444, ksh_616, ksh_617, \
                         ksh_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_6 * ksg0_442[k]
                   - f_7 * ksg1_442[k]
                   + f_3 * pc_x[k] * ksh_616[k];

        t_820[k] = f_6 * ksg0_443[k]
                   - f_7 * ksg1_443[k]
                   + f_3 * pc_x[k] * ksh_617[k];

        t_821[k] = f_6 * ksg0_444[k]
                   - f_7 * ksg1_444[k]
                   + f_3 * pc_x[k] * ksh_618[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, pa_z, pc_x, pc_z, isi0_598, isi1_598, ksg0_446, \
                         ksg0_447, ksg1_446, ksg1_447, ksh_620, \
                         ksh_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = pa_z[k] * isi0_598[k]
                   - f_10 * pc_z[k] * isi1_598[k];

        t_823[k] = f_4 * ksg0_446[k]
                   - f_5 * ksg1_446[k]
                   + f_3 * pc_x[k] * ksh_620[k];

        t_824[k] = f_4 * ksg0_447[k]
                   - f_5 * ksg1_447[k]
                   + f_3 * pc_x[k] * ksh_621[k];
    }

#pragma omp simd aligned(t_825, t_826, t_827, t_828, t_829, pc_x, ksg0_448, ksg0_449, \
                         ksg1_448, ksg1_449, ksh_622, ksh_623, ksh_624, ksh_625, \
                         ksh_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_825[k] = f_4 * ksg0_448[k]
                   - f_5 * ksg1_448[k]
                   + f_3 * pc_x[k] * ksh_622[k];

        t_826[k] = f_4 * ksg0_449[k]
                   - f_5 * ksg1_449[k]
                   + f_3 * pc_x[k] * ksh_623[k];

        t_827[k] = f_3 * pc_x[k] * ksh_624[k];

        t_828[k] = f_3 * pc_x[k] * ksh_625[k];

        t_829[k] = f_3 * pc_x[k] * ksh_626[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, t_834, pa_z, pc_x, pc_z, isi0_609, \
                         ish_456, isi1_609, ksh_624, ksh_627, ksh_628, \
                         ksh_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_3 * pc_x[k] * ksh_627[k];

        t_831[k] = f_3 * pc_x[k] * ksh_628[k];

        t_832[k] = f_3 * pc_x[k] * ksh_629[k];

        t_833[k] = pa_z[k] * isi0_609[k]
                   - f_10 * pc_z[k] * isi1_609[k];

        t_834[k] = f_11 * ish_456[k]
                   + f_3 * pc_z[k] * ksh_624[k];
    }

#pragma omp simd aligned(t_835, t_836, t_837, pa_z, pc_z, isi0_611, isi0_612, isi0_613, \
                         ish_457, ish_458, ish_459, isi1_611, isi1_612, \
                         isi1_613 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_835[k] = pa_z[k] * isi0_611[k]
                   + f_12 * ish_457[k]
                   - f_10 * pc_z[k] * isi1_611[k];

        t_836[k] = pa_z[k] * isi0_612[k]
                   + f_13 * ish_458[k]
                   - f_10 * pc_z[k] * isi1_612[k];

        t_837[k] = pa_z[k] * isi0_613[k]
                   + f_14 * ish_459[k]
                   - f_10 * pc_z[k] * isi1_613[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pc_x, pc_y, pc_z, ish_461, ish_482, ksg0_449, \
                         ksg0_450, ksg1_449, ksg1_450, ksh_629, \
                         ksh_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_15 * ish_482[k]
                   + f_3 * pc_y[k] * ksh_629[k];

        t_839[k] = f_11 * ish_461[k]
                   + f_1 * ksg0_449[k]
                   - f_2 * ksg1_449[k]
                   + f_3 * pc_z[k] * ksh_629[k];

        t_840[k] = f_1 * ksg0_450[k]
                   - f_2 * ksg1_450[k]
                   + f_3 * pc_x[k] * ksh_630[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
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
    const auto f_12 = 1.0 / q;
    const auto f_13 = 1.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);
    const auto f_18 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isi0_756 = buffer.data(isi0 + 756);
    const auto *isi0_758 = buffer.data(isi0 + 758);

    const auto *ish_477 = buffer.data(ish + 477);
    const auto *ish_482 = buffer.data(ish + 482);
    const auto *ish_498 = buffer.data(ish + 498);
    const auto *ish_500 = buffer.data(ish + 500);
    const auto *ish_501 = buffer.data(ish + 501);
    const auto *ish_502 = buffer.data(ish + 502);
    const auto *ish_503 = buffer.data(ish + 503);
    const auto *ish_519 = buffer.data(ish + 519);
    const auto *ish_521 = buffer.data(ish + 521);
    const auto *ish_522 = buffer.data(ish + 522);
    const auto *ish_523 = buffer.data(ish + 523);
    const auto *ish_524 = buffer.data(ish + 524);
    const auto *ish_540 = buffer.data(ish + 540);
    const auto *ish_542 = buffer.data(ish + 542);
    const auto *ish_543 = buffer.data(ish + 543);
    const auto *ish_544 = buffer.data(ish + 544);
    const auto *ish_545 = buffer.data(ish + 545);
    const auto *ish_561 = buffer.data(ish + 561);
    const auto *ish_563 = buffer.data(ish + 563);
    const auto *ish_564 = buffer.data(ish + 564);
    const auto *ish_565 = buffer.data(ish + 565);
    const auto *ish_566 = buffer.data(ish + 566);

    const auto *isi1_756 = buffer.data(isi1 + 756);
    const auto *isi1_758 = buffer.data(isi1 + 758);

    const auto *ksg0_451 = buffer.data(ksg0 + 451);
    const auto *ksg0_452 = buffer.data(ksg0 + 452);
    const auto *ksg0_453 = buffer.data(ksg0 + 453);
    const auto *ksg0_454 = buffer.data(ksg0 + 454);
    const auto *ksg0_455 = buffer.data(ksg0 + 455);
    const auto *ksg0_456 = buffer.data(ksg0 + 456);
    const auto *ksg0_457 = buffer.data(ksg0 + 457);
    const auto *ksg0_458 = buffer.data(ksg0 + 458);
    const auto *ksg0_459 = buffer.data(ksg0 + 459);
    const auto *ksg0_460 = buffer.data(ksg0 + 460);
    const auto *ksg0_461 = buffer.data(ksg0 + 461);
    const auto *ksg0_462 = buffer.data(ksg0 + 462);
    const auto *ksg0_463 = buffer.data(ksg0 + 463);
    const auto *ksg0_464 = buffer.data(ksg0 + 464);
    const auto *ksg0_465 = buffer.data(ksg0 + 465);
    const auto *ksg0_466 = buffer.data(ksg0 + 466);
    const auto *ksg0_467 = buffer.data(ksg0 + 467);
    const auto *ksg0_468 = buffer.data(ksg0 + 468);
    const auto *ksg0_469 = buffer.data(ksg0 + 469);
    const auto *ksg0_470 = buffer.data(ksg0 + 470);
    const auto *ksg0_471 = buffer.data(ksg0 + 471);
    const auto *ksg0_472 = buffer.data(ksg0 + 472);
    const auto *ksg0_473 = buffer.data(ksg0 + 473);
    const auto *ksg0_474 = buffer.data(ksg0 + 474);
    const auto *ksg0_475 = buffer.data(ksg0 + 475);
    const auto *ksg0_476 = buffer.data(ksg0 + 476);
    const auto *ksg0_477 = buffer.data(ksg0 + 477);
    const auto *ksg0_478 = buffer.data(ksg0 + 478);
    const auto *ksg0_479 = buffer.data(ksg0 + 479);
    const auto *ksg0_480 = buffer.data(ksg0 + 480);
    const auto *ksg0_481 = buffer.data(ksg0 + 481);
    const auto *ksg0_482 = buffer.data(ksg0 + 482);
    const auto *ksg0_483 = buffer.data(ksg0 + 483);
    const auto *ksg0_484 = buffer.data(ksg0 + 484);
    const auto *ksg0_485 = buffer.data(ksg0 + 485);
    const auto *ksg0_486 = buffer.data(ksg0 + 486);
    const auto *ksg0_487 = buffer.data(ksg0 + 487);
    const auto *ksg0_488 = buffer.data(ksg0 + 488);
    const auto *ksg0_489 = buffer.data(ksg0 + 489);
    const auto *ksg0_490 = buffer.data(ksg0 + 490);
    const auto *ksg0_491 = buffer.data(ksg0 + 491);
    const auto *ksg0_492 = buffer.data(ksg0 + 492);
    const auto *ksg0_493 = buffer.data(ksg0 + 493);
    const auto *ksg0_494 = buffer.data(ksg0 + 494);
    const auto *ksg0_495 = buffer.data(ksg0 + 495);
    const auto *ksg0_496 = buffer.data(ksg0 + 496);
    const auto *ksg0_497 = buffer.data(ksg0 + 497);
    const auto *ksg0_498 = buffer.data(ksg0 + 498);
    const auto *ksg0_499 = buffer.data(ksg0 + 499);
    const auto *ksg0_500 = buffer.data(ksg0 + 500);
    const auto *ksg0_501 = buffer.data(ksg0 + 501);
    const auto *ksg0_502 = buffer.data(ksg0 + 502);
    const auto *ksg0_503 = buffer.data(ksg0 + 503);
    const auto *ksg0_504 = buffer.data(ksg0 + 504);
    const auto *ksg0_505 = buffer.data(ksg0 + 505);
    const auto *ksg0_506 = buffer.data(ksg0 + 506);
    const auto *ksg0_507 = buffer.data(ksg0 + 507);
    const auto *ksg0_508 = buffer.data(ksg0 + 508);
    const auto *ksg0_509 = buffer.data(ksg0 + 509);
    const auto *ksg0_511 = buffer.data(ksg0 + 511);
    const auto *ksg0_513 = buffer.data(ksg0 + 513);

    const auto *ksg1_451 = buffer.data(ksg1 + 451);
    const auto *ksg1_452 = buffer.data(ksg1 + 452);
    const auto *ksg1_453 = buffer.data(ksg1 + 453);
    const auto *ksg1_454 = buffer.data(ksg1 + 454);
    const auto *ksg1_455 = buffer.data(ksg1 + 455);
    const auto *ksg1_456 = buffer.data(ksg1 + 456);
    const auto *ksg1_457 = buffer.data(ksg1 + 457);
    const auto *ksg1_458 = buffer.data(ksg1 + 458);
    const auto *ksg1_459 = buffer.data(ksg1 + 459);
    const auto *ksg1_460 = buffer.data(ksg1 + 460);
    const auto *ksg1_461 = buffer.data(ksg1 + 461);
    const auto *ksg1_462 = buffer.data(ksg1 + 462);
    const auto *ksg1_463 = buffer.data(ksg1 + 463);
    const auto *ksg1_464 = buffer.data(ksg1 + 464);
    const auto *ksg1_465 = buffer.data(ksg1 + 465);
    const auto *ksg1_466 = buffer.data(ksg1 + 466);
    const auto *ksg1_467 = buffer.data(ksg1 + 467);
    const auto *ksg1_468 = buffer.data(ksg1 + 468);
    const auto *ksg1_469 = buffer.data(ksg1 + 469);
    const auto *ksg1_470 = buffer.data(ksg1 + 470);
    const auto *ksg1_471 = buffer.data(ksg1 + 471);
    const auto *ksg1_472 = buffer.data(ksg1 + 472);
    const auto *ksg1_473 = buffer.data(ksg1 + 473);
    const auto *ksg1_474 = buffer.data(ksg1 + 474);
    const auto *ksg1_475 = buffer.data(ksg1 + 475);
    const auto *ksg1_476 = buffer.data(ksg1 + 476);
    const auto *ksg1_477 = buffer.data(ksg1 + 477);
    const auto *ksg1_478 = buffer.data(ksg1 + 478);
    const auto *ksg1_479 = buffer.data(ksg1 + 479);
    const auto *ksg1_480 = buffer.data(ksg1 + 480);
    const auto *ksg1_481 = buffer.data(ksg1 + 481);
    const auto *ksg1_482 = buffer.data(ksg1 + 482);
    const auto *ksg1_483 = buffer.data(ksg1 + 483);
    const auto *ksg1_484 = buffer.data(ksg1 + 484);
    const auto *ksg1_485 = buffer.data(ksg1 + 485);
    const auto *ksg1_486 = buffer.data(ksg1 + 486);
    const auto *ksg1_487 = buffer.data(ksg1 + 487);
    const auto *ksg1_488 = buffer.data(ksg1 + 488);
    const auto *ksg1_489 = buffer.data(ksg1 + 489);
    const auto *ksg1_490 = buffer.data(ksg1 + 490);
    const auto *ksg1_491 = buffer.data(ksg1 + 491);
    const auto *ksg1_492 = buffer.data(ksg1 + 492);
    const auto *ksg1_493 = buffer.data(ksg1 + 493);
    const auto *ksg1_494 = buffer.data(ksg1 + 494);
    const auto *ksg1_495 = buffer.data(ksg1 + 495);
    const auto *ksg1_496 = buffer.data(ksg1 + 496);
    const auto *ksg1_497 = buffer.data(ksg1 + 497);
    const auto *ksg1_498 = buffer.data(ksg1 + 498);
    const auto *ksg1_499 = buffer.data(ksg1 + 499);
    const auto *ksg1_500 = buffer.data(ksg1 + 500);
    const auto *ksg1_501 = buffer.data(ksg1 + 501);
    const auto *ksg1_502 = buffer.data(ksg1 + 502);
    const auto *ksg1_503 = buffer.data(ksg1 + 503);
    const auto *ksg1_504 = buffer.data(ksg1 + 504);
    const auto *ksg1_505 = buffer.data(ksg1 + 505);
    const auto *ksg1_506 = buffer.data(ksg1 + 506);
    const auto *ksg1_507 = buffer.data(ksg1 + 507);
    const auto *ksg1_508 = buffer.data(ksg1 + 508);
    const auto *ksg1_509 = buffer.data(ksg1 + 509);
    const auto *ksg1_511 = buffer.data(ksg1 + 511);
    const auto *ksg1_513 = buffer.data(ksg1 + 513);

    const auto *ksh_631 = buffer.data(ksh + 631);
    const auto *ksh_632 = buffer.data(ksh + 632);
    const auto *ksh_633 = buffer.data(ksh + 633);
    const auto *ksh_634 = buffer.data(ksh + 634);
    const auto *ksh_635 = buffer.data(ksh + 635);
    const auto *ksh_636 = buffer.data(ksh + 636);
    const auto *ksh_637 = buffer.data(ksh + 637);
    const auto *ksh_638 = buffer.data(ksh + 638);
    const auto *ksh_639 = buffer.data(ksh + 639);
    const auto *ksh_640 = buffer.data(ksh + 640);
    const auto *ksh_641 = buffer.data(ksh + 641);
    const auto *ksh_642 = buffer.data(ksh + 642);
    const auto *ksh_643 = buffer.data(ksh + 643);
    const auto *ksh_644 = buffer.data(ksh + 644);
    const auto *ksh_645 = buffer.data(ksh + 645);
    const auto *ksh_646 = buffer.data(ksh + 646);
    const auto *ksh_647 = buffer.data(ksh + 647);
    const auto *ksh_648 = buffer.data(ksh + 648);
    const auto *ksh_649 = buffer.data(ksh + 649);
    const auto *ksh_650 = buffer.data(ksh + 650);
    const auto *ksh_651 = buffer.data(ksh + 651);
    const auto *ksh_652 = buffer.data(ksh + 652);
    const auto *ksh_653 = buffer.data(ksh + 653);
    const auto *ksh_654 = buffer.data(ksh + 654);
    const auto *ksh_655 = buffer.data(ksh + 655);
    const auto *ksh_656 = buffer.data(ksh + 656);
    const auto *ksh_657 = buffer.data(ksh + 657);
    const auto *ksh_658 = buffer.data(ksh + 658);
    const auto *ksh_659 = buffer.data(ksh + 659);
    const auto *ksh_660 = buffer.data(ksh + 660);
    const auto *ksh_661 = buffer.data(ksh + 661);
    const auto *ksh_662 = buffer.data(ksh + 662);
    const auto *ksh_663 = buffer.data(ksh + 663);
    const auto *ksh_664 = buffer.data(ksh + 664);
    const auto *ksh_665 = buffer.data(ksh + 665);
    const auto *ksh_666 = buffer.data(ksh + 666);
    const auto *ksh_667 = buffer.data(ksh + 667);
    const auto *ksh_668 = buffer.data(ksh + 668);
    const auto *ksh_669 = buffer.data(ksh + 669);
    const auto *ksh_670 = buffer.data(ksh + 670);
    const auto *ksh_671 = buffer.data(ksh + 671);
    const auto *ksh_672 = buffer.data(ksh + 672);
    const auto *ksh_673 = buffer.data(ksh + 673);
    const auto *ksh_674 = buffer.data(ksh + 674);
    const auto *ksh_675 = buffer.data(ksh + 675);
    const auto *ksh_676 = buffer.data(ksh + 676);
    const auto *ksh_677 = buffer.data(ksh + 677);
    const auto *ksh_678 = buffer.data(ksh + 678);
    const auto *ksh_679 = buffer.data(ksh + 679);
    const auto *ksh_680 = buffer.data(ksh + 680);
    const auto *ksh_681 = buffer.data(ksh + 681);
    const auto *ksh_682 = buffer.data(ksh + 682);
    const auto *ksh_683 = buffer.data(ksh + 683);
    const auto *ksh_684 = buffer.data(ksh + 684);
    const auto *ksh_685 = buffer.data(ksh + 685);
    const auto *ksh_686 = buffer.data(ksh + 686);
    const auto *ksh_687 = buffer.data(ksh + 687);
    const auto *ksh_688 = buffer.data(ksh + 688);
    const auto *ksh_689 = buffer.data(ksh + 689);
    const auto *ksh_690 = buffer.data(ksh + 690);
    const auto *ksh_691 = buffer.data(ksh + 691);
    const auto *ksh_692 = buffer.data(ksh + 692);
    const auto *ksh_693 = buffer.data(ksh + 693);
    const auto *ksh_694 = buffer.data(ksh + 694);
    const auto *ksh_695 = buffer.data(ksh + 695);
    const auto *ksh_696 = buffer.data(ksh + 696);
    const auto *ksh_697 = buffer.data(ksh + 697);
    const auto *ksh_698 = buffer.data(ksh + 698);
    const auto *ksh_699 = buffer.data(ksh + 699);
    const auto *ksh_700 = buffer.data(ksh + 700);
    const auto *ksh_701 = buffer.data(ksh + 701);
    const auto *ksh_702 = buffer.data(ksh + 702);
    const auto *ksh_703 = buffer.data(ksh + 703);
    const auto *ksh_704 = buffer.data(ksh + 704);
    const auto *ksh_705 = buffer.data(ksh + 705);
    const auto *ksh_706 = buffer.data(ksh + 706);
    const auto *ksh_707 = buffer.data(ksh + 707);
    const auto *ksh_708 = buffer.data(ksh + 708);
    const auto *ksh_709 = buffer.data(ksh + 709);
    const auto *ksh_710 = buffer.data(ksh + 710);
    const auto *ksh_711 = buffer.data(ksh + 711);
    const auto *ksh_712 = buffer.data(ksh + 712);
    const auto *ksh_713 = buffer.data(ksh + 713);
    const auto *ksh_715 = buffer.data(ksh + 715);
    const auto *ksh_717 = buffer.data(ksh + 717);

#pragma omp simd aligned(t_841, t_842, t_843, pc_x, ksg0_451, ksg0_452, ksg0_453, ksg1_451, \
                         ksg1_452, ksg1_453, ksh_631, ksh_632, \
                         ksh_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_16 * ksg0_451[k]
                   - f_17 * ksg1_451[k]
                   + f_3 * pc_x[k] * ksh_631[k];

        t_842[k] = f_16 * ksg0_452[k]
                   - f_17 * ksg1_452[k]
                   + f_3 * pc_x[k] * ksh_632[k];

        t_843[k] = f_8 * ksg0_453[k]
                   - f_9 * ksg1_453[k]
                   + f_3 * pc_x[k] * ksh_633[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, pc_x, ksg0_454, ksg0_455, ksg0_456, ksg1_454, \
                         ksg1_455, ksg1_456, ksh_634, ksh_635, \
                         ksh_636 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = f_8 * ksg0_454[k]
                   - f_9 * ksg1_454[k]
                   + f_3 * pc_x[k] * ksh_634[k];

        t_845[k] = f_8 * ksg0_455[k]
                   - f_9 * ksg1_455[k]
                   + f_3 * pc_x[k] * ksh_635[k];

        t_846[k] = f_6 * ksg0_456[k]
                   - f_7 * ksg1_456[k]
                   + f_3 * pc_x[k] * ksh_636[k];
    }

#pragma omp simd aligned(t_847, t_848, t_849, pc_x, ksg0_457, ksg0_458, ksg0_459, ksg1_457, \
                         ksg1_458, ksg1_459, ksh_637, ksh_638, \
                         ksh_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_847[k] = f_6 * ksg0_457[k]
                   - f_7 * ksg1_457[k]
                   + f_3 * pc_x[k] * ksh_637[k];

        t_848[k] = f_6 * ksg0_458[k]
                   - f_7 * ksg1_458[k]
                   + f_3 * pc_x[k] * ksh_638[k];

        t_849[k] = f_6 * ksg0_459[k]
                   - f_7 * ksg1_459[k]
                   + f_3 * pc_x[k] * ksh_639[k];
    }

#pragma omp simd aligned(t_850, t_851, t_852, pc_x, ksg0_460, ksg0_461, ksg0_462, ksg1_460, \
                         ksg1_461, ksg1_462, ksh_640, ksh_641, \
                         ksh_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_850[k] = f_4 * ksg0_460[k]
                   - f_5 * ksg1_460[k]
                   + f_3 * pc_x[k] * ksh_640[k];

        t_851[k] = f_4 * ksg0_461[k]
                   - f_5 * ksg1_461[k]
                   + f_3 * pc_x[k] * ksh_641[k];

        t_852[k] = f_4 * ksg0_462[k]
                   - f_5 * ksg1_462[k]
                   + f_3 * pc_x[k] * ksh_642[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, t_857, pc_x, ksg0_463, ksg0_464, \
                         ksg1_463, ksg1_464, ksh_643, ksh_644, ksh_645, ksh_646, \
                         ksh_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_4 * ksg0_463[k]
                   - f_5 * ksg1_463[k]
                   + f_3 * pc_x[k] * ksh_643[k];

        t_854[k] = f_4 * ksg0_464[k]
                   - f_5 * ksg1_464[k]
                   + f_3 * pc_x[k] * ksh_644[k];

        t_855[k] = f_3 * pc_x[k] * ksh_645[k];

        t_856[k] = f_3 * pc_x[k] * ksh_646[k];

        t_857[k] = f_3 * pc_x[k] * ksh_647[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, t_862, pc_x, pc_y, pc_z, ish_477, \
                         ish_498, ksg0_460, ksg1_460, ksh_645, ksh_648, ksh_649, \
                         ksh_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_3 * pc_x[k] * ksh_648[k];

        t_859[k] = f_3 * pc_x[k] * ksh_649[k];

        t_860[k] = f_3 * pc_x[k] * ksh_650[k];

        t_861[k] = f_18 * ish_498[k]
                   + f_1 * ksg0_460[k]
                   - f_2 * ksg1_460[k]
                   + f_3 * pc_y[k] * ksh_645[k];

        t_862[k] = f_12 * ish_477[k]
                   + f_3 * pc_z[k] * ksh_645[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pc_y, ish_500, ish_501, ish_502, ksg0_462, \
                         ksg0_463, ksg0_464, ksg1_462, ksg1_463, ksg1_464, ksh_647, ksh_648, \
                         ksh_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_18 * ish_500[k]
                   + f_8 * ksg0_462[k]
                   - f_9 * ksg1_462[k]
                   + f_3 * pc_y[k] * ksh_647[k];

        t_864[k] = f_18 * ish_501[k]
                   + f_6 * ksg0_463[k]
                   - f_7 * ksg1_463[k]
                   + f_3 * pc_y[k] * ksh_648[k];

        t_865[k] = f_18 * ish_502[k]
                   + f_4 * ksg0_464[k]
                   - f_5 * ksg1_464[k]
                   + f_3 * pc_y[k] * ksh_649[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pc_x, pc_y, pc_z, ish_482, ish_503, ksg0_464, \
                         ksg0_465, ksg1_464, ksg1_465, ksh_650, \
                         ksh_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * ish_503[k]
                   + f_3 * pc_y[k] * ksh_650[k];

        t_867[k] = f_12 * ish_482[k]
                   + f_1 * ksg0_464[k]
                   - f_2 * ksg1_464[k]
                   + f_3 * pc_z[k] * ksh_650[k];

        t_868[k] = f_1 * ksg0_465[k]
                   - f_2 * ksg1_465[k]
                   + f_3 * pc_x[k] * ksh_651[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pc_x, ksg0_466, ksg0_467, ksg0_468, ksg1_466, \
                         ksg1_467, ksg1_468, ksh_652, ksh_653, \
                         ksh_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = f_16 * ksg0_466[k]
                   - f_17 * ksg1_466[k]
                   + f_3 * pc_x[k] * ksh_652[k];

        t_870[k] = f_16 * ksg0_467[k]
                   - f_17 * ksg1_467[k]
                   + f_3 * pc_x[k] * ksh_653[k];

        t_871[k] = f_8 * ksg0_468[k]
                   - f_9 * ksg1_468[k]
                   + f_3 * pc_x[k] * ksh_654[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pc_x, ksg0_469, ksg0_470, ksg0_471, ksg1_469, \
                         ksg1_470, ksg1_471, ksh_655, ksh_656, \
                         ksh_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_8 * ksg0_469[k]
                   - f_9 * ksg1_469[k]
                   + f_3 * pc_x[k] * ksh_655[k];

        t_873[k] = f_8 * ksg0_470[k]
                   - f_9 * ksg1_470[k]
                   + f_3 * pc_x[k] * ksh_656[k];

        t_874[k] = f_6 * ksg0_471[k]
                   - f_7 * ksg1_471[k]
                   + f_3 * pc_x[k] * ksh_657[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pc_x, ksg0_472, ksg0_473, ksg0_474, ksg1_472, \
                         ksg1_473, ksg1_474, ksh_658, ksh_659, \
                         ksh_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = f_6 * ksg0_472[k]
                   - f_7 * ksg1_472[k]
                   + f_3 * pc_x[k] * ksh_658[k];

        t_876[k] = f_6 * ksg0_473[k]
                   - f_7 * ksg1_473[k]
                   + f_3 * pc_x[k] * ksh_659[k];

        t_877[k] = f_6 * ksg0_474[k]
                   - f_7 * ksg1_474[k]
                   + f_3 * pc_x[k] * ksh_660[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pc_x, ksg0_475, ksg0_476, ksg0_477, ksg1_475, \
                         ksg1_476, ksg1_477, ksh_661, ksh_662, \
                         ksh_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = f_4 * ksg0_475[k]
                   - f_5 * ksg1_475[k]
                   + f_3 * pc_x[k] * ksh_661[k];

        t_879[k] = f_4 * ksg0_476[k]
                   - f_5 * ksg1_476[k]
                   + f_3 * pc_x[k] * ksh_662[k];

        t_880[k] = f_4 * ksg0_477[k]
                   - f_5 * ksg1_477[k]
                   + f_3 * pc_x[k] * ksh_663[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, t_884, t_885, pc_x, ksg0_478, ksg0_479, \
                         ksg1_478, ksg1_479, ksh_664, ksh_665, ksh_666, ksh_667, \
                         ksh_668 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_4 * ksg0_478[k]
                   - f_5 * ksg1_478[k]
                   + f_3 * pc_x[k] * ksh_664[k];

        t_882[k] = f_4 * ksg0_479[k]
                   - f_5 * ksg1_479[k]
                   + f_3 * pc_x[k] * ksh_665[k];

        t_883[k] = f_3 * pc_x[k] * ksh_666[k];

        t_884[k] = f_3 * pc_x[k] * ksh_667[k];

        t_885[k] = f_3 * pc_x[k] * ksh_668[k];
    }

#pragma omp simd aligned(t_886, t_887, t_888, t_889, t_890, pc_x, pc_y, pc_z, ish_498, \
                         ish_519, ksg0_475, ksg1_475, ksh_666, ksh_669, ksh_670, \
                         ksh_671 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_886[k] = f_3 * pc_x[k] * ksh_669[k];

        t_887[k] = f_3 * pc_x[k] * ksh_670[k];

        t_888[k] = f_3 * pc_x[k] * ksh_671[k];

        t_889[k] = f_14 * ish_519[k]
                   + f_1 * ksg0_475[k]
                   - f_2 * ksg1_475[k]
                   + f_3 * pc_y[k] * ksh_666[k];

        t_890[k] = f_13 * ish_498[k]
                   + f_3 * pc_z[k] * ksh_666[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, pc_y, ish_521, ish_522, ish_523, ksg0_477, \
                         ksg0_478, ksg0_479, ksg1_477, ksg1_478, ksg1_479, ksh_668, ksh_669, \
                         ksh_670 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = f_14 * ish_521[k]
                   + f_8 * ksg0_477[k]
                   - f_9 * ksg1_477[k]
                   + f_3 * pc_y[k] * ksh_668[k];

        t_892[k] = f_14 * ish_522[k]
                   + f_6 * ksg0_478[k]
                   - f_7 * ksg1_478[k]
                   + f_3 * pc_y[k] * ksh_669[k];

        t_893[k] = f_14 * ish_523[k]
                   + f_4 * ksg0_479[k]
                   - f_5 * ksg1_479[k]
                   + f_3 * pc_y[k] * ksh_670[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, pc_x, pc_y, pc_z, ish_503, ish_524, ksg0_479, \
                         ksg0_480, ksg1_479, ksg1_480, ksh_671, \
                         ksh_672 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_14 * ish_524[k]
                   + f_3 * pc_y[k] * ksh_671[k];

        t_895[k] = f_13 * ish_503[k]
                   + f_1 * ksg0_479[k]
                   - f_2 * ksg1_479[k]
                   + f_3 * pc_z[k] * ksh_671[k];

        t_896[k] = f_1 * ksg0_480[k]
                   - f_2 * ksg1_480[k]
                   + f_3 * pc_x[k] * ksh_672[k];
    }

#pragma omp simd aligned(t_897, t_898, t_899, pc_x, ksg0_481, ksg0_482, ksg0_483, ksg1_481, \
                         ksg1_482, ksg1_483, ksh_673, ksh_674, \
                         ksh_675 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_897[k] = f_16 * ksg0_481[k]
                   - f_17 * ksg1_481[k]
                   + f_3 * pc_x[k] * ksh_673[k];

        t_898[k] = f_16 * ksg0_482[k]
                   - f_17 * ksg1_482[k]
                   + f_3 * pc_x[k] * ksh_674[k];

        t_899[k] = f_8 * ksg0_483[k]
                   - f_9 * ksg1_483[k]
                   + f_3 * pc_x[k] * ksh_675[k];
    }

#pragma omp simd aligned(t_900, t_901, t_902, pc_x, ksg0_484, ksg0_485, ksg0_486, ksg1_484, \
                         ksg1_485, ksg1_486, ksh_676, ksh_677, \
                         ksh_678 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_900[k] = f_8 * ksg0_484[k]
                   - f_9 * ksg1_484[k]
                   + f_3 * pc_x[k] * ksh_676[k];

        t_901[k] = f_8 * ksg0_485[k]
                   - f_9 * ksg1_485[k]
                   + f_3 * pc_x[k] * ksh_677[k];

        t_902[k] = f_6 * ksg0_486[k]
                   - f_7 * ksg1_486[k]
                   + f_3 * pc_x[k] * ksh_678[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, ksg0_487, ksg0_488, ksg0_489, ksg1_487, \
                         ksg1_488, ksg1_489, ksh_679, ksh_680, \
                         ksh_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_6 * ksg0_487[k]
                   - f_7 * ksg1_487[k]
                   + f_3 * pc_x[k] * ksh_679[k];

        t_904[k] = f_6 * ksg0_488[k]
                   - f_7 * ksg1_488[k]
                   + f_3 * pc_x[k] * ksh_680[k];

        t_905[k] = f_6 * ksg0_489[k]
                   - f_7 * ksg1_489[k]
                   + f_3 * pc_x[k] * ksh_681[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_x, ksg0_490, ksg0_491, ksg0_492, ksg1_490, \
                         ksg1_491, ksg1_492, ksh_682, ksh_683, \
                         ksh_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_4 * ksg0_490[k]
                   - f_5 * ksg1_490[k]
                   + f_3 * pc_x[k] * ksh_682[k];

        t_907[k] = f_4 * ksg0_491[k]
                   - f_5 * ksg1_491[k]
                   + f_3 * pc_x[k] * ksh_683[k];

        t_908[k] = f_4 * ksg0_492[k]
                   - f_5 * ksg1_492[k]
                   + f_3 * pc_x[k] * ksh_684[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, t_913, pc_x, ksg0_493, ksg0_494, \
                         ksg1_493, ksg1_494, ksh_685, ksh_686, ksh_687, ksh_688, \
                         ksh_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_4 * ksg0_493[k]
                   - f_5 * ksg1_493[k]
                   + f_3 * pc_x[k] * ksh_685[k];

        t_910[k] = f_4 * ksg0_494[k]
                   - f_5 * ksg1_494[k]
                   + f_3 * pc_x[k] * ksh_686[k];

        t_911[k] = f_3 * pc_x[k] * ksh_687[k];

        t_912[k] = f_3 * pc_x[k] * ksh_688[k];

        t_913[k] = f_3 * pc_x[k] * ksh_689[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, pc_x, pc_y, pc_z, ish_519, \
                         ish_540, ksg0_490, ksg1_490, ksh_687, ksh_690, ksh_691, \
                         ksh_692 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_3 * pc_x[k] * ksh_690[k];

        t_915[k] = f_3 * pc_x[k] * ksh_691[k];

        t_916[k] = f_3 * pc_x[k] * ksh_692[k];

        t_917[k] = f_13 * ish_540[k]
                   + f_1 * ksg0_490[k]
                   - f_2 * ksg1_490[k]
                   + f_3 * pc_y[k] * ksh_687[k];

        t_918[k] = f_14 * ish_519[k]
                   + f_3 * pc_z[k] * ksh_687[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pc_y, ish_542, ish_543, ish_544, ksg0_492, \
                         ksg0_493, ksg0_494, ksg1_492, ksg1_493, ksg1_494, ksh_689, ksh_690, \
                         ksh_691 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = f_13 * ish_542[k]
                   + f_8 * ksg0_492[k]
                   - f_9 * ksg1_492[k]
                   + f_3 * pc_y[k] * ksh_689[k];

        t_920[k] = f_13 * ish_543[k]
                   + f_6 * ksg0_493[k]
                   - f_7 * ksg1_493[k]
                   + f_3 * pc_y[k] * ksh_690[k];

        t_921[k] = f_13 * ish_544[k]
                   + f_4 * ksg0_494[k]
                   - f_5 * ksg1_494[k]
                   + f_3 * pc_y[k] * ksh_691[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, pc_x, pc_y, pc_z, ish_524, ish_545, ksg0_494, \
                         ksg0_495, ksg1_494, ksg1_495, ksh_692, \
                         ksh_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_13 * ish_545[k]
                   + f_3 * pc_y[k] * ksh_692[k];

        t_923[k] = f_14 * ish_524[k]
                   + f_1 * ksg0_494[k]
                   - f_2 * ksg1_494[k]
                   + f_3 * pc_z[k] * ksh_692[k];

        t_924[k] = f_1 * ksg0_495[k]
                   - f_2 * ksg1_495[k]
                   + f_3 * pc_x[k] * ksh_693[k];
    }

#pragma omp simd aligned(t_925, t_926, t_927, pc_x, ksg0_496, ksg0_497, ksg0_498, ksg1_496, \
                         ksg1_497, ksg1_498, ksh_694, ksh_695, \
                         ksh_696 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_925[k] = f_16 * ksg0_496[k]
                   - f_17 * ksg1_496[k]
                   + f_3 * pc_x[k] * ksh_694[k];

        t_926[k] = f_16 * ksg0_497[k]
                   - f_17 * ksg1_497[k]
                   + f_3 * pc_x[k] * ksh_695[k];

        t_927[k] = f_8 * ksg0_498[k]
                   - f_9 * ksg1_498[k]
                   + f_3 * pc_x[k] * ksh_696[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, pc_x, ksg0_499, ksg0_500, ksg0_501, ksg1_499, \
                         ksg1_500, ksg1_501, ksh_697, ksh_698, \
                         ksh_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_8 * ksg0_499[k]
                   - f_9 * ksg1_499[k]
                   + f_3 * pc_x[k] * ksh_697[k];

        t_929[k] = f_8 * ksg0_500[k]
                   - f_9 * ksg1_500[k]
                   + f_3 * pc_x[k] * ksh_698[k];

        t_930[k] = f_6 * ksg0_501[k]
                   - f_7 * ksg1_501[k]
                   + f_3 * pc_x[k] * ksh_699[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, pc_x, ksg0_502, ksg0_503, ksg0_504, ksg1_502, \
                         ksg1_503, ksg1_504, ksh_700, ksh_701, \
                         ksh_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_6 * ksg0_502[k]
                   - f_7 * ksg1_502[k]
                   + f_3 * pc_x[k] * ksh_700[k];

        t_932[k] = f_6 * ksg0_503[k]
                   - f_7 * ksg1_503[k]
                   + f_3 * pc_x[k] * ksh_701[k];

        t_933[k] = f_6 * ksg0_504[k]
                   - f_7 * ksg1_504[k]
                   + f_3 * pc_x[k] * ksh_702[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, pc_x, ksg0_505, ksg0_506, ksg0_507, ksg1_505, \
                         ksg1_506, ksg1_507, ksh_703, ksh_704, \
                         ksh_705 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_4 * ksg0_505[k]
                   - f_5 * ksg1_505[k]
                   + f_3 * pc_x[k] * ksh_703[k];

        t_935[k] = f_4 * ksg0_506[k]
                   - f_5 * ksg1_506[k]
                   + f_3 * pc_x[k] * ksh_704[k];

        t_936[k] = f_4 * ksg0_507[k]
                   - f_5 * ksg1_507[k]
                   + f_3 * pc_x[k] * ksh_705[k];
    }

#pragma omp simd aligned(t_937, t_938, t_939, t_940, t_941, pc_x, ksg0_508, ksg0_509, \
                         ksg1_508, ksg1_509, ksh_706, ksh_707, ksh_708, ksh_709, \
                         ksh_710 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_937[k] = f_4 * ksg0_508[k]
                   - f_5 * ksg1_508[k]
                   + f_3 * pc_x[k] * ksh_706[k];

        t_938[k] = f_4 * ksg0_509[k]
                   - f_5 * ksg1_509[k]
                   + f_3 * pc_x[k] * ksh_707[k];

        t_939[k] = f_3 * pc_x[k] * ksh_708[k];

        t_940[k] = f_3 * pc_x[k] * ksh_709[k];

        t_941[k] = f_3 * pc_x[k] * ksh_710[k];
    }

#pragma omp simd aligned(t_942, t_943, t_944, t_945, t_946, pc_x, pc_y, pc_z, ish_540, \
                         ish_561, ksg0_505, ksg1_505, ksh_708, ksh_711, ksh_712, \
                         ksh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_942[k] = f_3 * pc_x[k] * ksh_711[k];

        t_943[k] = f_3 * pc_x[k] * ksh_712[k];

        t_944[k] = f_3 * pc_x[k] * ksh_713[k];

        t_945[k] = f_12 * ish_561[k]
                   + f_1 * ksg0_505[k]
                   - f_2 * ksg1_505[k]
                   + f_3 * pc_y[k] * ksh_708[k];

        t_946[k] = f_18 * ish_540[k]
                   + f_3 * pc_z[k] * ksh_708[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pc_y, ish_563, ish_564, ish_565, ksg0_507, \
                         ksg0_508, ksg0_509, ksg1_507, ksg1_508, ksg1_509, ksh_710, ksh_711, \
                         ksh_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_12 * ish_563[k]
                   + f_8 * ksg0_507[k]
                   - f_9 * ksg1_507[k]
                   + f_3 * pc_y[k] * ksh_710[k];

        t_948[k] = f_12 * ish_564[k]
                   + f_6 * ksg0_508[k]
                   - f_7 * ksg1_508[k]
                   + f_3 * pc_y[k] * ksh_711[k];

        t_949[k] = f_12 * ish_565[k]
                   + f_4 * ksg0_509[k]
                   - f_5 * ksg1_509[k]
                   + f_3 * pc_y[k] * ksh_712[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pa_y, pc_y, pc_z, isi0_756, ish_545, ish_566, \
                         isi1_756, ksg0_509, ksg1_509, ksh_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = f_12 * ish_566[k]
                   + f_3 * pc_y[k] * ksh_713[k];

        t_951[k] = f_18 * ish_545[k]
                   + f_1 * ksg0_509[k]
                   - f_2 * ksg1_509[k]
                   + f_3 * pc_z[k] * ksh_713[k];

        t_952[k] = pa_y[k] * isi0_756[k]
                   - f_10 * pc_y[k] * isi1_756[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, pa_y, pc_x, pc_y, isi0_758, isi1_758, ksg0_511, \
                         ksg0_513, ksg1_511, ksg1_513, ksh_715, \
                         ksh_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_16 * ksg0_511[k]
                   - f_17 * ksg1_511[k]
                   + f_3 * pc_x[k] * ksh_715[k];

        t_954[k] = pa_y[k] * isi0_758[k]
                   - f_10 * pc_y[k] * isi1_758[k];

        t_955[k] = f_8 * ksg0_513[k]
                   - f_9 * ksg1_513[k]
                   + f_3 * pc_x[k] * ksh_717[k];
    }
}

static auto
compute_prim_ksi_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t isi0,
                                                          const size_t ish, const size_t isi1,
                                                          const size_t ksg0, const size_t ksg1,
                                                          const size_t ksh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
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
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.0 / gamma;
    const auto f_17 = 2.0 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *isi0_761 = buffer.data(isi0 + 761);
    const auto *isi0_765 = buffer.data(isi0 + 765);
    const auto *isi0_770 = buffer.data(isi0 + 770);
    const auto *isi0_777 = buffer.data(isi0 + 777);
    const auto *isi0_779 = buffer.data(isi0 + 779);
    const auto *isi0_780 = buffer.data(isi0 + 780);
    const auto *isi0_781 = buffer.data(isi0 + 781);
    const auto *isi0_783 = buffer.data(isi0 + 783);

    const auto *ish_561 = buffer.data(ish + 561);
    const auto *ish_582 = buffer.data(ish + 582);
    const auto *ish_584 = buffer.data(ish + 584);
    const auto *ish_585 = buffer.data(ish + 585);
    const auto *ish_586 = buffer.data(ish + 586);
    const auto *ish_587 = buffer.data(ish + 587);

    const auto *isi1_761 = buffer.data(isi1 + 761);
    const auto *isi1_765 = buffer.data(isi1 + 765);
    const auto *isi1_770 = buffer.data(isi1 + 770);
    const auto *isi1_777 = buffer.data(isi1 + 777);
    const auto *isi1_779 = buffer.data(isi1 + 779);
    const auto *isi1_780 = buffer.data(isi1 + 780);
    const auto *isi1_781 = buffer.data(isi1 + 781);
    const auto *isi1_783 = buffer.data(isi1 + 783);

    const auto *ksg0_514 = buffer.data(ksg0 + 514);
    const auto *ksg0_516 = buffer.data(ksg0 + 516);
    const auto *ksg0_517 = buffer.data(ksg0 + 517);
    const auto *ksg0_518 = buffer.data(ksg0 + 518);
    const auto *ksg0_520 = buffer.data(ksg0 + 520);
    const auto *ksg0_521 = buffer.data(ksg0 + 521);
    const auto *ksg0_522 = buffer.data(ksg0 + 522);
    const auto *ksg0_523 = buffer.data(ksg0 + 523);
    const auto *ksg0_525 = buffer.data(ksg0 + 525);
    const auto *ksg0_527 = buffer.data(ksg0 + 527);
    const auto *ksg0_528 = buffer.data(ksg0 + 528);
    const auto *ksg0_530 = buffer.data(ksg0 + 530);
    const auto *ksg0_531 = buffer.data(ksg0 + 531);
    const auto *ksg0_532 = buffer.data(ksg0 + 532);
    const auto *ksg0_534 = buffer.data(ksg0 + 534);
    const auto *ksg0_535 = buffer.data(ksg0 + 535);
    const auto *ksg0_536 = buffer.data(ksg0 + 536);
    const auto *ksg0_537 = buffer.data(ksg0 + 537);
    const auto *ksg0_538 = buffer.data(ksg0 + 538);
    const auto *ksg0_539 = buffer.data(ksg0 + 539);

    const auto *ksg1_514 = buffer.data(ksg1 + 514);
    const auto *ksg1_516 = buffer.data(ksg1 + 516);
    const auto *ksg1_517 = buffer.data(ksg1 + 517);
    const auto *ksg1_518 = buffer.data(ksg1 + 518);
    const auto *ksg1_520 = buffer.data(ksg1 + 520);
    const auto *ksg1_521 = buffer.data(ksg1 + 521);
    const auto *ksg1_522 = buffer.data(ksg1 + 522);
    const auto *ksg1_523 = buffer.data(ksg1 + 523);
    const auto *ksg1_525 = buffer.data(ksg1 + 525);
    const auto *ksg1_527 = buffer.data(ksg1 + 527);
    const auto *ksg1_528 = buffer.data(ksg1 + 528);
    const auto *ksg1_530 = buffer.data(ksg1 + 530);
    const auto *ksg1_531 = buffer.data(ksg1 + 531);
    const auto *ksg1_532 = buffer.data(ksg1 + 532);
    const auto *ksg1_534 = buffer.data(ksg1 + 534);
    const auto *ksg1_535 = buffer.data(ksg1 + 535);
    const auto *ksg1_536 = buffer.data(ksg1 + 536);
    const auto *ksg1_537 = buffer.data(ksg1 + 537);
    const auto *ksg1_538 = buffer.data(ksg1 + 538);
    const auto *ksg1_539 = buffer.data(ksg1 + 539);

    const auto *ksh_718 = buffer.data(ksh + 718);
    const auto *ksh_720 = buffer.data(ksh + 720);
    const auto *ksh_721 = buffer.data(ksh + 721);
    const auto *ksh_722 = buffer.data(ksh + 722);
    const auto *ksh_724 = buffer.data(ksh + 724);
    const auto *ksh_725 = buffer.data(ksh + 725);
    const auto *ksh_726 = buffer.data(ksh + 726);
    const auto *ksh_727 = buffer.data(ksh + 727);
    const auto *ksh_729 = buffer.data(ksh + 729);
    const auto *ksh_730 = buffer.data(ksh + 730);
    const auto *ksh_731 = buffer.data(ksh + 731);
    const auto *ksh_732 = buffer.data(ksh + 732);
    const auto *ksh_733 = buffer.data(ksh + 733);
    const auto *ksh_734 = buffer.data(ksh + 734);
    const auto *ksh_735 = buffer.data(ksh + 735);
    const auto *ksh_737 = buffer.data(ksh + 737);
    const auto *ksh_738 = buffer.data(ksh + 738);
    const auto *ksh_740 = buffer.data(ksh + 740);
    const auto *ksh_741 = buffer.data(ksh + 741);
    const auto *ksh_742 = buffer.data(ksh + 742);
    const auto *ksh_744 = buffer.data(ksh + 744);
    const auto *ksh_745 = buffer.data(ksh + 745);
    const auto *ksh_746 = buffer.data(ksh + 746);
    const auto *ksh_747 = buffer.data(ksh + 747);
    const auto *ksh_749 = buffer.data(ksh + 749);
    const auto *ksh_750 = buffer.data(ksh + 750);
    const auto *ksh_751 = buffer.data(ksh + 751);
    const auto *ksh_752 = buffer.data(ksh + 752);
    const auto *ksh_753 = buffer.data(ksh + 753);
    const auto *ksh_754 = buffer.data(ksh + 754);
    const auto *ksh_755 = buffer.data(ksh + 755);

#pragma omp simd aligned(t_956, t_957, t_958, pa_y, pc_x, pc_y, isi0_761, isi1_761, ksg0_514, \
                         ksg0_516, ksg1_514, ksg1_516, ksh_718, \
                         ksh_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = f_8 * ksg0_514[k]
                   - f_9 * ksg1_514[k]
                   + f_3 * pc_x[k] * ksh_718[k];

        t_957[k] = pa_y[k] * isi0_761[k]
                   - f_10 * pc_y[k] * isi1_761[k];

        t_958[k] = f_6 * ksg0_516[k]
                   - f_7 * ksg1_516[k]
                   + f_3 * pc_x[k] * ksh_720[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, pa_y, pc_x, pc_y, isi0_765, isi1_765, ksg0_517, \
                         ksg0_518, ksg1_517, ksg1_518, ksh_721, \
                         ksh_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = f_6 * ksg0_517[k]
                   - f_7 * ksg1_517[k]
                   + f_3 * pc_x[k] * ksh_721[k];

        t_960[k] = f_6 * ksg0_518[k]
                   - f_7 * ksg1_518[k]
                   + f_3 * pc_x[k] * ksh_722[k];

        t_961[k] = pa_y[k] * isi0_765[k]
                   - f_10 * pc_y[k] * isi1_765[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, pc_x, ksg0_520, ksg0_521, ksg0_522, ksg1_520, \
                         ksg1_521, ksg1_522, ksh_724, ksh_725, \
                         ksh_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_4 * ksg0_520[k]
                   - f_5 * ksg1_520[k]
                   + f_3 * pc_x[k] * ksh_724[k];

        t_963[k] = f_4 * ksg0_521[k]
                   - f_5 * ksg1_521[k]
                   + f_3 * pc_x[k] * ksh_725[k];

        t_964[k] = f_4 * ksg0_522[k]
                   - f_5 * ksg1_522[k]
                   + f_3 * pc_x[k] * ksh_726[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, t_969, pa_y, pc_x, pc_y, isi0_770, \
                         isi1_770, ksg0_523, ksg1_523, ksh_727, ksh_729, ksh_730, \
                         ksh_731 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = f_4 * ksg0_523[k]
                   - f_5 * ksg1_523[k]
                   + f_3 * pc_x[k] * ksh_727[k];

        t_966[k] = pa_y[k] * isi0_770[k]
                   - f_10 * pc_y[k] * isi1_770[k];

        t_967[k] = f_3 * pc_x[k] * ksh_729[k];

        t_968[k] = f_3 * pc_x[k] * ksh_730[k];

        t_969[k] = f_3 * pc_x[k] * ksh_731[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, pa_y, pc_x, pc_y, isi0_777, ish_582, \
                         isi1_777, ksh_732, ksh_733, ksh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_3 * pc_x[k] * ksh_732[k];

        t_971[k] = f_3 * pc_x[k] * ksh_733[k];

        t_972[k] = f_3 * pc_x[k] * ksh_734[k];

        t_973[k] = pa_y[k] * isi0_777[k]
                   + f_15 * ish_582[k]
                   - f_10 * pc_y[k] * isi1_777[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, pa_y, pc_y, pc_z, isi0_779, isi0_780, ish_561, \
                         ish_584, ish_585, isi1_779, isi1_780, \
                         ksh_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_15 * ish_561[k]
                   + f_3 * pc_z[k] * ksh_729[k];

        t_975[k] = pa_y[k] * isi0_779[k]
                   + f_14 * ish_584[k]
                   - f_10 * pc_y[k] * isi1_779[k];

        t_976[k] = pa_y[k] * isi0_780[k]
                   + f_13 * ish_585[k]
                   - f_10 * pc_y[k] * isi1_780[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, pa_y, pc_y, isi0_781, isi0_783, ish_586, \
                         ish_587, isi1_781, isi1_783, ksh_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = pa_y[k] * isi0_781[k]
                   + f_12 * ish_586[k]
                   - f_10 * pc_y[k] * isi1_781[k];

        t_978[k] = f_11 * ish_587[k]
                   + f_3 * pc_y[k] * ksh_734[k];

        t_979[k] = pa_y[k] * isi0_783[k]
                   - f_10 * pc_y[k] * isi1_783[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, t_984, pc_x, pc_y, ksg0_525, ksg0_527, \
                         ksg0_528, ksg1_525, ksg1_527, ksg1_528, ksh_735, ksh_737, \
                         ksh_738 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_1 * ksg0_525[k]
                   - f_2 * ksg1_525[k]
                   + f_3 * pc_x[k] * ksh_735[k];

        t_981[k] = f_3 * pc_y[k] * ksh_735[k];

        t_982[k] = f_16 * ksg0_527[k]
                   - f_17 * ksg1_527[k]
                   + f_3 * pc_x[k] * ksh_737[k];

        t_983[k] = f_8 * ksg0_528[k]
                   - f_9 * ksg1_528[k]
                   + f_3 * pc_x[k] * ksh_738[k];

        t_984[k] = f_3 * pc_y[k] * ksh_737[k];
    }

#pragma omp simd aligned(t_985, t_986, t_987, t_988, pc_x, pc_y, ksg0_530, ksg0_531, ksg0_532, \
                         ksg1_530, ksg1_531, ksg1_532, ksh_740, ksh_741, \
                         ksh_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_985[k] = f_8 * ksg0_530[k]
                   - f_9 * ksg1_530[k]
                   + f_3 * pc_x[k] * ksh_740[k];

        t_986[k] = f_6 * ksg0_531[k]
                   - f_7 * ksg1_531[k]
                   + f_3 * pc_x[k] * ksh_741[k];

        t_987[k] = f_6 * ksg0_532[k]
                   - f_7 * ksg1_532[k]
                   + f_3 * pc_x[k] * ksh_742[k];

        t_988[k] = f_3 * pc_y[k] * ksh_740[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, pc_x, ksg0_534, ksg0_535, ksg0_536, ksg1_534, \
                         ksg1_535, ksg1_536, ksh_744, ksh_745, \
                         ksh_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = f_6 * ksg0_534[k]
                   - f_7 * ksg1_534[k]
                   + f_3 * pc_x[k] * ksh_744[k];

        t_990[k] = f_4 * ksg0_535[k]
                   - f_5 * ksg1_535[k]
                   + f_3 * pc_x[k] * ksh_745[k];

        t_991[k] = f_4 * ksg0_536[k]
                   - f_5 * ksg1_536[k]
                   + f_3 * pc_x[k] * ksh_746[k];
    }

#pragma omp simd aligned(t_992, t_993, t_994, t_995, t_996, pc_x, pc_y, ksg0_537, ksg0_539, \
                         ksg1_537, ksg1_539, ksh_744, ksh_747, ksh_749, ksh_750, \
                         ksh_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_992[k] = f_4 * ksg0_537[k]
                   - f_5 * ksg1_537[k]
                   + f_3 * pc_x[k] * ksh_747[k];

        t_993[k] = f_3 * pc_y[k] * ksh_744[k];

        t_994[k] = f_4 * ksg0_539[k]
                   - f_5 * ksg1_539[k]
                   + f_3 * pc_x[k] * ksh_749[k];

        t_995[k] = f_3 * pc_x[k] * ksh_750[k];

        t_996[k] = f_3 * pc_x[k] * ksh_751[k];
    }

#pragma omp simd aligned(t_997, t_998, t_999, t_1000, t_1001, pc_x, pc_y, ksg0_535, ksg1_535, \
                         ksh_750, ksh_752, ksh_753, ksh_754, ksh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_997[k] = f_3 * pc_x[k] * ksh_752[k];

        t_998[k] = f_3 * pc_x[k] * ksh_753[k];

        t_999[k] = f_3 * pc_x[k] * ksh_754[k];

        t_1000[k] = f_3 * pc_x[k] * ksh_755[k];

        t_1001[k] = f_1 * ksg0_535[k]
                    - f_2 * ksg1_535[k]
                    + f_3 * pc_y[k] * ksh_750[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_y, ksg0_536, ksg0_537, ksg0_538, ksg1_536, \
                         ksg1_537, ksg1_538, ksh_751, ksh_752, \
                         ksh_753 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_16 * ksg0_536[k]
                    - f_17 * ksg1_536[k]
                    + f_3 * pc_y[k] * ksh_751[k];

        t_1003[k] = f_8 * ksg0_537[k]
                    - f_9 * ksg1_537[k]
                    + f_3 * pc_y[k] * ksh_752[k];

        t_1004[k] = f_6 * ksg0_538[k]
                    - f_7 * ksg1_538[k]
                    + f_3 * pc_y[k] * ksh_753[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pc_y, pc_z, ish_587, ksg0_539, ksg1_539, \
                         ksh_754, ksh_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = f_4 * ksg0_539[k]
                    - f_5 * ksg1_539[k]
                    + f_3 * pc_y[k] * ksh_754[k];

        t_1006[k] = f_3 * pc_y[k] * ksh_755[k];

        t_1007[k] = f_0 * ish_587[k]
                    + f_1 * ksg0_539[k]
                    - f_2 * ksg1_539[k]
                    + f_3 * pc_z[k] * ksh_755[k];
    }
}

auto
compute_prim_ksi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t isi0, const size_t ish,
                                                   const size_t isi1, const size_t ksg0,
                                                   const size_t ksg1, const size_t ksh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_ksi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);

    compute_prim_ksi_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, isi0, ish,
                                                              isi1, ksg0, ksg1, ksh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
