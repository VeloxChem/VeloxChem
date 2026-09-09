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


#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_hsi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t hsg0, const size_t hsg1,
                                                          const size_t hsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);

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

    const auto *gsi0_0 = buffer.data(gsi0 + 0);
    const auto *gsi0_3 = buffer.data(gsi0 + 3);
    const auto *gsi0_5 = buffer.data(gsi0 + 5);
    const auto *gsi0_6 = buffer.data(gsi0 + 6);
    const auto *gsi0_9 = buffer.data(gsi0 + 9);
    const auto *gsi0_10 = buffer.data(gsi0 + 10);
    const auto *gsi0_14 = buffer.data(gsi0 + 14);
    const auto *gsi0_21 = buffer.data(gsi0 + 21);
    const auto *gsi0_27 = buffer.data(gsi0 + 27);
    const auto *gsi0_31 = buffer.data(gsi0 + 31);
    const auto *gsi0_34 = buffer.data(gsi0 + 34);
    const auto *gsi0_38 = buffer.data(gsi0 + 38);
    const auto *gsi0_56 = buffer.data(gsi0 + 56);
    const auto *gsi0_61 = buffer.data(gsi0 + 61);
    const auto *gsi0_65 = buffer.data(gsi0 + 65);
    const auto *gsi0_68 = buffer.data(gsi0 + 68);
    const auto *gsi0_70 = buffer.data(gsi0 + 70);

    const auto *gsh_0 = buffer.data(gsh + 0);
    const auto *gsh_1 = buffer.data(gsh + 1);
    const auto *gsh_2 = buffer.data(gsh + 2);
    const auto *gsh_3 = buffer.data(gsh + 3);
    const auto *gsh_5 = buffer.data(gsh + 5);
    const auto *gsh_6 = buffer.data(gsh + 6);
    const auto *gsh_9 = buffer.data(gsh + 9);
    const auto *gsh_15 = buffer.data(gsh + 15);
    const auto *gsh_17 = buffer.data(gsh + 17);
    const auto *gsh_18 = buffer.data(gsh + 18);
    const auto *gsh_20 = buffer.data(gsh + 20);
    const auto *gsh_21 = buffer.data(gsh + 21);
    const auto *gsh_24 = buffer.data(gsh + 24);
    const auto *gsh_26 = buffer.data(gsh + 26);
    const auto *gsh_27 = buffer.data(gsh + 27);
    const auto *gsh_30 = buffer.data(gsh + 30);
    const auto *gsh_36 = buffer.data(gsh + 36);
    const auto *gsh_38 = buffer.data(gsh + 38);
    const auto *gsh_39 = buffer.data(gsh + 39);
    const auto *gsh_40 = buffer.data(gsh + 40);
    const auto *gsh_41 = buffer.data(gsh + 41);
    const auto *gsh_42 = buffer.data(gsh + 42);
    const auto *gsh_44 = buffer.data(gsh + 44);
    const auto *gsh_47 = buffer.data(gsh + 47);
    const auto *gsh_50 = buffer.data(gsh + 50);
    const auto *gsh_51 = buffer.data(gsh + 51);
    const auto *gsh_57 = buffer.data(gsh + 57);
    const auto *gsh_58 = buffer.data(gsh + 58);
    const auto *gsh_59 = buffer.data(gsh + 59);
    const auto *gsh_60 = buffer.data(gsh + 60);
    const auto *gsh_62 = buffer.data(gsh + 62);
    const auto *gsh_63 = buffer.data(gsh + 63);
    const auto *gsh_66 = buffer.data(gsh + 66);
    const auto *gsh_69 = buffer.data(gsh + 69);
    const auto *gsh_73 = buffer.data(gsh + 73);
    const auto *gsh_78 = buffer.data(gsh + 78);
    const auto *gsh_80 = buffer.data(gsh + 80);
    const auto *gsh_81 = buffer.data(gsh + 81);
    const auto *gsh_82 = buffer.data(gsh + 82);
    const auto *gsh_83 = buffer.data(gsh + 83);
    const auto *gsh_99 = buffer.data(gsh + 99);

    const auto *gsi1_0 = buffer.data(gsi1 + 0);
    const auto *gsi1_3 = buffer.data(gsi1 + 3);
    const auto *gsi1_5 = buffer.data(gsi1 + 5);
    const auto *gsi1_6 = buffer.data(gsi1 + 6);
    const auto *gsi1_9 = buffer.data(gsi1 + 9);
    const auto *gsi1_10 = buffer.data(gsi1 + 10);
    const auto *gsi1_14 = buffer.data(gsi1 + 14);
    const auto *gsi1_21 = buffer.data(gsi1 + 21);
    const auto *gsi1_27 = buffer.data(gsi1 + 27);
    const auto *gsi1_31 = buffer.data(gsi1 + 31);
    const auto *gsi1_34 = buffer.data(gsi1 + 34);
    const auto *gsi1_38 = buffer.data(gsi1 + 38);
    const auto *gsi1_56 = buffer.data(gsi1 + 56);
    const auto *gsi1_61 = buffer.data(gsi1 + 61);
    const auto *gsi1_65 = buffer.data(gsi1 + 65);
    const auto *gsi1_68 = buffer.data(gsi1 + 68);
    const auto *gsi1_70 = buffer.data(gsi1 + 70);

    const auto *hsg0_0 = buffer.data(hsg0 + 0);
    const auto *hsg0_1 = buffer.data(hsg0 + 1);
    const auto *hsg0_2 = buffer.data(hsg0 + 2);
    const auto *hsg0_3 = buffer.data(hsg0 + 3);
    const auto *hsg0_5 = buffer.data(hsg0 + 5);
    const auto *hsg0_10 = buffer.data(hsg0 + 10);
    const auto *hsg0_12 = buffer.data(hsg0 + 12);
    const auto *hsg0_13 = buffer.data(hsg0 + 13);
    const auto *hsg0_14 = buffer.data(hsg0 + 14);
    const auto *hsg0_18 = buffer.data(hsg0 + 18);
    const auto *hsg0_25 = buffer.data(hsg0 + 25);
    const auto *hsg0_26 = buffer.data(hsg0 + 26);
    const auto *hsg0_27 = buffer.data(hsg0 + 27);
    const auto *hsg0_32 = buffer.data(hsg0 + 32);
    const auto *hsg0_34 = buffer.data(hsg0 + 34);
    const auto *hsg0_35 = buffer.data(hsg0 + 35);
    const auto *hsg0_41 = buffer.data(hsg0 + 41);
    const auto *hsg0_42 = buffer.data(hsg0 + 42);
    const auto *hsg0_43 = buffer.data(hsg0 + 43);
    const auto *hsg0_44 = buffer.data(hsg0 + 44);
    const auto *hsg0_45 = buffer.data(hsg0 + 45);
    const auto *hsg0_47 = buffer.data(hsg0 + 47);
    const auto *hsg0_48 = buffer.data(hsg0 + 48);
    const auto *hsg0_50 = buffer.data(hsg0 + 50);
    const auto *hsg0_51 = buffer.data(hsg0 + 51);
    const auto *hsg0_55 = buffer.data(hsg0 + 55);
    const auto *hsg0_56 = buffer.data(hsg0 + 56);
    const auto *hsg0_57 = buffer.data(hsg0 + 57);
    const auto *hsg0_59 = buffer.data(hsg0 + 59);

    const auto *hsg1_0 = buffer.data(hsg1 + 0);
    const auto *hsg1_1 = buffer.data(hsg1 + 1);
    const auto *hsg1_2 = buffer.data(hsg1 + 2);
    const auto *hsg1_3 = buffer.data(hsg1 + 3);
    const auto *hsg1_5 = buffer.data(hsg1 + 5);
    const auto *hsg1_10 = buffer.data(hsg1 + 10);
    const auto *hsg1_12 = buffer.data(hsg1 + 12);
    const auto *hsg1_13 = buffer.data(hsg1 + 13);
    const auto *hsg1_14 = buffer.data(hsg1 + 14);
    const auto *hsg1_18 = buffer.data(hsg1 + 18);
    const auto *hsg1_25 = buffer.data(hsg1 + 25);
    const auto *hsg1_26 = buffer.data(hsg1 + 26);
    const auto *hsg1_27 = buffer.data(hsg1 + 27);
    const auto *hsg1_32 = buffer.data(hsg1 + 32);
    const auto *hsg1_34 = buffer.data(hsg1 + 34);
    const auto *hsg1_35 = buffer.data(hsg1 + 35);
    const auto *hsg1_41 = buffer.data(hsg1 + 41);
    const auto *hsg1_42 = buffer.data(hsg1 + 42);
    const auto *hsg1_43 = buffer.data(hsg1 + 43);
    const auto *hsg1_44 = buffer.data(hsg1 + 44);
    const auto *hsg1_45 = buffer.data(hsg1 + 45);
    const auto *hsg1_47 = buffer.data(hsg1 + 47);
    const auto *hsg1_48 = buffer.data(hsg1 + 48);
    const auto *hsg1_50 = buffer.data(hsg1 + 50);
    const auto *hsg1_51 = buffer.data(hsg1 + 51);
    const auto *hsg1_55 = buffer.data(hsg1 + 55);
    const auto *hsg1_56 = buffer.data(hsg1 + 56);
    const auto *hsg1_57 = buffer.data(hsg1 + 57);
    const auto *hsg1_59 = buffer.data(hsg1 + 59);

    const auto *hsh_0 = buffer.data(hsh + 0);
    const auto *hsh_1 = buffer.data(hsh + 1);
    const auto *hsh_2 = buffer.data(hsh + 2);
    const auto *hsh_3 = buffer.data(hsh + 3);
    const auto *hsh_5 = buffer.data(hsh + 5);
    const auto *hsh_6 = buffer.data(hsh + 6);
    const auto *hsh_8 = buffer.data(hsh + 8);
    const auto *hsh_9 = buffer.data(hsh + 9);
    const auto *hsh_10 = buffer.data(hsh + 10);
    const auto *hsh_14 = buffer.data(hsh + 14);
    const auto *hsh_15 = buffer.data(hsh + 15);
    const auto *hsh_17 = buffer.data(hsh + 17);
    const auto *hsh_18 = buffer.data(hsh + 18);
    const auto *hsh_19 = buffer.data(hsh + 19);
    const auto *hsh_20 = buffer.data(hsh + 20);
    const auto *hsh_21 = buffer.data(hsh + 21);
    const auto *hsh_22 = buffer.data(hsh + 22);
    const auto *hsh_24 = buffer.data(hsh + 24);
    const auto *hsh_26 = buffer.data(hsh + 26);
    const auto *hsh_27 = buffer.data(hsh + 27);
    const auto *hsh_28 = buffer.data(hsh + 28);
    const auto *hsh_30 = buffer.data(hsh + 30);
    const auto *hsh_31 = buffer.data(hsh + 31);
    const auto *hsh_36 = buffer.data(hsh + 36);
    const auto *hsh_37 = buffer.data(hsh + 37);
    const auto *hsh_38 = buffer.data(hsh + 38);
    const auto *hsh_39 = buffer.data(hsh + 39);
    const auto *hsh_40 = buffer.data(hsh + 40);
    const auto *hsh_41 = buffer.data(hsh + 41);
    const auto *hsh_42 = buffer.data(hsh + 42);
    const auto *hsh_44 = buffer.data(hsh + 44);
    const auto *hsh_46 = buffer.data(hsh + 46);
    const auto *hsh_47 = buffer.data(hsh + 47);
    const auto *hsh_49 = buffer.data(hsh + 49);
    const auto *hsh_50 = buffer.data(hsh + 50);
    const auto *hsh_51 = buffer.data(hsh + 51);
    const auto *hsh_56 = buffer.data(hsh + 56);
    const auto *hsh_57 = buffer.data(hsh + 57);
    const auto *hsh_58 = buffer.data(hsh + 58);
    const auto *hsh_59 = buffer.data(hsh + 59);
    const auto *hsh_60 = buffer.data(hsh + 60);
    const auto *hsh_61 = buffer.data(hsh + 61);
    const auto *hsh_62 = buffer.data(hsh + 62);
    const auto *hsh_63 = buffer.data(hsh + 63);
    const auto *hsh_64 = buffer.data(hsh + 64);
    const auto *hsh_65 = buffer.data(hsh + 65);
    const auto *hsh_66 = buffer.data(hsh + 66);
    const auto *hsh_68 = buffer.data(hsh + 68);
    const auto *hsh_69 = buffer.data(hsh + 69);
    const auto *hsh_70 = buffer.data(hsh + 70);
    const auto *hsh_72 = buffer.data(hsh + 72);
    const auto *hsh_73 = buffer.data(hsh + 73);
    const auto *hsh_78 = buffer.data(hsh + 78);
    const auto *hsh_79 = buffer.data(hsh + 79);
    const auto *hsh_80 = buffer.data(hsh + 80);
    const auto *hsh_81 = buffer.data(hsh + 81);
    const auto *hsh_82 = buffer.data(hsh + 82);
    const auto *hsh_83 = buffer.data(hsh + 83);
    const auto *hsh_84 = buffer.data(hsh + 84);
    const auto *hsh_86 = buffer.data(hsh + 86);
    const auto *hsh_87 = buffer.data(hsh + 87);
    const auto *hsh_89 = buffer.data(hsh + 89);
    const auto *hsh_90 = buffer.data(hsh + 90);
    const auto *hsh_93 = buffer.data(hsh + 93);
    const auto *hsh_99 = buffer.data(hsh + 99);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, gsh_0, hsg0_0, \
                         hsg1_0, hsh_0, hsh_1, hsh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gsh_0[k]
                 + f_1 * hsg0_0[k]
                 - f_2 * hsg1_0[k]
                 + f_3 * pc_x[k] * hsh_0[k];

        t_1[k] = f_3 * pc_y[k] * hsh_0[k];

        t_2[k] = f_3 * pc_z[k] * hsh_0[k];

        t_3[k] = f_4 * hsg0_0[k]
                 - f_5 * hsg1_0[k]
                 + f_3 * pc_y[k] * hsh_1[k];

        t_4[k] = f_3 * pc_y[k] * hsh_2[k];

        t_5[k] = f_4 * hsg0_0[k]
                 - f_5 * hsg1_0[k]
                 + f_3 * pc_z[k] * hsh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, hsg0_1, hsg0_2, hsg0_3, hsg1_1, \
                         hsg1_2, hsg1_3, hsh_3, hsh_5, hsh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * hsg0_1[k]
                 - f_7 * hsg1_1[k]
                 + f_3 * pc_y[k] * hsh_3[k];

        t_7[k] = f_3 * pc_z[k] * hsh_3[k];

        t_8[k] = f_3 * pc_y[k] * hsh_5[k];

        t_9[k] = f_6 * hsg0_2[k]
                 - f_7 * hsg1_2[k]
                 + f_3 * pc_z[k] * hsh_5[k];

        t_10[k] = f_8 * hsg0_3[k]
                  - f_9 * hsg1_3[k]
                  + f_3 * pc_y[k] * hsh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, gsh_15, hsg0_5, \
                         hsg1_5, hsh_6, hsh_8, hsh_9, hsh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * hsh_6[k];

        t_12[k] = f_4 * hsg0_5[k]
                  - f_5 * hsg1_5[k]
                  + f_3 * pc_y[k] * hsh_8[k];

        t_13[k] = f_3 * pc_y[k] * hsh_9[k];

        t_14[k] = f_8 * hsg0_5[k]
                  - f_9 * hsg1_5[k]
                  + f_3 * pc_z[k] * hsh_9[k];

        t_15[k] = f_0 * gsh_15[k]
                  + f_3 * pc_x[k] * hsh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, pc_y, pc_z, gsh_17, gsh_18, \
                         gsh_20, hsh_10, hsh_14, hsh_17, hsh_18, \
                         hsh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * hsh_10[k];

        t_17[k] = f_0 * gsh_17[k]
                  + f_3 * pc_x[k] * hsh_17[k];

        t_18[k] = f_0 * gsh_18[k]
                  + f_3 * pc_x[k] * hsh_18[k];

        t_19[k] = f_3 * pc_y[k] * hsh_14[k];

        t_20[k] = f_0 * gsh_20[k]
                  + f_3 * pc_x[k] * hsh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_y, pc_z, hsg0_10, hsg0_12, hsg0_13, \
                         hsg1_10, hsg1_12, hsg1_13, hsh_15, hsh_17, \
                         hsh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * hsg0_10[k]
                  - f_2 * hsg1_10[k]
                  + f_3 * pc_y[k] * hsh_15[k];

        t_22[k] = f_3 * pc_z[k] * hsh_15[k];

        t_23[k] = f_8 * hsg0_12[k]
                  - f_9 * hsg1_12[k]
                  + f_3 * pc_y[k] * hsh_17[k];

        t_24[k] = f_6 * hsg0_13[k]
                  - f_7 * hsg1_13[k]
                  + f_3 * pc_y[k] * hsh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pc_y, pc_z, gsi0_0, gsh_0, \
                         gsi1_0, hsg0_14, hsg1_14, hsh_19, hsh_20, \
                         hsh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_4 * hsg0_14[k]
                  - f_5 * hsg1_14[k]
                  + f_3 * pc_y[k] * hsh_19[k];

        t_26[k] = f_3 * pc_y[k] * hsh_20[k];

        t_27[k] = f_1 * hsg0_14[k]
                  - f_2 * hsg1_14[k]
                  + f_3 * pc_z[k] * hsh_20[k];

        t_28[k] = pa_y[k] * gsi0_0[k]
                  - f_10 * pc_y[k] * gsi1_0[k];

        t_29[k] = f_11 * gsh_0[k]
                  + f_3 * pc_y[k] * hsh_21[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_y, pc_z, gsi0_3, gsi0_5, gsh_1, \
                         gsi1_3, gsi1_5, hsh_21, hsh_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_3 * pc_z[k] * hsh_21[k];

        t_31[k] = pa_y[k] * gsi0_3[k]
                  + f_12 * gsh_1[k]
                  - f_10 * pc_y[k] * gsi1_3[k];

        t_32[k] = f_3 * pc_z[k] * hsh_22[k];

        t_33[k] = pa_y[k] * gsi0_5[k]
                  - f_10 * pc_y[k] * gsi1_5[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_y, pc_y, pc_z, gsi0_6, gsi0_9, gsh_3, \
                         gsh_5, gsi1_6, gsi1_9, hsh_24, hsh_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pa_y[k] * gsi0_6[k]
                  + f_13 * gsh_3[k]
                  - f_10 * pc_y[k] * gsi1_6[k];

        t_35[k] = f_3 * pc_z[k] * hsh_24[k];

        t_36[k] = f_11 * gsh_5[k]
                  + f_3 * pc_y[k] * hsh_26[k];

        t_37[k] = pa_y[k] * gsi0_9[k]
                  - f_10 * pc_y[k] * gsi1_9[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_y, pc_y, pc_z, gsi0_10, gsh_6, gsh_9, \
                         gsi1_10, hsg0_18, hsg1_18, hsh_27, hsh_28, \
                         hsh_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_y[k] * gsi0_10[k]
                  + f_14 * gsh_6[k]
                  - f_10 * pc_y[k] * gsi1_10[k];

        t_39[k] = f_3 * pc_z[k] * hsh_27[k];

        t_40[k] = f_4 * hsg0_18[k]
                  - f_5 * hsg1_18[k]
                  + f_3 * pc_z[k] * hsh_28[k];

        t_41[k] = f_11 * gsh_9[k]
                  + f_3 * pc_y[k] * hsh_30[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pc_x, pc_y, pc_z, gsi0_14, gsh_36, \
                         gsh_38, gsi1_14, hsh_31, hsh_36, hsh_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * gsi0_14[k]
                  - f_10 * pc_y[k] * gsi1_14[k];

        t_43[k] = f_14 * gsh_36[k]
                  + f_3 * pc_x[k] * hsh_36[k];

        t_44[k] = f_3 * pc_z[k] * hsh_31[k];

        t_45[k] = f_14 * gsh_38[k]
                  + f_3 * pc_x[k] * hsh_38[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, gsh_15, gsh_39, gsh_40, gsh_41, \
                         hsg0_25, hsg1_25, hsh_36, hsh_39, hsh_40, \
                         hsh_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_14 * gsh_39[k]
                  + f_3 * pc_x[k] * hsh_39[k];

        t_47[k] = f_14 * gsh_40[k]
                  + f_3 * pc_x[k] * hsh_40[k];

        t_48[k] = f_14 * gsh_41[k]
                  + f_3 * pc_x[k] * hsh_41[k];

        t_49[k] = f_11 * gsh_15[k]
                  + f_1 * hsg0_25[k]
                  - f_2 * hsg1_25[k]
                  + f_3 * pc_y[k] * hsh_36[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_z, hsg0_25, hsg0_26, hsg0_27, hsg1_25, \
                         hsg1_26, hsg1_27, hsh_36, hsh_37, hsh_38, \
                         hsh_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_z[k] * hsh_36[k];

        t_51[k] = f_4 * hsg0_25[k]
                  - f_5 * hsg1_25[k]
                  + f_3 * pc_z[k] * hsh_37[k];

        t_52[k] = f_6 * hsg0_26[k]
                  - f_7 * hsg1_26[k]
                  + f_3 * pc_z[k] * hsh_38[k];

        t_53[k] = f_8 * hsg0_27[k]
                  - f_9 * hsg1_27[k]
                  + f_3 * pc_z[k] * hsh_39[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_y, pa_z, pc_y, pc_z, gsi0_0, gsi0_27, \
                         gsh_20, gsi1_0, gsi1_27, hsh_41, hsh_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * gsh_20[k]
                  + f_3 * pc_y[k] * hsh_41[k];

        t_55[k] = pa_y[k] * gsi0_27[k]
                  - f_10 * pc_y[k] * gsi1_27[k];

        t_56[k] = pa_z[k] * gsi0_0[k]
                  - f_10 * pc_z[k] * gsi1_0[k];

        t_57[k] = f_3 * pc_y[k] * hsh_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_z, pc_y, pc_z, gsi0_3, gsi0_5, gsh_0, \
                         gsh_2, gsi1_3, gsi1_5, hsh_42, hsh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_11 * gsh_0[k]
                  + f_3 * pc_z[k] * hsh_42[k];

        t_59[k] = pa_z[k] * gsi0_3[k]
                  - f_10 * pc_z[k] * gsi1_3[k];

        t_60[k] = f_3 * pc_y[k] * hsh_44[k];

        t_61[k] = pa_z[k] * gsi0_5[k]
                  + f_12 * gsh_2[k]
                  - f_10 * pc_z[k] * gsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_z, pc_y, pc_z, gsi0_6, gsi0_9, gsh_5, \
                         gsi1_6, gsi1_9, hsg0_32, hsg1_32, hsh_46, \
                         hsh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_z[k] * gsi0_6[k]
                  - f_10 * pc_z[k] * gsi1_6[k];

        t_63[k] = f_4 * hsg0_32[k]
                  - f_5 * hsg1_32[k]
                  + f_3 * pc_y[k] * hsh_46[k];

        t_64[k] = f_3 * pc_y[k] * hsh_47[k];

        t_65[k] = pa_z[k] * gsi0_9[k]
                  + f_13 * gsh_5[k]
                  - f_10 * pc_z[k] * gsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_z, pc_y, pc_z, gsi0_10, gsi1_10, hsg0_34, \
                         hsg0_35, hsg1_34, hsg1_35, hsh_49, hsh_50, \
                         hsh_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_z[k] * gsi0_10[k]
                  - f_10 * pc_z[k] * gsi1_10[k];

        t_67[k] = f_6 * hsg0_34[k]
                  - f_7 * hsg1_34[k]
                  + f_3 * pc_y[k] * hsh_49[k];

        t_68[k] = f_4 * hsg0_35[k]
                  - f_5 * hsg1_35[k]
                  + f_3 * pc_y[k] * hsh_50[k];

        t_69[k] = f_3 * pc_y[k] * hsh_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_z, pc_x, pc_z, gsi0_14, gsh_9, gsh_57, \
                         gsh_58, gsh_59, gsi1_14, hsh_57, hsh_58, \
                         hsh_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * gsi0_14[k]
                  + f_14 * gsh_9[k]
                  - f_10 * pc_z[k] * gsi1_14[k];

        t_71[k] = f_14 * gsh_57[k]
                  + f_3 * pc_x[k] * hsh_57[k];

        t_72[k] = f_14 * gsh_58[k]
                  + f_3 * pc_x[k] * hsh_58[k];

        t_73[k] = f_14 * gsh_59[k]
                  + f_3 * pc_x[k] * hsh_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_z, pc_x, pc_y, pc_z, gsi0_21, gsh_60, \
                         gsh_62, gsi1_21, hsh_56, hsh_60, hsh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_14 * gsh_60[k]
                  + f_3 * pc_x[k] * hsh_60[k];

        t_75[k] = f_3 * pc_y[k] * hsh_56[k];

        t_76[k] = f_14 * gsh_62[k]
                  + f_3 * pc_x[k] * hsh_62[k];

        t_77[k] = pa_z[k] * gsi0_21[k]
                  - f_10 * pc_z[k] * gsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, hsg0_41, hsg0_42, hsg0_43, hsg1_41, hsg1_42, \
                         hsg1_43, hsh_58, hsh_59, hsh_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_15 * hsg0_41[k]
                  - f_16 * hsg1_41[k]
                  + f_3 * pc_y[k] * hsh_58[k];

        t_79[k] = f_8 * hsg0_42[k]
                  - f_9 * hsg1_42[k]
                  + f_3 * pc_y[k] * hsh_59[k];

        t_80[k] = f_6 * hsg0_43[k]
                  - f_7 * hsg1_43[k]
                  + f_3 * pc_y[k] * hsh_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, gsh_20, gsh_63, hsg0_44, \
                         hsg0_45, hsg1_44, hsg1_45, hsh_61, hsh_62, \
                         hsh_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * hsg0_44[k]
                  - f_5 * hsg1_44[k]
                  + f_3 * pc_y[k] * hsh_61[k];

        t_82[k] = f_3 * pc_y[k] * hsh_62[k];

        t_83[k] = f_11 * gsh_20[k]
                  + f_1 * hsg0_44[k]
                  - f_2 * hsg1_44[k]
                  + f_3 * pc_z[k] * hsh_62[k];

        t_84[k] = f_13 * gsh_63[k]
                  + f_1 * hsg0_45[k]
                  - f_2 * hsg1_45[k]
                  + f_3 * pc_x[k] * hsh_63[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, gsh_21, gsh_66, hsg0_48, \
                         hsg1_48, hsh_63, hsh_64, hsh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_12 * gsh_21[k]
                  + f_3 * pc_y[k] * hsh_63[k];

        t_86[k] = f_3 * pc_z[k] * hsh_63[k];

        t_87[k] = f_13 * gsh_66[k]
                  + f_8 * hsg0_48[k]
                  - f_9 * hsg1_48[k]
                  + f_3 * pc_x[k] * hsh_66[k];

        t_88[k] = f_3 * pc_z[k] * hsh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pc_x, pc_z, gsh_69, hsg0_45, hsg0_51, hsg1_45, \
                         hsg1_51, hsh_65, hsh_66, hsh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_4 * hsg0_45[k]
                  - f_5 * hsg1_45[k]
                  + f_3 * pc_z[k] * hsh_65[k];

        t_90[k] = f_13 * gsh_69[k]
                  + f_6 * hsg0_51[k]
                  - f_7 * hsg1_51[k]
                  + f_3 * pc_x[k] * hsh_69[k];

        t_91[k] = f_3 * pc_z[k] * hsh_66[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, gsh_26, gsh_73, hsg0_47, \
                         hsg0_55, hsg1_47, hsg1_55, hsh_68, hsh_69, \
                         hsh_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_12 * gsh_26[k]
                  + f_3 * pc_y[k] * hsh_68[k];

        t_93[k] = f_6 * hsg0_47[k]
                  - f_7 * hsg1_47[k]
                  + f_3 * pc_z[k] * hsh_68[k];

        t_94[k] = f_13 * gsh_73[k]
                  + f_4 * hsg0_55[k]
                  - f_5 * hsg1_55[k]
                  + f_3 * pc_x[k] * hsh_73[k];

        t_95[k] = f_3 * pc_z[k] * hsh_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, pc_y, pc_z, gsh_30, gsh_78, hsg0_48, \
                         hsg0_50, hsg1_48, hsg1_50, hsh_70, hsh_72, \
                         hsh_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_4 * hsg0_48[k]
                  - f_5 * hsg1_48[k]
                  + f_3 * pc_z[k] * hsh_70[k];

        t_97[k] = f_12 * gsh_30[k]
                  + f_3 * pc_y[k] * hsh_72[k];

        t_98[k] = f_8 * hsg0_50[k]
                  - f_9 * hsg1_50[k]
                  + f_3 * pc_z[k] * hsh_72[k];

        t_99[k] = f_13 * gsh_78[k]
                  + f_3 * pc_x[k] * hsh_78[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, gsh_80, gsh_81, \
                         gsh_82, gsh_83, hsh_73, hsh_80, hsh_81, hsh_82, \
                         hsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_3 * pc_z[k] * hsh_73[k];

        t_101[k] = f_13 * gsh_80[k]
                   + f_3 * pc_x[k] * hsh_80[k];

        t_102[k] = f_13 * gsh_81[k]
                   + f_3 * pc_x[k] * hsh_81[k];

        t_103[k] = f_13 * gsh_82[k]
                   + f_3 * pc_x[k] * hsh_82[k];

        t_104[k] = f_13 * gsh_83[k]
                   + f_3 * pc_x[k] * hsh_83[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pc_y, pc_z, gsh_36, hsg0_55, hsg0_56, \
                         hsg1_55, hsg1_56, hsh_78, hsh_79, hsh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_12 * gsh_36[k]
                   + f_1 * hsg0_55[k]
                   - f_2 * hsg1_55[k]
                   + f_3 * pc_y[k] * hsh_78[k];

        t_106[k] = f_3 * pc_z[k] * hsh_78[k];

        t_107[k] = f_4 * hsg0_55[k]
                   - f_5 * hsg1_55[k]
                   + f_3 * pc_z[k] * hsh_79[k];

        t_108[k] = f_6 * hsg0_56[k]
                   - f_7 * hsg1_56[k]
                   + f_3 * pc_z[k] * hsh_80[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_y, pc_y, pc_z, gsi0_56, gsh_41, \
                         gsi1_56, hsg0_57, hsg0_59, hsg1_57, hsg1_59, hsh_81, \
                         hsh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_8 * hsg0_57[k]
                   - f_9 * hsg1_57[k]
                   + f_3 * pc_z[k] * hsh_81[k];

        t_110[k] = f_12 * gsh_41[k]
                   + f_3 * pc_y[k] * hsh_83[k];

        t_111[k] = f_1 * hsg0_59[k]
                   - f_2 * hsg1_59[k]
                   + f_3 * pc_z[k] * hsh_83[k];

        t_112[k] = pa_y[k] * gsi0_56[k]
                   - f_10 * pc_y[k] * gsi1_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_z, pc_y, pc_z, gsi0_31, gsh_21, \
                         gsh_42, gsh_44, gsi1_31, hsh_84, hsh_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_11 * gsh_42[k]
                   + f_3 * pc_y[k] * hsh_84[k];

        t_114[k] = f_11 * gsh_21[k]
                   + f_3 * pc_z[k] * hsh_84[k];

        t_115[k] = pa_z[k] * gsi0_31[k]
                   - f_10 * pc_z[k] * gsi1_31[k];

        t_116[k] = f_11 * gsh_44[k]
                   + f_3 * pc_y[k] * hsh_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, gsi0_34, gsi0_61, \
                         gsh_24, gsh_47, gsi1_34, gsi1_61, hsh_87, \
                         hsh_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_y[k] * gsi0_61[k]
                   - f_10 * pc_y[k] * gsi1_61[k];

        t_118[k] = pa_z[k] * gsi0_34[k]
                   - f_10 * pc_z[k] * gsi1_34[k];

        t_119[k] = f_11 * gsh_24[k]
                   + f_3 * pc_z[k] * hsh_87[k];

        t_120[k] = f_11 * gsh_47[k]
                   + f_3 * pc_y[k] * hsh_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pa_z, pc_y, pc_z, gsi0_38, gsi0_65, \
                         gsh_27, gsi1_38, gsi1_65, hsh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pa_y[k] * gsi0_65[k]
                   - f_10 * pc_y[k] * gsi1_65[k];

        t_122[k] = pa_z[k] * gsi0_38[k]
                   - f_10 * pc_z[k] * gsi1_38[k];

        t_123[k] = f_11 * gsh_27[k]
                   + f_3 * pc_z[k] * hsh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pc_x, pc_y, gsi0_68, gsi0_70, \
                         gsh_50, gsh_51, gsh_99, gsi1_68, gsi1_70, hsh_93, \
                         hsh_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pa_y[k] * gsi0_68[k]
                   + f_12 * gsh_50[k]
                   - f_10 * pc_y[k] * gsi1_68[k];

        t_125[k] = f_11 * gsh_51[k]
                   + f_3 * pc_y[k] * hsh_93[k];

        t_126[k] = pa_y[k] * gsi0_70[k]
                   - f_10 * pc_y[k] * gsi1_70[k];

        t_127[k] = f_13 * gsh_99[k]
                   + f_3 * pc_x[k] * hsh_99[k];
    }
}

static auto
compute_prim_hsi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t hsg0, const size_t hsg1,
                                                          const size_t hsh, const size_t ncols,
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);

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

    const auto *gsi0_49 = buffer.data(gsi0 + 49);
    const auto *gsi0_83 = buffer.data(gsi0 + 83);
    const auto *gsi0_84 = buffer.data(gsi0 + 84);
    const auto *gsi0_87 = buffer.data(gsi0 + 87);
    const auto *gsi0_90 = buffer.data(gsi0 + 90);
    const auto *gsi0_94 = buffer.data(gsi0 + 94);
    const auto *gsi0_96 = buffer.data(gsi0 + 96);
    const auto *gsi0_105 = buffer.data(gsi0 + 105);
    const auto *gsi0_140 = buffer.data(gsi0 + 140);
    const auto *gsi0_143 = buffer.data(gsi0 + 143);
    const auto *gsi0_145 = buffer.data(gsi0 + 145);
    const auto *gsi0_146 = buffer.data(gsi0 + 146);
    const auto *gsi0_149 = buffer.data(gsi0 + 149);
    const auto *gsi0_150 = buffer.data(gsi0 + 150);
    const auto *gsi0_152 = buffer.data(gsi0 + 152);
    const auto *gsi0_154 = buffer.data(gsi0 + 154);

    const auto *gsh_36 = buffer.data(gsh + 36);
    const auto *gsh_42 = buffer.data(gsh + 42);
    const auto *gsh_59 = buffer.data(gsh + 59);
    const auto *gsh_60 = buffer.data(gsh + 60);
    const auto *gsh_61 = buffer.data(gsh + 61);
    const auto *gsh_62 = buffer.data(gsh + 62);
    const auto *gsh_63 = buffer.data(gsh + 63);
    const auto *gsh_66 = buffer.data(gsh + 66);
    const auto *gsh_68 = buffer.data(gsh + 68);
    const auto *gsh_69 = buffer.data(gsh + 69);
    const auto *gsh_70 = buffer.data(gsh + 70);
    const auto *gsh_72 = buffer.data(gsh + 72);
    const auto *gsh_78 = buffer.data(gsh + 78);
    const auto *gsh_83 = buffer.data(gsh + 83);
    const auto *gsh_84 = buffer.data(gsh + 84);
    const auto *gsh_86 = buffer.data(gsh + 86);
    const auto *gsh_87 = buffer.data(gsh + 87);
    const auto *gsh_89 = buffer.data(gsh + 89);
    const auto *gsh_90 = buffer.data(gsh + 90);
    const auto *gsh_93 = buffer.data(gsh + 93);
    const auto *gsh_99 = buffer.data(gsh + 99);
    const auto *gsh_100 = buffer.data(gsh + 100);
    const auto *gsh_101 = buffer.data(gsh + 101);
    const auto *gsh_102 = buffer.data(gsh + 102);
    const auto *gsh_103 = buffer.data(gsh + 103);
    const auto *gsh_104 = buffer.data(gsh + 104);
    const auto *gsh_105 = buffer.data(gsh + 105);
    const auto *gsh_106 = buffer.data(gsh + 106);
    const auto *gsh_107 = buffer.data(gsh + 107);
    const auto *gsh_108 = buffer.data(gsh + 108);
    const auto *gsh_110 = buffer.data(gsh + 110);
    const auto *gsh_111 = buffer.data(gsh + 111);
    const auto *gsh_113 = buffer.data(gsh + 113);
    const auto *gsh_114 = buffer.data(gsh + 114);
    const auto *gsh_119 = buffer.data(gsh + 119);
    const auto *gsh_120 = buffer.data(gsh + 120);
    const auto *gsh_121 = buffer.data(gsh + 121);
    const auto *gsh_122 = buffer.data(gsh + 122);
    const auto *gsh_123 = buffer.data(gsh + 123);
    const auto *gsh_125 = buffer.data(gsh + 125);
    const auto *gsh_126 = buffer.data(gsh + 126);
    const auto *gsh_129 = buffer.data(gsh + 129);
    const auto *gsh_132 = buffer.data(gsh + 132);
    const auto *gsh_136 = buffer.data(gsh + 136);
    const auto *gsh_141 = buffer.data(gsh + 141);
    const auto *gsh_143 = buffer.data(gsh + 143);
    const auto *gsh_144 = buffer.data(gsh + 144);
    const auto *gsh_145 = buffer.data(gsh + 145);
    const auto *gsh_146 = buffer.data(gsh + 146);
    const auto *gsh_152 = buffer.data(gsh + 152);
    const auto *gsh_156 = buffer.data(gsh + 156);
    const auto *gsh_161 = buffer.data(gsh + 161);
    const auto *gsh_162 = buffer.data(gsh + 162);
    const auto *gsh_163 = buffer.data(gsh + 163);
    const auto *gsh_164 = buffer.data(gsh + 164);
    const auto *gsh_165 = buffer.data(gsh + 165);
    const auto *gsh_166 = buffer.data(gsh + 166);
    const auto *gsh_167 = buffer.data(gsh + 167);
    const auto *gsh_183 = buffer.data(gsh + 183);
    const auto *gsh_184 = buffer.data(gsh + 184);
    const auto *gsh_185 = buffer.data(gsh + 185);
    const auto *gsh_186 = buffer.data(gsh + 186);
    const auto *gsh_187 = buffer.data(gsh + 187);
    const auto *gsh_188 = buffer.data(gsh + 188);

    const auto *gsi1_49 = buffer.data(gsi1 + 49);
    const auto *gsi1_83 = buffer.data(gsi1 + 83);
    const auto *gsi1_84 = buffer.data(gsi1 + 84);
    const auto *gsi1_87 = buffer.data(gsi1 + 87);
    const auto *gsi1_90 = buffer.data(gsi1 + 90);
    const auto *gsi1_94 = buffer.data(gsi1 + 94);
    const auto *gsi1_96 = buffer.data(gsi1 + 96);
    const auto *gsi1_105 = buffer.data(gsi1 + 105);
    const auto *gsi1_140 = buffer.data(gsi1 + 140);
    const auto *gsi1_143 = buffer.data(gsi1 + 143);
    const auto *gsi1_145 = buffer.data(gsi1 + 145);
    const auto *gsi1_146 = buffer.data(gsi1 + 146);
    const auto *gsi1_149 = buffer.data(gsi1 + 149);
    const auto *gsi1_150 = buffer.data(gsi1 + 150);
    const auto *gsi1_152 = buffer.data(gsi1 + 152);
    const auto *gsi1_154 = buffer.data(gsi1 + 154);

    const auto *hsg0_72 = buffer.data(hsg0 + 72);
    const auto *hsg0_73 = buffer.data(hsg0 + 73);
    const auto *hsg0_74 = buffer.data(hsg0 + 74);
    const auto *hsg0_75 = buffer.data(hsg0 + 75);
    const auto *hsg0_76 = buffer.data(hsg0 + 76);
    const auto *hsg0_77 = buffer.data(hsg0 + 77);
    const auto *hsg0_78 = buffer.data(hsg0 + 78);
    const auto *hsg0_79 = buffer.data(hsg0 + 79);
    const auto *hsg0_80 = buffer.data(hsg0 + 80);
    const auto *hsg0_84 = buffer.data(hsg0 + 84);
    const auto *hsg0_85 = buffer.data(hsg0 + 85);
    const auto *hsg0_86 = buffer.data(hsg0 + 86);
    const auto *hsg0_87 = buffer.data(hsg0 + 87);
    const auto *hsg0_88 = buffer.data(hsg0 + 88);
    const auto *hsg0_89 = buffer.data(hsg0 + 89);
    const auto *hsg0_90 = buffer.data(hsg0 + 90);
    const auto *hsg0_92 = buffer.data(hsg0 + 92);
    const auto *hsg0_93 = buffer.data(hsg0 + 93);
    const auto *hsg0_95 = buffer.data(hsg0 + 95);
    const auto *hsg0_96 = buffer.data(hsg0 + 96);
    const auto *hsg0_100 = buffer.data(hsg0 + 100);
    const auto *hsg0_101 = buffer.data(hsg0 + 101);
    const auto *hsg0_102 = buffer.data(hsg0 + 102);
    const auto *hsg0_104 = buffer.data(hsg0 + 104);
    const auto *hsg0_110 = buffer.data(hsg0 + 110);
    const auto *hsg0_114 = buffer.data(hsg0 + 114);
    const auto *hsg0_117 = buffer.data(hsg0 + 117);
    const auto *hsg0_118 = buffer.data(hsg0 + 118);
    const auto *hsg0_119 = buffer.data(hsg0 + 119);
    const auto *hsg0_130 = buffer.data(hsg0 + 130);
    const auto *hsg0_132 = buffer.data(hsg0 + 132);

    const auto *hsg1_72 = buffer.data(hsg1 + 72);
    const auto *hsg1_73 = buffer.data(hsg1 + 73);
    const auto *hsg1_74 = buffer.data(hsg1 + 74);
    const auto *hsg1_75 = buffer.data(hsg1 + 75);
    const auto *hsg1_76 = buffer.data(hsg1 + 76);
    const auto *hsg1_77 = buffer.data(hsg1 + 77);
    const auto *hsg1_78 = buffer.data(hsg1 + 78);
    const auto *hsg1_79 = buffer.data(hsg1 + 79);
    const auto *hsg1_80 = buffer.data(hsg1 + 80);
    const auto *hsg1_84 = buffer.data(hsg1 + 84);
    const auto *hsg1_85 = buffer.data(hsg1 + 85);
    const auto *hsg1_86 = buffer.data(hsg1 + 86);
    const auto *hsg1_87 = buffer.data(hsg1 + 87);
    const auto *hsg1_88 = buffer.data(hsg1 + 88);
    const auto *hsg1_89 = buffer.data(hsg1 + 89);
    const auto *hsg1_90 = buffer.data(hsg1 + 90);
    const auto *hsg1_92 = buffer.data(hsg1 + 92);
    const auto *hsg1_93 = buffer.data(hsg1 + 93);
    const auto *hsg1_95 = buffer.data(hsg1 + 95);
    const auto *hsg1_96 = buffer.data(hsg1 + 96);
    const auto *hsg1_100 = buffer.data(hsg1 + 100);
    const auto *hsg1_101 = buffer.data(hsg1 + 101);
    const auto *hsg1_102 = buffer.data(hsg1 + 102);
    const auto *hsg1_104 = buffer.data(hsg1 + 104);
    const auto *hsg1_110 = buffer.data(hsg1 + 110);
    const auto *hsg1_114 = buffer.data(hsg1 + 114);
    const auto *hsg1_117 = buffer.data(hsg1 + 117);
    const auto *hsg1_118 = buffer.data(hsg1 + 118);
    const auto *hsg1_119 = buffer.data(hsg1 + 119);
    const auto *hsg1_130 = buffer.data(hsg1 + 130);
    const auto *hsg1_132 = buffer.data(hsg1 + 132);

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
    const auto *hsh_109 = buffer.data(hsh + 109);
    const auto *hsh_110 = buffer.data(hsh + 110);
    const auto *hsh_111 = buffer.data(hsh + 111);
    const auto *hsh_112 = buffer.data(hsh + 112);
    const auto *hsh_113 = buffer.data(hsh + 113);
    const auto *hsh_114 = buffer.data(hsh + 114);
    const auto *hsh_119 = buffer.data(hsh + 119);
    const auto *hsh_120 = buffer.data(hsh + 120);
    const auto *hsh_121 = buffer.data(hsh + 121);
    const auto *hsh_122 = buffer.data(hsh + 122);
    const auto *hsh_123 = buffer.data(hsh + 123);
    const auto *hsh_124 = buffer.data(hsh + 124);
    const auto *hsh_125 = buffer.data(hsh + 125);
    const auto *hsh_126 = buffer.data(hsh + 126);
    const auto *hsh_127 = buffer.data(hsh + 127);
    const auto *hsh_128 = buffer.data(hsh + 128);
    const auto *hsh_129 = buffer.data(hsh + 129);
    const auto *hsh_131 = buffer.data(hsh + 131);
    const auto *hsh_132 = buffer.data(hsh + 132);
    const auto *hsh_133 = buffer.data(hsh + 133);
    const auto *hsh_135 = buffer.data(hsh + 135);
    const auto *hsh_136 = buffer.data(hsh + 136);
    const auto *hsh_141 = buffer.data(hsh + 141);
    const auto *hsh_142 = buffer.data(hsh + 142);
    const auto *hsh_143 = buffer.data(hsh + 143);
    const auto *hsh_144 = buffer.data(hsh + 144);
    const auto *hsh_145 = buffer.data(hsh + 145);
    const auto *hsh_146 = buffer.data(hsh + 146);
    const auto *hsh_147 = buffer.data(hsh + 147);
    const auto *hsh_149 = buffer.data(hsh + 149);
    const auto *hsh_150 = buffer.data(hsh + 150);
    const auto *hsh_152 = buffer.data(hsh + 152);
    const auto *hsh_153 = buffer.data(hsh + 153);
    const auto *hsh_156 = buffer.data(hsh + 156);
    const auto *hsh_161 = buffer.data(hsh + 161);
    const auto *hsh_162 = buffer.data(hsh + 162);
    const auto *hsh_163 = buffer.data(hsh + 163);
    const auto *hsh_164 = buffer.data(hsh + 164);
    const auto *hsh_165 = buffer.data(hsh + 165);
    const auto *hsh_166 = buffer.data(hsh + 166);
    const auto *hsh_167 = buffer.data(hsh + 167);
    const auto *hsh_168 = buffer.data(hsh + 168);
    const auto *hsh_170 = buffer.data(hsh + 170);
    const auto *hsh_171 = buffer.data(hsh + 171);
    const auto *hsh_173 = buffer.data(hsh + 173);
    const auto *hsh_174 = buffer.data(hsh + 174);
    const auto *hsh_177 = buffer.data(hsh + 177);
    const auto *hsh_183 = buffer.data(hsh + 183);
    const auto *hsh_184 = buffer.data(hsh + 184);
    const auto *hsh_185 = buffer.data(hsh + 185);
    const auto *hsh_186 = buffer.data(hsh + 186);
    const auto *hsh_187 = buffer.data(hsh + 187);
    const auto *hsh_188 = buffer.data(hsh + 188);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, t_132, pc_x, gsh_100, gsh_101, gsh_102, \
                         gsh_103, gsh_104, hsh_100, hsh_101, hsh_102, hsh_103, \
                         hsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_13 * gsh_100[k]
                   + f_3 * pc_x[k] * hsh_100[k];

        t_129[k] = f_13 * gsh_101[k]
                   + f_3 * pc_x[k] * hsh_101[k];

        t_130[k] = f_13 * gsh_102[k]
                   + f_3 * pc_x[k] * hsh_102[k];

        t_131[k] = f_13 * gsh_103[k]
                   + f_3 * pc_x[k] * hsh_103[k];

        t_132[k] = f_13 * gsh_104[k]
                   + f_3 * pc_x[k] * hsh_104[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_z, pc_y, pc_z, gsi0_49, gsh_36, gsh_59, \
                         gsi1_49, hsg0_72, hsg1_72, hsh_99, hsh_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = pa_z[k] * gsi0_49[k]
                   - f_10 * pc_z[k] * gsi1_49[k];

        t_134[k] = f_11 * gsh_36[k]
                   + f_3 * pc_z[k] * hsh_99[k];

        t_135[k] = f_11 * gsh_59[k]
                   + f_8 * hsg0_72[k]
                   - f_9 * hsg1_72[k]
                   + f_3 * pc_y[k] * hsh_101[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pc_y, gsh_60, gsh_61, gsh_62, hsg0_73, hsg0_74, \
                         hsg1_73, hsg1_74, hsh_102, hsh_103, hsh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_11 * gsh_60[k]
                   + f_6 * hsg0_73[k]
                   - f_7 * hsg1_73[k]
                   + f_3 * pc_y[k] * hsh_102[k];

        t_137[k] = f_11 * gsh_61[k]
                   + f_4 * hsg0_74[k]
                   - f_5 * hsg1_74[k]
                   + f_3 * pc_y[k] * hsh_103[k];

        t_138[k] = f_11 * gsh_62[k]
                   + f_3 * pc_y[k] * hsh_104[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_y, pc_x, pc_y, pc_z, gsi0_83, gsh_42, \
                         gsh_105, gsi1_83, hsg0_75, hsg1_75, hsh_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * gsi0_83[k]
                   - f_10 * pc_y[k] * gsi1_83[k];

        t_140[k] = f_13 * gsh_105[k]
                   + f_1 * hsg0_75[k]
                   - f_2 * hsg1_75[k]
                   + f_3 * pc_x[k] * hsh_105[k];

        t_141[k] = f_3 * pc_y[k] * hsh_105[k];

        t_142[k] = f_12 * gsh_42[k]
                   + f_3 * pc_z[k] * hsh_105[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pc_x, pc_y, gsh_110, hsg0_75, hsg0_80, hsg1_75, \
                         hsg1_80, hsh_106, hsh_107, hsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * hsg0_75[k]
                   - f_5 * hsg1_75[k]
                   + f_3 * pc_y[k] * hsh_106[k];

        t_144[k] = f_3 * pc_y[k] * hsh_107[k];

        t_145[k] = f_13 * gsh_110[k]
                   + f_8 * hsg0_80[k]
                   - f_9 * hsg1_80[k]
                   + f_3 * pc_x[k] * hsh_110[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pc_y, hsg0_76, hsg0_77, hsg1_76, hsg1_77, \
                         hsh_108, hsh_109, hsh_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_6 * hsg0_76[k]
                   - f_7 * hsg1_76[k]
                   + f_3 * pc_y[k] * hsh_108[k];

        t_147[k] = f_4 * hsg0_77[k]
                   - f_5 * hsg1_77[k]
                   + f_3 * pc_y[k] * hsh_109[k];

        t_148[k] = f_3 * pc_y[k] * hsh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pc_x, pc_y, gsh_114, hsg0_78, hsg0_79, hsg0_84, \
                         hsg1_78, hsg1_79, hsg1_84, hsh_111, hsh_112, \
                         hsh_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_13 * gsh_114[k]
                   + f_6 * hsg0_84[k]
                   - f_7 * hsg1_84[k]
                   + f_3 * pc_x[k] * hsh_114[k];

        t_150[k] = f_8 * hsg0_78[k]
                   - f_9 * hsg1_78[k]
                   + f_3 * pc_y[k] * hsh_111[k];

        t_151[k] = f_6 * hsg0_79[k]
                   - f_7 * hsg1_79[k]
                   + f_3 * pc_y[k] * hsh_112[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pc_x, pc_y, gsh_119, gsh_120, hsg0_80, \
                         hsg0_89, hsg1_80, hsg1_89, hsh_113, hsh_114, hsh_119, \
                         hsh_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * hsg0_80[k]
                   - f_5 * hsg1_80[k]
                   + f_3 * pc_y[k] * hsh_113[k];

        t_153[k] = f_3 * pc_y[k] * hsh_114[k];

        t_154[k] = f_13 * gsh_119[k]
                   + f_4 * hsg0_89[k]
                   - f_5 * hsg1_89[k]
                   + f_3 * pc_x[k] * hsh_119[k];

        t_155[k] = f_13 * gsh_120[k]
                   + f_3 * pc_x[k] * hsh_120[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pc_x, pc_y, gsh_121, gsh_122, \
                         gsh_123, gsh_125, hsh_119, hsh_121, hsh_122, hsh_123, \
                         hsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_13 * gsh_121[k]
                   + f_3 * pc_x[k] * hsh_121[k];

        t_157[k] = f_13 * gsh_122[k]
                   + f_3 * pc_x[k] * hsh_122[k];

        t_158[k] = f_13 * gsh_123[k]
                   + f_3 * pc_x[k] * hsh_123[k];

        t_159[k] = f_3 * pc_y[k] * hsh_119[k];

        t_160[k] = f_13 * gsh_125[k]
                   + f_3 * pc_x[k] * hsh_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pc_y, hsg0_85, hsg0_86, hsg0_87, hsg1_85, \
                         hsg1_86, hsg1_87, hsh_120, hsh_121, hsh_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_1 * hsg0_85[k]
                   - f_2 * hsg1_85[k]
                   + f_3 * pc_y[k] * hsh_120[k];

        t_162[k] = f_15 * hsg0_86[k]
                   - f_16 * hsg1_86[k]
                   + f_3 * pc_y[k] * hsh_121[k];

        t_163[k] = f_8 * hsg0_87[k]
                   - f_9 * hsg1_87[k]
                   + f_3 * pc_y[k] * hsh_122[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pc_y, pc_z, gsh_62, hsg0_88, hsg0_89, \
                         hsg1_88, hsg1_89, hsh_123, hsh_124, hsh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_6 * hsg0_88[k]
                   - f_7 * hsg1_88[k]
                   + f_3 * pc_y[k] * hsh_123[k];

        t_165[k] = f_4 * hsg0_89[k]
                   - f_5 * hsg1_89[k]
                   + f_3 * pc_y[k] * hsh_124[k];

        t_166[k] = f_3 * pc_y[k] * hsh_125[k];

        t_167[k] = f_12 * gsh_62[k]
                   + f_1 * hsg0_89[k]
                   - f_2 * hsg1_89[k]
                   + f_3 * pc_z[k] * hsh_125[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, pc_x, pc_y, pc_z, gsh_63, gsh_126, \
                         gsh_129, hsg0_90, hsg0_93, hsg1_90, hsg1_93, hsh_126, \
                         hsh_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_12 * gsh_126[k]
                   + f_1 * hsg0_90[k]
                   - f_2 * hsg1_90[k]
                   + f_3 * pc_x[k] * hsh_126[k];

        t_169[k] = f_13 * gsh_63[k]
                   + f_3 * pc_y[k] * hsh_126[k];

        t_170[k] = f_3 * pc_z[k] * hsh_126[k];

        t_171[k] = f_12 * gsh_129[k]
                   + f_8 * hsg0_93[k]
                   - f_9 * hsg1_93[k]
                   + f_3 * pc_x[k] * hsh_129[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, t_175, pc_x, pc_z, gsh_132, hsg0_90, hsg0_96, \
                         hsg1_90, hsg1_96, hsh_127, hsh_128, hsh_129, \
                         hsh_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_3 * pc_z[k] * hsh_127[k];

        t_173[k] = f_4 * hsg0_90[k]
                   - f_5 * hsg1_90[k]
                   + f_3 * pc_z[k] * hsh_128[k];

        t_174[k] = f_12 * gsh_132[k]
                   + f_6 * hsg0_96[k]
                   - f_7 * hsg1_96[k]
                   + f_3 * pc_x[k] * hsh_132[k];

        t_175[k] = f_3 * pc_z[k] * hsh_129[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pc_x, pc_y, pc_z, gsh_68, gsh_136, \
                         hsg0_92, hsg0_100, hsg1_92, hsg1_100, hsh_131, hsh_132, \
                         hsh_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = f_13 * gsh_68[k]
                   + f_3 * pc_y[k] * hsh_131[k];

        t_177[k] = f_6 * hsg0_92[k]
                   - f_7 * hsg1_92[k]
                   + f_3 * pc_z[k] * hsh_131[k];

        t_178[k] = f_12 * gsh_136[k]
                   + f_4 * hsg0_100[k]
                   - f_5 * hsg1_100[k]
                   + f_3 * pc_x[k] * hsh_136[k];

        t_179[k] = f_3 * pc_z[k] * hsh_132[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pc_x, pc_y, pc_z, gsh_72, gsh_141, \
                         hsg0_93, hsg0_95, hsg1_93, hsg1_95, hsh_133, hsh_135, \
                         hsh_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * hsg0_93[k]
                   - f_5 * hsg1_93[k]
                   + f_3 * pc_z[k] * hsh_133[k];

        t_181[k] = f_13 * gsh_72[k]
                   + f_3 * pc_y[k] * hsh_135[k];

        t_182[k] = f_8 * hsg0_95[k]
                   - f_9 * hsg1_95[k]
                   + f_3 * pc_z[k] * hsh_135[k];

        t_183[k] = f_12 * gsh_141[k]
                   + f_3 * pc_x[k] * hsh_141[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, t_188, pc_x, pc_z, gsh_143, gsh_144, \
                         gsh_145, gsh_146, hsh_136, hsh_143, hsh_144, hsh_145, \
                         hsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_3 * pc_z[k] * hsh_136[k];

        t_185[k] = f_12 * gsh_143[k]
                   + f_3 * pc_x[k] * hsh_143[k];

        t_186[k] = f_12 * gsh_144[k]
                   + f_3 * pc_x[k] * hsh_144[k];

        t_187[k] = f_12 * gsh_145[k]
                   + f_3 * pc_x[k] * hsh_145[k];

        t_188[k] = f_12 * gsh_146[k]
                   + f_3 * pc_x[k] * hsh_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_y, pc_z, gsh_78, hsg0_100, hsg0_101, \
                         hsg1_100, hsg1_101, hsh_141, hsh_142, \
                         hsh_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_13 * gsh_78[k]
                   + f_1 * hsg0_100[k]
                   - f_2 * hsg1_100[k]
                   + f_3 * pc_y[k] * hsh_141[k];

        t_190[k] = f_3 * pc_z[k] * hsh_141[k];

        t_191[k] = f_4 * hsg0_100[k]
                   - f_5 * hsg1_100[k]
                   + f_3 * pc_z[k] * hsh_142[k];

        t_192[k] = f_6 * hsg0_101[k]
                   - f_7 * hsg1_101[k]
                   + f_3 * pc_z[k] * hsh_143[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pa_z, pc_y, pc_z, gsi0_84, gsh_83, \
                         gsi1_84, hsg0_102, hsg0_104, hsg1_102, hsg1_104, hsh_144, \
                         hsh_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * hsg0_102[k]
                   - f_9 * hsg1_102[k]
                   + f_3 * pc_z[k] * hsh_144[k];

        t_194[k] = f_13 * gsh_83[k]
                   + f_3 * pc_y[k] * hsh_146[k];

        t_195[k] = f_1 * hsg0_104[k]
                   - f_2 * hsg1_104[k]
                   + f_3 * pc_z[k] * hsh_146[k];

        t_196[k] = pa_z[k] * gsi0_84[k]
                   - f_10 * pc_z[k] * gsi1_84[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_z, pc_y, pc_z, gsi0_87, gsh_63, \
                         gsh_84, gsh_86, gsi1_87, hsh_147, hsh_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * gsh_84[k]
                   + f_3 * pc_y[k] * hsh_147[k];

        t_198[k] = f_11 * gsh_63[k]
                   + f_3 * pc_z[k] * hsh_147[k];

        t_199[k] = pa_z[k] * gsi0_87[k]
                   - f_10 * pc_z[k] * gsi1_87[k];

        t_200[k] = f_12 * gsh_86[k]
                   + f_3 * pc_y[k] * hsh_149[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pa_z, pc_x, pc_z, gsi0_90, gsh_66, gsh_152, \
                         gsi1_90, hsg0_110, hsg1_110, hsh_150, \
                         hsh_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_12 * gsh_152[k]
                   + f_8 * hsg0_110[k]
                   - f_9 * hsg1_110[k]
                   + f_3 * pc_x[k] * hsh_152[k];

        t_202[k] = pa_z[k] * gsi0_90[k]
                   - f_10 * pc_z[k] * gsi1_90[k];

        t_203[k] = f_11 * gsh_66[k]
                   + f_3 * pc_z[k] * hsh_150[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_y, pc_z, gsi0_94, gsh_89, \
                         gsh_156, gsi1_94, hsg0_114, hsg1_114, hsh_152, \
                         hsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_12 * gsh_89[k]
                   + f_3 * pc_y[k] * hsh_152[k];

        t_205[k] = f_12 * gsh_156[k]
                   + f_6 * hsg0_114[k]
                   - f_7 * hsg1_114[k]
                   + f_3 * pc_x[k] * hsh_156[k];

        t_206[k] = pa_z[k] * gsi0_94[k]
                   - f_10 * pc_z[k] * gsi1_94[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_z, pc_y, pc_z, gsi0_96, gsh_69, gsh_70, \
                         gsh_93, gsi1_96, hsh_153, hsh_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_11 * gsh_69[k]
                   + f_3 * pc_z[k] * hsh_153[k];

        t_208[k] = pa_z[k] * gsi0_96[k]
                   + f_12 * gsh_70[k]
                   - f_10 * pc_z[k] * gsi1_96[k];

        t_209[k] = f_12 * gsh_93[k]
                   + f_3 * pc_y[k] * hsh_156[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pc_x, gsh_161, gsh_162, gsh_163, gsh_164, \
                         hsg0_119, hsg1_119, hsh_161, hsh_162, hsh_163, \
                         hsh_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_12 * gsh_161[k]
                   + f_4 * hsg0_119[k]
                   - f_5 * hsg1_119[k]
                   + f_3 * pc_x[k] * hsh_161[k];

        t_211[k] = f_12 * gsh_162[k]
                   + f_3 * pc_x[k] * hsh_162[k];

        t_212[k] = f_12 * gsh_163[k]
                   + f_3 * pc_x[k] * hsh_163[k];

        t_213[k] = f_12 * gsh_164[k]
                   + f_3 * pc_x[k] * hsh_164[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_z, pc_x, pc_z, gsi0_105, gsh_165, \
                         gsh_166, gsh_167, gsi1_105, hsh_165, hsh_166, \
                         hsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_12 * gsh_165[k]
                   + f_3 * pc_x[k] * hsh_165[k];

        t_215[k] = f_12 * gsh_166[k]
                   + f_3 * pc_x[k] * hsh_166[k];

        t_216[k] = f_12 * gsh_167[k]
                   + f_3 * pc_x[k] * hsh_167[k];

        t_217[k] = pa_z[k] * gsi0_105[k]
                   - f_10 * pc_z[k] * gsi1_105[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pc_y, pc_z, gsh_78, gsh_101, gsh_102, hsg0_117, \
                         hsg0_118, hsg1_117, hsg1_118, hsh_162, hsh_164, \
                         hsh_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_11 * gsh_78[k]
                   + f_3 * pc_z[k] * hsh_162[k];

        t_219[k] = f_12 * gsh_101[k]
                   + f_8 * hsg0_117[k]
                   - f_9 * hsg1_117[k]
                   + f_3 * pc_y[k] * hsh_164[k];

        t_220[k] = f_12 * gsh_102[k]
                   + f_6 * hsg0_118[k]
                   - f_7 * hsg1_118[k]
                   + f_3 * pc_y[k] * hsh_165[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, t_224, pa_y, pc_y, pc_z, gsi0_140, gsh_83, \
                         gsh_103, gsh_104, gsi1_140, hsg0_119, hsg1_119, hsh_166, \
                         hsh_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_12 * gsh_103[k]
                   + f_4 * hsg0_119[k]
                   - f_5 * hsg1_119[k]
                   + f_3 * pc_y[k] * hsh_166[k];

        t_222[k] = f_12 * gsh_104[k]
                   + f_3 * pc_y[k] * hsh_167[k];

        t_223[k] = f_11 * gsh_83[k]
                   + f_1 * hsg0_119[k]
                   - f_2 * hsg1_119[k]
                   + f_3 * pc_z[k] * hsh_167[k];

        t_224[k] = pa_y[k] * gsi0_140[k]
                   - f_10 * pc_y[k] * gsi1_140[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_y, pc_y, pc_z, gsi0_143, gsh_84, \
                         gsh_105, gsh_106, gsh_107, gsi1_143, hsh_168, \
                         hsh_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_11 * gsh_105[k]
                   + f_3 * pc_y[k] * hsh_168[k];

        t_226[k] = f_12 * gsh_84[k]
                   + f_3 * pc_z[k] * hsh_168[k];

        t_227[k] = pa_y[k] * gsi0_143[k]
                   + f_12 * gsh_106[k]
                   - f_10 * pc_y[k] * gsi1_143[k];

        t_228[k] = f_11 * gsh_107[k]
                   + f_3 * pc_y[k] * hsh_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_y, pc_z, gsi0_145, gsi0_146, \
                         gsh_87, gsh_108, gsh_110, gsi1_145, gsi1_146, hsh_171, \
                         hsh_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_y[k] * gsi0_145[k]
                   - f_10 * pc_y[k] * gsi1_145[k];

        t_230[k] = pa_y[k] * gsi0_146[k]
                   + f_13 * gsh_108[k]
                   - f_10 * pc_y[k] * gsi1_146[k];

        t_231[k] = f_12 * gsh_87[k]
                   + f_3 * pc_z[k] * hsh_171[k];

        t_232[k] = f_11 * gsh_110[k]
                   + f_3 * pc_y[k] * hsh_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pa_y, pc_y, pc_z, gsi0_149, gsi0_150, gsh_90, \
                         gsh_111, gsi1_149, gsi1_150, hsh_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = pa_y[k] * gsi0_149[k]
                   - f_10 * pc_y[k] * gsi1_149[k];

        t_234[k] = pa_y[k] * gsi0_150[k]
                   + f_14 * gsh_111[k]
                   - f_10 * pc_y[k] * gsi1_150[k];

        t_235[k] = f_12 * gsh_90[k]
                   + f_3 * pc_z[k] * hsh_174[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pa_y, pc_x, pc_y, gsi0_152, gsi0_154, \
                         gsh_113, gsh_114, gsh_183, gsi1_152, gsi1_154, hsh_177, \
                         hsh_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_y[k] * gsi0_152[k]
                   + f_12 * gsh_113[k]
                   - f_10 * pc_y[k] * gsi1_152[k];

        t_237[k] = f_11 * gsh_114[k]
                   + f_3 * pc_y[k] * hsh_177[k];

        t_238[k] = pa_y[k] * gsi0_154[k]
                   - f_10 * pc_y[k] * gsi1_154[k];

        t_239[k] = f_12 * gsh_183[k]
                   + f_3 * pc_x[k] * hsh_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, gsh_184, gsh_185, gsh_186, \
                         gsh_187, gsh_188, hsh_184, hsh_185, hsh_186, hsh_187, \
                         hsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_12 * gsh_184[k]
                   + f_3 * pc_x[k] * hsh_184[k];

        t_241[k] = f_12 * gsh_185[k]
                   + f_3 * pc_x[k] * hsh_185[k];

        t_242[k] = f_12 * gsh_186[k]
                   + f_3 * pc_x[k] * hsh_186[k];

        t_243[k] = f_12 * gsh_187[k]
                   + f_3 * pc_x[k] * hsh_187[k];

        t_244[k] = f_12 * gsh_188[k]
                   + f_3 * pc_x[k] * hsh_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, pc_z, gsh_99, gsh_120, gsh_122, hsg0_130, \
                         hsg0_132, hsg1_130, hsg1_132, hsh_183, \
                         hsh_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_11 * gsh_120[k]
                   + f_1 * hsg0_130[k]
                   - f_2 * hsg1_130[k]
                   + f_3 * pc_y[k] * hsh_183[k];

        t_246[k] = f_12 * gsh_99[k]
                   + f_3 * pc_z[k] * hsh_183[k];

        t_247[k] = f_11 * gsh_122[k]
                   + f_8 * hsg0_132[k]
                   - f_9 * hsg1_132[k]
                   + f_3 * pc_y[k] * hsh_185[k];
    }
}

static auto
compute_prim_hsi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t hsg0, const size_t hsg1,
                                                          const size_t hsh, const size_t ncols,
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 3.0 / q;

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
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsi0_167 = buffer.data(gsi0 + 167);
    const auto *gsi0_168 = buffer.data(gsi0 + 168);
    const auto *gsi0_171 = buffer.data(gsi0 + 171);
    const auto *gsi0_174 = buffer.data(gsi0 + 174);
    const auto *gsi0_178 = buffer.data(gsi0 + 178);
    const auto *gsi0_252 = buffer.data(gsi0 + 252);
    const auto *gsi0_280 = buffer.data(gsi0 + 280);
    const auto *gsi0_283 = buffer.data(gsi0 + 283);
    const auto *gsi0_286 = buffer.data(gsi0 + 286);
    const auto *gsi0_290 = buffer.data(gsi0 + 290);
    const auto *gsi0_301 = buffer.data(gsi0 + 301);
    const auto *gsi0_303 = buffer.data(gsi0 + 303);
    const auto *gsi0_304 = buffer.data(gsi0 + 304);
    const auto *gsi0_305 = buffer.data(gsi0 + 305);
    const auto *gsi0_307 = buffer.data(gsi0 + 307);
    const auto *gsi0_313 = buffer.data(gsi0 + 313);
    const auto *gsi0_317 = buffer.data(gsi0 + 317);
    const auto *gsi0_320 = buffer.data(gsi0 + 320);
    const auto *gsi0_322 = buffer.data(gsi0 + 322);
    const auto *gsi0_329 = buffer.data(gsi0 + 329);
    const auto *gsi0_331 = buffer.data(gsi0 + 331);
    const auto *gsi0_332 = buffer.data(gsi0 + 332);
    const auto *gsi0_333 = buffer.data(gsi0 + 333);
    const auto *gsi0_335 = buffer.data(gsi0 + 335);
    const auto *gsi0_336 = buffer.data(gsi0 + 336);
    const auto *gsi0_339 = buffer.data(gsi0 + 339);
    const auto *gsi0_341 = buffer.data(gsi0 + 341);
    const auto *gsi0_342 = buffer.data(gsi0 + 342);
    const auto *gsi0_345 = buffer.data(gsi0 + 345);
    const auto *gsi0_346 = buffer.data(gsi0 + 346);
    const auto *gsi0_348 = buffer.data(gsi0 + 348);
    const auto *gsi0_350 = buffer.data(gsi0 + 350);
    const auto *gsi0_357 = buffer.data(gsi0 + 357);
    const auto *gsi0_359 = buffer.data(gsi0 + 359);
    const auto *gsi0_360 = buffer.data(gsi0 + 360);
    const auto *gsi0_361 = buffer.data(gsi0 + 361);
    const auto *gsi0_363 = buffer.data(gsi0 + 363);
    const auto *gsi0_367 = buffer.data(gsi0 + 367);

    const auto *gsh_105 = buffer.data(gsh + 105);
    const auto *gsh_123 = buffer.data(gsh + 123);
    const auto *gsh_124 = buffer.data(gsh + 124);
    const auto *gsh_125 = buffer.data(gsh + 125);
    const auto *gsh_126 = buffer.data(gsh + 126);
    const auto *gsh_129 = buffer.data(gsh + 129);
    const auto *gsh_131 = buffer.data(gsh + 131);
    const auto *gsh_132 = buffer.data(gsh + 132);
    const auto *gsh_135 = buffer.data(gsh + 135);
    const auto *gsh_141 = buffer.data(gsh + 141);
    const auto *gsh_146 = buffer.data(gsh + 146);
    const auto *gsh_147 = buffer.data(gsh + 147);
    const auto *gsh_149 = buffer.data(gsh + 149);
    const auto *gsh_150 = buffer.data(gsh + 150);
    const auto *gsh_152 = buffer.data(gsh + 152);
    const auto *gsh_153 = buffer.data(gsh + 153);
    const auto *gsh_156 = buffer.data(gsh + 156);
    const auto *gsh_162 = buffer.data(gsh + 162);
    const auto *gsh_167 = buffer.data(gsh + 167);
    const auto *gsh_168 = buffer.data(gsh + 168);
    const auto *gsh_170 = buffer.data(gsh + 170);
    const auto *gsh_173 = buffer.data(gsh + 173);
    const auto *gsh_177 = buffer.data(gsh + 177);
    const auto *gsh_188 = buffer.data(gsh + 188);
    const auto *gsh_189 = buffer.data(gsh + 189);
    const auto *gsh_191 = buffer.data(gsh + 191);
    const auto *gsh_194 = buffer.data(gsh + 194);
    const auto *gsh_198 = buffer.data(gsh + 198);
    const auto *gsh_203 = buffer.data(gsh + 203);
    const auto *gsh_204 = buffer.data(gsh + 204);
    const auto *gsh_205 = buffer.data(gsh + 205);
    const auto *gsh_206 = buffer.data(gsh + 206);
    const auto *gsh_207 = buffer.data(gsh + 207);
    const auto *gsh_209 = buffer.data(gsh + 209);
    const auto *gsh_210 = buffer.data(gsh + 210);
    const auto *gsh_213 = buffer.data(gsh + 213);
    const auto *gsh_216 = buffer.data(gsh + 216);
    const auto *gsh_220 = buffer.data(gsh + 220);
    const auto *gsh_225 = buffer.data(gsh + 225);
    const auto *gsh_227 = buffer.data(gsh + 227);
    const auto *gsh_228 = buffer.data(gsh + 228);
    const auto *gsh_229 = buffer.data(gsh + 229);
    const auto *gsh_230 = buffer.data(gsh + 230);
    const auto *gsh_236 = buffer.data(gsh + 236);
    const auto *gsh_240 = buffer.data(gsh + 240);
    const auto *gsh_243 = buffer.data(gsh + 243);
    const auto *gsh_245 = buffer.data(gsh + 245);
    const auto *gsh_246 = buffer.data(gsh + 246);
    const auto *gsh_247 = buffer.data(gsh + 247);
    const auto *gsh_248 = buffer.data(gsh + 248);
    const auto *gsh_249 = buffer.data(gsh + 249);
    const auto *gsh_250 = buffer.data(gsh + 250);
    const auto *gsh_251 = buffer.data(gsh + 251);
    const auto *gsh_252 = buffer.data(gsh + 252);
    const auto *gsh_255 = buffer.data(gsh + 255);
    const auto *gsh_257 = buffer.data(gsh + 257);
    const auto *gsh_258 = buffer.data(gsh + 258);
    const auto *gsh_261 = buffer.data(gsh + 261);
    const auto *gsh_262 = buffer.data(gsh + 262);
    const auto *gsh_264 = buffer.data(gsh + 264);
    const auto *gsh_266 = buffer.data(gsh + 266);
    const auto *gsh_267 = buffer.data(gsh + 267);
    const auto *gsh_268 = buffer.data(gsh + 268);
    const auto *gsh_269 = buffer.data(gsh + 269);
    const auto *gsh_270 = buffer.data(gsh + 270);
    const auto *gsh_271 = buffer.data(gsh + 271);
    const auto *gsh_272 = buffer.data(gsh + 272);
    const auto *gsh_276 = buffer.data(gsh + 276);

    const auto *gsi1_167 = buffer.data(gsi1 + 167);
    const auto *gsi1_168 = buffer.data(gsi1 + 168);
    const auto *gsi1_171 = buffer.data(gsi1 + 171);
    const auto *gsi1_174 = buffer.data(gsi1 + 174);
    const auto *gsi1_178 = buffer.data(gsi1 + 178);
    const auto *gsi1_252 = buffer.data(gsi1 + 252);
    const auto *gsi1_280 = buffer.data(gsi1 + 280);
    const auto *gsi1_283 = buffer.data(gsi1 + 283);
    const auto *gsi1_286 = buffer.data(gsi1 + 286);
    const auto *gsi1_290 = buffer.data(gsi1 + 290);
    const auto *gsi1_301 = buffer.data(gsi1 + 301);
    const auto *gsi1_303 = buffer.data(gsi1 + 303);
    const auto *gsi1_304 = buffer.data(gsi1 + 304);
    const auto *gsi1_305 = buffer.data(gsi1 + 305);
    const auto *gsi1_307 = buffer.data(gsi1 + 307);
    const auto *gsi1_313 = buffer.data(gsi1 + 313);
    const auto *gsi1_317 = buffer.data(gsi1 + 317);
    const auto *gsi1_320 = buffer.data(gsi1 + 320);
    const auto *gsi1_322 = buffer.data(gsi1 + 322);
    const auto *gsi1_329 = buffer.data(gsi1 + 329);
    const auto *gsi1_331 = buffer.data(gsi1 + 331);
    const auto *gsi1_332 = buffer.data(gsi1 + 332);
    const auto *gsi1_333 = buffer.data(gsi1 + 333);
    const auto *gsi1_335 = buffer.data(gsi1 + 335);
    const auto *gsi1_336 = buffer.data(gsi1 + 336);
    const auto *gsi1_339 = buffer.data(gsi1 + 339);
    const auto *gsi1_341 = buffer.data(gsi1 + 341);
    const auto *gsi1_342 = buffer.data(gsi1 + 342);
    const auto *gsi1_345 = buffer.data(gsi1 + 345);
    const auto *gsi1_346 = buffer.data(gsi1 + 346);
    const auto *gsi1_348 = buffer.data(gsi1 + 348);
    const auto *gsi1_350 = buffer.data(gsi1 + 350);
    const auto *gsi1_357 = buffer.data(gsi1 + 357);
    const auto *gsi1_359 = buffer.data(gsi1 + 359);
    const auto *gsi1_360 = buffer.data(gsi1 + 360);
    const auto *gsi1_361 = buffer.data(gsi1 + 361);
    const auto *gsi1_363 = buffer.data(gsi1 + 363);
    const auto *gsi1_367 = buffer.data(gsi1 + 367);

    const auto *hsg0_133 = buffer.data(hsg0 + 133);
    const auto *hsg0_134 = buffer.data(hsg0 + 134);
    const auto *hsg0_135 = buffer.data(hsg0 + 135);
    const auto *hsg0_136 = buffer.data(hsg0 + 136);
    const auto *hsg0_137 = buffer.data(hsg0 + 137);
    const auto *hsg0_138 = buffer.data(hsg0 + 138);
    const auto *hsg0_139 = buffer.data(hsg0 + 139);
    const auto *hsg0_140 = buffer.data(hsg0 + 140);
    const auto *hsg0_144 = buffer.data(hsg0 + 144);
    const auto *hsg0_145 = buffer.data(hsg0 + 145);
    const auto *hsg0_146 = buffer.data(hsg0 + 146);
    const auto *hsg0_147 = buffer.data(hsg0 + 147);
    const auto *hsg0_148 = buffer.data(hsg0 + 148);
    const auto *hsg0_149 = buffer.data(hsg0 + 149);
    const auto *hsg0_150 = buffer.data(hsg0 + 150);
    const auto *hsg0_152 = buffer.data(hsg0 + 152);
    const auto *hsg0_153 = buffer.data(hsg0 + 153);
    const auto *hsg0_155 = buffer.data(hsg0 + 155);

    const auto *hsg1_133 = buffer.data(hsg1 + 133);
    const auto *hsg1_134 = buffer.data(hsg1 + 134);
    const auto *hsg1_135 = buffer.data(hsg1 + 135);
    const auto *hsg1_136 = buffer.data(hsg1 + 136);
    const auto *hsg1_137 = buffer.data(hsg1 + 137);
    const auto *hsg1_138 = buffer.data(hsg1 + 138);
    const auto *hsg1_139 = buffer.data(hsg1 + 139);
    const auto *hsg1_140 = buffer.data(hsg1 + 140);
    const auto *hsg1_144 = buffer.data(hsg1 + 144);
    const auto *hsg1_145 = buffer.data(hsg1 + 145);
    const auto *hsg1_146 = buffer.data(hsg1 + 146);
    const auto *hsg1_147 = buffer.data(hsg1 + 147);
    const auto *hsg1_148 = buffer.data(hsg1 + 148);
    const auto *hsg1_149 = buffer.data(hsg1 + 149);
    const auto *hsg1_150 = buffer.data(hsg1 + 150);
    const auto *hsg1_152 = buffer.data(hsg1 + 152);
    const auto *hsg1_153 = buffer.data(hsg1 + 153);
    const auto *hsg1_155 = buffer.data(hsg1 + 155);

    const auto *hsh_186 = buffer.data(hsh + 186);
    const auto *hsh_187 = buffer.data(hsh + 187);
    const auto *hsh_188 = buffer.data(hsh + 188);
    const auto *hsh_189 = buffer.data(hsh + 189);
    const auto *hsh_190 = buffer.data(hsh + 190);
    const auto *hsh_191 = buffer.data(hsh + 191);
    const auto *hsh_192 = buffer.data(hsh + 192);
    const auto *hsh_193 = buffer.data(hsh + 193);
    const auto *hsh_194 = buffer.data(hsh + 194);
    const auto *hsh_195 = buffer.data(hsh + 195);
    const auto *hsh_196 = buffer.data(hsh + 196);
    const auto *hsh_197 = buffer.data(hsh + 197);
    const auto *hsh_198 = buffer.data(hsh + 198);
    const auto *hsh_203 = buffer.data(hsh + 203);
    const auto *hsh_204 = buffer.data(hsh + 204);
    const auto *hsh_205 = buffer.data(hsh + 205);
    const auto *hsh_206 = buffer.data(hsh + 206);
    const auto *hsh_207 = buffer.data(hsh + 207);
    const auto *hsh_208 = buffer.data(hsh + 208);
    const auto *hsh_209 = buffer.data(hsh + 209);
    const auto *hsh_210 = buffer.data(hsh + 210);
    const auto *hsh_211 = buffer.data(hsh + 211);
    const auto *hsh_212 = buffer.data(hsh + 212);
    const auto *hsh_213 = buffer.data(hsh + 213);
    const auto *hsh_215 = buffer.data(hsh + 215);
    const auto *hsh_216 = buffer.data(hsh + 216);
    const auto *hsh_217 = buffer.data(hsh + 217);
    const auto *hsh_219 = buffer.data(hsh + 219);
    const auto *hsh_220 = buffer.data(hsh + 220);
    const auto *hsh_225 = buffer.data(hsh + 225);
    const auto *hsh_227 = buffer.data(hsh + 227);
    const auto *hsh_228 = buffer.data(hsh + 228);
    const auto *hsh_229 = buffer.data(hsh + 229);
    const auto *hsh_230 = buffer.data(hsh + 230);
    const auto *hsh_231 = buffer.data(hsh + 231);
    const auto *hsh_233 = buffer.data(hsh + 233);
    const auto *hsh_234 = buffer.data(hsh + 234);
    const auto *hsh_236 = buffer.data(hsh + 236);
    const auto *hsh_237 = buffer.data(hsh + 237);
    const auto *hsh_240 = buffer.data(hsh + 240);
    const auto *hsh_246 = buffer.data(hsh + 246);
    const auto *hsh_247 = buffer.data(hsh + 247);
    const auto *hsh_248 = buffer.data(hsh + 248);
    const auto *hsh_249 = buffer.data(hsh + 249);
    const auto *hsh_250 = buffer.data(hsh + 250);
    const auto *hsh_251 = buffer.data(hsh + 251);
    const auto *hsh_252 = buffer.data(hsh + 252);
    const auto *hsh_254 = buffer.data(hsh + 254);
    const auto *hsh_255 = buffer.data(hsh + 255);
    const auto *hsh_257 = buffer.data(hsh + 257);
    const auto *hsh_258 = buffer.data(hsh + 258);
    const auto *hsh_261 = buffer.data(hsh + 261);
    const auto *hsh_267 = buffer.data(hsh + 267);
    const auto *hsh_268 = buffer.data(hsh + 268);
    const auto *hsh_269 = buffer.data(hsh + 269);
    const auto *hsh_270 = buffer.data(hsh + 270);
    const auto *hsh_271 = buffer.data(hsh + 271);
    const auto *hsh_272 = buffer.data(hsh + 272);
    const auto *hsh_273 = buffer.data(hsh + 273);
    const auto *hsh_275 = buffer.data(hsh + 275);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, gsh_123, gsh_124, gsh_125, hsg0_133, \
                         hsg0_134, hsg1_133, hsg1_134, hsh_186, hsh_187, \
                         hsh_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_11 * gsh_123[k]
                   + f_6 * hsg0_133[k]
                   - f_7 * hsg1_133[k]
                   + f_3 * pc_y[k] * hsh_186[k];

        t_249[k] = f_11 * gsh_124[k]
                   + f_4 * hsg0_134[k]
                   - f_5 * hsg1_134[k]
                   + f_3 * pc_y[k] * hsh_187[k];

        t_250[k] = f_11 * gsh_125[k]
                   + f_3 * pc_y[k] * hsh_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pa_y, pc_x, pc_y, pc_z, gsi0_167, \
                         gsh_105, gsh_189, gsi1_167, hsg0_135, hsg1_135, \
                         hsh_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pa_y[k] * gsi0_167[k]
                   - f_10 * pc_y[k] * gsi1_167[k];

        t_252[k] = f_12 * gsh_189[k]
                   + f_1 * hsg0_135[k]
                   - f_2 * hsg1_135[k]
                   + f_3 * pc_x[k] * hsh_189[k];

        t_253[k] = f_3 * pc_y[k] * hsh_189[k];

        t_254[k] = f_13 * gsh_105[k]
                   + f_3 * pc_z[k] * hsh_189[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, gsh_194, hsg0_135, hsg0_140, \
                         hsg1_135, hsg1_140, hsh_190, hsh_191, \
                         hsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_4 * hsg0_135[k]
                   - f_5 * hsg1_135[k]
                   + f_3 * pc_y[k] * hsh_190[k];

        t_256[k] = f_3 * pc_y[k] * hsh_191[k];

        t_257[k] = f_12 * gsh_194[k]
                   + f_8 * hsg0_140[k]
                   - f_9 * hsg1_140[k]
                   + f_3 * pc_x[k] * hsh_194[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_y, hsg0_136, hsg0_137, hsg1_136, hsg1_137, \
                         hsh_192, hsh_193, hsh_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_6 * hsg0_136[k]
                   - f_7 * hsg1_136[k]
                   + f_3 * pc_y[k] * hsh_192[k];

        t_259[k] = f_4 * hsg0_137[k]
                   - f_5 * hsg1_137[k]
                   + f_3 * pc_y[k] * hsh_193[k];

        t_260[k] = f_3 * pc_y[k] * hsh_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_y, gsh_198, hsg0_138, hsg0_139, \
                         hsg0_144, hsg1_138, hsg1_139, hsg1_144, hsh_195, hsh_196, \
                         hsh_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_12 * gsh_198[k]
                   + f_6 * hsg0_144[k]
                   - f_7 * hsg1_144[k]
                   + f_3 * pc_x[k] * hsh_198[k];

        t_262[k] = f_8 * hsg0_138[k]
                   - f_9 * hsg1_138[k]
                   + f_3 * pc_y[k] * hsh_195[k];

        t_263[k] = f_6 * hsg0_139[k]
                   - f_7 * hsg1_139[k]
                   + f_3 * pc_y[k] * hsh_196[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, gsh_203, gsh_204, hsg0_140, \
                         hsg0_149, hsg1_140, hsg1_149, hsh_197, hsh_198, hsh_203, \
                         hsh_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * hsg0_140[k]
                   - f_5 * hsg1_140[k]
                   + f_3 * pc_y[k] * hsh_197[k];

        t_265[k] = f_3 * pc_y[k] * hsh_198[k];

        t_266[k] = f_12 * gsh_203[k]
                   + f_4 * hsg0_149[k]
                   - f_5 * hsg1_149[k]
                   + f_3 * pc_x[k] * hsh_203[k];

        t_267[k] = f_12 * gsh_204[k]
                   + f_3 * pc_x[k] * hsh_204[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, t_272, pc_x, pc_y, gsh_205, gsh_206, \
                         gsh_207, gsh_209, hsh_203, hsh_205, hsh_206, hsh_207, \
                         hsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * gsh_205[k]
                   + f_3 * pc_x[k] * hsh_205[k];

        t_269[k] = f_12 * gsh_206[k]
                   + f_3 * pc_x[k] * hsh_206[k];

        t_270[k] = f_12 * gsh_207[k]
                   + f_3 * pc_x[k] * hsh_207[k];

        t_271[k] = f_3 * pc_y[k] * hsh_203[k];

        t_272[k] = f_12 * gsh_209[k]
                   + f_3 * pc_x[k] * hsh_209[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_y, hsg0_145, hsg0_146, hsg0_147, hsg1_145, \
                         hsg1_146, hsg1_147, hsh_204, hsh_205, \
                         hsh_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * hsg0_145[k]
                   - f_2 * hsg1_145[k]
                   + f_3 * pc_y[k] * hsh_204[k];

        t_274[k] = f_15 * hsg0_146[k]
                   - f_16 * hsg1_146[k]
                   + f_3 * pc_y[k] * hsh_205[k];

        t_275[k] = f_8 * hsg0_147[k]
                   - f_9 * hsg1_147[k]
                   + f_3 * pc_y[k] * hsh_206[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, gsh_125, hsg0_148, hsg0_149, \
                         hsg1_148, hsg1_149, hsh_207, hsh_208, \
                         hsh_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * hsg0_148[k]
                   - f_7 * hsg1_148[k]
                   + f_3 * pc_y[k] * hsh_207[k];

        t_277[k] = f_4 * hsg0_149[k]
                   - f_5 * hsg1_149[k]
                   + f_3 * pc_y[k] * hsh_208[k];

        t_278[k] = f_3 * pc_y[k] * hsh_209[k];

        t_279[k] = f_13 * gsh_125[k]
                   + f_1 * hsg0_149[k]
                   - f_2 * hsg1_149[k]
                   + f_3 * pc_z[k] * hsh_209[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_x, pc_x, pc_y, pc_z, gsi0_280, \
                         gsi0_283, gsh_126, gsh_210, gsh_213, gsi1_280, gsi1_283, \
                         hsh_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = pa_x[k] * gsi0_280[k]
                   + f_17 * gsh_210[k]
                   - f_10 * pc_x[k] * gsi1_280[k];

        t_281[k] = f_14 * gsh_126[k]
                   + f_3 * pc_y[k] * hsh_210[k];

        t_282[k] = f_3 * pc_z[k] * hsh_210[k];

        t_283[k] = pa_x[k] * gsi0_283[k]
                   + f_14 * gsh_213[k]
                   - f_10 * pc_x[k] * gsi1_283[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_x, pc_x, pc_z, gsi0_286, gsh_216, \
                         gsi1_286, hsg0_150, hsg1_150, hsh_211, hsh_212, \
                         hsh_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_3 * pc_z[k] * hsh_211[k];

        t_285[k] = f_4 * hsg0_150[k]
                   - f_5 * hsg1_150[k]
                   + f_3 * pc_z[k] * hsh_212[k];

        t_286[k] = pa_x[k] * gsi0_286[k]
                   + f_13 * gsh_216[k]
                   - f_10 * pc_x[k] * gsi1_286[k];

        t_287[k] = f_3 * pc_z[k] * hsh_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, pc_y, pc_z, gsi0_290, \
                         gsh_131, gsh_220, gsi1_290, hsg0_152, hsg1_152, hsh_215, \
                         hsh_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_14 * gsh_131[k]
                   + f_3 * pc_y[k] * hsh_215[k];

        t_289[k] = f_6 * hsg0_152[k]
                   - f_7 * hsg1_152[k]
                   + f_3 * pc_z[k] * hsh_215[k];

        t_290[k] = pa_x[k] * gsi0_290[k]
                   + f_12 * gsh_220[k]
                   - f_10 * pc_x[k] * gsi1_290[k];

        t_291[k] = f_3 * pc_z[k] * hsh_216[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pc_x, pc_y, pc_z, gsh_135, gsh_225, \
                         hsg0_153, hsg0_155, hsg1_153, hsg1_155, hsh_217, hsh_219, \
                         hsh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * hsg0_153[k]
                   - f_5 * hsg1_153[k]
                   + f_3 * pc_z[k] * hsh_217[k];

        t_293[k] = f_14 * gsh_135[k]
                   + f_3 * pc_y[k] * hsh_219[k];

        t_294[k] = f_8 * hsg0_155[k]
                   - f_9 * hsg1_155[k]
                   + f_3 * pc_z[k] * hsh_219[k];

        t_295[k] = f_11 * gsh_225[k]
                   + f_3 * pc_x[k] * hsh_225[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, t_300, pc_x, pc_z, gsh_227, gsh_228, \
                         gsh_229, gsh_230, hsh_220, hsh_227, hsh_228, hsh_229, \
                         hsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_3 * pc_z[k] * hsh_220[k];

        t_297[k] = f_11 * gsh_227[k]
                   + f_3 * pc_x[k] * hsh_227[k];

        t_298[k] = f_11 * gsh_228[k]
                   + f_3 * pc_x[k] * hsh_228[k];

        t_299[k] = f_11 * gsh_229[k]
                   + f_3 * pc_x[k] * hsh_229[k];

        t_300[k] = f_11 * gsh_230[k]
                   + f_3 * pc_x[k] * hsh_230[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_x, pc_x, pc_z, gsi0_301, gsi0_303, \
                         gsi0_304, gsi1_301, gsi1_303, gsi1_304, \
                         hsh_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = pa_x[k] * gsi0_301[k]
                   - f_10 * pc_x[k] * gsi1_301[k];

        t_302[k] = f_3 * pc_z[k] * hsh_225[k];

        t_303[k] = pa_x[k] * gsi0_303[k]
                   - f_10 * pc_x[k] * gsi1_303[k];

        t_304[k] = pa_x[k] * gsi0_304[k]
                   - f_10 * pc_x[k] * gsi1_304[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pa_x, pc_x, pc_y, gsi0_305, gsi0_307, gsh_146, \
                         gsi1_305, gsi1_307, hsh_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_x[k] * gsi0_305[k]
                   - f_10 * pc_x[k] * gsi1_305[k];

        t_306[k] = f_14 * gsh_146[k]
                   + f_3 * pc_y[k] * hsh_230[k];

        t_307[k] = pa_x[k] * gsi0_307[k]
                   - f_10 * pc_x[k] * gsi1_307[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_z, pc_y, pc_z, gsi0_168, gsi0_171, \
                         gsh_126, gsh_147, gsi1_168, gsi1_171, \
                         hsh_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pa_z[k] * gsi0_168[k]
                   - f_10 * pc_z[k] * gsi1_168[k];

        t_309[k] = f_13 * gsh_147[k]
                   + f_3 * pc_y[k] * hsh_231[k];

        t_310[k] = f_11 * gsh_126[k]
                   + f_3 * pc_z[k] * hsh_231[k];

        t_311[k] = pa_z[k] * gsi0_171[k]
                   - f_10 * pc_z[k] * gsi1_171[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pa_z, pc_x, pc_y, pc_z, gsi0_174, \
                         gsi0_313, gsh_149, gsh_236, gsi1_174, gsi1_313, \
                         hsh_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_13 * gsh_149[k]
                   + f_3 * pc_y[k] * hsh_233[k];

        t_313[k] = pa_x[k] * gsi0_313[k]
                   + f_14 * gsh_236[k]
                   - f_10 * pc_x[k] * gsi1_313[k];

        t_314[k] = pa_z[k] * gsi0_174[k]
                   - f_10 * pc_z[k] * gsi1_174[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_x, pc_x, pc_y, pc_z, gsi0_317, gsh_129, \
                         gsh_152, gsh_240, gsi1_317, hsh_234, hsh_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_11 * gsh_129[k]
                   + f_3 * pc_z[k] * hsh_234[k];

        t_316[k] = f_13 * gsh_152[k]
                   + f_3 * pc_y[k] * hsh_236[k];

        t_317[k] = pa_x[k] * gsi0_317[k]
                   + f_13 * gsh_240[k]
                   - f_10 * pc_x[k] * gsi1_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_x, pa_z, pc_x, pc_z, gsi0_178, gsi0_320, \
                         gsh_132, gsh_243, gsi1_178, gsi1_320, \
                         hsh_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_z[k] * gsi0_178[k]
                   - f_10 * pc_z[k] * gsi1_178[k];

        t_319[k] = f_11 * gsh_132[k]
                   + f_3 * pc_z[k] * hsh_237[k];

        t_320[k] = pa_x[k] * gsi0_320[k]
                   + f_12 * gsh_243[k]
                   - f_10 * pc_x[k] * gsi1_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pa_x, pc_x, pc_y, gsi0_322, gsh_156, \
                         gsh_245, gsh_246, gsh_247, gsi1_322, hsh_240, hsh_246, \
                         hsh_247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_13 * gsh_156[k]
                   + f_3 * pc_y[k] * hsh_240[k];

        t_322[k] = pa_x[k] * gsi0_322[k]
                   + f_12 * gsh_245[k]
                   - f_10 * pc_x[k] * gsi1_322[k];

        t_323[k] = f_11 * gsh_246[k]
                   + f_3 * pc_x[k] * hsh_246[k];

        t_324[k] = f_11 * gsh_247[k]
                   + f_3 * pc_x[k] * hsh_247[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pc_x, gsh_248, gsh_249, gsh_250, gsh_251, \
                         hsh_248, hsh_249, hsh_250, hsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_11 * gsh_248[k]
                   + f_3 * pc_x[k] * hsh_248[k];

        t_326[k] = f_11 * gsh_249[k]
                   + f_3 * pc_x[k] * hsh_249[k];

        t_327[k] = f_11 * gsh_250[k]
                   + f_3 * pc_x[k] * hsh_250[k];

        t_328[k] = f_11 * gsh_251[k]
                   + f_3 * pc_x[k] * hsh_251[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_x, pc_x, pc_z, gsi0_329, gsi0_331, \
                         gsi0_332, gsh_141, gsi1_329, gsi1_331, gsi1_332, \
                         hsh_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = pa_x[k] * gsi0_329[k]
                   - f_10 * pc_x[k] * gsi1_329[k];

        t_330[k] = f_11 * gsh_141[k]
                   + f_3 * pc_z[k] * hsh_246[k];

        t_331[k] = pa_x[k] * gsi0_331[k]
                   - f_10 * pc_x[k] * gsi1_331[k];

        t_332[k] = pa_x[k] * gsi0_332[k]
                   - f_10 * pc_x[k] * gsi1_332[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pa_x, pc_x, pc_y, gsi0_333, gsi0_335, \
                         gsi0_336, gsh_167, gsh_252, gsi1_333, gsi1_335, gsi1_336, \
                         hsh_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_x[k] * gsi0_333[k]
                   - f_10 * pc_x[k] * gsi1_333[k];

        t_334[k] = f_13 * gsh_167[k]
                   + f_3 * pc_y[k] * hsh_251[k];

        t_335[k] = pa_x[k] * gsi0_335[k]
                   - f_10 * pc_x[k] * gsi1_335[k];

        t_336[k] = pa_x[k] * gsi0_336[k]
                   + f_17 * gsh_252[k]
                   - f_10 * pc_x[k] * gsi1_336[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_x, pc_x, pc_y, pc_z, gsi0_339, \
                         gsh_147, gsh_168, gsh_170, gsh_255, gsi1_339, hsh_252, \
                         hsh_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_12 * gsh_168[k]
                   + f_3 * pc_y[k] * hsh_252[k];

        t_338[k] = f_12 * gsh_147[k]
                   + f_3 * pc_z[k] * hsh_252[k];

        t_339[k] = pa_x[k] * gsi0_339[k]
                   + f_14 * gsh_255[k]
                   - f_10 * pc_x[k] * gsi1_339[k];

        t_340[k] = f_12 * gsh_170[k]
                   + f_3 * pc_y[k] * hsh_254[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pa_x, pc_x, pc_z, gsi0_341, gsi0_342, gsh_150, \
                         gsh_257, gsh_258, gsi1_341, gsi1_342, \
                         hsh_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pa_x[k] * gsi0_341[k]
                   + f_14 * gsh_257[k]
                   - f_10 * pc_x[k] * gsi1_341[k];

        t_342[k] = pa_x[k] * gsi0_342[k]
                   + f_13 * gsh_258[k]
                   - f_10 * pc_x[k] * gsi1_342[k];

        t_343[k] = f_12 * gsh_150[k]
                   + f_3 * pc_z[k] * hsh_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, pa_x, pc_x, pc_y, gsi0_345, gsi0_346, gsh_173, \
                         gsh_261, gsh_262, gsi1_345, gsi1_346, \
                         hsh_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_12 * gsh_173[k]
                   + f_3 * pc_y[k] * hsh_257[k];

        t_345[k] = pa_x[k] * gsi0_345[k]
                   + f_13 * gsh_261[k]
                   - f_10 * pc_x[k] * gsi1_345[k];

        t_346[k] = pa_x[k] * gsi0_346[k]
                   + f_12 * gsh_262[k]
                   - f_10 * pc_x[k] * gsi1_346[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pa_x, pc_x, pc_y, pc_z, gsi0_348, gsh_153, \
                         gsh_177, gsh_264, gsi1_348, hsh_258, hsh_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_12 * gsh_153[k]
                   + f_3 * pc_z[k] * hsh_258[k];

        t_348[k] = pa_x[k] * gsi0_348[k]
                   + f_12 * gsh_264[k]
                   - f_10 * pc_x[k] * gsi1_348[k];

        t_349[k] = f_12 * gsh_177[k]
                   + f_3 * pc_y[k] * hsh_261[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pa_x, pc_x, gsi0_350, gsh_266, gsh_267, \
                         gsh_268, gsh_269, gsi1_350, hsh_267, hsh_268, \
                         hsh_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pa_x[k] * gsi0_350[k]
                   + f_12 * gsh_266[k]
                   - f_10 * pc_x[k] * gsi1_350[k];

        t_351[k] = f_11 * gsh_267[k]
                   + f_3 * pc_x[k] * hsh_267[k];

        t_352[k] = f_11 * gsh_268[k]
                   + f_3 * pc_x[k] * hsh_268[k];

        t_353[k] = f_11 * gsh_269[k]
                   + f_3 * pc_x[k] * hsh_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pa_x, pc_x, gsi0_357, gsh_270, gsh_271, \
                         gsh_272, gsi1_357, hsh_270, hsh_271, hsh_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_11 * gsh_270[k]
                   + f_3 * pc_x[k] * hsh_270[k];

        t_355[k] = f_11 * gsh_271[k]
                   + f_3 * pc_x[k] * hsh_271[k];

        t_356[k] = f_11 * gsh_272[k]
                   + f_3 * pc_x[k] * hsh_272[k];

        t_357[k] = pa_x[k] * gsi0_357[k]
                   - f_10 * pc_x[k] * gsi1_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pa_x, pc_x, pc_z, gsi0_359, gsi0_360, \
                         gsi0_361, gsh_162, gsi1_359, gsi1_360, gsi1_361, \
                         hsh_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_12 * gsh_162[k]
                   + f_3 * pc_z[k] * hsh_267[k];

        t_359[k] = pa_x[k] * gsi0_359[k]
                   - f_10 * pc_x[k] * gsi1_359[k];

        t_360[k] = pa_x[k] * gsi0_360[k]
                   - f_10 * pc_x[k] * gsi1_360[k];

        t_361[k] = pa_x[k] * gsi0_361[k]
                   - f_10 * pc_x[k] * gsi1_361[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pa_x, pa_y, pc_x, pc_y, gsi0_252, \
                         gsi0_363, gsh_188, gsh_189, gsi1_252, gsi1_363, hsh_272, \
                         hsh_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * gsh_188[k]
                   + f_3 * pc_y[k] * hsh_272[k];

        t_363[k] = pa_x[k] * gsi0_363[k]
                   - f_10 * pc_x[k] * gsi1_363[k];

        t_364[k] = pa_y[k] * gsi0_252[k]
                   - f_10 * pc_y[k] * gsi1_252[k];

        t_365[k] = f_11 * gsh_189[k]
                   + f_3 * pc_y[k] * hsh_273[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, pa_x, pc_x, pc_y, pc_z, gsi0_367, gsh_168, \
                         gsh_191, gsh_276, gsi1_367, hsh_273, hsh_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_13 * gsh_168[k]
                   + f_3 * pc_z[k] * hsh_273[k];

        t_367[k] = pa_x[k] * gsi0_367[k]
                   + f_14 * gsh_276[k]
                   - f_10 * pc_x[k] * gsi1_367[k];

        t_368[k] = f_11 * gsh_191[k]
                   + f_3 * pc_y[k] * hsh_275[k];
    }
}

static auto
compute_prim_hsi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t hsg0, const size_t hsg1,
                                                          const size_t hsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 3.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsi0_257 = buffer.data(gsi0 + 257);
    const auto *gsi0_261 = buffer.data(gsi0 + 261);
    const auto *gsi0_266 = buffer.data(gsi0 + 266);
    const auto *gsi0_280 = buffer.data(gsi0 + 280);
    const auto *gsi0_281 = buffer.data(gsi0 + 281);
    const auto *gsi0_283 = buffer.data(gsi0 + 283);
    const auto *gsi0_286 = buffer.data(gsi0 + 286);
    const auto *gsi0_290 = buffer.data(gsi0 + 290);
    const auto *gsi0_301 = buffer.data(gsi0 + 301);
    const auto *gsi0_303 = buffer.data(gsi0 + 303);
    const auto *gsi0_304 = buffer.data(gsi0 + 304);
    const auto *gsi0_305 = buffer.data(gsi0 + 305);
    const auto *gsi0_370 = buffer.data(gsi0 + 370);
    const auto *gsi0_374 = buffer.data(gsi0 + 374);
    const auto *gsi0_376 = buffer.data(gsi0 + 376);
    const auto *gsi0_385 = buffer.data(gsi0 + 385);
    const auto *gsi0_387 = buffer.data(gsi0 + 387);
    const auto *gsi0_388 = buffer.data(gsi0 + 388);
    const auto *gsi0_389 = buffer.data(gsi0 + 389);
    const auto *gsi0_391 = buffer.data(gsi0 + 391);
    const auto *gsi0_392 = buffer.data(gsi0 + 392);
    const auto *gsi0_397 = buffer.data(gsi0 + 397);
    const auto *gsi0_401 = buffer.data(gsi0 + 401);
    const auto *gsi0_406 = buffer.data(gsi0 + 406);
    const auto *gsi0_413 = buffer.data(gsi0 + 413);
    const auto *gsi0_414 = buffer.data(gsi0 + 414);
    const auto *gsi0_415 = buffer.data(gsi0 + 415);
    const auto *gsi0_416 = buffer.data(gsi0 + 416);
    const auto *gsi0_417 = buffer.data(gsi0 + 417);
    const auto *gsi0_419 = buffer.data(gsi0 + 419);

    const auto *gsh_171 = buffer.data(gsh + 171);
    const auto *gsh_174 = buffer.data(gsh + 174);
    const auto *gsh_183 = buffer.data(gsh + 183);
    const auto *gsh_189 = buffer.data(gsh + 189);
    const auto *gsh_194 = buffer.data(gsh + 194);
    const auto *gsh_198 = buffer.data(gsh + 198);
    const auto *gsh_209 = buffer.data(gsh + 209);
    const auto *gsh_225 = buffer.data(gsh + 225);
    const auto *gsh_226 = buffer.data(gsh + 226);
    const auto *gsh_227 = buffer.data(gsh + 227);
    const auto *gsh_228 = buffer.data(gsh + 228);
    const auto *gsh_230 = buffer.data(gsh + 230);
    const auto *gsh_251 = buffer.data(gsh + 251);
    const auto *gsh_279 = buffer.data(gsh + 279);
    const auto *gsh_283 = buffer.data(gsh + 283);
    const auto *gsh_285 = buffer.data(gsh + 285);
    const auto *gsh_288 = buffer.data(gsh + 288);
    const auto *gsh_289 = buffer.data(gsh + 289);
    const auto *gsh_290 = buffer.data(gsh + 290);
    const auto *gsh_291 = buffer.data(gsh + 291);
    const auto *gsh_292 = buffer.data(gsh + 292);
    const auto *gsh_293 = buffer.data(gsh + 293);
    const auto *gsh_294 = buffer.data(gsh + 294);
    const auto *gsh_299 = buffer.data(gsh + 299);
    const auto *gsh_303 = buffer.data(gsh + 303);
    const auto *gsh_308 = buffer.data(gsh + 308);
    const auto *gsh_309 = buffer.data(gsh + 309);
    const auto *gsh_310 = buffer.data(gsh + 310);
    const auto *gsh_311 = buffer.data(gsh + 311);
    const auto *gsh_312 = buffer.data(gsh + 312);
    const auto *gsh_314 = buffer.data(gsh + 314);

    const auto *gsi1_257 = buffer.data(gsi1 + 257);
    const auto *gsi1_261 = buffer.data(gsi1 + 261);
    const auto *gsi1_266 = buffer.data(gsi1 + 266);
    const auto *gsi1_280 = buffer.data(gsi1 + 280);
    const auto *gsi1_281 = buffer.data(gsi1 + 281);
    const auto *gsi1_283 = buffer.data(gsi1 + 283);
    const auto *gsi1_286 = buffer.data(gsi1 + 286);
    const auto *gsi1_290 = buffer.data(gsi1 + 290);
    const auto *gsi1_301 = buffer.data(gsi1 + 301);
    const auto *gsi1_303 = buffer.data(gsi1 + 303);
    const auto *gsi1_304 = buffer.data(gsi1 + 304);
    const auto *gsi1_305 = buffer.data(gsi1 + 305);
    const auto *gsi1_370 = buffer.data(gsi1 + 370);
    const auto *gsi1_374 = buffer.data(gsi1 + 374);
    const auto *gsi1_376 = buffer.data(gsi1 + 376);
    const auto *gsi1_385 = buffer.data(gsi1 + 385);
    const auto *gsi1_387 = buffer.data(gsi1 + 387);
    const auto *gsi1_388 = buffer.data(gsi1 + 388);
    const auto *gsi1_389 = buffer.data(gsi1 + 389);
    const auto *gsi1_391 = buffer.data(gsi1 + 391);
    const auto *gsi1_392 = buffer.data(gsi1 + 392);
    const auto *gsi1_397 = buffer.data(gsi1 + 397);
    const auto *gsi1_401 = buffer.data(gsi1 + 401);
    const auto *gsi1_406 = buffer.data(gsi1 + 406);
    const auto *gsi1_413 = buffer.data(gsi1 + 413);
    const auto *gsi1_414 = buffer.data(gsi1 + 414);
    const auto *gsi1_415 = buffer.data(gsi1 + 415);
    const auto *gsi1_416 = buffer.data(gsi1 + 416);
    const auto *gsi1_417 = buffer.data(gsi1 + 417);
    const auto *gsi1_419 = buffer.data(gsi1 + 419);

    const auto *hsg0_210 = buffer.data(hsg0 + 210);
    const auto *hsg0_211 = buffer.data(hsg0 + 211);
    const auto *hsg0_212 = buffer.data(hsg0 + 212);
    const auto *hsg0_213 = buffer.data(hsg0 + 213);
    const auto *hsg0_214 = buffer.data(hsg0 + 214);
    const auto *hsg0_215 = buffer.data(hsg0 + 215);
    const auto *hsg0_225 = buffer.data(hsg0 + 225);
    const auto *hsg0_226 = buffer.data(hsg0 + 226);
    const auto *hsg0_228 = buffer.data(hsg0 + 228);
    const auto *hsg0_230 = buffer.data(hsg0 + 230);
    const auto *hsg0_231 = buffer.data(hsg0 + 231);
    const auto *hsg0_233 = buffer.data(hsg0 + 233);
    const auto *hsg0_234 = buffer.data(hsg0 + 234);
    const auto *hsg0_235 = buffer.data(hsg0 + 235);
    const auto *hsg0_236 = buffer.data(hsg0 + 236);
    const auto *hsg0_237 = buffer.data(hsg0 + 237);
    const auto *hsg0_238 = buffer.data(hsg0 + 238);
    const auto *hsg0_239 = buffer.data(hsg0 + 239);
    const auto *hsg0_242 = buffer.data(hsg0 + 242);
    const auto *hsg0_244 = buffer.data(hsg0 + 244);
    const auto *hsg0_245 = buffer.data(hsg0 + 245);
    const auto *hsg0_247 = buffer.data(hsg0 + 247);
    const auto *hsg0_248 = buffer.data(hsg0 + 248);
    const auto *hsg0_249 = buffer.data(hsg0 + 249);
    const auto *hsg0_251 = buffer.data(hsg0 + 251);
    const auto *hsg0_252 = buffer.data(hsg0 + 252);
    const auto *hsg0_253 = buffer.data(hsg0 + 253);
    const auto *hsg0_254 = buffer.data(hsg0 + 254);
    const auto *hsg0_255 = buffer.data(hsg0 + 255);
    const auto *hsg0_256 = buffer.data(hsg0 + 256);
    const auto *hsg0_257 = buffer.data(hsg0 + 257);
    const auto *hsg0_258 = buffer.data(hsg0 + 258);
    const auto *hsg0_259 = buffer.data(hsg0 + 259);
    const auto *hsg0_260 = buffer.data(hsg0 + 260);
    const auto *hsg0_261 = buffer.data(hsg0 + 261);
    const auto *hsg0_262 = buffer.data(hsg0 + 262);
    const auto *hsg0_263 = buffer.data(hsg0 + 263);
    const auto *hsg0_264 = buffer.data(hsg0 + 264);
    const auto *hsg0_265 = buffer.data(hsg0 + 265);
    const auto *hsg0_266 = buffer.data(hsg0 + 266);
    const auto *hsg0_267 = buffer.data(hsg0 + 267);
    const auto *hsg0_268 = buffer.data(hsg0 + 268);
    const auto *hsg0_269 = buffer.data(hsg0 + 269);

    const auto *hsg1_210 = buffer.data(hsg1 + 210);
    const auto *hsg1_211 = buffer.data(hsg1 + 211);
    const auto *hsg1_212 = buffer.data(hsg1 + 212);
    const auto *hsg1_213 = buffer.data(hsg1 + 213);
    const auto *hsg1_214 = buffer.data(hsg1 + 214);
    const auto *hsg1_215 = buffer.data(hsg1 + 215);
    const auto *hsg1_225 = buffer.data(hsg1 + 225);
    const auto *hsg1_226 = buffer.data(hsg1 + 226);
    const auto *hsg1_228 = buffer.data(hsg1 + 228);
    const auto *hsg1_230 = buffer.data(hsg1 + 230);
    const auto *hsg1_231 = buffer.data(hsg1 + 231);
    const auto *hsg1_233 = buffer.data(hsg1 + 233);
    const auto *hsg1_234 = buffer.data(hsg1 + 234);
    const auto *hsg1_235 = buffer.data(hsg1 + 235);
    const auto *hsg1_236 = buffer.data(hsg1 + 236);
    const auto *hsg1_237 = buffer.data(hsg1 + 237);
    const auto *hsg1_238 = buffer.data(hsg1 + 238);
    const auto *hsg1_239 = buffer.data(hsg1 + 239);
    const auto *hsg1_242 = buffer.data(hsg1 + 242);
    const auto *hsg1_244 = buffer.data(hsg1 + 244);
    const auto *hsg1_245 = buffer.data(hsg1 + 245);
    const auto *hsg1_247 = buffer.data(hsg1 + 247);
    const auto *hsg1_248 = buffer.data(hsg1 + 248);
    const auto *hsg1_249 = buffer.data(hsg1 + 249);
    const auto *hsg1_251 = buffer.data(hsg1 + 251);
    const auto *hsg1_252 = buffer.data(hsg1 + 252);
    const auto *hsg1_253 = buffer.data(hsg1 + 253);
    const auto *hsg1_254 = buffer.data(hsg1 + 254);
    const auto *hsg1_255 = buffer.data(hsg1 + 255);
    const auto *hsg1_256 = buffer.data(hsg1 + 256);
    const auto *hsg1_257 = buffer.data(hsg1 + 257);
    const auto *hsg1_258 = buffer.data(hsg1 + 258);
    const auto *hsg1_259 = buffer.data(hsg1 + 259);
    const auto *hsg1_260 = buffer.data(hsg1 + 260);
    const auto *hsg1_261 = buffer.data(hsg1 + 261);
    const auto *hsg1_262 = buffer.data(hsg1 + 262);
    const auto *hsg1_263 = buffer.data(hsg1 + 263);
    const auto *hsg1_264 = buffer.data(hsg1 + 264);
    const auto *hsg1_265 = buffer.data(hsg1 + 265);
    const auto *hsg1_266 = buffer.data(hsg1 + 266);
    const auto *hsg1_267 = buffer.data(hsg1 + 267);
    const auto *hsg1_268 = buffer.data(hsg1 + 268);
    const auto *hsg1_269 = buffer.data(hsg1 + 269);

    const auto *hsh_276 = buffer.data(hsh + 276);
    const auto *hsh_278 = buffer.data(hsh + 278);
    const auto *hsh_279 = buffer.data(hsh + 279);
    const auto *hsh_282 = buffer.data(hsh + 282);
    const auto *hsh_288 = buffer.data(hsh + 288);
    const auto *hsh_289 = buffer.data(hsh + 289);
    const auto *hsh_290 = buffer.data(hsh + 290);
    const auto *hsh_291 = buffer.data(hsh + 291);
    const auto *hsh_292 = buffer.data(hsh + 292);
    const auto *hsh_293 = buffer.data(hsh + 293);
    const auto *hsh_294 = buffer.data(hsh + 294);
    const auto *hsh_295 = buffer.data(hsh + 295);
    const auto *hsh_296 = buffer.data(hsh + 296);
    const auto *hsh_297 = buffer.data(hsh + 297);
    const auto *hsh_298 = buffer.data(hsh + 298);
    const auto *hsh_299 = buffer.data(hsh + 299);
    const auto *hsh_300 = buffer.data(hsh + 300);
    const auto *hsh_301 = buffer.data(hsh + 301);
    const auto *hsh_302 = buffer.data(hsh + 302);
    const auto *hsh_303 = buffer.data(hsh + 303);
    const auto *hsh_308 = buffer.data(hsh + 308);
    const auto *hsh_309 = buffer.data(hsh + 309);
    const auto *hsh_310 = buffer.data(hsh + 310);
    const auto *hsh_311 = buffer.data(hsh + 311);
    const auto *hsh_312 = buffer.data(hsh + 312);
    const auto *hsh_314 = buffer.data(hsh + 314);
    const auto *hsh_315 = buffer.data(hsh + 315);
    const auto *hsh_316 = buffer.data(hsh + 316);
    const auto *hsh_318 = buffer.data(hsh + 318);
    const auto *hsh_320 = buffer.data(hsh + 320);
    const auto *hsh_321 = buffer.data(hsh + 321);
    const auto *hsh_323 = buffer.data(hsh + 323);
    const auto *hsh_324 = buffer.data(hsh + 324);
    const auto *hsh_325 = buffer.data(hsh + 325);
    const auto *hsh_327 = buffer.data(hsh + 327);
    const auto *hsh_328 = buffer.data(hsh + 328);
    const auto *hsh_329 = buffer.data(hsh + 329);
    const auto *hsh_330 = buffer.data(hsh + 330);
    const auto *hsh_331 = buffer.data(hsh + 331);
    const auto *hsh_332 = buffer.data(hsh + 332);
    const auto *hsh_333 = buffer.data(hsh + 333);
    const auto *hsh_334 = buffer.data(hsh + 334);
    const auto *hsh_335 = buffer.data(hsh + 335);
    const auto *hsh_338 = buffer.data(hsh + 338);
    const auto *hsh_340 = buffer.data(hsh + 340);
    const auto *hsh_341 = buffer.data(hsh + 341);
    const auto *hsh_343 = buffer.data(hsh + 343);
    const auto *hsh_344 = buffer.data(hsh + 344);
    const auto *hsh_345 = buffer.data(hsh + 345);
    const auto *hsh_347 = buffer.data(hsh + 347);
    const auto *hsh_348 = buffer.data(hsh + 348);
    const auto *hsh_349 = buffer.data(hsh + 349);
    const auto *hsh_350 = buffer.data(hsh + 350);
    const auto *hsh_351 = buffer.data(hsh + 351);
    const auto *hsh_352 = buffer.data(hsh + 352);
    const auto *hsh_353 = buffer.data(hsh + 353);
    const auto *hsh_354 = buffer.data(hsh + 354);
    const auto *hsh_355 = buffer.data(hsh + 355);
    const auto *hsh_356 = buffer.data(hsh + 356);
    const auto *hsh_357 = buffer.data(hsh + 357);
    const auto *hsh_358 = buffer.data(hsh + 358);
    const auto *hsh_359 = buffer.data(hsh + 359);
    const auto *hsh_360 = buffer.data(hsh + 360);
    const auto *hsh_361 = buffer.data(hsh + 361);
    const auto *hsh_362 = buffer.data(hsh + 362);
    const auto *hsh_363 = buffer.data(hsh + 363);
    const auto *hsh_364 = buffer.data(hsh + 364);
    const auto *hsh_365 = buffer.data(hsh + 365);
    const auto *hsh_366 = buffer.data(hsh + 366);
    const auto *hsh_367 = buffer.data(hsh + 367);
    const auto *hsh_368 = buffer.data(hsh + 368);
    const auto *hsh_369 = buffer.data(hsh + 369);
    const auto *hsh_370 = buffer.data(hsh + 370);
    const auto *hsh_371 = buffer.data(hsh + 371);
    const auto *hsh_372 = buffer.data(hsh + 372);
    const auto *hsh_373 = buffer.data(hsh + 373);
    const auto *hsh_374 = buffer.data(hsh + 374);

#pragma omp simd aligned(t_369, t_370, t_371, pa_x, pa_y, pc_x, pc_y, pc_z, gsi0_257, \
                         gsi0_370, gsh_171, gsh_279, gsi1_257, gsi1_370, \
                         hsh_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * gsi0_257[k]
                   - f_10 * pc_y[k] * gsi1_257[k];

        t_370[k] = pa_x[k] * gsi0_370[k]
                   + f_13 * gsh_279[k]
                   - f_10 * pc_x[k] * gsi1_370[k];

        t_371[k] = f_13 * gsh_171[k]
                   + f_3 * pc_z[k] * hsh_276[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_x, pa_y, pc_x, pc_y, gsi0_261, gsi0_374, \
                         gsh_194, gsh_283, gsi1_261, gsi1_374, \
                         hsh_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_11 * gsh_194[k]
                   + f_3 * pc_y[k] * hsh_278[k];

        t_373[k] = pa_y[k] * gsi0_261[k]
                   - f_10 * pc_y[k] * gsi1_261[k];

        t_374[k] = pa_x[k] * gsi0_374[k]
                   + f_12 * gsh_283[k]
                   - f_10 * pc_x[k] * gsi1_374[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pc_x, pc_y, pc_z, gsi0_376, gsh_174, \
                         gsh_198, gsh_285, gsi1_376, hsh_279, hsh_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_13 * gsh_174[k]
                   + f_3 * pc_z[k] * hsh_279[k];

        t_376[k] = pa_x[k] * gsi0_376[k]
                   + f_12 * gsh_285[k]
                   - f_10 * pc_x[k] * gsi1_376[k];

        t_377[k] = f_11 * gsh_198[k]
                   + f_3 * pc_y[k] * hsh_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pa_y, pc_x, pc_y, gsi0_266, gsh_288, \
                         gsh_289, gsh_290, gsi1_266, hsh_288, hsh_289, \
                         hsh_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_y[k] * gsi0_266[k]
                   - f_10 * pc_y[k] * gsi1_266[k];

        t_379[k] = f_11 * gsh_288[k]
                   + f_3 * pc_x[k] * hsh_288[k];

        t_380[k] = f_11 * gsh_289[k]
                   + f_3 * pc_x[k] * hsh_289[k];

        t_381[k] = f_11 * gsh_290[k]
                   + f_3 * pc_x[k] * hsh_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_x, pc_x, gsi0_385, gsh_291, gsh_292, \
                         gsh_293, gsi1_385, hsh_291, hsh_292, hsh_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_11 * gsh_291[k]
                   + f_3 * pc_x[k] * hsh_291[k];

        t_383[k] = f_11 * gsh_292[k]
                   + f_3 * pc_x[k] * hsh_292[k];

        t_384[k] = f_11 * gsh_293[k]
                   + f_3 * pc_x[k] * hsh_293[k];

        t_385[k] = pa_x[k] * gsi0_385[k]
                   - f_10 * pc_x[k] * gsi1_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_x, pc_x, pc_z, gsi0_387, gsi0_388, \
                         gsi0_389, gsh_183, gsi1_387, gsi1_388, gsi1_389, \
                         hsh_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_13 * gsh_183[k]
                   + f_3 * pc_z[k] * hsh_288[k];

        t_387[k] = pa_x[k] * gsi0_387[k]
                   - f_10 * pc_x[k] * gsi1_387[k];

        t_388[k] = pa_x[k] * gsi0_388[k]
                   - f_10 * pc_x[k] * gsi1_388[k];

        t_389[k] = pa_x[k] * gsi0_389[k]
                   - f_10 * pc_x[k] * gsi1_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_x, pc_x, pc_y, gsi0_391, gsi0_392, \
                         gsh_209, gsh_294, gsi1_391, gsi1_392, hsh_293, \
                         hsh_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_11 * gsh_209[k]
                   + f_3 * pc_y[k] * hsh_293[k];

        t_391[k] = pa_x[k] * gsi0_391[k]
                   - f_10 * pc_x[k] * gsi1_391[k];

        t_392[k] = pa_x[k] * gsi0_392[k]
                   + f_17 * gsh_294[k]
                   - f_10 * pc_x[k] * gsi1_392[k];

        t_393[k] = f_3 * pc_y[k] * hsh_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, gsh_189, hsg0_210, hsg1_210, \
                         hsh_294, hsh_295, hsh_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_14 * gsh_189[k]
                   + f_3 * pc_z[k] * hsh_294[k];

        t_395[k] = f_4 * hsg0_210[k]
                   - f_5 * hsg1_210[k]
                   + f_3 * pc_y[k] * hsh_295[k];

        t_396[k] = f_3 * pc_y[k] * hsh_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pa_x, pc_x, pc_y, gsi0_397, gsh_299, gsi1_397, \
                         hsg0_211, hsg0_212, hsg1_211, hsg1_212, hsh_297, \
                         hsh_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pa_x[k] * gsi0_397[k]
                   + f_14 * gsh_299[k]
                   - f_10 * pc_x[k] * gsi1_397[k];

        t_398[k] = f_6 * hsg0_211[k]
                   - f_7 * hsg1_211[k]
                   + f_3 * pc_y[k] * hsh_297[k];

        t_399[k] = f_4 * hsg0_212[k]
                   - f_5 * hsg1_212[k]
                   + f_3 * pc_y[k] * hsh_298[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pa_x, pc_x, pc_y, gsi0_401, gsh_303, gsi1_401, \
                         hsg0_213, hsg1_213, hsh_299, hsh_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_3 * pc_y[k] * hsh_299[k];

        t_401[k] = pa_x[k] * gsi0_401[k]
                   + f_13 * gsh_303[k]
                   - f_10 * pc_x[k] * gsi1_401[k];

        t_402[k] = f_8 * hsg0_213[k]
                   - f_9 * hsg1_213[k]
                   + f_3 * pc_y[k] * hsh_300[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pc_y, hsg0_214, hsg0_215, hsg1_214, hsg1_215, \
                         hsh_301, hsh_302, hsh_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_6 * hsg0_214[k]
                   - f_7 * hsg1_214[k]
                   + f_3 * pc_y[k] * hsh_301[k];

        t_404[k] = f_4 * hsg0_215[k]
                   - f_5 * hsg1_215[k]
                   + f_3 * pc_y[k] * hsh_302[k];

        t_405[k] = f_3 * pc_y[k] * hsh_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_x, pc_x, gsi0_406, gsh_308, gsh_309, \
                         gsh_310, gsh_311, gsi1_406, hsh_309, hsh_310, \
                         hsh_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_x[k] * gsi0_406[k]
                   + f_12 * gsh_308[k]
                   - f_10 * pc_x[k] * gsi1_406[k];

        t_407[k] = f_11 * gsh_309[k]
                   + f_3 * pc_x[k] * hsh_309[k];

        t_408[k] = f_11 * gsh_310[k]
                   + f_3 * pc_x[k] * hsh_310[k];

        t_409[k] = f_11 * gsh_311[k]
                   + f_3 * pc_x[k] * hsh_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pc_x, pc_y, gsi0_413, gsh_312, \
                         gsh_314, gsi1_413, hsh_308, hsh_312, hsh_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_11 * gsh_312[k]
                   + f_3 * pc_x[k] * hsh_312[k];

        t_411[k] = f_3 * pc_y[k] * hsh_308[k];

        t_412[k] = f_11 * gsh_314[k]
                   + f_3 * pc_x[k] * hsh_314[k];

        t_413[k] = pa_x[k] * gsi0_413[k]
                   - f_10 * pc_x[k] * gsi1_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pa_x, pc_x, gsi0_414, gsi0_415, gsi0_416, \
                         gsi0_417, gsi1_414, gsi1_415, gsi1_416, \
                         gsi1_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pa_x[k] * gsi0_414[k]
                   - f_10 * pc_x[k] * gsi1_414[k];

        t_415[k] = pa_x[k] * gsi0_415[k]
                   - f_10 * pc_x[k] * gsi1_415[k];

        t_416[k] = pa_x[k] * gsi0_416[k]
                   - f_10 * pc_x[k] * gsi1_416[k];

        t_417[k] = pa_x[k] * gsi0_417[k]
                   - f_10 * pc_x[k] * gsi1_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_x, pc_x, pc_y, gsi0_419, gsi1_419, \
                         hsg0_225, hsg0_226, hsg1_225, hsg1_226, hsh_314, hsh_315, \
                         hsh_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * hsh_314[k];

        t_419[k] = pa_x[k] * gsi0_419[k]
                   - f_10 * pc_x[k] * gsi1_419[k];

        t_420[k] = f_1 * hsg0_225[k]
                   - f_2 * hsg1_225[k]
                   + f_3 * pc_x[k] * hsh_315[k];

        t_421[k] = f_15 * hsg0_226[k]
                   - f_16 * hsg1_226[k]
                   + f_3 * pc_x[k] * hsh_316[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, t_425, pc_x, pc_z, hsg0_228, hsg0_230, hsg1_228, \
                         hsg1_230, hsh_315, hsh_316, hsh_318, hsh_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_3 * pc_z[k] * hsh_315[k];

        t_423[k] = f_8 * hsg0_228[k]
                   - f_9 * hsg1_228[k]
                   + f_3 * pc_x[k] * hsh_318[k];

        t_424[k] = f_3 * pc_z[k] * hsh_316[k];

        t_425[k] = f_8 * hsg0_230[k]
                   - f_9 * hsg1_230[k]
                   + f_3 * pc_x[k] * hsh_320[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pc_x, pc_z, hsg0_231, hsg0_233, hsg0_234, \
                         hsg1_231, hsg1_233, hsg1_234, hsh_318, hsh_321, hsh_323, \
                         hsh_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_6 * hsg0_231[k]
                   - f_7 * hsg1_231[k]
                   + f_3 * pc_x[k] * hsh_321[k];

        t_427[k] = f_3 * pc_z[k] * hsh_318[k];

        t_428[k] = f_6 * hsg0_233[k]
                   - f_7 * hsg1_233[k]
                   + f_3 * pc_x[k] * hsh_323[k];

        t_429[k] = f_6 * hsg0_234[k]
                   - f_7 * hsg1_234[k]
                   + f_3 * pc_x[k] * hsh_324[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pc_x, pc_z, hsg0_235, hsg0_237, hsg0_238, \
                         hsg1_235, hsg1_237, hsg1_238, hsh_321, hsh_325, hsh_327, \
                         hsh_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_4 * hsg0_235[k]
                   - f_5 * hsg1_235[k]
                   + f_3 * pc_x[k] * hsh_325[k];

        t_431[k] = f_3 * pc_z[k] * hsh_321[k];

        t_432[k] = f_4 * hsg0_237[k]
                   - f_5 * hsg1_237[k]
                   + f_3 * pc_x[k] * hsh_327[k];

        t_433[k] = f_4 * hsg0_238[k]
                   - f_5 * hsg1_238[k]
                   + f_3 * pc_x[k] * hsh_328[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, t_438, t_439, pc_x, hsg0_239, hsg1_239, \
                         hsh_329, hsh_330, hsh_331, hsh_332, hsh_333, \
                         hsh_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_4 * hsg0_239[k]
                   - f_5 * hsg1_239[k]
                   + f_3 * pc_x[k] * hsh_329[k];

        t_435[k] = f_3 * pc_x[k] * hsh_330[k];

        t_436[k] = f_3 * pc_x[k] * hsh_331[k];

        t_437[k] = f_3 * pc_x[k] * hsh_332[k];

        t_438[k] = f_3 * pc_x[k] * hsh_333[k];

        t_439[k] = f_3 * pc_x[k] * hsh_334[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, pc_x, pc_y, pc_z, gsh_225, hsg0_235, \
                         hsg1_235, hsh_330, hsh_331, hsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_3 * pc_x[k] * hsh_335[k];

        t_441[k] = f_0 * gsh_225[k]
                   + f_1 * hsg0_235[k]
                   - f_2 * hsg1_235[k]
                   + f_3 * pc_y[k] * hsh_330[k];

        t_442[k] = f_3 * pc_z[k] * hsh_330[k];

        t_443[k] = f_4 * hsg0_235[k]
                   - f_5 * hsg1_235[k]
                   + f_3 * pc_z[k] * hsh_331[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, gsh_230, hsg0_236, hsg0_237, \
                         hsg0_239, hsg1_236, hsg1_237, hsg1_239, hsh_332, hsh_333, \
                         hsh_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_6 * hsg0_236[k]
                   - f_7 * hsg1_236[k]
                   + f_3 * pc_z[k] * hsh_332[k];

        t_445[k] = f_8 * hsg0_237[k]
                   - f_9 * hsg1_237[k]
                   + f_3 * pc_z[k] * hsh_333[k];

        t_446[k] = f_0 * gsh_230[k]
                   + f_3 * pc_y[k] * hsh_335[k];

        t_447[k] = f_1 * hsg0_239[k]
                   - f_2 * hsg1_239[k]
                   + f_3 * pc_z[k] * hsh_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, t_451, pa_z, pc_x, pc_z, gsi0_280, gsi0_281, \
                         gsi0_283, gsi1_280, gsi1_281, gsi1_283, hsg0_242, hsg1_242, \
                         hsh_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pa_z[k] * gsi0_280[k]
                   - f_10 * pc_z[k] * gsi1_280[k];

        t_449[k] = pa_z[k] * gsi0_281[k]
                   - f_10 * pc_z[k] * gsi1_281[k];

        t_450[k] = f_15 * hsg0_242[k]
                   - f_16 * hsg1_242[k]
                   + f_3 * pc_x[k] * hsh_338[k];

        t_451[k] = pa_z[k] * gsi0_283[k]
                   - f_10 * pc_z[k] * gsi1_283[k];
    }

#pragma omp simd aligned(t_452, t_453, t_454, pa_z, pc_x, pc_z, gsi0_286, gsi1_286, hsg0_244, \
                         hsg0_245, hsg1_244, hsg1_245, hsh_340, \
                         hsh_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_452[k] = f_8 * hsg0_244[k]
                   - f_9 * hsg1_244[k]
                   + f_3 * pc_x[k] * hsh_340[k];

        t_453[k] = f_8 * hsg0_245[k]
                   - f_9 * hsg1_245[k]
                   + f_3 * pc_x[k] * hsh_341[k];

        t_454[k] = pa_z[k] * gsi0_286[k]
                   - f_10 * pc_z[k] * gsi1_286[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, pc_x, hsg0_247, hsg0_248, hsg0_249, hsg1_247, \
                         hsg1_248, hsg1_249, hsh_343, hsh_344, \
                         hsh_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_6 * hsg0_247[k]
                   - f_7 * hsg1_247[k]
                   + f_3 * pc_x[k] * hsh_343[k];

        t_456[k] = f_6 * hsg0_248[k]
                   - f_7 * hsg1_248[k]
                   + f_3 * pc_x[k] * hsh_344[k];

        t_457[k] = f_6 * hsg0_249[k]
                   - f_7 * hsg1_249[k]
                   + f_3 * pc_x[k] * hsh_345[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, pa_z, pc_x, pc_z, gsi0_290, gsi1_290, hsg0_251, \
                         hsg0_252, hsg1_251, hsg1_252, hsh_347, \
                         hsh_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = pa_z[k] * gsi0_290[k]
                   - f_10 * pc_z[k] * gsi1_290[k];

        t_459[k] = f_4 * hsg0_251[k]
                   - f_5 * hsg1_251[k]
                   + f_3 * pc_x[k] * hsh_347[k];

        t_460[k] = f_4 * hsg0_252[k]
                   - f_5 * hsg1_252[k]
                   + f_3 * pc_x[k] * hsh_348[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, t_465, pc_x, hsg0_253, hsg0_254, \
                         hsg1_253, hsg1_254, hsh_349, hsh_350, hsh_351, hsh_352, \
                         hsh_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_4 * hsg0_253[k]
                   - f_5 * hsg1_253[k]
                   + f_3 * pc_x[k] * hsh_349[k];

        t_462[k] = f_4 * hsg0_254[k]
                   - f_5 * hsg1_254[k]
                   + f_3 * pc_x[k] * hsh_350[k];

        t_463[k] = f_3 * pc_x[k] * hsh_351[k];

        t_464[k] = f_3 * pc_x[k] * hsh_352[k];

        t_465[k] = f_3 * pc_x[k] * hsh_353[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, t_470, pa_z, pc_x, pc_z, gsi0_301, \
                         gsh_225, gsi1_301, hsh_351, hsh_354, hsh_355, \
                         hsh_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = f_3 * pc_x[k] * hsh_354[k];

        t_467[k] = f_3 * pc_x[k] * hsh_355[k];

        t_468[k] = f_3 * pc_x[k] * hsh_356[k];

        t_469[k] = pa_z[k] * gsi0_301[k]
                   - f_10 * pc_z[k] * gsi1_301[k];

        t_470[k] = f_11 * gsh_225[k]
                   + f_3 * pc_z[k] * hsh_351[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pa_z, pc_z, gsi0_303, gsi0_304, gsi0_305, \
                         gsh_226, gsh_227, gsh_228, gsi1_303, gsi1_304, \
                         gsi1_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pa_z[k] * gsi0_303[k]
                   + f_12 * gsh_226[k]
                   - f_10 * pc_z[k] * gsi1_303[k];

        t_472[k] = pa_z[k] * gsi0_304[k]
                   + f_13 * gsh_227[k]
                   - f_10 * pc_z[k] * gsi1_304[k];

        t_473[k] = pa_z[k] * gsi0_305[k]
                   + f_14 * gsh_228[k]
                   - f_10 * pc_z[k] * gsi1_305[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, gsh_230, gsh_251, hsg0_254, \
                         hsg0_255, hsg1_254, hsg1_255, hsh_356, \
                         hsh_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_14 * gsh_251[k]
                   + f_3 * pc_y[k] * hsh_356[k];

        t_475[k] = f_11 * gsh_230[k]
                   + f_1 * hsg0_254[k]
                   - f_2 * hsg1_254[k]
                   + f_3 * pc_z[k] * hsh_356[k];

        t_476[k] = f_1 * hsg0_255[k]
                   - f_2 * hsg1_255[k]
                   + f_3 * pc_x[k] * hsh_357[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pc_x, hsg0_256, hsg0_257, hsg0_258, hsg1_256, \
                         hsg1_257, hsg1_258, hsh_358, hsh_359, \
                         hsh_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_15 * hsg0_256[k]
                   - f_16 * hsg1_256[k]
                   + f_3 * pc_x[k] * hsh_358[k];

        t_478[k] = f_15 * hsg0_257[k]
                   - f_16 * hsg1_257[k]
                   + f_3 * pc_x[k] * hsh_359[k];

        t_479[k] = f_8 * hsg0_258[k]
                   - f_9 * hsg1_258[k]
                   + f_3 * pc_x[k] * hsh_360[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, pc_x, hsg0_259, hsg0_260, hsg0_261, hsg1_259, \
                         hsg1_260, hsg1_261, hsh_361, hsh_362, \
                         hsh_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_8 * hsg0_259[k]
                   - f_9 * hsg1_259[k]
                   + f_3 * pc_x[k] * hsh_361[k];

        t_481[k] = f_8 * hsg0_260[k]
                   - f_9 * hsg1_260[k]
                   + f_3 * pc_x[k] * hsh_362[k];

        t_482[k] = f_6 * hsg0_261[k]
                   - f_7 * hsg1_261[k]
                   + f_3 * pc_x[k] * hsh_363[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, pc_x, hsg0_262, hsg0_263, hsg0_264, hsg1_262, \
                         hsg1_263, hsg1_264, hsh_364, hsh_365, \
                         hsh_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = f_6 * hsg0_262[k]
                   - f_7 * hsg1_262[k]
                   + f_3 * pc_x[k] * hsh_364[k];

        t_484[k] = f_6 * hsg0_263[k]
                   - f_7 * hsg1_263[k]
                   + f_3 * pc_x[k] * hsh_365[k];

        t_485[k] = f_6 * hsg0_264[k]
                   - f_7 * hsg1_264[k]
                   + f_3 * pc_x[k] * hsh_366[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_x, hsg0_265, hsg0_266, hsg0_267, hsg1_265, \
                         hsg1_266, hsg1_267, hsh_367, hsh_368, \
                         hsh_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_4 * hsg0_265[k]
                   - f_5 * hsg1_265[k]
                   + f_3 * pc_x[k] * hsh_367[k];

        t_487[k] = f_4 * hsg0_266[k]
                   - f_5 * hsg1_266[k]
                   + f_3 * pc_x[k] * hsh_368[k];

        t_488[k] = f_4 * hsg0_267[k]
                   - f_5 * hsg1_267[k]
                   + f_3 * pc_x[k] * hsh_369[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, t_492, t_493, pc_x, hsg0_268, hsg0_269, \
                         hsg1_268, hsg1_269, hsh_370, hsh_371, hsh_372, hsh_373, \
                         hsh_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_4 * hsg0_268[k]
                   - f_5 * hsg1_268[k]
                   + f_3 * pc_x[k] * hsh_370[k];

        t_490[k] = f_4 * hsg0_269[k]
                   - f_5 * hsg1_269[k]
                   + f_3 * pc_x[k] * hsh_371[k];

        t_491[k] = f_3 * pc_x[k] * hsh_372[k];

        t_492[k] = f_3 * pc_x[k] * hsh_373[k];

        t_493[k] = f_3 * pc_x[k] * hsh_374[k];
    }
}

static auto
compute_prim_hsi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t gsi0,
                                                          const size_t gsh, const size_t gsi1,
                                                          const size_t hsg0, const size_t hsg1,
                                                          const size_t hsh, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 3.0 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *gsi0_392 = buffer.data(gsi0 + 392);
    const auto *gsi0_394 = buffer.data(gsi0 + 394);
    const auto *gsi0_397 = buffer.data(gsi0 + 397);
    const auto *gsi0_401 = buffer.data(gsi0 + 401);
    const auto *gsi0_406 = buffer.data(gsi0 + 406);
    const auto *gsi0_413 = buffer.data(gsi0 + 413);
    const auto *gsi0_415 = buffer.data(gsi0 + 415);
    const auto *gsi0_416 = buffer.data(gsi0 + 416);
    const auto *gsi0_417 = buffer.data(gsi0 + 417);
    const auto *gsi0_419 = buffer.data(gsi0 + 419);

    const auto *gsh_246 = buffer.data(gsh + 246);
    const auto *gsh_251 = buffer.data(gsh + 251);
    const auto *gsh_267 = buffer.data(gsh + 267);
    const auto *gsh_269 = buffer.data(gsh + 269);
    const auto *gsh_270 = buffer.data(gsh + 270);
    const auto *gsh_271 = buffer.data(gsh + 271);
    const auto *gsh_272 = buffer.data(gsh + 272);
    const auto *gsh_288 = buffer.data(gsh + 288);
    const auto *gsh_290 = buffer.data(gsh + 290);
    const auto *gsh_291 = buffer.data(gsh + 291);
    const auto *gsh_292 = buffer.data(gsh + 292);
    const auto *gsh_293 = buffer.data(gsh + 293);
    const auto *gsh_309 = buffer.data(gsh + 309);
    const auto *gsh_311 = buffer.data(gsh + 311);
    const auto *gsh_312 = buffer.data(gsh + 312);
    const auto *gsh_313 = buffer.data(gsh + 313);
    const auto *gsh_314 = buffer.data(gsh + 314);

    const auto *gsi1_392 = buffer.data(gsi1 + 392);
    const auto *gsi1_394 = buffer.data(gsi1 + 394);
    const auto *gsi1_397 = buffer.data(gsi1 + 397);
    const auto *gsi1_401 = buffer.data(gsi1 + 401);
    const auto *gsi1_406 = buffer.data(gsi1 + 406);
    const auto *gsi1_413 = buffer.data(gsi1 + 413);
    const auto *gsi1_415 = buffer.data(gsi1 + 415);
    const auto *gsi1_416 = buffer.data(gsi1 + 416);
    const auto *gsi1_417 = buffer.data(gsi1 + 417);
    const auto *gsi1_419 = buffer.data(gsi1 + 419);

    const auto *hsg0_265 = buffer.data(hsg0 + 265);
    const auto *hsg0_267 = buffer.data(hsg0 + 267);
    const auto *hsg0_268 = buffer.data(hsg0 + 268);
    const auto *hsg0_269 = buffer.data(hsg0 + 269);
    const auto *hsg0_270 = buffer.data(hsg0 + 270);
    const auto *hsg0_271 = buffer.data(hsg0 + 271);
    const auto *hsg0_272 = buffer.data(hsg0 + 272);
    const auto *hsg0_273 = buffer.data(hsg0 + 273);
    const auto *hsg0_274 = buffer.data(hsg0 + 274);
    const auto *hsg0_275 = buffer.data(hsg0 + 275);
    const auto *hsg0_276 = buffer.data(hsg0 + 276);
    const auto *hsg0_277 = buffer.data(hsg0 + 277);
    const auto *hsg0_278 = buffer.data(hsg0 + 278);
    const auto *hsg0_279 = buffer.data(hsg0 + 279);
    const auto *hsg0_280 = buffer.data(hsg0 + 280);
    const auto *hsg0_281 = buffer.data(hsg0 + 281);
    const auto *hsg0_282 = buffer.data(hsg0 + 282);
    const auto *hsg0_283 = buffer.data(hsg0 + 283);
    const auto *hsg0_284 = buffer.data(hsg0 + 284);
    const auto *hsg0_286 = buffer.data(hsg0 + 286);
    const auto *hsg0_288 = buffer.data(hsg0 + 288);
    const auto *hsg0_289 = buffer.data(hsg0 + 289);
    const auto *hsg0_291 = buffer.data(hsg0 + 291);
    const auto *hsg0_292 = buffer.data(hsg0 + 292);
    const auto *hsg0_293 = buffer.data(hsg0 + 293);
    const auto *hsg0_295 = buffer.data(hsg0 + 295);
    const auto *hsg0_296 = buffer.data(hsg0 + 296);
    const auto *hsg0_297 = buffer.data(hsg0 + 297);
    const auto *hsg0_298 = buffer.data(hsg0 + 298);
    const auto *hsg0_300 = buffer.data(hsg0 + 300);
    const auto *hsg0_302 = buffer.data(hsg0 + 302);
    const auto *hsg0_303 = buffer.data(hsg0 + 303);
    const auto *hsg0_305 = buffer.data(hsg0 + 305);
    const auto *hsg0_306 = buffer.data(hsg0 + 306);
    const auto *hsg0_307 = buffer.data(hsg0 + 307);
    const auto *hsg0_309 = buffer.data(hsg0 + 309);
    const auto *hsg0_310 = buffer.data(hsg0 + 310);
    const auto *hsg0_311 = buffer.data(hsg0 + 311);
    const auto *hsg0_312 = buffer.data(hsg0 + 312);
    const auto *hsg0_313 = buffer.data(hsg0 + 313);
    const auto *hsg0_314 = buffer.data(hsg0 + 314);

    const auto *hsg1_265 = buffer.data(hsg1 + 265);
    const auto *hsg1_267 = buffer.data(hsg1 + 267);
    const auto *hsg1_268 = buffer.data(hsg1 + 268);
    const auto *hsg1_269 = buffer.data(hsg1 + 269);
    const auto *hsg1_270 = buffer.data(hsg1 + 270);
    const auto *hsg1_271 = buffer.data(hsg1 + 271);
    const auto *hsg1_272 = buffer.data(hsg1 + 272);
    const auto *hsg1_273 = buffer.data(hsg1 + 273);
    const auto *hsg1_274 = buffer.data(hsg1 + 274);
    const auto *hsg1_275 = buffer.data(hsg1 + 275);
    const auto *hsg1_276 = buffer.data(hsg1 + 276);
    const auto *hsg1_277 = buffer.data(hsg1 + 277);
    const auto *hsg1_278 = buffer.data(hsg1 + 278);
    const auto *hsg1_279 = buffer.data(hsg1 + 279);
    const auto *hsg1_280 = buffer.data(hsg1 + 280);
    const auto *hsg1_281 = buffer.data(hsg1 + 281);
    const auto *hsg1_282 = buffer.data(hsg1 + 282);
    const auto *hsg1_283 = buffer.data(hsg1 + 283);
    const auto *hsg1_284 = buffer.data(hsg1 + 284);
    const auto *hsg1_286 = buffer.data(hsg1 + 286);
    const auto *hsg1_288 = buffer.data(hsg1 + 288);
    const auto *hsg1_289 = buffer.data(hsg1 + 289);
    const auto *hsg1_291 = buffer.data(hsg1 + 291);
    const auto *hsg1_292 = buffer.data(hsg1 + 292);
    const auto *hsg1_293 = buffer.data(hsg1 + 293);
    const auto *hsg1_295 = buffer.data(hsg1 + 295);
    const auto *hsg1_296 = buffer.data(hsg1 + 296);
    const auto *hsg1_297 = buffer.data(hsg1 + 297);
    const auto *hsg1_298 = buffer.data(hsg1 + 298);
    const auto *hsg1_300 = buffer.data(hsg1 + 300);
    const auto *hsg1_302 = buffer.data(hsg1 + 302);
    const auto *hsg1_303 = buffer.data(hsg1 + 303);
    const auto *hsg1_305 = buffer.data(hsg1 + 305);
    const auto *hsg1_306 = buffer.data(hsg1 + 306);
    const auto *hsg1_307 = buffer.data(hsg1 + 307);
    const auto *hsg1_309 = buffer.data(hsg1 + 309);
    const auto *hsg1_310 = buffer.data(hsg1 + 310);
    const auto *hsg1_311 = buffer.data(hsg1 + 311);
    const auto *hsg1_312 = buffer.data(hsg1 + 312);
    const auto *hsg1_313 = buffer.data(hsg1 + 313);
    const auto *hsg1_314 = buffer.data(hsg1 + 314);

    const auto *hsh_372 = buffer.data(hsh + 372);
    const auto *hsh_374 = buffer.data(hsh + 374);
    const auto *hsh_375 = buffer.data(hsh + 375);
    const auto *hsh_376 = buffer.data(hsh + 376);
    const auto *hsh_377 = buffer.data(hsh + 377);
    const auto *hsh_378 = buffer.data(hsh + 378);
    const auto *hsh_379 = buffer.data(hsh + 379);
    const auto *hsh_380 = buffer.data(hsh + 380);
    const auto *hsh_381 = buffer.data(hsh + 381);
    const auto *hsh_382 = buffer.data(hsh + 382);
    const auto *hsh_383 = buffer.data(hsh + 383);
    const auto *hsh_384 = buffer.data(hsh + 384);
    const auto *hsh_385 = buffer.data(hsh + 385);
    const auto *hsh_386 = buffer.data(hsh + 386);
    const auto *hsh_387 = buffer.data(hsh + 387);
    const auto *hsh_388 = buffer.data(hsh + 388);
    const auto *hsh_389 = buffer.data(hsh + 389);
    const auto *hsh_390 = buffer.data(hsh + 390);
    const auto *hsh_391 = buffer.data(hsh + 391);
    const auto *hsh_392 = buffer.data(hsh + 392);
    const auto *hsh_393 = buffer.data(hsh + 393);
    const auto *hsh_394 = buffer.data(hsh + 394);
    const auto *hsh_395 = buffer.data(hsh + 395);
    const auto *hsh_396 = buffer.data(hsh + 396);
    const auto *hsh_397 = buffer.data(hsh + 397);
    const auto *hsh_398 = buffer.data(hsh + 398);
    const auto *hsh_400 = buffer.data(hsh + 400);
    const auto *hsh_402 = buffer.data(hsh + 402);
    const auto *hsh_403 = buffer.data(hsh + 403);
    const auto *hsh_405 = buffer.data(hsh + 405);
    const auto *hsh_406 = buffer.data(hsh + 406);
    const auto *hsh_407 = buffer.data(hsh + 407);
    const auto *hsh_409 = buffer.data(hsh + 409);
    const auto *hsh_410 = buffer.data(hsh + 410);
    const auto *hsh_411 = buffer.data(hsh + 411);
    const auto *hsh_412 = buffer.data(hsh + 412);
    const auto *hsh_414 = buffer.data(hsh + 414);
    const auto *hsh_415 = buffer.data(hsh + 415);
    const auto *hsh_416 = buffer.data(hsh + 416);
    const auto *hsh_417 = buffer.data(hsh + 417);
    const auto *hsh_418 = buffer.data(hsh + 418);
    const auto *hsh_419 = buffer.data(hsh + 419);
    const auto *hsh_420 = buffer.data(hsh + 420);
    const auto *hsh_422 = buffer.data(hsh + 422);
    const auto *hsh_423 = buffer.data(hsh + 423);
    const auto *hsh_425 = buffer.data(hsh + 425);
    const auto *hsh_426 = buffer.data(hsh + 426);
    const auto *hsh_427 = buffer.data(hsh + 427);
    const auto *hsh_429 = buffer.data(hsh + 429);
    const auto *hsh_430 = buffer.data(hsh + 430);
    const auto *hsh_431 = buffer.data(hsh + 431);
    const auto *hsh_432 = buffer.data(hsh + 432);
    const auto *hsh_434 = buffer.data(hsh + 434);
    const auto *hsh_435 = buffer.data(hsh + 435);
    const auto *hsh_436 = buffer.data(hsh + 436);
    const auto *hsh_437 = buffer.data(hsh + 437);
    const auto *hsh_438 = buffer.data(hsh + 438);
    const auto *hsh_439 = buffer.data(hsh + 439);
    const auto *hsh_440 = buffer.data(hsh + 440);

#pragma omp simd aligned(t_494, t_495, t_496, t_497, t_498, pc_x, pc_y, pc_z, gsh_246, \
                         gsh_267, hsg0_265, hsg1_265, hsh_372, hsh_375, hsh_376, \
                         hsh_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_494[k] = f_3 * pc_x[k] * hsh_375[k];

        t_495[k] = f_3 * pc_x[k] * hsh_376[k];

        t_496[k] = f_3 * pc_x[k] * hsh_377[k];

        t_497[k] = f_13 * gsh_267[k]
                   + f_1 * hsg0_265[k]
                   - f_2 * hsg1_265[k]
                   + f_3 * pc_y[k] * hsh_372[k];

        t_498[k] = f_12 * gsh_246[k]
                   + f_3 * pc_z[k] * hsh_372[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, pc_y, gsh_269, gsh_270, gsh_271, hsg0_267, \
                         hsg0_268, hsg0_269, hsg1_267, hsg1_268, hsg1_269, hsh_374, hsh_375, \
                         hsh_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_13 * gsh_269[k]
                   + f_8 * hsg0_267[k]
                   - f_9 * hsg1_267[k]
                   + f_3 * pc_y[k] * hsh_374[k];

        t_500[k] = f_13 * gsh_270[k]
                   + f_6 * hsg0_268[k]
                   - f_7 * hsg1_268[k]
                   + f_3 * pc_y[k] * hsh_375[k];

        t_501[k] = f_13 * gsh_271[k]
                   + f_4 * hsg0_269[k]
                   - f_5 * hsg1_269[k]
                   + f_3 * pc_y[k] * hsh_376[k];
    }

#pragma omp simd aligned(t_502, t_503, t_504, pc_x, pc_y, pc_z, gsh_251, gsh_272, hsg0_269, \
                         hsg0_270, hsg1_269, hsg1_270, hsh_377, \
                         hsh_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_502[k] = f_13 * gsh_272[k]
                   + f_3 * pc_y[k] * hsh_377[k];

        t_503[k] = f_12 * gsh_251[k]
                   + f_1 * hsg0_269[k]
                   - f_2 * hsg1_269[k]
                   + f_3 * pc_z[k] * hsh_377[k];

        t_504[k] = f_1 * hsg0_270[k]
                   - f_2 * hsg1_270[k]
                   + f_3 * pc_x[k] * hsh_378[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pc_x, hsg0_271, hsg0_272, hsg0_273, hsg1_271, \
                         hsg1_272, hsg1_273, hsh_379, hsh_380, \
                         hsh_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_15 * hsg0_271[k]
                   - f_16 * hsg1_271[k]
                   + f_3 * pc_x[k] * hsh_379[k];

        t_506[k] = f_15 * hsg0_272[k]
                   - f_16 * hsg1_272[k]
                   + f_3 * pc_x[k] * hsh_380[k];

        t_507[k] = f_8 * hsg0_273[k]
                   - f_9 * hsg1_273[k]
                   + f_3 * pc_x[k] * hsh_381[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, hsg0_274, hsg0_275, hsg0_276, hsg1_274, \
                         hsg1_275, hsg1_276, hsh_382, hsh_383, \
                         hsh_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * hsg0_274[k]
                   - f_9 * hsg1_274[k]
                   + f_3 * pc_x[k] * hsh_382[k];

        t_509[k] = f_8 * hsg0_275[k]
                   - f_9 * hsg1_275[k]
                   + f_3 * pc_x[k] * hsh_383[k];

        t_510[k] = f_6 * hsg0_276[k]
                   - f_7 * hsg1_276[k]
                   + f_3 * pc_x[k] * hsh_384[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pc_x, hsg0_277, hsg0_278, hsg0_279, hsg1_277, \
                         hsg1_278, hsg1_279, hsh_385, hsh_386, \
                         hsh_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_6 * hsg0_277[k]
                   - f_7 * hsg1_277[k]
                   + f_3 * pc_x[k] * hsh_385[k];

        t_512[k] = f_6 * hsg0_278[k]
                   - f_7 * hsg1_278[k]
                   + f_3 * pc_x[k] * hsh_386[k];

        t_513[k] = f_6 * hsg0_279[k]
                   - f_7 * hsg1_279[k]
                   + f_3 * pc_x[k] * hsh_387[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, pc_x, hsg0_280, hsg0_281, hsg0_282, hsg1_280, \
                         hsg1_281, hsg1_282, hsh_388, hsh_389, \
                         hsh_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_4 * hsg0_280[k]
                   - f_5 * hsg1_280[k]
                   + f_3 * pc_x[k] * hsh_388[k];

        t_515[k] = f_4 * hsg0_281[k]
                   - f_5 * hsg1_281[k]
                   + f_3 * pc_x[k] * hsh_389[k];

        t_516[k] = f_4 * hsg0_282[k]
                   - f_5 * hsg1_282[k]
                   + f_3 * pc_x[k] * hsh_390[k];
    }

#pragma omp simd aligned(t_517, t_518, t_519, t_520, t_521, pc_x, hsg0_283, hsg0_284, \
                         hsg1_283, hsg1_284, hsh_391, hsh_392, hsh_393, hsh_394, \
                         hsh_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_517[k] = f_4 * hsg0_283[k]
                   - f_5 * hsg1_283[k]
                   + f_3 * pc_x[k] * hsh_391[k];

        t_518[k] = f_4 * hsg0_284[k]
                   - f_5 * hsg1_284[k]
                   + f_3 * pc_x[k] * hsh_392[k];

        t_519[k] = f_3 * pc_x[k] * hsh_393[k];

        t_520[k] = f_3 * pc_x[k] * hsh_394[k];

        t_521[k] = f_3 * pc_x[k] * hsh_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, t_526, pc_x, pc_y, pc_z, gsh_267, \
                         gsh_288, hsg0_280, hsg1_280, hsh_393, hsh_396, hsh_397, \
                         hsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_3 * pc_x[k] * hsh_396[k];

        t_523[k] = f_3 * pc_x[k] * hsh_397[k];

        t_524[k] = f_3 * pc_x[k] * hsh_398[k];

        t_525[k] = f_12 * gsh_288[k]
                   + f_1 * hsg0_280[k]
                   - f_2 * hsg1_280[k]
                   + f_3 * pc_y[k] * hsh_393[k];

        t_526[k] = f_13 * gsh_267[k]
                   + f_3 * pc_z[k] * hsh_393[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, gsh_290, gsh_291, gsh_292, hsg0_282, \
                         hsg0_283, hsg0_284, hsg1_282, hsg1_283, hsg1_284, hsh_395, hsh_396, \
                         hsh_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_12 * gsh_290[k]
                   + f_8 * hsg0_282[k]
                   - f_9 * hsg1_282[k]
                   + f_3 * pc_y[k] * hsh_395[k];

        t_528[k] = f_12 * gsh_291[k]
                   + f_6 * hsg0_283[k]
                   - f_7 * hsg1_283[k]
                   + f_3 * pc_y[k] * hsh_396[k];

        t_529[k] = f_12 * gsh_292[k]
                   + f_4 * hsg0_284[k]
                   - f_5 * hsg1_284[k]
                   + f_3 * pc_y[k] * hsh_397[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pa_y, pc_y, pc_z, gsi0_392, gsh_272, gsh_293, \
                         gsi1_392, hsg0_284, hsg1_284, hsh_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_12 * gsh_293[k]
                   + f_3 * pc_y[k] * hsh_398[k];

        t_531[k] = f_13 * gsh_272[k]
                   + f_1 * hsg0_284[k]
                   - f_2 * hsg1_284[k]
                   + f_3 * pc_z[k] * hsh_398[k];

        t_532[k] = pa_y[k] * gsi0_392[k]
                   - f_10 * pc_y[k] * gsi1_392[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pa_y, pc_x, pc_y, gsi0_394, gsi1_394, hsg0_286, \
                         hsg0_288, hsg1_286, hsg1_288, hsh_400, \
                         hsh_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_15 * hsg0_286[k]
                   - f_16 * hsg1_286[k]
                   + f_3 * pc_x[k] * hsh_400[k];

        t_534[k] = pa_y[k] * gsi0_394[k]
                   - f_10 * pc_y[k] * gsi1_394[k];

        t_535[k] = f_8 * hsg0_288[k]
                   - f_9 * hsg1_288[k]
                   + f_3 * pc_x[k] * hsh_402[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pa_y, pc_x, pc_y, gsi0_397, gsi1_397, hsg0_289, \
                         hsg0_291, hsg1_289, hsg1_291, hsh_403, \
                         hsh_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_8 * hsg0_289[k]
                   - f_9 * hsg1_289[k]
                   + f_3 * pc_x[k] * hsh_403[k];

        t_537[k] = pa_y[k] * gsi0_397[k]
                   - f_10 * pc_y[k] * gsi1_397[k];

        t_538[k] = f_6 * hsg0_291[k]
                   - f_7 * hsg1_291[k]
                   + f_3 * pc_x[k] * hsh_405[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pa_y, pc_x, pc_y, gsi0_401, gsi1_401, hsg0_292, \
                         hsg0_293, hsg1_292, hsg1_293, hsh_406, \
                         hsh_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_6 * hsg0_292[k]
                   - f_7 * hsg1_292[k]
                   + f_3 * pc_x[k] * hsh_406[k];

        t_540[k] = f_6 * hsg0_293[k]
                   - f_7 * hsg1_293[k]
                   + f_3 * pc_x[k] * hsh_407[k];

        t_541[k] = pa_y[k] * gsi0_401[k]
                   - f_10 * pc_y[k] * gsi1_401[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, hsg0_295, hsg0_296, hsg0_297, hsg1_295, \
                         hsg1_296, hsg1_297, hsh_409, hsh_410, \
                         hsh_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_4 * hsg0_295[k]
                   - f_5 * hsg1_295[k]
                   + f_3 * pc_x[k] * hsh_409[k];

        t_543[k] = f_4 * hsg0_296[k]
                   - f_5 * hsg1_296[k]
                   + f_3 * pc_x[k] * hsh_410[k];

        t_544[k] = f_4 * hsg0_297[k]
                   - f_5 * hsg1_297[k]
                   + f_3 * pc_x[k] * hsh_411[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, pa_y, pc_x, pc_y, gsi0_406, \
                         gsi1_406, hsg0_298, hsg1_298, hsh_412, hsh_414, hsh_415, \
                         hsh_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * hsg0_298[k]
                   - f_5 * hsg1_298[k]
                   + f_3 * pc_x[k] * hsh_412[k];

        t_546[k] = pa_y[k] * gsi0_406[k]
                   - f_10 * pc_y[k] * gsi1_406[k];

        t_547[k] = f_3 * pc_x[k] * hsh_414[k];

        t_548[k] = f_3 * pc_x[k] * hsh_415[k];

        t_549[k] = f_3 * pc_x[k] * hsh_416[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pa_y, pc_x, pc_y, gsi0_413, gsh_309, \
                         gsi1_413, hsh_417, hsh_418, hsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_3 * pc_x[k] * hsh_417[k];

        t_551[k] = f_3 * pc_x[k] * hsh_418[k];

        t_552[k] = f_3 * pc_x[k] * hsh_419[k];

        t_553[k] = pa_y[k] * gsi0_413[k]
                   + f_17 * gsh_309[k]
                   - f_10 * pc_y[k] * gsi1_413[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pa_y, pc_y, pc_z, gsi0_415, gsi0_416, gsh_288, \
                         gsh_311, gsh_312, gsi1_415, gsi1_416, \
                         hsh_414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_14 * gsh_288[k]
                   + f_3 * pc_z[k] * hsh_414[k];

        t_555[k] = pa_y[k] * gsi0_415[k]
                   + f_14 * gsh_311[k]
                   - f_10 * pc_y[k] * gsi1_415[k];

        t_556[k] = pa_y[k] * gsi0_416[k]
                   + f_13 * gsh_312[k]
                   - f_10 * pc_y[k] * gsi1_416[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pa_y, pc_y, gsi0_417, gsi0_419, gsh_313, \
                         gsh_314, gsi1_417, gsi1_419, hsh_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = pa_y[k] * gsi0_417[k]
                   + f_12 * gsh_313[k]
                   - f_10 * pc_y[k] * gsi1_417[k];

        t_558[k] = f_11 * gsh_314[k]
                   + f_3 * pc_y[k] * hsh_419[k];

        t_559[k] = pa_y[k] * gsi0_419[k]
                   - f_10 * pc_y[k] * gsi1_419[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, pc_x, pc_y, hsg0_300, hsg0_302, \
                         hsg0_303, hsg1_300, hsg1_302, hsg1_303, hsh_420, hsh_422, \
                         hsh_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_1 * hsg0_300[k]
                   - f_2 * hsg1_300[k]
                   + f_3 * pc_x[k] * hsh_420[k];

        t_561[k] = f_3 * pc_y[k] * hsh_420[k];

        t_562[k] = f_15 * hsg0_302[k]
                   - f_16 * hsg1_302[k]
                   + f_3 * pc_x[k] * hsh_422[k];

        t_563[k] = f_8 * hsg0_303[k]
                   - f_9 * hsg1_303[k]
                   + f_3 * pc_x[k] * hsh_423[k];

        t_564[k] = f_3 * pc_y[k] * hsh_422[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, t_568, pc_x, pc_y, hsg0_305, hsg0_306, hsg0_307, \
                         hsg1_305, hsg1_306, hsg1_307, hsh_425, hsh_426, \
                         hsh_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_8 * hsg0_305[k]
                   - f_9 * hsg1_305[k]
                   + f_3 * pc_x[k] * hsh_425[k];

        t_566[k] = f_6 * hsg0_306[k]
                   - f_7 * hsg1_306[k]
                   + f_3 * pc_x[k] * hsh_426[k];

        t_567[k] = f_6 * hsg0_307[k]
                   - f_7 * hsg1_307[k]
                   + f_3 * pc_x[k] * hsh_427[k];

        t_568[k] = f_3 * pc_y[k] * hsh_425[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pc_x, hsg0_309, hsg0_310, hsg0_311, hsg1_309, \
                         hsg1_310, hsg1_311, hsh_429, hsh_430, \
                         hsh_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_6 * hsg0_309[k]
                   - f_7 * hsg1_309[k]
                   + f_3 * pc_x[k] * hsh_429[k];

        t_570[k] = f_4 * hsg0_310[k]
                   - f_5 * hsg1_310[k]
                   + f_3 * pc_x[k] * hsh_430[k];

        t_571[k] = f_4 * hsg0_311[k]
                   - f_5 * hsg1_311[k]
                   + f_3 * pc_x[k] * hsh_431[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, t_576, pc_x, pc_y, hsg0_312, hsg0_314, \
                         hsg1_312, hsg1_314, hsh_429, hsh_432, hsh_434, hsh_435, \
                         hsh_436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_4 * hsg0_312[k]
                   - f_5 * hsg1_312[k]
                   + f_3 * pc_x[k] * hsh_432[k];

        t_573[k] = f_3 * pc_y[k] * hsh_429[k];

        t_574[k] = f_4 * hsg0_314[k]
                   - f_5 * hsg1_314[k]
                   + f_3 * pc_x[k] * hsh_434[k];

        t_575[k] = f_3 * pc_x[k] * hsh_435[k];

        t_576[k] = f_3 * pc_x[k] * hsh_436[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, t_581, pc_x, pc_y, hsg0_310, hsg1_310, \
                         hsh_435, hsh_437, hsh_438, hsh_439, hsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_3 * pc_x[k] * hsh_437[k];

        t_578[k] = f_3 * pc_x[k] * hsh_438[k];

        t_579[k] = f_3 * pc_x[k] * hsh_439[k];

        t_580[k] = f_3 * pc_x[k] * hsh_440[k];

        t_581[k] = f_1 * hsg0_310[k]
                   - f_2 * hsg1_310[k]
                   + f_3 * pc_y[k] * hsh_435[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, hsg0_311, hsg0_312, hsg0_313, hsg1_311, \
                         hsg1_312, hsg1_313, hsh_436, hsh_437, \
                         hsh_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_15 * hsg0_311[k]
                   - f_16 * hsg1_311[k]
                   + f_3 * pc_y[k] * hsh_436[k];

        t_583[k] = f_8 * hsg0_312[k]
                   - f_9 * hsg1_312[k]
                   + f_3 * pc_y[k] * hsh_437[k];

        t_584[k] = f_6 * hsg0_313[k]
                   - f_7 * hsg1_313[k]
                   + f_3 * pc_y[k] * hsh_438[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_y, pc_z, gsh_314, hsg0_314, hsg1_314, \
                         hsh_439, hsh_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_4 * hsg0_314[k]
                   - f_5 * hsg1_314[k]
                   + f_3 * pc_y[k] * hsh_439[k];

        t_586[k] = f_3 * pc_y[k] * hsh_440[k];

        t_587[k] = f_0 * gsh_314[k]
                   + f_1 * hsg0_314[k]
                   - f_2 * hsg1_314[k]
                   + f_3 * pc_z[k] * hsh_440[k];
    }
}

auto
compute_prim_hsi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t gsi0, const size_t gsh,
                                                   const size_t gsi1, const size_t hsg0,
                                                   const size_t hsg1, const size_t hsh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_hsi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, gsi0, gsh,
                                                              gsi1, hsg0, hsg1, hsh, ncols,
                                                              gamma, p, q);

    compute_prim_hsi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, gsi0, gsh,
                                                              gsi1, hsg0, hsg1, hsh, ncols,
                                                              gamma, p, q);

    compute_prim_hsi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, gsi0, gsh,
                                                              gsi1, hsg0, hsg1, hsh, ncols,
                                                              gamma, p, q);

    compute_prim_hsi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, gsi0, gsh,
                                                              gsi1, hsg0, hsg1, hsh, ncols,
                                                              gamma, p, q);

    compute_prim_hsi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, gsi0, gsh,
                                                              gsi1, hsg0, hsg1, hsh, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
