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


#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfh0,
                                                          const size_t sfg, const size_t sfh1,
                                                          const size_t sgf0, const size_t sgf1,
                                                          const size_t sgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_0 = buffer.data(sfh0 + 0);
    const auto *sfh0_3 = buffer.data(sfh0 + 3);
    const auto *sfh0_5 = buffer.data(sfh0 + 5);
    const auto *sfh0_6 = buffer.data(sfh0 + 6);
    const auto *sfh0_9 = buffer.data(sfh0 + 9);
    const auto *sfh0_15 = buffer.data(sfh0 + 15);
    const auto *sfh0_20 = buffer.data(sfh0 + 20);
    const auto *sfh0_24 = buffer.data(sfh0 + 24);
    const auto *sfh0_27 = buffer.data(sfh0 + 27);
    const auto *sfh0_36 = buffer.data(sfh0 + 36);
    const auto *sfh0_42 = buffer.data(sfh0 + 42);
    const auto *sfh0_47 = buffer.data(sfh0 + 47);
    const auto *sfh0_51 = buffer.data(sfh0 + 51);
    const auto *sfh0_62 = buffer.data(sfh0 + 62);

    const auto *sfg_0 = buffer.data(sfg + 0);
    const auto *sfg_1 = buffer.data(sfg + 1);
    const auto *sfg_2 = buffer.data(sfg + 2);
    const auto *sfg_3 = buffer.data(sfg + 3);
    const auto *sfg_5 = buffer.data(sfg + 5);
    const auto *sfg_6 = buffer.data(sfg + 6);
    const auto *sfg_9 = buffer.data(sfg + 9);
    const auto *sfg_10 = buffer.data(sfg + 10);
    const auto *sfg_11 = buffer.data(sfg + 11);
    const auto *sfg_12 = buffer.data(sfg + 12);
    const auto *sfg_13 = buffer.data(sfg + 13);
    const auto *sfg_14 = buffer.data(sfg + 14);
    const auto *sfg_15 = buffer.data(sfg + 15);
    const auto *sfg_17 = buffer.data(sfg + 17);
    const auto *sfg_18 = buffer.data(sfg + 18);
    const auto *sfg_20 = buffer.data(sfg + 20);
    const auto *sfg_25 = buffer.data(sfg + 25);
    const auto *sfg_26 = buffer.data(sfg + 26);
    const auto *sfg_27 = buffer.data(sfg + 27);
    const auto *sfg_28 = buffer.data(sfg + 28);
    const auto *sfg_29 = buffer.data(sfg + 29);
    const auto *sfg_30 = buffer.data(sfg + 30);
    const auto *sfg_32 = buffer.data(sfg + 32);
    const auto *sfg_33 = buffer.data(sfg + 33);
    const auto *sfg_35 = buffer.data(sfg + 35);
    const auto *sfg_40 = buffer.data(sfg + 40);
    const auto *sfg_41 = buffer.data(sfg + 41);
    const auto *sfg_42 = buffer.data(sfg + 42);
    const auto *sfg_43 = buffer.data(sfg + 43);
    const auto *sfg_44 = buffer.data(sfg + 44);
    const auto *sfg_45 = buffer.data(sfg + 45);
    const auto *sfg_48 = buffer.data(sfg + 48);
    const auto *sfg_50 = buffer.data(sfg + 50);
    const auto *sfg_51 = buffer.data(sfg + 51);
    const auto *sfg_54 = buffer.data(sfg + 54);
    const auto *sfg_55 = buffer.data(sfg + 55);
    const auto *sfg_56 = buffer.data(sfg + 56);
    const auto *sfg_57 = buffer.data(sfg + 57);
    const auto *sfg_58 = buffer.data(sfg + 58);
    const auto *sfg_59 = buffer.data(sfg + 59);
    const auto *sfg_70 = buffer.data(sfg + 70);
    const auto *sfg_71 = buffer.data(sfg + 71);
    const auto *sfg_72 = buffer.data(sfg + 72);
    const auto *sfg_73 = buffer.data(sfg + 73);
    const auto *sfg_74 = buffer.data(sfg + 74);
    const auto *sfg_75 = buffer.data(sfg + 75);
    const auto *sfg_78 = buffer.data(sfg + 78);
    const auto *sfg_80 = buffer.data(sfg + 80);
    const auto *sfg_81 = buffer.data(sfg + 81);
    const auto *sfg_84 = buffer.data(sfg + 84);
    const auto *sfg_85 = buffer.data(sfg + 85);
    const auto *sfg_86 = buffer.data(sfg + 86);
    const auto *sfg_87 = buffer.data(sfg + 87);
    const auto *sfg_88 = buffer.data(sfg + 88);
    const auto *sfg_89 = buffer.data(sfg + 89);

    const auto *sfh1_0 = buffer.data(sfh1 + 0);
    const auto *sfh1_3 = buffer.data(sfh1 + 3);
    const auto *sfh1_5 = buffer.data(sfh1 + 5);
    const auto *sfh1_6 = buffer.data(sfh1 + 6);
    const auto *sfh1_9 = buffer.data(sfh1 + 9);
    const auto *sfh1_15 = buffer.data(sfh1 + 15);
    const auto *sfh1_20 = buffer.data(sfh1 + 20);
    const auto *sfh1_24 = buffer.data(sfh1 + 24);
    const auto *sfh1_27 = buffer.data(sfh1 + 27);
    const auto *sfh1_36 = buffer.data(sfh1 + 36);
    const auto *sfh1_42 = buffer.data(sfh1 + 42);
    const auto *sfh1_47 = buffer.data(sfh1 + 47);
    const auto *sfh1_51 = buffer.data(sfh1 + 51);
    const auto *sfh1_62 = buffer.data(sfh1 + 62);

    const auto *sgf0_0 = buffer.data(sgf0 + 0);
    const auto *sgf0_3 = buffer.data(sgf0 + 3);
    const auto *sgf0_5 = buffer.data(sgf0 + 5);
    const auto *sgf0_6 = buffer.data(sgf0 + 6);
    const auto *sgf0_8 = buffer.data(sgf0 + 8);
    const auto *sgf0_9 = buffer.data(sgf0 + 9);
    const auto *sgf0_16 = buffer.data(sgf0 + 16);
    const auto *sgf0_18 = buffer.data(sgf0 + 18);
    const auto *sgf0_19 = buffer.data(sgf0 + 19);
    const auto *sgf0_28 = buffer.data(sgf0 + 28);
    const auto *sgf0_29 = buffer.data(sgf0 + 29);
    const auto *sgf0_30 = buffer.data(sgf0 + 30);
    const auto *sgf0_33 = buffer.data(sgf0 + 33);
    const auto *sgf0_35 = buffer.data(sgf0 + 35);
    const auto *sgf0_36 = buffer.data(sgf0 + 36);
    const auto *sgf0_38 = buffer.data(sgf0 + 38);
    const auto *sgf0_39 = buffer.data(sgf0 + 39);
    const auto *sgf0_48 = buffer.data(sgf0 + 48);
    const auto *sgf0_49 = buffer.data(sgf0 + 49);
    const auto *sgf0_50 = buffer.data(sgf0 + 50);
    const auto *sgf0_53 = buffer.data(sgf0 + 53);
    const auto *sgf0_55 = buffer.data(sgf0 + 55);
    const auto *sgf0_56 = buffer.data(sgf0 + 56);
    const auto *sgf0_58 = buffer.data(sgf0 + 58);
    const auto *sgf0_59 = buffer.data(sgf0 + 59);

    const auto *sgf1_0 = buffer.data(sgf1 + 0);
    const auto *sgf1_3 = buffer.data(sgf1 + 3);
    const auto *sgf1_5 = buffer.data(sgf1 + 5);
    const auto *sgf1_6 = buffer.data(sgf1 + 6);
    const auto *sgf1_8 = buffer.data(sgf1 + 8);
    const auto *sgf1_9 = buffer.data(sgf1 + 9);
    const auto *sgf1_16 = buffer.data(sgf1 + 16);
    const auto *sgf1_18 = buffer.data(sgf1 + 18);
    const auto *sgf1_19 = buffer.data(sgf1 + 19);
    const auto *sgf1_28 = buffer.data(sgf1 + 28);
    const auto *sgf1_29 = buffer.data(sgf1 + 29);
    const auto *sgf1_30 = buffer.data(sgf1 + 30);
    const auto *sgf1_33 = buffer.data(sgf1 + 33);
    const auto *sgf1_35 = buffer.data(sgf1 + 35);
    const auto *sgf1_36 = buffer.data(sgf1 + 36);
    const auto *sgf1_38 = buffer.data(sgf1 + 38);
    const auto *sgf1_39 = buffer.data(sgf1 + 39);
    const auto *sgf1_48 = buffer.data(sgf1 + 48);
    const auto *sgf1_49 = buffer.data(sgf1 + 49);
    const auto *sgf1_50 = buffer.data(sgf1 + 50);
    const auto *sgf1_53 = buffer.data(sgf1 + 53);
    const auto *sgf1_55 = buffer.data(sgf1 + 55);
    const auto *sgf1_56 = buffer.data(sgf1 + 56);
    const auto *sgf1_58 = buffer.data(sgf1 + 58);
    const auto *sgf1_59 = buffer.data(sgf1 + 59);

    const auto *sgg_0 = buffer.data(sgg + 0);
    const auto *sgg_2 = buffer.data(sgg + 2);
    const auto *sgg_3 = buffer.data(sgg + 3);
    const auto *sgg_5 = buffer.data(sgg + 5);
    const auto *sgg_6 = buffer.data(sgg + 6);
    const auto *sgg_9 = buffer.data(sgg + 9);
    const auto *sgg_10 = buffer.data(sgg + 10);
    const auto *sgg_11 = buffer.data(sgg + 11);
    const auto *sgg_12 = buffer.data(sgg + 12);
    const auto *sgg_13 = buffer.data(sgg + 13);
    const auto *sgg_14 = buffer.data(sgg + 14);
    const auto *sgg_15 = buffer.data(sgg + 15);
    const auto *sgg_17 = buffer.data(sgg + 17);
    const auto *sgg_18 = buffer.data(sgg + 18);
    const auto *sgg_20 = buffer.data(sgg + 20);
    const auto *sgg_25 = buffer.data(sgg + 25);
    const auto *sgg_26 = buffer.data(sgg + 26);
    const auto *sgg_27 = buffer.data(sgg + 27);
    const auto *sgg_28 = buffer.data(sgg + 28);
    const auto *sgg_29 = buffer.data(sgg + 29);
    const auto *sgg_30 = buffer.data(sgg + 30);
    const auto *sgg_32 = buffer.data(sgg + 32);
    const auto *sgg_33 = buffer.data(sgg + 33);
    const auto *sgg_35 = buffer.data(sgg + 35);
    const auto *sgg_40 = buffer.data(sgg + 40);
    const auto *sgg_41 = buffer.data(sgg + 41);
    const auto *sgg_42 = buffer.data(sgg + 42);
    const auto *sgg_43 = buffer.data(sgg + 43);
    const auto *sgg_44 = buffer.data(sgg + 44);
    const auto *sgg_45 = buffer.data(sgg + 45);
    const auto *sgg_47 = buffer.data(sgg + 47);
    const auto *sgg_48 = buffer.data(sgg + 48);
    const auto *sgg_50 = buffer.data(sgg + 50);
    const auto *sgg_51 = buffer.data(sgg + 51);
    const auto *sgg_54 = buffer.data(sgg + 54);
    const auto *sgg_55 = buffer.data(sgg + 55);
    const auto *sgg_56 = buffer.data(sgg + 56);
    const auto *sgg_57 = buffer.data(sgg + 57);
    const auto *sgg_58 = buffer.data(sgg + 58);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_60 = buffer.data(sgg + 60);
    const auto *sgg_62 = buffer.data(sgg + 62);
    const auto *sgg_63 = buffer.data(sgg + 63);
    const auto *sgg_65 = buffer.data(sgg + 65);
    const auto *sgg_70 = buffer.data(sgg + 70);
    const auto *sgg_71 = buffer.data(sgg + 71);
    const auto *sgg_72 = buffer.data(sgg + 72);
    const auto *sgg_73 = buffer.data(sgg + 73);
    const auto *sgg_74 = buffer.data(sgg + 74);
    const auto *sgg_75 = buffer.data(sgg + 75);
    const auto *sgg_77 = buffer.data(sgg + 77);
    const auto *sgg_78 = buffer.data(sgg + 78);
    const auto *sgg_80 = buffer.data(sgg + 80);
    const auto *sgg_81 = buffer.data(sgg + 81);
    const auto *sgg_84 = buffer.data(sgg + 84);
    const auto *sgg_85 = buffer.data(sgg + 85);
    const auto *sgg_86 = buffer.data(sgg + 86);
    const auto *sgg_87 = buffer.data(sgg + 87);
    const auto *sgg_88 = buffer.data(sgg + 88);
    const auto *sgg_89 = buffer.data(sgg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sfg_0, sfg_3, sgf0_0, sgf0_3, \
                         sgf1_0, sgf1_3, sgg_0, sgg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfg_0[k]
                 + f_1 * sgf0_0[k]
                 - f_2 * sgf1_0[k]
                 + f_3 * pc_x[k] * sgg_0[k];

        t_1[k] = f_3 * pc_y[k] * sgg_0[k];

        t_2[k] = f_3 * pc_z[k] * sgg_0[k];

        t_3[k] = f_0 * sfg_3[k]
                 + f_4 * sgf0_3[k]
                 - f_5 * sgf1_3[k]
                 + f_3 * pc_x[k] * sgg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sfg_5, sfg_6, sgf0_5, sgf0_6, sgf1_5, \
                         sgf1_6, sgg_2, sgg_5, sgg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sgg_2[k];

        t_5[k] = f_0 * sfg_5[k]
                 + f_4 * sgf0_5[k]
                 - f_5 * sgf1_5[k]
                 + f_3 * pc_x[k] * sgg_5[k];

        t_6[k] = f_0 * sfg_6[k]
                 + f_6 * sgf0_6[k]
                 - f_7 * sgf1_6[k]
                 + f_3 * pc_x[k] * sgg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, sfg_9, sfg_10, sgf0_9, sgf1_9, \
                         sgg_3, sgg_5, sgg_9, sgg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sgg_3[k];

        t_8[k] = f_3 * pc_y[k] * sgg_5[k];

        t_9[k] = f_0 * sfg_9[k]
                 + f_6 * sgf0_9[k]
                 - f_7 * sgf1_9[k]
                 + f_3 * pc_x[k] * sgg_9[k];

        t_10[k] = f_0 * sfg_10[k]
                  + f_3 * pc_x[k] * sgg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, sfg_11, sfg_12, sfg_13, sfg_14, sgg_11, \
                         sgg_12, sgg_13, sgg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sfg_11[k]
                  + f_3 * pc_x[k] * sgg_11[k];

        t_12[k] = f_0 * sfg_12[k]
                  + f_3 * pc_x[k] * sgg_12[k];

        t_13[k] = f_0 * sfg_13[k]
                  + f_3 * pc_x[k] * sgg_13[k];

        t_14[k] = f_0 * sfg_14[k]
                  + f_3 * pc_x[k] * sgg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, sgf0_6, sgf0_8, sgf0_9, sgf1_6, \
                         sgf1_8, sgf1_9, sgg_10, sgg_12, sgg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sgf0_6[k]
                  - f_2 * sgf1_6[k]
                  + f_3 * pc_y[k] * sgg_10[k];

        t_16[k] = f_3 * pc_z[k] * sgg_10[k];

        t_17[k] = f_4 * sgf0_8[k]
                  - f_5 * sgf1_8[k]
                  + f_3 * pc_y[k] * sgg_12[k];

        t_18[k] = f_6 * sgf0_9[k]
                  - f_7 * sgf1_9[k]
                  + f_3 * pc_y[k] * sgg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, sfh0_0, sfg_0, \
                         sfh1_0, sgf0_9, sgf1_9, sgg_14, sgg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sgg_14[k];

        t_20[k] = f_1 * sgf0_9[k]
                  - f_2 * sgf1_9[k]
                  + f_3 * pc_z[k] * sgg_14[k];

        t_21[k] = pb_y[k] * sfh0_0[k]
                  - f_8 * pc_y[k] * sfh1_0[k];

        t_22[k] = f_9 * sfg_0[k]
                  + f_3 * pc_y[k] * sgg_15[k];

        t_23[k] = f_3 * pc_z[k] * sgg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, sfh0_3, sfh0_5, sfh0_6, sfg_1, \
                         sfg_2, sfg_3, sfh1_3, sfh1_5, sfh1_6, sgg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sfh0_3[k]
                  + f_10 * sfg_1[k]
                  - f_8 * pc_y[k] * sfh1_3[k];

        t_25[k] = f_9 * sfg_2[k]
                  + f_3 * pc_y[k] * sgg_17[k];

        t_26[k] = pb_y[k] * sfh0_5[k]
                  - f_8 * pc_y[k] * sfh1_5[k];

        t_27[k] = pb_y[k] * sfh0_6[k]
                  + f_11 * sfg_3[k]
                  - f_8 * pc_y[k] * sfh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, sfh0_9, sfg_5, \
                         sfg_25, sfh1_9, sgg_18, sgg_20, sgg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * sgg_18[k];

        t_29[k] = f_9 * sfg_5[k]
                  + f_3 * pc_y[k] * sgg_20[k];

        t_30[k] = pb_y[k] * sfh0_9[k]
                  - f_8 * pc_y[k] * sfh1_9[k];

        t_31[k] = f_11 * sfg_25[k]
                  + f_3 * pc_x[k] * sgg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, sfg_26, sfg_27, sfg_28, sfg_29, sgg_26, \
                         sgg_27, sgg_28, sgg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_11 * sfg_26[k]
                  + f_3 * pc_x[k] * sgg_26[k];

        t_33[k] = f_11 * sfg_27[k]
                  + f_3 * pc_x[k] * sgg_27[k];

        t_34[k] = f_11 * sfg_28[k]
                  + f_3 * pc_x[k] * sgg_28[k];

        t_35[k] = f_11 * sfg_29[k]
                  + f_3 * pc_x[k] * sgg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, sfg_10, sfg_12, sgf0_16, sgf0_18, \
                         sgf1_16, sgf1_18, sgg_25, sgg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * sfg_10[k]
                  + f_1 * sgf0_16[k]
                  - f_2 * sgf1_16[k]
                  + f_3 * pc_y[k] * sgg_25[k];

        t_37[k] = f_3 * pc_z[k] * sgg_25[k];

        t_38[k] = f_9 * sfg_12[k]
                  + f_4 * sgf0_18[k]
                  - f_5 * sgf1_18[k]
                  + f_3 * pc_y[k] * sgg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, sfh0_20, sfg_13, sfg_14, sfh1_20, \
                         sgf0_19, sgf1_19, sgg_28, sgg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * sfg_13[k]
                  + f_6 * sgf0_19[k]
                  - f_7 * sgf1_19[k]
                  + f_3 * pc_y[k] * sgg_28[k];

        t_40[k] = f_9 * sfg_14[k]
                  + f_3 * pc_y[k] * sgg_29[k];

        t_41[k] = pb_y[k] * sfh0_20[k]
                  - f_8 * pc_y[k] * sfh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, sfh0_0, sfh0_3, \
                         sfg_0, sfh1_0, sfh1_3, sgg_30, sgg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sfh0_0[k]
                  - f_8 * pc_z[k] * sfh1_0[k];

        t_43[k] = f_3 * pc_y[k] * sgg_30[k];

        t_44[k] = f_9 * sfg_0[k]
                  + f_3 * pc_z[k] * sgg_30[k];

        t_45[k] = pb_z[k] * sfh0_3[k]
                  - f_8 * pc_z[k] * sfh1_3[k];

        t_46[k] = f_3 * pc_y[k] * sgg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, sfh0_5, sfh0_6, sfg_2, \
                         sfg_3, sfh1_5, sfh1_6, sgg_33, sgg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * sfh0_5[k]
                  + f_10 * sfg_2[k]
                  - f_8 * pc_z[k] * sfh1_5[k];

        t_48[k] = pb_z[k] * sfh0_6[k]
                  - f_8 * pc_z[k] * sfh1_6[k];

        t_49[k] = f_9 * sfg_3[k]
                  + f_3 * pc_z[k] * sgg_33[k];

        t_50[k] = f_3 * pc_y[k] * sgg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, sfh0_9, sfg_5, sfg_40, \
                         sfg_41, sfg_42, sfh1_9, sgg_40, sgg_41, \
                         sgg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * sfh0_9[k]
                  + f_11 * sfg_5[k]
                  - f_8 * pc_z[k] * sfh1_9[k];

        t_52[k] = f_11 * sfg_40[k]
                  + f_3 * pc_x[k] * sgg_40[k];

        t_53[k] = f_11 * sfg_41[k]
                  + f_3 * pc_x[k] * sgg_41[k];

        t_54[k] = f_11 * sfg_42[k]
                  + f_3 * pc_x[k] * sgg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, sfh0_15, sfg_10, sfg_43, \
                         sfg_44, sfh1_15, sgg_40, sgg_43, sgg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_11 * sfg_43[k]
                  + f_3 * pc_x[k] * sgg_43[k];

        t_56[k] = f_11 * sfg_44[k]
                  + f_3 * pc_x[k] * sgg_44[k];

        t_57[k] = pb_z[k] * sfh0_15[k]
                  - f_8 * pc_z[k] * sfh1_15[k];

        t_58[k] = f_9 * sfg_10[k]
                  + f_3 * pc_z[k] * sgg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, sfg_14, sgf0_28, sgf0_29, \
                         sgf1_28, sgf1_29, sgg_42, sgg_43, sgg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * sgf0_28[k]
                  - f_5 * sgf1_28[k]
                  + f_3 * pc_y[k] * sgg_42[k];

        t_60[k] = f_6 * sgf0_29[k]
                  - f_7 * sgf1_29[k]
                  + f_3 * pc_y[k] * sgg_43[k];

        t_61[k] = f_3 * pc_y[k] * sgg_44[k];

        t_62[k] = f_9 * sfg_14[k]
                  + f_1 * sgf0_29[k]
                  - f_2 * sgf1_29[k]
                  + f_3 * pc_z[k] * sgg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, sfg_15, sfg_45, sfg_48, \
                         sgf0_30, sgf0_33, sgf1_30, sgf1_33, sgg_45, \
                         sgg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_10 * sfg_45[k]
                  + f_1 * sgf0_30[k]
                  - f_2 * sgf1_30[k]
                  + f_3 * pc_x[k] * sgg_45[k];

        t_64[k] = f_10 * sfg_15[k]
                  + f_3 * pc_y[k] * sgg_45[k];

        t_65[k] = f_3 * pc_z[k] * sgg_45[k];

        t_66[k] = f_10 * sfg_48[k]
                  + f_4 * sgf0_33[k]
                  - f_5 * sgf1_33[k]
                  + f_3 * pc_x[k] * sgg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, sfg_17, sfg_50, sfg_51, sgf0_35, \
                         sgf0_36, sgf1_35, sgf1_36, sgg_47, sgg_50, \
                         sgg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * sfg_17[k]
                  + f_3 * pc_y[k] * sgg_47[k];

        t_68[k] = f_10 * sfg_50[k]
                  + f_4 * sgf0_35[k]
                  - f_5 * sgf1_35[k]
                  + f_3 * pc_x[k] * sgg_50[k];

        t_69[k] = f_10 * sfg_51[k]
                  + f_6 * sgf0_36[k]
                  - f_7 * sgf1_36[k]
                  + f_3 * pc_x[k] * sgg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, sfg_20, sfg_54, sfg_55, \
                         sgf0_39, sgf1_39, sgg_48, sgg_50, sgg_54, \
                         sgg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * sgg_48[k];

        t_71[k] = f_10 * sfg_20[k]
                  + f_3 * pc_y[k] * sgg_50[k];

        t_72[k] = f_10 * sfg_54[k]
                  + f_6 * sgf0_39[k]
                  - f_7 * sgf1_39[k]
                  + f_3 * pc_x[k] * sgg_54[k];

        t_73[k] = f_10 * sfg_55[k]
                  + f_3 * pc_x[k] * sgg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, sfg_56, sfg_57, sfg_58, sfg_59, sgg_56, \
                         sgg_57, sgg_58, sgg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_10 * sfg_56[k]
                  + f_3 * pc_x[k] * sgg_56[k];

        t_75[k] = f_10 * sfg_57[k]
                  + f_3 * pc_x[k] * sgg_57[k];

        t_76[k] = f_10 * sfg_58[k]
                  + f_3 * pc_x[k] * sgg_58[k];

        t_77[k] = f_10 * sfg_59[k]
                  + f_3 * pc_x[k] * sgg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, sfg_25, sfg_27, sgf0_36, sgf0_38, \
                         sgf1_36, sgf1_38, sgg_55, sgg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * sfg_25[k]
                  + f_1 * sgf0_36[k]
                  - f_2 * sgf1_36[k]
                  + f_3 * pc_y[k] * sgg_55[k];

        t_79[k] = f_3 * pc_z[k] * sgg_55[k];

        t_80[k] = f_10 * sfg_27[k]
                  + f_4 * sgf0_38[k]
                  - f_5 * sgf1_38[k]
                  + f_3 * pc_y[k] * sgg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, sfh0_42, sfg_28, sfg_29, \
                         sfh1_42, sgf0_39, sgf1_39, sgg_58, sgg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * sfg_28[k]
                  + f_6 * sgf0_39[k]
                  - f_7 * sgf1_39[k]
                  + f_3 * pc_y[k] * sgg_58[k];

        t_82[k] = f_10 * sfg_29[k]
                  + f_3 * pc_y[k] * sgg_59[k];

        t_83[k] = f_1 * sgf0_39[k]
                  - f_2 * sgf1_39[k]
                  + f_3 * pc_z[k] * sgg_59[k];

        t_84[k] = pb_y[k] * sfh0_42[k]
                  - f_8 * pc_y[k] * sfh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, sfh0_24, sfg_15, sfg_30, \
                         sfg_32, sfh1_24, sgg_60, sgg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * sfg_30[k]
                  + f_3 * pc_y[k] * sgg_60[k];

        t_86[k] = f_9 * sfg_15[k]
                  + f_3 * pc_z[k] * sgg_60[k];

        t_87[k] = pb_z[k] * sfh0_24[k]
                  - f_8 * pc_z[k] * sfh1_24[k];

        t_88[k] = f_9 * sfg_32[k]
                  + f_3 * pc_y[k] * sgg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, sfh0_27, sfh0_47, \
                         sfg_18, sfg_35, sfh1_27, sfh1_47, sgg_63, \
                         sgg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * sfh0_47[k]
                  - f_8 * pc_y[k] * sfh1_47[k];

        t_90[k] = pb_z[k] * sfh0_27[k]
                  - f_8 * pc_z[k] * sfh1_27[k];

        t_91[k] = f_9 * sfg_18[k]
                  + f_3 * pc_z[k] * sgg_63[k];

        t_92[k] = f_9 * sfg_35[k]
                  + f_3 * pc_y[k] * sgg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, sfh0_51, sfg_70, sfg_71, \
                         sfg_72, sfh1_51, sgg_70, sgg_71, sgg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * sfh0_51[k]
                  - f_8 * pc_y[k] * sfh1_51[k];

        t_94[k] = f_10 * sfg_70[k]
                  + f_3 * pc_x[k] * sgg_70[k];

        t_95[k] = f_10 * sfg_71[k]
                  + f_3 * pc_x[k] * sgg_71[k];

        t_96[k] = f_10 * sfg_72[k]
                  + f_3 * pc_x[k] * sgg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, sfh0_36, sfg_25, sfg_73, \
                         sfg_74, sfh1_36, sgg_70, sgg_73, sgg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_10 * sfg_73[k]
                  + f_3 * pc_x[k] * sgg_73[k];

        t_98[k] = f_10 * sfg_74[k]
                  + f_3 * pc_x[k] * sgg_74[k];

        t_99[k] = pb_z[k] * sfh0_36[k]
                  - f_8 * pc_z[k] * sfh1_36[k];

        t_100[k] = f_9 * sfg_25[k]
                   + f_3 * pc_z[k] * sgg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, sfg_42, sfg_43, sfg_44, sgf0_48, sgf0_49, \
                         sgf1_48, sgf1_49, sgg_72, sgg_73, sgg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * sfg_42[k]
                   + f_4 * sgf0_48[k]
                   - f_5 * sgf1_48[k]
                   + f_3 * pc_y[k] * sgg_72[k];

        t_102[k] = f_9 * sfg_43[k]
                   + f_6 * sgf0_49[k]
                   - f_7 * sgf1_49[k]
                   + f_3 * pc_y[k] * sgg_73[k];

        t_103[k] = f_9 * sfg_44[k]
                   + f_3 * pc_y[k] * sgg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, sfh0_62, sfg_30, \
                         sfg_75, sfh1_62, sgf0_50, sgf1_50, sgg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * sfh0_62[k]
                   - f_8 * pc_y[k] * sfh1_62[k];

        t_105[k] = f_10 * sfg_75[k]
                   + f_1 * sgf0_50[k]
                   - f_2 * sgf1_50[k]
                   + f_3 * pc_x[k] * sgg_75[k];

        t_106[k] = f_3 * pc_y[k] * sgg_75[k];

        t_107[k] = f_10 * sfg_30[k]
                   + f_3 * pc_z[k] * sgg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, sfg_78, sfg_80, sgf0_53, sgf0_55, \
                         sgf1_53, sgf1_55, sgg_77, sgg_78, sgg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_10 * sfg_78[k]
                   + f_4 * sgf0_53[k]
                   - f_5 * sgf1_53[k]
                   + f_3 * pc_x[k] * sgg_78[k];

        t_109[k] = f_3 * pc_y[k] * sgg_77[k];

        t_110[k] = f_10 * sfg_80[k]
                   + f_4 * sgf0_55[k]
                   - f_5 * sgf1_55[k]
                   + f_3 * pc_x[k] * sgg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, sfg_33, sfg_81, sgf0_56, \
                         sgf1_56, sgg_78, sgg_80, sgg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_10 * sfg_81[k]
                   + f_6 * sgf0_56[k]
                   - f_7 * sgf1_56[k]
                   + f_3 * pc_x[k] * sgg_81[k];

        t_112[k] = f_10 * sfg_33[k]
                   + f_3 * pc_z[k] * sgg_78[k];

        t_113[k] = f_3 * pc_y[k] * sgg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, sfg_84, sfg_85, sfg_86, sfg_87, \
                         sgf0_59, sgf1_59, sgg_84, sgg_85, sgg_86, \
                         sgg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_10 * sfg_84[k]
                   + f_6 * sgf0_59[k]
                   - f_7 * sgf1_59[k]
                   + f_3 * pc_x[k] * sgg_84[k];

        t_115[k] = f_10 * sfg_85[k]
                   + f_3 * pc_x[k] * sgg_85[k];

        t_116[k] = f_10 * sfg_86[k]
                   + f_3 * pc_x[k] * sgg_86[k];

        t_117[k] = f_10 * sfg_87[k]
                   + f_3 * pc_x[k] * sgg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, sfg_40, sfg_88, sfg_89, \
                         sgf0_56, sgf1_56, sgg_85, sgg_88, sgg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_10 * sfg_88[k]
                   + f_3 * pc_x[k] * sgg_88[k];

        t_119[k] = f_10 * sfg_89[k]
                   + f_3 * pc_x[k] * sgg_89[k];

        t_120[k] = f_1 * sgf0_56[k]
                   - f_2 * sgf1_56[k]
                   + f_3 * pc_y[k] * sgg_85[k];

        t_121[k] = f_10 * sfg_40[k]
                   + f_3 * pc_z[k] * sgg_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, sfg_44, sgf0_58, sgf0_59, \
                         sgf1_58, sgf1_59, sgg_87, sgg_88, sgg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * sgf0_58[k]
                   - f_5 * sgf1_58[k]
                   + f_3 * pc_y[k] * sgg_87[k];

        t_123[k] = f_6 * sgf0_59[k]
                   - f_7 * sgf1_59[k]
                   + f_3 * pc_y[k] * sgg_88[k];

        t_124[k] = f_3 * pc_y[k] * sgg_89[k];

        t_125[k] = f_10 * sfg_44[k]
                   + f_1 * sgf0_59[k]
                   - f_2 * sgf1_59[k]
                   + f_3 * pc_z[k] * sgg_89[k];
    }
}

static auto
compute_prim_sgh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfh0,
                                                          const size_t sfg, const size_t sfh1,
                                                          const size_t sgf0, const size_t sgf1,
                                                          const size_t sgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;

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
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_63 = buffer.data(sfh0 + 63);
    const auto *sfh0_66 = buffer.data(sfh0 + 66);
    const auto *sfh0_69 = buffer.data(sfh0 + 69);
    const auto *sfh0_105 = buffer.data(sfh0 + 105);
    const auto *sfh0_110 = buffer.data(sfh0 + 110);
    const auto *sfh0_114 = buffer.data(sfh0 + 114);
    const auto *sfh0_126 = buffer.data(sfh0 + 126);
    const auto *sfh0_129 = buffer.data(sfh0 + 129);
    const auto *sfh0_131 = buffer.data(sfh0 + 131);
    const auto *sfh0_132 = buffer.data(sfh0 + 132);
    const auto *sfh0_135 = buffer.data(sfh0 + 135);
    const auto *sfh0_141 = buffer.data(sfh0 + 141);
    const auto *sfh0_143 = buffer.data(sfh0 + 143);
    const auto *sfh0_144 = buffer.data(sfh0 + 144);
    const auto *sfh0_146 = buffer.data(sfh0 + 146);
    const auto *sfh0_152 = buffer.data(sfh0 + 152);
    const auto *sfh0_156 = buffer.data(sfh0 + 156);
    const auto *sfh0_162 = buffer.data(sfh0 + 162);
    const auto *sfh0_164 = buffer.data(sfh0 + 164);
    const auto *sfh0_165 = buffer.data(sfh0 + 165);
    const auto *sfh0_167 = buffer.data(sfh0 + 167);
    const auto *sfh0_171 = buffer.data(sfh0 + 171);
    const auto *sfh0_174 = buffer.data(sfh0 + 174);
    const auto *sfh0_183 = buffer.data(sfh0 + 183);
    const auto *sfh0_185 = buffer.data(sfh0 + 185);
    const auto *sfh0_186 = buffer.data(sfh0 + 186);
    const auto *sfh0_188 = buffer.data(sfh0 + 188);
    const auto *sfh0_189 = buffer.data(sfh0 + 189);
    const auto *sfh0_192 = buffer.data(sfh0 + 192);
    const auto *sfh0_194 = buffer.data(sfh0 + 194);
    const auto *sfh0_195 = buffer.data(sfh0 + 195);
    const auto *sfh0_198 = buffer.data(sfh0 + 198);
    const auto *sfh0_204 = buffer.data(sfh0 + 204);
    const auto *sfh0_206 = buffer.data(sfh0 + 206);
    const auto *sfh0_207 = buffer.data(sfh0 + 207);
    const auto *sfh0_209 = buffer.data(sfh0 + 209);

    const auto *sfg_45 = buffer.data(sfg + 45);
    const auto *sfg_47 = buffer.data(sfg + 47);
    const auto *sfg_48 = buffer.data(sfg + 48);
    const auto *sfg_50 = buffer.data(sfg + 50);
    const auto *sfg_55 = buffer.data(sfg + 55);
    const auto *sfg_59 = buffer.data(sfg + 59);
    const auto *sfg_60 = buffer.data(sfg + 60);
    const auto *sfg_62 = buffer.data(sfg + 62);
    const auto *sfg_63 = buffer.data(sfg + 63);
    const auto *sfg_65 = buffer.data(sfg + 65);
    const auto *sfg_70 = buffer.data(sfg + 70);
    const auto *sfg_74 = buffer.data(sfg + 74);
    const auto *sfg_75 = buffer.data(sfg + 75);
    const auto *sfg_77 = buffer.data(sfg + 77);
    const auto *sfg_78 = buffer.data(sfg + 78);
    const auto *sfg_80 = buffer.data(sfg + 80);
    const auto *sfg_85 = buffer.data(sfg + 85);
    const auto *sfg_89 = buffer.data(sfg + 89);
    const auto *sfg_90 = buffer.data(sfg + 90);
    const auto *sfg_92 = buffer.data(sfg + 92);
    const auto *sfg_93 = buffer.data(sfg + 93);
    const auto *sfg_95 = buffer.data(sfg + 95);
    const auto *sfg_96 = buffer.data(sfg + 96);
    const auto *sfg_99 = buffer.data(sfg + 99);
    const auto *sfg_100 = buffer.data(sfg + 100);
    const auto *sfg_101 = buffer.data(sfg + 101);
    const auto *sfg_102 = buffer.data(sfg + 102);
    const auto *sfg_103 = buffer.data(sfg + 103);
    const auto *sfg_104 = buffer.data(sfg + 104);
    const auto *sfg_105 = buffer.data(sfg + 105);
    const auto *sfg_107 = buffer.data(sfg + 107);
    const auto *sfg_110 = buffer.data(sfg + 110);
    const auto *sfg_114 = buffer.data(sfg + 114);
    const auto *sfg_115 = buffer.data(sfg + 115);
    const auto *sfg_116 = buffer.data(sfg + 116);
    const auto *sfg_117 = buffer.data(sfg + 117);
    const auto *sfg_118 = buffer.data(sfg + 118);
    const auto *sfg_119 = buffer.data(sfg + 119);
    const auto *sfg_120 = buffer.data(sfg + 120);
    const auto *sfg_123 = buffer.data(sfg + 123);
    const auto *sfg_126 = buffer.data(sfg + 126);
    const auto *sfg_130 = buffer.data(sfg + 130);
    const auto *sfg_131 = buffer.data(sfg + 131);
    const auto *sfg_132 = buffer.data(sfg + 132);
    const auto *sfg_133 = buffer.data(sfg + 133);
    const auto *sfg_134 = buffer.data(sfg + 134);
    const auto *sfg_135 = buffer.data(sfg + 135);
    const auto *sfg_138 = buffer.data(sfg + 138);
    const auto *sfg_140 = buffer.data(sfg + 140);
    const auto *sfg_141 = buffer.data(sfg + 141);
    const auto *sfg_144 = buffer.data(sfg + 144);
    const auto *sfg_145 = buffer.data(sfg + 145);
    const auto *sfg_146 = buffer.data(sfg + 146);
    const auto *sfg_147 = buffer.data(sfg + 147);
    const auto *sfg_148 = buffer.data(sfg + 148);
    const auto *sfg_149 = buffer.data(sfg + 149);

    const auto *sfh1_63 = buffer.data(sfh1 + 63);
    const auto *sfh1_66 = buffer.data(sfh1 + 66);
    const auto *sfh1_69 = buffer.data(sfh1 + 69);
    const auto *sfh1_105 = buffer.data(sfh1 + 105);
    const auto *sfh1_110 = buffer.data(sfh1 + 110);
    const auto *sfh1_114 = buffer.data(sfh1 + 114);
    const auto *sfh1_126 = buffer.data(sfh1 + 126);
    const auto *sfh1_129 = buffer.data(sfh1 + 129);
    const auto *sfh1_131 = buffer.data(sfh1 + 131);
    const auto *sfh1_132 = buffer.data(sfh1 + 132);
    const auto *sfh1_135 = buffer.data(sfh1 + 135);
    const auto *sfh1_141 = buffer.data(sfh1 + 141);
    const auto *sfh1_143 = buffer.data(sfh1 + 143);
    const auto *sfh1_144 = buffer.data(sfh1 + 144);
    const auto *sfh1_146 = buffer.data(sfh1 + 146);
    const auto *sfh1_152 = buffer.data(sfh1 + 152);
    const auto *sfh1_156 = buffer.data(sfh1 + 156);
    const auto *sfh1_162 = buffer.data(sfh1 + 162);
    const auto *sfh1_164 = buffer.data(sfh1 + 164);
    const auto *sfh1_165 = buffer.data(sfh1 + 165);
    const auto *sfh1_167 = buffer.data(sfh1 + 167);
    const auto *sfh1_171 = buffer.data(sfh1 + 171);
    const auto *sfh1_174 = buffer.data(sfh1 + 174);
    const auto *sfh1_183 = buffer.data(sfh1 + 183);
    const auto *sfh1_185 = buffer.data(sfh1 + 185);
    const auto *sfh1_186 = buffer.data(sfh1 + 186);
    const auto *sfh1_188 = buffer.data(sfh1 + 188);
    const auto *sfh1_189 = buffer.data(sfh1 + 189);
    const auto *sfh1_192 = buffer.data(sfh1 + 192);
    const auto *sfh1_194 = buffer.data(sfh1 + 194);
    const auto *sfh1_195 = buffer.data(sfh1 + 195);
    const auto *sfh1_198 = buffer.data(sfh1 + 198);
    const auto *sfh1_204 = buffer.data(sfh1 + 204);
    const auto *sfh1_206 = buffer.data(sfh1 + 206);
    const auto *sfh1_207 = buffer.data(sfh1 + 207);
    const auto *sfh1_209 = buffer.data(sfh1 + 209);

    const auto *sgf0_100 = buffer.data(sgf0 + 100);
    const auto *sgf0_103 = buffer.data(sgf0 + 103);
    const auto *sgf0_105 = buffer.data(sgf0 + 105);
    const auto *sgf0_106 = buffer.data(sgf0 + 106);
    const auto *sgf0_108 = buffer.data(sgf0 + 108);
    const auto *sgf0_109 = buffer.data(sgf0 + 109);
    const auto *sgf0_115 = buffer.data(sgf0 + 115);
    const auto *sgf0_119 = buffer.data(sgf0 + 119);
    const auto *sgf0_120 = buffer.data(sgf0 + 120);

    const auto *sgf1_100 = buffer.data(sgf1 + 100);
    const auto *sgf1_103 = buffer.data(sgf1 + 103);
    const auto *sgf1_105 = buffer.data(sgf1 + 105);
    const auto *sgf1_106 = buffer.data(sgf1 + 106);
    const auto *sgf1_108 = buffer.data(sgf1 + 108);
    const auto *sgf1_109 = buffer.data(sgf1 + 109);
    const auto *sgf1_115 = buffer.data(sgf1 + 115);
    const auto *sgf1_119 = buffer.data(sgf1 + 119);
    const auto *sgf1_120 = buffer.data(sgf1 + 120);

    const auto *sgg_90 = buffer.data(sgg + 90);
    const auto *sgg_92 = buffer.data(sgg + 92);
    const auto *sgg_93 = buffer.data(sgg + 93);
    const auto *sgg_95 = buffer.data(sgg + 95);
    const auto *sgg_100 = buffer.data(sgg + 100);
    const auto *sgg_101 = buffer.data(sgg + 101);
    const auto *sgg_102 = buffer.data(sgg + 102);
    const auto *sgg_103 = buffer.data(sgg + 103);
    const auto *sgg_104 = buffer.data(sgg + 104);
    const auto *sgg_105 = buffer.data(sgg + 105);
    const auto *sgg_107 = buffer.data(sgg + 107);
    const auto *sgg_108 = buffer.data(sgg + 108);
    const auto *sgg_110 = buffer.data(sgg + 110);
    const auto *sgg_115 = buffer.data(sgg + 115);
    const auto *sgg_116 = buffer.data(sgg + 116);
    const auto *sgg_117 = buffer.data(sgg + 117);
    const auto *sgg_118 = buffer.data(sgg + 118);
    const auto *sgg_119 = buffer.data(sgg + 119);
    const auto *sgg_120 = buffer.data(sgg + 120);
    const auto *sgg_122 = buffer.data(sgg + 122);
    const auto *sgg_123 = buffer.data(sgg + 123);
    const auto *sgg_125 = buffer.data(sgg + 125);
    const auto *sgg_130 = buffer.data(sgg + 130);
    const auto *sgg_131 = buffer.data(sgg + 131);
    const auto *sgg_132 = buffer.data(sgg + 132);
    const auto *sgg_133 = buffer.data(sgg + 133);
    const auto *sgg_134 = buffer.data(sgg + 134);
    const auto *sgg_135 = buffer.data(sgg + 135);
    const auto *sgg_137 = buffer.data(sgg + 137);
    const auto *sgg_138 = buffer.data(sgg + 138);
    const auto *sgg_140 = buffer.data(sgg + 140);
    const auto *sgg_145 = buffer.data(sgg + 145);
    const auto *sgg_146 = buffer.data(sgg + 146);
    const auto *sgg_147 = buffer.data(sgg + 147);
    const auto *sgg_148 = buffer.data(sgg + 148);
    const auto *sgg_149 = buffer.data(sgg + 149);
    const auto *sgg_150 = buffer.data(sgg + 150);
    const auto *sgg_152 = buffer.data(sgg + 152);
    const auto *sgg_153 = buffer.data(sgg + 153);
    const auto *sgg_155 = buffer.data(sgg + 155);
    const auto *sgg_156 = buffer.data(sgg + 156);
    const auto *sgg_159 = buffer.data(sgg + 159);
    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_161 = buffer.data(sgg + 161);
    const auto *sgg_162 = buffer.data(sgg + 162);
    const auto *sgg_163 = buffer.data(sgg + 163);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_165 = buffer.data(sgg + 165);
    const auto *sgg_167 = buffer.data(sgg + 167);
    const auto *sgg_168 = buffer.data(sgg + 168);
    const auto *sgg_170 = buffer.data(sgg + 170);
    const auto *sgg_174 = buffer.data(sgg + 174);
    const auto *sgg_175 = buffer.data(sgg + 175);
    const auto *sgg_176 = buffer.data(sgg + 176);
    const auto *sgg_177 = buffer.data(sgg + 177);
    const auto *sgg_178 = buffer.data(sgg + 178);
    const auto *sgg_179 = buffer.data(sgg + 179);
    const auto *sgg_180 = buffer.data(sgg + 180);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pb_x, pc_x, pc_y, pc_z, sfh0_126, \
                         sfh0_129, sfg_45, sfg_90, sfg_93, sfh1_126, sfh1_129, \
                         sgg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pb_x[k] * sfh0_126[k]
                   + f_12 * sfg_90[k]
                   - f_8 * pc_x[k] * sfh1_126[k];

        t_127[k] = f_11 * sfg_45[k]
                   + f_3 * pc_y[k] * sgg_90[k];

        t_128[k] = f_3 * pc_z[k] * sgg_90[k];

        t_129[k] = pb_x[k] * sfh0_129[k]
                   + f_11 * sfg_93[k]
                   - f_8 * pc_x[k] * sfh1_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, pc_x, pc_y, sfh0_131, sfh0_132, sfg_47, \
                         sfg_95, sfg_96, sfh1_131, sfh1_132, sgg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * sfg_47[k]
                   + f_3 * pc_y[k] * sgg_92[k];

        t_131[k] = pb_x[k] * sfh0_131[k]
                   + f_11 * sfg_95[k]
                   - f_8 * pc_x[k] * sfh1_131[k];

        t_132[k] = pb_x[k] * sfh0_132[k]
                   + f_10 * sfg_96[k]
                   - f_8 * pc_x[k] * sfh1_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pb_x, pc_x, pc_y, pc_z, sfh0_135, sfg_50, \
                         sfg_99, sfg_100, sfh1_135, sgg_93, sgg_95, \
                         sgg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * sgg_93[k];

        t_134[k] = f_11 * sfg_50[k]
                   + f_3 * pc_y[k] * sgg_95[k];

        t_135[k] = pb_x[k] * sfh0_135[k]
                   + f_10 * sfg_99[k]
                   - f_8 * pc_x[k] * sfh1_135[k];

        t_136[k] = f_9 * sfg_100[k]
                   + f_3 * pc_x[k] * sgg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, sfg_101, sfg_102, sfg_103, sfg_104, \
                         sgg_101, sgg_102, sgg_103, sgg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_9 * sfg_101[k]
                   + f_3 * pc_x[k] * sgg_101[k];

        t_138[k] = f_9 * sfg_102[k]
                   + f_3 * pc_x[k] * sgg_102[k];

        t_139[k] = f_9 * sfg_103[k]
                   + f_3 * pc_x[k] * sgg_103[k];

        t_140[k] = f_9 * sfg_104[k]
                   + f_3 * pc_x[k] * sgg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pb_x, pc_x, pc_z, sfh0_141, sfh0_143, \
                         sfh0_144, sfh1_141, sfh1_143, sfh1_144, \
                         sgg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = pb_x[k] * sfh0_141[k]
                   - f_8 * pc_x[k] * sfh1_141[k];

        t_142[k] = f_3 * pc_z[k] * sgg_100[k];

        t_143[k] = pb_x[k] * sfh0_143[k]
                   - f_8 * pc_x[k] * sfh1_143[k];

        t_144[k] = pb_x[k] * sfh0_144[k]
                   - f_8 * pc_x[k] * sfh1_144[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_z, pc_x, pc_y, pc_z, sfh0_63, sfh0_146, \
                         sfg_59, sfh1_63, sfh1_146, sgg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_11 * sfg_59[k]
                   + f_3 * pc_y[k] * sgg_104[k];

        t_146[k] = pb_x[k] * sfh0_146[k]
                   - f_8 * pc_x[k] * sfh1_146[k];

        t_147[k] = pb_z[k] * sfh0_63[k]
                   - f_8 * pc_z[k] * sfh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, sfh0_66, sfg_45, \
                         sfg_60, sfg_62, sfh1_66, sgg_105, sgg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * sfg_60[k]
                   + f_3 * pc_y[k] * sgg_105[k];

        t_149[k] = f_9 * sfg_45[k]
                   + f_3 * pc_z[k] * sgg_105[k];

        t_150[k] = pb_z[k] * sfh0_66[k]
                   - f_8 * pc_z[k] * sfh1_66[k];

        t_151[k] = f_10 * sfg_62[k]
                   + f_3 * pc_y[k] * sgg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_x, pb_z, pc_x, pc_z, sfh0_69, sfh0_152, \
                         sfg_48, sfg_110, sfh1_69, sfh1_152, sgg_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pb_x[k] * sfh0_152[k]
                   + f_11 * sfg_110[k]
                   - f_8 * pc_x[k] * sfh1_152[k];

        t_153[k] = pb_z[k] * sfh0_69[k]
                   - f_8 * pc_z[k] * sfh1_69[k];

        t_154[k] = f_9 * sfg_48[k]
                   + f_3 * pc_z[k] * sgg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pb_x, pc_x, pc_y, sfh0_156, sfg_65, \
                         sfg_114, sfg_115, sfg_116, sfh1_156, sgg_110, sgg_115, \
                         sgg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * sfg_65[k]
                   + f_3 * pc_y[k] * sgg_110[k];

        t_156[k] = pb_x[k] * sfh0_156[k]
                   + f_10 * sfg_114[k]
                   - f_8 * pc_x[k] * sfh1_156[k];

        t_157[k] = f_9 * sfg_115[k]
                   + f_3 * pc_x[k] * sgg_115[k];

        t_158[k] = f_9 * sfg_116[k]
                   + f_3 * pc_x[k] * sgg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_x, pc_x, sfh0_162, sfg_117, sfg_118, \
                         sfg_119, sfh1_162, sgg_117, sgg_118, sgg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_9 * sfg_117[k]
                   + f_3 * pc_x[k] * sgg_117[k];

        t_160[k] = f_9 * sfg_118[k]
                   + f_3 * pc_x[k] * sgg_118[k];

        t_161[k] = f_9 * sfg_119[k]
                   + f_3 * pc_x[k] * sgg_119[k];

        t_162[k] = pb_x[k] * sfh0_162[k]
                   - f_8 * pc_x[k] * sfh1_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pb_x, pc_x, pc_y, pc_z, sfh0_164, \
                         sfh0_165, sfg_55, sfg_74, sfh1_164, sfh1_165, sgg_115, \
                         sgg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * sfg_55[k]
                   + f_3 * pc_z[k] * sgg_115[k];

        t_164[k] = pb_x[k] * sfh0_164[k]
                   - f_8 * pc_x[k] * sfh1_164[k];

        t_165[k] = pb_x[k] * sfh0_165[k]
                   - f_8 * pc_x[k] * sfh1_165[k];

        t_166[k] = f_10 * sfg_74[k]
                   + f_3 * pc_y[k] * sgg_119[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pb_x, pb_y, pc_x, pc_y, pc_z, sfh0_105, \
                         sfh0_167, sfg_60, sfg_75, sfh1_105, sfh1_167, \
                         sgg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_x[k] * sfh0_167[k]
                   - f_8 * pc_x[k] * sfh1_167[k];

        t_168[k] = pb_y[k] * sfh0_105[k]
                   - f_8 * pc_y[k] * sfh1_105[k];

        t_169[k] = f_9 * sfg_75[k]
                   + f_3 * pc_y[k] * sgg_120[k];

        t_170[k] = f_10 * sfg_60[k]
                   + f_3 * pc_z[k] * sgg_120[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pb_x, pb_y, pc_x, pc_y, sfh0_110, sfh0_171, \
                         sfg_77, sfg_123, sfh1_110, sfh1_171, sgg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pb_x[k] * sfh0_171[k]
                   + f_11 * sfg_123[k]
                   - f_8 * pc_x[k] * sfh1_171[k];

        t_172[k] = f_9 * sfg_77[k]
                   + f_3 * pc_y[k] * sgg_122[k];

        t_173[k] = pb_y[k] * sfh0_110[k]
                   - f_8 * pc_y[k] * sfh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_x, pc_x, pc_y, pc_z, sfh0_174, sfg_63, \
                         sfg_80, sfg_126, sfh1_174, sgg_123, sgg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_x[k] * sfh0_174[k]
                   + f_10 * sfg_126[k]
                   - f_8 * pc_x[k] * sfh1_174[k];

        t_175[k] = f_10 * sfg_63[k]
                   + f_3 * pc_z[k] * sgg_123[k];

        t_176[k] = f_9 * sfg_80[k]
                   + f_3 * pc_y[k] * sgg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_x, pc_y, sfh0_114, sfg_130, \
                         sfg_131, sfg_132, sfh1_114, sgg_130, sgg_131, \
                         sgg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_y[k] * sfh0_114[k]
                   - f_8 * pc_y[k] * sfh1_114[k];

        t_178[k] = f_9 * sfg_130[k]
                   + f_3 * pc_x[k] * sgg_130[k];

        t_179[k] = f_9 * sfg_131[k]
                   + f_3 * pc_x[k] * sgg_131[k];

        t_180[k] = f_9 * sfg_132[k]
                   + f_3 * pc_x[k] * sgg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_x, pc_x, pc_z, sfh0_183, sfg_70, \
                         sfg_133, sfg_134, sfh1_183, sgg_130, sgg_133, \
                         sgg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_9 * sfg_133[k]
                   + f_3 * pc_x[k] * sgg_133[k];

        t_182[k] = f_9 * sfg_134[k]
                   + f_3 * pc_x[k] * sgg_134[k];

        t_183[k] = pb_x[k] * sfh0_183[k]
                   - f_8 * pc_x[k] * sfh1_183[k];

        t_184[k] = f_10 * sfg_70[k]
                   + f_3 * pc_z[k] * sgg_130[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, pc_x, pc_y, sfh0_185, sfh0_186, \
                         sfh0_188, sfg_89, sfh1_185, sfh1_186, sfh1_188, \
                         sgg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_x[k] * sfh0_185[k]
                   - f_8 * pc_x[k] * sfh1_185[k];

        t_186[k] = pb_x[k] * sfh0_186[k]
                   - f_8 * pc_x[k] * sfh1_186[k];

        t_187[k] = f_9 * sfg_89[k]
                   + f_3 * pc_y[k] * sgg_134[k];

        t_188[k] = pb_x[k] * sfh0_188[k]
                   - f_8 * pc_x[k] * sfh1_188[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, pc_x, pc_y, pc_z, sfh0_189, \
                         sfh0_192, sfg_75, sfg_135, sfg_138, sfh1_189, sfh1_192, \
                         sgg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_x[k] * sfh0_189[k]
                   + f_12 * sfg_135[k]
                   - f_8 * pc_x[k] * sfh1_189[k];

        t_190[k] = f_3 * pc_y[k] * sgg_135[k];

        t_191[k] = f_11 * sfg_75[k]
                   + f_3 * pc_z[k] * sgg_135[k];

        t_192[k] = pb_x[k] * sfh0_192[k]
                   + f_11 * sfg_138[k]
                   - f_8 * pc_x[k] * sfh1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pb_x, pc_x, pc_y, sfh0_194, sfh0_195, sfg_140, \
                         sfg_141, sfh1_194, sfh1_195, sgg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * sgg_137[k];

        t_194[k] = pb_x[k] * sfh0_194[k]
                   + f_11 * sfg_140[k]
                   - f_8 * pc_x[k] * sfh1_194[k];

        t_195[k] = pb_x[k] * sfh0_195[k]
                   + f_10 * sfg_141[k]
                   - f_8 * pc_x[k] * sfh1_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_x, pc_x, pc_y, pc_z, sfh0_198, sfg_78, \
                         sfg_144, sfg_145, sfh1_198, sgg_138, sgg_140, \
                         sgg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * sfg_78[k]
                   + f_3 * pc_z[k] * sgg_138[k];

        t_197[k] = f_3 * pc_y[k] * sgg_140[k];

        t_198[k] = pb_x[k] * sfh0_198[k]
                   + f_10 * sfg_144[k]
                   - f_8 * pc_x[k] * sfh1_198[k];

        t_199[k] = f_9 * sfg_145[k]
                   + f_3 * pc_x[k] * sgg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, sfg_146, sfg_147, sfg_148, sfg_149, \
                         sgg_146, sgg_147, sgg_148, sgg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_9 * sfg_146[k]
                   + f_3 * pc_x[k] * sgg_146[k];

        t_201[k] = f_9 * sfg_147[k]
                   + f_3 * pc_x[k] * sgg_147[k];

        t_202[k] = f_9 * sfg_148[k]
                   + f_3 * pc_x[k] * sgg_148[k];

        t_203[k] = f_9 * sfg_149[k]
                   + f_3 * pc_x[k] * sgg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_x, pc_x, pc_z, sfh0_204, sfh0_206, \
                         sfh0_207, sfg_85, sfh1_204, sfh1_206, sfh1_207, \
                         sgg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_x[k] * sfh0_204[k]
                   - f_8 * pc_x[k] * sfh1_204[k];

        t_205[k] = f_11 * sfg_85[k]
                   + f_3 * pc_z[k] * sgg_145[k];

        t_206[k] = pb_x[k] * sfh0_206[k]
                   - f_8 * pc_x[k] * sfh1_206[k];

        t_207[k] = pb_x[k] * sfh0_207[k]
                   - f_8 * pc_x[k] * sfh1_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pb_x, pc_x, pc_y, pc_z, sfh0_209, \
                         sfg_90, sfh1_209, sgf0_100, sgf1_100, sgg_149, \
                         sgg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * sgg_149[k];

        t_209[k] = pb_x[k] * sfh0_209[k]
                   - f_8 * pc_x[k] * sfh1_209[k];

        t_210[k] = f_1 * sgf0_100[k]
                   - f_2 * sgf1_100[k]
                   + f_3 * pc_x[k] * sgg_150[k];

        t_211[k] = f_0 * sfg_90[k]
                   + f_3 * pc_y[k] * sgg_150[k];

        t_212[k] = f_3 * pc_z[k] * sgg_150[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pc_x, pc_y, sfg_92, sgf0_103, sgf0_105, \
                         sgf1_103, sgf1_105, sgg_152, sgg_153, \
                         sgg_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_4 * sgf0_103[k]
                   - f_5 * sgf1_103[k]
                   + f_3 * pc_x[k] * sgg_153[k];

        t_214[k] = f_0 * sfg_92[k]
                   + f_3 * pc_y[k] * sgg_152[k];

        t_215[k] = f_4 * sgf0_105[k]
                   - f_5 * sgf1_105[k]
                   + f_3 * pc_x[k] * sgg_155[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, pc_x, pc_y, pc_z, sfg_95, sgf0_106, \
                         sgf0_109, sgf1_106, sgf1_109, sgg_153, sgg_155, sgg_156, \
                         sgg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * sgf0_106[k]
                   - f_7 * sgf1_106[k]
                   + f_3 * pc_x[k] * sgg_156[k];

        t_217[k] = f_3 * pc_z[k] * sgg_153[k];

        t_218[k] = f_0 * sfg_95[k]
                   + f_3 * pc_y[k] * sgg_155[k];

        t_219[k] = f_6 * sgf0_109[k]
                   - f_7 * sgf1_109[k]
                   + f_3 * pc_x[k] * sgg_159[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, t_225, pc_x, pc_y, sfg_100, \
                         sgf0_106, sgf1_106, sgg_160, sgg_161, sgg_162, sgg_163, \
                         sgg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_3 * pc_x[k] * sgg_160[k];

        t_221[k] = f_3 * pc_x[k] * sgg_161[k];

        t_222[k] = f_3 * pc_x[k] * sgg_162[k];

        t_223[k] = f_3 * pc_x[k] * sgg_163[k];

        t_224[k] = f_3 * pc_x[k] * sgg_164[k];

        t_225[k] = f_0 * sfg_100[k]
                   + f_1 * sgf0_106[k]
                   - f_2 * sgf1_106[k]
                   + f_3 * pc_y[k] * sgg_160[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pc_z, sfg_102, sfg_103, sgf0_108, \
                         sgf0_109, sgf1_108, sgf1_109, sgg_160, sgg_162, \
                         sgg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * sgg_160[k];

        t_227[k] = f_0 * sfg_102[k]
                   + f_4 * sgf0_108[k]
                   - f_5 * sgf1_108[k]
                   + f_3 * pc_y[k] * sgg_162[k];

        t_228[k] = f_0 * sfg_103[k]
                   + f_6 * sgf0_109[k]
                   - f_7 * sgf1_109[k]
                   + f_3 * pc_y[k] * sgg_163[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pb_z, pc_y, pc_z, sfh0_126, sfg_104, \
                         sfg_105, sfh1_126, sgf0_109, sgf1_109, sgg_164, \
                         sgg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_0 * sfg_104[k]
                   + f_3 * pc_y[k] * sgg_164[k];

        t_230[k] = f_1 * sgf0_109[k]
                   - f_2 * sgf1_109[k]
                   + f_3 * pc_z[k] * sgg_164[k];

        t_231[k] = pb_z[k] * sfh0_126[k]
                   - f_8 * pc_z[k] * sfh1_126[k];

        t_232[k] = f_11 * sfg_105[k]
                   + f_3 * pc_y[k] * sgg_165[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pb_z, pc_y, pc_z, sfh0_129, sfg_90, sfg_107, \
                         sfh1_129, sgg_165, sgg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_9 * sfg_90[k]
                   + f_3 * pc_z[k] * sgg_165[k];

        t_234[k] = pb_z[k] * sfh0_129[k]
                   - f_8 * pc_z[k] * sfh1_129[k];

        t_235[k] = f_11 * sfg_107[k]
                   + f_3 * pc_y[k] * sgg_167[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pb_z, pc_x, pc_y, pc_z, sfh0_132, sfg_93, \
                         sfg_110, sfh1_132, sgf0_115, sgf1_115, sgg_168, \
                         sgg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_4 * sgf0_115[k]
                   - f_5 * sgf1_115[k]
                   + f_3 * pc_x[k] * sgg_170[k];

        t_237[k] = pb_z[k] * sfh0_132[k]
                   - f_8 * pc_z[k] * sfh1_132[k];

        t_238[k] = f_9 * sfg_93[k]
                   + f_3 * pc_z[k] * sgg_168[k];

        t_239[k] = f_11 * sfg_110[k]
                   + f_3 * pc_y[k] * sgg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, t_245, pc_x, sgf0_119, sgf1_119, \
                         sgg_174, sgg_175, sgg_176, sgg_177, sgg_178, \
                         sgg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_6 * sgf0_119[k]
                   - f_7 * sgf1_119[k]
                   + f_3 * pc_x[k] * sgg_174[k];

        t_241[k] = f_3 * pc_x[k] * sgg_175[k];

        t_242[k] = f_3 * pc_x[k] * sgg_176[k];

        t_243[k] = f_3 * pc_x[k] * sgg_177[k];

        t_244[k] = f_3 * pc_x[k] * sgg_178[k];

        t_245[k] = f_3 * pc_x[k] * sgg_179[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pb_z, pc_z, sfh0_141, sfh0_143, sfh0_144, \
                         sfg_100, sfg_101, sfg_102, sfh1_141, sfh1_143, sfh1_144, \
                         sgg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pb_z[k] * sfh0_141[k]
                   - f_8 * pc_z[k] * sfh1_141[k];

        t_247[k] = f_9 * sfg_100[k]
                   + f_3 * pc_z[k] * sgg_175[k];

        t_248[k] = pb_z[k] * sfh0_143[k]
                   + f_10 * sfg_101[k]
                   - f_8 * pc_z[k] * sfh1_143[k];

        t_249[k] = pb_z[k] * sfh0_144[k]
                   + f_11 * sfg_102[k]
                   - f_8 * pc_z[k] * sfh1_144[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pc_x, pc_y, pc_z, sfg_104, sfg_119, \
                         sfg_120, sgf0_119, sgf0_120, sgf1_119, sgf1_120, sgg_179, \
                         sgg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_11 * sfg_119[k]
                   + f_3 * pc_y[k] * sgg_179[k];

        t_251[k] = f_9 * sfg_104[k]
                   + f_1 * sgf0_119[k]
                   - f_2 * sgf1_119[k]
                   + f_3 * pc_z[k] * sgg_179[k];

        t_252[k] = f_1 * sgf0_120[k]
                   - f_2 * sgf1_120[k]
                   + f_3 * pc_x[k] * sgg_180[k];

        t_253[k] = f_10 * sfg_120[k]
                   + f_3 * pc_y[k] * sgg_180[k];
    }
}

static auto
compute_prim_sgh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfh0,
                                                          const size_t sfg, const size_t sfh1,
                                                          const size_t sgf0, const size_t sgf1,
                                                          const size_t sgg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfh0_189 = buffer.data(sfh0 + 189);
    const auto *sfh0_194 = buffer.data(sfh0 + 194);
    const auto *sfh0_198 = buffer.data(sfh0 + 198);
    const auto *sfh0_204 = buffer.data(sfh0 + 204);
    const auto *sfh0_206 = buffer.data(sfh0 + 206);
    const auto *sfh0_207 = buffer.data(sfh0 + 207);
    const auto *sfh0_209 = buffer.data(sfh0 + 209);

    const auto *sfg_105 = buffer.data(sfg + 105);
    const auto *sfg_108 = buffer.data(sfg + 108);
    const auto *sfg_115 = buffer.data(sfg + 115);
    const auto *sfg_119 = buffer.data(sfg + 119);
    const auto *sfg_120 = buffer.data(sfg + 120);
    const auto *sfg_122 = buffer.data(sfg + 122);
    const auto *sfg_123 = buffer.data(sfg + 123);
    const auto *sfg_125 = buffer.data(sfg + 125);
    const auto *sfg_130 = buffer.data(sfg + 130);
    const auto *sfg_132 = buffer.data(sfg + 132);
    const auto *sfg_133 = buffer.data(sfg + 133);
    const auto *sfg_134 = buffer.data(sfg + 134);
    const auto *sfg_135 = buffer.data(sfg + 135);
    const auto *sfg_137 = buffer.data(sfg + 137);
    const auto *sfg_138 = buffer.data(sfg + 138);
    const auto *sfg_140 = buffer.data(sfg + 140);
    const auto *sfg_145 = buffer.data(sfg + 145);
    const auto *sfg_147 = buffer.data(sfg + 147);
    const auto *sfg_148 = buffer.data(sfg + 148);
    const auto *sfg_149 = buffer.data(sfg + 149);

    const auto *sfh1_189 = buffer.data(sfh1 + 189);
    const auto *sfh1_194 = buffer.data(sfh1 + 194);
    const auto *sfh1_198 = buffer.data(sfh1 + 198);
    const auto *sfh1_204 = buffer.data(sfh1 + 204);
    const auto *sfh1_206 = buffer.data(sfh1 + 206);
    const auto *sfh1_207 = buffer.data(sfh1 + 207);
    const auto *sfh1_209 = buffer.data(sfh1 + 209);

    const auto *sgf0_123 = buffer.data(sgf0 + 123);
    const auto *sgf0_125 = buffer.data(sgf0 + 125);
    const auto *sgf0_126 = buffer.data(sgf0 + 126);
    const auto *sgf0_128 = buffer.data(sgf0 + 128);
    const auto *sgf0_129 = buffer.data(sgf0 + 129);
    const auto *sgf0_133 = buffer.data(sgf0 + 133);
    const auto *sgf0_136 = buffer.data(sgf0 + 136);
    const auto *sgf0_140 = buffer.data(sgf0 + 140);
    const auto *sgf0_143 = buffer.data(sgf0 + 143);
    const auto *sgf0_145 = buffer.data(sgf0 + 145);
    const auto *sgf0_146 = buffer.data(sgf0 + 146);
    const auto *sgf0_148 = buffer.data(sgf0 + 148);
    const auto *sgf0_149 = buffer.data(sgf0 + 149);

    const auto *sgf1_123 = buffer.data(sgf1 + 123);
    const auto *sgf1_125 = buffer.data(sgf1 + 125);
    const auto *sgf1_126 = buffer.data(sgf1 + 126);
    const auto *sgf1_128 = buffer.data(sgf1 + 128);
    const auto *sgf1_129 = buffer.data(sgf1 + 129);
    const auto *sgf1_133 = buffer.data(sgf1 + 133);
    const auto *sgf1_136 = buffer.data(sgf1 + 136);
    const auto *sgf1_140 = buffer.data(sgf1 + 140);
    const auto *sgf1_143 = buffer.data(sgf1 + 143);
    const auto *sgf1_145 = buffer.data(sgf1 + 145);
    const auto *sgf1_146 = buffer.data(sgf1 + 146);
    const auto *sgf1_148 = buffer.data(sgf1 + 148);
    const auto *sgf1_149 = buffer.data(sgf1 + 149);

    const auto *sgg_180 = buffer.data(sgg + 180);
    const auto *sgg_182 = buffer.data(sgg + 182);
    const auto *sgg_183 = buffer.data(sgg + 183);
    const auto *sgg_185 = buffer.data(sgg + 185);
    const auto *sgg_186 = buffer.data(sgg + 186);
    const auto *sgg_189 = buffer.data(sgg + 189);
    const auto *sgg_190 = buffer.data(sgg + 190);
    const auto *sgg_191 = buffer.data(sgg + 191);
    const auto *sgg_192 = buffer.data(sgg + 192);
    const auto *sgg_193 = buffer.data(sgg + 193);
    const auto *sgg_194 = buffer.data(sgg + 194);
    const auto *sgg_195 = buffer.data(sgg + 195);
    const auto *sgg_197 = buffer.data(sgg + 197);
    const auto *sgg_198 = buffer.data(sgg + 198);
    const auto *sgg_200 = buffer.data(sgg + 200);
    const auto *sgg_201 = buffer.data(sgg + 201);
    const auto *sgg_205 = buffer.data(sgg + 205);
    const auto *sgg_206 = buffer.data(sgg + 206);
    const auto *sgg_207 = buffer.data(sgg + 207);
    const auto *sgg_208 = buffer.data(sgg + 208);
    const auto *sgg_209 = buffer.data(sgg + 209);
    const auto *sgg_210 = buffer.data(sgg + 210);
    const auto *sgg_212 = buffer.data(sgg + 212);
    const auto *sgg_213 = buffer.data(sgg + 213);
    const auto *sgg_215 = buffer.data(sgg + 215);
    const auto *sgg_216 = buffer.data(sgg + 216);
    const auto *sgg_219 = buffer.data(sgg + 219);
    const auto *sgg_220 = buffer.data(sgg + 220);
    const auto *sgg_221 = buffer.data(sgg + 221);
    const auto *sgg_222 = buffer.data(sgg + 222);
    const auto *sgg_223 = buffer.data(sgg + 223);
    const auto *sgg_224 = buffer.data(sgg + 224);

#pragma omp simd aligned(t_254, t_255, t_256, pc_x, pc_y, pc_z, sfg_105, sfg_122, sgf0_123, \
                         sgf1_123, sgg_180, sgg_182, sgg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * sfg_105[k]
                   + f_3 * pc_z[k] * sgg_180[k];

        t_255[k] = f_4 * sgf0_123[k]
                   - f_5 * sgf1_123[k]
                   + f_3 * pc_x[k] * sgg_183[k];

        t_256[k] = f_10 * sfg_122[k]
                   + f_3 * pc_y[k] * sgg_182[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, t_260, pc_x, pc_y, pc_z, sfg_108, sfg_125, \
                         sgf0_125, sgf0_126, sgf1_125, sgf1_126, sgg_183, sgg_185, \
                         sgg_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_4 * sgf0_125[k]
                   - f_5 * sgf1_125[k]
                   + f_3 * pc_x[k] * sgg_185[k];

        t_258[k] = f_6 * sgf0_126[k]
                   - f_7 * sgf1_126[k]
                   + f_3 * pc_x[k] * sgg_186[k];

        t_259[k] = f_10 * sfg_108[k]
                   + f_3 * pc_z[k] * sgg_183[k];

        t_260[k] = f_10 * sfg_125[k]
                   + f_3 * pc_y[k] * sgg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, t_265, t_266, pc_x, sgf0_129, sgf1_129, \
                         sgg_189, sgg_190, sgg_191, sgg_192, sgg_193, \
                         sgg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_6 * sgf0_129[k]
                   - f_7 * sgf1_129[k]
                   + f_3 * pc_x[k] * sgg_189[k];

        t_262[k] = f_3 * pc_x[k] * sgg_190[k];

        t_263[k] = f_3 * pc_x[k] * sgg_191[k];

        t_264[k] = f_3 * pc_x[k] * sgg_192[k];

        t_265[k] = f_3 * pc_x[k] * sgg_193[k];

        t_266[k] = f_3 * pc_x[k] * sgg_194[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, sfg_115, sfg_130, sfg_132, sgf0_126, \
                         sgf0_128, sgf1_126, sgf1_128, sgg_190, \
                         sgg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_10 * sfg_130[k]
                   + f_1 * sgf0_126[k]
                   - f_2 * sgf1_126[k]
                   + f_3 * pc_y[k] * sgg_190[k];

        t_268[k] = f_10 * sfg_115[k]
                   + f_3 * pc_z[k] * sgg_190[k];

        t_269[k] = f_10 * sfg_132[k]
                   + f_4 * sgf0_128[k]
                   - f_5 * sgf1_128[k]
                   + f_3 * pc_y[k] * sgg_192[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_y, pc_y, pc_z, sfh0_189, sfg_119, \
                         sfg_133, sfg_134, sfh1_189, sgf0_129, sgf1_129, sgg_193, \
                         sgg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_10 * sfg_133[k]
                   + f_6 * sgf0_129[k]
                   - f_7 * sgf1_129[k]
                   + f_3 * pc_y[k] * sgg_193[k];

        t_271[k] = f_10 * sfg_134[k]
                   + f_3 * pc_y[k] * sgg_194[k];

        t_272[k] = f_10 * sfg_119[k]
                   + f_1 * sgf0_129[k]
                   - f_2 * sgf1_129[k]
                   + f_3 * pc_z[k] * sgg_194[k];

        t_273[k] = pb_y[k] * sfh0_189[k]
                   - f_8 * pc_y[k] * sfh1_189[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pc_x, pc_y, pc_z, sfg_120, sfg_135, \
                         sfg_137, sgf0_133, sgf1_133, sgg_195, sgg_197, \
                         sgg_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_9 * sfg_135[k]
                   + f_3 * pc_y[k] * sgg_195[k];

        t_275[k] = f_11 * sfg_120[k]
                   + f_3 * pc_z[k] * sgg_195[k];

        t_276[k] = f_4 * sgf0_133[k]
                   - f_5 * sgf1_133[k]
                   + f_3 * pc_x[k] * sgg_198[k];

        t_277[k] = f_9 * sfg_137[k]
                   + f_3 * pc_y[k] * sgg_197[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pb_y, pc_x, pc_y, pc_z, sfh0_194, sfg_123, \
                         sfh1_194, sgf0_136, sgf1_136, sgg_198, \
                         sgg_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = pb_y[k] * sfh0_194[k]
                   - f_8 * pc_y[k] * sfh1_194[k];

        t_279[k] = f_6 * sgf0_136[k]
                   - f_7 * sgf1_136[k]
                   + f_3 * pc_x[k] * sgg_201[k];

        t_280[k] = f_11 * sfg_123[k]
                   + f_3 * pc_z[k] * sgg_198[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, t_285, pb_y, pc_x, pc_y, sfh0_198, \
                         sfg_140, sfh1_198, sgg_200, sgg_205, sgg_206, \
                         sgg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_9 * sfg_140[k]
                   + f_3 * pc_y[k] * sgg_200[k];

        t_282[k] = pb_y[k] * sfh0_198[k]
                   - f_8 * pc_y[k] * sfh1_198[k];

        t_283[k] = f_3 * pc_x[k] * sgg_205[k];

        t_284[k] = f_3 * pc_x[k] * sgg_206[k];

        t_285[k] = f_3 * pc_x[k] * sgg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_x, pc_y, pc_z, sfh0_204, \
                         sfg_130, sfg_145, sfh1_204, sgg_205, sgg_208, \
                         sgg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_3 * pc_x[k] * sgg_208[k];

        t_287[k] = f_3 * pc_x[k] * sgg_209[k];

        t_288[k] = pb_y[k] * sfh0_204[k]
                   + f_12 * sfg_145[k]
                   - f_8 * pc_y[k] * sfh1_204[k];

        t_289[k] = f_11 * sfg_130[k]
                   + f_3 * pc_z[k] * sgg_205[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_y, sfh0_206, sfh0_207, sfh0_209, \
                         sfg_147, sfg_148, sfg_149, sfh1_206, sfh1_207, sfh1_209, \
                         sgg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * sfh0_206[k]
                   + f_11 * sfg_147[k]
                   - f_8 * pc_y[k] * sfh1_206[k];

        t_291[k] = pb_y[k] * sfh0_207[k]
                   + f_10 * sfg_148[k]
                   - f_8 * pc_y[k] * sfh1_207[k];

        t_292[k] = f_9 * sfg_149[k]
                   + f_3 * pc_y[k] * sgg_209[k];

        t_293[k] = pb_y[k] * sfh0_209[k]
                   - f_8 * pc_y[k] * sfh1_209[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, t_298, pc_x, pc_y, pc_z, sfg_135, \
                         sgf0_140, sgf0_143, sgf1_140, sgf1_143, sgg_210, sgg_212, \
                         sgg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_1 * sgf0_140[k]
                   - f_2 * sgf1_140[k]
                   + f_3 * pc_x[k] * sgg_210[k];

        t_295[k] = f_3 * pc_y[k] * sgg_210[k];

        t_296[k] = f_0 * sfg_135[k]
                   + f_3 * pc_z[k] * sgg_210[k];

        t_297[k] = f_4 * sgf0_143[k]
                   - f_5 * sgf1_143[k]
                   + f_3 * pc_x[k] * sgg_213[k];

        t_298[k] = f_3 * pc_y[k] * sgg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, t_302, pc_x, pc_y, pc_z, sfg_138, sgf0_145, \
                         sgf0_146, sgf1_145, sgf1_146, sgg_213, sgg_215, \
                         sgg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_4 * sgf0_145[k]
                   - f_5 * sgf1_145[k]
                   + f_3 * pc_x[k] * sgg_215[k];

        t_300[k] = f_6 * sgf0_146[k]
                   - f_7 * sgf1_146[k]
                   + f_3 * pc_x[k] * sgg_216[k];

        t_301[k] = f_0 * sfg_138[k]
                   + f_3 * pc_z[k] * sgg_213[k];

        t_302[k] = f_3 * pc_y[k] * sgg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, t_307, t_308, pc_x, sgf0_149, sgf1_149, \
                         sgg_219, sgg_220, sgg_221, sgg_222, sgg_223, \
                         sgg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_6 * sgf0_149[k]
                   - f_7 * sgf1_149[k]
                   + f_3 * pc_x[k] * sgg_219[k];

        t_304[k] = f_3 * pc_x[k] * sgg_220[k];

        t_305[k] = f_3 * pc_x[k] * sgg_221[k];

        t_306[k] = f_3 * pc_x[k] * sgg_222[k];

        t_307[k] = f_3 * pc_x[k] * sgg_223[k];

        t_308[k] = f_3 * pc_x[k] * sgg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pc_y, pc_z, sfg_145, sgf0_146, sgf0_148, \
                         sgf0_149, sgf1_146, sgf1_148, sgf1_149, sgg_220, sgg_222, \
                         sgg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * sgf0_146[k]
                   - f_2 * sgf1_146[k]
                   + f_3 * pc_y[k] * sgg_220[k];

        t_310[k] = f_0 * sfg_145[k]
                   + f_3 * pc_z[k] * sgg_220[k];

        t_311[k] = f_4 * sgf0_148[k]
                   - f_5 * sgf1_148[k]
                   + f_3 * pc_y[k] * sgg_222[k];

        t_312[k] = f_6 * sgf0_149[k]
                   - f_7 * sgf1_149[k]
                   + f_3 * pc_y[k] * sgg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, pc_y, pc_z, sfg_149, sgf0_149, sgf1_149, \
                         sgg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * sgg_224[k];

        t_314[k] = f_0 * sfg_149[k]
                   + f_1 * sgf0_149[k]
                   - f_2 * sgf1_149[k]
                   + f_3 * pc_z[k] * sgg_224[k];
    }
}

auto
compute_prim_sgh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfh0, const size_t sfg,
                                                   const size_t sfh1, const size_t sgf0,
                                                   const size_t sgf1, const size_t sgg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfh0, sfg,
                                                              sfh1, sgf0, sgf1, sgg, ncols,
                                                              gamma, p, q);

    compute_prim_sgh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sfh0, sfg,
                                                              sfh1, sgf0, sgf1, sgg, ncols,
                                                              gamma, p, q);

    compute_prim_sgh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sfh0, sfg,
                                                              sfh1, sgf0, sgf1, sgg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
