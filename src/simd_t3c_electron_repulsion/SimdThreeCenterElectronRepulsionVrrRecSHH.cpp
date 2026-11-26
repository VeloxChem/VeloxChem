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


#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shh_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh0,
                                                          const size_t sgg, const size_t sgh1,
                                                          const size_t shf0, const size_t shf1,
                                                          const size_t shg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    auto *t_125 = buffer.data(target + 125);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_0 = buffer.data(sgh0 + 0);
    const auto *sgh0_3 = buffer.data(sgh0 + 3);
    const auto *sgh0_5 = buffer.data(sgh0 + 5);
    const auto *sgh0_6 = buffer.data(sgh0 + 6);
    const auto *sgh0_9 = buffer.data(sgh0 + 9);
    const auto *sgh0_15 = buffer.data(sgh0 + 15);
    const auto *sgh0_20 = buffer.data(sgh0 + 20);
    const auto *sgh0_24 = buffer.data(sgh0 + 24);
    const auto *sgh0_27 = buffer.data(sgh0 + 27);
    const auto *sgh0_36 = buffer.data(sgh0 + 36);
    const auto *sgh0_42 = buffer.data(sgh0 + 42);
    const auto *sgh0_47 = buffer.data(sgh0 + 47);
    const auto *sgh0_51 = buffer.data(sgh0 + 51);
    const auto *sgh0_62 = buffer.data(sgh0 + 62);

    const auto *sgg_0 = buffer.data(sgg + 0);
    const auto *sgg_1 = buffer.data(sgg + 1);
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
    const auto *sgg_48 = buffer.data(sgg + 48);
    const auto *sgg_50 = buffer.data(sgg + 50);
    const auto *sgg_51 = buffer.data(sgg + 51);
    const auto *sgg_54 = buffer.data(sgg + 54);
    const auto *sgg_55 = buffer.data(sgg + 55);
    const auto *sgg_56 = buffer.data(sgg + 56);
    const auto *sgg_57 = buffer.data(sgg + 57);
    const auto *sgg_58 = buffer.data(sgg + 58);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_70 = buffer.data(sgg + 70);
    const auto *sgg_71 = buffer.data(sgg + 71);
    const auto *sgg_72 = buffer.data(sgg + 72);
    const auto *sgg_73 = buffer.data(sgg + 73);
    const auto *sgg_74 = buffer.data(sgg + 74);
    const auto *sgg_75 = buffer.data(sgg + 75);
    const auto *sgg_78 = buffer.data(sgg + 78);
    const auto *sgg_80 = buffer.data(sgg + 80);
    const auto *sgg_81 = buffer.data(sgg + 81);
    const auto *sgg_84 = buffer.data(sgg + 84);
    const auto *sgg_85 = buffer.data(sgg + 85);
    const auto *sgg_86 = buffer.data(sgg + 86);
    const auto *sgg_87 = buffer.data(sgg + 87);
    const auto *sgg_88 = buffer.data(sgg + 88);
    const auto *sgg_89 = buffer.data(sgg + 89);

    const auto *sgh1_0 = buffer.data(sgh1 + 0);
    const auto *sgh1_3 = buffer.data(sgh1 + 3);
    const auto *sgh1_5 = buffer.data(sgh1 + 5);
    const auto *sgh1_6 = buffer.data(sgh1 + 6);
    const auto *sgh1_9 = buffer.data(sgh1 + 9);
    const auto *sgh1_15 = buffer.data(sgh1 + 15);
    const auto *sgh1_20 = buffer.data(sgh1 + 20);
    const auto *sgh1_24 = buffer.data(sgh1 + 24);
    const auto *sgh1_27 = buffer.data(sgh1 + 27);
    const auto *sgh1_36 = buffer.data(sgh1 + 36);
    const auto *sgh1_42 = buffer.data(sgh1 + 42);
    const auto *sgh1_47 = buffer.data(sgh1 + 47);
    const auto *sgh1_51 = buffer.data(sgh1 + 51);
    const auto *sgh1_62 = buffer.data(sgh1 + 62);

    const auto *shf0_0 = buffer.data(shf0 + 0);
    const auto *shf0_3 = buffer.data(shf0 + 3);
    const auto *shf0_5 = buffer.data(shf0 + 5);
    const auto *shf0_6 = buffer.data(shf0 + 6);
    const auto *shf0_8 = buffer.data(shf0 + 8);
    const auto *shf0_9 = buffer.data(shf0 + 9);
    const auto *shf0_16 = buffer.data(shf0 + 16);
    const auto *shf0_18 = buffer.data(shf0 + 18);
    const auto *shf0_19 = buffer.data(shf0 + 19);
    const auto *shf0_28 = buffer.data(shf0 + 28);
    const auto *shf0_29 = buffer.data(shf0 + 29);
    const auto *shf0_30 = buffer.data(shf0 + 30);
    const auto *shf0_33 = buffer.data(shf0 + 33);
    const auto *shf0_35 = buffer.data(shf0 + 35);
    const auto *shf0_36 = buffer.data(shf0 + 36);
    const auto *shf0_38 = buffer.data(shf0 + 38);
    const auto *shf0_39 = buffer.data(shf0 + 39);
    const auto *shf0_48 = buffer.data(shf0 + 48);
    const auto *shf0_49 = buffer.data(shf0 + 49);
    const auto *shf0_50 = buffer.data(shf0 + 50);
    const auto *shf0_53 = buffer.data(shf0 + 53);
    const auto *shf0_55 = buffer.data(shf0 + 55);
    const auto *shf0_56 = buffer.data(shf0 + 56);
    const auto *shf0_58 = buffer.data(shf0 + 58);
    const auto *shf0_59 = buffer.data(shf0 + 59);

    const auto *shf1_0 = buffer.data(shf1 + 0);
    const auto *shf1_3 = buffer.data(shf1 + 3);
    const auto *shf1_5 = buffer.data(shf1 + 5);
    const auto *shf1_6 = buffer.data(shf1 + 6);
    const auto *shf1_8 = buffer.data(shf1 + 8);
    const auto *shf1_9 = buffer.data(shf1 + 9);
    const auto *shf1_16 = buffer.data(shf1 + 16);
    const auto *shf1_18 = buffer.data(shf1 + 18);
    const auto *shf1_19 = buffer.data(shf1 + 19);
    const auto *shf1_28 = buffer.data(shf1 + 28);
    const auto *shf1_29 = buffer.data(shf1 + 29);
    const auto *shf1_30 = buffer.data(shf1 + 30);
    const auto *shf1_33 = buffer.data(shf1 + 33);
    const auto *shf1_35 = buffer.data(shf1 + 35);
    const auto *shf1_36 = buffer.data(shf1 + 36);
    const auto *shf1_38 = buffer.data(shf1 + 38);
    const auto *shf1_39 = buffer.data(shf1 + 39);
    const auto *shf1_48 = buffer.data(shf1 + 48);
    const auto *shf1_49 = buffer.data(shf1 + 49);
    const auto *shf1_50 = buffer.data(shf1 + 50);
    const auto *shf1_53 = buffer.data(shf1 + 53);
    const auto *shf1_55 = buffer.data(shf1 + 55);
    const auto *shf1_56 = buffer.data(shf1 + 56);
    const auto *shf1_58 = buffer.data(shf1 + 58);
    const auto *shf1_59 = buffer.data(shf1 + 59);

    const auto *shg_0 = buffer.data(shg + 0);
    const auto *shg_2 = buffer.data(shg + 2);
    const auto *shg_3 = buffer.data(shg + 3);
    const auto *shg_5 = buffer.data(shg + 5);
    const auto *shg_6 = buffer.data(shg + 6);
    const auto *shg_9 = buffer.data(shg + 9);
    const auto *shg_10 = buffer.data(shg + 10);
    const auto *shg_11 = buffer.data(shg + 11);
    const auto *shg_12 = buffer.data(shg + 12);
    const auto *shg_13 = buffer.data(shg + 13);
    const auto *shg_14 = buffer.data(shg + 14);
    const auto *shg_15 = buffer.data(shg + 15);
    const auto *shg_17 = buffer.data(shg + 17);
    const auto *shg_18 = buffer.data(shg + 18);
    const auto *shg_20 = buffer.data(shg + 20);
    const auto *shg_25 = buffer.data(shg + 25);
    const auto *shg_26 = buffer.data(shg + 26);
    const auto *shg_27 = buffer.data(shg + 27);
    const auto *shg_28 = buffer.data(shg + 28);
    const auto *shg_29 = buffer.data(shg + 29);
    const auto *shg_30 = buffer.data(shg + 30);
    const auto *shg_32 = buffer.data(shg + 32);
    const auto *shg_33 = buffer.data(shg + 33);
    const auto *shg_35 = buffer.data(shg + 35);
    const auto *shg_40 = buffer.data(shg + 40);
    const auto *shg_41 = buffer.data(shg + 41);
    const auto *shg_42 = buffer.data(shg + 42);
    const auto *shg_43 = buffer.data(shg + 43);
    const auto *shg_44 = buffer.data(shg + 44);
    const auto *shg_45 = buffer.data(shg + 45);
    const auto *shg_47 = buffer.data(shg + 47);
    const auto *shg_48 = buffer.data(shg + 48);
    const auto *shg_50 = buffer.data(shg + 50);
    const auto *shg_51 = buffer.data(shg + 51);
    const auto *shg_54 = buffer.data(shg + 54);
    const auto *shg_55 = buffer.data(shg + 55);
    const auto *shg_56 = buffer.data(shg + 56);
    const auto *shg_57 = buffer.data(shg + 57);
    const auto *shg_58 = buffer.data(shg + 58);
    const auto *shg_59 = buffer.data(shg + 59);
    const auto *shg_60 = buffer.data(shg + 60);
    const auto *shg_62 = buffer.data(shg + 62);
    const auto *shg_63 = buffer.data(shg + 63);
    const auto *shg_65 = buffer.data(shg + 65);
    const auto *shg_70 = buffer.data(shg + 70);
    const auto *shg_71 = buffer.data(shg + 71);
    const auto *shg_72 = buffer.data(shg + 72);
    const auto *shg_73 = buffer.data(shg + 73);
    const auto *shg_74 = buffer.data(shg + 74);
    const auto *shg_75 = buffer.data(shg + 75);
    const auto *shg_77 = buffer.data(shg + 77);
    const auto *shg_78 = buffer.data(shg + 78);
    const auto *shg_80 = buffer.data(shg + 80);
    const auto *shg_81 = buffer.data(shg + 81);
    const auto *shg_84 = buffer.data(shg + 84);
    const auto *shg_85 = buffer.data(shg + 85);
    const auto *shg_86 = buffer.data(shg + 86);
    const auto *shg_87 = buffer.data(shg + 87);
    const auto *shg_88 = buffer.data(shg + 88);
    const auto *shg_89 = buffer.data(shg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgg_0, sgg_3, shf0_0, shf0_3, \
                         shf1_0, shf1_3, shg_0, shg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgg_0[k]
                 + f_1 * shf0_0[k]
                 - f_2 * shf1_0[k]
                 + f_3 * pc_x[k] * shg_0[k];

        t_1[k] = f_3 * pc_y[k] * shg_0[k];

        t_2[k] = f_3 * pc_z[k] * shg_0[k];

        t_3[k] = f_0 * sgg_3[k]
                 + f_4 * shf0_3[k]
                 - f_5 * shf1_3[k]
                 + f_3 * pc_x[k] * shg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sgg_5, sgg_6, shf0_5, shf0_6, shf1_5, \
                         shf1_6, shg_2, shg_5, shg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * shg_2[k];

        t_5[k] = f_0 * sgg_5[k]
                 + f_4 * shf0_5[k]
                 - f_5 * shf1_5[k]
                 + f_3 * pc_x[k] * shg_5[k];

        t_6[k] = f_0 * sgg_6[k]
                 + f_6 * shf0_6[k]
                 - f_7 * shf1_6[k]
                 + f_3 * pc_x[k] * shg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, sgg_9, sgg_10, shf0_9, shf1_9, \
                         shg_3, shg_5, shg_9, shg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * shg_3[k];

        t_8[k] = f_3 * pc_y[k] * shg_5[k];

        t_9[k] = f_0 * sgg_9[k]
                 + f_6 * shf0_9[k]
                 - f_7 * shf1_9[k]
                 + f_3 * pc_x[k] * shg_9[k];

        t_10[k] = f_0 * sgg_10[k]
                  + f_3 * pc_x[k] * shg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, sgg_11, sgg_12, sgg_13, sgg_14, shg_11, \
                         shg_12, shg_13, shg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * sgg_11[k]
                  + f_3 * pc_x[k] * shg_11[k];

        t_12[k] = f_0 * sgg_12[k]
                  + f_3 * pc_x[k] * shg_12[k];

        t_13[k] = f_0 * sgg_13[k]
                  + f_3 * pc_x[k] * shg_13[k];

        t_14[k] = f_0 * sgg_14[k]
                  + f_3 * pc_x[k] * shg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, shf0_6, shf0_8, shf0_9, shf1_6, \
                         shf1_8, shf1_9, shg_10, shg_12, shg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * shf0_6[k]
                  - f_2 * shf1_6[k]
                  + f_3 * pc_y[k] * shg_10[k];

        t_16[k] = f_3 * pc_z[k] * shg_10[k];

        t_17[k] = f_4 * shf0_8[k]
                  - f_5 * shf1_8[k]
                  + f_3 * pc_y[k] * shg_12[k];

        t_18[k] = f_6 * shf0_9[k]
                  - f_7 * shf1_9[k]
                  + f_3 * pc_y[k] * shg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, sgh0_0, sgg_0, \
                         sgh1_0, shf0_9, shf1_9, shg_14, shg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * shg_14[k];

        t_20[k] = f_1 * shf0_9[k]
                  - f_2 * shf1_9[k]
                  + f_3 * pc_z[k] * shg_14[k];

        t_21[k] = pb_y[k] * sgh0_0[k]
                  - f_8 * pc_y[k] * sgh1_0[k];

        t_22[k] = f_9 * sgg_0[k]
                  + f_3 * pc_y[k] * shg_15[k];

        t_23[k] = f_3 * pc_z[k] * shg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_y, sgh0_3, sgh0_5, sgh0_6, sgg_1, \
                         sgg_2, sgg_3, sgh1_3, sgh1_5, sgh1_6, shg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sgh0_3[k]
                  + f_10 * sgg_1[k]
                  - f_8 * pc_y[k] * sgh1_3[k];

        t_25[k] = f_9 * sgg_2[k]
                  + f_3 * pc_y[k] * shg_17[k];

        t_26[k] = pb_y[k] * sgh0_5[k]
                  - f_8 * pc_y[k] * sgh1_5[k];

        t_27[k] = pb_y[k] * sgh0_6[k]
                  + f_11 * sgg_3[k]
                  - f_8 * pc_y[k] * sgh1_6[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_x, pc_y, pc_z, sgh0_9, sgg_5, \
                         sgg_25, sgh1_9, shg_18, shg_20, shg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_3 * pc_z[k] * shg_18[k];

        t_29[k] = f_9 * sgg_5[k]
                  + f_3 * pc_y[k] * shg_20[k];

        t_30[k] = pb_y[k] * sgh0_9[k]
                  - f_8 * pc_y[k] * sgh1_9[k];

        t_31[k] = f_12 * sgg_25[k]
                  + f_3 * pc_x[k] * shg_25[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, sgg_26, sgg_27, sgg_28, sgg_29, shg_26, \
                         shg_27, shg_28, shg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_12 * sgg_26[k]
                  + f_3 * pc_x[k] * shg_26[k];

        t_33[k] = f_12 * sgg_27[k]
                  + f_3 * pc_x[k] * shg_27[k];

        t_34[k] = f_12 * sgg_28[k]
                  + f_3 * pc_x[k] * shg_28[k];

        t_35[k] = f_12 * sgg_29[k]
                  + f_3 * pc_x[k] * shg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pc_y, pc_z, sgg_10, sgg_12, shf0_16, shf0_18, \
                         shf1_16, shf1_18, shg_25, shg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * sgg_10[k]
                  + f_1 * shf0_16[k]
                  - f_2 * shf1_16[k]
                  + f_3 * pc_y[k] * shg_25[k];

        t_37[k] = f_3 * pc_z[k] * shg_25[k];

        t_38[k] = f_9 * sgg_12[k]
                  + f_4 * shf0_18[k]
                  - f_5 * shf1_18[k]
                  + f_3 * pc_y[k] * shg_27[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, sgh0_20, sgg_13, sgg_14, sgh1_20, \
                         shf0_19, shf1_19, shg_28, shg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * sgg_13[k]
                  + f_6 * shf0_19[k]
                  - f_7 * shf1_19[k]
                  + f_3 * pc_y[k] * shg_28[k];

        t_40[k] = f_9 * sgg_14[k]
                  + f_3 * pc_y[k] * shg_29[k];

        t_41[k] = pb_y[k] * sgh0_20[k]
                  - f_8 * pc_y[k] * sgh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, sgh0_0, sgh0_3, \
                         sgg_0, sgh1_0, sgh1_3, shg_30, shg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sgh0_0[k]
                  - f_8 * pc_z[k] * sgh1_0[k];

        t_43[k] = f_3 * pc_y[k] * shg_30[k];

        t_44[k] = f_9 * sgg_0[k]
                  + f_3 * pc_z[k] * shg_30[k];

        t_45[k] = pb_z[k] * sgh0_3[k]
                  - f_8 * pc_z[k] * sgh1_3[k];

        t_46[k] = f_3 * pc_y[k] * shg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, sgh0_5, sgh0_6, sgg_2, \
                         sgg_3, sgh1_5, sgh1_6, shg_33, shg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * sgh0_5[k]
                  + f_10 * sgg_2[k]
                  - f_8 * pc_z[k] * sgh1_5[k];

        t_48[k] = pb_z[k] * sgh0_6[k]
                  - f_8 * pc_z[k] * sgh1_6[k];

        t_49[k] = f_9 * sgg_3[k]
                  + f_3 * pc_z[k] * shg_33[k];

        t_50[k] = f_3 * pc_y[k] * shg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, sgh0_9, sgg_5, sgg_40, \
                         sgg_41, sgg_42, sgh1_9, shg_40, shg_41, \
                         shg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * sgh0_9[k]
                  + f_11 * sgg_5[k]
                  - f_8 * pc_z[k] * sgh1_9[k];

        t_52[k] = f_12 * sgg_40[k]
                  + f_3 * pc_x[k] * shg_40[k];

        t_53[k] = f_12 * sgg_41[k]
                  + f_3 * pc_x[k] * shg_41[k];

        t_54[k] = f_12 * sgg_42[k]
                  + f_3 * pc_x[k] * shg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_z, sgh0_15, sgg_10, sgg_43, \
                         sgg_44, sgh1_15, shg_40, shg_43, shg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_12 * sgg_43[k]
                  + f_3 * pc_x[k] * shg_43[k];

        t_56[k] = f_12 * sgg_44[k]
                  + f_3 * pc_x[k] * shg_44[k];

        t_57[k] = pb_z[k] * sgh0_15[k]
                  - f_8 * pc_z[k] * sgh1_15[k];

        t_58[k] = f_9 * sgg_10[k]
                  + f_3 * pc_z[k] * shg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_y, pc_z, sgg_14, shf0_28, shf0_29, \
                         shf1_28, shf1_29, shg_42, shg_43, shg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_4 * shf0_28[k]
                  - f_5 * shf1_28[k]
                  + f_3 * pc_y[k] * shg_42[k];

        t_60[k] = f_6 * shf0_29[k]
                  - f_7 * shf1_29[k]
                  + f_3 * pc_y[k] * shg_43[k];

        t_61[k] = f_3 * pc_y[k] * shg_44[k];

        t_62[k] = f_9 * sgg_14[k]
                  + f_1 * shf0_29[k]
                  - f_2 * shf1_29[k]
                  + f_3 * pc_z[k] * shg_44[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, pc_z, sgg_15, sgg_45, sgg_48, \
                         shf0_30, shf0_33, shf1_30, shf1_33, shg_45, \
                         shg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_11 * sgg_45[k]
                  + f_1 * shf0_30[k]
                  - f_2 * shf1_30[k]
                  + f_3 * pc_x[k] * shg_45[k];

        t_64[k] = f_10 * sgg_15[k]
                  + f_3 * pc_y[k] * shg_45[k];

        t_65[k] = f_3 * pc_z[k] * shg_45[k];

        t_66[k] = f_11 * sgg_48[k]
                  + f_4 * shf0_33[k]
                  - f_5 * shf1_33[k]
                  + f_3 * pc_x[k] * shg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pc_x, pc_y, sgg_17, sgg_50, sgg_51, shf0_35, \
                         shf0_36, shf1_35, shf1_36, shg_47, shg_50, \
                         shg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_10 * sgg_17[k]
                  + f_3 * pc_y[k] * shg_47[k];

        t_68[k] = f_11 * sgg_50[k]
                  + f_4 * shf0_35[k]
                  - f_5 * shf1_35[k]
                  + f_3 * pc_x[k] * shg_50[k];

        t_69[k] = f_11 * sgg_51[k]
                  + f_6 * shf0_36[k]
                  - f_7 * shf1_36[k]
                  + f_3 * pc_x[k] * shg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pc_x, pc_y, pc_z, sgg_20, sgg_54, sgg_55, \
                         shf0_39, shf1_39, shg_48, shg_50, shg_54, \
                         shg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * shg_48[k];

        t_71[k] = f_10 * sgg_20[k]
                  + f_3 * pc_y[k] * shg_50[k];

        t_72[k] = f_11 * sgg_54[k]
                  + f_6 * shf0_39[k]
                  - f_7 * shf1_39[k]
                  + f_3 * pc_x[k] * shg_54[k];

        t_73[k] = f_11 * sgg_55[k]
                  + f_3 * pc_x[k] * shg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, sgg_56, sgg_57, sgg_58, sgg_59, shg_56, \
                         shg_57, shg_58, shg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * sgg_56[k]
                  + f_3 * pc_x[k] * shg_56[k];

        t_75[k] = f_11 * sgg_57[k]
                  + f_3 * pc_x[k] * shg_57[k];

        t_76[k] = f_11 * sgg_58[k]
                  + f_3 * pc_x[k] * shg_58[k];

        t_77[k] = f_11 * sgg_59[k]
                  + f_3 * pc_x[k] * shg_59[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, sgg_25, sgg_27, shf0_36, shf0_38, \
                         shf1_36, shf1_38, shg_55, shg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_10 * sgg_25[k]
                  + f_1 * shf0_36[k]
                  - f_2 * shf1_36[k]
                  + f_3 * pc_y[k] * shg_55[k];

        t_79[k] = f_3 * pc_z[k] * shg_55[k];

        t_80[k] = f_10 * sgg_27[k]
                  + f_4 * shf0_38[k]
                  - f_5 * shf1_38[k]
                  + f_3 * pc_y[k] * shg_57[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pb_y, pc_y, pc_z, sgh0_42, sgg_28, sgg_29, \
                         sgh1_42, shf0_39, shf1_39, shg_58, shg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * sgg_28[k]
                  + f_6 * shf0_39[k]
                  - f_7 * shf1_39[k]
                  + f_3 * pc_y[k] * shg_58[k];

        t_82[k] = f_10 * sgg_29[k]
                  + f_3 * pc_y[k] * shg_59[k];

        t_83[k] = f_1 * shf0_39[k]
                  - f_2 * shf1_39[k]
                  + f_3 * pc_z[k] * shg_59[k];

        t_84[k] = pb_y[k] * sgh0_42[k]
                  - f_8 * pc_y[k] * sgh1_42[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_z, pc_y, pc_z, sgh0_24, sgg_15, sgg_30, \
                         sgg_32, sgh1_24, shg_60, shg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_9 * sgg_30[k]
                  + f_3 * pc_y[k] * shg_60[k];

        t_86[k] = f_9 * sgg_15[k]
                  + f_3 * pc_z[k] * shg_60[k];

        t_87[k] = pb_z[k] * sgh0_24[k]
                  - f_8 * pc_z[k] * sgh1_24[k];

        t_88[k] = f_9 * sgg_32[k]
                  + f_3 * pc_y[k] * shg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_y, pb_z, pc_y, pc_z, sgh0_27, sgh0_47, \
                         sgg_18, sgg_35, sgh1_27, sgh1_47, shg_63, \
                         shg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pb_y[k] * sgh0_47[k]
                  - f_8 * pc_y[k] * sgh1_47[k];

        t_90[k] = pb_z[k] * sgh0_27[k]
                  - f_8 * pc_z[k] * sgh1_27[k];

        t_91[k] = f_9 * sgg_18[k]
                  + f_3 * pc_z[k] * shg_63[k];

        t_92[k] = f_9 * sgg_35[k]
                  + f_3 * pc_y[k] * shg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_y, pc_x, pc_y, sgh0_51, sgg_70, sgg_71, \
                         sgg_72, sgh1_51, shg_70, shg_71, shg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pb_y[k] * sgh0_51[k]
                  - f_8 * pc_y[k] * sgh1_51[k];

        t_94[k] = f_11 * sgg_70[k]
                  + f_3 * pc_x[k] * shg_70[k];

        t_95[k] = f_11 * sgg_71[k]
                  + f_3 * pc_x[k] * shg_71[k];

        t_96[k] = f_11 * sgg_72[k]
                  + f_3 * pc_x[k] * shg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, sgh0_36, sgg_25, sgg_73, \
                         sgg_74, sgh1_36, shg_70, shg_73, shg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_11 * sgg_73[k]
                  + f_3 * pc_x[k] * shg_73[k];

        t_98[k] = f_11 * sgg_74[k]
                  + f_3 * pc_x[k] * shg_74[k];

        t_99[k] = pb_z[k] * sgh0_36[k]
                  - f_8 * pc_z[k] * sgh1_36[k];

        t_100[k] = f_9 * sgg_25[k]
                   + f_3 * pc_z[k] * shg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, sgg_42, sgg_43, sgg_44, shf0_48, shf0_49, \
                         shf1_48, shf1_49, shg_72, shg_73, shg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_9 * sgg_42[k]
                   + f_4 * shf0_48[k]
                   - f_5 * shf1_48[k]
                   + f_3 * pc_y[k] * shg_72[k];

        t_102[k] = f_9 * sgg_43[k]
                   + f_6 * shf0_49[k]
                   - f_7 * shf1_49[k]
                   + f_3 * pc_y[k] * shg_73[k];

        t_103[k] = f_9 * sgg_44[k]
                   + f_3 * pc_y[k] * shg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, sgh0_62, sgg_30, \
                         sgg_75, sgh1_62, shf0_50, shf1_50, shg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * sgh0_62[k]
                   - f_8 * pc_y[k] * sgh1_62[k];

        t_105[k] = f_11 * sgg_75[k]
                   + f_1 * shf0_50[k]
                   - f_2 * shf1_50[k]
                   + f_3 * pc_x[k] * shg_75[k];

        t_106[k] = f_3 * pc_y[k] * shg_75[k];

        t_107[k] = f_10 * sgg_30[k]
                   + f_3 * pc_z[k] * shg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_x, pc_y, sgg_78, sgg_80, shf0_53, shf0_55, \
                         shf1_53, shf1_55, shg_77, shg_78, shg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_11 * sgg_78[k]
                   + f_4 * shf0_53[k]
                   - f_5 * shf1_53[k]
                   + f_3 * pc_x[k] * shg_78[k];

        t_109[k] = f_3 * pc_y[k] * shg_77[k];

        t_110[k] = f_11 * sgg_80[k]
                   + f_4 * shf0_55[k]
                   - f_5 * shf1_55[k]
                   + f_3 * pc_x[k] * shg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_x, pc_y, pc_z, sgg_33, sgg_81, shf0_56, \
                         shf1_56, shg_78, shg_80, shg_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_11 * sgg_81[k]
                   + f_6 * shf0_56[k]
                   - f_7 * shf1_56[k]
                   + f_3 * pc_x[k] * shg_81[k];

        t_112[k] = f_10 * sgg_33[k]
                   + f_3 * pc_z[k] * shg_78[k];

        t_113[k] = f_3 * pc_y[k] * shg_80[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pc_x, sgg_84, sgg_85, sgg_86, sgg_87, \
                         shf0_59, shf1_59, shg_84, shg_85, shg_86, \
                         shg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * sgg_84[k]
                   + f_6 * shf0_59[k]
                   - f_7 * shf1_59[k]
                   + f_3 * pc_x[k] * shg_84[k];

        t_115[k] = f_11 * sgg_85[k]
                   + f_3 * pc_x[k] * shg_85[k];

        t_116[k] = f_11 * sgg_86[k]
                   + f_3 * pc_x[k] * shg_86[k];

        t_117[k] = f_11 * sgg_87[k]
                   + f_3 * pc_x[k] * shg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pc_x, pc_y, pc_z, sgg_40, sgg_88, sgg_89, \
                         shf0_56, shf1_56, shg_85, shg_88, shg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_11 * sgg_88[k]
                   + f_3 * pc_x[k] * shg_88[k];

        t_119[k] = f_11 * sgg_89[k]
                   + f_3 * pc_x[k] * shg_89[k];

        t_120[k] = f_1 * shf0_56[k]
                   - f_2 * shf1_56[k]
                   + f_3 * pc_y[k] * shg_85[k];

        t_121[k] = f_10 * sgg_40[k]
                   + f_3 * pc_z[k] * shg_85[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_y, pc_z, sgg_44, shf0_58, shf0_59, \
                         shf1_58, shf1_59, shg_87, shg_88, shg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_4 * shf0_58[k]
                   - f_5 * shf1_58[k]
                   + f_3 * pc_y[k] * shg_87[k];

        t_123[k] = f_6 * shf0_59[k]
                   - f_7 * shf1_59[k]
                   + f_3 * pc_y[k] * shg_88[k];

        t_124[k] = f_3 * pc_y[k] * shg_89[k];

        t_125[k] = f_10 * sgg_44[k]
                   + f_1 * shf0_59[k]
                   - f_2 * shf1_59[k]
                   + f_3 * pc_z[k] * shg_89[k];
    }
}

static auto
compute_prim_shh_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh0,
                                                          const size_t sgg, const size_t sgh1,
                                                          const size_t shf0, const size_t shf1,
                                                          const size_t shg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_63 = buffer.data(sgh0 + 63);
    const auto *sgh0_66 = buffer.data(sgh0 + 66);
    const auto *sgh0_69 = buffer.data(sgh0 + 69);
    const auto *sgh0_78 = buffer.data(sgh0 + 78);
    const auto *sgh0_105 = buffer.data(sgh0 + 105);
    const auto *sgh0_108 = buffer.data(sgh0 + 108);
    const auto *sgh0_110 = buffer.data(sgh0 + 110);
    const auto *sgh0_111 = buffer.data(sgh0 + 111);
    const auto *sgh0_114 = buffer.data(sgh0 + 114);
    const auto *sgh0_125 = buffer.data(sgh0 + 125);
    const auto *sgh0_126 = buffer.data(sgh0 + 126);
    const auto *sgh0_129 = buffer.data(sgh0 + 129);
    const auto *sgh0_132 = buffer.data(sgh0 + 132);
    const auto *sgh0_210 = buffer.data(sgh0 + 210);
    const auto *sgh0_213 = buffer.data(sgh0 + 213);
    const auto *sgh0_215 = buffer.data(sgh0 + 215);
    const auto *sgh0_216 = buffer.data(sgh0 + 216);
    const auto *sgh0_219 = buffer.data(sgh0 + 219);
    const auto *sgh0_225 = buffer.data(sgh0 + 225);
    const auto *sgh0_227 = buffer.data(sgh0 + 227);
    const auto *sgh0_228 = buffer.data(sgh0 + 228);
    const auto *sgh0_230 = buffer.data(sgh0 + 230);
    const auto *sgh0_236 = buffer.data(sgh0 + 236);
    const auto *sgh0_240 = buffer.data(sgh0 + 240);
    const auto *sgh0_246 = buffer.data(sgh0 + 246);

    const auto *sgg_45 = buffer.data(sgg + 45);
    const auto *sgg_47 = buffer.data(sgg + 47);
    const auto *sgg_48 = buffer.data(sgg + 48);
    const auto *sgg_50 = buffer.data(sgg + 50);
    const auto *sgg_55 = buffer.data(sgg + 55);
    const auto *sgg_57 = buffer.data(sgg + 57);
    const auto *sgg_58 = buffer.data(sgg + 58);
    const auto *sgg_59 = buffer.data(sgg + 59);
    const auto *sgg_60 = buffer.data(sgg + 60);
    const auto *sgg_62 = buffer.data(sgg + 62);
    const auto *sgg_63 = buffer.data(sgg + 63);
    const auto *sgg_65 = buffer.data(sgg + 65);
    const auto *sgg_70 = buffer.data(sgg + 70);
    const auto *sgg_72 = buffer.data(sgg + 72);
    const auto *sgg_73 = buffer.data(sgg + 73);
    const auto *sgg_74 = buffer.data(sgg + 74);
    const auto *sgg_75 = buffer.data(sgg + 75);
    const auto *sgg_76 = buffer.data(sgg + 76);
    const auto *sgg_77 = buffer.data(sgg + 77);
    const auto *sgg_78 = buffer.data(sgg + 78);
    const auto *sgg_80 = buffer.data(sgg + 80);
    const auto *sgg_85 = buffer.data(sgg + 85);
    const auto *sgg_87 = buffer.data(sgg + 87);
    const auto *sgg_88 = buffer.data(sgg + 88);
    const auto *sgg_89 = buffer.data(sgg + 89);
    const auto *sgg_90 = buffer.data(sgg + 90);
    const auto *sgg_92 = buffer.data(sgg + 92);
    const auto *sgg_93 = buffer.data(sgg + 93);
    const auto *sgg_95 = buffer.data(sgg + 95);
    const auto *sgg_96 = buffer.data(sgg + 96);
    const auto *sgg_99 = buffer.data(sgg + 99);
    const auto *sgg_100 = buffer.data(sgg + 100);
    const auto *sgg_101 = buffer.data(sgg + 101);
    const auto *sgg_102 = buffer.data(sgg + 102);
    const auto *sgg_103 = buffer.data(sgg + 103);
    const auto *sgg_104 = buffer.data(sgg + 104);
    const auto *sgg_105 = buffer.data(sgg + 105);
    const auto *sgg_107 = buffer.data(sgg + 107);
    const auto *sgg_110 = buffer.data(sgg + 110);
    const auto *sgg_114 = buffer.data(sgg + 114);
    const auto *sgg_115 = buffer.data(sgg + 115);
    const auto *sgg_116 = buffer.data(sgg + 116);
    const auto *sgg_117 = buffer.data(sgg + 117);
    const auto *sgg_118 = buffer.data(sgg + 118);
    const auto *sgg_119 = buffer.data(sgg + 119);
    const auto *sgg_130 = buffer.data(sgg + 130);
    const auto *sgg_131 = buffer.data(sgg + 131);
    const auto *sgg_132 = buffer.data(sgg + 132);
    const auto *sgg_133 = buffer.data(sgg + 133);
    const auto *sgg_134 = buffer.data(sgg + 134);
    const auto *sgg_135 = buffer.data(sgg + 135);
    const auto *sgg_138 = buffer.data(sgg + 138);
    const auto *sgg_140 = buffer.data(sgg + 140);
    const auto *sgg_141 = buffer.data(sgg + 141);
    const auto *sgg_144 = buffer.data(sgg + 144);
    const auto *sgg_145 = buffer.data(sgg + 145);
    const auto *sgg_146 = buffer.data(sgg + 146);
    const auto *sgg_147 = buffer.data(sgg + 147);
    const auto *sgg_148 = buffer.data(sgg + 148);
    const auto *sgg_149 = buffer.data(sgg + 149);
    const auto *sgg_150 = buffer.data(sgg + 150);
    const auto *sgg_153 = buffer.data(sgg + 153);
    const auto *sgg_155 = buffer.data(sgg + 155);
    const auto *sgg_156 = buffer.data(sgg + 156);
    const auto *sgg_159 = buffer.data(sgg + 159);
    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_161 = buffer.data(sgg + 161);
    const auto *sgg_162 = buffer.data(sgg + 162);
    const auto *sgg_163 = buffer.data(sgg + 163);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_170 = buffer.data(sgg + 170);
    const auto *sgg_174 = buffer.data(sgg + 174);
    const auto *sgg_175 = buffer.data(sgg + 175);
    const auto *sgg_176 = buffer.data(sgg + 176);
    const auto *sgg_177 = buffer.data(sgg + 177);
    const auto *sgg_178 = buffer.data(sgg + 178);
    const auto *sgg_179 = buffer.data(sgg + 179);

    const auto *sgh1_63 = buffer.data(sgh1 + 63);
    const auto *sgh1_66 = buffer.data(sgh1 + 66);
    const auto *sgh1_69 = buffer.data(sgh1 + 69);
    const auto *sgh1_78 = buffer.data(sgh1 + 78);
    const auto *sgh1_105 = buffer.data(sgh1 + 105);
    const auto *sgh1_108 = buffer.data(sgh1 + 108);
    const auto *sgh1_110 = buffer.data(sgh1 + 110);
    const auto *sgh1_111 = buffer.data(sgh1 + 111);
    const auto *sgh1_114 = buffer.data(sgh1 + 114);
    const auto *sgh1_125 = buffer.data(sgh1 + 125);
    const auto *sgh1_126 = buffer.data(sgh1 + 126);
    const auto *sgh1_129 = buffer.data(sgh1 + 129);
    const auto *sgh1_132 = buffer.data(sgh1 + 132);
    const auto *sgh1_210 = buffer.data(sgh1 + 210);
    const auto *sgh1_213 = buffer.data(sgh1 + 213);
    const auto *sgh1_215 = buffer.data(sgh1 + 215);
    const auto *sgh1_216 = buffer.data(sgh1 + 216);
    const auto *sgh1_219 = buffer.data(sgh1 + 219);
    const auto *sgh1_225 = buffer.data(sgh1 + 225);
    const auto *sgh1_227 = buffer.data(sgh1 + 227);
    const auto *sgh1_228 = buffer.data(sgh1 + 228);
    const auto *sgh1_230 = buffer.data(sgh1 + 230);
    const auto *sgh1_236 = buffer.data(sgh1 + 236);
    const auto *sgh1_240 = buffer.data(sgh1 + 240);
    const auto *sgh1_246 = buffer.data(sgh1 + 246);

    const auto *shf0_60 = buffer.data(shf0 + 60);
    const auto *shf0_63 = buffer.data(shf0 + 63);
    const auto *shf0_65 = buffer.data(shf0 + 65);
    const auto *shf0_66 = buffer.data(shf0 + 66);
    const auto *shf0_68 = buffer.data(shf0 + 68);
    const auto *shf0_69 = buffer.data(shf0 + 69);
    const auto *shf0_75 = buffer.data(shf0 + 75);
    const auto *shf0_78 = buffer.data(shf0 + 78);
    const auto *shf0_79 = buffer.data(shf0 + 79);
    const auto *shf0_86 = buffer.data(shf0 + 86);
    const auto *shf0_88 = buffer.data(shf0 + 88);
    const auto *shf0_89 = buffer.data(shf0 + 89);
    const auto *shf0_90 = buffer.data(shf0 + 90);
    const auto *shf0_93 = buffer.data(shf0 + 93);
    const auto *shf0_95 = buffer.data(shf0 + 95);
    const auto *shf0_96 = buffer.data(shf0 + 96);
    const auto *shf0_98 = buffer.data(shf0 + 98);
    const auto *shf0_99 = buffer.data(shf0 + 99);

    const auto *shf1_60 = buffer.data(shf1 + 60);
    const auto *shf1_63 = buffer.data(shf1 + 63);
    const auto *shf1_65 = buffer.data(shf1 + 65);
    const auto *shf1_66 = buffer.data(shf1 + 66);
    const auto *shf1_68 = buffer.data(shf1 + 68);
    const auto *shf1_69 = buffer.data(shf1 + 69);
    const auto *shf1_75 = buffer.data(shf1 + 75);
    const auto *shf1_78 = buffer.data(shf1 + 78);
    const auto *shf1_79 = buffer.data(shf1 + 79);
    const auto *shf1_86 = buffer.data(shf1 + 86);
    const auto *shf1_88 = buffer.data(shf1 + 88);
    const auto *shf1_89 = buffer.data(shf1 + 89);
    const auto *shf1_90 = buffer.data(shf1 + 90);
    const auto *shf1_93 = buffer.data(shf1 + 93);
    const auto *shf1_95 = buffer.data(shf1 + 95);
    const auto *shf1_96 = buffer.data(shf1 + 96);
    const auto *shf1_98 = buffer.data(shf1 + 98);
    const auto *shf1_99 = buffer.data(shf1 + 99);

    const auto *shg_90 = buffer.data(shg + 90);
    const auto *shg_92 = buffer.data(shg + 92);
    const auto *shg_93 = buffer.data(shg + 93);
    const auto *shg_95 = buffer.data(shg + 95);
    const auto *shg_96 = buffer.data(shg + 96);
    const auto *shg_99 = buffer.data(shg + 99);
    const auto *shg_100 = buffer.data(shg + 100);
    const auto *shg_101 = buffer.data(shg + 101);
    const auto *shg_102 = buffer.data(shg + 102);
    const auto *shg_103 = buffer.data(shg + 103);
    const auto *shg_104 = buffer.data(shg + 104);
    const auto *shg_105 = buffer.data(shg + 105);
    const auto *shg_107 = buffer.data(shg + 107);
    const auto *shg_108 = buffer.data(shg + 108);
    const auto *shg_110 = buffer.data(shg + 110);
    const auto *shg_114 = buffer.data(shg + 114);
    const auto *shg_115 = buffer.data(shg + 115);
    const auto *shg_116 = buffer.data(shg + 116);
    const auto *shg_117 = buffer.data(shg + 117);
    const auto *shg_118 = buffer.data(shg + 118);
    const auto *shg_119 = buffer.data(shg + 119);
    const auto *shg_120 = buffer.data(shg + 120);
    const auto *shg_122 = buffer.data(shg + 122);
    const auto *shg_123 = buffer.data(shg + 123);
    const auto *shg_125 = buffer.data(shg + 125);
    const auto *shg_130 = buffer.data(shg + 130);
    const auto *shg_131 = buffer.data(shg + 131);
    const auto *shg_132 = buffer.data(shg + 132);
    const auto *shg_133 = buffer.data(shg + 133);
    const auto *shg_134 = buffer.data(shg + 134);
    const auto *shg_135 = buffer.data(shg + 135);
    const auto *shg_137 = buffer.data(shg + 137);
    const auto *shg_138 = buffer.data(shg + 138);
    const auto *shg_140 = buffer.data(shg + 140);
    const auto *shg_141 = buffer.data(shg + 141);
    const auto *shg_144 = buffer.data(shg + 144);
    const auto *shg_145 = buffer.data(shg + 145);
    const auto *shg_146 = buffer.data(shg + 146);
    const auto *shg_147 = buffer.data(shg + 147);
    const auto *shg_148 = buffer.data(shg + 148);
    const auto *shg_149 = buffer.data(shg + 149);
    const auto *shg_150 = buffer.data(shg + 150);
    const auto *shg_152 = buffer.data(shg + 152);
    const auto *shg_153 = buffer.data(shg + 153);
    const auto *shg_155 = buffer.data(shg + 155);
    const auto *shg_160 = buffer.data(shg + 160);
    const auto *shg_161 = buffer.data(shg + 161);
    const auto *shg_162 = buffer.data(shg + 162);
    const auto *shg_163 = buffer.data(shg + 163);
    const auto *shg_164 = buffer.data(shg + 164);
    const auto *shg_165 = buffer.data(shg + 165);
    const auto *shg_167 = buffer.data(shg + 167);
    const auto *shg_168 = buffer.data(shg + 168);
    const auto *shg_170 = buffer.data(shg + 170);
    const auto *shg_175 = buffer.data(shg + 175);
    const auto *shg_176 = buffer.data(shg + 176);
    const auto *shg_177 = buffer.data(shg + 177);
    const auto *shg_178 = buffer.data(shg + 178);
    const auto *shg_179 = buffer.data(shg + 179);

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, pc_y, pc_z, sgg_45, sgg_90, sgg_93, \
                         shf0_60, shf0_63, shf1_60, shf1_63, shg_90, \
                         shg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_10 * sgg_90[k]
                   + f_1 * shf0_60[k]
                   - f_2 * shf1_60[k]
                   + f_3 * pc_x[k] * shg_90[k];

        t_127[k] = f_11 * sgg_45[k]
                   + f_3 * pc_y[k] * shg_90[k];

        t_128[k] = f_3 * pc_z[k] * shg_90[k];

        t_129[k] = f_10 * sgg_93[k]
                   + f_4 * shf0_63[k]
                   - f_5 * shf1_63[k]
                   + f_3 * pc_x[k] * shg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pc_x, pc_y, sgg_47, sgg_95, sgg_96, shf0_65, \
                         shf0_66, shf1_65, shf1_66, shg_92, shg_95, \
                         shg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_11 * sgg_47[k]
                   + f_3 * pc_y[k] * shg_92[k];

        t_131[k] = f_10 * sgg_95[k]
                   + f_4 * shf0_65[k]
                   - f_5 * shf1_65[k]
                   + f_3 * pc_x[k] * shg_95[k];

        t_132[k] = f_10 * sgg_96[k]
                   + f_6 * shf0_66[k]
                   - f_7 * shf1_66[k]
                   + f_3 * pc_x[k] * shg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, pc_z, sgg_50, sgg_99, \
                         sgg_100, shf0_69, shf1_69, shg_93, shg_95, shg_99, \
                         shg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_3 * pc_z[k] * shg_93[k];

        t_134[k] = f_11 * sgg_50[k]
                   + f_3 * pc_y[k] * shg_95[k];

        t_135[k] = f_10 * sgg_99[k]
                   + f_6 * shf0_69[k]
                   - f_7 * shf1_69[k]
                   + f_3 * pc_x[k] * shg_99[k];

        t_136[k] = f_10 * sgg_100[k]
                   + f_3 * pc_x[k] * shg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, sgg_101, sgg_102, sgg_103, sgg_104, \
                         shg_101, shg_102, shg_103, shg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_10 * sgg_101[k]
                   + f_3 * pc_x[k] * shg_101[k];

        t_138[k] = f_10 * sgg_102[k]
                   + f_3 * pc_x[k] * shg_102[k];

        t_139[k] = f_10 * sgg_103[k]
                   + f_3 * pc_x[k] * shg_103[k];

        t_140[k] = f_10 * sgg_104[k]
                   + f_3 * pc_x[k] * shg_104[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, sgg_55, sgg_57, shf0_66, shf0_68, \
                         shf1_66, shf1_68, shg_100, shg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_11 * sgg_55[k]
                   + f_1 * shf0_66[k]
                   - f_2 * shf1_66[k]
                   + f_3 * pc_y[k] * shg_100[k];

        t_142[k] = f_3 * pc_z[k] * shg_100[k];

        t_143[k] = f_11 * sgg_57[k]
                   + f_4 * shf0_68[k]
                   - f_5 * shf1_68[k]
                   + f_3 * pc_y[k] * shg_102[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_z, pc_y, pc_z, sgh0_63, sgg_58, \
                         sgg_59, sgh1_63, shf0_69, shf1_69, shg_103, \
                         shg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_11 * sgg_58[k]
                   + f_6 * shf0_69[k]
                   - f_7 * shf1_69[k]
                   + f_3 * pc_y[k] * shg_103[k];

        t_145[k] = f_11 * sgg_59[k]
                   + f_3 * pc_y[k] * shg_104[k];

        t_146[k] = f_1 * shf0_69[k]
                   - f_2 * shf1_69[k]
                   + f_3 * pc_z[k] * shg_104[k];

        t_147[k] = pb_z[k] * sgh0_63[k]
                   - f_8 * pc_z[k] * sgh1_63[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_z, pc_y, pc_z, sgh0_66, sgg_45, \
                         sgg_60, sgg_62, sgh1_66, shg_105, shg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_10 * sgg_60[k]
                   + f_3 * pc_y[k] * shg_105[k];

        t_149[k] = f_9 * sgg_45[k]
                   + f_3 * pc_z[k] * shg_105[k];

        t_150[k] = pb_z[k] * sgh0_66[k]
                   - f_8 * pc_z[k] * sgh1_66[k];

        t_151[k] = f_10 * sgg_62[k]
                   + f_3 * pc_y[k] * shg_107[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pb_z, pc_x, pc_z, sgh0_69, sgg_48, sgg_110, \
                         sgh1_69, shf0_75, shf1_75, shg_108, shg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_10 * sgg_110[k]
                   + f_4 * shf0_75[k]
                   - f_5 * shf1_75[k]
                   + f_3 * pc_x[k] * shg_110[k];

        t_153[k] = pb_z[k] * sgh0_69[k]
                   - f_8 * pc_z[k] * sgh1_69[k];

        t_154[k] = f_9 * sgg_48[k]
                   + f_3 * pc_z[k] * shg_108[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_y, sgg_65, sgg_114, sgg_115, \
                         sgg_116, shf0_79, shf1_79, shg_110, shg_114, shg_115, \
                         shg_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_10 * sgg_65[k]
                   + f_3 * pc_y[k] * shg_110[k];

        t_156[k] = f_10 * sgg_114[k]
                   + f_6 * shf0_79[k]
                   - f_7 * shf1_79[k]
                   + f_3 * pc_x[k] * shg_114[k];

        t_157[k] = f_10 * sgg_115[k]
                   + f_3 * pc_x[k] * shg_115[k];

        t_158[k] = f_10 * sgg_116[k]
                   + f_3 * pc_x[k] * shg_116[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_z, pc_x, pc_z, sgh0_78, sgg_117, \
                         sgg_118, sgg_119, sgh1_78, shg_117, shg_118, \
                         shg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_10 * sgg_117[k]
                   + f_3 * pc_x[k] * shg_117[k];

        t_160[k] = f_10 * sgg_118[k]
                   + f_3 * pc_x[k] * shg_118[k];

        t_161[k] = f_10 * sgg_119[k]
                   + f_3 * pc_x[k] * shg_119[k];

        t_162[k] = pb_z[k] * sgh0_78[k]
                   - f_8 * pc_z[k] * sgh1_78[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pc_y, pc_z, sgg_55, sgg_72, sgg_73, shf0_78, \
                         shf0_79, shf1_78, shf1_79, shg_115, shg_117, \
                         shg_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * sgg_55[k]
                   + f_3 * pc_z[k] * shg_115[k];

        t_164[k] = f_10 * sgg_72[k]
                   + f_4 * shf0_78[k]
                   - f_5 * shf1_78[k]
                   + f_3 * pc_y[k] * shg_117[k];

        t_165[k] = f_10 * sgg_73[k]
                   + f_6 * shf0_79[k]
                   - f_7 * shf1_79[k]
                   + f_3 * pc_y[k] * shg_118[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_y, pc_y, pc_z, sgh0_105, sgg_59, \
                         sgg_74, sgg_75, sgh1_105, shf0_79, shf1_79, shg_119, \
                         shg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * sgg_74[k]
                   + f_3 * pc_y[k] * shg_119[k];

        t_167[k] = f_9 * sgg_59[k]
                   + f_1 * shf0_79[k]
                   - f_2 * shf1_79[k]
                   + f_3 * pc_z[k] * shg_119[k];

        t_168[k] = pb_y[k] * sgh0_105[k]
                   - f_8 * pc_y[k] * sgh1_105[k];

        t_169[k] = f_9 * sgg_75[k]
                   + f_3 * pc_y[k] * shg_120[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_y, pc_y, pc_z, sgh0_108, sgh0_110, \
                         sgg_60, sgg_76, sgg_77, sgh1_108, sgh1_110, shg_120, \
                         shg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_10 * sgg_60[k]
                   + f_3 * pc_z[k] * shg_120[k];

        t_171[k] = pb_y[k] * sgh0_108[k]
                   + f_10 * sgg_76[k]
                   - f_8 * pc_y[k] * sgh1_108[k];

        t_172[k] = f_9 * sgg_77[k]
                   + f_3 * pc_y[k] * shg_122[k];

        t_173[k] = pb_y[k] * sgh0_110[k]
                   - f_8 * pc_y[k] * sgh1_110[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_y, pc_y, pc_z, sgh0_111, sgh0_114, \
                         sgg_63, sgg_78, sgg_80, sgh1_111, sgh1_114, shg_123, \
                         shg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pb_y[k] * sgh0_111[k]
                   + f_11 * sgg_78[k]
                   - f_8 * pc_y[k] * sgh1_111[k];

        t_175[k] = f_10 * sgg_63[k]
                   + f_3 * pc_z[k] * shg_123[k];

        t_176[k] = f_9 * sgg_80[k]
                   + f_3 * pc_y[k] * shg_125[k];

        t_177[k] = pb_y[k] * sgh0_114[k]
                   - f_8 * pc_y[k] * sgh1_114[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, t_182, pc_x, sgg_130, sgg_131, sgg_132, \
                         sgg_133, sgg_134, shg_130, shg_131, shg_132, shg_133, \
                         shg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * sgg_130[k]
                   + f_3 * pc_x[k] * shg_130[k];

        t_179[k] = f_10 * sgg_131[k]
                   + f_3 * pc_x[k] * shg_131[k];

        t_180[k] = f_10 * sgg_132[k]
                   + f_3 * pc_x[k] * shg_132[k];

        t_181[k] = f_10 * sgg_133[k]
                   + f_3 * pc_x[k] * shg_133[k];

        t_182[k] = f_10 * sgg_134[k]
                   + f_3 * pc_x[k] * shg_134[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_y, pc_z, sgg_70, sgg_85, sgg_87, shf0_86, \
                         shf0_88, shf1_86, shf1_88, shg_130, shg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_9 * sgg_85[k]
                   + f_1 * shf0_86[k]
                   - f_2 * shf1_86[k]
                   + f_3 * pc_y[k] * shg_130[k];

        t_184[k] = f_10 * sgg_70[k]
                   + f_3 * pc_z[k] * shg_130[k];

        t_185[k] = f_9 * sgg_87[k]
                   + f_4 * shf0_88[k]
                   - f_5 * shf1_88[k]
                   + f_3 * pc_y[k] * shg_132[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pb_y, pc_y, sgh0_125, sgg_88, sgg_89, sgh1_125, \
                         shf0_89, shf1_89, shg_133, shg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_9 * sgg_88[k]
                   + f_6 * shf0_89[k]
                   - f_7 * shf1_89[k]
                   + f_3 * pc_y[k] * shg_133[k];

        t_187[k] = f_9 * sgg_89[k]
                   + f_3 * pc_y[k] * shg_134[k];

        t_188[k] = pb_y[k] * sgh0_125[k]
                   - f_8 * pc_y[k] * sgh1_125[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pc_x, pc_y, pc_z, sgg_75, sgg_135, \
                         sgg_138, shf0_90, shf0_93, shf1_90, shf1_93, shg_135, \
                         shg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_10 * sgg_135[k]
                   + f_1 * shf0_90[k]
                   - f_2 * shf1_90[k]
                   + f_3 * pc_x[k] * shg_135[k];

        t_190[k] = f_3 * pc_y[k] * shg_135[k];

        t_191[k] = f_11 * sgg_75[k]
                   + f_3 * pc_z[k] * shg_135[k];

        t_192[k] = f_10 * sgg_138[k]
                   + f_4 * shf0_93[k]
                   - f_5 * shf1_93[k]
                   + f_3 * pc_x[k] * shg_138[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pc_x, pc_y, sgg_140, sgg_141, shf0_95, shf0_96, \
                         shf1_95, shf1_96, shg_137, shg_140, shg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_3 * pc_y[k] * shg_137[k];

        t_194[k] = f_10 * sgg_140[k]
                   + f_4 * shf0_95[k]
                   - f_5 * shf1_95[k]
                   + f_3 * pc_x[k] * shg_140[k];

        t_195[k] = f_10 * sgg_141[k]
                   + f_6 * shf0_96[k]
                   - f_7 * shf1_96[k]
                   + f_3 * pc_x[k] * shg_141[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, sgg_78, sgg_144, \
                         sgg_145, shf0_99, shf1_99, shg_138, shg_140, shg_144, \
                         shg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_11 * sgg_78[k]
                   + f_3 * pc_z[k] * shg_138[k];

        t_197[k] = f_3 * pc_y[k] * shg_140[k];

        t_198[k] = f_10 * sgg_144[k]
                   + f_6 * shf0_99[k]
                   - f_7 * shf1_99[k]
                   + f_3 * pc_x[k] * shg_144[k];

        t_199[k] = f_10 * sgg_145[k]
                   + f_3 * pc_x[k] * shg_145[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, sgg_146, sgg_147, sgg_148, sgg_149, \
                         shg_146, shg_147, shg_148, shg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_10 * sgg_146[k]
                   + f_3 * pc_x[k] * shg_146[k];

        t_201[k] = f_10 * sgg_147[k]
                   + f_3 * pc_x[k] * shg_147[k];

        t_202[k] = f_10 * sgg_148[k]
                   + f_3 * pc_x[k] * shg_148[k];

        t_203[k] = f_10 * sgg_149[k]
                   + f_3 * pc_x[k] * shg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, sgg_85, shf0_96, shf0_98, \
                         shf0_99, shf1_96, shf1_98, shf1_99, shg_145, shg_147, \
                         shg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * shf0_96[k]
                   - f_2 * shf1_96[k]
                   + f_3 * pc_y[k] * shg_145[k];

        t_205[k] = f_11 * sgg_85[k]
                   + f_3 * pc_z[k] * shg_145[k];

        t_206[k] = f_4 * shf0_98[k]
                   - f_5 * shf1_98[k]
                   + f_3 * pc_y[k] * shg_147[k];

        t_207[k] = f_6 * shf0_99[k]
                   - f_7 * shf1_99[k]
                   + f_3 * pc_y[k] * shg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pb_x, pc_x, pc_y, pc_z, sgh0_210, sgg_89, \
                         sgg_150, sgh1_210, shf0_99, shf1_99, shg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_3 * pc_y[k] * shg_149[k];

        t_209[k] = f_11 * sgg_89[k]
                   + f_1 * shf0_99[k]
                   - f_2 * shf1_99[k]
                   + f_3 * pc_z[k] * shg_149[k];

        t_210[k] = pb_x[k] * sgh0_210[k]
                   + f_0 * sgg_150[k]
                   - f_8 * pc_x[k] * sgh1_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pb_x, pc_x, pc_y, pc_z, sgh0_213, sgg_90, \
                         sgg_92, sgg_153, sgh1_213, shg_150, shg_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_12 * sgg_90[k]
                   + f_3 * pc_y[k] * shg_150[k];

        t_212[k] = f_3 * pc_z[k] * shg_150[k];

        t_213[k] = pb_x[k] * sgh0_213[k]
                   + f_11 * sgg_153[k]
                   - f_8 * pc_x[k] * sgh1_213[k];

        t_214[k] = f_12 * sgg_92[k]
                   + f_3 * pc_y[k] * shg_152[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pb_x, pc_x, pc_z, sgh0_215, sgh0_216, sgg_155, \
                         sgg_156, sgh1_215, sgh1_216, shg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * sgh0_215[k]
                   + f_11 * sgg_155[k]
                   - f_8 * pc_x[k] * sgh1_215[k];

        t_216[k] = pb_x[k] * sgh0_216[k]
                   + f_10 * sgg_156[k]
                   - f_8 * pc_x[k] * sgh1_216[k];

        t_217[k] = f_3 * pc_z[k] * shg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pb_x, pc_x, pc_y, sgh0_219, sgg_95, \
                         sgg_159, sgg_160, sgg_161, sgh1_219, shg_155, shg_160, \
                         shg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_12 * sgg_95[k]
                   + f_3 * pc_y[k] * shg_155[k];

        t_219[k] = pb_x[k] * sgh0_219[k]
                   + f_10 * sgg_159[k]
                   - f_8 * pc_x[k] * sgh1_219[k];

        t_220[k] = f_9 * sgg_160[k]
                   + f_3 * pc_x[k] * shg_160[k];

        t_221[k] = f_9 * sgg_161[k]
                   + f_3 * pc_x[k] * shg_161[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pb_x, pc_x, sgh0_225, sgg_162, sgg_163, \
                         sgg_164, sgh1_225, shg_162, shg_163, shg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_9 * sgg_162[k]
                   + f_3 * pc_x[k] * shg_162[k];

        t_223[k] = f_9 * sgg_163[k]
                   + f_3 * pc_x[k] * shg_163[k];

        t_224[k] = f_9 * sgg_164[k]
                   + f_3 * pc_x[k] * shg_164[k];

        t_225[k] = pb_x[k] * sgh0_225[k]
                   - f_8 * pc_x[k] * sgh1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pb_x, pc_x, pc_y, pc_z, sgh0_227, \
                         sgh0_228, sgg_104, sgh1_227, sgh1_228, shg_160, \
                         shg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_3 * pc_z[k] * shg_160[k];

        t_227[k] = pb_x[k] * sgh0_227[k]
                   - f_8 * pc_x[k] * sgh1_227[k];

        t_228[k] = pb_x[k] * sgh0_228[k]
                   - f_8 * pc_x[k] * sgh1_228[k];

        t_229[k] = f_12 * sgg_104[k]
                   + f_3 * pc_y[k] * shg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pb_x, pb_z, pc_x, pc_y, pc_z, sgh0_126, \
                         sgh0_230, sgg_90, sgg_105, sgh1_126, sgh1_230, \
                         shg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = pb_x[k] * sgh0_230[k]
                   - f_8 * pc_x[k] * sgh1_230[k];

        t_231[k] = pb_z[k] * sgh0_126[k]
                   - f_8 * pc_z[k] * sgh1_126[k];

        t_232[k] = f_11 * sgg_105[k]
                   + f_3 * pc_y[k] * shg_165[k];

        t_233[k] = f_9 * sgg_90[k]
                   + f_3 * pc_z[k] * shg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pb_x, pb_z, pc_x, pc_y, pc_z, sgh0_129, \
                         sgh0_236, sgg_107, sgg_170, sgh1_129, sgh1_236, \
                         shg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pb_z[k] * sgh0_129[k]
                   - f_8 * pc_z[k] * sgh1_129[k];

        t_235[k] = f_11 * sgg_107[k]
                   + f_3 * pc_y[k] * shg_167[k];

        t_236[k] = pb_x[k] * sgh0_236[k]
                   + f_11 * sgg_170[k]
                   - f_8 * pc_x[k] * sgh1_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pb_z, pc_y, pc_z, sgh0_132, sgg_93, sgg_110, \
                         sgh1_132, shg_168, shg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pb_z[k] * sgh0_132[k]
                   - f_8 * pc_z[k] * sgh1_132[k];

        t_238[k] = f_9 * sgg_93[k]
                   + f_3 * pc_z[k] * shg_168[k];

        t_239[k] = f_11 * sgg_110[k]
                   + f_3 * pc_y[k] * shg_170[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_x, pc_x, sgh0_240, sgg_174, sgg_175, \
                         sgg_176, sgg_177, sgh1_240, shg_175, shg_176, \
                         shg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pb_x[k] * sgh0_240[k]
                   + f_10 * sgg_174[k]
                   - f_8 * pc_x[k] * sgh1_240[k];

        t_241[k] = f_9 * sgg_175[k]
                   + f_3 * pc_x[k] * shg_175[k];

        t_242[k] = f_9 * sgg_176[k]
                   + f_3 * pc_x[k] * shg_176[k];

        t_243[k] = f_9 * sgg_177[k]
                   + f_3 * pc_x[k] * shg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_x, pc_x, pc_z, sgh0_246, sgg_100, \
                         sgg_178, sgg_179, sgh1_246, shg_175, shg_178, \
                         shg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_9 * sgg_178[k]
                   + f_3 * pc_x[k] * shg_178[k];

        t_245[k] = f_9 * sgg_179[k]
                   + f_3 * pc_x[k] * shg_179[k];

        t_246[k] = pb_x[k] * sgh0_246[k]
                   - f_8 * pc_x[k] * sgh1_246[k];

        t_247[k] = f_9 * sgg_100[k]
                   + f_3 * pc_z[k] * shg_175[k];
    }
}

static auto
compute_prim_shh_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh0,
                                                          const size_t sgg, const size_t sgh1,
                                                          const size_t shf0, const size_t shf1,
                                                          const size_t shg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;

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
    auto *t_369 = buffer.data(target + 369);
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_189 = buffer.data(sgh0 + 189);
    const auto *sgh0_194 = buffer.data(sgh0 + 194);
    const auto *sgh0_198 = buffer.data(sgh0 + 198);
    const auto *sgh0_210 = buffer.data(sgh0 + 210);
    const auto *sgh0_213 = buffer.data(sgh0 + 213);
    const auto *sgh0_216 = buffer.data(sgh0 + 216);
    const auto *sgh0_225 = buffer.data(sgh0 + 225);
    const auto *sgh0_227 = buffer.data(sgh0 + 227);
    const auto *sgh0_228 = buffer.data(sgh0 + 228);
    const auto *sgh0_248 = buffer.data(sgh0 + 248);
    const auto *sgh0_249 = buffer.data(sgh0 + 249);
    const auto *sgh0_251 = buffer.data(sgh0 + 251);
    const auto *sgh0_252 = buffer.data(sgh0 + 252);
    const auto *sgh0_255 = buffer.data(sgh0 + 255);
    const auto *sgh0_257 = buffer.data(sgh0 + 257);
    const auto *sgh0_258 = buffer.data(sgh0 + 258);
    const auto *sgh0_261 = buffer.data(sgh0 + 261);
    const auto *sgh0_267 = buffer.data(sgh0 + 267);
    const auto *sgh0_269 = buffer.data(sgh0 + 269);
    const auto *sgh0_270 = buffer.data(sgh0 + 270);
    const auto *sgh0_272 = buffer.data(sgh0 + 272);
    const auto *sgh0_276 = buffer.data(sgh0 + 276);
    const auto *sgh0_279 = buffer.data(sgh0 + 279);
    const auto *sgh0_288 = buffer.data(sgh0 + 288);
    const auto *sgh0_290 = buffer.data(sgh0 + 290);
    const auto *sgh0_291 = buffer.data(sgh0 + 291);
    const auto *sgh0_293 = buffer.data(sgh0 + 293);
    const auto *sgh0_294 = buffer.data(sgh0 + 294);
    const auto *sgh0_297 = buffer.data(sgh0 + 297);
    const auto *sgh0_299 = buffer.data(sgh0 + 299);
    const auto *sgh0_300 = buffer.data(sgh0 + 300);
    const auto *sgh0_303 = buffer.data(sgh0 + 303);
    const auto *sgh0_309 = buffer.data(sgh0 + 309);
    const auto *sgh0_311 = buffer.data(sgh0 + 311);
    const auto *sgh0_312 = buffer.data(sgh0 + 312);
    const auto *sgh0_314 = buffer.data(sgh0 + 314);

    const auto *sgg_105 = buffer.data(sgg + 105);
    const auto *sgg_108 = buffer.data(sgg + 108);
    const auto *sgg_115 = buffer.data(sgg + 115);
    const auto *sgg_119 = buffer.data(sgg + 119);
    const auto *sgg_120 = buffer.data(sgg + 120);
    const auto *sgg_122 = buffer.data(sgg + 122);
    const auto *sgg_123 = buffer.data(sgg + 123);
    const auto *sgg_125 = buffer.data(sgg + 125);
    const auto *sgg_130 = buffer.data(sgg + 130);
    const auto *sgg_134 = buffer.data(sgg + 134);
    const auto *sgg_135 = buffer.data(sgg + 135);
    const auto *sgg_137 = buffer.data(sgg + 137);
    const auto *sgg_138 = buffer.data(sgg + 138);
    const auto *sgg_140 = buffer.data(sgg + 140);
    const auto *sgg_145 = buffer.data(sgg + 145);
    const auto *sgg_149 = buffer.data(sgg + 149);
    const auto *sgg_150 = buffer.data(sgg + 150);
    const auto *sgg_152 = buffer.data(sgg + 152);
    const auto *sgg_153 = buffer.data(sgg + 153);
    const auto *sgg_155 = buffer.data(sgg + 155);
    const auto *sgg_160 = buffer.data(sgg + 160);
    const auto *sgg_161 = buffer.data(sgg + 161);
    const auto *sgg_162 = buffer.data(sgg + 162);
    const auto *sgg_163 = buffer.data(sgg + 163);
    const auto *sgg_164 = buffer.data(sgg + 164);
    const auto *sgg_165 = buffer.data(sgg + 165);
    const auto *sgg_167 = buffer.data(sgg + 167);
    const auto *sgg_168 = buffer.data(sgg + 168);
    const auto *sgg_170 = buffer.data(sgg + 170);
    const auto *sgg_175 = buffer.data(sgg + 175);
    const auto *sgg_179 = buffer.data(sgg + 179);
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
    const auto *sgg_198 = buffer.data(sgg + 198);
    const auto *sgg_201 = buffer.data(sgg + 201);
    const auto *sgg_205 = buffer.data(sgg + 205);
    const auto *sgg_206 = buffer.data(sgg + 206);
    const auto *sgg_207 = buffer.data(sgg + 207);
    const auto *sgg_208 = buffer.data(sgg + 208);
    const auto *sgg_209 = buffer.data(sgg + 209);
    const auto *sgg_210 = buffer.data(sgg + 210);
    const auto *sgg_213 = buffer.data(sgg + 213);
    const auto *sgg_215 = buffer.data(sgg + 215);
    const auto *sgg_216 = buffer.data(sgg + 216);
    const auto *sgg_219 = buffer.data(sgg + 219);
    const auto *sgg_220 = buffer.data(sgg + 220);
    const auto *sgg_221 = buffer.data(sgg + 221);
    const auto *sgg_222 = buffer.data(sgg + 222);
    const auto *sgg_223 = buffer.data(sgg + 223);
    const auto *sgg_224 = buffer.data(sgg + 224);

    const auto *sgh1_189 = buffer.data(sgh1 + 189);
    const auto *sgh1_194 = buffer.data(sgh1 + 194);
    const auto *sgh1_198 = buffer.data(sgh1 + 198);
    const auto *sgh1_210 = buffer.data(sgh1 + 210);
    const auto *sgh1_213 = buffer.data(sgh1 + 213);
    const auto *sgh1_216 = buffer.data(sgh1 + 216);
    const auto *sgh1_225 = buffer.data(sgh1 + 225);
    const auto *sgh1_227 = buffer.data(sgh1 + 227);
    const auto *sgh1_228 = buffer.data(sgh1 + 228);
    const auto *sgh1_248 = buffer.data(sgh1 + 248);
    const auto *sgh1_249 = buffer.data(sgh1 + 249);
    const auto *sgh1_251 = buffer.data(sgh1 + 251);
    const auto *sgh1_252 = buffer.data(sgh1 + 252);
    const auto *sgh1_255 = buffer.data(sgh1 + 255);
    const auto *sgh1_257 = buffer.data(sgh1 + 257);
    const auto *sgh1_258 = buffer.data(sgh1 + 258);
    const auto *sgh1_261 = buffer.data(sgh1 + 261);
    const auto *sgh1_267 = buffer.data(sgh1 + 267);
    const auto *sgh1_269 = buffer.data(sgh1 + 269);
    const auto *sgh1_270 = buffer.data(sgh1 + 270);
    const auto *sgh1_272 = buffer.data(sgh1 + 272);
    const auto *sgh1_276 = buffer.data(sgh1 + 276);
    const auto *sgh1_279 = buffer.data(sgh1 + 279);
    const auto *sgh1_288 = buffer.data(sgh1 + 288);
    const auto *sgh1_290 = buffer.data(sgh1 + 290);
    const auto *sgh1_291 = buffer.data(sgh1 + 291);
    const auto *sgh1_293 = buffer.data(sgh1 + 293);
    const auto *sgh1_294 = buffer.data(sgh1 + 294);
    const auto *sgh1_297 = buffer.data(sgh1 + 297);
    const auto *sgh1_299 = buffer.data(sgh1 + 299);
    const auto *sgh1_300 = buffer.data(sgh1 + 300);
    const auto *sgh1_303 = buffer.data(sgh1 + 303);
    const auto *sgh1_309 = buffer.data(sgh1 + 309);
    const auto *sgh1_311 = buffer.data(sgh1 + 311);
    const auto *sgh1_312 = buffer.data(sgh1 + 312);
    const auto *sgh1_314 = buffer.data(sgh1 + 314);

    const auto *shf0_150 = buffer.data(shf0 + 150);
    const auto *shf0_153 = buffer.data(shf0 + 153);
    const auto *shf0_155 = buffer.data(shf0 + 155);
    const auto *shf0_156 = buffer.data(shf0 + 156);
    const auto *shf0_158 = buffer.data(shf0 + 158);
    const auto *shf0_159 = buffer.data(shf0 + 159);
    const auto *shf0_165 = buffer.data(shf0 + 165);
    const auto *shf0_169 = buffer.data(shf0 + 169);
    const auto *shf0_170 = buffer.data(shf0 + 170);
    const auto *shf0_173 = buffer.data(shf0 + 173);
    const auto *shf0_175 = buffer.data(shf0 + 175);
    const auto *shf0_176 = buffer.data(shf0 + 176);
    const auto *shf0_178 = buffer.data(shf0 + 178);
    const auto *shf0_179 = buffer.data(shf0 + 179);

    const auto *shf1_150 = buffer.data(shf1 + 150);
    const auto *shf1_153 = buffer.data(shf1 + 153);
    const auto *shf1_155 = buffer.data(shf1 + 155);
    const auto *shf1_156 = buffer.data(shf1 + 156);
    const auto *shf1_158 = buffer.data(shf1 + 158);
    const auto *shf1_159 = buffer.data(shf1 + 159);
    const auto *shf1_165 = buffer.data(shf1 + 165);
    const auto *shf1_169 = buffer.data(shf1 + 169);
    const auto *shf1_170 = buffer.data(shf1 + 170);
    const auto *shf1_173 = buffer.data(shf1 + 173);
    const auto *shf1_175 = buffer.data(shf1 + 175);
    const auto *shf1_176 = buffer.data(shf1 + 176);
    const auto *shf1_178 = buffer.data(shf1 + 178);
    const auto *shf1_179 = buffer.data(shf1 + 179);

    const auto *shg_179 = buffer.data(shg + 179);
    const auto *shg_180 = buffer.data(shg + 180);
    const auto *shg_182 = buffer.data(shg + 182);
    const auto *shg_183 = buffer.data(shg + 183);
    const auto *shg_185 = buffer.data(shg + 185);
    const auto *shg_190 = buffer.data(shg + 190);
    const auto *shg_191 = buffer.data(shg + 191);
    const auto *shg_192 = buffer.data(shg + 192);
    const auto *shg_193 = buffer.data(shg + 193);
    const auto *shg_194 = buffer.data(shg + 194);
    const auto *shg_195 = buffer.data(shg + 195);
    const auto *shg_197 = buffer.data(shg + 197);
    const auto *shg_198 = buffer.data(shg + 198);
    const auto *shg_200 = buffer.data(shg + 200);
    const auto *shg_205 = buffer.data(shg + 205);
    const auto *shg_206 = buffer.data(shg + 206);
    const auto *shg_207 = buffer.data(shg + 207);
    const auto *shg_208 = buffer.data(shg + 208);
    const auto *shg_209 = buffer.data(shg + 209);
    const auto *shg_210 = buffer.data(shg + 210);
    const auto *shg_212 = buffer.data(shg + 212);
    const auto *shg_213 = buffer.data(shg + 213);
    const auto *shg_215 = buffer.data(shg + 215);
    const auto *shg_220 = buffer.data(shg + 220);
    const auto *shg_221 = buffer.data(shg + 221);
    const auto *shg_222 = buffer.data(shg + 222);
    const auto *shg_223 = buffer.data(shg + 223);
    const auto *shg_224 = buffer.data(shg + 224);
    const auto *shg_225 = buffer.data(shg + 225);
    const auto *shg_227 = buffer.data(shg + 227);
    const auto *shg_228 = buffer.data(shg + 228);
    const auto *shg_230 = buffer.data(shg + 230);
    const auto *shg_231 = buffer.data(shg + 231);
    const auto *shg_234 = buffer.data(shg + 234);
    const auto *shg_235 = buffer.data(shg + 235);
    const auto *shg_236 = buffer.data(shg + 236);
    const auto *shg_237 = buffer.data(shg + 237);
    const auto *shg_238 = buffer.data(shg + 238);
    const auto *shg_239 = buffer.data(shg + 239);
    const auto *shg_240 = buffer.data(shg + 240);
    const auto *shg_242 = buffer.data(shg + 242);
    const auto *shg_243 = buffer.data(shg + 243);
    const auto *shg_245 = buffer.data(shg + 245);
    const auto *shg_249 = buffer.data(shg + 249);
    const auto *shg_250 = buffer.data(shg + 250);
    const auto *shg_251 = buffer.data(shg + 251);
    const auto *shg_252 = buffer.data(shg + 252);
    const auto *shg_253 = buffer.data(shg + 253);
    const auto *shg_254 = buffer.data(shg + 254);
    const auto *shg_255 = buffer.data(shg + 255);
    const auto *shg_257 = buffer.data(shg + 257);
    const auto *shg_258 = buffer.data(shg + 258);
    const auto *shg_260 = buffer.data(shg + 260);
    const auto *shg_261 = buffer.data(shg + 261);
    const auto *shg_264 = buffer.data(shg + 264);
    const auto *shg_265 = buffer.data(shg + 265);
    const auto *shg_266 = buffer.data(shg + 266);
    const auto *shg_267 = buffer.data(shg + 267);
    const auto *shg_268 = buffer.data(shg + 268);
    const auto *shg_269 = buffer.data(shg + 269);

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pb_x, pc_x, pc_y, sgh0_248, sgh0_249, \
                         sgh0_251, sgg_119, sgh1_248, sgh1_249, sgh1_251, \
                         shg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pb_x[k] * sgh0_248[k]
                   - f_8 * pc_x[k] * sgh1_248[k];

        t_249[k] = pb_x[k] * sgh0_249[k]
                   - f_8 * pc_x[k] * sgh1_249[k];

        t_250[k] = f_11 * sgg_119[k]
                   + f_3 * pc_y[k] * shg_179[k];

        t_251[k] = pb_x[k] * sgh0_251[k]
                   - f_8 * pc_x[k] * sgh1_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pb_x, pc_x, pc_y, pc_z, sgh0_252, sgg_105, \
                         sgg_120, sgg_180, sgh1_252, shg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pb_x[k] * sgh0_252[k]
                   + f_0 * sgg_180[k]
                   - f_8 * pc_x[k] * sgh1_252[k];

        t_253[k] = f_10 * sgg_120[k]
                   + f_3 * pc_y[k] * shg_180[k];

        t_254[k] = f_10 * sgg_105[k]
                   + f_3 * pc_z[k] * shg_180[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_x, pc_x, pc_y, sgh0_255, sgh0_257, sgg_122, \
                         sgg_183, sgg_185, sgh1_255, sgh1_257, \
                         shg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = pb_x[k] * sgh0_255[k]
                   + f_11 * sgg_183[k]
                   - f_8 * pc_x[k] * sgh1_255[k];

        t_256[k] = f_10 * sgg_122[k]
                   + f_3 * pc_y[k] * shg_182[k];

        t_257[k] = pb_x[k] * sgh0_257[k]
                   + f_11 * sgg_185[k]
                   - f_8 * pc_x[k] * sgh1_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pb_x, pc_x, pc_y, pc_z, sgh0_258, sgg_108, \
                         sgg_125, sgg_186, sgh1_258, shg_183, shg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pb_x[k] * sgh0_258[k]
                   + f_10 * sgg_186[k]
                   - f_8 * pc_x[k] * sgh1_258[k];

        t_259[k] = f_10 * sgg_108[k]
                   + f_3 * pc_z[k] * shg_183[k];

        t_260[k] = f_10 * sgg_125[k]
                   + f_3 * pc_y[k] * shg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_x, pc_x, sgh0_261, sgg_189, sgg_190, \
                         sgg_191, sgg_192, sgh1_261, shg_190, shg_191, \
                         shg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = pb_x[k] * sgh0_261[k]
                   + f_10 * sgg_189[k]
                   - f_8 * pc_x[k] * sgh1_261[k];

        t_262[k] = f_9 * sgg_190[k]
                   + f_3 * pc_x[k] * shg_190[k];

        t_263[k] = f_9 * sgg_191[k]
                   + f_3 * pc_x[k] * shg_191[k];

        t_264[k] = f_9 * sgg_192[k]
                   + f_3 * pc_x[k] * shg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, pb_x, pc_x, pc_z, sgh0_267, sgg_115, \
                         sgg_193, sgg_194, sgh1_267, shg_190, shg_193, \
                         shg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_9 * sgg_193[k]
                   + f_3 * pc_x[k] * shg_193[k];

        t_266[k] = f_9 * sgg_194[k]
                   + f_3 * pc_x[k] * shg_194[k];

        t_267[k] = pb_x[k] * sgh0_267[k]
                   - f_8 * pc_x[k] * sgh1_267[k];

        t_268[k] = f_10 * sgg_115[k]
                   + f_3 * pc_z[k] * shg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, t_272, pb_x, pc_x, pc_y, sgh0_269, sgh0_270, \
                         sgh0_272, sgg_134, sgh1_269, sgh1_270, sgh1_272, \
                         shg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = pb_x[k] * sgh0_269[k]
                   - f_8 * pc_x[k] * sgh1_269[k];

        t_270[k] = pb_x[k] * sgh0_270[k]
                   - f_8 * pc_x[k] * sgh1_270[k];

        t_271[k] = f_10 * sgg_134[k]
                   + f_3 * pc_y[k] * shg_194[k];

        t_272[k] = pb_x[k] * sgh0_272[k]
                   - f_8 * pc_x[k] * sgh1_272[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pb_y, pc_y, pc_z, sgh0_189, sgg_120, sgg_135, \
                         sgh1_189, shg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = pb_y[k] * sgh0_189[k]
                   - f_8 * pc_y[k] * sgh1_189[k];

        t_274[k] = f_9 * sgg_135[k]
                   + f_3 * pc_y[k] * shg_195[k];

        t_275[k] = f_11 * sgg_120[k]
                   + f_3 * pc_z[k] * shg_195[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pb_x, pb_y, pc_x, pc_y, sgh0_194, sgh0_276, \
                         sgg_137, sgg_198, sgh1_194, sgh1_276, \
                         shg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pb_x[k] * sgh0_276[k]
                   + f_11 * sgg_198[k]
                   - f_8 * pc_x[k] * sgh1_276[k];

        t_277[k] = f_9 * sgg_137[k]
                   + f_3 * pc_y[k] * shg_197[k];

        t_278[k] = pb_y[k] * sgh0_194[k]
                   - f_8 * pc_y[k] * sgh1_194[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pb_x, pc_x, pc_y, pc_z, sgh0_279, sgg_123, \
                         sgg_140, sgg_201, sgh1_279, shg_198, shg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = pb_x[k] * sgh0_279[k]
                   + f_10 * sgg_201[k]
                   - f_8 * pc_x[k] * sgh1_279[k];

        t_280[k] = f_11 * sgg_123[k]
                   + f_3 * pc_z[k] * shg_198[k];

        t_281[k] = f_9 * sgg_140[k]
                   + f_3 * pc_y[k] * shg_200[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_x, pc_y, sgh0_198, sgg_205, \
                         sgg_206, sgg_207, sgh1_198, shg_205, shg_206, \
                         shg_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pb_y[k] * sgh0_198[k]
                   - f_8 * pc_y[k] * sgh1_198[k];

        t_283[k] = f_9 * sgg_205[k]
                   + f_3 * pc_x[k] * shg_205[k];

        t_284[k] = f_9 * sgg_206[k]
                   + f_3 * pc_x[k] * shg_206[k];

        t_285[k] = f_9 * sgg_207[k]
                   + f_3 * pc_x[k] * shg_207[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_x, pc_x, pc_z, sgh0_288, sgg_130, \
                         sgg_208, sgg_209, sgh1_288, shg_205, shg_208, \
                         shg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_9 * sgg_208[k]
                   + f_3 * pc_x[k] * shg_208[k];

        t_287[k] = f_9 * sgg_209[k]
                   + f_3 * pc_x[k] * shg_209[k];

        t_288[k] = pb_x[k] * sgh0_288[k]
                   - f_8 * pc_x[k] * sgh1_288[k];

        t_289[k] = f_11 * sgg_130[k]
                   + f_3 * pc_z[k] * shg_205[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_x, pc_x, pc_y, sgh0_290, sgh0_291, \
                         sgh0_293, sgg_149, sgh1_290, sgh1_291, sgh1_293, \
                         shg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * sgh0_290[k]
                   - f_8 * pc_x[k] * sgh1_290[k];

        t_291[k] = pb_x[k] * sgh0_291[k]
                   - f_8 * pc_x[k] * sgh1_291[k];

        t_292[k] = f_9 * sgg_149[k]
                   + f_3 * pc_y[k] * shg_209[k];

        t_293[k] = pb_x[k] * sgh0_293[k]
                   - f_8 * pc_x[k] * sgh1_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pb_x, pc_x, pc_y, pc_z, sgh0_294, \
                         sgh0_297, sgg_135, sgg_210, sgg_213, sgh1_294, sgh1_297, \
                         shg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = pb_x[k] * sgh0_294[k]
                   + f_0 * sgg_210[k]
                   - f_8 * pc_x[k] * sgh1_294[k];

        t_295[k] = f_3 * pc_y[k] * shg_210[k];

        t_296[k] = f_12 * sgg_135[k]
                   + f_3 * pc_z[k] * shg_210[k];

        t_297[k] = pb_x[k] * sgh0_297[k]
                   + f_11 * sgg_213[k]
                   - f_8 * pc_x[k] * sgh1_297[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, pb_x, pc_x, pc_y, sgh0_299, sgh0_300, sgg_215, \
                         sgg_216, sgh1_299, sgh1_300, shg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_3 * pc_y[k] * shg_212[k];

        t_299[k] = pb_x[k] * sgh0_299[k]
                   + f_11 * sgg_215[k]
                   - f_8 * pc_x[k] * sgh1_299[k];

        t_300[k] = pb_x[k] * sgh0_300[k]
                   + f_10 * sgg_216[k]
                   - f_8 * pc_x[k] * sgh1_300[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_x, pc_x, pc_y, pc_z, sgh0_303, \
                         sgg_138, sgg_219, sgg_220, sgh1_303, shg_213, shg_215, \
                         shg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_12 * sgg_138[k]
                   + f_3 * pc_z[k] * shg_213[k];

        t_302[k] = f_3 * pc_y[k] * shg_215[k];

        t_303[k] = pb_x[k] * sgh0_303[k]
                   + f_10 * sgg_219[k]
                   - f_8 * pc_x[k] * sgh1_303[k];

        t_304[k] = f_9 * sgg_220[k]
                   + f_3 * pc_x[k] * shg_220[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, sgg_221, sgg_222, sgg_223, sgg_224, \
                         shg_221, shg_222, shg_223, shg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_9 * sgg_221[k]
                   + f_3 * pc_x[k] * shg_221[k];

        t_306[k] = f_9 * sgg_222[k]
                   + f_3 * pc_x[k] * shg_222[k];

        t_307[k] = f_9 * sgg_223[k]
                   + f_3 * pc_x[k] * shg_223[k];

        t_308[k] = f_9 * sgg_224[k]
                   + f_3 * pc_x[k] * shg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pb_x, pc_x, pc_z, sgh0_309, sgh0_311, \
                         sgh0_312, sgg_145, sgh1_309, sgh1_311, sgh1_312, \
                         shg_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pb_x[k] * sgh0_309[k]
                   - f_8 * pc_x[k] * sgh1_309[k];

        t_310[k] = f_12 * sgg_145[k]
                   + f_3 * pc_z[k] * shg_220[k];

        t_311[k] = pb_x[k] * sgh0_311[k]
                   - f_8 * pc_x[k] * sgh1_311[k];

        t_312[k] = pb_x[k] * sgh0_312[k]
                   - f_8 * pc_x[k] * sgh1_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, t_316, t_317, pb_x, pc_x, pc_y, pc_z, sgh0_314, \
                         sgg_150, sgh1_314, shf0_150, shf1_150, shg_224, \
                         shg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_3 * pc_y[k] * shg_224[k];

        t_314[k] = pb_x[k] * sgh0_314[k]
                   - f_8 * pc_x[k] * sgh1_314[k];

        t_315[k] = f_1 * shf0_150[k]
                   - f_2 * shf1_150[k]
                   + f_3 * pc_x[k] * shg_225[k];

        t_316[k] = f_0 * sgg_150[k]
                   + f_3 * pc_y[k] * shg_225[k];

        t_317[k] = f_3 * pc_z[k] * shg_225[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pc_x, pc_y, sgg_152, shf0_153, shf0_155, \
                         shf1_153, shf1_155, shg_227, shg_228, \
                         shg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = f_4 * shf0_153[k]
                   - f_5 * shf1_153[k]
                   + f_3 * pc_x[k] * shg_228[k];

        t_319[k] = f_0 * sgg_152[k]
                   + f_3 * pc_y[k] * shg_227[k];

        t_320[k] = f_4 * shf0_155[k]
                   - f_5 * shf1_155[k]
                   + f_3 * pc_x[k] * shg_230[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pc_y, pc_z, sgg_155, shf0_156, \
                         shf0_159, shf1_156, shf1_159, shg_228, shg_230, shg_231, \
                         shg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_6 * shf0_156[k]
                   - f_7 * shf1_156[k]
                   + f_3 * pc_x[k] * shg_231[k];

        t_322[k] = f_3 * pc_z[k] * shg_228[k];

        t_323[k] = f_0 * sgg_155[k]
                   + f_3 * pc_y[k] * shg_230[k];

        t_324[k] = f_6 * shf0_159[k]
                   - f_7 * shf1_159[k]
                   + f_3 * pc_x[k] * shg_234[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, t_330, pc_x, pc_y, sgg_160, \
                         shf0_156, shf1_156, shg_235, shg_236, shg_237, shg_238, \
                         shg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_3 * pc_x[k] * shg_235[k];

        t_326[k] = f_3 * pc_x[k] * shg_236[k];

        t_327[k] = f_3 * pc_x[k] * shg_237[k];

        t_328[k] = f_3 * pc_x[k] * shg_238[k];

        t_329[k] = f_3 * pc_x[k] * shg_239[k];

        t_330[k] = f_0 * sgg_160[k]
                   + f_1 * shf0_156[k]
                   - f_2 * shf1_156[k]
                   + f_3 * pc_y[k] * shg_235[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pc_y, pc_z, sgg_162, sgg_163, shf0_158, \
                         shf0_159, shf1_158, shf1_159, shg_235, shg_237, \
                         shg_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_3 * pc_z[k] * shg_235[k];

        t_332[k] = f_0 * sgg_162[k]
                   + f_4 * shf0_158[k]
                   - f_5 * shf1_158[k]
                   + f_3 * pc_y[k] * shg_237[k];

        t_333[k] = f_0 * sgg_163[k]
                   + f_6 * shf0_159[k]
                   - f_7 * shf1_159[k]
                   + f_3 * pc_y[k] * shg_238[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pb_z, pc_y, pc_z, sgh0_210, sgg_164, \
                         sgg_165, sgh1_210, shf0_159, shf1_159, shg_239, \
                         shg_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_0 * sgg_164[k]
                   + f_3 * pc_y[k] * shg_239[k];

        t_335[k] = f_1 * shf0_159[k]
                   - f_2 * shf1_159[k]
                   + f_3 * pc_z[k] * shg_239[k];

        t_336[k] = pb_z[k] * sgh0_210[k]
                   - f_8 * pc_z[k] * sgh1_210[k];

        t_337[k] = f_12 * sgg_165[k]
                   + f_3 * pc_y[k] * shg_240[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_y, pc_z, sgh0_213, sgg_150, sgg_167, \
                         sgh1_213, shg_240, shg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_9 * sgg_150[k]
                   + f_3 * pc_z[k] * shg_240[k];

        t_339[k] = pb_z[k] * sgh0_213[k]
                   - f_8 * pc_z[k] * sgh1_213[k];

        t_340[k] = f_12 * sgg_167[k]
                   + f_3 * pc_y[k] * shg_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pb_z, pc_x, pc_y, pc_z, sgh0_216, \
                         sgg_153, sgg_170, sgh1_216, shf0_165, shf1_165, shg_243, \
                         shg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_4 * shf0_165[k]
                   - f_5 * shf1_165[k]
                   + f_3 * pc_x[k] * shg_245[k];

        t_342[k] = pb_z[k] * sgh0_216[k]
                   - f_8 * pc_z[k] * sgh1_216[k];

        t_343[k] = f_9 * sgg_153[k]
                   + f_3 * pc_z[k] * shg_243[k];

        t_344[k] = f_12 * sgg_170[k]
                   + f_3 * pc_y[k] * shg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, t_350, pc_x, shf0_169, shf1_169, \
                         shg_249, shg_250, shg_251, shg_252, shg_253, \
                         shg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_6 * shf0_169[k]
                   - f_7 * shf1_169[k]
                   + f_3 * pc_x[k] * shg_249[k];

        t_346[k] = f_3 * pc_x[k] * shg_250[k];

        t_347[k] = f_3 * pc_x[k] * shg_251[k];

        t_348[k] = f_3 * pc_x[k] * shg_252[k];

        t_349[k] = f_3 * pc_x[k] * shg_253[k];

        t_350[k] = f_3 * pc_x[k] * shg_254[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pb_z, pc_z, sgh0_225, sgh0_227, sgh0_228, \
                         sgg_160, sgg_161, sgg_162, sgh1_225, sgh1_227, sgh1_228, \
                         shg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = pb_z[k] * sgh0_225[k]
                   - f_8 * pc_z[k] * sgh1_225[k];

        t_352[k] = f_9 * sgg_160[k]
                   + f_3 * pc_z[k] * shg_250[k];

        t_353[k] = pb_z[k] * sgh0_227[k]
                   + f_10 * sgg_161[k]
                   - f_8 * pc_z[k] * sgh1_227[k];

        t_354[k] = pb_z[k] * sgh0_228[k]
                   + f_11 * sgg_162[k]
                   - f_8 * pc_z[k] * sgh1_228[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, pc_z, sgg_164, sgg_179, \
                         sgg_180, shf0_169, shf0_170, shf1_169, shf1_170, shg_254, \
                         shg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_12 * sgg_179[k]
                   + f_3 * pc_y[k] * shg_254[k];

        t_356[k] = f_9 * sgg_164[k]
                   + f_1 * shf0_169[k]
                   - f_2 * shf1_169[k]
                   + f_3 * pc_z[k] * shg_254[k];

        t_357[k] = f_1 * shf0_170[k]
                   - f_2 * shf1_170[k]
                   + f_3 * pc_x[k] * shg_255[k];

        t_358[k] = f_11 * sgg_180[k]
                   + f_3 * pc_y[k] * shg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_x, pc_y, pc_z, sgg_165, sgg_182, shf0_173, \
                         shf1_173, shg_255, shg_257, shg_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * sgg_165[k]
                   + f_3 * pc_z[k] * shg_255[k];

        t_360[k] = f_4 * shf0_173[k]
                   - f_5 * shf1_173[k]
                   + f_3 * pc_x[k] * shg_258[k];

        t_361[k] = f_11 * sgg_182[k]
                   + f_3 * pc_y[k] * shg_257[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pc_x, pc_y, pc_z, sgg_168, sgg_185, \
                         shf0_175, shf0_176, shf1_175, shf1_176, shg_258, shg_260, \
                         shg_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_4 * shf0_175[k]
                   - f_5 * shf1_175[k]
                   + f_3 * pc_x[k] * shg_260[k];

        t_363[k] = f_6 * shf0_176[k]
                   - f_7 * shf1_176[k]
                   + f_3 * pc_x[k] * shg_261[k];

        t_364[k] = f_10 * sgg_168[k]
                   + f_3 * pc_z[k] * shg_258[k];

        t_365[k] = f_11 * sgg_185[k]
                   + f_3 * pc_y[k] * shg_260[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, t_370, t_371, pc_x, shf0_179, shf1_179, \
                         shg_264, shg_265, shg_266, shg_267, shg_268, \
                         shg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_6 * shf0_179[k]
                   - f_7 * shf1_179[k]
                   + f_3 * pc_x[k] * shg_264[k];

        t_367[k] = f_3 * pc_x[k] * shg_265[k];

        t_368[k] = f_3 * pc_x[k] * shg_266[k];

        t_369[k] = f_3 * pc_x[k] * shg_267[k];

        t_370[k] = f_3 * pc_x[k] * shg_268[k];

        t_371[k] = f_3 * pc_x[k] * shg_269[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, sgg_175, sgg_190, sgg_192, shf0_176, \
                         shf0_178, shf1_176, shf1_178, shg_265, \
                         shg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_11 * sgg_190[k]
                   + f_1 * shf0_176[k]
                   - f_2 * shf1_176[k]
                   + f_3 * pc_y[k] * shg_265[k];

        t_373[k] = f_10 * sgg_175[k]
                   + f_3 * pc_z[k] * shg_265[k];

        t_374[k] = f_11 * sgg_192[k]
                   + f_4 * shf0_178[k]
                   - f_5 * shf1_178[k]
                   + f_3 * pc_y[k] * shg_267[k];
    }
}

static auto
compute_prim_shh_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgh0,
                                                          const size_t sgg, const size_t sgh1,
                                                          const size_t shf0, const size_t shf1,
                                                          const size_t shg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
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
    const auto f_12 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgh0_294 = buffer.data(sgh0 + 294);
    const auto *sgh0_299 = buffer.data(sgh0 + 299);
    const auto *sgh0_303 = buffer.data(sgh0 + 303);
    const auto *sgh0_309 = buffer.data(sgh0 + 309);
    const auto *sgh0_311 = buffer.data(sgh0 + 311);
    const auto *sgh0_312 = buffer.data(sgh0 + 312);
    const auto *sgh0_314 = buffer.data(sgh0 + 314);

    const auto *sgg_179 = buffer.data(sgg + 179);
    const auto *sgg_180 = buffer.data(sgg + 180);
    const auto *sgg_183 = buffer.data(sgg + 183);
    const auto *sgg_190 = buffer.data(sgg + 190);
    const auto *sgg_193 = buffer.data(sgg + 193);
    const auto *sgg_194 = buffer.data(sgg + 194);
    const auto *sgg_195 = buffer.data(sgg + 195);
    const auto *sgg_197 = buffer.data(sgg + 197);
    const auto *sgg_198 = buffer.data(sgg + 198);
    const auto *sgg_200 = buffer.data(sgg + 200);
    const auto *sgg_205 = buffer.data(sgg + 205);
    const auto *sgg_207 = buffer.data(sgg + 207);
    const auto *sgg_208 = buffer.data(sgg + 208);
    const auto *sgg_209 = buffer.data(sgg + 209);
    const auto *sgg_210 = buffer.data(sgg + 210);
    const auto *sgg_212 = buffer.data(sgg + 212);
    const auto *sgg_213 = buffer.data(sgg + 213);
    const auto *sgg_215 = buffer.data(sgg + 215);
    const auto *sgg_220 = buffer.data(sgg + 220);
    const auto *sgg_222 = buffer.data(sgg + 222);
    const auto *sgg_223 = buffer.data(sgg + 223);
    const auto *sgg_224 = buffer.data(sgg + 224);

    const auto *sgh1_294 = buffer.data(sgh1 + 294);
    const auto *sgh1_299 = buffer.data(sgh1 + 299);
    const auto *sgh1_303 = buffer.data(sgh1 + 303);
    const auto *sgh1_309 = buffer.data(sgh1 + 309);
    const auto *sgh1_311 = buffer.data(sgh1 + 311);
    const auto *sgh1_312 = buffer.data(sgh1 + 312);
    const auto *sgh1_314 = buffer.data(sgh1 + 314);

    const auto *shf0_179 = buffer.data(shf0 + 179);
    const auto *shf0_180 = buffer.data(shf0 + 180);
    const auto *shf0_183 = buffer.data(shf0 + 183);
    const auto *shf0_185 = buffer.data(shf0 + 185);
    const auto *shf0_186 = buffer.data(shf0 + 186);
    const auto *shf0_188 = buffer.data(shf0 + 188);
    const auto *shf0_189 = buffer.data(shf0 + 189);
    const auto *shf0_193 = buffer.data(shf0 + 193);
    const auto *shf0_196 = buffer.data(shf0 + 196);
    const auto *shf0_200 = buffer.data(shf0 + 200);
    const auto *shf0_203 = buffer.data(shf0 + 203);
    const auto *shf0_205 = buffer.data(shf0 + 205);
    const auto *shf0_206 = buffer.data(shf0 + 206);
    const auto *shf0_208 = buffer.data(shf0 + 208);
    const auto *shf0_209 = buffer.data(shf0 + 209);

    const auto *shf1_179 = buffer.data(shf1 + 179);
    const auto *shf1_180 = buffer.data(shf1 + 180);
    const auto *shf1_183 = buffer.data(shf1 + 183);
    const auto *shf1_185 = buffer.data(shf1 + 185);
    const auto *shf1_186 = buffer.data(shf1 + 186);
    const auto *shf1_188 = buffer.data(shf1 + 188);
    const auto *shf1_189 = buffer.data(shf1 + 189);
    const auto *shf1_193 = buffer.data(shf1 + 193);
    const auto *shf1_196 = buffer.data(shf1 + 196);
    const auto *shf1_200 = buffer.data(shf1 + 200);
    const auto *shf1_203 = buffer.data(shf1 + 203);
    const auto *shf1_205 = buffer.data(shf1 + 205);
    const auto *shf1_206 = buffer.data(shf1 + 206);
    const auto *shf1_208 = buffer.data(shf1 + 208);
    const auto *shf1_209 = buffer.data(shf1 + 209);

    const auto *shg_268 = buffer.data(shg + 268);
    const auto *shg_269 = buffer.data(shg + 269);
    const auto *shg_270 = buffer.data(shg + 270);
    const auto *shg_272 = buffer.data(shg + 272);
    const auto *shg_273 = buffer.data(shg + 273);
    const auto *shg_275 = buffer.data(shg + 275);
    const auto *shg_276 = buffer.data(shg + 276);
    const auto *shg_279 = buffer.data(shg + 279);
    const auto *shg_280 = buffer.data(shg + 280);
    const auto *shg_281 = buffer.data(shg + 281);
    const auto *shg_282 = buffer.data(shg + 282);
    const auto *shg_283 = buffer.data(shg + 283);
    const auto *shg_284 = buffer.data(shg + 284);
    const auto *shg_285 = buffer.data(shg + 285);
    const auto *shg_287 = buffer.data(shg + 287);
    const auto *shg_288 = buffer.data(shg + 288);
    const auto *shg_290 = buffer.data(shg + 290);
    const auto *shg_291 = buffer.data(shg + 291);
    const auto *shg_295 = buffer.data(shg + 295);
    const auto *shg_296 = buffer.data(shg + 296);
    const auto *shg_297 = buffer.data(shg + 297);
    const auto *shg_298 = buffer.data(shg + 298);
    const auto *shg_299 = buffer.data(shg + 299);
    const auto *shg_300 = buffer.data(shg + 300);
    const auto *shg_302 = buffer.data(shg + 302);
    const auto *shg_303 = buffer.data(shg + 303);
    const auto *shg_305 = buffer.data(shg + 305);
    const auto *shg_306 = buffer.data(shg + 306);
    const auto *shg_309 = buffer.data(shg + 309);
    const auto *shg_310 = buffer.data(shg + 310);
    const auto *shg_311 = buffer.data(shg + 311);
    const auto *shg_312 = buffer.data(shg + 312);
    const auto *shg_313 = buffer.data(shg + 313);
    const auto *shg_314 = buffer.data(shg + 314);

#pragma omp simd aligned(t_375, t_376, t_377, pc_y, pc_z, sgg_179, sgg_193, sgg_194, shf0_179, \
                         shf1_179, shg_268, shg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_11 * sgg_193[k]
                   + f_6 * shf0_179[k]
                   - f_7 * shf1_179[k]
                   + f_3 * pc_y[k] * shg_268[k];

        t_376[k] = f_11 * sgg_194[k]
                   + f_3 * pc_y[k] * shg_269[k];

        t_377[k] = f_10 * sgg_179[k]
                   + f_1 * shf0_179[k]
                   - f_2 * shf1_179[k]
                   + f_3 * pc_z[k] * shg_269[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, pc_y, pc_z, sgg_180, sgg_195, \
                         shf0_180, shf0_183, shf1_180, shf1_183, shg_270, \
                         shg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_1 * shf0_180[k]
                   - f_2 * shf1_180[k]
                   + f_3 * pc_x[k] * shg_270[k];

        t_379[k] = f_10 * sgg_195[k]
                   + f_3 * pc_y[k] * shg_270[k];

        t_380[k] = f_11 * sgg_180[k]
                   + f_3 * pc_z[k] * shg_270[k];

        t_381[k] = f_4 * shf0_183[k]
                   - f_5 * shf1_183[k]
                   + f_3 * pc_x[k] * shg_273[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pc_x, pc_y, sgg_197, shf0_185, shf0_186, \
                         shf1_185, shf1_186, shg_272, shg_275, \
                         shg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_10 * sgg_197[k]
                   + f_3 * pc_y[k] * shg_272[k];

        t_383[k] = f_4 * shf0_185[k]
                   - f_5 * shf1_185[k]
                   + f_3 * pc_x[k] * shg_275[k];

        t_384[k] = f_6 * shf0_186[k]
                   - f_7 * shf1_186[k]
                   + f_3 * pc_x[k] * shg_276[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pc_x, pc_y, pc_z, sgg_183, sgg_200, \
                         shf0_189, shf1_189, shg_273, shg_275, shg_279, \
                         shg_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = f_11 * sgg_183[k]
                   + f_3 * pc_z[k] * shg_273[k];

        t_386[k] = f_10 * sgg_200[k]
                   + f_3 * pc_y[k] * shg_275[k];

        t_387[k] = f_6 * shf0_189[k]
                   - f_7 * shf1_189[k]
                   + f_3 * pc_x[k] * shg_279[k];

        t_388[k] = f_3 * pc_x[k] * shg_280[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, pc_y, sgg_205, shf0_186, \
                         shf1_186, shg_280, shg_281, shg_282, shg_283, \
                         shg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_3 * pc_x[k] * shg_281[k];

        t_390[k] = f_3 * pc_x[k] * shg_282[k];

        t_391[k] = f_3 * pc_x[k] * shg_283[k];

        t_392[k] = f_3 * pc_x[k] * shg_284[k];

        t_393[k] = f_10 * sgg_205[k]
                   + f_1 * shf0_186[k]
                   - f_2 * shf1_186[k]
                   + f_3 * pc_y[k] * shg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pc_y, pc_z, sgg_190, sgg_207, sgg_208, shf0_188, \
                         shf0_189, shf1_188, shf1_189, shg_280, shg_282, \
                         shg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_11 * sgg_190[k]
                   + f_3 * pc_z[k] * shg_280[k];

        t_395[k] = f_10 * sgg_207[k]
                   + f_4 * shf0_188[k]
                   - f_5 * shf1_188[k]
                   + f_3 * pc_y[k] * shg_282[k];

        t_396[k] = f_10 * sgg_208[k]
                   + f_6 * shf0_189[k]
                   - f_7 * shf1_189[k]
                   + f_3 * pc_y[k] * shg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pb_y, pc_y, pc_z, sgh0_294, sgg_194, \
                         sgg_209, sgg_210, sgh1_294, shf0_189, shf1_189, shg_284, \
                         shg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_10 * sgg_209[k]
                   + f_3 * pc_y[k] * shg_284[k];

        t_398[k] = f_11 * sgg_194[k]
                   + f_1 * shf0_189[k]
                   - f_2 * shf1_189[k]
                   + f_3 * pc_z[k] * shg_284[k];

        t_399[k] = pb_y[k] * sgh0_294[k]
                   - f_8 * pc_y[k] * sgh1_294[k];

        t_400[k] = f_9 * sgg_210[k]
                   + f_3 * pc_y[k] * shg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_x, pc_y, pc_z, sgg_195, sgg_212, shf0_193, \
                         shf1_193, shg_285, shg_287, shg_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_12 * sgg_195[k]
                   + f_3 * pc_z[k] * shg_285[k];

        t_402[k] = f_4 * shf0_193[k]
                   - f_5 * shf1_193[k]
                   + f_3 * pc_x[k] * shg_288[k];

        t_403[k] = f_9 * sgg_212[k]
                   + f_3 * pc_y[k] * shg_287[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pb_y, pc_x, pc_y, pc_z, sgh0_299, sgg_198, \
                         sgh1_299, shf0_196, shf1_196, shg_288, \
                         shg_291 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * sgh0_299[k]
                   - f_8 * pc_y[k] * sgh1_299[k];

        t_405[k] = f_6 * shf0_196[k]
                   - f_7 * shf1_196[k]
                   + f_3 * pc_x[k] * shg_291[k];

        t_406[k] = f_12 * sgg_198[k]
                   + f_3 * pc_z[k] * shg_288[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, t_410, t_411, pb_y, pc_x, pc_y, sgh0_303, \
                         sgg_215, sgh1_303, shg_290, shg_295, shg_296, \
                         shg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_9 * sgg_215[k]
                   + f_3 * pc_y[k] * shg_290[k];

        t_408[k] = pb_y[k] * sgh0_303[k]
                   - f_8 * pc_y[k] * sgh1_303[k];

        t_409[k] = f_3 * pc_x[k] * shg_295[k];

        t_410[k] = f_3 * pc_x[k] * shg_296[k];

        t_411[k] = f_3 * pc_x[k] * shg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pb_y, pc_x, pc_y, pc_z, sgh0_309, \
                         sgg_205, sgg_220, sgh1_309, shg_295, shg_298, \
                         shg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_3 * pc_x[k] * shg_298[k];

        t_413[k] = f_3 * pc_x[k] * shg_299[k];

        t_414[k] = pb_y[k] * sgh0_309[k]
                   + f_0 * sgg_220[k]
                   - f_8 * pc_y[k] * sgh1_309[k];

        t_415[k] = f_12 * sgg_205[k]
                   + f_3 * pc_z[k] * shg_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pb_y, pc_y, sgh0_311, sgh0_312, sgh0_314, \
                         sgg_222, sgg_223, sgg_224, sgh1_311, sgh1_312, sgh1_314, \
                         shg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pb_y[k] * sgh0_311[k]
                   + f_11 * sgg_222[k]
                   - f_8 * pc_y[k] * sgh1_311[k];

        t_417[k] = pb_y[k] * sgh0_312[k]
                   + f_10 * sgg_223[k]
                   - f_8 * pc_y[k] * sgh1_312[k];

        t_418[k] = f_9 * sgg_224[k]
                   + f_3 * pc_y[k] * shg_299[k];

        t_419[k] = pb_y[k] * sgh0_314[k]
                   - f_8 * pc_y[k] * sgh1_314[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, sgg_210, \
                         shf0_200, shf0_203, shf1_200, shf1_203, shg_300, shg_302, \
                         shg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_1 * shf0_200[k]
                   - f_2 * shf1_200[k]
                   + f_3 * pc_x[k] * shg_300[k];

        t_421[k] = f_3 * pc_y[k] * shg_300[k];

        t_422[k] = f_0 * sgg_210[k]
                   + f_3 * pc_z[k] * shg_300[k];

        t_423[k] = f_4 * shf0_203[k]
                   - f_5 * shf1_203[k]
                   + f_3 * pc_x[k] * shg_303[k];

        t_424[k] = f_3 * pc_y[k] * shg_302[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, sgg_213, shf0_205, \
                         shf0_206, shf1_205, shf1_206, shg_303, shg_305, \
                         shg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * shf0_205[k]
                   - f_5 * shf1_205[k]
                   + f_3 * pc_x[k] * shg_305[k];

        t_426[k] = f_6 * shf0_206[k]
                   - f_7 * shf1_206[k]
                   + f_3 * pc_x[k] * shg_306[k];

        t_427[k] = f_0 * sgg_213[k]
                   + f_3 * pc_z[k] * shg_303[k];

        t_428[k] = f_3 * pc_y[k] * shg_305[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, t_434, pc_x, shf0_209, shf1_209, \
                         shg_309, shg_310, shg_311, shg_312, shg_313, \
                         shg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_6 * shf0_209[k]
                   - f_7 * shf1_209[k]
                   + f_3 * pc_x[k] * shg_309[k];

        t_430[k] = f_3 * pc_x[k] * shg_310[k];

        t_431[k] = f_3 * pc_x[k] * shg_311[k];

        t_432[k] = f_3 * pc_x[k] * shg_312[k];

        t_433[k] = f_3 * pc_x[k] * shg_313[k];

        t_434[k] = f_3 * pc_x[k] * shg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, pc_y, pc_z, sgg_220, shf0_206, shf0_208, \
                         shf0_209, shf1_206, shf1_208, shf1_209, shg_310, shg_312, \
                         shg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_1 * shf0_206[k]
                   - f_2 * shf1_206[k]
                   + f_3 * pc_y[k] * shg_310[k];

        t_436[k] = f_0 * sgg_220[k]
                   + f_3 * pc_z[k] * shg_310[k];

        t_437[k] = f_4 * shf0_208[k]
                   - f_5 * shf1_208[k]
                   + f_3 * pc_y[k] * shg_312[k];

        t_438[k] = f_6 * shf0_209[k]
                   - f_7 * shf1_209[k]
                   + f_3 * pc_y[k] * shg_313[k];
    }

#pragma omp simd aligned(t_439, t_440, pc_y, pc_z, sgg_224, shf0_209, shf1_209, \
                         shg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_439[k] = f_3 * pc_y[k] * shg_314[k];

        t_440[k] = f_0 * sgg_224[k]
                   + f_1 * shf0_209[k]
                   - f_2 * shf1_209[k]
                   + f_3 * pc_z[k] * shg_314[k];
    }
}

auto
compute_prim_shh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgh0, const size_t sgg,
                                                   const size_t sgh1, const size_t shf0,
                                                   const size_t shf1, const size_t shg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shh_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgh0, sgg,
                                                              sgh1, shf0, shf1, shg, ncols,
                                                              gamma, p, q);

    compute_prim_shh_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgh0, sgg,
                                                              sgh1, shf0, shf1, shg, ncols,
                                                              gamma, p, q);

    compute_prim_shh_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgh0, sgg,
                                                              sgh1, shf0, shf1, shg, ncols,
                                                              gamma, p, q);

    compute_prim_shh_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sgh0, sgg,
                                                              sgh1, shf0, shf1, shg, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
