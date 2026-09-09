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


#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgg0,
                                                          const size_t sgf, const size_t sgg1,
                                                          const size_t shd0, const size_t shd1,
                                                          const size_t shf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_0 = buffer.data(sgg0 + 0);
    const auto *sgg0_3 = buffer.data(sgg0 + 3);
    const auto *sgg0_5 = buffer.data(sgg0 + 5);
    const auto *sgg0_10 = buffer.data(sgg0 + 10);
    const auto *sgg0_14 = buffer.data(sgg0 + 14);
    const auto *sgg0_18 = buffer.data(sgg0 + 18);
    const auto *sgg0_25 = buffer.data(sgg0 + 25);
    const auto *sgg0_30 = buffer.data(sgg0 + 30);
    const auto *sgg0_35 = buffer.data(sgg0 + 35);
    const auto *sgg0_44 = buffer.data(sgg0 + 44);
    const auto *sgg0_45 = buffer.data(sgg0 + 45);
    const auto *sgg0_48 = buffer.data(sgg0 + 48);
    const auto *sgg0_55 = buffer.data(sgg0 + 55);
    const auto *sgg0_75 = buffer.data(sgg0 + 75);
    const auto *sgg0_78 = buffer.data(sgg0 + 78);

    const auto *sgf_0 = buffer.data(sgf + 0);
    const auto *sgf_1 = buffer.data(sgf + 1);
    const auto *sgf_2 = buffer.data(sgf + 2);
    const auto *sgf_3 = buffer.data(sgf + 3);
    const auto *sgf_5 = buffer.data(sgf + 5);
    const auto *sgf_6 = buffer.data(sgf + 6);
    const auto *sgf_7 = buffer.data(sgf + 7);
    const auto *sgf_8 = buffer.data(sgf + 8);
    const auto *sgf_9 = buffer.data(sgf + 9);
    const auto *sgf_10 = buffer.data(sgf + 10);
    const auto *sgf_12 = buffer.data(sgf + 12);
    const auto *sgf_16 = buffer.data(sgf + 16);
    const auto *sgf_17 = buffer.data(sgf + 17);
    const auto *sgf_18 = buffer.data(sgf + 18);
    const auto *sgf_19 = buffer.data(sgf + 19);
    const auto *sgf_20 = buffer.data(sgf + 20);
    const auto *sgf_22 = buffer.data(sgf + 22);
    const auto *sgf_26 = buffer.data(sgf + 26);
    const auto *sgf_27 = buffer.data(sgf + 27);
    const auto *sgf_28 = buffer.data(sgf + 28);
    const auto *sgf_29 = buffer.data(sgf + 29);
    const auto *sgf_30 = buffer.data(sgf + 30);
    const auto *sgf_32 = buffer.data(sgf + 32);
    const auto *sgf_33 = buffer.data(sgf + 33);
    const auto *sgf_35 = buffer.data(sgf + 35);
    const auto *sgf_36 = buffer.data(sgf + 36);
    const auto *sgf_37 = buffer.data(sgf + 37);
    const auto *sgf_38 = buffer.data(sgf + 38);
    const auto *sgf_39 = buffer.data(sgf + 39);
    const auto *sgf_40 = buffer.data(sgf + 40);
    const auto *sgf_42 = buffer.data(sgf + 42);
    const auto *sgf_46 = buffer.data(sgf + 46);
    const auto *sgf_47 = buffer.data(sgf + 47);
    const auto *sgf_48 = buffer.data(sgf + 48);
    const auto *sgf_49 = buffer.data(sgf + 49);
    const auto *sgf_50 = buffer.data(sgf + 50);
    const auto *sgf_51 = buffer.data(sgf + 51);
    const auto *sgf_52 = buffer.data(sgf + 52);
    const auto *sgf_53 = buffer.data(sgf + 53);
    const auto *sgf_55 = buffer.data(sgf + 55);
    const auto *sgf_56 = buffer.data(sgf + 56);
    const auto *sgf_57 = buffer.data(sgf + 57);
    const auto *sgf_58 = buffer.data(sgf + 58);
    const auto *sgf_59 = buffer.data(sgf + 59);
    const auto *sgf_60 = buffer.data(sgf + 60);
    const auto *sgf_63 = buffer.data(sgf + 63);
    const auto *sgf_65 = buffer.data(sgf + 65);
    const auto *sgf_66 = buffer.data(sgf + 66);
    const auto *sgf_67 = buffer.data(sgf + 67);
    const auto *sgf_68 = buffer.data(sgf + 68);
    const auto *sgf_69 = buffer.data(sgf + 69);
    const auto *sgf_75 = buffer.data(sgf + 75);
    const auto *sgf_76 = buffer.data(sgf + 76);
    const auto *sgf_77 = buffer.data(sgf + 77);
    const auto *sgf_78 = buffer.data(sgf + 78);
    const auto *sgf_79 = buffer.data(sgf + 79);

    const auto *sgg1_0 = buffer.data(sgg1 + 0);
    const auto *sgg1_3 = buffer.data(sgg1 + 3);
    const auto *sgg1_5 = buffer.data(sgg1 + 5);
    const auto *sgg1_10 = buffer.data(sgg1 + 10);
    const auto *sgg1_14 = buffer.data(sgg1 + 14);
    const auto *sgg1_18 = buffer.data(sgg1 + 18);
    const auto *sgg1_25 = buffer.data(sgg1 + 25);
    const auto *sgg1_30 = buffer.data(sgg1 + 30);
    const auto *sgg1_35 = buffer.data(sgg1 + 35);
    const auto *sgg1_44 = buffer.data(sgg1 + 44);
    const auto *sgg1_45 = buffer.data(sgg1 + 45);
    const auto *sgg1_48 = buffer.data(sgg1 + 48);
    const auto *sgg1_55 = buffer.data(sgg1 + 55);
    const auto *sgg1_75 = buffer.data(sgg1 + 75);
    const auto *sgg1_78 = buffer.data(sgg1 + 78);

    const auto *shd0_0 = buffer.data(shd0 + 0);
    const auto *shd0_3 = buffer.data(shd0 + 3);
    const auto *shd0_5 = buffer.data(shd0 + 5);
    const auto *shd0_9 = buffer.data(shd0 + 9);
    const auto *shd0_11 = buffer.data(shd0 + 11);
    const auto *shd0_17 = buffer.data(shd0 + 17);
    const auto *shd0_18 = buffer.data(shd0 + 18);
    const auto *shd0_21 = buffer.data(shd0 + 21);
    const auto *shd0_23 = buffer.data(shd0 + 23);
    const auto *shd0_29 = buffer.data(shd0 + 29);
    const auto *shd0_30 = buffer.data(shd0 + 30);
    const auto *shd0_33 = buffer.data(shd0 + 33);
    const auto *shd0_35 = buffer.data(shd0 + 35);
    const auto *shd0_36 = buffer.data(shd0 + 36);
    const auto *shd0_39 = buffer.data(shd0 + 39);
    const auto *shd0_41 = buffer.data(shd0 + 41);
    const auto *shd0_47 = buffer.data(shd0 + 47);

    const auto *shd1_0 = buffer.data(shd1 + 0);
    const auto *shd1_3 = buffer.data(shd1 + 3);
    const auto *shd1_5 = buffer.data(shd1 + 5);
    const auto *shd1_9 = buffer.data(shd1 + 9);
    const auto *shd1_11 = buffer.data(shd1 + 11);
    const auto *shd1_17 = buffer.data(shd1 + 17);
    const auto *shd1_18 = buffer.data(shd1 + 18);
    const auto *shd1_21 = buffer.data(shd1 + 21);
    const auto *shd1_23 = buffer.data(shd1 + 23);
    const auto *shd1_29 = buffer.data(shd1 + 29);
    const auto *shd1_30 = buffer.data(shd1 + 30);
    const auto *shd1_33 = buffer.data(shd1 + 33);
    const auto *shd1_35 = buffer.data(shd1 + 35);
    const auto *shd1_36 = buffer.data(shd1 + 36);
    const auto *shd1_39 = buffer.data(shd1 + 39);
    const auto *shd1_41 = buffer.data(shd1 + 41);
    const auto *shd1_47 = buffer.data(shd1 + 47);

    const auto *shf_0 = buffer.data(shf + 0);
    const auto *shf_2 = buffer.data(shf + 2);
    const auto *shf_3 = buffer.data(shf + 3);
    const auto *shf_5 = buffer.data(shf + 5);
    const auto *shf_6 = buffer.data(shf + 6);
    const auto *shf_7 = buffer.data(shf + 7);
    const auto *shf_8 = buffer.data(shf + 8);
    const auto *shf_9 = buffer.data(shf + 9);
    const auto *shf_10 = buffer.data(shf + 10);
    const auto *shf_12 = buffer.data(shf + 12);
    const auto *shf_16 = buffer.data(shf + 16);
    const auto *shf_17 = buffer.data(shf + 17);
    const auto *shf_18 = buffer.data(shf + 18);
    const auto *shf_19 = buffer.data(shf + 19);
    const auto *shf_20 = buffer.data(shf + 20);
    const auto *shf_22 = buffer.data(shf + 22);
    const auto *shf_26 = buffer.data(shf + 26);
    const auto *shf_27 = buffer.data(shf + 27);
    const auto *shf_28 = buffer.data(shf + 28);
    const auto *shf_29 = buffer.data(shf + 29);
    const auto *shf_30 = buffer.data(shf + 30);
    const auto *shf_32 = buffer.data(shf + 32);
    const auto *shf_33 = buffer.data(shf + 33);
    const auto *shf_35 = buffer.data(shf + 35);
    const auto *shf_36 = buffer.data(shf + 36);
    const auto *shf_37 = buffer.data(shf + 37);
    const auto *shf_38 = buffer.data(shf + 38);
    const auto *shf_39 = buffer.data(shf + 39);
    const auto *shf_40 = buffer.data(shf + 40);
    const auto *shf_42 = buffer.data(shf + 42);
    const auto *shf_46 = buffer.data(shf + 46);
    const auto *shf_47 = buffer.data(shf + 47);
    const auto *shf_48 = buffer.data(shf + 48);
    const auto *shf_49 = buffer.data(shf + 49);
    const auto *shf_50 = buffer.data(shf + 50);
    const auto *shf_52 = buffer.data(shf + 52);
    const auto *shf_53 = buffer.data(shf + 53);
    const auto *shf_55 = buffer.data(shf + 55);
    const auto *shf_56 = buffer.data(shf + 56);
    const auto *shf_57 = buffer.data(shf + 57);
    const auto *shf_58 = buffer.data(shf + 58);
    const auto *shf_59 = buffer.data(shf + 59);
    const auto *shf_60 = buffer.data(shf + 60);
    const auto *shf_62 = buffer.data(shf + 62);
    const auto *shf_63 = buffer.data(shf + 63);
    const auto *shf_65 = buffer.data(shf + 65);
    const auto *shf_66 = buffer.data(shf + 66);
    const auto *shf_67 = buffer.data(shf + 67);
    const auto *shf_68 = buffer.data(shf + 68);
    const auto *shf_69 = buffer.data(shf + 69);
    const auto *shf_70 = buffer.data(shf + 70);
    const auto *shf_72 = buffer.data(shf + 72);
    const auto *shf_75 = buffer.data(shf + 75);
    const auto *shf_76 = buffer.data(shf + 76);
    const auto *shf_77 = buffer.data(shf + 77);
    const auto *shf_78 = buffer.data(shf + 78);
    const auto *shf_79 = buffer.data(shf + 79);
    const auto *shf_80 = buffer.data(shf + 80);
    const auto *shf_82 = buffer.data(shf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sgf_0, sgf_3, shd0_0, shd0_3, \
                         shd1_0, shd1_3, shf_0, shf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgf_0[k]
                 + f_1 * shd0_0[k]
                 - f_2 * shd1_0[k]
                 + f_3 * pc_x[k] * shf_0[k];

        t_1[k] = f_3 * pc_y[k] * shf_0[k];

        t_2[k] = f_3 * pc_z[k] * shf_0[k];

        t_3[k] = f_0 * sgf_3[k]
                 + f_4 * shd0_3[k]
                 - f_5 * shd1_3[k]
                 + f_3 * pc_x[k] * shf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, sgf_5, sgf_6, sgf_7, shd0_5, shd1_5, \
                         shf_2, shf_5, shf_6, shf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * shf_2[k];

        t_5[k] = f_0 * sgf_5[k]
                 + f_4 * shd0_5[k]
                 - f_5 * shd1_5[k]
                 + f_3 * pc_x[k] * shf_5[k];

        t_6[k] = f_0 * sgf_6[k]
                 + f_3 * pc_x[k] * shf_6[k];

        t_7[k] = f_0 * sgf_7[k]
                 + f_3 * pc_x[k] * shf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, sgf_8, sgf_9, shd0_3, shd1_3, \
                         shf_6, shf_8, shf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sgf_8[k]
                 + f_3 * pc_x[k] * shf_8[k];

        t_9[k] = f_0 * sgf_9[k]
                 + f_3 * pc_x[k] * shf_9[k];

        t_10[k] = f_1 * shd0_3[k]
                  - f_2 * shd1_3[k]
                  + f_3 * pc_y[k] * shf_6[k];

        t_11[k] = f_3 * pc_z[k] * shf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, sgg0_0, sgf_0, \
                         sgg1_0, shd0_5, shd1_5, shf_8, shf_9, shf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * shd0_5[k]
                  - f_5 * shd1_5[k]
                  + f_3 * pc_y[k] * shf_8[k];

        t_13[k] = f_3 * pc_y[k] * shf_9[k];

        t_14[k] = f_1 * shd0_5[k]
                  - f_2 * shd1_5[k]
                  + f_3 * pc_z[k] * shf_9[k];

        t_15[k] = pb_y[k] * sgg0_0[k]
                  - f_6 * pc_y[k] * sgg1_0[k];

        t_16[k] = f_7 * sgf_0[k]
                  + f_3 * pc_y[k] * shf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, sgg0_3, sgg0_5, sgf_1, \
                         sgf_2, sgg1_3, sgg1_5, shf_10, shf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * shf_10[k];

        t_18[k] = pb_y[k] * sgg0_3[k]
                  + f_8 * sgf_1[k]
                  - f_6 * pc_y[k] * sgg1_3[k];

        t_19[k] = f_7 * sgf_2[k]
                  + f_3 * pc_y[k] * shf_12[k];

        t_20[k] = pb_y[k] * sgg0_5[k]
                  - f_6 * pc_y[k] * sgg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, sgf_16, sgf_17, sgf_18, sgf_19, shf_16, \
                         shf_17, shf_18, shf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * sgf_16[k]
                  + f_3 * pc_x[k] * shf_16[k];

        t_22[k] = f_9 * sgf_17[k]
                  + f_3 * pc_x[k] * shf_17[k];

        t_23[k] = f_9 * sgf_18[k]
                  + f_3 * pc_x[k] * shf_18[k];

        t_24[k] = f_9 * sgf_19[k]
                  + f_3 * pc_x[k] * shf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, sgf_6, sgf_8, sgf_9, shd0_9, \
                         shd0_11, shd1_9, shd1_11, shf_16, shf_18, \
                         shf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * sgf_6[k]
                  + f_1 * shd0_9[k]
                  - f_2 * shd1_9[k]
                  + f_3 * pc_y[k] * shf_16[k];

        t_26[k] = f_3 * pc_z[k] * shf_16[k];

        t_27[k] = f_7 * sgf_8[k]
                  + f_4 * shd0_11[k]
                  - f_5 * shd1_11[k]
                  + f_3 * pc_y[k] * shf_18[k];

        t_28[k] = f_7 * sgf_9[k]
                  + f_3 * pc_y[k] * shf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, sgg0_0, sgg0_14, \
                         sgf_0, sgg1_0, sgg1_14, shf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * sgg0_14[k]
                  - f_6 * pc_y[k] * sgg1_14[k];

        t_30[k] = pb_z[k] * sgg0_0[k]
                  - f_6 * pc_z[k] * sgg1_0[k];

        t_31[k] = f_3 * pc_y[k] * shf_20[k];

        t_32[k] = f_7 * sgf_0[k]
                  + f_3 * pc_z[k] * shf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, sgg0_3, sgg0_5, \
                         sgf_2, sgf_26, sgg1_3, sgg1_5, shf_22, \
                         shf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * sgg0_3[k]
                  - f_6 * pc_z[k] * sgg1_3[k];

        t_34[k] = f_3 * pc_y[k] * shf_22[k];

        t_35[k] = pb_z[k] * sgg0_5[k]
                  + f_8 * sgf_2[k]
                  - f_6 * pc_z[k] * sgg1_5[k];

        t_36[k] = f_9 * sgf_26[k]
                  + f_3 * pc_x[k] * shf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, sgg0_10, sgf_27, sgf_28, \
                         sgf_29, sgg1_10, shf_27, shf_28, shf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * sgf_27[k]
                  + f_3 * pc_x[k] * shf_27[k];

        t_38[k] = f_9 * sgf_28[k]
                  + f_3 * pc_x[k] * shf_28[k];

        t_39[k] = f_9 * sgf_29[k]
                  + f_3 * pc_x[k] * shf_29[k];

        t_40[k] = pb_z[k] * sgg0_10[k]
                  - f_6 * pc_z[k] * sgg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sgf_6, sgf_9, shd0_17, shd1_17, \
                         shf_26, shf_28, shf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sgf_6[k]
                  + f_3 * pc_z[k] * shf_26[k];

        t_42[k] = f_4 * shd0_17[k]
                  - f_5 * shd1_17[k]
                  + f_3 * pc_y[k] * shf_28[k];

        t_43[k] = f_3 * pc_y[k] * shf_29[k];

        t_44[k] = f_7 * sgf_9[k]
                  + f_1 * shd0_17[k]
                  - f_2 * shd1_17[k]
                  + f_3 * pc_z[k] * shf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, sgf_10, sgf_30, sgf_33, \
                         shd0_18, shd0_21, shd1_18, shd1_21, shf_30, \
                         shf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * sgf_30[k]
                  + f_1 * shd0_18[k]
                  - f_2 * shd1_18[k]
                  + f_3 * pc_x[k] * shf_30[k];

        t_46[k] = f_8 * sgf_10[k]
                  + f_3 * pc_y[k] * shf_30[k];

        t_47[k] = f_3 * pc_z[k] * shf_30[k];

        t_48[k] = f_10 * sgf_33[k]
                  + f_4 * shd0_21[k]
                  - f_5 * shd1_21[k]
                  + f_3 * pc_x[k] * shf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, sgf_12, sgf_35, sgf_36, sgf_37, \
                         shd0_23, shd1_23, shf_32, shf_35, shf_36, \
                         shf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sgf_12[k]
                  + f_3 * pc_y[k] * shf_32[k];

        t_50[k] = f_10 * sgf_35[k]
                  + f_4 * shd0_23[k]
                  - f_5 * shd1_23[k]
                  + f_3 * pc_x[k] * shf_35[k];

        t_51[k] = f_10 * sgf_36[k]
                  + f_3 * pc_x[k] * shf_36[k];

        t_52[k] = f_10 * sgf_37[k]
                  + f_3 * pc_x[k] * shf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sgf_16, sgf_38, sgf_39, \
                         shd0_21, shd1_21, shf_36, shf_38, shf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * sgf_38[k]
                  + f_3 * pc_x[k] * shf_38[k];

        t_54[k] = f_10 * sgf_39[k]
                  + f_3 * pc_x[k] * shf_39[k];

        t_55[k] = f_8 * sgf_16[k]
                  + f_1 * shd0_21[k]
                  - f_2 * shd1_21[k]
                  + f_3 * pc_y[k] * shf_36[k];

        t_56[k] = f_3 * pc_z[k] * shf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, sgg0_30, sgf_18, sgf_19, \
                         sgg1_30, shd0_23, shd1_23, shf_38, shf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * sgf_18[k]
                  + f_4 * shd0_23[k]
                  - f_5 * shd1_23[k]
                  + f_3 * pc_y[k] * shf_38[k];

        t_58[k] = f_8 * sgf_19[k]
                  + f_3 * pc_y[k] * shf_39[k];

        t_59[k] = f_1 * shd0_23[k]
                  - f_2 * shd1_23[k]
                  + f_3 * pc_z[k] * shf_39[k];

        t_60[k] = pb_y[k] * sgg0_30[k]
                  - f_6 * pc_y[k] * sgg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sgg0_18, sgf_10, sgf_20, \
                         sgf_22, sgg1_18, shf_40, shf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * sgf_20[k]
                  + f_3 * pc_y[k] * shf_40[k];

        t_62[k] = f_7 * sgf_10[k]
                  + f_3 * pc_z[k] * shf_40[k];

        t_63[k] = pb_z[k] * sgg0_18[k]
                  - f_6 * pc_z[k] * sgg1_18[k];

        t_64[k] = f_7 * sgf_22[k]
                  + f_3 * pc_y[k] * shf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, sgg0_35, sgf_46, sgf_47, \
                         sgf_48, sgg1_35, shf_46, shf_47, shf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sgg0_35[k]
                  - f_6 * pc_y[k] * sgg1_35[k];

        t_66[k] = f_10 * sgf_46[k]
                  + f_3 * pc_x[k] * shf_46[k];

        t_67[k] = f_10 * sgf_47[k]
                  + f_3 * pc_x[k] * shf_47[k];

        t_68[k] = f_10 * sgf_48[k]
                  + f_3 * pc_x[k] * shf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, sgg0_25, sgf_16, sgf_49, sgg1_25, \
                         shf_46, shf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * sgf_49[k]
                  + f_3 * pc_x[k] * shf_49[k];

        t_70[k] = pb_z[k] * sgg0_25[k]
                  - f_6 * pc_z[k] * sgg1_25[k];

        t_71[k] = f_7 * sgf_16[k]
                  + f_3 * pc_z[k] * shf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, sgg0_44, sgf_28, sgf_29, sgg1_44, \
                         shd0_29, shd1_29, shf_48, shf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * sgf_28[k]
                  + f_4 * shd0_29[k]
                  - f_5 * shd1_29[k]
                  + f_3 * pc_y[k] * shf_48[k];

        t_73[k] = f_7 * sgf_29[k]
                  + f_3 * pc_y[k] * shf_49[k];

        t_74[k] = pb_y[k] * sgg0_44[k]
                  - f_6 * pc_y[k] * sgg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, sgf_20, sgf_50, sgf_53, \
                         shd0_30, shd0_33, shd1_30, shd1_33, shf_50, \
                         shf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * sgf_50[k]
                  + f_1 * shd0_30[k]
                  - f_2 * shd1_30[k]
                  + f_3 * pc_x[k] * shf_50[k];

        t_76[k] = f_3 * pc_y[k] * shf_50[k];

        t_77[k] = f_8 * sgf_20[k]
                  + f_3 * pc_z[k] * shf_50[k];

        t_78[k] = f_10 * sgf_53[k]
                  + f_4 * shd0_33[k]
                  - f_5 * shd1_33[k]
                  + f_3 * pc_x[k] * shf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, sgf_55, sgf_56, sgf_57, shd0_35, \
                         shd1_35, shf_52, shf_55, shf_56, shf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * shf_52[k];

        t_80[k] = f_10 * sgf_55[k]
                  + f_4 * shd0_35[k]
                  - f_5 * shd1_35[k]
                  + f_3 * pc_x[k] * shf_55[k];

        t_81[k] = f_10 * sgf_56[k]
                  + f_3 * pc_x[k] * shf_56[k];

        t_82[k] = f_10 * sgf_57[k]
                  + f_3 * pc_x[k] * shf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, sgf_26, sgf_58, sgf_59, \
                         shd0_33, shd1_33, shf_56, shf_58, shf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * sgf_58[k]
                  + f_3 * pc_x[k] * shf_58[k];

        t_84[k] = f_10 * sgf_59[k]
                  + f_3 * pc_x[k] * shf_59[k];

        t_85[k] = f_1 * shd0_33[k]
                  - f_2 * shd1_33[k]
                  + f_3 * pc_y[k] * shf_56[k];

        t_86[k] = f_8 * sgf_26[k]
                  + f_3 * pc_z[k] * shf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, sgf_29, sgf_60, shd0_35, \
                         shd0_36, shd1_35, shd1_36, shf_58, shf_59, \
                         shf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * shd0_35[k]
                  - f_5 * shd1_35[k]
                  + f_3 * pc_y[k] * shf_58[k];

        t_88[k] = f_3 * pc_y[k] * shf_59[k];

        t_89[k] = f_8 * sgf_29[k]
                  + f_1 * shd0_35[k]
                  - f_2 * shd1_35[k]
                  + f_3 * pc_z[k] * shf_59[k];

        t_90[k] = f_8 * sgf_60[k]
                  + f_1 * shd0_36[k]
                  - f_2 * shd1_36[k]
                  + f_3 * pc_x[k] * shf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, sgf_30, sgf_32, sgf_63, \
                         shd0_39, shd1_39, shf_60, shf_62, shf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_10 * sgf_30[k]
                  + f_3 * pc_y[k] * shf_60[k];

        t_92[k] = f_3 * pc_z[k] * shf_60[k];

        t_93[k] = f_8 * sgf_63[k]
                  + f_4 * shd0_39[k]
                  - f_5 * shd1_39[k]
                  + f_3 * pc_x[k] * shf_63[k];

        t_94[k] = f_10 * sgf_32[k]
                  + f_3 * pc_y[k] * shf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, sgf_65, sgf_66, sgf_67, sgf_68, \
                         shd0_41, shd1_41, shf_65, shf_66, shf_67, \
                         shf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_8 * sgf_65[k]
                  + f_4 * shd0_41[k]
                  - f_5 * shd1_41[k]
                  + f_3 * pc_x[k] * shf_65[k];

        t_96[k] = f_8 * sgf_66[k]
                  + f_3 * pc_x[k] * shf_66[k];

        t_97[k] = f_8 * sgf_67[k]
                  + f_3 * pc_x[k] * shf_67[k];

        t_98[k] = f_8 * sgf_68[k]
                  + f_3 * pc_x[k] * shf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, sgf_36, sgf_69, shd0_39, \
                         shd1_39, shf_66, shf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_8 * sgf_69[k]
                  + f_3 * pc_x[k] * shf_69[k];

        t_100[k] = f_10 * sgf_36[k]
                   + f_1 * shd0_39[k]
                   - f_2 * shd1_39[k]
                   + f_3 * pc_y[k] * shf_66[k];

        t_101[k] = f_3 * pc_z[k] * shf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, sgg0_45, sgf_38, \
                         sgf_39, sgg1_45, shd0_41, shd1_41, shf_68, \
                         shf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_10 * sgf_38[k]
                   + f_4 * shd0_41[k]
                   - f_5 * shd1_41[k]
                   + f_3 * pc_y[k] * shf_68[k];

        t_103[k] = f_10 * sgf_39[k]
                   + f_3 * pc_y[k] * shf_69[k];

        t_104[k] = f_1 * shd0_41[k]
                   - f_2 * shd1_41[k]
                   + f_3 * pc_z[k] * shf_69[k];

        t_105[k] = pb_z[k] * sgg0_45[k]
                   - f_6 * pc_z[k] * sgg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, sgg0_48, sgf_30, \
                         sgf_40, sgf_42, sgg1_48, shf_70, shf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * sgf_40[k]
                   + f_3 * pc_y[k] * shf_70[k];

        t_107[k] = f_7 * sgf_30[k]
                   + f_3 * pc_z[k] * shf_70[k];

        t_108[k] = pb_z[k] * sgg0_48[k]
                   - f_6 * pc_z[k] * sgg1_48[k];

        t_109[k] = f_8 * sgf_42[k]
                   + f_3 * pc_y[k] * shf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, sgf_75, sgf_76, sgf_77, sgf_78, \
                         shd0_47, shd1_47, shf_75, shf_76, shf_77, \
                         shf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_8 * sgf_75[k]
                   + f_4 * shd0_47[k]
                   - f_5 * shd1_47[k]
                   + f_3 * pc_x[k] * shf_75[k];

        t_111[k] = f_8 * sgf_76[k]
                   + f_3 * pc_x[k] * shf_76[k];

        t_112[k] = f_8 * sgf_77[k]
                   + f_3 * pc_x[k] * shf_77[k];

        t_113[k] = f_8 * sgf_78[k]
                   + f_3 * pc_x[k] * shf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, sgg0_55, sgf_36, sgf_79, \
                         sgg1_55, shf_76, shf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_8 * sgf_79[k]
                   + f_3 * pc_x[k] * shf_79[k];

        t_115[k] = pb_z[k] * sgg0_55[k]
                   - f_6 * pc_z[k] * sgg1_55[k];

        t_116[k] = f_7 * sgf_36[k]
                   + f_3 * pc_z[k] * shf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, sgg0_75, sgf_39, \
                         sgf_48, sgf_49, sgg1_75, shd0_47, shd1_47, shf_78, \
                         shf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * sgf_48[k]
                   + f_4 * shd0_47[k]
                   - f_5 * shd1_47[k]
                   + f_3 * pc_y[k] * shf_78[k];

        t_118[k] = f_8 * sgf_49[k]
                   + f_3 * pc_y[k] * shf_79[k];

        t_119[k] = f_7 * sgf_39[k]
                   + f_1 * shd0_47[k]
                   - f_2 * shd1_47[k]
                   + f_3 * pc_z[k] * shf_79[k];

        t_120[k] = pb_y[k] * sgg0_75[k]
                   - f_6 * pc_y[k] * sgg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, sgg0_78, sgf_40, \
                         sgf_50, sgf_51, sgf_52, sgg1_78, shf_80, \
                         shf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * sgf_50[k]
                   + f_3 * pc_y[k] * shf_80[k];

        t_122[k] = f_8 * sgf_40[k]
                   + f_3 * pc_z[k] * shf_80[k];

        t_123[k] = pb_y[k] * sgg0_78[k]
                   + f_8 * sgf_51[k]
                   - f_6 * pc_y[k] * sgg1_78[k];

        t_124[k] = f_7 * sgf_52[k]
                   + f_3 * pc_y[k] * shf_82[k];
    }
}

static auto
compute_prim_shg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgg0,
                                                          const size_t sgf, const size_t sgg1,
                                                          const size_t shd0, const size_t shd1,
                                                          const size_t shf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgg0_80 = buffer.data(sgg0 + 80);
    const auto *sgg0_89 = buffer.data(sgg0 + 89);
    const auto *sgg0_90 = buffer.data(sgg0 + 90);
    const auto *sgg0_93 = buffer.data(sgg0 + 93);
    const auto *sgg0_135 = buffer.data(sgg0 + 135);
    const auto *sgg0_140 = buffer.data(sgg0 + 140);
    const auto *sgg0_150 = buffer.data(sgg0 + 150);
    const auto *sgg0_153 = buffer.data(sgg0 + 153);
    const auto *sgg0_155 = buffer.data(sgg0 + 155);
    const auto *sgg0_160 = buffer.data(sgg0 + 160);
    const auto *sgg0_162 = buffer.data(sgg0 + 162);
    const auto *sgg0_164 = buffer.data(sgg0 + 164);
    const auto *sgg0_170 = buffer.data(sgg0 + 170);
    const auto *sgg0_175 = buffer.data(sgg0 + 175);
    const auto *sgg0_177 = buffer.data(sgg0 + 177);
    const auto *sgg0_179 = buffer.data(sgg0 + 179);
    const auto *sgg0_180 = buffer.data(sgg0 + 180);
    const auto *sgg0_183 = buffer.data(sgg0 + 183);
    const auto *sgg0_185 = buffer.data(sgg0 + 185);
    const auto *sgg0_190 = buffer.data(sgg0 + 190);
    const auto *sgg0_192 = buffer.data(sgg0 + 192);
    const auto *sgg0_194 = buffer.data(sgg0 + 194);
    const auto *sgg0_198 = buffer.data(sgg0 + 198);
    const auto *sgg0_205 = buffer.data(sgg0 + 205);
    const auto *sgg0_207 = buffer.data(sgg0 + 207);
    const auto *sgg0_209 = buffer.data(sgg0 + 209);
    const auto *sgg0_210 = buffer.data(sgg0 + 210);
    const auto *sgg0_213 = buffer.data(sgg0 + 213);
    const auto *sgg0_215 = buffer.data(sgg0 + 215);
    const auto *sgg0_220 = buffer.data(sgg0 + 220);
    const auto *sgg0_222 = buffer.data(sgg0 + 222);
    const auto *sgg0_224 = buffer.data(sgg0 + 224);

    const auto *sgf_46 = buffer.data(sgf + 46);
    const auto *sgf_50 = buffer.data(sgf + 50);
    const auto *sgf_56 = buffer.data(sgf + 56);
    const auto *sgf_58 = buffer.data(sgf + 58);
    const auto *sgf_59 = buffer.data(sgf + 59);
    const auto *sgf_60 = buffer.data(sgf + 60);
    const auto *sgf_62 = buffer.data(sgf + 62);
    const auto *sgf_66 = buffer.data(sgf + 66);
    const auto *sgf_69 = buffer.data(sgf + 69);
    const auto *sgf_70 = buffer.data(sgf + 70);
    const auto *sgf_72 = buffer.data(sgf + 72);
    const auto *sgf_76 = buffer.data(sgf + 76);
    const auto *sgf_79 = buffer.data(sgf + 79);
    const auto *sgf_80 = buffer.data(sgf + 80);
    const auto *sgf_82 = buffer.data(sgf + 82);
    const auto *sgf_86 = buffer.data(sgf + 86);
    const auto *sgf_87 = buffer.data(sgf + 87);
    const auto *sgf_88 = buffer.data(sgf + 88);
    const auto *sgf_89 = buffer.data(sgf + 89);
    const auto *sgf_90 = buffer.data(sgf + 90);
    const auto *sgf_92 = buffer.data(sgf + 92);
    const auto *sgf_93 = buffer.data(sgf + 93);
    const auto *sgf_95 = buffer.data(sgf + 95);
    const auto *sgf_96 = buffer.data(sgf + 96);
    const auto *sgf_97 = buffer.data(sgf + 97);
    const auto *sgf_98 = buffer.data(sgf + 98);
    const auto *sgf_99 = buffer.data(sgf + 99);
    const auto *sgf_100 = buffer.data(sgf + 100);
    const auto *sgf_102 = buffer.data(sgf + 102);
    const auto *sgf_103 = buffer.data(sgf + 103);
    const auto *sgf_105 = buffer.data(sgf + 105);
    const auto *sgf_106 = buffer.data(sgf + 106);
    const auto *sgf_107 = buffer.data(sgf + 107);
    const auto *sgf_108 = buffer.data(sgf + 108);
    const auto *sgf_109 = buffer.data(sgf + 109);
    const auto *sgf_110 = buffer.data(sgf + 110);
    const auto *sgf_112 = buffer.data(sgf + 112);
    const auto *sgf_115 = buffer.data(sgf + 115);
    const auto *sgf_116 = buffer.data(sgf + 116);
    const auto *sgf_117 = buffer.data(sgf + 117);
    const auto *sgf_118 = buffer.data(sgf + 118);
    const auto *sgf_119 = buffer.data(sgf + 119);
    const auto *sgf_120 = buffer.data(sgf + 120);
    const auto *sgf_123 = buffer.data(sgf + 123);
    const auto *sgf_125 = buffer.data(sgf + 125);
    const auto *sgf_126 = buffer.data(sgf + 126);
    const auto *sgf_127 = buffer.data(sgf + 127);
    const auto *sgf_128 = buffer.data(sgf + 128);
    const auto *sgf_129 = buffer.data(sgf + 129);
    const auto *sgf_133 = buffer.data(sgf + 133);
    const auto *sgf_136 = buffer.data(sgf + 136);
    const auto *sgf_137 = buffer.data(sgf + 137);
    const auto *sgf_138 = buffer.data(sgf + 138);
    const auto *sgf_139 = buffer.data(sgf + 139);
    const auto *sgf_140 = buffer.data(sgf + 140);
    const auto *sgf_143 = buffer.data(sgf + 143);
    const auto *sgf_145 = buffer.data(sgf + 145);
    const auto *sgf_146 = buffer.data(sgf + 146);
    const auto *sgf_147 = buffer.data(sgf + 147);
    const auto *sgf_148 = buffer.data(sgf + 148);
    const auto *sgf_149 = buffer.data(sgf + 149);

    const auto *sgg1_80 = buffer.data(sgg1 + 80);
    const auto *sgg1_89 = buffer.data(sgg1 + 89);
    const auto *sgg1_90 = buffer.data(sgg1 + 90);
    const auto *sgg1_93 = buffer.data(sgg1 + 93);
    const auto *sgg1_135 = buffer.data(sgg1 + 135);
    const auto *sgg1_140 = buffer.data(sgg1 + 140);
    const auto *sgg1_150 = buffer.data(sgg1 + 150);
    const auto *sgg1_153 = buffer.data(sgg1 + 153);
    const auto *sgg1_155 = buffer.data(sgg1 + 155);
    const auto *sgg1_160 = buffer.data(sgg1 + 160);
    const auto *sgg1_162 = buffer.data(sgg1 + 162);
    const auto *sgg1_164 = buffer.data(sgg1 + 164);
    const auto *sgg1_170 = buffer.data(sgg1 + 170);
    const auto *sgg1_175 = buffer.data(sgg1 + 175);
    const auto *sgg1_177 = buffer.data(sgg1 + 177);
    const auto *sgg1_179 = buffer.data(sgg1 + 179);
    const auto *sgg1_180 = buffer.data(sgg1 + 180);
    const auto *sgg1_183 = buffer.data(sgg1 + 183);
    const auto *sgg1_185 = buffer.data(sgg1 + 185);
    const auto *sgg1_190 = buffer.data(sgg1 + 190);
    const auto *sgg1_192 = buffer.data(sgg1 + 192);
    const auto *sgg1_194 = buffer.data(sgg1 + 194);
    const auto *sgg1_198 = buffer.data(sgg1 + 198);
    const auto *sgg1_205 = buffer.data(sgg1 + 205);
    const auto *sgg1_207 = buffer.data(sgg1 + 207);
    const auto *sgg1_209 = buffer.data(sgg1 + 209);
    const auto *sgg1_210 = buffer.data(sgg1 + 210);
    const auto *sgg1_213 = buffer.data(sgg1 + 213);
    const auto *sgg1_215 = buffer.data(sgg1 + 215);
    const auto *sgg1_220 = buffer.data(sgg1 + 220);
    const auto *sgg1_222 = buffer.data(sgg1 + 222);
    const auto *sgg1_224 = buffer.data(sgg1 + 224);

    const auto *shd0_51 = buffer.data(shd0 + 51);
    const auto *shd0_53 = buffer.data(shd0 + 53);
    const auto *shd0_54 = buffer.data(shd0 + 54);
    const auto *shd0_57 = buffer.data(shd0 + 57);
    const auto *shd0_59 = buffer.data(shd0 + 59);
    const auto *shd0_90 = buffer.data(shd0 + 90);
    const auto *shd0_93 = buffer.data(shd0 + 93);
    const auto *shd0_95 = buffer.data(shd0 + 95);
    const auto *shd0_101 = buffer.data(shd0 + 101);

    const auto *shd1_51 = buffer.data(shd1 + 51);
    const auto *shd1_53 = buffer.data(shd1 + 53);
    const auto *shd1_54 = buffer.data(shd1 + 54);
    const auto *shd1_57 = buffer.data(shd1 + 57);
    const auto *shd1_59 = buffer.data(shd1 + 59);
    const auto *shd1_90 = buffer.data(shd1 + 90);
    const auto *shd1_93 = buffer.data(shd1 + 93);
    const auto *shd1_95 = buffer.data(shd1 + 95);
    const auto *shd1_101 = buffer.data(shd1 + 101);

    const auto *shf_86 = buffer.data(shf + 86);
    const auto *shf_87 = buffer.data(shf + 87);
    const auto *shf_88 = buffer.data(shf + 88);
    const auto *shf_89 = buffer.data(shf + 89);
    const auto *shf_90 = buffer.data(shf + 90);
    const auto *shf_92 = buffer.data(shf + 92);
    const auto *shf_93 = buffer.data(shf + 93);
    const auto *shf_95 = buffer.data(shf + 95);
    const auto *shf_96 = buffer.data(shf + 96);
    const auto *shf_97 = buffer.data(shf + 97);
    const auto *shf_98 = buffer.data(shf + 98);
    const auto *shf_99 = buffer.data(shf + 99);
    const auto *shf_100 = buffer.data(shf + 100);
    const auto *shf_102 = buffer.data(shf + 102);
    const auto *shf_106 = buffer.data(shf + 106);
    const auto *shf_107 = buffer.data(shf + 107);
    const auto *shf_108 = buffer.data(shf + 108);
    const auto *shf_109 = buffer.data(shf + 109);
    const auto *shf_110 = buffer.data(shf + 110);
    const auto *shf_112 = buffer.data(shf + 112);
    const auto *shf_116 = buffer.data(shf + 116);
    const auto *shf_117 = buffer.data(shf + 117);
    const auto *shf_118 = buffer.data(shf + 118);
    const auto *shf_119 = buffer.data(shf + 119);
    const auto *shf_120 = buffer.data(shf + 120);
    const auto *shf_122 = buffer.data(shf + 122);
    const auto *shf_126 = buffer.data(shf + 126);
    const auto *shf_127 = buffer.data(shf + 127);
    const auto *shf_128 = buffer.data(shf + 128);
    const auto *shf_129 = buffer.data(shf + 129);
    const auto *shf_130 = buffer.data(shf + 130);
    const auto *shf_132 = buffer.data(shf + 132);
    const auto *shf_136 = buffer.data(shf + 136);
    const auto *shf_137 = buffer.data(shf + 137);
    const auto *shf_138 = buffer.data(shf + 138);
    const auto *shf_139 = buffer.data(shf + 139);
    const auto *shf_140 = buffer.data(shf + 140);
    const auto *shf_142 = buffer.data(shf + 142);
    const auto *shf_146 = buffer.data(shf + 146);
    const auto *shf_147 = buffer.data(shf + 147);
    const auto *shf_148 = buffer.data(shf + 148);
    const auto *shf_149 = buffer.data(shf + 149);
    const auto *shf_150 = buffer.data(shf + 150);
    const auto *shf_152 = buffer.data(shf + 152);
    const auto *shf_153 = buffer.data(shf + 153);
    const auto *shf_155 = buffer.data(shf + 155);
    const auto *shf_156 = buffer.data(shf + 156);
    const auto *shf_157 = buffer.data(shf + 157);
    const auto *shf_158 = buffer.data(shf + 158);
    const auto *shf_159 = buffer.data(shf + 159);
    const auto *shf_160 = buffer.data(shf + 160);
    const auto *shf_162 = buffer.data(shf + 162);
    const auto *shf_165 = buffer.data(shf + 165);
    const auto *shf_166 = buffer.data(shf + 166);
    const auto *shf_167 = buffer.data(shf + 167);
    const auto *shf_168 = buffer.data(shf + 168);
    const auto *shf_169 = buffer.data(shf + 169);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, sgg0_80, sgf_86, \
                         sgf_87, sgf_88, sgg1_80, shf_86, shf_87, \
                         shf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * sgg0_80[k]
                   - f_6 * pc_y[k] * sgg1_80[k];

        t_126[k] = f_8 * sgf_86[k]
                   + f_3 * pc_x[k] * shf_86[k];

        t_127[k] = f_8 * sgf_87[k]
                   + f_3 * pc_x[k] * shf_87[k];

        t_128[k] = f_8 * sgf_88[k]
                   + f_3 * pc_x[k] * shf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, sgf_46, sgf_56, sgf_89, \
                         shd0_51, shd1_51, shf_86, shf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_8 * sgf_89[k]
                   + f_3 * pc_x[k] * shf_89[k];

        t_130[k] = f_7 * sgf_56[k]
                   + f_1 * shd0_51[k]
                   - f_2 * shd1_51[k]
                   + f_3 * pc_y[k] * shf_86[k];

        t_131[k] = f_8 * sgf_46[k]
                   + f_3 * pc_z[k] * shf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, sgg0_89, sgf_58, sgf_59, sgg1_89, \
                         shd0_53, shd1_53, shf_88, shf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * sgf_58[k]
                   + f_4 * shd0_53[k]
                   - f_5 * shd1_53[k]
                   + f_3 * pc_y[k] * shf_88[k];

        t_133[k] = f_7 * sgf_59[k]
                   + f_3 * pc_y[k] * shf_89[k];

        t_134[k] = pb_y[k] * sgg0_89[k]
                   - f_6 * pc_y[k] * sgg1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, sgf_50, sgf_90, sgf_93, \
                         shd0_54, shd0_57, shd1_54, shd1_57, shf_90, \
                         shf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_8 * sgf_90[k]
                   + f_1 * shd0_54[k]
                   - f_2 * shd1_54[k]
                   + f_3 * pc_x[k] * shf_90[k];

        t_136[k] = f_3 * pc_y[k] * shf_90[k];

        t_137[k] = f_10 * sgf_50[k]
                   + f_3 * pc_z[k] * shf_90[k];

        t_138[k] = f_8 * sgf_93[k]
                   + f_4 * shd0_57[k]
                   - f_5 * shd1_57[k]
                   + f_3 * pc_x[k] * shf_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, sgf_95, sgf_96, sgf_97, \
                         shd0_59, shd1_59, shf_92, shf_95, shf_96, \
                         shf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * shf_92[k];

        t_140[k] = f_8 * sgf_95[k]
                   + f_4 * shd0_59[k]
                   - f_5 * shd1_59[k]
                   + f_3 * pc_x[k] * shf_95[k];

        t_141[k] = f_8 * sgf_96[k]
                   + f_3 * pc_x[k] * shf_96[k];

        t_142[k] = f_8 * sgf_97[k]
                   + f_3 * pc_x[k] * shf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, sgf_56, sgf_98, sgf_99, \
                         shd0_57, shd1_57, shf_96, shf_98, shf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_8 * sgf_98[k]
                   + f_3 * pc_x[k] * shf_98[k];

        t_144[k] = f_8 * sgf_99[k]
                   + f_3 * pc_x[k] * shf_99[k];

        t_145[k] = f_1 * shd0_57[k]
                   - f_2 * shd1_57[k]
                   + f_3 * pc_y[k] * shf_96[k];

        t_146[k] = f_10 * sgf_56[k]
                   + f_3 * pc_z[k] * shf_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pb_x, pc_x, pc_y, pc_z, sgg0_150, sgf_59, \
                         sgf_100, sgg1_150, shd0_59, shd1_59, shf_98, \
                         shf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * shd0_59[k]
                   - f_5 * shd1_59[k]
                   + f_3 * pc_y[k] * shf_98[k];

        t_148[k] = f_3 * pc_y[k] * shf_99[k];

        t_149[k] = f_10 * sgf_59[k]
                   + f_1 * shd0_59[k]
                   - f_2 * shd1_59[k]
                   + f_3 * pc_z[k] * shf_99[k];

        t_150[k] = pb_x[k] * sgg0_150[k]
                   + f_9 * sgf_100[k]
                   - f_6 * pc_x[k] * sgg1_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pb_x, pc_x, pc_y, pc_z, sgg0_153, sgf_60, \
                         sgf_62, sgf_103, sgg1_153, shf_100, shf_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * sgf_60[k]
                   + f_3 * pc_y[k] * shf_100[k];

        t_152[k] = f_3 * pc_z[k] * shf_100[k];

        t_153[k] = pb_x[k] * sgg0_153[k]
                   + f_8 * sgf_103[k]
                   - f_6 * pc_x[k] * sgg1_153[k];

        t_154[k] = f_9 * sgf_62[k]
                   + f_3 * pc_y[k] * shf_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pb_x, pc_x, sgg0_155, sgf_105, sgf_106, \
                         sgf_107, sgf_108, sgg1_155, shf_106, shf_107, \
                         shf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = pb_x[k] * sgg0_155[k]
                   + f_8 * sgf_105[k]
                   - f_6 * pc_x[k] * sgg1_155[k];

        t_156[k] = f_7 * sgf_106[k]
                   + f_3 * pc_x[k] * shf_106[k];

        t_157[k] = f_7 * sgf_107[k]
                   + f_3 * pc_x[k] * shf_107[k];

        t_158[k] = f_7 * sgf_108[k]
                   + f_3 * pc_x[k] * shf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pb_x, pc_x, pc_z, sgg0_160, sgg0_162, \
                         sgf_109, sgg1_160, sgg1_162, shf_106, \
                         shf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_7 * sgf_109[k]
                   + f_3 * pc_x[k] * shf_109[k];

        t_160[k] = pb_x[k] * sgg0_160[k]
                   - f_6 * pc_x[k] * sgg1_160[k];

        t_161[k] = f_3 * pc_z[k] * shf_106[k];

        t_162[k] = pb_x[k] * sgg0_162[k]
                   - f_6 * pc_x[k] * sgg1_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pb_x, pb_z, pc_x, pc_y, pc_z, sgg0_90, sgg0_164, \
                         sgf_69, sgg1_90, sgg1_164, shf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * sgf_69[k]
                   + f_3 * pc_y[k] * shf_109[k];

        t_164[k] = pb_x[k] * sgg0_164[k]
                   - f_6 * pc_x[k] * sgg1_164[k];

        t_165[k] = pb_z[k] * sgg0_90[k]
                   - f_6 * pc_z[k] * sgg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, sgg0_93, sgf_60, \
                         sgf_70, sgf_72, sgg1_93, shf_110, shf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_10 * sgf_70[k]
                   + f_3 * pc_y[k] * shf_110[k];

        t_167[k] = f_7 * sgf_60[k]
                   + f_3 * pc_z[k] * shf_110[k];

        t_168[k] = pb_z[k] * sgg0_93[k]
                   - f_6 * pc_z[k] * sgg1_93[k];

        t_169[k] = f_10 * sgf_72[k]
                   + f_3 * pc_y[k] * shf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_x, pc_x, sgg0_170, sgf_115, sgf_116, \
                         sgf_117, sgf_118, sgg1_170, shf_116, shf_117, \
                         shf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_x[k] * sgg0_170[k]
                   + f_8 * sgf_115[k]
                   - f_6 * pc_x[k] * sgg1_170[k];

        t_171[k] = f_7 * sgf_116[k]
                   + f_3 * pc_x[k] * shf_116[k];

        t_172[k] = f_7 * sgf_117[k]
                   + f_3 * pc_x[k] * shf_117[k];

        t_173[k] = f_7 * sgf_118[k]
                   + f_3 * pc_x[k] * shf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pb_x, pc_x, pc_z, sgg0_175, sgg0_177, \
                         sgf_66, sgf_119, sgg1_175, sgg1_177, shf_116, \
                         shf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_7 * sgf_119[k]
                   + f_3 * pc_x[k] * shf_119[k];

        t_175[k] = pb_x[k] * sgg0_175[k]
                   - f_6 * pc_x[k] * sgg1_175[k];

        t_176[k] = f_7 * sgf_66[k]
                   + f_3 * pc_z[k] * shf_116[k];

        t_177[k] = pb_x[k] * sgg0_177[k]
                   - f_6 * pc_x[k] * sgg1_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pb_x, pc_x, pc_y, sgg0_179, sgg0_180, \
                         sgf_79, sgf_80, sgf_120, sgg1_179, sgg1_180, shf_119, \
                         shf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_10 * sgf_79[k]
                   + f_3 * pc_y[k] * shf_119[k];

        t_179[k] = pb_x[k] * sgg0_179[k]
                   - f_6 * pc_x[k] * sgg1_179[k];

        t_180[k] = pb_x[k] * sgg0_180[k]
                   + f_9 * sgf_120[k]
                   - f_6 * pc_x[k] * sgg1_180[k];

        t_181[k] = f_8 * sgf_80[k]
                   + f_3 * pc_y[k] * shf_120[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pb_x, pc_x, pc_y, pc_z, sgg0_183, sgf_70, \
                         sgf_82, sgf_123, sgg1_183, shf_120, shf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_8 * sgf_70[k]
                   + f_3 * pc_z[k] * shf_120[k];

        t_183[k] = pb_x[k] * sgg0_183[k]
                   + f_8 * sgf_123[k]
                   - f_6 * pc_x[k] * sgg1_183[k];

        t_184[k] = f_8 * sgf_82[k]
                   + f_3 * pc_y[k] * shf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_x, pc_x, sgg0_185, sgf_125, sgf_126, \
                         sgf_127, sgf_128, sgg1_185, shf_126, shf_127, \
                         shf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_x[k] * sgg0_185[k]
                   + f_8 * sgf_125[k]
                   - f_6 * pc_x[k] * sgg1_185[k];

        t_186[k] = f_7 * sgf_126[k]
                   + f_3 * pc_x[k] * shf_126[k];

        t_187[k] = f_7 * sgf_127[k]
                   + f_3 * pc_x[k] * shf_127[k];

        t_188[k] = f_7 * sgf_128[k]
                   + f_3 * pc_x[k] * shf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pb_x, pc_x, pc_z, sgg0_190, sgg0_192, \
                         sgf_76, sgf_129, sgg1_190, sgg1_192, shf_126, \
                         shf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_7 * sgf_129[k]
                   + f_3 * pc_x[k] * shf_129[k];

        t_190[k] = pb_x[k] * sgg0_190[k]
                   - f_6 * pc_x[k] * sgg1_190[k];

        t_191[k] = f_8 * sgf_76[k]
                   + f_3 * pc_z[k] * shf_126[k];

        t_192[k] = pb_x[k] * sgg0_192[k]
                   - f_6 * pc_x[k] * sgg1_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_x, pb_y, pc_x, pc_y, sgg0_135, \
                         sgg0_194, sgf_89, sgf_90, sgg1_135, sgg1_194, shf_129, \
                         shf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * sgf_89[k]
                   + f_3 * pc_y[k] * shf_129[k];

        t_194[k] = pb_x[k] * sgg0_194[k]
                   - f_6 * pc_x[k] * sgg1_194[k];

        t_195[k] = pb_y[k] * sgg0_135[k]
                   - f_6 * pc_y[k] * sgg1_135[k];

        t_196[k] = f_7 * sgf_90[k]
                   + f_3 * pc_y[k] * shf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pb_x, pc_x, pc_y, pc_z, sgg0_198, sgf_80, \
                         sgf_92, sgf_133, sgg1_198, shf_130, shf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_10 * sgf_80[k]
                   + f_3 * pc_z[k] * shf_130[k];

        t_198[k] = pb_x[k] * sgg0_198[k]
                   + f_8 * sgf_133[k]
                   - f_6 * pc_x[k] * sgg1_198[k];

        t_199[k] = f_7 * sgf_92[k]
                   + f_3 * pc_y[k] * shf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pb_y, pc_x, pc_y, sgg0_140, sgf_136, \
                         sgf_137, sgf_138, sgg1_140, shf_136, shf_137, \
                         shf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_y[k] * sgg0_140[k]
                   - f_6 * pc_y[k] * sgg1_140[k];

        t_201[k] = f_7 * sgf_136[k]
                   + f_3 * pc_x[k] * shf_136[k];

        t_202[k] = f_7 * sgf_137[k]
                   + f_3 * pc_x[k] * shf_137[k];

        t_203[k] = f_7 * sgf_138[k]
                   + f_3 * pc_x[k] * shf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_x, pc_x, pc_z, sgg0_205, sgg0_207, \
                         sgf_86, sgf_139, sgg1_205, sgg1_207, shf_136, \
                         shf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_7 * sgf_139[k]
                   + f_3 * pc_x[k] * shf_139[k];

        t_205[k] = pb_x[k] * sgg0_205[k]
                   - f_6 * pc_x[k] * sgg1_205[k];

        t_206[k] = f_10 * sgf_86[k]
                   + f_3 * pc_z[k] * shf_136[k];

        t_207[k] = pb_x[k] * sgg0_207[k]
                   - f_6 * pc_x[k] * sgg1_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_x, pc_x, pc_y, sgg0_209, sgg0_210, \
                         sgf_99, sgf_140, sgg1_209, sgg1_210, shf_139, \
                         shf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * sgf_99[k]
                   + f_3 * pc_y[k] * shf_139[k];

        t_209[k] = pb_x[k] * sgg0_209[k]
                   - f_6 * pc_x[k] * sgg1_209[k];

        t_210[k] = pb_x[k] * sgg0_210[k]
                   + f_9 * sgf_140[k]
                   - f_6 * pc_x[k] * sgg1_210[k];

        t_211[k] = f_3 * pc_y[k] * shf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pb_x, pc_x, pc_y, pc_z, sgg0_213, sgf_90, \
                         sgf_143, sgg1_213, shf_140, shf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_9 * sgf_90[k]
                   + f_3 * pc_z[k] * shf_140[k];

        t_213[k] = pb_x[k] * sgg0_213[k]
                   + f_8 * sgf_143[k]
                   - f_6 * pc_x[k] * sgg1_213[k];

        t_214[k] = f_3 * pc_y[k] * shf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pb_x, pc_x, sgg0_215, sgf_145, sgf_146, \
                         sgf_147, sgf_148, sgg1_215, shf_146, shf_147, \
                         shf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pb_x[k] * sgg0_215[k]
                   + f_8 * sgf_145[k]
                   - f_6 * pc_x[k] * sgg1_215[k];

        t_216[k] = f_7 * sgf_146[k]
                   + f_3 * pc_x[k] * shf_146[k];

        t_217[k] = f_7 * sgf_147[k]
                   + f_3 * pc_x[k] * shf_147[k];

        t_218[k] = f_7 * sgf_148[k]
                   + f_3 * pc_x[k] * shf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pb_x, pc_x, pc_z, sgg0_220, sgg0_222, \
                         sgf_96, sgf_149, sgg1_220, sgg1_222, shf_146, \
                         shf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_7 * sgf_149[k]
                   + f_3 * pc_x[k] * shf_149[k];

        t_220[k] = pb_x[k] * sgg0_220[k]
                   - f_6 * pc_x[k] * sgg1_220[k];

        t_221[k] = f_9 * sgf_96[k]
                   + f_3 * pc_z[k] * shf_146[k];

        t_222[k] = pb_x[k] * sgg0_222[k]
                   - f_6 * pc_x[k] * sgg1_222[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, t_227, pb_x, pc_x, pc_y, pc_z, sgg0_224, \
                         sgf_100, sgg1_224, shd0_90, shd1_90, shf_149, \
                         shf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * shf_149[k];

        t_224[k] = pb_x[k] * sgg0_224[k]
                   - f_6 * pc_x[k] * sgg1_224[k];

        t_225[k] = f_1 * shd0_90[k]
                   - f_2 * shd1_90[k]
                   + f_3 * pc_x[k] * shf_150[k];

        t_226[k] = f_0 * sgf_100[k]
                   + f_3 * pc_y[k] * shf_150[k];

        t_227[k] = f_3 * pc_z[k] * shf_150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pc_x, pc_y, sgf_102, shd0_93, shd0_95, \
                         shd1_93, shd1_95, shf_152, shf_153, shf_155, \
                         shf_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_4 * shd0_93[k]
                   - f_5 * shd1_93[k]
                   + f_3 * pc_x[k] * shf_153[k];

        t_229[k] = f_0 * sgf_102[k]
                   + f_3 * pc_y[k] * shf_152[k];

        t_230[k] = f_4 * shd0_95[k]
                   - f_5 * shd1_95[k]
                   + f_3 * pc_x[k] * shf_155[k];

        t_231[k] = f_3 * pc_x[k] * shf_156[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, t_236, pc_x, pc_y, pc_z, sgf_106, \
                         shd0_93, shd1_93, shf_156, shf_157, shf_158, \
                         shf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_3 * pc_x[k] * shf_157[k];

        t_233[k] = f_3 * pc_x[k] * shf_158[k];

        t_234[k] = f_3 * pc_x[k] * shf_159[k];

        t_235[k] = f_0 * sgf_106[k]
                   + f_1 * shd0_93[k]
                   - f_2 * shd1_93[k]
                   + f_3 * pc_y[k] * shf_156[k];

        t_236[k] = f_3 * pc_z[k] * shf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, sgg0_150, sgf_108, \
                         sgf_109, sgg1_150, shd0_95, shd1_95, shf_158, \
                         shf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_0 * sgf_108[k]
                   + f_4 * shd0_95[k]
                   - f_5 * shd1_95[k]
                   + f_3 * pc_y[k] * shf_158[k];

        t_238[k] = f_0 * sgf_109[k]
                   + f_3 * pc_y[k] * shf_159[k];

        t_239[k] = f_1 * shd0_95[k]
                   - f_2 * shd1_95[k]
                   + f_3 * pc_z[k] * shf_159[k];

        t_240[k] = pb_z[k] * sgg0_150[k]
                   - f_6 * pc_z[k] * sgg1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, sgg0_153, sgf_100, \
                         sgf_110, sgf_112, sgg1_153, shf_160, shf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_9 * sgf_110[k]
                   + f_3 * pc_y[k] * shf_160[k];

        t_242[k] = f_7 * sgf_100[k]
                   + f_3 * pc_z[k] * shf_160[k];

        t_243[k] = pb_z[k] * sgg0_153[k]
                   - f_6 * pc_z[k] * sgg1_153[k];

        t_244[k] = f_9 * sgf_112[k]
                   + f_3 * pc_y[k] * shf_162[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, pc_x, shd0_101, shd1_101, shf_165, \
                         shf_166, shf_167, shf_168, shf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_4 * shd0_101[k]
                   - f_5 * shd1_101[k]
                   + f_3 * pc_x[k] * shf_165[k];

        t_246[k] = f_3 * pc_x[k] * shf_166[k];

        t_247[k] = f_3 * pc_x[k] * shf_167[k];

        t_248[k] = f_3 * pc_x[k] * shf_168[k];

        t_249[k] = f_3 * pc_x[k] * shf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pb_z, pc_y, pc_z, sgg0_160, sgg0_162, \
                         sgf_106, sgf_107, sgf_119, sgg1_160, sgg1_162, shf_166, \
                         shf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pb_z[k] * sgg0_160[k]
                   - f_6 * pc_z[k] * sgg1_160[k];

        t_251[k] = f_7 * sgf_106[k]
                   + f_3 * pc_z[k] * shf_166[k];

        t_252[k] = pb_z[k] * sgg0_162[k]
                   + f_8 * sgf_107[k]
                   - f_6 * pc_z[k] * sgg1_162[k];

        t_253[k] = f_9 * sgf_119[k]
                   + f_3 * pc_y[k] * shf_169[k];
    }
}

static auto
compute_prim_shg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgg0,
                                                          const size_t sgf, const size_t sgg1,
                                                          const size_t shd0, const size_t shd1,
                                                          const size_t shf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *sgg0_210 = buffer.data(sgg0 + 210);
    const auto *sgg0_215 = buffer.data(sgg0 + 215);
    const auto *sgg0_220 = buffer.data(sgg0 + 220);
    const auto *sgg0_222 = buffer.data(sgg0 + 222);
    const auto *sgg0_224 = buffer.data(sgg0 + 224);

    const auto *sgf_109 = buffer.data(sgf + 109);
    const auto *sgf_110 = buffer.data(sgf + 110);
    const auto *sgf_116 = buffer.data(sgf + 116);
    const auto *sgf_119 = buffer.data(sgf + 119);
    const auto *sgf_120 = buffer.data(sgf + 120);
    const auto *sgf_122 = buffer.data(sgf + 122);
    const auto *sgf_126 = buffer.data(sgf + 126);
    const auto *sgf_128 = buffer.data(sgf + 128);
    const auto *sgf_129 = buffer.data(sgf + 129);
    const auto *sgf_130 = buffer.data(sgf + 130);
    const auto *sgf_132 = buffer.data(sgf + 132);
    const auto *sgf_136 = buffer.data(sgf + 136);
    const auto *sgf_138 = buffer.data(sgf + 138);
    const auto *sgf_139 = buffer.data(sgf + 139);
    const auto *sgf_140 = buffer.data(sgf + 140);
    const auto *sgf_142 = buffer.data(sgf + 142);
    const auto *sgf_146 = buffer.data(sgf + 146);
    const auto *sgf_148 = buffer.data(sgf + 148);
    const auto *sgf_149 = buffer.data(sgf + 149);

    const auto *sgg1_210 = buffer.data(sgg1 + 210);
    const auto *sgg1_215 = buffer.data(sgg1 + 215);
    const auto *sgg1_220 = buffer.data(sgg1 + 220);
    const auto *sgg1_222 = buffer.data(sgg1 + 222);
    const auto *sgg1_224 = buffer.data(sgg1 + 224);

    const auto *shd0_101 = buffer.data(shd0 + 101);
    const auto *shd0_102 = buffer.data(shd0 + 102);
    const auto *shd0_105 = buffer.data(shd0 + 105);
    const auto *shd0_107 = buffer.data(shd0 + 107);
    const auto *shd0_108 = buffer.data(shd0 + 108);
    const auto *shd0_111 = buffer.data(shd0 + 111);
    const auto *shd0_113 = buffer.data(shd0 + 113);
    const auto *shd0_117 = buffer.data(shd0 + 117);
    const auto *shd0_120 = buffer.data(shd0 + 120);
    const auto *shd0_123 = buffer.data(shd0 + 123);
    const auto *shd0_125 = buffer.data(shd0 + 125);

    const auto *shd1_101 = buffer.data(shd1 + 101);
    const auto *shd1_102 = buffer.data(shd1 + 102);
    const auto *shd1_105 = buffer.data(shd1 + 105);
    const auto *shd1_107 = buffer.data(shd1 + 107);
    const auto *shd1_108 = buffer.data(shd1 + 108);
    const auto *shd1_111 = buffer.data(shd1 + 111);
    const auto *shd1_113 = buffer.data(shd1 + 113);
    const auto *shd1_117 = buffer.data(shd1 + 117);
    const auto *shd1_120 = buffer.data(shd1 + 120);
    const auto *shd1_123 = buffer.data(shd1 + 123);
    const auto *shd1_125 = buffer.data(shd1 + 125);

    const auto *shf_169 = buffer.data(shf + 169);
    const auto *shf_170 = buffer.data(shf + 170);
    const auto *shf_172 = buffer.data(shf + 172);
    const auto *shf_173 = buffer.data(shf + 173);
    const auto *shf_175 = buffer.data(shf + 175);
    const auto *shf_176 = buffer.data(shf + 176);
    const auto *shf_177 = buffer.data(shf + 177);
    const auto *shf_178 = buffer.data(shf + 178);
    const auto *shf_179 = buffer.data(shf + 179);
    const auto *shf_180 = buffer.data(shf + 180);
    const auto *shf_182 = buffer.data(shf + 182);
    const auto *shf_183 = buffer.data(shf + 183);
    const auto *shf_185 = buffer.data(shf + 185);
    const auto *shf_186 = buffer.data(shf + 186);
    const auto *shf_187 = buffer.data(shf + 187);
    const auto *shf_188 = buffer.data(shf + 188);
    const auto *shf_189 = buffer.data(shf + 189);
    const auto *shf_190 = buffer.data(shf + 190);
    const auto *shf_192 = buffer.data(shf + 192);
    const auto *shf_193 = buffer.data(shf + 193);
    const auto *shf_196 = buffer.data(shf + 196);
    const auto *shf_197 = buffer.data(shf + 197);
    const auto *shf_198 = buffer.data(shf + 198);
    const auto *shf_199 = buffer.data(shf + 199);
    const auto *shf_200 = buffer.data(shf + 200);
    const auto *shf_202 = buffer.data(shf + 202);
    const auto *shf_203 = buffer.data(shf + 203);
    const auto *shf_205 = buffer.data(shf + 205);
    const auto *shf_206 = buffer.data(shf + 206);
    const auto *shf_207 = buffer.data(shf + 207);
    const auto *shf_208 = buffer.data(shf + 208);
    const auto *shf_209 = buffer.data(shf + 209);

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pc_x, pc_y, pc_z, sgf_109, sgf_110, \
                         sgf_120, shd0_101, shd0_102, shd1_101, shd1_102, shf_169, \
                         shf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_7 * sgf_109[k]
                   + f_1 * shd0_101[k]
                   - f_2 * shd1_101[k]
                   + f_3 * pc_z[k] * shf_169[k];

        t_255[k] = f_1 * shd0_102[k]
                   - f_2 * shd1_102[k]
                   + f_3 * pc_x[k] * shf_170[k];

        t_256[k] = f_10 * sgf_120[k]
                   + f_3 * pc_y[k] * shf_170[k];

        t_257[k] = f_8 * sgf_110[k]
                   + f_3 * pc_z[k] * shf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pc_x, pc_y, sgf_122, shd0_105, shd0_107, \
                         shd1_105, shd1_107, shf_172, shf_173, shf_175, \
                         shf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_4 * shd0_105[k]
                   - f_5 * shd1_105[k]
                   + f_3 * pc_x[k] * shf_173[k];

        t_259[k] = f_10 * sgf_122[k]
                   + f_3 * pc_y[k] * shf_172[k];

        t_260[k] = f_4 * shd0_107[k]
                   - f_5 * shd1_107[k]
                   + f_3 * pc_x[k] * shf_175[k];

        t_261[k] = f_3 * pc_x[k] * shf_176[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pc_x, pc_y, pc_z, sgf_116, \
                         sgf_126, shd0_105, shd1_105, shf_176, shf_177, shf_178, \
                         shf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_3 * pc_x[k] * shf_177[k];

        t_263[k] = f_3 * pc_x[k] * shf_178[k];

        t_264[k] = f_3 * pc_x[k] * shf_179[k];

        t_265[k] = f_10 * sgf_126[k]
                   + f_1 * shd0_105[k]
                   - f_2 * shd1_105[k]
                   + f_3 * pc_y[k] * shf_176[k];

        t_266[k] = f_8 * sgf_116[k]
                   + f_3 * pc_z[k] * shf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, sgf_119, sgf_128, sgf_129, shd0_107, \
                         shd1_107, shf_178, shf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_10 * sgf_128[k]
                   + f_4 * shd0_107[k]
                   - f_5 * shd1_107[k]
                   + f_3 * pc_y[k] * shf_178[k];

        t_268[k] = f_10 * sgf_129[k]
                   + f_3 * pc_y[k] * shf_179[k];

        t_269[k] = f_8 * sgf_119[k]
                   + f_1 * shd0_107[k]
                   - f_2 * shd1_107[k]
                   + f_3 * pc_z[k] * shf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_y, pc_z, sgf_120, sgf_130, \
                         shd0_108, shd0_111, shd1_108, shd1_111, shf_180, \
                         shf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * shd0_108[k]
                   - f_2 * shd1_108[k]
                   + f_3 * pc_x[k] * shf_180[k];

        t_271[k] = f_8 * sgf_130[k]
                   + f_3 * pc_y[k] * shf_180[k];

        t_272[k] = f_10 * sgf_120[k]
                   + f_3 * pc_z[k] * shf_180[k];

        t_273[k] = f_4 * shd0_111[k]
                   - f_5 * shd1_111[k]
                   + f_3 * pc_x[k] * shf_183[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_y, sgf_132, shd0_113, \
                         shd1_113, shf_182, shf_185, shf_186, shf_187, \
                         shf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_8 * sgf_132[k]
                   + f_3 * pc_y[k] * shf_182[k];

        t_275[k] = f_4 * shd0_113[k]
                   - f_5 * shd1_113[k]
                   + f_3 * pc_x[k] * shf_185[k];

        t_276[k] = f_3 * pc_x[k] * shf_186[k];

        t_277[k] = f_3 * pc_x[k] * shf_187[k];

        t_278[k] = f_3 * pc_x[k] * shf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, sgf_126, sgf_136, shd0_111, \
                         shd1_111, shf_186, shf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_3 * pc_x[k] * shf_189[k];

        t_280[k] = f_8 * sgf_136[k]
                   + f_1 * shd0_111[k]
                   - f_2 * shd1_111[k]
                   + f_3 * pc_y[k] * shf_186[k];

        t_281[k] = f_10 * sgf_126[k]
                   + f_3 * pc_z[k] * shf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, sgg0_210, sgf_129, \
                         sgf_138, sgf_139, sgg1_210, shd0_113, shd1_113, shf_188, \
                         shf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * sgf_138[k]
                   + f_4 * shd0_113[k]
                   - f_5 * shd1_113[k]
                   + f_3 * pc_y[k] * shf_188[k];

        t_283[k] = f_8 * sgf_139[k]
                   + f_3 * pc_y[k] * shf_189[k];

        t_284[k] = f_10 * sgf_129[k]
                   + f_1 * shd0_113[k]
                   - f_2 * shd1_113[k]
                   + f_3 * pc_z[k] * shf_189[k];

        t_285[k] = pb_y[k] * sgg0_210[k]
                   - f_6 * pc_y[k] * sgg1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pc_x, pc_y, pc_z, sgf_130, sgf_140, \
                         sgf_142, shd0_117, shd1_117, shf_190, shf_192, \
                         shf_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * sgf_140[k]
                   + f_3 * pc_y[k] * shf_190[k];

        t_287[k] = f_9 * sgf_130[k]
                   + f_3 * pc_z[k] * shf_190[k];

        t_288[k] = f_4 * shd0_117[k]
                   - f_5 * shd1_117[k]
                   + f_3 * pc_x[k] * shf_193[k];

        t_289[k] = f_7 * sgf_142[k]
                   + f_3 * pc_y[k] * shf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pb_y, pc_x, pc_y, sgg0_215, \
                         sgg1_215, shf_196, shf_197, shf_198, shf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * sgg0_215[k]
                   - f_6 * pc_y[k] * sgg1_215[k];

        t_291[k] = f_3 * pc_x[k] * shf_196[k];

        t_292[k] = f_3 * pc_x[k] * shf_197[k];

        t_293[k] = f_3 * pc_x[k] * shf_198[k];

        t_294[k] = f_3 * pc_x[k] * shf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pb_y, pc_y, pc_z, sgg0_220, sgg0_222, sgf_136, \
                         sgf_146, sgf_148, sgg1_220, sgg1_222, \
                         shf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = pb_y[k] * sgg0_220[k]
                   + f_9 * sgf_146[k]
                   - f_6 * pc_y[k] * sgg1_220[k];

        t_296[k] = f_9 * sgf_136[k]
                   + f_3 * pc_z[k] * shf_196[k];

        t_297[k] = pb_y[k] * sgg0_222[k]
                   + f_8 * sgf_148[k]
                   - f_6 * pc_y[k] * sgg1_222[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pb_y, pc_x, pc_y, sgg0_224, sgf_149, \
                         sgg1_224, shd0_120, shd1_120, shf_199, \
                         shf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * sgf_149[k]
                   + f_3 * pc_y[k] * shf_199[k];

        t_299[k] = pb_y[k] * sgg0_224[k]
                   - f_6 * pc_y[k] * sgg1_224[k];

        t_300[k] = f_1 * shd0_120[k]
                   - f_2 * shd1_120[k]
                   + f_3 * pc_x[k] * shf_200[k];

        t_301[k] = f_3 * pc_y[k] * shf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pc_x, pc_y, pc_z, sgf_140, shd0_123, \
                         shd0_125, shd1_123, shd1_125, shf_200, shf_202, shf_203, \
                         shf_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_0 * sgf_140[k]
                   + f_3 * pc_z[k] * shf_200[k];

        t_303[k] = f_4 * shd0_123[k]
                   - f_5 * shd1_123[k]
                   + f_3 * pc_x[k] * shf_203[k];

        t_304[k] = f_3 * pc_y[k] * shf_202[k];

        t_305[k] = f_4 * shd0_125[k]
                   - f_5 * shd1_125[k]
                   + f_3 * pc_x[k] * shf_205[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, sgf_146, \
                         shd0_123, shd1_123, shf_206, shf_207, shf_208, \
                         shf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_3 * pc_x[k] * shf_206[k];

        t_307[k] = f_3 * pc_x[k] * shf_207[k];

        t_308[k] = f_3 * pc_x[k] * shf_208[k];

        t_309[k] = f_3 * pc_x[k] * shf_209[k];

        t_310[k] = f_1 * shd0_123[k]
                   - f_2 * shd1_123[k]
                   + f_3 * pc_y[k] * shf_206[k];

        t_311[k] = f_0 * sgf_146[k]
                   + f_3 * pc_z[k] * shf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pc_y, pc_z, sgf_149, shd0_125, shd1_125, \
                         shf_208, shf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * shd0_125[k]
                   - f_5 * shd1_125[k]
                   + f_3 * pc_y[k] * shf_208[k];

        t_313[k] = f_3 * pc_y[k] * shf_209[k];

        t_314[k] = f_0 * sgf_149[k]
                   + f_1 * shd0_125[k]
                   - f_2 * shd1_125[k]
                   + f_3 * pc_z[k] * shf_209[k];
    }
}

auto
compute_prim_shg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgg0, const size_t sgf,
                                                   const size_t sgg1, const size_t shd0,
                                                   const size_t shd1, const size_t shf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgg0, sgf,
                                                              sgg1, shd0, shd1, shf, ncols,
                                                              gamma, p, q);

    compute_prim_shg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgg0, sgf,
                                                              sgg1, shd0, shd1, shf, ncols,
                                                              gamma, p, q);

    compute_prim_shg_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sgg0, sgf,
                                                              sgg1, shd0, shd1, shf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
