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


#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_shf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgf0,
                                                          const size_t sgd, const size_t sgf1,
                                                          const size_t shp0, const size_t shp1,
                                                          const size_t shd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 1.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_0 = buffer.data(sgf0 + 0);
    const auto *sgf0_6 = buffer.data(sgf0 + 6);
    const auto *sgf0_9 = buffer.data(sgf0 + 9);
    const auto *sgf0_16 = buffer.data(sgf0 + 16);
    const auto *sgf0_20 = buffer.data(sgf0 + 20);
    const auto *sgf0_29 = buffer.data(sgf0 + 29);
    const auto *sgf0_30 = buffer.data(sgf0 + 30);
    const auto *sgf0_36 = buffer.data(sgf0 + 36);
    const auto *sgf0_50 = buffer.data(sgf0 + 50);
    const auto *sgf0_59 = buffer.data(sgf0 + 59);
    const auto *sgf0_60 = buffer.data(sgf0 + 60);
    const auto *sgf0_100 = buffer.data(sgf0 + 100);
    const auto *sgf0_106 = buffer.data(sgf0 + 106);
    const auto *sgf0_109 = buffer.data(sgf0 + 109);
    const auto *sgf0_116 = buffer.data(sgf0 + 116);
    const auto *sgf0_119 = buffer.data(sgf0 + 119);
    const auto *sgf0_120 = buffer.data(sgf0 + 120);
    const auto *sgf0_126 = buffer.data(sgf0 + 126);
    const auto *sgf0_129 = buffer.data(sgf0 + 129);

    const auto *sgd_0 = buffer.data(sgd + 0);
    const auto *sgd_3 = buffer.data(sgd + 3);
    const auto *sgd_4 = buffer.data(sgd + 4);
    const auto *sgd_5 = buffer.data(sgd + 5);
    const auto *sgd_6 = buffer.data(sgd + 6);
    const auto *sgd_9 = buffer.data(sgd + 9);
    const auto *sgd_10 = buffer.data(sgd + 10);
    const auto *sgd_11 = buffer.data(sgd + 11);
    const auto *sgd_12 = buffer.data(sgd + 12);
    const auto *sgd_15 = buffer.data(sgd + 15);
    const auto *sgd_16 = buffer.data(sgd + 16);
    const auto *sgd_17 = buffer.data(sgd + 17);
    const auto *sgd_18 = buffer.data(sgd + 18);
    const auto *sgd_21 = buffer.data(sgd + 21);
    const auto *sgd_22 = buffer.data(sgd + 22);
    const auto *sgd_23 = buffer.data(sgd + 23);
    const auto *sgd_24 = buffer.data(sgd + 24);
    const auto *sgd_27 = buffer.data(sgd + 27);
    const auto *sgd_28 = buffer.data(sgd + 28);
    const auto *sgd_29 = buffer.data(sgd + 29);
    const auto *sgd_30 = buffer.data(sgd + 30);
    const auto *sgd_33 = buffer.data(sgd + 33);
    const auto *sgd_34 = buffer.data(sgd + 34);
    const auto *sgd_35 = buffer.data(sgd + 35);
    const auto *sgd_36 = buffer.data(sgd + 36);
    const auto *sgd_39 = buffer.data(sgd + 39);
    const auto *sgd_40 = buffer.data(sgd + 40);
    const auto *sgd_41 = buffer.data(sgd + 41);
    const auto *sgd_42 = buffer.data(sgd + 42);
    const auto *sgd_45 = buffer.data(sgd + 45);
    const auto *sgd_46 = buffer.data(sgd + 46);
    const auto *sgd_47 = buffer.data(sgd + 47);
    const auto *sgd_48 = buffer.data(sgd + 48);
    const auto *sgd_51 = buffer.data(sgd + 51);
    const auto *sgd_52 = buffer.data(sgd + 52);
    const auto *sgd_53 = buffer.data(sgd + 53);
    const auto *sgd_54 = buffer.data(sgd + 54);
    const auto *sgd_57 = buffer.data(sgd + 57);
    const auto *sgd_58 = buffer.data(sgd + 58);
    const auto *sgd_59 = buffer.data(sgd + 59);
    const auto *sgd_60 = buffer.data(sgd + 60);
    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_64 = buffer.data(sgd + 64);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_69 = buffer.data(sgd + 69);
    const auto *sgd_70 = buffer.data(sgd + 70);
    const auto *sgd_71 = buffer.data(sgd + 71);
    const auto *sgd_72 = buffer.data(sgd + 72);
    const auto *sgd_75 = buffer.data(sgd + 75);
    const auto *sgd_76 = buffer.data(sgd + 76);
    const auto *sgd_77 = buffer.data(sgd + 77);

    const auto *sgf1_0 = buffer.data(sgf1 + 0);
    const auto *sgf1_6 = buffer.data(sgf1 + 6);
    const auto *sgf1_9 = buffer.data(sgf1 + 9);
    const auto *sgf1_16 = buffer.data(sgf1 + 16);
    const auto *sgf1_20 = buffer.data(sgf1 + 20);
    const auto *sgf1_29 = buffer.data(sgf1 + 29);
    const auto *sgf1_30 = buffer.data(sgf1 + 30);
    const auto *sgf1_36 = buffer.data(sgf1 + 36);
    const auto *sgf1_50 = buffer.data(sgf1 + 50);
    const auto *sgf1_59 = buffer.data(sgf1 + 59);
    const auto *sgf1_60 = buffer.data(sgf1 + 60);
    const auto *sgf1_100 = buffer.data(sgf1 + 100);
    const auto *sgf1_106 = buffer.data(sgf1 + 106);
    const auto *sgf1_109 = buffer.data(sgf1 + 109);
    const auto *sgf1_116 = buffer.data(sgf1 + 116);
    const auto *sgf1_119 = buffer.data(sgf1 + 119);
    const auto *sgf1_120 = buffer.data(sgf1 + 120);
    const auto *sgf1_126 = buffer.data(sgf1 + 126);
    const auto *sgf1_129 = buffer.data(sgf1 + 129);

    const auto *shp0_0 = buffer.data(shp0 + 0);
    const auto *shp0_1 = buffer.data(shp0 + 1);
    const auto *shp0_2 = buffer.data(shp0 + 2);
    const auto *shp0_4 = buffer.data(shp0 + 4);
    const auto *shp0_8 = buffer.data(shp0 + 8);
    const auto *shp0_9 = buffer.data(shp0 + 9);
    const auto *shp0_10 = buffer.data(shp0 + 10);
    const auto *shp0_11 = buffer.data(shp0 + 11);
    const auto *shp0_15 = buffer.data(shp0 + 15);
    const auto *shp0_16 = buffer.data(shp0 + 16);
    const auto *shp0_17 = buffer.data(shp0 + 17);
    const auto *shp0_18 = buffer.data(shp0 + 18);
    const auto *shp0_19 = buffer.data(shp0 + 19);
    const auto *shp0_20 = buffer.data(shp0 + 20);
    const auto *shp0_23 = buffer.data(shp0 + 23);
    const auto *shp0_25 = buffer.data(shp0 + 25);
    const auto *shp0_27 = buffer.data(shp0 + 27);
    const auto *shp0_28 = buffer.data(shp0 + 28);
    const auto *shp0_29 = buffer.data(shp0 + 29);

    const auto *shp1_0 = buffer.data(shp1 + 0);
    const auto *shp1_1 = buffer.data(shp1 + 1);
    const auto *shp1_2 = buffer.data(shp1 + 2);
    const auto *shp1_4 = buffer.data(shp1 + 4);
    const auto *shp1_8 = buffer.data(shp1 + 8);
    const auto *shp1_9 = buffer.data(shp1 + 9);
    const auto *shp1_10 = buffer.data(shp1 + 10);
    const auto *shp1_11 = buffer.data(shp1 + 11);
    const auto *shp1_15 = buffer.data(shp1 + 15);
    const auto *shp1_16 = buffer.data(shp1 + 16);
    const auto *shp1_17 = buffer.data(shp1 + 17);
    const auto *shp1_18 = buffer.data(shp1 + 18);
    const auto *shp1_19 = buffer.data(shp1 + 19);
    const auto *shp1_20 = buffer.data(shp1 + 20);
    const auto *shp1_23 = buffer.data(shp1 + 23);
    const auto *shp1_25 = buffer.data(shp1 + 25);
    const auto *shp1_27 = buffer.data(shp1 + 27);
    const auto *shp1_28 = buffer.data(shp1 + 28);
    const auto *shp1_29 = buffer.data(shp1 + 29);

    const auto *shd_0 = buffer.data(shd + 0);
    const auto *shd_3 = buffer.data(shd + 3);
    const auto *shd_4 = buffer.data(shd + 4);
    const auto *shd_5 = buffer.data(shd + 5);
    const auto *shd_6 = buffer.data(shd + 6);
    const auto *shd_9 = buffer.data(shd + 9);
    const auto *shd_10 = buffer.data(shd + 10);
    const auto *shd_11 = buffer.data(shd + 11);
    const auto *shd_12 = buffer.data(shd + 12);
    const auto *shd_15 = buffer.data(shd + 15);
    const auto *shd_16 = buffer.data(shd + 16);
    const auto *shd_17 = buffer.data(shd + 17);
    const auto *shd_18 = buffer.data(shd + 18);
    const auto *shd_21 = buffer.data(shd + 21);
    const auto *shd_22 = buffer.data(shd + 22);
    const auto *shd_23 = buffer.data(shd + 23);
    const auto *shd_24 = buffer.data(shd + 24);
    const auto *shd_27 = buffer.data(shd + 27);
    const auto *shd_28 = buffer.data(shd + 28);
    const auto *shd_29 = buffer.data(shd + 29);
    const auto *shd_30 = buffer.data(shd + 30);
    const auto *shd_33 = buffer.data(shd + 33);
    const auto *shd_34 = buffer.data(shd + 34);
    const auto *shd_35 = buffer.data(shd + 35);
    const auto *shd_36 = buffer.data(shd + 36);
    const auto *shd_39 = buffer.data(shd + 39);
    const auto *shd_40 = buffer.data(shd + 40);
    const auto *shd_41 = buffer.data(shd + 41);
    const auto *shd_42 = buffer.data(shd + 42);
    const auto *shd_45 = buffer.data(shd + 45);
    const auto *shd_46 = buffer.data(shd + 46);
    const auto *shd_47 = buffer.data(shd + 47);
    const auto *shd_48 = buffer.data(shd + 48);
    const auto *shd_51 = buffer.data(shd + 51);
    const auto *shd_52 = buffer.data(shd + 52);
    const auto *shd_53 = buffer.data(shd + 53);
    const auto *shd_54 = buffer.data(shd + 54);
    const auto *shd_57 = buffer.data(shd + 57);
    const auto *shd_58 = buffer.data(shd + 58);
    const auto *shd_59 = buffer.data(shd + 59);
    const auto *shd_60 = buffer.data(shd + 60);
    const auto *shd_63 = buffer.data(shd + 63);
    const auto *shd_64 = buffer.data(shd + 64);
    const auto *shd_65 = buffer.data(shd + 65);
    const auto *shd_66 = buffer.data(shd + 66);
    const auto *shd_69 = buffer.data(shd + 69);
    const auto *shd_70 = buffer.data(shd + 70);
    const auto *shd_71 = buffer.data(shd + 71);
    const auto *shd_72 = buffer.data(shd + 72);
    const auto *shd_75 = buffer.data(shd + 75);
    const auto *shd_76 = buffer.data(shd + 76);
    const auto *shd_77 = buffer.data(shd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, sgd_0, sgd_3, sgd_4, \
                         shp0_0, shp1_0, shd_0, shd_3, shd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgd_0[k]
                 + f_1 * shp0_0[k]
                 - f_2 * shp1_0[k]
                 + f_3 * pc_x[k] * shd_0[k];

        t_1[k] = f_3 * pc_y[k] * shd_0[k];

        t_2[k] = f_3 * pc_z[k] * shd_0[k];

        t_3[k] = f_0 * sgd_3[k]
                 + f_3 * pc_x[k] * shd_3[k];

        t_4[k] = f_0 * sgd_4[k]
                 + f_3 * pc_x[k] * shd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sgd_5, shp0_1, shp0_2, \
                         shp1_1, shp1_2, shd_3, shd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * sgd_5[k]
                 + f_3 * pc_x[k] * shd_5[k];

        t_6[k] = f_1 * shp0_1[k]
                 - f_2 * shp1_1[k]
                 + f_3 * pc_y[k] * shd_3[k];

        t_7[k] = f_3 * pc_z[k] * shd_3[k];

        t_8[k] = f_3 * pc_y[k] * shd_5[k];

        t_9[k] = f_1 * shp0_2[k]
                 - f_2 * shp1_2[k]
                 + f_3 * pc_z[k] * shd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, sgf0_0, sgd_0, sgd_9, \
                         sgf1_0, shd_6, shd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * sgf0_0[k]
                  - f_4 * pc_y[k] * sgf1_0[k];

        t_11[k] = f_5 * sgd_0[k]
                  + f_3 * pc_y[k] * shd_6[k];

        t_12[k] = f_3 * pc_z[k] * shd_6[k];

        t_13[k] = f_6 * sgd_9[k]
                  + f_3 * pc_x[k] * shd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, sgd_3, sgd_10, sgd_11, \
                         shp0_4, shp1_4, shd_9, shd_10, shd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_6 * sgd_10[k]
                  + f_3 * pc_x[k] * shd_10[k];

        t_15[k] = f_6 * sgd_11[k]
                  + f_3 * pc_x[k] * shd_11[k];

        t_16[k] = f_5 * sgd_3[k]
                  + f_1 * shp0_4[k]
                  - f_2 * shp1_4[k]
                  + f_3 * pc_y[k] * shd_9[k];

        t_17[k] = f_3 * pc_z[k] * shd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_y, pb_z, pc_y, pc_z, sgf0_0, sgf0_9, \
                         sgd_5, sgf1_0, sgf1_9, shd_11, shd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * sgd_5[k]
                  + f_3 * pc_y[k] * shd_11[k];

        t_19[k] = pb_y[k] * sgf0_9[k]
                  - f_4 * pc_y[k] * sgf1_9[k];

        t_20[k] = pb_z[k] * sgf0_0[k]
                  - f_4 * pc_z[k] * sgf1_0[k];

        t_21[k] = f_3 * pc_y[k] * shd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, sgd_0, sgd_15, sgd_16, sgd_17, \
                         shd_12, shd_15, shd_16, shd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * sgd_0[k]
                  + f_3 * pc_z[k] * shd_12[k];

        t_23[k] = f_6 * sgd_15[k]
                  + f_3 * pc_x[k] * shd_15[k];

        t_24[k] = f_6 * sgd_16[k]
                  + f_3 * pc_x[k] * shd_16[k];

        t_25[k] = f_6 * sgd_17[k]
                  + f_3 * pc_x[k] * shd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_z, pc_y, pc_z, sgf0_6, sgd_3, sgd_5, \
                         sgf1_6, shp0_8, shp1_8, shd_15, shd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_z[k] * sgf0_6[k]
                  - f_4 * pc_z[k] * sgf1_6[k];

        t_27[k] = f_5 * sgd_3[k]
                  + f_3 * pc_z[k] * shd_15[k];

        t_28[k] = f_3 * pc_y[k] * shd_17[k];

        t_29[k] = f_5 * sgd_5[k]
                  + f_1 * shp0_8[k]
                  - f_2 * shp1_8[k]
                  + f_3 * pc_z[k] * shd_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, sgd_6, sgd_18, sgd_21, \
                         shp0_9, shp1_9, shd_18, shd_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sgd_18[k]
                  + f_1 * shp0_9[k]
                  - f_2 * shp1_9[k]
                  + f_3 * pc_x[k] * shd_18[k];

        t_31[k] = f_8 * sgd_6[k]
                  + f_3 * pc_y[k] * shd_18[k];

        t_32[k] = f_3 * pc_z[k] * shd_18[k];

        t_33[k] = f_7 * sgd_21[k]
                  + f_3 * pc_x[k] * shd_21[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pc_x, pc_y, pc_z, sgd_9, sgd_22, sgd_23, \
                         shp0_10, shp1_10, shd_21, shd_22, shd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_7 * sgd_22[k]
                  + f_3 * pc_x[k] * shd_22[k];

        t_35[k] = f_7 * sgd_23[k]
                  + f_3 * pc_x[k] * shd_23[k];

        t_36[k] = f_8 * sgd_9[k]
                  + f_1 * shp0_10[k]
                  - f_2 * shp1_10[k]
                  + f_3 * pc_y[k] * shd_21[k];

        t_37[k] = f_3 * pc_z[k] * shd_21[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, pc_z, sgf0_20, sgd_11, sgd_12, \
                         sgf1_20, shp0_11, shp1_11, shd_23, shd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * sgd_11[k]
                  + f_3 * pc_y[k] * shd_23[k];

        t_39[k] = f_1 * shp0_11[k]
                  - f_2 * shp1_11[k]
                  + f_3 * pc_z[k] * shd_23[k];

        t_40[k] = pb_y[k] * sgf0_20[k]
                  - f_4 * pc_y[k] * sgf1_20[k];

        t_41[k] = f_5 * sgd_12[k]
                  + f_3 * pc_y[k] * shd_24[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pc_x, pc_z, sgd_6, sgd_27, sgd_28, sgd_29, \
                         shd_24, shd_27, shd_28, shd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * sgd_6[k]
                  + f_3 * pc_z[k] * shd_24[k];

        t_43[k] = f_7 * sgd_27[k]
                  + f_3 * pc_x[k] * shd_27[k];

        t_44[k] = f_7 * sgd_28[k]
                  + f_3 * pc_x[k] * shd_28[k];

        t_45[k] = f_7 * sgd_29[k]
                  + f_3 * pc_x[k] * shd_29[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pb_y, pb_z, pc_y, pc_z, sgf0_16, sgf0_29, \
                         sgd_9, sgd_17, sgf1_16, sgf1_29, shd_27, \
                         shd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pb_z[k] * sgf0_16[k]
                  - f_4 * pc_z[k] * sgf1_16[k];

        t_47[k] = f_5 * sgd_9[k]
                  + f_3 * pc_z[k] * shd_27[k];

        t_48[k] = f_5 * sgd_17[k]
                  + f_3 * pc_y[k] * shd_29[k];

        t_49[k] = pb_y[k] * sgf0_29[k]
                  - f_4 * pc_y[k] * sgf1_29[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pc_x, pc_y, pc_z, sgd_12, sgd_30, sgd_33, \
                         shp0_15, shp1_15, shd_30, shd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_7 * sgd_30[k]
                  + f_1 * shp0_15[k]
                  - f_2 * shp1_15[k]
                  + f_3 * pc_x[k] * shd_30[k];

        t_51[k] = f_3 * pc_y[k] * shd_30[k];

        t_52[k] = f_8 * sgd_12[k]
                  + f_3 * pc_z[k] * shd_30[k];

        t_53[k] = f_7 * sgd_33[k]
                  + f_3 * pc_x[k] * shd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sgd_15, sgd_34, \
                         sgd_35, shp0_16, shp1_16, shd_33, shd_34, \
                         shd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * sgd_34[k]
                  + f_3 * pc_x[k] * shd_34[k];

        t_55[k] = f_7 * sgd_35[k]
                  + f_3 * pc_x[k] * shd_35[k];

        t_56[k] = f_1 * shp0_16[k]
                  - f_2 * shp1_16[k]
                  + f_3 * pc_y[k] * shd_33[k];

        t_57[k] = f_8 * sgd_15[k]
                  + f_3 * pc_z[k] * shd_33[k];

        t_58[k] = f_3 * pc_y[k] * shd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_y, pc_z, sgd_17, sgd_18, sgd_36, \
                         shp0_17, shp0_18, shp1_17, shp1_18, shd_35, \
                         shd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_8 * sgd_17[k]
                  + f_1 * shp0_17[k]
                  - f_2 * shp1_17[k]
                  + f_3 * pc_z[k] * shd_35[k];

        t_60[k] = f_8 * sgd_36[k]
                  + f_1 * shp0_18[k]
                  - f_2 * shp1_18[k]
                  + f_3 * pc_x[k] * shd_36[k];

        t_61[k] = f_7 * sgd_18[k]
                  + f_3 * pc_y[k] * shd_36[k];

        t_62[k] = f_3 * pc_z[k] * shd_36[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_y, sgd_21, sgd_39, sgd_40, sgd_41, \
                         shp0_19, shp1_19, shd_39, shd_40, shd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_8 * sgd_39[k]
                  + f_3 * pc_x[k] * shd_39[k];

        t_64[k] = f_8 * sgd_40[k]
                  + f_3 * pc_x[k] * shd_40[k];

        t_65[k] = f_8 * sgd_41[k]
                  + f_3 * pc_x[k] * shd_41[k];

        t_66[k] = f_7 * sgd_21[k]
                  + f_1 * shp0_19[k]
                  - f_2 * shp1_19[k]
                  + f_3 * pc_y[k] * shd_39[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_z, pc_y, pc_z, sgf0_30, sgd_23, sgf1_30, \
                         shp0_20, shp1_20, shd_39, shd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_3 * pc_z[k] * shd_39[k];

        t_68[k] = f_7 * sgd_23[k]
                  + f_3 * pc_y[k] * shd_41[k];

        t_69[k] = f_1 * shp0_20[k]
                  - f_2 * shp1_20[k]
                  + f_3 * pc_z[k] * shd_41[k];

        t_70[k] = pb_z[k] * sgf0_30[k]
                  - f_4 * pc_z[k] * sgf1_30[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_x, pc_y, pc_z, sgd_18, sgd_24, sgd_45, \
                         sgd_46, shd_42, shd_45, shd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_8 * sgd_24[k]
                  + f_3 * pc_y[k] * shd_42[k];

        t_72[k] = f_5 * sgd_18[k]
                  + f_3 * pc_z[k] * shd_42[k];

        t_73[k] = f_8 * sgd_45[k]
                  + f_3 * pc_x[k] * shd_45[k];

        t_74[k] = f_8 * sgd_46[k]
                  + f_3 * pc_x[k] * shd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_z, pc_x, pc_y, pc_z, sgf0_36, sgd_21, \
                         sgd_29, sgd_47, sgf1_36, shd_45, shd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * sgd_47[k]
                  + f_3 * pc_x[k] * shd_47[k];

        t_76[k] = pb_z[k] * sgf0_36[k]
                  - f_4 * pc_z[k] * sgf1_36[k];

        t_77[k] = f_5 * sgd_21[k]
                  + f_3 * pc_z[k] * shd_45[k];

        t_78[k] = f_8 * sgd_29[k]
                  + f_3 * pc_y[k] * shd_47[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_y, pc_y, pc_z, sgf0_50, sgd_23, sgd_24, \
                         sgd_30, sgf1_50, shp0_23, shp1_23, shd_47, \
                         shd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_5 * sgd_23[k]
                  + f_1 * shp0_23[k]
                  - f_2 * shp1_23[k]
                  + f_3 * pc_z[k] * shd_47[k];

        t_80[k] = pb_y[k] * sgf0_50[k]
                  - f_4 * pc_y[k] * sgf1_50[k];

        t_81[k] = f_5 * sgd_30[k]
                  + f_3 * pc_y[k] * shd_48[k];

        t_82[k] = f_8 * sgd_24[k]
                  + f_3 * pc_z[k] * shd_48[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, sgd_33, sgd_51, sgd_52, sgd_53, \
                         shp0_25, shp1_25, shd_51, shd_52, shd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_8 * sgd_51[k]
                  + f_3 * pc_x[k] * shd_51[k];

        t_84[k] = f_8 * sgd_52[k]
                  + f_3 * pc_x[k] * shd_52[k];

        t_85[k] = f_8 * sgd_53[k]
                  + f_3 * pc_x[k] * shd_53[k];

        t_86[k] = f_5 * sgd_33[k]
                  + f_1 * shp0_25[k]
                  - f_2 * shp1_25[k]
                  + f_3 * pc_y[k] * shd_51[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_y, pc_y, pc_z, sgf0_59, sgd_27, sgd_35, sgf1_59, \
                         shd_51, shd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_8 * sgd_27[k]
                  + f_3 * pc_z[k] * shd_51[k];

        t_88[k] = f_5 * sgd_35[k]
                  + f_3 * pc_y[k] * shd_53[k];

        t_89[k] = pb_y[k] * sgf0_59[k]
                  - f_4 * pc_y[k] * sgf1_59[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pc_x, pc_y, pc_z, sgd_30, sgd_54, sgd_57, \
                         shp0_27, shp1_27, shd_54, shd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_8 * sgd_54[k]
                  + f_1 * shp0_27[k]
                  - f_2 * shp1_27[k]
                  + f_3 * pc_x[k] * shd_54[k];

        t_91[k] = f_3 * pc_y[k] * shd_54[k];

        t_92[k] = f_7 * sgd_30[k]
                  + f_3 * pc_z[k] * shd_54[k];

        t_93[k] = f_8 * sgd_57[k]
                  + f_3 * pc_x[k] * shd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, sgd_33, sgd_58, \
                         sgd_59, shp0_28, shp1_28, shd_57, shd_58, \
                         shd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_8 * sgd_58[k]
                  + f_3 * pc_x[k] * shd_58[k];

        t_95[k] = f_8 * sgd_59[k]
                  + f_3 * pc_x[k] * shd_59[k];

        t_96[k] = f_1 * shp0_28[k]
                  - f_2 * shp1_28[k]
                  + f_3 * pc_y[k] * shd_57[k];

        t_97[k] = f_7 * sgd_33[k]
                  + f_3 * pc_z[k] * shd_57[k];

        t_98[k] = f_3 * pc_y[k] * shd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pb_x, pc_x, pc_y, pc_z, sgf0_100, sgd_35, sgd_36, \
                         sgd_60, sgf1_100, shp0_29, shp1_29, shd_59, \
                         shd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_7 * sgd_35[k]
                  + f_1 * shp0_29[k]
                  - f_2 * shp1_29[k]
                  + f_3 * pc_z[k] * shd_59[k];

        t_100[k] = pb_x[k] * sgf0_100[k]
                   + f_7 * sgd_60[k]
                   - f_4 * pc_x[k] * sgf1_100[k];

        t_101[k] = f_6 * sgd_36[k]
                   + f_3 * pc_y[k] * shd_60[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_z, sgd_63, sgd_64, sgd_65, \
                         shd_60, shd_63, shd_64, shd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_3 * pc_z[k] * shd_60[k];

        t_103[k] = f_5 * sgd_63[k]
                   + f_3 * pc_x[k] * shd_63[k];

        t_104[k] = f_5 * sgd_64[k]
                   + f_3 * pc_x[k] * shd_64[k];

        t_105[k] = f_5 * sgd_65[k]
                   + f_3 * pc_x[k] * shd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_x, pc_x, pc_y, pc_z, sgf0_106, \
                         sgf0_109, sgd_41, sgf1_106, sgf1_109, shd_63, \
                         shd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_x[k] * sgf0_106[k]
                   - f_4 * pc_x[k] * sgf1_106[k];

        t_107[k] = f_3 * pc_z[k] * shd_63[k];

        t_108[k] = f_6 * sgd_41[k]
                   + f_3 * pc_y[k] * shd_65[k];

        t_109[k] = pb_x[k] * sgf0_109[k]
                   - f_4 * pc_x[k] * sgf1_109[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_z, pc_x, pc_y, pc_z, sgf0_60, sgd_36, \
                         sgd_42, sgd_69, sgf1_60, shd_66, shd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_z[k] * sgf0_60[k]
                   - f_4 * pc_z[k] * sgf1_60[k];

        t_111[k] = f_7 * sgd_42[k]
                   + f_3 * pc_y[k] * shd_66[k];

        t_112[k] = f_5 * sgd_36[k]
                   + f_3 * pc_z[k] * shd_66[k];

        t_113[k] = f_5 * sgd_69[k]
                   + f_3 * pc_x[k] * shd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pc_x, pc_z, sgf0_116, sgd_39, \
                         sgd_70, sgd_71, sgf1_116, shd_69, shd_70, \
                         shd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_5 * sgd_70[k]
                   + f_3 * pc_x[k] * shd_70[k];

        t_115[k] = f_5 * sgd_71[k]
                   + f_3 * pc_x[k] * shd_71[k];

        t_116[k] = pb_x[k] * sgf0_116[k]
                   - f_4 * pc_x[k] * sgf1_116[k];

        t_117[k] = f_5 * sgd_39[k]
                   + f_3 * pc_z[k] * shd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pb_x, pc_x, pc_y, sgf0_119, sgf0_120, \
                         sgd_47, sgd_48, sgd_72, sgf1_119, sgf1_120, shd_71, \
                         shd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_7 * sgd_47[k]
                   + f_3 * pc_y[k] * shd_71[k];

        t_119[k] = pb_x[k] * sgf0_119[k]
                   - f_4 * pc_x[k] * sgf1_119[k];

        t_120[k] = pb_x[k] * sgf0_120[k]
                   + f_7 * sgd_72[k]
                   - f_4 * pc_x[k] * sgf1_120[k];

        t_121[k] = f_8 * sgd_48[k]
                   + f_3 * pc_y[k] * shd_72[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pc_x, pc_z, sgd_42, sgd_75, sgd_76, \
                         sgd_77, shd_72, shd_75, shd_76, shd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * sgd_42[k]
                   + f_3 * pc_z[k] * shd_72[k];

        t_123[k] = f_5 * sgd_75[k]
                   + f_3 * pc_x[k] * shd_75[k];

        t_124[k] = f_5 * sgd_76[k]
                   + f_3 * pc_x[k] * shd_76[k];

        t_125[k] = f_5 * sgd_77[k]
                   + f_3 * pc_x[k] * shd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pb_x, pc_x, pc_y, pc_z, sgf0_126, \
                         sgf0_129, sgd_45, sgd_53, sgf1_126, sgf1_129, shd_75, \
                         shd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pb_x[k] * sgf0_126[k]
                   - f_4 * pc_x[k] * sgf1_126[k];

        t_127[k] = f_8 * sgd_45[k]
                   + f_3 * pc_z[k] * shd_75[k];

        t_128[k] = f_8 * sgd_53[k]
                   + f_3 * pc_y[k] * shd_77[k];

        t_129[k] = pb_x[k] * sgf0_129[k]
                   - f_4 * pc_x[k] * sgf1_129[k];
    }
}

static auto
compute_prim_shf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sgf0,
                                                          const size_t sgd, const size_t sgf1,
                                                          const size_t shp0, const size_t shp1,
                                                          const size_t shd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 2.0 / q;
    const auto f_7 = 1.5 / q;
    const auto f_8 = 1.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgf0_90 = buffer.data(sgf0 + 90);
    const auto *sgf0_100 = buffer.data(sgf0 + 100);
    const auto *sgf0_106 = buffer.data(sgf0 + 106);
    const auto *sgf0_136 = buffer.data(sgf0 + 136);
    const auto *sgf0_139 = buffer.data(sgf0 + 139);
    const auto *sgf0_140 = buffer.data(sgf0 + 140);
    const auto *sgf0_146 = buffer.data(sgf0 + 146);
    const auto *sgf0_149 = buffer.data(sgf0 + 149);

    const auto *sgd_48 = buffer.data(sgd + 48);
    const auto *sgd_51 = buffer.data(sgd + 51);
    const auto *sgd_54 = buffer.data(sgd + 54);
    const auto *sgd_57 = buffer.data(sgd + 57);
    const auto *sgd_59 = buffer.data(sgd + 59);
    const auto *sgd_60 = buffer.data(sgd + 60);
    const auto *sgd_63 = buffer.data(sgd + 63);
    const auto *sgd_65 = buffer.data(sgd + 65);
    const auto *sgd_66 = buffer.data(sgd + 66);
    const auto *sgd_69 = buffer.data(sgd + 69);
    const auto *sgd_71 = buffer.data(sgd + 71);
    const auto *sgd_72 = buffer.data(sgd + 72);
    const auto *sgd_75 = buffer.data(sgd + 75);
    const auto *sgd_77 = buffer.data(sgd + 77);
    const auto *sgd_78 = buffer.data(sgd + 78);
    const auto *sgd_81 = buffer.data(sgd + 81);
    const auto *sgd_82 = buffer.data(sgd + 82);
    const auto *sgd_83 = buffer.data(sgd + 83);
    const auto *sgd_84 = buffer.data(sgd + 84);
    const auto *sgd_87 = buffer.data(sgd + 87);
    const auto *sgd_88 = buffer.data(sgd + 88);
    const auto *sgd_89 = buffer.data(sgd + 89);

    const auto *sgf1_90 = buffer.data(sgf1 + 90);
    const auto *sgf1_100 = buffer.data(sgf1 + 100);
    const auto *sgf1_106 = buffer.data(sgf1 + 106);
    const auto *sgf1_136 = buffer.data(sgf1 + 136);
    const auto *sgf1_139 = buffer.data(sgf1 + 139);
    const auto *sgf1_140 = buffer.data(sgf1 + 140);
    const auto *sgf1_146 = buffer.data(sgf1 + 146);
    const auto *sgf1_149 = buffer.data(sgf1 + 149);

    const auto *shp0_45 = buffer.data(shp0 + 45);
    const auto *shp0_46 = buffer.data(shp0 + 46);
    const auto *shp0_47 = buffer.data(shp0 + 47);
    const auto *shp0_50 = buffer.data(shp0 + 50);
    const auto *shp0_51 = buffer.data(shp0 + 51);
    const auto *shp0_52 = buffer.data(shp0 + 52);
    const auto *shp0_53 = buffer.data(shp0 + 53);
    const auto *shp0_54 = buffer.data(shp0 + 54);
    const auto *shp0_55 = buffer.data(shp0 + 55);
    const auto *shp0_56 = buffer.data(shp0 + 56);
    const auto *shp0_60 = buffer.data(shp0 + 60);
    const auto *shp0_61 = buffer.data(shp0 + 61);
    const auto *shp0_62 = buffer.data(shp0 + 62);

    const auto *shp1_45 = buffer.data(shp1 + 45);
    const auto *shp1_46 = buffer.data(shp1 + 46);
    const auto *shp1_47 = buffer.data(shp1 + 47);
    const auto *shp1_50 = buffer.data(shp1 + 50);
    const auto *shp1_51 = buffer.data(shp1 + 51);
    const auto *shp1_52 = buffer.data(shp1 + 52);
    const auto *shp1_53 = buffer.data(shp1 + 53);
    const auto *shp1_54 = buffer.data(shp1 + 54);
    const auto *shp1_55 = buffer.data(shp1 + 55);
    const auto *shp1_56 = buffer.data(shp1 + 56);
    const auto *shp1_60 = buffer.data(shp1 + 60);
    const auto *shp1_61 = buffer.data(shp1 + 61);
    const auto *shp1_62 = buffer.data(shp1 + 62);

    const auto *shd_78 = buffer.data(shd + 78);
    const auto *shd_81 = buffer.data(shd + 81);
    const auto *shd_82 = buffer.data(shd + 82);
    const auto *shd_83 = buffer.data(shd + 83);
    const auto *shd_84 = buffer.data(shd + 84);
    const auto *shd_87 = buffer.data(shd + 87);
    const auto *shd_88 = buffer.data(shd + 88);
    const auto *shd_89 = buffer.data(shd + 89);
    const auto *shd_90 = buffer.data(shd + 90);
    const auto *shd_93 = buffer.data(shd + 93);
    const auto *shd_94 = buffer.data(shd + 94);
    const auto *shd_95 = buffer.data(shd + 95);
    const auto *shd_96 = buffer.data(shd + 96);
    const auto *shd_99 = buffer.data(shd + 99);
    const auto *shd_100 = buffer.data(shd + 100);
    const auto *shd_101 = buffer.data(shd + 101);
    const auto *shd_102 = buffer.data(shd + 102);
    const auto *shd_105 = buffer.data(shd + 105);
    const auto *shd_106 = buffer.data(shd + 106);
    const auto *shd_107 = buffer.data(shd + 107);
    const auto *shd_108 = buffer.data(shd + 108);
    const auto *shd_111 = buffer.data(shd + 111);
    const auto *shd_112 = buffer.data(shd + 112);
    const auto *shd_113 = buffer.data(shd + 113);
    const auto *shd_114 = buffer.data(shd + 114);
    const auto *shd_117 = buffer.data(shd + 117);
    const auto *shd_118 = buffer.data(shd + 118);
    const auto *shd_119 = buffer.data(shd + 119);
    const auto *shd_120 = buffer.data(shd + 120);
    const auto *shd_123 = buffer.data(shd + 123);
    const auto *shd_124 = buffer.data(shd + 124);
    const auto *shd_125 = buffer.data(shd + 125);

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pb_y, pc_x, pc_y, pc_z, sgf0_90, sgd_48, \
                         sgd_54, sgd_81, sgf1_90, shd_78, shd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pb_y[k] * sgf0_90[k]
                   - f_4 * pc_y[k] * sgf1_90[k];

        t_131[k] = f_5 * sgd_54[k]
                   + f_3 * pc_y[k] * shd_78[k];

        t_132[k] = f_7 * sgd_48[k]
                   + f_3 * pc_z[k] * shd_78[k];

        t_133[k] = f_5 * sgd_81[k]
                   + f_3 * pc_x[k] * shd_81[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pb_x, pc_x, pc_z, sgf0_136, sgd_51, \
                         sgd_82, sgd_83, sgf1_136, shd_81, shd_82, \
                         shd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_5 * sgd_82[k]
                   + f_3 * pc_x[k] * shd_82[k];

        t_135[k] = f_5 * sgd_83[k]
                   + f_3 * pc_x[k] * shd_83[k];

        t_136[k] = pb_x[k] * sgf0_136[k]
                   - f_4 * pc_x[k] * sgf1_136[k];

        t_137[k] = f_7 * sgd_51[k]
                   + f_3 * pc_z[k] * shd_81[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pb_x, pc_x, pc_y, sgf0_139, sgf0_140, \
                         sgd_59, sgd_84, sgf1_139, sgf1_140, shd_83, \
                         shd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_5 * sgd_59[k]
                   + f_3 * pc_y[k] * shd_83[k];

        t_139[k] = pb_x[k] * sgf0_139[k]
                   - f_4 * pc_x[k] * sgf1_139[k];

        t_140[k] = pb_x[k] * sgf0_140[k]
                   + f_7 * sgd_84[k]
                   - f_4 * pc_x[k] * sgf1_140[k];

        t_141[k] = f_3 * pc_y[k] * shd_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_x, pc_z, sgd_54, sgd_87, sgd_88, \
                         sgd_89, shd_84, shd_87, shd_88, shd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_6 * sgd_54[k]
                   + f_3 * pc_z[k] * shd_84[k];

        t_143[k] = f_5 * sgd_87[k]
                   + f_3 * pc_x[k] * shd_87[k];

        t_144[k] = f_5 * sgd_88[k]
                   + f_3 * pc_x[k] * shd_88[k];

        t_145[k] = f_5 * sgd_89[k]
                   + f_3 * pc_x[k] * shd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pb_x, pc_x, pc_y, pc_z, sgf0_146, \
                         sgf0_149, sgd_57, sgf1_146, sgf1_149, shd_87, \
                         shd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pb_x[k] * sgf0_146[k]
                   - f_4 * pc_x[k] * sgf1_146[k];

        t_147[k] = f_6 * sgd_57[k]
                   + f_3 * pc_z[k] * shd_87[k];

        t_148[k] = f_3 * pc_y[k] * shd_89[k];

        t_149[k] = pb_x[k] * sgf0_149[k]
                   - f_4 * pc_x[k] * sgf1_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, pc_z, sgd_60, \
                         shp0_45, shp1_45, shd_90, shd_93, shd_94, \
                         shd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_1 * shp0_45[k]
                   - f_2 * shp1_45[k]
                   + f_3 * pc_x[k] * shd_90[k];

        t_151[k] = f_0 * sgd_60[k]
                   + f_3 * pc_y[k] * shd_90[k];

        t_152[k] = f_3 * pc_z[k] * shd_90[k];

        t_153[k] = f_3 * pc_x[k] * shd_93[k];

        t_154[k] = f_3 * pc_x[k] * shd_94[k];

        t_155[k] = f_3 * pc_x[k] * shd_95[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_y, pc_z, sgd_63, sgd_65, shp0_46, \
                         shp0_47, shp1_46, shp1_47, shd_93, shd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_0 * sgd_63[k]
                   + f_1 * shp0_46[k]
                   - f_2 * shp1_46[k]
                   + f_3 * pc_y[k] * shd_93[k];

        t_157[k] = f_3 * pc_z[k] * shd_93[k];

        t_158[k] = f_0 * sgd_65[k]
                   + f_3 * pc_y[k] * shd_95[k];

        t_159[k] = f_1 * shp0_47[k]
                   - f_2 * shp1_47[k]
                   + f_3 * pc_z[k] * shd_95[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, t_164, pb_z, pc_x, pc_y, pc_z, sgf0_100, \
                         sgd_60, sgd_66, sgf1_100, shd_96, shd_99, \
                         shd_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pb_z[k] * sgf0_100[k]
                   - f_4 * pc_z[k] * sgf1_100[k];

        t_161[k] = f_6 * sgd_66[k]
                   + f_3 * pc_y[k] * shd_96[k];

        t_162[k] = f_5 * sgd_60[k]
                   + f_3 * pc_z[k] * shd_96[k];

        t_163[k] = f_3 * pc_x[k] * shd_99[k];

        t_164[k] = f_3 * pc_x[k] * shd_100[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_z, pc_x, pc_y, pc_z, sgf0_106, sgd_63, \
                         sgd_71, sgf1_106, shd_99, shd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_3 * pc_x[k] * shd_101[k];

        t_166[k] = pb_z[k] * sgf0_106[k]
                   - f_4 * pc_z[k] * sgf1_106[k];

        t_167[k] = f_5 * sgd_63[k]
                   + f_3 * pc_z[k] * shd_99[k];

        t_168[k] = f_6 * sgd_71[k]
                   + f_3 * pc_y[k] * shd_101[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, sgd_65, sgd_66, sgd_72, \
                         shp0_50, shp0_51, shp1_50, shp1_51, shd_101, \
                         shd_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_5 * sgd_65[k]
                   + f_1 * shp0_50[k]
                   - f_2 * shp1_50[k]
                   + f_3 * pc_z[k] * shd_101[k];

        t_170[k] = f_1 * shp0_51[k]
                   - f_2 * shp1_51[k]
                   + f_3 * pc_x[k] * shd_102[k];

        t_171[k] = f_7 * sgd_72[k]
                   + f_3 * pc_y[k] * shd_102[k];

        t_172[k] = f_8 * sgd_66[k]
                   + f_3 * pc_z[k] * shd_102[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, sgd_69, sgd_75, \
                         shp0_52, shp1_52, shd_105, shd_106, shd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_3 * pc_x[k] * shd_105[k];

        t_174[k] = f_3 * pc_x[k] * shd_106[k];

        t_175[k] = f_3 * pc_x[k] * shd_107[k];

        t_176[k] = f_7 * sgd_75[k]
                   + f_1 * shp0_52[k]
                   - f_2 * shp1_52[k]
                   + f_3 * pc_y[k] * shd_105[k];

        t_177[k] = f_8 * sgd_69[k]
                   + f_3 * pc_z[k] * shd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, t_181, pc_x, pc_y, pc_z, sgd_71, sgd_77, sgd_78, \
                         shp0_53, shp0_54, shp1_53, shp1_54, shd_107, \
                         shd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_7 * sgd_77[k]
                   + f_3 * pc_y[k] * shd_107[k];

        t_179[k] = f_8 * sgd_71[k]
                   + f_1 * shp0_53[k]
                   - f_2 * shp1_53[k]
                   + f_3 * pc_z[k] * shd_107[k];

        t_180[k] = f_1 * shp0_54[k]
                   - f_2 * shp1_54[k]
                   + f_3 * pc_x[k] * shd_108[k];

        t_181[k] = f_8 * sgd_78[k]
                   + f_3 * pc_y[k] * shd_108[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, t_186, pc_x, pc_y, pc_z, sgd_72, sgd_81, \
                         shp0_55, shp1_55, shd_108, shd_111, shd_112, \
                         shd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_7 * sgd_72[k]
                   + f_3 * pc_z[k] * shd_108[k];

        t_183[k] = f_3 * pc_x[k] * shd_111[k];

        t_184[k] = f_3 * pc_x[k] * shd_112[k];

        t_185[k] = f_3 * pc_x[k] * shd_113[k];

        t_186[k] = f_8 * sgd_81[k]
                   + f_1 * shp0_55[k]
                   - f_2 * shp1_55[k]
                   + f_3 * pc_y[k] * shd_111[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, pb_y, pc_y, pc_z, sgf0_140, sgd_75, \
                         sgd_77, sgd_83, sgf1_140, shp0_56, shp1_56, shd_111, \
                         shd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_7 * sgd_75[k]
                   + f_3 * pc_z[k] * shd_111[k];

        t_188[k] = f_8 * sgd_83[k]
                   + f_3 * pc_y[k] * shd_113[k];

        t_189[k] = f_7 * sgd_77[k]
                   + f_1 * shp0_56[k]
                   - f_2 * shp1_56[k]
                   + f_3 * pc_z[k] * shd_113[k];

        t_190[k] = pb_y[k] * sgf0_140[k]
                   - f_4 * pc_y[k] * sgf1_140[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, t_195, pc_x, pc_y, pc_z, sgd_78, sgd_84, \
                         shd_114, shd_117, shd_118, shd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_5 * sgd_84[k]
                   + f_3 * pc_y[k] * shd_114[k];

        t_192[k] = f_6 * sgd_78[k]
                   + f_3 * pc_z[k] * shd_114[k];

        t_193[k] = f_3 * pc_x[k] * shd_117[k];

        t_194[k] = f_3 * pc_x[k] * shd_118[k];

        t_195[k] = f_3 * pc_x[k] * shd_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pb_y, pc_y, pc_z, sgf0_146, sgf0_149, \
                         sgd_81, sgd_87, sgd_89, sgf1_146, sgf1_149, shd_117, \
                         shd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pb_y[k] * sgf0_146[k]
                   + f_7 * sgd_87[k]
                   - f_4 * pc_y[k] * sgf1_146[k];

        t_197[k] = f_6 * sgd_81[k]
                   + f_3 * pc_z[k] * shd_117[k];

        t_198[k] = f_5 * sgd_89[k]
                   + f_3 * pc_y[k] * shd_119[k];

        t_199[k] = pb_y[k] * sgf0_149[k]
                   - f_4 * pc_y[k] * sgf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, t_205, pc_x, pc_y, pc_z, sgd_84, \
                         shp0_60, shp1_60, shd_120, shd_123, shd_124, \
                         shd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * shp0_60[k]
                   - f_2 * shp1_60[k]
                   + f_3 * pc_x[k] * shd_120[k];

        t_201[k] = f_3 * pc_y[k] * shd_120[k];

        t_202[k] = f_0 * sgd_84[k]
                   + f_3 * pc_z[k] * shd_120[k];

        t_203[k] = f_3 * pc_x[k] * shd_123[k];

        t_204[k] = f_3 * pc_x[k] * shd_124[k];

        t_205[k] = f_3 * pc_x[k] * shd_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pc_y, pc_z, sgd_87, sgd_89, shp0_61, \
                         shp0_62, shp1_61, shp1_62, shd_123, shd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_1 * shp0_61[k]
                   - f_2 * shp1_61[k]
                   + f_3 * pc_y[k] * shd_123[k];

        t_207[k] = f_0 * sgd_87[k]
                   + f_3 * pc_z[k] * shd_123[k];

        t_208[k] = f_3 * pc_y[k] * shd_125[k];

        t_209[k] = f_0 * sgd_89[k]
                   + f_1 * shp0_62[k]
                   - f_2 * shp1_62[k]
                   + f_3 * pc_z[k] * shd_125[k];
    }
}

auto
compute_prim_shf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgf0, const size_t sgd,
                                                   const size_t sgf1, const size_t shp0,
                                                   const size_t shp1, const size_t shd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_shf_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sgf0, sgd,
                                                              sgf1, shp0, shp1, shd, ncols,
                                                              gamma, p, q);

    compute_prim_shf_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sgf0, sgd,
                                                              sgf1, shp0, shp1, shd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
