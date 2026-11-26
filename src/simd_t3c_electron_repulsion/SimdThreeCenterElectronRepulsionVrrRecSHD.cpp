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


#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_shd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgd0, const size_t sgp,
                                                   const size_t sgd1, const size_t shs0,
                                                   const size_t shs1, const size_t shp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 2.0 / q;
    const auto f_6 = 0.5 / q;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgd0_0 = buffer.data(sgd0 + 0);
    const auto *sgd0_3 = buffer.data(sgd0 + 3);
    const auto *sgd0_5 = buffer.data(sgd0 + 5);
    const auto *sgd0_9 = buffer.data(sgd0 + 9);
    const auto *sgd0_12 = buffer.data(sgd0 + 12);
    const auto *sgd0_17 = buffer.data(sgd0 + 17);
    const auto *sgd0_18 = buffer.data(sgd0 + 18);
    const auto *sgd0_21 = buffer.data(sgd0 + 21);
    const auto *sgd0_30 = buffer.data(sgd0 + 30);
    const auto *sgd0_35 = buffer.data(sgd0 + 35);
    const auto *sgd0_36 = buffer.data(sgd0 + 36);
    const auto *sgd0_54 = buffer.data(sgd0 + 54);
    const auto *sgd0_60 = buffer.data(sgd0 + 60);
    const auto *sgd0_63 = buffer.data(sgd0 + 63);
    const auto *sgd0_65 = buffer.data(sgd0 + 65);
    const auto *sgd0_69 = buffer.data(sgd0 + 69);
    const auto *sgd0_71 = buffer.data(sgd0 + 71);
    const auto *sgd0_72 = buffer.data(sgd0 + 72);
    const auto *sgd0_75 = buffer.data(sgd0 + 75);
    const auto *sgd0_77 = buffer.data(sgd0 + 77);
    const auto *sgd0_81 = buffer.data(sgd0 + 81);
    const auto *sgd0_83 = buffer.data(sgd0 + 83);
    const auto *sgd0_84 = buffer.data(sgd0 + 84);
    const auto *sgd0_87 = buffer.data(sgd0 + 87);
    const auto *sgd0_89 = buffer.data(sgd0 + 89);

    const auto *sgp_0 = buffer.data(sgp + 0);
    const auto *sgp_1 = buffer.data(sgp + 1);
    const auto *sgp_2 = buffer.data(sgp + 2);
    const auto *sgp_4 = buffer.data(sgp + 4);
    const auto *sgp_5 = buffer.data(sgp + 5);
    const auto *sgp_7 = buffer.data(sgp + 7);
    const auto *sgp_8 = buffer.data(sgp + 8);
    const auto *sgp_9 = buffer.data(sgp + 9);
    const auto *sgp_10 = buffer.data(sgp + 10);
    const auto *sgp_11 = buffer.data(sgp + 11);
    const auto *sgp_13 = buffer.data(sgp + 13);
    const auto *sgp_14 = buffer.data(sgp + 14);
    const auto *sgp_15 = buffer.data(sgp + 15);
    const auto *sgp_16 = buffer.data(sgp + 16);
    const auto *sgp_17 = buffer.data(sgp + 17);
    const auto *sgp_18 = buffer.data(sgp + 18);
    const auto *sgp_19 = buffer.data(sgp + 19);
    const auto *sgp_20 = buffer.data(sgp + 20);
    const auto *sgp_22 = buffer.data(sgp + 22);
    const auto *sgp_23 = buffer.data(sgp + 23);
    const auto *sgp_25 = buffer.data(sgp + 25);
    const auto *sgp_26 = buffer.data(sgp + 26);
    const auto *sgp_27 = buffer.data(sgp + 27);
    const auto *sgp_28 = buffer.data(sgp + 28);
    const auto *sgp_29 = buffer.data(sgp + 29);
    const auto *sgp_30 = buffer.data(sgp + 30);
    const auto *sgp_31 = buffer.data(sgp + 31);
    const auto *sgp_32 = buffer.data(sgp + 32);
    const auto *sgp_34 = buffer.data(sgp + 34);
    const auto *sgp_35 = buffer.data(sgp + 35);
    const auto *sgp_36 = buffer.data(sgp + 36);
    const auto *sgp_37 = buffer.data(sgp + 37);
    const auto *sgp_38 = buffer.data(sgp + 38);
    const auto *sgp_40 = buffer.data(sgp + 40);
    const auto *sgp_41 = buffer.data(sgp + 41);
    const auto *sgp_42 = buffer.data(sgp + 42);
    const auto *sgp_43 = buffer.data(sgp + 43);
    const auto *sgp_44 = buffer.data(sgp + 44);

    const auto *sgd1_0 = buffer.data(sgd1 + 0);
    const auto *sgd1_3 = buffer.data(sgd1 + 3);
    const auto *sgd1_5 = buffer.data(sgd1 + 5);
    const auto *sgd1_9 = buffer.data(sgd1 + 9);
    const auto *sgd1_12 = buffer.data(sgd1 + 12);
    const auto *sgd1_17 = buffer.data(sgd1 + 17);
    const auto *sgd1_18 = buffer.data(sgd1 + 18);
    const auto *sgd1_21 = buffer.data(sgd1 + 21);
    const auto *sgd1_30 = buffer.data(sgd1 + 30);
    const auto *sgd1_35 = buffer.data(sgd1 + 35);
    const auto *sgd1_36 = buffer.data(sgd1 + 36);
    const auto *sgd1_54 = buffer.data(sgd1 + 54);
    const auto *sgd1_60 = buffer.data(sgd1 + 60);
    const auto *sgd1_63 = buffer.data(sgd1 + 63);
    const auto *sgd1_65 = buffer.data(sgd1 + 65);
    const auto *sgd1_69 = buffer.data(sgd1 + 69);
    const auto *sgd1_71 = buffer.data(sgd1 + 71);
    const auto *sgd1_72 = buffer.data(sgd1 + 72);
    const auto *sgd1_75 = buffer.data(sgd1 + 75);
    const auto *sgd1_77 = buffer.data(sgd1 + 77);
    const auto *sgd1_81 = buffer.data(sgd1 + 81);
    const auto *sgd1_83 = buffer.data(sgd1 + 83);
    const auto *sgd1_84 = buffer.data(sgd1 + 84);
    const auto *sgd1_87 = buffer.data(sgd1 + 87);
    const auto *sgd1_89 = buffer.data(sgd1 + 89);

    const auto *shs0_0 = buffer.data(shs0 + 0);
    const auto *shs0_1 = buffer.data(shs0 + 1);
    const auto *shs0_2 = buffer.data(shs0 + 2);
    const auto *shs0_3 = buffer.data(shs0 + 3);
    const auto *shs0_5 = buffer.data(shs0 + 5);
    const auto *shs0_6 = buffer.data(shs0 + 6);
    const auto *shs0_7 = buffer.data(shs0 + 7);
    const auto *shs0_8 = buffer.data(shs0 + 8);
    const auto *shs0_9 = buffer.data(shs0 + 9);
    const auto *shs0_15 = buffer.data(shs0 + 15);
    const auto *shs0_16 = buffer.data(shs0 + 16);
    const auto *shs0_17 = buffer.data(shs0 + 17);
    const auto *shs0_18 = buffer.data(shs0 + 18);
    const auto *shs0_20 = buffer.data(shs0 + 20);

    const auto *shs1_0 = buffer.data(shs1 + 0);
    const auto *shs1_1 = buffer.data(shs1 + 1);
    const auto *shs1_2 = buffer.data(shs1 + 2);
    const auto *shs1_3 = buffer.data(shs1 + 3);
    const auto *shs1_5 = buffer.data(shs1 + 5);
    const auto *shs1_6 = buffer.data(shs1 + 6);
    const auto *shs1_7 = buffer.data(shs1 + 7);
    const auto *shs1_8 = buffer.data(shs1 + 8);
    const auto *shs1_9 = buffer.data(shs1 + 9);
    const auto *shs1_15 = buffer.data(shs1 + 15);
    const auto *shs1_16 = buffer.data(shs1 + 16);
    const auto *shs1_17 = buffer.data(shs1 + 17);
    const auto *shs1_18 = buffer.data(shs1 + 18);
    const auto *shs1_20 = buffer.data(shs1 + 20);

    const auto *shp_0 = buffer.data(shp + 0);
    const auto *shp_1 = buffer.data(shp + 1);
    const auto *shp_2 = buffer.data(shp + 2);
    const auto *shp_4 = buffer.data(shp + 4);
    const auto *shp_5 = buffer.data(shp + 5);
    const auto *shp_7 = buffer.data(shp + 7);
    const auto *shp_8 = buffer.data(shp + 8);
    const auto *shp_9 = buffer.data(shp + 9);
    const auto *shp_10 = buffer.data(shp + 10);
    const auto *shp_11 = buffer.data(shp + 11);
    const auto *shp_13 = buffer.data(shp + 13);
    const auto *shp_14 = buffer.data(shp + 14);
    const auto *shp_15 = buffer.data(shp + 15);
    const auto *shp_16 = buffer.data(shp + 16);
    const auto *shp_17 = buffer.data(shp + 17);
    const auto *shp_18 = buffer.data(shp + 18);
    const auto *shp_19 = buffer.data(shp + 19);
    const auto *shp_20 = buffer.data(shp + 20);
    const auto *shp_22 = buffer.data(shp + 22);
    const auto *shp_23 = buffer.data(shp + 23);
    const auto *shp_25 = buffer.data(shp + 25);
    const auto *shp_26 = buffer.data(shp + 26);
    const auto *shp_27 = buffer.data(shp + 27);
    const auto *shp_28 = buffer.data(shp + 28);
    const auto *shp_29 = buffer.data(shp + 29);
    const auto *shp_31 = buffer.data(shp + 31);
    const auto *shp_32 = buffer.data(shp + 32);
    const auto *shp_34 = buffer.data(shp + 34);
    const auto *shp_35 = buffer.data(shp + 35);
    const auto *shp_37 = buffer.data(shp + 37);
    const auto *shp_38 = buffer.data(shp + 38);
    const auto *shp_40 = buffer.data(shp + 40);
    const auto *shp_41 = buffer.data(shp + 41);
    const auto *shp_43 = buffer.data(shp + 43);
    const auto *shp_44 = buffer.data(shp + 44);
    const auto *shp_45 = buffer.data(shp + 45);
    const auto *shp_46 = buffer.data(shp + 46);
    const auto *shp_47 = buffer.data(shp + 47);
    const auto *shp_49 = buffer.data(shp + 49);
    const auto *shp_50 = buffer.data(shp + 50);
    const auto *shp_51 = buffer.data(shp + 51);
    const auto *shp_52 = buffer.data(shp + 52);
    const auto *shp_53 = buffer.data(shp + 53);
    const auto *shp_54 = buffer.data(shp + 54);
    const auto *shp_55 = buffer.data(shp + 55);
    const auto *shp_56 = buffer.data(shp + 56);
    const auto *shp_58 = buffer.data(shp + 58);
    const auto *shp_59 = buffer.data(shp + 59);
    const auto *shp_60 = buffer.data(shp + 60);
    const auto *shp_61 = buffer.data(shp + 61);
    const auto *shp_62 = buffer.data(shp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, sgp_0, sgp_1, sgp_2, shs0_0, \
                         shs1_0, shp_0, shp_1, shp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgp_0[k]
                 + f_1 * shs0_0[k]
                 - f_2 * shs1_0[k]
                 + f_3 * pc_x[k] * shp_0[k];

        t_1[k] = f_0 * sgp_1[k]
                 + f_3 * pc_x[k] * shp_1[k];

        t_2[k] = f_0 * sgp_2[k]
                 + f_3 * pc_x[k] * shp_2[k];

        t_3[k] = f_1 * shs0_0[k]
                 - f_2 * shs1_0[k]
                 + f_3 * pc_y[k] * shp_1[k];

        t_4[k] = f_3 * pc_y[k] * shp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, sgd0_0, sgp_4, sgd1_0, shs0_0, \
                         shs1_0, shp_2, shp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * shs0_0[k]
                 - f_2 * shs1_0[k]
                 + f_3 * pc_z[k] * shp_2[k];

        t_6[k] = pb_y[k] * sgd0_0[k]
                 - f_4 * pc_y[k] * sgd1_0[k];

        t_7[k] = f_5 * sgp_4[k]
                 + f_3 * pc_x[k] * shp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, sgd0_5, sgp_1, sgp_2, sgp_5, \
                         sgd1_5, shs0_1, shs1_1, shp_4, shp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * sgp_5[k]
                 + f_3 * pc_x[k] * shp_5[k];

        t_9[k] = f_6 * sgp_1[k]
                 + f_1 * shs0_1[k]
                 - f_2 * shs1_1[k]
                 + f_3 * pc_y[k] * shp_4[k];

        t_10[k] = f_6 * sgp_2[k]
                  + f_3 * pc_y[k] * shp_5[k];

        t_11[k] = pb_y[k] * sgd0_5[k]
                  - f_4 * pc_y[k] * sgd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, sgd0_0, sgd0_3, sgp_7, \
                         sgp_8, sgd1_0, sgd1_3, shp_7, shp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * sgd0_0[k]
                  - f_4 * pc_z[k] * sgd1_0[k];

        t_13[k] = f_5 * sgp_7[k]
                  + f_3 * pc_x[k] * shp_7[k];

        t_14[k] = f_5 * sgp_8[k]
                  + f_3 * pc_x[k] * shp_8[k];

        t_15[k] = pb_z[k] * sgd0_3[k]
                  - f_4 * pc_z[k] * sgd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, sgp_2, sgp_9, shs0_2, shs0_3, \
                         shs1_2, shs1_3, shp_8, shp_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * shp_8[k];

        t_17[k] = f_6 * sgp_2[k]
                  + f_1 * shs0_2[k]
                  - f_2 * shs1_2[k]
                  + f_3 * pc_z[k] * shp_8[k];

        t_18[k] = f_7 * sgp_9[k]
                  + f_1 * shs0_3[k]
                  - f_2 * shs1_3[k]
                  + f_3 * pc_x[k] * shp_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, sgp_4, sgp_5, sgp_10, \
                         sgp_11, shs0_3, shs1_3, shp_10, shp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * sgp_10[k]
                  + f_3 * pc_x[k] * shp_10[k];

        t_20[k] = f_7 * sgp_11[k]
                  + f_3 * pc_x[k] * shp_11[k];

        t_21[k] = f_8 * sgp_4[k]
                  + f_1 * shs0_3[k]
                  - f_2 * shs1_3[k]
                  + f_3 * pc_y[k] * shp_10[k];

        t_22[k] = f_8 * sgp_5[k]
                  + f_3 * pc_y[k] * shp_11[k];

        t_23[k] = f_1 * shs0_3[k]
                  - f_2 * shs1_3[k]
                  + f_3 * pc_z[k] * shp_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, sgd0_12, sgp_13, sgp_14, sgd1_12, \
                         shp_13, shp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * sgd0_12[k]
                  - f_4 * pc_y[k] * sgd1_12[k];

        t_25[k] = f_7 * sgp_13[k]
                  + f_3 * pc_x[k] * shp_13[k];

        t_26[k] = f_7 * sgp_14[k]
                  + f_3 * pc_x[k] * shp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, sgd0_9, sgd0_17, sgp_8, \
                         sgd1_9, sgd1_17, shp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * sgd0_9[k]
                  - f_4 * pc_z[k] * sgd1_9[k];

        t_28[k] = f_6 * sgp_8[k]
                  + f_3 * pc_y[k] * shp_14[k];

        t_29[k] = pb_y[k] * sgd0_17[k]
                  - f_4 * pc_y[k] * sgd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, sgp_15, sgp_16, sgp_17, \
                         shs0_5, shs1_5, shp_15, shp_16, shp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * sgp_15[k]
                  + f_1 * shs0_5[k]
                  - f_2 * shs1_5[k]
                  + f_3 * pc_x[k] * shp_15[k];

        t_31[k] = f_7 * sgp_16[k]
                  + f_3 * pc_x[k] * shp_16[k];

        t_32[k] = f_7 * sgp_17[k]
                  + f_3 * pc_x[k] * shp_17[k];

        t_33[k] = f_1 * shs0_5[k]
                  - f_2 * shs1_5[k]
                  + f_3 * pc_y[k] * shp_16[k];

        t_34[k] = f_3 * pc_y[k] * shp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, sgp_8, sgp_18, sgp_19, shs0_5, shs0_6, \
                         shs1_5, shs1_6, shp_17, shp_18, shp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * sgp_8[k]
                  + f_1 * shs0_5[k]
                  - f_2 * shs1_5[k]
                  + f_3 * pc_z[k] * shp_17[k];

        t_36[k] = f_8 * sgp_18[k]
                  + f_1 * shs0_6[k]
                  - f_2 * shs1_6[k]
                  + f_3 * pc_x[k] * shp_18[k];

        t_37[k] = f_8 * sgp_19[k]
                  + f_3 * pc_x[k] * shp_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, sgp_10, sgp_11, sgp_20, \
                         shs0_6, shs1_6, shp_19, shp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_8 * sgp_20[k]
                  + f_3 * pc_x[k] * shp_20[k];

        t_39[k] = f_7 * sgp_10[k]
                  + f_1 * shs0_6[k]
                  - f_2 * shs1_6[k]
                  + f_3 * pc_y[k] * shp_19[k];

        t_40[k] = f_7 * sgp_11[k]
                  + f_3 * pc_y[k] * shp_20[k];

        t_41[k] = f_1 * shs0_6[k]
                  - f_2 * shs1_6[k]
                  + f_3 * pc_z[k] * shp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, sgd0_18, sgd0_21, sgp_22, \
                         sgp_23, sgd1_18, sgd1_21, shp_22, shp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sgd0_18[k]
                  - f_4 * pc_z[k] * sgd1_18[k];

        t_43[k] = f_8 * sgp_22[k]
                  + f_3 * pc_x[k] * shp_22[k];

        t_44[k] = f_8 * sgp_23[k]
                  + f_3 * pc_x[k] * shp_23[k];

        t_45[k] = pb_z[k] * sgd0_21[k]
                  - f_4 * pc_z[k] * sgd1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, sgd0_30, sgp_11, sgp_14, sgd1_30, \
                         shs0_7, shs1_7, shp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * sgp_14[k]
                  + f_3 * pc_y[k] * shp_23[k];

        t_47[k] = f_6 * sgp_11[k]
                  + f_1 * shs0_7[k]
                  - f_2 * shs1_7[k]
                  + f_3 * pc_z[k] * shp_23[k];

        t_48[k] = pb_y[k] * sgd0_30[k]
                  - f_4 * pc_y[k] * sgd1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, sgp_16, sgp_17, sgp_25, sgp_26, \
                         shs0_8, shs1_8, shp_25, shp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sgp_25[k]
                  + f_3 * pc_x[k] * shp_25[k];

        t_50[k] = f_8 * sgp_26[k]
                  + f_3 * pc_x[k] * shp_26[k];

        t_51[k] = f_6 * sgp_16[k]
                  + f_1 * shs0_8[k]
                  - f_2 * shs1_8[k]
                  + f_3 * pc_y[k] * shp_25[k];

        t_52[k] = f_6 * sgp_17[k]
                  + f_3 * pc_y[k] * shp_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, sgd0_35, sgp_27, sgp_28, sgd1_35, \
                         shs0_9, shs1_9, shp_27, shp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * sgd0_35[k]
                  - f_4 * pc_y[k] * sgd1_35[k];

        t_54[k] = f_8 * sgp_27[k]
                  + f_1 * shs0_9[k]
                  - f_2 * shs1_9[k]
                  + f_3 * pc_x[k] * shp_27[k];

        t_55[k] = f_8 * sgp_28[k]
                  + f_3 * pc_x[k] * shp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sgp_17, sgp_29, shs0_9, \
                         shs1_9, shp_28, shp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_8 * sgp_29[k]
                  + f_3 * pc_x[k] * shp_29[k];

        t_57[k] = f_1 * shs0_9[k]
                  - f_2 * shs1_9[k]
                  + f_3 * pc_y[k] * shp_28[k];

        t_58[k] = f_3 * pc_y[k] * shp_29[k];

        t_59[k] = f_7 * sgp_17[k]
                  + f_1 * shs0_9[k]
                  - f_2 * shs1_9[k]
                  + f_3 * pc_z[k] * shp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pb_x, pc_x, sgd0_60, sgd0_63, sgp_30, sgp_31, \
                         sgp_32, sgd1_60, sgd1_63, shp_31, shp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pb_x[k] * sgd0_60[k]
                  + f_8 * sgp_30[k]
                  - f_4 * pc_x[k] * sgd1_60[k];

        t_61[k] = f_6 * sgp_31[k]
                  + f_3 * pc_x[k] * shp_31[k];

        t_62[k] = f_6 * sgp_32[k]
                  + f_3 * pc_x[k] * shp_32[k];

        t_63[k] = pb_x[k] * sgd0_63[k]
                  - f_4 * pc_x[k] * sgd1_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_z, pc_x, pc_y, pc_z, sgd0_36, sgd0_65, \
                         sgp_20, sgd1_36, sgd1_65, shp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_5 * sgp_20[k]
                  + f_3 * pc_y[k] * shp_32[k];

        t_65[k] = pb_x[k] * sgd0_65[k]
                  - f_4 * pc_x[k] * sgd1_65[k];

        t_66[k] = pb_z[k] * sgd0_36[k]
                  - f_4 * pc_z[k] * sgd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pb_x, pc_x, pc_y, sgd0_69, sgp_23, sgp_34, \
                         sgp_35, sgd1_69, shp_34, shp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_6 * sgp_34[k]
                  + f_3 * pc_x[k] * shp_34[k];

        t_68[k] = f_6 * sgp_35[k]
                  + f_3 * pc_x[k] * shp_35[k];

        t_69[k] = pb_x[k] * sgd0_69[k]
                  - f_4 * pc_x[k] * sgd1_69[k];

        t_70[k] = f_7 * sgp_23[k]
                  + f_3 * pc_y[k] * shp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pb_x, pc_x, sgd0_71, sgd0_72, sgp_36, sgp_37, \
                         sgp_38, sgd1_71, sgd1_72, shp_37, shp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_x[k] * sgd0_71[k]
                  - f_4 * pc_x[k] * sgd1_71[k];

        t_72[k] = pb_x[k] * sgd0_72[k]
                  + f_8 * sgp_36[k]
                  - f_4 * pc_x[k] * sgd1_72[k];

        t_73[k] = f_6 * sgp_37[k]
                  + f_3 * pc_x[k] * shp_37[k];

        t_74[k] = f_6 * sgp_38[k]
                  + f_3 * pc_x[k] * shp_38[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_x, pb_y, pc_x, pc_y, sgd0_54, sgd0_75, \
                         sgd0_77, sgp_26, sgd1_54, sgd1_75, sgd1_77, \
                         shp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_x[k] * sgd0_75[k]
                  - f_4 * pc_x[k] * sgd1_75[k];

        t_76[k] = f_8 * sgp_26[k]
                  + f_3 * pc_y[k] * shp_38[k];

        t_77[k] = pb_x[k] * sgd0_77[k]
                  - f_4 * pc_x[k] * sgd1_77[k];

        t_78[k] = pb_y[k] * sgd0_54[k]
                  - f_4 * pc_y[k] * sgd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_x, pc_x, pc_y, sgd0_81, sgp_29, sgp_40, \
                         sgp_41, sgd1_81, shp_40, shp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_6 * sgp_40[k]
                  + f_3 * pc_x[k] * shp_40[k];

        t_80[k] = f_6 * sgp_41[k]
                  + f_3 * pc_x[k] * shp_41[k];

        t_81[k] = pb_x[k] * sgd0_81[k]
                  - f_4 * pc_x[k] * sgd1_81[k];

        t_82[k] = f_6 * sgp_29[k]
                  + f_3 * pc_y[k] * shp_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_x, pc_x, sgd0_83, sgd0_84, sgp_42, sgp_43, \
                         sgp_44, sgd1_83, sgd1_84, shp_43, shp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_x[k] * sgd0_83[k]
                  - f_4 * pc_x[k] * sgd1_83[k];

        t_84[k] = pb_x[k] * sgd0_84[k]
                  + f_8 * sgp_42[k]
                  - f_4 * pc_x[k] * sgd1_84[k];

        t_85[k] = f_6 * sgp_43[k]
                  + f_3 * pc_x[k] * shp_43[k];

        t_86[k] = f_6 * sgp_44[k]
                  + f_3 * pc_x[k] * shp_44[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pb_x, pc_x, pc_y, sgd0_87, sgd0_89, sgd1_87, \
                         sgd1_89, shs0_15, shs1_15, shp_44, shp_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_x[k] * sgd0_87[k]
                  - f_4 * pc_x[k] * sgd1_87[k];

        t_88[k] = f_3 * pc_y[k] * shp_44[k];

        t_89[k] = pb_x[k] * sgd0_89[k]
                  - f_4 * pc_x[k] * sgd1_89[k];

        t_90[k] = f_1 * shs0_15[k]
                  - f_2 * shs1_15[k]
                  + f_3 * pc_x[k] * shp_45[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, sgp_31, sgp_32, \
                         shs0_15, shs1_15, shp_46, shp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * pc_x[k] * shp_46[k];

        t_92[k] = f_3 * pc_x[k] * shp_47[k];

        t_93[k] = f_0 * sgp_31[k]
                  + f_1 * shs0_15[k]
                  - f_2 * shs1_15[k]
                  + f_3 * pc_y[k] * shp_46[k];

        t_94[k] = f_0 * sgp_32[k]
                  + f_3 * pc_y[k] * shp_47[k];

        t_95[k] = f_1 * shs0_15[k]
                  - f_2 * shs1_15[k]
                  + f_3 * pc_z[k] * shp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_z, pc_x, pc_y, pc_z, sgd0_60, \
                         sgd0_63, sgp_35, sgd1_60, sgd1_63, shp_49, \
                         shp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sgd0_60[k]
                  - f_4 * pc_z[k] * sgd1_60[k];

        t_97[k] = f_3 * pc_x[k] * shp_49[k];

        t_98[k] = f_3 * pc_x[k] * shp_50[k];

        t_99[k] = pb_z[k] * sgd0_63[k]
                  - f_4 * pc_z[k] * sgd1_63[k];

        t_100[k] = f_5 * sgp_35[k]
                   + f_3 * pc_y[k] * shp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_x, pc_z, sgp_32, shs0_16, shs0_17, \
                         shs1_16, shs1_17, shp_50, shp_51, shp_52, \
                         shp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_6 * sgp_32[k]
                   + f_1 * shs0_16[k]
                   - f_2 * shs1_16[k]
                   + f_3 * pc_z[k] * shp_50[k];

        t_102[k] = f_1 * shs0_17[k]
                   - f_2 * shs1_17[k]
                   + f_3 * pc_x[k] * shp_51[k];

        t_103[k] = f_3 * pc_x[k] * shp_52[k];

        t_104[k] = f_3 * pc_x[k] * shp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pc_y, pc_z, sgp_35, sgp_37, sgp_38, shs0_17, \
                         shs1_17, shp_52, shp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_7 * sgp_37[k]
                   + f_1 * shs0_17[k]
                   - f_2 * shs1_17[k]
                   + f_3 * pc_y[k] * shp_52[k];

        t_106[k] = f_7 * sgp_38[k]
                   + f_3 * pc_y[k] * shp_53[k];

        t_107[k] = f_8 * sgp_35[k]
                   + f_1 * shs0_17[k]
                   - f_2 * shs1_17[k]
                   + f_3 * pc_z[k] * shp_53[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pc_x, pc_y, sgp_40, sgp_41, \
                         shs0_18, shs1_18, shp_54, shp_55, shp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_1 * shs0_18[k]
                   - f_2 * shs1_18[k]
                   + f_3 * pc_x[k] * shp_54[k];

        t_109[k] = f_3 * pc_x[k] * shp_55[k];

        t_110[k] = f_3 * pc_x[k] * shp_56[k];

        t_111[k] = f_8 * sgp_40[k]
                   + f_1 * shs0_18[k]
                   - f_2 * shs1_18[k]
                   + f_3 * pc_y[k] * shp_55[k];

        t_112[k] = f_8 * sgp_41[k]
                   + f_3 * pc_y[k] * shp_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_y, pc_x, pc_y, pc_z, sgd0_84, sgp_38, \
                         sgd1_84, shs0_18, shs1_18, shp_56, shp_58, \
                         shp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_7 * sgp_38[k]
                   + f_1 * shs0_18[k]
                   - f_2 * shs1_18[k]
                   + f_3 * pc_z[k] * shp_56[k];

        t_114[k] = pb_y[k] * sgd0_84[k]
                   - f_4 * pc_y[k] * sgd1_84[k];

        t_115[k] = f_3 * pc_x[k] * shp_58[k];

        t_116[k] = f_3 * pc_x[k] * shp_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pb_y, pc_y, sgd0_87, sgd0_89, sgp_43, sgp_44, \
                         sgd1_87, sgd1_89, shp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_y[k] * sgd0_87[k]
                   + f_8 * sgp_43[k]
                   - f_4 * pc_y[k] * sgd1_87[k];

        t_118[k] = f_6 * sgp_44[k]
                   + f_3 * pc_y[k] * shp_59[k];

        t_119[k] = pb_y[k] * sgd0_89[k]
                   - f_4 * pc_y[k] * sgd1_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, t_125, pc_x, pc_y, pc_z, sgp_44, \
                         shs0_20, shs1_20, shp_60, shp_61, shp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * shs0_20[k]
                   - f_2 * shs1_20[k]
                   + f_3 * pc_x[k] * shp_60[k];

        t_121[k] = f_3 * pc_x[k] * shp_61[k];

        t_122[k] = f_3 * pc_x[k] * shp_62[k];

        t_123[k] = f_1 * shs0_20[k]
                   - f_2 * shs1_20[k]
                   + f_3 * pc_y[k] * shp_61[k];

        t_124[k] = f_3 * pc_y[k] * shp_62[k];

        t_125[k] = f_0 * sgp_44[k]
                   + f_1 * shs0_20[k]
                   - f_2 * shs1_20[k]
                   + f_3 * pc_z[k] * shp_62[k];
    }
}

}  // namespace simdt3ceri
