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


#include "SimdThreeCenterElectronRepulsionVrrRecSKL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_0 = buffer.data(sil0 + 0);
    const auto *sil0_3 = buffer.data(sil0 + 3);
    const auto *sil0_5 = buffer.data(sil0 + 5);
    const auto *sil0_6 = buffer.data(sil0 + 6);
    const auto *sil0_9 = buffer.data(sil0 + 9);
    const auto *sil0_10 = buffer.data(sil0 + 10);
    const auto *sil0_12 = buffer.data(sil0 + 12);
    const auto *sil0_14 = buffer.data(sil0 + 14);
    const auto *sil0_15 = buffer.data(sil0 + 15);
    const auto *sil0_17 = buffer.data(sil0 + 17);
    const auto *sil0_18 = buffer.data(sil0 + 18);
    const auto *sil0_20 = buffer.data(sil0 + 20);
    const auto *sil0_21 = buffer.data(sil0 + 21);
    const auto *sil0_23 = buffer.data(sil0 + 23);
    const auto *sil0_24 = buffer.data(sil0 + 24);
    const auto *sil0_25 = buffer.data(sil0 + 25);
    const auto *sil0_27 = buffer.data(sil0 + 27);
    const auto *sil0_44 = buffer.data(sil0 + 44);

    const auto *sik_0 = buffer.data(sik + 0);
    const auto *sik_1 = buffer.data(sik + 1);
    const auto *sik_2 = buffer.data(sik + 2);
    const auto *sik_3 = buffer.data(sik + 3);
    const auto *sik_5 = buffer.data(sik + 5);
    const auto *sik_6 = buffer.data(sik + 6);
    const auto *sik_7 = buffer.data(sik + 7);
    const auto *sik_8 = buffer.data(sik + 8);
    const auto *sik_9 = buffer.data(sik + 9);
    const auto *sik_10 = buffer.data(sik + 10);
    const auto *sik_11 = buffer.data(sik + 11);
    const auto *sik_12 = buffer.data(sik + 12);
    const auto *sik_13 = buffer.data(sik + 13);
    const auto *sik_14 = buffer.data(sik + 14);
    const auto *sik_15 = buffer.data(sik + 15);
    const auto *sik_16 = buffer.data(sik + 16);
    const auto *sik_17 = buffer.data(sik + 17);
    const auto *sik_18 = buffer.data(sik + 18);
    const auto *sik_19 = buffer.data(sik + 19);
    const auto *sik_20 = buffer.data(sik + 20);
    const auto *sik_21 = buffer.data(sik + 21);
    const auto *sik_23 = buffer.data(sik + 23);
    const auto *sik_24 = buffer.data(sik + 24);
    const auto *sik_25 = buffer.data(sik + 25);
    const auto *sik_27 = buffer.data(sik + 27);
    const auto *sik_28 = buffer.data(sik + 28);
    const auto *sik_29 = buffer.data(sik + 29);
    const auto *sik_30 = buffer.data(sik + 30);
    const auto *sik_31 = buffer.data(sik + 31);
    const auto *sik_32 = buffer.data(sik + 32);
    const auto *sik_33 = buffer.data(sik + 33);
    const auto *sik_34 = buffer.data(sik + 34);
    const auto *sik_35 = buffer.data(sik + 35);
    const auto *sik_64 = buffer.data(sik + 64);
    const auto *sik_65 = buffer.data(sik + 65);
    const auto *sik_66 = buffer.data(sik + 66);
    const auto *sik_67 = buffer.data(sik + 67);
    const auto *sik_68 = buffer.data(sik + 68);
    const auto *sik_69 = buffer.data(sik + 69);
    const auto *sik_70 = buffer.data(sik + 70);
    const auto *sik_71 = buffer.data(sik + 71);

    const auto *sil1_0 = buffer.data(sil1 + 0);
    const auto *sil1_3 = buffer.data(sil1 + 3);
    const auto *sil1_5 = buffer.data(sil1 + 5);
    const auto *sil1_6 = buffer.data(sil1 + 6);
    const auto *sil1_9 = buffer.data(sil1 + 9);
    const auto *sil1_10 = buffer.data(sil1 + 10);
    const auto *sil1_12 = buffer.data(sil1 + 12);
    const auto *sil1_14 = buffer.data(sil1 + 14);
    const auto *sil1_15 = buffer.data(sil1 + 15);
    const auto *sil1_17 = buffer.data(sil1 + 17);
    const auto *sil1_18 = buffer.data(sil1 + 18);
    const auto *sil1_20 = buffer.data(sil1 + 20);
    const auto *sil1_21 = buffer.data(sil1 + 21);
    const auto *sil1_23 = buffer.data(sil1 + 23);
    const auto *sil1_24 = buffer.data(sil1 + 24);
    const auto *sil1_25 = buffer.data(sil1 + 25);
    const auto *sil1_27 = buffer.data(sil1 + 27);
    const auto *sil1_44 = buffer.data(sil1 + 44);

    const auto *ski0_0 = buffer.data(ski0 + 0);
    const auto *ski0_3 = buffer.data(ski0 + 3);
    const auto *ski0_5 = buffer.data(ski0 + 5);
    const auto *ski0_6 = buffer.data(ski0 + 6);
    const auto *ski0_9 = buffer.data(ski0 + 9);
    const auto *ski0_10 = buffer.data(ski0 + 10);
    const auto *ski0_12 = buffer.data(ski0 + 12);
    const auto *ski0_14 = buffer.data(ski0 + 14);
    const auto *ski0_15 = buffer.data(ski0 + 15);
    const auto *ski0_17 = buffer.data(ski0 + 17);
    const auto *ski0_18 = buffer.data(ski0 + 18);
    const auto *ski0_20 = buffer.data(ski0 + 20);
    const auto *ski0_21 = buffer.data(ski0 + 21);
    const auto *ski0_23 = buffer.data(ski0 + 23);
    const auto *ski0_24 = buffer.data(ski0 + 24);
    const auto *ski0_25 = buffer.data(ski0 + 25);
    const auto *ski0_26 = buffer.data(ski0 + 26);
    const auto *ski0_27 = buffer.data(ski0 + 27);
    const auto *ski0_49 = buffer.data(ski0 + 49);
    const auto *ski0_51 = buffer.data(ski0 + 51);
    const auto *ski0_52 = buffer.data(ski0 + 52);
    const auto *ski0_53 = buffer.data(ski0 + 53);
    const auto *ski0_54 = buffer.data(ski0 + 54);
    const auto *ski0_55 = buffer.data(ski0 + 55);

    const auto *ski1_0 = buffer.data(ski1 + 0);
    const auto *ski1_3 = buffer.data(ski1 + 3);
    const auto *ski1_5 = buffer.data(ski1 + 5);
    const auto *ski1_6 = buffer.data(ski1 + 6);
    const auto *ski1_9 = buffer.data(ski1 + 9);
    const auto *ski1_10 = buffer.data(ski1 + 10);
    const auto *ski1_12 = buffer.data(ski1 + 12);
    const auto *ski1_14 = buffer.data(ski1 + 14);
    const auto *ski1_15 = buffer.data(ski1 + 15);
    const auto *ski1_17 = buffer.data(ski1 + 17);
    const auto *ski1_18 = buffer.data(ski1 + 18);
    const auto *ski1_20 = buffer.data(ski1 + 20);
    const auto *ski1_21 = buffer.data(ski1 + 21);
    const auto *ski1_23 = buffer.data(ski1 + 23);
    const auto *ski1_24 = buffer.data(ski1 + 24);
    const auto *ski1_25 = buffer.data(ski1 + 25);
    const auto *ski1_26 = buffer.data(ski1 + 26);
    const auto *ski1_27 = buffer.data(ski1 + 27);
    const auto *ski1_49 = buffer.data(ski1 + 49);
    const auto *ski1_51 = buffer.data(ski1 + 51);
    const auto *ski1_52 = buffer.data(ski1 + 52);
    const auto *ski1_53 = buffer.data(ski1 + 53);
    const auto *ski1_54 = buffer.data(ski1 + 54);
    const auto *ski1_55 = buffer.data(ski1 + 55);

    const auto *skk_0 = buffer.data(skk + 0);
    const auto *skk_2 = buffer.data(skk + 2);
    const auto *skk_3 = buffer.data(skk + 3);
    const auto *skk_5 = buffer.data(skk + 5);
    const auto *skk_6 = buffer.data(skk + 6);
    const auto *skk_9 = buffer.data(skk + 9);
    const auto *skk_10 = buffer.data(skk + 10);
    const auto *skk_12 = buffer.data(skk + 12);
    const auto *skk_14 = buffer.data(skk + 14);
    const auto *skk_15 = buffer.data(skk + 15);
    const auto *skk_17 = buffer.data(skk + 17);
    const auto *skk_18 = buffer.data(skk + 18);
    const auto *skk_20 = buffer.data(skk + 20);
    const auto *skk_21 = buffer.data(skk + 21);
    const auto *skk_23 = buffer.data(skk + 23);
    const auto *skk_24 = buffer.data(skk + 24);
    const auto *skk_25 = buffer.data(skk + 25);
    const auto *skk_27 = buffer.data(skk + 27);
    const auto *skk_28 = buffer.data(skk + 28);
    const auto *skk_29 = buffer.data(skk + 29);
    const auto *skk_30 = buffer.data(skk + 30);
    const auto *skk_31 = buffer.data(skk + 31);
    const auto *skk_32 = buffer.data(skk + 32);
    const auto *skk_33 = buffer.data(skk + 33);
    const auto *skk_34 = buffer.data(skk + 34);
    const auto *skk_35 = buffer.data(skk + 35);
    const auto *skk_36 = buffer.data(skk + 36);
    const auto *skk_38 = buffer.data(skk + 38);
    const auto *skk_39 = buffer.data(skk + 39);
    const auto *skk_41 = buffer.data(skk + 41);
    const auto *skk_42 = buffer.data(skk + 42);
    const auto *skk_45 = buffer.data(skk + 45);
    const auto *skk_46 = buffer.data(skk + 46);
    const auto *skk_50 = buffer.data(skk + 50);
    const auto *skk_51 = buffer.data(skk + 51);
    const auto *skk_56 = buffer.data(skk + 56);
    const auto *skk_64 = buffer.data(skk + 64);
    const auto *skk_65 = buffer.data(skk + 65);
    const auto *skk_66 = buffer.data(skk + 66);
    const auto *skk_67 = buffer.data(skk + 67);
    const auto *skk_68 = buffer.data(skk + 68);
    const auto *skk_69 = buffer.data(skk + 69);
    const auto *skk_70 = buffer.data(skk + 70);
    const auto *skk_71 = buffer.data(skk + 71);
    const auto *skk_72 = buffer.data(skk + 72);
    const auto *skk_74 = buffer.data(skk + 74);
    const auto *skk_75 = buffer.data(skk + 75);
    const auto *skk_77 = buffer.data(skk + 77);
    const auto *skk_78 = buffer.data(skk + 78);
    const auto *skk_81 = buffer.data(skk + 81);
    const auto *skk_82 = buffer.data(skk + 82);
    const auto *skk_86 = buffer.data(skk + 86);
    const auto *skk_87 = buffer.data(skk + 87);
    const auto *skk_92 = buffer.data(skk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sik_0, sik_3, ski0_0, ski0_3, \
                         ski1_0, ski1_3, skk_0, skk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sik_0[k]
                 + f_1 * ski0_0[k]
                 - f_2 * ski1_0[k]
                 + f_3 * pc_x[k] * skk_0[k];

        t_1[k] = f_3 * pc_y[k] * skk_0[k];

        t_2[k] = f_3 * pc_z[k] * skk_0[k];

        t_3[k] = f_0 * sik_3[k]
                 + f_4 * ski0_3[k]
                 - f_5 * ski1_3[k]
                 + f_3 * pc_x[k] * skk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, sik_5, sik_6, ski0_5, ski0_6, ski1_5, \
                         ski1_6, skk_2, skk_5, skk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * skk_2[k];

        t_5[k] = f_0 * sik_5[k]
                 + f_4 * ski0_5[k]
                 - f_5 * ski1_5[k]
                 + f_3 * pc_x[k] * skk_5[k];

        t_6[k] = f_0 * sik_6[k]
                 + f_6 * ski0_6[k]
                 - f_7 * ski1_6[k]
                 + f_3 * pc_x[k] * skk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, sik_9, ski0_9, ski1_9, skk_3, skk_5, \
                         skk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * skk_3[k];

        t_8[k] = f_3 * pc_y[k] * skk_5[k];

        t_9[k] = f_0 * sik_9[k]
                 + f_6 * ski0_9[k]
                 - f_7 * ski1_9[k]
                 + f_3 * pc_x[k] * skk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, sik_10, sik_12, ski0_10, ski0_12, \
                         ski1_10, ski1_12, skk_6, skk_10, skk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * sik_10[k]
                  + f_8 * ski0_10[k]
                  - f_9 * ski1_10[k]
                  + f_3 * pc_x[k] * skk_10[k];

        t_11[k] = f_3 * pc_z[k] * skk_6[k];

        t_12[k] = f_0 * sik_12[k]
                  + f_8 * ski0_12[k]
                  - f_9 * ski1_12[k]
                  + f_3 * pc_x[k] * skk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, sik_14, sik_15, ski0_14, ski0_15, \
                         ski1_14, ski1_15, skk_9, skk_14, skk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * skk_9[k];

        t_14[k] = f_0 * sik_14[k]
                  + f_8 * ski0_14[k]
                  - f_9 * ski1_14[k]
                  + f_3 * pc_x[k] * skk_14[k];

        t_15[k] = f_0 * sik_15[k]
                  + f_10 * ski0_15[k]
                  - f_11 * ski1_15[k]
                  + f_3 * pc_x[k] * skk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, sik_17, sik_18, ski0_17, ski0_18, \
                         ski1_17, ski1_18, skk_10, skk_17, skk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * skk_10[k];

        t_17[k] = f_0 * sik_17[k]
                  + f_10 * ski0_17[k]
                  - f_11 * ski1_17[k]
                  + f_3 * pc_x[k] * skk_17[k];

        t_18[k] = f_0 * sik_18[k]
                  + f_10 * ski0_18[k]
                  - f_11 * ski1_18[k]
                  + f_3 * pc_x[k] * skk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, sik_20, sik_21, ski0_20, ski0_21, \
                         ski1_20, ski1_21, skk_14, skk_20, skk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * skk_14[k];

        t_20[k] = f_0 * sik_20[k]
                  + f_10 * ski0_20[k]
                  - f_11 * ski1_20[k]
                  + f_3 * pc_x[k] * skk_20[k];

        t_21[k] = f_0 * sik_21[k]
                  + f_12 * ski0_21[k]
                  - f_13 * ski1_21[k]
                  + f_3 * pc_x[k] * skk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, sik_23, sik_24, ski0_23, ski0_24, \
                         ski1_23, ski1_24, skk_15, skk_23, skk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * skk_15[k];

        t_23[k] = f_0 * sik_23[k]
                  + f_12 * ski0_23[k]
                  - f_13 * ski1_23[k]
                  + f_3 * pc_x[k] * skk_23[k];

        t_24[k] = f_0 * sik_24[k]
                  + f_12 * ski0_24[k]
                  - f_13 * ski1_24[k]
                  + f_3 * pc_x[k] * skk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, sik_25, sik_27, ski0_25, ski0_27, \
                         ski1_25, ski1_27, skk_20, skk_25, skk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sik_25[k]
                  + f_12 * ski0_25[k]
                  - f_13 * ski1_25[k]
                  + f_3 * pc_x[k] * skk_25[k];

        t_26[k] = f_3 * pc_y[k] * skk_20[k];

        t_27[k] = f_0 * sik_27[k]
                  + f_12 * ski0_27[k]
                  - f_13 * ski1_27[k]
                  + f_3 * pc_x[k] * skk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, sik_28, sik_29, sik_30, sik_31, \
                         sik_32, skk_28, skk_29, skk_30, skk_31, \
                         skk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * sik_28[k]
                  + f_3 * pc_x[k] * skk_28[k];

        t_29[k] = f_0 * sik_29[k]
                  + f_3 * pc_x[k] * skk_29[k];

        t_30[k] = f_0 * sik_30[k]
                  + f_3 * pc_x[k] * skk_30[k];

        t_31[k] = f_0 * sik_31[k]
                  + f_3 * pc_x[k] * skk_31[k];

        t_32[k] = f_0 * sik_32[k]
                  + f_3 * pc_x[k] * skk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, sik_33, sik_34, sik_35, ski0_21, \
                         ski1_21, skk_28, skk_33, skk_34, skk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * sik_33[k]
                  + f_3 * pc_x[k] * skk_33[k];

        t_34[k] = f_0 * sik_34[k]
                  + f_3 * pc_x[k] * skk_34[k];

        t_35[k] = f_0 * sik_35[k]
                  + f_3 * pc_x[k] * skk_35[k];

        t_36[k] = f_1 * ski0_21[k]
                  - f_2 * ski1_21[k]
                  + f_3 * pc_y[k] * skk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, ski0_23, ski0_24, ski0_25, \
                         ski1_23, ski1_24, ski1_25, skk_28, skk_30, skk_31, \
                         skk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * skk_28[k];

        t_38[k] = f_4 * ski0_23[k]
                  - f_5 * ski1_23[k]
                  + f_3 * pc_y[k] * skk_30[k];

        t_39[k] = f_6 * ski0_24[k]
                  - f_7 * ski1_24[k]
                  + f_3 * pc_y[k] * skk_31[k];

        t_40[k] = f_8 * ski0_25[k]
                  - f_9 * ski1_25[k]
                  + f_3 * pc_y[k] * skk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, ski0_26, ski0_27, ski1_26, \
                         ski1_27, skk_33, skk_34, skk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * ski0_26[k]
                  - f_11 * ski1_26[k]
                  + f_3 * pc_y[k] * skk_33[k];

        t_42[k] = f_12 * ski0_27[k]
                  - f_13 * ski1_27[k]
                  + f_3 * pc_y[k] * skk_34[k];

        t_43[k] = f_3 * pc_y[k] * skk_35[k];

        t_44[k] = f_1 * ski0_27[k]
                  - f_2 * ski1_27[k]
                  + f_3 * pc_z[k] * skk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, sil0_0, sil0_3, sik_0, \
                         sik_1, sil1_0, sil1_3, skk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * sil0_0[k]
                  - f_14 * pc_y[k] * sil1_0[k];

        t_46[k] = f_15 * sik_0[k]
                  + f_3 * pc_y[k] * skk_36[k];

        t_47[k] = f_3 * pc_z[k] * skk_36[k];

        t_48[k] = pb_y[k] * sil0_3[k]
                  + f_16 * sik_1[k]
                  - f_14 * pc_y[k] * sil1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, sil0_5, sil0_6, sik_2, \
                         sik_3, sil1_5, sil1_6, skk_38, skk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * sik_2[k]
                  + f_3 * pc_y[k] * skk_38[k];

        t_50[k] = pb_y[k] * sil0_5[k]
                  - f_14 * pc_y[k] * sil1_5[k];

        t_51[k] = pb_y[k] * sil0_6[k]
                  + f_17 * sik_3[k]
                  - f_14 * pc_y[k] * sil1_6[k];

        t_52[k] = f_3 * pc_z[k] * skk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, sil0_9, sil0_10, sik_5, \
                         sik_6, sil1_9, sil1_10, skk_41, skk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * sik_5[k]
                  + f_3 * pc_y[k] * skk_41[k];

        t_54[k] = pb_y[k] * sil0_9[k]
                  - f_14 * pc_y[k] * sil1_9[k];

        t_55[k] = pb_y[k] * sil0_10[k]
                  + f_18 * sik_6[k]
                  - f_14 * pc_y[k] * sil1_10[k];

        t_56[k] = f_3 * pc_z[k] * skk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, sil0_12, sil0_14, sil0_15, sik_8, \
                         sik_9, sik_10, sil1_12, sil1_14, sil1_15, \
                         skk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * sil0_12[k]
                  + f_16 * sik_8[k]
                  - f_14 * pc_y[k] * sil1_12[k];

        t_58[k] = f_15 * sik_9[k]
                  + f_3 * pc_y[k] * skk_45[k];

        t_59[k] = pb_y[k] * sil0_14[k]
                  - f_14 * pc_y[k] * sil1_14[k];

        t_60[k] = pb_y[k] * sil0_15[k]
                  + f_19 * sik_10[k]
                  - f_14 * pc_y[k] * sil1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, sil0_17, sil0_18, sik_12, \
                         sik_13, sik_14, sil1_17, sil1_18, skk_46, \
                         skk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * skk_46[k];

        t_62[k] = pb_y[k] * sil0_17[k]
                  + f_17 * sik_12[k]
                  - f_14 * pc_y[k] * sil1_17[k];

        t_63[k] = pb_y[k] * sil0_18[k]
                  + f_16 * sik_13[k]
                  - f_14 * pc_y[k] * sil1_18[k];

        t_64[k] = f_15 * sik_14[k]
                  + f_3 * pc_y[k] * skk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, sil0_20, sil0_21, sil0_23, \
                         sik_15, sik_17, sil1_20, sil1_21, sil1_23, \
                         skk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sil0_20[k]
                  - f_14 * pc_y[k] * sil1_20[k];

        t_66[k] = pb_y[k] * sil0_21[k]
                  + f_20 * sik_15[k]
                  - f_14 * pc_y[k] * sil1_21[k];

        t_67[k] = f_3 * pc_z[k] * skk_51[k];

        t_68[k] = pb_y[k] * sil0_23[k]
                  + f_18 * sik_17[k]
                  - f_14 * pc_y[k] * sil1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, sil0_24, sil0_25, sil0_27, \
                         sik_18, sik_19, sik_20, sil1_24, sil1_25, sil1_27, \
                         skk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * sil0_24[k]
                  + f_17 * sik_18[k]
                  - f_14 * pc_y[k] * sil1_24[k];

        t_70[k] = pb_y[k] * sil0_25[k]
                  + f_16 * sik_19[k]
                  - f_14 * pc_y[k] * sil1_25[k];

        t_71[k] = f_15 * sik_20[k]
                  + f_3 * pc_y[k] * skk_56[k];

        t_72[k] = pb_y[k] * sil0_27[k]
                  - f_14 * pc_y[k] * sil1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, sik_64, sik_65, sik_66, sik_67, \
                         sik_68, skk_64, skk_65, skk_66, skk_67, \
                         skk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_20 * sik_64[k]
                  + f_3 * pc_x[k] * skk_64[k];

        t_74[k] = f_20 * sik_65[k]
                  + f_3 * pc_x[k] * skk_65[k];

        t_75[k] = f_20 * sik_66[k]
                  + f_3 * pc_x[k] * skk_66[k];

        t_76[k] = f_20 * sik_67[k]
                  + f_3 * pc_x[k] * skk_67[k];

        t_77[k] = f_20 * sik_68[k]
                  + f_3 * pc_x[k] * skk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, sik_28, sik_69, sik_70, sik_71, \
                         ski0_49, ski1_49, skk_64, skk_69, skk_70, \
                         skk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_20 * sik_69[k]
                  + f_3 * pc_x[k] * skk_69[k];

        t_79[k] = f_20 * sik_70[k]
                  + f_3 * pc_x[k] * skk_70[k];

        t_80[k] = f_20 * sik_71[k]
                  + f_3 * pc_x[k] * skk_71[k];

        t_81[k] = f_15 * sik_28[k]
                  + f_1 * ski0_49[k]
                  - f_2 * ski1_49[k]
                  + f_3 * pc_y[k] * skk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, sik_30, sik_31, ski0_51, ski0_52, \
                         ski1_51, ski1_52, skk_64, skk_66, skk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * skk_64[k];

        t_83[k] = f_15 * sik_30[k]
                  + f_4 * ski0_51[k]
                  - f_5 * ski1_51[k]
                  + f_3 * pc_y[k] * skk_66[k];

        t_84[k] = f_15 * sik_31[k]
                  + f_6 * ski0_52[k]
                  - f_7 * ski1_52[k]
                  + f_3 * pc_y[k] * skk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, sik_32, sik_33, sik_34, ski0_53, ski0_54, \
                         ski0_55, ski1_53, ski1_54, ski1_55, skk_68, skk_69, \
                         skk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * sik_32[k]
                  + f_8 * ski0_53[k]
                  - f_9 * ski1_53[k]
                  + f_3 * pc_y[k] * skk_68[k];

        t_86[k] = f_15 * sik_33[k]
                  + f_10 * ski0_54[k]
                  - f_11 * ski1_54[k]
                  + f_3 * pc_y[k] * skk_69[k];

        t_87[k] = f_15 * sik_34[k]
                  + f_12 * ski0_55[k]
                  - f_13 * ski1_55[k]
                  + f_3 * pc_y[k] * skk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, sil0_0, sil0_44, \
                         sik_35, sil1_0, sil1_44, skk_71, skk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * sik_35[k]
                  + f_3 * pc_y[k] * skk_71[k];

        t_89[k] = pb_y[k] * sil0_44[k]
                  - f_14 * pc_y[k] * sil1_44[k];

        t_90[k] = pb_z[k] * sil0_0[k]
                  - f_14 * pc_z[k] * sil1_0[k];

        t_91[k] = f_3 * pc_y[k] * skk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, sil0_3, sil0_5, sik_0, \
                         sik_2, sil1_3, sil1_5, skk_72, skk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * sik_0[k]
                  + f_3 * pc_z[k] * skk_72[k];

        t_93[k] = pb_z[k] * sil0_3[k]
                  - f_14 * pc_z[k] * sil1_3[k];

        t_94[k] = f_3 * pc_y[k] * skk_74[k];

        t_95[k] = pb_z[k] * sil0_5[k]
                  + f_16 * sik_2[k]
                  - f_14 * pc_z[k] * sil1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, sil0_6, sil0_9, sik_3, \
                         sik_5, sil1_6, sil1_9, skk_75, skk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * sil0_6[k]
                  - f_14 * pc_z[k] * sil1_6[k];

        t_97[k] = f_15 * sik_3[k]
                  + f_3 * pc_z[k] * skk_75[k];

        t_98[k] = f_3 * pc_y[k] * skk_77[k];

        t_99[k] = pb_z[k] * sil0_9[k]
                  + f_17 * sik_5[k]
                  - f_14 * pc_z[k] * sil1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, sil0_10, sil0_12, \
                         sik_6, sik_7, sil1_10, sil1_12, skk_78, \
                         skk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * sil0_10[k]
                   - f_14 * pc_z[k] * sil1_10[k];

        t_101[k] = f_15 * sik_6[k]
                   + f_3 * pc_z[k] * skk_78[k];

        t_102[k] = pb_z[k] * sil0_12[k]
                   + f_16 * sik_7[k]
                   - f_14 * pc_z[k] * sil1_12[k];

        t_103[k] = f_3 * pc_y[k] * skk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, sil0_14, sil0_15, sil0_17, \
                         sik_9, sik_10, sik_11, sil1_14, sil1_15, sil1_17, \
                         skk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * sil0_14[k]
                   + f_18 * sik_9[k]
                   - f_14 * pc_z[k] * sil1_14[k];

        t_105[k] = pb_z[k] * sil0_15[k]
                   - f_14 * pc_z[k] * sil1_15[k];

        t_106[k] = f_15 * sik_10[k]
                   + f_3 * pc_z[k] * skk_82[k];

        t_107[k] = pb_z[k] * sil0_17[k]
                   + f_16 * sik_11[k]
                   - f_14 * pc_z[k] * sil1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, sil0_18, sil0_20, \
                         sil0_21, sik_12, sik_14, sil1_18, sil1_20, sil1_21, \
                         skk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * sil0_18[k]
                   + f_17 * sik_12[k]
                   - f_14 * pc_z[k] * sil1_18[k];

        t_109[k] = f_3 * pc_y[k] * skk_86[k];

        t_110[k] = pb_z[k] * sil0_20[k]
                   + f_19 * sik_14[k]
                   - f_14 * pc_z[k] * sil1_20[k];

        t_111[k] = pb_z[k] * sil0_21[k]
                   - f_14 * pc_z[k] * sil1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, sil0_23, sil0_24, sik_15, sik_16, \
                         sik_17, sil1_23, sil1_24, skk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * sik_15[k]
                   + f_3 * pc_z[k] * skk_87[k];

        t_113[k] = pb_z[k] * sil0_23[k]
                   + f_16 * sik_16[k]
                   - f_14 * pc_z[k] * sil1_23[k];

        t_114[k] = pb_z[k] * sil0_24[k]
                   + f_17 * sik_17[k]
                   - f_14 * pc_z[k] * sil1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, sil0_25, sil0_27, sik_18, \
                         sik_20, sil1_25, sil1_27, skk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * sil0_25[k]
                   + f_18 * sik_18[k]
                   - f_14 * pc_z[k] * sil1_25[k];

        t_116[k] = f_3 * pc_y[k] * skk_92[k];

        t_117[k] = pb_z[k] * sil0_27[k]
                   + f_20 * sik_20[k]
                   - f_14 * pc_z[k] * sil1_27[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_36 = buffer.data(sil0 + 36);
    const auto *sil0_48 = buffer.data(sil0 + 48);
    const auto *sil0_51 = buffer.data(sil0 + 51);
    const auto *sil0_55 = buffer.data(sil0 + 55);
    const auto *sil0_60 = buffer.data(sil0 + 60);
    const auto *sil0_66 = buffer.data(sil0 + 66);
    const auto *sil0_81 = buffer.data(sil0 + 81);
    const auto *sil0_90 = buffer.data(sil0 + 90);
    const auto *sil0_95 = buffer.data(sil0 + 95);
    const auto *sil0_99 = buffer.data(sil0 + 99);
    const auto *sil0_102 = buffer.data(sil0 + 102);
    const auto *sil0_104 = buffer.data(sil0 + 104);
    const auto *sil0_107 = buffer.data(sil0 + 107);
    const auto *sil0_108 = buffer.data(sil0 + 108);
    const auto *sil0_110 = buffer.data(sil0 + 110);
    const auto *sil0_113 = buffer.data(sil0 + 113);
    const auto *sil0_114 = buffer.data(sil0 + 114);
    const auto *sil0_115 = buffer.data(sil0 + 115);
    const auto *sil0_117 = buffer.data(sil0 + 117);
    const auto *sil0_134 = buffer.data(sil0 + 134);

    const auto *sik_28 = buffer.data(sik + 28);
    const auto *sik_35 = buffer.data(sik + 35);
    const auto *sik_36 = buffer.data(sik + 36);
    const auto *sik_38 = buffer.data(sik + 38);
    const auto *sik_39 = buffer.data(sik + 39);
    const auto *sik_41 = buffer.data(sik + 41);
    const auto *sik_42 = buffer.data(sik + 42);
    const auto *sik_45 = buffer.data(sik + 45);
    const auto *sik_46 = buffer.data(sik + 46);
    const auto *sik_50 = buffer.data(sik + 50);
    const auto *sik_51 = buffer.data(sik + 51);
    const auto *sik_56 = buffer.data(sik + 56);
    const auto *sik_64 = buffer.data(sik + 64);
    const auto *sik_66 = buffer.data(sik + 66);
    const auto *sik_67 = buffer.data(sik + 67);
    const auto *sik_68 = buffer.data(sik + 68);
    const auto *sik_69 = buffer.data(sik + 69);
    const auto *sik_70 = buffer.data(sik + 70);
    const auto *sik_71 = buffer.data(sik + 71);
    const auto *sik_72 = buffer.data(sik + 72);
    const auto *sik_74 = buffer.data(sik + 74);
    const auto *sik_75 = buffer.data(sik + 75);
    const auto *sik_77 = buffer.data(sik + 77);
    const auto *sik_80 = buffer.data(sik + 80);
    const auto *sik_81 = buffer.data(sik + 81);
    const auto *sik_84 = buffer.data(sik + 84);
    const auto *sik_85 = buffer.data(sik + 85);
    const auto *sik_86 = buffer.data(sik + 86);
    const auto *sik_89 = buffer.data(sik + 89);
    const auto *sik_90 = buffer.data(sik + 90);
    const auto *sik_91 = buffer.data(sik + 91);
    const auto *sik_92 = buffer.data(sik + 92);
    const auto *sik_100 = buffer.data(sik + 100);
    const auto *sik_101 = buffer.data(sik + 101);
    const auto *sik_102 = buffer.data(sik + 102);
    const auto *sik_103 = buffer.data(sik + 103);
    const auto *sik_104 = buffer.data(sik + 104);
    const auto *sik_105 = buffer.data(sik + 105);
    const auto *sik_106 = buffer.data(sik + 106);
    const auto *sik_107 = buffer.data(sik + 107);
    const auto *sik_108 = buffer.data(sik + 108);
    const auto *sik_111 = buffer.data(sik + 111);
    const auto *sik_113 = buffer.data(sik + 113);
    const auto *sik_114 = buffer.data(sik + 114);
    const auto *sik_117 = buffer.data(sik + 117);
    const auto *sik_118 = buffer.data(sik + 118);
    const auto *sik_120 = buffer.data(sik + 120);
    const auto *sik_122 = buffer.data(sik + 122);
    const auto *sik_123 = buffer.data(sik + 123);
    const auto *sik_125 = buffer.data(sik + 125);
    const auto *sik_126 = buffer.data(sik + 126);
    const auto *sik_128 = buffer.data(sik + 128);
    const auto *sik_129 = buffer.data(sik + 129);
    const auto *sik_131 = buffer.data(sik + 131);
    const auto *sik_132 = buffer.data(sik + 132);
    const auto *sik_133 = buffer.data(sik + 133);
    const auto *sik_135 = buffer.data(sik + 135);
    const auto *sik_136 = buffer.data(sik + 136);
    const auto *sik_137 = buffer.data(sik + 137);
    const auto *sik_138 = buffer.data(sik + 138);
    const auto *sik_139 = buffer.data(sik + 139);
    const auto *sik_140 = buffer.data(sik + 140);
    const auto *sik_141 = buffer.data(sik + 141);
    const auto *sik_142 = buffer.data(sik + 142);
    const auto *sik_143 = buffer.data(sik + 143);
    const auto *sik_172 = buffer.data(sik + 172);
    const auto *sik_173 = buffer.data(sik + 173);
    const auto *sik_174 = buffer.data(sik + 174);
    const auto *sik_175 = buffer.data(sik + 175);
    const auto *sik_176 = buffer.data(sik + 176);
    const auto *sik_177 = buffer.data(sik + 177);
    const auto *sik_178 = buffer.data(sik + 178);
    const auto *sik_179 = buffer.data(sik + 179);
    const auto *sik_180 = buffer.data(sik + 180);
    const auto *sik_183 = buffer.data(sik + 183);
    const auto *sik_185 = buffer.data(sik + 185);
    const auto *sik_186 = buffer.data(sik + 186);

    const auto *sil1_36 = buffer.data(sil1 + 36);
    const auto *sil1_48 = buffer.data(sil1 + 48);
    const auto *sil1_51 = buffer.data(sil1 + 51);
    const auto *sil1_55 = buffer.data(sil1 + 55);
    const auto *sil1_60 = buffer.data(sil1 + 60);
    const auto *sil1_66 = buffer.data(sil1 + 66);
    const auto *sil1_81 = buffer.data(sil1 + 81);
    const auto *sil1_90 = buffer.data(sil1 + 90);
    const auto *sil1_95 = buffer.data(sil1 + 95);
    const auto *sil1_99 = buffer.data(sil1 + 99);
    const auto *sil1_102 = buffer.data(sil1 + 102);
    const auto *sil1_104 = buffer.data(sil1 + 104);
    const auto *sil1_107 = buffer.data(sil1 + 107);
    const auto *sil1_108 = buffer.data(sil1 + 108);
    const auto *sil1_110 = buffer.data(sil1 + 110);
    const auto *sil1_113 = buffer.data(sil1 + 113);
    const auto *sil1_114 = buffer.data(sil1 + 114);
    const auto *sil1_115 = buffer.data(sil1 + 115);
    const auto *sil1_117 = buffer.data(sil1 + 117);
    const auto *sil1_134 = buffer.data(sil1 + 134);

    const auto *ski0_79 = buffer.data(ski0 + 79);
    const auto *ski0_80 = buffer.data(ski0 + 80);
    const auto *ski0_81 = buffer.data(ski0 + 81);
    const auto *ski0_82 = buffer.data(ski0 + 82);
    const auto *ski0_83 = buffer.data(ski0 + 83);
    const auto *ski0_84 = buffer.data(ski0 + 84);
    const auto *ski0_87 = buffer.data(ski0 + 87);
    const auto *ski0_89 = buffer.data(ski0 + 89);
    const auto *ski0_90 = buffer.data(ski0 + 90);
    const auto *ski0_93 = buffer.data(ski0 + 93);
    const auto *ski0_94 = buffer.data(ski0 + 94);
    const auto *ski0_96 = buffer.data(ski0 + 96);
    const auto *ski0_98 = buffer.data(ski0 + 98);
    const auto *ski0_99 = buffer.data(ski0 + 99);
    const auto *ski0_101 = buffer.data(ski0 + 101);
    const auto *ski0_102 = buffer.data(ski0 + 102);
    const auto *ski0_104 = buffer.data(ski0 + 104);
    const auto *ski0_105 = buffer.data(ski0 + 105);
    const auto *ski0_107 = buffer.data(ski0 + 107);
    const auto *ski0_108 = buffer.data(ski0 + 108);
    const auto *ski0_109 = buffer.data(ski0 + 109);
    const auto *ski0_110 = buffer.data(ski0 + 110);
    const auto *ski0_111 = buffer.data(ski0 + 111);
    const auto *ski0_135 = buffer.data(ski0 + 135);
    const auto *ski0_136 = buffer.data(ski0 + 136);
    const auto *ski0_137 = buffer.data(ski0 + 137);
    const auto *ski0_138 = buffer.data(ski0 + 138);
    const auto *ski0_139 = buffer.data(ski0 + 139);
    const auto *ski0_140 = buffer.data(ski0 + 140);
    const auto *ski0_143 = buffer.data(ski0 + 143);
    const auto *ski0_145 = buffer.data(ski0 + 145);
    const auto *ski0_146 = buffer.data(ski0 + 146);

    const auto *ski1_79 = buffer.data(ski1 + 79);
    const auto *ski1_80 = buffer.data(ski1 + 80);
    const auto *ski1_81 = buffer.data(ski1 + 81);
    const auto *ski1_82 = buffer.data(ski1 + 82);
    const auto *ski1_83 = buffer.data(ski1 + 83);
    const auto *ski1_84 = buffer.data(ski1 + 84);
    const auto *ski1_87 = buffer.data(ski1 + 87);
    const auto *ski1_89 = buffer.data(ski1 + 89);
    const auto *ski1_90 = buffer.data(ski1 + 90);
    const auto *ski1_93 = buffer.data(ski1 + 93);
    const auto *ski1_94 = buffer.data(ski1 + 94);
    const auto *ski1_96 = buffer.data(ski1 + 96);
    const auto *ski1_98 = buffer.data(ski1 + 98);
    const auto *ski1_99 = buffer.data(ski1 + 99);
    const auto *ski1_101 = buffer.data(ski1 + 101);
    const auto *ski1_102 = buffer.data(ski1 + 102);
    const auto *ski1_104 = buffer.data(ski1 + 104);
    const auto *ski1_105 = buffer.data(ski1 + 105);
    const auto *ski1_107 = buffer.data(ski1 + 107);
    const auto *ski1_108 = buffer.data(ski1 + 108);
    const auto *ski1_109 = buffer.data(ski1 + 109);
    const auto *ski1_110 = buffer.data(ski1 + 110);
    const auto *ski1_111 = buffer.data(ski1 + 111);
    const auto *ski1_135 = buffer.data(ski1 + 135);
    const auto *ski1_136 = buffer.data(ski1 + 136);
    const auto *ski1_137 = buffer.data(ski1 + 137);
    const auto *ski1_138 = buffer.data(ski1 + 138);
    const auto *ski1_139 = buffer.data(ski1 + 139);
    const auto *ski1_140 = buffer.data(ski1 + 140);
    const auto *ski1_143 = buffer.data(ski1 + 143);
    const auto *ski1_145 = buffer.data(ski1 + 145);
    const auto *ski1_146 = buffer.data(ski1 + 146);

    const auto *skk_100 = buffer.data(skk + 100);
    const auto *skk_101 = buffer.data(skk + 101);
    const auto *skk_102 = buffer.data(skk + 102);
    const auto *skk_103 = buffer.data(skk + 103);
    const auto *skk_104 = buffer.data(skk + 104);
    const auto *skk_105 = buffer.data(skk + 105);
    const auto *skk_106 = buffer.data(skk + 106);
    const auto *skk_107 = buffer.data(skk + 107);
    const auto *skk_108 = buffer.data(skk + 108);
    const auto *skk_110 = buffer.data(skk + 110);
    const auto *skk_111 = buffer.data(skk + 111);
    const auto *skk_113 = buffer.data(skk + 113);
    const auto *skk_114 = buffer.data(skk + 114);
    const auto *skk_117 = buffer.data(skk + 117);
    const auto *skk_118 = buffer.data(skk + 118);
    const auto *skk_120 = buffer.data(skk + 120);
    const auto *skk_122 = buffer.data(skk + 122);
    const auto *skk_123 = buffer.data(skk + 123);
    const auto *skk_125 = buffer.data(skk + 125);
    const auto *skk_126 = buffer.data(skk + 126);
    const auto *skk_128 = buffer.data(skk + 128);
    const auto *skk_129 = buffer.data(skk + 129);
    const auto *skk_131 = buffer.data(skk + 131);
    const auto *skk_132 = buffer.data(skk + 132);
    const auto *skk_133 = buffer.data(skk + 133);
    const auto *skk_135 = buffer.data(skk + 135);
    const auto *skk_136 = buffer.data(skk + 136);
    const auto *skk_137 = buffer.data(skk + 137);
    const auto *skk_138 = buffer.data(skk + 138);
    const auto *skk_139 = buffer.data(skk + 139);
    const auto *skk_140 = buffer.data(skk + 140);
    const auto *skk_141 = buffer.data(skk + 141);
    const auto *skk_142 = buffer.data(skk + 142);
    const auto *skk_143 = buffer.data(skk + 143);
    const auto *skk_144 = buffer.data(skk + 144);
    const auto *skk_146 = buffer.data(skk + 146);
    const auto *skk_147 = buffer.data(skk + 147);
    const auto *skk_149 = buffer.data(skk + 149);
    const auto *skk_150 = buffer.data(skk + 150);
    const auto *skk_153 = buffer.data(skk + 153);
    const auto *skk_154 = buffer.data(skk + 154);
    const auto *skk_158 = buffer.data(skk + 158);
    const auto *skk_159 = buffer.data(skk + 159);
    const auto *skk_164 = buffer.data(skk + 164);
    const auto *skk_172 = buffer.data(skk + 172);
    const auto *skk_173 = buffer.data(skk + 173);
    const auto *skk_174 = buffer.data(skk + 174);
    const auto *skk_175 = buffer.data(skk + 175);
    const auto *skk_176 = buffer.data(skk + 176);
    const auto *skk_177 = buffer.data(skk + 177);
    const auto *skk_178 = buffer.data(skk + 178);
    const auto *skk_179 = buffer.data(skk + 179);
    const auto *skk_180 = buffer.data(skk + 180);
    const auto *skk_182 = buffer.data(skk + 182);
    const auto *skk_183 = buffer.data(skk + 183);
    const auto *skk_185 = buffer.data(skk + 185);
    const auto *skk_186 = buffer.data(skk + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, sik_100, sik_101, sik_102, \
                         sik_103, sik_104, skk_100, skk_101, skk_102, skk_103, \
                         skk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_20 * sik_100[k]
                   + f_3 * pc_x[k] * skk_100[k];

        t_119[k] = f_20 * sik_101[k]
                   + f_3 * pc_x[k] * skk_101[k];

        t_120[k] = f_20 * sik_102[k]
                   + f_3 * pc_x[k] * skk_102[k];

        t_121[k] = f_20 * sik_103[k]
                   + f_3 * pc_x[k] * skk_103[k];

        t_122[k] = f_20 * sik_104[k]
                   + f_3 * pc_x[k] * skk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, sil0_36, sik_105, \
                         sik_106, sik_107, sil1_36, skk_105, skk_106, \
                         skk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_20 * sik_105[k]
                   + f_3 * pc_x[k] * skk_105[k];

        t_124[k] = f_20 * sik_106[k]
                   + f_3 * pc_x[k] * skk_106[k];

        t_125[k] = f_20 * sik_107[k]
                   + f_3 * pc_x[k] * skk_107[k];

        t_126[k] = pb_z[k] * sil0_36[k]
                   - f_14 * pc_z[k] * sil1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, sik_28, ski0_79, ski0_80, ski1_79, \
                         ski1_80, skk_100, skk_102, skk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * sik_28[k]
                   + f_3 * pc_z[k] * skk_100[k];

        t_128[k] = f_4 * ski0_79[k]
                   - f_5 * ski1_79[k]
                   + f_3 * pc_y[k] * skk_102[k];

        t_129[k] = f_6 * ski0_80[k]
                   - f_7 * ski1_80[k]
                   + f_3 * pc_y[k] * skk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, ski0_81, ski0_82, ski0_83, ski1_81, \
                         ski1_82, ski1_83, skk_104, skk_105, skk_106, \
                         skk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * ski0_81[k]
                   - f_9 * ski1_81[k]
                   + f_3 * pc_y[k] * skk_104[k];

        t_131[k] = f_10 * ski0_82[k]
                   - f_11 * ski1_82[k]
                   + f_3 * pc_y[k] * skk_105[k];

        t_132[k] = f_12 * ski0_83[k]
                   - f_13 * ski1_83[k]
                   + f_3 * pc_y[k] * skk_106[k];

        t_133[k] = f_3 * pc_y[k] * skk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, sik_35, sik_36, \
                         sik_108, ski0_83, ski0_84, ski1_83, ski1_84, skk_107, \
                         skk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * sik_35[k]
                   + f_1 * ski0_83[k]
                   - f_2 * ski1_83[k]
                   + f_3 * pc_z[k] * skk_107[k];

        t_135[k] = f_19 * sik_108[k]
                   + f_1 * ski0_84[k]
                   - f_2 * ski1_84[k]
                   + f_3 * pc_x[k] * skk_108[k];

        t_136[k] = f_16 * sik_36[k]
                   + f_3 * pc_y[k] * skk_108[k];

        t_137[k] = f_3 * pc_z[k] * skk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, sik_38, sik_111, sik_113, ski0_87, \
                         ski0_89, ski1_87, ski1_89, skk_110, skk_111, \
                         skk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_19 * sik_111[k]
                   + f_4 * ski0_87[k]
                   - f_5 * ski1_87[k]
                   + f_3 * pc_x[k] * skk_111[k];

        t_139[k] = f_16 * sik_38[k]
                   + f_3 * pc_y[k] * skk_110[k];

        t_140[k] = f_19 * sik_113[k]
                   + f_4 * ski0_89[k]
                   - f_5 * ski1_89[k]
                   + f_3 * pc_x[k] * skk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, sik_41, sik_114, ski0_90, \
                         ski1_90, skk_111, skk_113, skk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_19 * sik_114[k]
                   + f_6 * ski0_90[k]
                   - f_7 * ski1_90[k]
                   + f_3 * pc_x[k] * skk_114[k];

        t_142[k] = f_3 * pc_z[k] * skk_111[k];

        t_143[k] = f_16 * sik_41[k]
                   + f_3 * pc_y[k] * skk_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, sik_117, sik_118, ski0_93, ski0_94, \
                         ski1_93, ski1_94, skk_114, skk_117, skk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_19 * sik_117[k]
                   + f_6 * ski0_93[k]
                   - f_7 * ski1_93[k]
                   + f_3 * pc_x[k] * skk_117[k];

        t_145[k] = f_19 * sik_118[k]
                   + f_8 * ski0_94[k]
                   - f_9 * ski1_94[k]
                   + f_3 * pc_x[k] * skk_118[k];

        t_146[k] = f_3 * pc_z[k] * skk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, sik_45, sik_120, sik_122, ski0_96, \
                         ski0_98, ski1_96, ski1_98, skk_117, skk_120, \
                         skk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_19 * sik_120[k]
                   + f_8 * ski0_96[k]
                   - f_9 * ski1_96[k]
                   + f_3 * pc_x[k] * skk_120[k];

        t_148[k] = f_16 * sik_45[k]
                   + f_3 * pc_y[k] * skk_117[k];

        t_149[k] = f_19 * sik_122[k]
                   + f_8 * ski0_98[k]
                   - f_9 * ski1_98[k]
                   + f_3 * pc_x[k] * skk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, sik_123, sik_125, ski0_99, ski0_101, \
                         ski1_99, ski1_101, skk_118, skk_123, skk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_19 * sik_123[k]
                   + f_10 * ski0_99[k]
                   - f_11 * ski1_99[k]
                   + f_3 * pc_x[k] * skk_123[k];

        t_151[k] = f_3 * pc_z[k] * skk_118[k];

        t_152[k] = f_19 * sik_125[k]
                   + f_10 * ski0_101[k]
                   - f_11 * ski1_101[k]
                   + f_3 * pc_x[k] * skk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, sik_50, sik_126, sik_128, ski0_102, \
                         ski0_104, ski1_102, ski1_104, skk_122, skk_126, \
                         skk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_19 * sik_126[k]
                   + f_10 * ski0_102[k]
                   - f_11 * ski1_102[k]
                   + f_3 * pc_x[k] * skk_126[k];

        t_154[k] = f_16 * sik_50[k]
                   + f_3 * pc_y[k] * skk_122[k];

        t_155[k] = f_19 * sik_128[k]
                   + f_10 * ski0_104[k]
                   - f_11 * ski1_104[k]
                   + f_3 * pc_x[k] * skk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, sik_129, sik_131, ski0_105, \
                         ski0_107, ski1_105, ski1_107, skk_123, skk_129, \
                         skk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_19 * sik_129[k]
                   + f_12 * ski0_105[k]
                   - f_13 * ski1_105[k]
                   + f_3 * pc_x[k] * skk_129[k];

        t_157[k] = f_3 * pc_z[k] * skk_123[k];

        t_158[k] = f_19 * sik_131[k]
                   + f_12 * ski0_107[k]
                   - f_13 * ski1_107[k]
                   + f_3 * pc_x[k] * skk_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, sik_56, sik_132, sik_133, ski0_108, \
                         ski0_109, ski1_108, ski1_109, skk_128, skk_132, \
                         skk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_19 * sik_132[k]
                   + f_12 * ski0_108[k]
                   - f_13 * ski1_108[k]
                   + f_3 * pc_x[k] * skk_132[k];

        t_160[k] = f_19 * sik_133[k]
                   + f_12 * ski0_109[k]
                   - f_13 * ski1_109[k]
                   + f_3 * pc_x[k] * skk_133[k];

        t_161[k] = f_16 * sik_56[k]
                   + f_3 * pc_y[k] * skk_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, sik_135, sik_136, sik_137, sik_138, \
                         ski0_111, ski1_111, skk_135, skk_136, skk_137, \
                         skk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_19 * sik_135[k]
                   + f_12 * ski0_111[k]
                   - f_13 * ski1_111[k]
                   + f_3 * pc_x[k] * skk_135[k];

        t_163[k] = f_19 * sik_136[k]
                   + f_3 * pc_x[k] * skk_136[k];

        t_164[k] = f_19 * sik_137[k]
                   + f_3 * pc_x[k] * skk_137[k];

        t_165[k] = f_19 * sik_138[k]
                   + f_3 * pc_x[k] * skk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, sik_139, sik_140, sik_141, \
                         sik_142, sik_143, skk_139, skk_140, skk_141, skk_142, \
                         skk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_19 * sik_139[k]
                   + f_3 * pc_x[k] * skk_139[k];

        t_167[k] = f_19 * sik_140[k]
                   + f_3 * pc_x[k] * skk_140[k];

        t_168[k] = f_19 * sik_141[k]
                   + f_3 * pc_x[k] * skk_141[k];

        t_169[k] = f_19 * sik_142[k]
                   + f_3 * pc_x[k] * skk_142[k];

        t_170[k] = f_19 * sik_143[k]
                   + f_3 * pc_x[k] * skk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, sik_64, sik_66, ski0_105, ski0_107, \
                         ski1_105, ski1_107, skk_136, skk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * sik_64[k]
                   + f_1 * ski0_105[k]
                   - f_2 * ski1_105[k]
                   + f_3 * pc_y[k] * skk_136[k];

        t_172[k] = f_3 * pc_z[k] * skk_136[k];

        t_173[k] = f_16 * sik_66[k]
                   + f_4 * ski0_107[k]
                   - f_5 * ski1_107[k]
                   + f_3 * pc_y[k] * skk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, sik_67, sik_68, sik_69, ski0_108, \
                         ski0_109, ski0_110, ski1_108, ski1_109, ski1_110, skk_139, skk_140, \
                         skk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * sik_67[k]
                   + f_6 * ski0_108[k]
                   - f_7 * ski1_108[k]
                   + f_3 * pc_y[k] * skk_139[k];

        t_175[k] = f_16 * sik_68[k]
                   + f_8 * ski0_109[k]
                   - f_9 * ski1_109[k]
                   + f_3 * pc_y[k] * skk_140[k];

        t_176[k] = f_16 * sik_69[k]
                   + f_10 * ski0_110[k]
                   - f_11 * ski1_110[k]
                   + f_3 * pc_y[k] * skk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, sil0_90, sik_70, \
                         sik_71, sil1_90, ski0_111, ski1_111, skk_142, \
                         skk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * sik_70[k]
                   + f_12 * ski0_111[k]
                   - f_13 * ski1_111[k]
                   + f_3 * pc_y[k] * skk_142[k];

        t_178[k] = f_16 * sik_71[k]
                   + f_3 * pc_y[k] * skk_143[k];

        t_179[k] = f_1 * ski0_111[k]
                   - f_2 * ski1_111[k]
                   + f_3 * pc_z[k] * skk_143[k];

        t_180[k] = pb_y[k] * sil0_90[k]
                   - f_14 * pc_y[k] * sil1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, sil0_48, sik_36, \
                         sik_72, sik_74, sil1_48, skk_144, skk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * sik_72[k]
                   + f_3 * pc_y[k] * skk_144[k];

        t_182[k] = f_15 * sik_36[k]
                   + f_3 * pc_z[k] * skk_144[k];

        t_183[k] = pb_z[k] * sil0_48[k]
                   - f_14 * pc_z[k] * sil1_48[k];

        t_184[k] = f_15 * sik_74[k]
                   + f_3 * pc_y[k] * skk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, sil0_51, sil0_95, \
                         sik_39, sik_77, sil1_51, sil1_95, skk_147, \
                         skk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * sil0_95[k]
                   - f_14 * pc_y[k] * sil1_95[k];

        t_186[k] = pb_z[k] * sil0_51[k]
                   - f_14 * pc_z[k] * sil1_51[k];

        t_187[k] = f_15 * sik_39[k]
                   + f_3 * pc_z[k] * skk_147[k];

        t_188[k] = f_15 * sik_77[k]
                   + f_3 * pc_y[k] * skk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, sil0_55, sil0_99, \
                         sik_42, sil1_55, sil1_99, skk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * sil0_99[k]
                   - f_14 * pc_y[k] * sil1_99[k];

        t_190[k] = pb_z[k] * sil0_55[k]
                   - f_14 * pc_z[k] * sil1_55[k];

        t_191[k] = f_15 * sik_42[k]
                   + f_3 * pc_z[k] * skk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, sil0_102, sil0_104, sik_80, sik_81, \
                         sil1_102, sil1_104, skk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * sil0_102[k]
                   + f_16 * sik_80[k]
                   - f_14 * pc_y[k] * sil1_102[k];

        t_193[k] = f_15 * sik_81[k]
                   + f_3 * pc_y[k] * skk_153[k];

        t_194[k] = pb_y[k] * sil0_104[k]
                   - f_14 * pc_y[k] * sil1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, sil0_60, sil0_107, \
                         sik_46, sik_84, sil1_60, sil1_107, skk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * sil0_60[k]
                   - f_14 * pc_z[k] * sil1_60[k];

        t_196[k] = f_15 * sik_46[k]
                   + f_3 * pc_z[k] * skk_154[k];

        t_197[k] = pb_y[k] * sil0_107[k]
                   + f_17 * sik_84[k]
                   - f_14 * pc_y[k] * sil1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, sil0_108, sil0_110, sik_85, sik_86, \
                         sil1_108, sil1_110, skk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * sil0_108[k]
                   + f_16 * sik_85[k]
                   - f_14 * pc_y[k] * sil1_108[k];

        t_199[k] = f_15 * sik_86[k]
                   + f_3 * pc_y[k] * skk_158[k];

        t_200[k] = pb_y[k] * sil0_110[k]
                   - f_14 * pc_y[k] * sil1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, sil0_66, sil0_113, \
                         sik_51, sik_89, sil1_66, sil1_113, skk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * sil0_66[k]
                   - f_14 * pc_z[k] * sil1_66[k];

        t_202[k] = f_15 * sik_51[k]
                   + f_3 * pc_z[k] * skk_159[k];

        t_203[k] = pb_y[k] * sil0_113[k]
                   + f_18 * sik_89[k]
                   - f_14 * pc_y[k] * sil1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, sil0_114, sil0_115, sil0_117, \
                         sik_90, sik_91, sik_92, sil1_114, sil1_115, sil1_117, \
                         skk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * sil0_114[k]
                   + f_17 * sik_90[k]
                   - f_14 * pc_y[k] * sil1_114[k];

        t_205[k] = pb_y[k] * sil0_115[k]
                   + f_16 * sik_91[k]
                   - f_14 * pc_y[k] * sil1_115[k];

        t_206[k] = f_15 * sik_92[k]
                   + f_3 * pc_y[k] * skk_164[k];

        t_207[k] = pb_y[k] * sil0_117[k]
                   - f_14 * pc_y[k] * sil1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, sik_172, sik_173, sik_174, \
                         sik_175, sik_176, skk_172, skk_173, skk_174, skk_175, \
                         skk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_19 * sik_172[k]
                   + f_3 * pc_x[k] * skk_172[k];

        t_209[k] = f_19 * sik_173[k]
                   + f_3 * pc_x[k] * skk_173[k];

        t_210[k] = f_19 * sik_174[k]
                   + f_3 * pc_x[k] * skk_174[k];

        t_211[k] = f_19 * sik_175[k]
                   + f_3 * pc_x[k] * skk_175[k];

        t_212[k] = f_19 * sik_176[k]
                   + f_3 * pc_x[k] * skk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, sil0_81, sik_177, \
                         sik_178, sik_179, sil1_81, skk_177, skk_178, \
                         skk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_19 * sik_177[k]
                   + f_3 * pc_x[k] * skk_177[k];

        t_214[k] = f_19 * sik_178[k]
                   + f_3 * pc_x[k] * skk_178[k];

        t_215[k] = f_19 * sik_179[k]
                   + f_3 * pc_x[k] * skk_179[k];

        t_216[k] = pb_z[k] * sil0_81[k]
                   - f_14 * pc_z[k] * sil1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, sik_64, sik_102, sik_103, ski0_135, \
                         ski0_136, ski1_135, ski1_136, skk_172, skk_174, \
                         skk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * sik_64[k]
                   + f_3 * pc_z[k] * skk_172[k];

        t_218[k] = f_15 * sik_102[k]
                   + f_4 * ski0_135[k]
                   - f_5 * ski1_135[k]
                   + f_3 * pc_y[k] * skk_174[k];

        t_219[k] = f_15 * sik_103[k]
                   + f_6 * ski0_136[k]
                   - f_7 * ski1_136[k]
                   + f_3 * pc_y[k] * skk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, sik_104, sik_105, sik_106, ski0_137, \
                         ski0_138, ski0_139, ski1_137, ski1_138, ski1_139, skk_176, skk_177, \
                         skk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * sik_104[k]
                   + f_8 * ski0_137[k]
                   - f_9 * ski1_137[k]
                   + f_3 * pc_y[k] * skk_176[k];

        t_221[k] = f_15 * sik_105[k]
                   + f_10 * ski0_138[k]
                   - f_11 * ski1_138[k]
                   + f_3 * pc_y[k] * skk_177[k];

        t_222[k] = f_15 * sik_106[k]
                   + f_12 * ski0_139[k]
                   - f_13 * ski1_139[k]
                   + f_3 * pc_y[k] * skk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, sil0_134, sik_107, \
                         sik_180, sil1_134, ski0_140, ski1_140, skk_179, \
                         skk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * sik_107[k]
                   + f_3 * pc_y[k] * skk_179[k];

        t_224[k] = pb_y[k] * sil0_134[k]
                   - f_14 * pc_y[k] * sil1_134[k];

        t_225[k] = f_19 * sik_180[k]
                   + f_1 * ski0_140[k]
                   - f_2 * ski1_140[k]
                   + f_3 * pc_x[k] * skk_180[k];

        t_226[k] = f_3 * pc_y[k] * skk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, sik_72, sik_183, ski0_143, \
                         ski1_143, skk_180, skk_182, skk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * sik_72[k]
                   + f_3 * pc_z[k] * skk_180[k];

        t_228[k] = f_19 * sik_183[k]
                   + f_4 * ski0_143[k]
                   - f_5 * ski1_143[k]
                   + f_3 * pc_x[k] * skk_183[k];

        t_229[k] = f_3 * pc_y[k] * skk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, sik_75, sik_185, sik_186, ski0_145, \
                         ski0_146, ski1_145, ski1_146, skk_183, skk_185, \
                         skk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_19 * sik_185[k]
                   + f_4 * ski0_145[k]
                   - f_5 * ski1_145[k]
                   + f_3 * pc_x[k] * skk_185[k];

        t_231[k] = f_19 * sik_186[k]
                   + f_6 * ski0_146[k]
                   - f_7 * ski1_146[k]
                   + f_3 * pc_x[k] * skk_186[k];

        t_232[k] = f_16 * sik_75[k]
                   + f_3 * pc_z[k] * skk_183[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_135 = buffer.data(sil0 + 135);
    const auto *sil0_138 = buffer.data(sil0 + 138);
    const auto *sil0_141 = buffer.data(sil0 + 141);
    const auto *sil0_145 = buffer.data(sil0 + 145);
    const auto *sil0_147 = buffer.data(sil0 + 147);
    const auto *sil0_150 = buffer.data(sil0 + 150);
    const auto *sil0_152 = buffer.data(sil0 + 152);
    const auto *sil0_153 = buffer.data(sil0 + 153);
    const auto *sil0_156 = buffer.data(sil0 + 156);
    const auto *sil0_158 = buffer.data(sil0 + 158);
    const auto *sil0_159 = buffer.data(sil0 + 159);
    const auto *sil0_160 = buffer.data(sil0 + 160);

    const auto *sik_78 = buffer.data(sik + 78);
    const auto *sik_82 = buffer.data(sik + 82);
    const auto *sik_87 = buffer.data(sik + 87);
    const auto *sik_100 = buffer.data(sik + 100);
    const auto *sik_107 = buffer.data(sik + 107);
    const auto *sik_108 = buffer.data(sik + 108);
    const auto *sik_110 = buffer.data(sik + 110);
    const auto *sik_111 = buffer.data(sik + 111);
    const auto *sik_113 = buffer.data(sik + 113);
    const auto *sik_114 = buffer.data(sik + 114);
    const auto *sik_115 = buffer.data(sik + 115);
    const auto *sik_117 = buffer.data(sik + 117);
    const auto *sik_118 = buffer.data(sik + 118);
    const auto *sik_119 = buffer.data(sik + 119);
    const auto *sik_120 = buffer.data(sik + 120);
    const auto *sik_122 = buffer.data(sik + 122);
    const auto *sik_123 = buffer.data(sik + 123);
    const auto *sik_124 = buffer.data(sik + 124);
    const auto *sik_125 = buffer.data(sik + 125);
    const auto *sik_126 = buffer.data(sik + 126);
    const auto *sik_128 = buffer.data(sik + 128);
    const auto *sik_136 = buffer.data(sik + 136);
    const auto *sik_138 = buffer.data(sik + 138);
    const auto *sik_139 = buffer.data(sik + 139);
    const auto *sik_140 = buffer.data(sik + 140);
    const auto *sik_141 = buffer.data(sik + 141);
    const auto *sik_142 = buffer.data(sik + 142);
    const auto *sik_143 = buffer.data(sik + 143);
    const auto *sik_144 = buffer.data(sik + 144);
    const auto *sik_146 = buffer.data(sik + 146);
    const auto *sik_149 = buffer.data(sik + 149);
    const auto *sik_153 = buffer.data(sik + 153);
    const auto *sik_158 = buffer.data(sik + 158);
    const auto *sik_189 = buffer.data(sik + 189);
    const auto *sik_190 = buffer.data(sik + 190);
    const auto *sik_192 = buffer.data(sik + 192);
    const auto *sik_194 = buffer.data(sik + 194);
    const auto *sik_195 = buffer.data(sik + 195);
    const auto *sik_197 = buffer.data(sik + 197);
    const auto *sik_198 = buffer.data(sik + 198);
    const auto *sik_200 = buffer.data(sik + 200);
    const auto *sik_201 = buffer.data(sik + 201);
    const auto *sik_203 = buffer.data(sik + 203);
    const auto *sik_204 = buffer.data(sik + 204);
    const auto *sik_205 = buffer.data(sik + 205);
    const auto *sik_207 = buffer.data(sik + 207);
    const auto *sik_208 = buffer.data(sik + 208);
    const auto *sik_209 = buffer.data(sik + 209);
    const auto *sik_210 = buffer.data(sik + 210);
    const auto *sik_211 = buffer.data(sik + 211);
    const auto *sik_212 = buffer.data(sik + 212);
    const auto *sik_213 = buffer.data(sik + 213);
    const auto *sik_214 = buffer.data(sik + 214);
    const auto *sik_215 = buffer.data(sik + 215);
    const auto *sik_216 = buffer.data(sik + 216);
    const auto *sik_219 = buffer.data(sik + 219);
    const auto *sik_221 = buffer.data(sik + 221);
    const auto *sik_222 = buffer.data(sik + 222);
    const auto *sik_225 = buffer.data(sik + 225);
    const auto *sik_226 = buffer.data(sik + 226);
    const auto *sik_228 = buffer.data(sik + 228);
    const auto *sik_230 = buffer.data(sik + 230);
    const auto *sik_231 = buffer.data(sik + 231);
    const auto *sik_233 = buffer.data(sik + 233);
    const auto *sik_234 = buffer.data(sik + 234);
    const auto *sik_236 = buffer.data(sik + 236);
    const auto *sik_237 = buffer.data(sik + 237);
    const auto *sik_239 = buffer.data(sik + 239);
    const auto *sik_240 = buffer.data(sik + 240);
    const auto *sik_241 = buffer.data(sik + 241);
    const auto *sik_243 = buffer.data(sik + 243);
    const auto *sik_244 = buffer.data(sik + 244);
    const auto *sik_245 = buffer.data(sik + 245);
    const auto *sik_246 = buffer.data(sik + 246);
    const auto *sik_247 = buffer.data(sik + 247);
    const auto *sik_248 = buffer.data(sik + 248);
    const auto *sik_249 = buffer.data(sik + 249);
    const auto *sik_250 = buffer.data(sik + 250);
    const auto *sik_251 = buffer.data(sik + 251);
    const auto *sik_257 = buffer.data(sik + 257);
    const auto *sik_261 = buffer.data(sik + 261);
    const auto *sik_266 = buffer.data(sik + 266);
    const auto *sik_272 = buffer.data(sik + 272);

    const auto *sil1_135 = buffer.data(sil1 + 135);
    const auto *sil1_138 = buffer.data(sil1 + 138);
    const auto *sil1_141 = buffer.data(sil1 + 141);
    const auto *sil1_145 = buffer.data(sil1 + 145);
    const auto *sil1_147 = buffer.data(sil1 + 147);
    const auto *sil1_150 = buffer.data(sil1 + 150);
    const auto *sil1_152 = buffer.data(sil1 + 152);
    const auto *sil1_153 = buffer.data(sil1 + 153);
    const auto *sil1_156 = buffer.data(sil1 + 156);
    const auto *sil1_158 = buffer.data(sil1 + 158);
    const auto *sil1_159 = buffer.data(sil1 + 159);
    const auto *sil1_160 = buffer.data(sil1 + 160);

    const auto *ski0_149 = buffer.data(ski0 + 149);
    const auto *ski0_150 = buffer.data(ski0 + 150);
    const auto *ski0_152 = buffer.data(ski0 + 152);
    const auto *ski0_154 = buffer.data(ski0 + 154);
    const auto *ski0_155 = buffer.data(ski0 + 155);
    const auto *ski0_157 = buffer.data(ski0 + 157);
    const auto *ski0_158 = buffer.data(ski0 + 158);
    const auto *ski0_160 = buffer.data(ski0 + 160);
    const auto *ski0_161 = buffer.data(ski0 + 161);
    const auto *ski0_163 = buffer.data(ski0 + 163);
    const auto *ski0_164 = buffer.data(ski0 + 164);
    const auto *ski0_165 = buffer.data(ski0 + 165);
    const auto *ski0_166 = buffer.data(ski0 + 166);
    const auto *ski0_167 = buffer.data(ski0 + 167);
    const auto *ski0_168 = buffer.data(ski0 + 168);
    const auto *ski0_171 = buffer.data(ski0 + 171);
    const auto *ski0_173 = buffer.data(ski0 + 173);
    const auto *ski0_174 = buffer.data(ski0 + 174);
    const auto *ski0_177 = buffer.data(ski0 + 177);
    const auto *ski0_178 = buffer.data(ski0 + 178);
    const auto *ski0_180 = buffer.data(ski0 + 180);
    const auto *ski0_182 = buffer.data(ski0 + 182);
    const auto *ski0_183 = buffer.data(ski0 + 183);
    const auto *ski0_185 = buffer.data(ski0 + 185);
    const auto *ski0_186 = buffer.data(ski0 + 186);
    const auto *ski0_188 = buffer.data(ski0 + 188);
    const auto *ski0_189 = buffer.data(ski0 + 189);
    const auto *ski0_191 = buffer.data(ski0 + 191);
    const auto *ski0_192 = buffer.data(ski0 + 192);
    const auto *ski0_193 = buffer.data(ski0 + 193);
    const auto *ski0_194 = buffer.data(ski0 + 194);
    const auto *ski0_195 = buffer.data(ski0 + 195);
    const auto *ski0_201 = buffer.data(ski0 + 201);
    const auto *ski0_205 = buffer.data(ski0 + 205);
    const auto *ski0_210 = buffer.data(ski0 + 210);
    const auto *ski0_216 = buffer.data(ski0 + 216);

    const auto *ski1_149 = buffer.data(ski1 + 149);
    const auto *ski1_150 = buffer.data(ski1 + 150);
    const auto *ski1_152 = buffer.data(ski1 + 152);
    const auto *ski1_154 = buffer.data(ski1 + 154);
    const auto *ski1_155 = buffer.data(ski1 + 155);
    const auto *ski1_157 = buffer.data(ski1 + 157);
    const auto *ski1_158 = buffer.data(ski1 + 158);
    const auto *ski1_160 = buffer.data(ski1 + 160);
    const auto *ski1_161 = buffer.data(ski1 + 161);
    const auto *ski1_163 = buffer.data(ski1 + 163);
    const auto *ski1_164 = buffer.data(ski1 + 164);
    const auto *ski1_165 = buffer.data(ski1 + 165);
    const auto *ski1_166 = buffer.data(ski1 + 166);
    const auto *ski1_167 = buffer.data(ski1 + 167);
    const auto *ski1_168 = buffer.data(ski1 + 168);
    const auto *ski1_171 = buffer.data(ski1 + 171);
    const auto *ski1_173 = buffer.data(ski1 + 173);
    const auto *ski1_174 = buffer.data(ski1 + 174);
    const auto *ski1_177 = buffer.data(ski1 + 177);
    const auto *ski1_178 = buffer.data(ski1 + 178);
    const auto *ski1_180 = buffer.data(ski1 + 180);
    const auto *ski1_182 = buffer.data(ski1 + 182);
    const auto *ski1_183 = buffer.data(ski1 + 183);
    const auto *ski1_185 = buffer.data(ski1 + 185);
    const auto *ski1_186 = buffer.data(ski1 + 186);
    const auto *ski1_188 = buffer.data(ski1 + 188);
    const auto *ski1_189 = buffer.data(ski1 + 189);
    const auto *ski1_191 = buffer.data(ski1 + 191);
    const auto *ski1_192 = buffer.data(ski1 + 192);
    const auto *ski1_193 = buffer.data(ski1 + 193);
    const auto *ski1_194 = buffer.data(ski1 + 194);
    const auto *ski1_195 = buffer.data(ski1 + 195);
    const auto *ski1_201 = buffer.data(ski1 + 201);
    const auto *ski1_205 = buffer.data(ski1 + 205);
    const auto *ski1_210 = buffer.data(ski1 + 210);
    const auto *ski1_216 = buffer.data(ski1 + 216);

    const auto *skk_185 = buffer.data(skk + 185);
    const auto *skk_186 = buffer.data(skk + 186);
    const auto *skk_189 = buffer.data(skk + 189);
    const auto *skk_190 = buffer.data(skk + 190);
    const auto *skk_192 = buffer.data(skk + 192);
    const auto *skk_194 = buffer.data(skk + 194);
    const auto *skk_195 = buffer.data(skk + 195);
    const auto *skk_197 = buffer.data(skk + 197);
    const auto *skk_198 = buffer.data(skk + 198);
    const auto *skk_200 = buffer.data(skk + 200);
    const auto *skk_201 = buffer.data(skk + 201);
    const auto *skk_203 = buffer.data(skk + 203);
    const auto *skk_204 = buffer.data(skk + 204);
    const auto *skk_205 = buffer.data(skk + 205);
    const auto *skk_207 = buffer.data(skk + 207);
    const auto *skk_208 = buffer.data(skk + 208);
    const auto *skk_209 = buffer.data(skk + 209);
    const auto *skk_210 = buffer.data(skk + 210);
    const auto *skk_211 = buffer.data(skk + 211);
    const auto *skk_212 = buffer.data(skk + 212);
    const auto *skk_213 = buffer.data(skk + 213);
    const auto *skk_214 = buffer.data(skk + 214);
    const auto *skk_215 = buffer.data(skk + 215);
    const auto *skk_216 = buffer.data(skk + 216);
    const auto *skk_218 = buffer.data(skk + 218);
    const auto *skk_219 = buffer.data(skk + 219);
    const auto *skk_221 = buffer.data(skk + 221);
    const auto *skk_222 = buffer.data(skk + 222);
    const auto *skk_225 = buffer.data(skk + 225);
    const auto *skk_226 = buffer.data(skk + 226);
    const auto *skk_228 = buffer.data(skk + 228);
    const auto *skk_230 = buffer.data(skk + 230);
    const auto *skk_231 = buffer.data(skk + 231);
    const auto *skk_233 = buffer.data(skk + 233);
    const auto *skk_234 = buffer.data(skk + 234);
    const auto *skk_236 = buffer.data(skk + 236);
    const auto *skk_237 = buffer.data(skk + 237);
    const auto *skk_239 = buffer.data(skk + 239);
    const auto *skk_240 = buffer.data(skk + 240);
    const auto *skk_241 = buffer.data(skk + 241);
    const auto *skk_243 = buffer.data(skk + 243);
    const auto *skk_244 = buffer.data(skk + 244);
    const auto *skk_245 = buffer.data(skk + 245);
    const auto *skk_246 = buffer.data(skk + 246);
    const auto *skk_247 = buffer.data(skk + 247);
    const auto *skk_248 = buffer.data(skk + 248);
    const auto *skk_249 = buffer.data(skk + 249);
    const auto *skk_250 = buffer.data(skk + 250);
    const auto *skk_251 = buffer.data(skk + 251);
    const auto *skk_252 = buffer.data(skk + 252);
    const auto *skk_254 = buffer.data(skk + 254);
    const auto *skk_255 = buffer.data(skk + 255);
    const auto *skk_257 = buffer.data(skk + 257);
    const auto *skk_258 = buffer.data(skk + 258);
    const auto *skk_261 = buffer.data(skk + 261);
    const auto *skk_262 = buffer.data(skk + 262);
    const auto *skk_266 = buffer.data(skk + 266);
    const auto *skk_267 = buffer.data(skk + 267);
    const auto *skk_272 = buffer.data(skk + 272);

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, sik_189, sik_190, ski0_149, \
                         ski0_150, ski1_149, ski1_150, skk_185, skk_189, \
                         skk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * skk_185[k];

        t_234[k] = f_19 * sik_189[k]
                   + f_6 * ski0_149[k]
                   - f_7 * ski1_149[k]
                   + f_3 * pc_x[k] * skk_189[k];

        t_235[k] = f_19 * sik_190[k]
                   + f_8 * ski0_150[k]
                   - f_9 * ski1_150[k]
                   + f_3 * pc_x[k] * skk_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, sik_78, sik_192, ski0_152, \
                         ski1_152, skk_186, skk_189, skk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * sik_78[k]
                   + f_3 * pc_z[k] * skk_186[k];

        t_237[k] = f_19 * sik_192[k]
                   + f_8 * ski0_152[k]
                   - f_9 * ski1_152[k]
                   + f_3 * pc_x[k] * skk_192[k];

        t_238[k] = f_3 * pc_y[k] * skk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, sik_82, sik_194, sik_195, ski0_154, \
                         ski0_155, ski1_154, ski1_155, skk_190, skk_194, \
                         skk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_19 * sik_194[k]
                   + f_8 * ski0_154[k]
                   - f_9 * ski1_154[k]
                   + f_3 * pc_x[k] * skk_194[k];

        t_240[k] = f_19 * sik_195[k]
                   + f_10 * ski0_155[k]
                   - f_11 * ski1_155[k]
                   + f_3 * pc_x[k] * skk_195[k];

        t_241[k] = f_16 * sik_82[k]
                   + f_3 * pc_z[k] * skk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, sik_197, sik_198, ski0_157, \
                         ski0_158, ski1_157, ski1_158, skk_194, skk_197, \
                         skk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_19 * sik_197[k]
                   + f_10 * ski0_157[k]
                   - f_11 * ski1_157[k]
                   + f_3 * pc_x[k] * skk_197[k];

        t_243[k] = f_19 * sik_198[k]
                   + f_10 * ski0_158[k]
                   - f_11 * ski1_158[k]
                   + f_3 * pc_x[k] * skk_198[k];

        t_244[k] = f_3 * pc_y[k] * skk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, sik_87, sik_200, sik_201, ski0_160, \
                         ski0_161, ski1_160, ski1_161, skk_195, skk_200, \
                         skk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_19 * sik_200[k]
                   + f_10 * ski0_160[k]
                   - f_11 * ski1_160[k]
                   + f_3 * pc_x[k] * skk_200[k];

        t_246[k] = f_19 * sik_201[k]
                   + f_12 * ski0_161[k]
                   - f_13 * ski1_161[k]
                   + f_3 * pc_x[k] * skk_201[k];

        t_247[k] = f_16 * sik_87[k]
                   + f_3 * pc_z[k] * skk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, sik_203, sik_204, sik_205, ski0_163, \
                         ski0_164, ski0_165, ski1_163, ski1_164, ski1_165, skk_203, skk_204, \
                         skk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_19 * sik_203[k]
                   + f_12 * ski0_163[k]
                   - f_13 * ski1_163[k]
                   + f_3 * pc_x[k] * skk_203[k];

        t_249[k] = f_19 * sik_204[k]
                   + f_12 * ski0_164[k]
                   - f_13 * ski1_164[k]
                   + f_3 * pc_x[k] * skk_204[k];

        t_250[k] = f_19 * sik_205[k]
                   + f_12 * ski0_165[k]
                   - f_13 * ski1_165[k]
                   + f_3 * pc_x[k] * skk_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, sik_207, sik_208, sik_209, \
                         ski0_167, ski1_167, skk_200, skk_207, skk_208, \
                         skk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * skk_200[k];

        t_252[k] = f_19 * sik_207[k]
                   + f_12 * ski0_167[k]
                   - f_13 * ski1_167[k]
                   + f_3 * pc_x[k] * skk_207[k];

        t_253[k] = f_19 * sik_208[k]
                   + f_3 * pc_x[k] * skk_208[k];

        t_254[k] = f_19 * sik_209[k]
                   + f_3 * pc_x[k] * skk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, sik_210, sik_211, sik_212, \
                         sik_213, sik_214, skk_210, skk_211, skk_212, skk_213, \
                         skk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_19 * sik_210[k]
                   + f_3 * pc_x[k] * skk_210[k];

        t_256[k] = f_19 * sik_211[k]
                   + f_3 * pc_x[k] * skk_211[k];

        t_257[k] = f_19 * sik_212[k]
                   + f_3 * pc_x[k] * skk_212[k];

        t_258[k] = f_19 * sik_213[k]
                   + f_3 * pc_x[k] * skk_213[k];

        t_259[k] = f_19 * sik_214[k]
                   + f_3 * pc_x[k] * skk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, sik_100, sik_215, \
                         ski0_161, ski0_163, ski1_161, ski1_163, skk_208, skk_210, \
                         skk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_19 * sik_215[k]
                   + f_3 * pc_x[k] * skk_215[k];

        t_261[k] = f_1 * ski0_161[k]
                   - f_2 * ski1_161[k]
                   + f_3 * pc_y[k] * skk_208[k];

        t_262[k] = f_16 * sik_100[k]
                   + f_3 * pc_z[k] * skk_208[k];

        t_263[k] = f_4 * ski0_163[k]
                   - f_5 * ski1_163[k]
                   + f_3 * pc_y[k] * skk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, ski0_164, ski0_165, ski0_166, ski1_164, \
                         ski1_165, ski1_166, skk_211, skk_212, \
                         skk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * ski0_164[k]
                   - f_7 * ski1_164[k]
                   + f_3 * pc_y[k] * skk_211[k];

        t_265[k] = f_8 * ski0_165[k]
                   - f_9 * ski1_165[k]
                   + f_3 * pc_y[k] * skk_212[k];

        t_266[k] = f_10 * ski0_166[k]
                   - f_11 * ski1_166[k]
                   + f_3 * pc_y[k] * skk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, sik_107, sik_216, \
                         ski0_167, ski0_168, ski1_167, ski1_168, skk_214, skk_215, \
                         skk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * ski0_167[k]
                   - f_13 * ski1_167[k]
                   + f_3 * pc_y[k] * skk_214[k];

        t_268[k] = f_3 * pc_y[k] * skk_215[k];

        t_269[k] = f_16 * sik_107[k]
                   + f_1 * ski0_167[k]
                   - f_2 * ski1_167[k]
                   + f_3 * pc_z[k] * skk_215[k];

        t_270[k] = f_18 * sik_216[k]
                   + f_1 * ski0_168[k]
                   - f_2 * ski1_168[k]
                   + f_3 * pc_x[k] * skk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, sik_108, sik_110, \
                         sik_219, ski0_171, ski1_171, skk_216, skk_218, \
                         skk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * sik_108[k]
                   + f_3 * pc_y[k] * skk_216[k];

        t_272[k] = f_3 * pc_z[k] * skk_216[k];

        t_273[k] = f_18 * sik_219[k]
                   + f_4 * ski0_171[k]
                   - f_5 * ski1_171[k]
                   + f_3 * pc_x[k] * skk_219[k];

        t_274[k] = f_17 * sik_110[k]
                   + f_3 * pc_y[k] * skk_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, sik_221, sik_222, ski0_173, \
                         ski0_174, ski1_173, ski1_174, skk_219, skk_221, \
                         skk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_18 * sik_221[k]
                   + f_4 * ski0_173[k]
                   - f_5 * ski1_173[k]
                   + f_3 * pc_x[k] * skk_221[k];

        t_276[k] = f_18 * sik_222[k]
                   + f_6 * ski0_174[k]
                   - f_7 * ski1_174[k]
                   + f_3 * pc_x[k] * skk_222[k];

        t_277[k] = f_3 * pc_z[k] * skk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, sik_113, sik_225, sik_226, ski0_177, \
                         ski0_178, ski1_177, ski1_178, skk_221, skk_225, \
                         skk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * sik_113[k]
                   + f_3 * pc_y[k] * skk_221[k];

        t_279[k] = f_18 * sik_225[k]
                   + f_6 * ski0_177[k]
                   - f_7 * ski1_177[k]
                   + f_3 * pc_x[k] * skk_225[k];

        t_280[k] = f_18 * sik_226[k]
                   + f_8 * ski0_178[k]
                   - f_9 * ski1_178[k]
                   + f_3 * pc_x[k] * skk_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, pc_z, sik_117, sik_228, ski0_180, \
                         ski1_180, skk_222, skk_225, skk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * skk_222[k];

        t_282[k] = f_18 * sik_228[k]
                   + f_8 * ski0_180[k]
                   - f_9 * ski1_180[k]
                   + f_3 * pc_x[k] * skk_228[k];

        t_283[k] = f_17 * sik_117[k]
                   + f_3 * pc_y[k] * skk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_z, sik_230, sik_231, ski0_182, \
                         ski0_183, ski1_182, ski1_183, skk_226, skk_230, \
                         skk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_18 * sik_230[k]
                   + f_8 * ski0_182[k]
                   - f_9 * ski1_182[k]
                   + f_3 * pc_x[k] * skk_230[k];

        t_285[k] = f_18 * sik_231[k]
                   + f_10 * ski0_183[k]
                   - f_11 * ski1_183[k]
                   + f_3 * pc_x[k] * skk_231[k];

        t_286[k] = f_3 * pc_z[k] * skk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, sik_122, sik_233, sik_234, ski0_185, \
                         ski0_186, ski1_185, ski1_186, skk_230, skk_233, \
                         skk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_18 * sik_233[k]
                   + f_10 * ski0_185[k]
                   - f_11 * ski1_185[k]
                   + f_3 * pc_x[k] * skk_233[k];

        t_288[k] = f_18 * sik_234[k]
                   + f_10 * ski0_186[k]
                   - f_11 * ski1_186[k]
                   + f_3 * pc_x[k] * skk_234[k];

        t_289[k] = f_17 * sik_122[k]
                   + f_3 * pc_y[k] * skk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, sik_236, sik_237, ski0_188, \
                         ski0_189, ski1_188, ski1_189, skk_231, skk_236, \
                         skk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_18 * sik_236[k]
                   + f_10 * ski0_188[k]
                   - f_11 * ski1_188[k]
                   + f_3 * pc_x[k] * skk_236[k];

        t_291[k] = f_18 * sik_237[k]
                   + f_12 * ski0_189[k]
                   - f_13 * ski1_189[k]
                   + f_3 * pc_x[k] * skk_237[k];

        t_292[k] = f_3 * pc_z[k] * skk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, sik_239, sik_240, sik_241, ski0_191, \
                         ski0_192, ski0_193, ski1_191, ski1_192, ski1_193, skk_239, skk_240, \
                         skk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_18 * sik_239[k]
                   + f_12 * ski0_191[k]
                   - f_13 * ski1_191[k]
                   + f_3 * pc_x[k] * skk_239[k];

        t_294[k] = f_18 * sik_240[k]
                   + f_12 * ski0_192[k]
                   - f_13 * ski1_192[k]
                   + f_3 * pc_x[k] * skk_240[k];

        t_295[k] = f_18 * sik_241[k]
                   + f_12 * ski0_193[k]
                   - f_13 * ski1_193[k]
                   + f_3 * pc_x[k] * skk_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pc_x, pc_y, sik_128, sik_243, sik_244, \
                         sik_245, ski0_195, ski1_195, skk_236, skk_243, skk_244, \
                         skk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * sik_128[k]
                   + f_3 * pc_y[k] * skk_236[k];

        t_297[k] = f_18 * sik_243[k]
                   + f_12 * ski0_195[k]
                   - f_13 * ski1_195[k]
                   + f_3 * pc_x[k] * skk_243[k];

        t_298[k] = f_18 * sik_244[k]
                   + f_3 * pc_x[k] * skk_244[k];

        t_299[k] = f_18 * sik_245[k]
                   + f_3 * pc_x[k] * skk_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, sik_246, sik_247, sik_248, \
                         sik_249, sik_250, skk_246, skk_247, skk_248, skk_249, \
                         skk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_18 * sik_246[k]
                   + f_3 * pc_x[k] * skk_246[k];

        t_301[k] = f_18 * sik_247[k]
                   + f_3 * pc_x[k] * skk_247[k];

        t_302[k] = f_18 * sik_248[k]
                   + f_3 * pc_x[k] * skk_248[k];

        t_303[k] = f_18 * sik_249[k]
                   + f_3 * pc_x[k] * skk_249[k];

        t_304[k] = f_18 * sik_250[k]
                   + f_3 * pc_x[k] * skk_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_y, pc_z, sik_136, sik_251, ski0_189, \
                         ski1_189, skk_244, skk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_18 * sik_251[k]
                   + f_3 * pc_x[k] * skk_251[k];

        t_306[k] = f_17 * sik_136[k]
                   + f_1 * ski0_189[k]
                   - f_2 * ski1_189[k]
                   + f_3 * pc_y[k] * skk_244[k];

        t_307[k] = f_3 * pc_z[k] * skk_244[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_y, sik_138, sik_139, sik_140, ski0_191, \
                         ski0_192, ski0_193, ski1_191, ski1_192, ski1_193, skk_246, skk_247, \
                         skk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_17 * sik_138[k]
                   + f_4 * ski0_191[k]
                   - f_5 * ski1_191[k]
                   + f_3 * pc_y[k] * skk_246[k];

        t_309[k] = f_17 * sik_139[k]
                   + f_6 * ski0_192[k]
                   - f_7 * ski1_192[k]
                   + f_3 * pc_y[k] * skk_247[k];

        t_310[k] = f_17 * sik_140[k]
                   + f_8 * ski0_193[k]
                   - f_9 * ski1_193[k]
                   + f_3 * pc_y[k] * skk_248[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, sik_141, sik_142, sik_143, \
                         ski0_194, ski0_195, ski1_194, ski1_195, skk_249, skk_250, \
                         skk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_17 * sik_141[k]
                   + f_10 * ski0_194[k]
                   - f_11 * ski1_194[k]
                   + f_3 * pc_y[k] * skk_249[k];

        t_312[k] = f_17 * sik_142[k]
                   + f_12 * ski0_195[k]
                   - f_13 * ski1_195[k]
                   + f_3 * pc_y[k] * skk_250[k];

        t_313[k] = f_17 * sik_143[k]
                   + f_3 * pc_y[k] * skk_251[k];

        t_314[k] = f_1 * ski0_195[k]
                   - f_2 * ski1_195[k]
                   + f_3 * pc_z[k] * skk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pb_z, pc_y, pc_z, sil0_135, sil0_138, \
                         sik_108, sik_144, sil1_135, sil1_138, \
                         skk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_z[k] * sil0_135[k]
                   - f_14 * pc_z[k] * sil1_135[k];

        t_316[k] = f_16 * sik_144[k]
                   + f_3 * pc_y[k] * skk_252[k];

        t_317[k] = f_15 * sik_108[k]
                   + f_3 * pc_z[k] * skk_252[k];

        t_318[k] = pb_z[k] * sil0_138[k]
                   - f_14 * pc_z[k] * sil1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_z, pc_x, pc_y, pc_z, sil0_141, sik_146, \
                         sik_257, sil1_141, ski0_201, ski1_201, skk_254, \
                         skk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * sik_146[k]
                   + f_3 * pc_y[k] * skk_254[k];

        t_320[k] = f_18 * sik_257[k]
                   + f_4 * ski0_201[k]
                   - f_5 * ski1_201[k]
                   + f_3 * pc_x[k] * skk_257[k];

        t_321[k] = pb_z[k] * sil0_141[k]
                   - f_14 * pc_z[k] * sil1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, sik_111, sik_149, sik_261, \
                         ski0_205, ski1_205, skk_255, skk_257, \
                         skk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * sik_111[k]
                   + f_3 * pc_z[k] * skk_255[k];

        t_323[k] = f_16 * sik_149[k]
                   + f_3 * pc_y[k] * skk_257[k];

        t_324[k] = f_18 * sik_261[k]
                   + f_6 * ski0_205[k]
                   - f_7 * ski1_205[k]
                   + f_3 * pc_x[k] * skk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_y, pc_z, sil0_145, sil0_147, \
                         sik_114, sik_115, sik_153, sil1_145, sil1_147, skk_258, \
                         skk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pb_z[k] * sil0_145[k]
                   - f_14 * pc_z[k] * sil1_145[k];

        t_326[k] = f_15 * sik_114[k]
                   + f_3 * pc_z[k] * skk_258[k];

        t_327[k] = pb_z[k] * sil0_147[k]
                   + f_16 * sik_115[k]
                   - f_14 * pc_z[k] * sil1_147[k];

        t_328[k] = f_16 * sik_153[k]
                   + f_3 * pc_y[k] * skk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_z, pc_x, pc_z, sil0_150, sik_118, sik_266, \
                         sil1_150, ski0_210, ski1_210, skk_262, \
                         skk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_18 * sik_266[k]
                   + f_8 * ski0_210[k]
                   - f_9 * ski1_210[k]
                   + f_3 * pc_x[k] * skk_266[k];

        t_330[k] = pb_z[k] * sil0_150[k]
                   - f_14 * pc_z[k] * sil1_150[k];

        t_331[k] = f_15 * sik_118[k]
                   + f_3 * pc_z[k] * skk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_y, pc_z, sil0_152, sil0_153, sik_119, \
                         sik_120, sik_158, sil1_152, sil1_153, \
                         skk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_z[k] * sil0_152[k]
                   + f_16 * sik_119[k]
                   - f_14 * pc_z[k] * sil1_152[k];

        t_333[k] = pb_z[k] * sil0_153[k]
                   + f_17 * sik_120[k]
                   - f_14 * pc_z[k] * sil1_153[k];

        t_334[k] = f_16 * sik_158[k]
                   + f_3 * pc_y[k] * skk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_z, pc_x, pc_z, sil0_156, sik_123, sik_272, \
                         sil1_156, ski0_216, ski1_216, skk_267, \
                         skk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_18 * sik_272[k]
                   + f_10 * ski0_216[k]
                   - f_11 * ski1_216[k]
                   + f_3 * pc_x[k] * skk_272[k];

        t_336[k] = pb_z[k] * sil0_156[k]
                   - f_14 * pc_z[k] * sil1_156[k];

        t_337[k] = f_15 * sik_123[k]
                   + f_3 * pc_z[k] * skk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_z, sil0_158, sil0_159, sil0_160, \
                         sik_124, sik_125, sik_126, sil1_158, sil1_159, \
                         sil1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_z[k] * sil0_158[k]
                   + f_16 * sik_124[k]
                   - f_14 * pc_z[k] * sil1_158[k];

        t_339[k] = pb_z[k] * sil0_159[k]
                   + f_17 * sik_125[k]
                   - f_14 * pc_z[k] * sil1_159[k];

        t_340[k] = pb_z[k] * sil0_160[k]
                   + f_18 * sik_126[k]
                   - f_14 * pc_z[k] * sil1_160[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_171 = buffer.data(sil0 + 171);
    const auto *sil0_225 = buffer.data(sil0 + 225);
    const auto *sil0_228 = buffer.data(sil0 + 228);
    const auto *sil0_230 = buffer.data(sil0 + 230);
    const auto *sil0_231 = buffer.data(sil0 + 231);
    const auto *sil0_234 = buffer.data(sil0 + 234);
    const auto *sil0_235 = buffer.data(sil0 + 235);
    const auto *sil0_237 = buffer.data(sil0 + 237);
    const auto *sil0_239 = buffer.data(sil0 + 239);
    const auto *sil0_240 = buffer.data(sil0 + 240);
    const auto *sil0_242 = buffer.data(sil0 + 242);
    const auto *sil0_243 = buffer.data(sil0 + 243);
    const auto *sil0_245 = buffer.data(sil0 + 245);
    const auto *sil0_246 = buffer.data(sil0 + 246);
    const auto *sil0_248 = buffer.data(sil0 + 248);
    const auto *sil0_249 = buffer.data(sil0 + 249);
    const auto *sil0_250 = buffer.data(sil0 + 250);
    const auto *sil0_252 = buffer.data(sil0 + 252);
    const auto *sil0_269 = buffer.data(sil0 + 269);

    const auto *sik_136 = buffer.data(sik + 136);
    const auto *sik_143 = buffer.data(sik + 143);
    const auto *sik_144 = buffer.data(sik + 144);
    const auto *sik_147 = buffer.data(sik + 147);
    const auto *sik_150 = buffer.data(sik + 150);
    const auto *sik_154 = buffer.data(sik + 154);
    const auto *sik_159 = buffer.data(sik + 159);
    const auto *sik_164 = buffer.data(sik + 164);
    const auto *sik_172 = buffer.data(sik + 172);
    const auto *sik_174 = buffer.data(sik + 174);
    const auto *sik_175 = buffer.data(sik + 175);
    const auto *sik_176 = buffer.data(sik + 176);
    const auto *sik_177 = buffer.data(sik + 177);
    const auto *sik_178 = buffer.data(sik + 178);
    const auto *sik_179 = buffer.data(sik + 179);
    const auto *sik_180 = buffer.data(sik + 180);
    const auto *sik_181 = buffer.data(sik + 181);
    const auto *sik_182 = buffer.data(sik + 182);
    const auto *sik_183 = buffer.data(sik + 183);
    const auto *sik_185 = buffer.data(sik + 185);
    const auto *sik_186 = buffer.data(sik + 186);
    const auto *sik_188 = buffer.data(sik + 188);
    const auto *sik_189 = buffer.data(sik + 189);
    const auto *sik_190 = buffer.data(sik + 190);
    const auto *sik_192 = buffer.data(sik + 192);
    const auto *sik_193 = buffer.data(sik + 193);
    const auto *sik_194 = buffer.data(sik + 194);
    const auto *sik_195 = buffer.data(sik + 195);
    const auto *sik_197 = buffer.data(sik + 197);
    const auto *sik_198 = buffer.data(sik + 198);
    const auto *sik_199 = buffer.data(sik + 199);
    const auto *sik_200 = buffer.data(sik + 200);
    const auto *sik_208 = buffer.data(sik + 208);
    const auto *sik_210 = buffer.data(sik + 210);
    const auto *sik_211 = buffer.data(sik + 211);
    const auto *sik_212 = buffer.data(sik + 212);
    const auto *sik_213 = buffer.data(sik + 213);
    const auto *sik_214 = buffer.data(sik + 214);
    const auto *sik_215 = buffer.data(sik + 215);
    const auto *sik_216 = buffer.data(sik + 216);
    const auto *sik_218 = buffer.data(sik + 218);
    const auto *sik_279 = buffer.data(sik + 279);
    const auto *sik_280 = buffer.data(sik + 280);
    const auto *sik_281 = buffer.data(sik + 281);
    const auto *sik_282 = buffer.data(sik + 282);
    const auto *sik_283 = buffer.data(sik + 283);
    const auto *sik_284 = buffer.data(sik + 284);
    const auto *sik_285 = buffer.data(sik + 285);
    const auto *sik_286 = buffer.data(sik + 286);
    const auto *sik_287 = buffer.data(sik + 287);
    const auto *sik_316 = buffer.data(sik + 316);
    const auto *sik_317 = buffer.data(sik + 317);
    const auto *sik_318 = buffer.data(sik + 318);
    const auto *sik_319 = buffer.data(sik + 319);
    const auto *sik_320 = buffer.data(sik + 320);
    const auto *sik_321 = buffer.data(sik + 321);
    const auto *sik_322 = buffer.data(sik + 322);
    const auto *sik_323 = buffer.data(sik + 323);
    const auto *sik_324 = buffer.data(sik + 324);
    const auto *sik_327 = buffer.data(sik + 327);
    const auto *sik_329 = buffer.data(sik + 329);
    const auto *sik_330 = buffer.data(sik + 330);
    const auto *sik_333 = buffer.data(sik + 333);
    const auto *sik_334 = buffer.data(sik + 334);
    const auto *sik_336 = buffer.data(sik + 336);
    const auto *sik_338 = buffer.data(sik + 338);
    const auto *sik_339 = buffer.data(sik + 339);
    const auto *sik_341 = buffer.data(sik + 341);
    const auto *sik_342 = buffer.data(sik + 342);
    const auto *sik_344 = buffer.data(sik + 344);
    const auto *sik_345 = buffer.data(sik + 345);
    const auto *sik_347 = buffer.data(sik + 347);
    const auto *sik_348 = buffer.data(sik + 348);
    const auto *sik_349 = buffer.data(sik + 349);
    const auto *sik_351 = buffer.data(sik + 351);
    const auto *sik_352 = buffer.data(sik + 352);
    const auto *sik_353 = buffer.data(sik + 353);
    const auto *sik_354 = buffer.data(sik + 354);
    const auto *sik_355 = buffer.data(sik + 355);
    const auto *sik_356 = buffer.data(sik + 356);
    const auto *sik_357 = buffer.data(sik + 357);
    const auto *sik_358 = buffer.data(sik + 358);
    const auto *sik_359 = buffer.data(sik + 359);
    const auto *sik_360 = buffer.data(sik + 360);
    const auto *sik_363 = buffer.data(sik + 363);
    const auto *sik_365 = buffer.data(sik + 365);

    const auto *sil1_171 = buffer.data(sil1 + 171);
    const auto *sil1_225 = buffer.data(sil1 + 225);
    const auto *sil1_228 = buffer.data(sil1 + 228);
    const auto *sil1_230 = buffer.data(sil1 + 230);
    const auto *sil1_231 = buffer.data(sil1 + 231);
    const auto *sil1_234 = buffer.data(sil1 + 234);
    const auto *sil1_235 = buffer.data(sil1 + 235);
    const auto *sil1_237 = buffer.data(sil1 + 237);
    const auto *sil1_239 = buffer.data(sil1 + 239);
    const auto *sil1_240 = buffer.data(sil1 + 240);
    const auto *sil1_242 = buffer.data(sil1 + 242);
    const auto *sil1_243 = buffer.data(sil1 + 243);
    const auto *sil1_245 = buffer.data(sil1 + 245);
    const auto *sil1_246 = buffer.data(sil1 + 246);
    const auto *sil1_248 = buffer.data(sil1 + 248);
    const auto *sil1_249 = buffer.data(sil1 + 249);
    const auto *sil1_250 = buffer.data(sil1 + 250);
    const auto *sil1_252 = buffer.data(sil1 + 252);
    const auto *sil1_269 = buffer.data(sil1 + 269);

    const auto *ski0_219 = buffer.data(ski0 + 219);
    const auto *ski0_220 = buffer.data(ski0 + 220);
    const auto *ski0_221 = buffer.data(ski0 + 221);
    const auto *ski0_222 = buffer.data(ski0 + 222);
    const auto *ski0_223 = buffer.data(ski0 + 223);
    const auto *ski0_245 = buffer.data(ski0 + 245);
    const auto *ski0_247 = buffer.data(ski0 + 247);
    const auto *ski0_248 = buffer.data(ski0 + 248);
    const auto *ski0_249 = buffer.data(ski0 + 249);
    const auto *ski0_250 = buffer.data(ski0 + 250);
    const auto *ski0_251 = buffer.data(ski0 + 251);
    const auto *ski0_252 = buffer.data(ski0 + 252);
    const auto *ski0_255 = buffer.data(ski0 + 255);
    const auto *ski0_257 = buffer.data(ski0 + 257);
    const auto *ski0_258 = buffer.data(ski0 + 258);
    const auto *ski0_261 = buffer.data(ski0 + 261);
    const auto *ski0_262 = buffer.data(ski0 + 262);
    const auto *ski0_264 = buffer.data(ski0 + 264);
    const auto *ski0_266 = buffer.data(ski0 + 266);
    const auto *ski0_267 = buffer.data(ski0 + 267);
    const auto *ski0_269 = buffer.data(ski0 + 269);
    const auto *ski0_270 = buffer.data(ski0 + 270);
    const auto *ski0_272 = buffer.data(ski0 + 272);
    const auto *ski0_273 = buffer.data(ski0 + 273);
    const auto *ski0_275 = buffer.data(ski0 + 275);
    const auto *ski0_276 = buffer.data(ski0 + 276);
    const auto *ski0_277 = buffer.data(ski0 + 277);
    const auto *ski0_278 = buffer.data(ski0 + 278);
    const auto *ski0_279 = buffer.data(ski0 + 279);
    const auto *ski0_280 = buffer.data(ski0 + 280);
    const auto *ski0_283 = buffer.data(ski0 + 283);
    const auto *ski0_285 = buffer.data(ski0 + 285);

    const auto *ski1_219 = buffer.data(ski1 + 219);
    const auto *ski1_220 = buffer.data(ski1 + 220);
    const auto *ski1_221 = buffer.data(ski1 + 221);
    const auto *ski1_222 = buffer.data(ski1 + 222);
    const auto *ski1_223 = buffer.data(ski1 + 223);
    const auto *ski1_245 = buffer.data(ski1 + 245);
    const auto *ski1_247 = buffer.data(ski1 + 247);
    const auto *ski1_248 = buffer.data(ski1 + 248);
    const auto *ski1_249 = buffer.data(ski1 + 249);
    const auto *ski1_250 = buffer.data(ski1 + 250);
    const auto *ski1_251 = buffer.data(ski1 + 251);
    const auto *ski1_252 = buffer.data(ski1 + 252);
    const auto *ski1_255 = buffer.data(ski1 + 255);
    const auto *ski1_257 = buffer.data(ski1 + 257);
    const auto *ski1_258 = buffer.data(ski1 + 258);
    const auto *ski1_261 = buffer.data(ski1 + 261);
    const auto *ski1_262 = buffer.data(ski1 + 262);
    const auto *ski1_264 = buffer.data(ski1 + 264);
    const auto *ski1_266 = buffer.data(ski1 + 266);
    const auto *ski1_267 = buffer.data(ski1 + 267);
    const auto *ski1_269 = buffer.data(ski1 + 269);
    const auto *ski1_270 = buffer.data(ski1 + 270);
    const auto *ski1_272 = buffer.data(ski1 + 272);
    const auto *ski1_273 = buffer.data(ski1 + 273);
    const auto *ski1_275 = buffer.data(ski1 + 275);
    const auto *ski1_276 = buffer.data(ski1 + 276);
    const auto *ski1_277 = buffer.data(ski1 + 277);
    const auto *ski1_278 = buffer.data(ski1 + 278);
    const auto *ski1_279 = buffer.data(ski1 + 279);
    const auto *ski1_280 = buffer.data(ski1 + 280);
    const auto *ski1_283 = buffer.data(ski1 + 283);
    const auto *ski1_285 = buffer.data(ski1 + 285);

    const auto *skk_272 = buffer.data(skk + 272);
    const auto *skk_279 = buffer.data(skk + 279);
    const auto *skk_280 = buffer.data(skk + 280);
    const auto *skk_281 = buffer.data(skk + 281);
    const auto *skk_282 = buffer.data(skk + 282);
    const auto *skk_283 = buffer.data(skk + 283);
    const auto *skk_284 = buffer.data(skk + 284);
    const auto *skk_285 = buffer.data(skk + 285);
    const auto *skk_286 = buffer.data(skk + 286);
    const auto *skk_287 = buffer.data(skk + 287);
    const auto *skk_288 = buffer.data(skk + 288);
    const auto *skk_290 = buffer.data(skk + 290);
    const auto *skk_291 = buffer.data(skk + 291);
    const auto *skk_293 = buffer.data(skk + 293);
    const auto *skk_294 = buffer.data(skk + 294);
    const auto *skk_297 = buffer.data(skk + 297);
    const auto *skk_298 = buffer.data(skk + 298);
    const auto *skk_302 = buffer.data(skk + 302);
    const auto *skk_303 = buffer.data(skk + 303);
    const auto *skk_308 = buffer.data(skk + 308);
    const auto *skk_316 = buffer.data(skk + 316);
    const auto *skk_317 = buffer.data(skk + 317);
    const auto *skk_318 = buffer.data(skk + 318);
    const auto *skk_319 = buffer.data(skk + 319);
    const auto *skk_320 = buffer.data(skk + 320);
    const auto *skk_321 = buffer.data(skk + 321);
    const auto *skk_322 = buffer.data(skk + 322);
    const auto *skk_323 = buffer.data(skk + 323);
    const auto *skk_324 = buffer.data(skk + 324);
    const auto *skk_326 = buffer.data(skk + 326);
    const auto *skk_327 = buffer.data(skk + 327);
    const auto *skk_329 = buffer.data(skk + 329);
    const auto *skk_330 = buffer.data(skk + 330);
    const auto *skk_333 = buffer.data(skk + 333);
    const auto *skk_334 = buffer.data(skk + 334);
    const auto *skk_336 = buffer.data(skk + 336);
    const auto *skk_338 = buffer.data(skk + 338);
    const auto *skk_339 = buffer.data(skk + 339);
    const auto *skk_341 = buffer.data(skk + 341);
    const auto *skk_342 = buffer.data(skk + 342);
    const auto *skk_344 = buffer.data(skk + 344);
    const auto *skk_345 = buffer.data(skk + 345);
    const auto *skk_347 = buffer.data(skk + 347);
    const auto *skk_348 = buffer.data(skk + 348);
    const auto *skk_349 = buffer.data(skk + 349);
    const auto *skk_351 = buffer.data(skk + 351);
    const auto *skk_352 = buffer.data(skk + 352);
    const auto *skk_353 = buffer.data(skk + 353);
    const auto *skk_354 = buffer.data(skk + 354);
    const auto *skk_355 = buffer.data(skk + 355);
    const auto *skk_356 = buffer.data(skk + 356);
    const auto *skk_357 = buffer.data(skk + 357);
    const auto *skk_358 = buffer.data(skk + 358);
    const auto *skk_359 = buffer.data(skk + 359);
    const auto *skk_360 = buffer.data(skk + 360);
    const auto *skk_362 = buffer.data(skk + 362);
    const auto *skk_363 = buffer.data(skk + 363);
    const auto *skk_365 = buffer.data(skk + 365);

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, sik_164, sik_279, sik_280, \
                         sik_281, ski0_223, ski1_223, skk_272, skk_279, skk_280, \
                         skk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * sik_164[k]
                   + f_3 * pc_y[k] * skk_272[k];

        t_342[k] = f_18 * sik_279[k]
                   + f_12 * ski0_223[k]
                   - f_13 * ski1_223[k]
                   + f_3 * pc_x[k] * skk_279[k];

        t_343[k] = f_18 * sik_280[k]
                   + f_3 * pc_x[k] * skk_280[k];

        t_344[k] = f_18 * sik_281[k]
                   + f_3 * pc_x[k] * skk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, sik_282, sik_283, sik_284, \
                         sik_285, sik_286, skk_282, skk_283, skk_284, skk_285, \
                         skk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_18 * sik_282[k]
                   + f_3 * pc_x[k] * skk_282[k];

        t_346[k] = f_18 * sik_283[k]
                   + f_3 * pc_x[k] * skk_283[k];

        t_347[k] = f_18 * sik_284[k]
                   + f_3 * pc_x[k] * skk_284[k];

        t_348[k] = f_18 * sik_285[k]
                   + f_3 * pc_x[k] * skk_285[k];

        t_349[k] = f_18 * sik_286[k]
                   + f_3 * pc_x[k] * skk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pb_z, pc_x, pc_z, sil0_171, sik_136, sik_287, \
                         sil1_171, skk_280, skk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_18 * sik_287[k]
                   + f_3 * pc_x[k] * skk_287[k];

        t_351[k] = pb_z[k] * sil0_171[k]
                   - f_14 * pc_z[k] * sil1_171[k];

        t_352[k] = f_15 * sik_136[k]
                   + f_3 * pc_z[k] * skk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, sik_174, sik_175, sik_176, ski0_219, \
                         ski0_220, ski0_221, ski1_219, ski1_220, ski1_221, skk_282, skk_283, \
                         skk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * sik_174[k]
                   + f_4 * ski0_219[k]
                   - f_5 * ski1_219[k]
                   + f_3 * pc_y[k] * skk_282[k];

        t_354[k] = f_16 * sik_175[k]
                   + f_6 * ski0_220[k]
                   - f_7 * ski1_220[k]
                   + f_3 * pc_y[k] * skk_283[k];

        t_355[k] = f_16 * sik_176[k]
                   + f_8 * ski0_221[k]
                   - f_9 * ski1_221[k]
                   + f_3 * pc_y[k] * skk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, sik_177, sik_178, sik_179, ski0_222, \
                         ski0_223, ski1_222, ski1_223, skk_285, skk_286, \
                         skk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * sik_177[k]
                   + f_10 * ski0_222[k]
                   - f_11 * ski1_222[k]
                   + f_3 * pc_y[k] * skk_285[k];

        t_357[k] = f_16 * sik_178[k]
                   + f_12 * ski0_223[k]
                   - f_13 * ski1_223[k]
                   + f_3 * pc_y[k] * skk_286[k];

        t_358[k] = f_16 * sik_179[k]
                   + f_3 * pc_y[k] * skk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, sil0_225, sik_143, \
                         sik_144, sik_180, sil1_225, ski0_223, ski1_223, skk_287, \
                         skk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * sik_143[k]
                   + f_1 * ski0_223[k]
                   - f_2 * ski1_223[k]
                   + f_3 * pc_z[k] * skk_287[k];

        t_360[k] = pb_y[k] * sil0_225[k]
                   - f_14 * pc_y[k] * sil1_225[k];

        t_361[k] = f_15 * sik_180[k]
                   + f_3 * pc_y[k] * skk_288[k];

        t_362[k] = f_16 * sik_144[k]
                   + f_3 * pc_z[k] * skk_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, sil0_228, sil0_230, sil0_231, \
                         sik_181, sik_182, sik_183, sil1_228, sil1_230, sil1_231, \
                         skk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_y[k] * sil0_228[k]
                   + f_16 * sik_181[k]
                   - f_14 * pc_y[k] * sil1_228[k];

        t_364[k] = f_15 * sik_182[k]
                   + f_3 * pc_y[k] * skk_290[k];

        t_365[k] = pb_y[k] * sil0_230[k]
                   - f_14 * pc_y[k] * sil1_230[k];

        t_366[k] = pb_y[k] * sil0_231[k]
                   + f_17 * sik_183[k]
                   - f_14 * pc_y[k] * sil1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pb_y, pc_y, pc_z, sil0_234, sil0_235, \
                         sik_147, sik_185, sik_186, sil1_234, sil1_235, skk_291, \
                         skk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * sik_147[k]
                   + f_3 * pc_z[k] * skk_291[k];

        t_368[k] = f_15 * sik_185[k]
                   + f_3 * pc_y[k] * skk_293[k];

        t_369[k] = pb_y[k] * sil0_234[k]
                   - f_14 * pc_y[k] * sil1_234[k];

        t_370[k] = pb_y[k] * sil0_235[k]
                   + f_18 * sik_186[k]
                   - f_14 * pc_y[k] * sil1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pb_y, pc_y, pc_z, sil0_237, sil0_239, \
                         sik_150, sik_188, sik_189, sil1_237, sil1_239, skk_294, \
                         skk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * sik_150[k]
                   + f_3 * pc_z[k] * skk_294[k];

        t_372[k] = pb_y[k] * sil0_237[k]
                   + f_16 * sik_188[k]
                   - f_14 * pc_y[k] * sil1_237[k];

        t_373[k] = f_15 * sik_189[k]
                   + f_3 * pc_y[k] * skk_297[k];

        t_374[k] = pb_y[k] * sil0_239[k]
                   - f_14 * pc_y[k] * sil1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_y, pc_y, pc_z, sil0_240, sil0_242, sik_154, \
                         sik_190, sik_192, sil1_240, sil1_242, \
                         skk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pb_y[k] * sil0_240[k]
                   + f_19 * sik_190[k]
                   - f_14 * pc_y[k] * sil1_240[k];

        t_376[k] = f_16 * sik_154[k]
                   + f_3 * pc_z[k] * skk_298[k];

        t_377[k] = pb_y[k] * sil0_242[k]
                   + f_17 * sik_192[k]
                   - f_14 * pc_y[k] * sil1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_y, sil0_243, sil0_245, sil0_246, \
                         sik_193, sik_194, sik_195, sil1_243, sil1_245, sil1_246, \
                         skk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * sil0_243[k]
                   + f_16 * sik_193[k]
                   - f_14 * pc_y[k] * sil1_243[k];

        t_379[k] = f_15 * sik_194[k]
                   + f_3 * pc_y[k] * skk_302[k];

        t_380[k] = pb_y[k] * sil0_245[k]
                   - f_14 * pc_y[k] * sil1_245[k];

        t_381[k] = pb_y[k] * sil0_246[k]
                   + f_20 * sik_195[k]
                   - f_14 * pc_y[k] * sil1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_y, pc_y, pc_z, sil0_248, sil0_249, sik_159, \
                         sik_197, sik_198, sil1_248, sil1_249, \
                         skk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * sik_159[k]
                   + f_3 * pc_z[k] * skk_303[k];

        t_383[k] = pb_y[k] * sil0_248[k]
                   + f_18 * sik_197[k]
                   - f_14 * pc_y[k] * sil1_248[k];

        t_384[k] = pb_y[k] * sil0_249[k]
                   + f_17 * sik_198[k]
                   - f_14 * pc_y[k] * sil1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_y, pc_x, pc_y, sil0_250, sil0_252, \
                         sik_199, sik_200, sik_316, sil1_250, sil1_252, skk_308, \
                         skk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * sil0_250[k]
                   + f_16 * sik_199[k]
                   - f_14 * pc_y[k] * sil1_250[k];

        t_386[k] = f_15 * sik_200[k]
                   + f_3 * pc_y[k] * skk_308[k];

        t_387[k] = pb_y[k] * sil0_252[k]
                   - f_14 * pc_y[k] * sil1_252[k];

        t_388[k] = f_18 * sik_316[k]
                   + f_3 * pc_x[k] * skk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, sik_317, sik_318, sik_319, \
                         sik_320, sik_321, skk_317, skk_318, skk_319, skk_320, \
                         skk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_18 * sik_317[k]
                   + f_3 * pc_x[k] * skk_317[k];

        t_390[k] = f_18 * sik_318[k]
                   + f_3 * pc_x[k] * skk_318[k];

        t_391[k] = f_18 * sik_319[k]
                   + f_3 * pc_x[k] * skk_319[k];

        t_392[k] = f_18 * sik_320[k]
                   + f_3 * pc_x[k] * skk_320[k];

        t_393[k] = f_18 * sik_321[k]
                   + f_3 * pc_x[k] * skk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, sik_172, sik_208, \
                         sik_322, sik_323, ski0_245, ski1_245, skk_316, skk_322, \
                         skk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_18 * sik_322[k]
                   + f_3 * pc_x[k] * skk_322[k];

        t_395[k] = f_18 * sik_323[k]
                   + f_3 * pc_x[k] * skk_323[k];

        t_396[k] = f_15 * sik_208[k]
                   + f_1 * ski0_245[k]
                   - f_2 * ski1_245[k]
                   + f_3 * pc_y[k] * skk_316[k];

        t_397[k] = f_16 * sik_172[k]
                   + f_3 * pc_z[k] * skk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, sik_210, sik_211, sik_212, ski0_247, \
                         ski0_248, ski0_249, ski1_247, ski1_248, ski1_249, skk_318, skk_319, \
                         skk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * sik_210[k]
                   + f_4 * ski0_247[k]
                   - f_5 * ski1_247[k]
                   + f_3 * pc_y[k] * skk_318[k];

        t_399[k] = f_15 * sik_211[k]
                   + f_6 * ski0_248[k]
                   - f_7 * ski1_248[k]
                   + f_3 * pc_y[k] * skk_319[k];

        t_400[k] = f_15 * sik_212[k]
                   + f_8 * ski0_249[k]
                   - f_9 * ski1_249[k]
                   + f_3 * pc_y[k] * skk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, sik_213, sik_214, sik_215, ski0_250, \
                         ski0_251, ski1_250, ski1_251, skk_321, skk_322, \
                         skk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * sik_213[k]
                   + f_10 * ski0_250[k]
                   - f_11 * ski1_250[k]
                   + f_3 * pc_y[k] * skk_321[k];

        t_402[k] = f_15 * sik_214[k]
                   + f_12 * ski0_251[k]
                   - f_13 * ski1_251[k]
                   + f_3 * pc_y[k] * skk_322[k];

        t_403[k] = f_15 * sik_215[k]
                   + f_3 * pc_y[k] * skk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_y, pc_x, pc_y, pc_z, sil0_269, \
                         sik_180, sik_324, sil1_269, ski0_252, ski1_252, \
                         skk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * sil0_269[k]
                   - f_14 * pc_y[k] * sil1_269[k];

        t_405[k] = f_18 * sik_324[k]
                   + f_1 * ski0_252[k]
                   - f_2 * ski1_252[k]
                   + f_3 * pc_x[k] * skk_324[k];

        t_406[k] = f_3 * pc_y[k] * skk_324[k];

        t_407[k] = f_17 * sik_180[k]
                   + f_3 * pc_z[k] * skk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, sik_327, sik_329, ski0_255, \
                         ski0_257, ski1_255, ski1_257, skk_326, skk_327, \
                         skk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_18 * sik_327[k]
                   + f_4 * ski0_255[k]
                   - f_5 * ski1_255[k]
                   + f_3 * pc_x[k] * skk_327[k];

        t_409[k] = f_3 * pc_y[k] * skk_326[k];

        t_410[k] = f_18 * sik_329[k]
                   + f_4 * ski0_257[k]
                   - f_5 * ski1_257[k]
                   + f_3 * pc_x[k] * skk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_x, pc_y, pc_z, sik_183, sik_330, ski0_258, \
                         ski1_258, skk_327, skk_329, skk_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_18 * sik_330[k]
                   + f_6 * ski0_258[k]
                   - f_7 * ski1_258[k]
                   + f_3 * pc_x[k] * skk_330[k];

        t_412[k] = f_17 * sik_183[k]
                   + f_3 * pc_z[k] * skk_327[k];

        t_413[k] = f_3 * pc_y[k] * skk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_z, sik_186, sik_333, sik_334, ski0_261, \
                         ski0_262, ski1_261, ski1_262, skk_330, skk_333, \
                         skk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_18 * sik_333[k]
                   + f_6 * ski0_261[k]
                   - f_7 * ski1_261[k]
                   + f_3 * pc_x[k] * skk_333[k];

        t_415[k] = f_18 * sik_334[k]
                   + f_8 * ski0_262[k]
                   - f_9 * ski1_262[k]
                   + f_3 * pc_x[k] * skk_334[k];

        t_416[k] = f_17 * sik_186[k]
                   + f_3 * pc_z[k] * skk_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, sik_336, sik_338, ski0_264, \
                         ski0_266, ski1_264, ski1_266, skk_333, skk_336, \
                         skk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_18 * sik_336[k]
                   + f_8 * ski0_264[k]
                   - f_9 * ski1_264[k]
                   + f_3 * pc_x[k] * skk_336[k];

        t_418[k] = f_3 * pc_y[k] * skk_333[k];

        t_419[k] = f_18 * sik_338[k]
                   + f_8 * ski0_266[k]
                   - f_9 * ski1_266[k]
                   + f_3 * pc_x[k] * skk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_x, pc_z, sik_190, sik_339, sik_341, ski0_267, \
                         ski0_269, ski1_267, ski1_269, skk_334, skk_339, \
                         skk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_18 * sik_339[k]
                   + f_10 * ski0_267[k]
                   - f_11 * ski1_267[k]
                   + f_3 * pc_x[k] * skk_339[k];

        t_421[k] = f_17 * sik_190[k]
                   + f_3 * pc_z[k] * skk_334[k];

        t_422[k] = f_18 * sik_341[k]
                   + f_10 * ski0_269[k]
                   - f_11 * ski1_269[k]
                   + f_3 * pc_x[k] * skk_341[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, sik_342, sik_344, ski0_270, \
                         ski0_272, ski1_270, ski1_272, skk_338, skk_342, \
                         skk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_18 * sik_342[k]
                   + f_10 * ski0_270[k]
                   - f_11 * ski1_270[k]
                   + f_3 * pc_x[k] * skk_342[k];

        t_424[k] = f_3 * pc_y[k] * skk_338[k];

        t_425[k] = f_18 * sik_344[k]
                   + f_10 * ski0_272[k]
                   - f_11 * ski1_272[k]
                   + f_3 * pc_x[k] * skk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, sik_195, sik_345, sik_347, ski0_273, \
                         ski0_275, ski1_273, ski1_275, skk_339, skk_345, \
                         skk_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_18 * sik_345[k]
                   + f_12 * ski0_273[k]
                   - f_13 * ski1_273[k]
                   + f_3 * pc_x[k] * skk_345[k];

        t_427[k] = f_17 * sik_195[k]
                   + f_3 * pc_z[k] * skk_339[k];

        t_428[k] = f_18 * sik_347[k]
                   + f_12 * ski0_275[k]
                   - f_13 * ski1_275[k]
                   + f_3 * pc_x[k] * skk_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, sik_348, sik_349, ski0_276, \
                         ski0_277, ski1_276, ski1_277, skk_344, skk_348, \
                         skk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_18 * sik_348[k]
                   + f_12 * ski0_276[k]
                   - f_13 * ski1_276[k]
                   + f_3 * pc_x[k] * skk_348[k];

        t_430[k] = f_18 * sik_349[k]
                   + f_12 * ski0_277[k]
                   - f_13 * ski1_277[k]
                   + f_3 * pc_x[k] * skk_349[k];

        t_431[k] = f_3 * pc_y[k] * skk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, sik_351, sik_352, sik_353, sik_354, \
                         ski0_279, ski1_279, skk_351, skk_352, skk_353, \
                         skk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_18 * sik_351[k]
                   + f_12 * ski0_279[k]
                   - f_13 * ski1_279[k]
                   + f_3 * pc_x[k] * skk_351[k];

        t_433[k] = f_18 * sik_352[k]
                   + f_3 * pc_x[k] * skk_352[k];

        t_434[k] = f_18 * sik_353[k]
                   + f_3 * pc_x[k] * skk_353[k];

        t_435[k] = f_18 * sik_354[k]
                   + f_3 * pc_x[k] * skk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, sik_355, sik_356, sik_357, \
                         sik_358, sik_359, skk_355, skk_356, skk_357, skk_358, \
                         skk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_18 * sik_355[k]
                   + f_3 * pc_x[k] * skk_355[k];

        t_437[k] = f_18 * sik_356[k]
                   + f_3 * pc_x[k] * skk_356[k];

        t_438[k] = f_18 * sik_357[k]
                   + f_3 * pc_x[k] * skk_357[k];

        t_439[k] = f_18 * sik_358[k]
                   + f_3 * pc_x[k] * skk_358[k];

        t_440[k] = f_18 * sik_359[k]
                   + f_3 * pc_x[k] * skk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, sik_208, ski0_273, ski0_275, \
                         ski0_276, ski1_273, ski1_275, ski1_276, skk_352, skk_354, \
                         skk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * ski0_273[k]
                   - f_2 * ski1_273[k]
                   + f_3 * pc_y[k] * skk_352[k];

        t_442[k] = f_17 * sik_208[k]
                   + f_3 * pc_z[k] * skk_352[k];

        t_443[k] = f_4 * ski0_275[k]
                   - f_5 * ski1_275[k]
                   + f_3 * pc_y[k] * skk_354[k];

        t_444[k] = f_6 * ski0_276[k]
                   - f_7 * ski1_276[k]
                   + f_3 * pc_y[k] * skk_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, ski0_277, ski0_278, ski0_279, \
                         ski1_277, ski1_278, ski1_279, skk_356, skk_357, skk_358, \
                         skk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * ski0_277[k]
                   - f_9 * ski1_277[k]
                   + f_3 * pc_y[k] * skk_356[k];

        t_446[k] = f_10 * ski0_278[k]
                   - f_11 * ski1_278[k]
                   + f_3 * pc_y[k] * skk_357[k];

        t_447[k] = f_12 * ski0_279[k]
                   - f_13 * ski1_279[k]
                   + f_3 * pc_y[k] * skk_358[k];

        t_448[k] = f_3 * pc_y[k] * skk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, sik_215, sik_216, \
                         sik_360, ski0_279, ski0_280, ski1_279, ski1_280, skk_359, \
                         skk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * sik_215[k]
                   + f_1 * ski0_279[k]
                   - f_2 * ski1_279[k]
                   + f_3 * pc_z[k] * skk_359[k];

        t_450[k] = f_17 * sik_360[k]
                   + f_1 * ski0_280[k]
                   - f_2 * ski1_280[k]
                   + f_3 * pc_x[k] * skk_360[k];

        t_451[k] = f_18 * sik_216[k]
                   + f_3 * pc_y[k] * skk_360[k];

        t_452[k] = f_3 * pc_z[k] * skk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, sik_218, sik_363, sik_365, ski0_283, \
                         ski0_285, ski1_283, ski1_285, skk_362, skk_363, \
                         skk_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_17 * sik_363[k]
                   + f_4 * ski0_283[k]
                   - f_5 * ski1_283[k]
                   + f_3 * pc_x[k] * skk_363[k];

        t_454[k] = f_18 * sik_218[k]
                   + f_3 * pc_y[k] * skk_362[k];

        t_455[k] = f_17 * sik_365[k]
                   + f_4 * ski0_285[k]
                   - f_5 * ski1_285[k]
                   + f_3 * pc_x[k] * skk_365[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_270 = buffer.data(sil0 + 270);
    const auto *sil0_273 = buffer.data(sil0 + 273);
    const auto *sil0_276 = buffer.data(sil0 + 276);
    const auto *sil0_280 = buffer.data(sil0 + 280);
    const auto *sil0_282 = buffer.data(sil0 + 282);
    const auto *sil0_285 = buffer.data(sil0 + 285);
    const auto *sil0_287 = buffer.data(sil0 + 287);
    const auto *sil0_288 = buffer.data(sil0 + 288);
    const auto *sil0_291 = buffer.data(sil0 + 291);
    const auto *sil0_293 = buffer.data(sil0 + 293);
    const auto *sil0_294 = buffer.data(sil0 + 294);
    const auto *sil0_295 = buffer.data(sil0 + 295);
    const auto *sil0_306 = buffer.data(sil0 + 306);

    const auto *sik_216 = buffer.data(sik + 216);
    const auto *sik_219 = buffer.data(sik + 219);
    const auto *sik_221 = buffer.data(sik + 221);
    const auto *sik_222 = buffer.data(sik + 222);
    const auto *sik_223 = buffer.data(sik + 223);
    const auto *sik_225 = buffer.data(sik + 225);
    const auto *sik_226 = buffer.data(sik + 226);
    const auto *sik_227 = buffer.data(sik + 227);
    const auto *sik_228 = buffer.data(sik + 228);
    const auto *sik_230 = buffer.data(sik + 230);
    const auto *sik_231 = buffer.data(sik + 231);
    const auto *sik_232 = buffer.data(sik + 232);
    const auto *sik_233 = buffer.data(sik + 233);
    const auto *sik_234 = buffer.data(sik + 234);
    const auto *sik_236 = buffer.data(sik + 236);
    const auto *sik_244 = buffer.data(sik + 244);
    const auto *sik_246 = buffer.data(sik + 246);
    const auto *sik_247 = buffer.data(sik + 247);
    const auto *sik_248 = buffer.data(sik + 248);
    const auto *sik_249 = buffer.data(sik + 249);
    const auto *sik_250 = buffer.data(sik + 250);
    const auto *sik_251 = buffer.data(sik + 251);
    const auto *sik_252 = buffer.data(sik + 252);
    const auto *sik_254 = buffer.data(sik + 254);
    const auto *sik_255 = buffer.data(sik + 255);
    const auto *sik_257 = buffer.data(sik + 257);
    const auto *sik_258 = buffer.data(sik + 258);
    const auto *sik_261 = buffer.data(sik + 261);
    const auto *sik_262 = buffer.data(sik + 262);
    const auto *sik_266 = buffer.data(sik + 266);
    const auto *sik_267 = buffer.data(sik + 267);
    const auto *sik_272 = buffer.data(sik + 272);
    const auto *sik_282 = buffer.data(sik + 282);
    const auto *sik_283 = buffer.data(sik + 283);
    const auto *sik_284 = buffer.data(sik + 284);
    const auto *sik_285 = buffer.data(sik + 285);
    const auto *sik_286 = buffer.data(sik + 286);
    const auto *sik_287 = buffer.data(sik + 287);
    const auto *sik_288 = buffer.data(sik + 288);
    const auto *sik_290 = buffer.data(sik + 290);
    const auto *sik_293 = buffer.data(sik + 293);
    const auto *sik_297 = buffer.data(sik + 297);
    const auto *sik_302 = buffer.data(sik + 302);
    const auto *sik_366 = buffer.data(sik + 366);
    const auto *sik_369 = buffer.data(sik + 369);
    const auto *sik_370 = buffer.data(sik + 370);
    const auto *sik_372 = buffer.data(sik + 372);
    const auto *sik_374 = buffer.data(sik + 374);
    const auto *sik_375 = buffer.data(sik + 375);
    const auto *sik_377 = buffer.data(sik + 377);
    const auto *sik_378 = buffer.data(sik + 378);
    const auto *sik_380 = buffer.data(sik + 380);
    const auto *sik_381 = buffer.data(sik + 381);
    const auto *sik_383 = buffer.data(sik + 383);
    const auto *sik_384 = buffer.data(sik + 384);
    const auto *sik_385 = buffer.data(sik + 385);
    const auto *sik_387 = buffer.data(sik + 387);
    const auto *sik_388 = buffer.data(sik + 388);
    const auto *sik_389 = buffer.data(sik + 389);
    const auto *sik_390 = buffer.data(sik + 390);
    const auto *sik_391 = buffer.data(sik + 391);
    const auto *sik_392 = buffer.data(sik + 392);
    const auto *sik_393 = buffer.data(sik + 393);
    const auto *sik_394 = buffer.data(sik + 394);
    const auto *sik_395 = buffer.data(sik + 395);
    const auto *sik_401 = buffer.data(sik + 401);
    const auto *sik_405 = buffer.data(sik + 405);
    const auto *sik_410 = buffer.data(sik + 410);
    const auto *sik_416 = buffer.data(sik + 416);
    const auto *sik_423 = buffer.data(sik + 423);
    const auto *sik_424 = buffer.data(sik + 424);
    const auto *sik_425 = buffer.data(sik + 425);
    const auto *sik_426 = buffer.data(sik + 426);
    const auto *sik_427 = buffer.data(sik + 427);
    const auto *sik_428 = buffer.data(sik + 428);
    const auto *sik_429 = buffer.data(sik + 429);
    const auto *sik_430 = buffer.data(sik + 430);
    const auto *sik_431 = buffer.data(sik + 431);
    const auto *sik_432 = buffer.data(sik + 432);
    const auto *sik_435 = buffer.data(sik + 435);
    const auto *sik_437 = buffer.data(sik + 437);
    const auto *sik_438 = buffer.data(sik + 438);
    const auto *sik_441 = buffer.data(sik + 441);
    const auto *sik_442 = buffer.data(sik + 442);
    const auto *sik_444 = buffer.data(sik + 444);
    const auto *sik_446 = buffer.data(sik + 446);
    const auto *sik_447 = buffer.data(sik + 447);
    const auto *sik_449 = buffer.data(sik + 449);
    const auto *sik_450 = buffer.data(sik + 450);
    const auto *sik_452 = buffer.data(sik + 452);
    const auto *sik_453 = buffer.data(sik + 453);

    const auto *sil1_270 = buffer.data(sil1 + 270);
    const auto *sil1_273 = buffer.data(sil1 + 273);
    const auto *sil1_276 = buffer.data(sil1 + 276);
    const auto *sil1_280 = buffer.data(sil1 + 280);
    const auto *sil1_282 = buffer.data(sil1 + 282);
    const auto *sil1_285 = buffer.data(sil1 + 285);
    const auto *sil1_287 = buffer.data(sil1 + 287);
    const auto *sil1_288 = buffer.data(sil1 + 288);
    const auto *sil1_291 = buffer.data(sil1 + 291);
    const auto *sil1_293 = buffer.data(sil1 + 293);
    const auto *sil1_294 = buffer.data(sil1 + 294);
    const auto *sil1_295 = buffer.data(sil1 + 295);
    const auto *sil1_306 = buffer.data(sil1 + 306);

    const auto *ski0_286 = buffer.data(ski0 + 286);
    const auto *ski0_289 = buffer.data(ski0 + 289);
    const auto *ski0_290 = buffer.data(ski0 + 290);
    const auto *ski0_292 = buffer.data(ski0 + 292);
    const auto *ski0_294 = buffer.data(ski0 + 294);
    const auto *ski0_295 = buffer.data(ski0 + 295);
    const auto *ski0_297 = buffer.data(ski0 + 297);
    const auto *ski0_298 = buffer.data(ski0 + 298);
    const auto *ski0_300 = buffer.data(ski0 + 300);
    const auto *ski0_301 = buffer.data(ski0 + 301);
    const auto *ski0_303 = buffer.data(ski0 + 303);
    const auto *ski0_304 = buffer.data(ski0 + 304);
    const auto *ski0_305 = buffer.data(ski0 + 305);
    const auto *ski0_306 = buffer.data(ski0 + 306);
    const auto *ski0_307 = buffer.data(ski0 + 307);
    const auto *ski0_313 = buffer.data(ski0 + 313);
    const auto *ski0_317 = buffer.data(ski0 + 317);
    const auto *ski0_322 = buffer.data(ski0 + 322);
    const auto *ski0_328 = buffer.data(ski0 + 328);
    const auto *ski0_331 = buffer.data(ski0 + 331);
    const auto *ski0_332 = buffer.data(ski0 + 332);
    const auto *ski0_333 = buffer.data(ski0 + 333);
    const auto *ski0_334 = buffer.data(ski0 + 334);
    const auto *ski0_335 = buffer.data(ski0 + 335);
    const auto *ski0_336 = buffer.data(ski0 + 336);
    const auto *ski0_339 = buffer.data(ski0 + 339);
    const auto *ski0_341 = buffer.data(ski0 + 341);
    const auto *ski0_342 = buffer.data(ski0 + 342);
    const auto *ski0_345 = buffer.data(ski0 + 345);
    const auto *ski0_346 = buffer.data(ski0 + 346);
    const auto *ski0_348 = buffer.data(ski0 + 348);
    const auto *ski0_350 = buffer.data(ski0 + 350);
    const auto *ski0_351 = buffer.data(ski0 + 351);
    const auto *ski0_353 = buffer.data(ski0 + 353);
    const auto *ski0_354 = buffer.data(ski0 + 354);
    const auto *ski0_356 = buffer.data(ski0 + 356);
    const auto *ski0_357 = buffer.data(ski0 + 357);

    const auto *ski1_286 = buffer.data(ski1 + 286);
    const auto *ski1_289 = buffer.data(ski1 + 289);
    const auto *ski1_290 = buffer.data(ski1 + 290);
    const auto *ski1_292 = buffer.data(ski1 + 292);
    const auto *ski1_294 = buffer.data(ski1 + 294);
    const auto *ski1_295 = buffer.data(ski1 + 295);
    const auto *ski1_297 = buffer.data(ski1 + 297);
    const auto *ski1_298 = buffer.data(ski1 + 298);
    const auto *ski1_300 = buffer.data(ski1 + 300);
    const auto *ski1_301 = buffer.data(ski1 + 301);
    const auto *ski1_303 = buffer.data(ski1 + 303);
    const auto *ski1_304 = buffer.data(ski1 + 304);
    const auto *ski1_305 = buffer.data(ski1 + 305);
    const auto *ski1_306 = buffer.data(ski1 + 306);
    const auto *ski1_307 = buffer.data(ski1 + 307);
    const auto *ski1_313 = buffer.data(ski1 + 313);
    const auto *ski1_317 = buffer.data(ski1 + 317);
    const auto *ski1_322 = buffer.data(ski1 + 322);
    const auto *ski1_328 = buffer.data(ski1 + 328);
    const auto *ski1_331 = buffer.data(ski1 + 331);
    const auto *ski1_332 = buffer.data(ski1 + 332);
    const auto *ski1_333 = buffer.data(ski1 + 333);
    const auto *ski1_334 = buffer.data(ski1 + 334);
    const auto *ski1_335 = buffer.data(ski1 + 335);
    const auto *ski1_336 = buffer.data(ski1 + 336);
    const auto *ski1_339 = buffer.data(ski1 + 339);
    const auto *ski1_341 = buffer.data(ski1 + 341);
    const auto *ski1_342 = buffer.data(ski1 + 342);
    const auto *ski1_345 = buffer.data(ski1 + 345);
    const auto *ski1_346 = buffer.data(ski1 + 346);
    const auto *ski1_348 = buffer.data(ski1 + 348);
    const auto *ski1_350 = buffer.data(ski1 + 350);
    const auto *ski1_351 = buffer.data(ski1 + 351);
    const auto *ski1_353 = buffer.data(ski1 + 353);
    const auto *ski1_354 = buffer.data(ski1 + 354);
    const auto *ski1_356 = buffer.data(ski1 + 356);
    const auto *ski1_357 = buffer.data(ski1 + 357);

    const auto *skk_363 = buffer.data(skk + 363);
    const auto *skk_365 = buffer.data(skk + 365);
    const auto *skk_366 = buffer.data(skk + 366);
    const auto *skk_369 = buffer.data(skk + 369);
    const auto *skk_370 = buffer.data(skk + 370);
    const auto *skk_372 = buffer.data(skk + 372);
    const auto *skk_374 = buffer.data(skk + 374);
    const auto *skk_375 = buffer.data(skk + 375);
    const auto *skk_377 = buffer.data(skk + 377);
    const auto *skk_378 = buffer.data(skk + 378);
    const auto *skk_380 = buffer.data(skk + 380);
    const auto *skk_381 = buffer.data(skk + 381);
    const auto *skk_383 = buffer.data(skk + 383);
    const auto *skk_384 = buffer.data(skk + 384);
    const auto *skk_385 = buffer.data(skk + 385);
    const auto *skk_387 = buffer.data(skk + 387);
    const auto *skk_388 = buffer.data(skk + 388);
    const auto *skk_389 = buffer.data(skk + 389);
    const auto *skk_390 = buffer.data(skk + 390);
    const auto *skk_391 = buffer.data(skk + 391);
    const auto *skk_392 = buffer.data(skk + 392);
    const auto *skk_393 = buffer.data(skk + 393);
    const auto *skk_394 = buffer.data(skk + 394);
    const auto *skk_395 = buffer.data(skk + 395);
    const auto *skk_396 = buffer.data(skk + 396);
    const auto *skk_398 = buffer.data(skk + 398);
    const auto *skk_399 = buffer.data(skk + 399);
    const auto *skk_401 = buffer.data(skk + 401);
    const auto *skk_402 = buffer.data(skk + 402);
    const auto *skk_405 = buffer.data(skk + 405);
    const auto *skk_406 = buffer.data(skk + 406);
    const auto *skk_410 = buffer.data(skk + 410);
    const auto *skk_411 = buffer.data(skk + 411);
    const auto *skk_416 = buffer.data(skk + 416);
    const auto *skk_423 = buffer.data(skk + 423);
    const auto *skk_424 = buffer.data(skk + 424);
    const auto *skk_425 = buffer.data(skk + 425);
    const auto *skk_426 = buffer.data(skk + 426);
    const auto *skk_427 = buffer.data(skk + 427);
    const auto *skk_428 = buffer.data(skk + 428);
    const auto *skk_429 = buffer.data(skk + 429);
    const auto *skk_430 = buffer.data(skk + 430);
    const auto *skk_431 = buffer.data(skk + 431);
    const auto *skk_432 = buffer.data(skk + 432);
    const auto *skk_434 = buffer.data(skk + 434);
    const auto *skk_435 = buffer.data(skk + 435);
    const auto *skk_437 = buffer.data(skk + 437);
    const auto *skk_438 = buffer.data(skk + 438);
    const auto *skk_441 = buffer.data(skk + 441);
    const auto *skk_442 = buffer.data(skk + 442);
    const auto *skk_444 = buffer.data(skk + 444);
    const auto *skk_446 = buffer.data(skk + 446);
    const auto *skk_447 = buffer.data(skk + 447);
    const auto *skk_449 = buffer.data(skk + 449);
    const auto *skk_450 = buffer.data(skk + 450);
    const auto *skk_452 = buffer.data(skk + 452);
    const auto *skk_453 = buffer.data(skk + 453);

#pragma omp simd aligned(t_456, t_457, t_458, pc_x, pc_y, pc_z, sik_221, sik_366, ski0_286, \
                         ski1_286, skk_363, skk_365, skk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_17 * sik_366[k]
                   + f_6 * ski0_286[k]
                   - f_7 * ski1_286[k]
                   + f_3 * pc_x[k] * skk_366[k];

        t_457[k] = f_3 * pc_z[k] * skk_363[k];

        t_458[k] = f_18 * sik_221[k]
                   + f_3 * pc_y[k] * skk_365[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_z, sik_369, sik_370, ski0_289, \
                         ski0_290, ski1_289, ski1_290, skk_366, skk_369, \
                         skk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_17 * sik_369[k]
                   + f_6 * ski0_289[k]
                   - f_7 * ski1_289[k]
                   + f_3 * pc_x[k] * skk_369[k];

        t_460[k] = f_17 * sik_370[k]
                   + f_8 * ski0_290[k]
                   - f_9 * ski1_290[k]
                   + f_3 * pc_x[k] * skk_370[k];

        t_461[k] = f_3 * pc_z[k] * skk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_x, pc_y, sik_225, sik_372, sik_374, ski0_292, \
                         ski0_294, ski1_292, ski1_294, skk_369, skk_372, \
                         skk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_17 * sik_372[k]
                   + f_8 * ski0_292[k]
                   - f_9 * ski1_292[k]
                   + f_3 * pc_x[k] * skk_372[k];

        t_463[k] = f_18 * sik_225[k]
                   + f_3 * pc_y[k] * skk_369[k];

        t_464[k] = f_17 * sik_374[k]
                   + f_8 * ski0_294[k]
                   - f_9 * ski1_294[k]
                   + f_3 * pc_x[k] * skk_374[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, sik_375, sik_377, ski0_295, \
                         ski0_297, ski1_295, ski1_297, skk_370, skk_375, \
                         skk_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_17 * sik_375[k]
                   + f_10 * ski0_295[k]
                   - f_11 * ski1_295[k]
                   + f_3 * pc_x[k] * skk_375[k];

        t_466[k] = f_3 * pc_z[k] * skk_370[k];

        t_467[k] = f_17 * sik_377[k]
                   + f_10 * ski0_297[k]
                   - f_11 * ski1_297[k]
                   + f_3 * pc_x[k] * skk_377[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, sik_230, sik_378, sik_380, ski0_298, \
                         ski0_300, ski1_298, ski1_300, skk_374, skk_378, \
                         skk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_17 * sik_378[k]
                   + f_10 * ski0_298[k]
                   - f_11 * ski1_298[k]
                   + f_3 * pc_x[k] * skk_378[k];

        t_469[k] = f_18 * sik_230[k]
                   + f_3 * pc_y[k] * skk_374[k];

        t_470[k] = f_17 * sik_380[k]
                   + f_10 * ski0_300[k]
                   - f_11 * ski1_300[k]
                   + f_3 * pc_x[k] * skk_380[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, sik_381, sik_383, ski0_301, \
                         ski0_303, ski1_301, ski1_303, skk_375, skk_381, \
                         skk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_17 * sik_381[k]
                   + f_12 * ski0_301[k]
                   - f_13 * ski1_301[k]
                   + f_3 * pc_x[k] * skk_381[k];

        t_472[k] = f_3 * pc_z[k] * skk_375[k];

        t_473[k] = f_17 * sik_383[k]
                   + f_12 * ski0_303[k]
                   - f_13 * ski1_303[k]
                   + f_3 * pc_x[k] * skk_383[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, sik_236, sik_384, sik_385, ski0_304, \
                         ski0_305, ski1_304, ski1_305, skk_380, skk_384, \
                         skk_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_17 * sik_384[k]
                   + f_12 * ski0_304[k]
                   - f_13 * ski1_304[k]
                   + f_3 * pc_x[k] * skk_384[k];

        t_475[k] = f_17 * sik_385[k]
                   + f_12 * ski0_305[k]
                   - f_13 * ski1_305[k]
                   + f_3 * pc_x[k] * skk_385[k];

        t_476[k] = f_18 * sik_236[k]
                   + f_3 * pc_y[k] * skk_380[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pc_x, sik_387, sik_388, sik_389, sik_390, \
                         ski0_307, ski1_307, skk_387, skk_388, skk_389, \
                         skk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_17 * sik_387[k]
                   + f_12 * ski0_307[k]
                   - f_13 * ski1_307[k]
                   + f_3 * pc_x[k] * skk_387[k];

        t_478[k] = f_17 * sik_388[k]
                   + f_3 * pc_x[k] * skk_388[k];

        t_479[k] = f_17 * sik_389[k]
                   + f_3 * pc_x[k] * skk_389[k];

        t_480[k] = f_17 * sik_390[k]
                   + f_3 * pc_x[k] * skk_390[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pc_x, sik_391, sik_392, sik_393, \
                         sik_394, sik_395, skk_391, skk_392, skk_393, skk_394, \
                         skk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_17 * sik_391[k]
                   + f_3 * pc_x[k] * skk_391[k];

        t_482[k] = f_17 * sik_392[k]
                   + f_3 * pc_x[k] * skk_392[k];

        t_483[k] = f_17 * sik_393[k]
                   + f_3 * pc_x[k] * skk_393[k];

        t_484[k] = f_17 * sik_394[k]
                   + f_3 * pc_x[k] * skk_394[k];

        t_485[k] = f_17 * sik_395[k]
                   + f_3 * pc_x[k] * skk_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_y, pc_z, sik_244, sik_246, ski0_301, \
                         ski0_303, ski1_301, ski1_303, skk_388, \
                         skk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_18 * sik_244[k]
                   + f_1 * ski0_301[k]
                   - f_2 * ski1_301[k]
                   + f_3 * pc_y[k] * skk_388[k];

        t_487[k] = f_3 * pc_z[k] * skk_388[k];

        t_488[k] = f_18 * sik_246[k]
                   + f_4 * ski0_303[k]
                   - f_5 * ski1_303[k]
                   + f_3 * pc_y[k] * skk_390[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_y, sik_247, sik_248, sik_249, ski0_304, \
                         ski0_305, ski0_306, ski1_304, ski1_305, ski1_306, skk_391, skk_392, \
                         skk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_18 * sik_247[k]
                   + f_6 * ski0_304[k]
                   - f_7 * ski1_304[k]
                   + f_3 * pc_y[k] * skk_391[k];

        t_490[k] = f_18 * sik_248[k]
                   + f_8 * ski0_305[k]
                   - f_9 * ski1_305[k]
                   + f_3 * pc_y[k] * skk_392[k];

        t_491[k] = f_18 * sik_249[k]
                   + f_10 * ski0_306[k]
                   - f_11 * ski1_306[k]
                   + f_3 * pc_y[k] * skk_393[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_z, pc_y, pc_z, sil0_270, sik_250, \
                         sik_251, sil1_270, ski0_307, ski1_307, skk_394, \
                         skk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_18 * sik_250[k]
                   + f_12 * ski0_307[k]
                   - f_13 * ski1_307[k]
                   + f_3 * pc_y[k] * skk_394[k];

        t_493[k] = f_18 * sik_251[k]
                   + f_3 * pc_y[k] * skk_395[k];

        t_494[k] = f_1 * ski0_307[k]
                   - f_2 * ski1_307[k]
                   + f_3 * pc_z[k] * skk_395[k];

        t_495[k] = pb_z[k] * sil0_270[k]
                   - f_14 * pc_z[k] * sil1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, sil0_273, sik_216, \
                         sik_252, sik_254, sil1_273, skk_396, skk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * sik_252[k]
                   + f_3 * pc_y[k] * skk_396[k];

        t_497[k] = f_15 * sik_216[k]
                   + f_3 * pc_z[k] * skk_396[k];

        t_498[k] = pb_z[k] * sil0_273[k]
                   - f_14 * pc_z[k] * sil1_273[k];

        t_499[k] = f_17 * sik_254[k]
                   + f_3 * pc_y[k] * skk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pb_z, pc_x, pc_z, sil0_276, sik_219, sik_401, \
                         sil1_276, ski0_313, ski1_313, skk_399, \
                         skk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_17 * sik_401[k]
                   + f_4 * ski0_313[k]
                   - f_5 * ski1_313[k]
                   + f_3 * pc_x[k] * skk_401[k];

        t_501[k] = pb_z[k] * sil0_276[k]
                   - f_14 * pc_z[k] * sil1_276[k];

        t_502[k] = f_15 * sik_219[k]
                   + f_3 * pc_z[k] * skk_399[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pb_z, pc_x, pc_y, pc_z, sil0_280, sik_257, \
                         sik_405, sil1_280, ski0_317, ski1_317, skk_401, \
                         skk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * sik_257[k]
                   + f_3 * pc_y[k] * skk_401[k];

        t_504[k] = f_17 * sik_405[k]
                   + f_6 * ski0_317[k]
                   - f_7 * ski1_317[k]
                   + f_3 * pc_x[k] * skk_405[k];

        t_505[k] = pb_z[k] * sil0_280[k]
                   - f_14 * pc_z[k] * sil1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_z, pc_y, pc_z, sil0_282, sik_222, sik_223, \
                         sik_261, sil1_282, skk_402, skk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_15 * sik_222[k]
                   + f_3 * pc_z[k] * skk_402[k];

        t_507[k] = pb_z[k] * sil0_282[k]
                   + f_16 * sik_223[k]
                   - f_14 * pc_z[k] * sil1_282[k];

        t_508[k] = f_17 * sik_261[k]
                   + f_3 * pc_y[k] * skk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_z, pc_x, pc_z, sil0_285, sik_226, sik_410, \
                         sil1_285, ski0_322, ski1_322, skk_406, \
                         skk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_17 * sik_410[k]
                   + f_8 * ski0_322[k]
                   - f_9 * ski1_322[k]
                   + f_3 * pc_x[k] * skk_410[k];

        t_510[k] = pb_z[k] * sil0_285[k]
                   - f_14 * pc_z[k] * sil1_285[k];

        t_511[k] = f_15 * sik_226[k]
                   + f_3 * pc_z[k] * skk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_z, pc_y, pc_z, sil0_287, sil0_288, sik_227, \
                         sik_228, sik_266, sil1_287, sil1_288, \
                         skk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_z[k] * sil0_287[k]
                   + f_16 * sik_227[k]
                   - f_14 * pc_z[k] * sil1_287[k];

        t_513[k] = pb_z[k] * sil0_288[k]
                   + f_17 * sik_228[k]
                   - f_14 * pc_z[k] * sil1_288[k];

        t_514[k] = f_17 * sik_266[k]
                   + f_3 * pc_y[k] * skk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_z, pc_x, pc_z, sil0_291, sik_231, sik_416, \
                         sil1_291, ski0_328, ski1_328, skk_411, \
                         skk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_17 * sik_416[k]
                   + f_10 * ski0_328[k]
                   - f_11 * ski1_328[k]
                   + f_3 * pc_x[k] * skk_416[k];

        t_516[k] = pb_z[k] * sil0_291[k]
                   - f_14 * pc_z[k] * sil1_291[k];

        t_517[k] = f_15 * sik_231[k]
                   + f_3 * pc_z[k] * skk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pb_z, pc_z, sil0_293, sil0_294, sil0_295, \
                         sik_232, sik_233, sik_234, sil1_293, sil1_294, \
                         sil1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_z[k] * sil0_293[k]
                   + f_16 * sik_232[k]
                   - f_14 * pc_z[k] * sil1_293[k];

        t_519[k] = pb_z[k] * sil0_294[k]
                   + f_17 * sik_233[k]
                   - f_14 * pc_z[k] * sil1_294[k];

        t_520[k] = pb_z[k] * sil0_295[k]
                   + f_18 * sik_234[k]
                   - f_14 * pc_z[k] * sil1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, sik_272, sik_423, sik_424, \
                         sik_425, ski0_335, ski1_335, skk_416, skk_423, skk_424, \
                         skk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * sik_272[k]
                   + f_3 * pc_y[k] * skk_416[k];

        t_522[k] = f_17 * sik_423[k]
                   + f_12 * ski0_335[k]
                   - f_13 * ski1_335[k]
                   + f_3 * pc_x[k] * skk_423[k];

        t_523[k] = f_17 * sik_424[k]
                   + f_3 * pc_x[k] * skk_424[k];

        t_524[k] = f_17 * sik_425[k]
                   + f_3 * pc_x[k] * skk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, sik_426, sik_427, sik_428, \
                         sik_429, sik_430, skk_426, skk_427, skk_428, skk_429, \
                         skk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_17 * sik_426[k]
                   + f_3 * pc_x[k] * skk_426[k];

        t_526[k] = f_17 * sik_427[k]
                   + f_3 * pc_x[k] * skk_427[k];

        t_527[k] = f_17 * sik_428[k]
                   + f_3 * pc_x[k] * skk_428[k];

        t_528[k] = f_17 * sik_429[k]
                   + f_3 * pc_x[k] * skk_429[k];

        t_529[k] = f_17 * sik_430[k]
                   + f_3 * pc_x[k] * skk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pb_z, pc_x, pc_z, sil0_306, sik_244, sik_431, \
                         sil1_306, skk_424, skk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_17 * sik_431[k]
                   + f_3 * pc_x[k] * skk_431[k];

        t_531[k] = pb_z[k] * sil0_306[k]
                   - f_14 * pc_z[k] * sil1_306[k];

        t_532[k] = f_15 * sik_244[k]
                   + f_3 * pc_z[k] * skk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, sik_282, sik_283, sik_284, ski0_331, \
                         ski0_332, ski0_333, ski1_331, ski1_332, ski1_333, skk_426, skk_427, \
                         skk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * sik_282[k]
                   + f_4 * ski0_331[k]
                   - f_5 * ski1_331[k]
                   + f_3 * pc_y[k] * skk_426[k];

        t_534[k] = f_17 * sik_283[k]
                   + f_6 * ski0_332[k]
                   - f_7 * ski1_332[k]
                   + f_3 * pc_y[k] * skk_427[k];

        t_535[k] = f_17 * sik_284[k]
                   + f_8 * ski0_333[k]
                   - f_9 * ski1_333[k]
                   + f_3 * pc_y[k] * skk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, sik_285, sik_286, sik_287, ski0_334, \
                         ski0_335, ski1_334, ski1_335, skk_429, skk_430, \
                         skk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * sik_285[k]
                   + f_10 * ski0_334[k]
                   - f_11 * ski1_334[k]
                   + f_3 * pc_y[k] * skk_429[k];

        t_537[k] = f_17 * sik_286[k]
                   + f_12 * ski0_335[k]
                   - f_13 * ski1_335[k]
                   + f_3 * pc_y[k] * skk_430[k];

        t_538[k] = f_17 * sik_287[k]
                   + f_3 * pc_y[k] * skk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, sik_251, sik_288, sik_432, \
                         ski0_335, ski0_336, ski1_335, ski1_336, skk_431, \
                         skk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * sik_251[k]
                   + f_1 * ski0_335[k]
                   - f_2 * ski1_335[k]
                   + f_3 * pc_z[k] * skk_431[k];

        t_540[k] = f_17 * sik_432[k]
                   + f_1 * ski0_336[k]
                   - f_2 * ski1_336[k]
                   + f_3 * pc_x[k] * skk_432[k];

        t_541[k] = f_16 * sik_288[k]
                   + f_3 * pc_y[k] * skk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, sik_252, sik_290, sik_435, \
                         ski0_339, ski1_339, skk_432, skk_434, \
                         skk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * sik_252[k]
                   + f_3 * pc_z[k] * skk_432[k];

        t_543[k] = f_17 * sik_435[k]
                   + f_4 * ski0_339[k]
                   - f_5 * ski1_339[k]
                   + f_3 * pc_x[k] * skk_435[k];

        t_544[k] = f_16 * sik_290[k]
                   + f_3 * pc_y[k] * skk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, sik_255, sik_437, sik_438, ski0_341, \
                         ski0_342, ski1_341, ski1_342, skk_435, skk_437, \
                         skk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_17 * sik_437[k]
                   + f_4 * ski0_341[k]
                   - f_5 * ski1_341[k]
                   + f_3 * pc_x[k] * skk_437[k];

        t_546[k] = f_17 * sik_438[k]
                   + f_6 * ski0_342[k]
                   - f_7 * ski1_342[k]
                   + f_3 * pc_x[k] * skk_438[k];

        t_547[k] = f_16 * sik_255[k]
                   + f_3 * pc_z[k] * skk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, sik_293, sik_441, sik_442, ski0_345, \
                         ski0_346, ski1_345, ski1_346, skk_437, skk_441, \
                         skk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * sik_293[k]
                   + f_3 * pc_y[k] * skk_437[k];

        t_549[k] = f_17 * sik_441[k]
                   + f_6 * ski0_345[k]
                   - f_7 * ski1_345[k]
                   + f_3 * pc_x[k] * skk_441[k];

        t_550[k] = f_17 * sik_442[k]
                   + f_8 * ski0_346[k]
                   - f_9 * ski1_346[k]
                   + f_3 * pc_x[k] * skk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, sik_258, sik_297, sik_444, \
                         ski0_348, ski1_348, skk_438, skk_441, \
                         skk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * sik_258[k]
                   + f_3 * pc_z[k] * skk_438[k];

        t_552[k] = f_17 * sik_444[k]
                   + f_8 * ski0_348[k]
                   - f_9 * ski1_348[k]
                   + f_3 * pc_x[k] * skk_444[k];

        t_553[k] = f_16 * sik_297[k]
                   + f_3 * pc_y[k] * skk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, sik_262, sik_446, sik_447, ski0_350, \
                         ski0_351, ski1_350, ski1_351, skk_442, skk_446, \
                         skk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_17 * sik_446[k]
                   + f_8 * ski0_350[k]
                   - f_9 * ski1_350[k]
                   + f_3 * pc_x[k] * skk_446[k];

        t_555[k] = f_17 * sik_447[k]
                   + f_10 * ski0_351[k]
                   - f_11 * ski1_351[k]
                   + f_3 * pc_x[k] * skk_447[k];

        t_556[k] = f_16 * sik_262[k]
                   + f_3 * pc_z[k] * skk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, sik_302, sik_449, sik_450, ski0_353, \
                         ski0_354, ski1_353, ski1_354, skk_446, skk_449, \
                         skk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_17 * sik_449[k]
                   + f_10 * ski0_353[k]
                   - f_11 * ski1_353[k]
                   + f_3 * pc_x[k] * skk_449[k];

        t_558[k] = f_17 * sik_450[k]
                   + f_10 * ski0_354[k]
                   - f_11 * ski1_354[k]
                   + f_3 * pc_x[k] * skk_450[k];

        t_559[k] = f_16 * sik_302[k]
                   + f_3 * pc_y[k] * skk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, sik_267, sik_452, sik_453, ski0_356, \
                         ski0_357, ski1_356, ski1_357, skk_447, skk_452, \
                         skk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_17 * sik_452[k]
                   + f_10 * ski0_356[k]
                   - f_11 * ski1_356[k]
                   + f_3 * pc_x[k] * skk_452[k];

        t_561[k] = f_17 * sik_453[k]
                   + f_12 * ski0_357[k]
                   - f_13 * ski1_357[k]
                   + f_3 * pc_x[k] * skk_453[k];

        t_562[k] = f_16 * sik_267[k]
                   + f_3 * pc_z[k] * skk_447[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_405 = buffer.data(sil0 + 405);
    const auto *sil0_408 = buffer.data(sil0 + 408);
    const auto *sil0_410 = buffer.data(sil0 + 410);
    const auto *sil0_411 = buffer.data(sil0 + 411);
    const auto *sil0_414 = buffer.data(sil0 + 414);
    const auto *sil0_415 = buffer.data(sil0 + 415);
    const auto *sil0_417 = buffer.data(sil0 + 417);
    const auto *sil0_419 = buffer.data(sil0 + 419);
    const auto *sil0_420 = buffer.data(sil0 + 420);
    const auto *sil0_422 = buffer.data(sil0 + 422);
    const auto *sil0_423 = buffer.data(sil0 + 423);
    const auto *sil0_425 = buffer.data(sil0 + 425);
    const auto *sil0_426 = buffer.data(sil0 + 426);
    const auto *sil0_428 = buffer.data(sil0 + 428);
    const auto *sil0_429 = buffer.data(sil0 + 429);
    const auto *sil0_430 = buffer.data(sil0 + 430);
    const auto *sil0_432 = buffer.data(sil0 + 432);
    const auto *sil0_449 = buffer.data(sil0 + 449);

    const auto *sik_280 = buffer.data(sik + 280);
    const auto *sik_287 = buffer.data(sik + 287);
    const auto *sik_288 = buffer.data(sik + 288);
    const auto *sik_291 = buffer.data(sik + 291);
    const auto *sik_294 = buffer.data(sik + 294);
    const auto *sik_298 = buffer.data(sik + 298);
    const auto *sik_303 = buffer.data(sik + 303);
    const auto *sik_308 = buffer.data(sik + 308);
    const auto *sik_316 = buffer.data(sik + 316);
    const auto *sik_318 = buffer.data(sik + 318);
    const auto *sik_319 = buffer.data(sik + 319);
    const auto *sik_320 = buffer.data(sik + 320);
    const auto *sik_321 = buffer.data(sik + 321);
    const auto *sik_322 = buffer.data(sik + 322);
    const auto *sik_323 = buffer.data(sik + 323);
    const auto *sik_324 = buffer.data(sik + 324);
    const auto *sik_325 = buffer.data(sik + 325);
    const auto *sik_326 = buffer.data(sik + 326);
    const auto *sik_327 = buffer.data(sik + 327);
    const auto *sik_329 = buffer.data(sik + 329);
    const auto *sik_330 = buffer.data(sik + 330);
    const auto *sik_332 = buffer.data(sik + 332);
    const auto *sik_333 = buffer.data(sik + 333);
    const auto *sik_334 = buffer.data(sik + 334);
    const auto *sik_336 = buffer.data(sik + 336);
    const auto *sik_337 = buffer.data(sik + 337);
    const auto *sik_338 = buffer.data(sik + 338);
    const auto *sik_339 = buffer.data(sik + 339);
    const auto *sik_341 = buffer.data(sik + 341);
    const auto *sik_342 = buffer.data(sik + 342);
    const auto *sik_343 = buffer.data(sik + 343);
    const auto *sik_344 = buffer.data(sik + 344);
    const auto *sik_352 = buffer.data(sik + 352);
    const auto *sik_354 = buffer.data(sik + 354);
    const auto *sik_355 = buffer.data(sik + 355);
    const auto *sik_356 = buffer.data(sik + 356);
    const auto *sik_357 = buffer.data(sik + 357);
    const auto *sik_358 = buffer.data(sik + 358);
    const auto *sik_359 = buffer.data(sik + 359);
    const auto *sik_455 = buffer.data(sik + 455);
    const auto *sik_456 = buffer.data(sik + 456);
    const auto *sik_457 = buffer.data(sik + 457);
    const auto *sik_459 = buffer.data(sik + 459);
    const auto *sik_460 = buffer.data(sik + 460);
    const auto *sik_461 = buffer.data(sik + 461);
    const auto *sik_462 = buffer.data(sik + 462);
    const auto *sik_463 = buffer.data(sik + 463);
    const auto *sik_464 = buffer.data(sik + 464);
    const auto *sik_465 = buffer.data(sik + 465);
    const auto *sik_466 = buffer.data(sik + 466);
    const auto *sik_467 = buffer.data(sik + 467);
    const auto *sik_496 = buffer.data(sik + 496);
    const auto *sik_497 = buffer.data(sik + 497);
    const auto *sik_498 = buffer.data(sik + 498);
    const auto *sik_499 = buffer.data(sik + 499);
    const auto *sik_500 = buffer.data(sik + 500);
    const auto *sik_501 = buffer.data(sik + 501);
    const auto *sik_502 = buffer.data(sik + 502);
    const auto *sik_503 = buffer.data(sik + 503);
    const auto *sik_504 = buffer.data(sik + 504);
    const auto *sik_507 = buffer.data(sik + 507);
    const auto *sik_509 = buffer.data(sik + 509);
    const auto *sik_510 = buffer.data(sik + 510);
    const auto *sik_513 = buffer.data(sik + 513);
    const auto *sik_514 = buffer.data(sik + 514);
    const auto *sik_516 = buffer.data(sik + 516);
    const auto *sik_518 = buffer.data(sik + 518);
    const auto *sik_519 = buffer.data(sik + 519);
    const auto *sik_521 = buffer.data(sik + 521);
    const auto *sik_522 = buffer.data(sik + 522);
    const auto *sik_524 = buffer.data(sik + 524);
    const auto *sik_525 = buffer.data(sik + 525);
    const auto *sik_527 = buffer.data(sik + 527);
    const auto *sik_528 = buffer.data(sik + 528);
    const auto *sik_529 = buffer.data(sik + 529);
    const auto *sik_531 = buffer.data(sik + 531);
    const auto *sik_532 = buffer.data(sik + 532);
    const auto *sik_533 = buffer.data(sik + 533);
    const auto *sik_534 = buffer.data(sik + 534);
    const auto *sik_535 = buffer.data(sik + 535);
    const auto *sik_536 = buffer.data(sik + 536);
    const auto *sik_537 = buffer.data(sik + 537);
    const auto *sik_538 = buffer.data(sik + 538);
    const auto *sik_539 = buffer.data(sik + 539);

    const auto *sil1_405 = buffer.data(sil1 + 405);
    const auto *sil1_408 = buffer.data(sil1 + 408);
    const auto *sil1_410 = buffer.data(sil1 + 410);
    const auto *sil1_411 = buffer.data(sil1 + 411);
    const auto *sil1_414 = buffer.data(sil1 + 414);
    const auto *sil1_415 = buffer.data(sil1 + 415);
    const auto *sil1_417 = buffer.data(sil1 + 417);
    const auto *sil1_419 = buffer.data(sil1 + 419);
    const auto *sil1_420 = buffer.data(sil1 + 420);
    const auto *sil1_422 = buffer.data(sil1 + 422);
    const auto *sil1_423 = buffer.data(sil1 + 423);
    const auto *sil1_425 = buffer.data(sil1 + 425);
    const auto *sil1_426 = buffer.data(sil1 + 426);
    const auto *sil1_428 = buffer.data(sil1 + 428);
    const auto *sil1_429 = buffer.data(sil1 + 429);
    const auto *sil1_430 = buffer.data(sil1 + 430);
    const auto *sil1_432 = buffer.data(sil1 + 432);
    const auto *sil1_449 = buffer.data(sil1 + 449);

    const auto *ski0_357 = buffer.data(ski0 + 357);
    const auto *ski0_359 = buffer.data(ski0 + 359);
    const auto *ski0_360 = buffer.data(ski0 + 360);
    const auto *ski0_361 = buffer.data(ski0 + 361);
    const auto *ski0_362 = buffer.data(ski0 + 362);
    const auto *ski0_363 = buffer.data(ski0 + 363);
    const auto *ski0_385 = buffer.data(ski0 + 385);
    const auto *ski0_387 = buffer.data(ski0 + 387);
    const auto *ski0_388 = buffer.data(ski0 + 388);
    const auto *ski0_389 = buffer.data(ski0 + 389);
    const auto *ski0_390 = buffer.data(ski0 + 390);
    const auto *ski0_391 = buffer.data(ski0 + 391);
    const auto *ski0_392 = buffer.data(ski0 + 392);
    const auto *ski0_395 = buffer.data(ski0 + 395);
    const auto *ski0_397 = buffer.data(ski0 + 397);
    const auto *ski0_398 = buffer.data(ski0 + 398);
    const auto *ski0_401 = buffer.data(ski0 + 401);
    const auto *ski0_402 = buffer.data(ski0 + 402);
    const auto *ski0_404 = buffer.data(ski0 + 404);
    const auto *ski0_406 = buffer.data(ski0 + 406);
    const auto *ski0_407 = buffer.data(ski0 + 407);
    const auto *ski0_409 = buffer.data(ski0 + 409);
    const auto *ski0_410 = buffer.data(ski0 + 410);
    const auto *ski0_412 = buffer.data(ski0 + 412);
    const auto *ski0_413 = buffer.data(ski0 + 413);
    const auto *ski0_415 = buffer.data(ski0 + 415);
    const auto *ski0_416 = buffer.data(ski0 + 416);
    const auto *ski0_417 = buffer.data(ski0 + 417);
    const auto *ski0_418 = buffer.data(ski0 + 418);
    const auto *ski0_419 = buffer.data(ski0 + 419);

    const auto *ski1_357 = buffer.data(ski1 + 357);
    const auto *ski1_359 = buffer.data(ski1 + 359);
    const auto *ski1_360 = buffer.data(ski1 + 360);
    const auto *ski1_361 = buffer.data(ski1 + 361);
    const auto *ski1_362 = buffer.data(ski1 + 362);
    const auto *ski1_363 = buffer.data(ski1 + 363);
    const auto *ski1_385 = buffer.data(ski1 + 385);
    const auto *ski1_387 = buffer.data(ski1 + 387);
    const auto *ski1_388 = buffer.data(ski1 + 388);
    const auto *ski1_389 = buffer.data(ski1 + 389);
    const auto *ski1_390 = buffer.data(ski1 + 390);
    const auto *ski1_391 = buffer.data(ski1 + 391);
    const auto *ski1_392 = buffer.data(ski1 + 392);
    const auto *ski1_395 = buffer.data(ski1 + 395);
    const auto *ski1_397 = buffer.data(ski1 + 397);
    const auto *ski1_398 = buffer.data(ski1 + 398);
    const auto *ski1_401 = buffer.data(ski1 + 401);
    const auto *ski1_402 = buffer.data(ski1 + 402);
    const auto *ski1_404 = buffer.data(ski1 + 404);
    const auto *ski1_406 = buffer.data(ski1 + 406);
    const auto *ski1_407 = buffer.data(ski1 + 407);
    const auto *ski1_409 = buffer.data(ski1 + 409);
    const auto *ski1_410 = buffer.data(ski1 + 410);
    const auto *ski1_412 = buffer.data(ski1 + 412);
    const auto *ski1_413 = buffer.data(ski1 + 413);
    const auto *ski1_415 = buffer.data(ski1 + 415);
    const auto *ski1_416 = buffer.data(ski1 + 416);
    const auto *ski1_417 = buffer.data(ski1 + 417);
    const auto *ski1_418 = buffer.data(ski1 + 418);
    const auto *ski1_419 = buffer.data(ski1 + 419);

    const auto *skk_452 = buffer.data(skk + 452);
    const auto *skk_455 = buffer.data(skk + 455);
    const auto *skk_456 = buffer.data(skk + 456);
    const auto *skk_457 = buffer.data(skk + 457);
    const auto *skk_459 = buffer.data(skk + 459);
    const auto *skk_460 = buffer.data(skk + 460);
    const auto *skk_461 = buffer.data(skk + 461);
    const auto *skk_462 = buffer.data(skk + 462);
    const auto *skk_463 = buffer.data(skk + 463);
    const auto *skk_464 = buffer.data(skk + 464);
    const auto *skk_465 = buffer.data(skk + 465);
    const auto *skk_466 = buffer.data(skk + 466);
    const auto *skk_467 = buffer.data(skk + 467);
    const auto *skk_468 = buffer.data(skk + 468);
    const auto *skk_470 = buffer.data(skk + 470);
    const auto *skk_471 = buffer.data(skk + 471);
    const auto *skk_473 = buffer.data(skk + 473);
    const auto *skk_474 = buffer.data(skk + 474);
    const auto *skk_477 = buffer.data(skk + 477);
    const auto *skk_478 = buffer.data(skk + 478);
    const auto *skk_482 = buffer.data(skk + 482);
    const auto *skk_483 = buffer.data(skk + 483);
    const auto *skk_488 = buffer.data(skk + 488);
    const auto *skk_496 = buffer.data(skk + 496);
    const auto *skk_497 = buffer.data(skk + 497);
    const auto *skk_498 = buffer.data(skk + 498);
    const auto *skk_499 = buffer.data(skk + 499);
    const auto *skk_500 = buffer.data(skk + 500);
    const auto *skk_501 = buffer.data(skk + 501);
    const auto *skk_502 = buffer.data(skk + 502);
    const auto *skk_503 = buffer.data(skk + 503);
    const auto *skk_504 = buffer.data(skk + 504);
    const auto *skk_506 = buffer.data(skk + 506);
    const auto *skk_507 = buffer.data(skk + 507);
    const auto *skk_509 = buffer.data(skk + 509);
    const auto *skk_510 = buffer.data(skk + 510);
    const auto *skk_513 = buffer.data(skk + 513);
    const auto *skk_514 = buffer.data(skk + 514);
    const auto *skk_516 = buffer.data(skk + 516);
    const auto *skk_518 = buffer.data(skk + 518);
    const auto *skk_519 = buffer.data(skk + 519);
    const auto *skk_521 = buffer.data(skk + 521);
    const auto *skk_522 = buffer.data(skk + 522);
    const auto *skk_524 = buffer.data(skk + 524);
    const auto *skk_525 = buffer.data(skk + 525);
    const auto *skk_527 = buffer.data(skk + 527);
    const auto *skk_528 = buffer.data(skk + 528);
    const auto *skk_529 = buffer.data(skk + 529);
    const auto *skk_531 = buffer.data(skk + 531);
    const auto *skk_532 = buffer.data(skk + 532);
    const auto *skk_533 = buffer.data(skk + 533);
    const auto *skk_534 = buffer.data(skk + 534);
    const auto *skk_535 = buffer.data(skk + 535);
    const auto *skk_536 = buffer.data(skk + 536);
    const auto *skk_537 = buffer.data(skk + 537);
    const auto *skk_538 = buffer.data(skk + 538);
    const auto *skk_539 = buffer.data(skk + 539);

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, sik_455, sik_456, sik_457, ski0_359, \
                         ski0_360, ski0_361, ski1_359, ski1_360, ski1_361, skk_455, skk_456, \
                         skk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_17 * sik_455[k]
                   + f_12 * ski0_359[k]
                   - f_13 * ski1_359[k]
                   + f_3 * pc_x[k] * skk_455[k];

        t_564[k] = f_17 * sik_456[k]
                   + f_12 * ski0_360[k]
                   - f_13 * ski1_360[k]
                   + f_3 * pc_x[k] * skk_456[k];

        t_565[k] = f_17 * sik_457[k]
                   + f_12 * ski0_361[k]
                   - f_13 * ski1_361[k]
                   + f_3 * pc_x[k] * skk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, sik_308, sik_459, sik_460, \
                         sik_461, ski0_363, ski1_363, skk_452, skk_459, skk_460, \
                         skk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * sik_308[k]
                   + f_3 * pc_y[k] * skk_452[k];

        t_567[k] = f_17 * sik_459[k]
                   + f_12 * ski0_363[k]
                   - f_13 * ski1_363[k]
                   + f_3 * pc_x[k] * skk_459[k];

        t_568[k] = f_17 * sik_460[k]
                   + f_3 * pc_x[k] * skk_460[k];

        t_569[k] = f_17 * sik_461[k]
                   + f_3 * pc_x[k] * skk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, sik_462, sik_463, sik_464, \
                         sik_465, sik_466, skk_462, skk_463, skk_464, skk_465, \
                         skk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_17 * sik_462[k]
                   + f_3 * pc_x[k] * skk_462[k];

        t_571[k] = f_17 * sik_463[k]
                   + f_3 * pc_x[k] * skk_463[k];

        t_572[k] = f_17 * sik_464[k]
                   + f_3 * pc_x[k] * skk_464[k];

        t_573[k] = f_17 * sik_465[k]
                   + f_3 * pc_x[k] * skk_465[k];

        t_574[k] = f_17 * sik_466[k]
                   + f_3 * pc_x[k] * skk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, sik_280, sik_316, sik_467, \
                         ski0_357, ski1_357, skk_460, skk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_17 * sik_467[k]
                   + f_3 * pc_x[k] * skk_467[k];

        t_576[k] = f_16 * sik_316[k]
                   + f_1 * ski0_357[k]
                   - f_2 * ski1_357[k]
                   + f_3 * pc_y[k] * skk_460[k];

        t_577[k] = f_16 * sik_280[k]
                   + f_3 * pc_z[k] * skk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, sik_318, sik_319, sik_320, ski0_359, \
                         ski0_360, ski0_361, ski1_359, ski1_360, ski1_361, skk_462, skk_463, \
                         skk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * sik_318[k]
                   + f_4 * ski0_359[k]
                   - f_5 * ski1_359[k]
                   + f_3 * pc_y[k] * skk_462[k];

        t_579[k] = f_16 * sik_319[k]
                   + f_6 * ski0_360[k]
                   - f_7 * ski1_360[k]
                   + f_3 * pc_y[k] * skk_463[k];

        t_580[k] = f_16 * sik_320[k]
                   + f_8 * ski0_361[k]
                   - f_9 * ski1_361[k]
                   + f_3 * pc_y[k] * skk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, sik_321, sik_322, sik_323, ski0_362, \
                         ski0_363, ski1_362, ski1_363, skk_465, skk_466, \
                         skk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * sik_321[k]
                   + f_10 * ski0_362[k]
                   - f_11 * ski1_362[k]
                   + f_3 * pc_y[k] * skk_465[k];

        t_582[k] = f_16 * sik_322[k]
                   + f_12 * ski0_363[k]
                   - f_13 * ski1_363[k]
                   + f_3 * pc_y[k] * skk_466[k];

        t_583[k] = f_16 * sik_323[k]
                   + f_3 * pc_y[k] * skk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pc_y, pc_z, sil0_405, sik_287, \
                         sik_288, sik_324, sil1_405, ski0_363, ski1_363, skk_467, \
                         skk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * sik_287[k]
                   + f_1 * ski0_363[k]
                   - f_2 * ski1_363[k]
                   + f_3 * pc_z[k] * skk_467[k];

        t_585[k] = pb_y[k] * sil0_405[k]
                   - f_14 * pc_y[k] * sil1_405[k];

        t_586[k] = f_15 * sik_324[k]
                   + f_3 * pc_y[k] * skk_468[k];

        t_587[k] = f_17 * sik_288[k]
                   + f_3 * pc_z[k] * skk_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_y, pc_y, sil0_408, sil0_410, sil0_411, \
                         sik_325, sik_326, sik_327, sil1_408, sil1_410, sil1_411, \
                         skk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_y[k] * sil0_408[k]
                   + f_16 * sik_325[k]
                   - f_14 * pc_y[k] * sil1_408[k];

        t_589[k] = f_15 * sik_326[k]
                   + f_3 * pc_y[k] * skk_470[k];

        t_590[k] = pb_y[k] * sil0_410[k]
                   - f_14 * pc_y[k] * sil1_410[k];

        t_591[k] = pb_y[k] * sil0_411[k]
                   + f_17 * sik_327[k]
                   - f_14 * pc_y[k] * sil1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pb_y, pc_y, pc_z, sil0_414, sil0_415, \
                         sik_291, sik_329, sik_330, sil1_414, sil1_415, skk_471, \
                         skk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * sik_291[k]
                   + f_3 * pc_z[k] * skk_471[k];

        t_593[k] = f_15 * sik_329[k]
                   + f_3 * pc_y[k] * skk_473[k];

        t_594[k] = pb_y[k] * sil0_414[k]
                   - f_14 * pc_y[k] * sil1_414[k];

        t_595[k] = pb_y[k] * sil0_415[k]
                   + f_18 * sik_330[k]
                   - f_14 * pc_y[k] * sil1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_y, pc_y, pc_z, sil0_417, sil0_419, \
                         sik_294, sik_332, sik_333, sil1_417, sil1_419, skk_474, \
                         skk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * sik_294[k]
                   + f_3 * pc_z[k] * skk_474[k];

        t_597[k] = pb_y[k] * sil0_417[k]
                   + f_16 * sik_332[k]
                   - f_14 * pc_y[k] * sil1_417[k];

        t_598[k] = f_15 * sik_333[k]
                   + f_3 * pc_y[k] * skk_477[k];

        t_599[k] = pb_y[k] * sil0_419[k]
                   - f_14 * pc_y[k] * sil1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pb_y, pc_y, pc_z, sil0_420, sil0_422, sik_298, \
                         sik_334, sik_336, sil1_420, sil1_422, \
                         skk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pb_y[k] * sil0_420[k]
                   + f_19 * sik_334[k]
                   - f_14 * pc_y[k] * sil1_420[k];

        t_601[k] = f_17 * sik_298[k]
                   + f_3 * pc_z[k] * skk_478[k];

        t_602[k] = pb_y[k] * sil0_422[k]
                   + f_17 * sik_336[k]
                   - f_14 * pc_y[k] * sil1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pb_y, pc_y, sil0_423, sil0_425, sil0_426, \
                         sik_337, sik_338, sik_339, sil1_423, sil1_425, sil1_426, \
                         skk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pb_y[k] * sil0_423[k]
                   + f_16 * sik_337[k]
                   - f_14 * pc_y[k] * sil1_423[k];

        t_604[k] = f_15 * sik_338[k]
                   + f_3 * pc_y[k] * skk_482[k];

        t_605[k] = pb_y[k] * sil0_425[k]
                   - f_14 * pc_y[k] * sil1_425[k];

        t_606[k] = pb_y[k] * sil0_426[k]
                   + f_20 * sik_339[k]
                   - f_14 * pc_y[k] * sil1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pb_y, pc_y, pc_z, sil0_428, sil0_429, sik_303, \
                         sik_341, sik_342, sil1_428, sil1_429, \
                         skk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * sik_303[k]
                   + f_3 * pc_z[k] * skk_483[k];

        t_608[k] = pb_y[k] * sil0_428[k]
                   + f_18 * sik_341[k]
                   - f_14 * pc_y[k] * sil1_428[k];

        t_609[k] = pb_y[k] * sil0_429[k]
                   + f_17 * sik_342[k]
                   - f_14 * pc_y[k] * sil1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pb_y, pc_x, pc_y, sil0_430, sil0_432, \
                         sik_343, sik_344, sik_496, sil1_430, sil1_432, skk_488, \
                         skk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pb_y[k] * sil0_430[k]
                   + f_16 * sik_343[k]
                   - f_14 * pc_y[k] * sil1_430[k];

        t_611[k] = f_15 * sik_344[k]
                   + f_3 * pc_y[k] * skk_488[k];

        t_612[k] = pb_y[k] * sil0_432[k]
                   - f_14 * pc_y[k] * sil1_432[k];

        t_613[k] = f_17 * sik_496[k]
                   + f_3 * pc_x[k] * skk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, sik_497, sik_498, sik_499, \
                         sik_500, sik_501, skk_497, skk_498, skk_499, skk_500, \
                         skk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_17 * sik_497[k]
                   + f_3 * pc_x[k] * skk_497[k];

        t_615[k] = f_17 * sik_498[k]
                   + f_3 * pc_x[k] * skk_498[k];

        t_616[k] = f_17 * sik_499[k]
                   + f_3 * pc_x[k] * skk_499[k];

        t_617[k] = f_17 * sik_500[k]
                   + f_3 * pc_x[k] * skk_500[k];

        t_618[k] = f_17 * sik_501[k]
                   + f_3 * pc_x[k] * skk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, sik_316, sik_352, \
                         sik_502, sik_503, ski0_385, ski1_385, skk_496, skk_502, \
                         skk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_17 * sik_502[k]
                   + f_3 * pc_x[k] * skk_502[k];

        t_620[k] = f_17 * sik_503[k]
                   + f_3 * pc_x[k] * skk_503[k];

        t_621[k] = f_15 * sik_352[k]
                   + f_1 * ski0_385[k]
                   - f_2 * ski1_385[k]
                   + f_3 * pc_y[k] * skk_496[k];

        t_622[k] = f_17 * sik_316[k]
                   + f_3 * pc_z[k] * skk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, sik_354, sik_355, sik_356, ski0_387, \
                         ski0_388, ski0_389, ski1_387, ski1_388, ski1_389, skk_498, skk_499, \
                         skk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * sik_354[k]
                   + f_4 * ski0_387[k]
                   - f_5 * ski1_387[k]
                   + f_3 * pc_y[k] * skk_498[k];

        t_624[k] = f_15 * sik_355[k]
                   + f_6 * ski0_388[k]
                   - f_7 * ski1_388[k]
                   + f_3 * pc_y[k] * skk_499[k];

        t_625[k] = f_15 * sik_356[k]
                   + f_8 * ski0_389[k]
                   - f_9 * ski1_389[k]
                   + f_3 * pc_y[k] * skk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, sik_357, sik_358, sik_359, ski0_390, \
                         ski0_391, ski1_390, ski1_391, skk_501, skk_502, \
                         skk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * sik_357[k]
                   + f_10 * ski0_390[k]
                   - f_11 * ski1_390[k]
                   + f_3 * pc_y[k] * skk_501[k];

        t_627[k] = f_15 * sik_358[k]
                   + f_12 * ski0_391[k]
                   - f_13 * ski1_391[k]
                   + f_3 * pc_y[k] * skk_502[k];

        t_628[k] = f_15 * sik_359[k]
                   + f_3 * pc_y[k] * skk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pb_y, pc_x, pc_y, pc_z, sil0_449, \
                         sik_324, sik_504, sil1_449, ski0_392, ski1_392, \
                         skk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pb_y[k] * sil0_449[k]
                   - f_14 * pc_y[k] * sil1_449[k];

        t_630[k] = f_17 * sik_504[k]
                   + f_1 * ski0_392[k]
                   - f_2 * ski1_392[k]
                   + f_3 * pc_x[k] * skk_504[k];

        t_631[k] = f_3 * pc_y[k] * skk_504[k];

        t_632[k] = f_18 * sik_324[k]
                   + f_3 * pc_z[k] * skk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, sik_507, sik_509, ski0_395, \
                         ski0_397, ski1_395, ski1_397, skk_506, skk_507, \
                         skk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_17 * sik_507[k]
                   + f_4 * ski0_395[k]
                   - f_5 * ski1_395[k]
                   + f_3 * pc_x[k] * skk_507[k];

        t_634[k] = f_3 * pc_y[k] * skk_506[k];

        t_635[k] = f_17 * sik_509[k]
                   + f_4 * ski0_397[k]
                   - f_5 * ski1_397[k]
                   + f_3 * pc_x[k] * skk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_x, pc_y, pc_z, sik_327, sik_510, ski0_398, \
                         ski1_398, skk_507, skk_509, skk_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_17 * sik_510[k]
                   + f_6 * ski0_398[k]
                   - f_7 * ski1_398[k]
                   + f_3 * pc_x[k] * skk_510[k];

        t_637[k] = f_18 * sik_327[k]
                   + f_3 * pc_z[k] * skk_507[k];

        t_638[k] = f_3 * pc_y[k] * skk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_z, sik_330, sik_513, sik_514, ski0_401, \
                         ski0_402, ski1_401, ski1_402, skk_510, skk_513, \
                         skk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_17 * sik_513[k]
                   + f_6 * ski0_401[k]
                   - f_7 * ski1_401[k]
                   + f_3 * pc_x[k] * skk_513[k];

        t_640[k] = f_17 * sik_514[k]
                   + f_8 * ski0_402[k]
                   - f_9 * ski1_402[k]
                   + f_3 * pc_x[k] * skk_514[k];

        t_641[k] = f_18 * sik_330[k]
                   + f_3 * pc_z[k] * skk_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, sik_516, sik_518, ski0_404, \
                         ski0_406, ski1_404, ski1_406, skk_513, skk_516, \
                         skk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_17 * sik_516[k]
                   + f_8 * ski0_404[k]
                   - f_9 * ski1_404[k]
                   + f_3 * pc_x[k] * skk_516[k];

        t_643[k] = f_3 * pc_y[k] * skk_513[k];

        t_644[k] = f_17 * sik_518[k]
                   + f_8 * ski0_406[k]
                   - f_9 * ski1_406[k]
                   + f_3 * pc_x[k] * skk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, pc_z, sik_334, sik_519, sik_521, ski0_407, \
                         ski0_409, ski1_407, ski1_409, skk_514, skk_519, \
                         skk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_17 * sik_519[k]
                   + f_10 * ski0_407[k]
                   - f_11 * ski1_407[k]
                   + f_3 * pc_x[k] * skk_519[k];

        t_646[k] = f_18 * sik_334[k]
                   + f_3 * pc_z[k] * skk_514[k];

        t_647[k] = f_17 * sik_521[k]
                   + f_10 * ski0_409[k]
                   - f_11 * ski1_409[k]
                   + f_3 * pc_x[k] * skk_521[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, sik_522, sik_524, ski0_410, \
                         ski0_412, ski1_410, ski1_412, skk_518, skk_522, \
                         skk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_17 * sik_522[k]
                   + f_10 * ski0_410[k]
                   - f_11 * ski1_410[k]
                   + f_3 * pc_x[k] * skk_522[k];

        t_649[k] = f_3 * pc_y[k] * skk_518[k];

        t_650[k] = f_17 * sik_524[k]
                   + f_10 * ski0_412[k]
                   - f_11 * ski1_412[k]
                   + f_3 * pc_x[k] * skk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, pc_z, sik_339, sik_525, sik_527, ski0_413, \
                         ski0_415, ski1_413, ski1_415, skk_519, skk_525, \
                         skk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_17 * sik_525[k]
                   + f_12 * ski0_413[k]
                   - f_13 * ski1_413[k]
                   + f_3 * pc_x[k] * skk_525[k];

        t_652[k] = f_18 * sik_339[k]
                   + f_3 * pc_z[k] * skk_519[k];

        t_653[k] = f_17 * sik_527[k]
                   + f_12 * ski0_415[k]
                   - f_13 * ski1_415[k]
                   + f_3 * pc_x[k] * skk_527[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, sik_528, sik_529, ski0_416, \
                         ski0_417, ski1_416, ski1_417, skk_524, skk_528, \
                         skk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_17 * sik_528[k]
                   + f_12 * ski0_416[k]
                   - f_13 * ski1_416[k]
                   + f_3 * pc_x[k] * skk_528[k];

        t_655[k] = f_17 * sik_529[k]
                   + f_12 * ski0_417[k]
                   - f_13 * ski1_417[k]
                   + f_3 * pc_x[k] * skk_529[k];

        t_656[k] = f_3 * pc_y[k] * skk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, sik_531, sik_532, sik_533, sik_534, \
                         ski0_419, ski1_419, skk_531, skk_532, skk_533, \
                         skk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_17 * sik_531[k]
                   + f_12 * ski0_419[k]
                   - f_13 * ski1_419[k]
                   + f_3 * pc_x[k] * skk_531[k];

        t_658[k] = f_17 * sik_532[k]
                   + f_3 * pc_x[k] * skk_532[k];

        t_659[k] = f_17 * sik_533[k]
                   + f_3 * pc_x[k] * skk_533[k];

        t_660[k] = f_17 * sik_534[k]
                   + f_3 * pc_x[k] * skk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, sik_535, sik_536, sik_537, \
                         sik_538, sik_539, skk_535, skk_536, skk_537, skk_538, \
                         skk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_17 * sik_535[k]
                   + f_3 * pc_x[k] * skk_535[k];

        t_662[k] = f_17 * sik_536[k]
                   + f_3 * pc_x[k] * skk_536[k];

        t_663[k] = f_17 * sik_537[k]
                   + f_3 * pc_x[k] * skk_537[k];

        t_664[k] = f_17 * sik_538[k]
                   + f_3 * pc_x[k] * skk_538[k];

        t_665[k] = f_17 * sik_539[k]
                   + f_3 * pc_x[k] * skk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, sik_352, ski0_413, ski0_415, \
                         ski0_416, ski1_413, ski1_415, ski1_416, skk_532, skk_534, \
                         skk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * ski0_413[k]
                   - f_2 * ski1_413[k]
                   + f_3 * pc_y[k] * skk_532[k];

        t_667[k] = f_18 * sik_352[k]
                   + f_3 * pc_z[k] * skk_532[k];

        t_668[k] = f_4 * ski0_415[k]
                   - f_5 * ski1_415[k]
                   + f_3 * pc_y[k] * skk_534[k];

        t_669[k] = f_6 * ski0_416[k]
                   - f_7 * ski1_416[k]
                   + f_3 * pc_y[k] * skk_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, ski0_417, ski0_418, ski0_419, \
                         ski1_417, ski1_418, ski1_419, skk_536, skk_537, skk_538, \
                         skk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_8 * ski0_417[k]
                   - f_9 * ski1_417[k]
                   + f_3 * pc_y[k] * skk_536[k];

        t_671[k] = f_10 * ski0_418[k]
                   - f_11 * ski1_418[k]
                   + f_3 * pc_y[k] * skk_537[k];

        t_672[k] = f_12 * ski0_419[k]
                   - f_13 * ski1_419[k]
                   + f_3 * pc_y[k] * skk_538[k];

        t_673[k] = f_3 * pc_y[k] * skk_539[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_450 = buffer.data(sil0 + 450);
    const auto *sil0_453 = buffer.data(sil0 + 453);
    const auto *sil0_456 = buffer.data(sil0 + 456);
    const auto *sil0_460 = buffer.data(sil0 + 460);
    const auto *sil0_462 = buffer.data(sil0 + 462);
    const auto *sil0_465 = buffer.data(sil0 + 465);
    const auto *sil0_467 = buffer.data(sil0 + 467);
    const auto *sil0_468 = buffer.data(sil0 + 468);
    const auto *sil0_471 = buffer.data(sil0 + 471);
    const auto *sil0_473 = buffer.data(sil0 + 473);
    const auto *sil0_474 = buffer.data(sil0 + 474);
    const auto *sil0_475 = buffer.data(sil0 + 475);
    const auto *sil0_486 = buffer.data(sil0 + 486);

    const auto *sik_359 = buffer.data(sik + 359);
    const auto *sik_360 = buffer.data(sik + 360);
    const auto *sik_362 = buffer.data(sik + 362);
    const auto *sik_363 = buffer.data(sik + 363);
    const auto *sik_365 = buffer.data(sik + 365);
    const auto *sik_366 = buffer.data(sik + 366);
    const auto *sik_367 = buffer.data(sik + 367);
    const auto *sik_369 = buffer.data(sik + 369);
    const auto *sik_370 = buffer.data(sik + 370);
    const auto *sik_371 = buffer.data(sik + 371);
    const auto *sik_372 = buffer.data(sik + 372);
    const auto *sik_374 = buffer.data(sik + 374);
    const auto *sik_375 = buffer.data(sik + 375);
    const auto *sik_376 = buffer.data(sik + 376);
    const auto *sik_377 = buffer.data(sik + 377);
    const auto *sik_378 = buffer.data(sik + 378);
    const auto *sik_380 = buffer.data(sik + 380);
    const auto *sik_388 = buffer.data(sik + 388);
    const auto *sik_390 = buffer.data(sik + 390);
    const auto *sik_391 = buffer.data(sik + 391);
    const auto *sik_392 = buffer.data(sik + 392);
    const auto *sik_393 = buffer.data(sik + 393);
    const auto *sik_394 = buffer.data(sik + 394);
    const auto *sik_395 = buffer.data(sik + 395);
    const auto *sik_396 = buffer.data(sik + 396);
    const auto *sik_398 = buffer.data(sik + 398);
    const auto *sik_399 = buffer.data(sik + 399);
    const auto *sik_401 = buffer.data(sik + 401);
    const auto *sik_402 = buffer.data(sik + 402);
    const auto *sik_405 = buffer.data(sik + 405);
    const auto *sik_406 = buffer.data(sik + 406);
    const auto *sik_410 = buffer.data(sik + 410);
    const auto *sik_416 = buffer.data(sik + 416);
    const auto *sik_426 = buffer.data(sik + 426);
    const auto *sik_427 = buffer.data(sik + 427);
    const auto *sik_428 = buffer.data(sik + 428);
    const auto *sik_429 = buffer.data(sik + 429);
    const auto *sik_430 = buffer.data(sik + 430);
    const auto *sik_431 = buffer.data(sik + 431);
    const auto *sik_432 = buffer.data(sik + 432);
    const auto *sik_434 = buffer.data(sik + 434);
    const auto *sik_437 = buffer.data(sik + 437);
    const auto *sik_441 = buffer.data(sik + 441);
    const auto *sik_540 = buffer.data(sik + 540);
    const auto *sik_543 = buffer.data(sik + 543);
    const auto *sik_545 = buffer.data(sik + 545);
    const auto *sik_546 = buffer.data(sik + 546);
    const auto *sik_549 = buffer.data(sik + 549);
    const auto *sik_550 = buffer.data(sik + 550);
    const auto *sik_552 = buffer.data(sik + 552);
    const auto *sik_554 = buffer.data(sik + 554);
    const auto *sik_555 = buffer.data(sik + 555);
    const auto *sik_557 = buffer.data(sik + 557);
    const auto *sik_558 = buffer.data(sik + 558);
    const auto *sik_560 = buffer.data(sik + 560);
    const auto *sik_561 = buffer.data(sik + 561);
    const auto *sik_563 = buffer.data(sik + 563);
    const auto *sik_564 = buffer.data(sik + 564);
    const auto *sik_565 = buffer.data(sik + 565);
    const auto *sik_567 = buffer.data(sik + 567);
    const auto *sik_568 = buffer.data(sik + 568);
    const auto *sik_569 = buffer.data(sik + 569);
    const auto *sik_570 = buffer.data(sik + 570);
    const auto *sik_571 = buffer.data(sik + 571);
    const auto *sik_572 = buffer.data(sik + 572);
    const auto *sik_573 = buffer.data(sik + 573);
    const auto *sik_574 = buffer.data(sik + 574);
    const auto *sik_575 = buffer.data(sik + 575);
    const auto *sik_581 = buffer.data(sik + 581);
    const auto *sik_585 = buffer.data(sik + 585);
    const auto *sik_590 = buffer.data(sik + 590);
    const auto *sik_596 = buffer.data(sik + 596);
    const auto *sik_603 = buffer.data(sik + 603);
    const auto *sik_604 = buffer.data(sik + 604);
    const auto *sik_605 = buffer.data(sik + 605);
    const auto *sik_606 = buffer.data(sik + 606);
    const auto *sik_607 = buffer.data(sik + 607);
    const auto *sik_608 = buffer.data(sik + 608);
    const auto *sik_609 = buffer.data(sik + 609);
    const auto *sik_610 = buffer.data(sik + 610);
    const auto *sik_611 = buffer.data(sik + 611);
    const auto *sik_612 = buffer.data(sik + 612);
    const auto *sik_615 = buffer.data(sik + 615);
    const auto *sik_617 = buffer.data(sik + 617);
    const auto *sik_618 = buffer.data(sik + 618);
    const auto *sik_621 = buffer.data(sik + 621);
    const auto *sik_622 = buffer.data(sik + 622);
    const auto *sik_624 = buffer.data(sik + 624);
    const auto *sik_626 = buffer.data(sik + 626);
    const auto *sik_627 = buffer.data(sik + 627);

    const auto *sil1_450 = buffer.data(sil1 + 450);
    const auto *sil1_453 = buffer.data(sil1 + 453);
    const auto *sil1_456 = buffer.data(sil1 + 456);
    const auto *sil1_460 = buffer.data(sil1 + 460);
    const auto *sil1_462 = buffer.data(sil1 + 462);
    const auto *sil1_465 = buffer.data(sil1 + 465);
    const auto *sil1_467 = buffer.data(sil1 + 467);
    const auto *sil1_468 = buffer.data(sil1 + 468);
    const auto *sil1_471 = buffer.data(sil1 + 471);
    const auto *sil1_473 = buffer.data(sil1 + 473);
    const auto *sil1_474 = buffer.data(sil1 + 474);
    const auto *sil1_475 = buffer.data(sil1 + 475);
    const auto *sil1_486 = buffer.data(sil1 + 486);

    const auto *ski0_419 = buffer.data(ski0 + 419);
    const auto *ski0_420 = buffer.data(ski0 + 420);
    const auto *ski0_423 = buffer.data(ski0 + 423);
    const auto *ski0_425 = buffer.data(ski0 + 425);
    const auto *ski0_426 = buffer.data(ski0 + 426);
    const auto *ski0_429 = buffer.data(ski0 + 429);
    const auto *ski0_430 = buffer.data(ski0 + 430);
    const auto *ski0_432 = buffer.data(ski0 + 432);
    const auto *ski0_434 = buffer.data(ski0 + 434);
    const auto *ski0_435 = buffer.data(ski0 + 435);
    const auto *ski0_437 = buffer.data(ski0 + 437);
    const auto *ski0_438 = buffer.data(ski0 + 438);
    const auto *ski0_440 = buffer.data(ski0 + 440);
    const auto *ski0_441 = buffer.data(ski0 + 441);
    const auto *ski0_443 = buffer.data(ski0 + 443);
    const auto *ski0_444 = buffer.data(ski0 + 444);
    const auto *ski0_445 = buffer.data(ski0 + 445);
    const auto *ski0_446 = buffer.data(ski0 + 446);
    const auto *ski0_447 = buffer.data(ski0 + 447);
    const auto *ski0_453 = buffer.data(ski0 + 453);
    const auto *ski0_457 = buffer.data(ski0 + 457);
    const auto *ski0_462 = buffer.data(ski0 + 462);
    const auto *ski0_468 = buffer.data(ski0 + 468);
    const auto *ski0_471 = buffer.data(ski0 + 471);
    const auto *ski0_472 = buffer.data(ski0 + 472);
    const auto *ski0_473 = buffer.data(ski0 + 473);
    const auto *ski0_474 = buffer.data(ski0 + 474);
    const auto *ski0_475 = buffer.data(ski0 + 475);
    const auto *ski0_476 = buffer.data(ski0 + 476);
    const auto *ski0_479 = buffer.data(ski0 + 479);
    const auto *ski0_481 = buffer.data(ski0 + 481);
    const auto *ski0_482 = buffer.data(ski0 + 482);
    const auto *ski0_485 = buffer.data(ski0 + 485);
    const auto *ski0_486 = buffer.data(ski0 + 486);
    const auto *ski0_488 = buffer.data(ski0 + 488);
    const auto *ski0_490 = buffer.data(ski0 + 490);
    const auto *ski0_491 = buffer.data(ski0 + 491);

    const auto *ski1_419 = buffer.data(ski1 + 419);
    const auto *ski1_420 = buffer.data(ski1 + 420);
    const auto *ski1_423 = buffer.data(ski1 + 423);
    const auto *ski1_425 = buffer.data(ski1 + 425);
    const auto *ski1_426 = buffer.data(ski1 + 426);
    const auto *ski1_429 = buffer.data(ski1 + 429);
    const auto *ski1_430 = buffer.data(ski1 + 430);
    const auto *ski1_432 = buffer.data(ski1 + 432);
    const auto *ski1_434 = buffer.data(ski1 + 434);
    const auto *ski1_435 = buffer.data(ski1 + 435);
    const auto *ski1_437 = buffer.data(ski1 + 437);
    const auto *ski1_438 = buffer.data(ski1 + 438);
    const auto *ski1_440 = buffer.data(ski1 + 440);
    const auto *ski1_441 = buffer.data(ski1 + 441);
    const auto *ski1_443 = buffer.data(ski1 + 443);
    const auto *ski1_444 = buffer.data(ski1 + 444);
    const auto *ski1_445 = buffer.data(ski1 + 445);
    const auto *ski1_446 = buffer.data(ski1 + 446);
    const auto *ski1_447 = buffer.data(ski1 + 447);
    const auto *ski1_453 = buffer.data(ski1 + 453);
    const auto *ski1_457 = buffer.data(ski1 + 457);
    const auto *ski1_462 = buffer.data(ski1 + 462);
    const auto *ski1_468 = buffer.data(ski1 + 468);
    const auto *ski1_471 = buffer.data(ski1 + 471);
    const auto *ski1_472 = buffer.data(ski1 + 472);
    const auto *ski1_473 = buffer.data(ski1 + 473);
    const auto *ski1_474 = buffer.data(ski1 + 474);
    const auto *ski1_475 = buffer.data(ski1 + 475);
    const auto *ski1_476 = buffer.data(ski1 + 476);
    const auto *ski1_479 = buffer.data(ski1 + 479);
    const auto *ski1_481 = buffer.data(ski1 + 481);
    const auto *ski1_482 = buffer.data(ski1 + 482);
    const auto *ski1_485 = buffer.data(ski1 + 485);
    const auto *ski1_486 = buffer.data(ski1 + 486);
    const auto *ski1_488 = buffer.data(ski1 + 488);
    const auto *ski1_490 = buffer.data(ski1 + 490);
    const auto *ski1_491 = buffer.data(ski1 + 491);

    const auto *skk_539 = buffer.data(skk + 539);
    const auto *skk_540 = buffer.data(skk + 540);
    const auto *skk_542 = buffer.data(skk + 542);
    const auto *skk_543 = buffer.data(skk + 543);
    const auto *skk_545 = buffer.data(skk + 545);
    const auto *skk_546 = buffer.data(skk + 546);
    const auto *skk_549 = buffer.data(skk + 549);
    const auto *skk_550 = buffer.data(skk + 550);
    const auto *skk_552 = buffer.data(skk + 552);
    const auto *skk_554 = buffer.data(skk + 554);
    const auto *skk_555 = buffer.data(skk + 555);
    const auto *skk_557 = buffer.data(skk + 557);
    const auto *skk_558 = buffer.data(skk + 558);
    const auto *skk_560 = buffer.data(skk + 560);
    const auto *skk_561 = buffer.data(skk + 561);
    const auto *skk_563 = buffer.data(skk + 563);
    const auto *skk_564 = buffer.data(skk + 564);
    const auto *skk_565 = buffer.data(skk + 565);
    const auto *skk_567 = buffer.data(skk + 567);
    const auto *skk_568 = buffer.data(skk + 568);
    const auto *skk_569 = buffer.data(skk + 569);
    const auto *skk_570 = buffer.data(skk + 570);
    const auto *skk_571 = buffer.data(skk + 571);
    const auto *skk_572 = buffer.data(skk + 572);
    const auto *skk_573 = buffer.data(skk + 573);
    const auto *skk_574 = buffer.data(skk + 574);
    const auto *skk_575 = buffer.data(skk + 575);
    const auto *skk_576 = buffer.data(skk + 576);
    const auto *skk_578 = buffer.data(skk + 578);
    const auto *skk_579 = buffer.data(skk + 579);
    const auto *skk_581 = buffer.data(skk + 581);
    const auto *skk_582 = buffer.data(skk + 582);
    const auto *skk_585 = buffer.data(skk + 585);
    const auto *skk_586 = buffer.data(skk + 586);
    const auto *skk_590 = buffer.data(skk + 590);
    const auto *skk_591 = buffer.data(skk + 591);
    const auto *skk_596 = buffer.data(skk + 596);
    const auto *skk_603 = buffer.data(skk + 603);
    const auto *skk_604 = buffer.data(skk + 604);
    const auto *skk_605 = buffer.data(skk + 605);
    const auto *skk_606 = buffer.data(skk + 606);
    const auto *skk_607 = buffer.data(skk + 607);
    const auto *skk_608 = buffer.data(skk + 608);
    const auto *skk_609 = buffer.data(skk + 609);
    const auto *skk_610 = buffer.data(skk + 610);
    const auto *skk_611 = buffer.data(skk + 611);
    const auto *skk_612 = buffer.data(skk + 612);
    const auto *skk_614 = buffer.data(skk + 614);
    const auto *skk_615 = buffer.data(skk + 615);
    const auto *skk_617 = buffer.data(skk + 617);
    const auto *skk_618 = buffer.data(skk + 618);
    const auto *skk_621 = buffer.data(skk + 621);
    const auto *skk_622 = buffer.data(skk + 622);
    const auto *skk_624 = buffer.data(skk + 624);
    const auto *skk_626 = buffer.data(skk + 626);
    const auto *skk_627 = buffer.data(skk + 627);

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pc_x, pc_y, pc_z, sik_359, sik_360, \
                         sik_540, ski0_419, ski0_420, ski1_419, ski1_420, skk_539, \
                         skk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_18 * sik_359[k]
                   + f_1 * ski0_419[k]
                   - f_2 * ski1_419[k]
                   + f_3 * pc_z[k] * skk_539[k];

        t_675[k] = f_16 * sik_540[k]
                   + f_1 * ski0_420[k]
                   - f_2 * ski1_420[k]
                   + f_3 * pc_x[k] * skk_540[k];

        t_676[k] = f_19 * sik_360[k]
                   + f_3 * pc_y[k] * skk_540[k];

        t_677[k] = f_3 * pc_z[k] * skk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pc_x, pc_y, sik_362, sik_543, sik_545, ski0_423, \
                         ski0_425, ski1_423, ski1_425, skk_542, skk_543, \
                         skk_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_16 * sik_543[k]
                   + f_4 * ski0_423[k]
                   - f_5 * ski1_423[k]
                   + f_3 * pc_x[k] * skk_543[k];

        t_679[k] = f_19 * sik_362[k]
                   + f_3 * pc_y[k] * skk_542[k];

        t_680[k] = f_16 * sik_545[k]
                   + f_4 * ski0_425[k]
                   - f_5 * ski1_425[k]
                   + f_3 * pc_x[k] * skk_545[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, pc_x, pc_y, pc_z, sik_365, sik_546, ski0_426, \
                         ski1_426, skk_543, skk_545, skk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_16 * sik_546[k]
                   + f_6 * ski0_426[k]
                   - f_7 * ski1_426[k]
                   + f_3 * pc_x[k] * skk_546[k];

        t_682[k] = f_3 * pc_z[k] * skk_543[k];

        t_683[k] = f_19 * sik_365[k]
                   + f_3 * pc_y[k] * skk_545[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, pc_x, pc_z, sik_549, sik_550, ski0_429, \
                         ski0_430, ski1_429, ski1_430, skk_546, skk_549, \
                         skk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_16 * sik_549[k]
                   + f_6 * ski0_429[k]
                   - f_7 * ski1_429[k]
                   + f_3 * pc_x[k] * skk_549[k];

        t_685[k] = f_16 * sik_550[k]
                   + f_8 * ski0_430[k]
                   - f_9 * ski1_430[k]
                   + f_3 * pc_x[k] * skk_550[k];

        t_686[k] = f_3 * pc_z[k] * skk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_x, pc_y, sik_369, sik_552, sik_554, ski0_432, \
                         ski0_434, ski1_432, ski1_434, skk_549, skk_552, \
                         skk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_16 * sik_552[k]
                   + f_8 * ski0_432[k]
                   - f_9 * ski1_432[k]
                   + f_3 * pc_x[k] * skk_552[k];

        t_688[k] = f_19 * sik_369[k]
                   + f_3 * pc_y[k] * skk_549[k];

        t_689[k] = f_16 * sik_554[k]
                   + f_8 * ski0_434[k]
                   - f_9 * ski1_434[k]
                   + f_3 * pc_x[k] * skk_554[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, sik_555, sik_557, ski0_435, \
                         ski0_437, ski1_435, ski1_437, skk_550, skk_555, \
                         skk_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_16 * sik_555[k]
                   + f_10 * ski0_435[k]
                   - f_11 * ski1_435[k]
                   + f_3 * pc_x[k] * skk_555[k];

        t_691[k] = f_3 * pc_z[k] * skk_550[k];

        t_692[k] = f_16 * sik_557[k]
                   + f_10 * ski0_437[k]
                   - f_11 * ski1_437[k]
                   + f_3 * pc_x[k] * skk_557[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_x, pc_y, sik_374, sik_558, sik_560, ski0_438, \
                         ski0_440, ski1_438, ski1_440, skk_554, skk_558, \
                         skk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_16 * sik_558[k]
                   + f_10 * ski0_438[k]
                   - f_11 * ski1_438[k]
                   + f_3 * pc_x[k] * skk_558[k];

        t_694[k] = f_19 * sik_374[k]
                   + f_3 * pc_y[k] * skk_554[k];

        t_695[k] = f_16 * sik_560[k]
                   + f_10 * ski0_440[k]
                   - f_11 * ski1_440[k]
                   + f_3 * pc_x[k] * skk_560[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, sik_561, sik_563, ski0_441, \
                         ski0_443, ski1_441, ski1_443, skk_555, skk_561, \
                         skk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_16 * sik_561[k]
                   + f_12 * ski0_441[k]
                   - f_13 * ski1_441[k]
                   + f_3 * pc_x[k] * skk_561[k];

        t_697[k] = f_3 * pc_z[k] * skk_555[k];

        t_698[k] = f_16 * sik_563[k]
                   + f_12 * ski0_443[k]
                   - f_13 * ski1_443[k]
                   + f_3 * pc_x[k] * skk_563[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pc_x, pc_y, sik_380, sik_564, sik_565, ski0_444, \
                         ski0_445, ski1_444, ski1_445, skk_560, skk_564, \
                         skk_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_16 * sik_564[k]
                   + f_12 * ski0_444[k]
                   - f_13 * ski1_444[k]
                   + f_3 * pc_x[k] * skk_564[k];

        t_700[k] = f_16 * sik_565[k]
                   + f_12 * ski0_445[k]
                   - f_13 * ski1_445[k]
                   + f_3 * pc_x[k] * skk_565[k];

        t_701[k] = f_19 * sik_380[k]
                   + f_3 * pc_y[k] * skk_560[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pc_x, sik_567, sik_568, sik_569, sik_570, \
                         ski0_447, ski1_447, skk_567, skk_568, skk_569, \
                         skk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_16 * sik_567[k]
                   + f_12 * ski0_447[k]
                   - f_13 * ski1_447[k]
                   + f_3 * pc_x[k] * skk_567[k];

        t_703[k] = f_16 * sik_568[k]
                   + f_3 * pc_x[k] * skk_568[k];

        t_704[k] = f_16 * sik_569[k]
                   + f_3 * pc_x[k] * skk_569[k];

        t_705[k] = f_16 * sik_570[k]
                   + f_3 * pc_x[k] * skk_570[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, pc_x, sik_571, sik_572, sik_573, \
                         sik_574, sik_575, skk_571, skk_572, skk_573, skk_574, \
                         skk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_16 * sik_571[k]
                   + f_3 * pc_x[k] * skk_571[k];

        t_707[k] = f_16 * sik_572[k]
                   + f_3 * pc_x[k] * skk_572[k];

        t_708[k] = f_16 * sik_573[k]
                   + f_3 * pc_x[k] * skk_573[k];

        t_709[k] = f_16 * sik_574[k]
                   + f_3 * pc_x[k] * skk_574[k];

        t_710[k] = f_16 * sik_575[k]
                   + f_3 * pc_x[k] * skk_575[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pc_y, pc_z, sik_388, sik_390, ski0_441, \
                         ski0_443, ski1_441, ski1_443, skk_568, \
                         skk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_19 * sik_388[k]
                   + f_1 * ski0_441[k]
                   - f_2 * ski1_441[k]
                   + f_3 * pc_y[k] * skk_568[k];

        t_712[k] = f_3 * pc_z[k] * skk_568[k];

        t_713[k] = f_19 * sik_390[k]
                   + f_4 * ski0_443[k]
                   - f_5 * ski1_443[k]
                   + f_3 * pc_y[k] * skk_570[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, sik_391, sik_392, sik_393, ski0_444, \
                         ski0_445, ski0_446, ski1_444, ski1_445, ski1_446, skk_571, skk_572, \
                         skk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_19 * sik_391[k]
                   + f_6 * ski0_444[k]
                   - f_7 * ski1_444[k]
                   + f_3 * pc_y[k] * skk_571[k];

        t_715[k] = f_19 * sik_392[k]
                   + f_8 * ski0_445[k]
                   - f_9 * ski1_445[k]
                   + f_3 * pc_y[k] * skk_572[k];

        t_716[k] = f_19 * sik_393[k]
                   + f_10 * ski0_446[k]
                   - f_11 * ski1_446[k]
                   + f_3 * pc_y[k] * skk_573[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, pb_z, pc_y, pc_z, sil0_450, sik_394, \
                         sik_395, sil1_450, ski0_447, ski1_447, skk_574, \
                         skk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_19 * sik_394[k]
                   + f_12 * ski0_447[k]
                   - f_13 * ski1_447[k]
                   + f_3 * pc_y[k] * skk_574[k];

        t_718[k] = f_19 * sik_395[k]
                   + f_3 * pc_y[k] * skk_575[k];

        t_719[k] = f_1 * ski0_447[k]
                   - f_2 * ski1_447[k]
                   + f_3 * pc_z[k] * skk_575[k];

        t_720[k] = pb_z[k] * sil0_450[k]
                   - f_14 * pc_z[k] * sil1_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pb_z, pc_y, pc_z, sil0_453, sik_360, \
                         sik_396, sik_398, sil1_453, skk_576, skk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_18 * sik_396[k]
                   + f_3 * pc_y[k] * skk_576[k];

        t_722[k] = f_15 * sik_360[k]
                   + f_3 * pc_z[k] * skk_576[k];

        t_723[k] = pb_z[k] * sil0_453[k]
                   - f_14 * pc_z[k] * sil1_453[k];

        t_724[k] = f_18 * sik_398[k]
                   + f_3 * pc_y[k] * skk_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pb_z, pc_x, pc_z, sil0_456, sik_363, sik_581, \
                         sil1_456, ski0_453, ski1_453, skk_579, \
                         skk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_16 * sik_581[k]
                   + f_4 * ski0_453[k]
                   - f_5 * ski1_453[k]
                   + f_3 * pc_x[k] * skk_581[k];

        t_726[k] = pb_z[k] * sil0_456[k]
                   - f_14 * pc_z[k] * sil1_456[k];

        t_727[k] = f_15 * sik_363[k]
                   + f_3 * pc_z[k] * skk_579[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pb_z, pc_x, pc_y, pc_z, sil0_460, sik_401, \
                         sik_585, sil1_460, ski0_457, ski1_457, skk_581, \
                         skk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_18 * sik_401[k]
                   + f_3 * pc_y[k] * skk_581[k];

        t_729[k] = f_16 * sik_585[k]
                   + f_6 * ski0_457[k]
                   - f_7 * ski1_457[k]
                   + f_3 * pc_x[k] * skk_585[k];

        t_730[k] = pb_z[k] * sil0_460[k]
                   - f_14 * pc_z[k] * sil1_460[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_z, pc_y, pc_z, sil0_462, sik_366, sik_367, \
                         sik_405, sil1_462, skk_582, skk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_15 * sik_366[k]
                   + f_3 * pc_z[k] * skk_582[k];

        t_732[k] = pb_z[k] * sil0_462[k]
                   + f_16 * sik_367[k]
                   - f_14 * pc_z[k] * sil1_462[k];

        t_733[k] = f_18 * sik_405[k]
                   + f_3 * pc_y[k] * skk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_z, pc_x, pc_z, sil0_465, sik_370, sik_590, \
                         sil1_465, ski0_462, ski1_462, skk_586, \
                         skk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_16 * sik_590[k]
                   + f_8 * ski0_462[k]
                   - f_9 * ski1_462[k]
                   + f_3 * pc_x[k] * skk_590[k];

        t_735[k] = pb_z[k] * sil0_465[k]
                   - f_14 * pc_z[k] * sil1_465[k];

        t_736[k] = f_15 * sik_370[k]
                   + f_3 * pc_z[k] * skk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_z, pc_y, pc_z, sil0_467, sil0_468, sik_371, \
                         sik_372, sik_410, sil1_467, sil1_468, \
                         skk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_z[k] * sil0_467[k]
                   + f_16 * sik_371[k]
                   - f_14 * pc_z[k] * sil1_467[k];

        t_738[k] = pb_z[k] * sil0_468[k]
                   + f_17 * sik_372[k]
                   - f_14 * pc_z[k] * sil1_468[k];

        t_739[k] = f_18 * sik_410[k]
                   + f_3 * pc_y[k] * skk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pb_z, pc_x, pc_z, sil0_471, sik_375, sik_596, \
                         sil1_471, ski0_468, ski1_468, skk_591, \
                         skk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_16 * sik_596[k]
                   + f_10 * ski0_468[k]
                   - f_11 * ski1_468[k]
                   + f_3 * pc_x[k] * skk_596[k];

        t_741[k] = pb_z[k] * sil0_471[k]
                   - f_14 * pc_z[k] * sil1_471[k];

        t_742[k] = f_15 * sik_375[k]
                   + f_3 * pc_z[k] * skk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pb_z, pc_z, sil0_473, sil0_474, sil0_475, \
                         sik_376, sik_377, sik_378, sil1_473, sil1_474, \
                         sil1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pb_z[k] * sil0_473[k]
                   + f_16 * sik_376[k]
                   - f_14 * pc_z[k] * sil1_473[k];

        t_744[k] = pb_z[k] * sil0_474[k]
                   + f_17 * sik_377[k]
                   - f_14 * pc_z[k] * sil1_474[k];

        t_745[k] = pb_z[k] * sil0_475[k]
                   + f_18 * sik_378[k]
                   - f_14 * pc_z[k] * sil1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, sik_416, sik_603, sik_604, \
                         sik_605, ski0_475, ski1_475, skk_596, skk_603, skk_604, \
                         skk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * sik_416[k]
                   + f_3 * pc_y[k] * skk_596[k];

        t_747[k] = f_16 * sik_603[k]
                   + f_12 * ski0_475[k]
                   - f_13 * ski1_475[k]
                   + f_3 * pc_x[k] * skk_603[k];

        t_748[k] = f_16 * sik_604[k]
                   + f_3 * pc_x[k] * skk_604[k];

        t_749[k] = f_16 * sik_605[k]
                   + f_3 * pc_x[k] * skk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, sik_606, sik_607, sik_608, \
                         sik_609, sik_610, skk_606, skk_607, skk_608, skk_609, \
                         skk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_16 * sik_606[k]
                   + f_3 * pc_x[k] * skk_606[k];

        t_751[k] = f_16 * sik_607[k]
                   + f_3 * pc_x[k] * skk_607[k];

        t_752[k] = f_16 * sik_608[k]
                   + f_3 * pc_x[k] * skk_608[k];

        t_753[k] = f_16 * sik_609[k]
                   + f_3 * pc_x[k] * skk_609[k];

        t_754[k] = f_16 * sik_610[k]
                   + f_3 * pc_x[k] * skk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pb_z, pc_x, pc_z, sil0_486, sik_388, sik_611, \
                         sil1_486, skk_604, skk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_16 * sik_611[k]
                   + f_3 * pc_x[k] * skk_611[k];

        t_756[k] = pb_z[k] * sil0_486[k]
                   - f_14 * pc_z[k] * sil1_486[k];

        t_757[k] = f_15 * sik_388[k]
                   + f_3 * pc_z[k] * skk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, sik_426, sik_427, sik_428, ski0_471, \
                         ski0_472, ski0_473, ski1_471, ski1_472, ski1_473, skk_606, skk_607, \
                         skk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * sik_426[k]
                   + f_4 * ski0_471[k]
                   - f_5 * ski1_471[k]
                   + f_3 * pc_y[k] * skk_606[k];

        t_759[k] = f_18 * sik_427[k]
                   + f_6 * ski0_472[k]
                   - f_7 * ski1_472[k]
                   + f_3 * pc_y[k] * skk_607[k];

        t_760[k] = f_18 * sik_428[k]
                   + f_8 * ski0_473[k]
                   - f_9 * ski1_473[k]
                   + f_3 * pc_y[k] * skk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, sik_429, sik_430, sik_431, ski0_474, \
                         ski0_475, ski1_474, ski1_475, skk_609, skk_610, \
                         skk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * sik_429[k]
                   + f_10 * ski0_474[k]
                   - f_11 * ski1_474[k]
                   + f_3 * pc_y[k] * skk_609[k];

        t_762[k] = f_18 * sik_430[k]
                   + f_12 * ski0_475[k]
                   - f_13 * ski1_475[k]
                   + f_3 * pc_y[k] * skk_610[k];

        t_763[k] = f_18 * sik_431[k]
                   + f_3 * pc_y[k] * skk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, sik_395, sik_432, sik_612, \
                         ski0_475, ski0_476, ski1_475, ski1_476, skk_611, \
                         skk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * sik_395[k]
                   + f_1 * ski0_475[k]
                   - f_2 * ski1_475[k]
                   + f_3 * pc_z[k] * skk_611[k];

        t_765[k] = f_16 * sik_612[k]
                   + f_1 * ski0_476[k]
                   - f_2 * ski1_476[k]
                   + f_3 * pc_x[k] * skk_612[k];

        t_766[k] = f_17 * sik_432[k]
                   + f_3 * pc_y[k] * skk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, sik_396, sik_434, sik_615, \
                         ski0_479, ski1_479, skk_612, skk_614, \
                         skk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * sik_396[k]
                   + f_3 * pc_z[k] * skk_612[k];

        t_768[k] = f_16 * sik_615[k]
                   + f_4 * ski0_479[k]
                   - f_5 * ski1_479[k]
                   + f_3 * pc_x[k] * skk_615[k];

        t_769[k] = f_17 * sik_434[k]
                   + f_3 * pc_y[k] * skk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, sik_399, sik_617, sik_618, ski0_481, \
                         ski0_482, ski1_481, ski1_482, skk_615, skk_617, \
                         skk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_16 * sik_617[k]
                   + f_4 * ski0_481[k]
                   - f_5 * ski1_481[k]
                   + f_3 * pc_x[k] * skk_617[k];

        t_771[k] = f_16 * sik_618[k]
                   + f_6 * ski0_482[k]
                   - f_7 * ski1_482[k]
                   + f_3 * pc_x[k] * skk_618[k];

        t_772[k] = f_16 * sik_399[k]
                   + f_3 * pc_z[k] * skk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, sik_437, sik_621, sik_622, ski0_485, \
                         ski0_486, ski1_485, ski1_486, skk_617, skk_621, \
                         skk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * sik_437[k]
                   + f_3 * pc_y[k] * skk_617[k];

        t_774[k] = f_16 * sik_621[k]
                   + f_6 * ski0_485[k]
                   - f_7 * ski1_485[k]
                   + f_3 * pc_x[k] * skk_621[k];

        t_775[k] = f_16 * sik_622[k]
                   + f_8 * ski0_486[k]
                   - f_9 * ski1_486[k]
                   + f_3 * pc_x[k] * skk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, sik_402, sik_441, sik_624, \
                         ski0_488, ski1_488, skk_618, skk_621, \
                         skk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * sik_402[k]
                   + f_3 * pc_z[k] * skk_618[k];

        t_777[k] = f_16 * sik_624[k]
                   + f_8 * ski0_488[k]
                   - f_9 * ski1_488[k]
                   + f_3 * pc_x[k] * skk_624[k];

        t_778[k] = f_17 * sik_441[k]
                   + f_3 * pc_y[k] * skk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, sik_406, sik_626, sik_627, ski0_490, \
                         ski0_491, ski1_490, ski1_491, skk_622, skk_626, \
                         skk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_16 * sik_626[k]
                   + f_8 * ski0_490[k]
                   - f_9 * ski1_490[k]
                   + f_3 * pc_x[k] * skk_626[k];

        t_780[k] = f_16 * sik_627[k]
                   + f_10 * ski0_491[k]
                   - f_11 * ski1_491[k]
                   + f_3 * pc_x[k] * skk_627[k];

        t_781[k] = f_16 * sik_406[k]
                   + f_3 * pc_z[k] * skk_622[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_630 = buffer.data(sil0 + 630);
    const auto *sil0_633 = buffer.data(sil0 + 633);
    const auto *sil0_635 = buffer.data(sil0 + 635);
    const auto *sil0_636 = buffer.data(sil0 + 636);
    const auto *sil0_639 = buffer.data(sil0 + 639);
    const auto *sil0_640 = buffer.data(sil0 + 640);
    const auto *sil0_642 = buffer.data(sil0 + 642);
    const auto *sil0_644 = buffer.data(sil0 + 644);
    const auto *sil0_645 = buffer.data(sil0 + 645);
    const auto *sil0_647 = buffer.data(sil0 + 647);
    const auto *sil0_648 = buffer.data(sil0 + 648);
    const auto *sil0_650 = buffer.data(sil0 + 650);
    const auto *sil0_651 = buffer.data(sil0 + 651);
    const auto *sil0_653 = buffer.data(sil0 + 653);
    const auto *sil0_654 = buffer.data(sil0 + 654);
    const auto *sil0_655 = buffer.data(sil0 + 655);
    const auto *sil0_657 = buffer.data(sil0 + 657);

    const auto *sik_411 = buffer.data(sik + 411);
    const auto *sik_424 = buffer.data(sik + 424);
    const auto *sik_431 = buffer.data(sik + 431);
    const auto *sik_432 = buffer.data(sik + 432);
    const auto *sik_435 = buffer.data(sik + 435);
    const auto *sik_438 = buffer.data(sik + 438);
    const auto *sik_442 = buffer.data(sik + 442);
    const auto *sik_446 = buffer.data(sik + 446);
    const auto *sik_447 = buffer.data(sik + 447);
    const auto *sik_452 = buffer.data(sik + 452);
    const auto *sik_460 = buffer.data(sik + 460);
    const auto *sik_462 = buffer.data(sik + 462);
    const auto *sik_463 = buffer.data(sik + 463);
    const auto *sik_464 = buffer.data(sik + 464);
    const auto *sik_465 = buffer.data(sik + 465);
    const auto *sik_466 = buffer.data(sik + 466);
    const auto *sik_467 = buffer.data(sik + 467);
    const auto *sik_468 = buffer.data(sik + 468);
    const auto *sik_470 = buffer.data(sik + 470);
    const auto *sik_471 = buffer.data(sik + 471);
    const auto *sik_473 = buffer.data(sik + 473);
    const auto *sik_474 = buffer.data(sik + 474);
    const auto *sik_477 = buffer.data(sik + 477);
    const auto *sik_478 = buffer.data(sik + 478);
    const auto *sik_482 = buffer.data(sik + 482);
    const auto *sik_483 = buffer.data(sik + 483);
    const auto *sik_488 = buffer.data(sik + 488);
    const auto *sik_496 = buffer.data(sik + 496);
    const auto *sik_498 = buffer.data(sik + 498);
    const auto *sik_499 = buffer.data(sik + 499);
    const auto *sik_500 = buffer.data(sik + 500);
    const auto *sik_501 = buffer.data(sik + 501);
    const auto *sik_502 = buffer.data(sik + 502);
    const auto *sik_503 = buffer.data(sik + 503);
    const auto *sik_504 = buffer.data(sik + 504);
    const auto *sik_505 = buffer.data(sik + 505);
    const auto *sik_506 = buffer.data(sik + 506);
    const auto *sik_507 = buffer.data(sik + 507);
    const auto *sik_509 = buffer.data(sik + 509);
    const auto *sik_510 = buffer.data(sik + 510);
    const auto *sik_512 = buffer.data(sik + 512);
    const auto *sik_513 = buffer.data(sik + 513);
    const auto *sik_514 = buffer.data(sik + 514);
    const auto *sik_516 = buffer.data(sik + 516);
    const auto *sik_517 = buffer.data(sik + 517);
    const auto *sik_518 = buffer.data(sik + 518);
    const auto *sik_519 = buffer.data(sik + 519);
    const auto *sik_521 = buffer.data(sik + 521);
    const auto *sik_522 = buffer.data(sik + 522);
    const auto *sik_523 = buffer.data(sik + 523);
    const auto *sik_524 = buffer.data(sik + 524);
    const auto *sik_629 = buffer.data(sik + 629);
    const auto *sik_630 = buffer.data(sik + 630);
    const auto *sik_632 = buffer.data(sik + 632);
    const auto *sik_633 = buffer.data(sik + 633);
    const auto *sik_635 = buffer.data(sik + 635);
    const auto *sik_636 = buffer.data(sik + 636);
    const auto *sik_637 = buffer.data(sik + 637);
    const auto *sik_639 = buffer.data(sik + 639);
    const auto *sik_640 = buffer.data(sik + 640);
    const auto *sik_641 = buffer.data(sik + 641);
    const auto *sik_642 = buffer.data(sik + 642);
    const auto *sik_643 = buffer.data(sik + 643);
    const auto *sik_644 = buffer.data(sik + 644);
    const auto *sik_645 = buffer.data(sik + 645);
    const auto *sik_646 = buffer.data(sik + 646);
    const auto *sik_647 = buffer.data(sik + 647);
    const auto *sik_648 = buffer.data(sik + 648);
    const auto *sik_651 = buffer.data(sik + 651);
    const auto *sik_653 = buffer.data(sik + 653);
    const auto *sik_654 = buffer.data(sik + 654);
    const auto *sik_657 = buffer.data(sik + 657);
    const auto *sik_658 = buffer.data(sik + 658);
    const auto *sik_660 = buffer.data(sik + 660);
    const auto *sik_662 = buffer.data(sik + 662);
    const auto *sik_663 = buffer.data(sik + 663);
    const auto *sik_665 = buffer.data(sik + 665);
    const auto *sik_666 = buffer.data(sik + 666);
    const auto *sik_668 = buffer.data(sik + 668);
    const auto *sik_669 = buffer.data(sik + 669);
    const auto *sik_671 = buffer.data(sik + 671);
    const auto *sik_672 = buffer.data(sik + 672);
    const auto *sik_673 = buffer.data(sik + 673);
    const auto *sik_675 = buffer.data(sik + 675);
    const auto *sik_676 = buffer.data(sik + 676);
    const auto *sik_677 = buffer.data(sik + 677);
    const auto *sik_678 = buffer.data(sik + 678);
    const auto *sik_679 = buffer.data(sik + 679);
    const auto *sik_680 = buffer.data(sik + 680);
    const auto *sik_681 = buffer.data(sik + 681);
    const auto *sik_682 = buffer.data(sik + 682);
    const auto *sik_683 = buffer.data(sik + 683);
    const auto *sik_712 = buffer.data(sik + 712);
    const auto *sik_713 = buffer.data(sik + 713);
    const auto *sik_714 = buffer.data(sik + 714);
    const auto *sik_715 = buffer.data(sik + 715);
    const auto *sik_716 = buffer.data(sik + 716);
    const auto *sik_717 = buffer.data(sik + 717);

    const auto *sil1_630 = buffer.data(sil1 + 630);
    const auto *sil1_633 = buffer.data(sil1 + 633);
    const auto *sil1_635 = buffer.data(sil1 + 635);
    const auto *sil1_636 = buffer.data(sil1 + 636);
    const auto *sil1_639 = buffer.data(sil1 + 639);
    const auto *sil1_640 = buffer.data(sil1 + 640);
    const auto *sil1_642 = buffer.data(sil1 + 642);
    const auto *sil1_644 = buffer.data(sil1 + 644);
    const auto *sil1_645 = buffer.data(sil1 + 645);
    const auto *sil1_647 = buffer.data(sil1 + 647);
    const auto *sil1_648 = buffer.data(sil1 + 648);
    const auto *sil1_650 = buffer.data(sil1 + 650);
    const auto *sil1_651 = buffer.data(sil1 + 651);
    const auto *sil1_653 = buffer.data(sil1 + 653);
    const auto *sil1_654 = buffer.data(sil1 + 654);
    const auto *sil1_655 = buffer.data(sil1 + 655);
    const auto *sil1_657 = buffer.data(sil1 + 657);

    const auto *ski0_493 = buffer.data(ski0 + 493);
    const auto *ski0_494 = buffer.data(ski0 + 494);
    const auto *ski0_496 = buffer.data(ski0 + 496);
    const auto *ski0_497 = buffer.data(ski0 + 497);
    const auto *ski0_499 = buffer.data(ski0 + 499);
    const auto *ski0_500 = buffer.data(ski0 + 500);
    const auto *ski0_501 = buffer.data(ski0 + 501);
    const auto *ski0_502 = buffer.data(ski0 + 502);
    const auto *ski0_503 = buffer.data(ski0 + 503);
    const auto *ski0_504 = buffer.data(ski0 + 504);
    const auto *ski0_507 = buffer.data(ski0 + 507);
    const auto *ski0_509 = buffer.data(ski0 + 509);
    const auto *ski0_510 = buffer.data(ski0 + 510);
    const auto *ski0_513 = buffer.data(ski0 + 513);
    const auto *ski0_514 = buffer.data(ski0 + 514);
    const auto *ski0_516 = buffer.data(ski0 + 516);
    const auto *ski0_518 = buffer.data(ski0 + 518);
    const auto *ski0_519 = buffer.data(ski0 + 519);
    const auto *ski0_521 = buffer.data(ski0 + 521);
    const auto *ski0_522 = buffer.data(ski0 + 522);
    const auto *ski0_524 = buffer.data(ski0 + 524);
    const auto *ski0_525 = buffer.data(ski0 + 525);
    const auto *ski0_527 = buffer.data(ski0 + 527);
    const auto *ski0_528 = buffer.data(ski0 + 528);
    const auto *ski0_529 = buffer.data(ski0 + 529);
    const auto *ski0_530 = buffer.data(ski0 + 530);
    const auto *ski0_531 = buffer.data(ski0 + 531);

    const auto *ski1_493 = buffer.data(ski1 + 493);
    const auto *ski1_494 = buffer.data(ski1 + 494);
    const auto *ski1_496 = buffer.data(ski1 + 496);
    const auto *ski1_497 = buffer.data(ski1 + 497);
    const auto *ski1_499 = buffer.data(ski1 + 499);
    const auto *ski1_500 = buffer.data(ski1 + 500);
    const auto *ski1_501 = buffer.data(ski1 + 501);
    const auto *ski1_502 = buffer.data(ski1 + 502);
    const auto *ski1_503 = buffer.data(ski1 + 503);
    const auto *ski1_504 = buffer.data(ski1 + 504);
    const auto *ski1_507 = buffer.data(ski1 + 507);
    const auto *ski1_509 = buffer.data(ski1 + 509);
    const auto *ski1_510 = buffer.data(ski1 + 510);
    const auto *ski1_513 = buffer.data(ski1 + 513);
    const auto *ski1_514 = buffer.data(ski1 + 514);
    const auto *ski1_516 = buffer.data(ski1 + 516);
    const auto *ski1_518 = buffer.data(ski1 + 518);
    const auto *ski1_519 = buffer.data(ski1 + 519);
    const auto *ski1_521 = buffer.data(ski1 + 521);
    const auto *ski1_522 = buffer.data(ski1 + 522);
    const auto *ski1_524 = buffer.data(ski1 + 524);
    const auto *ski1_525 = buffer.data(ski1 + 525);
    const auto *ski1_527 = buffer.data(ski1 + 527);
    const auto *ski1_528 = buffer.data(ski1 + 528);
    const auto *ski1_529 = buffer.data(ski1 + 529);
    const auto *ski1_530 = buffer.data(ski1 + 530);
    const auto *ski1_531 = buffer.data(ski1 + 531);

    const auto *skk_626 = buffer.data(skk + 626);
    const auto *skk_627 = buffer.data(skk + 627);
    const auto *skk_629 = buffer.data(skk + 629);
    const auto *skk_630 = buffer.data(skk + 630);
    const auto *skk_632 = buffer.data(skk + 632);
    const auto *skk_633 = buffer.data(skk + 633);
    const auto *skk_635 = buffer.data(skk + 635);
    const auto *skk_636 = buffer.data(skk + 636);
    const auto *skk_637 = buffer.data(skk + 637);
    const auto *skk_639 = buffer.data(skk + 639);
    const auto *skk_640 = buffer.data(skk + 640);
    const auto *skk_641 = buffer.data(skk + 641);
    const auto *skk_642 = buffer.data(skk + 642);
    const auto *skk_643 = buffer.data(skk + 643);
    const auto *skk_644 = buffer.data(skk + 644);
    const auto *skk_645 = buffer.data(skk + 645);
    const auto *skk_646 = buffer.data(skk + 646);
    const auto *skk_647 = buffer.data(skk + 647);
    const auto *skk_648 = buffer.data(skk + 648);
    const auto *skk_650 = buffer.data(skk + 650);
    const auto *skk_651 = buffer.data(skk + 651);
    const auto *skk_653 = buffer.data(skk + 653);
    const auto *skk_654 = buffer.data(skk + 654);
    const auto *skk_657 = buffer.data(skk + 657);
    const auto *skk_658 = buffer.data(skk + 658);
    const auto *skk_660 = buffer.data(skk + 660);
    const auto *skk_662 = buffer.data(skk + 662);
    const auto *skk_663 = buffer.data(skk + 663);
    const auto *skk_665 = buffer.data(skk + 665);
    const auto *skk_666 = buffer.data(skk + 666);
    const auto *skk_668 = buffer.data(skk + 668);
    const auto *skk_669 = buffer.data(skk + 669);
    const auto *skk_671 = buffer.data(skk + 671);
    const auto *skk_672 = buffer.data(skk + 672);
    const auto *skk_673 = buffer.data(skk + 673);
    const auto *skk_675 = buffer.data(skk + 675);
    const auto *skk_676 = buffer.data(skk + 676);
    const auto *skk_677 = buffer.data(skk + 677);
    const auto *skk_678 = buffer.data(skk + 678);
    const auto *skk_679 = buffer.data(skk + 679);
    const auto *skk_680 = buffer.data(skk + 680);
    const auto *skk_681 = buffer.data(skk + 681);
    const auto *skk_682 = buffer.data(skk + 682);
    const auto *skk_683 = buffer.data(skk + 683);
    const auto *skk_684 = buffer.data(skk + 684);
    const auto *skk_686 = buffer.data(skk + 686);
    const auto *skk_687 = buffer.data(skk + 687);
    const auto *skk_689 = buffer.data(skk + 689);
    const auto *skk_690 = buffer.data(skk + 690);
    const auto *skk_693 = buffer.data(skk + 693);
    const auto *skk_694 = buffer.data(skk + 694);
    const auto *skk_698 = buffer.data(skk + 698);
    const auto *skk_699 = buffer.data(skk + 699);
    const auto *skk_704 = buffer.data(skk + 704);
    const auto *skk_712 = buffer.data(skk + 712);
    const auto *skk_713 = buffer.data(skk + 713);
    const auto *skk_714 = buffer.data(skk + 714);
    const auto *skk_715 = buffer.data(skk + 715);
    const auto *skk_716 = buffer.data(skk + 716);
    const auto *skk_717 = buffer.data(skk + 717);

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, sik_446, sik_629, sik_630, ski0_493, \
                         ski0_494, ski1_493, ski1_494, skk_626, skk_629, \
                         skk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_16 * sik_629[k]
                   + f_10 * ski0_493[k]
                   - f_11 * ski1_493[k]
                   + f_3 * pc_x[k] * skk_629[k];

        t_783[k] = f_16 * sik_630[k]
                   + f_10 * ski0_494[k]
                   - f_11 * ski1_494[k]
                   + f_3 * pc_x[k] * skk_630[k];

        t_784[k] = f_17 * sik_446[k]
                   + f_3 * pc_y[k] * skk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, sik_411, sik_632, sik_633, ski0_496, \
                         ski0_497, ski1_496, ski1_497, skk_627, skk_632, \
                         skk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_16 * sik_632[k]
                   + f_10 * ski0_496[k]
                   - f_11 * ski1_496[k]
                   + f_3 * pc_x[k] * skk_632[k];

        t_786[k] = f_16 * sik_633[k]
                   + f_12 * ski0_497[k]
                   - f_13 * ski1_497[k]
                   + f_3 * pc_x[k] * skk_633[k];

        t_787[k] = f_16 * sik_411[k]
                   + f_3 * pc_z[k] * skk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, sik_635, sik_636, sik_637, ski0_499, \
                         ski0_500, ski0_501, ski1_499, ski1_500, ski1_501, skk_635, skk_636, \
                         skk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_16 * sik_635[k]
                   + f_12 * ski0_499[k]
                   - f_13 * ski1_499[k]
                   + f_3 * pc_x[k] * skk_635[k];

        t_789[k] = f_16 * sik_636[k]
                   + f_12 * ski0_500[k]
                   - f_13 * ski1_500[k]
                   + f_3 * pc_x[k] * skk_636[k];

        t_790[k] = f_16 * sik_637[k]
                   + f_12 * ski0_501[k]
                   - f_13 * ski1_501[k]
                   + f_3 * pc_x[k] * skk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, sik_452, sik_639, sik_640, \
                         sik_641, ski0_503, ski1_503, skk_632, skk_639, skk_640, \
                         skk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * sik_452[k]
                   + f_3 * pc_y[k] * skk_632[k];

        t_792[k] = f_16 * sik_639[k]
                   + f_12 * ski0_503[k]
                   - f_13 * ski1_503[k]
                   + f_3 * pc_x[k] * skk_639[k];

        t_793[k] = f_16 * sik_640[k]
                   + f_3 * pc_x[k] * skk_640[k];

        t_794[k] = f_16 * sik_641[k]
                   + f_3 * pc_x[k] * skk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, sik_642, sik_643, sik_644, \
                         sik_645, sik_646, skk_642, skk_643, skk_644, skk_645, \
                         skk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_16 * sik_642[k]
                   + f_3 * pc_x[k] * skk_642[k];

        t_796[k] = f_16 * sik_643[k]
                   + f_3 * pc_x[k] * skk_643[k];

        t_797[k] = f_16 * sik_644[k]
                   + f_3 * pc_x[k] * skk_644[k];

        t_798[k] = f_16 * sik_645[k]
                   + f_3 * pc_x[k] * skk_645[k];

        t_799[k] = f_16 * sik_646[k]
                   + f_3 * pc_x[k] * skk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, sik_424, sik_460, sik_647, \
                         ski0_497, ski1_497, skk_640, skk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_16 * sik_647[k]
                   + f_3 * pc_x[k] * skk_647[k];

        t_801[k] = f_17 * sik_460[k]
                   + f_1 * ski0_497[k]
                   - f_2 * ski1_497[k]
                   + f_3 * pc_y[k] * skk_640[k];

        t_802[k] = f_16 * sik_424[k]
                   + f_3 * pc_z[k] * skk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, sik_462, sik_463, sik_464, ski0_499, \
                         ski0_500, ski0_501, ski1_499, ski1_500, ski1_501, skk_642, skk_643, \
                         skk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * sik_462[k]
                   + f_4 * ski0_499[k]
                   - f_5 * ski1_499[k]
                   + f_3 * pc_y[k] * skk_642[k];

        t_804[k] = f_17 * sik_463[k]
                   + f_6 * ski0_500[k]
                   - f_7 * ski1_500[k]
                   + f_3 * pc_y[k] * skk_643[k];

        t_805[k] = f_17 * sik_464[k]
                   + f_8 * ski0_501[k]
                   - f_9 * ski1_501[k]
                   + f_3 * pc_y[k] * skk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, sik_465, sik_466, sik_467, ski0_502, \
                         ski0_503, ski1_502, ski1_503, skk_645, skk_646, \
                         skk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * sik_465[k]
                   + f_10 * ski0_502[k]
                   - f_11 * ski1_502[k]
                   + f_3 * pc_y[k] * skk_645[k];

        t_807[k] = f_17 * sik_466[k]
                   + f_12 * ski0_503[k]
                   - f_13 * ski1_503[k]
                   + f_3 * pc_y[k] * skk_646[k];

        t_808[k] = f_17 * sik_467[k]
                   + f_3 * pc_y[k] * skk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, sik_431, sik_468, sik_648, \
                         ski0_503, ski0_504, ski1_503, ski1_504, skk_647, \
                         skk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * sik_431[k]
                   + f_1 * ski0_503[k]
                   - f_2 * ski1_503[k]
                   + f_3 * pc_z[k] * skk_647[k];

        t_810[k] = f_16 * sik_648[k]
                   + f_1 * ski0_504[k]
                   - f_2 * ski1_504[k]
                   + f_3 * pc_x[k] * skk_648[k];

        t_811[k] = f_16 * sik_468[k]
                   + f_3 * pc_y[k] * skk_648[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, sik_432, sik_470, sik_651, \
                         ski0_507, ski1_507, skk_648, skk_650, \
                         skk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * sik_432[k]
                   + f_3 * pc_z[k] * skk_648[k];

        t_813[k] = f_16 * sik_651[k]
                   + f_4 * ski0_507[k]
                   - f_5 * ski1_507[k]
                   + f_3 * pc_x[k] * skk_651[k];

        t_814[k] = f_16 * sik_470[k]
                   + f_3 * pc_y[k] * skk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, sik_435, sik_653, sik_654, ski0_509, \
                         ski0_510, ski1_509, ski1_510, skk_651, skk_653, \
                         skk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_16 * sik_653[k]
                   + f_4 * ski0_509[k]
                   - f_5 * ski1_509[k]
                   + f_3 * pc_x[k] * skk_653[k];

        t_816[k] = f_16 * sik_654[k]
                   + f_6 * ski0_510[k]
                   - f_7 * ski1_510[k]
                   + f_3 * pc_x[k] * skk_654[k];

        t_817[k] = f_17 * sik_435[k]
                   + f_3 * pc_z[k] * skk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, sik_473, sik_657, sik_658, ski0_513, \
                         ski0_514, ski1_513, ski1_514, skk_653, skk_657, \
                         skk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * sik_473[k]
                   + f_3 * pc_y[k] * skk_653[k];

        t_819[k] = f_16 * sik_657[k]
                   + f_6 * ski0_513[k]
                   - f_7 * ski1_513[k]
                   + f_3 * pc_x[k] * skk_657[k];

        t_820[k] = f_16 * sik_658[k]
                   + f_8 * ski0_514[k]
                   - f_9 * ski1_514[k]
                   + f_3 * pc_x[k] * skk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, sik_438, sik_477, sik_660, \
                         ski0_516, ski1_516, skk_654, skk_657, \
                         skk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * sik_438[k]
                   + f_3 * pc_z[k] * skk_654[k];

        t_822[k] = f_16 * sik_660[k]
                   + f_8 * ski0_516[k]
                   - f_9 * ski1_516[k]
                   + f_3 * pc_x[k] * skk_660[k];

        t_823[k] = f_16 * sik_477[k]
                   + f_3 * pc_y[k] * skk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, sik_442, sik_662, sik_663, ski0_518, \
                         ski0_519, ski1_518, ski1_519, skk_658, skk_662, \
                         skk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_16 * sik_662[k]
                   + f_8 * ski0_518[k]
                   - f_9 * ski1_518[k]
                   + f_3 * pc_x[k] * skk_662[k];

        t_825[k] = f_16 * sik_663[k]
                   + f_10 * ski0_519[k]
                   - f_11 * ski1_519[k]
                   + f_3 * pc_x[k] * skk_663[k];

        t_826[k] = f_17 * sik_442[k]
                   + f_3 * pc_z[k] * skk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, sik_482, sik_665, sik_666, ski0_521, \
                         ski0_522, ski1_521, ski1_522, skk_662, skk_665, \
                         skk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_16 * sik_665[k]
                   + f_10 * ski0_521[k]
                   - f_11 * ski1_521[k]
                   + f_3 * pc_x[k] * skk_665[k];

        t_828[k] = f_16 * sik_666[k]
                   + f_10 * ski0_522[k]
                   - f_11 * ski1_522[k]
                   + f_3 * pc_x[k] * skk_666[k];

        t_829[k] = f_16 * sik_482[k]
                   + f_3 * pc_y[k] * skk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, sik_447, sik_668, sik_669, ski0_524, \
                         ski0_525, ski1_524, ski1_525, skk_663, skk_668, \
                         skk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_16 * sik_668[k]
                   + f_10 * ski0_524[k]
                   - f_11 * ski1_524[k]
                   + f_3 * pc_x[k] * skk_668[k];

        t_831[k] = f_16 * sik_669[k]
                   + f_12 * ski0_525[k]
                   - f_13 * ski1_525[k]
                   + f_3 * pc_x[k] * skk_669[k];

        t_832[k] = f_17 * sik_447[k]
                   + f_3 * pc_z[k] * skk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, sik_671, sik_672, sik_673, ski0_527, \
                         ski0_528, ski0_529, ski1_527, ski1_528, ski1_529, skk_671, skk_672, \
                         skk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_16 * sik_671[k]
                   + f_12 * ski0_527[k]
                   - f_13 * ski1_527[k]
                   + f_3 * pc_x[k] * skk_671[k];

        t_834[k] = f_16 * sik_672[k]
                   + f_12 * ski0_528[k]
                   - f_13 * ski1_528[k]
                   + f_3 * pc_x[k] * skk_672[k];

        t_835[k] = f_16 * sik_673[k]
                   + f_12 * ski0_529[k]
                   - f_13 * ski1_529[k]
                   + f_3 * pc_x[k] * skk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, sik_488, sik_675, sik_676, \
                         sik_677, ski0_531, ski1_531, skk_668, skk_675, skk_676, \
                         skk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * sik_488[k]
                   + f_3 * pc_y[k] * skk_668[k];

        t_837[k] = f_16 * sik_675[k]
                   + f_12 * ski0_531[k]
                   - f_13 * ski1_531[k]
                   + f_3 * pc_x[k] * skk_675[k];

        t_838[k] = f_16 * sik_676[k]
                   + f_3 * pc_x[k] * skk_676[k];

        t_839[k] = f_16 * sik_677[k]
                   + f_3 * pc_x[k] * skk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, sik_678, sik_679, sik_680, \
                         sik_681, sik_682, skk_678, skk_679, skk_680, skk_681, \
                         skk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_16 * sik_678[k]
                   + f_3 * pc_x[k] * skk_678[k];

        t_841[k] = f_16 * sik_679[k]
                   + f_3 * pc_x[k] * skk_679[k];

        t_842[k] = f_16 * sik_680[k]
                   + f_3 * pc_x[k] * skk_680[k];

        t_843[k] = f_16 * sik_681[k]
                   + f_3 * pc_x[k] * skk_681[k];

        t_844[k] = f_16 * sik_682[k]
                   + f_3 * pc_x[k] * skk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, sik_460, sik_496, sik_683, \
                         ski0_525, ski1_525, skk_676, skk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_16 * sik_683[k]
                   + f_3 * pc_x[k] * skk_683[k];

        t_846[k] = f_16 * sik_496[k]
                   + f_1 * ski0_525[k]
                   - f_2 * ski1_525[k]
                   + f_3 * pc_y[k] * skk_676[k];

        t_847[k] = f_17 * sik_460[k]
                   + f_3 * pc_z[k] * skk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, sik_498, sik_499, sik_500, ski0_527, \
                         ski0_528, ski0_529, ski1_527, ski1_528, ski1_529, skk_678, skk_679, \
                         skk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * sik_498[k]
                   + f_4 * ski0_527[k]
                   - f_5 * ski1_527[k]
                   + f_3 * pc_y[k] * skk_678[k];

        t_849[k] = f_16 * sik_499[k]
                   + f_6 * ski0_528[k]
                   - f_7 * ski1_528[k]
                   + f_3 * pc_y[k] * skk_679[k];

        t_850[k] = f_16 * sik_500[k]
                   + f_8 * ski0_529[k]
                   - f_9 * ski1_529[k]
                   + f_3 * pc_y[k] * skk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, sik_501, sik_502, sik_503, ski0_530, \
                         ski0_531, ski1_530, ski1_531, skk_681, skk_682, \
                         skk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * sik_501[k]
                   + f_10 * ski0_530[k]
                   - f_11 * ski1_530[k]
                   + f_3 * pc_y[k] * skk_681[k];

        t_852[k] = f_16 * sik_502[k]
                   + f_12 * ski0_531[k]
                   - f_13 * ski1_531[k]
                   + f_3 * pc_y[k] * skk_682[k];

        t_853[k] = f_16 * sik_503[k]
                   + f_3 * pc_y[k] * skk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_y, pc_y, pc_z, sil0_630, sik_467, \
                         sik_468, sik_504, sil1_630, ski0_531, ski1_531, skk_683, \
                         skk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * sik_467[k]
                   + f_1 * ski0_531[k]
                   - f_2 * ski1_531[k]
                   + f_3 * pc_z[k] * skk_683[k];

        t_855[k] = pb_y[k] * sil0_630[k]
                   - f_14 * pc_y[k] * sil1_630[k];

        t_856[k] = f_15 * sik_504[k]
                   + f_3 * pc_y[k] * skk_684[k];

        t_857[k] = f_18 * sik_468[k]
                   + f_3 * pc_z[k] * skk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_y, pc_y, sil0_633, sil0_635, sil0_636, \
                         sik_505, sik_506, sik_507, sil1_633, sil1_635, sil1_636, \
                         skk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pb_y[k] * sil0_633[k]
                   + f_16 * sik_505[k]
                   - f_14 * pc_y[k] * sil1_633[k];

        t_859[k] = f_15 * sik_506[k]
                   + f_3 * pc_y[k] * skk_686[k];

        t_860[k] = pb_y[k] * sil0_635[k]
                   - f_14 * pc_y[k] * sil1_635[k];

        t_861[k] = pb_y[k] * sil0_636[k]
                   + f_17 * sik_507[k]
                   - f_14 * pc_y[k] * sil1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_y, pc_y, pc_z, sil0_639, sil0_640, \
                         sik_471, sik_509, sik_510, sil1_639, sil1_640, skk_687, \
                         skk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * sik_471[k]
                   + f_3 * pc_z[k] * skk_687[k];

        t_863[k] = f_15 * sik_509[k]
                   + f_3 * pc_y[k] * skk_689[k];

        t_864[k] = pb_y[k] * sil0_639[k]
                   - f_14 * pc_y[k] * sil1_639[k];

        t_865[k] = pb_y[k] * sil0_640[k]
                   + f_18 * sik_510[k]
                   - f_14 * pc_y[k] * sil1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pb_y, pc_y, pc_z, sil0_642, sil0_644, \
                         sik_474, sik_512, sik_513, sil1_642, sil1_644, skk_690, \
                         skk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * sik_474[k]
                   + f_3 * pc_z[k] * skk_690[k];

        t_867[k] = pb_y[k] * sil0_642[k]
                   + f_16 * sik_512[k]
                   - f_14 * pc_y[k] * sil1_642[k];

        t_868[k] = f_15 * sik_513[k]
                   + f_3 * pc_y[k] * skk_693[k];

        t_869[k] = pb_y[k] * sil0_644[k]
                   - f_14 * pc_y[k] * sil1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pb_y, pc_y, pc_z, sil0_645, sil0_647, sik_478, \
                         sik_514, sik_516, sil1_645, sil1_647, \
                         skk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pb_y[k] * sil0_645[k]
                   + f_19 * sik_514[k]
                   - f_14 * pc_y[k] * sil1_645[k];

        t_871[k] = f_18 * sik_478[k]
                   + f_3 * pc_z[k] * skk_694[k];

        t_872[k] = pb_y[k] * sil0_647[k]
                   + f_17 * sik_516[k]
                   - f_14 * pc_y[k] * sil1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pb_y, pc_y, sil0_648, sil0_650, sil0_651, \
                         sik_517, sik_518, sik_519, sil1_648, sil1_650, sil1_651, \
                         skk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pb_y[k] * sil0_648[k]
                   + f_16 * sik_517[k]
                   - f_14 * pc_y[k] * sil1_648[k];

        t_874[k] = f_15 * sik_518[k]
                   + f_3 * pc_y[k] * skk_698[k];

        t_875[k] = pb_y[k] * sil0_650[k]
                   - f_14 * pc_y[k] * sil1_650[k];

        t_876[k] = pb_y[k] * sil0_651[k]
                   + f_20 * sik_519[k]
                   - f_14 * pc_y[k] * sil1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pb_y, pc_y, pc_z, sil0_653, sil0_654, sik_483, \
                         sik_521, sik_522, sil1_653, sil1_654, \
                         skk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * sik_483[k]
                   + f_3 * pc_z[k] * skk_699[k];

        t_878[k] = pb_y[k] * sil0_653[k]
                   + f_18 * sik_521[k]
                   - f_14 * pc_y[k] * sil1_653[k];

        t_879[k] = pb_y[k] * sil0_654[k]
                   + f_17 * sik_522[k]
                   - f_14 * pc_y[k] * sil1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pb_y, pc_x, pc_y, sil0_655, sil0_657, \
                         sik_523, sik_524, sik_712, sil1_655, sil1_657, skk_704, \
                         skk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pb_y[k] * sil0_655[k]
                   + f_16 * sik_523[k]
                   - f_14 * pc_y[k] * sil1_655[k];

        t_881[k] = f_15 * sik_524[k]
                   + f_3 * pc_y[k] * skk_704[k];

        t_882[k] = pb_y[k] * sil0_657[k]
                   - f_14 * pc_y[k] * sil1_657[k];

        t_883[k] = f_16 * sik_712[k]
                   + f_3 * pc_x[k] * skk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, sik_713, sik_714, sik_715, \
                         sik_716, sik_717, skk_713, skk_714, skk_715, skk_716, \
                         skk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_16 * sik_713[k]
                   + f_3 * pc_x[k] * skk_713[k];

        t_885[k] = f_16 * sik_714[k]
                   + f_3 * pc_x[k] * skk_714[k];

        t_886[k] = f_16 * sik_715[k]
                   + f_3 * pc_x[k] * skk_715[k];

        t_887[k] = f_16 * sik_716[k]
                   + f_3 * pc_x[k] * skk_716[k];

        t_888[k] = f_16 * sik_717[k]
                   + f_3 * pc_x[k] * skk_717[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t ski0, const size_t ski1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_674 = buffer.data(sil0 + 674);
    const auto *sil0_675 = buffer.data(sil0 + 675);
    const auto *sil0_678 = buffer.data(sil0 + 678);
    const auto *sil0_681 = buffer.data(sil0 + 681);
    const auto *sil0_685 = buffer.data(sil0 + 685);
    const auto *sil0_945 = buffer.data(sil0 + 945);
    const auto *sil0_948 = buffer.data(sil0 + 948);
    const auto *sil0_950 = buffer.data(sil0 + 950);
    const auto *sil0_951 = buffer.data(sil0 + 951);
    const auto *sil0_954 = buffer.data(sil0 + 954);
    const auto *sil0_955 = buffer.data(sil0 + 955);
    const auto *sil0_957 = buffer.data(sil0 + 957);
    const auto *sil0_959 = buffer.data(sil0 + 959);
    const auto *sil0_960 = buffer.data(sil0 + 960);
    const auto *sil0_962 = buffer.data(sil0 + 962);
    const auto *sil0_963 = buffer.data(sil0 + 963);
    const auto *sil0_965 = buffer.data(sil0 + 965);
    const auto *sil0_966 = buffer.data(sil0 + 966);
    const auto *sil0_968 = buffer.data(sil0 + 968);
    const auto *sil0_969 = buffer.data(sil0 + 969);
    const auto *sil0_970 = buffer.data(sil0 + 970);
    const auto *sil0_972 = buffer.data(sil0 + 972);
    const auto *sil0_981 = buffer.data(sil0 + 981);
    const auto *sil0_983 = buffer.data(sil0 + 983);
    const auto *sil0_984 = buffer.data(sil0 + 984);
    const auto *sil0_985 = buffer.data(sil0 + 985);
    const auto *sil0_986 = buffer.data(sil0 + 986);
    const auto *sil0_987 = buffer.data(sil0 + 987);
    const auto *sil0_989 = buffer.data(sil0 + 989);
    const auto *sil0_995 = buffer.data(sil0 + 995);
    const auto *sil0_999 = buffer.data(sil0 + 999);
    const auto *sil0_1002 = buffer.data(sil0 + 1002);

    const auto *sik_496 = buffer.data(sik + 496);
    const auto *sik_504 = buffer.data(sik + 504);
    const auto *sik_507 = buffer.data(sik + 507);
    const auto *sik_510 = buffer.data(sik + 510);
    const auto *sik_514 = buffer.data(sik + 514);
    const auto *sik_519 = buffer.data(sik + 519);
    const auto *sik_532 = buffer.data(sik + 532);
    const auto *sik_534 = buffer.data(sik + 534);
    const auto *sik_535 = buffer.data(sik + 535);
    const auto *sik_536 = buffer.data(sik + 536);
    const auto *sik_537 = buffer.data(sik + 537);
    const auto *sik_538 = buffer.data(sik + 538);
    const auto *sik_539 = buffer.data(sik + 539);
    const auto *sik_540 = buffer.data(sik + 540);
    const auto *sik_542 = buffer.data(sik + 542);
    const auto *sik_543 = buffer.data(sik + 543);
    const auto *sik_545 = buffer.data(sik + 545);
    const auto *sik_546 = buffer.data(sik + 546);
    const auto *sik_549 = buffer.data(sik + 549);
    const auto *sik_554 = buffer.data(sik + 554);
    const auto *sik_560 = buffer.data(sik + 560);
    const auto *sik_575 = buffer.data(sik + 575);
    const auto *sik_576 = buffer.data(sik + 576);
    const auto *sik_578 = buffer.data(sik + 578);
    const auto *sik_581 = buffer.data(sik + 581);
    const auto *sik_585 = buffer.data(sik + 585);
    const auto *sik_718 = buffer.data(sik + 718);
    const auto *sik_719 = buffer.data(sik + 719);
    const auto *sik_720 = buffer.data(sik + 720);
    const auto *sik_723 = buffer.data(sik + 723);
    const auto *sik_725 = buffer.data(sik + 725);
    const auto *sik_726 = buffer.data(sik + 726);
    const auto *sik_729 = buffer.data(sik + 729);
    const auto *sik_730 = buffer.data(sik + 730);
    const auto *sik_732 = buffer.data(sik + 732);
    const auto *sik_734 = buffer.data(sik + 734);
    const auto *sik_735 = buffer.data(sik + 735);
    const auto *sik_737 = buffer.data(sik + 737);
    const auto *sik_738 = buffer.data(sik + 738);
    const auto *sik_740 = buffer.data(sik + 740);
    const auto *sik_741 = buffer.data(sik + 741);
    const auto *sik_743 = buffer.data(sik + 743);
    const auto *sik_744 = buffer.data(sik + 744);
    const auto *sik_745 = buffer.data(sik + 745);
    const auto *sik_747 = buffer.data(sik + 747);
    const auto *sik_748 = buffer.data(sik + 748);
    const auto *sik_749 = buffer.data(sik + 749);
    const auto *sik_750 = buffer.data(sik + 750);
    const auto *sik_751 = buffer.data(sik + 751);
    const auto *sik_752 = buffer.data(sik + 752);
    const auto *sik_753 = buffer.data(sik + 753);
    const auto *sik_754 = buffer.data(sik + 754);
    const auto *sik_755 = buffer.data(sik + 755);
    const auto *sik_756 = buffer.data(sik + 756);
    const auto *sik_759 = buffer.data(sik + 759);
    const auto *sik_761 = buffer.data(sik + 761);
    const auto *sik_762 = buffer.data(sik + 762);
    const auto *sik_765 = buffer.data(sik + 765);
    const auto *sik_766 = buffer.data(sik + 766);
    const auto *sik_768 = buffer.data(sik + 768);
    const auto *sik_770 = buffer.data(sik + 770);
    const auto *sik_771 = buffer.data(sik + 771);
    const auto *sik_773 = buffer.data(sik + 773);
    const auto *sik_774 = buffer.data(sik + 774);
    const auto *sik_776 = buffer.data(sik + 776);
    const auto *sik_777 = buffer.data(sik + 777);
    const auto *sik_779 = buffer.data(sik + 779);
    const auto *sik_780 = buffer.data(sik + 780);
    const auto *sik_781 = buffer.data(sik + 781);
    const auto *sik_783 = buffer.data(sik + 783);
    const auto *sik_784 = buffer.data(sik + 784);
    const auto *sik_785 = buffer.data(sik + 785);
    const auto *sik_786 = buffer.data(sik + 786);
    const auto *sik_787 = buffer.data(sik + 787);
    const auto *sik_788 = buffer.data(sik + 788);
    const auto *sik_789 = buffer.data(sik + 789);
    const auto *sik_790 = buffer.data(sik + 790);
    const auto *sik_791 = buffer.data(sik + 791);
    const auto *sik_797 = buffer.data(sik + 797);
    const auto *sik_801 = buffer.data(sik + 801);
    const auto *sik_804 = buffer.data(sik + 804);

    const auto *sil1_674 = buffer.data(sil1 + 674);
    const auto *sil1_675 = buffer.data(sil1 + 675);
    const auto *sil1_678 = buffer.data(sil1 + 678);
    const auto *sil1_681 = buffer.data(sil1 + 681);
    const auto *sil1_685 = buffer.data(sil1 + 685);
    const auto *sil1_945 = buffer.data(sil1 + 945);
    const auto *sil1_948 = buffer.data(sil1 + 948);
    const auto *sil1_950 = buffer.data(sil1 + 950);
    const auto *sil1_951 = buffer.data(sil1 + 951);
    const auto *sil1_954 = buffer.data(sil1 + 954);
    const auto *sil1_955 = buffer.data(sil1 + 955);
    const auto *sil1_957 = buffer.data(sil1 + 957);
    const auto *sil1_959 = buffer.data(sil1 + 959);
    const auto *sil1_960 = buffer.data(sil1 + 960);
    const auto *sil1_962 = buffer.data(sil1 + 962);
    const auto *sil1_963 = buffer.data(sil1 + 963);
    const auto *sil1_965 = buffer.data(sil1 + 965);
    const auto *sil1_966 = buffer.data(sil1 + 966);
    const auto *sil1_968 = buffer.data(sil1 + 968);
    const auto *sil1_969 = buffer.data(sil1 + 969);
    const auto *sil1_970 = buffer.data(sil1 + 970);
    const auto *sil1_972 = buffer.data(sil1 + 972);
    const auto *sil1_981 = buffer.data(sil1 + 981);
    const auto *sil1_983 = buffer.data(sil1 + 983);
    const auto *sil1_984 = buffer.data(sil1 + 984);
    const auto *sil1_985 = buffer.data(sil1 + 985);
    const auto *sil1_986 = buffer.data(sil1 + 986);
    const auto *sil1_987 = buffer.data(sil1 + 987);
    const auto *sil1_989 = buffer.data(sil1 + 989);
    const auto *sil1_995 = buffer.data(sil1 + 995);
    const auto *sil1_999 = buffer.data(sil1 + 999);
    const auto *sil1_1002 = buffer.data(sil1 + 1002);

    const auto *ski0_553 = buffer.data(ski0 + 553);
    const auto *ski0_555 = buffer.data(ski0 + 555);
    const auto *ski0_556 = buffer.data(ski0 + 556);
    const auto *ski0_557 = buffer.data(ski0 + 557);
    const auto *ski0_558 = buffer.data(ski0 + 558);
    const auto *ski0_559 = buffer.data(ski0 + 559);
    const auto *ski0_560 = buffer.data(ski0 + 560);
    const auto *ski0_563 = buffer.data(ski0 + 563);
    const auto *ski0_565 = buffer.data(ski0 + 565);
    const auto *ski0_566 = buffer.data(ski0 + 566);
    const auto *ski0_569 = buffer.data(ski0 + 569);
    const auto *ski0_570 = buffer.data(ski0 + 570);
    const auto *ski0_572 = buffer.data(ski0 + 572);
    const auto *ski0_574 = buffer.data(ski0 + 574);
    const auto *ski0_575 = buffer.data(ski0 + 575);
    const auto *ski0_577 = buffer.data(ski0 + 577);
    const auto *ski0_578 = buffer.data(ski0 + 578);
    const auto *ski0_580 = buffer.data(ski0 + 580);
    const auto *ski0_581 = buffer.data(ski0 + 581);
    const auto *ski0_583 = buffer.data(ski0 + 583);
    const auto *ski0_584 = buffer.data(ski0 + 584);
    const auto *ski0_585 = buffer.data(ski0 + 585);
    const auto *ski0_586 = buffer.data(ski0 + 586);
    const auto *ski0_587 = buffer.data(ski0 + 587);

    const auto *ski1_553 = buffer.data(ski1 + 553);
    const auto *ski1_555 = buffer.data(ski1 + 555);
    const auto *ski1_556 = buffer.data(ski1 + 556);
    const auto *ski1_557 = buffer.data(ski1 + 557);
    const auto *ski1_558 = buffer.data(ski1 + 558);
    const auto *ski1_559 = buffer.data(ski1 + 559);
    const auto *ski1_560 = buffer.data(ski1 + 560);
    const auto *ski1_563 = buffer.data(ski1 + 563);
    const auto *ski1_565 = buffer.data(ski1 + 565);
    const auto *ski1_566 = buffer.data(ski1 + 566);
    const auto *ski1_569 = buffer.data(ski1 + 569);
    const auto *ski1_570 = buffer.data(ski1 + 570);
    const auto *ski1_572 = buffer.data(ski1 + 572);
    const auto *ski1_574 = buffer.data(ski1 + 574);
    const auto *ski1_575 = buffer.data(ski1 + 575);
    const auto *ski1_577 = buffer.data(ski1 + 577);
    const auto *ski1_578 = buffer.data(ski1 + 578);
    const auto *ski1_580 = buffer.data(ski1 + 580);
    const auto *ski1_581 = buffer.data(ski1 + 581);
    const auto *ski1_583 = buffer.data(ski1 + 583);
    const auto *ski1_584 = buffer.data(ski1 + 584);
    const auto *ski1_585 = buffer.data(ski1 + 585);
    const auto *ski1_586 = buffer.data(ski1 + 586);
    const auto *ski1_587 = buffer.data(ski1 + 587);

    const auto *skk_712 = buffer.data(skk + 712);
    const auto *skk_714 = buffer.data(skk + 714);
    const auto *skk_715 = buffer.data(skk + 715);
    const auto *skk_716 = buffer.data(skk + 716);
    const auto *skk_717 = buffer.data(skk + 717);
    const auto *skk_718 = buffer.data(skk + 718);
    const auto *skk_719 = buffer.data(skk + 719);
    const auto *skk_720 = buffer.data(skk + 720);
    const auto *skk_722 = buffer.data(skk + 722);
    const auto *skk_723 = buffer.data(skk + 723);
    const auto *skk_725 = buffer.data(skk + 725);
    const auto *skk_726 = buffer.data(skk + 726);
    const auto *skk_729 = buffer.data(skk + 729);
    const auto *skk_730 = buffer.data(skk + 730);
    const auto *skk_732 = buffer.data(skk + 732);
    const auto *skk_734 = buffer.data(skk + 734);
    const auto *skk_735 = buffer.data(skk + 735);
    const auto *skk_737 = buffer.data(skk + 737);
    const auto *skk_738 = buffer.data(skk + 738);
    const auto *skk_740 = buffer.data(skk + 740);
    const auto *skk_741 = buffer.data(skk + 741);
    const auto *skk_743 = buffer.data(skk + 743);
    const auto *skk_744 = buffer.data(skk + 744);
    const auto *skk_745 = buffer.data(skk + 745);
    const auto *skk_747 = buffer.data(skk + 747);
    const auto *skk_748 = buffer.data(skk + 748);
    const auto *skk_749 = buffer.data(skk + 749);
    const auto *skk_750 = buffer.data(skk + 750);
    const auto *skk_751 = buffer.data(skk + 751);
    const auto *skk_752 = buffer.data(skk + 752);
    const auto *skk_753 = buffer.data(skk + 753);
    const auto *skk_754 = buffer.data(skk + 754);
    const auto *skk_755 = buffer.data(skk + 755);
    const auto *skk_756 = buffer.data(skk + 756);
    const auto *skk_758 = buffer.data(skk + 758);
    const auto *skk_759 = buffer.data(skk + 759);
    const auto *skk_761 = buffer.data(skk + 761);
    const auto *skk_762 = buffer.data(skk + 762);
    const auto *skk_765 = buffer.data(skk + 765);
    const auto *skk_766 = buffer.data(skk + 766);
    const auto *skk_770 = buffer.data(skk + 770);
    const auto *skk_771 = buffer.data(skk + 771);
    const auto *skk_776 = buffer.data(skk + 776);
    const auto *skk_784 = buffer.data(skk + 784);
    const auto *skk_785 = buffer.data(skk + 785);
    const auto *skk_786 = buffer.data(skk + 786);
    const auto *skk_787 = buffer.data(skk + 787);
    const auto *skk_788 = buffer.data(skk + 788);
    const auto *skk_789 = buffer.data(skk + 789);
    const auto *skk_790 = buffer.data(skk + 790);
    const auto *skk_791 = buffer.data(skk + 791);
    const auto *skk_792 = buffer.data(skk + 792);
    const auto *skk_794 = buffer.data(skk + 794);
    const auto *skk_795 = buffer.data(skk + 795);
    const auto *skk_797 = buffer.data(skk + 797);
    const auto *skk_798 = buffer.data(skk + 798);
    const auto *skk_801 = buffer.data(skk + 801);

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, sik_496, sik_532, \
                         sik_718, sik_719, ski0_553, ski1_553, skk_712, skk_718, \
                         skk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_16 * sik_718[k]
                   + f_3 * pc_x[k] * skk_718[k];

        t_890[k] = f_16 * sik_719[k]
                   + f_3 * pc_x[k] * skk_719[k];

        t_891[k] = f_15 * sik_532[k]
                   + f_1 * ski0_553[k]
                   - f_2 * ski1_553[k]
                   + f_3 * pc_y[k] * skk_712[k];

        t_892[k] = f_18 * sik_496[k]
                   + f_3 * pc_z[k] * skk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, sik_534, sik_535, sik_536, ski0_555, \
                         ski0_556, ski0_557, ski1_555, ski1_556, ski1_557, skk_714, skk_715, \
                         skk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * sik_534[k]
                   + f_4 * ski0_555[k]
                   - f_5 * ski1_555[k]
                   + f_3 * pc_y[k] * skk_714[k];

        t_894[k] = f_15 * sik_535[k]
                   + f_6 * ski0_556[k]
                   - f_7 * ski1_556[k]
                   + f_3 * pc_y[k] * skk_715[k];

        t_895[k] = f_15 * sik_536[k]
                   + f_8 * ski0_557[k]
                   - f_9 * ski1_557[k]
                   + f_3 * pc_y[k] * skk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, sik_537, sik_538, sik_539, ski0_558, \
                         ski0_559, ski1_558, ski1_559, skk_717, skk_718, \
                         skk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * sik_537[k]
                   + f_10 * ski0_558[k]
                   - f_11 * ski1_558[k]
                   + f_3 * pc_y[k] * skk_717[k];

        t_897[k] = f_15 * sik_538[k]
                   + f_12 * ski0_559[k]
                   - f_13 * ski1_559[k]
                   + f_3 * pc_y[k] * skk_718[k];

        t_898[k] = f_15 * sik_539[k]
                   + f_3 * pc_y[k] * skk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pb_y, pc_x, pc_y, pc_z, sil0_674, \
                         sik_504, sik_720, sil1_674, ski0_560, ski1_560, \
                         skk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pb_y[k] * sil0_674[k]
                   - f_14 * pc_y[k] * sil1_674[k];

        t_900[k] = f_16 * sik_720[k]
                   + f_1 * ski0_560[k]
                   - f_2 * ski1_560[k]
                   + f_3 * pc_x[k] * skk_720[k];

        t_901[k] = f_3 * pc_y[k] * skk_720[k];

        t_902[k] = f_19 * sik_504[k]
                   + f_3 * pc_z[k] * skk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, sik_723, sik_725, ski0_563, \
                         ski0_565, ski1_563, ski1_565, skk_722, skk_723, \
                         skk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_16 * sik_723[k]
                   + f_4 * ski0_563[k]
                   - f_5 * ski1_563[k]
                   + f_3 * pc_x[k] * skk_723[k];

        t_904[k] = f_3 * pc_y[k] * skk_722[k];

        t_905[k] = f_16 * sik_725[k]
                   + f_4 * ski0_565[k]
                   - f_5 * ski1_565[k]
                   + f_3 * pc_x[k] * skk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_x, pc_y, pc_z, sik_507, sik_726, ski0_566, \
                         ski1_566, skk_723, skk_725, skk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_16 * sik_726[k]
                   + f_6 * ski0_566[k]
                   - f_7 * ski1_566[k]
                   + f_3 * pc_x[k] * skk_726[k];

        t_907[k] = f_19 * sik_507[k]
                   + f_3 * pc_z[k] * skk_723[k];

        t_908[k] = f_3 * pc_y[k] * skk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_z, sik_510, sik_729, sik_730, ski0_569, \
                         ski0_570, ski1_569, ski1_570, skk_726, skk_729, \
                         skk_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_16 * sik_729[k]
                   + f_6 * ski0_569[k]
                   - f_7 * ski1_569[k]
                   + f_3 * pc_x[k] * skk_729[k];

        t_910[k] = f_16 * sik_730[k]
                   + f_8 * ski0_570[k]
                   - f_9 * ski1_570[k]
                   + f_3 * pc_x[k] * skk_730[k];

        t_911[k] = f_19 * sik_510[k]
                   + f_3 * pc_z[k] * skk_726[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, sik_732, sik_734, ski0_572, \
                         ski0_574, ski1_572, ski1_574, skk_729, skk_732, \
                         skk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_16 * sik_732[k]
                   + f_8 * ski0_572[k]
                   - f_9 * ski1_572[k]
                   + f_3 * pc_x[k] * skk_732[k];

        t_913[k] = f_3 * pc_y[k] * skk_729[k];

        t_914[k] = f_16 * sik_734[k]
                   + f_8 * ski0_574[k]
                   - f_9 * ski1_574[k]
                   + f_3 * pc_x[k] * skk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_x, pc_z, sik_514, sik_735, sik_737, ski0_575, \
                         ski0_577, ski1_575, ski1_577, skk_730, skk_735, \
                         skk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_16 * sik_735[k]
                   + f_10 * ski0_575[k]
                   - f_11 * ski1_575[k]
                   + f_3 * pc_x[k] * skk_735[k];

        t_916[k] = f_19 * sik_514[k]
                   + f_3 * pc_z[k] * skk_730[k];

        t_917[k] = f_16 * sik_737[k]
                   + f_10 * ski0_577[k]
                   - f_11 * ski1_577[k]
                   + f_3 * pc_x[k] * skk_737[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, sik_738, sik_740, ski0_578, \
                         ski0_580, ski1_578, ski1_580, skk_734, skk_738, \
                         skk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_16 * sik_738[k]
                   + f_10 * ski0_578[k]
                   - f_11 * ski1_578[k]
                   + f_3 * pc_x[k] * skk_738[k];

        t_919[k] = f_3 * pc_y[k] * skk_734[k];

        t_920[k] = f_16 * sik_740[k]
                   + f_10 * ski0_580[k]
                   - f_11 * ski1_580[k]
                   + f_3 * pc_x[k] * skk_740[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pc_x, pc_z, sik_519, sik_741, sik_743, ski0_581, \
                         ski0_583, ski1_581, ski1_583, skk_735, skk_741, \
                         skk_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_16 * sik_741[k]
                   + f_12 * ski0_581[k]
                   - f_13 * ski1_581[k]
                   + f_3 * pc_x[k] * skk_741[k];

        t_922[k] = f_19 * sik_519[k]
                   + f_3 * pc_z[k] * skk_735[k];

        t_923[k] = f_16 * sik_743[k]
                   + f_12 * ski0_583[k]
                   - f_13 * ski1_583[k]
                   + f_3 * pc_x[k] * skk_743[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_x, pc_y, sik_744, sik_745, ski0_584, \
                         ski0_585, ski1_584, ski1_585, skk_740, skk_744, \
                         skk_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_16 * sik_744[k]
                   + f_12 * ski0_584[k]
                   - f_13 * ski1_584[k]
                   + f_3 * pc_x[k] * skk_744[k];

        t_925[k] = f_16 * sik_745[k]
                   + f_12 * ski0_585[k]
                   - f_13 * ski1_585[k]
                   + f_3 * pc_x[k] * skk_745[k];

        t_926[k] = f_3 * pc_y[k] * skk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, sik_747, sik_748, sik_749, sik_750, \
                         ski0_587, ski1_587, skk_747, skk_748, skk_749, \
                         skk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_16 * sik_747[k]
                   + f_12 * ski0_587[k]
                   - f_13 * ski1_587[k]
                   + f_3 * pc_x[k] * skk_747[k];

        t_928[k] = f_16 * sik_748[k]
                   + f_3 * pc_x[k] * skk_748[k];

        t_929[k] = f_16 * sik_749[k]
                   + f_3 * pc_x[k] * skk_749[k];

        t_930[k] = f_16 * sik_750[k]
                   + f_3 * pc_x[k] * skk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, sik_751, sik_752, sik_753, \
                         sik_754, sik_755, skk_751, skk_752, skk_753, skk_754, \
                         skk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_16 * sik_751[k]
                   + f_3 * pc_x[k] * skk_751[k];

        t_932[k] = f_16 * sik_752[k]
                   + f_3 * pc_x[k] * skk_752[k];

        t_933[k] = f_16 * sik_753[k]
                   + f_3 * pc_x[k] * skk_753[k];

        t_934[k] = f_16 * sik_754[k]
                   + f_3 * pc_x[k] * skk_754[k];

        t_935[k] = f_16 * sik_755[k]
                   + f_3 * pc_x[k] * skk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pc_y, pc_z, sik_532, ski0_581, ski0_583, \
                         ski0_584, ski1_581, ski1_583, ski1_584, skk_748, skk_750, \
                         skk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * ski0_581[k]
                   - f_2 * ski1_581[k]
                   + f_3 * pc_y[k] * skk_748[k];

        t_937[k] = f_19 * sik_532[k]
                   + f_3 * pc_z[k] * skk_748[k];

        t_938[k] = f_4 * ski0_583[k]
                   - f_5 * ski1_583[k]
                   + f_3 * pc_y[k] * skk_750[k];

        t_939[k] = f_6 * ski0_584[k]
                   - f_7 * ski1_584[k]
                   + f_3 * pc_y[k] * skk_751[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pc_y, ski0_585, ski0_586, ski0_587, \
                         ski1_585, ski1_586, ski1_587, skk_752, skk_753, skk_754, \
                         skk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_8 * ski0_585[k]
                   - f_9 * ski1_585[k]
                   + f_3 * pc_y[k] * skk_752[k];

        t_941[k] = f_10 * ski0_586[k]
                   - f_11 * ski1_586[k]
                   + f_3 * pc_y[k] * skk_753[k];

        t_942[k] = f_12 * ski0_587[k]
                   - f_13 * ski1_587[k]
                   + f_3 * pc_y[k] * skk_754[k];

        t_943[k] = f_3 * pc_y[k] * skk_755[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, pb_x, pc_x, pc_y, pc_z, sil0_945, sik_539, \
                         sik_540, sik_756, sil1_945, ski0_587, ski1_587, skk_755, \
                         skk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_19 * sik_539[k]
                   + f_1 * ski0_587[k]
                   - f_2 * ski1_587[k]
                   + f_3 * pc_z[k] * skk_755[k];

        t_945[k] = pb_x[k] * sil0_945[k]
                   + f_21 * sik_756[k]
                   - f_14 * pc_x[k] * sil1_945[k];

        t_946[k] = f_20 * sik_540[k]
                   + f_3 * pc_y[k] * skk_756[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pb_x, pc_x, pc_y, pc_z, sil0_948, sik_542, \
                         sik_759, sil1_948, skk_756, skk_758 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_3 * pc_z[k] * skk_756[k];

        t_948[k] = pb_x[k] * sil0_948[k]
                   + f_20 * sik_759[k]
                   - f_14 * pc_x[k] * sil1_948[k];

        t_949[k] = f_20 * sik_542[k]
                   + f_3 * pc_y[k] * skk_758[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, pb_x, pc_x, pc_z, sil0_950, sil0_951, sik_761, \
                         sik_762, sil1_950, sil1_951, skk_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pb_x[k] * sil0_950[k]
                   + f_20 * sik_761[k]
                   - f_14 * pc_x[k] * sil1_950[k];

        t_951[k] = pb_x[k] * sil0_951[k]
                   + f_19 * sik_762[k]
                   - f_14 * pc_x[k] * sil1_951[k];

        t_952[k] = f_3 * pc_z[k] * skk_759[k];
    }

#pragma omp simd aligned(t_953, t_954, t_955, pb_x, pc_x, pc_y, sil0_954, sil0_955, sik_545, \
                         sik_765, sik_766, sil1_954, sil1_955, \
                         skk_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_953[k] = f_20 * sik_545[k]
                   + f_3 * pc_y[k] * skk_761[k];

        t_954[k] = pb_x[k] * sil0_954[k]
                   + f_19 * sik_765[k]
                   - f_14 * pc_x[k] * sil1_954[k];

        t_955[k] = pb_x[k] * sil0_955[k]
                   + f_18 * sik_766[k]
                   - f_14 * pc_x[k] * sil1_955[k];
    }

#pragma omp simd aligned(t_956, t_957, t_958, pb_x, pc_x, pc_y, pc_z, sil0_957, sik_549, \
                         sik_768, sil1_957, skk_762, skk_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_956[k] = f_3 * pc_z[k] * skk_762[k];

        t_957[k] = pb_x[k] * sil0_957[k]
                   + f_18 * sik_768[k]
                   - f_14 * pc_x[k] * sil1_957[k];

        t_958[k] = f_20 * sik_549[k]
                   + f_3 * pc_y[k] * skk_765[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, pb_x, pc_x, pc_z, sil0_959, sil0_960, sik_770, \
                         sik_771, sil1_959, sil1_960, skk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = pb_x[k] * sil0_959[k]
                   + f_18 * sik_770[k]
                   - f_14 * pc_x[k] * sil1_959[k];

        t_960[k] = pb_x[k] * sil0_960[k]
                   + f_17 * sik_771[k]
                   - f_14 * pc_x[k] * sil1_960[k];

        t_961[k] = f_3 * pc_z[k] * skk_766[k];
    }

#pragma omp simd aligned(t_962, t_963, t_964, pb_x, pc_x, pc_y, sil0_962, sil0_963, sik_554, \
                         sik_773, sik_774, sil1_962, sil1_963, \
                         skk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = pb_x[k] * sil0_962[k]
                   + f_17 * sik_773[k]
                   - f_14 * pc_x[k] * sil1_962[k];

        t_963[k] = pb_x[k] * sil0_963[k]
                   + f_17 * sik_774[k]
                   - f_14 * pc_x[k] * sil1_963[k];

        t_964[k] = f_20 * sik_554[k]
                   + f_3 * pc_y[k] * skk_770[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, pb_x, pc_x, pc_z, sil0_965, sil0_966, sik_776, \
                         sik_777, sil1_965, sil1_966, skk_771 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pb_x[k] * sil0_965[k]
                   + f_17 * sik_776[k]
                   - f_14 * pc_x[k] * sil1_965[k];

        t_966[k] = pb_x[k] * sil0_966[k]
                   + f_16 * sik_777[k]
                   - f_14 * pc_x[k] * sil1_966[k];

        t_967[k] = f_3 * pc_z[k] * skk_771[k];
    }

#pragma omp simd aligned(t_968, t_969, t_970, pb_x, pc_x, sil0_968, sil0_969, sil0_970, \
                         sik_779, sik_780, sik_781, sil1_968, sil1_969, \
                         sil1_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_968[k] = pb_x[k] * sil0_968[k]
                   + f_16 * sik_779[k]
                   - f_14 * pc_x[k] * sil1_968[k];

        t_969[k] = pb_x[k] * sil0_969[k]
                   + f_16 * sik_780[k]
                   - f_14 * pc_x[k] * sil1_969[k];

        t_970[k] = pb_x[k] * sil0_970[k]
                   + f_16 * sik_781[k]
                   - f_14 * pc_x[k] * sil1_970[k];
    }

#pragma omp simd aligned(t_971, t_972, t_973, t_974, pb_x, pc_x, pc_y, sil0_972, sik_560, \
                         sik_783, sik_784, sik_785, sil1_972, skk_776, skk_784, \
                         skk_785 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_971[k] = f_20 * sik_560[k]
                   + f_3 * pc_y[k] * skk_776[k];

        t_972[k] = pb_x[k] * sil0_972[k]
                   + f_16 * sik_783[k]
                   - f_14 * pc_x[k] * sil1_972[k];

        t_973[k] = f_15 * sik_784[k]
                   + f_3 * pc_x[k] * skk_784[k];

        t_974[k] = f_15 * sik_785[k]
                   + f_3 * pc_x[k] * skk_785[k];
    }

#pragma omp simd aligned(t_975, t_976, t_977, t_978, t_979, pc_x, sik_786, sik_787, sik_788, \
                         sik_789, sik_790, skk_786, skk_787, skk_788, skk_789, \
                         skk_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_975[k] = f_15 * sik_786[k]
                   + f_3 * pc_x[k] * skk_786[k];

        t_976[k] = f_15 * sik_787[k]
                   + f_3 * pc_x[k] * skk_787[k];

        t_977[k] = f_15 * sik_788[k]
                   + f_3 * pc_x[k] * skk_788[k];

        t_978[k] = f_15 * sik_789[k]
                   + f_3 * pc_x[k] * skk_789[k];

        t_979[k] = f_15 * sik_790[k]
                   + f_3 * pc_x[k] * skk_790[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pb_x, pc_x, pc_z, sil0_981, sil0_983, \
                         sik_791, sil1_981, sil1_983, skk_784, \
                         skk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = f_15 * sik_791[k]
                   + f_3 * pc_x[k] * skk_791[k];

        t_981[k] = pb_x[k] * sil0_981[k]
                   - f_14 * pc_x[k] * sil1_981[k];

        t_982[k] = f_3 * pc_z[k] * skk_784[k];

        t_983[k] = pb_x[k] * sil0_983[k]
                   - f_14 * pc_x[k] * sil1_983[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, pb_x, pc_x, sil0_984, sil0_985, sil0_986, \
                         sil0_987, sil1_984, sil1_985, sil1_986, \
                         sil1_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = pb_x[k] * sil0_984[k]
                   - f_14 * pc_x[k] * sil1_984[k];

        t_985[k] = pb_x[k] * sil0_985[k]
                   - f_14 * pc_x[k] * sil1_985[k];

        t_986[k] = pb_x[k] * sil0_986[k]
                   - f_14 * pc_x[k] * sil1_986[k];

        t_987[k] = pb_x[k] * sil0_987[k]
                   - f_14 * pc_x[k] * sil1_987[k];
    }

#pragma omp simd aligned(t_988, t_989, t_990, pb_x, pb_z, pc_x, pc_y, pc_z, sil0_675, \
                         sil0_989, sik_575, sil1_675, sil1_989, \
                         skk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_988[k] = f_20 * sik_575[k]
                   + f_3 * pc_y[k] * skk_791[k];

        t_989[k] = pb_x[k] * sil0_989[k]
                   - f_14 * pc_x[k] * sil1_989[k];

        t_990[k] = pb_z[k] * sil0_675[k]
                   - f_14 * pc_z[k] * sil1_675[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_z, pc_y, pc_z, sil0_678, sik_540, \
                         sik_576, sik_578, sil1_678, skk_792, skk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_19 * sik_576[k]
                   + f_3 * pc_y[k] * skk_792[k];

        t_992[k] = f_15 * sik_540[k]
                   + f_3 * pc_z[k] * skk_792[k];

        t_993[k] = pb_z[k] * sil0_678[k]
                   - f_14 * pc_z[k] * sil1_678[k];

        t_994[k] = f_19 * sik_578[k]
                   + f_3 * pc_y[k] * skk_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, pb_x, pb_z, pc_x, pc_z, sil0_681, sil0_995, \
                         sik_543, sik_797, sil1_681, sil1_995, \
                         skk_795 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = pb_x[k] * sil0_995[k]
                   + f_20 * sik_797[k]
                   - f_14 * pc_x[k] * sil1_995[k];

        t_996[k] = pb_z[k] * sil0_681[k]
                   - f_14 * pc_z[k] * sil1_681[k];

        t_997[k] = f_15 * sik_543[k]
                   + f_3 * pc_z[k] * skk_795[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, pb_x, pb_z, pc_x, pc_y, pc_z, sil0_685, \
                         sil0_999, sik_581, sik_801, sil1_685, sil1_999, \
                         skk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_19 * sik_581[k]
                   + f_3 * pc_y[k] * skk_797[k];

        t_999[k] = pb_x[k] * sil0_999[k]
                   + f_19 * sik_801[k]
                   - f_14 * pc_x[k] * sil1_999[k];

        t_1000[k] = pb_z[k] * sil0_685[k]
                    - f_14 * pc_z[k] * sil1_685[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pb_x, pc_x, pc_y, pc_z, sil0_1002, sik_546, \
                         sik_585, sik_804, sil1_1002, skk_798, \
                         skk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_15 * sik_546[k]
                    + f_3 * pc_z[k] * skk_798[k];

        t_1002[k] = pb_x[k] * sil0_1002[k]
                    + f_18 * sik_804[k]
                    - f_14 * pc_x[k] * sil1_1002[k];

        t_1003[k] = f_19 * sik_585[k]
                    + f_3 * pc_y[k] * skk_801[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sil0,
                                                          const size_t sik, const size_t sil1,
                                                          const size_t skk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_690 = buffer.data(sil0 + 690);
    const auto *sil0_696 = buffer.data(sil0 + 696);
    const auto *sil0_1004 = buffer.data(sil0 + 1004);
    const auto *sil0_1007 = buffer.data(sil0 + 1007);
    const auto *sil0_1008 = buffer.data(sil0 + 1008);
    const auto *sil0_1010 = buffer.data(sil0 + 1010);
    const auto *sil0_1013 = buffer.data(sil0 + 1013);
    const auto *sil0_1014 = buffer.data(sil0 + 1014);
    const auto *sil0_1015 = buffer.data(sil0 + 1015);
    const auto *sil0_1017 = buffer.data(sil0 + 1017);
    const auto *sil0_1026 = buffer.data(sil0 + 1026);
    const auto *sil0_1028 = buffer.data(sil0 + 1028);
    const auto *sil0_1029 = buffer.data(sil0 + 1029);
    const auto *sil0_1030 = buffer.data(sil0 + 1030);
    const auto *sil0_1031 = buffer.data(sil0 + 1031);
    const auto *sil0_1032 = buffer.data(sil0 + 1032);
    const auto *sil0_1034 = buffer.data(sil0 + 1034);
    const auto *sil0_1035 = buffer.data(sil0 + 1035);
    const auto *sil0_1038 = buffer.data(sil0 + 1038);
    const auto *sil0_1040 = buffer.data(sil0 + 1040);
    const auto *sil0_1041 = buffer.data(sil0 + 1041);
    const auto *sil0_1044 = buffer.data(sil0 + 1044);
    const auto *sil0_1045 = buffer.data(sil0 + 1045);
    const auto *sil0_1047 = buffer.data(sil0 + 1047);
    const auto *sil0_1049 = buffer.data(sil0 + 1049);
    const auto *sil0_1050 = buffer.data(sil0 + 1050);
    const auto *sil0_1052 = buffer.data(sil0 + 1052);
    const auto *sil0_1053 = buffer.data(sil0 + 1053);
    const auto *sil0_1055 = buffer.data(sil0 + 1055);
    const auto *sil0_1056 = buffer.data(sil0 + 1056);
    const auto *sil0_1058 = buffer.data(sil0 + 1058);
    const auto *sil0_1059 = buffer.data(sil0 + 1059);
    const auto *sil0_1060 = buffer.data(sil0 + 1060);
    const auto *sil0_1062 = buffer.data(sil0 + 1062);
    const auto *sil0_1071 = buffer.data(sil0 + 1071);
    const auto *sil0_1073 = buffer.data(sil0 + 1073);
    const auto *sil0_1074 = buffer.data(sil0 + 1074);
    const auto *sil0_1075 = buffer.data(sil0 + 1075);
    const auto *sil0_1076 = buffer.data(sil0 + 1076);
    const auto *sil0_1077 = buffer.data(sil0 + 1077);
    const auto *sil0_1079 = buffer.data(sil0 + 1079);
    const auto *sil0_1080 = buffer.data(sil0 + 1080);
    const auto *sil0_1083 = buffer.data(sil0 + 1083);
    const auto *sil0_1085 = buffer.data(sil0 + 1085);
    const auto *sil0_1086 = buffer.data(sil0 + 1086);
    const auto *sil0_1089 = buffer.data(sil0 + 1089);
    const auto *sil0_1090 = buffer.data(sil0 + 1090);
    const auto *sil0_1092 = buffer.data(sil0 + 1092);
    const auto *sil0_1094 = buffer.data(sil0 + 1094);
    const auto *sil0_1095 = buffer.data(sil0 + 1095);
    const auto *sil0_1097 = buffer.data(sil0 + 1097);
    const auto *sil0_1098 = buffer.data(sil0 + 1098);
    const auto *sil0_1100 = buffer.data(sil0 + 1100);
    const auto *sil0_1101 = buffer.data(sil0 + 1101);
    const auto *sil0_1103 = buffer.data(sil0 + 1103);
    const auto *sil0_1104 = buffer.data(sil0 + 1104);
    const auto *sil0_1105 = buffer.data(sil0 + 1105);
    const auto *sil0_1107 = buffer.data(sil0 + 1107);
    const auto *sil0_1116 = buffer.data(sil0 + 1116);
    const auto *sil0_1118 = buffer.data(sil0 + 1118);
    const auto *sil0_1119 = buffer.data(sil0 + 1119);
    const auto *sil0_1120 = buffer.data(sil0 + 1120);
    const auto *sil0_1121 = buffer.data(sil0 + 1121);
    const auto *sil0_1122 = buffer.data(sil0 + 1122);

    const auto *sik_550 = buffer.data(sik + 550);
    const auto *sik_555 = buffer.data(sik + 555);
    const auto *sik_568 = buffer.data(sik + 568);
    const auto *sik_576 = buffer.data(sik + 576);
    const auto *sik_579 = buffer.data(sik + 579);
    const auto *sik_582 = buffer.data(sik + 582);
    const auto *sik_586 = buffer.data(sik + 586);
    const auto *sik_590 = buffer.data(sik + 590);
    const auto *sik_591 = buffer.data(sik + 591);
    const auto *sik_596 = buffer.data(sik + 596);
    const auto *sik_604 = buffer.data(sik + 604);
    const auto *sik_611 = buffer.data(sik + 611);
    const auto *sik_612 = buffer.data(sik + 612);
    const auto *sik_614 = buffer.data(sik + 614);
    const auto *sik_615 = buffer.data(sik + 615);
    const auto *sik_617 = buffer.data(sik + 617);
    const auto *sik_618 = buffer.data(sik + 618);
    const auto *sik_621 = buffer.data(sik + 621);
    const auto *sik_622 = buffer.data(sik + 622);
    const auto *sik_626 = buffer.data(sik + 626);
    const auto *sik_627 = buffer.data(sik + 627);
    const auto *sik_632 = buffer.data(sik + 632);
    const auto *sik_640 = buffer.data(sik + 640);
    const auto *sik_647 = buffer.data(sik + 647);
    const auto *sik_648 = buffer.data(sik + 648);
    const auto *sik_650 = buffer.data(sik + 650);
    const auto *sik_653 = buffer.data(sik + 653);
    const auto *sik_657 = buffer.data(sik + 657);
    const auto *sik_662 = buffer.data(sik + 662);
    const auto *sik_668 = buffer.data(sik + 668);
    const auto *sik_806 = buffer.data(sik + 806);
    const auto *sik_809 = buffer.data(sik + 809);
    const auto *sik_810 = buffer.data(sik + 810);
    const auto *sik_812 = buffer.data(sik + 812);
    const auto *sik_815 = buffer.data(sik + 815);
    const auto *sik_816 = buffer.data(sik + 816);
    const auto *sik_817 = buffer.data(sik + 817);
    const auto *sik_819 = buffer.data(sik + 819);
    const auto *sik_820 = buffer.data(sik + 820);
    const auto *sik_821 = buffer.data(sik + 821);
    const auto *sik_822 = buffer.data(sik + 822);
    const auto *sik_823 = buffer.data(sik + 823);
    const auto *sik_824 = buffer.data(sik + 824);
    const auto *sik_825 = buffer.data(sik + 825);
    const auto *sik_826 = buffer.data(sik + 826);
    const auto *sik_827 = buffer.data(sik + 827);
    const auto *sik_828 = buffer.data(sik + 828);
    const auto *sik_831 = buffer.data(sik + 831);
    const auto *sik_833 = buffer.data(sik + 833);
    const auto *sik_834 = buffer.data(sik + 834);
    const auto *sik_837 = buffer.data(sik + 837);
    const auto *sik_838 = buffer.data(sik + 838);
    const auto *sik_840 = buffer.data(sik + 840);
    const auto *sik_842 = buffer.data(sik + 842);
    const auto *sik_843 = buffer.data(sik + 843);
    const auto *sik_845 = buffer.data(sik + 845);
    const auto *sik_846 = buffer.data(sik + 846);
    const auto *sik_848 = buffer.data(sik + 848);
    const auto *sik_849 = buffer.data(sik + 849);
    const auto *sik_851 = buffer.data(sik + 851);
    const auto *sik_852 = buffer.data(sik + 852);
    const auto *sik_853 = buffer.data(sik + 853);
    const auto *sik_855 = buffer.data(sik + 855);
    const auto *sik_856 = buffer.data(sik + 856);
    const auto *sik_857 = buffer.data(sik + 857);
    const auto *sik_858 = buffer.data(sik + 858);
    const auto *sik_859 = buffer.data(sik + 859);
    const auto *sik_860 = buffer.data(sik + 860);
    const auto *sik_861 = buffer.data(sik + 861);
    const auto *sik_862 = buffer.data(sik + 862);
    const auto *sik_863 = buffer.data(sik + 863);
    const auto *sik_864 = buffer.data(sik + 864);
    const auto *sik_867 = buffer.data(sik + 867);
    const auto *sik_869 = buffer.data(sik + 869);
    const auto *sik_870 = buffer.data(sik + 870);
    const auto *sik_873 = buffer.data(sik + 873);
    const auto *sik_874 = buffer.data(sik + 874);
    const auto *sik_876 = buffer.data(sik + 876);
    const auto *sik_878 = buffer.data(sik + 878);
    const auto *sik_879 = buffer.data(sik + 879);
    const auto *sik_881 = buffer.data(sik + 881);
    const auto *sik_882 = buffer.data(sik + 882);
    const auto *sik_884 = buffer.data(sik + 884);
    const auto *sik_885 = buffer.data(sik + 885);
    const auto *sik_887 = buffer.data(sik + 887);
    const auto *sik_888 = buffer.data(sik + 888);
    const auto *sik_889 = buffer.data(sik + 889);
    const auto *sik_891 = buffer.data(sik + 891);
    const auto *sik_892 = buffer.data(sik + 892);
    const auto *sik_893 = buffer.data(sik + 893);
    const auto *sik_894 = buffer.data(sik + 894);
    const auto *sik_895 = buffer.data(sik + 895);
    const auto *sik_896 = buffer.data(sik + 896);
    const auto *sik_897 = buffer.data(sik + 897);
    const auto *sik_898 = buffer.data(sik + 898);
    const auto *sik_899 = buffer.data(sik + 899);

    const auto *sil1_690 = buffer.data(sil1 + 690);
    const auto *sil1_696 = buffer.data(sil1 + 696);
    const auto *sil1_1004 = buffer.data(sil1 + 1004);
    const auto *sil1_1007 = buffer.data(sil1 + 1007);
    const auto *sil1_1008 = buffer.data(sil1 + 1008);
    const auto *sil1_1010 = buffer.data(sil1 + 1010);
    const auto *sil1_1013 = buffer.data(sil1 + 1013);
    const auto *sil1_1014 = buffer.data(sil1 + 1014);
    const auto *sil1_1015 = buffer.data(sil1 + 1015);
    const auto *sil1_1017 = buffer.data(sil1 + 1017);
    const auto *sil1_1026 = buffer.data(sil1 + 1026);
    const auto *sil1_1028 = buffer.data(sil1 + 1028);
    const auto *sil1_1029 = buffer.data(sil1 + 1029);
    const auto *sil1_1030 = buffer.data(sil1 + 1030);
    const auto *sil1_1031 = buffer.data(sil1 + 1031);
    const auto *sil1_1032 = buffer.data(sil1 + 1032);
    const auto *sil1_1034 = buffer.data(sil1 + 1034);
    const auto *sil1_1035 = buffer.data(sil1 + 1035);
    const auto *sil1_1038 = buffer.data(sil1 + 1038);
    const auto *sil1_1040 = buffer.data(sil1 + 1040);
    const auto *sil1_1041 = buffer.data(sil1 + 1041);
    const auto *sil1_1044 = buffer.data(sil1 + 1044);
    const auto *sil1_1045 = buffer.data(sil1 + 1045);
    const auto *sil1_1047 = buffer.data(sil1 + 1047);
    const auto *sil1_1049 = buffer.data(sil1 + 1049);
    const auto *sil1_1050 = buffer.data(sil1 + 1050);
    const auto *sil1_1052 = buffer.data(sil1 + 1052);
    const auto *sil1_1053 = buffer.data(sil1 + 1053);
    const auto *sil1_1055 = buffer.data(sil1 + 1055);
    const auto *sil1_1056 = buffer.data(sil1 + 1056);
    const auto *sil1_1058 = buffer.data(sil1 + 1058);
    const auto *sil1_1059 = buffer.data(sil1 + 1059);
    const auto *sil1_1060 = buffer.data(sil1 + 1060);
    const auto *sil1_1062 = buffer.data(sil1 + 1062);
    const auto *sil1_1071 = buffer.data(sil1 + 1071);
    const auto *sil1_1073 = buffer.data(sil1 + 1073);
    const auto *sil1_1074 = buffer.data(sil1 + 1074);
    const auto *sil1_1075 = buffer.data(sil1 + 1075);
    const auto *sil1_1076 = buffer.data(sil1 + 1076);
    const auto *sil1_1077 = buffer.data(sil1 + 1077);
    const auto *sil1_1079 = buffer.data(sil1 + 1079);
    const auto *sil1_1080 = buffer.data(sil1 + 1080);
    const auto *sil1_1083 = buffer.data(sil1 + 1083);
    const auto *sil1_1085 = buffer.data(sil1 + 1085);
    const auto *sil1_1086 = buffer.data(sil1 + 1086);
    const auto *sil1_1089 = buffer.data(sil1 + 1089);
    const auto *sil1_1090 = buffer.data(sil1 + 1090);
    const auto *sil1_1092 = buffer.data(sil1 + 1092);
    const auto *sil1_1094 = buffer.data(sil1 + 1094);
    const auto *sil1_1095 = buffer.data(sil1 + 1095);
    const auto *sil1_1097 = buffer.data(sil1 + 1097);
    const auto *sil1_1098 = buffer.data(sil1 + 1098);
    const auto *sil1_1100 = buffer.data(sil1 + 1100);
    const auto *sil1_1101 = buffer.data(sil1 + 1101);
    const auto *sil1_1103 = buffer.data(sil1 + 1103);
    const auto *sil1_1104 = buffer.data(sil1 + 1104);
    const auto *sil1_1105 = buffer.data(sil1 + 1105);
    const auto *sil1_1107 = buffer.data(sil1 + 1107);
    const auto *sil1_1116 = buffer.data(sil1 + 1116);
    const auto *sil1_1118 = buffer.data(sil1 + 1118);
    const auto *sil1_1119 = buffer.data(sil1 + 1119);
    const auto *sil1_1120 = buffer.data(sil1 + 1120);
    const auto *sil1_1121 = buffer.data(sil1 + 1121);
    const auto *sil1_1122 = buffer.data(sil1 + 1122);

    const auto *skk_802 = buffer.data(skk + 802);
    const auto *skk_806 = buffer.data(skk + 806);
    const auto *skk_807 = buffer.data(skk + 807);
    const auto *skk_812 = buffer.data(skk + 812);
    const auto *skk_820 = buffer.data(skk + 820);
    const auto *skk_821 = buffer.data(skk + 821);
    const auto *skk_822 = buffer.data(skk + 822);
    const auto *skk_823 = buffer.data(skk + 823);
    const auto *skk_824 = buffer.data(skk + 824);
    const auto *skk_825 = buffer.data(skk + 825);
    const auto *skk_826 = buffer.data(skk + 826);
    const auto *skk_827 = buffer.data(skk + 827);
    const auto *skk_828 = buffer.data(skk + 828);
    const auto *skk_830 = buffer.data(skk + 830);
    const auto *skk_831 = buffer.data(skk + 831);
    const auto *skk_833 = buffer.data(skk + 833);
    const auto *skk_834 = buffer.data(skk + 834);
    const auto *skk_837 = buffer.data(skk + 837);
    const auto *skk_838 = buffer.data(skk + 838);
    const auto *skk_842 = buffer.data(skk + 842);
    const auto *skk_843 = buffer.data(skk + 843);
    const auto *skk_848 = buffer.data(skk + 848);
    const auto *skk_856 = buffer.data(skk + 856);
    const auto *skk_857 = buffer.data(skk + 857);
    const auto *skk_858 = buffer.data(skk + 858);
    const auto *skk_859 = buffer.data(skk + 859);
    const auto *skk_860 = buffer.data(skk + 860);
    const auto *skk_861 = buffer.data(skk + 861);
    const auto *skk_862 = buffer.data(skk + 862);
    const auto *skk_863 = buffer.data(skk + 863);
    const auto *skk_864 = buffer.data(skk + 864);
    const auto *skk_866 = buffer.data(skk + 866);
    const auto *skk_867 = buffer.data(skk + 867);
    const auto *skk_869 = buffer.data(skk + 869);
    const auto *skk_870 = buffer.data(skk + 870);
    const auto *skk_873 = buffer.data(skk + 873);
    const auto *skk_874 = buffer.data(skk + 874);
    const auto *skk_878 = buffer.data(skk + 878);
    const auto *skk_879 = buffer.data(skk + 879);
    const auto *skk_884 = buffer.data(skk + 884);
    const auto *skk_892 = buffer.data(skk + 892);
    const auto *skk_893 = buffer.data(skk + 893);
    const auto *skk_894 = buffer.data(skk + 894);
    const auto *skk_895 = buffer.data(skk + 895);
    const auto *skk_896 = buffer.data(skk + 896);
    const auto *skk_897 = buffer.data(skk + 897);
    const auto *skk_898 = buffer.data(skk + 898);
    const auto *skk_899 = buffer.data(skk + 899);

#pragma omp simd aligned(t_1004, t_1005, t_1006, pb_x, pb_z, pc_x, pc_z, sil0_690, sil0_1004, \
                         sik_550, sik_806, sil1_690, sil1_1004, \
                         skk_802 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = pb_x[k] * sil0_1004[k]
                    + f_18 * sik_806[k]
                    - f_14 * pc_x[k] * sil1_1004[k];

        t_1005[k] = pb_z[k] * sil0_690[k]
                    - f_14 * pc_z[k] * sil1_690[k];

        t_1006[k] = f_15 * sik_550[k]
                    + f_3 * pc_z[k] * skk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pb_x, pc_x, pc_y, sil0_1007, sil0_1008, \
                         sik_590, sik_809, sik_810, sil1_1007, sil1_1008, \
                         skk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = pb_x[k] * sil0_1007[k]
                    + f_17 * sik_809[k]
                    - f_14 * pc_x[k] * sil1_1007[k];

        t_1008[k] = pb_x[k] * sil0_1008[k]
                    + f_17 * sik_810[k]
                    - f_14 * pc_x[k] * sil1_1008[k];

        t_1009[k] = f_19 * sik_590[k]
                    + f_3 * pc_y[k] * skk_806[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pb_x, pb_z, pc_x, pc_z, sil0_696, sil0_1010, \
                         sik_555, sik_812, sil1_696, sil1_1010, \
                         skk_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = pb_x[k] * sil0_1010[k]
                    + f_17 * sik_812[k]
                    - f_14 * pc_x[k] * sil1_1010[k];

        t_1011[k] = pb_z[k] * sil0_696[k]
                    - f_14 * pc_z[k] * sil1_696[k];

        t_1012[k] = f_15 * sik_555[k]
                    + f_3 * pc_z[k] * skk_807[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pb_x, pc_x, sil0_1013, sil0_1014, sil0_1015, \
                         sik_815, sik_816, sik_817, sil1_1013, sil1_1014, \
                         sil1_1015 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = pb_x[k] * sil0_1013[k]
                    + f_16 * sik_815[k]
                    - f_14 * pc_x[k] * sil1_1013[k];

        t_1014[k] = pb_x[k] * sil0_1014[k]
                    + f_16 * sik_816[k]
                    - f_14 * pc_x[k] * sil1_1014[k];

        t_1015[k] = pb_x[k] * sil0_1015[k]
                    + f_16 * sik_817[k]
                    - f_14 * pc_x[k] * sil1_1015[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pb_x, pc_x, pc_y, sil0_1017, sik_596, \
                         sik_819, sik_820, sik_821, sil1_1017, skk_812, skk_820, \
                         skk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * sik_596[k]
                    + f_3 * pc_y[k] * skk_812[k];

        t_1017[k] = pb_x[k] * sil0_1017[k]
                    + f_16 * sik_819[k]
                    - f_14 * pc_x[k] * sil1_1017[k];

        t_1018[k] = f_15 * sik_820[k]
                    + f_3 * pc_x[k] * skk_820[k];

        t_1019[k] = f_15 * sik_821[k]
                    + f_3 * pc_x[k] * skk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pc_x, sik_822, sik_823, \
                         sik_824, sik_825, sik_826, skk_822, skk_823, skk_824, skk_825, \
                         skk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_15 * sik_822[k]
                    + f_3 * pc_x[k] * skk_822[k];

        t_1021[k] = f_15 * sik_823[k]
                    + f_3 * pc_x[k] * skk_823[k];

        t_1022[k] = f_15 * sik_824[k]
                    + f_3 * pc_x[k] * skk_824[k];

        t_1023[k] = f_15 * sik_825[k]
                    + f_3 * pc_x[k] * skk_825[k];

        t_1024[k] = f_15 * sik_826[k]
                    + f_3 * pc_x[k] * skk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, pb_x, pc_x, pc_z, sil0_1026, \
                         sil0_1028, sik_568, sik_827, sil1_1026, sil1_1028, skk_820, \
                         skk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_15 * sik_827[k]
                    + f_3 * pc_x[k] * skk_827[k];

        t_1026[k] = pb_x[k] * sil0_1026[k]
                    - f_14 * pc_x[k] * sil1_1026[k];

        t_1027[k] = f_15 * sik_568[k]
                    + f_3 * pc_z[k] * skk_820[k];

        t_1028[k] = pb_x[k] * sil0_1028[k]
                    - f_14 * pc_x[k] * sil1_1028[k];
    }

#pragma omp simd aligned(t_1029, t_1030, t_1031, t_1032, pb_x, pc_x, sil0_1029, sil0_1030, \
                         sil0_1031, sil0_1032, sil1_1029, sil1_1030, sil1_1031, \
                         sil1_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1029[k] = pb_x[k] * sil0_1029[k]
                    - f_14 * pc_x[k] * sil1_1029[k];

        t_1030[k] = pb_x[k] * sil0_1030[k]
                    - f_14 * pc_x[k] * sil1_1030[k];

        t_1031[k] = pb_x[k] * sil0_1031[k]
                    - f_14 * pc_x[k] * sil1_1031[k];

        t_1032[k] = pb_x[k] * sil0_1032[k]
                    - f_14 * pc_x[k] * sil1_1032[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, t_1036, pb_x, pc_x, pc_y, sil0_1034, \
                         sil0_1035, sik_611, sik_612, sik_828, sil1_1034, sil1_1035, skk_827, \
                         skk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_19 * sik_611[k]
                    + f_3 * pc_y[k] * skk_827[k];

        t_1034[k] = pb_x[k] * sil0_1034[k]
                    - f_14 * pc_x[k] * sil1_1034[k];

        t_1035[k] = pb_x[k] * sil0_1035[k]
                    + f_21 * sik_828[k]
                    - f_14 * pc_x[k] * sil1_1035[k];

        t_1036[k] = f_18 * sik_612[k]
                    + f_3 * pc_y[k] * skk_828[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, pb_x, pc_x, pc_y, pc_z, sil0_1038, sik_576, \
                         sik_614, sik_831, sil1_1038, skk_828, \
                         skk_830 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_16 * sik_576[k]
                    + f_3 * pc_z[k] * skk_828[k];

        t_1038[k] = pb_x[k] * sil0_1038[k]
                    + f_20 * sik_831[k]
                    - f_14 * pc_x[k] * sil1_1038[k];

        t_1039[k] = f_18 * sik_614[k]
                    + f_3 * pc_y[k] * skk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pb_x, pc_x, pc_z, sil0_1040, sil0_1041, \
                         sik_579, sik_833, sik_834, sil1_1040, sil1_1041, \
                         skk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = pb_x[k] * sil0_1040[k]
                    + f_20 * sik_833[k]
                    - f_14 * pc_x[k] * sil1_1040[k];

        t_1041[k] = pb_x[k] * sil0_1041[k]
                    + f_19 * sik_834[k]
                    - f_14 * pc_x[k] * sil1_1041[k];

        t_1042[k] = f_16 * sik_579[k]
                    + f_3 * pc_z[k] * skk_831[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pb_x, pc_x, pc_y, sil0_1044, sil0_1045, \
                         sik_617, sik_837, sik_838, sil1_1044, sil1_1045, \
                         skk_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_18 * sik_617[k]
                    + f_3 * pc_y[k] * skk_833[k];

        t_1044[k] = pb_x[k] * sil0_1044[k]
                    + f_19 * sik_837[k]
                    - f_14 * pc_x[k] * sil1_1044[k];

        t_1045[k] = pb_x[k] * sil0_1045[k]
                    + f_18 * sik_838[k]
                    - f_14 * pc_x[k] * sil1_1045[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pb_x, pc_x, pc_y, pc_z, sil0_1047, sik_582, \
                         sik_621, sik_840, sil1_1047, skk_834, \
                         skk_837 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_16 * sik_582[k]
                    + f_3 * pc_z[k] * skk_834[k];

        t_1047[k] = pb_x[k] * sil0_1047[k]
                    + f_18 * sik_840[k]
                    - f_14 * pc_x[k] * sil1_1047[k];

        t_1048[k] = f_18 * sik_621[k]
                    + f_3 * pc_y[k] * skk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pb_x, pc_x, pc_z, sil0_1049, sil0_1050, \
                         sik_586, sik_842, sik_843, sil1_1049, sil1_1050, \
                         skk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = pb_x[k] * sil0_1049[k]
                    + f_18 * sik_842[k]
                    - f_14 * pc_x[k] * sil1_1049[k];

        t_1050[k] = pb_x[k] * sil0_1050[k]
                    + f_17 * sik_843[k]
                    - f_14 * pc_x[k] * sil1_1050[k];

        t_1051[k] = f_16 * sik_586[k]
                    + f_3 * pc_z[k] * skk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pb_x, pc_x, pc_y, sil0_1052, sil0_1053, \
                         sik_626, sik_845, sik_846, sil1_1052, sil1_1053, \
                         skk_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = pb_x[k] * sil0_1052[k]
                    + f_17 * sik_845[k]
                    - f_14 * pc_x[k] * sil1_1052[k];

        t_1053[k] = pb_x[k] * sil0_1053[k]
                    + f_17 * sik_846[k]
                    - f_14 * pc_x[k] * sil1_1053[k];

        t_1054[k] = f_18 * sik_626[k]
                    + f_3 * pc_y[k] * skk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pb_x, pc_x, pc_z, sil0_1055, sil0_1056, \
                         sik_591, sik_848, sik_849, sil1_1055, sil1_1056, \
                         skk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = pb_x[k] * sil0_1055[k]
                    + f_17 * sik_848[k]
                    - f_14 * pc_x[k] * sil1_1055[k];

        t_1056[k] = pb_x[k] * sil0_1056[k]
                    + f_16 * sik_849[k]
                    - f_14 * pc_x[k] * sil1_1056[k];

        t_1057[k] = f_16 * sik_591[k]
                    + f_3 * pc_z[k] * skk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pb_x, pc_x, sil0_1058, sil0_1059, sil0_1060, \
                         sik_851, sik_852, sik_853, sil1_1058, sil1_1059, \
                         sil1_1060 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = pb_x[k] * sil0_1058[k]
                    + f_16 * sik_851[k]
                    - f_14 * pc_x[k] * sil1_1058[k];

        t_1059[k] = pb_x[k] * sil0_1059[k]
                    + f_16 * sik_852[k]
                    - f_14 * pc_x[k] * sil1_1059[k];

        t_1060[k] = pb_x[k] * sil0_1060[k]
                    + f_16 * sik_853[k]
                    - f_14 * pc_x[k] * sil1_1060[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pb_x, pc_x, pc_y, sil0_1062, sik_632, \
                         sik_855, sik_856, sik_857, sil1_1062, skk_848, skk_856, \
                         skk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * sik_632[k]
                    + f_3 * pc_y[k] * skk_848[k];

        t_1062[k] = pb_x[k] * sil0_1062[k]
                    + f_16 * sik_855[k]
                    - f_14 * pc_x[k] * sil1_1062[k];

        t_1063[k] = f_15 * sik_856[k]
                    + f_3 * pc_x[k] * skk_856[k];

        t_1064[k] = f_15 * sik_857[k]
                    + f_3 * pc_x[k] * skk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, sik_858, sik_859, \
                         sik_860, sik_861, sik_862, skk_858, skk_859, skk_860, skk_861, \
                         skk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_15 * sik_858[k]
                    + f_3 * pc_x[k] * skk_858[k];

        t_1066[k] = f_15 * sik_859[k]
                    + f_3 * pc_x[k] * skk_859[k];

        t_1067[k] = f_15 * sik_860[k]
                    + f_3 * pc_x[k] * skk_860[k];

        t_1068[k] = f_15 * sik_861[k]
                    + f_3 * pc_x[k] * skk_861[k];

        t_1069[k] = f_15 * sik_862[k]
                    + f_3 * pc_x[k] * skk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, t_1073, pb_x, pc_x, pc_z, sil0_1071, \
                         sil0_1073, sik_604, sik_863, sil1_1071, sil1_1073, skk_856, \
                         skk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_15 * sik_863[k]
                    + f_3 * pc_x[k] * skk_863[k];

        t_1071[k] = pb_x[k] * sil0_1071[k]
                    - f_14 * pc_x[k] * sil1_1071[k];

        t_1072[k] = f_16 * sik_604[k]
                    + f_3 * pc_z[k] * skk_856[k];

        t_1073[k] = pb_x[k] * sil0_1073[k]
                    - f_14 * pc_x[k] * sil1_1073[k];
    }

#pragma omp simd aligned(t_1074, t_1075, t_1076, t_1077, pb_x, pc_x, sil0_1074, sil0_1075, \
                         sil0_1076, sil0_1077, sil1_1074, sil1_1075, sil1_1076, \
                         sil1_1077 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1074[k] = pb_x[k] * sil0_1074[k]
                    - f_14 * pc_x[k] * sil1_1074[k];

        t_1075[k] = pb_x[k] * sil0_1075[k]
                    - f_14 * pc_x[k] * sil1_1075[k];

        t_1076[k] = pb_x[k] * sil0_1076[k]
                    - f_14 * pc_x[k] * sil1_1076[k];

        t_1077[k] = pb_x[k] * sil0_1077[k]
                    - f_14 * pc_x[k] * sil1_1077[k];
    }

#pragma omp simd aligned(t_1078, t_1079, t_1080, t_1081, pb_x, pc_x, pc_y, sil0_1079, \
                         sil0_1080, sik_647, sik_648, sik_864, sil1_1079, sil1_1080, skk_863, \
                         skk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1078[k] = f_18 * sik_647[k]
                    + f_3 * pc_y[k] * skk_863[k];

        t_1079[k] = pb_x[k] * sil0_1079[k]
                    - f_14 * pc_x[k] * sil1_1079[k];

        t_1080[k] = pb_x[k] * sil0_1080[k]
                    + f_21 * sik_864[k]
                    - f_14 * pc_x[k] * sil1_1080[k];

        t_1081[k] = f_17 * sik_648[k]
                    + f_3 * pc_y[k] * skk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pb_x, pc_x, pc_y, pc_z, sil0_1083, sik_612, \
                         sik_650, sik_867, sil1_1083, skk_864, \
                         skk_866 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_17 * sik_612[k]
                    + f_3 * pc_z[k] * skk_864[k];

        t_1083[k] = pb_x[k] * sil0_1083[k]
                    + f_20 * sik_867[k]
                    - f_14 * pc_x[k] * sil1_1083[k];

        t_1084[k] = f_17 * sik_650[k]
                    + f_3 * pc_y[k] * skk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pb_x, pc_x, pc_z, sil0_1085, sil0_1086, \
                         sik_615, sik_869, sik_870, sil1_1085, sil1_1086, \
                         skk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = pb_x[k] * sil0_1085[k]
                    + f_20 * sik_869[k]
                    - f_14 * pc_x[k] * sil1_1085[k];

        t_1086[k] = pb_x[k] * sil0_1086[k]
                    + f_19 * sik_870[k]
                    - f_14 * pc_x[k] * sil1_1086[k];

        t_1087[k] = f_17 * sik_615[k]
                    + f_3 * pc_z[k] * skk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pb_x, pc_x, pc_y, sil0_1089, sil0_1090, \
                         sik_653, sik_873, sik_874, sil1_1089, sil1_1090, \
                         skk_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * sik_653[k]
                    + f_3 * pc_y[k] * skk_869[k];

        t_1089[k] = pb_x[k] * sil0_1089[k]
                    + f_19 * sik_873[k]
                    - f_14 * pc_x[k] * sil1_1089[k];

        t_1090[k] = pb_x[k] * sil0_1090[k]
                    + f_18 * sik_874[k]
                    - f_14 * pc_x[k] * sil1_1090[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pb_x, pc_x, pc_y, pc_z, sil0_1092, sik_618, \
                         sik_657, sik_876, sil1_1092, skk_870, \
                         skk_873 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_17 * sik_618[k]
                    + f_3 * pc_z[k] * skk_870[k];

        t_1092[k] = pb_x[k] * sil0_1092[k]
                    + f_18 * sik_876[k]
                    - f_14 * pc_x[k] * sil1_1092[k];

        t_1093[k] = f_17 * sik_657[k]
                    + f_3 * pc_y[k] * skk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pb_x, pc_x, pc_z, sil0_1094, sil0_1095, \
                         sik_622, sik_878, sik_879, sil1_1094, sil1_1095, \
                         skk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = pb_x[k] * sil0_1094[k]
                    + f_18 * sik_878[k]
                    - f_14 * pc_x[k] * sil1_1094[k];

        t_1095[k] = pb_x[k] * sil0_1095[k]
                    + f_17 * sik_879[k]
                    - f_14 * pc_x[k] * sil1_1095[k];

        t_1096[k] = f_17 * sik_622[k]
                    + f_3 * pc_z[k] * skk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pb_x, pc_x, pc_y, sil0_1097, sil0_1098, \
                         sik_662, sik_881, sik_882, sil1_1097, sil1_1098, \
                         skk_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = pb_x[k] * sil0_1097[k]
                    + f_17 * sik_881[k]
                    - f_14 * pc_x[k] * sil1_1097[k];

        t_1098[k] = pb_x[k] * sil0_1098[k]
                    + f_17 * sik_882[k]
                    - f_14 * pc_x[k] * sil1_1098[k];

        t_1099[k] = f_17 * sik_662[k]
                    + f_3 * pc_y[k] * skk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pb_x, pc_x, pc_z, sil0_1100, sil0_1101, \
                         sik_627, sik_884, sik_885, sil1_1100, sil1_1101, \
                         skk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = pb_x[k] * sil0_1100[k]
                    + f_17 * sik_884[k]
                    - f_14 * pc_x[k] * sil1_1100[k];

        t_1101[k] = pb_x[k] * sil0_1101[k]
                    + f_16 * sik_885[k]
                    - f_14 * pc_x[k] * sil1_1101[k];

        t_1102[k] = f_17 * sik_627[k]
                    + f_3 * pc_z[k] * skk_879[k];
    }

#pragma omp simd aligned(t_1103, t_1104, t_1105, pb_x, pc_x, sil0_1103, sil0_1104, sil0_1105, \
                         sik_887, sik_888, sik_889, sil1_1103, sil1_1104, \
                         sil1_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = pb_x[k] * sil0_1103[k]
                    + f_16 * sik_887[k]
                    - f_14 * pc_x[k] * sil1_1103[k];

        t_1104[k] = pb_x[k] * sil0_1104[k]
                    + f_16 * sik_888[k]
                    - f_14 * pc_x[k] * sil1_1104[k];

        t_1105[k] = pb_x[k] * sil0_1105[k]
                    + f_16 * sik_889[k]
                    - f_14 * pc_x[k] * sil1_1105[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pb_x, pc_x, pc_y, sil0_1107, sik_668, \
                         sik_891, sik_892, sik_893, sil1_1107, skk_884, skk_892, \
                         skk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_17 * sik_668[k]
                    + f_3 * pc_y[k] * skk_884[k];

        t_1107[k] = pb_x[k] * sil0_1107[k]
                    + f_16 * sik_891[k]
                    - f_14 * pc_x[k] * sil1_1107[k];

        t_1108[k] = f_15 * sik_892[k]
                    + f_3 * pc_x[k] * skk_892[k];

        t_1109[k] = f_15 * sik_893[k]
                    + f_3 * pc_x[k] * skk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, sik_894, sik_895, \
                         sik_896, sik_897, sik_898, skk_894, skk_895, skk_896, skk_897, \
                         skk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_15 * sik_894[k]
                    + f_3 * pc_x[k] * skk_894[k];

        t_1111[k] = f_15 * sik_895[k]
                    + f_3 * pc_x[k] * skk_895[k];

        t_1112[k] = f_15 * sik_896[k]
                    + f_3 * pc_x[k] * skk_896[k];

        t_1113[k] = f_15 * sik_897[k]
                    + f_3 * pc_x[k] * skk_897[k];

        t_1114[k] = f_15 * sik_898[k]
                    + f_3 * pc_x[k] * skk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, t_1118, pb_x, pc_x, pc_z, sil0_1116, \
                         sil0_1118, sik_640, sik_899, sil1_1116, sil1_1118, skk_892, \
                         skk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_15 * sik_899[k]
                    + f_3 * pc_x[k] * skk_899[k];

        t_1116[k] = pb_x[k] * sil0_1116[k]
                    - f_14 * pc_x[k] * sil1_1116[k];

        t_1117[k] = f_17 * sik_640[k]
                    + f_3 * pc_z[k] * skk_892[k];

        t_1118[k] = pb_x[k] * sil0_1118[k]
                    - f_14 * pc_x[k] * sil1_1118[k];
    }

#pragma omp simd aligned(t_1119, t_1120, t_1121, t_1122, pb_x, pc_x, sil0_1119, sil0_1120, \
                         sil0_1121, sil0_1122, sil1_1119, sil1_1120, sil1_1121, \
                         sil1_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1119[k] = pb_x[k] * sil0_1119[k]
                    - f_14 * pc_x[k] * sil1_1119[k];

        t_1120[k] = pb_x[k] * sil0_1120[k]
                    - f_14 * pc_x[k] * sil1_1120[k];

        t_1121[k] = pb_x[k] * sil0_1121[k]
                    - f_14 * pc_x[k] * sil1_1121[k];

        t_1122[k] = pb_x[k] * sil0_1122[k]
                    - f_14 * pc_x[k] * sil1_1122[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sil0,
                                                           const size_t sik, const size_t sil1,
                                                           const size_t skk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_900 = buffer.data(sil0 + 900);
    const auto *sil0_905 = buffer.data(sil0 + 905);
    const auto *sil0_909 = buffer.data(sil0 + 909);
    const auto *sil0_914 = buffer.data(sil0 + 914);
    const auto *sil0_920 = buffer.data(sil0 + 920);
    const auto *sil0_927 = buffer.data(sil0 + 927);
    const auto *sil0_1124 = buffer.data(sil0 + 1124);
    const auto *sil0_1125 = buffer.data(sil0 + 1125);
    const auto *sil0_1128 = buffer.data(sil0 + 1128);
    const auto *sil0_1130 = buffer.data(sil0 + 1130);
    const auto *sil0_1131 = buffer.data(sil0 + 1131);
    const auto *sil0_1134 = buffer.data(sil0 + 1134);
    const auto *sil0_1135 = buffer.data(sil0 + 1135);
    const auto *sil0_1137 = buffer.data(sil0 + 1137);
    const auto *sil0_1139 = buffer.data(sil0 + 1139);
    const auto *sil0_1140 = buffer.data(sil0 + 1140);
    const auto *sil0_1142 = buffer.data(sil0 + 1142);
    const auto *sil0_1143 = buffer.data(sil0 + 1143);
    const auto *sil0_1145 = buffer.data(sil0 + 1145);
    const auto *sil0_1146 = buffer.data(sil0 + 1146);
    const auto *sil0_1148 = buffer.data(sil0 + 1148);
    const auto *sil0_1149 = buffer.data(sil0 + 1149);
    const auto *sil0_1150 = buffer.data(sil0 + 1150);
    const auto *sil0_1152 = buffer.data(sil0 + 1152);
    const auto *sil0_1161 = buffer.data(sil0 + 1161);
    const auto *sil0_1163 = buffer.data(sil0 + 1163);
    const auto *sil0_1164 = buffer.data(sil0 + 1164);
    const auto *sil0_1165 = buffer.data(sil0 + 1165);
    const auto *sil0_1166 = buffer.data(sil0 + 1166);
    const auto *sil0_1167 = buffer.data(sil0 + 1167);
    const auto *sil0_1169 = buffer.data(sil0 + 1169);
    const auto *sil0_1173 = buffer.data(sil0 + 1173);
    const auto *sil0_1176 = buffer.data(sil0 + 1176);
    const auto *sil0_1180 = buffer.data(sil0 + 1180);
    const auto *sil0_1182 = buffer.data(sil0 + 1182);
    const auto *sil0_1185 = buffer.data(sil0 + 1185);
    const auto *sil0_1187 = buffer.data(sil0 + 1187);
    const auto *sil0_1188 = buffer.data(sil0 + 1188);
    const auto *sil0_1191 = buffer.data(sil0 + 1191);
    const auto *sil0_1193 = buffer.data(sil0 + 1193);
    const auto *sil0_1194 = buffer.data(sil0 + 1194);
    const auto *sil0_1195 = buffer.data(sil0 + 1195);
    const auto *sil0_1206 = buffer.data(sil0 + 1206);
    const auto *sil0_1208 = buffer.data(sil0 + 1208);
    const auto *sil0_1209 = buffer.data(sil0 + 1209);
    const auto *sil0_1210 = buffer.data(sil0 + 1210);
    const auto *sil0_1211 = buffer.data(sil0 + 1211);
    const auto *sil0_1212 = buffer.data(sil0 + 1212);
    const auto *sil0_1214 = buffer.data(sil0 + 1214);
    const auto *sil0_1215 = buffer.data(sil0 + 1215);
    const auto *sil0_1218 = buffer.data(sil0 + 1218);
    const auto *sil0_1220 = buffer.data(sil0 + 1220);
    const auto *sil0_1221 = buffer.data(sil0 + 1221);
    const auto *sil0_1224 = buffer.data(sil0 + 1224);
    const auto *sil0_1225 = buffer.data(sil0 + 1225);
    const auto *sil0_1227 = buffer.data(sil0 + 1227);
    const auto *sil0_1229 = buffer.data(sil0 + 1229);
    const auto *sil0_1230 = buffer.data(sil0 + 1230);
    const auto *sil0_1232 = buffer.data(sil0 + 1232);
    const auto *sil0_1233 = buffer.data(sil0 + 1233);
    const auto *sil0_1235 = buffer.data(sil0 + 1235);
    const auto *sil0_1236 = buffer.data(sil0 + 1236);

    const auto *sik_648 = buffer.data(sik + 648);
    const auto *sik_651 = buffer.data(sik + 651);
    const auto *sik_654 = buffer.data(sik + 654);
    const auto *sik_658 = buffer.data(sik + 658);
    const auto *sik_663 = buffer.data(sik + 663);
    const auto *sik_676 = buffer.data(sik + 676);
    const auto *sik_683 = buffer.data(sik + 683);
    const auto *sik_684 = buffer.data(sik + 684);
    const auto *sik_686 = buffer.data(sik + 686);
    const auto *sik_687 = buffer.data(sik + 687);
    const auto *sik_689 = buffer.data(sik + 689);
    const auto *sik_690 = buffer.data(sik + 690);
    const auto *sik_693 = buffer.data(sik + 693);
    const auto *sik_694 = buffer.data(sik + 694);
    const auto *sik_698 = buffer.data(sik + 698);
    const auto *sik_699 = buffer.data(sik + 699);
    const auto *sik_704 = buffer.data(sik + 704);
    const auto *sik_712 = buffer.data(sik + 712);
    const auto *sik_719 = buffer.data(sik + 719);
    const auto *sik_720 = buffer.data(sik + 720);
    const auto *sik_722 = buffer.data(sik + 722);
    const auto *sik_723 = buffer.data(sik + 723);
    const auto *sik_725 = buffer.data(sik + 725);
    const auto *sik_726 = buffer.data(sik + 726);
    const auto *sik_729 = buffer.data(sik + 729);
    const auto *sik_730 = buffer.data(sik + 730);
    const auto *sik_734 = buffer.data(sik + 734);
    const auto *sik_735 = buffer.data(sik + 735);
    const auto *sik_740 = buffer.data(sik + 740);
    const auto *sik_755 = buffer.data(sik + 755);
    const auto *sik_900 = buffer.data(sik + 900);
    const auto *sik_903 = buffer.data(sik + 903);
    const auto *sik_905 = buffer.data(sik + 905);
    const auto *sik_906 = buffer.data(sik + 906);
    const auto *sik_909 = buffer.data(sik + 909);
    const auto *sik_910 = buffer.data(sik + 910);
    const auto *sik_912 = buffer.data(sik + 912);
    const auto *sik_914 = buffer.data(sik + 914);
    const auto *sik_915 = buffer.data(sik + 915);
    const auto *sik_917 = buffer.data(sik + 917);
    const auto *sik_918 = buffer.data(sik + 918);
    const auto *sik_920 = buffer.data(sik + 920);
    const auto *sik_921 = buffer.data(sik + 921);
    const auto *sik_923 = buffer.data(sik + 923);
    const auto *sik_924 = buffer.data(sik + 924);
    const auto *sik_925 = buffer.data(sik + 925);
    const auto *sik_927 = buffer.data(sik + 927);
    const auto *sik_928 = buffer.data(sik + 928);
    const auto *sik_929 = buffer.data(sik + 929);
    const auto *sik_930 = buffer.data(sik + 930);
    const auto *sik_931 = buffer.data(sik + 931);
    const auto *sik_932 = buffer.data(sik + 932);
    const auto *sik_933 = buffer.data(sik + 933);
    const auto *sik_934 = buffer.data(sik + 934);
    const auto *sik_935 = buffer.data(sik + 935);
    const auto *sik_939 = buffer.data(sik + 939);
    const auto *sik_942 = buffer.data(sik + 942);
    const auto *sik_946 = buffer.data(sik + 946);
    const auto *sik_948 = buffer.data(sik + 948);
    const auto *sik_951 = buffer.data(sik + 951);
    const auto *sik_953 = buffer.data(sik + 953);
    const auto *sik_954 = buffer.data(sik + 954);
    const auto *sik_957 = buffer.data(sik + 957);
    const auto *sik_959 = buffer.data(sik + 959);
    const auto *sik_960 = buffer.data(sik + 960);
    const auto *sik_961 = buffer.data(sik + 961);
    const auto *sik_964 = buffer.data(sik + 964);
    const auto *sik_965 = buffer.data(sik + 965);
    const auto *sik_966 = buffer.data(sik + 966);
    const auto *sik_967 = buffer.data(sik + 967);
    const auto *sik_968 = buffer.data(sik + 968);
    const auto *sik_969 = buffer.data(sik + 969);
    const auto *sik_970 = buffer.data(sik + 970);
    const auto *sik_971 = buffer.data(sik + 971);
    const auto *sik_972 = buffer.data(sik + 972);
    const auto *sik_975 = buffer.data(sik + 975);
    const auto *sik_977 = buffer.data(sik + 977);
    const auto *sik_978 = buffer.data(sik + 978);
    const auto *sik_981 = buffer.data(sik + 981);
    const auto *sik_982 = buffer.data(sik + 982);
    const auto *sik_984 = buffer.data(sik + 984);
    const auto *sik_986 = buffer.data(sik + 986);
    const auto *sik_987 = buffer.data(sik + 987);
    const auto *sik_989 = buffer.data(sik + 989);
    const auto *sik_990 = buffer.data(sik + 990);
    const auto *sik_992 = buffer.data(sik + 992);
    const auto *sik_993 = buffer.data(sik + 993);

    const auto *sil1_900 = buffer.data(sil1 + 900);
    const auto *sil1_905 = buffer.data(sil1 + 905);
    const auto *sil1_909 = buffer.data(sil1 + 909);
    const auto *sil1_914 = buffer.data(sil1 + 914);
    const auto *sil1_920 = buffer.data(sil1 + 920);
    const auto *sil1_927 = buffer.data(sil1 + 927);
    const auto *sil1_1124 = buffer.data(sil1 + 1124);
    const auto *sil1_1125 = buffer.data(sil1 + 1125);
    const auto *sil1_1128 = buffer.data(sil1 + 1128);
    const auto *sil1_1130 = buffer.data(sil1 + 1130);
    const auto *sil1_1131 = buffer.data(sil1 + 1131);
    const auto *sil1_1134 = buffer.data(sil1 + 1134);
    const auto *sil1_1135 = buffer.data(sil1 + 1135);
    const auto *sil1_1137 = buffer.data(sil1 + 1137);
    const auto *sil1_1139 = buffer.data(sil1 + 1139);
    const auto *sil1_1140 = buffer.data(sil1 + 1140);
    const auto *sil1_1142 = buffer.data(sil1 + 1142);
    const auto *sil1_1143 = buffer.data(sil1 + 1143);
    const auto *sil1_1145 = buffer.data(sil1 + 1145);
    const auto *sil1_1146 = buffer.data(sil1 + 1146);
    const auto *sil1_1148 = buffer.data(sil1 + 1148);
    const auto *sil1_1149 = buffer.data(sil1 + 1149);
    const auto *sil1_1150 = buffer.data(sil1 + 1150);
    const auto *sil1_1152 = buffer.data(sil1 + 1152);
    const auto *sil1_1161 = buffer.data(sil1 + 1161);
    const auto *sil1_1163 = buffer.data(sil1 + 1163);
    const auto *sil1_1164 = buffer.data(sil1 + 1164);
    const auto *sil1_1165 = buffer.data(sil1 + 1165);
    const auto *sil1_1166 = buffer.data(sil1 + 1166);
    const auto *sil1_1167 = buffer.data(sil1 + 1167);
    const auto *sil1_1169 = buffer.data(sil1 + 1169);
    const auto *sil1_1173 = buffer.data(sil1 + 1173);
    const auto *sil1_1176 = buffer.data(sil1 + 1176);
    const auto *sil1_1180 = buffer.data(sil1 + 1180);
    const auto *sil1_1182 = buffer.data(sil1 + 1182);
    const auto *sil1_1185 = buffer.data(sil1 + 1185);
    const auto *sil1_1187 = buffer.data(sil1 + 1187);
    const auto *sil1_1188 = buffer.data(sil1 + 1188);
    const auto *sil1_1191 = buffer.data(sil1 + 1191);
    const auto *sil1_1193 = buffer.data(sil1 + 1193);
    const auto *sil1_1194 = buffer.data(sil1 + 1194);
    const auto *sil1_1195 = buffer.data(sil1 + 1195);
    const auto *sil1_1206 = buffer.data(sil1 + 1206);
    const auto *sil1_1208 = buffer.data(sil1 + 1208);
    const auto *sil1_1209 = buffer.data(sil1 + 1209);
    const auto *sil1_1210 = buffer.data(sil1 + 1210);
    const auto *sil1_1211 = buffer.data(sil1 + 1211);
    const auto *sil1_1212 = buffer.data(sil1 + 1212);
    const auto *sil1_1214 = buffer.data(sil1 + 1214);
    const auto *sil1_1215 = buffer.data(sil1 + 1215);
    const auto *sil1_1218 = buffer.data(sil1 + 1218);
    const auto *sil1_1220 = buffer.data(sil1 + 1220);
    const auto *sil1_1221 = buffer.data(sil1 + 1221);
    const auto *sil1_1224 = buffer.data(sil1 + 1224);
    const auto *sil1_1225 = buffer.data(sil1 + 1225);
    const auto *sil1_1227 = buffer.data(sil1 + 1227);
    const auto *sil1_1229 = buffer.data(sil1 + 1229);
    const auto *sil1_1230 = buffer.data(sil1 + 1230);
    const auto *sil1_1232 = buffer.data(sil1 + 1232);
    const auto *sil1_1233 = buffer.data(sil1 + 1233);
    const auto *sil1_1235 = buffer.data(sil1 + 1235);
    const auto *sil1_1236 = buffer.data(sil1 + 1236);

    const auto *skk_899 = buffer.data(skk + 899);
    const auto *skk_900 = buffer.data(skk + 900);
    const auto *skk_902 = buffer.data(skk + 902);
    const auto *skk_903 = buffer.data(skk + 903);
    const auto *skk_905 = buffer.data(skk + 905);
    const auto *skk_906 = buffer.data(skk + 906);
    const auto *skk_909 = buffer.data(skk + 909);
    const auto *skk_910 = buffer.data(skk + 910);
    const auto *skk_914 = buffer.data(skk + 914);
    const auto *skk_915 = buffer.data(skk + 915);
    const auto *skk_920 = buffer.data(skk + 920);
    const auto *skk_928 = buffer.data(skk + 928);
    const auto *skk_929 = buffer.data(skk + 929);
    const auto *skk_930 = buffer.data(skk + 930);
    const auto *skk_931 = buffer.data(skk + 931);
    const auto *skk_932 = buffer.data(skk + 932);
    const auto *skk_933 = buffer.data(skk + 933);
    const auto *skk_934 = buffer.data(skk + 934);
    const auto *skk_935 = buffer.data(skk + 935);
    const auto *skk_936 = buffer.data(skk + 936);
    const auto *skk_938 = buffer.data(skk + 938);
    const auto *skk_939 = buffer.data(skk + 939);
    const auto *skk_941 = buffer.data(skk + 941);
    const auto *skk_942 = buffer.data(skk + 942);
    const auto *skk_945 = buffer.data(skk + 945);
    const auto *skk_946 = buffer.data(skk + 946);
    const auto *skk_950 = buffer.data(skk + 950);
    const auto *skk_951 = buffer.data(skk + 951);
    const auto *skk_956 = buffer.data(skk + 956);
    const auto *skk_964 = buffer.data(skk + 964);
    const auto *skk_965 = buffer.data(skk + 965);
    const auto *skk_966 = buffer.data(skk + 966);
    const auto *skk_967 = buffer.data(skk + 967);
    const auto *skk_968 = buffer.data(skk + 968);
    const auto *skk_969 = buffer.data(skk + 969);
    const auto *skk_970 = buffer.data(skk + 970);
    const auto *skk_971 = buffer.data(skk + 971);
    const auto *skk_972 = buffer.data(skk + 972);
    const auto *skk_974 = buffer.data(skk + 974);
    const auto *skk_975 = buffer.data(skk + 975);
    const auto *skk_977 = buffer.data(skk + 977);
    const auto *skk_978 = buffer.data(skk + 978);
    const auto *skk_981 = buffer.data(skk + 981);
    const auto *skk_982 = buffer.data(skk + 982);
    const auto *skk_986 = buffer.data(skk + 986);
    const auto *skk_987 = buffer.data(skk + 987);

#pragma omp simd aligned(t_1123, t_1124, t_1125, t_1126, pb_x, pc_x, pc_y, sil0_1124, \
                         sil0_1125, sik_683, sik_684, sik_900, sil1_1124, sil1_1125, skk_899, \
                         skk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1123[k] = f_17 * sik_683[k]
                    + f_3 * pc_y[k] * skk_899[k];

        t_1124[k] = pb_x[k] * sil0_1124[k]
                    - f_14 * pc_x[k] * sil1_1124[k];

        t_1125[k] = pb_x[k] * sil0_1125[k]
                    + f_21 * sik_900[k]
                    - f_14 * pc_x[k] * sil1_1125[k];

        t_1126[k] = f_16 * sik_684[k]
                    + f_3 * pc_y[k] * skk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pb_x, pc_x, pc_y, pc_z, sil0_1128, sik_648, \
                         sik_686, sik_903, sil1_1128, skk_900, \
                         skk_902 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_18 * sik_648[k]
                    + f_3 * pc_z[k] * skk_900[k];

        t_1128[k] = pb_x[k] * sil0_1128[k]
                    + f_20 * sik_903[k]
                    - f_14 * pc_x[k] * sil1_1128[k];

        t_1129[k] = f_16 * sik_686[k]
                    + f_3 * pc_y[k] * skk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pb_x, pc_x, pc_z, sil0_1130, sil0_1131, \
                         sik_651, sik_905, sik_906, sil1_1130, sil1_1131, \
                         skk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = pb_x[k] * sil0_1130[k]
                    + f_20 * sik_905[k]
                    - f_14 * pc_x[k] * sil1_1130[k];

        t_1131[k] = pb_x[k] * sil0_1131[k]
                    + f_19 * sik_906[k]
                    - f_14 * pc_x[k] * sil1_1131[k];

        t_1132[k] = f_18 * sik_651[k]
                    + f_3 * pc_z[k] * skk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pb_x, pc_x, pc_y, sil0_1134, sil0_1135, \
                         sik_689, sik_909, sik_910, sil1_1134, sil1_1135, \
                         skk_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_16 * sik_689[k]
                    + f_3 * pc_y[k] * skk_905[k];

        t_1134[k] = pb_x[k] * sil0_1134[k]
                    + f_19 * sik_909[k]
                    - f_14 * pc_x[k] * sil1_1134[k];

        t_1135[k] = pb_x[k] * sil0_1135[k]
                    + f_18 * sik_910[k]
                    - f_14 * pc_x[k] * sil1_1135[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pb_x, pc_x, pc_y, pc_z, sil0_1137, sik_654, \
                         sik_693, sik_912, sil1_1137, skk_906, \
                         skk_909 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_18 * sik_654[k]
                    + f_3 * pc_z[k] * skk_906[k];

        t_1137[k] = pb_x[k] * sil0_1137[k]
                    + f_18 * sik_912[k]
                    - f_14 * pc_x[k] * sil1_1137[k];

        t_1138[k] = f_16 * sik_693[k]
                    + f_3 * pc_y[k] * skk_909[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, pb_x, pc_x, pc_z, sil0_1139, sil0_1140, \
                         sik_658, sik_914, sik_915, sil1_1139, sil1_1140, \
                         skk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = pb_x[k] * sil0_1139[k]
                    + f_18 * sik_914[k]
                    - f_14 * pc_x[k] * sil1_1139[k];

        t_1140[k] = pb_x[k] * sil0_1140[k]
                    + f_17 * sik_915[k]
                    - f_14 * pc_x[k] * sil1_1140[k];

        t_1141[k] = f_18 * sik_658[k]
                    + f_3 * pc_z[k] * skk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pb_x, pc_x, pc_y, sil0_1142, sil0_1143, \
                         sik_698, sik_917, sik_918, sil1_1142, sil1_1143, \
                         skk_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = pb_x[k] * sil0_1142[k]
                    + f_17 * sik_917[k]
                    - f_14 * pc_x[k] * sil1_1142[k];

        t_1143[k] = pb_x[k] * sil0_1143[k]
                    + f_17 * sik_918[k]
                    - f_14 * pc_x[k] * sil1_1143[k];

        t_1144[k] = f_16 * sik_698[k]
                    + f_3 * pc_y[k] * skk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pb_x, pc_x, pc_z, sil0_1145, sil0_1146, \
                         sik_663, sik_920, sik_921, sil1_1145, sil1_1146, \
                         skk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = pb_x[k] * sil0_1145[k]
                    + f_17 * sik_920[k]
                    - f_14 * pc_x[k] * sil1_1145[k];

        t_1146[k] = pb_x[k] * sil0_1146[k]
                    + f_16 * sik_921[k]
                    - f_14 * pc_x[k] * sil1_1146[k];

        t_1147[k] = f_18 * sik_663[k]
                    + f_3 * pc_z[k] * skk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pb_x, pc_x, sil0_1148, sil0_1149, sil0_1150, \
                         sik_923, sik_924, sik_925, sil1_1148, sil1_1149, \
                         sil1_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = pb_x[k] * sil0_1148[k]
                    + f_16 * sik_923[k]
                    - f_14 * pc_x[k] * sil1_1148[k];

        t_1149[k] = pb_x[k] * sil0_1149[k]
                    + f_16 * sik_924[k]
                    - f_14 * pc_x[k] * sil1_1149[k];

        t_1150[k] = pb_x[k] * sil0_1150[k]
                    + f_16 * sik_925[k]
                    - f_14 * pc_x[k] * sil1_1150[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pb_x, pc_x, pc_y, sil0_1152, sik_704, \
                         sik_927, sik_928, sik_929, sil1_1152, skk_920, skk_928, \
                         skk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_16 * sik_704[k]
                    + f_3 * pc_y[k] * skk_920[k];

        t_1152[k] = pb_x[k] * sil0_1152[k]
                    + f_16 * sik_927[k]
                    - f_14 * pc_x[k] * sil1_1152[k];

        t_1153[k] = f_15 * sik_928[k]
                    + f_3 * pc_x[k] * skk_928[k];

        t_1154[k] = f_15 * sik_929[k]
                    + f_3 * pc_x[k] * skk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, sik_930, sik_931, \
                         sik_932, sik_933, sik_934, skk_930, skk_931, skk_932, skk_933, \
                         skk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_15 * sik_930[k]
                    + f_3 * pc_x[k] * skk_930[k];

        t_1156[k] = f_15 * sik_931[k]
                    + f_3 * pc_x[k] * skk_931[k];

        t_1157[k] = f_15 * sik_932[k]
                    + f_3 * pc_x[k] * skk_932[k];

        t_1158[k] = f_15 * sik_933[k]
                    + f_3 * pc_x[k] * skk_933[k];

        t_1159[k] = f_15 * sik_934[k]
                    + f_3 * pc_x[k] * skk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, pb_x, pc_x, pc_z, sil0_1161, \
                         sil0_1163, sik_676, sik_935, sil1_1161, sil1_1163, skk_928, \
                         skk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_15 * sik_935[k]
                    + f_3 * pc_x[k] * skk_935[k];

        t_1161[k] = pb_x[k] * sil0_1161[k]
                    - f_14 * pc_x[k] * sil1_1161[k];

        t_1162[k] = f_18 * sik_676[k]
                    + f_3 * pc_z[k] * skk_928[k];

        t_1163[k] = pb_x[k] * sil0_1163[k]
                    - f_14 * pc_x[k] * sil1_1163[k];
    }

#pragma omp simd aligned(t_1164, t_1165, t_1166, t_1167, pb_x, pc_x, sil0_1164, sil0_1165, \
                         sil0_1166, sil0_1167, sil1_1164, sil1_1165, sil1_1166, \
                         sil1_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1164[k] = pb_x[k] * sil0_1164[k]
                    - f_14 * pc_x[k] * sil1_1164[k];

        t_1165[k] = pb_x[k] * sil0_1165[k]
                    - f_14 * pc_x[k] * sil1_1165[k];

        t_1166[k] = pb_x[k] * sil0_1166[k]
                    - f_14 * pc_x[k] * sil1_1166[k];

        t_1167[k] = pb_x[k] * sil0_1167[k]
                    - f_14 * pc_x[k] * sil1_1167[k];
    }

#pragma omp simd aligned(t_1168, t_1169, t_1170, t_1171, pb_x, pb_y, pc_x, pc_y, sil0_900, \
                         sil0_1169, sik_719, sik_720, sil1_900, sil1_1169, skk_935, \
                         skk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1168[k] = f_16 * sik_719[k]
                    + f_3 * pc_y[k] * skk_935[k];

        t_1169[k] = pb_x[k] * sil0_1169[k]
                    - f_14 * pc_x[k] * sil1_1169[k];

        t_1170[k] = pb_y[k] * sil0_900[k]
                    - f_14 * pc_y[k] * sil1_900[k];

        t_1171[k] = f_15 * sik_720[k]
                    + f_3 * pc_y[k] * skk_936[k];
    }

#pragma omp simd aligned(t_1172, t_1173, t_1174, pb_x, pc_x, pc_y, pc_z, sil0_1173, sik_684, \
                         sik_722, sik_939, sil1_1173, skk_936, \
                         skk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1172[k] = f_19 * sik_684[k]
                    + f_3 * pc_z[k] * skk_936[k];

        t_1173[k] = pb_x[k] * sil0_1173[k]
                    + f_20 * sik_939[k]
                    - f_14 * pc_x[k] * sil1_1173[k];

        t_1174[k] = f_15 * sik_722[k]
                    + f_3 * pc_y[k] * skk_938[k];
    }

#pragma omp simd aligned(t_1175, t_1176, t_1177, pb_x, pb_y, pc_x, pc_y, pc_z, sil0_905, \
                         sil0_1176, sik_687, sik_942, sil1_905, sil1_1176, \
                         skk_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1175[k] = pb_y[k] * sil0_905[k]
                    - f_14 * pc_y[k] * sil1_905[k];

        t_1176[k] = pb_x[k] * sil0_1176[k]
                    + f_19 * sik_942[k]
                    - f_14 * pc_x[k] * sil1_1176[k];

        t_1177[k] = f_19 * sik_687[k]
                    + f_3 * pc_z[k] * skk_939[k];
    }

#pragma omp simd aligned(t_1178, t_1179, t_1180, pb_x, pb_y, pc_x, pc_y, sil0_909, sil0_1180, \
                         sik_725, sik_946, sil1_909, sil1_1180, \
                         skk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1178[k] = f_15 * sik_725[k]
                    + f_3 * pc_y[k] * skk_941[k];

        t_1179[k] = pb_y[k] * sil0_909[k]
                    - f_14 * pc_y[k] * sil1_909[k];

        t_1180[k] = pb_x[k] * sil0_1180[k]
                    + f_18 * sik_946[k]
                    - f_14 * pc_x[k] * sil1_1180[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, pb_x, pc_x, pc_y, pc_z, sil0_1182, sik_690, \
                         sik_729, sik_948, sil1_1182, skk_942, \
                         skk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_19 * sik_690[k]
                    + f_3 * pc_z[k] * skk_942[k];

        t_1182[k] = pb_x[k] * sil0_1182[k]
                    + f_18 * sik_948[k]
                    - f_14 * pc_x[k] * sil1_1182[k];

        t_1183[k] = f_15 * sik_729[k]
                    + f_3 * pc_y[k] * skk_945[k];
    }

#pragma omp simd aligned(t_1184, t_1185, t_1186, pb_x, pb_y, pc_x, pc_y, pc_z, sil0_914, \
                         sil0_1185, sik_694, sik_951, sil1_914, sil1_1185, \
                         skk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1184[k] = pb_y[k] * sil0_914[k]
                    - f_14 * pc_y[k] * sil1_914[k];

        t_1185[k] = pb_x[k] * sil0_1185[k]
                    + f_17 * sik_951[k]
                    - f_14 * pc_x[k] * sil1_1185[k];

        t_1186[k] = f_19 * sik_694[k]
                    + f_3 * pc_z[k] * skk_946[k];
    }

#pragma omp simd aligned(t_1187, t_1188, t_1189, pb_x, pc_x, pc_y, sil0_1187, sil0_1188, \
                         sik_734, sik_953, sik_954, sil1_1187, sil1_1188, \
                         skk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1187[k] = pb_x[k] * sil0_1187[k]
                    + f_17 * sik_953[k]
                    - f_14 * pc_x[k] * sil1_1187[k];

        t_1188[k] = pb_x[k] * sil0_1188[k]
                    + f_17 * sik_954[k]
                    - f_14 * pc_x[k] * sil1_1188[k];

        t_1189[k] = f_15 * sik_734[k]
                    + f_3 * pc_y[k] * skk_950[k];
    }

#pragma omp simd aligned(t_1190, t_1191, t_1192, pb_x, pb_y, pc_x, pc_y, pc_z, sil0_920, \
                         sil0_1191, sik_699, sik_957, sil1_920, sil1_1191, \
                         skk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1190[k] = pb_y[k] * sil0_920[k]
                    - f_14 * pc_y[k] * sil1_920[k];

        t_1191[k] = pb_x[k] * sil0_1191[k]
                    + f_16 * sik_957[k]
                    - f_14 * pc_x[k] * sil1_1191[k];

        t_1192[k] = f_19 * sik_699[k]
                    + f_3 * pc_z[k] * skk_951[k];
    }

#pragma omp simd aligned(t_1193, t_1194, t_1195, pb_x, pc_x, sil0_1193, sil0_1194, sil0_1195, \
                         sik_959, sik_960, sik_961, sil1_1193, sil1_1194, \
                         sil1_1195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1193[k] = pb_x[k] * sil0_1193[k]
                    + f_16 * sik_959[k]
                    - f_14 * pc_x[k] * sil1_1193[k];

        t_1194[k] = pb_x[k] * sil0_1194[k]
                    + f_16 * sik_960[k]
                    - f_14 * pc_x[k] * sil1_1194[k];

        t_1195[k] = pb_x[k] * sil0_1195[k]
                    + f_16 * sik_961[k]
                    - f_14 * pc_x[k] * sil1_1195[k];
    }

#pragma omp simd aligned(t_1196, t_1197, t_1198, t_1199, pb_y, pc_x, pc_y, sil0_927, sik_740, \
                         sik_964, sik_965, sil1_927, skk_956, skk_964, \
                         skk_965 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1196[k] = f_15 * sik_740[k]
                    + f_3 * pc_y[k] * skk_956[k];

        t_1197[k] = pb_y[k] * sil0_927[k]
                    - f_14 * pc_y[k] * sil1_927[k];

        t_1198[k] = f_15 * sik_964[k]
                    + f_3 * pc_x[k] * skk_964[k];

        t_1199[k] = f_15 * sik_965[k]
                    + f_3 * pc_x[k] * skk_965[k];
    }

#pragma omp simd aligned(t_1200, t_1201, t_1202, t_1203, t_1204, pc_x, sik_966, sik_967, \
                         sik_968, sik_969, sik_970, skk_966, skk_967, skk_968, skk_969, \
                         skk_970 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1200[k] = f_15 * sik_966[k]
                    + f_3 * pc_x[k] * skk_966[k];

        t_1201[k] = f_15 * sik_967[k]
                    + f_3 * pc_x[k] * skk_967[k];

        t_1202[k] = f_15 * sik_968[k]
                    + f_3 * pc_x[k] * skk_968[k];

        t_1203[k] = f_15 * sik_969[k]
                    + f_3 * pc_x[k] * skk_969[k];

        t_1204[k] = f_15 * sik_970[k]
                    + f_3 * pc_x[k] * skk_970[k];
    }

#pragma omp simd aligned(t_1205, t_1206, t_1207, t_1208, pb_x, pc_x, pc_z, sil0_1206, \
                         sil0_1208, sik_712, sik_971, sil1_1206, sil1_1208, skk_964, \
                         skk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1205[k] = f_15 * sik_971[k]
                    + f_3 * pc_x[k] * skk_971[k];

        t_1206[k] = pb_x[k] * sil0_1206[k]
                    - f_14 * pc_x[k] * sil1_1206[k];

        t_1207[k] = f_19 * sik_712[k]
                    + f_3 * pc_z[k] * skk_964[k];

        t_1208[k] = pb_x[k] * sil0_1208[k]
                    - f_14 * pc_x[k] * sil1_1208[k];
    }

#pragma omp simd aligned(t_1209, t_1210, t_1211, t_1212, pb_x, pc_x, sil0_1209, sil0_1210, \
                         sil0_1211, sil0_1212, sil1_1209, sil1_1210, sil1_1211, \
                         sil1_1212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1209[k] = pb_x[k] * sil0_1209[k]
                    - f_14 * pc_x[k] * sil1_1209[k];

        t_1210[k] = pb_x[k] * sil0_1210[k]
                    - f_14 * pc_x[k] * sil1_1210[k];

        t_1211[k] = pb_x[k] * sil0_1211[k]
                    - f_14 * pc_x[k] * sil1_1211[k];

        t_1212[k] = pb_x[k] * sil0_1212[k]
                    - f_14 * pc_x[k] * sil1_1212[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, t_1216, pb_x, pc_x, pc_y, sil0_1214, \
                         sil0_1215, sik_755, sik_972, sil1_1214, sil1_1215, skk_971, \
                         skk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = f_15 * sik_755[k]
                    + f_3 * pc_y[k] * skk_971[k];

        t_1214[k] = pb_x[k] * sil0_1214[k]
                    - f_14 * pc_x[k] * sil1_1214[k];

        t_1215[k] = pb_x[k] * sil0_1215[k]
                    + f_21 * sik_972[k]
                    - f_14 * pc_x[k] * sil1_1215[k];

        t_1216[k] = f_3 * pc_y[k] * skk_972[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, pb_x, pc_x, pc_y, pc_z, sil0_1218, sik_720, \
                         sik_975, sil1_1218, skk_972, skk_974 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_20 * sik_720[k]
                    + f_3 * pc_z[k] * skk_972[k];

        t_1218[k] = pb_x[k] * sil0_1218[k]
                    + f_20 * sik_975[k]
                    - f_14 * pc_x[k] * sil1_1218[k];

        t_1219[k] = f_3 * pc_y[k] * skk_974[k];
    }

#pragma omp simd aligned(t_1220, t_1221, t_1222, pb_x, pc_x, pc_z, sil0_1220, sil0_1221, \
                         sik_723, sik_977, sik_978, sil1_1220, sil1_1221, \
                         skk_975 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1220[k] = pb_x[k] * sil0_1220[k]
                    + f_20 * sik_977[k]
                    - f_14 * pc_x[k] * sil1_1220[k];

        t_1221[k] = pb_x[k] * sil0_1221[k]
                    + f_19 * sik_978[k]
                    - f_14 * pc_x[k] * sil1_1221[k];

        t_1222[k] = f_20 * sik_723[k]
                    + f_3 * pc_z[k] * skk_975[k];
    }

#pragma omp simd aligned(t_1223, t_1224, t_1225, pb_x, pc_x, pc_y, sil0_1224, sil0_1225, \
                         sik_981, sik_982, sil1_1224, sil1_1225, \
                         skk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1223[k] = f_3 * pc_y[k] * skk_977[k];

        t_1224[k] = pb_x[k] * sil0_1224[k]
                    + f_19 * sik_981[k]
                    - f_14 * pc_x[k] * sil1_1224[k];

        t_1225[k] = pb_x[k] * sil0_1225[k]
                    + f_18 * sik_982[k]
                    - f_14 * pc_x[k] * sil1_1225[k];
    }

#pragma omp simd aligned(t_1226, t_1227, t_1228, pb_x, pc_x, pc_y, pc_z, sil0_1227, sik_726, \
                         sik_984, sil1_1227, skk_978, skk_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1226[k] = f_20 * sik_726[k]
                    + f_3 * pc_z[k] * skk_978[k];

        t_1227[k] = pb_x[k] * sil0_1227[k]
                    + f_18 * sik_984[k]
                    - f_14 * pc_x[k] * sil1_1227[k];

        t_1228[k] = f_3 * pc_y[k] * skk_981[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, pb_x, pc_x, pc_z, sil0_1229, sil0_1230, \
                         sik_730, sik_986, sik_987, sil1_1229, sil1_1230, \
                         skk_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = pb_x[k] * sil0_1229[k]
                    + f_18 * sik_986[k]
                    - f_14 * pc_x[k] * sil1_1229[k];

        t_1230[k] = pb_x[k] * sil0_1230[k]
                    + f_17 * sik_987[k]
                    - f_14 * pc_x[k] * sil1_1230[k];

        t_1231[k] = f_20 * sik_730[k]
                    + f_3 * pc_z[k] * skk_982[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, pb_x, pc_x, pc_y, sil0_1232, sil0_1233, \
                         sik_989, sik_990, sil1_1232, sil1_1233, \
                         skk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = pb_x[k] * sil0_1232[k]
                    + f_17 * sik_989[k]
                    - f_14 * pc_x[k] * sil1_1232[k];

        t_1233[k] = pb_x[k] * sil0_1233[k]
                    + f_17 * sik_990[k]
                    - f_14 * pc_x[k] * sil1_1233[k];

        t_1234[k] = f_3 * pc_y[k] * skk_986[k];
    }

#pragma omp simd aligned(t_1235, t_1236, t_1237, pb_x, pc_x, pc_z, sil0_1235, sil0_1236, \
                         sik_735, sik_992, sik_993, sil1_1235, sil1_1236, \
                         skk_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1235[k] = pb_x[k] * sil0_1235[k]
                    + f_17 * sik_992[k]
                    - f_14 * pc_x[k] * sil1_1235[k];

        t_1236[k] = pb_x[k] * sil0_1236[k]
                    + f_16 * sik_993[k]
                    - f_14 * pc_x[k] * sil1_1236[k];

        t_1237[k] = f_20 * sik_735[k]
                    + f_3 * pc_z[k] * skk_987[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sil0,
                                                           const size_t sik, const size_t sil1,
                                                           const size_t ski0, const size_t ski1,
                                                           const size_t skk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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
    auto *t_1260 = buffer.data(target + 1260);
    auto *t_1261 = buffer.data(target + 1261);
    auto *t_1262 = buffer.data(target + 1262);
    auto *t_1263 = buffer.data(target + 1263);
    auto *t_1264 = buffer.data(target + 1264);
    auto *t_1265 = buffer.data(target + 1265);
    auto *t_1266 = buffer.data(target + 1266);
    auto *t_1267 = buffer.data(target + 1267);
    auto *t_1268 = buffer.data(target + 1268);
    auto *t_1269 = buffer.data(target + 1269);
    auto *t_1270 = buffer.data(target + 1270);
    auto *t_1271 = buffer.data(target + 1271);
    auto *t_1272 = buffer.data(target + 1272);
    auto *t_1273 = buffer.data(target + 1273);
    auto *t_1274 = buffer.data(target + 1274);
    auto *t_1275 = buffer.data(target + 1275);
    auto *t_1276 = buffer.data(target + 1276);
    auto *t_1277 = buffer.data(target + 1277);
    auto *t_1278 = buffer.data(target + 1278);
    auto *t_1279 = buffer.data(target + 1279);
    auto *t_1280 = buffer.data(target + 1280);
    auto *t_1281 = buffer.data(target + 1281);
    auto *t_1282 = buffer.data(target + 1282);
    auto *t_1283 = buffer.data(target + 1283);
    auto *t_1284 = buffer.data(target + 1284);
    auto *t_1285 = buffer.data(target + 1285);
    auto *t_1286 = buffer.data(target + 1286);
    auto *t_1287 = buffer.data(target + 1287);
    auto *t_1288 = buffer.data(target + 1288);
    auto *t_1289 = buffer.data(target + 1289);
    auto *t_1290 = buffer.data(target + 1290);
    auto *t_1291 = buffer.data(target + 1291);
    auto *t_1292 = buffer.data(target + 1292);
    auto *t_1293 = buffer.data(target + 1293);
    auto *t_1294 = buffer.data(target + 1294);
    auto *t_1295 = buffer.data(target + 1295);
    auto *t_1296 = buffer.data(target + 1296);
    auto *t_1297 = buffer.data(target + 1297);
    auto *t_1298 = buffer.data(target + 1298);
    auto *t_1299 = buffer.data(target + 1299);
    auto *t_1300 = buffer.data(target + 1300);
    auto *t_1301 = buffer.data(target + 1301);
    auto *t_1302 = buffer.data(target + 1302);
    auto *t_1303 = buffer.data(target + 1303);
    auto *t_1304 = buffer.data(target + 1304);
    auto *t_1305 = buffer.data(target + 1305);
    auto *t_1306 = buffer.data(target + 1306);
    auto *t_1307 = buffer.data(target + 1307);
    auto *t_1308 = buffer.data(target + 1308);
    auto *t_1309 = buffer.data(target + 1309);
    auto *t_1310 = buffer.data(target + 1310);
    auto *t_1311 = buffer.data(target + 1311);
    auto *t_1312 = buffer.data(target + 1312);
    auto *t_1313 = buffer.data(target + 1313);
    auto *t_1314 = buffer.data(target + 1314);
    auto *t_1315 = buffer.data(target + 1315);
    auto *t_1316 = buffer.data(target + 1316);
    auto *t_1317 = buffer.data(target + 1317);
    auto *t_1318 = buffer.data(target + 1318);
    auto *t_1319 = buffer.data(target + 1319);
    auto *t_1320 = buffer.data(target + 1320);
    auto *t_1321 = buffer.data(target + 1321);
    auto *t_1322 = buffer.data(target + 1322);
    auto *t_1323 = buffer.data(target + 1323);
    auto *t_1324 = buffer.data(target + 1324);
    auto *t_1325 = buffer.data(target + 1325);
    auto *t_1326 = buffer.data(target + 1326);
    auto *t_1327 = buffer.data(target + 1327);
    auto *t_1328 = buffer.data(target + 1328);
    auto *t_1329 = buffer.data(target + 1329);
    auto *t_1330 = buffer.data(target + 1330);
    auto *t_1331 = buffer.data(target + 1331);
    auto *t_1332 = buffer.data(target + 1332);
    auto *t_1333 = buffer.data(target + 1333);
    auto *t_1334 = buffer.data(target + 1334);
    auto *t_1335 = buffer.data(target + 1335);
    auto *t_1336 = buffer.data(target + 1336);
    auto *t_1337 = buffer.data(target + 1337);
    auto *t_1338 = buffer.data(target + 1338);
    auto *t_1339 = buffer.data(target + 1339);
    auto *t_1340 = buffer.data(target + 1340);
    auto *t_1341 = buffer.data(target + 1341);
    auto *t_1342 = buffer.data(target + 1342);
    auto *t_1343 = buffer.data(target + 1343);
    auto *t_1344 = buffer.data(target + 1344);
    auto *t_1345 = buffer.data(target + 1345);
    auto *t_1346 = buffer.data(target + 1346);
    auto *t_1347 = buffer.data(target + 1347);
    auto *t_1348 = buffer.data(target + 1348);
    auto *t_1349 = buffer.data(target + 1349);
    auto *t_1350 = buffer.data(target + 1350);
    auto *t_1351 = buffer.data(target + 1351);
    auto *t_1352 = buffer.data(target + 1352);
    auto *t_1353 = buffer.data(target + 1353);
    auto *t_1354 = buffer.data(target + 1354);
    auto *t_1355 = buffer.data(target + 1355);
    auto *t_1356 = buffer.data(target + 1356);
    auto *t_1357 = buffer.data(target + 1357);
    auto *t_1358 = buffer.data(target + 1358);
    auto *t_1359 = buffer.data(target + 1359);
    auto *t_1360 = buffer.data(target + 1360);
    auto *t_1361 = buffer.data(target + 1361);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_945 = buffer.data(sil0 + 945);
    const auto *sil0_948 = buffer.data(sil0 + 948);
    const auto *sil0_951 = buffer.data(sil0 + 951);
    const auto *sil0_955 = buffer.data(sil0 + 955);
    const auto *sil0_960 = buffer.data(sil0 + 960);
    const auto *sil0_966 = buffer.data(sil0 + 966);
    const auto *sil0_981 = buffer.data(sil0 + 981);
    const auto *sil0_983 = buffer.data(sil0 + 983);
    const auto *sil0_984 = buffer.data(sil0 + 984);
    const auto *sil0_985 = buffer.data(sil0 + 985);
    const auto *sil0_986 = buffer.data(sil0 + 986);
    const auto *sil0_987 = buffer.data(sil0 + 987);
    const auto *sil0_1238 = buffer.data(sil0 + 1238);
    const auto *sil0_1239 = buffer.data(sil0 + 1239);
    const auto *sil0_1240 = buffer.data(sil0 + 1240);
    const auto *sil0_1242 = buffer.data(sil0 + 1242);
    const auto *sil0_1251 = buffer.data(sil0 + 1251);
    const auto *sil0_1253 = buffer.data(sil0 + 1253);
    const auto *sil0_1254 = buffer.data(sil0 + 1254);
    const auto *sil0_1255 = buffer.data(sil0 + 1255);
    const auto *sil0_1256 = buffer.data(sil0 + 1256);
    const auto *sil0_1257 = buffer.data(sil0 + 1257);
    const auto *sil0_1259 = buffer.data(sil0 + 1259);

    const auto *sik_748 = buffer.data(sik + 748);
    const auto *sik_756 = buffer.data(sik + 756);
    const auto *sik_758 = buffer.data(sik + 758);
    const auto *sik_759 = buffer.data(sik + 759);
    const auto *sik_761 = buffer.data(sik + 761);
    const auto *sik_762 = buffer.data(sik + 762);
    const auto *sik_765 = buffer.data(sik + 765);
    const auto *sik_766 = buffer.data(sik + 766);
    const auto *sik_770 = buffer.data(sik + 770);
    const auto *sik_771 = buffer.data(sik + 771);
    const auto *sik_776 = buffer.data(sik + 776);
    const auto *sik_784 = buffer.data(sik + 784);
    const auto *sik_785 = buffer.data(sik + 785);
    const auto *sik_786 = buffer.data(sik + 786);
    const auto *sik_787 = buffer.data(sik + 787);
    const auto *sik_788 = buffer.data(sik + 788);
    const auto *sik_789 = buffer.data(sik + 789);
    const auto *sik_790 = buffer.data(sik + 790);
    const auto *sik_791 = buffer.data(sik + 791);
    const auto *sik_792 = buffer.data(sik + 792);
    const auto *sik_794 = buffer.data(sik + 794);
    const auto *sik_795 = buffer.data(sik + 795);
    const auto *sik_797 = buffer.data(sik + 797);
    const auto *sik_798 = buffer.data(sik + 798);
    const auto *sik_801 = buffer.data(sik + 801);
    const auto *sik_806 = buffer.data(sik + 806);
    const auto *sik_812 = buffer.data(sik + 812);
    const auto *sik_827 = buffer.data(sik + 827);
    const auto *sik_828 = buffer.data(sik + 828);
    const auto *sik_830 = buffer.data(sik + 830);
    const auto *sik_833 = buffer.data(sik + 833);
    const auto *sik_995 = buffer.data(sik + 995);
    const auto *sik_996 = buffer.data(sik + 996);
    const auto *sik_997 = buffer.data(sik + 997);
    const auto *sik_999 = buffer.data(sik + 999);
    const auto *sik_1000 = buffer.data(sik + 1000);
    const auto *sik_1001 = buffer.data(sik + 1001);
    const auto *sik_1002 = buffer.data(sik + 1002);
    const auto *sik_1003 = buffer.data(sik + 1003);
    const auto *sik_1004 = buffer.data(sik + 1004);
    const auto *sik_1005 = buffer.data(sik + 1005);
    const auto *sik_1006 = buffer.data(sik + 1006);
    const auto *sik_1007 = buffer.data(sik + 1007);

    const auto *sil1_945 = buffer.data(sil1 + 945);
    const auto *sil1_948 = buffer.data(sil1 + 948);
    const auto *sil1_951 = buffer.data(sil1 + 951);
    const auto *sil1_955 = buffer.data(sil1 + 955);
    const auto *sil1_960 = buffer.data(sil1 + 960);
    const auto *sil1_966 = buffer.data(sil1 + 966);
    const auto *sil1_981 = buffer.data(sil1 + 981);
    const auto *sil1_983 = buffer.data(sil1 + 983);
    const auto *sil1_984 = buffer.data(sil1 + 984);
    const auto *sil1_985 = buffer.data(sil1 + 985);
    const auto *sil1_986 = buffer.data(sil1 + 986);
    const auto *sil1_987 = buffer.data(sil1 + 987);
    const auto *sil1_1238 = buffer.data(sil1 + 1238);
    const auto *sil1_1239 = buffer.data(sil1 + 1239);
    const auto *sil1_1240 = buffer.data(sil1 + 1240);
    const auto *sil1_1242 = buffer.data(sil1 + 1242);
    const auto *sil1_1251 = buffer.data(sil1 + 1251);
    const auto *sil1_1253 = buffer.data(sil1 + 1253);
    const auto *sil1_1254 = buffer.data(sil1 + 1254);
    const auto *sil1_1255 = buffer.data(sil1 + 1255);
    const auto *sil1_1256 = buffer.data(sil1 + 1256);
    const auto *sil1_1257 = buffer.data(sil1 + 1257);
    const auto *sil1_1259 = buffer.data(sil1 + 1259);

    const auto *ski0_784 = buffer.data(ski0 + 784);
    const auto *ski0_787 = buffer.data(ski0 + 787);
    const auto *ski0_789 = buffer.data(ski0 + 789);
    const auto *ski0_790 = buffer.data(ski0 + 790);
    const auto *ski0_793 = buffer.data(ski0 + 793);
    const auto *ski0_794 = buffer.data(ski0 + 794);
    const auto *ski0_796 = buffer.data(ski0 + 796);
    const auto *ski0_798 = buffer.data(ski0 + 798);
    const auto *ski0_799 = buffer.data(ski0 + 799);
    const auto *ski0_801 = buffer.data(ski0 + 801);
    const auto *ski0_802 = buffer.data(ski0 + 802);
    const auto *ski0_804 = buffer.data(ski0 + 804);
    const auto *ski0_805 = buffer.data(ski0 + 805);
    const auto *ski0_807 = buffer.data(ski0 + 807);
    const auto *ski0_808 = buffer.data(ski0 + 808);
    const auto *ski0_809 = buffer.data(ski0 + 809);
    const auto *ski0_810 = buffer.data(ski0 + 810);
    const auto *ski0_811 = buffer.data(ski0 + 811);
    const auto *ski0_817 = buffer.data(ski0 + 817);
    const auto *ski0_821 = buffer.data(ski0 + 821);
    const auto *ski0_824 = buffer.data(ski0 + 824);
    const auto *ski0_826 = buffer.data(ski0 + 826);
    const auto *ski0_829 = buffer.data(ski0 + 829);
    const auto *ski0_830 = buffer.data(ski0 + 830);
    const auto *ski0_832 = buffer.data(ski0 + 832);
    const auto *ski0_835 = buffer.data(ski0 + 835);
    const auto *ski0_836 = buffer.data(ski0 + 836);
    const auto *ski0_837 = buffer.data(ski0 + 837);
    const auto *ski0_839 = buffer.data(ski0 + 839);
    const auto *ski0_840 = buffer.data(ski0 + 840);
    const auto *ski0_843 = buffer.data(ski0 + 843);
    const auto *ski0_845 = buffer.data(ski0 + 845);
    const auto *ski0_846 = buffer.data(ski0 + 846);
    const auto *ski0_849 = buffer.data(ski0 + 849);
    const auto *ski0_850 = buffer.data(ski0 + 850);

    const auto *ski1_784 = buffer.data(ski1 + 784);
    const auto *ski1_787 = buffer.data(ski1 + 787);
    const auto *ski1_789 = buffer.data(ski1 + 789);
    const auto *ski1_790 = buffer.data(ski1 + 790);
    const auto *ski1_793 = buffer.data(ski1 + 793);
    const auto *ski1_794 = buffer.data(ski1 + 794);
    const auto *ski1_796 = buffer.data(ski1 + 796);
    const auto *ski1_798 = buffer.data(ski1 + 798);
    const auto *ski1_799 = buffer.data(ski1 + 799);
    const auto *ski1_801 = buffer.data(ski1 + 801);
    const auto *ski1_802 = buffer.data(ski1 + 802);
    const auto *ski1_804 = buffer.data(ski1 + 804);
    const auto *ski1_805 = buffer.data(ski1 + 805);
    const auto *ski1_807 = buffer.data(ski1 + 807);
    const auto *ski1_808 = buffer.data(ski1 + 808);
    const auto *ski1_809 = buffer.data(ski1 + 809);
    const auto *ski1_810 = buffer.data(ski1 + 810);
    const auto *ski1_811 = buffer.data(ski1 + 811);
    const auto *ski1_817 = buffer.data(ski1 + 817);
    const auto *ski1_821 = buffer.data(ski1 + 821);
    const auto *ski1_824 = buffer.data(ski1 + 824);
    const auto *ski1_826 = buffer.data(ski1 + 826);
    const auto *ski1_829 = buffer.data(ski1 + 829);
    const auto *ski1_830 = buffer.data(ski1 + 830);
    const auto *ski1_832 = buffer.data(ski1 + 832);
    const auto *ski1_835 = buffer.data(ski1 + 835);
    const auto *ski1_836 = buffer.data(ski1 + 836);
    const auto *ski1_837 = buffer.data(ski1 + 837);
    const auto *ski1_839 = buffer.data(ski1 + 839);
    const auto *ski1_840 = buffer.data(ski1 + 840);
    const auto *ski1_843 = buffer.data(ski1 + 843);
    const auto *ski1_845 = buffer.data(ski1 + 845);
    const auto *ski1_846 = buffer.data(ski1 + 846);
    const auto *ski1_849 = buffer.data(ski1 + 849);
    const auto *ski1_850 = buffer.data(ski1 + 850);

    const auto *skk_992 = buffer.data(skk + 992);
    const auto *skk_1000 = buffer.data(skk + 1000);
    const auto *skk_1001 = buffer.data(skk + 1001);
    const auto *skk_1002 = buffer.data(skk + 1002);
    const auto *skk_1003 = buffer.data(skk + 1003);
    const auto *skk_1004 = buffer.data(skk + 1004);
    const auto *skk_1005 = buffer.data(skk + 1005);
    const auto *skk_1006 = buffer.data(skk + 1006);
    const auto *skk_1007 = buffer.data(skk + 1007);
    const auto *skk_1008 = buffer.data(skk + 1008);
    const auto *skk_1010 = buffer.data(skk + 1010);
    const auto *skk_1011 = buffer.data(skk + 1011);
    const auto *skk_1013 = buffer.data(skk + 1013);
    const auto *skk_1014 = buffer.data(skk + 1014);
    const auto *skk_1017 = buffer.data(skk + 1017);
    const auto *skk_1018 = buffer.data(skk + 1018);
    const auto *skk_1020 = buffer.data(skk + 1020);
    const auto *skk_1022 = buffer.data(skk + 1022);
    const auto *skk_1023 = buffer.data(skk + 1023);
    const auto *skk_1025 = buffer.data(skk + 1025);
    const auto *skk_1026 = buffer.data(skk + 1026);
    const auto *skk_1028 = buffer.data(skk + 1028);
    const auto *skk_1029 = buffer.data(skk + 1029);
    const auto *skk_1031 = buffer.data(skk + 1031);
    const auto *skk_1032 = buffer.data(skk + 1032);
    const auto *skk_1033 = buffer.data(skk + 1033);
    const auto *skk_1035 = buffer.data(skk + 1035);
    const auto *skk_1036 = buffer.data(skk + 1036);
    const auto *skk_1037 = buffer.data(skk + 1037);
    const auto *skk_1038 = buffer.data(skk + 1038);
    const auto *skk_1039 = buffer.data(skk + 1039);
    const auto *skk_1040 = buffer.data(skk + 1040);
    const auto *skk_1041 = buffer.data(skk + 1041);
    const auto *skk_1042 = buffer.data(skk + 1042);
    const auto *skk_1043 = buffer.data(skk + 1043);
    const auto *skk_1044 = buffer.data(skk + 1044);
    const auto *skk_1046 = buffer.data(skk + 1046);
    const auto *skk_1047 = buffer.data(skk + 1047);
    const auto *skk_1049 = buffer.data(skk + 1049);
    const auto *skk_1050 = buffer.data(skk + 1050);
    const auto *skk_1053 = buffer.data(skk + 1053);
    const auto *skk_1054 = buffer.data(skk + 1054);
    const auto *skk_1056 = buffer.data(skk + 1056);
    const auto *skk_1058 = buffer.data(skk + 1058);
    const auto *skk_1059 = buffer.data(skk + 1059);
    const auto *skk_1061 = buffer.data(skk + 1061);
    const auto *skk_1062 = buffer.data(skk + 1062);
    const auto *skk_1064 = buffer.data(skk + 1064);
    const auto *skk_1067 = buffer.data(skk + 1067);
    const auto *skk_1068 = buffer.data(skk + 1068);
    const auto *skk_1069 = buffer.data(skk + 1069);
    const auto *skk_1071 = buffer.data(skk + 1071);
    const auto *skk_1072 = buffer.data(skk + 1072);
    const auto *skk_1073 = buffer.data(skk + 1073);
    const auto *skk_1074 = buffer.data(skk + 1074);
    const auto *skk_1075 = buffer.data(skk + 1075);
    const auto *skk_1076 = buffer.data(skk + 1076);
    const auto *skk_1077 = buffer.data(skk + 1077);
    const auto *skk_1078 = buffer.data(skk + 1078);
    const auto *skk_1079 = buffer.data(skk + 1079);
    const auto *skk_1080 = buffer.data(skk + 1080);
    const auto *skk_1082 = buffer.data(skk + 1082);
    const auto *skk_1083 = buffer.data(skk + 1083);
    const auto *skk_1085 = buffer.data(skk + 1085);
    const auto *skk_1086 = buffer.data(skk + 1086);
    const auto *skk_1089 = buffer.data(skk + 1089);
    const auto *skk_1090 = buffer.data(skk + 1090);

#pragma omp simd aligned(t_1238, t_1239, t_1240, pb_x, pc_x, sil0_1238, sil0_1239, sil0_1240, \
                         sik_995, sik_996, sik_997, sil1_1238, sil1_1239, \
                         sil1_1240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1238[k] = pb_x[k] * sil0_1238[k]
                    + f_16 * sik_995[k]
                    - f_14 * pc_x[k] * sil1_1238[k];

        t_1239[k] = pb_x[k] * sil0_1239[k]
                    + f_16 * sik_996[k]
                    - f_14 * pc_x[k] * sil1_1239[k];

        t_1240[k] = pb_x[k] * sil0_1240[k]
                    + f_16 * sik_997[k]
                    - f_14 * pc_x[k] * sil1_1240[k];
    }

#pragma omp simd aligned(t_1241, t_1242, t_1243, t_1244, pb_x, pc_x, pc_y, sil0_1242, sik_999, \
                         sik_1000, sik_1001, sil1_1242, skk_992, skk_1000, \
                         skk_1001 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1241[k] = f_3 * pc_y[k] * skk_992[k];

        t_1242[k] = pb_x[k] * sil0_1242[k]
                    + f_16 * sik_999[k]
                    - f_14 * pc_x[k] * sil1_1242[k];

        t_1243[k] = f_15 * sik_1000[k]
                    + f_3 * pc_x[k] * skk_1000[k];

        t_1244[k] = f_15 * sik_1001[k]
                    + f_3 * pc_x[k] * skk_1001[k];
    }

#pragma omp simd aligned(t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, sik_1002, sik_1003, \
                         sik_1004, sik_1005, sik_1006, skk_1002, skk_1003, skk_1004, skk_1005, \
                         skk_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1245[k] = f_15 * sik_1002[k]
                    + f_3 * pc_x[k] * skk_1002[k];

        t_1246[k] = f_15 * sik_1003[k]
                    + f_3 * pc_x[k] * skk_1003[k];

        t_1247[k] = f_15 * sik_1004[k]
                    + f_3 * pc_x[k] * skk_1004[k];

        t_1248[k] = f_15 * sik_1005[k]
                    + f_3 * pc_x[k] * skk_1005[k];

        t_1249[k] = f_15 * sik_1006[k]
                    + f_3 * pc_x[k] * skk_1006[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pb_x, pc_x, pc_z, sil0_1251, \
                         sil0_1253, sik_748, sik_1007, sil1_1251, sil1_1253, skk_1000, \
                         skk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_15 * sik_1007[k]
                    + f_3 * pc_x[k] * skk_1007[k];

        t_1251[k] = pb_x[k] * sil0_1251[k]
                    - f_14 * pc_x[k] * sil1_1251[k];

        t_1252[k] = f_20 * sik_748[k]
                    + f_3 * pc_z[k] * skk_1000[k];

        t_1253[k] = pb_x[k] * sil0_1253[k]
                    - f_14 * pc_x[k] * sil1_1253[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, t_1257, pb_x, pc_x, sil0_1254, sil0_1255, \
                         sil0_1256, sil0_1257, sil1_1254, sil1_1255, sil1_1256, \
                         sil1_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = pb_x[k] * sil0_1254[k]
                    - f_14 * pc_x[k] * sil1_1254[k];

        t_1255[k] = pb_x[k] * sil0_1255[k]
                    - f_14 * pc_x[k] * sil1_1255[k];

        t_1256[k] = pb_x[k] * sil0_1256[k]
                    - f_14 * pc_x[k] * sil1_1256[k];

        t_1257[k] = pb_x[k] * sil0_1257[k]
                    - f_14 * pc_x[k] * sil1_1257[k];
    }

#pragma omp simd aligned(t_1258, t_1259, t_1260, t_1261, t_1262, pb_x, pc_x, pc_y, pc_z, \
                         sil0_1259, sik_756, sil1_1259, ski0_784, ski1_784, skk_1007, \
                         skk_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1258[k] = f_3 * pc_y[k] * skk_1007[k];

        t_1259[k] = pb_x[k] * sil0_1259[k]
                    - f_14 * pc_x[k] * sil1_1259[k];

        t_1260[k] = f_1 * ski0_784[k]
                    - f_2 * ski1_784[k]
                    + f_3 * pc_x[k] * skk_1008[k];

        t_1261[k] = f_0 * sik_756[k]
                    + f_3 * pc_y[k] * skk_1008[k];

        t_1262[k] = f_3 * pc_z[k] * skk_1008[k];
    }

#pragma omp simd aligned(t_1263, t_1264, t_1265, pc_x, pc_y, sik_758, ski0_787, ski0_789, \
                         ski1_787, ski1_789, skk_1010, skk_1011, \
                         skk_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1263[k] = f_4 * ski0_787[k]
                    - f_5 * ski1_787[k]
                    + f_3 * pc_x[k] * skk_1011[k];

        t_1264[k] = f_0 * sik_758[k]
                    + f_3 * pc_y[k] * skk_1010[k];

        t_1265[k] = f_4 * ski0_789[k]
                    - f_5 * ski1_789[k]
                    + f_3 * pc_x[k] * skk_1013[k];
    }

#pragma omp simd aligned(t_1266, t_1267, t_1268, t_1269, pc_x, pc_y, pc_z, sik_761, ski0_790, \
                         ski0_793, ski1_790, ski1_793, skk_1011, skk_1013, skk_1014, \
                         skk_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1266[k] = f_6 * ski0_790[k]
                    - f_7 * ski1_790[k]
                    + f_3 * pc_x[k] * skk_1014[k];

        t_1267[k] = f_3 * pc_z[k] * skk_1011[k];

        t_1268[k] = f_0 * sik_761[k]
                    + f_3 * pc_y[k] * skk_1013[k];

        t_1269[k] = f_6 * ski0_793[k]
                    - f_7 * ski1_793[k]
                    + f_3 * pc_x[k] * skk_1017[k];
    }

#pragma omp simd aligned(t_1270, t_1271, t_1272, t_1273, pc_x, pc_y, pc_z, sik_765, ski0_794, \
                         ski0_796, ski1_794, ski1_796, skk_1014, skk_1017, skk_1018, \
                         skk_1020 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1270[k] = f_8 * ski0_794[k]
                    - f_9 * ski1_794[k]
                    + f_3 * pc_x[k] * skk_1018[k];

        t_1271[k] = f_3 * pc_z[k] * skk_1014[k];

        t_1272[k] = f_8 * ski0_796[k]
                    - f_9 * ski1_796[k]
                    + f_3 * pc_x[k] * skk_1020[k];

        t_1273[k] = f_0 * sik_765[k]
                    + f_3 * pc_y[k] * skk_1017[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, t_1277, pc_x, pc_z, ski0_798, ski0_799, \
                         ski0_801, ski1_798, ski1_799, ski1_801, skk_1018, skk_1022, skk_1023, \
                         skk_1025 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = f_8 * ski0_798[k]
                    - f_9 * ski1_798[k]
                    + f_3 * pc_x[k] * skk_1022[k];

        t_1275[k] = f_10 * ski0_799[k]
                    - f_11 * ski1_799[k]
                    + f_3 * pc_x[k] * skk_1023[k];

        t_1276[k] = f_3 * pc_z[k] * skk_1018[k];

        t_1277[k] = f_10 * ski0_801[k]
                    - f_11 * ski1_801[k]
                    + f_3 * pc_x[k] * skk_1025[k];
    }

#pragma omp simd aligned(t_1278, t_1279, t_1280, pc_x, pc_y, sik_770, ski0_802, ski0_804, \
                         ski1_802, ski1_804, skk_1022, skk_1026, \
                         skk_1028 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1278[k] = f_10 * ski0_802[k]
                    - f_11 * ski1_802[k]
                    + f_3 * pc_x[k] * skk_1026[k];

        t_1279[k] = f_0 * sik_770[k]
                    + f_3 * pc_y[k] * skk_1022[k];

        t_1280[k] = f_10 * ski0_804[k]
                    - f_11 * ski1_804[k]
                    + f_3 * pc_x[k] * skk_1028[k];
    }

#pragma omp simd aligned(t_1281, t_1282, t_1283, t_1284, pc_x, pc_z, ski0_805, ski0_807, \
                         ski0_808, ski1_805, ski1_807, ski1_808, skk_1023, skk_1029, skk_1031, \
                         skk_1032 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1281[k] = f_12 * ski0_805[k]
                    - f_13 * ski1_805[k]
                    + f_3 * pc_x[k] * skk_1029[k];

        t_1282[k] = f_3 * pc_z[k] * skk_1023[k];

        t_1283[k] = f_12 * ski0_807[k]
                    - f_13 * ski1_807[k]
                    + f_3 * pc_x[k] * skk_1031[k];

        t_1284[k] = f_12 * ski0_808[k]
                    - f_13 * ski1_808[k]
                    + f_3 * pc_x[k] * skk_1032[k];
    }

#pragma omp simd aligned(t_1285, t_1286, t_1287, t_1288, pc_x, pc_y, sik_776, ski0_809, \
                         ski0_811, ski1_809, ski1_811, skk_1028, skk_1033, skk_1035, \
                         skk_1036 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1285[k] = f_12 * ski0_809[k]
                    - f_13 * ski1_809[k]
                    + f_3 * pc_x[k] * skk_1033[k];

        t_1286[k] = f_0 * sik_776[k]
                    + f_3 * pc_y[k] * skk_1028[k];

        t_1287[k] = f_12 * ski0_811[k]
                    - f_13 * ski1_811[k]
                    + f_3 * pc_x[k] * skk_1035[k];

        t_1288[k] = f_3 * pc_x[k] * skk_1036[k];
    }

#pragma omp simd aligned(t_1289, t_1290, t_1291, t_1292, t_1293, t_1294, t_1295, pc_x, \
                         skk_1037, skk_1038, skk_1039, skk_1040, skk_1041, skk_1042, \
                         skk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1289[k] = f_3 * pc_x[k] * skk_1037[k];

        t_1290[k] = f_3 * pc_x[k] * skk_1038[k];

        t_1291[k] = f_3 * pc_x[k] * skk_1039[k];

        t_1292[k] = f_3 * pc_x[k] * skk_1040[k];

        t_1293[k] = f_3 * pc_x[k] * skk_1041[k];

        t_1294[k] = f_3 * pc_x[k] * skk_1042[k];

        t_1295[k] = f_3 * pc_x[k] * skk_1043[k];
    }

#pragma omp simd aligned(t_1296, t_1297, t_1298, pc_y, pc_z, sik_784, sik_786, ski0_805, \
                         ski0_807, ski1_805, ski1_807, skk_1036, \
                         skk_1038 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1296[k] = f_0 * sik_784[k]
                    + f_1 * ski0_805[k]
                    - f_2 * ski1_805[k]
                    + f_3 * pc_y[k] * skk_1036[k];

        t_1297[k] = f_3 * pc_z[k] * skk_1036[k];

        t_1298[k] = f_0 * sik_786[k]
                    + f_4 * ski0_807[k]
                    - f_5 * ski1_807[k]
                    + f_3 * pc_y[k] * skk_1038[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, pc_y, sik_787, sik_788, sik_789, ski0_808, \
                         ski0_809, ski0_810, ski1_808, ski1_809, ski1_810, skk_1039, skk_1040, \
                         skk_1041 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = f_0 * sik_787[k]
                    + f_6 * ski0_808[k]
                    - f_7 * ski1_808[k]
                    + f_3 * pc_y[k] * skk_1039[k];

        t_1300[k] = f_0 * sik_788[k]
                    + f_8 * ski0_809[k]
                    - f_9 * ski1_809[k]
                    + f_3 * pc_y[k] * skk_1040[k];

        t_1301[k] = f_0 * sik_789[k]
                    + f_10 * ski0_810[k]
                    - f_11 * ski1_810[k]
                    + f_3 * pc_y[k] * skk_1041[k];
    }

#pragma omp simd aligned(t_1302, t_1303, t_1304, t_1305, pb_z, pc_y, pc_z, sil0_945, sik_790, \
                         sik_791, sil1_945, ski0_811, ski1_811, skk_1042, \
                         skk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1302[k] = f_0 * sik_790[k]
                    + f_12 * ski0_811[k]
                    - f_13 * ski1_811[k]
                    + f_3 * pc_y[k] * skk_1042[k];

        t_1303[k] = f_0 * sik_791[k]
                    + f_3 * pc_y[k] * skk_1043[k];

        t_1304[k] = f_1 * ski0_811[k]
                    - f_2 * ski1_811[k]
                    + f_3 * pc_z[k] * skk_1043[k];

        t_1305[k] = pb_z[k] * sil0_945[k]
                    - f_14 * pc_z[k] * sil1_945[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pb_z, pc_y, pc_z, sil0_948, sik_756, \
                         sik_792, sik_794, sil1_948, skk_1044, \
                         skk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_20 * sik_792[k]
                    + f_3 * pc_y[k] * skk_1044[k];

        t_1307[k] = f_15 * sik_756[k]
                    + f_3 * pc_z[k] * skk_1044[k];

        t_1308[k] = pb_z[k] * sil0_948[k]
                    - f_14 * pc_z[k] * sil1_948[k];

        t_1309[k] = f_20 * sik_794[k]
                    + f_3 * pc_y[k] * skk_1046[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, t_1313, pb_z, pc_x, pc_y, pc_z, sil0_951, \
                         sik_759, sik_797, sil1_951, ski0_817, ski1_817, skk_1047, \
                         skk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = f_4 * ski0_817[k]
                    - f_5 * ski1_817[k]
                    + f_3 * pc_x[k] * skk_1049[k];

        t_1311[k] = pb_z[k] * sil0_951[k]
                    - f_14 * pc_z[k] * sil1_951[k];

        t_1312[k] = f_15 * sik_759[k]
                    + f_3 * pc_z[k] * skk_1047[k];

        t_1313[k] = f_20 * sik_797[k]
                    + f_3 * pc_y[k] * skk_1049[k];
    }

#pragma omp simd aligned(t_1314, t_1315, t_1316, pb_z, pc_x, pc_z, sil0_955, sik_762, \
                         sil1_955, ski0_821, ski1_821, skk_1050, \
                         skk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1314[k] = f_6 * ski0_821[k]
                    - f_7 * ski1_821[k]
                    + f_3 * pc_x[k] * skk_1053[k];

        t_1315[k] = pb_z[k] * sil0_955[k]
                    - f_14 * pc_z[k] * sil1_955[k];

        t_1316[k] = f_15 * sik_762[k]
                    + f_3 * pc_z[k] * skk_1050[k];
    }

#pragma omp simd aligned(t_1317, t_1318, t_1319, pc_x, pc_y, sik_801, ski0_824, ski0_826, \
                         ski1_824, ski1_826, skk_1053, skk_1056, \
                         skk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1317[k] = f_8 * ski0_824[k]
                    - f_9 * ski1_824[k]
                    + f_3 * pc_x[k] * skk_1056[k];

        t_1318[k] = f_20 * sik_801[k]
                    + f_3 * pc_y[k] * skk_1053[k];

        t_1319[k] = f_8 * ski0_826[k]
                    - f_9 * ski1_826[k]
                    + f_3 * pc_x[k] * skk_1058[k];
    }

#pragma omp simd aligned(t_1320, t_1321, t_1322, pb_z, pc_x, pc_z, sil0_960, sik_766, \
                         sil1_960, ski0_829, ski1_829, skk_1054, \
                         skk_1061 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1320[k] = pb_z[k] * sil0_960[k]
                    - f_14 * pc_z[k] * sil1_960[k];

        t_1321[k] = f_15 * sik_766[k]
                    + f_3 * pc_z[k] * skk_1054[k];

        t_1322[k] = f_10 * ski0_829[k]
                    - f_11 * ski1_829[k]
                    + f_3 * pc_x[k] * skk_1061[k];
    }

#pragma omp simd aligned(t_1323, t_1324, t_1325, pc_x, pc_y, sik_806, ski0_830, ski0_832, \
                         ski1_830, ski1_832, skk_1058, skk_1062, \
                         skk_1064 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1323[k] = f_10 * ski0_830[k]
                    - f_11 * ski1_830[k]
                    + f_3 * pc_x[k] * skk_1062[k];

        t_1324[k] = f_20 * sik_806[k]
                    + f_3 * pc_y[k] * skk_1058[k];

        t_1325[k] = f_10 * ski0_832[k]
                    - f_11 * ski1_832[k]
                    + f_3 * pc_x[k] * skk_1064[k];
    }

#pragma omp simd aligned(t_1326, t_1327, t_1328, pb_z, pc_x, pc_z, sil0_966, sik_771, \
                         sil1_966, ski0_835, ski1_835, skk_1059, \
                         skk_1067 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1326[k] = pb_z[k] * sil0_966[k]
                    - f_14 * pc_z[k] * sil1_966[k];

        t_1327[k] = f_15 * sik_771[k]
                    + f_3 * pc_z[k] * skk_1059[k];

        t_1328[k] = f_12 * ski0_835[k]
                    - f_13 * ski1_835[k]
                    + f_3 * pc_x[k] * skk_1067[k];
    }

#pragma omp simd aligned(t_1329, t_1330, t_1331, pc_x, pc_y, sik_812, ski0_836, ski0_837, \
                         ski1_836, ski1_837, skk_1064, skk_1068, \
                         skk_1069 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1329[k] = f_12 * ski0_836[k]
                    - f_13 * ski1_836[k]
                    + f_3 * pc_x[k] * skk_1068[k];

        t_1330[k] = f_12 * ski0_837[k]
                    - f_13 * ski1_837[k]
                    + f_3 * pc_x[k] * skk_1069[k];

        t_1331[k] = f_20 * sik_812[k]
                    + f_3 * pc_y[k] * skk_1064[k];
    }

#pragma omp simd aligned(t_1332, t_1333, t_1334, t_1335, t_1336, t_1337, pc_x, ski0_839, \
                         ski1_839, skk_1071, skk_1072, skk_1073, skk_1074, skk_1075, \
                         skk_1076 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1332[k] = f_12 * ski0_839[k]
                    - f_13 * ski1_839[k]
                    + f_3 * pc_x[k] * skk_1071[k];

        t_1333[k] = f_3 * pc_x[k] * skk_1072[k];

        t_1334[k] = f_3 * pc_x[k] * skk_1073[k];

        t_1335[k] = f_3 * pc_x[k] * skk_1074[k];

        t_1336[k] = f_3 * pc_x[k] * skk_1075[k];

        t_1337[k] = f_3 * pc_x[k] * skk_1076[k];
    }

#pragma omp simd aligned(t_1338, t_1339, t_1340, t_1341, t_1342, pb_z, pc_x, pc_z, sil0_981, \
                         sik_784, sil1_981, skk_1072, skk_1077, skk_1078, \
                         skk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1338[k] = f_3 * pc_x[k] * skk_1077[k];

        t_1339[k] = f_3 * pc_x[k] * skk_1078[k];

        t_1340[k] = f_3 * pc_x[k] * skk_1079[k];

        t_1341[k] = pb_z[k] * sil0_981[k]
                    - f_14 * pc_z[k] * sil1_981[k];

        t_1342[k] = f_15 * sik_784[k]
                    + f_3 * pc_z[k] * skk_1072[k];
    }

#pragma omp simd aligned(t_1343, t_1344, t_1345, pb_z, pc_z, sil0_983, sil0_984, sil0_985, \
                         sik_785, sik_786, sik_787, sil1_983, sil1_984, \
                         sil1_985 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1343[k] = pb_z[k] * sil0_983[k]
                    + f_16 * sik_785[k]
                    - f_14 * pc_z[k] * sil1_983[k];

        t_1344[k] = pb_z[k] * sil0_984[k]
                    + f_17 * sik_786[k]
                    - f_14 * pc_z[k] * sil1_984[k];

        t_1345[k] = pb_z[k] * sil0_985[k]
                    + f_18 * sik_787[k]
                    - f_14 * pc_z[k] * sil1_985[k];
    }

#pragma omp simd aligned(t_1346, t_1347, t_1348, pb_z, pc_y, pc_z, sil0_986, sil0_987, \
                         sik_788, sik_789, sik_827, sil1_986, sil1_987, \
                         skk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1346[k] = pb_z[k] * sil0_986[k]
                    + f_19 * sik_788[k]
                    - f_14 * pc_z[k] * sil1_986[k];

        t_1347[k] = pb_z[k] * sil0_987[k]
                    + f_20 * sik_789[k]
                    - f_14 * pc_z[k] * sil1_987[k];

        t_1348[k] = f_20 * sik_827[k]
                    + f_3 * pc_y[k] * skk_1079[k];
    }

#pragma omp simd aligned(t_1349, t_1350, t_1351, t_1352, pc_x, pc_y, pc_z, sik_791, sik_792, \
                         sik_828, ski0_839, ski0_840, ski1_839, ski1_840, skk_1079, \
                         skk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1349[k] = f_15 * sik_791[k]
                    + f_1 * ski0_839[k]
                    - f_2 * ski1_839[k]
                    + f_3 * pc_z[k] * skk_1079[k];

        t_1350[k] = f_1 * ski0_840[k]
                    - f_2 * ski1_840[k]
                    + f_3 * pc_x[k] * skk_1080[k];

        t_1351[k] = f_19 * sik_828[k]
                    + f_3 * pc_y[k] * skk_1080[k];

        t_1352[k] = f_16 * sik_792[k]
                    + f_3 * pc_z[k] * skk_1080[k];
    }

#pragma omp simd aligned(t_1353, t_1354, t_1355, pc_x, pc_y, sik_830, ski0_843, ski0_845, \
                         ski1_843, ski1_845, skk_1082, skk_1083, \
                         skk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1353[k] = f_4 * ski0_843[k]
                    - f_5 * ski1_843[k]
                    + f_3 * pc_x[k] * skk_1083[k];

        t_1354[k] = f_19 * sik_830[k]
                    + f_3 * pc_y[k] * skk_1082[k];

        t_1355[k] = f_4 * ski0_845[k]
                    - f_5 * ski1_845[k]
                    + f_3 * pc_x[k] * skk_1085[k];
    }

#pragma omp simd aligned(t_1356, t_1357, t_1358, pc_x, pc_y, pc_z, sik_795, sik_833, ski0_846, \
                         ski1_846, skk_1083, skk_1085, skk_1086 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1356[k] = f_6 * ski0_846[k]
                    - f_7 * ski1_846[k]
                    + f_3 * pc_x[k] * skk_1086[k];

        t_1357[k] = f_16 * sik_795[k]
                    + f_3 * pc_z[k] * skk_1083[k];

        t_1358[k] = f_19 * sik_833[k]
                    + f_3 * pc_y[k] * skk_1085[k];
    }

#pragma omp simd aligned(t_1359, t_1360, t_1361, pc_x, pc_z, sik_798, ski0_849, ski0_850, \
                         ski1_849, ski1_850, skk_1086, skk_1089, \
                         skk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1359[k] = f_6 * ski0_849[k]
                    - f_7 * ski1_849[k]
                    + f_3 * pc_x[k] * skk_1089[k];

        t_1360[k] = f_8 * ski0_850[k]
                    - f_9 * ski1_850[k]
                    + f_3 * pc_x[k] * skk_1090[k];

        t_1361[k] = f_16 * sik_798[k]
                    + f_3 * pc_z[k] * skk_1086[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t sik, const size_t ski0,
                                                           const size_t ski1, const size_t skk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;

    auto *t_1362 = buffer.data(target + 1362);
    auto *t_1363 = buffer.data(target + 1363);
    auto *t_1364 = buffer.data(target + 1364);
    auto *t_1365 = buffer.data(target + 1365);
    auto *t_1366 = buffer.data(target + 1366);
    auto *t_1367 = buffer.data(target + 1367);
    auto *t_1368 = buffer.data(target + 1368);
    auto *t_1369 = buffer.data(target + 1369);
    auto *t_1370 = buffer.data(target + 1370);
    auto *t_1371 = buffer.data(target + 1371);
    auto *t_1372 = buffer.data(target + 1372);
    auto *t_1373 = buffer.data(target + 1373);
    auto *t_1374 = buffer.data(target + 1374);
    auto *t_1375 = buffer.data(target + 1375);
    auto *t_1376 = buffer.data(target + 1376);
    auto *t_1377 = buffer.data(target + 1377);
    auto *t_1378 = buffer.data(target + 1378);
    auto *t_1379 = buffer.data(target + 1379);
    auto *t_1380 = buffer.data(target + 1380);
    auto *t_1381 = buffer.data(target + 1381);
    auto *t_1382 = buffer.data(target + 1382);
    auto *t_1383 = buffer.data(target + 1383);
    auto *t_1384 = buffer.data(target + 1384);
    auto *t_1385 = buffer.data(target + 1385);
    auto *t_1386 = buffer.data(target + 1386);
    auto *t_1387 = buffer.data(target + 1387);
    auto *t_1388 = buffer.data(target + 1388);
    auto *t_1389 = buffer.data(target + 1389);
    auto *t_1390 = buffer.data(target + 1390);
    auto *t_1391 = buffer.data(target + 1391);
    auto *t_1392 = buffer.data(target + 1392);
    auto *t_1393 = buffer.data(target + 1393);
    auto *t_1394 = buffer.data(target + 1394);
    auto *t_1395 = buffer.data(target + 1395);
    auto *t_1396 = buffer.data(target + 1396);
    auto *t_1397 = buffer.data(target + 1397);
    auto *t_1398 = buffer.data(target + 1398);
    auto *t_1399 = buffer.data(target + 1399);
    auto *t_1400 = buffer.data(target + 1400);
    auto *t_1401 = buffer.data(target + 1401);
    auto *t_1402 = buffer.data(target + 1402);
    auto *t_1403 = buffer.data(target + 1403);
    auto *t_1404 = buffer.data(target + 1404);
    auto *t_1405 = buffer.data(target + 1405);
    auto *t_1406 = buffer.data(target + 1406);
    auto *t_1407 = buffer.data(target + 1407);
    auto *t_1408 = buffer.data(target + 1408);
    auto *t_1409 = buffer.data(target + 1409);
    auto *t_1410 = buffer.data(target + 1410);
    auto *t_1411 = buffer.data(target + 1411);
    auto *t_1412 = buffer.data(target + 1412);
    auto *t_1413 = buffer.data(target + 1413);
    auto *t_1414 = buffer.data(target + 1414);
    auto *t_1415 = buffer.data(target + 1415);
    auto *t_1416 = buffer.data(target + 1416);
    auto *t_1417 = buffer.data(target + 1417);
    auto *t_1418 = buffer.data(target + 1418);
    auto *t_1419 = buffer.data(target + 1419);
    auto *t_1420 = buffer.data(target + 1420);
    auto *t_1421 = buffer.data(target + 1421);
    auto *t_1422 = buffer.data(target + 1422);
    auto *t_1423 = buffer.data(target + 1423);
    auto *t_1424 = buffer.data(target + 1424);
    auto *t_1425 = buffer.data(target + 1425);
    auto *t_1426 = buffer.data(target + 1426);
    auto *t_1427 = buffer.data(target + 1427);
    auto *t_1428 = buffer.data(target + 1428);
    auto *t_1429 = buffer.data(target + 1429);
    auto *t_1430 = buffer.data(target + 1430);
    auto *t_1431 = buffer.data(target + 1431);
    auto *t_1432 = buffer.data(target + 1432);
    auto *t_1433 = buffer.data(target + 1433);
    auto *t_1434 = buffer.data(target + 1434);
    auto *t_1435 = buffer.data(target + 1435);
    auto *t_1436 = buffer.data(target + 1436);
    auto *t_1437 = buffer.data(target + 1437);
    auto *t_1438 = buffer.data(target + 1438);
    auto *t_1439 = buffer.data(target + 1439);
    auto *t_1440 = buffer.data(target + 1440);
    auto *t_1441 = buffer.data(target + 1441);
    auto *t_1442 = buffer.data(target + 1442);
    auto *t_1443 = buffer.data(target + 1443);
    auto *t_1444 = buffer.data(target + 1444);
    auto *t_1445 = buffer.data(target + 1445);
    auto *t_1446 = buffer.data(target + 1446);
    auto *t_1447 = buffer.data(target + 1447);
    auto *t_1448 = buffer.data(target + 1448);
    auto *t_1449 = buffer.data(target + 1449);
    auto *t_1450 = buffer.data(target + 1450);
    auto *t_1451 = buffer.data(target + 1451);
    auto *t_1452 = buffer.data(target + 1452);
    auto *t_1453 = buffer.data(target + 1453);
    auto *t_1454 = buffer.data(target + 1454);
    auto *t_1455 = buffer.data(target + 1455);
    auto *t_1456 = buffer.data(target + 1456);
    auto *t_1457 = buffer.data(target + 1457);
    auto *t_1458 = buffer.data(target + 1458);
    auto *t_1459 = buffer.data(target + 1459);
    auto *t_1460 = buffer.data(target + 1460);
    auto *t_1461 = buffer.data(target + 1461);
    auto *t_1462 = buffer.data(target + 1462);
    auto *t_1463 = buffer.data(target + 1463);
    auto *t_1464 = buffer.data(target + 1464);
    auto *t_1465 = buffer.data(target + 1465);
    auto *t_1466 = buffer.data(target + 1466);
    auto *t_1467 = buffer.data(target + 1467);
    auto *t_1468 = buffer.data(target + 1468);
    auto *t_1469 = buffer.data(target + 1469);
    auto *t_1470 = buffer.data(target + 1470);
    auto *t_1471 = buffer.data(target + 1471);
    auto *t_1472 = buffer.data(target + 1472);
    auto *t_1473 = buffer.data(target + 1473);
    auto *t_1474 = buffer.data(target + 1474);
    auto *t_1475 = buffer.data(target + 1475);
    auto *t_1476 = buffer.data(target + 1476);
    auto *t_1477 = buffer.data(target + 1477);
    auto *t_1478 = buffer.data(target + 1478);
    auto *t_1479 = buffer.data(target + 1479);
    auto *t_1480 = buffer.data(target + 1480);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik_802 = buffer.data(sik + 802);
    const auto *sik_807 = buffer.data(sik + 807);
    const auto *sik_820 = buffer.data(sik + 820);
    const auto *sik_827 = buffer.data(sik + 827);
    const auto *sik_828 = buffer.data(sik + 828);
    const auto *sik_831 = buffer.data(sik + 831);
    const auto *sik_834 = buffer.data(sik + 834);
    const auto *sik_837 = buffer.data(sik + 837);
    const auto *sik_838 = buffer.data(sik + 838);
    const auto *sik_842 = buffer.data(sik + 842);
    const auto *sik_843 = buffer.data(sik + 843);
    const auto *sik_848 = buffer.data(sik + 848);
    const auto *sik_856 = buffer.data(sik + 856);
    const auto *sik_858 = buffer.data(sik + 858);
    const auto *sik_859 = buffer.data(sik + 859);
    const auto *sik_860 = buffer.data(sik + 860);
    const auto *sik_861 = buffer.data(sik + 861);
    const auto *sik_862 = buffer.data(sik + 862);
    const auto *sik_863 = buffer.data(sik + 863);
    const auto *sik_864 = buffer.data(sik + 864);
    const auto *sik_866 = buffer.data(sik + 866);
    const auto *sik_867 = buffer.data(sik + 867);
    const auto *sik_869 = buffer.data(sik + 869);
    const auto *sik_870 = buffer.data(sik + 870);
    const auto *sik_873 = buffer.data(sik + 873);
    const auto *sik_874 = buffer.data(sik + 874);
    const auto *sik_878 = buffer.data(sik + 878);
    const auto *sik_879 = buffer.data(sik + 879);
    const auto *sik_884 = buffer.data(sik + 884);
    const auto *sik_892 = buffer.data(sik + 892);
    const auto *sik_894 = buffer.data(sik + 894);
    const auto *sik_895 = buffer.data(sik + 895);
    const auto *sik_896 = buffer.data(sik + 896);
    const auto *sik_897 = buffer.data(sik + 897);
    const auto *sik_898 = buffer.data(sik + 898);
    const auto *sik_899 = buffer.data(sik + 899);
    const auto *sik_900 = buffer.data(sik + 900);
    const auto *sik_902 = buffer.data(sik + 902);
    const auto *sik_905 = buffer.data(sik + 905);
    const auto *sik_909 = buffer.data(sik + 909);
    const auto *sik_914 = buffer.data(sik + 914);
    const auto *sik_920 = buffer.data(sik + 920);
    const auto *sik_928 = buffer.data(sik + 928);
    const auto *sik_930 = buffer.data(sik + 930);
    const auto *sik_931 = buffer.data(sik + 931);
    const auto *sik_932 = buffer.data(sik + 932);

    const auto *ski0_852 = buffer.data(ski0 + 852);
    const auto *ski0_854 = buffer.data(ski0 + 854);
    const auto *ski0_855 = buffer.data(ski0 + 855);
    const auto *ski0_857 = buffer.data(ski0 + 857);
    const auto *ski0_858 = buffer.data(ski0 + 858);
    const auto *ski0_860 = buffer.data(ski0 + 860);
    const auto *ski0_861 = buffer.data(ski0 + 861);
    const auto *ski0_863 = buffer.data(ski0 + 863);
    const auto *ski0_864 = buffer.data(ski0 + 864);
    const auto *ski0_865 = buffer.data(ski0 + 865);
    const auto *ski0_866 = buffer.data(ski0 + 866);
    const auto *ski0_867 = buffer.data(ski0 + 867);
    const auto *ski0_868 = buffer.data(ski0 + 868);
    const auto *ski0_871 = buffer.data(ski0 + 871);
    const auto *ski0_873 = buffer.data(ski0 + 873);
    const auto *ski0_874 = buffer.data(ski0 + 874);
    const auto *ski0_877 = buffer.data(ski0 + 877);
    const auto *ski0_878 = buffer.data(ski0 + 878);
    const auto *ski0_880 = buffer.data(ski0 + 880);
    const auto *ski0_882 = buffer.data(ski0 + 882);
    const auto *ski0_883 = buffer.data(ski0 + 883);
    const auto *ski0_885 = buffer.data(ski0 + 885);
    const auto *ski0_886 = buffer.data(ski0 + 886);
    const auto *ski0_888 = buffer.data(ski0 + 888);
    const auto *ski0_889 = buffer.data(ski0 + 889);
    const auto *ski0_891 = buffer.data(ski0 + 891);
    const auto *ski0_892 = buffer.data(ski0 + 892);
    const auto *ski0_893 = buffer.data(ski0 + 893);
    const auto *ski0_894 = buffer.data(ski0 + 894);
    const auto *ski0_895 = buffer.data(ski0 + 895);
    const auto *ski0_896 = buffer.data(ski0 + 896);
    const auto *ski0_899 = buffer.data(ski0 + 899);
    const auto *ski0_901 = buffer.data(ski0 + 901);
    const auto *ski0_902 = buffer.data(ski0 + 902);
    const auto *ski0_905 = buffer.data(ski0 + 905);
    const auto *ski0_906 = buffer.data(ski0 + 906);
    const auto *ski0_908 = buffer.data(ski0 + 908);
    const auto *ski0_910 = buffer.data(ski0 + 910);
    const auto *ski0_911 = buffer.data(ski0 + 911);
    const auto *ski0_913 = buffer.data(ski0 + 913);
    const auto *ski0_914 = buffer.data(ski0 + 914);
    const auto *ski0_916 = buffer.data(ski0 + 916);
    const auto *ski0_917 = buffer.data(ski0 + 917);
    const auto *ski0_919 = buffer.data(ski0 + 919);
    const auto *ski0_920 = buffer.data(ski0 + 920);
    const auto *ski0_921 = buffer.data(ski0 + 921);
    const auto *ski0_923 = buffer.data(ski0 + 923);

    const auto *ski1_852 = buffer.data(ski1 + 852);
    const auto *ski1_854 = buffer.data(ski1 + 854);
    const auto *ski1_855 = buffer.data(ski1 + 855);
    const auto *ski1_857 = buffer.data(ski1 + 857);
    const auto *ski1_858 = buffer.data(ski1 + 858);
    const auto *ski1_860 = buffer.data(ski1 + 860);
    const auto *ski1_861 = buffer.data(ski1 + 861);
    const auto *ski1_863 = buffer.data(ski1 + 863);
    const auto *ski1_864 = buffer.data(ski1 + 864);
    const auto *ski1_865 = buffer.data(ski1 + 865);
    const auto *ski1_866 = buffer.data(ski1 + 866);
    const auto *ski1_867 = buffer.data(ski1 + 867);
    const auto *ski1_868 = buffer.data(ski1 + 868);
    const auto *ski1_871 = buffer.data(ski1 + 871);
    const auto *ski1_873 = buffer.data(ski1 + 873);
    const auto *ski1_874 = buffer.data(ski1 + 874);
    const auto *ski1_877 = buffer.data(ski1 + 877);
    const auto *ski1_878 = buffer.data(ski1 + 878);
    const auto *ski1_880 = buffer.data(ski1 + 880);
    const auto *ski1_882 = buffer.data(ski1 + 882);
    const auto *ski1_883 = buffer.data(ski1 + 883);
    const auto *ski1_885 = buffer.data(ski1 + 885);
    const auto *ski1_886 = buffer.data(ski1 + 886);
    const auto *ski1_888 = buffer.data(ski1 + 888);
    const auto *ski1_889 = buffer.data(ski1 + 889);
    const auto *ski1_891 = buffer.data(ski1 + 891);
    const auto *ski1_892 = buffer.data(ski1 + 892);
    const auto *ski1_893 = buffer.data(ski1 + 893);
    const auto *ski1_894 = buffer.data(ski1 + 894);
    const auto *ski1_895 = buffer.data(ski1 + 895);
    const auto *ski1_896 = buffer.data(ski1 + 896);
    const auto *ski1_899 = buffer.data(ski1 + 899);
    const auto *ski1_901 = buffer.data(ski1 + 901);
    const auto *ski1_902 = buffer.data(ski1 + 902);
    const auto *ski1_905 = buffer.data(ski1 + 905);
    const auto *ski1_906 = buffer.data(ski1 + 906);
    const auto *ski1_908 = buffer.data(ski1 + 908);
    const auto *ski1_910 = buffer.data(ski1 + 910);
    const auto *ski1_911 = buffer.data(ski1 + 911);
    const auto *ski1_913 = buffer.data(ski1 + 913);
    const auto *ski1_914 = buffer.data(ski1 + 914);
    const auto *ski1_916 = buffer.data(ski1 + 916);
    const auto *ski1_917 = buffer.data(ski1 + 917);
    const auto *ski1_919 = buffer.data(ski1 + 919);
    const auto *ski1_920 = buffer.data(ski1 + 920);
    const auto *ski1_921 = buffer.data(ski1 + 921);
    const auto *ski1_923 = buffer.data(ski1 + 923);

    const auto *skk_1089 = buffer.data(skk + 1089);
    const auto *skk_1090 = buffer.data(skk + 1090);
    const auto *skk_1092 = buffer.data(skk + 1092);
    const auto *skk_1094 = buffer.data(skk + 1094);
    const auto *skk_1095 = buffer.data(skk + 1095);
    const auto *skk_1097 = buffer.data(skk + 1097);
    const auto *skk_1098 = buffer.data(skk + 1098);
    const auto *skk_1100 = buffer.data(skk + 1100);
    const auto *skk_1101 = buffer.data(skk + 1101);
    const auto *skk_1103 = buffer.data(skk + 1103);
    const auto *skk_1104 = buffer.data(skk + 1104);
    const auto *skk_1105 = buffer.data(skk + 1105);
    const auto *skk_1107 = buffer.data(skk + 1107);
    const auto *skk_1108 = buffer.data(skk + 1108);
    const auto *skk_1109 = buffer.data(skk + 1109);
    const auto *skk_1110 = buffer.data(skk + 1110);
    const auto *skk_1111 = buffer.data(skk + 1111);
    const auto *skk_1112 = buffer.data(skk + 1112);
    const auto *skk_1113 = buffer.data(skk + 1113);
    const auto *skk_1114 = buffer.data(skk + 1114);
    const auto *skk_1115 = buffer.data(skk + 1115);
    const auto *skk_1116 = buffer.data(skk + 1116);
    const auto *skk_1118 = buffer.data(skk + 1118);
    const auto *skk_1119 = buffer.data(skk + 1119);
    const auto *skk_1121 = buffer.data(skk + 1121);
    const auto *skk_1122 = buffer.data(skk + 1122);
    const auto *skk_1125 = buffer.data(skk + 1125);
    const auto *skk_1126 = buffer.data(skk + 1126);
    const auto *skk_1128 = buffer.data(skk + 1128);
    const auto *skk_1130 = buffer.data(skk + 1130);
    const auto *skk_1131 = buffer.data(skk + 1131);
    const auto *skk_1133 = buffer.data(skk + 1133);
    const auto *skk_1134 = buffer.data(skk + 1134);
    const auto *skk_1136 = buffer.data(skk + 1136);
    const auto *skk_1137 = buffer.data(skk + 1137);
    const auto *skk_1139 = buffer.data(skk + 1139);
    const auto *skk_1140 = buffer.data(skk + 1140);
    const auto *skk_1141 = buffer.data(skk + 1141);
    const auto *skk_1143 = buffer.data(skk + 1143);
    const auto *skk_1144 = buffer.data(skk + 1144);
    const auto *skk_1145 = buffer.data(skk + 1145);
    const auto *skk_1146 = buffer.data(skk + 1146);
    const auto *skk_1147 = buffer.data(skk + 1147);
    const auto *skk_1148 = buffer.data(skk + 1148);
    const auto *skk_1149 = buffer.data(skk + 1149);
    const auto *skk_1150 = buffer.data(skk + 1150);
    const auto *skk_1151 = buffer.data(skk + 1151);
    const auto *skk_1152 = buffer.data(skk + 1152);
    const auto *skk_1154 = buffer.data(skk + 1154);
    const auto *skk_1155 = buffer.data(skk + 1155);
    const auto *skk_1157 = buffer.data(skk + 1157);
    const auto *skk_1158 = buffer.data(skk + 1158);
    const auto *skk_1161 = buffer.data(skk + 1161);
    const auto *skk_1162 = buffer.data(skk + 1162);
    const auto *skk_1164 = buffer.data(skk + 1164);
    const auto *skk_1166 = buffer.data(skk + 1166);
    const auto *skk_1167 = buffer.data(skk + 1167);
    const auto *skk_1169 = buffer.data(skk + 1169);
    const auto *skk_1170 = buffer.data(skk + 1170);
    const auto *skk_1172 = buffer.data(skk + 1172);
    const auto *skk_1173 = buffer.data(skk + 1173);
    const auto *skk_1175 = buffer.data(skk + 1175);
    const auto *skk_1176 = buffer.data(skk + 1176);
    const auto *skk_1177 = buffer.data(skk + 1177);
    const auto *skk_1179 = buffer.data(skk + 1179);
    const auto *skk_1180 = buffer.data(skk + 1180);
    const auto *skk_1181 = buffer.data(skk + 1181);
    const auto *skk_1182 = buffer.data(skk + 1182);
    const auto *skk_1183 = buffer.data(skk + 1183);
    const auto *skk_1184 = buffer.data(skk + 1184);
    const auto *skk_1185 = buffer.data(skk + 1185);
    const auto *skk_1186 = buffer.data(skk + 1186);
    const auto *skk_1187 = buffer.data(skk + 1187);

#pragma omp simd aligned(t_1362, t_1363, t_1364, pc_x, pc_y, sik_837, ski0_852, ski0_854, \
                         ski1_852, ski1_854, skk_1089, skk_1092, \
                         skk_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1362[k] = f_8 * ski0_852[k]
                    - f_9 * ski1_852[k]
                    + f_3 * pc_x[k] * skk_1092[k];

        t_1363[k] = f_19 * sik_837[k]
                    + f_3 * pc_y[k] * skk_1089[k];

        t_1364[k] = f_8 * ski0_854[k]
                    - f_9 * ski1_854[k]
                    + f_3 * pc_x[k] * skk_1094[k];
    }

#pragma omp simd aligned(t_1365, t_1366, t_1367, pc_x, pc_z, sik_802, ski0_855, ski0_857, \
                         ski1_855, ski1_857, skk_1090, skk_1095, \
                         skk_1097 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1365[k] = f_10 * ski0_855[k]
                    - f_11 * ski1_855[k]
                    + f_3 * pc_x[k] * skk_1095[k];

        t_1366[k] = f_16 * sik_802[k]
                    + f_3 * pc_z[k] * skk_1090[k];

        t_1367[k] = f_10 * ski0_857[k]
                    - f_11 * ski1_857[k]
                    + f_3 * pc_x[k] * skk_1097[k];
    }

#pragma omp simd aligned(t_1368, t_1369, t_1370, pc_x, pc_y, sik_842, ski0_858, ski0_860, \
                         ski1_858, ski1_860, skk_1094, skk_1098, \
                         skk_1100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1368[k] = f_10 * ski0_858[k]
                    - f_11 * ski1_858[k]
                    + f_3 * pc_x[k] * skk_1098[k];

        t_1369[k] = f_19 * sik_842[k]
                    + f_3 * pc_y[k] * skk_1094[k];

        t_1370[k] = f_10 * ski0_860[k]
                    - f_11 * ski1_860[k]
                    + f_3 * pc_x[k] * skk_1100[k];
    }

#pragma omp simd aligned(t_1371, t_1372, t_1373, pc_x, pc_z, sik_807, ski0_861, ski0_863, \
                         ski1_861, ski1_863, skk_1095, skk_1101, \
                         skk_1103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1371[k] = f_12 * ski0_861[k]
                    - f_13 * ski1_861[k]
                    + f_3 * pc_x[k] * skk_1101[k];

        t_1372[k] = f_16 * sik_807[k]
                    + f_3 * pc_z[k] * skk_1095[k];

        t_1373[k] = f_12 * ski0_863[k]
                    - f_13 * ski1_863[k]
                    + f_3 * pc_x[k] * skk_1103[k];
    }

#pragma omp simd aligned(t_1374, t_1375, t_1376, pc_x, pc_y, sik_848, ski0_864, ski0_865, \
                         ski1_864, ski1_865, skk_1100, skk_1104, \
                         skk_1105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1374[k] = f_12 * ski0_864[k]
                    - f_13 * ski1_864[k]
                    + f_3 * pc_x[k] * skk_1104[k];

        t_1375[k] = f_12 * ski0_865[k]
                    - f_13 * ski1_865[k]
                    + f_3 * pc_x[k] * skk_1105[k];

        t_1376[k] = f_19 * sik_848[k]
                    + f_3 * pc_y[k] * skk_1100[k];
    }

#pragma omp simd aligned(t_1377, t_1378, t_1379, t_1380, t_1381, t_1382, pc_x, ski0_867, \
                         ski1_867, skk_1107, skk_1108, skk_1109, skk_1110, skk_1111, \
                         skk_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1377[k] = f_12 * ski0_867[k]
                    - f_13 * ski1_867[k]
                    + f_3 * pc_x[k] * skk_1107[k];

        t_1378[k] = f_3 * pc_x[k] * skk_1108[k];

        t_1379[k] = f_3 * pc_x[k] * skk_1109[k];

        t_1380[k] = f_3 * pc_x[k] * skk_1110[k];

        t_1381[k] = f_3 * pc_x[k] * skk_1111[k];

        t_1382[k] = f_3 * pc_x[k] * skk_1112[k];
    }

#pragma omp simd aligned(t_1383, t_1384, t_1385, t_1386, t_1387, pc_x, pc_y, pc_z, sik_820, \
                         sik_856, ski0_861, ski1_861, skk_1108, skk_1113, skk_1114, \
                         skk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1383[k] = f_3 * pc_x[k] * skk_1113[k];

        t_1384[k] = f_3 * pc_x[k] * skk_1114[k];

        t_1385[k] = f_3 * pc_x[k] * skk_1115[k];

        t_1386[k] = f_19 * sik_856[k]
                    + f_1 * ski0_861[k]
                    - f_2 * ski1_861[k]
                    + f_3 * pc_y[k] * skk_1108[k];

        t_1387[k] = f_16 * sik_820[k]
                    + f_3 * pc_z[k] * skk_1108[k];
    }

#pragma omp simd aligned(t_1388, t_1389, t_1390, pc_y, sik_858, sik_859, sik_860, ski0_863, \
                         ski0_864, ski0_865, ski1_863, ski1_864, ski1_865, skk_1110, skk_1111, \
                         skk_1112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1388[k] = f_19 * sik_858[k]
                    + f_4 * ski0_863[k]
                    - f_5 * ski1_863[k]
                    + f_3 * pc_y[k] * skk_1110[k];

        t_1389[k] = f_19 * sik_859[k]
                    + f_6 * ski0_864[k]
                    - f_7 * ski1_864[k]
                    + f_3 * pc_y[k] * skk_1111[k];

        t_1390[k] = f_19 * sik_860[k]
                    + f_8 * ski0_865[k]
                    - f_9 * ski1_865[k]
                    + f_3 * pc_y[k] * skk_1112[k];
    }

#pragma omp simd aligned(t_1391, t_1392, t_1393, pc_y, sik_861, sik_862, sik_863, ski0_866, \
                         ski0_867, ski1_866, ski1_867, skk_1113, skk_1114, \
                         skk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1391[k] = f_19 * sik_861[k]
                    + f_10 * ski0_866[k]
                    - f_11 * ski1_866[k]
                    + f_3 * pc_y[k] * skk_1113[k];

        t_1392[k] = f_19 * sik_862[k]
                    + f_12 * ski0_867[k]
                    - f_13 * ski1_867[k]
                    + f_3 * pc_y[k] * skk_1114[k];

        t_1393[k] = f_19 * sik_863[k]
                    + f_3 * pc_y[k] * skk_1115[k];
    }

#pragma omp simd aligned(t_1394, t_1395, t_1396, t_1397, pc_x, pc_y, pc_z, sik_827, sik_828, \
                         sik_864, ski0_867, ski0_868, ski1_867, ski1_868, skk_1115, \
                         skk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1394[k] = f_16 * sik_827[k]
                    + f_1 * ski0_867[k]
                    - f_2 * ski1_867[k]
                    + f_3 * pc_z[k] * skk_1115[k];

        t_1395[k] = f_1 * ski0_868[k]
                    - f_2 * ski1_868[k]
                    + f_3 * pc_x[k] * skk_1116[k];

        t_1396[k] = f_18 * sik_864[k]
                    + f_3 * pc_y[k] * skk_1116[k];

        t_1397[k] = f_17 * sik_828[k]
                    + f_3 * pc_z[k] * skk_1116[k];
    }

#pragma omp simd aligned(t_1398, t_1399, t_1400, pc_x, pc_y, sik_866, ski0_871, ski0_873, \
                         ski1_871, ski1_873, skk_1118, skk_1119, \
                         skk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1398[k] = f_4 * ski0_871[k]
                    - f_5 * ski1_871[k]
                    + f_3 * pc_x[k] * skk_1119[k];

        t_1399[k] = f_18 * sik_866[k]
                    + f_3 * pc_y[k] * skk_1118[k];

        t_1400[k] = f_4 * ski0_873[k]
                    - f_5 * ski1_873[k]
                    + f_3 * pc_x[k] * skk_1121[k];
    }

#pragma omp simd aligned(t_1401, t_1402, t_1403, pc_x, pc_y, pc_z, sik_831, sik_869, ski0_874, \
                         ski1_874, skk_1119, skk_1121, skk_1122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1401[k] = f_6 * ski0_874[k]
                    - f_7 * ski1_874[k]
                    + f_3 * pc_x[k] * skk_1122[k];

        t_1402[k] = f_17 * sik_831[k]
                    + f_3 * pc_z[k] * skk_1119[k];

        t_1403[k] = f_18 * sik_869[k]
                    + f_3 * pc_y[k] * skk_1121[k];
    }

#pragma omp simd aligned(t_1404, t_1405, t_1406, pc_x, pc_z, sik_834, ski0_877, ski0_878, \
                         ski1_877, ski1_878, skk_1122, skk_1125, \
                         skk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1404[k] = f_6 * ski0_877[k]
                    - f_7 * ski1_877[k]
                    + f_3 * pc_x[k] * skk_1125[k];

        t_1405[k] = f_8 * ski0_878[k]
                    - f_9 * ski1_878[k]
                    + f_3 * pc_x[k] * skk_1126[k];

        t_1406[k] = f_17 * sik_834[k]
                    + f_3 * pc_z[k] * skk_1122[k];
    }

#pragma omp simd aligned(t_1407, t_1408, t_1409, pc_x, pc_y, sik_873, ski0_880, ski0_882, \
                         ski1_880, ski1_882, skk_1125, skk_1128, \
                         skk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1407[k] = f_8 * ski0_880[k]
                    - f_9 * ski1_880[k]
                    + f_3 * pc_x[k] * skk_1128[k];

        t_1408[k] = f_18 * sik_873[k]
                    + f_3 * pc_y[k] * skk_1125[k];

        t_1409[k] = f_8 * ski0_882[k]
                    - f_9 * ski1_882[k]
                    + f_3 * pc_x[k] * skk_1130[k];
    }

#pragma omp simd aligned(t_1410, t_1411, t_1412, pc_x, pc_z, sik_838, ski0_883, ski0_885, \
                         ski1_883, ski1_885, skk_1126, skk_1131, \
                         skk_1133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1410[k] = f_10 * ski0_883[k]
                    - f_11 * ski1_883[k]
                    + f_3 * pc_x[k] * skk_1131[k];

        t_1411[k] = f_17 * sik_838[k]
                    + f_3 * pc_z[k] * skk_1126[k];

        t_1412[k] = f_10 * ski0_885[k]
                    - f_11 * ski1_885[k]
                    + f_3 * pc_x[k] * skk_1133[k];
    }

#pragma omp simd aligned(t_1413, t_1414, t_1415, pc_x, pc_y, sik_878, ski0_886, ski0_888, \
                         ski1_886, ski1_888, skk_1130, skk_1134, \
                         skk_1136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1413[k] = f_10 * ski0_886[k]
                    - f_11 * ski1_886[k]
                    + f_3 * pc_x[k] * skk_1134[k];

        t_1414[k] = f_18 * sik_878[k]
                    + f_3 * pc_y[k] * skk_1130[k];

        t_1415[k] = f_10 * ski0_888[k]
                    - f_11 * ski1_888[k]
                    + f_3 * pc_x[k] * skk_1136[k];
    }

#pragma omp simd aligned(t_1416, t_1417, t_1418, pc_x, pc_z, sik_843, ski0_889, ski0_891, \
                         ski1_889, ski1_891, skk_1131, skk_1137, \
                         skk_1139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1416[k] = f_12 * ski0_889[k]
                    - f_13 * ski1_889[k]
                    + f_3 * pc_x[k] * skk_1137[k];

        t_1417[k] = f_17 * sik_843[k]
                    + f_3 * pc_z[k] * skk_1131[k];

        t_1418[k] = f_12 * ski0_891[k]
                    - f_13 * ski1_891[k]
                    + f_3 * pc_x[k] * skk_1139[k];
    }

#pragma omp simd aligned(t_1419, t_1420, t_1421, pc_x, pc_y, sik_884, ski0_892, ski0_893, \
                         ski1_892, ski1_893, skk_1136, skk_1140, \
                         skk_1141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1419[k] = f_12 * ski0_892[k]
                    - f_13 * ski1_892[k]
                    + f_3 * pc_x[k] * skk_1140[k];

        t_1420[k] = f_12 * ski0_893[k]
                    - f_13 * ski1_893[k]
                    + f_3 * pc_x[k] * skk_1141[k];

        t_1421[k] = f_18 * sik_884[k]
                    + f_3 * pc_y[k] * skk_1136[k];
    }

#pragma omp simd aligned(t_1422, t_1423, t_1424, t_1425, t_1426, t_1427, pc_x, ski0_895, \
                         ski1_895, skk_1143, skk_1144, skk_1145, skk_1146, skk_1147, \
                         skk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1422[k] = f_12 * ski0_895[k]
                    - f_13 * ski1_895[k]
                    + f_3 * pc_x[k] * skk_1143[k];

        t_1423[k] = f_3 * pc_x[k] * skk_1144[k];

        t_1424[k] = f_3 * pc_x[k] * skk_1145[k];

        t_1425[k] = f_3 * pc_x[k] * skk_1146[k];

        t_1426[k] = f_3 * pc_x[k] * skk_1147[k];

        t_1427[k] = f_3 * pc_x[k] * skk_1148[k];
    }

#pragma omp simd aligned(t_1428, t_1429, t_1430, t_1431, t_1432, pc_x, pc_y, pc_z, sik_856, \
                         sik_892, ski0_889, ski1_889, skk_1144, skk_1149, skk_1150, \
                         skk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1428[k] = f_3 * pc_x[k] * skk_1149[k];

        t_1429[k] = f_3 * pc_x[k] * skk_1150[k];

        t_1430[k] = f_3 * pc_x[k] * skk_1151[k];

        t_1431[k] = f_18 * sik_892[k]
                    + f_1 * ski0_889[k]
                    - f_2 * ski1_889[k]
                    + f_3 * pc_y[k] * skk_1144[k];

        t_1432[k] = f_17 * sik_856[k]
                    + f_3 * pc_z[k] * skk_1144[k];
    }

#pragma omp simd aligned(t_1433, t_1434, t_1435, pc_y, sik_894, sik_895, sik_896, ski0_891, \
                         ski0_892, ski0_893, ski1_891, ski1_892, ski1_893, skk_1146, skk_1147, \
                         skk_1148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1433[k] = f_18 * sik_894[k]
                    + f_4 * ski0_891[k]
                    - f_5 * ski1_891[k]
                    + f_3 * pc_y[k] * skk_1146[k];

        t_1434[k] = f_18 * sik_895[k]
                    + f_6 * ski0_892[k]
                    - f_7 * ski1_892[k]
                    + f_3 * pc_y[k] * skk_1147[k];

        t_1435[k] = f_18 * sik_896[k]
                    + f_8 * ski0_893[k]
                    - f_9 * ski1_893[k]
                    + f_3 * pc_y[k] * skk_1148[k];
    }

#pragma omp simd aligned(t_1436, t_1437, t_1438, pc_y, sik_897, sik_898, sik_899, ski0_894, \
                         ski0_895, ski1_894, ski1_895, skk_1149, skk_1150, \
                         skk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1436[k] = f_18 * sik_897[k]
                    + f_10 * ski0_894[k]
                    - f_11 * ski1_894[k]
                    + f_3 * pc_y[k] * skk_1149[k];

        t_1437[k] = f_18 * sik_898[k]
                    + f_12 * ski0_895[k]
                    - f_13 * ski1_895[k]
                    + f_3 * pc_y[k] * skk_1150[k];

        t_1438[k] = f_18 * sik_899[k]
                    + f_3 * pc_y[k] * skk_1151[k];
    }

#pragma omp simd aligned(t_1439, t_1440, t_1441, t_1442, pc_x, pc_y, pc_z, sik_863, sik_864, \
                         sik_900, ski0_895, ski0_896, ski1_895, ski1_896, skk_1151, \
                         skk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1439[k] = f_17 * sik_863[k]
                    + f_1 * ski0_895[k]
                    - f_2 * ski1_895[k]
                    + f_3 * pc_z[k] * skk_1151[k];

        t_1440[k] = f_1 * ski0_896[k]
                    - f_2 * ski1_896[k]
                    + f_3 * pc_x[k] * skk_1152[k];

        t_1441[k] = f_17 * sik_900[k]
                    + f_3 * pc_y[k] * skk_1152[k];

        t_1442[k] = f_18 * sik_864[k]
                    + f_3 * pc_z[k] * skk_1152[k];
    }

#pragma omp simd aligned(t_1443, t_1444, t_1445, pc_x, pc_y, sik_902, ski0_899, ski0_901, \
                         ski1_899, ski1_901, skk_1154, skk_1155, \
                         skk_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1443[k] = f_4 * ski0_899[k]
                    - f_5 * ski1_899[k]
                    + f_3 * pc_x[k] * skk_1155[k];

        t_1444[k] = f_17 * sik_902[k]
                    + f_3 * pc_y[k] * skk_1154[k];

        t_1445[k] = f_4 * ski0_901[k]
                    - f_5 * ski1_901[k]
                    + f_3 * pc_x[k] * skk_1157[k];
    }

#pragma omp simd aligned(t_1446, t_1447, t_1448, pc_x, pc_y, pc_z, sik_867, sik_905, ski0_902, \
                         ski1_902, skk_1155, skk_1157, skk_1158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1446[k] = f_6 * ski0_902[k]
                    - f_7 * ski1_902[k]
                    + f_3 * pc_x[k] * skk_1158[k];

        t_1447[k] = f_18 * sik_867[k]
                    + f_3 * pc_z[k] * skk_1155[k];

        t_1448[k] = f_17 * sik_905[k]
                    + f_3 * pc_y[k] * skk_1157[k];
    }

#pragma omp simd aligned(t_1449, t_1450, t_1451, pc_x, pc_z, sik_870, ski0_905, ski0_906, \
                         ski1_905, ski1_906, skk_1158, skk_1161, \
                         skk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1449[k] = f_6 * ski0_905[k]
                    - f_7 * ski1_905[k]
                    + f_3 * pc_x[k] * skk_1161[k];

        t_1450[k] = f_8 * ski0_906[k]
                    - f_9 * ski1_906[k]
                    + f_3 * pc_x[k] * skk_1162[k];

        t_1451[k] = f_18 * sik_870[k]
                    + f_3 * pc_z[k] * skk_1158[k];
    }

#pragma omp simd aligned(t_1452, t_1453, t_1454, pc_x, pc_y, sik_909, ski0_908, ski0_910, \
                         ski1_908, ski1_910, skk_1161, skk_1164, \
                         skk_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1452[k] = f_8 * ski0_908[k]
                    - f_9 * ski1_908[k]
                    + f_3 * pc_x[k] * skk_1164[k];

        t_1453[k] = f_17 * sik_909[k]
                    + f_3 * pc_y[k] * skk_1161[k];

        t_1454[k] = f_8 * ski0_910[k]
                    - f_9 * ski1_910[k]
                    + f_3 * pc_x[k] * skk_1166[k];
    }

#pragma omp simd aligned(t_1455, t_1456, t_1457, pc_x, pc_z, sik_874, ski0_911, ski0_913, \
                         ski1_911, ski1_913, skk_1162, skk_1167, \
                         skk_1169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1455[k] = f_10 * ski0_911[k]
                    - f_11 * ski1_911[k]
                    + f_3 * pc_x[k] * skk_1167[k];

        t_1456[k] = f_18 * sik_874[k]
                    + f_3 * pc_z[k] * skk_1162[k];

        t_1457[k] = f_10 * ski0_913[k]
                    - f_11 * ski1_913[k]
                    + f_3 * pc_x[k] * skk_1169[k];
    }

#pragma omp simd aligned(t_1458, t_1459, t_1460, pc_x, pc_y, sik_914, ski0_914, ski0_916, \
                         ski1_914, ski1_916, skk_1166, skk_1170, \
                         skk_1172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1458[k] = f_10 * ski0_914[k]
                    - f_11 * ski1_914[k]
                    + f_3 * pc_x[k] * skk_1170[k];

        t_1459[k] = f_17 * sik_914[k]
                    + f_3 * pc_y[k] * skk_1166[k];

        t_1460[k] = f_10 * ski0_916[k]
                    - f_11 * ski1_916[k]
                    + f_3 * pc_x[k] * skk_1172[k];
    }

#pragma omp simd aligned(t_1461, t_1462, t_1463, pc_x, pc_z, sik_879, ski0_917, ski0_919, \
                         ski1_917, ski1_919, skk_1167, skk_1173, \
                         skk_1175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1461[k] = f_12 * ski0_917[k]
                    - f_13 * ski1_917[k]
                    + f_3 * pc_x[k] * skk_1173[k];

        t_1462[k] = f_18 * sik_879[k]
                    + f_3 * pc_z[k] * skk_1167[k];

        t_1463[k] = f_12 * ski0_919[k]
                    - f_13 * ski1_919[k]
                    + f_3 * pc_x[k] * skk_1175[k];
    }

#pragma omp simd aligned(t_1464, t_1465, t_1466, pc_x, pc_y, sik_920, ski0_920, ski0_921, \
                         ski1_920, ski1_921, skk_1172, skk_1176, \
                         skk_1177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1464[k] = f_12 * ski0_920[k]
                    - f_13 * ski1_920[k]
                    + f_3 * pc_x[k] * skk_1176[k];

        t_1465[k] = f_12 * ski0_921[k]
                    - f_13 * ski1_921[k]
                    + f_3 * pc_x[k] * skk_1177[k];

        t_1466[k] = f_17 * sik_920[k]
                    + f_3 * pc_y[k] * skk_1172[k];
    }

#pragma omp simd aligned(t_1467, t_1468, t_1469, t_1470, t_1471, t_1472, pc_x, ski0_923, \
                         ski1_923, skk_1179, skk_1180, skk_1181, skk_1182, skk_1183, \
                         skk_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1467[k] = f_12 * ski0_923[k]
                    - f_13 * ski1_923[k]
                    + f_3 * pc_x[k] * skk_1179[k];

        t_1468[k] = f_3 * pc_x[k] * skk_1180[k];

        t_1469[k] = f_3 * pc_x[k] * skk_1181[k];

        t_1470[k] = f_3 * pc_x[k] * skk_1182[k];

        t_1471[k] = f_3 * pc_x[k] * skk_1183[k];

        t_1472[k] = f_3 * pc_x[k] * skk_1184[k];
    }

#pragma omp simd aligned(t_1473, t_1474, t_1475, t_1476, t_1477, pc_x, pc_y, pc_z, sik_892, \
                         sik_928, ski0_917, ski1_917, skk_1180, skk_1185, skk_1186, \
                         skk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1473[k] = f_3 * pc_x[k] * skk_1185[k];

        t_1474[k] = f_3 * pc_x[k] * skk_1186[k];

        t_1475[k] = f_3 * pc_x[k] * skk_1187[k];

        t_1476[k] = f_17 * sik_928[k]
                    + f_1 * ski0_917[k]
                    - f_2 * ski1_917[k]
                    + f_3 * pc_y[k] * skk_1180[k];

        t_1477[k] = f_18 * sik_892[k]
                    + f_3 * pc_z[k] * skk_1180[k];
    }

#pragma omp simd aligned(t_1478, t_1479, t_1480, pc_y, sik_930, sik_931, sik_932, ski0_919, \
                         ski0_920, ski0_921, ski1_919, ski1_920, ski1_921, skk_1182, skk_1183, \
                         skk_1184 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1478[k] = f_17 * sik_930[k]
                    + f_4 * ski0_919[k]
                    - f_5 * ski1_919[k]
                    + f_3 * pc_y[k] * skk_1182[k];

        t_1479[k] = f_17 * sik_931[k]
                    + f_6 * ski0_920[k]
                    - f_7 * ski1_920[k]
                    + f_3 * pc_y[k] * skk_1183[k];

        t_1480[k] = f_17 * sik_932[k]
                    + f_8 * ski0_921[k]
                    - f_9 * ski1_921[k]
                    + f_3 * pc_y[k] * skk_1184[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t sil0,
                                                           const size_t sik, const size_t sil1,
                                                           const size_t ski0, const size_t ski1,
                                                           const size_t skk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;
    const auto f_21 = 4.0 / q;

    auto *t_1481 = buffer.data(target + 1481);
    auto *t_1482 = buffer.data(target + 1482);
    auto *t_1483 = buffer.data(target + 1483);
    auto *t_1484 = buffer.data(target + 1484);
    auto *t_1485 = buffer.data(target + 1485);
    auto *t_1486 = buffer.data(target + 1486);
    auto *t_1487 = buffer.data(target + 1487);
    auto *t_1488 = buffer.data(target + 1488);
    auto *t_1489 = buffer.data(target + 1489);
    auto *t_1490 = buffer.data(target + 1490);
    auto *t_1491 = buffer.data(target + 1491);
    auto *t_1492 = buffer.data(target + 1492);
    auto *t_1493 = buffer.data(target + 1493);
    auto *t_1494 = buffer.data(target + 1494);
    auto *t_1495 = buffer.data(target + 1495);
    auto *t_1496 = buffer.data(target + 1496);
    auto *t_1497 = buffer.data(target + 1497);
    auto *t_1498 = buffer.data(target + 1498);
    auto *t_1499 = buffer.data(target + 1499);
    auto *t_1500 = buffer.data(target + 1500);
    auto *t_1501 = buffer.data(target + 1501);
    auto *t_1502 = buffer.data(target + 1502);
    auto *t_1503 = buffer.data(target + 1503);
    auto *t_1504 = buffer.data(target + 1504);
    auto *t_1505 = buffer.data(target + 1505);
    auto *t_1506 = buffer.data(target + 1506);
    auto *t_1507 = buffer.data(target + 1507);
    auto *t_1508 = buffer.data(target + 1508);
    auto *t_1509 = buffer.data(target + 1509);
    auto *t_1510 = buffer.data(target + 1510);
    auto *t_1511 = buffer.data(target + 1511);
    auto *t_1512 = buffer.data(target + 1512);
    auto *t_1513 = buffer.data(target + 1513);
    auto *t_1514 = buffer.data(target + 1514);
    auto *t_1515 = buffer.data(target + 1515);
    auto *t_1516 = buffer.data(target + 1516);
    auto *t_1517 = buffer.data(target + 1517);
    auto *t_1518 = buffer.data(target + 1518);
    auto *t_1519 = buffer.data(target + 1519);
    auto *t_1520 = buffer.data(target + 1520);
    auto *t_1521 = buffer.data(target + 1521);
    auto *t_1522 = buffer.data(target + 1522);
    auto *t_1523 = buffer.data(target + 1523);
    auto *t_1524 = buffer.data(target + 1524);
    auto *t_1525 = buffer.data(target + 1525);
    auto *t_1526 = buffer.data(target + 1526);
    auto *t_1527 = buffer.data(target + 1527);
    auto *t_1528 = buffer.data(target + 1528);
    auto *t_1529 = buffer.data(target + 1529);
    auto *t_1530 = buffer.data(target + 1530);
    auto *t_1531 = buffer.data(target + 1531);
    auto *t_1532 = buffer.data(target + 1532);
    auto *t_1533 = buffer.data(target + 1533);
    auto *t_1534 = buffer.data(target + 1534);
    auto *t_1535 = buffer.data(target + 1535);
    auto *t_1536 = buffer.data(target + 1536);
    auto *t_1537 = buffer.data(target + 1537);
    auto *t_1538 = buffer.data(target + 1538);
    auto *t_1539 = buffer.data(target + 1539);
    auto *t_1540 = buffer.data(target + 1540);
    auto *t_1541 = buffer.data(target + 1541);
    auto *t_1542 = buffer.data(target + 1542);
    auto *t_1543 = buffer.data(target + 1543);
    auto *t_1544 = buffer.data(target + 1544);
    auto *t_1545 = buffer.data(target + 1545);
    auto *t_1546 = buffer.data(target + 1546);
    auto *t_1547 = buffer.data(target + 1547);
    auto *t_1548 = buffer.data(target + 1548);
    auto *t_1549 = buffer.data(target + 1549);
    auto *t_1550 = buffer.data(target + 1550);
    auto *t_1551 = buffer.data(target + 1551);
    auto *t_1552 = buffer.data(target + 1552);
    auto *t_1553 = buffer.data(target + 1553);
    auto *t_1554 = buffer.data(target + 1554);
    auto *t_1555 = buffer.data(target + 1555);
    auto *t_1556 = buffer.data(target + 1556);
    auto *t_1557 = buffer.data(target + 1557);
    auto *t_1558 = buffer.data(target + 1558);
    auto *t_1559 = buffer.data(target + 1559);
    auto *t_1560 = buffer.data(target + 1560);
    auto *t_1561 = buffer.data(target + 1561);
    auto *t_1562 = buffer.data(target + 1562);
    auto *t_1563 = buffer.data(target + 1563);
    auto *t_1564 = buffer.data(target + 1564);
    auto *t_1565 = buffer.data(target + 1565);
    auto *t_1566 = buffer.data(target + 1566);
    auto *t_1567 = buffer.data(target + 1567);
    auto *t_1568 = buffer.data(target + 1568);
    auto *t_1569 = buffer.data(target + 1569);
    auto *t_1570 = buffer.data(target + 1570);
    auto *t_1571 = buffer.data(target + 1571);
    auto *t_1572 = buffer.data(target + 1572);
    auto *t_1573 = buffer.data(target + 1573);
    auto *t_1574 = buffer.data(target + 1574);
    auto *t_1575 = buffer.data(target + 1575);
    auto *t_1576 = buffer.data(target + 1576);
    auto *t_1577 = buffer.data(target + 1577);
    auto *t_1578 = buffer.data(target + 1578);
    auto *t_1579 = buffer.data(target + 1579);
    auto *t_1580 = buffer.data(target + 1580);
    auto *t_1581 = buffer.data(target + 1581);
    auto *t_1582 = buffer.data(target + 1582);
    auto *t_1583 = buffer.data(target + 1583);
    auto *t_1584 = buffer.data(target + 1584);
    auto *t_1585 = buffer.data(target + 1585);
    auto *t_1586 = buffer.data(target + 1586);
    auto *t_1587 = buffer.data(target + 1587);
    auto *t_1588 = buffer.data(target + 1588);
    auto *t_1589 = buffer.data(target + 1589);
    auto *t_1590 = buffer.data(target + 1590);
    auto *t_1591 = buffer.data(target + 1591);
    auto *t_1592 = buffer.data(target + 1592);
    auto *t_1593 = buffer.data(target + 1593);
    auto *t_1594 = buffer.data(target + 1594);
    auto *t_1595 = buffer.data(target + 1595);
    auto *t_1596 = buffer.data(target + 1596);
    auto *t_1597 = buffer.data(target + 1597);
    auto *t_1598 = buffer.data(target + 1598);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sil0_1215 = buffer.data(sil0 + 1215);
    const auto *sil0_1220 = buffer.data(sil0 + 1220);
    const auto *sil0_1224 = buffer.data(sil0 + 1224);
    const auto *sil0_1229 = buffer.data(sil0 + 1229);
    const auto *sil0_1235 = buffer.data(sil0 + 1235);
    const auto *sil0_1242 = buffer.data(sil0 + 1242);
    const auto *sil0_1251 = buffer.data(sil0 + 1251);
    const auto *sil0_1253 = buffer.data(sil0 + 1253);
    const auto *sil0_1254 = buffer.data(sil0 + 1254);
    const auto *sil0_1255 = buffer.data(sil0 + 1255);
    const auto *sil0_1256 = buffer.data(sil0 + 1256);
    const auto *sil0_1257 = buffer.data(sil0 + 1257);
    const auto *sil0_1259 = buffer.data(sil0 + 1259);

    const auto *sik_899 = buffer.data(sik + 899);
    const auto *sik_900 = buffer.data(sik + 900);
    const auto *sik_903 = buffer.data(sik + 903);
    const auto *sik_906 = buffer.data(sik + 906);
    const auto *sik_910 = buffer.data(sik + 910);
    const auto *sik_915 = buffer.data(sik + 915);
    const auto *sik_928 = buffer.data(sik + 928);
    const auto *sik_933 = buffer.data(sik + 933);
    const auto *sik_934 = buffer.data(sik + 934);
    const auto *sik_935 = buffer.data(sik + 935);
    const auto *sik_936 = buffer.data(sik + 936);
    const auto *sik_938 = buffer.data(sik + 938);
    const auto *sik_939 = buffer.data(sik + 939);
    const auto *sik_941 = buffer.data(sik + 941);
    const auto *sik_942 = buffer.data(sik + 942);
    const auto *sik_945 = buffer.data(sik + 945);
    const auto *sik_946 = buffer.data(sik + 946);
    const auto *sik_950 = buffer.data(sik + 950);
    const auto *sik_951 = buffer.data(sik + 951);
    const auto *sik_956 = buffer.data(sik + 956);
    const auto *sik_964 = buffer.data(sik + 964);
    const auto *sik_966 = buffer.data(sik + 966);
    const auto *sik_967 = buffer.data(sik + 967);
    const auto *sik_968 = buffer.data(sik + 968);
    const auto *sik_969 = buffer.data(sik + 969);
    const auto *sik_970 = buffer.data(sik + 970);
    const auto *sik_971 = buffer.data(sik + 971);
    const auto *sik_972 = buffer.data(sik + 972);
    const auto *sik_974 = buffer.data(sik + 974);
    const auto *sik_975 = buffer.data(sik + 975);
    const auto *sik_977 = buffer.data(sik + 977);
    const auto *sik_978 = buffer.data(sik + 978);
    const auto *sik_981 = buffer.data(sik + 981);
    const auto *sik_982 = buffer.data(sik + 982);
    const auto *sik_986 = buffer.data(sik + 986);
    const auto *sik_987 = buffer.data(sik + 987);
    const auto *sik_992 = buffer.data(sik + 992);
    const auto *sik_1000 = buffer.data(sik + 1000);
    const auto *sik_1002 = buffer.data(sik + 1002);
    const auto *sik_1003 = buffer.data(sik + 1003);
    const auto *sik_1004 = buffer.data(sik + 1004);
    const auto *sik_1005 = buffer.data(sik + 1005);
    const auto *sik_1006 = buffer.data(sik + 1006);
    const auto *sik_1007 = buffer.data(sik + 1007);

    const auto *sil1_1215 = buffer.data(sil1 + 1215);
    const auto *sil1_1220 = buffer.data(sil1 + 1220);
    const auto *sil1_1224 = buffer.data(sil1 + 1224);
    const auto *sil1_1229 = buffer.data(sil1 + 1229);
    const auto *sil1_1235 = buffer.data(sil1 + 1235);
    const auto *sil1_1242 = buffer.data(sil1 + 1242);
    const auto *sil1_1251 = buffer.data(sil1 + 1251);
    const auto *sil1_1253 = buffer.data(sil1 + 1253);
    const auto *sil1_1254 = buffer.data(sil1 + 1254);
    const auto *sil1_1255 = buffer.data(sil1 + 1255);
    const auto *sil1_1256 = buffer.data(sil1 + 1256);
    const auto *sil1_1257 = buffer.data(sil1 + 1257);
    const auto *sil1_1259 = buffer.data(sil1 + 1259);

    const auto *ski0_922 = buffer.data(ski0 + 922);
    const auto *ski0_923 = buffer.data(ski0 + 923);
    const auto *ski0_924 = buffer.data(ski0 + 924);
    const auto *ski0_927 = buffer.data(ski0 + 927);
    const auto *ski0_929 = buffer.data(ski0 + 929);
    const auto *ski0_930 = buffer.data(ski0 + 930);
    const auto *ski0_933 = buffer.data(ski0 + 933);
    const auto *ski0_934 = buffer.data(ski0 + 934);
    const auto *ski0_936 = buffer.data(ski0 + 936);
    const auto *ski0_938 = buffer.data(ski0 + 938);
    const auto *ski0_939 = buffer.data(ski0 + 939);
    const auto *ski0_941 = buffer.data(ski0 + 941);
    const auto *ski0_942 = buffer.data(ski0 + 942);
    const auto *ski0_944 = buffer.data(ski0 + 944);
    const auto *ski0_945 = buffer.data(ski0 + 945);
    const auto *ski0_947 = buffer.data(ski0 + 947);
    const auto *ski0_948 = buffer.data(ski0 + 948);
    const auto *ski0_949 = buffer.data(ski0 + 949);
    const auto *ski0_950 = buffer.data(ski0 + 950);
    const auto *ski0_951 = buffer.data(ski0 + 951);
    const auto *ski0_955 = buffer.data(ski0 + 955);
    const auto *ski0_958 = buffer.data(ski0 + 958);
    const auto *ski0_962 = buffer.data(ski0 + 962);
    const auto *ski0_964 = buffer.data(ski0 + 964);
    const auto *ski0_967 = buffer.data(ski0 + 967);
    const auto *ski0_969 = buffer.data(ski0 + 969);
    const auto *ski0_970 = buffer.data(ski0 + 970);
    const auto *ski0_973 = buffer.data(ski0 + 973);
    const auto *ski0_975 = buffer.data(ski0 + 975);
    const auto *ski0_976 = buffer.data(ski0 + 976);
    const auto *ski0_977 = buffer.data(ski0 + 977);
    const auto *ski0_980 = buffer.data(ski0 + 980);
    const auto *ski0_983 = buffer.data(ski0 + 983);
    const auto *ski0_985 = buffer.data(ski0 + 985);
    const auto *ski0_986 = buffer.data(ski0 + 986);
    const auto *ski0_989 = buffer.data(ski0 + 989);
    const auto *ski0_990 = buffer.data(ski0 + 990);
    const auto *ski0_992 = buffer.data(ski0 + 992);
    const auto *ski0_994 = buffer.data(ski0 + 994);
    const auto *ski0_995 = buffer.data(ski0 + 995);
    const auto *ski0_997 = buffer.data(ski0 + 997);
    const auto *ski0_998 = buffer.data(ski0 + 998);
    const auto *ski0_1000 = buffer.data(ski0 + 1000);
    const auto *ski0_1001 = buffer.data(ski0 + 1001);
    const auto *ski0_1003 = buffer.data(ski0 + 1003);

    const auto *ski1_922 = buffer.data(ski1 + 922);
    const auto *ski1_923 = buffer.data(ski1 + 923);
    const auto *ski1_924 = buffer.data(ski1 + 924);
    const auto *ski1_927 = buffer.data(ski1 + 927);
    const auto *ski1_929 = buffer.data(ski1 + 929);
    const auto *ski1_930 = buffer.data(ski1 + 930);
    const auto *ski1_933 = buffer.data(ski1 + 933);
    const auto *ski1_934 = buffer.data(ski1 + 934);
    const auto *ski1_936 = buffer.data(ski1 + 936);
    const auto *ski1_938 = buffer.data(ski1 + 938);
    const auto *ski1_939 = buffer.data(ski1 + 939);
    const auto *ski1_941 = buffer.data(ski1 + 941);
    const auto *ski1_942 = buffer.data(ski1 + 942);
    const auto *ski1_944 = buffer.data(ski1 + 944);
    const auto *ski1_945 = buffer.data(ski1 + 945);
    const auto *ski1_947 = buffer.data(ski1 + 947);
    const auto *ski1_948 = buffer.data(ski1 + 948);
    const auto *ski1_949 = buffer.data(ski1 + 949);
    const auto *ski1_950 = buffer.data(ski1 + 950);
    const auto *ski1_951 = buffer.data(ski1 + 951);
    const auto *ski1_955 = buffer.data(ski1 + 955);
    const auto *ski1_958 = buffer.data(ski1 + 958);
    const auto *ski1_962 = buffer.data(ski1 + 962);
    const auto *ski1_964 = buffer.data(ski1 + 964);
    const auto *ski1_967 = buffer.data(ski1 + 967);
    const auto *ski1_969 = buffer.data(ski1 + 969);
    const auto *ski1_970 = buffer.data(ski1 + 970);
    const auto *ski1_973 = buffer.data(ski1 + 973);
    const auto *ski1_975 = buffer.data(ski1 + 975);
    const auto *ski1_976 = buffer.data(ski1 + 976);
    const auto *ski1_977 = buffer.data(ski1 + 977);
    const auto *ski1_980 = buffer.data(ski1 + 980);
    const auto *ski1_983 = buffer.data(ski1 + 983);
    const auto *ski1_985 = buffer.data(ski1 + 985);
    const auto *ski1_986 = buffer.data(ski1 + 986);
    const auto *ski1_989 = buffer.data(ski1 + 989);
    const auto *ski1_990 = buffer.data(ski1 + 990);
    const auto *ski1_992 = buffer.data(ski1 + 992);
    const auto *ski1_994 = buffer.data(ski1 + 994);
    const auto *ski1_995 = buffer.data(ski1 + 995);
    const auto *ski1_997 = buffer.data(ski1 + 997);
    const auto *ski1_998 = buffer.data(ski1 + 998);
    const auto *ski1_1000 = buffer.data(ski1 + 1000);
    const auto *ski1_1001 = buffer.data(ski1 + 1001);
    const auto *ski1_1003 = buffer.data(ski1 + 1003);

    const auto *skk_1185 = buffer.data(skk + 1185);
    const auto *skk_1186 = buffer.data(skk + 1186);
    const auto *skk_1187 = buffer.data(skk + 1187);
    const auto *skk_1188 = buffer.data(skk + 1188);
    const auto *skk_1190 = buffer.data(skk + 1190);
    const auto *skk_1191 = buffer.data(skk + 1191);
    const auto *skk_1193 = buffer.data(skk + 1193);
    const auto *skk_1194 = buffer.data(skk + 1194);
    const auto *skk_1197 = buffer.data(skk + 1197);
    const auto *skk_1198 = buffer.data(skk + 1198);
    const auto *skk_1200 = buffer.data(skk + 1200);
    const auto *skk_1202 = buffer.data(skk + 1202);
    const auto *skk_1203 = buffer.data(skk + 1203);
    const auto *skk_1205 = buffer.data(skk + 1205);
    const auto *skk_1206 = buffer.data(skk + 1206);
    const auto *skk_1208 = buffer.data(skk + 1208);
    const auto *skk_1209 = buffer.data(skk + 1209);
    const auto *skk_1211 = buffer.data(skk + 1211);
    const auto *skk_1212 = buffer.data(skk + 1212);
    const auto *skk_1213 = buffer.data(skk + 1213);
    const auto *skk_1215 = buffer.data(skk + 1215);
    const auto *skk_1216 = buffer.data(skk + 1216);
    const auto *skk_1217 = buffer.data(skk + 1217);
    const auto *skk_1218 = buffer.data(skk + 1218);
    const auto *skk_1219 = buffer.data(skk + 1219);
    const auto *skk_1220 = buffer.data(skk + 1220);
    const auto *skk_1221 = buffer.data(skk + 1221);
    const auto *skk_1222 = buffer.data(skk + 1222);
    const auto *skk_1223 = buffer.data(skk + 1223);
    const auto *skk_1224 = buffer.data(skk + 1224);
    const auto *skk_1226 = buffer.data(skk + 1226);
    const auto *skk_1227 = buffer.data(skk + 1227);
    const auto *skk_1229 = buffer.data(skk + 1229);
    const auto *skk_1230 = buffer.data(skk + 1230);
    const auto *skk_1233 = buffer.data(skk + 1233);
    const auto *skk_1234 = buffer.data(skk + 1234);
    const auto *skk_1236 = buffer.data(skk + 1236);
    const auto *skk_1238 = buffer.data(skk + 1238);
    const auto *skk_1239 = buffer.data(skk + 1239);
    const auto *skk_1241 = buffer.data(skk + 1241);
    const auto *skk_1242 = buffer.data(skk + 1242);
    const auto *skk_1244 = buffer.data(skk + 1244);
    const auto *skk_1245 = buffer.data(skk + 1245);
    const auto *skk_1247 = buffer.data(skk + 1247);
    const auto *skk_1248 = buffer.data(skk + 1248);
    const auto *skk_1249 = buffer.data(skk + 1249);
    const auto *skk_1252 = buffer.data(skk + 1252);
    const auto *skk_1253 = buffer.data(skk + 1253);
    const auto *skk_1254 = buffer.data(skk + 1254);
    const auto *skk_1255 = buffer.data(skk + 1255);
    const auto *skk_1256 = buffer.data(skk + 1256);
    const auto *skk_1257 = buffer.data(skk + 1257);
    const auto *skk_1258 = buffer.data(skk + 1258);
    const auto *skk_1259 = buffer.data(skk + 1259);
    const auto *skk_1260 = buffer.data(skk + 1260);
    const auto *skk_1262 = buffer.data(skk + 1262);
    const auto *skk_1263 = buffer.data(skk + 1263);
    const auto *skk_1265 = buffer.data(skk + 1265);
    const auto *skk_1266 = buffer.data(skk + 1266);
    const auto *skk_1269 = buffer.data(skk + 1269);
    const auto *skk_1270 = buffer.data(skk + 1270);
    const auto *skk_1272 = buffer.data(skk + 1272);
    const auto *skk_1274 = buffer.data(skk + 1274);
    const auto *skk_1275 = buffer.data(skk + 1275);
    const auto *skk_1277 = buffer.data(skk + 1277);
    const auto *skk_1278 = buffer.data(skk + 1278);
    const auto *skk_1280 = buffer.data(skk + 1280);
    const auto *skk_1281 = buffer.data(skk + 1281);
    const auto *skk_1283 = buffer.data(skk + 1283);

#pragma omp simd aligned(t_1481, t_1482, t_1483, pc_y, sik_933, sik_934, sik_935, ski0_922, \
                         ski0_923, ski1_922, ski1_923, skk_1185, skk_1186, \
                         skk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1481[k] = f_17 * sik_933[k]
                    + f_10 * ski0_922[k]
                    - f_11 * ski1_922[k]
                    + f_3 * pc_y[k] * skk_1185[k];

        t_1482[k] = f_17 * sik_934[k]
                    + f_12 * ski0_923[k]
                    - f_13 * ski1_923[k]
                    + f_3 * pc_y[k] * skk_1186[k];

        t_1483[k] = f_17 * sik_935[k]
                    + f_3 * pc_y[k] * skk_1187[k];
    }

#pragma omp simd aligned(t_1484, t_1485, t_1486, t_1487, pc_x, pc_y, pc_z, sik_899, sik_900, \
                         sik_936, ski0_923, ski0_924, ski1_923, ski1_924, skk_1187, \
                         skk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1484[k] = f_18 * sik_899[k]
                    + f_1 * ski0_923[k]
                    - f_2 * ski1_923[k]
                    + f_3 * pc_z[k] * skk_1187[k];

        t_1485[k] = f_1 * ski0_924[k]
                    - f_2 * ski1_924[k]
                    + f_3 * pc_x[k] * skk_1188[k];

        t_1486[k] = f_16 * sik_936[k]
                    + f_3 * pc_y[k] * skk_1188[k];

        t_1487[k] = f_19 * sik_900[k]
                    + f_3 * pc_z[k] * skk_1188[k];
    }

#pragma omp simd aligned(t_1488, t_1489, t_1490, pc_x, pc_y, sik_938, ski0_927, ski0_929, \
                         ski1_927, ski1_929, skk_1190, skk_1191, \
                         skk_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1488[k] = f_4 * ski0_927[k]
                    - f_5 * ski1_927[k]
                    + f_3 * pc_x[k] * skk_1191[k];

        t_1489[k] = f_16 * sik_938[k]
                    + f_3 * pc_y[k] * skk_1190[k];

        t_1490[k] = f_4 * ski0_929[k]
                    - f_5 * ski1_929[k]
                    + f_3 * pc_x[k] * skk_1193[k];
    }

#pragma omp simd aligned(t_1491, t_1492, t_1493, pc_x, pc_y, pc_z, sik_903, sik_941, ski0_930, \
                         ski1_930, skk_1191, skk_1193, skk_1194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1491[k] = f_6 * ski0_930[k]
                    - f_7 * ski1_930[k]
                    + f_3 * pc_x[k] * skk_1194[k];

        t_1492[k] = f_19 * sik_903[k]
                    + f_3 * pc_z[k] * skk_1191[k];

        t_1493[k] = f_16 * sik_941[k]
                    + f_3 * pc_y[k] * skk_1193[k];
    }

#pragma omp simd aligned(t_1494, t_1495, t_1496, pc_x, pc_z, sik_906, ski0_933, ski0_934, \
                         ski1_933, ski1_934, skk_1194, skk_1197, \
                         skk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1494[k] = f_6 * ski0_933[k]
                    - f_7 * ski1_933[k]
                    + f_3 * pc_x[k] * skk_1197[k];

        t_1495[k] = f_8 * ski0_934[k]
                    - f_9 * ski1_934[k]
                    + f_3 * pc_x[k] * skk_1198[k];

        t_1496[k] = f_19 * sik_906[k]
                    + f_3 * pc_z[k] * skk_1194[k];
    }

#pragma omp simd aligned(t_1497, t_1498, t_1499, pc_x, pc_y, sik_945, ski0_936, ski0_938, \
                         ski1_936, ski1_938, skk_1197, skk_1200, \
                         skk_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1497[k] = f_8 * ski0_936[k]
                    - f_9 * ski1_936[k]
                    + f_3 * pc_x[k] * skk_1200[k];

        t_1498[k] = f_16 * sik_945[k]
                    + f_3 * pc_y[k] * skk_1197[k];

        t_1499[k] = f_8 * ski0_938[k]
                    - f_9 * ski1_938[k]
                    + f_3 * pc_x[k] * skk_1202[k];
    }

#pragma omp simd aligned(t_1500, t_1501, t_1502, pc_x, pc_z, sik_910, ski0_939, ski0_941, \
                         ski1_939, ski1_941, skk_1198, skk_1203, \
                         skk_1205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1500[k] = f_10 * ski0_939[k]
                    - f_11 * ski1_939[k]
                    + f_3 * pc_x[k] * skk_1203[k];

        t_1501[k] = f_19 * sik_910[k]
                    + f_3 * pc_z[k] * skk_1198[k];

        t_1502[k] = f_10 * ski0_941[k]
                    - f_11 * ski1_941[k]
                    + f_3 * pc_x[k] * skk_1205[k];
    }

#pragma omp simd aligned(t_1503, t_1504, t_1505, pc_x, pc_y, sik_950, ski0_942, ski0_944, \
                         ski1_942, ski1_944, skk_1202, skk_1206, \
                         skk_1208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1503[k] = f_10 * ski0_942[k]
                    - f_11 * ski1_942[k]
                    + f_3 * pc_x[k] * skk_1206[k];

        t_1504[k] = f_16 * sik_950[k]
                    + f_3 * pc_y[k] * skk_1202[k];

        t_1505[k] = f_10 * ski0_944[k]
                    - f_11 * ski1_944[k]
                    + f_3 * pc_x[k] * skk_1208[k];
    }

#pragma omp simd aligned(t_1506, t_1507, t_1508, pc_x, pc_z, sik_915, ski0_945, ski0_947, \
                         ski1_945, ski1_947, skk_1203, skk_1209, \
                         skk_1211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1506[k] = f_12 * ski0_945[k]
                    - f_13 * ski1_945[k]
                    + f_3 * pc_x[k] * skk_1209[k];

        t_1507[k] = f_19 * sik_915[k]
                    + f_3 * pc_z[k] * skk_1203[k];

        t_1508[k] = f_12 * ski0_947[k]
                    - f_13 * ski1_947[k]
                    + f_3 * pc_x[k] * skk_1211[k];
    }

#pragma omp simd aligned(t_1509, t_1510, t_1511, pc_x, pc_y, sik_956, ski0_948, ski0_949, \
                         ski1_948, ski1_949, skk_1208, skk_1212, \
                         skk_1213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1509[k] = f_12 * ski0_948[k]
                    - f_13 * ski1_948[k]
                    + f_3 * pc_x[k] * skk_1212[k];

        t_1510[k] = f_12 * ski0_949[k]
                    - f_13 * ski1_949[k]
                    + f_3 * pc_x[k] * skk_1213[k];

        t_1511[k] = f_16 * sik_956[k]
                    + f_3 * pc_y[k] * skk_1208[k];
    }

#pragma omp simd aligned(t_1512, t_1513, t_1514, t_1515, t_1516, t_1517, pc_x, ski0_951, \
                         ski1_951, skk_1215, skk_1216, skk_1217, skk_1218, skk_1219, \
                         skk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1512[k] = f_12 * ski0_951[k]
                    - f_13 * ski1_951[k]
                    + f_3 * pc_x[k] * skk_1215[k];

        t_1513[k] = f_3 * pc_x[k] * skk_1216[k];

        t_1514[k] = f_3 * pc_x[k] * skk_1217[k];

        t_1515[k] = f_3 * pc_x[k] * skk_1218[k];

        t_1516[k] = f_3 * pc_x[k] * skk_1219[k];

        t_1517[k] = f_3 * pc_x[k] * skk_1220[k];
    }

#pragma omp simd aligned(t_1518, t_1519, t_1520, t_1521, t_1522, pc_x, pc_y, pc_z, sik_928, \
                         sik_964, ski0_945, ski1_945, skk_1216, skk_1221, skk_1222, \
                         skk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1518[k] = f_3 * pc_x[k] * skk_1221[k];

        t_1519[k] = f_3 * pc_x[k] * skk_1222[k];

        t_1520[k] = f_3 * pc_x[k] * skk_1223[k];

        t_1521[k] = f_16 * sik_964[k]
                    + f_1 * ski0_945[k]
                    - f_2 * ski1_945[k]
                    + f_3 * pc_y[k] * skk_1216[k];

        t_1522[k] = f_19 * sik_928[k]
                    + f_3 * pc_z[k] * skk_1216[k];
    }

#pragma omp simd aligned(t_1523, t_1524, t_1525, pc_y, sik_966, sik_967, sik_968, ski0_947, \
                         ski0_948, ski0_949, ski1_947, ski1_948, ski1_949, skk_1218, skk_1219, \
                         skk_1220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1523[k] = f_16 * sik_966[k]
                    + f_4 * ski0_947[k]
                    - f_5 * ski1_947[k]
                    + f_3 * pc_y[k] * skk_1218[k];

        t_1524[k] = f_16 * sik_967[k]
                    + f_6 * ski0_948[k]
                    - f_7 * ski1_948[k]
                    + f_3 * pc_y[k] * skk_1219[k];

        t_1525[k] = f_16 * sik_968[k]
                    + f_8 * ski0_949[k]
                    - f_9 * ski1_949[k]
                    + f_3 * pc_y[k] * skk_1220[k];
    }

#pragma omp simd aligned(t_1526, t_1527, t_1528, pc_y, sik_969, sik_970, sik_971, ski0_950, \
                         ski0_951, ski1_950, ski1_951, skk_1221, skk_1222, \
                         skk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1526[k] = f_16 * sik_969[k]
                    + f_10 * ski0_950[k]
                    - f_11 * ski1_950[k]
                    + f_3 * pc_y[k] * skk_1221[k];

        t_1527[k] = f_16 * sik_970[k]
                    + f_12 * ski0_951[k]
                    - f_13 * ski1_951[k]
                    + f_3 * pc_y[k] * skk_1222[k];

        t_1528[k] = f_16 * sik_971[k]
                    + f_3 * pc_y[k] * skk_1223[k];
    }

#pragma omp simd aligned(t_1529, t_1530, t_1531, t_1532, pb_y, pc_y, pc_z, sil0_1215, sik_935, \
                         sik_936, sik_972, sil1_1215, ski0_951, ski1_951, skk_1223, \
                         skk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1529[k] = f_19 * sik_935[k]
                    + f_1 * ski0_951[k]
                    - f_2 * ski1_951[k]
                    + f_3 * pc_z[k] * skk_1223[k];

        t_1530[k] = pb_y[k] * sil0_1215[k]
                    - f_14 * pc_y[k] * sil1_1215[k];

        t_1531[k] = f_15 * sik_972[k]
                    + f_3 * pc_y[k] * skk_1224[k];

        t_1532[k] = f_20 * sik_936[k]
                    + f_3 * pc_z[k] * skk_1224[k];
    }

#pragma omp simd aligned(t_1533, t_1534, t_1535, pb_y, pc_x, pc_y, sil0_1220, sik_974, \
                         sil1_1220, ski0_955, ski1_955, skk_1226, \
                         skk_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1533[k] = f_4 * ski0_955[k]
                    - f_5 * ski1_955[k]
                    + f_3 * pc_x[k] * skk_1227[k];

        t_1534[k] = f_15 * sik_974[k]
                    + f_3 * pc_y[k] * skk_1226[k];

        t_1535[k] = pb_y[k] * sil0_1220[k]
                    - f_14 * pc_y[k] * sil1_1220[k];
    }

#pragma omp simd aligned(t_1536, t_1537, t_1538, pc_x, pc_y, pc_z, sik_939, sik_977, ski0_958, \
                         ski1_958, skk_1227, skk_1229, skk_1230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1536[k] = f_6 * ski0_958[k]
                    - f_7 * ski1_958[k]
                    + f_3 * pc_x[k] * skk_1230[k];

        t_1537[k] = f_20 * sik_939[k]
                    + f_3 * pc_z[k] * skk_1227[k];

        t_1538[k] = f_15 * sik_977[k]
                    + f_3 * pc_y[k] * skk_1229[k];
    }

#pragma omp simd aligned(t_1539, t_1540, t_1541, pb_y, pc_x, pc_y, pc_z, sil0_1224, sik_942, \
                         sil1_1224, ski0_962, ski1_962, skk_1230, \
                         skk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1539[k] = pb_y[k] * sil0_1224[k]
                    - f_14 * pc_y[k] * sil1_1224[k];

        t_1540[k] = f_8 * ski0_962[k]
                    - f_9 * ski1_962[k]
                    + f_3 * pc_x[k] * skk_1234[k];

        t_1541[k] = f_20 * sik_942[k]
                    + f_3 * pc_z[k] * skk_1230[k];
    }

#pragma omp simd aligned(t_1542, t_1543, t_1544, pb_y, pc_x, pc_y, sil0_1229, sik_981, \
                         sil1_1229, ski0_964, ski1_964, skk_1233, \
                         skk_1236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1542[k] = f_8 * ski0_964[k]
                    - f_9 * ski1_964[k]
                    + f_3 * pc_x[k] * skk_1236[k];

        t_1543[k] = f_15 * sik_981[k]
                    + f_3 * pc_y[k] * skk_1233[k];

        t_1544[k] = pb_y[k] * sil0_1229[k]
                    - f_14 * pc_y[k] * sil1_1229[k];
    }

#pragma omp simd aligned(t_1545, t_1546, t_1547, pc_x, pc_z, sik_946, ski0_967, ski0_969, \
                         ski1_967, ski1_969, skk_1234, skk_1239, \
                         skk_1241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1545[k] = f_10 * ski0_967[k]
                    - f_11 * ski1_967[k]
                    + f_3 * pc_x[k] * skk_1239[k];

        t_1546[k] = f_20 * sik_946[k]
                    + f_3 * pc_z[k] * skk_1234[k];

        t_1547[k] = f_10 * ski0_969[k]
                    - f_11 * ski1_969[k]
                    + f_3 * pc_x[k] * skk_1241[k];
    }

#pragma omp simd aligned(t_1548, t_1549, t_1550, pb_y, pc_x, pc_y, sil0_1235, sik_986, \
                         sil1_1235, ski0_970, ski1_970, skk_1238, \
                         skk_1242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1548[k] = f_10 * ski0_970[k]
                    - f_11 * ski1_970[k]
                    + f_3 * pc_x[k] * skk_1242[k];

        t_1549[k] = f_15 * sik_986[k]
                    + f_3 * pc_y[k] * skk_1238[k];

        t_1550[k] = pb_y[k] * sil0_1235[k]
                    - f_14 * pc_y[k] * sil1_1235[k];
    }

#pragma omp simd aligned(t_1551, t_1552, t_1553, pc_x, pc_z, sik_951, ski0_973, ski0_975, \
                         ski1_973, ski1_975, skk_1239, skk_1245, \
                         skk_1247 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1551[k] = f_12 * ski0_973[k]
                    - f_13 * ski1_973[k]
                    + f_3 * pc_x[k] * skk_1245[k];

        t_1552[k] = f_20 * sik_951[k]
                    + f_3 * pc_z[k] * skk_1239[k];

        t_1553[k] = f_12 * ski0_975[k]
                    - f_13 * ski1_975[k]
                    + f_3 * pc_x[k] * skk_1247[k];
    }

#pragma omp simd aligned(t_1554, t_1555, t_1556, pc_x, pc_y, sik_992, ski0_976, ski0_977, \
                         ski1_976, ski1_977, skk_1244, skk_1248, \
                         skk_1249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1554[k] = f_12 * ski0_976[k]
                    - f_13 * ski1_976[k]
                    + f_3 * pc_x[k] * skk_1248[k];

        t_1555[k] = f_12 * ski0_977[k]
                    - f_13 * ski1_977[k]
                    + f_3 * pc_x[k] * skk_1249[k];

        t_1556[k] = f_15 * sik_992[k]
                    + f_3 * pc_y[k] * skk_1244[k];
    }

#pragma omp simd aligned(t_1557, t_1558, t_1559, t_1560, t_1561, t_1562, pb_y, pc_x, pc_y, \
                         sil0_1242, sil1_1242, skk_1252, skk_1253, skk_1254, skk_1255, \
                         skk_1256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1557[k] = pb_y[k] * sil0_1242[k]
                    - f_14 * pc_y[k] * sil1_1242[k];

        t_1558[k] = f_3 * pc_x[k] * skk_1252[k];

        t_1559[k] = f_3 * pc_x[k] * skk_1253[k];

        t_1560[k] = f_3 * pc_x[k] * skk_1254[k];

        t_1561[k] = f_3 * pc_x[k] * skk_1255[k];

        t_1562[k] = f_3 * pc_x[k] * skk_1256[k];
    }

#pragma omp simd aligned(t_1563, t_1564, t_1565, t_1566, pb_y, pc_x, pc_y, sil0_1251, \
                         sik_1000, sil1_1251, skk_1257, skk_1258, \
                         skk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1563[k] = f_3 * pc_x[k] * skk_1257[k];

        t_1564[k] = f_3 * pc_x[k] * skk_1258[k];

        t_1565[k] = f_3 * pc_x[k] * skk_1259[k];

        t_1566[k] = pb_y[k] * sil0_1251[k]
                    + f_21 * sik_1000[k]
                    - f_14 * pc_y[k] * sil1_1251[k];
    }

#pragma omp simd aligned(t_1567, t_1568, t_1569, pb_y, pc_y, pc_z, sil0_1253, sil0_1254, \
                         sik_964, sik_1002, sik_1003, sil1_1253, sil1_1254, \
                         skk_1252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1567[k] = f_20 * sik_964[k]
                    + f_3 * pc_z[k] * skk_1252[k];

        t_1568[k] = pb_y[k] * sil0_1253[k]
                    + f_20 * sik_1002[k]
                    - f_14 * pc_y[k] * sil1_1253[k];

        t_1569[k] = pb_y[k] * sil0_1254[k]
                    + f_19 * sik_1003[k]
                    - f_14 * pc_y[k] * sil1_1254[k];
    }

#pragma omp simd aligned(t_1570, t_1571, t_1572, pb_y, pc_y, sil0_1255, sil0_1256, sil0_1257, \
                         sik_1004, sik_1005, sik_1006, sil1_1255, sil1_1256, \
                         sil1_1257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1570[k] = pb_y[k] * sil0_1255[k]
                    + f_18 * sik_1004[k]
                    - f_14 * pc_y[k] * sil1_1255[k];

        t_1571[k] = pb_y[k] * sil0_1256[k]
                    + f_17 * sik_1005[k]
                    - f_14 * pc_y[k] * sil1_1256[k];

        t_1572[k] = pb_y[k] * sil0_1257[k]
                    + f_16 * sik_1006[k]
                    - f_14 * pc_y[k] * sil1_1257[k];
    }

#pragma omp simd aligned(t_1573, t_1574, t_1575, t_1576, pb_y, pc_x, pc_y, sil0_1259, \
                         sik_1007, sil1_1259, ski0_980, ski1_980, skk_1259, \
                         skk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1573[k] = f_15 * sik_1007[k]
                    + f_3 * pc_y[k] * skk_1259[k];

        t_1574[k] = pb_y[k] * sil0_1259[k]
                    - f_14 * pc_y[k] * sil1_1259[k];

        t_1575[k] = f_1 * ski0_980[k]
                    - f_2 * ski1_980[k]
                    + f_3 * pc_x[k] * skk_1260[k];

        t_1576[k] = f_3 * pc_y[k] * skk_1260[k];
    }

#pragma omp simd aligned(t_1577, t_1578, t_1579, t_1580, pc_x, pc_y, pc_z, sik_972, ski0_983, \
                         ski0_985, ski1_983, ski1_985, skk_1260, skk_1262, skk_1263, \
                         skk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1577[k] = f_0 * sik_972[k]
                    + f_3 * pc_z[k] * skk_1260[k];

        t_1578[k] = f_4 * ski0_983[k]
                    - f_5 * ski1_983[k]
                    + f_3 * pc_x[k] * skk_1263[k];

        t_1579[k] = f_3 * pc_y[k] * skk_1262[k];

        t_1580[k] = f_4 * ski0_985[k]
                    - f_5 * ski1_985[k]
                    + f_3 * pc_x[k] * skk_1265[k];
    }

#pragma omp simd aligned(t_1581, t_1582, t_1583, t_1584, pc_x, pc_y, pc_z, sik_975, ski0_986, \
                         ski0_989, ski1_986, ski1_989, skk_1263, skk_1265, skk_1266, \
                         skk_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1581[k] = f_6 * ski0_986[k]
                    - f_7 * ski1_986[k]
                    + f_3 * pc_x[k] * skk_1266[k];

        t_1582[k] = f_0 * sik_975[k]
                    + f_3 * pc_z[k] * skk_1263[k];

        t_1583[k] = f_3 * pc_y[k] * skk_1265[k];

        t_1584[k] = f_6 * ski0_989[k]
                    - f_7 * ski1_989[k]
                    + f_3 * pc_x[k] * skk_1269[k];
    }

#pragma omp simd aligned(t_1585, t_1586, t_1587, t_1588, pc_x, pc_y, pc_z, sik_978, ski0_990, \
                         ski0_992, ski1_990, ski1_992, skk_1266, skk_1269, skk_1270, \
                         skk_1272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1585[k] = f_8 * ski0_990[k]
                    - f_9 * ski1_990[k]
                    + f_3 * pc_x[k] * skk_1270[k];

        t_1586[k] = f_0 * sik_978[k]
                    + f_3 * pc_z[k] * skk_1266[k];

        t_1587[k] = f_8 * ski0_992[k]
                    - f_9 * ski1_992[k]
                    + f_3 * pc_x[k] * skk_1272[k];

        t_1588[k] = f_3 * pc_y[k] * skk_1269[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, pc_x, pc_z, sik_982, ski0_994, ski0_995, \
                         ski1_994, ski1_995, skk_1270, skk_1274, \
                         skk_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = f_8 * ski0_994[k]
                    - f_9 * ski1_994[k]
                    + f_3 * pc_x[k] * skk_1274[k];

        t_1590[k] = f_10 * ski0_995[k]
                    - f_11 * ski1_995[k]
                    + f_3 * pc_x[k] * skk_1275[k];

        t_1591[k] = f_0 * sik_982[k]
                    + f_3 * pc_z[k] * skk_1270[k];
    }

#pragma omp simd aligned(t_1592, t_1593, t_1594, t_1595, pc_x, pc_y, ski0_997, ski0_998, \
                         ski0_1000, ski1_997, ski1_998, ski1_1000, skk_1274, skk_1277, \
                         skk_1278, skk_1280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1592[k] = f_10 * ski0_997[k]
                    - f_11 * ski1_997[k]
                    + f_3 * pc_x[k] * skk_1277[k];

        t_1593[k] = f_10 * ski0_998[k]
                    - f_11 * ski1_998[k]
                    + f_3 * pc_x[k] * skk_1278[k];

        t_1594[k] = f_3 * pc_y[k] * skk_1274[k];

        t_1595[k] = f_10 * ski0_1000[k]
                    - f_11 * ski1_1000[k]
                    + f_3 * pc_x[k] * skk_1280[k];
    }

#pragma omp simd aligned(t_1596, t_1597, t_1598, pc_x, pc_z, sik_987, ski0_1001, ski0_1003, \
                         ski1_1001, ski1_1003, skk_1275, skk_1281, \
                         skk_1283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1596[k] = f_12 * ski0_1001[k]
                    - f_13 * ski1_1001[k]
                    + f_3 * pc_x[k] * skk_1281[k];

        t_1597[k] = f_0 * sik_987[k]
                    + f_3 * pc_z[k] * skk_1275[k];

        t_1598[k] = f_12 * ski0_1003[k]
                    - f_13 * ski1_1003[k]
                    + f_3 * pc_x[k] * skk_1283[k];
    }
}

static auto
compute_prim_skl_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t sik, const size_t ski0,
                                                           const size_t ski1, const size_t skk,
                                                           const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 2.5 / gamma;
    const auto f_5 = 2.5 * p / (gamma * q);
    const auto f_6 = 2.0 / gamma;
    const auto f_7 = 2.0 * p / (gamma * q);
    const auto f_8 = 1.5 / gamma;
    const auto f_9 = 1.5 * p / (gamma * q);
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / gamma;
    const auto f_13 = 0.5 * p / (gamma * q);

    auto *t_1599 = buffer.data(target + 1599);
    auto *t_1600 = buffer.data(target + 1600);
    auto *t_1601 = buffer.data(target + 1601);
    auto *t_1602 = buffer.data(target + 1602);
    auto *t_1603 = buffer.data(target + 1603);
    auto *t_1604 = buffer.data(target + 1604);
    auto *t_1605 = buffer.data(target + 1605);
    auto *t_1606 = buffer.data(target + 1606);
    auto *t_1607 = buffer.data(target + 1607);
    auto *t_1608 = buffer.data(target + 1608);
    auto *t_1609 = buffer.data(target + 1609);
    auto *t_1610 = buffer.data(target + 1610);
    auto *t_1611 = buffer.data(target + 1611);
    auto *t_1612 = buffer.data(target + 1612);
    auto *t_1613 = buffer.data(target + 1613);
    auto *t_1614 = buffer.data(target + 1614);
    auto *t_1615 = buffer.data(target + 1615);
    auto *t_1616 = buffer.data(target + 1616);
    auto *t_1617 = buffer.data(target + 1617);
    auto *t_1618 = buffer.data(target + 1618);
    auto *t_1619 = buffer.data(target + 1619);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sik_1000 = buffer.data(sik + 1000);
    const auto *sik_1007 = buffer.data(sik + 1007);

    const auto *ski0_1001 = buffer.data(ski0 + 1001);
    const auto *ski0_1003 = buffer.data(ski0 + 1003);
    const auto *ski0_1004 = buffer.data(ski0 + 1004);
    const auto *ski0_1005 = buffer.data(ski0 + 1005);
    const auto *ski0_1006 = buffer.data(ski0 + 1006);
    const auto *ski0_1007 = buffer.data(ski0 + 1007);

    const auto *ski1_1001 = buffer.data(ski1 + 1001);
    const auto *ski1_1003 = buffer.data(ski1 + 1003);
    const auto *ski1_1004 = buffer.data(ski1 + 1004);
    const auto *ski1_1005 = buffer.data(ski1 + 1005);
    const auto *ski1_1006 = buffer.data(ski1 + 1006);
    const auto *ski1_1007 = buffer.data(ski1 + 1007);

    const auto *skk_1280 = buffer.data(skk + 1280);
    const auto *skk_1284 = buffer.data(skk + 1284);
    const auto *skk_1285 = buffer.data(skk + 1285);
    const auto *skk_1287 = buffer.data(skk + 1287);
    const auto *skk_1288 = buffer.data(skk + 1288);
    const auto *skk_1289 = buffer.data(skk + 1289);
    const auto *skk_1290 = buffer.data(skk + 1290);
    const auto *skk_1291 = buffer.data(skk + 1291);
    const auto *skk_1292 = buffer.data(skk + 1292);
    const auto *skk_1293 = buffer.data(skk + 1293);
    const auto *skk_1294 = buffer.data(skk + 1294);
    const auto *skk_1295 = buffer.data(skk + 1295);

#pragma omp simd aligned(t_1599, t_1600, t_1601, t_1602, pc_x, pc_y, ski0_1004, ski0_1005, \
                         ski0_1007, ski1_1004, ski1_1005, ski1_1007, skk_1280, skk_1284, \
                         skk_1285, skk_1287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1599[k] = f_12 * ski0_1004[k]
                    - f_13 * ski1_1004[k]
                    + f_3 * pc_x[k] * skk_1284[k];

        t_1600[k] = f_12 * ski0_1005[k]
                    - f_13 * ski1_1005[k]
                    + f_3 * pc_x[k] * skk_1285[k];

        t_1601[k] = f_3 * pc_y[k] * skk_1280[k];

        t_1602[k] = f_12 * ski0_1007[k]
                    - f_13 * ski1_1007[k]
                    + f_3 * pc_x[k] * skk_1287[k];
    }

#pragma omp simd aligned(t_1603, t_1604, t_1605, t_1606, t_1607, t_1608, t_1609, pc_x, \
                         skk_1288, skk_1289, skk_1290, skk_1291, skk_1292, skk_1293, \
                         skk_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1603[k] = f_3 * pc_x[k] * skk_1288[k];

        t_1604[k] = f_3 * pc_x[k] * skk_1289[k];

        t_1605[k] = f_3 * pc_x[k] * skk_1290[k];

        t_1606[k] = f_3 * pc_x[k] * skk_1291[k];

        t_1607[k] = f_3 * pc_x[k] * skk_1292[k];

        t_1608[k] = f_3 * pc_x[k] * skk_1293[k];

        t_1609[k] = f_3 * pc_x[k] * skk_1294[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pc_x, pc_y, pc_z, sik_1000, \
                         ski0_1001, ski0_1003, ski1_1001, ski1_1003, skk_1288, skk_1290, \
                         skk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = f_3 * pc_x[k] * skk_1295[k];

        t_1611[k] = f_1 * ski0_1001[k]
                    - f_2 * ski1_1001[k]
                    + f_3 * pc_y[k] * skk_1288[k];

        t_1612[k] = f_0 * sik_1000[k]
                    + f_3 * pc_z[k] * skk_1288[k];

        t_1613[k] = f_4 * ski0_1003[k]
                    - f_5 * ski1_1003[k]
                    + f_3 * pc_y[k] * skk_1290[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, pc_y, ski0_1004, ski0_1005, ski0_1006, \
                         ski1_1004, ski1_1005, ski1_1006, skk_1291, skk_1292, \
                         skk_1293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = f_6 * ski0_1004[k]
                    - f_7 * ski1_1004[k]
                    + f_3 * pc_y[k] * skk_1291[k];

        t_1615[k] = f_8 * ski0_1005[k]
                    - f_9 * ski1_1005[k]
                    + f_3 * pc_y[k] * skk_1292[k];

        t_1616[k] = f_10 * ski0_1006[k]
                    - f_11 * ski1_1006[k]
                    + f_3 * pc_y[k] * skk_1293[k];
    }

#pragma omp simd aligned(t_1617, t_1618, t_1619, pc_y, pc_z, sik_1007, ski0_1007, ski1_1007, \
                         skk_1294, skk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1617[k] = f_12 * ski0_1007[k]
                    - f_13 * ski1_1007[k]
                    + f_3 * pc_y[k] * skk_1294[k];

        t_1618[k] = f_3 * pc_y[k] * skk_1295[k];

        t_1619[k] = f_0 * sik_1007[k]
                    + f_1 * ski0_1007[k]
                    - f_2 * ski1_1007[k]
                    + f_3 * pc_z[k] * skk_1295[k];
    }
}

auto
compute_prim_skl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sil0, const size_t sik,
                                                   const size_t sil1, const size_t ski0,
                                                   const size_t ski1, const size_t skk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skl_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, sil0, sik,
                                                              sil1, ski0, ski1, skk, ncols,
                                                              gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, sil0, sik,
                                                              sil1, skk, ncols, gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, sil0,
                                                               sik, sil1, skk, ncols, gamma, p,
                                                               q);

    compute_prim_skl_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, sil0,
                                                               sik, sil1, ski0, ski1, skk,
                                                               ncols, gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece12(buffer, target, pc, sik, ski0,
                                                               ski1, skk, ncols, gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, sil0,
                                                               sik, sil1, ski0, ski1, skk,
                                                               ncols, gamma, p, q);

    compute_prim_skl_three_center_electron_repulsion_0_piece14(buffer, target, pc, sik, ski0,
                                                               ski1, skk, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
