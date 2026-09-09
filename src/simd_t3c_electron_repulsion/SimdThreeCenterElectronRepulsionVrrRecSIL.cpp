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


#include "SimdThreeCenterElectronRepulsionVrrRecSIL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *shl0_0 = buffer.data(shl0 + 0);
    const auto *shl0_3 = buffer.data(shl0 + 3);
    const auto *shl0_5 = buffer.data(shl0 + 5);
    const auto *shl0_6 = buffer.data(shl0 + 6);
    const auto *shl0_9 = buffer.data(shl0 + 9);
    const auto *shl0_10 = buffer.data(shl0 + 10);
    const auto *shl0_12 = buffer.data(shl0 + 12);
    const auto *shl0_14 = buffer.data(shl0 + 14);
    const auto *shl0_15 = buffer.data(shl0 + 15);
    const auto *shl0_17 = buffer.data(shl0 + 17);
    const auto *shl0_18 = buffer.data(shl0 + 18);
    const auto *shl0_20 = buffer.data(shl0 + 20);
    const auto *shl0_21 = buffer.data(shl0 + 21);
    const auto *shl0_23 = buffer.data(shl0 + 23);
    const auto *shl0_24 = buffer.data(shl0 + 24);
    const auto *shl0_25 = buffer.data(shl0 + 25);
    const auto *shl0_27 = buffer.data(shl0 + 27);
    const auto *shl0_44 = buffer.data(shl0 + 44);

    const auto *shk_0 = buffer.data(shk + 0);
    const auto *shk_1 = buffer.data(shk + 1);
    const auto *shk_2 = buffer.data(shk + 2);
    const auto *shk_3 = buffer.data(shk + 3);
    const auto *shk_5 = buffer.data(shk + 5);
    const auto *shk_6 = buffer.data(shk + 6);
    const auto *shk_7 = buffer.data(shk + 7);
    const auto *shk_8 = buffer.data(shk + 8);
    const auto *shk_9 = buffer.data(shk + 9);
    const auto *shk_10 = buffer.data(shk + 10);
    const auto *shk_11 = buffer.data(shk + 11);
    const auto *shk_12 = buffer.data(shk + 12);
    const auto *shk_13 = buffer.data(shk + 13);
    const auto *shk_14 = buffer.data(shk + 14);
    const auto *shk_15 = buffer.data(shk + 15);
    const auto *shk_16 = buffer.data(shk + 16);
    const auto *shk_17 = buffer.data(shk + 17);
    const auto *shk_18 = buffer.data(shk + 18);
    const auto *shk_19 = buffer.data(shk + 19);
    const auto *shk_20 = buffer.data(shk + 20);
    const auto *shk_21 = buffer.data(shk + 21);
    const auto *shk_23 = buffer.data(shk + 23);
    const auto *shk_24 = buffer.data(shk + 24);
    const auto *shk_25 = buffer.data(shk + 25);
    const auto *shk_27 = buffer.data(shk + 27);
    const auto *shk_28 = buffer.data(shk + 28);
    const auto *shk_29 = buffer.data(shk + 29);
    const auto *shk_30 = buffer.data(shk + 30);
    const auto *shk_31 = buffer.data(shk + 31);
    const auto *shk_32 = buffer.data(shk + 32);
    const auto *shk_33 = buffer.data(shk + 33);
    const auto *shk_34 = buffer.data(shk + 34);
    const auto *shk_35 = buffer.data(shk + 35);
    const auto *shk_64 = buffer.data(shk + 64);
    const auto *shk_65 = buffer.data(shk + 65);
    const auto *shk_66 = buffer.data(shk + 66);
    const auto *shk_67 = buffer.data(shk + 67);
    const auto *shk_68 = buffer.data(shk + 68);
    const auto *shk_69 = buffer.data(shk + 69);
    const auto *shk_70 = buffer.data(shk + 70);
    const auto *shk_71 = buffer.data(shk + 71);

    const auto *shl1_0 = buffer.data(shl1 + 0);
    const auto *shl1_3 = buffer.data(shl1 + 3);
    const auto *shl1_5 = buffer.data(shl1 + 5);
    const auto *shl1_6 = buffer.data(shl1 + 6);
    const auto *shl1_9 = buffer.data(shl1 + 9);
    const auto *shl1_10 = buffer.data(shl1 + 10);
    const auto *shl1_12 = buffer.data(shl1 + 12);
    const auto *shl1_14 = buffer.data(shl1 + 14);
    const auto *shl1_15 = buffer.data(shl1 + 15);
    const auto *shl1_17 = buffer.data(shl1 + 17);
    const auto *shl1_18 = buffer.data(shl1 + 18);
    const auto *shl1_20 = buffer.data(shl1 + 20);
    const auto *shl1_21 = buffer.data(shl1 + 21);
    const auto *shl1_23 = buffer.data(shl1 + 23);
    const auto *shl1_24 = buffer.data(shl1 + 24);
    const auto *shl1_25 = buffer.data(shl1 + 25);
    const auto *shl1_27 = buffer.data(shl1 + 27);
    const auto *shl1_44 = buffer.data(shl1 + 44);

    const auto *sii0_0 = buffer.data(sii0 + 0);
    const auto *sii0_3 = buffer.data(sii0 + 3);
    const auto *sii0_5 = buffer.data(sii0 + 5);
    const auto *sii0_6 = buffer.data(sii0 + 6);
    const auto *sii0_9 = buffer.data(sii0 + 9);
    const auto *sii0_10 = buffer.data(sii0 + 10);
    const auto *sii0_12 = buffer.data(sii0 + 12);
    const auto *sii0_14 = buffer.data(sii0 + 14);
    const auto *sii0_15 = buffer.data(sii0 + 15);
    const auto *sii0_17 = buffer.data(sii0 + 17);
    const auto *sii0_18 = buffer.data(sii0 + 18);
    const auto *sii0_20 = buffer.data(sii0 + 20);
    const auto *sii0_21 = buffer.data(sii0 + 21);
    const auto *sii0_23 = buffer.data(sii0 + 23);
    const auto *sii0_24 = buffer.data(sii0 + 24);
    const auto *sii0_25 = buffer.data(sii0 + 25);
    const auto *sii0_26 = buffer.data(sii0 + 26);
    const auto *sii0_27 = buffer.data(sii0 + 27);
    const auto *sii0_49 = buffer.data(sii0 + 49);
    const auto *sii0_51 = buffer.data(sii0 + 51);
    const auto *sii0_52 = buffer.data(sii0 + 52);
    const auto *sii0_53 = buffer.data(sii0 + 53);
    const auto *sii0_54 = buffer.data(sii0 + 54);
    const auto *sii0_55 = buffer.data(sii0 + 55);

    const auto *sii1_0 = buffer.data(sii1 + 0);
    const auto *sii1_3 = buffer.data(sii1 + 3);
    const auto *sii1_5 = buffer.data(sii1 + 5);
    const auto *sii1_6 = buffer.data(sii1 + 6);
    const auto *sii1_9 = buffer.data(sii1 + 9);
    const auto *sii1_10 = buffer.data(sii1 + 10);
    const auto *sii1_12 = buffer.data(sii1 + 12);
    const auto *sii1_14 = buffer.data(sii1 + 14);
    const auto *sii1_15 = buffer.data(sii1 + 15);
    const auto *sii1_17 = buffer.data(sii1 + 17);
    const auto *sii1_18 = buffer.data(sii1 + 18);
    const auto *sii1_20 = buffer.data(sii1 + 20);
    const auto *sii1_21 = buffer.data(sii1 + 21);
    const auto *sii1_23 = buffer.data(sii1 + 23);
    const auto *sii1_24 = buffer.data(sii1 + 24);
    const auto *sii1_25 = buffer.data(sii1 + 25);
    const auto *sii1_26 = buffer.data(sii1 + 26);
    const auto *sii1_27 = buffer.data(sii1 + 27);
    const auto *sii1_49 = buffer.data(sii1 + 49);
    const auto *sii1_51 = buffer.data(sii1 + 51);
    const auto *sii1_52 = buffer.data(sii1 + 52);
    const auto *sii1_53 = buffer.data(sii1 + 53);
    const auto *sii1_54 = buffer.data(sii1 + 54);
    const auto *sii1_55 = buffer.data(sii1 + 55);

    const auto *sik_0 = buffer.data(sik + 0);
    const auto *sik_2 = buffer.data(sik + 2);
    const auto *sik_3 = buffer.data(sik + 3);
    const auto *sik_5 = buffer.data(sik + 5);
    const auto *sik_6 = buffer.data(sik + 6);
    const auto *sik_9 = buffer.data(sik + 9);
    const auto *sik_10 = buffer.data(sik + 10);
    const auto *sik_12 = buffer.data(sik + 12);
    const auto *sik_14 = buffer.data(sik + 14);
    const auto *sik_15 = buffer.data(sik + 15);
    const auto *sik_17 = buffer.data(sik + 17);
    const auto *sik_18 = buffer.data(sik + 18);
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
    const auto *sik_65 = buffer.data(sik + 65);
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
    const auto *sik_78 = buffer.data(sik + 78);
    const auto *sik_81 = buffer.data(sik + 81);
    const auto *sik_82 = buffer.data(sik + 82);
    const auto *sik_86 = buffer.data(sik + 86);
    const auto *sik_87 = buffer.data(sik + 87);
    const auto *sik_92 = buffer.data(sik + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, shk_0, shk_3, sii0_0, sii0_3, \
                         sii1_0, sii1_3, sik_0, sik_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shk_0[k]
                 + f_1 * sii0_0[k]
                 - f_2 * sii1_0[k]
                 + f_3 * pc_x[k] * sik_0[k];

        t_1[k] = f_3 * pc_y[k] * sik_0[k];

        t_2[k] = f_3 * pc_z[k] * sik_0[k];

        t_3[k] = f_0 * shk_3[k]
                 + f_4 * sii0_3[k]
                 - f_5 * sii1_3[k]
                 + f_3 * pc_x[k] * sik_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, shk_5, shk_6, sii0_5, sii0_6, sii1_5, \
                         sii1_6, sik_2, sik_5, sik_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sik_2[k];

        t_5[k] = f_0 * shk_5[k]
                 + f_4 * sii0_5[k]
                 - f_5 * sii1_5[k]
                 + f_3 * pc_x[k] * sik_5[k];

        t_6[k] = f_0 * shk_6[k]
                 + f_6 * sii0_6[k]
                 - f_7 * sii1_6[k]
                 + f_3 * pc_x[k] * sik_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, shk_9, sii0_9, sii1_9, sik_3, sik_5, \
                         sik_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sik_3[k];

        t_8[k] = f_3 * pc_y[k] * sik_5[k];

        t_9[k] = f_0 * shk_9[k]
                 + f_6 * sii0_9[k]
                 - f_7 * sii1_9[k]
                 + f_3 * pc_x[k] * sik_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, shk_10, shk_12, sii0_10, sii0_12, \
                         sii1_10, sii1_12, sik_6, sik_10, sik_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * shk_10[k]
                  + f_8 * sii0_10[k]
                  - f_9 * sii1_10[k]
                  + f_3 * pc_x[k] * sik_10[k];

        t_11[k] = f_3 * pc_z[k] * sik_6[k];

        t_12[k] = f_0 * shk_12[k]
                  + f_8 * sii0_12[k]
                  - f_9 * sii1_12[k]
                  + f_3 * pc_x[k] * sik_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, shk_14, shk_15, sii0_14, sii0_15, \
                         sii1_14, sii1_15, sik_9, sik_14, sik_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * sik_9[k];

        t_14[k] = f_0 * shk_14[k]
                  + f_8 * sii0_14[k]
                  - f_9 * sii1_14[k]
                  + f_3 * pc_x[k] * sik_14[k];

        t_15[k] = f_0 * shk_15[k]
                  + f_10 * sii0_15[k]
                  - f_11 * sii1_15[k]
                  + f_3 * pc_x[k] * sik_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, shk_17, shk_18, sii0_17, sii0_18, \
                         sii1_17, sii1_18, sik_10, sik_17, sik_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * sik_10[k];

        t_17[k] = f_0 * shk_17[k]
                  + f_10 * sii0_17[k]
                  - f_11 * sii1_17[k]
                  + f_3 * pc_x[k] * sik_17[k];

        t_18[k] = f_0 * shk_18[k]
                  + f_10 * sii0_18[k]
                  - f_11 * sii1_18[k]
                  + f_3 * pc_x[k] * sik_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, shk_20, shk_21, sii0_20, sii0_21, \
                         sii1_20, sii1_21, sik_14, sik_20, sik_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sik_14[k];

        t_20[k] = f_0 * shk_20[k]
                  + f_10 * sii0_20[k]
                  - f_11 * sii1_20[k]
                  + f_3 * pc_x[k] * sik_20[k];

        t_21[k] = f_0 * shk_21[k]
                  + f_12 * sii0_21[k]
                  - f_13 * sii1_21[k]
                  + f_3 * pc_x[k] * sik_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, shk_23, shk_24, sii0_23, sii0_24, \
                         sii1_23, sii1_24, sik_15, sik_23, sik_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * sik_15[k];

        t_23[k] = f_0 * shk_23[k]
                  + f_12 * sii0_23[k]
                  - f_13 * sii1_23[k]
                  + f_3 * pc_x[k] * sik_23[k];

        t_24[k] = f_0 * shk_24[k]
                  + f_12 * sii0_24[k]
                  - f_13 * sii1_24[k]
                  + f_3 * pc_x[k] * sik_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, shk_25, shk_27, sii0_25, sii0_27, \
                         sii1_25, sii1_27, sik_20, sik_25, sik_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * shk_25[k]
                  + f_12 * sii0_25[k]
                  - f_13 * sii1_25[k]
                  + f_3 * pc_x[k] * sik_25[k];

        t_26[k] = f_3 * pc_y[k] * sik_20[k];

        t_27[k] = f_0 * shk_27[k]
                  + f_12 * sii0_27[k]
                  - f_13 * sii1_27[k]
                  + f_3 * pc_x[k] * sik_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, shk_28, shk_29, shk_30, shk_31, \
                         shk_32, sik_28, sik_29, sik_30, sik_31, \
                         sik_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * shk_28[k]
                  + f_3 * pc_x[k] * sik_28[k];

        t_29[k] = f_0 * shk_29[k]
                  + f_3 * pc_x[k] * sik_29[k];

        t_30[k] = f_0 * shk_30[k]
                  + f_3 * pc_x[k] * sik_30[k];

        t_31[k] = f_0 * shk_31[k]
                  + f_3 * pc_x[k] * sik_31[k];

        t_32[k] = f_0 * shk_32[k]
                  + f_3 * pc_x[k] * sik_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, shk_33, shk_34, shk_35, sii0_21, \
                         sii1_21, sik_28, sik_33, sik_34, sik_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * shk_33[k]
                  + f_3 * pc_x[k] * sik_33[k];

        t_34[k] = f_0 * shk_34[k]
                  + f_3 * pc_x[k] * sik_34[k];

        t_35[k] = f_0 * shk_35[k]
                  + f_3 * pc_x[k] * sik_35[k];

        t_36[k] = f_1 * sii0_21[k]
                  - f_2 * sii1_21[k]
                  + f_3 * pc_y[k] * sik_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, sii0_23, sii0_24, sii0_25, \
                         sii1_23, sii1_24, sii1_25, sik_28, sik_30, sik_31, \
                         sik_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * sik_28[k];

        t_38[k] = f_4 * sii0_23[k]
                  - f_5 * sii1_23[k]
                  + f_3 * pc_y[k] * sik_30[k];

        t_39[k] = f_6 * sii0_24[k]
                  - f_7 * sii1_24[k]
                  + f_3 * pc_y[k] * sik_31[k];

        t_40[k] = f_8 * sii0_25[k]
                  - f_9 * sii1_25[k]
                  + f_3 * pc_y[k] * sik_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sii0_26, sii0_27, sii1_26, \
                         sii1_27, sik_33, sik_34, sik_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * sii0_26[k]
                  - f_11 * sii1_26[k]
                  + f_3 * pc_y[k] * sik_33[k];

        t_42[k] = f_12 * sii0_27[k]
                  - f_13 * sii1_27[k]
                  + f_3 * pc_y[k] * sik_34[k];

        t_43[k] = f_3 * pc_y[k] * sik_35[k];

        t_44[k] = f_1 * sii0_27[k]
                  - f_2 * sii1_27[k]
                  + f_3 * pc_z[k] * sik_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, shl0_0, shl0_3, shk_0, \
                         shk_1, shl1_0, shl1_3, sik_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * shl0_0[k]
                  - f_14 * pc_y[k] * shl1_0[k];

        t_46[k] = f_15 * shk_0[k]
                  + f_3 * pc_y[k] * sik_36[k];

        t_47[k] = f_3 * pc_z[k] * sik_36[k];

        t_48[k] = pb_y[k] * shl0_3[k]
                  + f_16 * shk_1[k]
                  - f_14 * pc_y[k] * shl1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, shl0_5, shl0_6, shk_2, \
                         shk_3, shl1_5, shl1_6, sik_38, sik_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * shk_2[k]
                  + f_3 * pc_y[k] * sik_38[k];

        t_50[k] = pb_y[k] * shl0_5[k]
                  - f_14 * pc_y[k] * shl1_5[k];

        t_51[k] = pb_y[k] * shl0_6[k]
                  + f_17 * shk_3[k]
                  - f_14 * pc_y[k] * shl1_6[k];

        t_52[k] = f_3 * pc_z[k] * sik_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, shl0_9, shl0_10, shk_5, \
                         shk_6, shl1_9, shl1_10, sik_41, sik_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * shk_5[k]
                  + f_3 * pc_y[k] * sik_41[k];

        t_54[k] = pb_y[k] * shl0_9[k]
                  - f_14 * pc_y[k] * shl1_9[k];

        t_55[k] = pb_y[k] * shl0_10[k]
                  + f_18 * shk_6[k]
                  - f_14 * pc_y[k] * shl1_10[k];

        t_56[k] = f_3 * pc_z[k] * sik_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, shl0_12, shl0_14, shl0_15, shk_8, \
                         shk_9, shk_10, shl1_12, shl1_14, shl1_15, \
                         sik_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * shl0_12[k]
                  + f_16 * shk_8[k]
                  - f_14 * pc_y[k] * shl1_12[k];

        t_58[k] = f_15 * shk_9[k]
                  + f_3 * pc_y[k] * sik_45[k];

        t_59[k] = pb_y[k] * shl0_14[k]
                  - f_14 * pc_y[k] * shl1_14[k];

        t_60[k] = pb_y[k] * shl0_15[k]
                  + f_19 * shk_10[k]
                  - f_14 * pc_y[k] * shl1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, shl0_17, shl0_18, shk_12, \
                         shk_13, shk_14, shl1_17, shl1_18, sik_46, \
                         sik_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * sik_46[k];

        t_62[k] = pb_y[k] * shl0_17[k]
                  + f_17 * shk_12[k]
                  - f_14 * pc_y[k] * shl1_17[k];

        t_63[k] = pb_y[k] * shl0_18[k]
                  + f_16 * shk_13[k]
                  - f_14 * pc_y[k] * shl1_18[k];

        t_64[k] = f_15 * shk_14[k]
                  + f_3 * pc_y[k] * sik_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, shl0_20, shl0_21, shl0_23, \
                         shk_15, shk_17, shl1_20, shl1_21, shl1_23, \
                         sik_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * shl0_20[k]
                  - f_14 * pc_y[k] * shl1_20[k];

        t_66[k] = pb_y[k] * shl0_21[k]
                  + f_0 * shk_15[k]
                  - f_14 * pc_y[k] * shl1_21[k];

        t_67[k] = f_3 * pc_z[k] * sik_51[k];

        t_68[k] = pb_y[k] * shl0_23[k]
                  + f_18 * shk_17[k]
                  - f_14 * pc_y[k] * shl1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, shl0_24, shl0_25, shl0_27, \
                         shk_18, shk_19, shk_20, shl1_24, shl1_25, shl1_27, \
                         sik_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * shl0_24[k]
                  + f_17 * shk_18[k]
                  - f_14 * pc_y[k] * shl1_24[k];

        t_70[k] = pb_y[k] * shl0_25[k]
                  + f_16 * shk_19[k]
                  - f_14 * pc_y[k] * shl1_25[k];

        t_71[k] = f_15 * shk_20[k]
                  + f_3 * pc_y[k] * sik_56[k];

        t_72[k] = pb_y[k] * shl0_27[k]
                  - f_14 * pc_y[k] * shl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, shk_64, shk_65, shk_66, shk_67, \
                         shk_68, sik_64, sik_65, sik_66, sik_67, \
                         sik_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_19 * shk_64[k]
                  + f_3 * pc_x[k] * sik_64[k];

        t_74[k] = f_19 * shk_65[k]
                  + f_3 * pc_x[k] * sik_65[k];

        t_75[k] = f_19 * shk_66[k]
                  + f_3 * pc_x[k] * sik_66[k];

        t_76[k] = f_19 * shk_67[k]
                  + f_3 * pc_x[k] * sik_67[k];

        t_77[k] = f_19 * shk_68[k]
                  + f_3 * pc_x[k] * sik_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, shk_28, shk_69, shk_70, shk_71, \
                         sii0_49, sii1_49, sik_64, sik_69, sik_70, \
                         sik_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_19 * shk_69[k]
                  + f_3 * pc_x[k] * sik_69[k];

        t_79[k] = f_19 * shk_70[k]
                  + f_3 * pc_x[k] * sik_70[k];

        t_80[k] = f_19 * shk_71[k]
                  + f_3 * pc_x[k] * sik_71[k];

        t_81[k] = f_15 * shk_28[k]
                  + f_1 * sii0_49[k]
                  - f_2 * sii1_49[k]
                  + f_3 * pc_y[k] * sik_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, shk_30, shk_31, sii0_51, sii0_52, \
                         sii1_51, sii1_52, sik_64, sik_66, sik_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * sik_64[k];

        t_83[k] = f_15 * shk_30[k]
                  + f_4 * sii0_51[k]
                  - f_5 * sii1_51[k]
                  + f_3 * pc_y[k] * sik_66[k];

        t_84[k] = f_15 * shk_31[k]
                  + f_6 * sii0_52[k]
                  - f_7 * sii1_52[k]
                  + f_3 * pc_y[k] * sik_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, shk_32, shk_33, shk_34, sii0_53, sii0_54, \
                         sii0_55, sii1_53, sii1_54, sii1_55, sik_68, sik_69, \
                         sik_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * shk_32[k]
                  + f_8 * sii0_53[k]
                  - f_9 * sii1_53[k]
                  + f_3 * pc_y[k] * sik_68[k];

        t_86[k] = f_15 * shk_33[k]
                  + f_10 * sii0_54[k]
                  - f_11 * sii1_54[k]
                  + f_3 * pc_y[k] * sik_69[k];

        t_87[k] = f_15 * shk_34[k]
                  + f_12 * sii0_55[k]
                  - f_13 * sii1_55[k]
                  + f_3 * pc_y[k] * sik_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, shl0_0, shl0_44, \
                         shk_35, shl1_0, shl1_44, sik_71, sik_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * shk_35[k]
                  + f_3 * pc_y[k] * sik_71[k];

        t_89[k] = pb_y[k] * shl0_44[k]
                  - f_14 * pc_y[k] * shl1_44[k];

        t_90[k] = pb_z[k] * shl0_0[k]
                  - f_14 * pc_z[k] * shl1_0[k];

        t_91[k] = f_3 * pc_y[k] * sik_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, shl0_3, shl0_5, shk_0, \
                         shk_2, shl1_3, shl1_5, sik_72, sik_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * shk_0[k]
                  + f_3 * pc_z[k] * sik_72[k];

        t_93[k] = pb_z[k] * shl0_3[k]
                  - f_14 * pc_z[k] * shl1_3[k];

        t_94[k] = f_3 * pc_y[k] * sik_74[k];

        t_95[k] = pb_z[k] * shl0_5[k]
                  + f_16 * shk_2[k]
                  - f_14 * pc_z[k] * shl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, shl0_6, shl0_9, shk_3, \
                         shk_5, shl1_6, shl1_9, sik_75, sik_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * shl0_6[k]
                  - f_14 * pc_z[k] * shl1_6[k];

        t_97[k] = f_15 * shk_3[k]
                  + f_3 * pc_z[k] * sik_75[k];

        t_98[k] = f_3 * pc_y[k] * sik_77[k];

        t_99[k] = pb_z[k] * shl0_9[k]
                  + f_17 * shk_5[k]
                  - f_14 * pc_z[k] * shl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, shl0_10, shl0_12, \
                         shk_6, shk_7, shl1_10, shl1_12, sik_78, \
                         sik_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * shl0_10[k]
                   - f_14 * pc_z[k] * shl1_10[k];

        t_101[k] = f_15 * shk_6[k]
                   + f_3 * pc_z[k] * sik_78[k];

        t_102[k] = pb_z[k] * shl0_12[k]
                   + f_16 * shk_7[k]
                   - f_14 * pc_z[k] * shl1_12[k];

        t_103[k] = f_3 * pc_y[k] * sik_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, shl0_14, shl0_15, shl0_17, \
                         shk_9, shk_10, shk_11, shl1_14, shl1_15, shl1_17, \
                         sik_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * shl0_14[k]
                   + f_18 * shk_9[k]
                   - f_14 * pc_z[k] * shl1_14[k];

        t_105[k] = pb_z[k] * shl0_15[k]
                   - f_14 * pc_z[k] * shl1_15[k];

        t_106[k] = f_15 * shk_10[k]
                   + f_3 * pc_z[k] * sik_82[k];

        t_107[k] = pb_z[k] * shl0_17[k]
                   + f_16 * shk_11[k]
                   - f_14 * pc_z[k] * shl1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, shl0_18, shl0_20, \
                         shl0_21, shk_12, shk_14, shl1_18, shl1_20, shl1_21, \
                         sik_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * shl0_18[k]
                   + f_17 * shk_12[k]
                   - f_14 * pc_z[k] * shl1_18[k];

        t_109[k] = f_3 * pc_y[k] * sik_86[k];

        t_110[k] = pb_z[k] * shl0_20[k]
                   + f_19 * shk_14[k]
                   - f_14 * pc_z[k] * shl1_20[k];

        t_111[k] = pb_z[k] * shl0_21[k]
                   - f_14 * pc_z[k] * shl1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, shl0_23, shl0_24, shk_15, shk_16, \
                         shk_17, shl1_23, shl1_24, sik_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * shk_15[k]
                   + f_3 * pc_z[k] * sik_87[k];

        t_113[k] = pb_z[k] * shl0_23[k]
                   + f_16 * shk_16[k]
                   - f_14 * pc_z[k] * shl1_23[k];

        t_114[k] = pb_z[k] * shl0_24[k]
                   + f_17 * shk_17[k]
                   - f_14 * pc_z[k] * shl1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, shl0_25, shl0_27, shk_18, \
                         shk_20, shl1_25, shl1_27, sik_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * shl0_25[k]
                   + f_18 * shk_18[k]
                   - f_14 * pc_z[k] * shl1_25[k];

        t_116[k] = f_3 * pc_y[k] * sik_92[k];

        t_117[k] = pb_z[k] * shl0_27[k]
                   + f_0 * shk_20[k]
                   - f_14 * pc_z[k] * shl1_27[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
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

    const auto *shl0_36 = buffer.data(shl0 + 36);
    const auto *shl0_48 = buffer.data(shl0 + 48);
    const auto *shl0_51 = buffer.data(shl0 + 51);
    const auto *shl0_55 = buffer.data(shl0 + 55);
    const auto *shl0_60 = buffer.data(shl0 + 60);
    const auto *shl0_66 = buffer.data(shl0 + 66);
    const auto *shl0_81 = buffer.data(shl0 + 81);
    const auto *shl0_90 = buffer.data(shl0 + 90);
    const auto *shl0_95 = buffer.data(shl0 + 95);
    const auto *shl0_99 = buffer.data(shl0 + 99);
    const auto *shl0_102 = buffer.data(shl0 + 102);
    const auto *shl0_104 = buffer.data(shl0 + 104);
    const auto *shl0_107 = buffer.data(shl0 + 107);
    const auto *shl0_108 = buffer.data(shl0 + 108);
    const auto *shl0_110 = buffer.data(shl0 + 110);
    const auto *shl0_113 = buffer.data(shl0 + 113);
    const auto *shl0_114 = buffer.data(shl0 + 114);
    const auto *shl0_115 = buffer.data(shl0 + 115);
    const auto *shl0_117 = buffer.data(shl0 + 117);
    const auto *shl0_134 = buffer.data(shl0 + 134);

    const auto *shk_28 = buffer.data(shk + 28);
    const auto *shk_35 = buffer.data(shk + 35);
    const auto *shk_36 = buffer.data(shk + 36);
    const auto *shk_38 = buffer.data(shk + 38);
    const auto *shk_39 = buffer.data(shk + 39);
    const auto *shk_41 = buffer.data(shk + 41);
    const auto *shk_42 = buffer.data(shk + 42);
    const auto *shk_45 = buffer.data(shk + 45);
    const auto *shk_46 = buffer.data(shk + 46);
    const auto *shk_50 = buffer.data(shk + 50);
    const auto *shk_51 = buffer.data(shk + 51);
    const auto *shk_56 = buffer.data(shk + 56);
    const auto *shk_64 = buffer.data(shk + 64);
    const auto *shk_66 = buffer.data(shk + 66);
    const auto *shk_67 = buffer.data(shk + 67);
    const auto *shk_68 = buffer.data(shk + 68);
    const auto *shk_69 = buffer.data(shk + 69);
    const auto *shk_70 = buffer.data(shk + 70);
    const auto *shk_71 = buffer.data(shk + 71);
    const auto *shk_72 = buffer.data(shk + 72);
    const auto *shk_74 = buffer.data(shk + 74);
    const auto *shk_75 = buffer.data(shk + 75);
    const auto *shk_77 = buffer.data(shk + 77);
    const auto *shk_80 = buffer.data(shk + 80);
    const auto *shk_81 = buffer.data(shk + 81);
    const auto *shk_84 = buffer.data(shk + 84);
    const auto *shk_85 = buffer.data(shk + 85);
    const auto *shk_86 = buffer.data(shk + 86);
    const auto *shk_89 = buffer.data(shk + 89);
    const auto *shk_90 = buffer.data(shk + 90);
    const auto *shk_91 = buffer.data(shk + 91);
    const auto *shk_92 = buffer.data(shk + 92);
    const auto *shk_100 = buffer.data(shk + 100);
    const auto *shk_101 = buffer.data(shk + 101);
    const auto *shk_102 = buffer.data(shk + 102);
    const auto *shk_103 = buffer.data(shk + 103);
    const auto *shk_104 = buffer.data(shk + 104);
    const auto *shk_105 = buffer.data(shk + 105);
    const auto *shk_106 = buffer.data(shk + 106);
    const auto *shk_107 = buffer.data(shk + 107);
    const auto *shk_108 = buffer.data(shk + 108);
    const auto *shk_111 = buffer.data(shk + 111);
    const auto *shk_113 = buffer.data(shk + 113);
    const auto *shk_114 = buffer.data(shk + 114);
    const auto *shk_117 = buffer.data(shk + 117);
    const auto *shk_118 = buffer.data(shk + 118);
    const auto *shk_120 = buffer.data(shk + 120);
    const auto *shk_122 = buffer.data(shk + 122);
    const auto *shk_123 = buffer.data(shk + 123);
    const auto *shk_125 = buffer.data(shk + 125);
    const auto *shk_126 = buffer.data(shk + 126);
    const auto *shk_128 = buffer.data(shk + 128);
    const auto *shk_129 = buffer.data(shk + 129);
    const auto *shk_131 = buffer.data(shk + 131);
    const auto *shk_132 = buffer.data(shk + 132);
    const auto *shk_133 = buffer.data(shk + 133);
    const auto *shk_135 = buffer.data(shk + 135);
    const auto *shk_136 = buffer.data(shk + 136);
    const auto *shk_137 = buffer.data(shk + 137);
    const auto *shk_138 = buffer.data(shk + 138);
    const auto *shk_139 = buffer.data(shk + 139);
    const auto *shk_140 = buffer.data(shk + 140);
    const auto *shk_141 = buffer.data(shk + 141);
    const auto *shk_142 = buffer.data(shk + 142);
    const auto *shk_143 = buffer.data(shk + 143);
    const auto *shk_172 = buffer.data(shk + 172);
    const auto *shk_173 = buffer.data(shk + 173);
    const auto *shk_174 = buffer.data(shk + 174);
    const auto *shk_175 = buffer.data(shk + 175);
    const auto *shk_176 = buffer.data(shk + 176);
    const auto *shk_177 = buffer.data(shk + 177);
    const auto *shk_178 = buffer.data(shk + 178);
    const auto *shk_179 = buffer.data(shk + 179);
    const auto *shk_180 = buffer.data(shk + 180);
    const auto *shk_183 = buffer.data(shk + 183);
    const auto *shk_185 = buffer.data(shk + 185);
    const auto *shk_186 = buffer.data(shk + 186);

    const auto *shl1_36 = buffer.data(shl1 + 36);
    const auto *shl1_48 = buffer.data(shl1 + 48);
    const auto *shl1_51 = buffer.data(shl1 + 51);
    const auto *shl1_55 = buffer.data(shl1 + 55);
    const auto *shl1_60 = buffer.data(shl1 + 60);
    const auto *shl1_66 = buffer.data(shl1 + 66);
    const auto *shl1_81 = buffer.data(shl1 + 81);
    const auto *shl1_90 = buffer.data(shl1 + 90);
    const auto *shl1_95 = buffer.data(shl1 + 95);
    const auto *shl1_99 = buffer.data(shl1 + 99);
    const auto *shl1_102 = buffer.data(shl1 + 102);
    const auto *shl1_104 = buffer.data(shl1 + 104);
    const auto *shl1_107 = buffer.data(shl1 + 107);
    const auto *shl1_108 = buffer.data(shl1 + 108);
    const auto *shl1_110 = buffer.data(shl1 + 110);
    const auto *shl1_113 = buffer.data(shl1 + 113);
    const auto *shl1_114 = buffer.data(shl1 + 114);
    const auto *shl1_115 = buffer.data(shl1 + 115);
    const auto *shl1_117 = buffer.data(shl1 + 117);
    const auto *shl1_134 = buffer.data(shl1 + 134);

    const auto *sii0_79 = buffer.data(sii0 + 79);
    const auto *sii0_80 = buffer.data(sii0 + 80);
    const auto *sii0_81 = buffer.data(sii0 + 81);
    const auto *sii0_82 = buffer.data(sii0 + 82);
    const auto *sii0_83 = buffer.data(sii0 + 83);
    const auto *sii0_84 = buffer.data(sii0 + 84);
    const auto *sii0_87 = buffer.data(sii0 + 87);
    const auto *sii0_89 = buffer.data(sii0 + 89);
    const auto *sii0_90 = buffer.data(sii0 + 90);
    const auto *sii0_93 = buffer.data(sii0 + 93);
    const auto *sii0_94 = buffer.data(sii0 + 94);
    const auto *sii0_96 = buffer.data(sii0 + 96);
    const auto *sii0_98 = buffer.data(sii0 + 98);
    const auto *sii0_99 = buffer.data(sii0 + 99);
    const auto *sii0_101 = buffer.data(sii0 + 101);
    const auto *sii0_102 = buffer.data(sii0 + 102);
    const auto *sii0_104 = buffer.data(sii0 + 104);
    const auto *sii0_105 = buffer.data(sii0 + 105);
    const auto *sii0_107 = buffer.data(sii0 + 107);
    const auto *sii0_108 = buffer.data(sii0 + 108);
    const auto *sii0_109 = buffer.data(sii0 + 109);
    const auto *sii0_110 = buffer.data(sii0 + 110);
    const auto *sii0_111 = buffer.data(sii0 + 111);
    const auto *sii0_135 = buffer.data(sii0 + 135);
    const auto *sii0_136 = buffer.data(sii0 + 136);
    const auto *sii0_137 = buffer.data(sii0 + 137);
    const auto *sii0_138 = buffer.data(sii0 + 138);
    const auto *sii0_139 = buffer.data(sii0 + 139);
    const auto *sii0_140 = buffer.data(sii0 + 140);
    const auto *sii0_143 = buffer.data(sii0 + 143);
    const auto *sii0_145 = buffer.data(sii0 + 145);
    const auto *sii0_146 = buffer.data(sii0 + 146);

    const auto *sii1_79 = buffer.data(sii1 + 79);
    const auto *sii1_80 = buffer.data(sii1 + 80);
    const auto *sii1_81 = buffer.data(sii1 + 81);
    const auto *sii1_82 = buffer.data(sii1 + 82);
    const auto *sii1_83 = buffer.data(sii1 + 83);
    const auto *sii1_84 = buffer.data(sii1 + 84);
    const auto *sii1_87 = buffer.data(sii1 + 87);
    const auto *sii1_89 = buffer.data(sii1 + 89);
    const auto *sii1_90 = buffer.data(sii1 + 90);
    const auto *sii1_93 = buffer.data(sii1 + 93);
    const auto *sii1_94 = buffer.data(sii1 + 94);
    const auto *sii1_96 = buffer.data(sii1 + 96);
    const auto *sii1_98 = buffer.data(sii1 + 98);
    const auto *sii1_99 = buffer.data(sii1 + 99);
    const auto *sii1_101 = buffer.data(sii1 + 101);
    const auto *sii1_102 = buffer.data(sii1 + 102);
    const auto *sii1_104 = buffer.data(sii1 + 104);
    const auto *sii1_105 = buffer.data(sii1 + 105);
    const auto *sii1_107 = buffer.data(sii1 + 107);
    const auto *sii1_108 = buffer.data(sii1 + 108);
    const auto *sii1_109 = buffer.data(sii1 + 109);
    const auto *sii1_110 = buffer.data(sii1 + 110);
    const auto *sii1_111 = buffer.data(sii1 + 111);
    const auto *sii1_135 = buffer.data(sii1 + 135);
    const auto *sii1_136 = buffer.data(sii1 + 136);
    const auto *sii1_137 = buffer.data(sii1 + 137);
    const auto *sii1_138 = buffer.data(sii1 + 138);
    const auto *sii1_139 = buffer.data(sii1 + 139);
    const auto *sii1_140 = buffer.data(sii1 + 140);
    const auto *sii1_143 = buffer.data(sii1 + 143);
    const auto *sii1_145 = buffer.data(sii1 + 145);
    const auto *sii1_146 = buffer.data(sii1 + 146);

    const auto *sik_100 = buffer.data(sik + 100);
    const auto *sik_101 = buffer.data(sik + 101);
    const auto *sik_102 = buffer.data(sik + 102);
    const auto *sik_103 = buffer.data(sik + 103);
    const auto *sik_104 = buffer.data(sik + 104);
    const auto *sik_105 = buffer.data(sik + 105);
    const auto *sik_106 = buffer.data(sik + 106);
    const auto *sik_107 = buffer.data(sik + 107);
    const auto *sik_108 = buffer.data(sik + 108);
    const auto *sik_110 = buffer.data(sik + 110);
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
    const auto *sik_144 = buffer.data(sik + 144);
    const auto *sik_146 = buffer.data(sik + 146);
    const auto *sik_147 = buffer.data(sik + 147);
    const auto *sik_149 = buffer.data(sik + 149);
    const auto *sik_150 = buffer.data(sik + 150);
    const auto *sik_153 = buffer.data(sik + 153);
    const auto *sik_154 = buffer.data(sik + 154);
    const auto *sik_158 = buffer.data(sik + 158);
    const auto *sik_159 = buffer.data(sik + 159);
    const auto *sik_164 = buffer.data(sik + 164);
    const auto *sik_172 = buffer.data(sik + 172);
    const auto *sik_173 = buffer.data(sik + 173);
    const auto *sik_174 = buffer.data(sik + 174);
    const auto *sik_175 = buffer.data(sik + 175);
    const auto *sik_176 = buffer.data(sik + 176);
    const auto *sik_177 = buffer.data(sik + 177);
    const auto *sik_178 = buffer.data(sik + 178);
    const auto *sik_179 = buffer.data(sik + 179);
    const auto *sik_180 = buffer.data(sik + 180);
    const auto *sik_182 = buffer.data(sik + 182);
    const auto *sik_183 = buffer.data(sik + 183);
    const auto *sik_185 = buffer.data(sik + 185);
    const auto *sik_186 = buffer.data(sik + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, shk_100, shk_101, shk_102, \
                         shk_103, shk_104, sik_100, sik_101, sik_102, sik_103, \
                         sik_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_19 * shk_100[k]
                   + f_3 * pc_x[k] * sik_100[k];

        t_119[k] = f_19 * shk_101[k]
                   + f_3 * pc_x[k] * sik_101[k];

        t_120[k] = f_19 * shk_102[k]
                   + f_3 * pc_x[k] * sik_102[k];

        t_121[k] = f_19 * shk_103[k]
                   + f_3 * pc_x[k] * sik_103[k];

        t_122[k] = f_19 * shk_104[k]
                   + f_3 * pc_x[k] * sik_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, shl0_36, shk_105, \
                         shk_106, shk_107, shl1_36, sik_105, sik_106, \
                         sik_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_19 * shk_105[k]
                   + f_3 * pc_x[k] * sik_105[k];

        t_124[k] = f_19 * shk_106[k]
                   + f_3 * pc_x[k] * sik_106[k];

        t_125[k] = f_19 * shk_107[k]
                   + f_3 * pc_x[k] * sik_107[k];

        t_126[k] = pb_z[k] * shl0_36[k]
                   - f_14 * pc_z[k] * shl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, shk_28, sii0_79, sii0_80, sii1_79, \
                         sii1_80, sik_100, sik_102, sik_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * shk_28[k]
                   + f_3 * pc_z[k] * sik_100[k];

        t_128[k] = f_4 * sii0_79[k]
                   - f_5 * sii1_79[k]
                   + f_3 * pc_y[k] * sik_102[k];

        t_129[k] = f_6 * sii0_80[k]
                   - f_7 * sii1_80[k]
                   + f_3 * pc_y[k] * sik_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, sii0_81, sii0_82, sii0_83, sii1_81, \
                         sii1_82, sii1_83, sik_104, sik_105, sik_106, \
                         sik_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * sii0_81[k]
                   - f_9 * sii1_81[k]
                   + f_3 * pc_y[k] * sik_104[k];

        t_131[k] = f_10 * sii0_82[k]
                   - f_11 * sii1_82[k]
                   + f_3 * pc_y[k] * sik_105[k];

        t_132[k] = f_12 * sii0_83[k]
                   - f_13 * sii1_83[k]
                   + f_3 * pc_y[k] * sik_106[k];

        t_133[k] = f_3 * pc_y[k] * sik_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, shk_35, shk_36, \
                         shk_108, sii0_83, sii0_84, sii1_83, sii1_84, sik_107, \
                         sik_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * shk_35[k]
                   + f_1 * sii0_83[k]
                   - f_2 * sii1_83[k]
                   + f_3 * pc_z[k] * sik_107[k];

        t_135[k] = f_18 * shk_108[k]
                   + f_1 * sii0_84[k]
                   - f_2 * sii1_84[k]
                   + f_3 * pc_x[k] * sik_108[k];

        t_136[k] = f_16 * shk_36[k]
                   + f_3 * pc_y[k] * sik_108[k];

        t_137[k] = f_3 * pc_z[k] * sik_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, shk_38, shk_111, shk_113, sii0_87, \
                         sii0_89, sii1_87, sii1_89, sik_110, sik_111, \
                         sik_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_18 * shk_111[k]
                   + f_4 * sii0_87[k]
                   - f_5 * sii1_87[k]
                   + f_3 * pc_x[k] * sik_111[k];

        t_139[k] = f_16 * shk_38[k]
                   + f_3 * pc_y[k] * sik_110[k];

        t_140[k] = f_18 * shk_113[k]
                   + f_4 * sii0_89[k]
                   - f_5 * sii1_89[k]
                   + f_3 * pc_x[k] * sik_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, shk_41, shk_114, sii0_90, \
                         sii1_90, sik_111, sik_113, sik_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_18 * shk_114[k]
                   + f_6 * sii0_90[k]
                   - f_7 * sii1_90[k]
                   + f_3 * pc_x[k] * sik_114[k];

        t_142[k] = f_3 * pc_z[k] * sik_111[k];

        t_143[k] = f_16 * shk_41[k]
                   + f_3 * pc_y[k] * sik_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, shk_117, shk_118, sii0_93, sii0_94, \
                         sii1_93, sii1_94, sik_114, sik_117, sik_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_18 * shk_117[k]
                   + f_6 * sii0_93[k]
                   - f_7 * sii1_93[k]
                   + f_3 * pc_x[k] * sik_117[k];

        t_145[k] = f_18 * shk_118[k]
                   + f_8 * sii0_94[k]
                   - f_9 * sii1_94[k]
                   + f_3 * pc_x[k] * sik_118[k];

        t_146[k] = f_3 * pc_z[k] * sik_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, shk_45, shk_120, shk_122, sii0_96, \
                         sii0_98, sii1_96, sii1_98, sik_117, sik_120, \
                         sik_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_18 * shk_120[k]
                   + f_8 * sii0_96[k]
                   - f_9 * sii1_96[k]
                   + f_3 * pc_x[k] * sik_120[k];

        t_148[k] = f_16 * shk_45[k]
                   + f_3 * pc_y[k] * sik_117[k];

        t_149[k] = f_18 * shk_122[k]
                   + f_8 * sii0_98[k]
                   - f_9 * sii1_98[k]
                   + f_3 * pc_x[k] * sik_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, shk_123, shk_125, sii0_99, sii0_101, \
                         sii1_99, sii1_101, sik_118, sik_123, sik_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_18 * shk_123[k]
                   + f_10 * sii0_99[k]
                   - f_11 * sii1_99[k]
                   + f_3 * pc_x[k] * sik_123[k];

        t_151[k] = f_3 * pc_z[k] * sik_118[k];

        t_152[k] = f_18 * shk_125[k]
                   + f_10 * sii0_101[k]
                   - f_11 * sii1_101[k]
                   + f_3 * pc_x[k] * sik_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, shk_50, shk_126, shk_128, sii0_102, \
                         sii0_104, sii1_102, sii1_104, sik_122, sik_126, \
                         sik_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_18 * shk_126[k]
                   + f_10 * sii0_102[k]
                   - f_11 * sii1_102[k]
                   + f_3 * pc_x[k] * sik_126[k];

        t_154[k] = f_16 * shk_50[k]
                   + f_3 * pc_y[k] * sik_122[k];

        t_155[k] = f_18 * shk_128[k]
                   + f_10 * sii0_104[k]
                   - f_11 * sii1_104[k]
                   + f_3 * pc_x[k] * sik_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, shk_129, shk_131, sii0_105, \
                         sii0_107, sii1_105, sii1_107, sik_123, sik_129, \
                         sik_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_18 * shk_129[k]
                   + f_12 * sii0_105[k]
                   - f_13 * sii1_105[k]
                   + f_3 * pc_x[k] * sik_129[k];

        t_157[k] = f_3 * pc_z[k] * sik_123[k];

        t_158[k] = f_18 * shk_131[k]
                   + f_12 * sii0_107[k]
                   - f_13 * sii1_107[k]
                   + f_3 * pc_x[k] * sik_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, shk_56, shk_132, shk_133, sii0_108, \
                         sii0_109, sii1_108, sii1_109, sik_128, sik_132, \
                         sik_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_18 * shk_132[k]
                   + f_12 * sii0_108[k]
                   - f_13 * sii1_108[k]
                   + f_3 * pc_x[k] * sik_132[k];

        t_160[k] = f_18 * shk_133[k]
                   + f_12 * sii0_109[k]
                   - f_13 * sii1_109[k]
                   + f_3 * pc_x[k] * sik_133[k];

        t_161[k] = f_16 * shk_56[k]
                   + f_3 * pc_y[k] * sik_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, shk_135, shk_136, shk_137, shk_138, \
                         sii0_111, sii1_111, sik_135, sik_136, sik_137, \
                         sik_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_18 * shk_135[k]
                   + f_12 * sii0_111[k]
                   - f_13 * sii1_111[k]
                   + f_3 * pc_x[k] * sik_135[k];

        t_163[k] = f_18 * shk_136[k]
                   + f_3 * pc_x[k] * sik_136[k];

        t_164[k] = f_18 * shk_137[k]
                   + f_3 * pc_x[k] * sik_137[k];

        t_165[k] = f_18 * shk_138[k]
                   + f_3 * pc_x[k] * sik_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, shk_139, shk_140, shk_141, \
                         shk_142, shk_143, sik_139, sik_140, sik_141, sik_142, \
                         sik_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_18 * shk_139[k]
                   + f_3 * pc_x[k] * sik_139[k];

        t_167[k] = f_18 * shk_140[k]
                   + f_3 * pc_x[k] * sik_140[k];

        t_168[k] = f_18 * shk_141[k]
                   + f_3 * pc_x[k] * sik_141[k];

        t_169[k] = f_18 * shk_142[k]
                   + f_3 * pc_x[k] * sik_142[k];

        t_170[k] = f_18 * shk_143[k]
                   + f_3 * pc_x[k] * sik_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, shk_64, shk_66, sii0_105, sii0_107, \
                         sii1_105, sii1_107, sik_136, sik_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * shk_64[k]
                   + f_1 * sii0_105[k]
                   - f_2 * sii1_105[k]
                   + f_3 * pc_y[k] * sik_136[k];

        t_172[k] = f_3 * pc_z[k] * sik_136[k];

        t_173[k] = f_16 * shk_66[k]
                   + f_4 * sii0_107[k]
                   - f_5 * sii1_107[k]
                   + f_3 * pc_y[k] * sik_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, shk_67, shk_68, shk_69, sii0_108, \
                         sii0_109, sii0_110, sii1_108, sii1_109, sii1_110, sik_139, sik_140, \
                         sik_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * shk_67[k]
                   + f_6 * sii0_108[k]
                   - f_7 * sii1_108[k]
                   + f_3 * pc_y[k] * sik_139[k];

        t_175[k] = f_16 * shk_68[k]
                   + f_8 * sii0_109[k]
                   - f_9 * sii1_109[k]
                   + f_3 * pc_y[k] * sik_140[k];

        t_176[k] = f_16 * shk_69[k]
                   + f_10 * sii0_110[k]
                   - f_11 * sii1_110[k]
                   + f_3 * pc_y[k] * sik_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, shl0_90, shk_70, \
                         shk_71, shl1_90, sii0_111, sii1_111, sik_142, \
                         sik_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * shk_70[k]
                   + f_12 * sii0_111[k]
                   - f_13 * sii1_111[k]
                   + f_3 * pc_y[k] * sik_142[k];

        t_178[k] = f_16 * shk_71[k]
                   + f_3 * pc_y[k] * sik_143[k];

        t_179[k] = f_1 * sii0_111[k]
                   - f_2 * sii1_111[k]
                   + f_3 * pc_z[k] * sik_143[k];

        t_180[k] = pb_y[k] * shl0_90[k]
                   - f_14 * pc_y[k] * shl1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, shl0_48, shk_36, \
                         shk_72, shk_74, shl1_48, sik_144, sik_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * shk_72[k]
                   + f_3 * pc_y[k] * sik_144[k];

        t_182[k] = f_15 * shk_36[k]
                   + f_3 * pc_z[k] * sik_144[k];

        t_183[k] = pb_z[k] * shl0_48[k]
                   - f_14 * pc_z[k] * shl1_48[k];

        t_184[k] = f_15 * shk_74[k]
                   + f_3 * pc_y[k] * sik_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, shl0_51, shl0_95, \
                         shk_39, shk_77, shl1_51, shl1_95, sik_147, \
                         sik_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * shl0_95[k]
                   - f_14 * pc_y[k] * shl1_95[k];

        t_186[k] = pb_z[k] * shl0_51[k]
                   - f_14 * pc_z[k] * shl1_51[k];

        t_187[k] = f_15 * shk_39[k]
                   + f_3 * pc_z[k] * sik_147[k];

        t_188[k] = f_15 * shk_77[k]
                   + f_3 * pc_y[k] * sik_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, shl0_55, shl0_99, \
                         shk_42, shl1_55, shl1_99, sik_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * shl0_99[k]
                   - f_14 * pc_y[k] * shl1_99[k];

        t_190[k] = pb_z[k] * shl0_55[k]
                   - f_14 * pc_z[k] * shl1_55[k];

        t_191[k] = f_15 * shk_42[k]
                   + f_3 * pc_z[k] * sik_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, shl0_102, shl0_104, shk_80, shk_81, \
                         shl1_102, shl1_104, sik_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * shl0_102[k]
                   + f_16 * shk_80[k]
                   - f_14 * pc_y[k] * shl1_102[k];

        t_193[k] = f_15 * shk_81[k]
                   + f_3 * pc_y[k] * sik_153[k];

        t_194[k] = pb_y[k] * shl0_104[k]
                   - f_14 * pc_y[k] * shl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, shl0_60, shl0_107, \
                         shk_46, shk_84, shl1_60, shl1_107, sik_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * shl0_60[k]
                   - f_14 * pc_z[k] * shl1_60[k];

        t_196[k] = f_15 * shk_46[k]
                   + f_3 * pc_z[k] * sik_154[k];

        t_197[k] = pb_y[k] * shl0_107[k]
                   + f_17 * shk_84[k]
                   - f_14 * pc_y[k] * shl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, shl0_108, shl0_110, shk_85, shk_86, \
                         shl1_108, shl1_110, sik_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * shl0_108[k]
                   + f_16 * shk_85[k]
                   - f_14 * pc_y[k] * shl1_108[k];

        t_199[k] = f_15 * shk_86[k]
                   + f_3 * pc_y[k] * sik_158[k];

        t_200[k] = pb_y[k] * shl0_110[k]
                   - f_14 * pc_y[k] * shl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, shl0_66, shl0_113, \
                         shk_51, shk_89, shl1_66, shl1_113, sik_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * shl0_66[k]
                   - f_14 * pc_z[k] * shl1_66[k];

        t_202[k] = f_15 * shk_51[k]
                   + f_3 * pc_z[k] * sik_159[k];

        t_203[k] = pb_y[k] * shl0_113[k]
                   + f_18 * shk_89[k]
                   - f_14 * pc_y[k] * shl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, shl0_114, shl0_115, shl0_117, \
                         shk_90, shk_91, shk_92, shl1_114, shl1_115, shl1_117, \
                         sik_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * shl0_114[k]
                   + f_17 * shk_90[k]
                   - f_14 * pc_y[k] * shl1_114[k];

        t_205[k] = pb_y[k] * shl0_115[k]
                   + f_16 * shk_91[k]
                   - f_14 * pc_y[k] * shl1_115[k];

        t_206[k] = f_15 * shk_92[k]
                   + f_3 * pc_y[k] * sik_164[k];

        t_207[k] = pb_y[k] * shl0_117[k]
                   - f_14 * pc_y[k] * shl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, shk_172, shk_173, shk_174, \
                         shk_175, shk_176, sik_172, sik_173, sik_174, sik_175, \
                         sik_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_18 * shk_172[k]
                   + f_3 * pc_x[k] * sik_172[k];

        t_209[k] = f_18 * shk_173[k]
                   + f_3 * pc_x[k] * sik_173[k];

        t_210[k] = f_18 * shk_174[k]
                   + f_3 * pc_x[k] * sik_174[k];

        t_211[k] = f_18 * shk_175[k]
                   + f_3 * pc_x[k] * sik_175[k];

        t_212[k] = f_18 * shk_176[k]
                   + f_3 * pc_x[k] * sik_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, shl0_81, shk_177, \
                         shk_178, shk_179, shl1_81, sik_177, sik_178, \
                         sik_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_18 * shk_177[k]
                   + f_3 * pc_x[k] * sik_177[k];

        t_214[k] = f_18 * shk_178[k]
                   + f_3 * pc_x[k] * sik_178[k];

        t_215[k] = f_18 * shk_179[k]
                   + f_3 * pc_x[k] * sik_179[k];

        t_216[k] = pb_z[k] * shl0_81[k]
                   - f_14 * pc_z[k] * shl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, shk_64, shk_102, shk_103, sii0_135, \
                         sii0_136, sii1_135, sii1_136, sik_172, sik_174, \
                         sik_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * shk_64[k]
                   + f_3 * pc_z[k] * sik_172[k];

        t_218[k] = f_15 * shk_102[k]
                   + f_4 * sii0_135[k]
                   - f_5 * sii1_135[k]
                   + f_3 * pc_y[k] * sik_174[k];

        t_219[k] = f_15 * shk_103[k]
                   + f_6 * sii0_136[k]
                   - f_7 * sii1_136[k]
                   + f_3 * pc_y[k] * sik_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, shk_104, shk_105, shk_106, sii0_137, \
                         sii0_138, sii0_139, sii1_137, sii1_138, sii1_139, sik_176, sik_177, \
                         sik_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * shk_104[k]
                   + f_8 * sii0_137[k]
                   - f_9 * sii1_137[k]
                   + f_3 * pc_y[k] * sik_176[k];

        t_221[k] = f_15 * shk_105[k]
                   + f_10 * sii0_138[k]
                   - f_11 * sii1_138[k]
                   + f_3 * pc_y[k] * sik_177[k];

        t_222[k] = f_15 * shk_106[k]
                   + f_12 * sii0_139[k]
                   - f_13 * sii1_139[k]
                   + f_3 * pc_y[k] * sik_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, shl0_134, shk_107, \
                         shk_180, shl1_134, sii0_140, sii1_140, sik_179, \
                         sik_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * shk_107[k]
                   + f_3 * pc_y[k] * sik_179[k];

        t_224[k] = pb_y[k] * shl0_134[k]
                   - f_14 * pc_y[k] * shl1_134[k];

        t_225[k] = f_18 * shk_180[k]
                   + f_1 * sii0_140[k]
                   - f_2 * sii1_140[k]
                   + f_3 * pc_x[k] * sik_180[k];

        t_226[k] = f_3 * pc_y[k] * sik_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, shk_72, shk_183, sii0_143, \
                         sii1_143, sik_180, sik_182, sik_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * shk_72[k]
                   + f_3 * pc_z[k] * sik_180[k];

        t_228[k] = f_18 * shk_183[k]
                   + f_4 * sii0_143[k]
                   - f_5 * sii1_143[k]
                   + f_3 * pc_x[k] * sik_183[k];

        t_229[k] = f_3 * pc_y[k] * sik_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, shk_75, shk_185, shk_186, sii0_145, \
                         sii0_146, sii1_145, sii1_146, sik_183, sik_185, \
                         sik_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_18 * shk_185[k]
                   + f_4 * sii0_145[k]
                   - f_5 * sii1_145[k]
                   + f_3 * pc_x[k] * sik_185[k];

        t_231[k] = f_18 * shk_186[k]
                   + f_6 * sii0_146[k]
                   - f_7 * sii1_146[k]
                   + f_3 * pc_x[k] * sik_186[k];

        t_232[k] = f_16 * shk_75[k]
                   + f_3 * pc_z[k] * sik_183[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
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

    const auto *shl0_135 = buffer.data(shl0 + 135);
    const auto *shl0_138 = buffer.data(shl0 + 138);
    const auto *shl0_141 = buffer.data(shl0 + 141);
    const auto *shl0_145 = buffer.data(shl0 + 145);
    const auto *shl0_147 = buffer.data(shl0 + 147);
    const auto *shl0_150 = buffer.data(shl0 + 150);
    const auto *shl0_152 = buffer.data(shl0 + 152);
    const auto *shl0_153 = buffer.data(shl0 + 153);
    const auto *shl0_156 = buffer.data(shl0 + 156);
    const auto *shl0_158 = buffer.data(shl0 + 158);
    const auto *shl0_159 = buffer.data(shl0 + 159);
    const auto *shl0_160 = buffer.data(shl0 + 160);

    const auto *shk_78 = buffer.data(shk + 78);
    const auto *shk_82 = buffer.data(shk + 82);
    const auto *shk_87 = buffer.data(shk + 87);
    const auto *shk_100 = buffer.data(shk + 100);
    const auto *shk_107 = buffer.data(shk + 107);
    const auto *shk_108 = buffer.data(shk + 108);
    const auto *shk_110 = buffer.data(shk + 110);
    const auto *shk_111 = buffer.data(shk + 111);
    const auto *shk_113 = buffer.data(shk + 113);
    const auto *shk_114 = buffer.data(shk + 114);
    const auto *shk_115 = buffer.data(shk + 115);
    const auto *shk_117 = buffer.data(shk + 117);
    const auto *shk_118 = buffer.data(shk + 118);
    const auto *shk_119 = buffer.data(shk + 119);
    const auto *shk_120 = buffer.data(shk + 120);
    const auto *shk_122 = buffer.data(shk + 122);
    const auto *shk_123 = buffer.data(shk + 123);
    const auto *shk_124 = buffer.data(shk + 124);
    const auto *shk_125 = buffer.data(shk + 125);
    const auto *shk_126 = buffer.data(shk + 126);
    const auto *shk_128 = buffer.data(shk + 128);
    const auto *shk_136 = buffer.data(shk + 136);
    const auto *shk_138 = buffer.data(shk + 138);
    const auto *shk_139 = buffer.data(shk + 139);
    const auto *shk_140 = buffer.data(shk + 140);
    const auto *shk_141 = buffer.data(shk + 141);
    const auto *shk_142 = buffer.data(shk + 142);
    const auto *shk_143 = buffer.data(shk + 143);
    const auto *shk_144 = buffer.data(shk + 144);
    const auto *shk_146 = buffer.data(shk + 146);
    const auto *shk_149 = buffer.data(shk + 149);
    const auto *shk_153 = buffer.data(shk + 153);
    const auto *shk_158 = buffer.data(shk + 158);
    const auto *shk_189 = buffer.data(shk + 189);
    const auto *shk_190 = buffer.data(shk + 190);
    const auto *shk_192 = buffer.data(shk + 192);
    const auto *shk_194 = buffer.data(shk + 194);
    const auto *shk_195 = buffer.data(shk + 195);
    const auto *shk_197 = buffer.data(shk + 197);
    const auto *shk_198 = buffer.data(shk + 198);
    const auto *shk_200 = buffer.data(shk + 200);
    const auto *shk_201 = buffer.data(shk + 201);
    const auto *shk_203 = buffer.data(shk + 203);
    const auto *shk_204 = buffer.data(shk + 204);
    const auto *shk_205 = buffer.data(shk + 205);
    const auto *shk_207 = buffer.data(shk + 207);
    const auto *shk_208 = buffer.data(shk + 208);
    const auto *shk_209 = buffer.data(shk + 209);
    const auto *shk_210 = buffer.data(shk + 210);
    const auto *shk_211 = buffer.data(shk + 211);
    const auto *shk_212 = buffer.data(shk + 212);
    const auto *shk_213 = buffer.data(shk + 213);
    const auto *shk_214 = buffer.data(shk + 214);
    const auto *shk_215 = buffer.data(shk + 215);
    const auto *shk_216 = buffer.data(shk + 216);
    const auto *shk_219 = buffer.data(shk + 219);
    const auto *shk_221 = buffer.data(shk + 221);
    const auto *shk_222 = buffer.data(shk + 222);
    const auto *shk_225 = buffer.data(shk + 225);
    const auto *shk_226 = buffer.data(shk + 226);
    const auto *shk_228 = buffer.data(shk + 228);
    const auto *shk_230 = buffer.data(shk + 230);
    const auto *shk_231 = buffer.data(shk + 231);
    const auto *shk_233 = buffer.data(shk + 233);
    const auto *shk_234 = buffer.data(shk + 234);
    const auto *shk_236 = buffer.data(shk + 236);
    const auto *shk_237 = buffer.data(shk + 237);
    const auto *shk_239 = buffer.data(shk + 239);
    const auto *shk_240 = buffer.data(shk + 240);
    const auto *shk_241 = buffer.data(shk + 241);
    const auto *shk_243 = buffer.data(shk + 243);
    const auto *shk_244 = buffer.data(shk + 244);
    const auto *shk_245 = buffer.data(shk + 245);
    const auto *shk_246 = buffer.data(shk + 246);
    const auto *shk_247 = buffer.data(shk + 247);
    const auto *shk_248 = buffer.data(shk + 248);
    const auto *shk_249 = buffer.data(shk + 249);
    const auto *shk_250 = buffer.data(shk + 250);
    const auto *shk_251 = buffer.data(shk + 251);
    const auto *shk_257 = buffer.data(shk + 257);
    const auto *shk_261 = buffer.data(shk + 261);
    const auto *shk_266 = buffer.data(shk + 266);
    const auto *shk_272 = buffer.data(shk + 272);

    const auto *shl1_135 = buffer.data(shl1 + 135);
    const auto *shl1_138 = buffer.data(shl1 + 138);
    const auto *shl1_141 = buffer.data(shl1 + 141);
    const auto *shl1_145 = buffer.data(shl1 + 145);
    const auto *shl1_147 = buffer.data(shl1 + 147);
    const auto *shl1_150 = buffer.data(shl1 + 150);
    const auto *shl1_152 = buffer.data(shl1 + 152);
    const auto *shl1_153 = buffer.data(shl1 + 153);
    const auto *shl1_156 = buffer.data(shl1 + 156);
    const auto *shl1_158 = buffer.data(shl1 + 158);
    const auto *shl1_159 = buffer.data(shl1 + 159);
    const auto *shl1_160 = buffer.data(shl1 + 160);

    const auto *sii0_149 = buffer.data(sii0 + 149);
    const auto *sii0_150 = buffer.data(sii0 + 150);
    const auto *sii0_152 = buffer.data(sii0 + 152);
    const auto *sii0_154 = buffer.data(sii0 + 154);
    const auto *sii0_155 = buffer.data(sii0 + 155);
    const auto *sii0_157 = buffer.data(sii0 + 157);
    const auto *sii0_158 = buffer.data(sii0 + 158);
    const auto *sii0_160 = buffer.data(sii0 + 160);
    const auto *sii0_161 = buffer.data(sii0 + 161);
    const auto *sii0_163 = buffer.data(sii0 + 163);
    const auto *sii0_164 = buffer.data(sii0 + 164);
    const auto *sii0_165 = buffer.data(sii0 + 165);
    const auto *sii0_166 = buffer.data(sii0 + 166);
    const auto *sii0_167 = buffer.data(sii0 + 167);
    const auto *sii0_168 = buffer.data(sii0 + 168);
    const auto *sii0_171 = buffer.data(sii0 + 171);
    const auto *sii0_173 = buffer.data(sii0 + 173);
    const auto *sii0_174 = buffer.data(sii0 + 174);
    const auto *sii0_177 = buffer.data(sii0 + 177);
    const auto *sii0_178 = buffer.data(sii0 + 178);
    const auto *sii0_180 = buffer.data(sii0 + 180);
    const auto *sii0_182 = buffer.data(sii0 + 182);
    const auto *sii0_183 = buffer.data(sii0 + 183);
    const auto *sii0_185 = buffer.data(sii0 + 185);
    const auto *sii0_186 = buffer.data(sii0 + 186);
    const auto *sii0_188 = buffer.data(sii0 + 188);
    const auto *sii0_189 = buffer.data(sii0 + 189);
    const auto *sii0_191 = buffer.data(sii0 + 191);
    const auto *sii0_192 = buffer.data(sii0 + 192);
    const auto *sii0_193 = buffer.data(sii0 + 193);
    const auto *sii0_194 = buffer.data(sii0 + 194);
    const auto *sii0_195 = buffer.data(sii0 + 195);
    const auto *sii0_201 = buffer.data(sii0 + 201);
    const auto *sii0_205 = buffer.data(sii0 + 205);
    const auto *sii0_210 = buffer.data(sii0 + 210);
    const auto *sii0_216 = buffer.data(sii0 + 216);

    const auto *sii1_149 = buffer.data(sii1 + 149);
    const auto *sii1_150 = buffer.data(sii1 + 150);
    const auto *sii1_152 = buffer.data(sii1 + 152);
    const auto *sii1_154 = buffer.data(sii1 + 154);
    const auto *sii1_155 = buffer.data(sii1 + 155);
    const auto *sii1_157 = buffer.data(sii1 + 157);
    const auto *sii1_158 = buffer.data(sii1 + 158);
    const auto *sii1_160 = buffer.data(sii1 + 160);
    const auto *sii1_161 = buffer.data(sii1 + 161);
    const auto *sii1_163 = buffer.data(sii1 + 163);
    const auto *sii1_164 = buffer.data(sii1 + 164);
    const auto *sii1_165 = buffer.data(sii1 + 165);
    const auto *sii1_166 = buffer.data(sii1 + 166);
    const auto *sii1_167 = buffer.data(sii1 + 167);
    const auto *sii1_168 = buffer.data(sii1 + 168);
    const auto *sii1_171 = buffer.data(sii1 + 171);
    const auto *sii1_173 = buffer.data(sii1 + 173);
    const auto *sii1_174 = buffer.data(sii1 + 174);
    const auto *sii1_177 = buffer.data(sii1 + 177);
    const auto *sii1_178 = buffer.data(sii1 + 178);
    const auto *sii1_180 = buffer.data(sii1 + 180);
    const auto *sii1_182 = buffer.data(sii1 + 182);
    const auto *sii1_183 = buffer.data(sii1 + 183);
    const auto *sii1_185 = buffer.data(sii1 + 185);
    const auto *sii1_186 = buffer.data(sii1 + 186);
    const auto *sii1_188 = buffer.data(sii1 + 188);
    const auto *sii1_189 = buffer.data(sii1 + 189);
    const auto *sii1_191 = buffer.data(sii1 + 191);
    const auto *sii1_192 = buffer.data(sii1 + 192);
    const auto *sii1_193 = buffer.data(sii1 + 193);
    const auto *sii1_194 = buffer.data(sii1 + 194);
    const auto *sii1_195 = buffer.data(sii1 + 195);
    const auto *sii1_201 = buffer.data(sii1 + 201);
    const auto *sii1_205 = buffer.data(sii1 + 205);
    const auto *sii1_210 = buffer.data(sii1 + 210);
    const auto *sii1_216 = buffer.data(sii1 + 216);

    const auto *sik_185 = buffer.data(sik + 185);
    const auto *sik_186 = buffer.data(sik + 186);
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
    const auto *sik_218 = buffer.data(sik + 218);
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

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, shk_189, shk_190, sii0_149, \
                         sii0_150, sii1_149, sii1_150, sik_185, sik_189, \
                         sik_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * sik_185[k];

        t_234[k] = f_18 * shk_189[k]
                   + f_6 * sii0_149[k]
                   - f_7 * sii1_149[k]
                   + f_3 * pc_x[k] * sik_189[k];

        t_235[k] = f_18 * shk_190[k]
                   + f_8 * sii0_150[k]
                   - f_9 * sii1_150[k]
                   + f_3 * pc_x[k] * sik_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, shk_78, shk_192, sii0_152, \
                         sii1_152, sik_186, sik_189, sik_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * shk_78[k]
                   + f_3 * pc_z[k] * sik_186[k];

        t_237[k] = f_18 * shk_192[k]
                   + f_8 * sii0_152[k]
                   - f_9 * sii1_152[k]
                   + f_3 * pc_x[k] * sik_192[k];

        t_238[k] = f_3 * pc_y[k] * sik_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, shk_82, shk_194, shk_195, sii0_154, \
                         sii0_155, sii1_154, sii1_155, sik_190, sik_194, \
                         sik_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_18 * shk_194[k]
                   + f_8 * sii0_154[k]
                   - f_9 * sii1_154[k]
                   + f_3 * pc_x[k] * sik_194[k];

        t_240[k] = f_18 * shk_195[k]
                   + f_10 * sii0_155[k]
                   - f_11 * sii1_155[k]
                   + f_3 * pc_x[k] * sik_195[k];

        t_241[k] = f_16 * shk_82[k]
                   + f_3 * pc_z[k] * sik_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, shk_197, shk_198, sii0_157, \
                         sii0_158, sii1_157, sii1_158, sik_194, sik_197, \
                         sik_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_18 * shk_197[k]
                   + f_10 * sii0_157[k]
                   - f_11 * sii1_157[k]
                   + f_3 * pc_x[k] * sik_197[k];

        t_243[k] = f_18 * shk_198[k]
                   + f_10 * sii0_158[k]
                   - f_11 * sii1_158[k]
                   + f_3 * pc_x[k] * sik_198[k];

        t_244[k] = f_3 * pc_y[k] * sik_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, shk_87, shk_200, shk_201, sii0_160, \
                         sii0_161, sii1_160, sii1_161, sik_195, sik_200, \
                         sik_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_18 * shk_200[k]
                   + f_10 * sii0_160[k]
                   - f_11 * sii1_160[k]
                   + f_3 * pc_x[k] * sik_200[k];

        t_246[k] = f_18 * shk_201[k]
                   + f_12 * sii0_161[k]
                   - f_13 * sii1_161[k]
                   + f_3 * pc_x[k] * sik_201[k];

        t_247[k] = f_16 * shk_87[k]
                   + f_3 * pc_z[k] * sik_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, shk_203, shk_204, shk_205, sii0_163, \
                         sii0_164, sii0_165, sii1_163, sii1_164, sii1_165, sik_203, sik_204, \
                         sik_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_18 * shk_203[k]
                   + f_12 * sii0_163[k]
                   - f_13 * sii1_163[k]
                   + f_3 * pc_x[k] * sik_203[k];

        t_249[k] = f_18 * shk_204[k]
                   + f_12 * sii0_164[k]
                   - f_13 * sii1_164[k]
                   + f_3 * pc_x[k] * sik_204[k];

        t_250[k] = f_18 * shk_205[k]
                   + f_12 * sii0_165[k]
                   - f_13 * sii1_165[k]
                   + f_3 * pc_x[k] * sik_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, shk_207, shk_208, shk_209, \
                         sii0_167, sii1_167, sik_200, sik_207, sik_208, \
                         sik_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * sik_200[k];

        t_252[k] = f_18 * shk_207[k]
                   + f_12 * sii0_167[k]
                   - f_13 * sii1_167[k]
                   + f_3 * pc_x[k] * sik_207[k];

        t_253[k] = f_18 * shk_208[k]
                   + f_3 * pc_x[k] * sik_208[k];

        t_254[k] = f_18 * shk_209[k]
                   + f_3 * pc_x[k] * sik_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, shk_210, shk_211, shk_212, \
                         shk_213, shk_214, sik_210, sik_211, sik_212, sik_213, \
                         sik_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_18 * shk_210[k]
                   + f_3 * pc_x[k] * sik_210[k];

        t_256[k] = f_18 * shk_211[k]
                   + f_3 * pc_x[k] * sik_211[k];

        t_257[k] = f_18 * shk_212[k]
                   + f_3 * pc_x[k] * sik_212[k];

        t_258[k] = f_18 * shk_213[k]
                   + f_3 * pc_x[k] * sik_213[k];

        t_259[k] = f_18 * shk_214[k]
                   + f_3 * pc_x[k] * sik_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, shk_100, shk_215, \
                         sii0_161, sii0_163, sii1_161, sii1_163, sik_208, sik_210, \
                         sik_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_18 * shk_215[k]
                   + f_3 * pc_x[k] * sik_215[k];

        t_261[k] = f_1 * sii0_161[k]
                   - f_2 * sii1_161[k]
                   + f_3 * pc_y[k] * sik_208[k];

        t_262[k] = f_16 * shk_100[k]
                   + f_3 * pc_z[k] * sik_208[k];

        t_263[k] = f_4 * sii0_163[k]
                   - f_5 * sii1_163[k]
                   + f_3 * pc_y[k] * sik_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, sii0_164, sii0_165, sii0_166, sii1_164, \
                         sii1_165, sii1_166, sik_211, sik_212, \
                         sik_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * sii0_164[k]
                   - f_7 * sii1_164[k]
                   + f_3 * pc_y[k] * sik_211[k];

        t_265[k] = f_8 * sii0_165[k]
                   - f_9 * sii1_165[k]
                   + f_3 * pc_y[k] * sik_212[k];

        t_266[k] = f_10 * sii0_166[k]
                   - f_11 * sii1_166[k]
                   + f_3 * pc_y[k] * sik_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, shk_107, shk_216, \
                         sii0_167, sii0_168, sii1_167, sii1_168, sik_214, sik_215, \
                         sik_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * sii0_167[k]
                   - f_13 * sii1_167[k]
                   + f_3 * pc_y[k] * sik_214[k];

        t_268[k] = f_3 * pc_y[k] * sik_215[k];

        t_269[k] = f_16 * shk_107[k]
                   + f_1 * sii0_167[k]
                   - f_2 * sii1_167[k]
                   + f_3 * pc_z[k] * sik_215[k];

        t_270[k] = f_17 * shk_216[k]
                   + f_1 * sii0_168[k]
                   - f_2 * sii1_168[k]
                   + f_3 * pc_x[k] * sik_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, shk_108, shk_110, \
                         shk_219, sii0_171, sii1_171, sik_216, sik_218, \
                         sik_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * shk_108[k]
                   + f_3 * pc_y[k] * sik_216[k];

        t_272[k] = f_3 * pc_z[k] * sik_216[k];

        t_273[k] = f_17 * shk_219[k]
                   + f_4 * sii0_171[k]
                   - f_5 * sii1_171[k]
                   + f_3 * pc_x[k] * sik_219[k];

        t_274[k] = f_17 * shk_110[k]
                   + f_3 * pc_y[k] * sik_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, shk_221, shk_222, sii0_173, \
                         sii0_174, sii1_173, sii1_174, sik_219, sik_221, \
                         sik_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_17 * shk_221[k]
                   + f_4 * sii0_173[k]
                   - f_5 * sii1_173[k]
                   + f_3 * pc_x[k] * sik_221[k];

        t_276[k] = f_17 * shk_222[k]
                   + f_6 * sii0_174[k]
                   - f_7 * sii1_174[k]
                   + f_3 * pc_x[k] * sik_222[k];

        t_277[k] = f_3 * pc_z[k] * sik_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, shk_113, shk_225, shk_226, sii0_177, \
                         sii0_178, sii1_177, sii1_178, sik_221, sik_225, \
                         sik_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * shk_113[k]
                   + f_3 * pc_y[k] * sik_221[k];

        t_279[k] = f_17 * shk_225[k]
                   + f_6 * sii0_177[k]
                   - f_7 * sii1_177[k]
                   + f_3 * pc_x[k] * sik_225[k];

        t_280[k] = f_17 * shk_226[k]
                   + f_8 * sii0_178[k]
                   - f_9 * sii1_178[k]
                   + f_3 * pc_x[k] * sik_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, pc_z, shk_117, shk_228, sii0_180, \
                         sii1_180, sik_222, sik_225, sik_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * sik_222[k];

        t_282[k] = f_17 * shk_228[k]
                   + f_8 * sii0_180[k]
                   - f_9 * sii1_180[k]
                   + f_3 * pc_x[k] * sik_228[k];

        t_283[k] = f_17 * shk_117[k]
                   + f_3 * pc_y[k] * sik_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_z, shk_230, shk_231, sii0_182, \
                         sii0_183, sii1_182, sii1_183, sik_226, sik_230, \
                         sik_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_17 * shk_230[k]
                   + f_8 * sii0_182[k]
                   - f_9 * sii1_182[k]
                   + f_3 * pc_x[k] * sik_230[k];

        t_285[k] = f_17 * shk_231[k]
                   + f_10 * sii0_183[k]
                   - f_11 * sii1_183[k]
                   + f_3 * pc_x[k] * sik_231[k];

        t_286[k] = f_3 * pc_z[k] * sik_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, shk_122, shk_233, shk_234, sii0_185, \
                         sii0_186, sii1_185, sii1_186, sik_230, sik_233, \
                         sik_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_17 * shk_233[k]
                   + f_10 * sii0_185[k]
                   - f_11 * sii1_185[k]
                   + f_3 * pc_x[k] * sik_233[k];

        t_288[k] = f_17 * shk_234[k]
                   + f_10 * sii0_186[k]
                   - f_11 * sii1_186[k]
                   + f_3 * pc_x[k] * sik_234[k];

        t_289[k] = f_17 * shk_122[k]
                   + f_3 * pc_y[k] * sik_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, shk_236, shk_237, sii0_188, \
                         sii0_189, sii1_188, sii1_189, sik_231, sik_236, \
                         sik_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_17 * shk_236[k]
                   + f_10 * sii0_188[k]
                   - f_11 * sii1_188[k]
                   + f_3 * pc_x[k] * sik_236[k];

        t_291[k] = f_17 * shk_237[k]
                   + f_12 * sii0_189[k]
                   - f_13 * sii1_189[k]
                   + f_3 * pc_x[k] * sik_237[k];

        t_292[k] = f_3 * pc_z[k] * sik_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, shk_239, shk_240, shk_241, sii0_191, \
                         sii0_192, sii0_193, sii1_191, sii1_192, sii1_193, sik_239, sik_240, \
                         sik_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_17 * shk_239[k]
                   + f_12 * sii0_191[k]
                   - f_13 * sii1_191[k]
                   + f_3 * pc_x[k] * sik_239[k];

        t_294[k] = f_17 * shk_240[k]
                   + f_12 * sii0_192[k]
                   - f_13 * sii1_192[k]
                   + f_3 * pc_x[k] * sik_240[k];

        t_295[k] = f_17 * shk_241[k]
                   + f_12 * sii0_193[k]
                   - f_13 * sii1_193[k]
                   + f_3 * pc_x[k] * sik_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pc_x, pc_y, shk_128, shk_243, shk_244, \
                         shk_245, sii0_195, sii1_195, sik_236, sik_243, sik_244, \
                         sik_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * shk_128[k]
                   + f_3 * pc_y[k] * sik_236[k];

        t_297[k] = f_17 * shk_243[k]
                   + f_12 * sii0_195[k]
                   - f_13 * sii1_195[k]
                   + f_3 * pc_x[k] * sik_243[k];

        t_298[k] = f_17 * shk_244[k]
                   + f_3 * pc_x[k] * sik_244[k];

        t_299[k] = f_17 * shk_245[k]
                   + f_3 * pc_x[k] * sik_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, shk_246, shk_247, shk_248, \
                         shk_249, shk_250, sik_246, sik_247, sik_248, sik_249, \
                         sik_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_17 * shk_246[k]
                   + f_3 * pc_x[k] * sik_246[k];

        t_301[k] = f_17 * shk_247[k]
                   + f_3 * pc_x[k] * sik_247[k];

        t_302[k] = f_17 * shk_248[k]
                   + f_3 * pc_x[k] * sik_248[k];

        t_303[k] = f_17 * shk_249[k]
                   + f_3 * pc_x[k] * sik_249[k];

        t_304[k] = f_17 * shk_250[k]
                   + f_3 * pc_x[k] * sik_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_y, pc_z, shk_136, shk_251, sii0_189, \
                         sii1_189, sik_244, sik_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_17 * shk_251[k]
                   + f_3 * pc_x[k] * sik_251[k];

        t_306[k] = f_17 * shk_136[k]
                   + f_1 * sii0_189[k]
                   - f_2 * sii1_189[k]
                   + f_3 * pc_y[k] * sik_244[k];

        t_307[k] = f_3 * pc_z[k] * sik_244[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_y, shk_138, shk_139, shk_140, sii0_191, \
                         sii0_192, sii0_193, sii1_191, sii1_192, sii1_193, sik_246, sik_247, \
                         sik_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_17 * shk_138[k]
                   + f_4 * sii0_191[k]
                   - f_5 * sii1_191[k]
                   + f_3 * pc_y[k] * sik_246[k];

        t_309[k] = f_17 * shk_139[k]
                   + f_6 * sii0_192[k]
                   - f_7 * sii1_192[k]
                   + f_3 * pc_y[k] * sik_247[k];

        t_310[k] = f_17 * shk_140[k]
                   + f_8 * sii0_193[k]
                   - f_9 * sii1_193[k]
                   + f_3 * pc_y[k] * sik_248[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, shk_141, shk_142, shk_143, \
                         sii0_194, sii0_195, sii1_194, sii1_195, sik_249, sik_250, \
                         sik_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_17 * shk_141[k]
                   + f_10 * sii0_194[k]
                   - f_11 * sii1_194[k]
                   + f_3 * pc_y[k] * sik_249[k];

        t_312[k] = f_17 * shk_142[k]
                   + f_12 * sii0_195[k]
                   - f_13 * sii1_195[k]
                   + f_3 * pc_y[k] * sik_250[k];

        t_313[k] = f_17 * shk_143[k]
                   + f_3 * pc_y[k] * sik_251[k];

        t_314[k] = f_1 * sii0_195[k]
                   - f_2 * sii1_195[k]
                   + f_3 * pc_z[k] * sik_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pb_z, pc_y, pc_z, shl0_135, shl0_138, \
                         shk_108, shk_144, shl1_135, shl1_138, \
                         sik_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_z[k] * shl0_135[k]
                   - f_14 * pc_z[k] * shl1_135[k];

        t_316[k] = f_16 * shk_144[k]
                   + f_3 * pc_y[k] * sik_252[k];

        t_317[k] = f_15 * shk_108[k]
                   + f_3 * pc_z[k] * sik_252[k];

        t_318[k] = pb_z[k] * shl0_138[k]
                   - f_14 * pc_z[k] * shl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_z, pc_x, pc_y, pc_z, shl0_141, shk_146, \
                         shk_257, shl1_141, sii0_201, sii1_201, sik_254, \
                         sik_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * shk_146[k]
                   + f_3 * pc_y[k] * sik_254[k];

        t_320[k] = f_17 * shk_257[k]
                   + f_4 * sii0_201[k]
                   - f_5 * sii1_201[k]
                   + f_3 * pc_x[k] * sik_257[k];

        t_321[k] = pb_z[k] * shl0_141[k]
                   - f_14 * pc_z[k] * shl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, shk_111, shk_149, shk_261, \
                         sii0_205, sii1_205, sik_255, sik_257, \
                         sik_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * shk_111[k]
                   + f_3 * pc_z[k] * sik_255[k];

        t_323[k] = f_16 * shk_149[k]
                   + f_3 * pc_y[k] * sik_257[k];

        t_324[k] = f_17 * shk_261[k]
                   + f_6 * sii0_205[k]
                   - f_7 * sii1_205[k]
                   + f_3 * pc_x[k] * sik_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_y, pc_z, shl0_145, shl0_147, \
                         shk_114, shk_115, shk_153, shl1_145, shl1_147, sik_258, \
                         sik_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pb_z[k] * shl0_145[k]
                   - f_14 * pc_z[k] * shl1_145[k];

        t_326[k] = f_15 * shk_114[k]
                   + f_3 * pc_z[k] * sik_258[k];

        t_327[k] = pb_z[k] * shl0_147[k]
                   + f_16 * shk_115[k]
                   - f_14 * pc_z[k] * shl1_147[k];

        t_328[k] = f_16 * shk_153[k]
                   + f_3 * pc_y[k] * sik_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_z, pc_x, pc_z, shl0_150, shk_118, shk_266, \
                         shl1_150, sii0_210, sii1_210, sik_262, \
                         sik_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_17 * shk_266[k]
                   + f_8 * sii0_210[k]
                   - f_9 * sii1_210[k]
                   + f_3 * pc_x[k] * sik_266[k];

        t_330[k] = pb_z[k] * shl0_150[k]
                   - f_14 * pc_z[k] * shl1_150[k];

        t_331[k] = f_15 * shk_118[k]
                   + f_3 * pc_z[k] * sik_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_y, pc_z, shl0_152, shl0_153, shk_119, \
                         shk_120, shk_158, shl1_152, shl1_153, \
                         sik_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_z[k] * shl0_152[k]
                   + f_16 * shk_119[k]
                   - f_14 * pc_z[k] * shl1_152[k];

        t_333[k] = pb_z[k] * shl0_153[k]
                   + f_17 * shk_120[k]
                   - f_14 * pc_z[k] * shl1_153[k];

        t_334[k] = f_16 * shk_158[k]
                   + f_3 * pc_y[k] * sik_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_z, pc_x, pc_z, shl0_156, shk_123, shk_272, \
                         shl1_156, sii0_216, sii1_216, sik_267, \
                         sik_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_17 * shk_272[k]
                   + f_10 * sii0_216[k]
                   - f_11 * sii1_216[k]
                   + f_3 * pc_x[k] * sik_272[k];

        t_336[k] = pb_z[k] * shl0_156[k]
                   - f_14 * pc_z[k] * shl1_156[k];

        t_337[k] = f_15 * shk_123[k]
                   + f_3 * pc_z[k] * sik_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_z, shl0_158, shl0_159, shl0_160, \
                         shk_124, shk_125, shk_126, shl1_158, shl1_159, \
                         shl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_z[k] * shl0_158[k]
                   + f_16 * shk_124[k]
                   - f_14 * pc_z[k] * shl1_158[k];

        t_339[k] = pb_z[k] * shl0_159[k]
                   + f_17 * shk_125[k]
                   - f_14 * pc_z[k] * shl1_159[k];

        t_340[k] = pb_z[k] * shl0_160[k]
                   + f_18 * shk_126[k]
                   - f_14 * pc_z[k] * shl1_160[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *shl0_171 = buffer.data(shl0 + 171);
    const auto *shl0_225 = buffer.data(shl0 + 225);
    const auto *shl0_228 = buffer.data(shl0 + 228);
    const auto *shl0_230 = buffer.data(shl0 + 230);
    const auto *shl0_231 = buffer.data(shl0 + 231);
    const auto *shl0_234 = buffer.data(shl0 + 234);
    const auto *shl0_235 = buffer.data(shl0 + 235);
    const auto *shl0_237 = buffer.data(shl0 + 237);
    const auto *shl0_239 = buffer.data(shl0 + 239);
    const auto *shl0_240 = buffer.data(shl0 + 240);
    const auto *shl0_242 = buffer.data(shl0 + 242);
    const auto *shl0_243 = buffer.data(shl0 + 243);
    const auto *shl0_245 = buffer.data(shl0 + 245);
    const auto *shl0_246 = buffer.data(shl0 + 246);
    const auto *shl0_248 = buffer.data(shl0 + 248);
    const auto *shl0_249 = buffer.data(shl0 + 249);
    const auto *shl0_250 = buffer.data(shl0 + 250);
    const auto *shl0_252 = buffer.data(shl0 + 252);
    const auto *shl0_269 = buffer.data(shl0 + 269);

    const auto *shk_136 = buffer.data(shk + 136);
    const auto *shk_143 = buffer.data(shk + 143);
    const auto *shk_144 = buffer.data(shk + 144);
    const auto *shk_147 = buffer.data(shk + 147);
    const auto *shk_150 = buffer.data(shk + 150);
    const auto *shk_154 = buffer.data(shk + 154);
    const auto *shk_159 = buffer.data(shk + 159);
    const auto *shk_164 = buffer.data(shk + 164);
    const auto *shk_172 = buffer.data(shk + 172);
    const auto *shk_174 = buffer.data(shk + 174);
    const auto *shk_175 = buffer.data(shk + 175);
    const auto *shk_176 = buffer.data(shk + 176);
    const auto *shk_177 = buffer.data(shk + 177);
    const auto *shk_178 = buffer.data(shk + 178);
    const auto *shk_179 = buffer.data(shk + 179);
    const auto *shk_180 = buffer.data(shk + 180);
    const auto *shk_181 = buffer.data(shk + 181);
    const auto *shk_182 = buffer.data(shk + 182);
    const auto *shk_183 = buffer.data(shk + 183);
    const auto *shk_185 = buffer.data(shk + 185);
    const auto *shk_186 = buffer.data(shk + 186);
    const auto *shk_188 = buffer.data(shk + 188);
    const auto *shk_189 = buffer.data(shk + 189);
    const auto *shk_190 = buffer.data(shk + 190);
    const auto *shk_192 = buffer.data(shk + 192);
    const auto *shk_193 = buffer.data(shk + 193);
    const auto *shk_194 = buffer.data(shk + 194);
    const auto *shk_195 = buffer.data(shk + 195);
    const auto *shk_197 = buffer.data(shk + 197);
    const auto *shk_198 = buffer.data(shk + 198);
    const auto *shk_199 = buffer.data(shk + 199);
    const auto *shk_200 = buffer.data(shk + 200);
    const auto *shk_208 = buffer.data(shk + 208);
    const auto *shk_210 = buffer.data(shk + 210);
    const auto *shk_211 = buffer.data(shk + 211);
    const auto *shk_212 = buffer.data(shk + 212);
    const auto *shk_213 = buffer.data(shk + 213);
    const auto *shk_214 = buffer.data(shk + 214);
    const auto *shk_215 = buffer.data(shk + 215);
    const auto *shk_216 = buffer.data(shk + 216);
    const auto *shk_218 = buffer.data(shk + 218);
    const auto *shk_279 = buffer.data(shk + 279);
    const auto *shk_280 = buffer.data(shk + 280);
    const auto *shk_281 = buffer.data(shk + 281);
    const auto *shk_282 = buffer.data(shk + 282);
    const auto *shk_283 = buffer.data(shk + 283);
    const auto *shk_284 = buffer.data(shk + 284);
    const auto *shk_285 = buffer.data(shk + 285);
    const auto *shk_286 = buffer.data(shk + 286);
    const auto *shk_287 = buffer.data(shk + 287);
    const auto *shk_316 = buffer.data(shk + 316);
    const auto *shk_317 = buffer.data(shk + 317);
    const auto *shk_318 = buffer.data(shk + 318);
    const auto *shk_319 = buffer.data(shk + 319);
    const auto *shk_320 = buffer.data(shk + 320);
    const auto *shk_321 = buffer.data(shk + 321);
    const auto *shk_322 = buffer.data(shk + 322);
    const auto *shk_323 = buffer.data(shk + 323);
    const auto *shk_324 = buffer.data(shk + 324);
    const auto *shk_327 = buffer.data(shk + 327);
    const auto *shk_329 = buffer.data(shk + 329);
    const auto *shk_330 = buffer.data(shk + 330);
    const auto *shk_333 = buffer.data(shk + 333);
    const auto *shk_334 = buffer.data(shk + 334);
    const auto *shk_336 = buffer.data(shk + 336);
    const auto *shk_338 = buffer.data(shk + 338);
    const auto *shk_339 = buffer.data(shk + 339);
    const auto *shk_341 = buffer.data(shk + 341);
    const auto *shk_342 = buffer.data(shk + 342);
    const auto *shk_344 = buffer.data(shk + 344);
    const auto *shk_345 = buffer.data(shk + 345);
    const auto *shk_347 = buffer.data(shk + 347);
    const auto *shk_348 = buffer.data(shk + 348);
    const auto *shk_349 = buffer.data(shk + 349);
    const auto *shk_351 = buffer.data(shk + 351);
    const auto *shk_352 = buffer.data(shk + 352);
    const auto *shk_353 = buffer.data(shk + 353);
    const auto *shk_354 = buffer.data(shk + 354);
    const auto *shk_355 = buffer.data(shk + 355);
    const auto *shk_356 = buffer.data(shk + 356);
    const auto *shk_357 = buffer.data(shk + 357);
    const auto *shk_358 = buffer.data(shk + 358);
    const auto *shk_359 = buffer.data(shk + 359);
    const auto *shk_360 = buffer.data(shk + 360);
    const auto *shk_363 = buffer.data(shk + 363);
    const auto *shk_365 = buffer.data(shk + 365);

    const auto *shl1_171 = buffer.data(shl1 + 171);
    const auto *shl1_225 = buffer.data(shl1 + 225);
    const auto *shl1_228 = buffer.data(shl1 + 228);
    const auto *shl1_230 = buffer.data(shl1 + 230);
    const auto *shl1_231 = buffer.data(shl1 + 231);
    const auto *shl1_234 = buffer.data(shl1 + 234);
    const auto *shl1_235 = buffer.data(shl1 + 235);
    const auto *shl1_237 = buffer.data(shl1 + 237);
    const auto *shl1_239 = buffer.data(shl1 + 239);
    const auto *shl1_240 = buffer.data(shl1 + 240);
    const auto *shl1_242 = buffer.data(shl1 + 242);
    const auto *shl1_243 = buffer.data(shl1 + 243);
    const auto *shl1_245 = buffer.data(shl1 + 245);
    const auto *shl1_246 = buffer.data(shl1 + 246);
    const auto *shl1_248 = buffer.data(shl1 + 248);
    const auto *shl1_249 = buffer.data(shl1 + 249);
    const auto *shl1_250 = buffer.data(shl1 + 250);
    const auto *shl1_252 = buffer.data(shl1 + 252);
    const auto *shl1_269 = buffer.data(shl1 + 269);

    const auto *sii0_219 = buffer.data(sii0 + 219);
    const auto *sii0_220 = buffer.data(sii0 + 220);
    const auto *sii0_221 = buffer.data(sii0 + 221);
    const auto *sii0_222 = buffer.data(sii0 + 222);
    const auto *sii0_223 = buffer.data(sii0 + 223);
    const auto *sii0_245 = buffer.data(sii0 + 245);
    const auto *sii0_247 = buffer.data(sii0 + 247);
    const auto *sii0_248 = buffer.data(sii0 + 248);
    const auto *sii0_249 = buffer.data(sii0 + 249);
    const auto *sii0_250 = buffer.data(sii0 + 250);
    const auto *sii0_251 = buffer.data(sii0 + 251);
    const auto *sii0_252 = buffer.data(sii0 + 252);
    const auto *sii0_255 = buffer.data(sii0 + 255);
    const auto *sii0_257 = buffer.data(sii0 + 257);
    const auto *sii0_258 = buffer.data(sii0 + 258);
    const auto *sii0_261 = buffer.data(sii0 + 261);
    const auto *sii0_262 = buffer.data(sii0 + 262);
    const auto *sii0_264 = buffer.data(sii0 + 264);
    const auto *sii0_266 = buffer.data(sii0 + 266);
    const auto *sii0_267 = buffer.data(sii0 + 267);
    const auto *sii0_269 = buffer.data(sii0 + 269);
    const auto *sii0_270 = buffer.data(sii0 + 270);
    const auto *sii0_272 = buffer.data(sii0 + 272);
    const auto *sii0_273 = buffer.data(sii0 + 273);
    const auto *sii0_275 = buffer.data(sii0 + 275);
    const auto *sii0_276 = buffer.data(sii0 + 276);
    const auto *sii0_277 = buffer.data(sii0 + 277);
    const auto *sii0_278 = buffer.data(sii0 + 278);
    const auto *sii0_279 = buffer.data(sii0 + 279);
    const auto *sii0_280 = buffer.data(sii0 + 280);
    const auto *sii0_283 = buffer.data(sii0 + 283);
    const auto *sii0_285 = buffer.data(sii0 + 285);

    const auto *sii1_219 = buffer.data(sii1 + 219);
    const auto *sii1_220 = buffer.data(sii1 + 220);
    const auto *sii1_221 = buffer.data(sii1 + 221);
    const auto *sii1_222 = buffer.data(sii1 + 222);
    const auto *sii1_223 = buffer.data(sii1 + 223);
    const auto *sii1_245 = buffer.data(sii1 + 245);
    const auto *sii1_247 = buffer.data(sii1 + 247);
    const auto *sii1_248 = buffer.data(sii1 + 248);
    const auto *sii1_249 = buffer.data(sii1 + 249);
    const auto *sii1_250 = buffer.data(sii1 + 250);
    const auto *sii1_251 = buffer.data(sii1 + 251);
    const auto *sii1_252 = buffer.data(sii1 + 252);
    const auto *sii1_255 = buffer.data(sii1 + 255);
    const auto *sii1_257 = buffer.data(sii1 + 257);
    const auto *sii1_258 = buffer.data(sii1 + 258);
    const auto *sii1_261 = buffer.data(sii1 + 261);
    const auto *sii1_262 = buffer.data(sii1 + 262);
    const auto *sii1_264 = buffer.data(sii1 + 264);
    const auto *sii1_266 = buffer.data(sii1 + 266);
    const auto *sii1_267 = buffer.data(sii1 + 267);
    const auto *sii1_269 = buffer.data(sii1 + 269);
    const auto *sii1_270 = buffer.data(sii1 + 270);
    const auto *sii1_272 = buffer.data(sii1 + 272);
    const auto *sii1_273 = buffer.data(sii1 + 273);
    const auto *sii1_275 = buffer.data(sii1 + 275);
    const auto *sii1_276 = buffer.data(sii1 + 276);
    const auto *sii1_277 = buffer.data(sii1 + 277);
    const auto *sii1_278 = buffer.data(sii1 + 278);
    const auto *sii1_279 = buffer.data(sii1 + 279);
    const auto *sii1_280 = buffer.data(sii1 + 280);
    const auto *sii1_283 = buffer.data(sii1 + 283);
    const auto *sii1_285 = buffer.data(sii1 + 285);

    const auto *sik_272 = buffer.data(sik + 272);
    const auto *sik_279 = buffer.data(sik + 279);
    const auto *sik_280 = buffer.data(sik + 280);
    const auto *sik_281 = buffer.data(sik + 281);
    const auto *sik_282 = buffer.data(sik + 282);
    const auto *sik_283 = buffer.data(sik + 283);
    const auto *sik_284 = buffer.data(sik + 284);
    const auto *sik_285 = buffer.data(sik + 285);
    const auto *sik_286 = buffer.data(sik + 286);
    const auto *sik_287 = buffer.data(sik + 287);
    const auto *sik_288 = buffer.data(sik + 288);
    const auto *sik_290 = buffer.data(sik + 290);
    const auto *sik_291 = buffer.data(sik + 291);
    const auto *sik_293 = buffer.data(sik + 293);
    const auto *sik_294 = buffer.data(sik + 294);
    const auto *sik_297 = buffer.data(sik + 297);
    const auto *sik_298 = buffer.data(sik + 298);
    const auto *sik_302 = buffer.data(sik + 302);
    const auto *sik_303 = buffer.data(sik + 303);
    const auto *sik_308 = buffer.data(sik + 308);
    const auto *sik_316 = buffer.data(sik + 316);
    const auto *sik_317 = buffer.data(sik + 317);
    const auto *sik_318 = buffer.data(sik + 318);
    const auto *sik_319 = buffer.data(sik + 319);
    const auto *sik_320 = buffer.data(sik + 320);
    const auto *sik_321 = buffer.data(sik + 321);
    const auto *sik_322 = buffer.data(sik + 322);
    const auto *sik_323 = buffer.data(sik + 323);
    const auto *sik_324 = buffer.data(sik + 324);
    const auto *sik_326 = buffer.data(sik + 326);
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
    const auto *sik_362 = buffer.data(sik + 362);
    const auto *sik_363 = buffer.data(sik + 363);
    const auto *sik_365 = buffer.data(sik + 365);

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, shk_164, shk_279, shk_280, \
                         shk_281, sii0_223, sii1_223, sik_272, sik_279, sik_280, \
                         sik_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * shk_164[k]
                   + f_3 * pc_y[k] * sik_272[k];

        t_342[k] = f_17 * shk_279[k]
                   + f_12 * sii0_223[k]
                   - f_13 * sii1_223[k]
                   + f_3 * pc_x[k] * sik_279[k];

        t_343[k] = f_17 * shk_280[k]
                   + f_3 * pc_x[k] * sik_280[k];

        t_344[k] = f_17 * shk_281[k]
                   + f_3 * pc_x[k] * sik_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, shk_282, shk_283, shk_284, \
                         shk_285, shk_286, sik_282, sik_283, sik_284, sik_285, \
                         sik_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_17 * shk_282[k]
                   + f_3 * pc_x[k] * sik_282[k];

        t_346[k] = f_17 * shk_283[k]
                   + f_3 * pc_x[k] * sik_283[k];

        t_347[k] = f_17 * shk_284[k]
                   + f_3 * pc_x[k] * sik_284[k];

        t_348[k] = f_17 * shk_285[k]
                   + f_3 * pc_x[k] * sik_285[k];

        t_349[k] = f_17 * shk_286[k]
                   + f_3 * pc_x[k] * sik_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pb_z, pc_x, pc_z, shl0_171, shk_136, shk_287, \
                         shl1_171, sik_280, sik_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_17 * shk_287[k]
                   + f_3 * pc_x[k] * sik_287[k];

        t_351[k] = pb_z[k] * shl0_171[k]
                   - f_14 * pc_z[k] * shl1_171[k];

        t_352[k] = f_15 * shk_136[k]
                   + f_3 * pc_z[k] * sik_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, shk_174, shk_175, shk_176, sii0_219, \
                         sii0_220, sii0_221, sii1_219, sii1_220, sii1_221, sik_282, sik_283, \
                         sik_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * shk_174[k]
                   + f_4 * sii0_219[k]
                   - f_5 * sii1_219[k]
                   + f_3 * pc_y[k] * sik_282[k];

        t_354[k] = f_16 * shk_175[k]
                   + f_6 * sii0_220[k]
                   - f_7 * sii1_220[k]
                   + f_3 * pc_y[k] * sik_283[k];

        t_355[k] = f_16 * shk_176[k]
                   + f_8 * sii0_221[k]
                   - f_9 * sii1_221[k]
                   + f_3 * pc_y[k] * sik_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, shk_177, shk_178, shk_179, sii0_222, \
                         sii0_223, sii1_222, sii1_223, sik_285, sik_286, \
                         sik_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * shk_177[k]
                   + f_10 * sii0_222[k]
                   - f_11 * sii1_222[k]
                   + f_3 * pc_y[k] * sik_285[k];

        t_357[k] = f_16 * shk_178[k]
                   + f_12 * sii0_223[k]
                   - f_13 * sii1_223[k]
                   + f_3 * pc_y[k] * sik_286[k];

        t_358[k] = f_16 * shk_179[k]
                   + f_3 * pc_y[k] * sik_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, shl0_225, shk_143, \
                         shk_144, shk_180, shl1_225, sii0_223, sii1_223, sik_287, \
                         sik_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * shk_143[k]
                   + f_1 * sii0_223[k]
                   - f_2 * sii1_223[k]
                   + f_3 * pc_z[k] * sik_287[k];

        t_360[k] = pb_y[k] * shl0_225[k]
                   - f_14 * pc_y[k] * shl1_225[k];

        t_361[k] = f_15 * shk_180[k]
                   + f_3 * pc_y[k] * sik_288[k];

        t_362[k] = f_16 * shk_144[k]
                   + f_3 * pc_z[k] * sik_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, shl0_228, shl0_230, shl0_231, \
                         shk_181, shk_182, shk_183, shl1_228, shl1_230, shl1_231, \
                         sik_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_y[k] * shl0_228[k]
                   + f_16 * shk_181[k]
                   - f_14 * pc_y[k] * shl1_228[k];

        t_364[k] = f_15 * shk_182[k]
                   + f_3 * pc_y[k] * sik_290[k];

        t_365[k] = pb_y[k] * shl0_230[k]
                   - f_14 * pc_y[k] * shl1_230[k];

        t_366[k] = pb_y[k] * shl0_231[k]
                   + f_17 * shk_183[k]
                   - f_14 * pc_y[k] * shl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pb_y, pc_y, pc_z, shl0_234, shl0_235, \
                         shk_147, shk_185, shk_186, shl1_234, shl1_235, sik_291, \
                         sik_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * shk_147[k]
                   + f_3 * pc_z[k] * sik_291[k];

        t_368[k] = f_15 * shk_185[k]
                   + f_3 * pc_y[k] * sik_293[k];

        t_369[k] = pb_y[k] * shl0_234[k]
                   - f_14 * pc_y[k] * shl1_234[k];

        t_370[k] = pb_y[k] * shl0_235[k]
                   + f_18 * shk_186[k]
                   - f_14 * pc_y[k] * shl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pb_y, pc_y, pc_z, shl0_237, shl0_239, \
                         shk_150, shk_188, shk_189, shl1_237, shl1_239, sik_294, \
                         sik_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * shk_150[k]
                   + f_3 * pc_z[k] * sik_294[k];

        t_372[k] = pb_y[k] * shl0_237[k]
                   + f_16 * shk_188[k]
                   - f_14 * pc_y[k] * shl1_237[k];

        t_373[k] = f_15 * shk_189[k]
                   + f_3 * pc_y[k] * sik_297[k];

        t_374[k] = pb_y[k] * shl0_239[k]
                   - f_14 * pc_y[k] * shl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_y, pc_y, pc_z, shl0_240, shl0_242, shk_154, \
                         shk_190, shk_192, shl1_240, shl1_242, \
                         sik_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pb_y[k] * shl0_240[k]
                   + f_19 * shk_190[k]
                   - f_14 * pc_y[k] * shl1_240[k];

        t_376[k] = f_16 * shk_154[k]
                   + f_3 * pc_z[k] * sik_298[k];

        t_377[k] = pb_y[k] * shl0_242[k]
                   + f_17 * shk_192[k]
                   - f_14 * pc_y[k] * shl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_y, shl0_243, shl0_245, shl0_246, \
                         shk_193, shk_194, shk_195, shl1_243, shl1_245, shl1_246, \
                         sik_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * shl0_243[k]
                   + f_16 * shk_193[k]
                   - f_14 * pc_y[k] * shl1_243[k];

        t_379[k] = f_15 * shk_194[k]
                   + f_3 * pc_y[k] * sik_302[k];

        t_380[k] = pb_y[k] * shl0_245[k]
                   - f_14 * pc_y[k] * shl1_245[k];

        t_381[k] = pb_y[k] * shl0_246[k]
                   + f_0 * shk_195[k]
                   - f_14 * pc_y[k] * shl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_y, pc_y, pc_z, shl0_248, shl0_249, shk_159, \
                         shk_197, shk_198, shl1_248, shl1_249, \
                         sik_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * shk_159[k]
                   + f_3 * pc_z[k] * sik_303[k];

        t_383[k] = pb_y[k] * shl0_248[k]
                   + f_18 * shk_197[k]
                   - f_14 * pc_y[k] * shl1_248[k];

        t_384[k] = pb_y[k] * shl0_249[k]
                   + f_17 * shk_198[k]
                   - f_14 * pc_y[k] * shl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_y, pc_x, pc_y, shl0_250, shl0_252, \
                         shk_199, shk_200, shk_316, shl1_250, shl1_252, sik_308, \
                         sik_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * shl0_250[k]
                   + f_16 * shk_199[k]
                   - f_14 * pc_y[k] * shl1_250[k];

        t_386[k] = f_15 * shk_200[k]
                   + f_3 * pc_y[k] * sik_308[k];

        t_387[k] = pb_y[k] * shl0_252[k]
                   - f_14 * pc_y[k] * shl1_252[k];

        t_388[k] = f_17 * shk_316[k]
                   + f_3 * pc_x[k] * sik_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, shk_317, shk_318, shk_319, \
                         shk_320, shk_321, sik_317, sik_318, sik_319, sik_320, \
                         sik_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_17 * shk_317[k]
                   + f_3 * pc_x[k] * sik_317[k];

        t_390[k] = f_17 * shk_318[k]
                   + f_3 * pc_x[k] * sik_318[k];

        t_391[k] = f_17 * shk_319[k]
                   + f_3 * pc_x[k] * sik_319[k];

        t_392[k] = f_17 * shk_320[k]
                   + f_3 * pc_x[k] * sik_320[k];

        t_393[k] = f_17 * shk_321[k]
                   + f_3 * pc_x[k] * sik_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, shk_172, shk_208, \
                         shk_322, shk_323, sii0_245, sii1_245, sik_316, sik_322, \
                         sik_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_17 * shk_322[k]
                   + f_3 * pc_x[k] * sik_322[k];

        t_395[k] = f_17 * shk_323[k]
                   + f_3 * pc_x[k] * sik_323[k];

        t_396[k] = f_15 * shk_208[k]
                   + f_1 * sii0_245[k]
                   - f_2 * sii1_245[k]
                   + f_3 * pc_y[k] * sik_316[k];

        t_397[k] = f_16 * shk_172[k]
                   + f_3 * pc_z[k] * sik_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, shk_210, shk_211, shk_212, sii0_247, \
                         sii0_248, sii0_249, sii1_247, sii1_248, sii1_249, sik_318, sik_319, \
                         sik_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * shk_210[k]
                   + f_4 * sii0_247[k]
                   - f_5 * sii1_247[k]
                   + f_3 * pc_y[k] * sik_318[k];

        t_399[k] = f_15 * shk_211[k]
                   + f_6 * sii0_248[k]
                   - f_7 * sii1_248[k]
                   + f_3 * pc_y[k] * sik_319[k];

        t_400[k] = f_15 * shk_212[k]
                   + f_8 * sii0_249[k]
                   - f_9 * sii1_249[k]
                   + f_3 * pc_y[k] * sik_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, shk_213, shk_214, shk_215, sii0_250, \
                         sii0_251, sii1_250, sii1_251, sik_321, sik_322, \
                         sik_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * shk_213[k]
                   + f_10 * sii0_250[k]
                   - f_11 * sii1_250[k]
                   + f_3 * pc_y[k] * sik_321[k];

        t_402[k] = f_15 * shk_214[k]
                   + f_12 * sii0_251[k]
                   - f_13 * sii1_251[k]
                   + f_3 * pc_y[k] * sik_322[k];

        t_403[k] = f_15 * shk_215[k]
                   + f_3 * pc_y[k] * sik_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_y, pc_x, pc_y, pc_z, shl0_269, \
                         shk_180, shk_324, shl1_269, sii0_252, sii1_252, \
                         sik_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * shl0_269[k]
                   - f_14 * pc_y[k] * shl1_269[k];

        t_405[k] = f_17 * shk_324[k]
                   + f_1 * sii0_252[k]
                   - f_2 * sii1_252[k]
                   + f_3 * pc_x[k] * sik_324[k];

        t_406[k] = f_3 * pc_y[k] * sik_324[k];

        t_407[k] = f_17 * shk_180[k]
                   + f_3 * pc_z[k] * sik_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, shk_327, shk_329, sii0_255, \
                         sii0_257, sii1_255, sii1_257, sik_326, sik_327, \
                         sik_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_17 * shk_327[k]
                   + f_4 * sii0_255[k]
                   - f_5 * sii1_255[k]
                   + f_3 * pc_x[k] * sik_327[k];

        t_409[k] = f_3 * pc_y[k] * sik_326[k];

        t_410[k] = f_17 * shk_329[k]
                   + f_4 * sii0_257[k]
                   - f_5 * sii1_257[k]
                   + f_3 * pc_x[k] * sik_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_x, pc_y, pc_z, shk_183, shk_330, sii0_258, \
                         sii1_258, sik_327, sik_329, sik_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_17 * shk_330[k]
                   + f_6 * sii0_258[k]
                   - f_7 * sii1_258[k]
                   + f_3 * pc_x[k] * sik_330[k];

        t_412[k] = f_17 * shk_183[k]
                   + f_3 * pc_z[k] * sik_327[k];

        t_413[k] = f_3 * pc_y[k] * sik_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_z, shk_186, shk_333, shk_334, sii0_261, \
                         sii0_262, sii1_261, sii1_262, sik_330, sik_333, \
                         sik_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_17 * shk_333[k]
                   + f_6 * sii0_261[k]
                   - f_7 * sii1_261[k]
                   + f_3 * pc_x[k] * sik_333[k];

        t_415[k] = f_17 * shk_334[k]
                   + f_8 * sii0_262[k]
                   - f_9 * sii1_262[k]
                   + f_3 * pc_x[k] * sik_334[k];

        t_416[k] = f_17 * shk_186[k]
                   + f_3 * pc_z[k] * sik_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, shk_336, shk_338, sii0_264, \
                         sii0_266, sii1_264, sii1_266, sik_333, sik_336, \
                         sik_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_17 * shk_336[k]
                   + f_8 * sii0_264[k]
                   - f_9 * sii1_264[k]
                   + f_3 * pc_x[k] * sik_336[k];

        t_418[k] = f_3 * pc_y[k] * sik_333[k];

        t_419[k] = f_17 * shk_338[k]
                   + f_8 * sii0_266[k]
                   - f_9 * sii1_266[k]
                   + f_3 * pc_x[k] * sik_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_x, pc_z, shk_190, shk_339, shk_341, sii0_267, \
                         sii0_269, sii1_267, sii1_269, sik_334, sik_339, \
                         sik_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_17 * shk_339[k]
                   + f_10 * sii0_267[k]
                   - f_11 * sii1_267[k]
                   + f_3 * pc_x[k] * sik_339[k];

        t_421[k] = f_17 * shk_190[k]
                   + f_3 * pc_z[k] * sik_334[k];

        t_422[k] = f_17 * shk_341[k]
                   + f_10 * sii0_269[k]
                   - f_11 * sii1_269[k]
                   + f_3 * pc_x[k] * sik_341[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, shk_342, shk_344, sii0_270, \
                         sii0_272, sii1_270, sii1_272, sik_338, sik_342, \
                         sik_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_17 * shk_342[k]
                   + f_10 * sii0_270[k]
                   - f_11 * sii1_270[k]
                   + f_3 * pc_x[k] * sik_342[k];

        t_424[k] = f_3 * pc_y[k] * sik_338[k];

        t_425[k] = f_17 * shk_344[k]
                   + f_10 * sii0_272[k]
                   - f_11 * sii1_272[k]
                   + f_3 * pc_x[k] * sik_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, shk_195, shk_345, shk_347, sii0_273, \
                         sii0_275, sii1_273, sii1_275, sik_339, sik_345, \
                         sik_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_17 * shk_345[k]
                   + f_12 * sii0_273[k]
                   - f_13 * sii1_273[k]
                   + f_3 * pc_x[k] * sik_345[k];

        t_427[k] = f_17 * shk_195[k]
                   + f_3 * pc_z[k] * sik_339[k];

        t_428[k] = f_17 * shk_347[k]
                   + f_12 * sii0_275[k]
                   - f_13 * sii1_275[k]
                   + f_3 * pc_x[k] * sik_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, shk_348, shk_349, sii0_276, \
                         sii0_277, sii1_276, sii1_277, sik_344, sik_348, \
                         sik_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_17 * shk_348[k]
                   + f_12 * sii0_276[k]
                   - f_13 * sii1_276[k]
                   + f_3 * pc_x[k] * sik_348[k];

        t_430[k] = f_17 * shk_349[k]
                   + f_12 * sii0_277[k]
                   - f_13 * sii1_277[k]
                   + f_3 * pc_x[k] * sik_349[k];

        t_431[k] = f_3 * pc_y[k] * sik_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, shk_351, shk_352, shk_353, shk_354, \
                         sii0_279, sii1_279, sik_351, sik_352, sik_353, \
                         sik_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_17 * shk_351[k]
                   + f_12 * sii0_279[k]
                   - f_13 * sii1_279[k]
                   + f_3 * pc_x[k] * sik_351[k];

        t_433[k] = f_17 * shk_352[k]
                   + f_3 * pc_x[k] * sik_352[k];

        t_434[k] = f_17 * shk_353[k]
                   + f_3 * pc_x[k] * sik_353[k];

        t_435[k] = f_17 * shk_354[k]
                   + f_3 * pc_x[k] * sik_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, shk_355, shk_356, shk_357, \
                         shk_358, shk_359, sik_355, sik_356, sik_357, sik_358, \
                         sik_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_17 * shk_355[k]
                   + f_3 * pc_x[k] * sik_355[k];

        t_437[k] = f_17 * shk_356[k]
                   + f_3 * pc_x[k] * sik_356[k];

        t_438[k] = f_17 * shk_357[k]
                   + f_3 * pc_x[k] * sik_357[k];

        t_439[k] = f_17 * shk_358[k]
                   + f_3 * pc_x[k] * sik_358[k];

        t_440[k] = f_17 * shk_359[k]
                   + f_3 * pc_x[k] * sik_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, shk_208, sii0_273, sii0_275, \
                         sii0_276, sii1_273, sii1_275, sii1_276, sik_352, sik_354, \
                         sik_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * sii0_273[k]
                   - f_2 * sii1_273[k]
                   + f_3 * pc_y[k] * sik_352[k];

        t_442[k] = f_17 * shk_208[k]
                   + f_3 * pc_z[k] * sik_352[k];

        t_443[k] = f_4 * sii0_275[k]
                   - f_5 * sii1_275[k]
                   + f_3 * pc_y[k] * sik_354[k];

        t_444[k] = f_6 * sii0_276[k]
                   - f_7 * sii1_276[k]
                   + f_3 * pc_y[k] * sik_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, sii0_277, sii0_278, sii0_279, \
                         sii1_277, sii1_278, sii1_279, sik_356, sik_357, sik_358, \
                         sik_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * sii0_277[k]
                   - f_9 * sii1_277[k]
                   + f_3 * pc_y[k] * sik_356[k];

        t_446[k] = f_10 * sii0_278[k]
                   - f_11 * sii1_278[k]
                   + f_3 * pc_y[k] * sik_357[k];

        t_447[k] = f_12 * sii0_279[k]
                   - f_13 * sii1_279[k]
                   + f_3 * pc_y[k] * sik_358[k];

        t_448[k] = f_3 * pc_y[k] * sik_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, shk_215, shk_216, \
                         shk_360, sii0_279, sii0_280, sii1_279, sii1_280, sik_359, \
                         sik_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * shk_215[k]
                   + f_1 * sii0_279[k]
                   - f_2 * sii1_279[k]
                   + f_3 * pc_z[k] * sik_359[k];

        t_450[k] = f_16 * shk_360[k]
                   + f_1 * sii0_280[k]
                   - f_2 * sii1_280[k]
                   + f_3 * pc_x[k] * sik_360[k];

        t_451[k] = f_18 * shk_216[k]
                   + f_3 * pc_y[k] * sik_360[k];

        t_452[k] = f_3 * pc_z[k] * sik_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, shk_218, shk_363, shk_365, sii0_283, \
                         sii0_285, sii1_283, sii1_285, sik_362, sik_363, \
                         sik_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_16 * shk_363[k]
                   + f_4 * sii0_283[k]
                   - f_5 * sii1_283[k]
                   + f_3 * pc_x[k] * sik_363[k];

        t_454[k] = f_18 * shk_218[k]
                   + f_3 * pc_y[k] * sik_362[k];

        t_455[k] = f_16 * shk_365[k]
                   + f_4 * sii0_285[k]
                   - f_5 * sii1_285[k]
                   + f_3 * pc_x[k] * sik_365[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
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

    const auto *shl0_270 = buffer.data(shl0 + 270);
    const auto *shl0_273 = buffer.data(shl0 + 273);
    const auto *shl0_276 = buffer.data(shl0 + 276);
    const auto *shl0_280 = buffer.data(shl0 + 280);
    const auto *shl0_282 = buffer.data(shl0 + 282);
    const auto *shl0_285 = buffer.data(shl0 + 285);
    const auto *shl0_287 = buffer.data(shl0 + 287);
    const auto *shl0_288 = buffer.data(shl0 + 288);
    const auto *shl0_291 = buffer.data(shl0 + 291);
    const auto *shl0_293 = buffer.data(shl0 + 293);
    const auto *shl0_294 = buffer.data(shl0 + 294);
    const auto *shl0_295 = buffer.data(shl0 + 295);
    const auto *shl0_306 = buffer.data(shl0 + 306);

    const auto *shk_216 = buffer.data(shk + 216);
    const auto *shk_219 = buffer.data(shk + 219);
    const auto *shk_221 = buffer.data(shk + 221);
    const auto *shk_222 = buffer.data(shk + 222);
    const auto *shk_223 = buffer.data(shk + 223);
    const auto *shk_225 = buffer.data(shk + 225);
    const auto *shk_226 = buffer.data(shk + 226);
    const auto *shk_227 = buffer.data(shk + 227);
    const auto *shk_228 = buffer.data(shk + 228);
    const auto *shk_230 = buffer.data(shk + 230);
    const auto *shk_231 = buffer.data(shk + 231);
    const auto *shk_232 = buffer.data(shk + 232);
    const auto *shk_233 = buffer.data(shk + 233);
    const auto *shk_234 = buffer.data(shk + 234);
    const auto *shk_236 = buffer.data(shk + 236);
    const auto *shk_244 = buffer.data(shk + 244);
    const auto *shk_246 = buffer.data(shk + 246);
    const auto *shk_247 = buffer.data(shk + 247);
    const auto *shk_248 = buffer.data(shk + 248);
    const auto *shk_249 = buffer.data(shk + 249);
    const auto *shk_250 = buffer.data(shk + 250);
    const auto *shk_251 = buffer.data(shk + 251);
    const auto *shk_252 = buffer.data(shk + 252);
    const auto *shk_254 = buffer.data(shk + 254);
    const auto *shk_255 = buffer.data(shk + 255);
    const auto *shk_257 = buffer.data(shk + 257);
    const auto *shk_258 = buffer.data(shk + 258);
    const auto *shk_261 = buffer.data(shk + 261);
    const auto *shk_262 = buffer.data(shk + 262);
    const auto *shk_266 = buffer.data(shk + 266);
    const auto *shk_267 = buffer.data(shk + 267);
    const auto *shk_272 = buffer.data(shk + 272);
    const auto *shk_282 = buffer.data(shk + 282);
    const auto *shk_283 = buffer.data(shk + 283);
    const auto *shk_284 = buffer.data(shk + 284);
    const auto *shk_285 = buffer.data(shk + 285);
    const auto *shk_286 = buffer.data(shk + 286);
    const auto *shk_287 = buffer.data(shk + 287);
    const auto *shk_288 = buffer.data(shk + 288);
    const auto *shk_290 = buffer.data(shk + 290);
    const auto *shk_293 = buffer.data(shk + 293);
    const auto *shk_297 = buffer.data(shk + 297);
    const auto *shk_302 = buffer.data(shk + 302);
    const auto *shk_366 = buffer.data(shk + 366);
    const auto *shk_369 = buffer.data(shk + 369);
    const auto *shk_370 = buffer.data(shk + 370);
    const auto *shk_372 = buffer.data(shk + 372);
    const auto *shk_374 = buffer.data(shk + 374);
    const auto *shk_375 = buffer.data(shk + 375);
    const auto *shk_377 = buffer.data(shk + 377);
    const auto *shk_378 = buffer.data(shk + 378);
    const auto *shk_380 = buffer.data(shk + 380);
    const auto *shk_381 = buffer.data(shk + 381);
    const auto *shk_383 = buffer.data(shk + 383);
    const auto *shk_384 = buffer.data(shk + 384);
    const auto *shk_385 = buffer.data(shk + 385);
    const auto *shk_387 = buffer.data(shk + 387);
    const auto *shk_388 = buffer.data(shk + 388);
    const auto *shk_389 = buffer.data(shk + 389);
    const auto *shk_390 = buffer.data(shk + 390);
    const auto *shk_391 = buffer.data(shk + 391);
    const auto *shk_392 = buffer.data(shk + 392);
    const auto *shk_393 = buffer.data(shk + 393);
    const auto *shk_394 = buffer.data(shk + 394);
    const auto *shk_395 = buffer.data(shk + 395);
    const auto *shk_401 = buffer.data(shk + 401);
    const auto *shk_405 = buffer.data(shk + 405);
    const auto *shk_410 = buffer.data(shk + 410);
    const auto *shk_416 = buffer.data(shk + 416);
    const auto *shk_423 = buffer.data(shk + 423);
    const auto *shk_424 = buffer.data(shk + 424);
    const auto *shk_425 = buffer.data(shk + 425);
    const auto *shk_426 = buffer.data(shk + 426);
    const auto *shk_427 = buffer.data(shk + 427);
    const auto *shk_428 = buffer.data(shk + 428);
    const auto *shk_429 = buffer.data(shk + 429);
    const auto *shk_430 = buffer.data(shk + 430);
    const auto *shk_431 = buffer.data(shk + 431);
    const auto *shk_432 = buffer.data(shk + 432);
    const auto *shk_435 = buffer.data(shk + 435);
    const auto *shk_437 = buffer.data(shk + 437);
    const auto *shk_438 = buffer.data(shk + 438);
    const auto *shk_441 = buffer.data(shk + 441);
    const auto *shk_442 = buffer.data(shk + 442);
    const auto *shk_444 = buffer.data(shk + 444);
    const auto *shk_446 = buffer.data(shk + 446);
    const auto *shk_447 = buffer.data(shk + 447);
    const auto *shk_449 = buffer.data(shk + 449);
    const auto *shk_450 = buffer.data(shk + 450);
    const auto *shk_452 = buffer.data(shk + 452);
    const auto *shk_453 = buffer.data(shk + 453);

    const auto *shl1_270 = buffer.data(shl1 + 270);
    const auto *shl1_273 = buffer.data(shl1 + 273);
    const auto *shl1_276 = buffer.data(shl1 + 276);
    const auto *shl1_280 = buffer.data(shl1 + 280);
    const auto *shl1_282 = buffer.data(shl1 + 282);
    const auto *shl1_285 = buffer.data(shl1 + 285);
    const auto *shl1_287 = buffer.data(shl1 + 287);
    const auto *shl1_288 = buffer.data(shl1 + 288);
    const auto *shl1_291 = buffer.data(shl1 + 291);
    const auto *shl1_293 = buffer.data(shl1 + 293);
    const auto *shl1_294 = buffer.data(shl1 + 294);
    const auto *shl1_295 = buffer.data(shl1 + 295);
    const auto *shl1_306 = buffer.data(shl1 + 306);

    const auto *sii0_286 = buffer.data(sii0 + 286);
    const auto *sii0_289 = buffer.data(sii0 + 289);
    const auto *sii0_290 = buffer.data(sii0 + 290);
    const auto *sii0_292 = buffer.data(sii0 + 292);
    const auto *sii0_294 = buffer.data(sii0 + 294);
    const auto *sii0_295 = buffer.data(sii0 + 295);
    const auto *sii0_297 = buffer.data(sii0 + 297);
    const auto *sii0_298 = buffer.data(sii0 + 298);
    const auto *sii0_300 = buffer.data(sii0 + 300);
    const auto *sii0_301 = buffer.data(sii0 + 301);
    const auto *sii0_303 = buffer.data(sii0 + 303);
    const auto *sii0_304 = buffer.data(sii0 + 304);
    const auto *sii0_305 = buffer.data(sii0 + 305);
    const auto *sii0_306 = buffer.data(sii0 + 306);
    const auto *sii0_307 = buffer.data(sii0 + 307);
    const auto *sii0_313 = buffer.data(sii0 + 313);
    const auto *sii0_317 = buffer.data(sii0 + 317);
    const auto *sii0_322 = buffer.data(sii0 + 322);
    const auto *sii0_328 = buffer.data(sii0 + 328);
    const auto *sii0_331 = buffer.data(sii0 + 331);
    const auto *sii0_332 = buffer.data(sii0 + 332);
    const auto *sii0_333 = buffer.data(sii0 + 333);
    const auto *sii0_334 = buffer.data(sii0 + 334);
    const auto *sii0_335 = buffer.data(sii0 + 335);
    const auto *sii0_336 = buffer.data(sii0 + 336);
    const auto *sii0_339 = buffer.data(sii0 + 339);
    const auto *sii0_341 = buffer.data(sii0 + 341);
    const auto *sii0_342 = buffer.data(sii0 + 342);
    const auto *sii0_345 = buffer.data(sii0 + 345);
    const auto *sii0_346 = buffer.data(sii0 + 346);
    const auto *sii0_348 = buffer.data(sii0 + 348);
    const auto *sii0_350 = buffer.data(sii0 + 350);
    const auto *sii0_351 = buffer.data(sii0 + 351);
    const auto *sii0_353 = buffer.data(sii0 + 353);
    const auto *sii0_354 = buffer.data(sii0 + 354);
    const auto *sii0_356 = buffer.data(sii0 + 356);
    const auto *sii0_357 = buffer.data(sii0 + 357);

    const auto *sii1_286 = buffer.data(sii1 + 286);
    const auto *sii1_289 = buffer.data(sii1 + 289);
    const auto *sii1_290 = buffer.data(sii1 + 290);
    const auto *sii1_292 = buffer.data(sii1 + 292);
    const auto *sii1_294 = buffer.data(sii1 + 294);
    const auto *sii1_295 = buffer.data(sii1 + 295);
    const auto *sii1_297 = buffer.data(sii1 + 297);
    const auto *sii1_298 = buffer.data(sii1 + 298);
    const auto *sii1_300 = buffer.data(sii1 + 300);
    const auto *sii1_301 = buffer.data(sii1 + 301);
    const auto *sii1_303 = buffer.data(sii1 + 303);
    const auto *sii1_304 = buffer.data(sii1 + 304);
    const auto *sii1_305 = buffer.data(sii1 + 305);
    const auto *sii1_306 = buffer.data(sii1 + 306);
    const auto *sii1_307 = buffer.data(sii1 + 307);
    const auto *sii1_313 = buffer.data(sii1 + 313);
    const auto *sii1_317 = buffer.data(sii1 + 317);
    const auto *sii1_322 = buffer.data(sii1 + 322);
    const auto *sii1_328 = buffer.data(sii1 + 328);
    const auto *sii1_331 = buffer.data(sii1 + 331);
    const auto *sii1_332 = buffer.data(sii1 + 332);
    const auto *sii1_333 = buffer.data(sii1 + 333);
    const auto *sii1_334 = buffer.data(sii1 + 334);
    const auto *sii1_335 = buffer.data(sii1 + 335);
    const auto *sii1_336 = buffer.data(sii1 + 336);
    const auto *sii1_339 = buffer.data(sii1 + 339);
    const auto *sii1_341 = buffer.data(sii1 + 341);
    const auto *sii1_342 = buffer.data(sii1 + 342);
    const auto *sii1_345 = buffer.data(sii1 + 345);
    const auto *sii1_346 = buffer.data(sii1 + 346);
    const auto *sii1_348 = buffer.data(sii1 + 348);
    const auto *sii1_350 = buffer.data(sii1 + 350);
    const auto *sii1_351 = buffer.data(sii1 + 351);
    const auto *sii1_353 = buffer.data(sii1 + 353);
    const auto *sii1_354 = buffer.data(sii1 + 354);
    const auto *sii1_356 = buffer.data(sii1 + 356);
    const auto *sii1_357 = buffer.data(sii1 + 357);

    const auto *sik_363 = buffer.data(sik + 363);
    const auto *sik_365 = buffer.data(sik + 365);
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
    const auto *sik_396 = buffer.data(sik + 396);
    const auto *sik_398 = buffer.data(sik + 398);
    const auto *sik_399 = buffer.data(sik + 399);
    const auto *sik_401 = buffer.data(sik + 401);
    const auto *sik_402 = buffer.data(sik + 402);
    const auto *sik_405 = buffer.data(sik + 405);
    const auto *sik_406 = buffer.data(sik + 406);
    const auto *sik_410 = buffer.data(sik + 410);
    const auto *sik_411 = buffer.data(sik + 411);
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
    const auto *sik_434 = buffer.data(sik + 434);
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

#pragma omp simd aligned(t_456, t_457, t_458, pc_x, pc_y, pc_z, shk_221, shk_366, sii0_286, \
                         sii1_286, sik_363, sik_365, sik_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_16 * shk_366[k]
                   + f_6 * sii0_286[k]
                   - f_7 * sii1_286[k]
                   + f_3 * pc_x[k] * sik_366[k];

        t_457[k] = f_3 * pc_z[k] * sik_363[k];

        t_458[k] = f_18 * shk_221[k]
                   + f_3 * pc_y[k] * sik_365[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_z, shk_369, shk_370, sii0_289, \
                         sii0_290, sii1_289, sii1_290, sik_366, sik_369, \
                         sik_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_16 * shk_369[k]
                   + f_6 * sii0_289[k]
                   - f_7 * sii1_289[k]
                   + f_3 * pc_x[k] * sik_369[k];

        t_460[k] = f_16 * shk_370[k]
                   + f_8 * sii0_290[k]
                   - f_9 * sii1_290[k]
                   + f_3 * pc_x[k] * sik_370[k];

        t_461[k] = f_3 * pc_z[k] * sik_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_x, pc_y, shk_225, shk_372, shk_374, sii0_292, \
                         sii0_294, sii1_292, sii1_294, sik_369, sik_372, \
                         sik_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_16 * shk_372[k]
                   + f_8 * sii0_292[k]
                   - f_9 * sii1_292[k]
                   + f_3 * pc_x[k] * sik_372[k];

        t_463[k] = f_18 * shk_225[k]
                   + f_3 * pc_y[k] * sik_369[k];

        t_464[k] = f_16 * shk_374[k]
                   + f_8 * sii0_294[k]
                   - f_9 * sii1_294[k]
                   + f_3 * pc_x[k] * sik_374[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, shk_375, shk_377, sii0_295, \
                         sii0_297, sii1_295, sii1_297, sik_370, sik_375, \
                         sik_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_16 * shk_375[k]
                   + f_10 * sii0_295[k]
                   - f_11 * sii1_295[k]
                   + f_3 * pc_x[k] * sik_375[k];

        t_466[k] = f_3 * pc_z[k] * sik_370[k];

        t_467[k] = f_16 * shk_377[k]
                   + f_10 * sii0_297[k]
                   - f_11 * sii1_297[k]
                   + f_3 * pc_x[k] * sik_377[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, shk_230, shk_378, shk_380, sii0_298, \
                         sii0_300, sii1_298, sii1_300, sik_374, sik_378, \
                         sik_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_16 * shk_378[k]
                   + f_10 * sii0_298[k]
                   - f_11 * sii1_298[k]
                   + f_3 * pc_x[k] * sik_378[k];

        t_469[k] = f_18 * shk_230[k]
                   + f_3 * pc_y[k] * sik_374[k];

        t_470[k] = f_16 * shk_380[k]
                   + f_10 * sii0_300[k]
                   - f_11 * sii1_300[k]
                   + f_3 * pc_x[k] * sik_380[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, shk_381, shk_383, sii0_301, \
                         sii0_303, sii1_301, sii1_303, sik_375, sik_381, \
                         sik_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_16 * shk_381[k]
                   + f_12 * sii0_301[k]
                   - f_13 * sii1_301[k]
                   + f_3 * pc_x[k] * sik_381[k];

        t_472[k] = f_3 * pc_z[k] * sik_375[k];

        t_473[k] = f_16 * shk_383[k]
                   + f_12 * sii0_303[k]
                   - f_13 * sii1_303[k]
                   + f_3 * pc_x[k] * sik_383[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, shk_236, shk_384, shk_385, sii0_304, \
                         sii0_305, sii1_304, sii1_305, sik_380, sik_384, \
                         sik_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_16 * shk_384[k]
                   + f_12 * sii0_304[k]
                   - f_13 * sii1_304[k]
                   + f_3 * pc_x[k] * sik_384[k];

        t_475[k] = f_16 * shk_385[k]
                   + f_12 * sii0_305[k]
                   - f_13 * sii1_305[k]
                   + f_3 * pc_x[k] * sik_385[k];

        t_476[k] = f_18 * shk_236[k]
                   + f_3 * pc_y[k] * sik_380[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pc_x, shk_387, shk_388, shk_389, shk_390, \
                         sii0_307, sii1_307, sik_387, sik_388, sik_389, \
                         sik_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_16 * shk_387[k]
                   + f_12 * sii0_307[k]
                   - f_13 * sii1_307[k]
                   + f_3 * pc_x[k] * sik_387[k];

        t_478[k] = f_16 * shk_388[k]
                   + f_3 * pc_x[k] * sik_388[k];

        t_479[k] = f_16 * shk_389[k]
                   + f_3 * pc_x[k] * sik_389[k];

        t_480[k] = f_16 * shk_390[k]
                   + f_3 * pc_x[k] * sik_390[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pc_x, shk_391, shk_392, shk_393, \
                         shk_394, shk_395, sik_391, sik_392, sik_393, sik_394, \
                         sik_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_16 * shk_391[k]
                   + f_3 * pc_x[k] * sik_391[k];

        t_482[k] = f_16 * shk_392[k]
                   + f_3 * pc_x[k] * sik_392[k];

        t_483[k] = f_16 * shk_393[k]
                   + f_3 * pc_x[k] * sik_393[k];

        t_484[k] = f_16 * shk_394[k]
                   + f_3 * pc_x[k] * sik_394[k];

        t_485[k] = f_16 * shk_395[k]
                   + f_3 * pc_x[k] * sik_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_y, pc_z, shk_244, shk_246, sii0_301, \
                         sii0_303, sii1_301, sii1_303, sik_388, \
                         sik_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_18 * shk_244[k]
                   + f_1 * sii0_301[k]
                   - f_2 * sii1_301[k]
                   + f_3 * pc_y[k] * sik_388[k];

        t_487[k] = f_3 * pc_z[k] * sik_388[k];

        t_488[k] = f_18 * shk_246[k]
                   + f_4 * sii0_303[k]
                   - f_5 * sii1_303[k]
                   + f_3 * pc_y[k] * sik_390[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_y, shk_247, shk_248, shk_249, sii0_304, \
                         sii0_305, sii0_306, sii1_304, sii1_305, sii1_306, sik_391, sik_392, \
                         sik_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_18 * shk_247[k]
                   + f_6 * sii0_304[k]
                   - f_7 * sii1_304[k]
                   + f_3 * pc_y[k] * sik_391[k];

        t_490[k] = f_18 * shk_248[k]
                   + f_8 * sii0_305[k]
                   - f_9 * sii1_305[k]
                   + f_3 * pc_y[k] * sik_392[k];

        t_491[k] = f_18 * shk_249[k]
                   + f_10 * sii0_306[k]
                   - f_11 * sii1_306[k]
                   + f_3 * pc_y[k] * sik_393[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_z, pc_y, pc_z, shl0_270, shk_250, \
                         shk_251, shl1_270, sii0_307, sii1_307, sik_394, \
                         sik_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_18 * shk_250[k]
                   + f_12 * sii0_307[k]
                   - f_13 * sii1_307[k]
                   + f_3 * pc_y[k] * sik_394[k];

        t_493[k] = f_18 * shk_251[k]
                   + f_3 * pc_y[k] * sik_395[k];

        t_494[k] = f_1 * sii0_307[k]
                   - f_2 * sii1_307[k]
                   + f_3 * pc_z[k] * sik_395[k];

        t_495[k] = pb_z[k] * shl0_270[k]
                   - f_14 * pc_z[k] * shl1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, shl0_273, shk_216, \
                         shk_252, shk_254, shl1_273, sik_396, sik_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * shk_252[k]
                   + f_3 * pc_y[k] * sik_396[k];

        t_497[k] = f_15 * shk_216[k]
                   + f_3 * pc_z[k] * sik_396[k];

        t_498[k] = pb_z[k] * shl0_273[k]
                   - f_14 * pc_z[k] * shl1_273[k];

        t_499[k] = f_17 * shk_254[k]
                   + f_3 * pc_y[k] * sik_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pb_z, pc_x, pc_z, shl0_276, shk_219, shk_401, \
                         shl1_276, sii0_313, sii1_313, sik_399, \
                         sik_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_16 * shk_401[k]
                   + f_4 * sii0_313[k]
                   - f_5 * sii1_313[k]
                   + f_3 * pc_x[k] * sik_401[k];

        t_501[k] = pb_z[k] * shl0_276[k]
                   - f_14 * pc_z[k] * shl1_276[k];

        t_502[k] = f_15 * shk_219[k]
                   + f_3 * pc_z[k] * sik_399[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pb_z, pc_x, pc_y, pc_z, shl0_280, shk_257, \
                         shk_405, shl1_280, sii0_317, sii1_317, sik_401, \
                         sik_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * shk_257[k]
                   + f_3 * pc_y[k] * sik_401[k];

        t_504[k] = f_16 * shk_405[k]
                   + f_6 * sii0_317[k]
                   - f_7 * sii1_317[k]
                   + f_3 * pc_x[k] * sik_405[k];

        t_505[k] = pb_z[k] * shl0_280[k]
                   - f_14 * pc_z[k] * shl1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_z, pc_y, pc_z, shl0_282, shk_222, shk_223, \
                         shk_261, shl1_282, sik_402, sik_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_15 * shk_222[k]
                   + f_3 * pc_z[k] * sik_402[k];

        t_507[k] = pb_z[k] * shl0_282[k]
                   + f_16 * shk_223[k]
                   - f_14 * pc_z[k] * shl1_282[k];

        t_508[k] = f_17 * shk_261[k]
                   + f_3 * pc_y[k] * sik_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_z, pc_x, pc_z, shl0_285, shk_226, shk_410, \
                         shl1_285, sii0_322, sii1_322, sik_406, \
                         sik_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_16 * shk_410[k]
                   + f_8 * sii0_322[k]
                   - f_9 * sii1_322[k]
                   + f_3 * pc_x[k] * sik_410[k];

        t_510[k] = pb_z[k] * shl0_285[k]
                   - f_14 * pc_z[k] * shl1_285[k];

        t_511[k] = f_15 * shk_226[k]
                   + f_3 * pc_z[k] * sik_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_z, pc_y, pc_z, shl0_287, shl0_288, shk_227, \
                         shk_228, shk_266, shl1_287, shl1_288, \
                         sik_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_z[k] * shl0_287[k]
                   + f_16 * shk_227[k]
                   - f_14 * pc_z[k] * shl1_287[k];

        t_513[k] = pb_z[k] * shl0_288[k]
                   + f_17 * shk_228[k]
                   - f_14 * pc_z[k] * shl1_288[k];

        t_514[k] = f_17 * shk_266[k]
                   + f_3 * pc_y[k] * sik_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_z, pc_x, pc_z, shl0_291, shk_231, shk_416, \
                         shl1_291, sii0_328, sii1_328, sik_411, \
                         sik_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_16 * shk_416[k]
                   + f_10 * sii0_328[k]
                   - f_11 * sii1_328[k]
                   + f_3 * pc_x[k] * sik_416[k];

        t_516[k] = pb_z[k] * shl0_291[k]
                   - f_14 * pc_z[k] * shl1_291[k];

        t_517[k] = f_15 * shk_231[k]
                   + f_3 * pc_z[k] * sik_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pb_z, pc_z, shl0_293, shl0_294, shl0_295, \
                         shk_232, shk_233, shk_234, shl1_293, shl1_294, \
                         shl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_z[k] * shl0_293[k]
                   + f_16 * shk_232[k]
                   - f_14 * pc_z[k] * shl1_293[k];

        t_519[k] = pb_z[k] * shl0_294[k]
                   + f_17 * shk_233[k]
                   - f_14 * pc_z[k] * shl1_294[k];

        t_520[k] = pb_z[k] * shl0_295[k]
                   + f_18 * shk_234[k]
                   - f_14 * pc_z[k] * shl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, shk_272, shk_423, shk_424, \
                         shk_425, sii0_335, sii1_335, sik_416, sik_423, sik_424, \
                         sik_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * shk_272[k]
                   + f_3 * pc_y[k] * sik_416[k];

        t_522[k] = f_16 * shk_423[k]
                   + f_12 * sii0_335[k]
                   - f_13 * sii1_335[k]
                   + f_3 * pc_x[k] * sik_423[k];

        t_523[k] = f_16 * shk_424[k]
                   + f_3 * pc_x[k] * sik_424[k];

        t_524[k] = f_16 * shk_425[k]
                   + f_3 * pc_x[k] * sik_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, shk_426, shk_427, shk_428, \
                         shk_429, shk_430, sik_426, sik_427, sik_428, sik_429, \
                         sik_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_16 * shk_426[k]
                   + f_3 * pc_x[k] * sik_426[k];

        t_526[k] = f_16 * shk_427[k]
                   + f_3 * pc_x[k] * sik_427[k];

        t_527[k] = f_16 * shk_428[k]
                   + f_3 * pc_x[k] * sik_428[k];

        t_528[k] = f_16 * shk_429[k]
                   + f_3 * pc_x[k] * sik_429[k];

        t_529[k] = f_16 * shk_430[k]
                   + f_3 * pc_x[k] * sik_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pb_z, pc_x, pc_z, shl0_306, shk_244, shk_431, \
                         shl1_306, sik_424, sik_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_16 * shk_431[k]
                   + f_3 * pc_x[k] * sik_431[k];

        t_531[k] = pb_z[k] * shl0_306[k]
                   - f_14 * pc_z[k] * shl1_306[k];

        t_532[k] = f_15 * shk_244[k]
                   + f_3 * pc_z[k] * sik_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, shk_282, shk_283, shk_284, sii0_331, \
                         sii0_332, sii0_333, sii1_331, sii1_332, sii1_333, sik_426, sik_427, \
                         sik_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * shk_282[k]
                   + f_4 * sii0_331[k]
                   - f_5 * sii1_331[k]
                   + f_3 * pc_y[k] * sik_426[k];

        t_534[k] = f_17 * shk_283[k]
                   + f_6 * sii0_332[k]
                   - f_7 * sii1_332[k]
                   + f_3 * pc_y[k] * sik_427[k];

        t_535[k] = f_17 * shk_284[k]
                   + f_8 * sii0_333[k]
                   - f_9 * sii1_333[k]
                   + f_3 * pc_y[k] * sik_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, shk_285, shk_286, shk_287, sii0_334, \
                         sii0_335, sii1_334, sii1_335, sik_429, sik_430, \
                         sik_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * shk_285[k]
                   + f_10 * sii0_334[k]
                   - f_11 * sii1_334[k]
                   + f_3 * pc_y[k] * sik_429[k];

        t_537[k] = f_17 * shk_286[k]
                   + f_12 * sii0_335[k]
                   - f_13 * sii1_335[k]
                   + f_3 * pc_y[k] * sik_430[k];

        t_538[k] = f_17 * shk_287[k]
                   + f_3 * pc_y[k] * sik_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, shk_251, shk_288, shk_432, \
                         sii0_335, sii0_336, sii1_335, sii1_336, sik_431, \
                         sik_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * shk_251[k]
                   + f_1 * sii0_335[k]
                   - f_2 * sii1_335[k]
                   + f_3 * pc_z[k] * sik_431[k];

        t_540[k] = f_16 * shk_432[k]
                   + f_1 * sii0_336[k]
                   - f_2 * sii1_336[k]
                   + f_3 * pc_x[k] * sik_432[k];

        t_541[k] = f_16 * shk_288[k]
                   + f_3 * pc_y[k] * sik_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, shk_252, shk_290, shk_435, \
                         sii0_339, sii1_339, sik_432, sik_434, \
                         sik_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * shk_252[k]
                   + f_3 * pc_z[k] * sik_432[k];

        t_543[k] = f_16 * shk_435[k]
                   + f_4 * sii0_339[k]
                   - f_5 * sii1_339[k]
                   + f_3 * pc_x[k] * sik_435[k];

        t_544[k] = f_16 * shk_290[k]
                   + f_3 * pc_y[k] * sik_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, shk_255, shk_437, shk_438, sii0_341, \
                         sii0_342, sii1_341, sii1_342, sik_435, sik_437, \
                         sik_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_16 * shk_437[k]
                   + f_4 * sii0_341[k]
                   - f_5 * sii1_341[k]
                   + f_3 * pc_x[k] * sik_437[k];

        t_546[k] = f_16 * shk_438[k]
                   + f_6 * sii0_342[k]
                   - f_7 * sii1_342[k]
                   + f_3 * pc_x[k] * sik_438[k];

        t_547[k] = f_16 * shk_255[k]
                   + f_3 * pc_z[k] * sik_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, shk_293, shk_441, shk_442, sii0_345, \
                         sii0_346, sii1_345, sii1_346, sik_437, sik_441, \
                         sik_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * shk_293[k]
                   + f_3 * pc_y[k] * sik_437[k];

        t_549[k] = f_16 * shk_441[k]
                   + f_6 * sii0_345[k]
                   - f_7 * sii1_345[k]
                   + f_3 * pc_x[k] * sik_441[k];

        t_550[k] = f_16 * shk_442[k]
                   + f_8 * sii0_346[k]
                   - f_9 * sii1_346[k]
                   + f_3 * pc_x[k] * sik_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, shk_258, shk_297, shk_444, \
                         sii0_348, sii1_348, sik_438, sik_441, \
                         sik_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * shk_258[k]
                   + f_3 * pc_z[k] * sik_438[k];

        t_552[k] = f_16 * shk_444[k]
                   + f_8 * sii0_348[k]
                   - f_9 * sii1_348[k]
                   + f_3 * pc_x[k] * sik_444[k];

        t_553[k] = f_16 * shk_297[k]
                   + f_3 * pc_y[k] * sik_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, shk_262, shk_446, shk_447, sii0_350, \
                         sii0_351, sii1_350, sii1_351, sik_442, sik_446, \
                         sik_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_16 * shk_446[k]
                   + f_8 * sii0_350[k]
                   - f_9 * sii1_350[k]
                   + f_3 * pc_x[k] * sik_446[k];

        t_555[k] = f_16 * shk_447[k]
                   + f_10 * sii0_351[k]
                   - f_11 * sii1_351[k]
                   + f_3 * pc_x[k] * sik_447[k];

        t_556[k] = f_16 * shk_262[k]
                   + f_3 * pc_z[k] * sik_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, shk_302, shk_449, shk_450, sii0_353, \
                         sii0_354, sii1_353, sii1_354, sik_446, sik_449, \
                         sik_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_16 * shk_449[k]
                   + f_10 * sii0_353[k]
                   - f_11 * sii1_353[k]
                   + f_3 * pc_x[k] * sik_449[k];

        t_558[k] = f_16 * shk_450[k]
                   + f_10 * sii0_354[k]
                   - f_11 * sii1_354[k]
                   + f_3 * pc_x[k] * sik_450[k];

        t_559[k] = f_16 * shk_302[k]
                   + f_3 * pc_y[k] * sik_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, shk_267, shk_452, shk_453, sii0_356, \
                         sii0_357, sii1_356, sii1_357, sik_447, sik_452, \
                         sik_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_16 * shk_452[k]
                   + f_10 * sii0_356[k]
                   - f_11 * sii1_356[k]
                   + f_3 * pc_x[k] * sik_452[k];

        t_561[k] = f_16 * shk_453[k]
                   + f_12 * sii0_357[k]
                   - f_13 * sii1_357[k]
                   + f_3 * pc_x[k] * sik_453[k];

        t_562[k] = f_16 * shk_267[k]
                   + f_3 * pc_z[k] * sik_447[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *shl0_405 = buffer.data(shl0 + 405);
    const auto *shl0_408 = buffer.data(shl0 + 408);
    const auto *shl0_410 = buffer.data(shl0 + 410);
    const auto *shl0_411 = buffer.data(shl0 + 411);
    const auto *shl0_414 = buffer.data(shl0 + 414);
    const auto *shl0_415 = buffer.data(shl0 + 415);
    const auto *shl0_417 = buffer.data(shl0 + 417);
    const auto *shl0_419 = buffer.data(shl0 + 419);
    const auto *shl0_420 = buffer.data(shl0 + 420);
    const auto *shl0_422 = buffer.data(shl0 + 422);
    const auto *shl0_423 = buffer.data(shl0 + 423);
    const auto *shl0_425 = buffer.data(shl0 + 425);
    const auto *shl0_426 = buffer.data(shl0 + 426);
    const auto *shl0_428 = buffer.data(shl0 + 428);
    const auto *shl0_429 = buffer.data(shl0 + 429);
    const auto *shl0_430 = buffer.data(shl0 + 430);
    const auto *shl0_432 = buffer.data(shl0 + 432);
    const auto *shl0_449 = buffer.data(shl0 + 449);

    const auto *shk_280 = buffer.data(shk + 280);
    const auto *shk_287 = buffer.data(shk + 287);
    const auto *shk_288 = buffer.data(shk + 288);
    const auto *shk_291 = buffer.data(shk + 291);
    const auto *shk_294 = buffer.data(shk + 294);
    const auto *shk_298 = buffer.data(shk + 298);
    const auto *shk_303 = buffer.data(shk + 303);
    const auto *shk_308 = buffer.data(shk + 308);
    const auto *shk_316 = buffer.data(shk + 316);
    const auto *shk_318 = buffer.data(shk + 318);
    const auto *shk_319 = buffer.data(shk + 319);
    const auto *shk_320 = buffer.data(shk + 320);
    const auto *shk_321 = buffer.data(shk + 321);
    const auto *shk_322 = buffer.data(shk + 322);
    const auto *shk_323 = buffer.data(shk + 323);
    const auto *shk_324 = buffer.data(shk + 324);
    const auto *shk_325 = buffer.data(shk + 325);
    const auto *shk_326 = buffer.data(shk + 326);
    const auto *shk_327 = buffer.data(shk + 327);
    const auto *shk_329 = buffer.data(shk + 329);
    const auto *shk_330 = buffer.data(shk + 330);
    const auto *shk_332 = buffer.data(shk + 332);
    const auto *shk_333 = buffer.data(shk + 333);
    const auto *shk_334 = buffer.data(shk + 334);
    const auto *shk_336 = buffer.data(shk + 336);
    const auto *shk_337 = buffer.data(shk + 337);
    const auto *shk_338 = buffer.data(shk + 338);
    const auto *shk_339 = buffer.data(shk + 339);
    const auto *shk_341 = buffer.data(shk + 341);
    const auto *shk_342 = buffer.data(shk + 342);
    const auto *shk_343 = buffer.data(shk + 343);
    const auto *shk_344 = buffer.data(shk + 344);
    const auto *shk_352 = buffer.data(shk + 352);
    const auto *shk_354 = buffer.data(shk + 354);
    const auto *shk_355 = buffer.data(shk + 355);
    const auto *shk_356 = buffer.data(shk + 356);
    const auto *shk_357 = buffer.data(shk + 357);
    const auto *shk_358 = buffer.data(shk + 358);
    const auto *shk_359 = buffer.data(shk + 359);
    const auto *shk_455 = buffer.data(shk + 455);
    const auto *shk_456 = buffer.data(shk + 456);
    const auto *shk_457 = buffer.data(shk + 457);
    const auto *shk_459 = buffer.data(shk + 459);
    const auto *shk_460 = buffer.data(shk + 460);
    const auto *shk_461 = buffer.data(shk + 461);
    const auto *shk_462 = buffer.data(shk + 462);
    const auto *shk_463 = buffer.data(shk + 463);
    const auto *shk_464 = buffer.data(shk + 464);
    const auto *shk_465 = buffer.data(shk + 465);
    const auto *shk_466 = buffer.data(shk + 466);
    const auto *shk_467 = buffer.data(shk + 467);
    const auto *shk_496 = buffer.data(shk + 496);
    const auto *shk_497 = buffer.data(shk + 497);
    const auto *shk_498 = buffer.data(shk + 498);
    const auto *shk_499 = buffer.data(shk + 499);
    const auto *shk_500 = buffer.data(shk + 500);
    const auto *shk_501 = buffer.data(shk + 501);
    const auto *shk_502 = buffer.data(shk + 502);
    const auto *shk_503 = buffer.data(shk + 503);
    const auto *shk_504 = buffer.data(shk + 504);
    const auto *shk_507 = buffer.data(shk + 507);
    const auto *shk_509 = buffer.data(shk + 509);
    const auto *shk_510 = buffer.data(shk + 510);
    const auto *shk_513 = buffer.data(shk + 513);
    const auto *shk_514 = buffer.data(shk + 514);
    const auto *shk_516 = buffer.data(shk + 516);
    const auto *shk_518 = buffer.data(shk + 518);
    const auto *shk_519 = buffer.data(shk + 519);
    const auto *shk_521 = buffer.data(shk + 521);
    const auto *shk_522 = buffer.data(shk + 522);
    const auto *shk_524 = buffer.data(shk + 524);
    const auto *shk_525 = buffer.data(shk + 525);
    const auto *shk_527 = buffer.data(shk + 527);
    const auto *shk_528 = buffer.data(shk + 528);
    const auto *shk_529 = buffer.data(shk + 529);
    const auto *shk_531 = buffer.data(shk + 531);
    const auto *shk_532 = buffer.data(shk + 532);
    const auto *shk_533 = buffer.data(shk + 533);
    const auto *shk_534 = buffer.data(shk + 534);
    const auto *shk_535 = buffer.data(shk + 535);
    const auto *shk_536 = buffer.data(shk + 536);
    const auto *shk_537 = buffer.data(shk + 537);
    const auto *shk_538 = buffer.data(shk + 538);
    const auto *shk_539 = buffer.data(shk + 539);

    const auto *shl1_405 = buffer.data(shl1 + 405);
    const auto *shl1_408 = buffer.data(shl1 + 408);
    const auto *shl1_410 = buffer.data(shl1 + 410);
    const auto *shl1_411 = buffer.data(shl1 + 411);
    const auto *shl1_414 = buffer.data(shl1 + 414);
    const auto *shl1_415 = buffer.data(shl1 + 415);
    const auto *shl1_417 = buffer.data(shl1 + 417);
    const auto *shl1_419 = buffer.data(shl1 + 419);
    const auto *shl1_420 = buffer.data(shl1 + 420);
    const auto *shl1_422 = buffer.data(shl1 + 422);
    const auto *shl1_423 = buffer.data(shl1 + 423);
    const auto *shl1_425 = buffer.data(shl1 + 425);
    const auto *shl1_426 = buffer.data(shl1 + 426);
    const auto *shl1_428 = buffer.data(shl1 + 428);
    const auto *shl1_429 = buffer.data(shl1 + 429);
    const auto *shl1_430 = buffer.data(shl1 + 430);
    const auto *shl1_432 = buffer.data(shl1 + 432);
    const auto *shl1_449 = buffer.data(shl1 + 449);

    const auto *sii0_357 = buffer.data(sii0 + 357);
    const auto *sii0_359 = buffer.data(sii0 + 359);
    const auto *sii0_360 = buffer.data(sii0 + 360);
    const auto *sii0_361 = buffer.data(sii0 + 361);
    const auto *sii0_362 = buffer.data(sii0 + 362);
    const auto *sii0_363 = buffer.data(sii0 + 363);
    const auto *sii0_385 = buffer.data(sii0 + 385);
    const auto *sii0_387 = buffer.data(sii0 + 387);
    const auto *sii0_388 = buffer.data(sii0 + 388);
    const auto *sii0_389 = buffer.data(sii0 + 389);
    const auto *sii0_390 = buffer.data(sii0 + 390);
    const auto *sii0_391 = buffer.data(sii0 + 391);
    const auto *sii0_392 = buffer.data(sii0 + 392);
    const auto *sii0_395 = buffer.data(sii0 + 395);
    const auto *sii0_397 = buffer.data(sii0 + 397);
    const auto *sii0_398 = buffer.data(sii0 + 398);
    const auto *sii0_401 = buffer.data(sii0 + 401);
    const auto *sii0_402 = buffer.data(sii0 + 402);
    const auto *sii0_404 = buffer.data(sii0 + 404);
    const auto *sii0_406 = buffer.data(sii0 + 406);
    const auto *sii0_407 = buffer.data(sii0 + 407);
    const auto *sii0_409 = buffer.data(sii0 + 409);
    const auto *sii0_410 = buffer.data(sii0 + 410);
    const auto *sii0_412 = buffer.data(sii0 + 412);
    const auto *sii0_413 = buffer.data(sii0 + 413);
    const auto *sii0_415 = buffer.data(sii0 + 415);
    const auto *sii0_416 = buffer.data(sii0 + 416);
    const auto *sii0_417 = buffer.data(sii0 + 417);
    const auto *sii0_418 = buffer.data(sii0 + 418);
    const auto *sii0_419 = buffer.data(sii0 + 419);

    const auto *sii1_357 = buffer.data(sii1 + 357);
    const auto *sii1_359 = buffer.data(sii1 + 359);
    const auto *sii1_360 = buffer.data(sii1 + 360);
    const auto *sii1_361 = buffer.data(sii1 + 361);
    const auto *sii1_362 = buffer.data(sii1 + 362);
    const auto *sii1_363 = buffer.data(sii1 + 363);
    const auto *sii1_385 = buffer.data(sii1 + 385);
    const auto *sii1_387 = buffer.data(sii1 + 387);
    const auto *sii1_388 = buffer.data(sii1 + 388);
    const auto *sii1_389 = buffer.data(sii1 + 389);
    const auto *sii1_390 = buffer.data(sii1 + 390);
    const auto *sii1_391 = buffer.data(sii1 + 391);
    const auto *sii1_392 = buffer.data(sii1 + 392);
    const auto *sii1_395 = buffer.data(sii1 + 395);
    const auto *sii1_397 = buffer.data(sii1 + 397);
    const auto *sii1_398 = buffer.data(sii1 + 398);
    const auto *sii1_401 = buffer.data(sii1 + 401);
    const auto *sii1_402 = buffer.data(sii1 + 402);
    const auto *sii1_404 = buffer.data(sii1 + 404);
    const auto *sii1_406 = buffer.data(sii1 + 406);
    const auto *sii1_407 = buffer.data(sii1 + 407);
    const auto *sii1_409 = buffer.data(sii1 + 409);
    const auto *sii1_410 = buffer.data(sii1 + 410);
    const auto *sii1_412 = buffer.data(sii1 + 412);
    const auto *sii1_413 = buffer.data(sii1 + 413);
    const auto *sii1_415 = buffer.data(sii1 + 415);
    const auto *sii1_416 = buffer.data(sii1 + 416);
    const auto *sii1_417 = buffer.data(sii1 + 417);
    const auto *sii1_418 = buffer.data(sii1 + 418);
    const auto *sii1_419 = buffer.data(sii1 + 419);

    const auto *sik_452 = buffer.data(sik + 452);
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
    const auto *sik_497 = buffer.data(sik + 497);
    const auto *sik_498 = buffer.data(sik + 498);
    const auto *sik_499 = buffer.data(sik + 499);
    const auto *sik_500 = buffer.data(sik + 500);
    const auto *sik_501 = buffer.data(sik + 501);
    const auto *sik_502 = buffer.data(sik + 502);
    const auto *sik_503 = buffer.data(sik + 503);
    const auto *sik_504 = buffer.data(sik + 504);
    const auto *sik_506 = buffer.data(sik + 506);
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

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, shk_455, shk_456, shk_457, sii0_359, \
                         sii0_360, sii0_361, sii1_359, sii1_360, sii1_361, sik_455, sik_456, \
                         sik_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_16 * shk_455[k]
                   + f_12 * sii0_359[k]
                   - f_13 * sii1_359[k]
                   + f_3 * pc_x[k] * sik_455[k];

        t_564[k] = f_16 * shk_456[k]
                   + f_12 * sii0_360[k]
                   - f_13 * sii1_360[k]
                   + f_3 * pc_x[k] * sik_456[k];

        t_565[k] = f_16 * shk_457[k]
                   + f_12 * sii0_361[k]
                   - f_13 * sii1_361[k]
                   + f_3 * pc_x[k] * sik_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, shk_308, shk_459, shk_460, \
                         shk_461, sii0_363, sii1_363, sik_452, sik_459, sik_460, \
                         sik_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * shk_308[k]
                   + f_3 * pc_y[k] * sik_452[k];

        t_567[k] = f_16 * shk_459[k]
                   + f_12 * sii0_363[k]
                   - f_13 * sii1_363[k]
                   + f_3 * pc_x[k] * sik_459[k];

        t_568[k] = f_16 * shk_460[k]
                   + f_3 * pc_x[k] * sik_460[k];

        t_569[k] = f_16 * shk_461[k]
                   + f_3 * pc_x[k] * sik_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, shk_462, shk_463, shk_464, \
                         shk_465, shk_466, sik_462, sik_463, sik_464, sik_465, \
                         sik_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_16 * shk_462[k]
                   + f_3 * pc_x[k] * sik_462[k];

        t_571[k] = f_16 * shk_463[k]
                   + f_3 * pc_x[k] * sik_463[k];

        t_572[k] = f_16 * shk_464[k]
                   + f_3 * pc_x[k] * sik_464[k];

        t_573[k] = f_16 * shk_465[k]
                   + f_3 * pc_x[k] * sik_465[k];

        t_574[k] = f_16 * shk_466[k]
                   + f_3 * pc_x[k] * sik_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, shk_280, shk_316, shk_467, \
                         sii0_357, sii1_357, sik_460, sik_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_16 * shk_467[k]
                   + f_3 * pc_x[k] * sik_467[k];

        t_576[k] = f_16 * shk_316[k]
                   + f_1 * sii0_357[k]
                   - f_2 * sii1_357[k]
                   + f_3 * pc_y[k] * sik_460[k];

        t_577[k] = f_16 * shk_280[k]
                   + f_3 * pc_z[k] * sik_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, shk_318, shk_319, shk_320, sii0_359, \
                         sii0_360, sii0_361, sii1_359, sii1_360, sii1_361, sik_462, sik_463, \
                         sik_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * shk_318[k]
                   + f_4 * sii0_359[k]
                   - f_5 * sii1_359[k]
                   + f_3 * pc_y[k] * sik_462[k];

        t_579[k] = f_16 * shk_319[k]
                   + f_6 * sii0_360[k]
                   - f_7 * sii1_360[k]
                   + f_3 * pc_y[k] * sik_463[k];

        t_580[k] = f_16 * shk_320[k]
                   + f_8 * sii0_361[k]
                   - f_9 * sii1_361[k]
                   + f_3 * pc_y[k] * sik_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, shk_321, shk_322, shk_323, sii0_362, \
                         sii0_363, sii1_362, sii1_363, sik_465, sik_466, \
                         sik_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * shk_321[k]
                   + f_10 * sii0_362[k]
                   - f_11 * sii1_362[k]
                   + f_3 * pc_y[k] * sik_465[k];

        t_582[k] = f_16 * shk_322[k]
                   + f_12 * sii0_363[k]
                   - f_13 * sii1_363[k]
                   + f_3 * pc_y[k] * sik_466[k];

        t_583[k] = f_16 * shk_323[k]
                   + f_3 * pc_y[k] * sik_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pc_y, pc_z, shl0_405, shk_287, \
                         shk_288, shk_324, shl1_405, sii0_363, sii1_363, sik_467, \
                         sik_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * shk_287[k]
                   + f_1 * sii0_363[k]
                   - f_2 * sii1_363[k]
                   + f_3 * pc_z[k] * sik_467[k];

        t_585[k] = pb_y[k] * shl0_405[k]
                   - f_14 * pc_y[k] * shl1_405[k];

        t_586[k] = f_15 * shk_324[k]
                   + f_3 * pc_y[k] * sik_468[k];

        t_587[k] = f_17 * shk_288[k]
                   + f_3 * pc_z[k] * sik_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_y, pc_y, shl0_408, shl0_410, shl0_411, \
                         shk_325, shk_326, shk_327, shl1_408, shl1_410, shl1_411, \
                         sik_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_y[k] * shl0_408[k]
                   + f_16 * shk_325[k]
                   - f_14 * pc_y[k] * shl1_408[k];

        t_589[k] = f_15 * shk_326[k]
                   + f_3 * pc_y[k] * sik_470[k];

        t_590[k] = pb_y[k] * shl0_410[k]
                   - f_14 * pc_y[k] * shl1_410[k];

        t_591[k] = pb_y[k] * shl0_411[k]
                   + f_17 * shk_327[k]
                   - f_14 * pc_y[k] * shl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pb_y, pc_y, pc_z, shl0_414, shl0_415, \
                         shk_291, shk_329, shk_330, shl1_414, shl1_415, sik_471, \
                         sik_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * shk_291[k]
                   + f_3 * pc_z[k] * sik_471[k];

        t_593[k] = f_15 * shk_329[k]
                   + f_3 * pc_y[k] * sik_473[k];

        t_594[k] = pb_y[k] * shl0_414[k]
                   - f_14 * pc_y[k] * shl1_414[k];

        t_595[k] = pb_y[k] * shl0_415[k]
                   + f_18 * shk_330[k]
                   - f_14 * pc_y[k] * shl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_y, pc_y, pc_z, shl0_417, shl0_419, \
                         shk_294, shk_332, shk_333, shl1_417, shl1_419, sik_474, \
                         sik_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * shk_294[k]
                   + f_3 * pc_z[k] * sik_474[k];

        t_597[k] = pb_y[k] * shl0_417[k]
                   + f_16 * shk_332[k]
                   - f_14 * pc_y[k] * shl1_417[k];

        t_598[k] = f_15 * shk_333[k]
                   + f_3 * pc_y[k] * sik_477[k];

        t_599[k] = pb_y[k] * shl0_419[k]
                   - f_14 * pc_y[k] * shl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pb_y, pc_y, pc_z, shl0_420, shl0_422, shk_298, \
                         shk_334, shk_336, shl1_420, shl1_422, \
                         sik_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pb_y[k] * shl0_420[k]
                   + f_19 * shk_334[k]
                   - f_14 * pc_y[k] * shl1_420[k];

        t_601[k] = f_17 * shk_298[k]
                   + f_3 * pc_z[k] * sik_478[k];

        t_602[k] = pb_y[k] * shl0_422[k]
                   + f_17 * shk_336[k]
                   - f_14 * pc_y[k] * shl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pb_y, pc_y, shl0_423, shl0_425, shl0_426, \
                         shk_337, shk_338, shk_339, shl1_423, shl1_425, shl1_426, \
                         sik_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pb_y[k] * shl0_423[k]
                   + f_16 * shk_337[k]
                   - f_14 * pc_y[k] * shl1_423[k];

        t_604[k] = f_15 * shk_338[k]
                   + f_3 * pc_y[k] * sik_482[k];

        t_605[k] = pb_y[k] * shl0_425[k]
                   - f_14 * pc_y[k] * shl1_425[k];

        t_606[k] = pb_y[k] * shl0_426[k]
                   + f_0 * shk_339[k]
                   - f_14 * pc_y[k] * shl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pb_y, pc_y, pc_z, shl0_428, shl0_429, shk_303, \
                         shk_341, shk_342, shl1_428, shl1_429, \
                         sik_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * shk_303[k]
                   + f_3 * pc_z[k] * sik_483[k];

        t_608[k] = pb_y[k] * shl0_428[k]
                   + f_18 * shk_341[k]
                   - f_14 * pc_y[k] * shl1_428[k];

        t_609[k] = pb_y[k] * shl0_429[k]
                   + f_17 * shk_342[k]
                   - f_14 * pc_y[k] * shl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pb_y, pc_x, pc_y, shl0_430, shl0_432, \
                         shk_343, shk_344, shk_496, shl1_430, shl1_432, sik_488, \
                         sik_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pb_y[k] * shl0_430[k]
                   + f_16 * shk_343[k]
                   - f_14 * pc_y[k] * shl1_430[k];

        t_611[k] = f_15 * shk_344[k]
                   + f_3 * pc_y[k] * sik_488[k];

        t_612[k] = pb_y[k] * shl0_432[k]
                   - f_14 * pc_y[k] * shl1_432[k];

        t_613[k] = f_16 * shk_496[k]
                   + f_3 * pc_x[k] * sik_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, shk_497, shk_498, shk_499, \
                         shk_500, shk_501, sik_497, sik_498, sik_499, sik_500, \
                         sik_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_16 * shk_497[k]
                   + f_3 * pc_x[k] * sik_497[k];

        t_615[k] = f_16 * shk_498[k]
                   + f_3 * pc_x[k] * sik_498[k];

        t_616[k] = f_16 * shk_499[k]
                   + f_3 * pc_x[k] * sik_499[k];

        t_617[k] = f_16 * shk_500[k]
                   + f_3 * pc_x[k] * sik_500[k];

        t_618[k] = f_16 * shk_501[k]
                   + f_3 * pc_x[k] * sik_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, shk_316, shk_352, \
                         shk_502, shk_503, sii0_385, sii1_385, sik_496, sik_502, \
                         sik_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_16 * shk_502[k]
                   + f_3 * pc_x[k] * sik_502[k];

        t_620[k] = f_16 * shk_503[k]
                   + f_3 * pc_x[k] * sik_503[k];

        t_621[k] = f_15 * shk_352[k]
                   + f_1 * sii0_385[k]
                   - f_2 * sii1_385[k]
                   + f_3 * pc_y[k] * sik_496[k];

        t_622[k] = f_17 * shk_316[k]
                   + f_3 * pc_z[k] * sik_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, shk_354, shk_355, shk_356, sii0_387, \
                         sii0_388, sii0_389, sii1_387, sii1_388, sii1_389, sik_498, sik_499, \
                         sik_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * shk_354[k]
                   + f_4 * sii0_387[k]
                   - f_5 * sii1_387[k]
                   + f_3 * pc_y[k] * sik_498[k];

        t_624[k] = f_15 * shk_355[k]
                   + f_6 * sii0_388[k]
                   - f_7 * sii1_388[k]
                   + f_3 * pc_y[k] * sik_499[k];

        t_625[k] = f_15 * shk_356[k]
                   + f_8 * sii0_389[k]
                   - f_9 * sii1_389[k]
                   + f_3 * pc_y[k] * sik_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, shk_357, shk_358, shk_359, sii0_390, \
                         sii0_391, sii1_390, sii1_391, sik_501, sik_502, \
                         sik_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * shk_357[k]
                   + f_10 * sii0_390[k]
                   - f_11 * sii1_390[k]
                   + f_3 * pc_y[k] * sik_501[k];

        t_627[k] = f_15 * shk_358[k]
                   + f_12 * sii0_391[k]
                   - f_13 * sii1_391[k]
                   + f_3 * pc_y[k] * sik_502[k];

        t_628[k] = f_15 * shk_359[k]
                   + f_3 * pc_y[k] * sik_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pb_y, pc_x, pc_y, pc_z, shl0_449, \
                         shk_324, shk_504, shl1_449, sii0_392, sii1_392, \
                         sik_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pb_y[k] * shl0_449[k]
                   - f_14 * pc_y[k] * shl1_449[k];

        t_630[k] = f_16 * shk_504[k]
                   + f_1 * sii0_392[k]
                   - f_2 * sii1_392[k]
                   + f_3 * pc_x[k] * sik_504[k];

        t_631[k] = f_3 * pc_y[k] * sik_504[k];

        t_632[k] = f_18 * shk_324[k]
                   + f_3 * pc_z[k] * sik_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, shk_507, shk_509, sii0_395, \
                         sii0_397, sii1_395, sii1_397, sik_506, sik_507, \
                         sik_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_16 * shk_507[k]
                   + f_4 * sii0_395[k]
                   - f_5 * sii1_395[k]
                   + f_3 * pc_x[k] * sik_507[k];

        t_634[k] = f_3 * pc_y[k] * sik_506[k];

        t_635[k] = f_16 * shk_509[k]
                   + f_4 * sii0_397[k]
                   - f_5 * sii1_397[k]
                   + f_3 * pc_x[k] * sik_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_x, pc_y, pc_z, shk_327, shk_510, sii0_398, \
                         sii1_398, sik_507, sik_509, sik_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_16 * shk_510[k]
                   + f_6 * sii0_398[k]
                   - f_7 * sii1_398[k]
                   + f_3 * pc_x[k] * sik_510[k];

        t_637[k] = f_18 * shk_327[k]
                   + f_3 * pc_z[k] * sik_507[k];

        t_638[k] = f_3 * pc_y[k] * sik_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_z, shk_330, shk_513, shk_514, sii0_401, \
                         sii0_402, sii1_401, sii1_402, sik_510, sik_513, \
                         sik_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_16 * shk_513[k]
                   + f_6 * sii0_401[k]
                   - f_7 * sii1_401[k]
                   + f_3 * pc_x[k] * sik_513[k];

        t_640[k] = f_16 * shk_514[k]
                   + f_8 * sii0_402[k]
                   - f_9 * sii1_402[k]
                   + f_3 * pc_x[k] * sik_514[k];

        t_641[k] = f_18 * shk_330[k]
                   + f_3 * pc_z[k] * sik_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, shk_516, shk_518, sii0_404, \
                         sii0_406, sii1_404, sii1_406, sik_513, sik_516, \
                         sik_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_16 * shk_516[k]
                   + f_8 * sii0_404[k]
                   - f_9 * sii1_404[k]
                   + f_3 * pc_x[k] * sik_516[k];

        t_643[k] = f_3 * pc_y[k] * sik_513[k];

        t_644[k] = f_16 * shk_518[k]
                   + f_8 * sii0_406[k]
                   - f_9 * sii1_406[k]
                   + f_3 * pc_x[k] * sik_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, pc_z, shk_334, shk_519, shk_521, sii0_407, \
                         sii0_409, sii1_407, sii1_409, sik_514, sik_519, \
                         sik_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_16 * shk_519[k]
                   + f_10 * sii0_407[k]
                   - f_11 * sii1_407[k]
                   + f_3 * pc_x[k] * sik_519[k];

        t_646[k] = f_18 * shk_334[k]
                   + f_3 * pc_z[k] * sik_514[k];

        t_647[k] = f_16 * shk_521[k]
                   + f_10 * sii0_409[k]
                   - f_11 * sii1_409[k]
                   + f_3 * pc_x[k] * sik_521[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, shk_522, shk_524, sii0_410, \
                         sii0_412, sii1_410, sii1_412, sik_518, sik_522, \
                         sik_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_16 * shk_522[k]
                   + f_10 * sii0_410[k]
                   - f_11 * sii1_410[k]
                   + f_3 * pc_x[k] * sik_522[k];

        t_649[k] = f_3 * pc_y[k] * sik_518[k];

        t_650[k] = f_16 * shk_524[k]
                   + f_10 * sii0_412[k]
                   - f_11 * sii1_412[k]
                   + f_3 * pc_x[k] * sik_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, pc_z, shk_339, shk_525, shk_527, sii0_413, \
                         sii0_415, sii1_413, sii1_415, sik_519, sik_525, \
                         sik_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_16 * shk_525[k]
                   + f_12 * sii0_413[k]
                   - f_13 * sii1_413[k]
                   + f_3 * pc_x[k] * sik_525[k];

        t_652[k] = f_18 * shk_339[k]
                   + f_3 * pc_z[k] * sik_519[k];

        t_653[k] = f_16 * shk_527[k]
                   + f_12 * sii0_415[k]
                   - f_13 * sii1_415[k]
                   + f_3 * pc_x[k] * sik_527[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, shk_528, shk_529, sii0_416, \
                         sii0_417, sii1_416, sii1_417, sik_524, sik_528, \
                         sik_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_16 * shk_528[k]
                   + f_12 * sii0_416[k]
                   - f_13 * sii1_416[k]
                   + f_3 * pc_x[k] * sik_528[k];

        t_655[k] = f_16 * shk_529[k]
                   + f_12 * sii0_417[k]
                   - f_13 * sii1_417[k]
                   + f_3 * pc_x[k] * sik_529[k];

        t_656[k] = f_3 * pc_y[k] * sik_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, shk_531, shk_532, shk_533, shk_534, \
                         sii0_419, sii1_419, sik_531, sik_532, sik_533, \
                         sik_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_16 * shk_531[k]
                   + f_12 * sii0_419[k]
                   - f_13 * sii1_419[k]
                   + f_3 * pc_x[k] * sik_531[k];

        t_658[k] = f_16 * shk_532[k]
                   + f_3 * pc_x[k] * sik_532[k];

        t_659[k] = f_16 * shk_533[k]
                   + f_3 * pc_x[k] * sik_533[k];

        t_660[k] = f_16 * shk_534[k]
                   + f_3 * pc_x[k] * sik_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, shk_535, shk_536, shk_537, \
                         shk_538, shk_539, sik_535, sik_536, sik_537, sik_538, \
                         sik_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_16 * shk_535[k]
                   + f_3 * pc_x[k] * sik_535[k];

        t_662[k] = f_16 * shk_536[k]
                   + f_3 * pc_x[k] * sik_536[k];

        t_663[k] = f_16 * shk_537[k]
                   + f_3 * pc_x[k] * sik_537[k];

        t_664[k] = f_16 * shk_538[k]
                   + f_3 * pc_x[k] * sik_538[k];

        t_665[k] = f_16 * shk_539[k]
                   + f_3 * pc_x[k] * sik_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, shk_352, sii0_413, sii0_415, \
                         sii0_416, sii1_413, sii1_415, sii1_416, sik_532, sik_534, \
                         sik_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * sii0_413[k]
                   - f_2 * sii1_413[k]
                   + f_3 * pc_y[k] * sik_532[k];

        t_667[k] = f_18 * shk_352[k]
                   + f_3 * pc_z[k] * sik_532[k];

        t_668[k] = f_4 * sii0_415[k]
                   - f_5 * sii1_415[k]
                   + f_3 * pc_y[k] * sik_534[k];

        t_669[k] = f_6 * sii0_416[k]
                   - f_7 * sii1_416[k]
                   + f_3 * pc_y[k] * sik_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, sii0_417, sii0_418, sii0_419, \
                         sii1_417, sii1_418, sii1_419, sik_536, sik_537, sik_538, \
                         sik_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_8 * sii0_417[k]
                   - f_9 * sii1_417[k]
                   + f_3 * pc_y[k] * sik_536[k];

        t_671[k] = f_10 * sii0_418[k]
                   - f_11 * sii1_418[k]
                   + f_3 * pc_y[k] * sik_537[k];

        t_672[k] = f_12 * sii0_419[k]
                   - f_13 * sii1_419[k]
                   + f_3 * pc_y[k] * sik_538[k];

        t_673[k] = f_3 * pc_y[k] * sik_539[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 3.5 / gamma;
    const auto f_2 = 3.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 4.0 / q;

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
    auto *t_782 = buffer.data(target + 782);
    auto *t_783 = buffer.data(target + 783);
    auto *t_784 = buffer.data(target + 784);
    auto *t_785 = buffer.data(target + 785);
    auto *t_786 = buffer.data(target + 786);
    auto *t_787 = buffer.data(target + 787);
    auto *t_788 = buffer.data(target + 788);
    auto *t_789 = buffer.data(target + 789);
    auto *t_790 = buffer.data(target + 790);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shl0_450 = buffer.data(shl0 + 450);
    const auto *shl0_453 = buffer.data(shl0 + 453);
    const auto *shl0_456 = buffer.data(shl0 + 456);
    const auto *shl0_460 = buffer.data(shl0 + 460);
    const auto *shl0_465 = buffer.data(shl0 + 465);
    const auto *shl0_471 = buffer.data(shl0 + 471);
    const auto *shl0_675 = buffer.data(shl0 + 675);
    const auto *shl0_678 = buffer.data(shl0 + 678);
    const auto *shl0_680 = buffer.data(shl0 + 680);
    const auto *shl0_681 = buffer.data(shl0 + 681);
    const auto *shl0_684 = buffer.data(shl0 + 684);
    const auto *shl0_685 = buffer.data(shl0 + 685);
    const auto *shl0_687 = buffer.data(shl0 + 687);
    const auto *shl0_689 = buffer.data(shl0 + 689);
    const auto *shl0_690 = buffer.data(shl0 + 690);
    const auto *shl0_692 = buffer.data(shl0 + 692);
    const auto *shl0_693 = buffer.data(shl0 + 693);
    const auto *shl0_695 = buffer.data(shl0 + 695);
    const auto *shl0_696 = buffer.data(shl0 + 696);
    const auto *shl0_698 = buffer.data(shl0 + 698);
    const auto *shl0_699 = buffer.data(shl0 + 699);
    const auto *shl0_700 = buffer.data(shl0 + 700);
    const auto *shl0_702 = buffer.data(shl0 + 702);
    const auto *shl0_711 = buffer.data(shl0 + 711);
    const auto *shl0_713 = buffer.data(shl0 + 713);
    const auto *shl0_714 = buffer.data(shl0 + 714);
    const auto *shl0_715 = buffer.data(shl0 + 715);
    const auto *shl0_716 = buffer.data(shl0 + 716);
    const auto *shl0_717 = buffer.data(shl0 + 717);
    const auto *shl0_719 = buffer.data(shl0 + 719);
    const auto *shl0_725 = buffer.data(shl0 + 725);
    const auto *shl0_729 = buffer.data(shl0 + 729);
    const auto *shl0_732 = buffer.data(shl0 + 732);
    const auto *shl0_734 = buffer.data(shl0 + 734);
    const auto *shl0_737 = buffer.data(shl0 + 737);
    const auto *shl0_738 = buffer.data(shl0 + 738);
    const auto *shl0_740 = buffer.data(shl0 + 740);
    const auto *shl0_743 = buffer.data(shl0 + 743);
    const auto *shl0_744 = buffer.data(shl0 + 744);
    const auto *shl0_745 = buffer.data(shl0 + 745);
    const auto *shl0_747 = buffer.data(shl0 + 747);
    const auto *shl0_756 = buffer.data(shl0 + 756);
    const auto *shl0_758 = buffer.data(shl0 + 758);
    const auto *shl0_759 = buffer.data(shl0 + 759);
    const auto *shl0_760 = buffer.data(shl0 + 760);
    const auto *shl0_761 = buffer.data(shl0 + 761);
    const auto *shl0_762 = buffer.data(shl0 + 762);
    const auto *shl0_764 = buffer.data(shl0 + 764);
    const auto *shl0_765 = buffer.data(shl0 + 765);
    const auto *shl0_768 = buffer.data(shl0 + 768);
    const auto *shl0_770 = buffer.data(shl0 + 770);
    const auto *shl0_771 = buffer.data(shl0 + 771);
    const auto *shl0_774 = buffer.data(shl0 + 774);
    const auto *shl0_775 = buffer.data(shl0 + 775);
    const auto *shl0_777 = buffer.data(shl0 + 777);
    const auto *shl0_779 = buffer.data(shl0 + 779);
    const auto *shl0_780 = buffer.data(shl0 + 780);
    const auto *shl0_782 = buffer.data(shl0 + 782);
    const auto *shl0_783 = buffer.data(shl0 + 783);
    const auto *shl0_785 = buffer.data(shl0 + 785);
    const auto *shl0_786 = buffer.data(shl0 + 786);
    const auto *shl0_788 = buffer.data(shl0 + 788);
    const auto *shl0_789 = buffer.data(shl0 + 789);
    const auto *shl0_790 = buffer.data(shl0 + 790);

    const auto *shk_359 = buffer.data(shk + 359);
    const auto *shk_360 = buffer.data(shk + 360);
    const auto *shk_362 = buffer.data(shk + 362);
    const auto *shk_363 = buffer.data(shk + 363);
    const auto *shk_365 = buffer.data(shk + 365);
    const auto *shk_366 = buffer.data(shk + 366);
    const auto *shk_369 = buffer.data(shk + 369);
    const auto *shk_370 = buffer.data(shk + 370);
    const auto *shk_374 = buffer.data(shk + 374);
    const auto *shk_375 = buffer.data(shk + 375);
    const auto *shk_380 = buffer.data(shk + 380);
    const auto *shk_388 = buffer.data(shk + 388);
    const auto *shk_395 = buffer.data(shk + 395);
    const auto *shk_396 = buffer.data(shk + 396);
    const auto *shk_398 = buffer.data(shk + 398);
    const auto *shk_399 = buffer.data(shk + 399);
    const auto *shk_401 = buffer.data(shk + 401);
    const auto *shk_402 = buffer.data(shk + 402);
    const auto *shk_405 = buffer.data(shk + 405);
    const auto *shk_406 = buffer.data(shk + 406);
    const auto *shk_410 = buffer.data(shk + 410);
    const auto *shk_411 = buffer.data(shk + 411);
    const auto *shk_416 = buffer.data(shk + 416);
    const auto *shk_431 = buffer.data(shk + 431);
    const auto *shk_432 = buffer.data(shk + 432);
    const auto *shk_434 = buffer.data(shk + 434);
    const auto *shk_437 = buffer.data(shk + 437);
    const auto *shk_441 = buffer.data(shk + 441);
    const auto *shk_446 = buffer.data(shk + 446);
    const auto *shk_540 = buffer.data(shk + 540);
    const auto *shk_543 = buffer.data(shk + 543);
    const auto *shk_545 = buffer.data(shk + 545);
    const auto *shk_546 = buffer.data(shk + 546);
    const auto *shk_549 = buffer.data(shk + 549);
    const auto *shk_550 = buffer.data(shk + 550);
    const auto *shk_552 = buffer.data(shk + 552);
    const auto *shk_554 = buffer.data(shk + 554);
    const auto *shk_555 = buffer.data(shk + 555);
    const auto *shk_557 = buffer.data(shk + 557);
    const auto *shk_558 = buffer.data(shk + 558);
    const auto *shk_560 = buffer.data(shk + 560);
    const auto *shk_561 = buffer.data(shk + 561);
    const auto *shk_563 = buffer.data(shk + 563);
    const auto *shk_564 = buffer.data(shk + 564);
    const auto *shk_565 = buffer.data(shk + 565);
    const auto *shk_567 = buffer.data(shk + 567);
    const auto *shk_568 = buffer.data(shk + 568);
    const auto *shk_569 = buffer.data(shk + 569);
    const auto *shk_570 = buffer.data(shk + 570);
    const auto *shk_571 = buffer.data(shk + 571);
    const auto *shk_572 = buffer.data(shk + 572);
    const auto *shk_573 = buffer.data(shk + 573);
    const auto *shk_574 = buffer.data(shk + 574);
    const auto *shk_575 = buffer.data(shk + 575);
    const auto *shk_581 = buffer.data(shk + 581);
    const auto *shk_585 = buffer.data(shk + 585);
    const auto *shk_588 = buffer.data(shk + 588);
    const auto *shk_590 = buffer.data(shk + 590);
    const auto *shk_593 = buffer.data(shk + 593);
    const auto *shk_594 = buffer.data(shk + 594);
    const auto *shk_596 = buffer.data(shk + 596);
    const auto *shk_599 = buffer.data(shk + 599);
    const auto *shk_600 = buffer.data(shk + 600);
    const auto *shk_601 = buffer.data(shk + 601);
    const auto *shk_603 = buffer.data(shk + 603);
    const auto *shk_604 = buffer.data(shk + 604);
    const auto *shk_605 = buffer.data(shk + 605);
    const auto *shk_606 = buffer.data(shk + 606);
    const auto *shk_607 = buffer.data(shk + 607);
    const auto *shk_608 = buffer.data(shk + 608);
    const auto *shk_609 = buffer.data(shk + 609);
    const auto *shk_610 = buffer.data(shk + 610);
    const auto *shk_611 = buffer.data(shk + 611);
    const auto *shk_612 = buffer.data(shk + 612);
    const auto *shk_615 = buffer.data(shk + 615);
    const auto *shk_617 = buffer.data(shk + 617);
    const auto *shk_618 = buffer.data(shk + 618);
    const auto *shk_621 = buffer.data(shk + 621);
    const auto *shk_622 = buffer.data(shk + 622);
    const auto *shk_624 = buffer.data(shk + 624);
    const auto *shk_626 = buffer.data(shk + 626);
    const auto *shk_627 = buffer.data(shk + 627);
    const auto *shk_629 = buffer.data(shk + 629);
    const auto *shk_630 = buffer.data(shk + 630);
    const auto *shk_632 = buffer.data(shk + 632);
    const auto *shk_633 = buffer.data(shk + 633);
    const auto *shk_635 = buffer.data(shk + 635);
    const auto *shk_636 = buffer.data(shk + 636);
    const auto *shk_637 = buffer.data(shk + 637);

    const auto *shl1_450 = buffer.data(shl1 + 450);
    const auto *shl1_453 = buffer.data(shl1 + 453);
    const auto *shl1_456 = buffer.data(shl1 + 456);
    const auto *shl1_460 = buffer.data(shl1 + 460);
    const auto *shl1_465 = buffer.data(shl1 + 465);
    const auto *shl1_471 = buffer.data(shl1 + 471);
    const auto *shl1_675 = buffer.data(shl1 + 675);
    const auto *shl1_678 = buffer.data(shl1 + 678);
    const auto *shl1_680 = buffer.data(shl1 + 680);
    const auto *shl1_681 = buffer.data(shl1 + 681);
    const auto *shl1_684 = buffer.data(shl1 + 684);
    const auto *shl1_685 = buffer.data(shl1 + 685);
    const auto *shl1_687 = buffer.data(shl1 + 687);
    const auto *shl1_689 = buffer.data(shl1 + 689);
    const auto *shl1_690 = buffer.data(shl1 + 690);
    const auto *shl1_692 = buffer.data(shl1 + 692);
    const auto *shl1_693 = buffer.data(shl1 + 693);
    const auto *shl1_695 = buffer.data(shl1 + 695);
    const auto *shl1_696 = buffer.data(shl1 + 696);
    const auto *shl1_698 = buffer.data(shl1 + 698);
    const auto *shl1_699 = buffer.data(shl1 + 699);
    const auto *shl1_700 = buffer.data(shl1 + 700);
    const auto *shl1_702 = buffer.data(shl1 + 702);
    const auto *shl1_711 = buffer.data(shl1 + 711);
    const auto *shl1_713 = buffer.data(shl1 + 713);
    const auto *shl1_714 = buffer.data(shl1 + 714);
    const auto *shl1_715 = buffer.data(shl1 + 715);
    const auto *shl1_716 = buffer.data(shl1 + 716);
    const auto *shl1_717 = buffer.data(shl1 + 717);
    const auto *shl1_719 = buffer.data(shl1 + 719);
    const auto *shl1_725 = buffer.data(shl1 + 725);
    const auto *shl1_729 = buffer.data(shl1 + 729);
    const auto *shl1_732 = buffer.data(shl1 + 732);
    const auto *shl1_734 = buffer.data(shl1 + 734);
    const auto *shl1_737 = buffer.data(shl1 + 737);
    const auto *shl1_738 = buffer.data(shl1 + 738);
    const auto *shl1_740 = buffer.data(shl1 + 740);
    const auto *shl1_743 = buffer.data(shl1 + 743);
    const auto *shl1_744 = buffer.data(shl1 + 744);
    const auto *shl1_745 = buffer.data(shl1 + 745);
    const auto *shl1_747 = buffer.data(shl1 + 747);
    const auto *shl1_756 = buffer.data(shl1 + 756);
    const auto *shl1_758 = buffer.data(shl1 + 758);
    const auto *shl1_759 = buffer.data(shl1 + 759);
    const auto *shl1_760 = buffer.data(shl1 + 760);
    const auto *shl1_761 = buffer.data(shl1 + 761);
    const auto *shl1_762 = buffer.data(shl1 + 762);
    const auto *shl1_764 = buffer.data(shl1 + 764);
    const auto *shl1_765 = buffer.data(shl1 + 765);
    const auto *shl1_768 = buffer.data(shl1 + 768);
    const auto *shl1_770 = buffer.data(shl1 + 770);
    const auto *shl1_771 = buffer.data(shl1 + 771);
    const auto *shl1_774 = buffer.data(shl1 + 774);
    const auto *shl1_775 = buffer.data(shl1 + 775);
    const auto *shl1_777 = buffer.data(shl1 + 777);
    const auto *shl1_779 = buffer.data(shl1 + 779);
    const auto *shl1_780 = buffer.data(shl1 + 780);
    const auto *shl1_782 = buffer.data(shl1 + 782);
    const auto *shl1_783 = buffer.data(shl1 + 783);
    const auto *shl1_785 = buffer.data(shl1 + 785);
    const auto *shl1_786 = buffer.data(shl1 + 786);
    const auto *shl1_788 = buffer.data(shl1 + 788);
    const auto *shl1_789 = buffer.data(shl1 + 789);
    const auto *shl1_790 = buffer.data(shl1 + 790);

    const auto *sii0_419 = buffer.data(sii0 + 419);

    const auto *sii1_419 = buffer.data(sii1 + 419);

    const auto *sik_539 = buffer.data(sik + 539);
    const auto *sik_540 = buffer.data(sik + 540);
    const auto *sik_542 = buffer.data(sik + 542);
    const auto *sik_543 = buffer.data(sik + 543);
    const auto *sik_545 = buffer.data(sik + 545);
    const auto *sik_546 = buffer.data(sik + 546);
    const auto *sik_549 = buffer.data(sik + 549);
    const auto *sik_550 = buffer.data(sik + 550);
    const auto *sik_554 = buffer.data(sik + 554);
    const auto *sik_555 = buffer.data(sik + 555);
    const auto *sik_560 = buffer.data(sik + 560);
    const auto *sik_568 = buffer.data(sik + 568);
    const auto *sik_569 = buffer.data(sik + 569);
    const auto *sik_570 = buffer.data(sik + 570);
    const auto *sik_571 = buffer.data(sik + 571);
    const auto *sik_572 = buffer.data(sik + 572);
    const auto *sik_573 = buffer.data(sik + 573);
    const auto *sik_574 = buffer.data(sik + 574);
    const auto *sik_575 = buffer.data(sik + 575);
    const auto *sik_576 = buffer.data(sik + 576);
    const auto *sik_578 = buffer.data(sik + 578);
    const auto *sik_579 = buffer.data(sik + 579);
    const auto *sik_581 = buffer.data(sik + 581);
    const auto *sik_582 = buffer.data(sik + 582);
    const auto *sik_585 = buffer.data(sik + 585);
    const auto *sik_586 = buffer.data(sik + 586);
    const auto *sik_590 = buffer.data(sik + 590);
    const auto *sik_591 = buffer.data(sik + 591);
    const auto *sik_596 = buffer.data(sik + 596);
    const auto *sik_604 = buffer.data(sik + 604);
    const auto *sik_605 = buffer.data(sik + 605);
    const auto *sik_606 = buffer.data(sik + 606);
    const auto *sik_607 = buffer.data(sik + 607);
    const auto *sik_608 = buffer.data(sik + 608);
    const auto *sik_609 = buffer.data(sik + 609);
    const auto *sik_610 = buffer.data(sik + 610);
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

#pragma omp simd aligned(t_674, t_675, t_676, pb_x, pc_x, pc_y, pc_z, shl0_675, shk_359, \
                         shk_360, shk_540, shl1_675, sii0_419, sii1_419, sik_539, \
                         sik_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_18 * shk_359[k]
                   + f_1 * sii0_419[k]
                   - f_2 * sii1_419[k]
                   + f_3 * pc_z[k] * sik_539[k];

        t_675[k] = pb_x[k] * shl0_675[k]
                   + f_20 * shk_540[k]
                   - f_14 * pc_x[k] * shl1_675[k];

        t_676[k] = f_19 * shk_360[k]
                   + f_3 * pc_y[k] * sik_540[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pb_x, pc_x, pc_y, pc_z, shl0_678, shk_362, \
                         shk_543, shl1_678, sik_540, sik_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_3 * pc_z[k] * sik_540[k];

        t_678[k] = pb_x[k] * shl0_678[k]
                   + f_0 * shk_543[k]
                   - f_14 * pc_x[k] * shl1_678[k];

        t_679[k] = f_19 * shk_362[k]
                   + f_3 * pc_y[k] * sik_542[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, pb_x, pc_x, pc_z, shl0_680, shl0_681, shk_545, \
                         shk_546, shl1_680, shl1_681, sik_543 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = pb_x[k] * shl0_680[k]
                   + f_0 * shk_545[k]
                   - f_14 * pc_x[k] * shl1_680[k];

        t_681[k] = pb_x[k] * shl0_681[k]
                   + f_19 * shk_546[k]
                   - f_14 * pc_x[k] * shl1_681[k];

        t_682[k] = f_3 * pc_z[k] * sik_543[k];
    }

#pragma omp simd aligned(t_683, t_684, t_685, pb_x, pc_x, pc_y, shl0_684, shl0_685, shk_365, \
                         shk_549, shk_550, shl1_684, shl1_685, \
                         sik_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_683[k] = f_19 * shk_365[k]
                   + f_3 * pc_y[k] * sik_545[k];

        t_684[k] = pb_x[k] * shl0_684[k]
                   + f_19 * shk_549[k]
                   - f_14 * pc_x[k] * shl1_684[k];

        t_685[k] = pb_x[k] * shl0_685[k]
                   + f_18 * shk_550[k]
                   - f_14 * pc_x[k] * shl1_685[k];
    }

#pragma omp simd aligned(t_686, t_687, t_688, pb_x, pc_x, pc_y, pc_z, shl0_687, shk_369, \
                         shk_552, shl1_687, sik_546, sik_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_686[k] = f_3 * pc_z[k] * sik_546[k];

        t_687[k] = pb_x[k] * shl0_687[k]
                   + f_18 * shk_552[k]
                   - f_14 * pc_x[k] * shl1_687[k];

        t_688[k] = f_19 * shk_369[k]
                   + f_3 * pc_y[k] * sik_549[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pb_x, pc_x, pc_z, shl0_689, shl0_690, shk_554, \
                         shk_555, shl1_689, shl1_690, sik_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = pb_x[k] * shl0_689[k]
                   + f_18 * shk_554[k]
                   - f_14 * pc_x[k] * shl1_689[k];

        t_690[k] = pb_x[k] * shl0_690[k]
                   + f_17 * shk_555[k]
                   - f_14 * pc_x[k] * shl1_690[k];

        t_691[k] = f_3 * pc_z[k] * sik_550[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, pb_x, pc_x, pc_y, shl0_692, shl0_693, shk_374, \
                         shk_557, shk_558, shl1_692, shl1_693, \
                         sik_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = pb_x[k] * shl0_692[k]
                   + f_17 * shk_557[k]
                   - f_14 * pc_x[k] * shl1_692[k];

        t_693[k] = pb_x[k] * shl0_693[k]
                   + f_17 * shk_558[k]
                   - f_14 * pc_x[k] * shl1_693[k];

        t_694[k] = f_19 * shk_374[k]
                   + f_3 * pc_y[k] * sik_554[k];
    }

#pragma omp simd aligned(t_695, t_696, t_697, pb_x, pc_x, pc_z, shl0_695, shl0_696, shk_560, \
                         shk_561, shl1_695, shl1_696, sik_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_695[k] = pb_x[k] * shl0_695[k]
                   + f_17 * shk_560[k]
                   - f_14 * pc_x[k] * shl1_695[k];

        t_696[k] = pb_x[k] * shl0_696[k]
                   + f_16 * shk_561[k]
                   - f_14 * pc_x[k] * shl1_696[k];

        t_697[k] = f_3 * pc_z[k] * sik_555[k];
    }

#pragma omp simd aligned(t_698, t_699, t_700, pb_x, pc_x, shl0_698, shl0_699, shl0_700, \
                         shk_563, shk_564, shk_565, shl1_698, shl1_699, \
                         shl1_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_698[k] = pb_x[k] * shl0_698[k]
                   + f_16 * shk_563[k]
                   - f_14 * pc_x[k] * shl1_698[k];

        t_699[k] = pb_x[k] * shl0_699[k]
                   + f_16 * shk_564[k]
                   - f_14 * pc_x[k] * shl1_699[k];

        t_700[k] = pb_x[k] * shl0_700[k]
                   + f_16 * shk_565[k]
                   - f_14 * pc_x[k] * shl1_700[k];
    }

#pragma omp simd aligned(t_701, t_702, t_703, t_704, pb_x, pc_x, pc_y, shl0_702, shk_380, \
                         shk_567, shk_568, shk_569, shl1_702, sik_560, sik_568, \
                         sik_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_701[k] = f_19 * shk_380[k]
                   + f_3 * pc_y[k] * sik_560[k];

        t_702[k] = pb_x[k] * shl0_702[k]
                   + f_16 * shk_567[k]
                   - f_14 * pc_x[k] * shl1_702[k];

        t_703[k] = f_15 * shk_568[k]
                   + f_3 * pc_x[k] * sik_568[k];

        t_704[k] = f_15 * shk_569[k]
                   + f_3 * pc_x[k] * sik_569[k];
    }

#pragma omp simd aligned(t_705, t_706, t_707, t_708, t_709, pc_x, shk_570, shk_571, shk_572, \
                         shk_573, shk_574, sik_570, sik_571, sik_572, sik_573, \
                         sik_574 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_705[k] = f_15 * shk_570[k]
                   + f_3 * pc_x[k] * sik_570[k];

        t_706[k] = f_15 * shk_571[k]
                   + f_3 * pc_x[k] * sik_571[k];

        t_707[k] = f_15 * shk_572[k]
                   + f_3 * pc_x[k] * sik_572[k];

        t_708[k] = f_15 * shk_573[k]
                   + f_3 * pc_x[k] * sik_573[k];

        t_709[k] = f_15 * shk_574[k]
                   + f_3 * pc_x[k] * sik_574[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pb_x, pc_x, pc_z, shl0_711, shl0_713, \
                         shk_575, shl1_711, shl1_713, sik_568, \
                         sik_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_15 * shk_575[k]
                   + f_3 * pc_x[k] * sik_575[k];

        t_711[k] = pb_x[k] * shl0_711[k]
                   - f_14 * pc_x[k] * shl1_711[k];

        t_712[k] = f_3 * pc_z[k] * sik_568[k];

        t_713[k] = pb_x[k] * shl0_713[k]
                   - f_14 * pc_x[k] * shl1_713[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, t_717, pb_x, pc_x, shl0_714, shl0_715, shl0_716, \
                         shl0_717, shl1_714, shl1_715, shl1_716, \
                         shl1_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = pb_x[k] * shl0_714[k]
                   - f_14 * pc_x[k] * shl1_714[k];

        t_715[k] = pb_x[k] * shl0_715[k]
                   - f_14 * pc_x[k] * shl1_715[k];

        t_716[k] = pb_x[k] * shl0_716[k]
                   - f_14 * pc_x[k] * shl1_716[k];

        t_717[k] = pb_x[k] * shl0_717[k]
                   - f_14 * pc_x[k] * shl1_717[k];
    }

#pragma omp simd aligned(t_718, t_719, t_720, pb_x, pb_z, pc_x, pc_y, pc_z, shl0_450, \
                         shl0_719, shk_395, shl1_450, shl1_719, \
                         sik_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_19 * shk_395[k]
                   + f_3 * pc_y[k] * sik_575[k];

        t_719[k] = pb_x[k] * shl0_719[k]
                   - f_14 * pc_x[k] * shl1_719[k];

        t_720[k] = pb_z[k] * shl0_450[k]
                   - f_14 * pc_z[k] * shl1_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pb_z, pc_y, pc_z, shl0_453, shk_360, \
                         shk_396, shk_398, shl1_453, sik_576, sik_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_18 * shk_396[k]
                   + f_3 * pc_y[k] * sik_576[k];

        t_722[k] = f_15 * shk_360[k]
                   + f_3 * pc_z[k] * sik_576[k];

        t_723[k] = pb_z[k] * shl0_453[k]
                   - f_14 * pc_z[k] * shl1_453[k];

        t_724[k] = f_18 * shk_398[k]
                   + f_3 * pc_y[k] * sik_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pb_x, pb_z, pc_x, pc_z, shl0_456, shl0_725, \
                         shk_363, shk_581, shl1_456, shl1_725, \
                         sik_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = pb_x[k] * shl0_725[k]
                   + f_0 * shk_581[k]
                   - f_14 * pc_x[k] * shl1_725[k];

        t_726[k] = pb_z[k] * shl0_456[k]
                   - f_14 * pc_z[k] * shl1_456[k];

        t_727[k] = f_15 * shk_363[k]
                   + f_3 * pc_z[k] * sik_579[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pb_x, pb_z, pc_x, pc_y, pc_z, shl0_460, \
                         shl0_729, shk_401, shk_585, shl1_460, shl1_729, \
                         sik_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_18 * shk_401[k]
                   + f_3 * pc_y[k] * sik_581[k];

        t_729[k] = pb_x[k] * shl0_729[k]
                   + f_19 * shk_585[k]
                   - f_14 * pc_x[k] * shl1_729[k];

        t_730[k] = pb_z[k] * shl0_460[k]
                   - f_14 * pc_z[k] * shl1_460[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_x, pc_x, pc_y, pc_z, shl0_732, shk_366, \
                         shk_405, shk_588, shl1_732, sik_582, sik_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_15 * shk_366[k]
                   + f_3 * pc_z[k] * sik_582[k];

        t_732[k] = pb_x[k] * shl0_732[k]
                   + f_18 * shk_588[k]
                   - f_14 * pc_x[k] * shl1_732[k];

        t_733[k] = f_18 * shk_405[k]
                   + f_3 * pc_y[k] * sik_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_x, pb_z, pc_x, pc_z, shl0_465, shl0_734, \
                         shk_370, shk_590, shl1_465, shl1_734, \
                         sik_586 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = pb_x[k] * shl0_734[k]
                   + f_18 * shk_590[k]
                   - f_14 * pc_x[k] * shl1_734[k];

        t_735[k] = pb_z[k] * shl0_465[k]
                   - f_14 * pc_z[k] * shl1_465[k];

        t_736[k] = f_15 * shk_370[k]
                   + f_3 * pc_z[k] * sik_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_x, pc_x, pc_y, shl0_737, shl0_738, shk_410, \
                         shk_593, shk_594, shl1_737, shl1_738, \
                         sik_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_x[k] * shl0_737[k]
                   + f_17 * shk_593[k]
                   - f_14 * pc_x[k] * shl1_737[k];

        t_738[k] = pb_x[k] * shl0_738[k]
                   + f_17 * shk_594[k]
                   - f_14 * pc_x[k] * shl1_738[k];

        t_739[k] = f_18 * shk_410[k]
                   + f_3 * pc_y[k] * sik_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pb_x, pb_z, pc_x, pc_z, shl0_471, shl0_740, \
                         shk_375, shk_596, shl1_471, shl1_740, \
                         sik_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = pb_x[k] * shl0_740[k]
                   + f_17 * shk_596[k]
                   - f_14 * pc_x[k] * shl1_740[k];

        t_741[k] = pb_z[k] * shl0_471[k]
                   - f_14 * pc_z[k] * shl1_471[k];

        t_742[k] = f_15 * shk_375[k]
                   + f_3 * pc_z[k] * sik_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pb_x, pc_x, shl0_743, shl0_744, shl0_745, \
                         shk_599, shk_600, shk_601, shl1_743, shl1_744, \
                         shl1_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pb_x[k] * shl0_743[k]
                   + f_16 * shk_599[k]
                   - f_14 * pc_x[k] * shl1_743[k];

        t_744[k] = pb_x[k] * shl0_744[k]
                   + f_16 * shk_600[k]
                   - f_14 * pc_x[k] * shl1_744[k];

        t_745[k] = pb_x[k] * shl0_745[k]
                   + f_16 * shk_601[k]
                   - f_14 * pc_x[k] * shl1_745[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pb_x, pc_x, pc_y, shl0_747, shk_416, \
                         shk_603, shk_604, shk_605, shl1_747, sik_596, sik_604, \
                         sik_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * shk_416[k]
                   + f_3 * pc_y[k] * sik_596[k];

        t_747[k] = pb_x[k] * shl0_747[k]
                   + f_16 * shk_603[k]
                   - f_14 * pc_x[k] * shl1_747[k];

        t_748[k] = f_15 * shk_604[k]
                   + f_3 * pc_x[k] * sik_604[k];

        t_749[k] = f_15 * shk_605[k]
                   + f_3 * pc_x[k] * sik_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, shk_606, shk_607, shk_608, \
                         shk_609, shk_610, sik_606, sik_607, sik_608, sik_609, \
                         sik_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_15 * shk_606[k]
                   + f_3 * pc_x[k] * sik_606[k];

        t_751[k] = f_15 * shk_607[k]
                   + f_3 * pc_x[k] * sik_607[k];

        t_752[k] = f_15 * shk_608[k]
                   + f_3 * pc_x[k] * sik_608[k];

        t_753[k] = f_15 * shk_609[k]
                   + f_3 * pc_x[k] * sik_609[k];

        t_754[k] = f_15 * shk_610[k]
                   + f_3 * pc_x[k] * sik_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, t_758, pb_x, pc_x, pc_z, shl0_756, shl0_758, \
                         shk_388, shk_611, shl1_756, shl1_758, sik_604, \
                         sik_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_15 * shk_611[k]
                   + f_3 * pc_x[k] * sik_611[k];

        t_756[k] = pb_x[k] * shl0_756[k]
                   - f_14 * pc_x[k] * shl1_756[k];

        t_757[k] = f_15 * shk_388[k]
                   + f_3 * pc_z[k] * sik_604[k];

        t_758[k] = pb_x[k] * shl0_758[k]
                   - f_14 * pc_x[k] * shl1_758[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, t_762, pb_x, pc_x, shl0_759, shl0_760, shl0_761, \
                         shl0_762, shl1_759, shl1_760, shl1_761, \
                         shl1_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = pb_x[k] * shl0_759[k]
                   - f_14 * pc_x[k] * shl1_759[k];

        t_760[k] = pb_x[k] * shl0_760[k]
                   - f_14 * pc_x[k] * shl1_760[k];

        t_761[k] = pb_x[k] * shl0_761[k]
                   - f_14 * pc_x[k] * shl1_761[k];

        t_762[k] = pb_x[k] * shl0_762[k]
                   - f_14 * pc_x[k] * shl1_762[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, t_766, pb_x, pc_x, pc_y, shl0_764, shl0_765, \
                         shk_431, shk_432, shk_612, shl1_764, shl1_765, sik_611, \
                         sik_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_18 * shk_431[k]
                   + f_3 * pc_y[k] * sik_611[k];

        t_764[k] = pb_x[k] * shl0_764[k]
                   - f_14 * pc_x[k] * shl1_764[k];

        t_765[k] = pb_x[k] * shl0_765[k]
                   + f_20 * shk_612[k]
                   - f_14 * pc_x[k] * shl1_765[k];

        t_766[k] = f_17 * shk_432[k]
                   + f_3 * pc_y[k] * sik_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pb_x, pc_x, pc_y, pc_z, shl0_768, shk_396, \
                         shk_434, shk_615, shl1_768, sik_612, sik_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * shk_396[k]
                   + f_3 * pc_z[k] * sik_612[k];

        t_768[k] = pb_x[k] * shl0_768[k]
                   + f_0 * shk_615[k]
                   - f_14 * pc_x[k] * shl1_768[k];

        t_769[k] = f_17 * shk_434[k]
                   + f_3 * pc_y[k] * sik_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pb_x, pc_x, pc_z, shl0_770, shl0_771, shk_399, \
                         shk_617, shk_618, shl1_770, shl1_771, \
                         sik_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = pb_x[k] * shl0_770[k]
                   + f_0 * shk_617[k]
                   - f_14 * pc_x[k] * shl1_770[k];

        t_771[k] = pb_x[k] * shl0_771[k]
                   + f_19 * shk_618[k]
                   - f_14 * pc_x[k] * shl1_771[k];

        t_772[k] = f_16 * shk_399[k]
                   + f_3 * pc_z[k] * sik_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pb_x, pc_x, pc_y, shl0_774, shl0_775, shk_437, \
                         shk_621, shk_622, shl1_774, shl1_775, \
                         sik_617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * shk_437[k]
                   + f_3 * pc_y[k] * sik_617[k];

        t_774[k] = pb_x[k] * shl0_774[k]
                   + f_19 * shk_621[k]
                   - f_14 * pc_x[k] * shl1_774[k];

        t_775[k] = pb_x[k] * shl0_775[k]
                   + f_18 * shk_622[k]
                   - f_14 * pc_x[k] * shl1_775[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pb_x, pc_x, pc_y, pc_z, shl0_777, shk_402, \
                         shk_441, shk_624, shl1_777, sik_618, sik_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * shk_402[k]
                   + f_3 * pc_z[k] * sik_618[k];

        t_777[k] = pb_x[k] * shl0_777[k]
                   + f_18 * shk_624[k]
                   - f_14 * pc_x[k] * shl1_777[k];

        t_778[k] = f_17 * shk_441[k]
                   + f_3 * pc_y[k] * sik_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pb_x, pc_x, pc_z, shl0_779, shl0_780, shk_406, \
                         shk_626, shk_627, shl1_779, shl1_780, \
                         sik_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = pb_x[k] * shl0_779[k]
                   + f_18 * shk_626[k]
                   - f_14 * pc_x[k] * shl1_779[k];

        t_780[k] = pb_x[k] * shl0_780[k]
                   + f_17 * shk_627[k]
                   - f_14 * pc_x[k] * shl1_780[k];

        t_781[k] = f_16 * shk_406[k]
                   + f_3 * pc_z[k] * sik_622[k];
    }

#pragma omp simd aligned(t_782, t_783, t_784, pb_x, pc_x, pc_y, shl0_782, shl0_783, shk_446, \
                         shk_629, shk_630, shl1_782, shl1_783, \
                         sik_626 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = pb_x[k] * shl0_782[k]
                   + f_17 * shk_629[k]
                   - f_14 * pc_x[k] * shl1_782[k];

        t_783[k] = pb_x[k] * shl0_783[k]
                   + f_17 * shk_630[k]
                   - f_14 * pc_x[k] * shl1_783[k];

        t_784[k] = f_17 * shk_446[k]
                   + f_3 * pc_y[k] * sik_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pb_x, pc_x, pc_z, shl0_785, shl0_786, shk_411, \
                         shk_632, shk_633, shl1_785, shl1_786, \
                         sik_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = pb_x[k] * shl0_785[k]
                   + f_17 * shk_632[k]
                   - f_14 * pc_x[k] * shl1_785[k];

        t_786[k] = pb_x[k] * shl0_786[k]
                   + f_16 * shk_633[k]
                   - f_14 * pc_x[k] * shl1_786[k];

        t_787[k] = f_16 * shk_411[k]
                   + f_3 * pc_z[k] * sik_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pb_x, pc_x, shl0_788, shl0_789, shl0_790, \
                         shk_635, shk_636, shk_637, shl1_788, shl1_789, \
                         shl1_790 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = pb_x[k] * shl0_788[k]
                   + f_16 * shk_635[k]
                   - f_14 * pc_x[k] * shl1_788[k];

        t_789[k] = pb_x[k] * shl0_789[k]
                   + f_16 * shk_636[k]
                   - f_14 * pc_x[k] * shl1_789[k];

        t_790[k] = pb_x[k] * shl0_790[k]
                   + f_16 * shk_637[k]
                   - f_14 * pc_x[k] * shl1_790[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 4.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shl0_630 = buffer.data(shl0 + 630);
    const auto *shl0_635 = buffer.data(shl0 + 635);
    const auto *shl0_639 = buffer.data(shl0 + 639);
    const auto *shl0_644 = buffer.data(shl0 + 644);
    const auto *shl0_650 = buffer.data(shl0 + 650);
    const auto *shl0_657 = buffer.data(shl0 + 657);
    const auto *shl0_792 = buffer.data(shl0 + 792);
    const auto *shl0_801 = buffer.data(shl0 + 801);
    const auto *shl0_803 = buffer.data(shl0 + 803);
    const auto *shl0_804 = buffer.data(shl0 + 804);
    const auto *shl0_805 = buffer.data(shl0 + 805);
    const auto *shl0_806 = buffer.data(shl0 + 806);
    const auto *shl0_807 = buffer.data(shl0 + 807);
    const auto *shl0_809 = buffer.data(shl0 + 809);
    const auto *shl0_810 = buffer.data(shl0 + 810);
    const auto *shl0_813 = buffer.data(shl0 + 813);
    const auto *shl0_815 = buffer.data(shl0 + 815);
    const auto *shl0_816 = buffer.data(shl0 + 816);
    const auto *shl0_819 = buffer.data(shl0 + 819);
    const auto *shl0_820 = buffer.data(shl0 + 820);
    const auto *shl0_822 = buffer.data(shl0 + 822);
    const auto *shl0_824 = buffer.data(shl0 + 824);
    const auto *shl0_825 = buffer.data(shl0 + 825);
    const auto *shl0_827 = buffer.data(shl0 + 827);
    const auto *shl0_828 = buffer.data(shl0 + 828);
    const auto *shl0_830 = buffer.data(shl0 + 830);
    const auto *shl0_831 = buffer.data(shl0 + 831);
    const auto *shl0_833 = buffer.data(shl0 + 833);
    const auto *shl0_834 = buffer.data(shl0 + 834);
    const auto *shl0_835 = buffer.data(shl0 + 835);
    const auto *shl0_837 = buffer.data(shl0 + 837);
    const auto *shl0_846 = buffer.data(shl0 + 846);
    const auto *shl0_848 = buffer.data(shl0 + 848);
    const auto *shl0_849 = buffer.data(shl0 + 849);
    const auto *shl0_850 = buffer.data(shl0 + 850);
    const auto *shl0_851 = buffer.data(shl0 + 851);
    const auto *shl0_852 = buffer.data(shl0 + 852);
    const auto *shl0_854 = buffer.data(shl0 + 854);
    const auto *shl0_858 = buffer.data(shl0 + 858);
    const auto *shl0_861 = buffer.data(shl0 + 861);
    const auto *shl0_865 = buffer.data(shl0 + 865);
    const auto *shl0_867 = buffer.data(shl0 + 867);
    const auto *shl0_870 = buffer.data(shl0 + 870);
    const auto *shl0_872 = buffer.data(shl0 + 872);
    const auto *shl0_873 = buffer.data(shl0 + 873);
    const auto *shl0_876 = buffer.data(shl0 + 876);
    const auto *shl0_878 = buffer.data(shl0 + 878);
    const auto *shl0_879 = buffer.data(shl0 + 879);
    const auto *shl0_880 = buffer.data(shl0 + 880);
    const auto *shl0_891 = buffer.data(shl0 + 891);
    const auto *shl0_893 = buffer.data(shl0 + 893);
    const auto *shl0_894 = buffer.data(shl0 + 894);
    const auto *shl0_895 = buffer.data(shl0 + 895);
    const auto *shl0_896 = buffer.data(shl0 + 896);
    const auto *shl0_897 = buffer.data(shl0 + 897);
    const auto *shl0_899 = buffer.data(shl0 + 899);
    const auto *shl0_900 = buffer.data(shl0 + 900);
    const auto *shl0_903 = buffer.data(shl0 + 903);
    const auto *shl0_905 = buffer.data(shl0 + 905);
    const auto *shl0_906 = buffer.data(shl0 + 906);
    const auto *shl0_909 = buffer.data(shl0 + 909);
    const auto *shl0_910 = buffer.data(shl0 + 910);

    const auto *shk_424 = buffer.data(shk + 424);
    const auto *shk_432 = buffer.data(shk + 432);
    const auto *shk_435 = buffer.data(shk + 435);
    const auto *shk_438 = buffer.data(shk + 438);
    const auto *shk_442 = buffer.data(shk + 442);
    const auto *shk_447 = buffer.data(shk + 447);
    const auto *shk_452 = buffer.data(shk + 452);
    const auto *shk_460 = buffer.data(shk + 460);
    const auto *shk_467 = buffer.data(shk + 467);
    const auto *shk_468 = buffer.data(shk + 468);
    const auto *shk_470 = buffer.data(shk + 470);
    const auto *shk_471 = buffer.data(shk + 471);
    const auto *shk_473 = buffer.data(shk + 473);
    const auto *shk_474 = buffer.data(shk + 474);
    const auto *shk_477 = buffer.data(shk + 477);
    const auto *shk_478 = buffer.data(shk + 478);
    const auto *shk_482 = buffer.data(shk + 482);
    const auto *shk_483 = buffer.data(shk + 483);
    const auto *shk_488 = buffer.data(shk + 488);
    const auto *shk_496 = buffer.data(shk + 496);
    const auto *shk_503 = buffer.data(shk + 503);
    const auto *shk_504 = buffer.data(shk + 504);
    const auto *shk_506 = buffer.data(shk + 506);
    const auto *shk_507 = buffer.data(shk + 507);
    const auto *shk_509 = buffer.data(shk + 509);
    const auto *shk_513 = buffer.data(shk + 513);
    const auto *shk_518 = buffer.data(shk + 518);
    const auto *shk_524 = buffer.data(shk + 524);
    const auto *shk_539 = buffer.data(shk + 539);
    const auto *shk_639 = buffer.data(shk + 639);
    const auto *shk_640 = buffer.data(shk + 640);
    const auto *shk_641 = buffer.data(shk + 641);
    const auto *shk_642 = buffer.data(shk + 642);
    const auto *shk_643 = buffer.data(shk + 643);
    const auto *shk_644 = buffer.data(shk + 644);
    const auto *shk_645 = buffer.data(shk + 645);
    const auto *shk_646 = buffer.data(shk + 646);
    const auto *shk_647 = buffer.data(shk + 647);
    const auto *shk_648 = buffer.data(shk + 648);
    const auto *shk_651 = buffer.data(shk + 651);
    const auto *shk_653 = buffer.data(shk + 653);
    const auto *shk_654 = buffer.data(shk + 654);
    const auto *shk_657 = buffer.data(shk + 657);
    const auto *shk_658 = buffer.data(shk + 658);
    const auto *shk_660 = buffer.data(shk + 660);
    const auto *shk_662 = buffer.data(shk + 662);
    const auto *shk_663 = buffer.data(shk + 663);
    const auto *shk_665 = buffer.data(shk + 665);
    const auto *shk_666 = buffer.data(shk + 666);
    const auto *shk_668 = buffer.data(shk + 668);
    const auto *shk_669 = buffer.data(shk + 669);
    const auto *shk_671 = buffer.data(shk + 671);
    const auto *shk_672 = buffer.data(shk + 672);
    const auto *shk_673 = buffer.data(shk + 673);
    const auto *shk_675 = buffer.data(shk + 675);
    const auto *shk_676 = buffer.data(shk + 676);
    const auto *shk_677 = buffer.data(shk + 677);
    const auto *shk_678 = buffer.data(shk + 678);
    const auto *shk_679 = buffer.data(shk + 679);
    const auto *shk_680 = buffer.data(shk + 680);
    const auto *shk_681 = buffer.data(shk + 681);
    const auto *shk_682 = buffer.data(shk + 682);
    const auto *shk_683 = buffer.data(shk + 683);
    const auto *shk_687 = buffer.data(shk + 687);
    const auto *shk_690 = buffer.data(shk + 690);
    const auto *shk_694 = buffer.data(shk + 694);
    const auto *shk_696 = buffer.data(shk + 696);
    const auto *shk_699 = buffer.data(shk + 699);
    const auto *shk_701 = buffer.data(shk + 701);
    const auto *shk_702 = buffer.data(shk + 702);
    const auto *shk_705 = buffer.data(shk + 705);
    const auto *shk_707 = buffer.data(shk + 707);
    const auto *shk_708 = buffer.data(shk + 708);
    const auto *shk_709 = buffer.data(shk + 709);
    const auto *shk_712 = buffer.data(shk + 712);
    const auto *shk_713 = buffer.data(shk + 713);
    const auto *shk_714 = buffer.data(shk + 714);
    const auto *shk_715 = buffer.data(shk + 715);
    const auto *shk_716 = buffer.data(shk + 716);
    const auto *shk_717 = buffer.data(shk + 717);
    const auto *shk_718 = buffer.data(shk + 718);
    const auto *shk_719 = buffer.data(shk + 719);
    const auto *shk_720 = buffer.data(shk + 720);
    const auto *shk_723 = buffer.data(shk + 723);
    const auto *shk_725 = buffer.data(shk + 725);
    const auto *shk_726 = buffer.data(shk + 726);
    const auto *shk_729 = buffer.data(shk + 729);
    const auto *shk_730 = buffer.data(shk + 730);

    const auto *shl1_630 = buffer.data(shl1 + 630);
    const auto *shl1_635 = buffer.data(shl1 + 635);
    const auto *shl1_639 = buffer.data(shl1 + 639);
    const auto *shl1_644 = buffer.data(shl1 + 644);
    const auto *shl1_650 = buffer.data(shl1 + 650);
    const auto *shl1_657 = buffer.data(shl1 + 657);
    const auto *shl1_792 = buffer.data(shl1 + 792);
    const auto *shl1_801 = buffer.data(shl1 + 801);
    const auto *shl1_803 = buffer.data(shl1 + 803);
    const auto *shl1_804 = buffer.data(shl1 + 804);
    const auto *shl1_805 = buffer.data(shl1 + 805);
    const auto *shl1_806 = buffer.data(shl1 + 806);
    const auto *shl1_807 = buffer.data(shl1 + 807);
    const auto *shl1_809 = buffer.data(shl1 + 809);
    const auto *shl1_810 = buffer.data(shl1 + 810);
    const auto *shl1_813 = buffer.data(shl1 + 813);
    const auto *shl1_815 = buffer.data(shl1 + 815);
    const auto *shl1_816 = buffer.data(shl1 + 816);
    const auto *shl1_819 = buffer.data(shl1 + 819);
    const auto *shl1_820 = buffer.data(shl1 + 820);
    const auto *shl1_822 = buffer.data(shl1 + 822);
    const auto *shl1_824 = buffer.data(shl1 + 824);
    const auto *shl1_825 = buffer.data(shl1 + 825);
    const auto *shl1_827 = buffer.data(shl1 + 827);
    const auto *shl1_828 = buffer.data(shl1 + 828);
    const auto *shl1_830 = buffer.data(shl1 + 830);
    const auto *shl1_831 = buffer.data(shl1 + 831);
    const auto *shl1_833 = buffer.data(shl1 + 833);
    const auto *shl1_834 = buffer.data(shl1 + 834);
    const auto *shl1_835 = buffer.data(shl1 + 835);
    const auto *shl1_837 = buffer.data(shl1 + 837);
    const auto *shl1_846 = buffer.data(shl1 + 846);
    const auto *shl1_848 = buffer.data(shl1 + 848);
    const auto *shl1_849 = buffer.data(shl1 + 849);
    const auto *shl1_850 = buffer.data(shl1 + 850);
    const auto *shl1_851 = buffer.data(shl1 + 851);
    const auto *shl1_852 = buffer.data(shl1 + 852);
    const auto *shl1_854 = buffer.data(shl1 + 854);
    const auto *shl1_858 = buffer.data(shl1 + 858);
    const auto *shl1_861 = buffer.data(shl1 + 861);
    const auto *shl1_865 = buffer.data(shl1 + 865);
    const auto *shl1_867 = buffer.data(shl1 + 867);
    const auto *shl1_870 = buffer.data(shl1 + 870);
    const auto *shl1_872 = buffer.data(shl1 + 872);
    const auto *shl1_873 = buffer.data(shl1 + 873);
    const auto *shl1_876 = buffer.data(shl1 + 876);
    const auto *shl1_878 = buffer.data(shl1 + 878);
    const auto *shl1_879 = buffer.data(shl1 + 879);
    const auto *shl1_880 = buffer.data(shl1 + 880);
    const auto *shl1_891 = buffer.data(shl1 + 891);
    const auto *shl1_893 = buffer.data(shl1 + 893);
    const auto *shl1_894 = buffer.data(shl1 + 894);
    const auto *shl1_895 = buffer.data(shl1 + 895);
    const auto *shl1_896 = buffer.data(shl1 + 896);
    const auto *shl1_897 = buffer.data(shl1 + 897);
    const auto *shl1_899 = buffer.data(shl1 + 899);
    const auto *shl1_900 = buffer.data(shl1 + 900);
    const auto *shl1_903 = buffer.data(shl1 + 903);
    const auto *shl1_905 = buffer.data(shl1 + 905);
    const auto *shl1_906 = buffer.data(shl1 + 906);
    const auto *shl1_909 = buffer.data(shl1 + 909);
    const auto *shl1_910 = buffer.data(shl1 + 910);

    const auto *sik_632 = buffer.data(sik + 632);
    const auto *sik_640 = buffer.data(sik + 640);
    const auto *sik_641 = buffer.data(sik + 641);
    const auto *sik_642 = buffer.data(sik + 642);
    const auto *sik_643 = buffer.data(sik + 643);
    const auto *sik_644 = buffer.data(sik + 644);
    const auto *sik_645 = buffer.data(sik + 645);
    const auto *sik_646 = buffer.data(sik + 646);
    const auto *sik_647 = buffer.data(sik + 647);
    const auto *sik_648 = buffer.data(sik + 648);
    const auto *sik_650 = buffer.data(sik + 650);
    const auto *sik_651 = buffer.data(sik + 651);
    const auto *sik_653 = buffer.data(sik + 653);
    const auto *sik_654 = buffer.data(sik + 654);
    const auto *sik_657 = buffer.data(sik + 657);
    const auto *sik_658 = buffer.data(sik + 658);
    const auto *sik_662 = buffer.data(sik + 662);
    const auto *sik_663 = buffer.data(sik + 663);
    const auto *sik_668 = buffer.data(sik + 668);
    const auto *sik_676 = buffer.data(sik + 676);
    const auto *sik_677 = buffer.data(sik + 677);
    const auto *sik_678 = buffer.data(sik + 678);
    const auto *sik_679 = buffer.data(sik + 679);
    const auto *sik_680 = buffer.data(sik + 680);
    const auto *sik_681 = buffer.data(sik + 681);
    const auto *sik_682 = buffer.data(sik + 682);
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
    const auto *sik_713 = buffer.data(sik + 713);
    const auto *sik_714 = buffer.data(sik + 714);
    const auto *sik_715 = buffer.data(sik + 715);
    const auto *sik_716 = buffer.data(sik + 716);
    const auto *sik_717 = buffer.data(sik + 717);
    const auto *sik_718 = buffer.data(sik + 718);
    const auto *sik_719 = buffer.data(sik + 719);
    const auto *sik_720 = buffer.data(sik + 720);
    const auto *sik_722 = buffer.data(sik + 722);
    const auto *sik_723 = buffer.data(sik + 723);
    const auto *sik_725 = buffer.data(sik + 725);

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pb_x, pc_x, pc_y, shl0_792, shk_452, \
                         shk_639, shk_640, shk_641, shl1_792, sik_632, sik_640, \
                         sik_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * shk_452[k]
                   + f_3 * pc_y[k] * sik_632[k];

        t_792[k] = pb_x[k] * shl0_792[k]
                   + f_16 * shk_639[k]
                   - f_14 * pc_x[k] * shl1_792[k];

        t_793[k] = f_15 * shk_640[k]
                   + f_3 * pc_x[k] * sik_640[k];

        t_794[k] = f_15 * shk_641[k]
                   + f_3 * pc_x[k] * sik_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, shk_642, shk_643, shk_644, \
                         shk_645, shk_646, sik_642, sik_643, sik_644, sik_645, \
                         sik_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_15 * shk_642[k]
                   + f_3 * pc_x[k] * sik_642[k];

        t_796[k] = f_15 * shk_643[k]
                   + f_3 * pc_x[k] * sik_643[k];

        t_797[k] = f_15 * shk_644[k]
                   + f_3 * pc_x[k] * sik_644[k];

        t_798[k] = f_15 * shk_645[k]
                   + f_3 * pc_x[k] * sik_645[k];

        t_799[k] = f_15 * shk_646[k]
                   + f_3 * pc_x[k] * sik_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, t_803, pb_x, pc_x, pc_z, shl0_801, shl0_803, \
                         shk_424, shk_647, shl1_801, shl1_803, sik_640, \
                         sik_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_15 * shk_647[k]
                   + f_3 * pc_x[k] * sik_647[k];

        t_801[k] = pb_x[k] * shl0_801[k]
                   - f_14 * pc_x[k] * shl1_801[k];

        t_802[k] = f_16 * shk_424[k]
                   + f_3 * pc_z[k] * sik_640[k];

        t_803[k] = pb_x[k] * shl0_803[k]
                   - f_14 * pc_x[k] * shl1_803[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, t_807, pb_x, pc_x, shl0_804, shl0_805, shl0_806, \
                         shl0_807, shl1_804, shl1_805, shl1_806, \
                         shl1_807 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = pb_x[k] * shl0_804[k]
                   - f_14 * pc_x[k] * shl1_804[k];

        t_805[k] = pb_x[k] * shl0_805[k]
                   - f_14 * pc_x[k] * shl1_805[k];

        t_806[k] = pb_x[k] * shl0_806[k]
                   - f_14 * pc_x[k] * shl1_806[k];

        t_807[k] = pb_x[k] * shl0_807[k]
                   - f_14 * pc_x[k] * shl1_807[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pb_x, pc_x, pc_y, shl0_809, shl0_810, \
                         shk_467, shk_468, shk_648, shl1_809, shl1_810, sik_647, \
                         sik_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_17 * shk_467[k]
                   + f_3 * pc_y[k] * sik_647[k];

        t_809[k] = pb_x[k] * shl0_809[k]
                   - f_14 * pc_x[k] * shl1_809[k];

        t_810[k] = pb_x[k] * shl0_810[k]
                   + f_20 * shk_648[k]
                   - f_14 * pc_x[k] * shl1_810[k];

        t_811[k] = f_16 * shk_468[k]
                   + f_3 * pc_y[k] * sik_648[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pb_x, pc_x, pc_y, pc_z, shl0_813, shk_432, \
                         shk_470, shk_651, shl1_813, sik_648, sik_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * shk_432[k]
                   + f_3 * pc_z[k] * sik_648[k];

        t_813[k] = pb_x[k] * shl0_813[k]
                   + f_0 * shk_651[k]
                   - f_14 * pc_x[k] * shl1_813[k];

        t_814[k] = f_16 * shk_470[k]
                   + f_3 * pc_y[k] * sik_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pb_x, pc_x, pc_z, shl0_815, shl0_816, shk_435, \
                         shk_653, shk_654, shl1_815, shl1_816, \
                         sik_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = pb_x[k] * shl0_815[k]
                   + f_0 * shk_653[k]
                   - f_14 * pc_x[k] * shl1_815[k];

        t_816[k] = pb_x[k] * shl0_816[k]
                   + f_19 * shk_654[k]
                   - f_14 * pc_x[k] * shl1_816[k];

        t_817[k] = f_17 * shk_435[k]
                   + f_3 * pc_z[k] * sik_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pb_x, pc_x, pc_y, shl0_819, shl0_820, shk_473, \
                         shk_657, shk_658, shl1_819, shl1_820, \
                         sik_653 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * shk_473[k]
                   + f_3 * pc_y[k] * sik_653[k];

        t_819[k] = pb_x[k] * shl0_819[k]
                   + f_19 * shk_657[k]
                   - f_14 * pc_x[k] * shl1_819[k];

        t_820[k] = pb_x[k] * shl0_820[k]
                   + f_18 * shk_658[k]
                   - f_14 * pc_x[k] * shl1_820[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pb_x, pc_x, pc_y, pc_z, shl0_822, shk_438, \
                         shk_477, shk_660, shl1_822, sik_654, sik_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * shk_438[k]
                   + f_3 * pc_z[k] * sik_654[k];

        t_822[k] = pb_x[k] * shl0_822[k]
                   + f_18 * shk_660[k]
                   - f_14 * pc_x[k] * shl1_822[k];

        t_823[k] = f_16 * shk_477[k]
                   + f_3 * pc_y[k] * sik_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pb_x, pc_x, pc_z, shl0_824, shl0_825, shk_442, \
                         shk_662, shk_663, shl1_824, shl1_825, \
                         sik_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = pb_x[k] * shl0_824[k]
                   + f_18 * shk_662[k]
                   - f_14 * pc_x[k] * shl1_824[k];

        t_825[k] = pb_x[k] * shl0_825[k]
                   + f_17 * shk_663[k]
                   - f_14 * pc_x[k] * shl1_825[k];

        t_826[k] = f_17 * shk_442[k]
                   + f_3 * pc_z[k] * sik_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pb_x, pc_x, pc_y, shl0_827, shl0_828, shk_482, \
                         shk_665, shk_666, shl1_827, shl1_828, \
                         sik_662 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = pb_x[k] * shl0_827[k]
                   + f_17 * shk_665[k]
                   - f_14 * pc_x[k] * shl1_827[k];

        t_828[k] = pb_x[k] * shl0_828[k]
                   + f_17 * shk_666[k]
                   - f_14 * pc_x[k] * shl1_828[k];

        t_829[k] = f_16 * shk_482[k]
                   + f_3 * pc_y[k] * sik_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pb_x, pc_x, pc_z, shl0_830, shl0_831, shk_447, \
                         shk_668, shk_669, shl1_830, shl1_831, \
                         sik_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = pb_x[k] * shl0_830[k]
                   + f_17 * shk_668[k]
                   - f_14 * pc_x[k] * shl1_830[k];

        t_831[k] = pb_x[k] * shl0_831[k]
                   + f_16 * shk_669[k]
                   - f_14 * pc_x[k] * shl1_831[k];

        t_832[k] = f_17 * shk_447[k]
                   + f_3 * pc_z[k] * sik_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pb_x, pc_x, shl0_833, shl0_834, shl0_835, \
                         shk_671, shk_672, shk_673, shl1_833, shl1_834, \
                         shl1_835 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = pb_x[k] * shl0_833[k]
                   + f_16 * shk_671[k]
                   - f_14 * pc_x[k] * shl1_833[k];

        t_834[k] = pb_x[k] * shl0_834[k]
                   + f_16 * shk_672[k]
                   - f_14 * pc_x[k] * shl1_834[k];

        t_835[k] = pb_x[k] * shl0_835[k]
                   + f_16 * shk_673[k]
                   - f_14 * pc_x[k] * shl1_835[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pb_x, pc_x, pc_y, shl0_837, shk_488, \
                         shk_675, shk_676, shk_677, shl1_837, sik_668, sik_676, \
                         sik_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * shk_488[k]
                   + f_3 * pc_y[k] * sik_668[k];

        t_837[k] = pb_x[k] * shl0_837[k]
                   + f_16 * shk_675[k]
                   - f_14 * pc_x[k] * shl1_837[k];

        t_838[k] = f_15 * shk_676[k]
                   + f_3 * pc_x[k] * sik_676[k];

        t_839[k] = f_15 * shk_677[k]
                   + f_3 * pc_x[k] * sik_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, shk_678, shk_679, shk_680, \
                         shk_681, shk_682, sik_678, sik_679, sik_680, sik_681, \
                         sik_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_15 * shk_678[k]
                   + f_3 * pc_x[k] * sik_678[k];

        t_841[k] = f_15 * shk_679[k]
                   + f_3 * pc_x[k] * sik_679[k];

        t_842[k] = f_15 * shk_680[k]
                   + f_3 * pc_x[k] * sik_680[k];

        t_843[k] = f_15 * shk_681[k]
                   + f_3 * pc_x[k] * sik_681[k];

        t_844[k] = f_15 * shk_682[k]
                   + f_3 * pc_x[k] * sik_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pb_x, pc_x, pc_z, shl0_846, shl0_848, \
                         shk_460, shk_683, shl1_846, shl1_848, sik_676, \
                         sik_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_15 * shk_683[k]
                   + f_3 * pc_x[k] * sik_683[k];

        t_846[k] = pb_x[k] * shl0_846[k]
                   - f_14 * pc_x[k] * shl1_846[k];

        t_847[k] = f_17 * shk_460[k]
                   + f_3 * pc_z[k] * sik_676[k];

        t_848[k] = pb_x[k] * shl0_848[k]
                   - f_14 * pc_x[k] * shl1_848[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pb_x, pc_x, shl0_849, shl0_850, shl0_851, \
                         shl0_852, shl1_849, shl1_850, shl1_851, \
                         shl1_852 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = pb_x[k] * shl0_849[k]
                   - f_14 * pc_x[k] * shl1_849[k];

        t_850[k] = pb_x[k] * shl0_850[k]
                   - f_14 * pc_x[k] * shl1_850[k];

        t_851[k] = pb_x[k] * shl0_851[k]
                   - f_14 * pc_x[k] * shl1_851[k];

        t_852[k] = pb_x[k] * shl0_852[k]
                   - f_14 * pc_x[k] * shl1_852[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pb_x, pb_y, pc_x, pc_y, shl0_630, \
                         shl0_854, shk_503, shk_504, shl1_630, shl1_854, sik_683, \
                         sik_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_16 * shk_503[k]
                   + f_3 * pc_y[k] * sik_683[k];

        t_854[k] = pb_x[k] * shl0_854[k]
                   - f_14 * pc_x[k] * shl1_854[k];

        t_855[k] = pb_y[k] * shl0_630[k]
                   - f_14 * pc_y[k] * shl1_630[k];

        t_856[k] = f_15 * shk_504[k]
                   + f_3 * pc_y[k] * sik_684[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pb_x, pc_x, pc_y, pc_z, shl0_858, shk_468, \
                         shk_506, shk_687, shl1_858, sik_684, sik_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_18 * shk_468[k]
                   + f_3 * pc_z[k] * sik_684[k];

        t_858[k] = pb_x[k] * shl0_858[k]
                   + f_0 * shk_687[k]
                   - f_14 * pc_x[k] * shl1_858[k];

        t_859[k] = f_15 * shk_506[k]
                   + f_3 * pc_y[k] * sik_686[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, pb_x, pb_y, pc_x, pc_y, pc_z, shl0_635, \
                         shl0_861, shk_471, shk_690, shl1_635, shl1_861, \
                         sik_687 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = pb_y[k] * shl0_635[k]
                   - f_14 * pc_y[k] * shl1_635[k];

        t_861[k] = pb_x[k] * shl0_861[k]
                   + f_19 * shk_690[k]
                   - f_14 * pc_x[k] * shl1_861[k];

        t_862[k] = f_18 * shk_471[k]
                   + f_3 * pc_z[k] * sik_687[k];
    }

#pragma omp simd aligned(t_863, t_864, t_865, pb_x, pb_y, pc_x, pc_y, shl0_639, shl0_865, \
                         shk_509, shk_694, shl1_639, shl1_865, \
                         sik_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_863[k] = f_15 * shk_509[k]
                   + f_3 * pc_y[k] * sik_689[k];

        t_864[k] = pb_y[k] * shl0_639[k]
                   - f_14 * pc_y[k] * shl1_639[k];

        t_865[k] = pb_x[k] * shl0_865[k]
                   + f_18 * shk_694[k]
                   - f_14 * pc_x[k] * shl1_865[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, pb_x, pc_x, pc_y, pc_z, shl0_867, shk_474, \
                         shk_513, shk_696, shl1_867, sik_690, sik_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * shk_474[k]
                   + f_3 * pc_z[k] * sik_690[k];

        t_867[k] = pb_x[k] * shl0_867[k]
                   + f_18 * shk_696[k]
                   - f_14 * pc_x[k] * shl1_867[k];

        t_868[k] = f_15 * shk_513[k]
                   + f_3 * pc_y[k] * sik_693[k];
    }

#pragma omp simd aligned(t_869, t_870, t_871, pb_x, pb_y, pc_x, pc_y, pc_z, shl0_644, \
                         shl0_870, shk_478, shk_699, shl1_644, shl1_870, \
                         sik_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_869[k] = pb_y[k] * shl0_644[k]
                   - f_14 * pc_y[k] * shl1_644[k];

        t_870[k] = pb_x[k] * shl0_870[k]
                   + f_17 * shk_699[k]
                   - f_14 * pc_x[k] * shl1_870[k];

        t_871[k] = f_18 * shk_478[k]
                   + f_3 * pc_z[k] * sik_694[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pb_x, pc_x, pc_y, shl0_872, shl0_873, shk_518, \
                         shk_701, shk_702, shl1_872, shl1_873, \
                         sik_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = pb_x[k] * shl0_872[k]
                   + f_17 * shk_701[k]
                   - f_14 * pc_x[k] * shl1_872[k];

        t_873[k] = pb_x[k] * shl0_873[k]
                   + f_17 * shk_702[k]
                   - f_14 * pc_x[k] * shl1_873[k];

        t_874[k] = f_15 * shk_518[k]
                   + f_3 * pc_y[k] * sik_698[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, pb_x, pb_y, pc_x, pc_y, pc_z, shl0_650, \
                         shl0_876, shk_483, shk_705, shl1_650, shl1_876, \
                         sik_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = pb_y[k] * shl0_650[k]
                   - f_14 * pc_y[k] * shl1_650[k];

        t_876[k] = pb_x[k] * shl0_876[k]
                   + f_16 * shk_705[k]
                   - f_14 * pc_x[k] * shl1_876[k];

        t_877[k] = f_18 * shk_483[k]
                   + f_3 * pc_z[k] * sik_699[k];
    }

#pragma omp simd aligned(t_878, t_879, t_880, pb_x, pc_x, shl0_878, shl0_879, shl0_880, \
                         shk_707, shk_708, shk_709, shl1_878, shl1_879, \
                         shl1_880 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_878[k] = pb_x[k] * shl0_878[k]
                   + f_16 * shk_707[k]
                   - f_14 * pc_x[k] * shl1_878[k];

        t_879[k] = pb_x[k] * shl0_879[k]
                   + f_16 * shk_708[k]
                   - f_14 * pc_x[k] * shl1_879[k];

        t_880[k] = pb_x[k] * shl0_880[k]
                   + f_16 * shk_709[k]
                   - f_14 * pc_x[k] * shl1_880[k];
    }

#pragma omp simd aligned(t_881, t_882, t_883, t_884, pb_y, pc_x, pc_y, shl0_657, shk_524, \
                         shk_712, shk_713, shl1_657, sik_704, sik_712, \
                         sik_713 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_881[k] = f_15 * shk_524[k]
                   + f_3 * pc_y[k] * sik_704[k];

        t_882[k] = pb_y[k] * shl0_657[k]
                   - f_14 * pc_y[k] * shl1_657[k];

        t_883[k] = f_15 * shk_712[k]
                   + f_3 * pc_x[k] * sik_712[k];

        t_884[k] = f_15 * shk_713[k]
                   + f_3 * pc_x[k] * sik_713[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, t_888, t_889, pc_x, shk_714, shk_715, shk_716, \
                         shk_717, shk_718, sik_714, sik_715, sik_716, sik_717, \
                         sik_718 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = f_15 * shk_714[k]
                   + f_3 * pc_x[k] * sik_714[k];

        t_886[k] = f_15 * shk_715[k]
                   + f_3 * pc_x[k] * sik_715[k];

        t_887[k] = f_15 * shk_716[k]
                   + f_3 * pc_x[k] * sik_716[k];

        t_888[k] = f_15 * shk_717[k]
                   + f_3 * pc_x[k] * sik_717[k];

        t_889[k] = f_15 * shk_718[k]
                   + f_3 * pc_x[k] * sik_718[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pb_x, pc_x, pc_z, shl0_891, shl0_893, \
                         shk_496, shk_719, shl1_891, shl1_893, sik_712, \
                         sik_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = f_15 * shk_719[k]
                   + f_3 * pc_x[k] * sik_719[k];

        t_891[k] = pb_x[k] * shl0_891[k]
                   - f_14 * pc_x[k] * shl1_891[k];

        t_892[k] = f_18 * shk_496[k]
                   + f_3 * pc_z[k] * sik_712[k];

        t_893[k] = pb_x[k] * shl0_893[k]
                   - f_14 * pc_x[k] * shl1_893[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pb_x, pc_x, shl0_894, shl0_895, shl0_896, \
                         shl0_897, shl1_894, shl1_895, shl1_896, \
                         shl1_897 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = pb_x[k] * shl0_894[k]
                   - f_14 * pc_x[k] * shl1_894[k];

        t_895[k] = pb_x[k] * shl0_895[k]
                   - f_14 * pc_x[k] * shl1_895[k];

        t_896[k] = pb_x[k] * shl0_896[k]
                   - f_14 * pc_x[k] * shl1_896[k];

        t_897[k] = pb_x[k] * shl0_897[k]
                   - f_14 * pc_x[k] * shl1_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pb_x, pc_x, pc_y, shl0_899, shl0_900, \
                         shk_539, shk_720, shl1_899, shl1_900, sik_719, \
                         sik_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_15 * shk_539[k]
                   + f_3 * pc_y[k] * sik_719[k];

        t_899[k] = pb_x[k] * shl0_899[k]
                   - f_14 * pc_x[k] * shl1_899[k];

        t_900[k] = pb_x[k] * shl0_900[k]
                   + f_20 * shk_720[k]
                   - f_14 * pc_x[k] * shl1_900[k];

        t_901[k] = f_3 * pc_y[k] * sik_720[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pb_x, pc_x, pc_y, pc_z, shl0_903, shk_504, \
                         shk_723, shl1_903, sik_720, sik_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_19 * shk_504[k]
                   + f_3 * pc_z[k] * sik_720[k];

        t_903[k] = pb_x[k] * shl0_903[k]
                   + f_0 * shk_723[k]
                   - f_14 * pc_x[k] * shl1_903[k];

        t_904[k] = f_3 * pc_y[k] * sik_722[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, pb_x, pc_x, pc_z, shl0_905, shl0_906, shk_507, \
                         shk_725, shk_726, shl1_905, shl1_906, \
                         sik_723 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = pb_x[k] * shl0_905[k]
                   + f_0 * shk_725[k]
                   - f_14 * pc_x[k] * shl1_905[k];

        t_906[k] = pb_x[k] * shl0_906[k]
                   + f_19 * shk_726[k]
                   - f_14 * pc_x[k] * shl1_906[k];

        t_907[k] = f_19 * shk_507[k]
                   + f_3 * pc_z[k] * sik_723[k];
    }

#pragma omp simd aligned(t_908, t_909, t_910, pb_x, pc_x, pc_y, shl0_909, shl0_910, shk_729, \
                         shk_730, shl1_909, shl1_910, sik_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_908[k] = f_3 * pc_y[k] * sik_725[k];

        t_909[k] = pb_x[k] * shl0_909[k]
                   + f_19 * shk_729[k]
                   - f_14 * pc_x[k] * shl1_909[k];

        t_910[k] = pb_x[k] * shl0_910[k]
                   + f_18 * shk_730[k]
                   - f_14 * pc_x[k] * shl1_910[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shl0,
                                                          const size_t shk, const size_t shl1,
                                                          const size_t sii0, const size_t sii1,
                                                          const size_t sik, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shl0_675 = buffer.data(shl0 + 675);
    const auto *shl0_678 = buffer.data(shl0 + 678);
    const auto *shl0_681 = buffer.data(shl0 + 681);
    const auto *shl0_685 = buffer.data(shl0 + 685);
    const auto *shl0_690 = buffer.data(shl0 + 690);
    const auto *shl0_696 = buffer.data(shl0 + 696);
    const auto *shl0_711 = buffer.data(shl0 + 711);
    const auto *shl0_713 = buffer.data(shl0 + 713);
    const auto *shl0_714 = buffer.data(shl0 + 714);
    const auto *shl0_715 = buffer.data(shl0 + 715);
    const auto *shl0_716 = buffer.data(shl0 + 716);
    const auto *shl0_717 = buffer.data(shl0 + 717);
    const auto *shl0_912 = buffer.data(shl0 + 912);
    const auto *shl0_914 = buffer.data(shl0 + 914);
    const auto *shl0_915 = buffer.data(shl0 + 915);
    const auto *shl0_917 = buffer.data(shl0 + 917);
    const auto *shl0_918 = buffer.data(shl0 + 918);
    const auto *shl0_920 = buffer.data(shl0 + 920);
    const auto *shl0_921 = buffer.data(shl0 + 921);
    const auto *shl0_923 = buffer.data(shl0 + 923);
    const auto *shl0_924 = buffer.data(shl0 + 924);
    const auto *shl0_925 = buffer.data(shl0 + 925);
    const auto *shl0_927 = buffer.data(shl0 + 927);
    const auto *shl0_936 = buffer.data(shl0 + 936);
    const auto *shl0_938 = buffer.data(shl0 + 938);
    const auto *shl0_939 = buffer.data(shl0 + 939);
    const auto *shl0_940 = buffer.data(shl0 + 940);
    const auto *shl0_941 = buffer.data(shl0 + 941);
    const auto *shl0_942 = buffer.data(shl0 + 942);
    const auto *shl0_944 = buffer.data(shl0 + 944);

    const auto *shk_510 = buffer.data(shk + 510);
    const auto *shk_514 = buffer.data(shk + 514);
    const auto *shk_519 = buffer.data(shk + 519);
    const auto *shk_532 = buffer.data(shk + 532);
    const auto *shk_540 = buffer.data(shk + 540);
    const auto *shk_542 = buffer.data(shk + 542);
    const auto *shk_543 = buffer.data(shk + 543);
    const auto *shk_545 = buffer.data(shk + 545);
    const auto *shk_546 = buffer.data(shk + 546);
    const auto *shk_549 = buffer.data(shk + 549);
    const auto *shk_550 = buffer.data(shk + 550);
    const auto *shk_554 = buffer.data(shk + 554);
    const auto *shk_555 = buffer.data(shk + 555);
    const auto *shk_560 = buffer.data(shk + 560);
    const auto *shk_568 = buffer.data(shk + 568);
    const auto *shk_569 = buffer.data(shk + 569);
    const auto *shk_570 = buffer.data(shk + 570);
    const auto *shk_571 = buffer.data(shk + 571);
    const auto *shk_572 = buffer.data(shk + 572);
    const auto *shk_573 = buffer.data(shk + 573);
    const auto *shk_574 = buffer.data(shk + 574);
    const auto *shk_575 = buffer.data(shk + 575);
    const auto *shk_576 = buffer.data(shk + 576);
    const auto *shk_578 = buffer.data(shk + 578);
    const auto *shk_581 = buffer.data(shk + 581);
    const auto *shk_585 = buffer.data(shk + 585);
    const auto *shk_590 = buffer.data(shk + 590);
    const auto *shk_596 = buffer.data(shk + 596);
    const auto *shk_611 = buffer.data(shk + 611);
    const auto *shk_732 = buffer.data(shk + 732);
    const auto *shk_734 = buffer.data(shk + 734);
    const auto *shk_735 = buffer.data(shk + 735);
    const auto *shk_737 = buffer.data(shk + 737);
    const auto *shk_738 = buffer.data(shk + 738);
    const auto *shk_740 = buffer.data(shk + 740);
    const auto *shk_741 = buffer.data(shk + 741);
    const auto *shk_743 = buffer.data(shk + 743);
    const auto *shk_744 = buffer.data(shk + 744);
    const auto *shk_745 = buffer.data(shk + 745);
    const auto *shk_747 = buffer.data(shk + 747);
    const auto *shk_748 = buffer.data(shk + 748);
    const auto *shk_749 = buffer.data(shk + 749);
    const auto *shk_750 = buffer.data(shk + 750);
    const auto *shk_751 = buffer.data(shk + 751);
    const auto *shk_752 = buffer.data(shk + 752);
    const auto *shk_753 = buffer.data(shk + 753);
    const auto *shk_754 = buffer.data(shk + 754);
    const auto *shk_755 = buffer.data(shk + 755);

    const auto *shl1_675 = buffer.data(shl1 + 675);
    const auto *shl1_678 = buffer.data(shl1 + 678);
    const auto *shl1_681 = buffer.data(shl1 + 681);
    const auto *shl1_685 = buffer.data(shl1 + 685);
    const auto *shl1_690 = buffer.data(shl1 + 690);
    const auto *shl1_696 = buffer.data(shl1 + 696);
    const auto *shl1_711 = buffer.data(shl1 + 711);
    const auto *shl1_713 = buffer.data(shl1 + 713);
    const auto *shl1_714 = buffer.data(shl1 + 714);
    const auto *shl1_715 = buffer.data(shl1 + 715);
    const auto *shl1_716 = buffer.data(shl1 + 716);
    const auto *shl1_717 = buffer.data(shl1 + 717);
    const auto *shl1_912 = buffer.data(shl1 + 912);
    const auto *shl1_914 = buffer.data(shl1 + 914);
    const auto *shl1_915 = buffer.data(shl1 + 915);
    const auto *shl1_917 = buffer.data(shl1 + 917);
    const auto *shl1_918 = buffer.data(shl1 + 918);
    const auto *shl1_920 = buffer.data(shl1 + 920);
    const auto *shl1_921 = buffer.data(shl1 + 921);
    const auto *shl1_923 = buffer.data(shl1 + 923);
    const auto *shl1_924 = buffer.data(shl1 + 924);
    const auto *shl1_925 = buffer.data(shl1 + 925);
    const auto *shl1_927 = buffer.data(shl1 + 927);
    const auto *shl1_936 = buffer.data(shl1 + 936);
    const auto *shl1_938 = buffer.data(shl1 + 938);
    const auto *shl1_939 = buffer.data(shl1 + 939);
    const auto *shl1_940 = buffer.data(shl1 + 940);
    const auto *shl1_941 = buffer.data(shl1 + 941);
    const auto *shl1_942 = buffer.data(shl1 + 942);
    const auto *shl1_944 = buffer.data(shl1 + 944);

    const auto *sii0_588 = buffer.data(sii0 + 588);
    const auto *sii0_591 = buffer.data(sii0 + 591);
    const auto *sii0_593 = buffer.data(sii0 + 593);
    const auto *sii0_594 = buffer.data(sii0 + 594);
    const auto *sii0_597 = buffer.data(sii0 + 597);
    const auto *sii0_598 = buffer.data(sii0 + 598);
    const auto *sii0_600 = buffer.data(sii0 + 600);
    const auto *sii0_602 = buffer.data(sii0 + 602);
    const auto *sii0_603 = buffer.data(sii0 + 603);
    const auto *sii0_605 = buffer.data(sii0 + 605);
    const auto *sii0_606 = buffer.data(sii0 + 606);
    const auto *sii0_608 = buffer.data(sii0 + 608);
    const auto *sii0_609 = buffer.data(sii0 + 609);
    const auto *sii0_611 = buffer.data(sii0 + 611);
    const auto *sii0_612 = buffer.data(sii0 + 612);
    const auto *sii0_613 = buffer.data(sii0 + 613);
    const auto *sii0_614 = buffer.data(sii0 + 614);
    const auto *sii0_615 = buffer.data(sii0 + 615);
    const auto *sii0_621 = buffer.data(sii0 + 621);
    const auto *sii0_625 = buffer.data(sii0 + 625);
    const auto *sii0_628 = buffer.data(sii0 + 628);
    const auto *sii0_630 = buffer.data(sii0 + 630);
    const auto *sii0_633 = buffer.data(sii0 + 633);
    const auto *sii0_634 = buffer.data(sii0 + 634);
    const auto *sii0_636 = buffer.data(sii0 + 636);
    const auto *sii0_639 = buffer.data(sii0 + 639);
    const auto *sii0_640 = buffer.data(sii0 + 640);
    const auto *sii0_641 = buffer.data(sii0 + 641);
    const auto *sii0_643 = buffer.data(sii0 + 643);

    const auto *sii1_588 = buffer.data(sii1 + 588);
    const auto *sii1_591 = buffer.data(sii1 + 591);
    const auto *sii1_593 = buffer.data(sii1 + 593);
    const auto *sii1_594 = buffer.data(sii1 + 594);
    const auto *sii1_597 = buffer.data(sii1 + 597);
    const auto *sii1_598 = buffer.data(sii1 + 598);
    const auto *sii1_600 = buffer.data(sii1 + 600);
    const auto *sii1_602 = buffer.data(sii1 + 602);
    const auto *sii1_603 = buffer.data(sii1 + 603);
    const auto *sii1_605 = buffer.data(sii1 + 605);
    const auto *sii1_606 = buffer.data(sii1 + 606);
    const auto *sii1_608 = buffer.data(sii1 + 608);
    const auto *sii1_609 = buffer.data(sii1 + 609);
    const auto *sii1_611 = buffer.data(sii1 + 611);
    const auto *sii1_612 = buffer.data(sii1 + 612);
    const auto *sii1_613 = buffer.data(sii1 + 613);
    const auto *sii1_614 = buffer.data(sii1 + 614);
    const auto *sii1_615 = buffer.data(sii1 + 615);
    const auto *sii1_621 = buffer.data(sii1 + 621);
    const auto *sii1_625 = buffer.data(sii1 + 625);
    const auto *sii1_628 = buffer.data(sii1 + 628);
    const auto *sii1_630 = buffer.data(sii1 + 630);
    const auto *sii1_633 = buffer.data(sii1 + 633);
    const auto *sii1_634 = buffer.data(sii1 + 634);
    const auto *sii1_636 = buffer.data(sii1 + 636);
    const auto *sii1_639 = buffer.data(sii1 + 639);
    const auto *sii1_640 = buffer.data(sii1 + 640);
    const auto *sii1_641 = buffer.data(sii1 + 641);
    const auto *sii1_643 = buffer.data(sii1 + 643);

    const auto *sik_726 = buffer.data(sik + 726);
    const auto *sik_729 = buffer.data(sik + 729);
    const auto *sik_730 = buffer.data(sik + 730);
    const auto *sik_734 = buffer.data(sik + 734);
    const auto *sik_735 = buffer.data(sik + 735);
    const auto *sik_740 = buffer.data(sik + 740);
    const auto *sik_748 = buffer.data(sik + 748);
    const auto *sik_749 = buffer.data(sik + 749);
    const auto *sik_750 = buffer.data(sik + 750);
    const auto *sik_751 = buffer.data(sik + 751);
    const auto *sik_752 = buffer.data(sik + 752);
    const auto *sik_753 = buffer.data(sik + 753);
    const auto *sik_754 = buffer.data(sik + 754);
    const auto *sik_755 = buffer.data(sik + 755);
    const auto *sik_756 = buffer.data(sik + 756);
    const auto *sik_758 = buffer.data(sik + 758);
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
    const auto *sik_792 = buffer.data(sik + 792);
    const auto *sik_794 = buffer.data(sik + 794);
    const auto *sik_795 = buffer.data(sik + 795);
    const auto *sik_797 = buffer.data(sik + 797);
    const auto *sik_798 = buffer.data(sik + 798);
    const auto *sik_801 = buffer.data(sik + 801);
    const auto *sik_802 = buffer.data(sik + 802);
    const auto *sik_804 = buffer.data(sik + 804);
    const auto *sik_806 = buffer.data(sik + 806);
    const auto *sik_807 = buffer.data(sik + 807);
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

#pragma omp simd aligned(t_911, t_912, t_913, pb_x, pc_x, pc_y, pc_z, shl0_912, shk_510, \
                         shk_732, shl1_912, sik_726, sik_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_911[k] = f_19 * shk_510[k]
                   + f_3 * pc_z[k] * sik_726[k];

        t_912[k] = pb_x[k] * shl0_912[k]
                   + f_18 * shk_732[k]
                   - f_14 * pc_x[k] * shl1_912[k];

        t_913[k] = f_3 * pc_y[k] * sik_729[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, pb_x, pc_x, pc_z, shl0_914, shl0_915, shk_514, \
                         shk_734, shk_735, shl1_914, shl1_915, \
                         sik_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = pb_x[k] * shl0_914[k]
                   + f_18 * shk_734[k]
                   - f_14 * pc_x[k] * shl1_914[k];

        t_915[k] = pb_x[k] * shl0_915[k]
                   + f_17 * shk_735[k]
                   - f_14 * pc_x[k] * shl1_915[k];

        t_916[k] = f_19 * shk_514[k]
                   + f_3 * pc_z[k] * sik_730[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pb_x, pc_x, pc_y, shl0_917, shl0_918, shk_737, \
                         shk_738, shl1_917, shl1_918, sik_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = pb_x[k] * shl0_917[k]
                   + f_17 * shk_737[k]
                   - f_14 * pc_x[k] * shl1_917[k];

        t_918[k] = pb_x[k] * shl0_918[k]
                   + f_17 * shk_738[k]
                   - f_14 * pc_x[k] * shl1_918[k];

        t_919[k] = f_3 * pc_y[k] * sik_734[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, pb_x, pc_x, pc_z, shl0_920, shl0_921, shk_519, \
                         shk_740, shk_741, shl1_920, shl1_921, \
                         sik_735 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pb_x[k] * shl0_920[k]
                   + f_17 * shk_740[k]
                   - f_14 * pc_x[k] * shl1_920[k];

        t_921[k] = pb_x[k] * shl0_921[k]
                   + f_16 * shk_741[k]
                   - f_14 * pc_x[k] * shl1_921[k];

        t_922[k] = f_19 * shk_519[k]
                   + f_3 * pc_z[k] * sik_735[k];
    }

#pragma omp simd aligned(t_923, t_924, t_925, pb_x, pc_x, shl0_923, shl0_924, shl0_925, \
                         shk_743, shk_744, shk_745, shl1_923, shl1_924, \
                         shl1_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_923[k] = pb_x[k] * shl0_923[k]
                   + f_16 * shk_743[k]
                   - f_14 * pc_x[k] * shl1_923[k];

        t_924[k] = pb_x[k] * shl0_924[k]
                   + f_16 * shk_744[k]
                   - f_14 * pc_x[k] * shl1_924[k];

        t_925[k] = pb_x[k] * shl0_925[k]
                   + f_16 * shk_745[k]
                   - f_14 * pc_x[k] * shl1_925[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, t_929, pb_x, pc_x, pc_y, shl0_927, shk_747, \
                         shk_748, shk_749, shl1_927, sik_740, sik_748, \
                         sik_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_3 * pc_y[k] * sik_740[k];

        t_927[k] = pb_x[k] * shl0_927[k]
                   + f_16 * shk_747[k]
                   - f_14 * pc_x[k] * shl1_927[k];

        t_928[k] = f_15 * shk_748[k]
                   + f_3 * pc_x[k] * sik_748[k];

        t_929[k] = f_15 * shk_749[k]
                   + f_3 * pc_x[k] * sik_749[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, t_934, pc_x, shk_750, shk_751, shk_752, \
                         shk_753, shk_754, sik_750, sik_751, sik_752, sik_753, \
                         sik_754 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_15 * shk_750[k]
                   + f_3 * pc_x[k] * sik_750[k];

        t_931[k] = f_15 * shk_751[k]
                   + f_3 * pc_x[k] * sik_751[k];

        t_932[k] = f_15 * shk_752[k]
                   + f_3 * pc_x[k] * sik_752[k];

        t_933[k] = f_15 * shk_753[k]
                   + f_3 * pc_x[k] * sik_753[k];

        t_934[k] = f_15 * shk_754[k]
                   + f_3 * pc_x[k] * sik_754[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pb_x, pc_x, pc_z, shl0_936, shl0_938, \
                         shk_532, shk_755, shl1_936, shl1_938, sik_748, \
                         sik_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = f_15 * shk_755[k]
                   + f_3 * pc_x[k] * sik_755[k];

        t_936[k] = pb_x[k] * shl0_936[k]
                   - f_14 * pc_x[k] * shl1_936[k];

        t_937[k] = f_19 * shk_532[k]
                   + f_3 * pc_z[k] * sik_748[k];

        t_938[k] = pb_x[k] * shl0_938[k]
                   - f_14 * pc_x[k] * shl1_938[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pb_x, pc_x, shl0_939, shl0_940, shl0_941, \
                         shl0_942, shl1_939, shl1_940, shl1_941, \
                         shl1_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = pb_x[k] * shl0_939[k]
                   - f_14 * pc_x[k] * shl1_939[k];

        t_940[k] = pb_x[k] * shl0_940[k]
                   - f_14 * pc_x[k] * shl1_940[k];

        t_941[k] = pb_x[k] * shl0_941[k]
                   - f_14 * pc_x[k] * shl1_941[k];

        t_942[k] = pb_x[k] * shl0_942[k]
                   - f_14 * pc_x[k] * shl1_942[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, t_947, pb_x, pc_x, pc_y, pc_z, shl0_944, \
                         shk_540, shl1_944, sii0_588, sii1_588, sik_755, \
                         sik_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_3 * pc_y[k] * sik_755[k];

        t_944[k] = pb_x[k] * shl0_944[k]
                   - f_14 * pc_x[k] * shl1_944[k];

        t_945[k] = f_1 * sii0_588[k]
                   - f_2 * sii1_588[k]
                   + f_3 * pc_x[k] * sik_756[k];

        t_946[k] = f_0 * shk_540[k]
                   + f_3 * pc_y[k] * sik_756[k];

        t_947[k] = f_3 * pc_z[k] * sik_756[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pc_x, pc_y, shk_542, sii0_591, sii0_593, \
                         sii1_591, sii1_593, sik_758, sik_759, \
                         sik_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_4 * sii0_591[k]
                   - f_5 * sii1_591[k]
                   + f_3 * pc_x[k] * sik_759[k];

        t_949[k] = f_0 * shk_542[k]
                   + f_3 * pc_y[k] * sik_758[k];

        t_950[k] = f_4 * sii0_593[k]
                   - f_5 * sii1_593[k]
                   + f_3 * pc_x[k] * sik_761[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, t_954, pc_x, pc_y, pc_z, shk_545, sii0_594, \
                         sii0_597, sii1_594, sii1_597, sik_759, sik_761, sik_762, \
                         sik_765 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_6 * sii0_594[k]
                   - f_7 * sii1_594[k]
                   + f_3 * pc_x[k] * sik_762[k];

        t_952[k] = f_3 * pc_z[k] * sik_759[k];

        t_953[k] = f_0 * shk_545[k]
                   + f_3 * pc_y[k] * sik_761[k];

        t_954[k] = f_6 * sii0_597[k]
                   - f_7 * sii1_597[k]
                   + f_3 * pc_x[k] * sik_765[k];
    }

#pragma omp simd aligned(t_955, t_956, t_957, t_958, pc_x, pc_y, pc_z, shk_549, sii0_598, \
                         sii0_600, sii1_598, sii1_600, sik_762, sik_765, sik_766, \
                         sik_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_955[k] = f_8 * sii0_598[k]
                   - f_9 * sii1_598[k]
                   + f_3 * pc_x[k] * sik_766[k];

        t_956[k] = f_3 * pc_z[k] * sik_762[k];

        t_957[k] = f_8 * sii0_600[k]
                   - f_9 * sii1_600[k]
                   + f_3 * pc_x[k] * sik_768[k];

        t_958[k] = f_0 * shk_549[k]
                   + f_3 * pc_y[k] * sik_765[k];
    }

#pragma omp simd aligned(t_959, t_960, t_961, t_962, pc_x, pc_z, sii0_602, sii0_603, sii0_605, \
                         sii1_602, sii1_603, sii1_605, sik_766, sik_770, sik_771, \
                         sik_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_959[k] = f_8 * sii0_602[k]
                   - f_9 * sii1_602[k]
                   + f_3 * pc_x[k] * sik_770[k];

        t_960[k] = f_10 * sii0_603[k]
                   - f_11 * sii1_603[k]
                   + f_3 * pc_x[k] * sik_771[k];

        t_961[k] = f_3 * pc_z[k] * sik_766[k];

        t_962[k] = f_10 * sii0_605[k]
                   - f_11 * sii1_605[k]
                   + f_3 * pc_x[k] * sik_773[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_x, pc_y, shk_554, sii0_606, sii0_608, \
                         sii1_606, sii1_608, sik_770, sik_774, \
                         sik_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_10 * sii0_606[k]
                   - f_11 * sii1_606[k]
                   + f_3 * pc_x[k] * sik_774[k];

        t_964[k] = f_0 * shk_554[k]
                   + f_3 * pc_y[k] * sik_770[k];

        t_965[k] = f_10 * sii0_608[k]
                   - f_11 * sii1_608[k]
                   + f_3 * pc_x[k] * sik_776[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, t_969, pc_x, pc_z, sii0_609, sii0_611, sii0_612, \
                         sii1_609, sii1_611, sii1_612, sik_771, sik_777, sik_779, \
                         sik_780 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_12 * sii0_609[k]
                   - f_13 * sii1_609[k]
                   + f_3 * pc_x[k] * sik_777[k];

        t_967[k] = f_3 * pc_z[k] * sik_771[k];

        t_968[k] = f_12 * sii0_611[k]
                   - f_13 * sii1_611[k]
                   + f_3 * pc_x[k] * sik_779[k];

        t_969[k] = f_12 * sii0_612[k]
                   - f_13 * sii1_612[k]
                   + f_3 * pc_x[k] * sik_780[k];
    }

#pragma omp simd aligned(t_970, t_971, t_972, t_973, pc_x, pc_y, shk_560, sii0_613, sii0_615, \
                         sii1_613, sii1_615, sik_776, sik_781, sik_783, \
                         sik_784 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_970[k] = f_12 * sii0_613[k]
                   - f_13 * sii1_613[k]
                   + f_3 * pc_x[k] * sik_781[k];

        t_971[k] = f_0 * shk_560[k]
                   + f_3 * pc_y[k] * sik_776[k];

        t_972[k] = f_12 * sii0_615[k]
                   - f_13 * sii1_615[k]
                   + f_3 * pc_x[k] * sik_783[k];

        t_973[k] = f_3 * pc_x[k] * sik_784[k];
    }

#pragma omp simd aligned(t_974, t_975, t_976, t_977, t_978, t_979, t_980, pc_x, sik_785, \
                         sik_786, sik_787, sik_788, sik_789, sik_790, \
                         sik_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_974[k] = f_3 * pc_x[k] * sik_785[k];

        t_975[k] = f_3 * pc_x[k] * sik_786[k];

        t_976[k] = f_3 * pc_x[k] * sik_787[k];

        t_977[k] = f_3 * pc_x[k] * sik_788[k];

        t_978[k] = f_3 * pc_x[k] * sik_789[k];

        t_979[k] = f_3 * pc_x[k] * sik_790[k];

        t_980[k] = f_3 * pc_x[k] * sik_791[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_y, pc_z, shk_568, shk_570, sii0_609, \
                         sii0_611, sii1_609, sii1_611, sik_784, \
                         sik_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_0 * shk_568[k]
                   + f_1 * sii0_609[k]
                   - f_2 * sii1_609[k]
                   + f_3 * pc_y[k] * sik_784[k];

        t_982[k] = f_3 * pc_z[k] * sik_784[k];

        t_983[k] = f_0 * shk_570[k]
                   + f_4 * sii0_611[k]
                   - f_5 * sii1_611[k]
                   + f_3 * pc_y[k] * sik_786[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_y, shk_571, shk_572, shk_573, sii0_612, \
                         sii0_613, sii0_614, sii1_612, sii1_613, sii1_614, sik_787, sik_788, \
                         sik_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_0 * shk_571[k]
                   + f_6 * sii0_612[k]
                   - f_7 * sii1_612[k]
                   + f_3 * pc_y[k] * sik_787[k];

        t_985[k] = f_0 * shk_572[k]
                   + f_8 * sii0_613[k]
                   - f_9 * sii1_613[k]
                   + f_3 * pc_y[k] * sik_788[k];

        t_986[k] = f_0 * shk_573[k]
                   + f_10 * sii0_614[k]
                   - f_11 * sii1_614[k]
                   + f_3 * pc_y[k] * sik_789[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, pb_z, pc_y, pc_z, shl0_675, shk_574, \
                         shk_575, shl1_675, sii0_615, sii1_615, sik_790, \
                         sik_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_0 * shk_574[k]
                   + f_12 * sii0_615[k]
                   - f_13 * sii1_615[k]
                   + f_3 * pc_y[k] * sik_790[k];

        t_988[k] = f_0 * shk_575[k]
                   + f_3 * pc_y[k] * sik_791[k];

        t_989[k] = f_1 * sii0_615[k]
                   - f_2 * sii1_615[k]
                   + f_3 * pc_z[k] * sik_791[k];

        t_990[k] = pb_z[k] * shl0_675[k]
                   - f_14 * pc_z[k] * shl1_675[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_z, pc_y, pc_z, shl0_678, shk_540, \
                         shk_576, shk_578, shl1_678, sik_792, sik_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_19 * shk_576[k]
                   + f_3 * pc_y[k] * sik_792[k];

        t_992[k] = f_15 * shk_540[k]
                   + f_3 * pc_z[k] * sik_792[k];

        t_993[k] = pb_z[k] * shl0_678[k]
                   - f_14 * pc_z[k] * shl1_678[k];

        t_994[k] = f_19 * shk_578[k]
                   + f_3 * pc_y[k] * sik_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, t_998, pb_z, pc_x, pc_y, pc_z, shl0_681, \
                         shk_543, shk_581, shl1_681, sii0_621, sii1_621, sik_795, \
                         sik_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_4 * sii0_621[k]
                   - f_5 * sii1_621[k]
                   + f_3 * pc_x[k] * sik_797[k];

        t_996[k] = pb_z[k] * shl0_681[k]
                   - f_14 * pc_z[k] * shl1_681[k];

        t_997[k] = f_15 * shk_543[k]
                   + f_3 * pc_z[k] * sik_795[k];

        t_998[k] = f_19 * shk_581[k]
                   + f_3 * pc_y[k] * sik_797[k];
    }

#pragma omp simd aligned(t_999, t_1000, t_1001, pb_z, pc_x, pc_z, shl0_685, shk_546, shl1_685, \
                         sii0_625, sii1_625, sik_798, sik_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_999[k] = f_6 * sii0_625[k]
                   - f_7 * sii1_625[k]
                   + f_3 * pc_x[k] * sik_801[k];

        t_1000[k] = pb_z[k] * shl0_685[k]
                    - f_14 * pc_z[k] * shl1_685[k];

        t_1001[k] = f_15 * shk_546[k]
                    + f_3 * pc_z[k] * sik_798[k];
    }

#pragma omp simd aligned(t_1002, t_1003, t_1004, pc_x, pc_y, shk_585, sii0_628, sii0_630, \
                         sii1_628, sii1_630, sik_801, sik_804, \
                         sik_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1002[k] = f_8 * sii0_628[k]
                    - f_9 * sii1_628[k]
                    + f_3 * pc_x[k] * sik_804[k];

        t_1003[k] = f_19 * shk_585[k]
                    + f_3 * pc_y[k] * sik_801[k];

        t_1004[k] = f_8 * sii0_630[k]
                    - f_9 * sii1_630[k]
                    + f_3 * pc_x[k] * sik_806[k];
    }

#pragma omp simd aligned(t_1005, t_1006, t_1007, pb_z, pc_x, pc_z, shl0_690, shk_550, \
                         shl1_690, sii0_633, sii1_633, sik_802, \
                         sik_809 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1005[k] = pb_z[k] * shl0_690[k]
                    - f_14 * pc_z[k] * shl1_690[k];

        t_1006[k] = f_15 * shk_550[k]
                    + f_3 * pc_z[k] * sik_802[k];

        t_1007[k] = f_10 * sii0_633[k]
                    - f_11 * sii1_633[k]
                    + f_3 * pc_x[k] * sik_809[k];
    }

#pragma omp simd aligned(t_1008, t_1009, t_1010, pc_x, pc_y, shk_590, sii0_634, sii0_636, \
                         sii1_634, sii1_636, sik_806, sik_810, \
                         sik_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1008[k] = f_10 * sii0_634[k]
                    - f_11 * sii1_634[k]
                    + f_3 * pc_x[k] * sik_810[k];

        t_1009[k] = f_19 * shk_590[k]
                    + f_3 * pc_y[k] * sik_806[k];

        t_1010[k] = f_10 * sii0_636[k]
                    - f_11 * sii1_636[k]
                    + f_3 * pc_x[k] * sik_812[k];
    }

#pragma omp simd aligned(t_1011, t_1012, t_1013, pb_z, pc_x, pc_z, shl0_696, shk_555, \
                         shl1_696, sii0_639, sii1_639, sik_807, \
                         sik_815 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1011[k] = pb_z[k] * shl0_696[k]
                    - f_14 * pc_z[k] * shl1_696[k];

        t_1012[k] = f_15 * shk_555[k]
                    + f_3 * pc_z[k] * sik_807[k];

        t_1013[k] = f_12 * sii0_639[k]
                    - f_13 * sii1_639[k]
                    + f_3 * pc_x[k] * sik_815[k];
    }

#pragma omp simd aligned(t_1014, t_1015, t_1016, pc_x, pc_y, shk_596, sii0_640, sii0_641, \
                         sii1_640, sii1_641, sik_812, sik_816, \
                         sik_817 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1014[k] = f_12 * sii0_640[k]
                    - f_13 * sii1_640[k]
                    + f_3 * pc_x[k] * sik_816[k];

        t_1015[k] = f_12 * sii0_641[k]
                    - f_13 * sii1_641[k]
                    + f_3 * pc_x[k] * sik_817[k];

        t_1016[k] = f_19 * shk_596[k]
                    + f_3 * pc_y[k] * sik_812[k];
    }

#pragma omp simd aligned(t_1017, t_1018, t_1019, t_1020, t_1021, t_1022, pc_x, sii0_643, \
                         sii1_643, sik_819, sik_820, sik_821, sik_822, sik_823, \
                         sik_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1017[k] = f_12 * sii0_643[k]
                    - f_13 * sii1_643[k]
                    + f_3 * pc_x[k] * sik_819[k];

        t_1018[k] = f_3 * pc_x[k] * sik_820[k];

        t_1019[k] = f_3 * pc_x[k] * sik_821[k];

        t_1020[k] = f_3 * pc_x[k] * sik_822[k];

        t_1021[k] = f_3 * pc_x[k] * sik_823[k];

        t_1022[k] = f_3 * pc_x[k] * sik_824[k];
    }

#pragma omp simd aligned(t_1023, t_1024, t_1025, t_1026, t_1027, pb_z, pc_x, pc_z, shl0_711, \
                         shk_568, shl1_711, sik_820, sik_825, sik_826, \
                         sik_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1023[k] = f_3 * pc_x[k] * sik_825[k];

        t_1024[k] = f_3 * pc_x[k] * sik_826[k];

        t_1025[k] = f_3 * pc_x[k] * sik_827[k];

        t_1026[k] = pb_z[k] * shl0_711[k]
                    - f_14 * pc_z[k] * shl1_711[k];

        t_1027[k] = f_15 * shk_568[k]
                    + f_3 * pc_z[k] * sik_820[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pb_z, pc_z, shl0_713, shl0_714, shl0_715, \
                         shk_569, shk_570, shk_571, shl1_713, shl1_714, \
                         shl1_715 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = pb_z[k] * shl0_713[k]
                    + f_16 * shk_569[k]
                    - f_14 * pc_z[k] * shl1_713[k];

        t_1029[k] = pb_z[k] * shl0_714[k]
                    + f_17 * shk_570[k]
                    - f_14 * pc_z[k] * shl1_714[k];

        t_1030[k] = pb_z[k] * shl0_715[k]
                    + f_18 * shk_571[k]
                    - f_14 * pc_z[k] * shl1_715[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pb_z, pc_y, pc_z, shl0_716, shl0_717, \
                         shk_572, shk_573, shk_611, shl1_716, shl1_717, \
                         sik_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = pb_z[k] * shl0_716[k]
                    + f_19 * shk_572[k]
                    - f_14 * pc_z[k] * shl1_716[k];

        t_1032[k] = pb_z[k] * shl0_717[k]
                    + f_0 * shk_573[k]
                    - f_14 * pc_z[k] * shl1_717[k];

        t_1033[k] = f_19 * shk_611[k]
                    + f_3 * pc_y[k] * sik_827[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t shk, const size_t sii0,
                                                          const size_t sii1, const size_t sik,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
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
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shk_575 = buffer.data(shk + 575);
    const auto *shk_576 = buffer.data(shk + 576);
    const auto *shk_579 = buffer.data(shk + 579);
    const auto *shk_582 = buffer.data(shk + 582);
    const auto *shk_586 = buffer.data(shk + 586);
    const auto *shk_591 = buffer.data(shk + 591);
    const auto *shk_604 = buffer.data(shk + 604);
    const auto *shk_611 = buffer.data(shk + 611);
    const auto *shk_612 = buffer.data(shk + 612);
    const auto *shk_614 = buffer.data(shk + 614);
    const auto *shk_615 = buffer.data(shk + 615);
    const auto *shk_617 = buffer.data(shk + 617);
    const auto *shk_618 = buffer.data(shk + 618);
    const auto *shk_621 = buffer.data(shk + 621);
    const auto *shk_622 = buffer.data(shk + 622);
    const auto *shk_626 = buffer.data(shk + 626);
    const auto *shk_627 = buffer.data(shk + 627);
    const auto *shk_632 = buffer.data(shk + 632);
    const auto *shk_640 = buffer.data(shk + 640);
    const auto *shk_642 = buffer.data(shk + 642);
    const auto *shk_643 = buffer.data(shk + 643);
    const auto *shk_644 = buffer.data(shk + 644);
    const auto *shk_645 = buffer.data(shk + 645);
    const auto *shk_646 = buffer.data(shk + 646);
    const auto *shk_647 = buffer.data(shk + 647);
    const auto *shk_648 = buffer.data(shk + 648);
    const auto *shk_650 = buffer.data(shk + 650);
    const auto *shk_651 = buffer.data(shk + 651);
    const auto *shk_653 = buffer.data(shk + 653);
    const auto *shk_654 = buffer.data(shk + 654);
    const auto *shk_657 = buffer.data(shk + 657);
    const auto *shk_658 = buffer.data(shk + 658);
    const auto *shk_662 = buffer.data(shk + 662);
    const auto *shk_663 = buffer.data(shk + 663);
    const auto *shk_668 = buffer.data(shk + 668);
    const auto *shk_676 = buffer.data(shk + 676);
    const auto *shk_678 = buffer.data(shk + 678);
    const auto *shk_679 = buffer.data(shk + 679);
    const auto *shk_680 = buffer.data(shk + 680);
    const auto *shk_681 = buffer.data(shk + 681);
    const auto *shk_682 = buffer.data(shk + 682);
    const auto *shk_683 = buffer.data(shk + 683);
    const auto *shk_684 = buffer.data(shk + 684);
    const auto *shk_686 = buffer.data(shk + 686);
    const auto *shk_689 = buffer.data(shk + 689);
    const auto *shk_693 = buffer.data(shk + 693);
    const auto *shk_698 = buffer.data(shk + 698);

    const auto *sii0_643 = buffer.data(sii0 + 643);
    const auto *sii0_644 = buffer.data(sii0 + 644);
    const auto *sii0_647 = buffer.data(sii0 + 647);
    const auto *sii0_649 = buffer.data(sii0 + 649);
    const auto *sii0_650 = buffer.data(sii0 + 650);
    const auto *sii0_653 = buffer.data(sii0 + 653);
    const auto *sii0_654 = buffer.data(sii0 + 654);
    const auto *sii0_656 = buffer.data(sii0 + 656);
    const auto *sii0_658 = buffer.data(sii0 + 658);
    const auto *sii0_659 = buffer.data(sii0 + 659);
    const auto *sii0_661 = buffer.data(sii0 + 661);
    const auto *sii0_662 = buffer.data(sii0 + 662);
    const auto *sii0_664 = buffer.data(sii0 + 664);
    const auto *sii0_665 = buffer.data(sii0 + 665);
    const auto *sii0_667 = buffer.data(sii0 + 667);
    const auto *sii0_668 = buffer.data(sii0 + 668);
    const auto *sii0_669 = buffer.data(sii0 + 669);
    const auto *sii0_670 = buffer.data(sii0 + 670);
    const auto *sii0_671 = buffer.data(sii0 + 671);
    const auto *sii0_672 = buffer.data(sii0 + 672);
    const auto *sii0_675 = buffer.data(sii0 + 675);
    const auto *sii0_677 = buffer.data(sii0 + 677);
    const auto *sii0_678 = buffer.data(sii0 + 678);
    const auto *sii0_681 = buffer.data(sii0 + 681);
    const auto *sii0_682 = buffer.data(sii0 + 682);
    const auto *sii0_684 = buffer.data(sii0 + 684);
    const auto *sii0_686 = buffer.data(sii0 + 686);
    const auto *sii0_687 = buffer.data(sii0 + 687);
    const auto *sii0_689 = buffer.data(sii0 + 689);
    const auto *sii0_690 = buffer.data(sii0 + 690);
    const auto *sii0_692 = buffer.data(sii0 + 692);
    const auto *sii0_693 = buffer.data(sii0 + 693);
    const auto *sii0_695 = buffer.data(sii0 + 695);
    const auto *sii0_696 = buffer.data(sii0 + 696);
    const auto *sii0_697 = buffer.data(sii0 + 697);
    const auto *sii0_698 = buffer.data(sii0 + 698);
    const auto *sii0_699 = buffer.data(sii0 + 699);
    const auto *sii0_700 = buffer.data(sii0 + 700);
    const auto *sii0_703 = buffer.data(sii0 + 703);
    const auto *sii0_705 = buffer.data(sii0 + 705);
    const auto *sii0_706 = buffer.data(sii0 + 706);
    const auto *sii0_709 = buffer.data(sii0 + 709);
    const auto *sii0_710 = buffer.data(sii0 + 710);
    const auto *sii0_712 = buffer.data(sii0 + 712);
    const auto *sii0_714 = buffer.data(sii0 + 714);
    const auto *sii0_715 = buffer.data(sii0 + 715);
    const auto *sii0_717 = buffer.data(sii0 + 717);
    const auto *sii0_718 = buffer.data(sii0 + 718);
    const auto *sii0_720 = buffer.data(sii0 + 720);
    const auto *sii0_721 = buffer.data(sii0 + 721);
    const auto *sii0_723 = buffer.data(sii0 + 723);

    const auto *sii1_643 = buffer.data(sii1 + 643);
    const auto *sii1_644 = buffer.data(sii1 + 644);
    const auto *sii1_647 = buffer.data(sii1 + 647);
    const auto *sii1_649 = buffer.data(sii1 + 649);
    const auto *sii1_650 = buffer.data(sii1 + 650);
    const auto *sii1_653 = buffer.data(sii1 + 653);
    const auto *sii1_654 = buffer.data(sii1 + 654);
    const auto *sii1_656 = buffer.data(sii1 + 656);
    const auto *sii1_658 = buffer.data(sii1 + 658);
    const auto *sii1_659 = buffer.data(sii1 + 659);
    const auto *sii1_661 = buffer.data(sii1 + 661);
    const auto *sii1_662 = buffer.data(sii1 + 662);
    const auto *sii1_664 = buffer.data(sii1 + 664);
    const auto *sii1_665 = buffer.data(sii1 + 665);
    const auto *sii1_667 = buffer.data(sii1 + 667);
    const auto *sii1_668 = buffer.data(sii1 + 668);
    const auto *sii1_669 = buffer.data(sii1 + 669);
    const auto *sii1_670 = buffer.data(sii1 + 670);
    const auto *sii1_671 = buffer.data(sii1 + 671);
    const auto *sii1_672 = buffer.data(sii1 + 672);
    const auto *sii1_675 = buffer.data(sii1 + 675);
    const auto *sii1_677 = buffer.data(sii1 + 677);
    const auto *sii1_678 = buffer.data(sii1 + 678);
    const auto *sii1_681 = buffer.data(sii1 + 681);
    const auto *sii1_682 = buffer.data(sii1 + 682);
    const auto *sii1_684 = buffer.data(sii1 + 684);
    const auto *sii1_686 = buffer.data(sii1 + 686);
    const auto *sii1_687 = buffer.data(sii1 + 687);
    const auto *sii1_689 = buffer.data(sii1 + 689);
    const auto *sii1_690 = buffer.data(sii1 + 690);
    const auto *sii1_692 = buffer.data(sii1 + 692);
    const auto *sii1_693 = buffer.data(sii1 + 693);
    const auto *sii1_695 = buffer.data(sii1 + 695);
    const auto *sii1_696 = buffer.data(sii1 + 696);
    const auto *sii1_697 = buffer.data(sii1 + 697);
    const auto *sii1_698 = buffer.data(sii1 + 698);
    const auto *sii1_699 = buffer.data(sii1 + 699);
    const auto *sii1_700 = buffer.data(sii1 + 700);
    const auto *sii1_703 = buffer.data(sii1 + 703);
    const auto *sii1_705 = buffer.data(sii1 + 705);
    const auto *sii1_706 = buffer.data(sii1 + 706);
    const auto *sii1_709 = buffer.data(sii1 + 709);
    const auto *sii1_710 = buffer.data(sii1 + 710);
    const auto *sii1_712 = buffer.data(sii1 + 712);
    const auto *sii1_714 = buffer.data(sii1 + 714);
    const auto *sii1_715 = buffer.data(sii1 + 715);
    const auto *sii1_717 = buffer.data(sii1 + 717);
    const auto *sii1_718 = buffer.data(sii1 + 718);
    const auto *sii1_720 = buffer.data(sii1 + 720);
    const auto *sii1_721 = buffer.data(sii1 + 721);
    const auto *sii1_723 = buffer.data(sii1 + 723);

    const auto *sik_827 = buffer.data(sik + 827);
    const auto *sik_828 = buffer.data(sik + 828);
    const auto *sik_830 = buffer.data(sik + 830);
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
    const auto *sik_866 = buffer.data(sik + 866);
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
    const auto *sik_900 = buffer.data(sik + 900);
    const auto *sik_902 = buffer.data(sik + 902);
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

#pragma omp simd aligned(t_1034, t_1035, t_1036, t_1037, pc_x, pc_y, pc_z, shk_575, shk_576, \
                         shk_612, sii0_643, sii0_644, sii1_643, sii1_644, sik_827, \
                         sik_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_15 * shk_575[k]
                    + f_1 * sii0_643[k]
                    - f_2 * sii1_643[k]
                    + f_3 * pc_z[k] * sik_827[k];

        t_1035[k] = f_1 * sii0_644[k]
                    - f_2 * sii1_644[k]
                    + f_3 * pc_x[k] * sik_828[k];

        t_1036[k] = f_18 * shk_612[k]
                    + f_3 * pc_y[k] * sik_828[k];

        t_1037[k] = f_16 * shk_576[k]
                    + f_3 * pc_z[k] * sik_828[k];
    }

#pragma omp simd aligned(t_1038, t_1039, t_1040, pc_x, pc_y, shk_614, sii0_647, sii0_649, \
                         sii1_647, sii1_649, sik_830, sik_831, \
                         sik_833 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1038[k] = f_4 * sii0_647[k]
                    - f_5 * sii1_647[k]
                    + f_3 * pc_x[k] * sik_831[k];

        t_1039[k] = f_18 * shk_614[k]
                    + f_3 * pc_y[k] * sik_830[k];

        t_1040[k] = f_4 * sii0_649[k]
                    - f_5 * sii1_649[k]
                    + f_3 * pc_x[k] * sik_833[k];
    }

#pragma omp simd aligned(t_1041, t_1042, t_1043, pc_x, pc_y, pc_z, shk_579, shk_617, sii0_650, \
                         sii1_650, sik_831, sik_833, sik_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1041[k] = f_6 * sii0_650[k]
                    - f_7 * sii1_650[k]
                    + f_3 * pc_x[k] * sik_834[k];

        t_1042[k] = f_16 * shk_579[k]
                    + f_3 * pc_z[k] * sik_831[k];

        t_1043[k] = f_18 * shk_617[k]
                    + f_3 * pc_y[k] * sik_833[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pc_x, pc_z, shk_582, sii0_653, sii0_654, \
                         sii1_653, sii1_654, sik_834, sik_837, \
                         sik_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_6 * sii0_653[k]
                    - f_7 * sii1_653[k]
                    + f_3 * pc_x[k] * sik_837[k];

        t_1045[k] = f_8 * sii0_654[k]
                    - f_9 * sii1_654[k]
                    + f_3 * pc_x[k] * sik_838[k];

        t_1046[k] = f_16 * shk_582[k]
                    + f_3 * pc_z[k] * sik_834[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pc_x, pc_y, shk_621, sii0_656, sii0_658, \
                         sii1_656, sii1_658, sik_837, sik_840, \
                         sik_842 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_8 * sii0_656[k]
                    - f_9 * sii1_656[k]
                    + f_3 * pc_x[k] * sik_840[k];

        t_1048[k] = f_18 * shk_621[k]
                    + f_3 * pc_y[k] * sik_837[k];

        t_1049[k] = f_8 * sii0_658[k]
                    - f_9 * sii1_658[k]
                    + f_3 * pc_x[k] * sik_842[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, pc_x, pc_z, shk_586, sii0_659, sii0_661, \
                         sii1_659, sii1_661, sik_838, sik_843, \
                         sik_845 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_10 * sii0_659[k]
                    - f_11 * sii1_659[k]
                    + f_3 * pc_x[k] * sik_843[k];

        t_1051[k] = f_16 * shk_586[k]
                    + f_3 * pc_z[k] * sik_838[k];

        t_1052[k] = f_10 * sii0_661[k]
                    - f_11 * sii1_661[k]
                    + f_3 * pc_x[k] * sik_845[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, pc_x, pc_y, shk_626, sii0_662, sii0_664, \
                         sii1_662, sii1_664, sik_842, sik_846, \
                         sik_848 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = f_10 * sii0_662[k]
                    - f_11 * sii1_662[k]
                    + f_3 * pc_x[k] * sik_846[k];

        t_1054[k] = f_18 * shk_626[k]
                    + f_3 * pc_y[k] * sik_842[k];

        t_1055[k] = f_10 * sii0_664[k]
                    - f_11 * sii1_664[k]
                    + f_3 * pc_x[k] * sik_848[k];
    }

#pragma omp simd aligned(t_1056, t_1057, t_1058, pc_x, pc_z, shk_591, sii0_665, sii0_667, \
                         sii1_665, sii1_667, sik_843, sik_849, \
                         sik_851 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1056[k] = f_12 * sii0_665[k]
                    - f_13 * sii1_665[k]
                    + f_3 * pc_x[k] * sik_849[k];

        t_1057[k] = f_16 * shk_591[k]
                    + f_3 * pc_z[k] * sik_843[k];

        t_1058[k] = f_12 * sii0_667[k]
                    - f_13 * sii1_667[k]
                    + f_3 * pc_x[k] * sik_851[k];
    }

#pragma omp simd aligned(t_1059, t_1060, t_1061, pc_x, pc_y, shk_632, sii0_668, sii0_669, \
                         sii1_668, sii1_669, sik_848, sik_852, \
                         sik_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1059[k] = f_12 * sii0_668[k]
                    - f_13 * sii1_668[k]
                    + f_3 * pc_x[k] * sik_852[k];

        t_1060[k] = f_12 * sii0_669[k]
                    - f_13 * sii1_669[k]
                    + f_3 * pc_x[k] * sik_853[k];

        t_1061[k] = f_18 * shk_632[k]
                    + f_3 * pc_y[k] * sik_848[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, t_1065, t_1066, t_1067, pc_x, sii0_671, \
                         sii1_671, sik_855, sik_856, sik_857, sik_858, sik_859, \
                         sik_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_12 * sii0_671[k]
                    - f_13 * sii1_671[k]
                    + f_3 * pc_x[k] * sik_855[k];

        t_1063[k] = f_3 * pc_x[k] * sik_856[k];

        t_1064[k] = f_3 * pc_x[k] * sik_857[k];

        t_1065[k] = f_3 * pc_x[k] * sik_858[k];

        t_1066[k] = f_3 * pc_x[k] * sik_859[k];

        t_1067[k] = f_3 * pc_x[k] * sik_860[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, shk_604, \
                         shk_640, sii0_665, sii1_665, sik_856, sik_861, sik_862, \
                         sik_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_3 * pc_x[k] * sik_861[k];

        t_1069[k] = f_3 * pc_x[k] * sik_862[k];

        t_1070[k] = f_3 * pc_x[k] * sik_863[k];

        t_1071[k] = f_18 * shk_640[k]
                    + f_1 * sii0_665[k]
                    - f_2 * sii1_665[k]
                    + f_3 * pc_y[k] * sik_856[k];

        t_1072[k] = f_16 * shk_604[k]
                    + f_3 * pc_z[k] * sik_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_y, shk_642, shk_643, shk_644, sii0_667, \
                         sii0_668, sii0_669, sii1_667, sii1_668, sii1_669, sik_858, sik_859, \
                         sik_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_18 * shk_642[k]
                    + f_4 * sii0_667[k]
                    - f_5 * sii1_667[k]
                    + f_3 * pc_y[k] * sik_858[k];

        t_1074[k] = f_18 * shk_643[k]
                    + f_6 * sii0_668[k]
                    - f_7 * sii1_668[k]
                    + f_3 * pc_y[k] * sik_859[k];

        t_1075[k] = f_18 * shk_644[k]
                    + f_8 * sii0_669[k]
                    - f_9 * sii1_669[k]
                    + f_3 * pc_y[k] * sik_860[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_y, shk_645, shk_646, shk_647, sii0_670, \
                         sii0_671, sii1_670, sii1_671, sik_861, sik_862, \
                         sik_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_18 * shk_645[k]
                    + f_10 * sii0_670[k]
                    - f_11 * sii1_670[k]
                    + f_3 * pc_y[k] * sik_861[k];

        t_1077[k] = f_18 * shk_646[k]
                    + f_12 * sii0_671[k]
                    - f_13 * sii1_671[k]
                    + f_3 * pc_y[k] * sik_862[k];

        t_1078[k] = f_18 * shk_647[k]
                    + f_3 * pc_y[k] * sik_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, t_1082, pc_x, pc_y, pc_z, shk_611, shk_612, \
                         shk_648, sii0_671, sii0_672, sii1_671, sii1_672, sik_863, \
                         sik_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_16 * shk_611[k]
                    + f_1 * sii0_671[k]
                    - f_2 * sii1_671[k]
                    + f_3 * pc_z[k] * sik_863[k];

        t_1080[k] = f_1 * sii0_672[k]
                    - f_2 * sii1_672[k]
                    + f_3 * pc_x[k] * sik_864[k];

        t_1081[k] = f_17 * shk_648[k]
                    + f_3 * pc_y[k] * sik_864[k];

        t_1082[k] = f_17 * shk_612[k]
                    + f_3 * pc_z[k] * sik_864[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, pc_x, pc_y, shk_650, sii0_675, sii0_677, \
                         sii1_675, sii1_677, sik_866, sik_867, \
                         sik_869 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_4 * sii0_675[k]
                    - f_5 * sii1_675[k]
                    + f_3 * pc_x[k] * sik_867[k];

        t_1084[k] = f_17 * shk_650[k]
                    + f_3 * pc_y[k] * sik_866[k];

        t_1085[k] = f_4 * sii0_677[k]
                    - f_5 * sii1_677[k]
                    + f_3 * pc_x[k] * sik_869[k];
    }

#pragma omp simd aligned(t_1086, t_1087, t_1088, pc_x, pc_y, pc_z, shk_615, shk_653, sii0_678, \
                         sii1_678, sik_867, sik_869, sik_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1086[k] = f_6 * sii0_678[k]
                    - f_7 * sii1_678[k]
                    + f_3 * pc_x[k] * sik_870[k];

        t_1087[k] = f_17 * shk_615[k]
                    + f_3 * pc_z[k] * sik_867[k];

        t_1088[k] = f_17 * shk_653[k]
                    + f_3 * pc_y[k] * sik_869[k];
    }

#pragma omp simd aligned(t_1089, t_1090, t_1091, pc_x, pc_z, shk_618, sii0_681, sii0_682, \
                         sii1_681, sii1_682, sik_870, sik_873, \
                         sik_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1089[k] = f_6 * sii0_681[k]
                    - f_7 * sii1_681[k]
                    + f_3 * pc_x[k] * sik_873[k];

        t_1090[k] = f_8 * sii0_682[k]
                    - f_9 * sii1_682[k]
                    + f_3 * pc_x[k] * sik_874[k];

        t_1091[k] = f_17 * shk_618[k]
                    + f_3 * pc_z[k] * sik_870[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_x, pc_y, shk_657, sii0_684, sii0_686, \
                         sii1_684, sii1_686, sik_873, sik_876, \
                         sik_878 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_8 * sii0_684[k]
                    - f_9 * sii1_684[k]
                    + f_3 * pc_x[k] * sik_876[k];

        t_1093[k] = f_17 * shk_657[k]
                    + f_3 * pc_y[k] * sik_873[k];

        t_1094[k] = f_8 * sii0_686[k]
                    - f_9 * sii1_686[k]
                    + f_3 * pc_x[k] * sik_878[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, pc_z, shk_622, sii0_687, sii0_689, \
                         sii1_687, sii1_689, sik_874, sik_879, \
                         sik_881 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_10 * sii0_687[k]
                    - f_11 * sii1_687[k]
                    + f_3 * pc_x[k] * sik_879[k];

        t_1096[k] = f_17 * shk_622[k]
                    + f_3 * pc_z[k] * sik_874[k];

        t_1097[k] = f_10 * sii0_689[k]
                    - f_11 * sii1_689[k]
                    + f_3 * pc_x[k] * sik_881[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, pc_x, pc_y, shk_662, sii0_690, sii0_692, \
                         sii1_690, sii1_692, sik_878, sik_882, \
                         sik_884 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_10 * sii0_690[k]
                    - f_11 * sii1_690[k]
                    + f_3 * pc_x[k] * sik_882[k];

        t_1099[k] = f_17 * shk_662[k]
                    + f_3 * pc_y[k] * sik_878[k];

        t_1100[k] = f_10 * sii0_692[k]
                    - f_11 * sii1_692[k]
                    + f_3 * pc_x[k] * sik_884[k];
    }

#pragma omp simd aligned(t_1101, t_1102, t_1103, pc_x, pc_z, shk_627, sii0_693, sii0_695, \
                         sii1_693, sii1_695, sik_879, sik_885, \
                         sik_887 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1101[k] = f_12 * sii0_693[k]
                    - f_13 * sii1_693[k]
                    + f_3 * pc_x[k] * sik_885[k];

        t_1102[k] = f_17 * shk_627[k]
                    + f_3 * pc_z[k] * sik_879[k];

        t_1103[k] = f_12 * sii0_695[k]
                    - f_13 * sii1_695[k]
                    + f_3 * pc_x[k] * sik_887[k];
    }

#pragma omp simd aligned(t_1104, t_1105, t_1106, pc_x, pc_y, shk_668, sii0_696, sii0_697, \
                         sii1_696, sii1_697, sik_884, sik_888, \
                         sik_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1104[k] = f_12 * sii0_696[k]
                    - f_13 * sii1_696[k]
                    + f_3 * pc_x[k] * sik_888[k];

        t_1105[k] = f_12 * sii0_697[k]
                    - f_13 * sii1_697[k]
                    + f_3 * pc_x[k] * sik_889[k];

        t_1106[k] = f_17 * shk_668[k]
                    + f_3 * pc_y[k] * sik_884[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, t_1110, t_1111, t_1112, pc_x, sii0_699, \
                         sii1_699, sik_891, sik_892, sik_893, sik_894, sik_895, \
                         sik_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_12 * sii0_699[k]
                    - f_13 * sii1_699[k]
                    + f_3 * pc_x[k] * sik_891[k];

        t_1108[k] = f_3 * pc_x[k] * sik_892[k];

        t_1109[k] = f_3 * pc_x[k] * sik_893[k];

        t_1110[k] = f_3 * pc_x[k] * sik_894[k];

        t_1111[k] = f_3 * pc_x[k] * sik_895[k];

        t_1112[k] = f_3 * pc_x[k] * sik_896[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, shk_640, \
                         shk_676, sii0_693, sii1_693, sik_892, sik_897, sik_898, \
                         sik_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = f_3 * pc_x[k] * sik_897[k];

        t_1114[k] = f_3 * pc_x[k] * sik_898[k];

        t_1115[k] = f_3 * pc_x[k] * sik_899[k];

        t_1116[k] = f_17 * shk_676[k]
                    + f_1 * sii0_693[k]
                    - f_2 * sii1_693[k]
                    + f_3 * pc_y[k] * sik_892[k];

        t_1117[k] = f_17 * shk_640[k]
                    + f_3 * pc_z[k] * sik_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, shk_678, shk_679, shk_680, sii0_695, \
                         sii0_696, sii0_697, sii1_695, sii1_696, sii1_697, sik_894, sik_895, \
                         sik_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * shk_678[k]
                    + f_4 * sii0_695[k]
                    - f_5 * sii1_695[k]
                    + f_3 * pc_y[k] * sik_894[k];

        t_1119[k] = f_17 * shk_679[k]
                    + f_6 * sii0_696[k]
                    - f_7 * sii1_696[k]
                    + f_3 * pc_y[k] * sik_895[k];

        t_1120[k] = f_17 * shk_680[k]
                    + f_8 * sii0_697[k]
                    - f_9 * sii1_697[k]
                    + f_3 * pc_y[k] * sik_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, shk_681, shk_682, shk_683, sii0_698, \
                         sii0_699, sii1_698, sii1_699, sik_897, sik_898, \
                         sik_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * shk_681[k]
                    + f_10 * sii0_698[k]
                    - f_11 * sii1_698[k]
                    + f_3 * pc_y[k] * sik_897[k];

        t_1122[k] = f_17 * shk_682[k]
                    + f_12 * sii0_699[k]
                    - f_13 * sii1_699[k]
                    + f_3 * pc_y[k] * sik_898[k];

        t_1123[k] = f_17 * shk_683[k]
                    + f_3 * pc_y[k] * sik_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, t_1127, pc_x, pc_y, pc_z, shk_647, shk_648, \
                         shk_684, sii0_699, sii0_700, sii1_699, sii1_700, sik_899, \
                         sik_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * shk_647[k]
                    + f_1 * sii0_699[k]
                    - f_2 * sii1_699[k]
                    + f_3 * pc_z[k] * sik_899[k];

        t_1125[k] = f_1 * sii0_700[k]
                    - f_2 * sii1_700[k]
                    + f_3 * pc_x[k] * sik_900[k];

        t_1126[k] = f_16 * shk_684[k]
                    + f_3 * pc_y[k] * sik_900[k];

        t_1127[k] = f_18 * shk_648[k]
                    + f_3 * pc_z[k] * sik_900[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, pc_x, pc_y, shk_686, sii0_703, sii0_705, \
                         sii1_703, sii1_705, sik_902, sik_903, \
                         sik_905 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_4 * sii0_703[k]
                    - f_5 * sii1_703[k]
                    + f_3 * pc_x[k] * sik_903[k];

        t_1129[k] = f_16 * shk_686[k]
                    + f_3 * pc_y[k] * sik_902[k];

        t_1130[k] = f_4 * sii0_705[k]
                    - f_5 * sii1_705[k]
                    + f_3 * pc_x[k] * sik_905[k];
    }

#pragma omp simd aligned(t_1131, t_1132, t_1133, pc_x, pc_y, pc_z, shk_651, shk_689, sii0_706, \
                         sii1_706, sik_903, sik_905, sik_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1131[k] = f_6 * sii0_706[k]
                    - f_7 * sii1_706[k]
                    + f_3 * pc_x[k] * sik_906[k];

        t_1132[k] = f_18 * shk_651[k]
                    + f_3 * pc_z[k] * sik_903[k];

        t_1133[k] = f_16 * shk_689[k]
                    + f_3 * pc_y[k] * sik_905[k];
    }

#pragma omp simd aligned(t_1134, t_1135, t_1136, pc_x, pc_z, shk_654, sii0_709, sii0_710, \
                         sii1_709, sii1_710, sik_906, sik_909, \
                         sik_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1134[k] = f_6 * sii0_709[k]
                    - f_7 * sii1_709[k]
                    + f_3 * pc_x[k] * sik_909[k];

        t_1135[k] = f_8 * sii0_710[k]
                    - f_9 * sii1_710[k]
                    + f_3 * pc_x[k] * sik_910[k];

        t_1136[k] = f_18 * shk_654[k]
                    + f_3 * pc_z[k] * sik_906[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, pc_x, pc_y, shk_693, sii0_712, sii0_714, \
                         sii1_712, sii1_714, sik_909, sik_912, \
                         sik_914 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_8 * sii0_712[k]
                    - f_9 * sii1_712[k]
                    + f_3 * pc_x[k] * sik_912[k];

        t_1138[k] = f_16 * shk_693[k]
                    + f_3 * pc_y[k] * sik_909[k];

        t_1139[k] = f_8 * sii0_714[k]
                    - f_9 * sii1_714[k]
                    + f_3 * pc_x[k] * sik_914[k];
    }

#pragma omp simd aligned(t_1140, t_1141, t_1142, pc_x, pc_z, shk_658, sii0_715, sii0_717, \
                         sii1_715, sii1_717, sik_910, sik_915, \
                         sik_917 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1140[k] = f_10 * sii0_715[k]
                    - f_11 * sii1_715[k]
                    + f_3 * pc_x[k] * sik_915[k];

        t_1141[k] = f_18 * shk_658[k]
                    + f_3 * pc_z[k] * sik_910[k];

        t_1142[k] = f_10 * sii0_717[k]
                    - f_11 * sii1_717[k]
                    + f_3 * pc_x[k] * sik_917[k];
    }

#pragma omp simd aligned(t_1143, t_1144, t_1145, pc_x, pc_y, shk_698, sii0_718, sii0_720, \
                         sii1_718, sii1_720, sik_914, sik_918, \
                         sik_920 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1143[k] = f_10 * sii0_718[k]
                    - f_11 * sii1_718[k]
                    + f_3 * pc_x[k] * sik_918[k];

        t_1144[k] = f_16 * shk_698[k]
                    + f_3 * pc_y[k] * sik_914[k];

        t_1145[k] = f_10 * sii0_720[k]
                    - f_11 * sii1_720[k]
                    + f_3 * pc_x[k] * sik_920[k];
    }

#pragma omp simd aligned(t_1146, t_1147, t_1148, pc_x, pc_z, shk_663, sii0_721, sii0_723, \
                         sii1_721, sii1_723, sik_915, sik_921, \
                         sik_923 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1146[k] = f_12 * sii0_721[k]
                    - f_13 * sii1_721[k]
                    + f_3 * pc_x[k] * sik_921[k];

        t_1147[k] = f_18 * shk_663[k]
                    + f_3 * pc_z[k] * sik_915[k];

        t_1148[k] = f_12 * sii0_723[k]
                    - f_13 * sii1_723[k]
                    + f_3 * pc_x[k] * sik_923[k];
    }
}

static auto
compute_prim_sil_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t shl0,
                                                           const size_t shk, const size_t shl1,
                                                           const size_t sii0, const size_t sii1,
                                                           const size_t sik, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
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
    const auto f_20 = 4.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shl0_900 = buffer.data(shl0 + 900);
    const auto *shl0_905 = buffer.data(shl0 + 905);
    const auto *shl0_909 = buffer.data(shl0 + 909);
    const auto *shl0_914 = buffer.data(shl0 + 914);
    const auto *shl0_920 = buffer.data(shl0 + 920);
    const auto *shl0_927 = buffer.data(shl0 + 927);
    const auto *shl0_936 = buffer.data(shl0 + 936);
    const auto *shl0_938 = buffer.data(shl0 + 938);
    const auto *shl0_939 = buffer.data(shl0 + 939);
    const auto *shl0_940 = buffer.data(shl0 + 940);
    const auto *shl0_941 = buffer.data(shl0 + 941);
    const auto *shl0_942 = buffer.data(shl0 + 942);
    const auto *shl0_944 = buffer.data(shl0 + 944);

    const auto *shk_676 = buffer.data(shk + 676);
    const auto *shk_683 = buffer.data(shk + 683);
    const auto *shk_684 = buffer.data(shk + 684);
    const auto *shk_687 = buffer.data(shk + 687);
    const auto *shk_690 = buffer.data(shk + 690);
    const auto *shk_694 = buffer.data(shk + 694);
    const auto *shk_699 = buffer.data(shk + 699);
    const auto *shk_704 = buffer.data(shk + 704);
    const auto *shk_712 = buffer.data(shk + 712);
    const auto *shk_714 = buffer.data(shk + 714);
    const auto *shk_715 = buffer.data(shk + 715);
    const auto *shk_716 = buffer.data(shk + 716);
    const auto *shk_717 = buffer.data(shk + 717);
    const auto *shk_718 = buffer.data(shk + 718);
    const auto *shk_719 = buffer.data(shk + 719);
    const auto *shk_720 = buffer.data(shk + 720);
    const auto *shk_722 = buffer.data(shk + 722);
    const auto *shk_723 = buffer.data(shk + 723);
    const auto *shk_725 = buffer.data(shk + 725);
    const auto *shk_726 = buffer.data(shk + 726);
    const auto *shk_729 = buffer.data(shk + 729);
    const auto *shk_730 = buffer.data(shk + 730);
    const auto *shk_734 = buffer.data(shk + 734);
    const auto *shk_735 = buffer.data(shk + 735);
    const auto *shk_740 = buffer.data(shk + 740);
    const auto *shk_748 = buffer.data(shk + 748);
    const auto *shk_750 = buffer.data(shk + 750);
    const auto *shk_751 = buffer.data(shk + 751);
    const auto *shk_752 = buffer.data(shk + 752);
    const auto *shk_753 = buffer.data(shk + 753);
    const auto *shk_754 = buffer.data(shk + 754);
    const auto *shk_755 = buffer.data(shk + 755);

    const auto *shl1_900 = buffer.data(shl1 + 900);
    const auto *shl1_905 = buffer.data(shl1 + 905);
    const auto *shl1_909 = buffer.data(shl1 + 909);
    const auto *shl1_914 = buffer.data(shl1 + 914);
    const auto *shl1_920 = buffer.data(shl1 + 920);
    const auto *shl1_927 = buffer.data(shl1 + 927);
    const auto *shl1_936 = buffer.data(shl1 + 936);
    const auto *shl1_938 = buffer.data(shl1 + 938);
    const auto *shl1_939 = buffer.data(shl1 + 939);
    const auto *shl1_940 = buffer.data(shl1 + 940);
    const auto *shl1_941 = buffer.data(shl1 + 941);
    const auto *shl1_942 = buffer.data(shl1 + 942);
    const auto *shl1_944 = buffer.data(shl1 + 944);

    const auto *sii0_721 = buffer.data(sii0 + 721);
    const auto *sii0_723 = buffer.data(sii0 + 723);
    const auto *sii0_724 = buffer.data(sii0 + 724);
    const auto *sii0_725 = buffer.data(sii0 + 725);
    const auto *sii0_726 = buffer.data(sii0 + 726);
    const auto *sii0_727 = buffer.data(sii0 + 727);
    const auto *sii0_731 = buffer.data(sii0 + 731);
    const auto *sii0_734 = buffer.data(sii0 + 734);
    const auto *sii0_738 = buffer.data(sii0 + 738);
    const auto *sii0_740 = buffer.data(sii0 + 740);
    const auto *sii0_743 = buffer.data(sii0 + 743);
    const auto *sii0_745 = buffer.data(sii0 + 745);
    const auto *sii0_746 = buffer.data(sii0 + 746);
    const auto *sii0_749 = buffer.data(sii0 + 749);
    const auto *sii0_751 = buffer.data(sii0 + 751);
    const auto *sii0_752 = buffer.data(sii0 + 752);
    const auto *sii0_753 = buffer.data(sii0 + 753);
    const auto *sii0_756 = buffer.data(sii0 + 756);
    const auto *sii0_759 = buffer.data(sii0 + 759);
    const auto *sii0_761 = buffer.data(sii0 + 761);
    const auto *sii0_762 = buffer.data(sii0 + 762);
    const auto *sii0_765 = buffer.data(sii0 + 765);
    const auto *sii0_766 = buffer.data(sii0 + 766);
    const auto *sii0_768 = buffer.data(sii0 + 768);
    const auto *sii0_770 = buffer.data(sii0 + 770);
    const auto *sii0_771 = buffer.data(sii0 + 771);
    const auto *sii0_773 = buffer.data(sii0 + 773);
    const auto *sii0_774 = buffer.data(sii0 + 774);
    const auto *sii0_776 = buffer.data(sii0 + 776);
    const auto *sii0_777 = buffer.data(sii0 + 777);
    const auto *sii0_779 = buffer.data(sii0 + 779);
    const auto *sii0_780 = buffer.data(sii0 + 780);
    const auto *sii0_781 = buffer.data(sii0 + 781);
    const auto *sii0_782 = buffer.data(sii0 + 782);
    const auto *sii0_783 = buffer.data(sii0 + 783);

    const auto *sii1_721 = buffer.data(sii1 + 721);
    const auto *sii1_723 = buffer.data(sii1 + 723);
    const auto *sii1_724 = buffer.data(sii1 + 724);
    const auto *sii1_725 = buffer.data(sii1 + 725);
    const auto *sii1_726 = buffer.data(sii1 + 726);
    const auto *sii1_727 = buffer.data(sii1 + 727);
    const auto *sii1_731 = buffer.data(sii1 + 731);
    const auto *sii1_734 = buffer.data(sii1 + 734);
    const auto *sii1_738 = buffer.data(sii1 + 738);
    const auto *sii1_740 = buffer.data(sii1 + 740);
    const auto *sii1_743 = buffer.data(sii1 + 743);
    const auto *sii1_745 = buffer.data(sii1 + 745);
    const auto *sii1_746 = buffer.data(sii1 + 746);
    const auto *sii1_749 = buffer.data(sii1 + 749);
    const auto *sii1_751 = buffer.data(sii1 + 751);
    const auto *sii1_752 = buffer.data(sii1 + 752);
    const auto *sii1_753 = buffer.data(sii1 + 753);
    const auto *sii1_756 = buffer.data(sii1 + 756);
    const auto *sii1_759 = buffer.data(sii1 + 759);
    const auto *sii1_761 = buffer.data(sii1 + 761);
    const auto *sii1_762 = buffer.data(sii1 + 762);
    const auto *sii1_765 = buffer.data(sii1 + 765);
    const auto *sii1_766 = buffer.data(sii1 + 766);
    const auto *sii1_768 = buffer.data(sii1 + 768);
    const auto *sii1_770 = buffer.data(sii1 + 770);
    const auto *sii1_771 = buffer.data(sii1 + 771);
    const auto *sii1_773 = buffer.data(sii1 + 773);
    const auto *sii1_774 = buffer.data(sii1 + 774);
    const auto *sii1_776 = buffer.data(sii1 + 776);
    const auto *sii1_777 = buffer.data(sii1 + 777);
    const auto *sii1_779 = buffer.data(sii1 + 779);
    const auto *sii1_780 = buffer.data(sii1 + 780);
    const auto *sii1_781 = buffer.data(sii1 + 781);
    const auto *sii1_782 = buffer.data(sii1 + 782);
    const auto *sii1_783 = buffer.data(sii1 + 783);

    const auto *sik_920 = buffer.data(sik + 920);
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
    const auto *sik_936 = buffer.data(sik + 936);
    const auto *sik_938 = buffer.data(sik + 938);
    const auto *sik_939 = buffer.data(sik + 939);
    const auto *sik_941 = buffer.data(sik + 941);
    const auto *sik_942 = buffer.data(sik + 942);
    const auto *sik_945 = buffer.data(sik + 945);
    const auto *sik_946 = buffer.data(sik + 946);
    const auto *sik_948 = buffer.data(sik + 948);
    const auto *sik_950 = buffer.data(sik + 950);
    const auto *sik_951 = buffer.data(sik + 951);
    const auto *sik_953 = buffer.data(sik + 953);
    const auto *sik_954 = buffer.data(sik + 954);
    const auto *sik_956 = buffer.data(sik + 956);
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
    const auto *sik_974 = buffer.data(sik + 974);
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

#pragma omp simd aligned(t_1149, t_1150, t_1151, pc_x, pc_y, shk_704, sii0_724, sii0_725, \
                         sii1_724, sii1_725, sik_920, sik_924, \
                         sik_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_12 * sii0_724[k]
                    - f_13 * sii1_724[k]
                    + f_3 * pc_x[k] * sik_924[k];

        t_1150[k] = f_12 * sii0_725[k]
                    - f_13 * sii1_725[k]
                    + f_3 * pc_x[k] * sik_925[k];

        t_1151[k] = f_16 * shk_704[k]
                    + f_3 * pc_y[k] * sik_920[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, t_1155, t_1156, t_1157, pc_x, sii0_727, \
                         sii1_727, sik_927, sik_928, sik_929, sik_930, sik_931, \
                         sik_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = f_12 * sii0_727[k]
                    - f_13 * sii1_727[k]
                    + f_3 * pc_x[k] * sik_927[k];

        t_1153[k] = f_3 * pc_x[k] * sik_928[k];

        t_1154[k] = f_3 * pc_x[k] * sik_929[k];

        t_1155[k] = f_3 * pc_x[k] * sik_930[k];

        t_1156[k] = f_3 * pc_x[k] * sik_931[k];

        t_1157[k] = f_3 * pc_x[k] * sik_932[k];
    }

#pragma omp simd aligned(t_1158, t_1159, t_1160, t_1161, t_1162, pc_x, pc_y, pc_z, shk_676, \
                         shk_712, sii0_721, sii1_721, sik_928, sik_933, sik_934, \
                         sik_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1158[k] = f_3 * pc_x[k] * sik_933[k];

        t_1159[k] = f_3 * pc_x[k] * sik_934[k];

        t_1160[k] = f_3 * pc_x[k] * sik_935[k];

        t_1161[k] = f_16 * shk_712[k]
                    + f_1 * sii0_721[k]
                    - f_2 * sii1_721[k]
                    + f_3 * pc_y[k] * sik_928[k];

        t_1162[k] = f_18 * shk_676[k]
                    + f_3 * pc_z[k] * sik_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_y, shk_714, shk_715, shk_716, sii0_723, \
                         sii0_724, sii0_725, sii1_723, sii1_724, sii1_725, sik_930, sik_931, \
                         sik_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * shk_714[k]
                    + f_4 * sii0_723[k]
                    - f_5 * sii1_723[k]
                    + f_3 * pc_y[k] * sik_930[k];

        t_1164[k] = f_16 * shk_715[k]
                    + f_6 * sii0_724[k]
                    - f_7 * sii1_724[k]
                    + f_3 * pc_y[k] * sik_931[k];

        t_1165[k] = f_16 * shk_716[k]
                    + f_8 * sii0_725[k]
                    - f_9 * sii1_725[k]
                    + f_3 * pc_y[k] * sik_932[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_y, shk_717, shk_718, shk_719, sii0_726, \
                         sii0_727, sii1_726, sii1_727, sik_933, sik_934, \
                         sik_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_16 * shk_717[k]
                    + f_10 * sii0_726[k]
                    - f_11 * sii1_726[k]
                    + f_3 * pc_y[k] * sik_933[k];

        t_1167[k] = f_16 * shk_718[k]
                    + f_12 * sii0_727[k]
                    - f_13 * sii1_727[k]
                    + f_3 * pc_y[k] * sik_934[k];

        t_1168[k] = f_16 * shk_719[k]
                    + f_3 * pc_y[k] * sik_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pb_y, pc_y, pc_z, shl0_900, shk_683, \
                         shk_684, shk_720, shl1_900, sii0_727, sii1_727, sik_935, \
                         sik_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_18 * shk_683[k]
                    + f_1 * sii0_727[k]
                    - f_2 * sii1_727[k]
                    + f_3 * pc_z[k] * sik_935[k];

        t_1170[k] = pb_y[k] * shl0_900[k]
                    - f_14 * pc_y[k] * shl1_900[k];

        t_1171[k] = f_15 * shk_720[k]
                    + f_3 * pc_y[k] * sik_936[k];

        t_1172[k] = f_19 * shk_684[k]
                    + f_3 * pc_z[k] * sik_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, pb_y, pc_x, pc_y, shl0_905, shk_722, \
                         shl1_905, sii0_731, sii1_731, sik_938, \
                         sik_939 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = f_4 * sii0_731[k]
                    - f_5 * sii1_731[k]
                    + f_3 * pc_x[k] * sik_939[k];

        t_1174[k] = f_15 * shk_722[k]
                    + f_3 * pc_y[k] * sik_938[k];

        t_1175[k] = pb_y[k] * shl0_905[k]
                    - f_14 * pc_y[k] * shl1_905[k];
    }

#pragma omp simd aligned(t_1176, t_1177, t_1178, pc_x, pc_y, pc_z, shk_687, shk_725, sii0_734, \
                         sii1_734, sik_939, sik_941, sik_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1176[k] = f_6 * sii0_734[k]
                    - f_7 * sii1_734[k]
                    + f_3 * pc_x[k] * sik_942[k];

        t_1177[k] = f_19 * shk_687[k]
                    + f_3 * pc_z[k] * sik_939[k];

        t_1178[k] = f_15 * shk_725[k]
                    + f_3 * pc_y[k] * sik_941[k];
    }

#pragma omp simd aligned(t_1179, t_1180, t_1181, pb_y, pc_x, pc_y, pc_z, shl0_909, shk_690, \
                         shl1_909, sii0_738, sii1_738, sik_942, \
                         sik_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1179[k] = pb_y[k] * shl0_909[k]
                    - f_14 * pc_y[k] * shl1_909[k];

        t_1180[k] = f_8 * sii0_738[k]
                    - f_9 * sii1_738[k]
                    + f_3 * pc_x[k] * sik_946[k];

        t_1181[k] = f_19 * shk_690[k]
                    + f_3 * pc_z[k] * sik_942[k];
    }

#pragma omp simd aligned(t_1182, t_1183, t_1184, pb_y, pc_x, pc_y, shl0_914, shk_729, \
                         shl1_914, sii0_740, sii1_740, sik_945, \
                         sik_948 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1182[k] = f_8 * sii0_740[k]
                    - f_9 * sii1_740[k]
                    + f_3 * pc_x[k] * sik_948[k];

        t_1183[k] = f_15 * shk_729[k]
                    + f_3 * pc_y[k] * sik_945[k];

        t_1184[k] = pb_y[k] * shl0_914[k]
                    - f_14 * pc_y[k] * shl1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pc_x, pc_z, shk_694, sii0_743, sii0_745, \
                         sii1_743, sii1_745, sik_946, sik_951, \
                         sik_953 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = f_10 * sii0_743[k]
                    - f_11 * sii1_743[k]
                    + f_3 * pc_x[k] * sik_951[k];

        t_1186[k] = f_19 * shk_694[k]
                    + f_3 * pc_z[k] * sik_946[k];

        t_1187[k] = f_10 * sii0_745[k]
                    - f_11 * sii1_745[k]
                    + f_3 * pc_x[k] * sik_953[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, pb_y, pc_x, pc_y, shl0_920, shk_734, \
                         shl1_920, sii0_746, sii1_746, sik_950, \
                         sik_954 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = f_10 * sii0_746[k]
                    - f_11 * sii1_746[k]
                    + f_3 * pc_x[k] * sik_954[k];

        t_1189[k] = f_15 * shk_734[k]
                    + f_3 * pc_y[k] * sik_950[k];

        t_1190[k] = pb_y[k] * shl0_920[k]
                    - f_14 * pc_y[k] * shl1_920[k];
    }

#pragma omp simd aligned(t_1191, t_1192, t_1193, pc_x, pc_z, shk_699, sii0_749, sii0_751, \
                         sii1_749, sii1_751, sik_951, sik_957, \
                         sik_959 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1191[k] = f_12 * sii0_749[k]
                    - f_13 * sii1_749[k]
                    + f_3 * pc_x[k] * sik_957[k];

        t_1192[k] = f_19 * shk_699[k]
                    + f_3 * pc_z[k] * sik_951[k];

        t_1193[k] = f_12 * sii0_751[k]
                    - f_13 * sii1_751[k]
                    + f_3 * pc_x[k] * sik_959[k];
    }

#pragma omp simd aligned(t_1194, t_1195, t_1196, pc_x, pc_y, shk_740, sii0_752, sii0_753, \
                         sii1_752, sii1_753, sik_956, sik_960, \
                         sik_961 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1194[k] = f_12 * sii0_752[k]
                    - f_13 * sii1_752[k]
                    + f_3 * pc_x[k] * sik_960[k];

        t_1195[k] = f_12 * sii0_753[k]
                    - f_13 * sii1_753[k]
                    + f_3 * pc_x[k] * sik_961[k];

        t_1196[k] = f_15 * shk_740[k]
                    + f_3 * pc_y[k] * sik_956[k];
    }

#pragma omp simd aligned(t_1197, t_1198, t_1199, t_1200, t_1201, t_1202, pb_y, pc_x, pc_y, \
                         shl0_927, shl1_927, sik_964, sik_965, sik_966, sik_967, \
                         sik_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1197[k] = pb_y[k] * shl0_927[k]
                    - f_14 * pc_y[k] * shl1_927[k];

        t_1198[k] = f_3 * pc_x[k] * sik_964[k];

        t_1199[k] = f_3 * pc_x[k] * sik_965[k];

        t_1200[k] = f_3 * pc_x[k] * sik_966[k];

        t_1201[k] = f_3 * pc_x[k] * sik_967[k];

        t_1202[k] = f_3 * pc_x[k] * sik_968[k];
    }

#pragma omp simd aligned(t_1203, t_1204, t_1205, t_1206, pb_y, pc_x, pc_y, shl0_936, shk_748, \
                         shl1_936, sik_969, sik_970, sik_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1203[k] = f_3 * pc_x[k] * sik_969[k];

        t_1204[k] = f_3 * pc_x[k] * sik_970[k];

        t_1205[k] = f_3 * pc_x[k] * sik_971[k];

        t_1206[k] = pb_y[k] * shl0_936[k]
                    + f_20 * shk_748[k]
                    - f_14 * pc_y[k] * shl1_936[k];
    }

#pragma omp simd aligned(t_1207, t_1208, t_1209, pb_y, pc_y, pc_z, shl0_938, shl0_939, \
                         shk_712, shk_750, shk_751, shl1_938, shl1_939, \
                         sik_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1207[k] = f_19 * shk_712[k]
                    + f_3 * pc_z[k] * sik_964[k];

        t_1208[k] = pb_y[k] * shl0_938[k]
                    + f_0 * shk_750[k]
                    - f_14 * pc_y[k] * shl1_938[k];

        t_1209[k] = pb_y[k] * shl0_939[k]
                    + f_19 * shk_751[k]
                    - f_14 * pc_y[k] * shl1_939[k];
    }

#pragma omp simd aligned(t_1210, t_1211, t_1212, pb_y, pc_y, shl0_940, shl0_941, shl0_942, \
                         shk_752, shk_753, shk_754, shl1_940, shl1_941, \
                         shl1_942 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1210[k] = pb_y[k] * shl0_940[k]
                    + f_18 * shk_752[k]
                    - f_14 * pc_y[k] * shl1_940[k];

        t_1211[k] = pb_y[k] * shl0_941[k]
                    + f_17 * shk_753[k]
                    - f_14 * pc_y[k] * shl1_941[k];

        t_1212[k] = pb_y[k] * shl0_942[k]
                    + f_16 * shk_754[k]
                    - f_14 * pc_y[k] * shl1_942[k];
    }

#pragma omp simd aligned(t_1213, t_1214, t_1215, t_1216, pb_y, pc_x, pc_y, shl0_944, shk_755, \
                         shl1_944, sii0_756, sii1_756, sik_971, \
                         sik_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1213[k] = f_15 * shk_755[k]
                    + f_3 * pc_y[k] * sik_971[k];

        t_1214[k] = pb_y[k] * shl0_944[k]
                    - f_14 * pc_y[k] * shl1_944[k];

        t_1215[k] = f_1 * sii0_756[k]
                    - f_2 * sii1_756[k]
                    + f_3 * pc_x[k] * sik_972[k];

        t_1216[k] = f_3 * pc_y[k] * sik_972[k];
    }

#pragma omp simd aligned(t_1217, t_1218, t_1219, t_1220, pc_x, pc_y, pc_z, shk_720, sii0_759, \
                         sii0_761, sii1_759, sii1_761, sik_972, sik_974, sik_975, \
                         sik_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1217[k] = f_0 * shk_720[k]
                    + f_3 * pc_z[k] * sik_972[k];

        t_1218[k] = f_4 * sii0_759[k]
                    - f_5 * sii1_759[k]
                    + f_3 * pc_x[k] * sik_975[k];

        t_1219[k] = f_3 * pc_y[k] * sik_974[k];

        t_1220[k] = f_4 * sii0_761[k]
                    - f_5 * sii1_761[k]
                    + f_3 * pc_x[k] * sik_977[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, t_1224, pc_x, pc_y, pc_z, shk_723, sii0_762, \
                         sii0_765, sii1_762, sii1_765, sik_975, sik_977, sik_978, \
                         sik_981 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_6 * sii0_762[k]
                    - f_7 * sii1_762[k]
                    + f_3 * pc_x[k] * sik_978[k];

        t_1222[k] = f_0 * shk_723[k]
                    + f_3 * pc_z[k] * sik_975[k];

        t_1223[k] = f_3 * pc_y[k] * sik_977[k];

        t_1224[k] = f_6 * sii0_765[k]
                    - f_7 * sii1_765[k]
                    + f_3 * pc_x[k] * sik_981[k];
    }

#pragma omp simd aligned(t_1225, t_1226, t_1227, t_1228, pc_x, pc_y, pc_z, shk_726, sii0_766, \
                         sii0_768, sii1_766, sii1_768, sik_978, sik_981, sik_982, \
                         sik_984 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1225[k] = f_8 * sii0_766[k]
                    - f_9 * sii1_766[k]
                    + f_3 * pc_x[k] * sik_982[k];

        t_1226[k] = f_0 * shk_726[k]
                    + f_3 * pc_z[k] * sik_978[k];

        t_1227[k] = f_8 * sii0_768[k]
                    - f_9 * sii1_768[k]
                    + f_3 * pc_x[k] * sik_984[k];

        t_1228[k] = f_3 * pc_y[k] * sik_981[k];
    }

#pragma omp simd aligned(t_1229, t_1230, t_1231, pc_x, pc_z, shk_730, sii0_770, sii0_771, \
                         sii1_770, sii1_771, sik_982, sik_986, \
                         sik_987 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1229[k] = f_8 * sii0_770[k]
                    - f_9 * sii1_770[k]
                    + f_3 * pc_x[k] * sik_986[k];

        t_1230[k] = f_10 * sii0_771[k]
                    - f_11 * sii1_771[k]
                    + f_3 * pc_x[k] * sik_987[k];

        t_1231[k] = f_0 * shk_730[k]
                    + f_3 * pc_z[k] * sik_982[k];
    }

#pragma omp simd aligned(t_1232, t_1233, t_1234, t_1235, pc_x, pc_y, sii0_773, sii0_774, \
                         sii0_776, sii1_773, sii1_774, sii1_776, sik_986, sik_989, sik_990, \
                         sik_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1232[k] = f_10 * sii0_773[k]
                    - f_11 * sii1_773[k]
                    + f_3 * pc_x[k] * sik_989[k];

        t_1233[k] = f_10 * sii0_774[k]
                    - f_11 * sii1_774[k]
                    + f_3 * pc_x[k] * sik_990[k];

        t_1234[k] = f_3 * pc_y[k] * sik_986[k];

        t_1235[k] = f_10 * sii0_776[k]
                    - f_11 * sii1_776[k]
                    + f_3 * pc_x[k] * sik_992[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_x, pc_z, shk_735, sii0_777, sii0_779, \
                         sii1_777, sii1_779, sik_987, sik_993, \
                         sik_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_12 * sii0_777[k]
                    - f_13 * sii1_777[k]
                    + f_3 * pc_x[k] * sik_993[k];

        t_1237[k] = f_0 * shk_735[k]
                    + f_3 * pc_z[k] * sik_987[k];

        t_1238[k] = f_12 * sii0_779[k]
                    - f_13 * sii1_779[k]
                    + f_3 * pc_x[k] * sik_995[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, t_1242, pc_x, pc_y, sii0_780, sii0_781, \
                         sii0_783, sii1_780, sii1_781, sii1_783, sik_992, sik_996, sik_997, \
                         sik_999 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_12 * sii0_780[k]
                    - f_13 * sii1_780[k]
                    + f_3 * pc_x[k] * sik_996[k];

        t_1240[k] = f_12 * sii0_781[k]
                    - f_13 * sii1_781[k]
                    + f_3 * pc_x[k] * sik_997[k];

        t_1241[k] = f_3 * pc_y[k] * sik_992[k];

        t_1242[k] = f_12 * sii0_783[k]
                    - f_13 * sii1_783[k]
                    + f_3 * pc_x[k] * sik_999[k];
    }

#pragma omp simd aligned(t_1243, t_1244, t_1245, t_1246, t_1247, t_1248, t_1249, pc_x, \
                         sik_1000, sik_1001, sik_1002, sik_1003, sik_1004, sik_1005, \
                         sik_1006 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1243[k] = f_3 * pc_x[k] * sik_1000[k];

        t_1244[k] = f_3 * pc_x[k] * sik_1001[k];

        t_1245[k] = f_3 * pc_x[k] * sik_1002[k];

        t_1246[k] = f_3 * pc_x[k] * sik_1003[k];

        t_1247[k] = f_3 * pc_x[k] * sik_1004[k];

        t_1248[k] = f_3 * pc_x[k] * sik_1005[k];

        t_1249[k] = f_3 * pc_x[k] * sik_1006[k];
    }

#pragma omp simd aligned(t_1250, t_1251, t_1252, t_1253, pc_x, pc_y, pc_z, shk_748, sii0_777, \
                         sii0_779, sii1_777, sii1_779, sik_1000, sik_1002, \
                         sik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1250[k] = f_3 * pc_x[k] * sik_1007[k];

        t_1251[k] = f_1 * sii0_777[k]
                    - f_2 * sii1_777[k]
                    + f_3 * pc_y[k] * sik_1000[k];

        t_1252[k] = f_0 * shk_748[k]
                    + f_3 * pc_z[k] * sik_1000[k];

        t_1253[k] = f_4 * sii0_779[k]
                    - f_5 * sii1_779[k]
                    + f_3 * pc_y[k] * sik_1002[k];
    }

#pragma omp simd aligned(t_1254, t_1255, t_1256, pc_y, sii0_780, sii0_781, sii0_782, sii1_780, \
                         sii1_781, sii1_782, sik_1003, sik_1004, \
                         sik_1005 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1254[k] = f_6 * sii0_780[k]
                    - f_7 * sii1_780[k]
                    + f_3 * pc_y[k] * sik_1003[k];

        t_1255[k] = f_8 * sii0_781[k]
                    - f_9 * sii1_781[k]
                    + f_3 * pc_y[k] * sik_1004[k];

        t_1256[k] = f_10 * sii0_782[k]
                    - f_11 * sii1_782[k]
                    + f_3 * pc_y[k] * sik_1005[k];
    }

#pragma omp simd aligned(t_1257, t_1258, t_1259, pc_y, pc_z, shk_755, sii0_783, sii1_783, \
                         sik_1006, sik_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1257[k] = f_12 * sii0_783[k]
                    - f_13 * sii1_783[k]
                    + f_3 * pc_y[k] * sik_1006[k];

        t_1258[k] = f_3 * pc_y[k] * sik_1007[k];

        t_1259[k] = f_0 * shk_755[k]
                    + f_1 * sii0_783[k]
                    - f_2 * sii1_783[k]
                    + f_3 * pc_z[k] * sik_1007[k];
    }
}

auto
compute_prim_sil_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shl0, const size_t shk,
                                                   const size_t shl1, const size_t sii0,
                                                   const size_t sii1, const size_t sik,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sil_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sik, ncols, gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, shl0, shk,
                                                              shl1, sii0, sii1, sik, ncols,
                                                              gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece9(buffer, target, pc, shk, sii0,
                                                              sii1, sik, ncols, gamma, p, q);

    compute_prim_sil_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, shl0,
                                                               shk, shl1, sii0, sii1, sik,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
