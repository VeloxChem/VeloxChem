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


#include "SimdThreeCenterElectronRepulsionVrrRecSLL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_21 = 3.5 / q;

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

    const auto *skl0_0 = buffer.data(skl0 + 0);
    const auto *skl0_3 = buffer.data(skl0 + 3);
    const auto *skl0_5 = buffer.data(skl0 + 5);
    const auto *skl0_6 = buffer.data(skl0 + 6);
    const auto *skl0_9 = buffer.data(skl0 + 9);
    const auto *skl0_10 = buffer.data(skl0 + 10);
    const auto *skl0_12 = buffer.data(skl0 + 12);
    const auto *skl0_14 = buffer.data(skl0 + 14);
    const auto *skl0_15 = buffer.data(skl0 + 15);
    const auto *skl0_17 = buffer.data(skl0 + 17);
    const auto *skl0_18 = buffer.data(skl0 + 18);
    const auto *skl0_20 = buffer.data(skl0 + 20);
    const auto *skl0_21 = buffer.data(skl0 + 21);
    const auto *skl0_23 = buffer.data(skl0 + 23);
    const auto *skl0_24 = buffer.data(skl0 + 24);
    const auto *skl0_25 = buffer.data(skl0 + 25);
    const auto *skl0_27 = buffer.data(skl0 + 27);
    const auto *skl0_44 = buffer.data(skl0 + 44);

    const auto *skk_0 = buffer.data(skk + 0);
    const auto *skk_1 = buffer.data(skk + 1);
    const auto *skk_2 = buffer.data(skk + 2);
    const auto *skk_3 = buffer.data(skk + 3);
    const auto *skk_5 = buffer.data(skk + 5);
    const auto *skk_6 = buffer.data(skk + 6);
    const auto *skk_7 = buffer.data(skk + 7);
    const auto *skk_8 = buffer.data(skk + 8);
    const auto *skk_9 = buffer.data(skk + 9);
    const auto *skk_10 = buffer.data(skk + 10);
    const auto *skk_11 = buffer.data(skk + 11);
    const auto *skk_12 = buffer.data(skk + 12);
    const auto *skk_13 = buffer.data(skk + 13);
    const auto *skk_14 = buffer.data(skk + 14);
    const auto *skk_15 = buffer.data(skk + 15);
    const auto *skk_16 = buffer.data(skk + 16);
    const auto *skk_17 = buffer.data(skk + 17);
    const auto *skk_18 = buffer.data(skk + 18);
    const auto *skk_19 = buffer.data(skk + 19);
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
    const auto *skk_64 = buffer.data(skk + 64);
    const auto *skk_65 = buffer.data(skk + 65);
    const auto *skk_66 = buffer.data(skk + 66);
    const auto *skk_67 = buffer.data(skk + 67);
    const auto *skk_68 = buffer.data(skk + 68);
    const auto *skk_69 = buffer.data(skk + 69);
    const auto *skk_70 = buffer.data(skk + 70);
    const auto *skk_71 = buffer.data(skk + 71);

    const auto *skl1_0 = buffer.data(skl1 + 0);
    const auto *skl1_3 = buffer.data(skl1 + 3);
    const auto *skl1_5 = buffer.data(skl1 + 5);
    const auto *skl1_6 = buffer.data(skl1 + 6);
    const auto *skl1_9 = buffer.data(skl1 + 9);
    const auto *skl1_10 = buffer.data(skl1 + 10);
    const auto *skl1_12 = buffer.data(skl1 + 12);
    const auto *skl1_14 = buffer.data(skl1 + 14);
    const auto *skl1_15 = buffer.data(skl1 + 15);
    const auto *skl1_17 = buffer.data(skl1 + 17);
    const auto *skl1_18 = buffer.data(skl1 + 18);
    const auto *skl1_20 = buffer.data(skl1 + 20);
    const auto *skl1_21 = buffer.data(skl1 + 21);
    const auto *skl1_23 = buffer.data(skl1 + 23);
    const auto *skl1_24 = buffer.data(skl1 + 24);
    const auto *skl1_25 = buffer.data(skl1 + 25);
    const auto *skl1_27 = buffer.data(skl1 + 27);
    const auto *skl1_44 = buffer.data(skl1 + 44);

    const auto *sli0_0 = buffer.data(sli0 + 0);
    const auto *sli0_3 = buffer.data(sli0 + 3);
    const auto *sli0_5 = buffer.data(sli0 + 5);
    const auto *sli0_6 = buffer.data(sli0 + 6);
    const auto *sli0_9 = buffer.data(sli0 + 9);
    const auto *sli0_10 = buffer.data(sli0 + 10);
    const auto *sli0_12 = buffer.data(sli0 + 12);
    const auto *sli0_14 = buffer.data(sli0 + 14);
    const auto *sli0_15 = buffer.data(sli0 + 15);
    const auto *sli0_17 = buffer.data(sli0 + 17);
    const auto *sli0_18 = buffer.data(sli0 + 18);
    const auto *sli0_20 = buffer.data(sli0 + 20);
    const auto *sli0_21 = buffer.data(sli0 + 21);
    const auto *sli0_23 = buffer.data(sli0 + 23);
    const auto *sli0_24 = buffer.data(sli0 + 24);
    const auto *sli0_25 = buffer.data(sli0 + 25);
    const auto *sli0_26 = buffer.data(sli0 + 26);
    const auto *sli0_27 = buffer.data(sli0 + 27);
    const auto *sli0_49 = buffer.data(sli0 + 49);
    const auto *sli0_51 = buffer.data(sli0 + 51);
    const auto *sli0_52 = buffer.data(sli0 + 52);
    const auto *sli0_53 = buffer.data(sli0 + 53);
    const auto *sli0_54 = buffer.data(sli0 + 54);
    const auto *sli0_55 = buffer.data(sli0 + 55);

    const auto *sli1_0 = buffer.data(sli1 + 0);
    const auto *sli1_3 = buffer.data(sli1 + 3);
    const auto *sli1_5 = buffer.data(sli1 + 5);
    const auto *sli1_6 = buffer.data(sli1 + 6);
    const auto *sli1_9 = buffer.data(sli1 + 9);
    const auto *sli1_10 = buffer.data(sli1 + 10);
    const auto *sli1_12 = buffer.data(sli1 + 12);
    const auto *sli1_14 = buffer.data(sli1 + 14);
    const auto *sli1_15 = buffer.data(sli1 + 15);
    const auto *sli1_17 = buffer.data(sli1 + 17);
    const auto *sli1_18 = buffer.data(sli1 + 18);
    const auto *sli1_20 = buffer.data(sli1 + 20);
    const auto *sli1_21 = buffer.data(sli1 + 21);
    const auto *sli1_23 = buffer.data(sli1 + 23);
    const auto *sli1_24 = buffer.data(sli1 + 24);
    const auto *sli1_25 = buffer.data(sli1 + 25);
    const auto *sli1_26 = buffer.data(sli1 + 26);
    const auto *sli1_27 = buffer.data(sli1 + 27);
    const auto *sli1_49 = buffer.data(sli1 + 49);
    const auto *sli1_51 = buffer.data(sli1 + 51);
    const auto *sli1_52 = buffer.data(sli1 + 52);
    const auto *sli1_53 = buffer.data(sli1 + 53);
    const auto *sli1_54 = buffer.data(sli1 + 54);
    const auto *sli1_55 = buffer.data(sli1 + 55);

    const auto *slk_0 = buffer.data(slk + 0);
    const auto *slk_2 = buffer.data(slk + 2);
    const auto *slk_3 = buffer.data(slk + 3);
    const auto *slk_5 = buffer.data(slk + 5);
    const auto *slk_6 = buffer.data(slk + 6);
    const auto *slk_9 = buffer.data(slk + 9);
    const auto *slk_10 = buffer.data(slk + 10);
    const auto *slk_12 = buffer.data(slk + 12);
    const auto *slk_14 = buffer.data(slk + 14);
    const auto *slk_15 = buffer.data(slk + 15);
    const auto *slk_17 = buffer.data(slk + 17);
    const auto *slk_18 = buffer.data(slk + 18);
    const auto *slk_20 = buffer.data(slk + 20);
    const auto *slk_21 = buffer.data(slk + 21);
    const auto *slk_23 = buffer.data(slk + 23);
    const auto *slk_24 = buffer.data(slk + 24);
    const auto *slk_25 = buffer.data(slk + 25);
    const auto *slk_27 = buffer.data(slk + 27);
    const auto *slk_28 = buffer.data(slk + 28);
    const auto *slk_29 = buffer.data(slk + 29);
    const auto *slk_30 = buffer.data(slk + 30);
    const auto *slk_31 = buffer.data(slk + 31);
    const auto *slk_32 = buffer.data(slk + 32);
    const auto *slk_33 = buffer.data(slk + 33);
    const auto *slk_34 = buffer.data(slk + 34);
    const auto *slk_35 = buffer.data(slk + 35);
    const auto *slk_36 = buffer.data(slk + 36);
    const auto *slk_38 = buffer.data(slk + 38);
    const auto *slk_39 = buffer.data(slk + 39);
    const auto *slk_41 = buffer.data(slk + 41);
    const auto *slk_42 = buffer.data(slk + 42);
    const auto *slk_45 = buffer.data(slk + 45);
    const auto *slk_46 = buffer.data(slk + 46);
    const auto *slk_50 = buffer.data(slk + 50);
    const auto *slk_51 = buffer.data(slk + 51);
    const auto *slk_56 = buffer.data(slk + 56);
    const auto *slk_64 = buffer.data(slk + 64);
    const auto *slk_65 = buffer.data(slk + 65);
    const auto *slk_66 = buffer.data(slk + 66);
    const auto *slk_67 = buffer.data(slk + 67);
    const auto *slk_68 = buffer.data(slk + 68);
    const auto *slk_69 = buffer.data(slk + 69);
    const auto *slk_70 = buffer.data(slk + 70);
    const auto *slk_71 = buffer.data(slk + 71);
    const auto *slk_72 = buffer.data(slk + 72);
    const auto *slk_74 = buffer.data(slk + 74);
    const auto *slk_75 = buffer.data(slk + 75);
    const auto *slk_77 = buffer.data(slk + 77);
    const auto *slk_78 = buffer.data(slk + 78);
    const auto *slk_81 = buffer.data(slk + 81);
    const auto *slk_82 = buffer.data(slk + 82);
    const auto *slk_86 = buffer.data(slk + 86);
    const auto *slk_87 = buffer.data(slk + 87);
    const auto *slk_92 = buffer.data(slk + 92);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, skk_0, skk_3, sli0_0, sli0_3, \
                         sli1_0, sli1_3, slk_0, slk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * skk_0[k]
                 + f_1 * sli0_0[k]
                 - f_2 * sli1_0[k]
                 + f_3 * pc_x[k] * slk_0[k];

        t_1[k] = f_3 * pc_y[k] * slk_0[k];

        t_2[k] = f_3 * pc_z[k] * slk_0[k];

        t_3[k] = f_0 * skk_3[k]
                 + f_4 * sli0_3[k]
                 - f_5 * sli1_3[k]
                 + f_3 * pc_x[k] * slk_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, skk_5, skk_6, sli0_5, sli0_6, sli1_5, \
                         sli1_6, slk_2, slk_5, slk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * slk_2[k];

        t_5[k] = f_0 * skk_5[k]
                 + f_4 * sli0_5[k]
                 - f_5 * sli1_5[k]
                 + f_3 * pc_x[k] * slk_5[k];

        t_6[k] = f_0 * skk_6[k]
                 + f_6 * sli0_6[k]
                 - f_7 * sli1_6[k]
                 + f_3 * pc_x[k] * slk_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pc_x, pc_y, pc_z, skk_9, sli0_9, sli1_9, slk_3, slk_5, \
                         slk_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * slk_3[k];

        t_8[k] = f_3 * pc_y[k] * slk_5[k];

        t_9[k] = f_0 * skk_9[k]
                 + f_6 * sli0_9[k]
                 - f_7 * sli1_9[k]
                 + f_3 * pc_x[k] * slk_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pc_x, pc_z, skk_10, skk_12, sli0_10, sli0_12, \
                         sli1_10, sli1_12, slk_6, slk_10, slk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * skk_10[k]
                  + f_8 * sli0_10[k]
                  - f_9 * sli1_10[k]
                  + f_3 * pc_x[k] * slk_10[k];

        t_11[k] = f_3 * pc_z[k] * slk_6[k];

        t_12[k] = f_0 * skk_12[k]
                  + f_8 * sli0_12[k]
                  - f_9 * sli1_12[k]
                  + f_3 * pc_x[k] * slk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pc_x, pc_y, skk_14, skk_15, sli0_14, sli0_15, \
                         sli1_14, sli1_15, slk_9, slk_14, slk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_3 * pc_y[k] * slk_9[k];

        t_14[k] = f_0 * skk_14[k]
                  + f_8 * sli0_14[k]
                  - f_9 * sli1_14[k]
                  + f_3 * pc_x[k] * slk_14[k];

        t_15[k] = f_0 * skk_15[k]
                  + f_10 * sli0_15[k]
                  - f_11 * sli1_15[k]
                  + f_3 * pc_x[k] * slk_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_z, skk_17, skk_18, sli0_17, sli0_18, \
                         sli1_17, sli1_18, slk_10, slk_17, slk_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_z[k] * slk_10[k];

        t_17[k] = f_0 * skk_17[k]
                  + f_10 * sli0_17[k]
                  - f_11 * sli1_17[k]
                  + f_3 * pc_x[k] * slk_17[k];

        t_18[k] = f_0 * skk_18[k]
                  + f_10 * sli0_18[k]
                  - f_11 * sli1_18[k]
                  + f_3 * pc_x[k] * slk_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, pc_x, pc_y, skk_20, skk_21, sli0_20, sli0_21, \
                         sli1_20, sli1_21, slk_14, slk_20, slk_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * slk_14[k];

        t_20[k] = f_0 * skk_20[k]
                  + f_10 * sli0_20[k]
                  - f_11 * sli1_20[k]
                  + f_3 * pc_x[k] * slk_20[k];

        t_21[k] = f_0 * skk_21[k]
                  + f_12 * sli0_21[k]
                  - f_13 * sli1_21[k]
                  + f_3 * pc_x[k] * slk_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, pc_x, pc_z, skk_23, skk_24, sli0_23, sli0_24, \
                         sli1_23, sli1_24, slk_15, slk_23, slk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_3 * pc_z[k] * slk_15[k];

        t_23[k] = f_0 * skk_23[k]
                  + f_12 * sli0_23[k]
                  - f_13 * sli1_23[k]
                  + f_3 * pc_x[k] * slk_23[k];

        t_24[k] = f_0 * skk_24[k]
                  + f_12 * sli0_24[k]
                  - f_13 * sli1_24[k]
                  + f_3 * pc_x[k] * slk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, pc_x, pc_y, skk_25, skk_27, sli0_25, sli0_27, \
                         sli1_25, sli1_27, slk_20, slk_25, slk_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * skk_25[k]
                  + f_12 * sli0_25[k]
                  - f_13 * sli1_25[k]
                  + f_3 * pc_x[k] * slk_25[k];

        t_26[k] = f_3 * pc_y[k] * slk_20[k];

        t_27[k] = f_0 * skk_27[k]
                  + f_12 * sli0_27[k]
                  - f_13 * sli1_27[k]
                  + f_3 * pc_x[k] * slk_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, pc_x, skk_28, skk_29, skk_30, skk_31, \
                         skk_32, slk_28, slk_29, slk_30, slk_31, \
                         slk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * skk_28[k]
                  + f_3 * pc_x[k] * slk_28[k];

        t_29[k] = f_0 * skk_29[k]
                  + f_3 * pc_x[k] * slk_29[k];

        t_30[k] = f_0 * skk_30[k]
                  + f_3 * pc_x[k] * slk_30[k];

        t_31[k] = f_0 * skk_31[k]
                  + f_3 * pc_x[k] * slk_31[k];

        t_32[k] = f_0 * skk_32[k]
                  + f_3 * pc_x[k] * slk_32[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pc_x, pc_y, skk_33, skk_34, skk_35, sli0_21, \
                         sli1_21, slk_28, slk_33, slk_34, slk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_0 * skk_33[k]
                  + f_3 * pc_x[k] * slk_33[k];

        t_34[k] = f_0 * skk_34[k]
                  + f_3 * pc_x[k] * slk_34[k];

        t_35[k] = f_0 * skk_35[k]
                  + f_3 * pc_x[k] * slk_35[k];

        t_36[k] = f_1 * sli0_21[k]
                  - f_2 * sli1_21[k]
                  + f_3 * pc_y[k] * slk_28[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pc_y, pc_z, sli0_23, sli0_24, sli0_25, \
                         sli1_23, sli1_24, sli1_25, slk_28, slk_30, slk_31, \
                         slk_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * pc_z[k] * slk_28[k];

        t_38[k] = f_4 * sli0_23[k]
                  - f_5 * sli1_23[k]
                  + f_3 * pc_y[k] * slk_30[k];

        t_39[k] = f_6 * sli0_24[k]
                  - f_7 * sli1_24[k]
                  + f_3 * pc_y[k] * slk_31[k];

        t_40[k] = f_8 * sli0_25[k]
                  - f_9 * sli1_25[k]
                  + f_3 * pc_y[k] * slk_32[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sli0_26, sli0_27, sli1_26, \
                         sli1_27, slk_33, slk_34, slk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * sli0_26[k]
                  - f_11 * sli1_26[k]
                  + f_3 * pc_y[k] * slk_33[k];

        t_42[k] = f_12 * sli0_27[k]
                  - f_13 * sli1_27[k]
                  + f_3 * pc_y[k] * slk_34[k];

        t_43[k] = f_3 * pc_y[k] * slk_35[k];

        t_44[k] = f_1 * sli0_27[k]
                  - f_2 * sli1_27[k]
                  + f_3 * pc_z[k] * slk_35[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_y, pc_y, pc_z, skl0_0, skl0_3, skk_0, \
                         skk_1, skl1_0, skl1_3, slk_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pb_y[k] * skl0_0[k]
                  - f_14 * pc_y[k] * skl1_0[k];

        t_46[k] = f_15 * skk_0[k]
                  + f_3 * pc_y[k] * slk_36[k];

        t_47[k] = f_3 * pc_z[k] * slk_36[k];

        t_48[k] = pb_y[k] * skl0_3[k]
                  + f_16 * skk_1[k]
                  - f_14 * pc_y[k] * skl1_3[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pb_y, pc_y, pc_z, skl0_5, skl0_6, skk_2, \
                         skk_3, skl1_5, skl1_6, slk_38, slk_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_15 * skk_2[k]
                  + f_3 * pc_y[k] * slk_38[k];

        t_50[k] = pb_y[k] * skl0_5[k]
                  - f_14 * pc_y[k] * skl1_5[k];

        t_51[k] = pb_y[k] * skl0_6[k]
                  + f_17 * skk_3[k]
                  - f_14 * pc_y[k] * skl1_6[k];

        t_52[k] = f_3 * pc_z[k] * slk_39[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pb_y, pc_y, pc_z, skl0_9, skl0_10, skk_5, \
                         skk_6, skl1_9, skl1_10, slk_41, slk_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_15 * skk_5[k]
                  + f_3 * pc_y[k] * slk_41[k];

        t_54[k] = pb_y[k] * skl0_9[k]
                  - f_14 * pc_y[k] * skl1_9[k];

        t_55[k] = pb_y[k] * skl0_10[k]
                  + f_18 * skk_6[k]
                  - f_14 * pc_y[k] * skl1_10[k];

        t_56[k] = f_3 * pc_z[k] * slk_42[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, skl0_12, skl0_14, skl0_15, skk_8, \
                         skk_9, skk_10, skl1_12, skl1_14, skl1_15, \
                         slk_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_y[k] * skl0_12[k]
                  + f_16 * skk_8[k]
                  - f_14 * pc_y[k] * skl1_12[k];

        t_58[k] = f_15 * skk_9[k]
                  + f_3 * pc_y[k] * slk_45[k];

        t_59[k] = pb_y[k] * skl0_14[k]
                  - f_14 * pc_y[k] * skl1_14[k];

        t_60[k] = pb_y[k] * skl0_15[k]
                  + f_19 * skk_10[k]
                  - f_14 * pc_y[k] * skl1_15[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_y, pc_y, pc_z, skl0_17, skl0_18, skk_12, \
                         skk_13, skk_14, skl1_17, skl1_18, slk_46, \
                         slk_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_3 * pc_z[k] * slk_46[k];

        t_62[k] = pb_y[k] * skl0_17[k]
                  + f_17 * skk_12[k]
                  - f_14 * pc_y[k] * skl1_17[k];

        t_63[k] = pb_y[k] * skl0_18[k]
                  + f_16 * skk_13[k]
                  - f_14 * pc_y[k] * skl1_18[k];

        t_64[k] = f_15 * skk_14[k]
                  + f_3 * pc_y[k] * slk_50[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_y, pc_z, skl0_20, skl0_21, skl0_23, \
                         skk_15, skk_17, skl1_20, skl1_21, skl1_23, \
                         slk_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * skl0_20[k]
                  - f_14 * pc_y[k] * skl1_20[k];

        t_66[k] = pb_y[k] * skl0_21[k]
                  + f_20 * skk_15[k]
                  - f_14 * pc_y[k] * skl1_21[k];

        t_67[k] = f_3 * pc_z[k] * slk_51[k];

        t_68[k] = pb_y[k] * skl0_23[k]
                  + f_18 * skk_17[k]
                  - f_14 * pc_y[k] * skl1_23[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pb_y, pc_y, skl0_24, skl0_25, skl0_27, \
                         skk_18, skk_19, skk_20, skl1_24, skl1_25, skl1_27, \
                         slk_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * skl0_24[k]
                  + f_17 * skk_18[k]
                  - f_14 * pc_y[k] * skl1_24[k];

        t_70[k] = pb_y[k] * skl0_25[k]
                  + f_16 * skk_19[k]
                  - f_14 * pc_y[k] * skl1_25[k];

        t_71[k] = f_15 * skk_20[k]
                  + f_3 * pc_y[k] * slk_56[k];

        t_72[k] = pb_y[k] * skl0_27[k]
                  - f_14 * pc_y[k] * skl1_27[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pc_x, skk_64, skk_65, skk_66, skk_67, \
                         skk_68, slk_64, slk_65, slk_66, slk_67, \
                         slk_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_21 * skk_64[k]
                  + f_3 * pc_x[k] * slk_64[k];

        t_74[k] = f_21 * skk_65[k]
                  + f_3 * pc_x[k] * slk_65[k];

        t_75[k] = f_21 * skk_66[k]
                  + f_3 * pc_x[k] * slk_66[k];

        t_76[k] = f_21 * skk_67[k]
                  + f_3 * pc_x[k] * slk_67[k];

        t_77[k] = f_21 * skk_68[k]
                  + f_3 * pc_x[k] * slk_68[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pc_x, pc_y, skk_28, skk_69, skk_70, skk_71, \
                         sli0_49, sli1_49, slk_64, slk_69, slk_70, \
                         slk_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_21 * skk_69[k]
                  + f_3 * pc_x[k] * slk_69[k];

        t_79[k] = f_21 * skk_70[k]
                  + f_3 * pc_x[k] * slk_70[k];

        t_80[k] = f_21 * skk_71[k]
                  + f_3 * pc_x[k] * slk_71[k];

        t_81[k] = f_15 * skk_28[k]
                  + f_1 * sli0_49[k]
                  - f_2 * sli1_49[k]
                  + f_3 * pc_y[k] * slk_64[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pc_y, pc_z, skk_30, skk_31, sli0_51, sli0_52, \
                         sli1_51, sli1_52, slk_64, slk_66, slk_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_3 * pc_z[k] * slk_64[k];

        t_83[k] = f_15 * skk_30[k]
                  + f_4 * sli0_51[k]
                  - f_5 * sli1_51[k]
                  + f_3 * pc_y[k] * slk_66[k];

        t_84[k] = f_15 * skk_31[k]
                  + f_6 * sli0_52[k]
                  - f_7 * sli1_52[k]
                  + f_3 * pc_y[k] * slk_67[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pc_y, skk_32, skk_33, skk_34, sli0_53, sli0_54, \
                         sli0_55, sli1_53, sli1_54, sli1_55, slk_68, slk_69, \
                         slk_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_15 * skk_32[k]
                  + f_8 * sli0_53[k]
                  - f_9 * sli1_53[k]
                  + f_3 * pc_y[k] * slk_68[k];

        t_86[k] = f_15 * skk_33[k]
                  + f_10 * sli0_54[k]
                  - f_11 * sli1_54[k]
                  + f_3 * pc_y[k] * slk_69[k];

        t_87[k] = f_15 * skk_34[k]
                  + f_12 * sli0_55[k]
                  - f_13 * sli1_55[k]
                  + f_3 * pc_y[k] * slk_70[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pb_y, pb_z, pc_y, pc_z, skl0_0, skl0_44, \
                         skk_35, skl1_0, skl1_44, slk_71, slk_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_15 * skk_35[k]
                  + f_3 * pc_y[k] * slk_71[k];

        t_89[k] = pb_y[k] * skl0_44[k]
                  - f_14 * pc_y[k] * skl1_44[k];

        t_90[k] = pb_z[k] * skl0_0[k]
                  - f_14 * pc_z[k] * skl1_0[k];

        t_91[k] = f_3 * pc_y[k] * slk_72[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pb_z, pc_y, pc_z, skl0_3, skl0_5, skk_0, \
                         skk_2, skl1_3, skl1_5, slk_72, slk_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_15 * skk_0[k]
                  + f_3 * pc_z[k] * slk_72[k];

        t_93[k] = pb_z[k] * skl0_3[k]
                  - f_14 * pc_z[k] * skl1_3[k];

        t_94[k] = f_3 * pc_y[k] * slk_74[k];

        t_95[k] = pb_z[k] * skl0_5[k]
                  + f_16 * skk_2[k]
                  - f_14 * pc_z[k] * skl1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pb_z, pc_y, pc_z, skl0_6, skl0_9, skk_3, \
                         skk_5, skl1_6, skl1_9, slk_75, slk_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_z[k] * skl0_6[k]
                  - f_14 * pc_z[k] * skl1_6[k];

        t_97[k] = f_15 * skk_3[k]
                  + f_3 * pc_z[k] * slk_75[k];

        t_98[k] = f_3 * pc_y[k] * slk_77[k];

        t_99[k] = pb_z[k] * skl0_9[k]
                  + f_17 * skk_5[k]
                  - f_14 * pc_z[k] * skl1_9[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pb_z, pc_y, pc_z, skl0_10, skl0_12, \
                         skk_6, skk_7, skl1_10, skl1_12, slk_78, \
                         slk_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pb_z[k] * skl0_10[k]
                   - f_14 * pc_z[k] * skl1_10[k];

        t_101[k] = f_15 * skk_6[k]
                   + f_3 * pc_z[k] * slk_78[k];

        t_102[k] = pb_z[k] * skl0_12[k]
                   + f_16 * skk_7[k]
                   - f_14 * pc_z[k] * skl1_12[k];

        t_103[k] = f_3 * pc_y[k] * slk_81[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_z, pc_z, skl0_14, skl0_15, skl0_17, \
                         skk_9, skk_10, skk_11, skl1_14, skl1_15, skl1_17, \
                         slk_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_z[k] * skl0_14[k]
                   + f_18 * skk_9[k]
                   - f_14 * pc_z[k] * skl1_14[k];

        t_105[k] = pb_z[k] * skl0_15[k]
                   - f_14 * pc_z[k] * skl1_15[k];

        t_106[k] = f_15 * skk_10[k]
                   + f_3 * pc_z[k] * slk_82[k];

        t_107[k] = pb_z[k] * skl0_17[k]
                   + f_16 * skk_11[k]
                   - f_14 * pc_z[k] * skl1_17[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pb_z, pc_y, pc_z, skl0_18, skl0_20, \
                         skl0_21, skk_12, skk_14, skl1_18, skl1_20, skl1_21, \
                         slk_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * skl0_18[k]
                   + f_17 * skk_12[k]
                   - f_14 * pc_z[k] * skl1_18[k];

        t_109[k] = f_3 * pc_y[k] * slk_86[k];

        t_110[k] = pb_z[k] * skl0_20[k]
                   + f_19 * skk_14[k]
                   - f_14 * pc_z[k] * skl1_20[k];

        t_111[k] = pb_z[k] * skl0_21[k]
                   - f_14 * pc_z[k] * skl1_21[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_z, pc_z, skl0_23, skl0_24, skk_15, skk_16, \
                         skk_17, skl1_23, skl1_24, slk_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_15 * skk_15[k]
                   + f_3 * pc_z[k] * slk_87[k];

        t_113[k] = pb_z[k] * skl0_23[k]
                   + f_16 * skk_16[k]
                   - f_14 * pc_z[k] * skl1_23[k];

        t_114[k] = pb_z[k] * skl0_24[k]
                   + f_17 * skk_17[k]
                   - f_14 * pc_z[k] * skl1_24[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_z, pc_y, pc_z, skl0_25, skl0_27, skk_18, \
                         skk_20, skl1_25, skl1_27, slk_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pb_z[k] * skl0_25[k]
                   + f_18 * skk_18[k]
                   - f_14 * pc_z[k] * skl1_25[k];

        t_116[k] = f_3 * pc_y[k] * slk_92[k];

        t_117[k] = pb_z[k] * skl0_27[k]
                   + f_20 * skk_20[k]
                   - f_14 * pc_z[k] * skl1_27[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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
    const auto f_20 = 3.0 / q;
    const auto f_21 = 3.5 / q;

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

    const auto *skl0_36 = buffer.data(skl0 + 36);
    const auto *skl0_48 = buffer.data(skl0 + 48);
    const auto *skl0_51 = buffer.data(skl0 + 51);
    const auto *skl0_55 = buffer.data(skl0 + 55);
    const auto *skl0_60 = buffer.data(skl0 + 60);
    const auto *skl0_66 = buffer.data(skl0 + 66);
    const auto *skl0_81 = buffer.data(skl0 + 81);
    const auto *skl0_90 = buffer.data(skl0 + 90);
    const auto *skl0_95 = buffer.data(skl0 + 95);
    const auto *skl0_99 = buffer.data(skl0 + 99);
    const auto *skl0_102 = buffer.data(skl0 + 102);
    const auto *skl0_104 = buffer.data(skl0 + 104);
    const auto *skl0_107 = buffer.data(skl0 + 107);
    const auto *skl0_108 = buffer.data(skl0 + 108);
    const auto *skl0_110 = buffer.data(skl0 + 110);
    const auto *skl0_113 = buffer.data(skl0 + 113);
    const auto *skl0_114 = buffer.data(skl0 + 114);
    const auto *skl0_115 = buffer.data(skl0 + 115);
    const auto *skl0_117 = buffer.data(skl0 + 117);
    const auto *skl0_134 = buffer.data(skl0 + 134);

    const auto *skk_28 = buffer.data(skk + 28);
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
    const auto *skk_80 = buffer.data(skk + 80);
    const auto *skk_81 = buffer.data(skk + 81);
    const auto *skk_84 = buffer.data(skk + 84);
    const auto *skk_85 = buffer.data(skk + 85);
    const auto *skk_86 = buffer.data(skk + 86);
    const auto *skk_89 = buffer.data(skk + 89);
    const auto *skk_90 = buffer.data(skk + 90);
    const auto *skk_91 = buffer.data(skk + 91);
    const auto *skk_92 = buffer.data(skk + 92);
    const auto *skk_100 = buffer.data(skk + 100);
    const auto *skk_101 = buffer.data(skk + 101);
    const auto *skk_102 = buffer.data(skk + 102);
    const auto *skk_103 = buffer.data(skk + 103);
    const auto *skk_104 = buffer.data(skk + 104);
    const auto *skk_105 = buffer.data(skk + 105);
    const auto *skk_106 = buffer.data(skk + 106);
    const auto *skk_107 = buffer.data(skk + 107);
    const auto *skk_108 = buffer.data(skk + 108);
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
    const auto *skk_172 = buffer.data(skk + 172);
    const auto *skk_173 = buffer.data(skk + 173);
    const auto *skk_174 = buffer.data(skk + 174);
    const auto *skk_175 = buffer.data(skk + 175);
    const auto *skk_176 = buffer.data(skk + 176);
    const auto *skk_177 = buffer.data(skk + 177);
    const auto *skk_178 = buffer.data(skk + 178);
    const auto *skk_179 = buffer.data(skk + 179);
    const auto *skk_180 = buffer.data(skk + 180);
    const auto *skk_183 = buffer.data(skk + 183);
    const auto *skk_185 = buffer.data(skk + 185);
    const auto *skk_186 = buffer.data(skk + 186);

    const auto *skl1_36 = buffer.data(skl1 + 36);
    const auto *skl1_48 = buffer.data(skl1 + 48);
    const auto *skl1_51 = buffer.data(skl1 + 51);
    const auto *skl1_55 = buffer.data(skl1 + 55);
    const auto *skl1_60 = buffer.data(skl1 + 60);
    const auto *skl1_66 = buffer.data(skl1 + 66);
    const auto *skl1_81 = buffer.data(skl1 + 81);
    const auto *skl1_90 = buffer.data(skl1 + 90);
    const auto *skl1_95 = buffer.data(skl1 + 95);
    const auto *skl1_99 = buffer.data(skl1 + 99);
    const auto *skl1_102 = buffer.data(skl1 + 102);
    const auto *skl1_104 = buffer.data(skl1 + 104);
    const auto *skl1_107 = buffer.data(skl1 + 107);
    const auto *skl1_108 = buffer.data(skl1 + 108);
    const auto *skl1_110 = buffer.data(skl1 + 110);
    const auto *skl1_113 = buffer.data(skl1 + 113);
    const auto *skl1_114 = buffer.data(skl1 + 114);
    const auto *skl1_115 = buffer.data(skl1 + 115);
    const auto *skl1_117 = buffer.data(skl1 + 117);
    const auto *skl1_134 = buffer.data(skl1 + 134);

    const auto *sli0_79 = buffer.data(sli0 + 79);
    const auto *sli0_80 = buffer.data(sli0 + 80);
    const auto *sli0_81 = buffer.data(sli0 + 81);
    const auto *sli0_82 = buffer.data(sli0 + 82);
    const auto *sli0_83 = buffer.data(sli0 + 83);
    const auto *sli0_84 = buffer.data(sli0 + 84);
    const auto *sli0_87 = buffer.data(sli0 + 87);
    const auto *sli0_89 = buffer.data(sli0 + 89);
    const auto *sli0_90 = buffer.data(sli0 + 90);
    const auto *sli0_93 = buffer.data(sli0 + 93);
    const auto *sli0_94 = buffer.data(sli0 + 94);
    const auto *sli0_96 = buffer.data(sli0 + 96);
    const auto *sli0_98 = buffer.data(sli0 + 98);
    const auto *sli0_99 = buffer.data(sli0 + 99);
    const auto *sli0_101 = buffer.data(sli0 + 101);
    const auto *sli0_102 = buffer.data(sli0 + 102);
    const auto *sli0_104 = buffer.data(sli0 + 104);
    const auto *sli0_105 = buffer.data(sli0 + 105);
    const auto *sli0_107 = buffer.data(sli0 + 107);
    const auto *sli0_108 = buffer.data(sli0 + 108);
    const auto *sli0_109 = buffer.data(sli0 + 109);
    const auto *sli0_110 = buffer.data(sli0 + 110);
    const auto *sli0_111 = buffer.data(sli0 + 111);
    const auto *sli0_135 = buffer.data(sli0 + 135);
    const auto *sli0_136 = buffer.data(sli0 + 136);
    const auto *sli0_137 = buffer.data(sli0 + 137);
    const auto *sli0_138 = buffer.data(sli0 + 138);
    const auto *sli0_139 = buffer.data(sli0 + 139);
    const auto *sli0_140 = buffer.data(sli0 + 140);
    const auto *sli0_143 = buffer.data(sli0 + 143);
    const auto *sli0_145 = buffer.data(sli0 + 145);
    const auto *sli0_146 = buffer.data(sli0 + 146);

    const auto *sli1_79 = buffer.data(sli1 + 79);
    const auto *sli1_80 = buffer.data(sli1 + 80);
    const auto *sli1_81 = buffer.data(sli1 + 81);
    const auto *sli1_82 = buffer.data(sli1 + 82);
    const auto *sli1_83 = buffer.data(sli1 + 83);
    const auto *sli1_84 = buffer.data(sli1 + 84);
    const auto *sli1_87 = buffer.data(sli1 + 87);
    const auto *sli1_89 = buffer.data(sli1 + 89);
    const auto *sli1_90 = buffer.data(sli1 + 90);
    const auto *sli1_93 = buffer.data(sli1 + 93);
    const auto *sli1_94 = buffer.data(sli1 + 94);
    const auto *sli1_96 = buffer.data(sli1 + 96);
    const auto *sli1_98 = buffer.data(sli1 + 98);
    const auto *sli1_99 = buffer.data(sli1 + 99);
    const auto *sli1_101 = buffer.data(sli1 + 101);
    const auto *sli1_102 = buffer.data(sli1 + 102);
    const auto *sli1_104 = buffer.data(sli1 + 104);
    const auto *sli1_105 = buffer.data(sli1 + 105);
    const auto *sli1_107 = buffer.data(sli1 + 107);
    const auto *sli1_108 = buffer.data(sli1 + 108);
    const auto *sli1_109 = buffer.data(sli1 + 109);
    const auto *sli1_110 = buffer.data(sli1 + 110);
    const auto *sli1_111 = buffer.data(sli1 + 111);
    const auto *sli1_135 = buffer.data(sli1 + 135);
    const auto *sli1_136 = buffer.data(sli1 + 136);
    const auto *sli1_137 = buffer.data(sli1 + 137);
    const auto *sli1_138 = buffer.data(sli1 + 138);
    const auto *sli1_139 = buffer.data(sli1 + 139);
    const auto *sli1_140 = buffer.data(sli1 + 140);
    const auto *sli1_143 = buffer.data(sli1 + 143);
    const auto *sli1_145 = buffer.data(sli1 + 145);
    const auto *sli1_146 = buffer.data(sli1 + 146);

    const auto *slk_100 = buffer.data(slk + 100);
    const auto *slk_101 = buffer.data(slk + 101);
    const auto *slk_102 = buffer.data(slk + 102);
    const auto *slk_103 = buffer.data(slk + 103);
    const auto *slk_104 = buffer.data(slk + 104);
    const auto *slk_105 = buffer.data(slk + 105);
    const auto *slk_106 = buffer.data(slk + 106);
    const auto *slk_107 = buffer.data(slk + 107);
    const auto *slk_108 = buffer.data(slk + 108);
    const auto *slk_110 = buffer.data(slk + 110);
    const auto *slk_111 = buffer.data(slk + 111);
    const auto *slk_113 = buffer.data(slk + 113);
    const auto *slk_114 = buffer.data(slk + 114);
    const auto *slk_117 = buffer.data(slk + 117);
    const auto *slk_118 = buffer.data(slk + 118);
    const auto *slk_120 = buffer.data(slk + 120);
    const auto *slk_122 = buffer.data(slk + 122);
    const auto *slk_123 = buffer.data(slk + 123);
    const auto *slk_125 = buffer.data(slk + 125);
    const auto *slk_126 = buffer.data(slk + 126);
    const auto *slk_128 = buffer.data(slk + 128);
    const auto *slk_129 = buffer.data(slk + 129);
    const auto *slk_131 = buffer.data(slk + 131);
    const auto *slk_132 = buffer.data(slk + 132);
    const auto *slk_133 = buffer.data(slk + 133);
    const auto *slk_135 = buffer.data(slk + 135);
    const auto *slk_136 = buffer.data(slk + 136);
    const auto *slk_137 = buffer.data(slk + 137);
    const auto *slk_138 = buffer.data(slk + 138);
    const auto *slk_139 = buffer.data(slk + 139);
    const auto *slk_140 = buffer.data(slk + 140);
    const auto *slk_141 = buffer.data(slk + 141);
    const auto *slk_142 = buffer.data(slk + 142);
    const auto *slk_143 = buffer.data(slk + 143);
    const auto *slk_144 = buffer.data(slk + 144);
    const auto *slk_146 = buffer.data(slk + 146);
    const auto *slk_147 = buffer.data(slk + 147);
    const auto *slk_149 = buffer.data(slk + 149);
    const auto *slk_150 = buffer.data(slk + 150);
    const auto *slk_153 = buffer.data(slk + 153);
    const auto *slk_154 = buffer.data(slk + 154);
    const auto *slk_158 = buffer.data(slk + 158);
    const auto *slk_159 = buffer.data(slk + 159);
    const auto *slk_164 = buffer.data(slk + 164);
    const auto *slk_172 = buffer.data(slk + 172);
    const auto *slk_173 = buffer.data(slk + 173);
    const auto *slk_174 = buffer.data(slk + 174);
    const auto *slk_175 = buffer.data(slk + 175);
    const auto *slk_176 = buffer.data(slk + 176);
    const auto *slk_177 = buffer.data(slk + 177);
    const auto *slk_178 = buffer.data(slk + 178);
    const auto *slk_179 = buffer.data(slk + 179);
    const auto *slk_180 = buffer.data(slk + 180);
    const auto *slk_182 = buffer.data(slk + 182);
    const auto *slk_183 = buffer.data(slk + 183);
    const auto *slk_185 = buffer.data(slk + 185);
    const auto *slk_186 = buffer.data(slk + 186);

#pragma omp simd aligned(t_118, t_119, t_120, t_121, t_122, pc_x, skk_100, skk_101, skk_102, \
                         skk_103, skk_104, slk_100, slk_101, slk_102, slk_103, \
                         slk_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_21 * skk_100[k]
                   + f_3 * pc_x[k] * slk_100[k];

        t_119[k] = f_21 * skk_101[k]
                   + f_3 * pc_x[k] * slk_101[k];

        t_120[k] = f_21 * skk_102[k]
                   + f_3 * pc_x[k] * slk_102[k];

        t_121[k] = f_21 * skk_103[k]
                   + f_3 * pc_x[k] * slk_103[k];

        t_122[k] = f_21 * skk_104[k]
                   + f_3 * pc_x[k] * slk_104[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pb_z, pc_x, pc_z, skl0_36, skk_105, \
                         skk_106, skk_107, skl1_36, slk_105, slk_106, \
                         slk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_21 * skk_105[k]
                   + f_3 * pc_x[k] * slk_105[k];

        t_124[k] = f_21 * skk_106[k]
                   + f_3 * pc_x[k] * slk_106[k];

        t_125[k] = f_21 * skk_107[k]
                   + f_3 * pc_x[k] * slk_107[k];

        t_126[k] = pb_z[k] * skl0_36[k]
                   - f_14 * pc_z[k] * skl1_36[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pc_y, pc_z, skk_28, sli0_79, sli0_80, sli1_79, \
                         sli1_80, slk_100, slk_102, slk_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_15 * skk_28[k]
                   + f_3 * pc_z[k] * slk_100[k];

        t_128[k] = f_4 * sli0_79[k]
                   - f_5 * sli1_79[k]
                   + f_3 * pc_y[k] * slk_102[k];

        t_129[k] = f_6 * sli0_80[k]
                   - f_7 * sli1_80[k]
                   + f_3 * pc_y[k] * slk_103[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, pc_y, sli0_81, sli0_82, sli0_83, sli1_81, \
                         sli1_82, sli1_83, slk_104, slk_105, slk_106, \
                         slk_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_8 * sli0_81[k]
                   - f_9 * sli1_81[k]
                   + f_3 * pc_y[k] * slk_104[k];

        t_131[k] = f_10 * sli0_82[k]
                   - f_11 * sli1_82[k]
                   + f_3 * pc_y[k] * slk_105[k];

        t_132[k] = f_12 * sli0_83[k]
                   - f_13 * sli1_83[k]
                   + f_3 * pc_y[k] * slk_106[k];

        t_133[k] = f_3 * pc_y[k] * slk_107[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, t_137, pc_x, pc_y, pc_z, skk_35, skk_36, \
                         skk_108, sli0_83, sli0_84, sli1_83, sli1_84, slk_107, \
                         slk_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_15 * skk_35[k]
                   + f_1 * sli0_83[k]
                   - f_2 * sli1_83[k]
                   + f_3 * pc_z[k] * slk_107[k];

        t_135[k] = f_20 * skk_108[k]
                   + f_1 * sli0_84[k]
                   - f_2 * sli1_84[k]
                   + f_3 * pc_x[k] * slk_108[k];

        t_136[k] = f_16 * skk_36[k]
                   + f_3 * pc_y[k] * slk_108[k];

        t_137[k] = f_3 * pc_z[k] * slk_108[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_y, skk_38, skk_111, skk_113, sli0_87, \
                         sli0_89, sli1_87, sli1_89, slk_110, slk_111, \
                         slk_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_20 * skk_111[k]
                   + f_4 * sli0_87[k]
                   - f_5 * sli1_87[k]
                   + f_3 * pc_x[k] * slk_111[k];

        t_139[k] = f_16 * skk_38[k]
                   + f_3 * pc_y[k] * slk_110[k];

        t_140[k] = f_20 * skk_113[k]
                   + f_4 * sli0_89[k]
                   - f_5 * sli1_89[k]
                   + f_3 * pc_x[k] * slk_113[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_x, pc_y, pc_z, skk_41, skk_114, sli0_90, \
                         sli1_90, slk_111, slk_113, slk_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_20 * skk_114[k]
                   + f_6 * sli0_90[k]
                   - f_7 * sli1_90[k]
                   + f_3 * pc_x[k] * slk_114[k];

        t_142[k] = f_3 * pc_z[k] * slk_111[k];

        t_143[k] = f_16 * skk_41[k]
                   + f_3 * pc_y[k] * slk_113[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_z, skk_117, skk_118, sli0_93, sli0_94, \
                         sli1_93, sli1_94, slk_114, slk_117, slk_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_20 * skk_117[k]
                   + f_6 * sli0_93[k]
                   - f_7 * sli1_93[k]
                   + f_3 * pc_x[k] * slk_117[k];

        t_145[k] = f_20 * skk_118[k]
                   + f_8 * sli0_94[k]
                   - f_9 * sli1_94[k]
                   + f_3 * pc_x[k] * slk_118[k];

        t_146[k] = f_3 * pc_z[k] * slk_114[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, pc_x, pc_y, skk_45, skk_120, skk_122, sli0_96, \
                         sli0_98, sli1_96, sli1_98, slk_117, slk_120, \
                         slk_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_20 * skk_120[k]
                   + f_8 * sli0_96[k]
                   - f_9 * sli1_96[k]
                   + f_3 * pc_x[k] * slk_120[k];

        t_148[k] = f_16 * skk_45[k]
                   + f_3 * pc_y[k] * slk_117[k];

        t_149[k] = f_20 * skk_122[k]
                   + f_8 * sli0_98[k]
                   - f_9 * sli1_98[k]
                   + f_3 * pc_x[k] * slk_122[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_z, skk_123, skk_125, sli0_99, sli0_101, \
                         sli1_99, sli1_101, slk_118, slk_123, slk_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_20 * skk_123[k]
                   + f_10 * sli0_99[k]
                   - f_11 * sli1_99[k]
                   + f_3 * pc_x[k] * slk_123[k];

        t_151[k] = f_3 * pc_z[k] * slk_118[k];

        t_152[k] = f_20 * skk_125[k]
                   + f_10 * sli0_101[k]
                   - f_11 * sli1_101[k]
                   + f_3 * pc_x[k] * slk_125[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, pc_x, pc_y, skk_50, skk_126, skk_128, sli0_102, \
                         sli0_104, sli1_102, sli1_104, slk_122, slk_126, \
                         slk_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_20 * skk_126[k]
                   + f_10 * sli0_102[k]
                   - f_11 * sli1_102[k]
                   + f_3 * pc_x[k] * slk_126[k];

        t_154[k] = f_16 * skk_50[k]
                   + f_3 * pc_y[k] * slk_122[k];

        t_155[k] = f_20 * skk_128[k]
                   + f_10 * sli0_104[k]
                   - f_11 * sli1_104[k]
                   + f_3 * pc_x[k] * slk_128[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pc_x, pc_z, skk_129, skk_131, sli0_105, \
                         sli0_107, sli1_105, sli1_107, slk_123, slk_129, \
                         slk_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_20 * skk_129[k]
                   + f_12 * sli0_105[k]
                   - f_13 * sli1_105[k]
                   + f_3 * pc_x[k] * slk_129[k];

        t_157[k] = f_3 * pc_z[k] * slk_123[k];

        t_158[k] = f_20 * skk_131[k]
                   + f_12 * sli0_107[k]
                   - f_13 * sli1_107[k]
                   + f_3 * pc_x[k] * slk_131[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, skk_56, skk_132, skk_133, sli0_108, \
                         sli0_109, sli1_108, sli1_109, slk_128, slk_132, \
                         slk_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_20 * skk_132[k]
                   + f_12 * sli0_108[k]
                   - f_13 * sli1_108[k]
                   + f_3 * pc_x[k] * slk_132[k];

        t_160[k] = f_20 * skk_133[k]
                   + f_12 * sli0_109[k]
                   - f_13 * sli1_109[k]
                   + f_3 * pc_x[k] * slk_133[k];

        t_161[k] = f_16 * skk_56[k]
                   + f_3 * pc_y[k] * slk_128[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, skk_135, skk_136, skk_137, skk_138, \
                         sli0_111, sli1_111, slk_135, slk_136, slk_137, \
                         slk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_20 * skk_135[k]
                   + f_12 * sli0_111[k]
                   - f_13 * sli1_111[k]
                   + f_3 * pc_x[k] * slk_135[k];

        t_163[k] = f_20 * skk_136[k]
                   + f_3 * pc_x[k] * slk_136[k];

        t_164[k] = f_20 * skk_137[k]
                   + f_3 * pc_x[k] * slk_137[k];

        t_165[k] = f_20 * skk_138[k]
                   + f_3 * pc_x[k] * slk_138[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, t_170, pc_x, skk_139, skk_140, skk_141, \
                         skk_142, skk_143, slk_139, slk_140, slk_141, slk_142, \
                         slk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_20 * skk_139[k]
                   + f_3 * pc_x[k] * slk_139[k];

        t_167[k] = f_20 * skk_140[k]
                   + f_3 * pc_x[k] * slk_140[k];

        t_168[k] = f_20 * skk_141[k]
                   + f_3 * pc_x[k] * slk_141[k];

        t_169[k] = f_20 * skk_142[k]
                   + f_3 * pc_x[k] * slk_142[k];

        t_170[k] = f_20 * skk_143[k]
                   + f_3 * pc_x[k] * slk_143[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pc_y, pc_z, skk_64, skk_66, sli0_105, sli0_107, \
                         sli1_105, sli1_107, slk_136, slk_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_16 * skk_64[k]
                   + f_1 * sli0_105[k]
                   - f_2 * sli1_105[k]
                   + f_3 * pc_y[k] * slk_136[k];

        t_172[k] = f_3 * pc_z[k] * slk_136[k];

        t_173[k] = f_16 * skk_66[k]
                   + f_4 * sli0_107[k]
                   - f_5 * sli1_107[k]
                   + f_3 * pc_y[k] * slk_138[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pc_y, skk_67, skk_68, skk_69, sli0_108, \
                         sli0_109, sli0_110, sli1_108, sli1_109, sli1_110, slk_139, slk_140, \
                         slk_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_16 * skk_67[k]
                   + f_6 * sli0_108[k]
                   - f_7 * sli1_108[k]
                   + f_3 * pc_y[k] * slk_139[k];

        t_175[k] = f_16 * skk_68[k]
                   + f_8 * sli0_109[k]
                   - f_9 * sli1_109[k]
                   + f_3 * pc_y[k] * slk_140[k];

        t_176[k] = f_16 * skk_69[k]
                   + f_10 * sli0_110[k]
                   - f_11 * sli1_110[k]
                   + f_3 * pc_y[k] * slk_141[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pb_y, pc_y, pc_z, skl0_90, skk_70, \
                         skk_71, skl1_90, sli0_111, sli1_111, slk_142, \
                         slk_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_16 * skk_70[k]
                   + f_12 * sli0_111[k]
                   - f_13 * sli1_111[k]
                   + f_3 * pc_y[k] * slk_142[k];

        t_178[k] = f_16 * skk_71[k]
                   + f_3 * pc_y[k] * slk_143[k];

        t_179[k] = f_1 * sli0_111[k]
                   - f_2 * sli1_111[k]
                   + f_3 * pc_z[k] * slk_143[k];

        t_180[k] = pb_y[k] * skl0_90[k]
                   - f_14 * pc_y[k] * skl1_90[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_y, pc_z, skl0_48, skk_36, \
                         skk_72, skk_74, skl1_48, slk_144, slk_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_15 * skk_72[k]
                   + f_3 * pc_y[k] * slk_144[k];

        t_182[k] = f_15 * skk_36[k]
                   + f_3 * pc_z[k] * slk_144[k];

        t_183[k] = pb_z[k] * skl0_48[k]
                   - f_14 * pc_z[k] * skl1_48[k];

        t_184[k] = f_15 * skk_74[k]
                   + f_3 * pc_y[k] * slk_146[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pb_y, pb_z, pc_y, pc_z, skl0_51, skl0_95, \
                         skk_39, skk_77, skl1_51, skl1_95, slk_147, \
                         slk_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pb_y[k] * skl0_95[k]
                   - f_14 * pc_y[k] * skl1_95[k];

        t_186[k] = pb_z[k] * skl0_51[k]
                   - f_14 * pc_z[k] * skl1_51[k];

        t_187[k] = f_15 * skk_39[k]
                   + f_3 * pc_z[k] * slk_147[k];

        t_188[k] = f_15 * skk_77[k]
                   + f_3 * pc_y[k] * slk_149[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_y, pb_z, pc_y, pc_z, skl0_55, skl0_99, \
                         skk_42, skl1_55, skl1_99, slk_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = pb_y[k] * skl0_99[k]
                   - f_14 * pc_y[k] * skl1_99[k];

        t_190[k] = pb_z[k] * skl0_55[k]
                   - f_14 * pc_z[k] * skl1_55[k];

        t_191[k] = f_15 * skk_42[k]
                   + f_3 * pc_z[k] * slk_150[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pb_y, pc_y, skl0_102, skl0_104, skk_80, skk_81, \
                         skl1_102, skl1_104, slk_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pb_y[k] * skl0_102[k]
                   + f_16 * skk_80[k]
                   - f_14 * pc_y[k] * skl1_102[k];

        t_193[k] = f_15 * skk_81[k]
                   + f_3 * pc_y[k] * slk_153[k];

        t_194[k] = pb_y[k] * skl0_104[k]
                   - f_14 * pc_y[k] * skl1_104[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_y, pb_z, pc_y, pc_z, skl0_60, skl0_107, \
                         skk_46, skk_84, skl1_60, skl1_107, slk_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_z[k] * skl0_60[k]
                   - f_14 * pc_z[k] * skl1_60[k];

        t_196[k] = f_15 * skk_46[k]
                   + f_3 * pc_z[k] * slk_154[k];

        t_197[k] = pb_y[k] * skl0_107[k]
                   + f_17 * skk_84[k]
                   - f_14 * pc_y[k] * skl1_107[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pb_y, pc_y, skl0_108, skl0_110, skk_85, skk_86, \
                         skl1_108, skl1_110, slk_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_y[k] * skl0_108[k]
                   + f_16 * skk_85[k]
                   - f_14 * pc_y[k] * skl1_108[k];

        t_199[k] = f_15 * skk_86[k]
                   + f_3 * pc_y[k] * slk_158[k];

        t_200[k] = pb_y[k] * skl0_110[k]
                   - f_14 * pc_y[k] * skl1_110[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_y, pb_z, pc_y, pc_z, skl0_66, skl0_113, \
                         skk_51, skk_89, skl1_66, skl1_113, slk_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = pb_z[k] * skl0_66[k]
                   - f_14 * pc_z[k] * skl1_66[k];

        t_202[k] = f_15 * skk_51[k]
                   + f_3 * pc_z[k] * slk_159[k];

        t_203[k] = pb_y[k] * skl0_113[k]
                   + f_18 * skk_89[k]
                   - f_14 * pc_y[k] * skl1_113[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pb_y, pc_y, skl0_114, skl0_115, skl0_117, \
                         skk_90, skk_91, skk_92, skl1_114, skl1_115, skl1_117, \
                         slk_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pb_y[k] * skl0_114[k]
                   + f_17 * skk_90[k]
                   - f_14 * pc_y[k] * skl1_114[k];

        t_205[k] = pb_y[k] * skl0_115[k]
                   + f_16 * skk_91[k]
                   - f_14 * pc_y[k] * skl1_115[k];

        t_206[k] = f_15 * skk_92[k]
                   + f_3 * pc_y[k] * slk_164[k];

        t_207[k] = pb_y[k] * skl0_117[k]
                   - f_14 * pc_y[k] * skl1_117[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, t_212, pc_x, skk_172, skk_173, skk_174, \
                         skk_175, skk_176, slk_172, slk_173, slk_174, slk_175, \
                         slk_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_20 * skk_172[k]
                   + f_3 * pc_x[k] * slk_172[k];

        t_209[k] = f_20 * skk_173[k]
                   + f_3 * pc_x[k] * slk_173[k];

        t_210[k] = f_20 * skk_174[k]
                   + f_3 * pc_x[k] * slk_174[k];

        t_211[k] = f_20 * skk_175[k]
                   + f_3 * pc_x[k] * slk_175[k];

        t_212[k] = f_20 * skk_176[k]
                   + f_3 * pc_x[k] * slk_176[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_z, pc_x, pc_z, skl0_81, skk_177, \
                         skk_178, skk_179, skl1_81, slk_177, slk_178, \
                         slk_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_20 * skk_177[k]
                   + f_3 * pc_x[k] * slk_177[k];

        t_214[k] = f_20 * skk_178[k]
                   + f_3 * pc_x[k] * slk_178[k];

        t_215[k] = f_20 * skk_179[k]
                   + f_3 * pc_x[k] * slk_179[k];

        t_216[k] = pb_z[k] * skl0_81[k]
                   - f_14 * pc_z[k] * skl1_81[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pc_y, pc_z, skk_64, skk_102, skk_103, sli0_135, \
                         sli0_136, sli1_135, sli1_136, slk_172, slk_174, \
                         slk_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_15 * skk_64[k]
                   + f_3 * pc_z[k] * slk_172[k];

        t_218[k] = f_15 * skk_102[k]
                   + f_4 * sli0_135[k]
                   - f_5 * sli1_135[k]
                   + f_3 * pc_y[k] * slk_174[k];

        t_219[k] = f_15 * skk_103[k]
                   + f_6 * sli0_136[k]
                   - f_7 * sli1_136[k]
                   + f_3 * pc_y[k] * slk_175[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, pc_y, skk_104, skk_105, skk_106, sli0_137, \
                         sli0_138, sli0_139, sli1_137, sli1_138, sli1_139, slk_176, slk_177, \
                         slk_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_15 * skk_104[k]
                   + f_8 * sli0_137[k]
                   - f_9 * sli1_137[k]
                   + f_3 * pc_y[k] * slk_176[k];

        t_221[k] = f_15 * skk_105[k]
                   + f_10 * sli0_138[k]
                   - f_11 * sli1_138[k]
                   + f_3 * pc_y[k] * slk_177[k];

        t_222[k] = f_15 * skk_106[k]
                   + f_12 * sli0_139[k]
                   - f_13 * sli1_139[k]
                   + f_3 * pc_y[k] * slk_178[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pb_y, pc_x, pc_y, skl0_134, skk_107, \
                         skk_180, skl1_134, sli0_140, sli1_140, slk_179, \
                         slk_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_15 * skk_107[k]
                   + f_3 * pc_y[k] * slk_179[k];

        t_224[k] = pb_y[k] * skl0_134[k]
                   - f_14 * pc_y[k] * skl1_134[k];

        t_225[k] = f_20 * skk_180[k]
                   + f_1 * sli0_140[k]
                   - f_2 * sli1_140[k]
                   + f_3 * pc_x[k] * slk_180[k];

        t_226[k] = f_3 * pc_y[k] * slk_180[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, skk_72, skk_183, sli0_143, \
                         sli1_143, slk_180, slk_182, slk_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_16 * skk_72[k]
                   + f_3 * pc_z[k] * slk_180[k];

        t_228[k] = f_20 * skk_183[k]
                   + f_4 * sli0_143[k]
                   - f_5 * sli1_143[k]
                   + f_3 * pc_x[k] * slk_183[k];

        t_229[k] = f_3 * pc_y[k] * slk_182[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pc_x, pc_z, skk_75, skk_185, skk_186, sli0_145, \
                         sli0_146, sli1_145, sli1_146, slk_183, slk_185, \
                         slk_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_20 * skk_185[k]
                   + f_4 * sli0_145[k]
                   - f_5 * sli1_145[k]
                   + f_3 * pc_x[k] * slk_185[k];

        t_231[k] = f_20 * skk_186[k]
                   + f_6 * sli0_146[k]
                   - f_7 * sli1_146[k]
                   + f_3 * pc_x[k] * slk_186[k];

        t_232[k] = f_16 * skk_75[k]
                   + f_3 * pc_z[k] * slk_183[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_135 = buffer.data(skl0 + 135);
    const auto *skl0_138 = buffer.data(skl0 + 138);
    const auto *skl0_141 = buffer.data(skl0 + 141);
    const auto *skl0_145 = buffer.data(skl0 + 145);
    const auto *skl0_147 = buffer.data(skl0 + 147);
    const auto *skl0_150 = buffer.data(skl0 + 150);
    const auto *skl0_152 = buffer.data(skl0 + 152);
    const auto *skl0_153 = buffer.data(skl0 + 153);
    const auto *skl0_156 = buffer.data(skl0 + 156);
    const auto *skl0_158 = buffer.data(skl0 + 158);
    const auto *skl0_159 = buffer.data(skl0 + 159);
    const auto *skl0_160 = buffer.data(skl0 + 160);

    const auto *skk_78 = buffer.data(skk + 78);
    const auto *skk_82 = buffer.data(skk + 82);
    const auto *skk_87 = buffer.data(skk + 87);
    const auto *skk_100 = buffer.data(skk + 100);
    const auto *skk_107 = buffer.data(skk + 107);
    const auto *skk_108 = buffer.data(skk + 108);
    const auto *skk_110 = buffer.data(skk + 110);
    const auto *skk_111 = buffer.data(skk + 111);
    const auto *skk_113 = buffer.data(skk + 113);
    const auto *skk_114 = buffer.data(skk + 114);
    const auto *skk_115 = buffer.data(skk + 115);
    const auto *skk_117 = buffer.data(skk + 117);
    const auto *skk_118 = buffer.data(skk + 118);
    const auto *skk_119 = buffer.data(skk + 119);
    const auto *skk_120 = buffer.data(skk + 120);
    const auto *skk_122 = buffer.data(skk + 122);
    const auto *skk_123 = buffer.data(skk + 123);
    const auto *skk_124 = buffer.data(skk + 124);
    const auto *skk_125 = buffer.data(skk + 125);
    const auto *skk_126 = buffer.data(skk + 126);
    const auto *skk_128 = buffer.data(skk + 128);
    const auto *skk_136 = buffer.data(skk + 136);
    const auto *skk_138 = buffer.data(skk + 138);
    const auto *skk_139 = buffer.data(skk + 139);
    const auto *skk_140 = buffer.data(skk + 140);
    const auto *skk_141 = buffer.data(skk + 141);
    const auto *skk_142 = buffer.data(skk + 142);
    const auto *skk_143 = buffer.data(skk + 143);
    const auto *skk_144 = buffer.data(skk + 144);
    const auto *skk_146 = buffer.data(skk + 146);
    const auto *skk_149 = buffer.data(skk + 149);
    const auto *skk_153 = buffer.data(skk + 153);
    const auto *skk_158 = buffer.data(skk + 158);
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
    const auto *skk_257 = buffer.data(skk + 257);
    const auto *skk_261 = buffer.data(skk + 261);
    const auto *skk_266 = buffer.data(skk + 266);
    const auto *skk_272 = buffer.data(skk + 272);

    const auto *skl1_135 = buffer.data(skl1 + 135);
    const auto *skl1_138 = buffer.data(skl1 + 138);
    const auto *skl1_141 = buffer.data(skl1 + 141);
    const auto *skl1_145 = buffer.data(skl1 + 145);
    const auto *skl1_147 = buffer.data(skl1 + 147);
    const auto *skl1_150 = buffer.data(skl1 + 150);
    const auto *skl1_152 = buffer.data(skl1 + 152);
    const auto *skl1_153 = buffer.data(skl1 + 153);
    const auto *skl1_156 = buffer.data(skl1 + 156);
    const auto *skl1_158 = buffer.data(skl1 + 158);
    const auto *skl1_159 = buffer.data(skl1 + 159);
    const auto *skl1_160 = buffer.data(skl1 + 160);

    const auto *sli0_149 = buffer.data(sli0 + 149);
    const auto *sli0_150 = buffer.data(sli0 + 150);
    const auto *sli0_152 = buffer.data(sli0 + 152);
    const auto *sli0_154 = buffer.data(sli0 + 154);
    const auto *sli0_155 = buffer.data(sli0 + 155);
    const auto *sli0_157 = buffer.data(sli0 + 157);
    const auto *sli0_158 = buffer.data(sli0 + 158);
    const auto *sli0_160 = buffer.data(sli0 + 160);
    const auto *sli0_161 = buffer.data(sli0 + 161);
    const auto *sli0_163 = buffer.data(sli0 + 163);
    const auto *sli0_164 = buffer.data(sli0 + 164);
    const auto *sli0_165 = buffer.data(sli0 + 165);
    const auto *sli0_166 = buffer.data(sli0 + 166);
    const auto *sli0_167 = buffer.data(sli0 + 167);
    const auto *sli0_168 = buffer.data(sli0 + 168);
    const auto *sli0_171 = buffer.data(sli0 + 171);
    const auto *sli0_173 = buffer.data(sli0 + 173);
    const auto *sli0_174 = buffer.data(sli0 + 174);
    const auto *sli0_177 = buffer.data(sli0 + 177);
    const auto *sli0_178 = buffer.data(sli0 + 178);
    const auto *sli0_180 = buffer.data(sli0 + 180);
    const auto *sli0_182 = buffer.data(sli0 + 182);
    const auto *sli0_183 = buffer.data(sli0 + 183);
    const auto *sli0_185 = buffer.data(sli0 + 185);
    const auto *sli0_186 = buffer.data(sli0 + 186);
    const auto *sli0_188 = buffer.data(sli0 + 188);
    const auto *sli0_189 = buffer.data(sli0 + 189);
    const auto *sli0_191 = buffer.data(sli0 + 191);
    const auto *sli0_192 = buffer.data(sli0 + 192);
    const auto *sli0_193 = buffer.data(sli0 + 193);
    const auto *sli0_194 = buffer.data(sli0 + 194);
    const auto *sli0_195 = buffer.data(sli0 + 195);
    const auto *sli0_201 = buffer.data(sli0 + 201);
    const auto *sli0_205 = buffer.data(sli0 + 205);
    const auto *sli0_210 = buffer.data(sli0 + 210);
    const auto *sli0_216 = buffer.data(sli0 + 216);

    const auto *sli1_149 = buffer.data(sli1 + 149);
    const auto *sli1_150 = buffer.data(sli1 + 150);
    const auto *sli1_152 = buffer.data(sli1 + 152);
    const auto *sli1_154 = buffer.data(sli1 + 154);
    const auto *sli1_155 = buffer.data(sli1 + 155);
    const auto *sli1_157 = buffer.data(sli1 + 157);
    const auto *sli1_158 = buffer.data(sli1 + 158);
    const auto *sli1_160 = buffer.data(sli1 + 160);
    const auto *sli1_161 = buffer.data(sli1 + 161);
    const auto *sli1_163 = buffer.data(sli1 + 163);
    const auto *sli1_164 = buffer.data(sli1 + 164);
    const auto *sli1_165 = buffer.data(sli1 + 165);
    const auto *sli1_166 = buffer.data(sli1 + 166);
    const auto *sli1_167 = buffer.data(sli1 + 167);
    const auto *sli1_168 = buffer.data(sli1 + 168);
    const auto *sli1_171 = buffer.data(sli1 + 171);
    const auto *sli1_173 = buffer.data(sli1 + 173);
    const auto *sli1_174 = buffer.data(sli1 + 174);
    const auto *sli1_177 = buffer.data(sli1 + 177);
    const auto *sli1_178 = buffer.data(sli1 + 178);
    const auto *sli1_180 = buffer.data(sli1 + 180);
    const auto *sli1_182 = buffer.data(sli1 + 182);
    const auto *sli1_183 = buffer.data(sli1 + 183);
    const auto *sli1_185 = buffer.data(sli1 + 185);
    const auto *sli1_186 = buffer.data(sli1 + 186);
    const auto *sli1_188 = buffer.data(sli1 + 188);
    const auto *sli1_189 = buffer.data(sli1 + 189);
    const auto *sli1_191 = buffer.data(sli1 + 191);
    const auto *sli1_192 = buffer.data(sli1 + 192);
    const auto *sli1_193 = buffer.data(sli1 + 193);
    const auto *sli1_194 = buffer.data(sli1 + 194);
    const auto *sli1_195 = buffer.data(sli1 + 195);
    const auto *sli1_201 = buffer.data(sli1 + 201);
    const auto *sli1_205 = buffer.data(sli1 + 205);
    const auto *sli1_210 = buffer.data(sli1 + 210);
    const auto *sli1_216 = buffer.data(sli1 + 216);

    const auto *slk_185 = buffer.data(slk + 185);
    const auto *slk_186 = buffer.data(slk + 186);
    const auto *slk_189 = buffer.data(slk + 189);
    const auto *slk_190 = buffer.data(slk + 190);
    const auto *slk_192 = buffer.data(slk + 192);
    const auto *slk_194 = buffer.data(slk + 194);
    const auto *slk_195 = buffer.data(slk + 195);
    const auto *slk_197 = buffer.data(slk + 197);
    const auto *slk_198 = buffer.data(slk + 198);
    const auto *slk_200 = buffer.data(slk + 200);
    const auto *slk_201 = buffer.data(slk + 201);
    const auto *slk_203 = buffer.data(slk + 203);
    const auto *slk_204 = buffer.data(slk + 204);
    const auto *slk_205 = buffer.data(slk + 205);
    const auto *slk_207 = buffer.data(slk + 207);
    const auto *slk_208 = buffer.data(slk + 208);
    const auto *slk_209 = buffer.data(slk + 209);
    const auto *slk_210 = buffer.data(slk + 210);
    const auto *slk_211 = buffer.data(slk + 211);
    const auto *slk_212 = buffer.data(slk + 212);
    const auto *slk_213 = buffer.data(slk + 213);
    const auto *slk_214 = buffer.data(slk + 214);
    const auto *slk_215 = buffer.data(slk + 215);
    const auto *slk_216 = buffer.data(slk + 216);
    const auto *slk_218 = buffer.data(slk + 218);
    const auto *slk_219 = buffer.data(slk + 219);
    const auto *slk_221 = buffer.data(slk + 221);
    const auto *slk_222 = buffer.data(slk + 222);
    const auto *slk_225 = buffer.data(slk + 225);
    const auto *slk_226 = buffer.data(slk + 226);
    const auto *slk_228 = buffer.data(slk + 228);
    const auto *slk_230 = buffer.data(slk + 230);
    const auto *slk_231 = buffer.data(slk + 231);
    const auto *slk_233 = buffer.data(slk + 233);
    const auto *slk_234 = buffer.data(slk + 234);
    const auto *slk_236 = buffer.data(slk + 236);
    const auto *slk_237 = buffer.data(slk + 237);
    const auto *slk_239 = buffer.data(slk + 239);
    const auto *slk_240 = buffer.data(slk + 240);
    const auto *slk_241 = buffer.data(slk + 241);
    const auto *slk_243 = buffer.data(slk + 243);
    const auto *slk_244 = buffer.data(slk + 244);
    const auto *slk_245 = buffer.data(slk + 245);
    const auto *slk_246 = buffer.data(slk + 246);
    const auto *slk_247 = buffer.data(slk + 247);
    const auto *slk_248 = buffer.data(slk + 248);
    const auto *slk_249 = buffer.data(slk + 249);
    const auto *slk_250 = buffer.data(slk + 250);
    const auto *slk_251 = buffer.data(slk + 251);
    const auto *slk_252 = buffer.data(slk + 252);
    const auto *slk_254 = buffer.data(slk + 254);
    const auto *slk_255 = buffer.data(slk + 255);
    const auto *slk_257 = buffer.data(slk + 257);
    const auto *slk_258 = buffer.data(slk + 258);
    const auto *slk_261 = buffer.data(slk + 261);
    const auto *slk_262 = buffer.data(slk + 262);
    const auto *slk_266 = buffer.data(slk + 266);
    const auto *slk_267 = buffer.data(slk + 267);
    const auto *slk_272 = buffer.data(slk + 272);

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, skk_189, skk_190, sli0_149, \
                         sli0_150, sli1_149, sli1_150, slk_185, slk_189, \
                         slk_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_3 * pc_y[k] * slk_185[k];

        t_234[k] = f_20 * skk_189[k]
                   + f_6 * sli0_149[k]
                   - f_7 * sli1_149[k]
                   + f_3 * pc_x[k] * slk_189[k];

        t_235[k] = f_20 * skk_190[k]
                   + f_8 * sli0_150[k]
                   - f_9 * sli1_150[k]
                   + f_3 * pc_x[k] * slk_190[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pc_x, pc_y, pc_z, skk_78, skk_192, sli0_152, \
                         sli1_152, slk_186, slk_189, slk_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_16 * skk_78[k]
                   + f_3 * pc_z[k] * slk_186[k];

        t_237[k] = f_20 * skk_192[k]
                   + f_8 * sli0_152[k]
                   - f_9 * sli1_152[k]
                   + f_3 * pc_x[k] * slk_192[k];

        t_238[k] = f_3 * pc_y[k] * slk_189[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, pc_x, pc_z, skk_82, skk_194, skk_195, sli0_154, \
                         sli0_155, sli1_154, sli1_155, slk_190, slk_194, \
                         slk_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_20 * skk_194[k]
                   + f_8 * sli0_154[k]
                   - f_9 * sli1_154[k]
                   + f_3 * pc_x[k] * slk_194[k];

        t_240[k] = f_20 * skk_195[k]
                   + f_10 * sli0_155[k]
                   - f_11 * sli1_155[k]
                   + f_3 * pc_x[k] * slk_195[k];

        t_241[k] = f_16 * skk_82[k]
                   + f_3 * pc_z[k] * slk_190[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pc_x, pc_y, skk_197, skk_198, sli0_157, \
                         sli0_158, sli1_157, sli1_158, slk_194, slk_197, \
                         slk_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_20 * skk_197[k]
                   + f_10 * sli0_157[k]
                   - f_11 * sli1_157[k]
                   + f_3 * pc_x[k] * slk_197[k];

        t_243[k] = f_20 * skk_198[k]
                   + f_10 * sli0_158[k]
                   - f_11 * sli1_158[k]
                   + f_3 * pc_x[k] * slk_198[k];

        t_244[k] = f_3 * pc_y[k] * slk_194[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_x, pc_z, skk_87, skk_200, skk_201, sli0_160, \
                         sli0_161, sli1_160, sli1_161, slk_195, slk_200, \
                         slk_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_20 * skk_200[k]
                   + f_10 * sli0_160[k]
                   - f_11 * sli1_160[k]
                   + f_3 * pc_x[k] * slk_200[k];

        t_246[k] = f_20 * skk_201[k]
                   + f_12 * sli0_161[k]
                   - f_13 * sli1_161[k]
                   + f_3 * pc_x[k] * slk_201[k];

        t_247[k] = f_16 * skk_87[k]
                   + f_3 * pc_z[k] * slk_195[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pc_x, skk_203, skk_204, skk_205, sli0_163, \
                         sli0_164, sli0_165, sli1_163, sli1_164, sli1_165, slk_203, slk_204, \
                         slk_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_20 * skk_203[k]
                   + f_12 * sli0_163[k]
                   - f_13 * sli1_163[k]
                   + f_3 * pc_x[k] * slk_203[k];

        t_249[k] = f_20 * skk_204[k]
                   + f_12 * sli0_164[k]
                   - f_13 * sli1_164[k]
                   + f_3 * pc_x[k] * slk_204[k];

        t_250[k] = f_20 * skk_205[k]
                   + f_12 * sli0_165[k]
                   - f_13 * sli1_165[k]
                   + f_3 * pc_x[k] * slk_205[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pc_x, pc_y, skk_207, skk_208, skk_209, \
                         sli0_167, sli1_167, slk_200, slk_207, slk_208, \
                         slk_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_3 * pc_y[k] * slk_200[k];

        t_252[k] = f_20 * skk_207[k]
                   + f_12 * sli0_167[k]
                   - f_13 * sli1_167[k]
                   + f_3 * pc_x[k] * slk_207[k];

        t_253[k] = f_20 * skk_208[k]
                   + f_3 * pc_x[k] * slk_208[k];

        t_254[k] = f_20 * skk_209[k]
                   + f_3 * pc_x[k] * slk_209[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, skk_210, skk_211, skk_212, \
                         skk_213, skk_214, slk_210, slk_211, slk_212, slk_213, \
                         slk_214 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_20 * skk_210[k]
                   + f_3 * pc_x[k] * slk_210[k];

        t_256[k] = f_20 * skk_211[k]
                   + f_3 * pc_x[k] * slk_211[k];

        t_257[k] = f_20 * skk_212[k]
                   + f_3 * pc_x[k] * slk_212[k];

        t_258[k] = f_20 * skk_213[k]
                   + f_3 * pc_x[k] * slk_213[k];

        t_259[k] = f_20 * skk_214[k]
                   + f_3 * pc_x[k] * slk_214[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, pc_z, skk_100, skk_215, \
                         sli0_161, sli0_163, sli1_161, sli1_163, slk_208, slk_210, \
                         slk_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_20 * skk_215[k]
                   + f_3 * pc_x[k] * slk_215[k];

        t_261[k] = f_1 * sli0_161[k]
                   - f_2 * sli1_161[k]
                   + f_3 * pc_y[k] * slk_208[k];

        t_262[k] = f_16 * skk_100[k]
                   + f_3 * pc_z[k] * slk_208[k];

        t_263[k] = f_4 * sli0_163[k]
                   - f_5 * sli1_163[k]
                   + f_3 * pc_y[k] * slk_210[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, sli0_164, sli0_165, sli0_166, sli1_164, \
                         sli1_165, sli1_166, slk_211, slk_212, \
                         slk_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_6 * sli0_164[k]
                   - f_7 * sli1_164[k]
                   + f_3 * pc_y[k] * slk_211[k];

        t_265[k] = f_8 * sli0_165[k]
                   - f_9 * sli1_165[k]
                   + f_3 * pc_y[k] * slk_212[k];

        t_266[k] = f_10 * sli0_166[k]
                   - f_11 * sli1_166[k]
                   + f_3 * pc_y[k] * slk_213[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_y, pc_z, skk_107, skk_216, \
                         sli0_167, sli0_168, sli1_167, sli1_168, slk_214, slk_215, \
                         slk_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_12 * sli0_167[k]
                   - f_13 * sli1_167[k]
                   + f_3 * pc_y[k] * slk_214[k];

        t_268[k] = f_3 * pc_y[k] * slk_215[k];

        t_269[k] = f_16 * skk_107[k]
                   + f_1 * sli0_167[k]
                   - f_2 * sli1_167[k]
                   + f_3 * pc_z[k] * slk_215[k];

        t_270[k] = f_19 * skk_216[k]
                   + f_1 * sli0_168[k]
                   - f_2 * sli1_168[k]
                   + f_3 * pc_x[k] * slk_216[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, skk_108, skk_110, \
                         skk_219, sli0_171, sli1_171, slk_216, slk_218, \
                         slk_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_17 * skk_108[k]
                   + f_3 * pc_y[k] * slk_216[k];

        t_272[k] = f_3 * pc_z[k] * slk_216[k];

        t_273[k] = f_19 * skk_219[k]
                   + f_4 * sli0_171[k]
                   - f_5 * sli1_171[k]
                   + f_3 * pc_x[k] * slk_219[k];

        t_274[k] = f_17 * skk_110[k]
                   + f_3 * pc_y[k] * slk_218[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pc_x, pc_z, skk_221, skk_222, sli0_173, \
                         sli0_174, sli1_173, sli1_174, slk_219, slk_221, \
                         slk_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_19 * skk_221[k]
                   + f_4 * sli0_173[k]
                   - f_5 * sli1_173[k]
                   + f_3 * pc_x[k] * slk_221[k];

        t_276[k] = f_19 * skk_222[k]
                   + f_6 * sli0_174[k]
                   - f_7 * sli1_174[k]
                   + f_3 * pc_x[k] * slk_222[k];

        t_277[k] = f_3 * pc_z[k] * slk_219[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pc_x, pc_y, skk_113, skk_225, skk_226, sli0_177, \
                         sli0_178, sli1_177, sli1_178, slk_221, slk_225, \
                         slk_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_17 * skk_113[k]
                   + f_3 * pc_y[k] * slk_221[k];

        t_279[k] = f_19 * skk_225[k]
                   + f_6 * sli0_177[k]
                   - f_7 * sli1_177[k]
                   + f_3 * pc_x[k] * slk_225[k];

        t_280[k] = f_19 * skk_226[k]
                   + f_8 * sli0_178[k]
                   - f_9 * sli1_178[k]
                   + f_3 * pc_x[k] * slk_226[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, pc_z, skk_117, skk_228, sli0_180, \
                         sli1_180, slk_222, slk_225, slk_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * pc_z[k] * slk_222[k];

        t_282[k] = f_19 * skk_228[k]
                   + f_8 * sli0_180[k]
                   - f_9 * sli1_180[k]
                   + f_3 * pc_x[k] * slk_228[k];

        t_283[k] = f_17 * skk_117[k]
                   + f_3 * pc_y[k] * slk_225[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, pc_x, pc_z, skk_230, skk_231, sli0_182, \
                         sli0_183, sli1_182, sli1_183, slk_226, slk_230, \
                         slk_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_19 * skk_230[k]
                   + f_8 * sli0_182[k]
                   - f_9 * sli1_182[k]
                   + f_3 * pc_x[k] * slk_230[k];

        t_285[k] = f_19 * skk_231[k]
                   + f_10 * sli0_183[k]
                   - f_11 * sli1_183[k]
                   + f_3 * pc_x[k] * slk_231[k];

        t_286[k] = f_3 * pc_z[k] * slk_226[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pc_x, pc_y, skk_122, skk_233, skk_234, sli0_185, \
                         sli0_186, sli1_185, sli1_186, slk_230, slk_233, \
                         slk_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_19 * skk_233[k]
                   + f_10 * sli0_185[k]
                   - f_11 * sli1_185[k]
                   + f_3 * pc_x[k] * slk_233[k];

        t_288[k] = f_19 * skk_234[k]
                   + f_10 * sli0_186[k]
                   - f_11 * sli1_186[k]
                   + f_3 * pc_x[k] * slk_234[k];

        t_289[k] = f_17 * skk_122[k]
                   + f_3 * pc_y[k] * slk_230[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pc_x, pc_z, skk_236, skk_237, sli0_188, \
                         sli0_189, sli1_188, sli1_189, slk_231, slk_236, \
                         slk_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_19 * skk_236[k]
                   + f_10 * sli0_188[k]
                   - f_11 * sli1_188[k]
                   + f_3 * pc_x[k] * slk_236[k];

        t_291[k] = f_19 * skk_237[k]
                   + f_12 * sli0_189[k]
                   - f_13 * sli1_189[k]
                   + f_3 * pc_x[k] * slk_237[k];

        t_292[k] = f_3 * pc_z[k] * slk_231[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pc_x, skk_239, skk_240, skk_241, sli0_191, \
                         sli0_192, sli0_193, sli1_191, sli1_192, sli1_193, slk_239, slk_240, \
                         slk_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_19 * skk_239[k]
                   + f_12 * sli0_191[k]
                   - f_13 * sli1_191[k]
                   + f_3 * pc_x[k] * slk_239[k];

        t_294[k] = f_19 * skk_240[k]
                   + f_12 * sli0_192[k]
                   - f_13 * sli1_192[k]
                   + f_3 * pc_x[k] * slk_240[k];

        t_295[k] = f_19 * skk_241[k]
                   + f_12 * sli0_193[k]
                   - f_13 * sli1_193[k]
                   + f_3 * pc_x[k] * slk_241[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, t_299, pc_x, pc_y, skk_128, skk_243, skk_244, \
                         skk_245, sli0_195, sli1_195, slk_236, slk_243, slk_244, \
                         slk_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_17 * skk_128[k]
                   + f_3 * pc_y[k] * slk_236[k];

        t_297[k] = f_19 * skk_243[k]
                   + f_12 * sli0_195[k]
                   - f_13 * sli1_195[k]
                   + f_3 * pc_x[k] * slk_243[k];

        t_298[k] = f_19 * skk_244[k]
                   + f_3 * pc_x[k] * slk_244[k];

        t_299[k] = f_19 * skk_245[k]
                   + f_3 * pc_x[k] * slk_245[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, t_304, pc_x, skk_246, skk_247, skk_248, \
                         skk_249, skk_250, slk_246, slk_247, slk_248, slk_249, \
                         slk_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_19 * skk_246[k]
                   + f_3 * pc_x[k] * slk_246[k];

        t_301[k] = f_19 * skk_247[k]
                   + f_3 * pc_x[k] * slk_247[k];

        t_302[k] = f_19 * skk_248[k]
                   + f_3 * pc_x[k] * slk_248[k];

        t_303[k] = f_19 * skk_249[k]
                   + f_3 * pc_x[k] * slk_249[k];

        t_304[k] = f_19 * skk_250[k]
                   + f_3 * pc_x[k] * slk_250[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pc_x, pc_y, pc_z, skk_136, skk_251, sli0_189, \
                         sli1_189, slk_244, slk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_19 * skk_251[k]
                   + f_3 * pc_x[k] * slk_251[k];

        t_306[k] = f_17 * skk_136[k]
                   + f_1 * sli0_189[k]
                   - f_2 * sli1_189[k]
                   + f_3 * pc_y[k] * slk_244[k];

        t_307[k] = f_3 * pc_z[k] * slk_244[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_y, skk_138, skk_139, skk_140, sli0_191, \
                         sli0_192, sli0_193, sli1_191, sli1_192, sli1_193, slk_246, slk_247, \
                         slk_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_17 * skk_138[k]
                   + f_4 * sli0_191[k]
                   - f_5 * sli1_191[k]
                   + f_3 * pc_y[k] * slk_246[k];

        t_309[k] = f_17 * skk_139[k]
                   + f_6 * sli0_192[k]
                   - f_7 * sli1_192[k]
                   + f_3 * pc_y[k] * slk_247[k];

        t_310[k] = f_17 * skk_140[k]
                   + f_8 * sli0_193[k]
                   - f_9 * sli1_193[k]
                   + f_3 * pc_y[k] * slk_248[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_y, pc_z, skk_141, skk_142, skk_143, \
                         sli0_194, sli0_195, sli1_194, sli1_195, slk_249, slk_250, \
                         slk_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_17 * skk_141[k]
                   + f_10 * sli0_194[k]
                   - f_11 * sli1_194[k]
                   + f_3 * pc_y[k] * slk_249[k];

        t_312[k] = f_17 * skk_142[k]
                   + f_12 * sli0_195[k]
                   - f_13 * sli1_195[k]
                   + f_3 * pc_y[k] * slk_250[k];

        t_313[k] = f_17 * skk_143[k]
                   + f_3 * pc_y[k] * slk_251[k];

        t_314[k] = f_1 * sli0_195[k]
                   - f_2 * sli1_195[k]
                   + f_3 * pc_z[k] * slk_251[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pb_z, pc_y, pc_z, skl0_135, skl0_138, \
                         skk_108, skk_144, skl1_135, skl1_138, \
                         slk_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = pb_z[k] * skl0_135[k]
                   - f_14 * pc_z[k] * skl1_135[k];

        t_316[k] = f_16 * skk_144[k]
                   + f_3 * pc_y[k] * slk_252[k];

        t_317[k] = f_15 * skk_108[k]
                   + f_3 * pc_z[k] * slk_252[k];

        t_318[k] = pb_z[k] * skl0_138[k]
                   - f_14 * pc_z[k] * skl1_138[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pb_z, pc_x, pc_y, pc_z, skl0_141, skk_146, \
                         skk_257, skl1_141, sli0_201, sli1_201, slk_254, \
                         slk_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_16 * skk_146[k]
                   + f_3 * pc_y[k] * slk_254[k];

        t_320[k] = f_19 * skk_257[k]
                   + f_4 * sli0_201[k]
                   - f_5 * sli1_201[k]
                   + f_3 * pc_x[k] * slk_257[k];

        t_321[k] = pb_z[k] * skl0_141[k]
                   - f_14 * pc_z[k] * skl1_141[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pc_x, pc_y, pc_z, skk_111, skk_149, skk_261, \
                         sli0_205, sli1_205, slk_255, slk_257, \
                         slk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_15 * skk_111[k]
                   + f_3 * pc_z[k] * slk_255[k];

        t_323[k] = f_16 * skk_149[k]
                   + f_3 * pc_y[k] * slk_257[k];

        t_324[k] = f_19 * skk_261[k]
                   + f_6 * sli0_205[k]
                   - f_7 * sli1_205[k]
                   + f_3 * pc_x[k] * slk_261[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, pb_z, pc_y, pc_z, skl0_145, skl0_147, \
                         skk_114, skk_115, skk_153, skl1_145, skl1_147, slk_258, \
                         slk_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pb_z[k] * skl0_145[k]
                   - f_14 * pc_z[k] * skl1_145[k];

        t_326[k] = f_15 * skk_114[k]
                   + f_3 * pc_z[k] * slk_258[k];

        t_327[k] = pb_z[k] * skl0_147[k]
                   + f_16 * skk_115[k]
                   - f_14 * pc_z[k] * skl1_147[k];

        t_328[k] = f_16 * skk_153[k]
                   + f_3 * pc_y[k] * slk_261[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, pb_z, pc_x, pc_z, skl0_150, skk_118, skk_266, \
                         skl1_150, sli0_210, sli1_210, slk_262, \
                         slk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_19 * skk_266[k]
                   + f_8 * sli0_210[k]
                   - f_9 * sli1_210[k]
                   + f_3 * pc_x[k] * slk_266[k];

        t_330[k] = pb_z[k] * skl0_150[k]
                   - f_14 * pc_z[k] * skl1_150[k];

        t_331[k] = f_15 * skk_118[k]
                   + f_3 * pc_z[k] * slk_262[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pb_z, pc_y, pc_z, skl0_152, skl0_153, skk_119, \
                         skk_120, skk_158, skl1_152, skl1_153, \
                         slk_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pb_z[k] * skl0_152[k]
                   + f_16 * skk_119[k]
                   - f_14 * pc_z[k] * skl1_152[k];

        t_333[k] = pb_z[k] * skl0_153[k]
                   + f_17 * skk_120[k]
                   - f_14 * pc_z[k] * skl1_153[k];

        t_334[k] = f_16 * skk_158[k]
                   + f_3 * pc_y[k] * slk_266[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, pb_z, pc_x, pc_z, skl0_156, skk_123, skk_272, \
                         skl1_156, sli0_216, sli1_216, slk_267, \
                         slk_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_19 * skk_272[k]
                   + f_10 * sli0_216[k]
                   - f_11 * sli1_216[k]
                   + f_3 * pc_x[k] * slk_272[k];

        t_336[k] = pb_z[k] * skl0_156[k]
                   - f_14 * pc_z[k] * skl1_156[k];

        t_337[k] = f_15 * skk_123[k]
                   + f_3 * pc_z[k] * slk_267[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, pb_z, pc_z, skl0_158, skl0_159, skl0_160, \
                         skk_124, skk_125, skk_126, skl1_158, skl1_159, \
                         skl1_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = pb_z[k] * skl0_158[k]
                   + f_16 * skk_124[k]
                   - f_14 * pc_z[k] * skl1_158[k];

        t_339[k] = pb_z[k] * skl0_159[k]
                   + f_17 * skk_125[k]
                   - f_14 * pc_z[k] * skl1_159[k];

        t_340[k] = pb_z[k] * skl0_160[k]
                   + f_18 * skk_126[k]
                   - f_14 * pc_z[k] * skl1_160[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_171 = buffer.data(skl0 + 171);
    const auto *skl0_225 = buffer.data(skl0 + 225);
    const auto *skl0_228 = buffer.data(skl0 + 228);
    const auto *skl0_230 = buffer.data(skl0 + 230);
    const auto *skl0_231 = buffer.data(skl0 + 231);
    const auto *skl0_234 = buffer.data(skl0 + 234);
    const auto *skl0_235 = buffer.data(skl0 + 235);
    const auto *skl0_237 = buffer.data(skl0 + 237);
    const auto *skl0_239 = buffer.data(skl0 + 239);
    const auto *skl0_240 = buffer.data(skl0 + 240);
    const auto *skl0_242 = buffer.data(skl0 + 242);
    const auto *skl0_243 = buffer.data(skl0 + 243);
    const auto *skl0_245 = buffer.data(skl0 + 245);
    const auto *skl0_246 = buffer.data(skl0 + 246);
    const auto *skl0_248 = buffer.data(skl0 + 248);
    const auto *skl0_249 = buffer.data(skl0 + 249);
    const auto *skl0_250 = buffer.data(skl0 + 250);
    const auto *skl0_252 = buffer.data(skl0 + 252);
    const auto *skl0_269 = buffer.data(skl0 + 269);

    const auto *skk_136 = buffer.data(skk + 136);
    const auto *skk_143 = buffer.data(skk + 143);
    const auto *skk_144 = buffer.data(skk + 144);
    const auto *skk_147 = buffer.data(skk + 147);
    const auto *skk_150 = buffer.data(skk + 150);
    const auto *skk_154 = buffer.data(skk + 154);
    const auto *skk_159 = buffer.data(skk + 159);
    const auto *skk_164 = buffer.data(skk + 164);
    const auto *skk_172 = buffer.data(skk + 172);
    const auto *skk_174 = buffer.data(skk + 174);
    const auto *skk_175 = buffer.data(skk + 175);
    const auto *skk_176 = buffer.data(skk + 176);
    const auto *skk_177 = buffer.data(skk + 177);
    const auto *skk_178 = buffer.data(skk + 178);
    const auto *skk_179 = buffer.data(skk + 179);
    const auto *skk_180 = buffer.data(skk + 180);
    const auto *skk_181 = buffer.data(skk + 181);
    const auto *skk_182 = buffer.data(skk + 182);
    const auto *skk_183 = buffer.data(skk + 183);
    const auto *skk_185 = buffer.data(skk + 185);
    const auto *skk_186 = buffer.data(skk + 186);
    const auto *skk_188 = buffer.data(skk + 188);
    const auto *skk_189 = buffer.data(skk + 189);
    const auto *skk_190 = buffer.data(skk + 190);
    const auto *skk_192 = buffer.data(skk + 192);
    const auto *skk_193 = buffer.data(skk + 193);
    const auto *skk_194 = buffer.data(skk + 194);
    const auto *skk_195 = buffer.data(skk + 195);
    const auto *skk_197 = buffer.data(skk + 197);
    const auto *skk_198 = buffer.data(skk + 198);
    const auto *skk_199 = buffer.data(skk + 199);
    const auto *skk_200 = buffer.data(skk + 200);
    const auto *skk_208 = buffer.data(skk + 208);
    const auto *skk_210 = buffer.data(skk + 210);
    const auto *skk_211 = buffer.data(skk + 211);
    const auto *skk_212 = buffer.data(skk + 212);
    const auto *skk_213 = buffer.data(skk + 213);
    const auto *skk_214 = buffer.data(skk + 214);
    const auto *skk_215 = buffer.data(skk + 215);
    const auto *skk_216 = buffer.data(skk + 216);
    const auto *skk_218 = buffer.data(skk + 218);
    const auto *skk_279 = buffer.data(skk + 279);
    const auto *skk_280 = buffer.data(skk + 280);
    const auto *skk_281 = buffer.data(skk + 281);
    const auto *skk_282 = buffer.data(skk + 282);
    const auto *skk_283 = buffer.data(skk + 283);
    const auto *skk_284 = buffer.data(skk + 284);
    const auto *skk_285 = buffer.data(skk + 285);
    const auto *skk_286 = buffer.data(skk + 286);
    const auto *skk_287 = buffer.data(skk + 287);
    const auto *skk_316 = buffer.data(skk + 316);
    const auto *skk_317 = buffer.data(skk + 317);
    const auto *skk_318 = buffer.data(skk + 318);
    const auto *skk_319 = buffer.data(skk + 319);
    const auto *skk_320 = buffer.data(skk + 320);
    const auto *skk_321 = buffer.data(skk + 321);
    const auto *skk_322 = buffer.data(skk + 322);
    const auto *skk_323 = buffer.data(skk + 323);
    const auto *skk_324 = buffer.data(skk + 324);
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
    const auto *skk_363 = buffer.data(skk + 363);
    const auto *skk_365 = buffer.data(skk + 365);

    const auto *skl1_171 = buffer.data(skl1 + 171);
    const auto *skl1_225 = buffer.data(skl1 + 225);
    const auto *skl1_228 = buffer.data(skl1 + 228);
    const auto *skl1_230 = buffer.data(skl1 + 230);
    const auto *skl1_231 = buffer.data(skl1 + 231);
    const auto *skl1_234 = buffer.data(skl1 + 234);
    const auto *skl1_235 = buffer.data(skl1 + 235);
    const auto *skl1_237 = buffer.data(skl1 + 237);
    const auto *skl1_239 = buffer.data(skl1 + 239);
    const auto *skl1_240 = buffer.data(skl1 + 240);
    const auto *skl1_242 = buffer.data(skl1 + 242);
    const auto *skl1_243 = buffer.data(skl1 + 243);
    const auto *skl1_245 = buffer.data(skl1 + 245);
    const auto *skl1_246 = buffer.data(skl1 + 246);
    const auto *skl1_248 = buffer.data(skl1 + 248);
    const auto *skl1_249 = buffer.data(skl1 + 249);
    const auto *skl1_250 = buffer.data(skl1 + 250);
    const auto *skl1_252 = buffer.data(skl1 + 252);
    const auto *skl1_269 = buffer.data(skl1 + 269);

    const auto *sli0_219 = buffer.data(sli0 + 219);
    const auto *sli0_220 = buffer.data(sli0 + 220);
    const auto *sli0_221 = buffer.data(sli0 + 221);
    const auto *sli0_222 = buffer.data(sli0 + 222);
    const auto *sli0_223 = buffer.data(sli0 + 223);
    const auto *sli0_245 = buffer.data(sli0 + 245);
    const auto *sli0_247 = buffer.data(sli0 + 247);
    const auto *sli0_248 = buffer.data(sli0 + 248);
    const auto *sli0_249 = buffer.data(sli0 + 249);
    const auto *sli0_250 = buffer.data(sli0 + 250);
    const auto *sli0_251 = buffer.data(sli0 + 251);
    const auto *sli0_252 = buffer.data(sli0 + 252);
    const auto *sli0_255 = buffer.data(sli0 + 255);
    const auto *sli0_257 = buffer.data(sli0 + 257);
    const auto *sli0_258 = buffer.data(sli0 + 258);
    const auto *sli0_261 = buffer.data(sli0 + 261);
    const auto *sli0_262 = buffer.data(sli0 + 262);
    const auto *sli0_264 = buffer.data(sli0 + 264);
    const auto *sli0_266 = buffer.data(sli0 + 266);
    const auto *sli0_267 = buffer.data(sli0 + 267);
    const auto *sli0_269 = buffer.data(sli0 + 269);
    const auto *sli0_270 = buffer.data(sli0 + 270);
    const auto *sli0_272 = buffer.data(sli0 + 272);
    const auto *sli0_273 = buffer.data(sli0 + 273);
    const auto *sli0_275 = buffer.data(sli0 + 275);
    const auto *sli0_276 = buffer.data(sli0 + 276);
    const auto *sli0_277 = buffer.data(sli0 + 277);
    const auto *sli0_278 = buffer.data(sli0 + 278);
    const auto *sli0_279 = buffer.data(sli0 + 279);
    const auto *sli0_280 = buffer.data(sli0 + 280);
    const auto *sli0_283 = buffer.data(sli0 + 283);
    const auto *sli0_285 = buffer.data(sli0 + 285);

    const auto *sli1_219 = buffer.data(sli1 + 219);
    const auto *sli1_220 = buffer.data(sli1 + 220);
    const auto *sli1_221 = buffer.data(sli1 + 221);
    const auto *sli1_222 = buffer.data(sli1 + 222);
    const auto *sli1_223 = buffer.data(sli1 + 223);
    const auto *sli1_245 = buffer.data(sli1 + 245);
    const auto *sli1_247 = buffer.data(sli1 + 247);
    const auto *sli1_248 = buffer.data(sli1 + 248);
    const auto *sli1_249 = buffer.data(sli1 + 249);
    const auto *sli1_250 = buffer.data(sli1 + 250);
    const auto *sli1_251 = buffer.data(sli1 + 251);
    const auto *sli1_252 = buffer.data(sli1 + 252);
    const auto *sli1_255 = buffer.data(sli1 + 255);
    const auto *sli1_257 = buffer.data(sli1 + 257);
    const auto *sli1_258 = buffer.data(sli1 + 258);
    const auto *sli1_261 = buffer.data(sli1 + 261);
    const auto *sli1_262 = buffer.data(sli1 + 262);
    const auto *sli1_264 = buffer.data(sli1 + 264);
    const auto *sli1_266 = buffer.data(sli1 + 266);
    const auto *sli1_267 = buffer.data(sli1 + 267);
    const auto *sli1_269 = buffer.data(sli1 + 269);
    const auto *sli1_270 = buffer.data(sli1 + 270);
    const auto *sli1_272 = buffer.data(sli1 + 272);
    const auto *sli1_273 = buffer.data(sli1 + 273);
    const auto *sli1_275 = buffer.data(sli1 + 275);
    const auto *sli1_276 = buffer.data(sli1 + 276);
    const auto *sli1_277 = buffer.data(sli1 + 277);
    const auto *sli1_278 = buffer.data(sli1 + 278);
    const auto *sli1_279 = buffer.data(sli1 + 279);
    const auto *sli1_280 = buffer.data(sli1 + 280);
    const auto *sli1_283 = buffer.data(sli1 + 283);
    const auto *sli1_285 = buffer.data(sli1 + 285);

    const auto *slk_272 = buffer.data(slk + 272);
    const auto *slk_279 = buffer.data(slk + 279);
    const auto *slk_280 = buffer.data(slk + 280);
    const auto *slk_281 = buffer.data(slk + 281);
    const auto *slk_282 = buffer.data(slk + 282);
    const auto *slk_283 = buffer.data(slk + 283);
    const auto *slk_284 = buffer.data(slk + 284);
    const auto *slk_285 = buffer.data(slk + 285);
    const auto *slk_286 = buffer.data(slk + 286);
    const auto *slk_287 = buffer.data(slk + 287);
    const auto *slk_288 = buffer.data(slk + 288);
    const auto *slk_290 = buffer.data(slk + 290);
    const auto *slk_291 = buffer.data(slk + 291);
    const auto *slk_293 = buffer.data(slk + 293);
    const auto *slk_294 = buffer.data(slk + 294);
    const auto *slk_297 = buffer.data(slk + 297);
    const auto *slk_298 = buffer.data(slk + 298);
    const auto *slk_302 = buffer.data(slk + 302);
    const auto *slk_303 = buffer.data(slk + 303);
    const auto *slk_308 = buffer.data(slk + 308);
    const auto *slk_316 = buffer.data(slk + 316);
    const auto *slk_317 = buffer.data(slk + 317);
    const auto *slk_318 = buffer.data(slk + 318);
    const auto *slk_319 = buffer.data(slk + 319);
    const auto *slk_320 = buffer.data(slk + 320);
    const auto *slk_321 = buffer.data(slk + 321);
    const auto *slk_322 = buffer.data(slk + 322);
    const auto *slk_323 = buffer.data(slk + 323);
    const auto *slk_324 = buffer.data(slk + 324);
    const auto *slk_326 = buffer.data(slk + 326);
    const auto *slk_327 = buffer.data(slk + 327);
    const auto *slk_329 = buffer.data(slk + 329);
    const auto *slk_330 = buffer.data(slk + 330);
    const auto *slk_333 = buffer.data(slk + 333);
    const auto *slk_334 = buffer.data(slk + 334);
    const auto *slk_336 = buffer.data(slk + 336);
    const auto *slk_338 = buffer.data(slk + 338);
    const auto *slk_339 = buffer.data(slk + 339);
    const auto *slk_341 = buffer.data(slk + 341);
    const auto *slk_342 = buffer.data(slk + 342);
    const auto *slk_344 = buffer.data(slk + 344);
    const auto *slk_345 = buffer.data(slk + 345);
    const auto *slk_347 = buffer.data(slk + 347);
    const auto *slk_348 = buffer.data(slk + 348);
    const auto *slk_349 = buffer.data(slk + 349);
    const auto *slk_351 = buffer.data(slk + 351);
    const auto *slk_352 = buffer.data(slk + 352);
    const auto *slk_353 = buffer.data(slk + 353);
    const auto *slk_354 = buffer.data(slk + 354);
    const auto *slk_355 = buffer.data(slk + 355);
    const auto *slk_356 = buffer.data(slk + 356);
    const auto *slk_357 = buffer.data(slk + 357);
    const auto *slk_358 = buffer.data(slk + 358);
    const auto *slk_359 = buffer.data(slk + 359);
    const auto *slk_360 = buffer.data(slk + 360);
    const auto *slk_362 = buffer.data(slk + 362);
    const auto *slk_363 = buffer.data(slk + 363);
    const auto *slk_365 = buffer.data(slk + 365);

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pc_x, pc_y, skk_164, skk_279, skk_280, \
                         skk_281, sli0_223, sli1_223, slk_272, slk_279, slk_280, \
                         slk_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_16 * skk_164[k]
                   + f_3 * pc_y[k] * slk_272[k];

        t_342[k] = f_19 * skk_279[k]
                   + f_12 * sli0_223[k]
                   - f_13 * sli1_223[k]
                   + f_3 * pc_x[k] * slk_279[k];

        t_343[k] = f_19 * skk_280[k]
                   + f_3 * pc_x[k] * slk_280[k];

        t_344[k] = f_19 * skk_281[k]
                   + f_3 * pc_x[k] * slk_281[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, t_349, pc_x, skk_282, skk_283, skk_284, \
                         skk_285, skk_286, slk_282, slk_283, slk_284, slk_285, \
                         slk_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_19 * skk_282[k]
                   + f_3 * pc_x[k] * slk_282[k];

        t_346[k] = f_19 * skk_283[k]
                   + f_3 * pc_x[k] * slk_283[k];

        t_347[k] = f_19 * skk_284[k]
                   + f_3 * pc_x[k] * slk_284[k];

        t_348[k] = f_19 * skk_285[k]
                   + f_3 * pc_x[k] * slk_285[k];

        t_349[k] = f_19 * skk_286[k]
                   + f_3 * pc_x[k] * slk_286[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, pb_z, pc_x, pc_z, skl0_171, skk_136, skk_287, \
                         skl1_171, slk_280, slk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_19 * skk_287[k]
                   + f_3 * pc_x[k] * slk_287[k];

        t_351[k] = pb_z[k] * skl0_171[k]
                   - f_14 * pc_z[k] * skl1_171[k];

        t_352[k] = f_15 * skk_136[k]
                   + f_3 * pc_z[k] * slk_280[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, pc_y, skk_174, skk_175, skk_176, sli0_219, \
                         sli0_220, sli0_221, sli1_219, sli1_220, sli1_221, slk_282, slk_283, \
                         slk_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_16 * skk_174[k]
                   + f_4 * sli0_219[k]
                   - f_5 * sli1_219[k]
                   + f_3 * pc_y[k] * slk_282[k];

        t_354[k] = f_16 * skk_175[k]
                   + f_6 * sli0_220[k]
                   - f_7 * sli1_220[k]
                   + f_3 * pc_y[k] * slk_283[k];

        t_355[k] = f_16 * skk_176[k]
                   + f_8 * sli0_221[k]
                   - f_9 * sli1_221[k]
                   + f_3 * pc_y[k] * slk_284[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pc_y, skk_177, skk_178, skk_179, sli0_222, \
                         sli0_223, sli1_222, sli1_223, slk_285, slk_286, \
                         slk_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_16 * skk_177[k]
                   + f_10 * sli0_222[k]
                   - f_11 * sli1_222[k]
                   + f_3 * pc_y[k] * slk_285[k];

        t_357[k] = f_16 * skk_178[k]
                   + f_12 * sli0_223[k]
                   - f_13 * sli1_223[k]
                   + f_3 * pc_y[k] * slk_286[k];

        t_358[k] = f_16 * skk_179[k]
                   + f_3 * pc_y[k] * slk_287[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pb_y, pc_y, pc_z, skl0_225, skk_143, \
                         skk_144, skk_180, skl1_225, sli0_223, sli1_223, slk_287, \
                         slk_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_15 * skk_143[k]
                   + f_1 * sli0_223[k]
                   - f_2 * sli1_223[k]
                   + f_3 * pc_z[k] * slk_287[k];

        t_360[k] = pb_y[k] * skl0_225[k]
                   - f_14 * pc_y[k] * skl1_225[k];

        t_361[k] = f_15 * skk_180[k]
                   + f_3 * pc_y[k] * slk_288[k];

        t_362[k] = f_16 * skk_144[k]
                   + f_3 * pc_z[k] * slk_288[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, pb_y, pc_y, skl0_228, skl0_230, skl0_231, \
                         skk_181, skk_182, skk_183, skl1_228, skl1_230, skl1_231, \
                         slk_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = pb_y[k] * skl0_228[k]
                   + f_16 * skk_181[k]
                   - f_14 * pc_y[k] * skl1_228[k];

        t_364[k] = f_15 * skk_182[k]
                   + f_3 * pc_y[k] * slk_290[k];

        t_365[k] = pb_y[k] * skl0_230[k]
                   - f_14 * pc_y[k] * skl1_230[k];

        t_366[k] = pb_y[k] * skl0_231[k]
                   + f_17 * skk_183[k]
                   - f_14 * pc_y[k] * skl1_231[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pb_y, pc_y, pc_z, skl0_234, skl0_235, \
                         skk_147, skk_185, skk_186, skl1_234, skl1_235, slk_291, \
                         slk_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_16 * skk_147[k]
                   + f_3 * pc_z[k] * slk_291[k];

        t_368[k] = f_15 * skk_185[k]
                   + f_3 * pc_y[k] * slk_293[k];

        t_369[k] = pb_y[k] * skl0_234[k]
                   - f_14 * pc_y[k] * skl1_234[k];

        t_370[k] = pb_y[k] * skl0_235[k]
                   + f_18 * skk_186[k]
                   - f_14 * pc_y[k] * skl1_235[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, t_374, pb_y, pc_y, pc_z, skl0_237, skl0_239, \
                         skk_150, skk_188, skk_189, skl1_237, skl1_239, slk_294, \
                         slk_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_16 * skk_150[k]
                   + f_3 * pc_z[k] * slk_294[k];

        t_372[k] = pb_y[k] * skl0_237[k]
                   + f_16 * skk_188[k]
                   - f_14 * pc_y[k] * skl1_237[k];

        t_373[k] = f_15 * skk_189[k]
                   + f_3 * pc_y[k] * slk_297[k];

        t_374[k] = pb_y[k] * skl0_239[k]
                   - f_14 * pc_y[k] * skl1_239[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pb_y, pc_y, pc_z, skl0_240, skl0_242, skk_154, \
                         skk_190, skk_192, skl1_240, skl1_242, \
                         slk_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = pb_y[k] * skl0_240[k]
                   + f_19 * skk_190[k]
                   - f_14 * pc_y[k] * skl1_240[k];

        t_376[k] = f_16 * skk_154[k]
                   + f_3 * pc_z[k] * slk_298[k];

        t_377[k] = pb_y[k] * skl0_242[k]
                   + f_17 * skk_192[k]
                   - f_14 * pc_y[k] * skl1_242[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_y, pc_y, skl0_243, skl0_245, skl0_246, \
                         skk_193, skk_194, skk_195, skl1_243, skl1_245, skl1_246, \
                         slk_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_y[k] * skl0_243[k]
                   + f_16 * skk_193[k]
                   - f_14 * pc_y[k] * skl1_243[k];

        t_379[k] = f_15 * skk_194[k]
                   + f_3 * pc_y[k] * slk_302[k];

        t_380[k] = pb_y[k] * skl0_245[k]
                   - f_14 * pc_y[k] * skl1_245[k];

        t_381[k] = pb_y[k] * skl0_246[k]
                   + f_20 * skk_195[k]
                   - f_14 * pc_y[k] * skl1_246[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_y, pc_y, pc_z, skl0_248, skl0_249, skk_159, \
                         skk_197, skk_198, skl1_248, skl1_249, \
                         slk_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_16 * skk_159[k]
                   + f_3 * pc_z[k] * slk_303[k];

        t_383[k] = pb_y[k] * skl0_248[k]
                   + f_18 * skk_197[k]
                   - f_14 * pc_y[k] * skl1_248[k];

        t_384[k] = pb_y[k] * skl0_249[k]
                   + f_17 * skk_198[k]
                   - f_14 * pc_y[k] * skl1_249[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pb_y, pc_x, pc_y, skl0_250, skl0_252, \
                         skk_199, skk_200, skk_316, skl1_250, skl1_252, slk_308, \
                         slk_316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pb_y[k] * skl0_250[k]
                   + f_16 * skk_199[k]
                   - f_14 * pc_y[k] * skl1_250[k];

        t_386[k] = f_15 * skk_200[k]
                   + f_3 * pc_y[k] * slk_308[k];

        t_387[k] = pb_y[k] * skl0_252[k]
                   - f_14 * pc_y[k] * skl1_252[k];

        t_388[k] = f_19 * skk_316[k]
                   + f_3 * pc_x[k] * slk_316[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, t_393, pc_x, skk_317, skk_318, skk_319, \
                         skk_320, skk_321, slk_317, slk_318, slk_319, slk_320, \
                         slk_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = f_19 * skk_317[k]
                   + f_3 * pc_x[k] * slk_317[k];

        t_390[k] = f_19 * skk_318[k]
                   + f_3 * pc_x[k] * slk_318[k];

        t_391[k] = f_19 * skk_319[k]
                   + f_3 * pc_x[k] * slk_319[k];

        t_392[k] = f_19 * skk_320[k]
                   + f_3 * pc_x[k] * slk_320[k];

        t_393[k] = f_19 * skk_321[k]
                   + f_3 * pc_x[k] * slk_321[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_x, pc_y, pc_z, skk_172, skk_208, \
                         skk_322, skk_323, sli0_245, sli1_245, slk_316, slk_322, \
                         slk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_19 * skk_322[k]
                   + f_3 * pc_x[k] * slk_322[k];

        t_395[k] = f_19 * skk_323[k]
                   + f_3 * pc_x[k] * slk_323[k];

        t_396[k] = f_15 * skk_208[k]
                   + f_1 * sli0_245[k]
                   - f_2 * sli1_245[k]
                   + f_3 * pc_y[k] * slk_316[k];

        t_397[k] = f_16 * skk_172[k]
                   + f_3 * pc_z[k] * slk_316[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pc_y, skk_210, skk_211, skk_212, sli0_247, \
                         sli0_248, sli0_249, sli1_247, sli1_248, sli1_249, slk_318, slk_319, \
                         slk_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_15 * skk_210[k]
                   + f_4 * sli0_247[k]
                   - f_5 * sli1_247[k]
                   + f_3 * pc_y[k] * slk_318[k];

        t_399[k] = f_15 * skk_211[k]
                   + f_6 * sli0_248[k]
                   - f_7 * sli1_248[k]
                   + f_3 * pc_y[k] * slk_319[k];

        t_400[k] = f_15 * skk_212[k]
                   + f_8 * sli0_249[k]
                   - f_9 * sli1_249[k]
                   + f_3 * pc_y[k] * slk_320[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pc_y, skk_213, skk_214, skk_215, sli0_250, \
                         sli0_251, sli1_250, sli1_251, slk_321, slk_322, \
                         slk_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_15 * skk_213[k]
                   + f_10 * sli0_250[k]
                   - f_11 * sli1_250[k]
                   + f_3 * pc_y[k] * slk_321[k];

        t_402[k] = f_15 * skk_214[k]
                   + f_12 * sli0_251[k]
                   - f_13 * sli1_251[k]
                   + f_3 * pc_y[k] * slk_322[k];

        t_403[k] = f_15 * skk_215[k]
                   + f_3 * pc_y[k] * slk_323[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pb_y, pc_x, pc_y, pc_z, skl0_269, \
                         skk_180, skk_324, skl1_269, sli0_252, sli1_252, \
                         slk_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pb_y[k] * skl0_269[k]
                   - f_14 * pc_y[k] * skl1_269[k];

        t_405[k] = f_19 * skk_324[k]
                   + f_1 * sli0_252[k]
                   - f_2 * sli1_252[k]
                   + f_3 * pc_x[k] * slk_324[k];

        t_406[k] = f_3 * pc_y[k] * slk_324[k];

        t_407[k] = f_17 * skk_180[k]
                   + f_3 * pc_z[k] * slk_324[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, skk_327, skk_329, sli0_255, \
                         sli0_257, sli1_255, sli1_257, slk_326, slk_327, \
                         slk_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_19 * skk_327[k]
                   + f_4 * sli0_255[k]
                   - f_5 * sli1_255[k]
                   + f_3 * pc_x[k] * slk_327[k];

        t_409[k] = f_3 * pc_y[k] * slk_326[k];

        t_410[k] = f_19 * skk_329[k]
                   + f_4 * sli0_257[k]
                   - f_5 * sli1_257[k]
                   + f_3 * pc_x[k] * slk_329[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, pc_x, pc_y, pc_z, skk_183, skk_330, sli0_258, \
                         sli1_258, slk_327, slk_329, slk_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_19 * skk_330[k]
                   + f_6 * sli0_258[k]
                   - f_7 * sli1_258[k]
                   + f_3 * pc_x[k] * slk_330[k];

        t_412[k] = f_17 * skk_183[k]
                   + f_3 * pc_z[k] * slk_327[k];

        t_413[k] = f_3 * pc_y[k] * slk_329[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_z, skk_186, skk_333, skk_334, sli0_261, \
                         sli0_262, sli1_261, sli1_262, slk_330, slk_333, \
                         slk_334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_19 * skk_333[k]
                   + f_6 * sli0_261[k]
                   - f_7 * sli1_261[k]
                   + f_3 * pc_x[k] * slk_333[k];

        t_415[k] = f_19 * skk_334[k]
                   + f_8 * sli0_262[k]
                   - f_9 * sli1_262[k]
                   + f_3 * pc_x[k] * slk_334[k];

        t_416[k] = f_17 * skk_186[k]
                   + f_3 * pc_z[k] * slk_330[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pc_x, pc_y, skk_336, skk_338, sli0_264, \
                         sli0_266, sli1_264, sli1_266, slk_333, slk_336, \
                         slk_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_19 * skk_336[k]
                   + f_8 * sli0_264[k]
                   - f_9 * sli1_264[k]
                   + f_3 * pc_x[k] * slk_336[k];

        t_418[k] = f_3 * pc_y[k] * slk_333[k];

        t_419[k] = f_19 * skk_338[k]
                   + f_8 * sli0_266[k]
                   - f_9 * sli1_266[k]
                   + f_3 * pc_x[k] * slk_338[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pc_x, pc_z, skk_190, skk_339, skk_341, sli0_267, \
                         sli0_269, sli1_267, sli1_269, slk_334, slk_339, \
                         slk_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = f_19 * skk_339[k]
                   + f_10 * sli0_267[k]
                   - f_11 * sli1_267[k]
                   + f_3 * pc_x[k] * slk_339[k];

        t_421[k] = f_17 * skk_190[k]
                   + f_3 * pc_z[k] * slk_334[k];

        t_422[k] = f_19 * skk_341[k]
                   + f_10 * sli0_269[k]
                   - f_11 * sli1_269[k]
                   + f_3 * pc_x[k] * slk_341[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pc_x, pc_y, skk_342, skk_344, sli0_270, \
                         sli0_272, sli1_270, sli1_272, slk_338, slk_342, \
                         slk_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_19 * skk_342[k]
                   + f_10 * sli0_270[k]
                   - f_11 * sli1_270[k]
                   + f_3 * pc_x[k] * slk_342[k];

        t_424[k] = f_3 * pc_y[k] * slk_338[k];

        t_425[k] = f_19 * skk_344[k]
                   + f_10 * sli0_272[k]
                   - f_11 * sli1_272[k]
                   + f_3 * pc_x[k] * slk_344[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pc_x, pc_z, skk_195, skk_345, skk_347, sli0_273, \
                         sli0_275, sli1_273, sli1_275, slk_339, slk_345, \
                         slk_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_19 * skk_345[k]
                   + f_12 * sli0_273[k]
                   - f_13 * sli1_273[k]
                   + f_3 * pc_x[k] * slk_345[k];

        t_427[k] = f_17 * skk_195[k]
                   + f_3 * pc_z[k] * slk_339[k];

        t_428[k] = f_19 * skk_347[k]
                   + f_12 * sli0_275[k]
                   - f_13 * sli1_275[k]
                   + f_3 * pc_x[k] * slk_347[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pc_x, pc_y, skk_348, skk_349, sli0_276, \
                         sli0_277, sli1_276, sli1_277, slk_344, slk_348, \
                         slk_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_19 * skk_348[k]
                   + f_12 * sli0_276[k]
                   - f_13 * sli1_276[k]
                   + f_3 * pc_x[k] * slk_348[k];

        t_430[k] = f_19 * skk_349[k]
                   + f_12 * sli0_277[k]
                   - f_13 * sli1_277[k]
                   + f_3 * pc_x[k] * slk_349[k];

        t_431[k] = f_3 * pc_y[k] * slk_344[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pc_x, skk_351, skk_352, skk_353, skk_354, \
                         sli0_279, sli1_279, slk_351, slk_352, slk_353, \
                         slk_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_19 * skk_351[k]
                   + f_12 * sli0_279[k]
                   - f_13 * sli1_279[k]
                   + f_3 * pc_x[k] * slk_351[k];

        t_433[k] = f_19 * skk_352[k]
                   + f_3 * pc_x[k] * slk_352[k];

        t_434[k] = f_19 * skk_353[k]
                   + f_3 * pc_x[k] * slk_353[k];

        t_435[k] = f_19 * skk_354[k]
                   + f_3 * pc_x[k] * slk_354[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, skk_355, skk_356, skk_357, \
                         skk_358, skk_359, slk_355, slk_356, slk_357, slk_358, \
                         slk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_19 * skk_355[k]
                   + f_3 * pc_x[k] * slk_355[k];

        t_437[k] = f_19 * skk_356[k]
                   + f_3 * pc_x[k] * slk_356[k];

        t_438[k] = f_19 * skk_357[k]
                   + f_3 * pc_x[k] * slk_357[k];

        t_439[k] = f_19 * skk_358[k]
                   + f_3 * pc_x[k] * slk_358[k];

        t_440[k] = f_19 * skk_359[k]
                   + f_3 * pc_x[k] * slk_359[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_y, pc_z, skk_208, sli0_273, sli0_275, \
                         sli0_276, sli1_273, sli1_275, sli1_276, slk_352, slk_354, \
                         slk_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_1 * sli0_273[k]
                   - f_2 * sli1_273[k]
                   + f_3 * pc_y[k] * slk_352[k];

        t_442[k] = f_17 * skk_208[k]
                   + f_3 * pc_z[k] * slk_352[k];

        t_443[k] = f_4 * sli0_275[k]
                   - f_5 * sli1_275[k]
                   + f_3 * pc_y[k] * slk_354[k];

        t_444[k] = f_6 * sli0_276[k]
                   - f_7 * sli1_276[k]
                   + f_3 * pc_y[k] * slk_355[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, sli0_277, sli0_278, sli0_279, \
                         sli1_277, sli1_278, sli1_279, slk_356, slk_357, slk_358, \
                         slk_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_8 * sli0_277[k]
                   - f_9 * sli1_277[k]
                   + f_3 * pc_y[k] * slk_356[k];

        t_446[k] = f_10 * sli0_278[k]
                   - f_11 * sli1_278[k]
                   + f_3 * pc_y[k] * slk_357[k];

        t_447[k] = f_12 * sli0_279[k]
                   - f_13 * sli1_279[k]
                   + f_3 * pc_y[k] * slk_358[k];

        t_448[k] = f_3 * pc_y[k] * slk_359[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, skk_215, skk_216, \
                         skk_360, sli0_279, sli0_280, sli1_279, sli1_280, slk_359, \
                         slk_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_17 * skk_215[k]
                   + f_1 * sli0_279[k]
                   - f_2 * sli1_279[k]
                   + f_3 * pc_z[k] * slk_359[k];

        t_450[k] = f_18 * skk_360[k]
                   + f_1 * sli0_280[k]
                   - f_2 * sli1_280[k]
                   + f_3 * pc_x[k] * slk_360[k];

        t_451[k] = f_18 * skk_216[k]
                   + f_3 * pc_y[k] * slk_360[k];

        t_452[k] = f_3 * pc_z[k] * slk_360[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pc_x, pc_y, skk_218, skk_363, skk_365, sli0_283, \
                         sli0_285, sli1_283, sli1_285, slk_362, slk_363, \
                         slk_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_18 * skk_363[k]
                   + f_4 * sli0_283[k]
                   - f_5 * sli1_283[k]
                   + f_3 * pc_x[k] * slk_363[k];

        t_454[k] = f_18 * skk_218[k]
                   + f_3 * pc_y[k] * slk_362[k];

        t_455[k] = f_18 * skk_365[k]
                   + f_4 * sli0_285[k]
                   - f_5 * sli1_285[k]
                   + f_3 * pc_x[k] * slk_365[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_270 = buffer.data(skl0 + 270);
    const auto *skl0_273 = buffer.data(skl0 + 273);
    const auto *skl0_276 = buffer.data(skl0 + 276);
    const auto *skl0_280 = buffer.data(skl0 + 280);
    const auto *skl0_282 = buffer.data(skl0 + 282);
    const auto *skl0_285 = buffer.data(skl0 + 285);
    const auto *skl0_287 = buffer.data(skl0 + 287);
    const auto *skl0_288 = buffer.data(skl0 + 288);
    const auto *skl0_291 = buffer.data(skl0 + 291);
    const auto *skl0_293 = buffer.data(skl0 + 293);
    const auto *skl0_294 = buffer.data(skl0 + 294);
    const auto *skl0_295 = buffer.data(skl0 + 295);
    const auto *skl0_306 = buffer.data(skl0 + 306);

    const auto *skk_216 = buffer.data(skk + 216);
    const auto *skk_219 = buffer.data(skk + 219);
    const auto *skk_221 = buffer.data(skk + 221);
    const auto *skk_222 = buffer.data(skk + 222);
    const auto *skk_223 = buffer.data(skk + 223);
    const auto *skk_225 = buffer.data(skk + 225);
    const auto *skk_226 = buffer.data(skk + 226);
    const auto *skk_227 = buffer.data(skk + 227);
    const auto *skk_228 = buffer.data(skk + 228);
    const auto *skk_230 = buffer.data(skk + 230);
    const auto *skk_231 = buffer.data(skk + 231);
    const auto *skk_232 = buffer.data(skk + 232);
    const auto *skk_233 = buffer.data(skk + 233);
    const auto *skk_234 = buffer.data(skk + 234);
    const auto *skk_236 = buffer.data(skk + 236);
    const auto *skk_244 = buffer.data(skk + 244);
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
    const auto *skk_282 = buffer.data(skk + 282);
    const auto *skk_283 = buffer.data(skk + 283);
    const auto *skk_284 = buffer.data(skk + 284);
    const auto *skk_285 = buffer.data(skk + 285);
    const auto *skk_286 = buffer.data(skk + 286);
    const auto *skk_287 = buffer.data(skk + 287);
    const auto *skk_288 = buffer.data(skk + 288);
    const auto *skk_290 = buffer.data(skk + 290);
    const auto *skk_293 = buffer.data(skk + 293);
    const auto *skk_297 = buffer.data(skk + 297);
    const auto *skk_302 = buffer.data(skk + 302);
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
    const auto *skk_401 = buffer.data(skk + 401);
    const auto *skk_405 = buffer.data(skk + 405);
    const auto *skk_410 = buffer.data(skk + 410);
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

    const auto *skl1_270 = buffer.data(skl1 + 270);
    const auto *skl1_273 = buffer.data(skl1 + 273);
    const auto *skl1_276 = buffer.data(skl1 + 276);
    const auto *skl1_280 = buffer.data(skl1 + 280);
    const auto *skl1_282 = buffer.data(skl1 + 282);
    const auto *skl1_285 = buffer.data(skl1 + 285);
    const auto *skl1_287 = buffer.data(skl1 + 287);
    const auto *skl1_288 = buffer.data(skl1 + 288);
    const auto *skl1_291 = buffer.data(skl1 + 291);
    const auto *skl1_293 = buffer.data(skl1 + 293);
    const auto *skl1_294 = buffer.data(skl1 + 294);
    const auto *skl1_295 = buffer.data(skl1 + 295);
    const auto *skl1_306 = buffer.data(skl1 + 306);

    const auto *sli0_286 = buffer.data(sli0 + 286);
    const auto *sli0_289 = buffer.data(sli0 + 289);
    const auto *sli0_290 = buffer.data(sli0 + 290);
    const auto *sli0_292 = buffer.data(sli0 + 292);
    const auto *sli0_294 = buffer.data(sli0 + 294);
    const auto *sli0_295 = buffer.data(sli0 + 295);
    const auto *sli0_297 = buffer.data(sli0 + 297);
    const auto *sli0_298 = buffer.data(sli0 + 298);
    const auto *sli0_300 = buffer.data(sli0 + 300);
    const auto *sli0_301 = buffer.data(sli0 + 301);
    const auto *sli0_303 = buffer.data(sli0 + 303);
    const auto *sli0_304 = buffer.data(sli0 + 304);
    const auto *sli0_305 = buffer.data(sli0 + 305);
    const auto *sli0_306 = buffer.data(sli0 + 306);
    const auto *sli0_307 = buffer.data(sli0 + 307);
    const auto *sli0_313 = buffer.data(sli0 + 313);
    const auto *sli0_317 = buffer.data(sli0 + 317);
    const auto *sli0_322 = buffer.data(sli0 + 322);
    const auto *sli0_328 = buffer.data(sli0 + 328);
    const auto *sli0_331 = buffer.data(sli0 + 331);
    const auto *sli0_332 = buffer.data(sli0 + 332);
    const auto *sli0_333 = buffer.data(sli0 + 333);
    const auto *sli0_334 = buffer.data(sli0 + 334);
    const auto *sli0_335 = buffer.data(sli0 + 335);
    const auto *sli0_336 = buffer.data(sli0 + 336);
    const auto *sli0_339 = buffer.data(sli0 + 339);
    const auto *sli0_341 = buffer.data(sli0 + 341);
    const auto *sli0_342 = buffer.data(sli0 + 342);
    const auto *sli0_345 = buffer.data(sli0 + 345);
    const auto *sli0_346 = buffer.data(sli0 + 346);
    const auto *sli0_348 = buffer.data(sli0 + 348);
    const auto *sli0_350 = buffer.data(sli0 + 350);
    const auto *sli0_351 = buffer.data(sli0 + 351);
    const auto *sli0_353 = buffer.data(sli0 + 353);
    const auto *sli0_354 = buffer.data(sli0 + 354);
    const auto *sli0_356 = buffer.data(sli0 + 356);
    const auto *sli0_357 = buffer.data(sli0 + 357);

    const auto *sli1_286 = buffer.data(sli1 + 286);
    const auto *sli1_289 = buffer.data(sli1 + 289);
    const auto *sli1_290 = buffer.data(sli1 + 290);
    const auto *sli1_292 = buffer.data(sli1 + 292);
    const auto *sli1_294 = buffer.data(sli1 + 294);
    const auto *sli1_295 = buffer.data(sli1 + 295);
    const auto *sli1_297 = buffer.data(sli1 + 297);
    const auto *sli1_298 = buffer.data(sli1 + 298);
    const auto *sli1_300 = buffer.data(sli1 + 300);
    const auto *sli1_301 = buffer.data(sli1 + 301);
    const auto *sli1_303 = buffer.data(sli1 + 303);
    const auto *sli1_304 = buffer.data(sli1 + 304);
    const auto *sli1_305 = buffer.data(sli1 + 305);
    const auto *sli1_306 = buffer.data(sli1 + 306);
    const auto *sli1_307 = buffer.data(sli1 + 307);
    const auto *sli1_313 = buffer.data(sli1 + 313);
    const auto *sli1_317 = buffer.data(sli1 + 317);
    const auto *sli1_322 = buffer.data(sli1 + 322);
    const auto *sli1_328 = buffer.data(sli1 + 328);
    const auto *sli1_331 = buffer.data(sli1 + 331);
    const auto *sli1_332 = buffer.data(sli1 + 332);
    const auto *sli1_333 = buffer.data(sli1 + 333);
    const auto *sli1_334 = buffer.data(sli1 + 334);
    const auto *sli1_335 = buffer.data(sli1 + 335);
    const auto *sli1_336 = buffer.data(sli1 + 336);
    const auto *sli1_339 = buffer.data(sli1 + 339);
    const auto *sli1_341 = buffer.data(sli1 + 341);
    const auto *sli1_342 = buffer.data(sli1 + 342);
    const auto *sli1_345 = buffer.data(sli1 + 345);
    const auto *sli1_346 = buffer.data(sli1 + 346);
    const auto *sli1_348 = buffer.data(sli1 + 348);
    const auto *sli1_350 = buffer.data(sli1 + 350);
    const auto *sli1_351 = buffer.data(sli1 + 351);
    const auto *sli1_353 = buffer.data(sli1 + 353);
    const auto *sli1_354 = buffer.data(sli1 + 354);
    const auto *sli1_356 = buffer.data(sli1 + 356);
    const auto *sli1_357 = buffer.data(sli1 + 357);

    const auto *slk_363 = buffer.data(slk + 363);
    const auto *slk_365 = buffer.data(slk + 365);
    const auto *slk_366 = buffer.data(slk + 366);
    const auto *slk_369 = buffer.data(slk + 369);
    const auto *slk_370 = buffer.data(slk + 370);
    const auto *slk_372 = buffer.data(slk + 372);
    const auto *slk_374 = buffer.data(slk + 374);
    const auto *slk_375 = buffer.data(slk + 375);
    const auto *slk_377 = buffer.data(slk + 377);
    const auto *slk_378 = buffer.data(slk + 378);
    const auto *slk_380 = buffer.data(slk + 380);
    const auto *slk_381 = buffer.data(slk + 381);
    const auto *slk_383 = buffer.data(slk + 383);
    const auto *slk_384 = buffer.data(slk + 384);
    const auto *slk_385 = buffer.data(slk + 385);
    const auto *slk_387 = buffer.data(slk + 387);
    const auto *slk_388 = buffer.data(slk + 388);
    const auto *slk_389 = buffer.data(slk + 389);
    const auto *slk_390 = buffer.data(slk + 390);
    const auto *slk_391 = buffer.data(slk + 391);
    const auto *slk_392 = buffer.data(slk + 392);
    const auto *slk_393 = buffer.data(slk + 393);
    const auto *slk_394 = buffer.data(slk + 394);
    const auto *slk_395 = buffer.data(slk + 395);
    const auto *slk_396 = buffer.data(slk + 396);
    const auto *slk_398 = buffer.data(slk + 398);
    const auto *slk_399 = buffer.data(slk + 399);
    const auto *slk_401 = buffer.data(slk + 401);
    const auto *slk_402 = buffer.data(slk + 402);
    const auto *slk_405 = buffer.data(slk + 405);
    const auto *slk_406 = buffer.data(slk + 406);
    const auto *slk_410 = buffer.data(slk + 410);
    const auto *slk_411 = buffer.data(slk + 411);
    const auto *slk_416 = buffer.data(slk + 416);
    const auto *slk_423 = buffer.data(slk + 423);
    const auto *slk_424 = buffer.data(slk + 424);
    const auto *slk_425 = buffer.data(slk + 425);
    const auto *slk_426 = buffer.data(slk + 426);
    const auto *slk_427 = buffer.data(slk + 427);
    const auto *slk_428 = buffer.data(slk + 428);
    const auto *slk_429 = buffer.data(slk + 429);
    const auto *slk_430 = buffer.data(slk + 430);
    const auto *slk_431 = buffer.data(slk + 431);
    const auto *slk_432 = buffer.data(slk + 432);
    const auto *slk_434 = buffer.data(slk + 434);
    const auto *slk_435 = buffer.data(slk + 435);
    const auto *slk_437 = buffer.data(slk + 437);
    const auto *slk_438 = buffer.data(slk + 438);
    const auto *slk_441 = buffer.data(slk + 441);
    const auto *slk_442 = buffer.data(slk + 442);
    const auto *slk_444 = buffer.data(slk + 444);
    const auto *slk_446 = buffer.data(slk + 446);
    const auto *slk_447 = buffer.data(slk + 447);
    const auto *slk_449 = buffer.data(slk + 449);
    const auto *slk_450 = buffer.data(slk + 450);
    const auto *slk_452 = buffer.data(slk + 452);
    const auto *slk_453 = buffer.data(slk + 453);

#pragma omp simd aligned(t_456, t_457, t_458, pc_x, pc_y, pc_z, skk_221, skk_366, sli0_286, \
                         sli1_286, slk_363, slk_365, slk_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_18 * skk_366[k]
                   + f_6 * sli0_286[k]
                   - f_7 * sli1_286[k]
                   + f_3 * pc_x[k] * slk_366[k];

        t_457[k] = f_3 * pc_z[k] * slk_363[k];

        t_458[k] = f_18 * skk_221[k]
                   + f_3 * pc_y[k] * slk_365[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_z, skk_369, skk_370, sli0_289, \
                         sli0_290, sli1_289, sli1_290, slk_366, slk_369, \
                         slk_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_18 * skk_369[k]
                   + f_6 * sli0_289[k]
                   - f_7 * sli1_289[k]
                   + f_3 * pc_x[k] * slk_369[k];

        t_460[k] = f_18 * skk_370[k]
                   + f_8 * sli0_290[k]
                   - f_9 * sli1_290[k]
                   + f_3 * pc_x[k] * slk_370[k];

        t_461[k] = f_3 * pc_z[k] * slk_366[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_x, pc_y, skk_225, skk_372, skk_374, sli0_292, \
                         sli0_294, sli1_292, sli1_294, slk_369, slk_372, \
                         slk_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_18 * skk_372[k]
                   + f_8 * sli0_292[k]
                   - f_9 * sli1_292[k]
                   + f_3 * pc_x[k] * slk_372[k];

        t_463[k] = f_18 * skk_225[k]
                   + f_3 * pc_y[k] * slk_369[k];

        t_464[k] = f_18 * skk_374[k]
                   + f_8 * sli0_294[k]
                   - f_9 * sli1_294[k]
                   + f_3 * pc_x[k] * slk_374[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_z, skk_375, skk_377, sli0_295, \
                         sli0_297, sli1_295, sli1_297, slk_370, slk_375, \
                         slk_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_18 * skk_375[k]
                   + f_10 * sli0_295[k]
                   - f_11 * sli1_295[k]
                   + f_3 * pc_x[k] * slk_375[k];

        t_466[k] = f_3 * pc_z[k] * slk_370[k];

        t_467[k] = f_18 * skk_377[k]
                   + f_10 * sli0_297[k]
                   - f_11 * sli1_297[k]
                   + f_3 * pc_x[k] * slk_377[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, skk_230, skk_378, skk_380, sli0_298, \
                         sli0_300, sli1_298, sli1_300, slk_374, slk_378, \
                         slk_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_18 * skk_378[k]
                   + f_10 * sli0_298[k]
                   - f_11 * sli1_298[k]
                   + f_3 * pc_x[k] * slk_378[k];

        t_469[k] = f_18 * skk_230[k]
                   + f_3 * pc_y[k] * slk_374[k];

        t_470[k] = f_18 * skk_380[k]
                   + f_10 * sli0_300[k]
                   - f_11 * sli1_300[k]
                   + f_3 * pc_x[k] * slk_380[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pc_x, pc_z, skk_381, skk_383, sli0_301, \
                         sli0_303, sli1_301, sli1_303, slk_375, slk_381, \
                         slk_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_18 * skk_381[k]
                   + f_12 * sli0_301[k]
                   - f_13 * sli1_301[k]
                   + f_3 * pc_x[k] * slk_381[k];

        t_472[k] = f_3 * pc_z[k] * slk_375[k];

        t_473[k] = f_18 * skk_383[k]
                   + f_12 * sli0_303[k]
                   - f_13 * sli1_303[k]
                   + f_3 * pc_x[k] * slk_383[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, skk_236, skk_384, skk_385, sli0_304, \
                         sli0_305, sli1_304, sli1_305, slk_380, slk_384, \
                         slk_385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_18 * skk_384[k]
                   + f_12 * sli0_304[k]
                   - f_13 * sli1_304[k]
                   + f_3 * pc_x[k] * slk_384[k];

        t_475[k] = f_18 * skk_385[k]
                   + f_12 * sli0_305[k]
                   - f_13 * sli1_305[k]
                   + f_3 * pc_x[k] * slk_385[k];

        t_476[k] = f_18 * skk_236[k]
                   + f_3 * pc_y[k] * slk_380[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pc_x, skk_387, skk_388, skk_389, skk_390, \
                         sli0_307, sli1_307, slk_387, slk_388, slk_389, \
                         slk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_18 * skk_387[k]
                   + f_12 * sli0_307[k]
                   - f_13 * sli1_307[k]
                   + f_3 * pc_x[k] * slk_387[k];

        t_478[k] = f_18 * skk_388[k]
                   + f_3 * pc_x[k] * slk_388[k];

        t_479[k] = f_18 * skk_389[k]
                   + f_3 * pc_x[k] * slk_389[k];

        t_480[k] = f_18 * skk_390[k]
                   + f_3 * pc_x[k] * slk_390[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, t_485, pc_x, skk_391, skk_392, skk_393, \
                         skk_394, skk_395, slk_391, slk_392, slk_393, slk_394, \
                         slk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_18 * skk_391[k]
                   + f_3 * pc_x[k] * slk_391[k];

        t_482[k] = f_18 * skk_392[k]
                   + f_3 * pc_x[k] * slk_392[k];

        t_483[k] = f_18 * skk_393[k]
                   + f_3 * pc_x[k] * slk_393[k];

        t_484[k] = f_18 * skk_394[k]
                   + f_3 * pc_x[k] * slk_394[k];

        t_485[k] = f_18 * skk_395[k]
                   + f_3 * pc_x[k] * slk_395[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, pc_y, pc_z, skk_244, skk_246, sli0_301, \
                         sli0_303, sli1_301, sli1_303, slk_388, \
                         slk_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_18 * skk_244[k]
                   + f_1 * sli0_301[k]
                   - f_2 * sli1_301[k]
                   + f_3 * pc_y[k] * slk_388[k];

        t_487[k] = f_3 * pc_z[k] * slk_388[k];

        t_488[k] = f_18 * skk_246[k]
                   + f_4 * sli0_303[k]
                   - f_5 * sli1_303[k]
                   + f_3 * pc_y[k] * slk_390[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_y, skk_247, skk_248, skk_249, sli0_304, \
                         sli0_305, sli0_306, sli1_304, sli1_305, sli1_306, slk_391, slk_392, \
                         slk_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_18 * skk_247[k]
                   + f_6 * sli0_304[k]
                   - f_7 * sli1_304[k]
                   + f_3 * pc_y[k] * slk_391[k];

        t_490[k] = f_18 * skk_248[k]
                   + f_8 * sli0_305[k]
                   - f_9 * sli1_305[k]
                   + f_3 * pc_y[k] * slk_392[k];

        t_491[k] = f_18 * skk_249[k]
                   + f_10 * sli0_306[k]
                   - f_11 * sli1_306[k]
                   + f_3 * pc_y[k] * slk_393[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pb_z, pc_y, pc_z, skl0_270, skk_250, \
                         skk_251, skl1_270, sli0_307, sli1_307, slk_394, \
                         slk_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_18 * skk_250[k]
                   + f_12 * sli0_307[k]
                   - f_13 * sli1_307[k]
                   + f_3 * pc_y[k] * slk_394[k];

        t_493[k] = f_18 * skk_251[k]
                   + f_3 * pc_y[k] * slk_395[k];

        t_494[k] = f_1 * sli0_307[k]
                   - f_2 * sli1_307[k]
                   + f_3 * pc_z[k] * slk_395[k];

        t_495[k] = pb_z[k] * skl0_270[k]
                   - f_14 * pc_z[k] * skl1_270[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pb_z, pc_y, pc_z, skl0_273, skk_216, \
                         skk_252, skk_254, skl1_273, slk_396, slk_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_17 * skk_252[k]
                   + f_3 * pc_y[k] * slk_396[k];

        t_497[k] = f_15 * skk_216[k]
                   + f_3 * pc_z[k] * slk_396[k];

        t_498[k] = pb_z[k] * skl0_273[k]
                   - f_14 * pc_z[k] * skl1_273[k];

        t_499[k] = f_17 * skk_254[k]
                   + f_3 * pc_y[k] * slk_398[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, pb_z, pc_x, pc_z, skl0_276, skk_219, skk_401, \
                         skl1_276, sli0_313, sli1_313, slk_399, \
                         slk_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_18 * skk_401[k]
                   + f_4 * sli0_313[k]
                   - f_5 * sli1_313[k]
                   + f_3 * pc_x[k] * slk_401[k];

        t_501[k] = pb_z[k] * skl0_276[k]
                   - f_14 * pc_z[k] * skl1_276[k];

        t_502[k] = f_15 * skk_219[k]
                   + f_3 * pc_z[k] * slk_399[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, pb_z, pc_x, pc_y, pc_z, skl0_280, skk_257, \
                         skk_405, skl1_280, sli0_317, sli1_317, slk_401, \
                         slk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = f_17 * skk_257[k]
                   + f_3 * pc_y[k] * slk_401[k];

        t_504[k] = f_18 * skk_405[k]
                   + f_6 * sli0_317[k]
                   - f_7 * sli1_317[k]
                   + f_3 * pc_x[k] * slk_405[k];

        t_505[k] = pb_z[k] * skl0_280[k]
                   - f_14 * pc_z[k] * skl1_280[k];
    }

#pragma omp simd aligned(t_506, t_507, t_508, pb_z, pc_y, pc_z, skl0_282, skk_222, skk_223, \
                         skk_261, skl1_282, slk_402, slk_405 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_506[k] = f_15 * skk_222[k]
                   + f_3 * pc_z[k] * slk_402[k];

        t_507[k] = pb_z[k] * skl0_282[k]
                   + f_16 * skk_223[k]
                   - f_14 * pc_z[k] * skl1_282[k];

        t_508[k] = f_17 * skk_261[k]
                   + f_3 * pc_y[k] * slk_405[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_z, pc_x, pc_z, skl0_285, skk_226, skk_410, \
                         skl1_285, sli0_322, sli1_322, slk_406, \
                         slk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = f_18 * skk_410[k]
                   + f_8 * sli0_322[k]
                   - f_9 * sli1_322[k]
                   + f_3 * pc_x[k] * slk_410[k];

        t_510[k] = pb_z[k] * skl0_285[k]
                   - f_14 * pc_z[k] * skl1_285[k];

        t_511[k] = f_15 * skk_226[k]
                   + f_3 * pc_z[k] * slk_406[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_z, pc_y, pc_z, skl0_287, skl0_288, skk_227, \
                         skk_228, skk_266, skl1_287, skl1_288, \
                         slk_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_z[k] * skl0_287[k]
                   + f_16 * skk_227[k]
                   - f_14 * pc_z[k] * skl1_287[k];

        t_513[k] = pb_z[k] * skl0_288[k]
                   + f_17 * skk_228[k]
                   - f_14 * pc_z[k] * skl1_288[k];

        t_514[k] = f_17 * skk_266[k]
                   + f_3 * pc_y[k] * slk_410[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_z, pc_x, pc_z, skl0_291, skk_231, skk_416, \
                         skl1_291, sli0_328, sli1_328, slk_411, \
                         slk_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_18 * skk_416[k]
                   + f_10 * sli0_328[k]
                   - f_11 * sli1_328[k]
                   + f_3 * pc_x[k] * slk_416[k];

        t_516[k] = pb_z[k] * skl0_291[k]
                   - f_14 * pc_z[k] * skl1_291[k];

        t_517[k] = f_15 * skk_231[k]
                   + f_3 * pc_z[k] * slk_411[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, pb_z, pc_z, skl0_293, skl0_294, skl0_295, \
                         skk_232, skk_233, skk_234, skl1_293, skl1_294, \
                         skl1_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_z[k] * skl0_293[k]
                   + f_16 * skk_232[k]
                   - f_14 * pc_z[k] * skl1_293[k];

        t_519[k] = pb_z[k] * skl0_294[k]
                   + f_17 * skk_233[k]
                   - f_14 * pc_z[k] * skl1_294[k];

        t_520[k] = pb_z[k] * skl0_295[k]
                   + f_18 * skk_234[k]
                   - f_14 * pc_z[k] * skl1_295[k];
    }

#pragma omp simd aligned(t_521, t_522, t_523, t_524, pc_x, pc_y, skk_272, skk_423, skk_424, \
                         skk_425, sli0_335, sli1_335, slk_416, slk_423, slk_424, \
                         slk_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_521[k] = f_17 * skk_272[k]
                   + f_3 * pc_y[k] * slk_416[k];

        t_522[k] = f_18 * skk_423[k]
                   + f_12 * sli0_335[k]
                   - f_13 * sli1_335[k]
                   + f_3 * pc_x[k] * slk_423[k];

        t_523[k] = f_18 * skk_424[k]
                   + f_3 * pc_x[k] * slk_424[k];

        t_524[k] = f_18 * skk_425[k]
                   + f_3 * pc_x[k] * slk_425[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, t_528, t_529, pc_x, skk_426, skk_427, skk_428, \
                         skk_429, skk_430, slk_426, slk_427, slk_428, slk_429, \
                         slk_430 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_18 * skk_426[k]
                   + f_3 * pc_x[k] * slk_426[k];

        t_526[k] = f_18 * skk_427[k]
                   + f_3 * pc_x[k] * slk_427[k];

        t_527[k] = f_18 * skk_428[k]
                   + f_3 * pc_x[k] * slk_428[k];

        t_528[k] = f_18 * skk_429[k]
                   + f_3 * pc_x[k] * slk_429[k];

        t_529[k] = f_18 * skk_430[k]
                   + f_3 * pc_x[k] * slk_430[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pb_z, pc_x, pc_z, skl0_306, skk_244, skk_431, \
                         skl1_306, slk_424, slk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_18 * skk_431[k]
                   + f_3 * pc_x[k] * slk_431[k];

        t_531[k] = pb_z[k] * skl0_306[k]
                   - f_14 * pc_z[k] * skl1_306[k];

        t_532[k] = f_15 * skk_244[k]
                   + f_3 * pc_z[k] * slk_424[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, pc_y, skk_282, skk_283, skk_284, sli0_331, \
                         sli0_332, sli0_333, sli1_331, sli1_332, sli1_333, slk_426, slk_427, \
                         slk_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * skk_282[k]
                   + f_4 * sli0_331[k]
                   - f_5 * sli1_331[k]
                   + f_3 * pc_y[k] * slk_426[k];

        t_534[k] = f_17 * skk_283[k]
                   + f_6 * sli0_332[k]
                   - f_7 * sli1_332[k]
                   + f_3 * pc_y[k] * slk_427[k];

        t_535[k] = f_17 * skk_284[k]
                   + f_8 * sli0_333[k]
                   - f_9 * sli1_333[k]
                   + f_3 * pc_y[k] * slk_428[k];
    }

#pragma omp simd aligned(t_536, t_537, t_538, pc_y, skk_285, skk_286, skk_287, sli0_334, \
                         sli0_335, sli1_334, sli1_335, slk_429, slk_430, \
                         slk_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_536[k] = f_17 * skk_285[k]
                   + f_10 * sli0_334[k]
                   - f_11 * sli1_334[k]
                   + f_3 * pc_y[k] * slk_429[k];

        t_537[k] = f_17 * skk_286[k]
                   + f_12 * sli0_335[k]
                   - f_13 * sli1_335[k]
                   + f_3 * pc_y[k] * slk_430[k];

        t_538[k] = f_17 * skk_287[k]
                   + f_3 * pc_y[k] * slk_431[k];
    }

#pragma omp simd aligned(t_539, t_540, t_541, pc_x, pc_y, pc_z, skk_251, skk_288, skk_432, \
                         sli0_335, sli0_336, sli1_335, sli1_336, slk_431, \
                         slk_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_539[k] = f_15 * skk_251[k]
                   + f_1 * sli0_335[k]
                   - f_2 * sli1_335[k]
                   + f_3 * pc_z[k] * slk_431[k];

        t_540[k] = f_18 * skk_432[k]
                   + f_1 * sli0_336[k]
                   - f_2 * sli1_336[k]
                   + f_3 * pc_x[k] * slk_432[k];

        t_541[k] = f_16 * skk_288[k]
                   + f_3 * pc_y[k] * slk_432[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, pc_x, pc_y, pc_z, skk_252, skk_290, skk_435, \
                         sli0_339, sli1_339, slk_432, slk_434, \
                         slk_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_16 * skk_252[k]
                   + f_3 * pc_z[k] * slk_432[k];

        t_543[k] = f_18 * skk_435[k]
                   + f_4 * sli0_339[k]
                   - f_5 * sli1_339[k]
                   + f_3 * pc_x[k] * slk_435[k];

        t_544[k] = f_16 * skk_290[k]
                   + f_3 * pc_y[k] * slk_434[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, pc_x, pc_z, skk_255, skk_437, skk_438, sli0_341, \
                         sli0_342, sli1_341, sli1_342, slk_435, slk_437, \
                         slk_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_18 * skk_437[k]
                   + f_4 * sli0_341[k]
                   - f_5 * sli1_341[k]
                   + f_3 * pc_x[k] * slk_437[k];

        t_546[k] = f_18 * skk_438[k]
                   + f_6 * sli0_342[k]
                   - f_7 * sli1_342[k]
                   + f_3 * pc_x[k] * slk_438[k];

        t_547[k] = f_16 * skk_255[k]
                   + f_3 * pc_z[k] * slk_435[k];
    }

#pragma omp simd aligned(t_548, t_549, t_550, pc_x, pc_y, skk_293, skk_441, skk_442, sli0_345, \
                         sli0_346, sli1_345, sli1_346, slk_437, slk_441, \
                         slk_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_548[k] = f_16 * skk_293[k]
                   + f_3 * pc_y[k] * slk_437[k];

        t_549[k] = f_18 * skk_441[k]
                   + f_6 * sli0_345[k]
                   - f_7 * sli1_345[k]
                   + f_3 * pc_x[k] * slk_441[k];

        t_550[k] = f_18 * skk_442[k]
                   + f_8 * sli0_346[k]
                   - f_9 * sli1_346[k]
                   + f_3 * pc_x[k] * slk_442[k];
    }

#pragma omp simd aligned(t_551, t_552, t_553, pc_x, pc_y, pc_z, skk_258, skk_297, skk_444, \
                         sli0_348, sli1_348, slk_438, slk_441, \
                         slk_444 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_551[k] = f_16 * skk_258[k]
                   + f_3 * pc_z[k] * slk_438[k];

        t_552[k] = f_18 * skk_444[k]
                   + f_8 * sli0_348[k]
                   - f_9 * sli1_348[k]
                   + f_3 * pc_x[k] * slk_444[k];

        t_553[k] = f_16 * skk_297[k]
                   + f_3 * pc_y[k] * slk_441[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pc_x, pc_z, skk_262, skk_446, skk_447, sli0_350, \
                         sli0_351, sli1_350, sli1_351, slk_442, slk_446, \
                         slk_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_18 * skk_446[k]
                   + f_8 * sli0_350[k]
                   - f_9 * sli1_350[k]
                   + f_3 * pc_x[k] * slk_446[k];

        t_555[k] = f_18 * skk_447[k]
                   + f_10 * sli0_351[k]
                   - f_11 * sli1_351[k]
                   + f_3 * pc_x[k] * slk_447[k];

        t_556[k] = f_16 * skk_262[k]
                   + f_3 * pc_z[k] * slk_442[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, pc_y, skk_302, skk_449, skk_450, sli0_353, \
                         sli0_354, sli1_353, sli1_354, slk_446, slk_449, \
                         slk_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_18 * skk_449[k]
                   + f_10 * sli0_353[k]
                   - f_11 * sli1_353[k]
                   + f_3 * pc_x[k] * slk_449[k];

        t_558[k] = f_18 * skk_450[k]
                   + f_10 * sli0_354[k]
                   - f_11 * sli1_354[k]
                   + f_3 * pc_x[k] * slk_450[k];

        t_559[k] = f_16 * skk_302[k]
                   + f_3 * pc_y[k] * slk_446[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, pc_x, pc_z, skk_267, skk_452, skk_453, sli0_356, \
                         sli0_357, sli1_356, sli1_357, slk_447, slk_452, \
                         slk_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_18 * skk_452[k]
                   + f_10 * sli0_356[k]
                   - f_11 * sli1_356[k]
                   + f_3 * pc_x[k] * slk_452[k];

        t_561[k] = f_18 * skk_453[k]
                   + f_12 * sli0_357[k]
                   - f_13 * sli1_357[k]
                   + f_3 * pc_x[k] * slk_453[k];

        t_562[k] = f_16 * skk_267[k]
                   + f_3 * pc_z[k] * slk_447[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_405 = buffer.data(skl0 + 405);
    const auto *skl0_408 = buffer.data(skl0 + 408);
    const auto *skl0_410 = buffer.data(skl0 + 410);
    const auto *skl0_411 = buffer.data(skl0 + 411);
    const auto *skl0_414 = buffer.data(skl0 + 414);
    const auto *skl0_415 = buffer.data(skl0 + 415);
    const auto *skl0_417 = buffer.data(skl0 + 417);
    const auto *skl0_419 = buffer.data(skl0 + 419);
    const auto *skl0_420 = buffer.data(skl0 + 420);
    const auto *skl0_422 = buffer.data(skl0 + 422);
    const auto *skl0_423 = buffer.data(skl0 + 423);
    const auto *skl0_425 = buffer.data(skl0 + 425);
    const auto *skl0_426 = buffer.data(skl0 + 426);
    const auto *skl0_428 = buffer.data(skl0 + 428);
    const auto *skl0_429 = buffer.data(skl0 + 429);
    const auto *skl0_430 = buffer.data(skl0 + 430);
    const auto *skl0_432 = buffer.data(skl0 + 432);
    const auto *skl0_449 = buffer.data(skl0 + 449);

    const auto *skk_280 = buffer.data(skk + 280);
    const auto *skk_287 = buffer.data(skk + 287);
    const auto *skk_288 = buffer.data(skk + 288);
    const auto *skk_291 = buffer.data(skk + 291);
    const auto *skk_294 = buffer.data(skk + 294);
    const auto *skk_298 = buffer.data(skk + 298);
    const auto *skk_303 = buffer.data(skk + 303);
    const auto *skk_308 = buffer.data(skk + 308);
    const auto *skk_316 = buffer.data(skk + 316);
    const auto *skk_318 = buffer.data(skk + 318);
    const auto *skk_319 = buffer.data(skk + 319);
    const auto *skk_320 = buffer.data(skk + 320);
    const auto *skk_321 = buffer.data(skk + 321);
    const auto *skk_322 = buffer.data(skk + 322);
    const auto *skk_323 = buffer.data(skk + 323);
    const auto *skk_324 = buffer.data(skk + 324);
    const auto *skk_325 = buffer.data(skk + 325);
    const auto *skk_326 = buffer.data(skk + 326);
    const auto *skk_327 = buffer.data(skk + 327);
    const auto *skk_329 = buffer.data(skk + 329);
    const auto *skk_330 = buffer.data(skk + 330);
    const auto *skk_332 = buffer.data(skk + 332);
    const auto *skk_333 = buffer.data(skk + 333);
    const auto *skk_334 = buffer.data(skk + 334);
    const auto *skk_336 = buffer.data(skk + 336);
    const auto *skk_337 = buffer.data(skk + 337);
    const auto *skk_338 = buffer.data(skk + 338);
    const auto *skk_339 = buffer.data(skk + 339);
    const auto *skk_341 = buffer.data(skk + 341);
    const auto *skk_342 = buffer.data(skk + 342);
    const auto *skk_343 = buffer.data(skk + 343);
    const auto *skk_344 = buffer.data(skk + 344);
    const auto *skk_352 = buffer.data(skk + 352);
    const auto *skk_354 = buffer.data(skk + 354);
    const auto *skk_355 = buffer.data(skk + 355);
    const auto *skk_356 = buffer.data(skk + 356);
    const auto *skk_357 = buffer.data(skk + 357);
    const auto *skk_358 = buffer.data(skk + 358);
    const auto *skk_359 = buffer.data(skk + 359);
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
    const auto *skk_496 = buffer.data(skk + 496);
    const auto *skk_497 = buffer.data(skk + 497);
    const auto *skk_498 = buffer.data(skk + 498);
    const auto *skk_499 = buffer.data(skk + 499);
    const auto *skk_500 = buffer.data(skk + 500);
    const auto *skk_501 = buffer.data(skk + 501);
    const auto *skk_502 = buffer.data(skk + 502);
    const auto *skk_503 = buffer.data(skk + 503);
    const auto *skk_504 = buffer.data(skk + 504);
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

    const auto *skl1_405 = buffer.data(skl1 + 405);
    const auto *skl1_408 = buffer.data(skl1 + 408);
    const auto *skl1_410 = buffer.data(skl1 + 410);
    const auto *skl1_411 = buffer.data(skl1 + 411);
    const auto *skl1_414 = buffer.data(skl1 + 414);
    const auto *skl1_415 = buffer.data(skl1 + 415);
    const auto *skl1_417 = buffer.data(skl1 + 417);
    const auto *skl1_419 = buffer.data(skl1 + 419);
    const auto *skl1_420 = buffer.data(skl1 + 420);
    const auto *skl1_422 = buffer.data(skl1 + 422);
    const auto *skl1_423 = buffer.data(skl1 + 423);
    const auto *skl1_425 = buffer.data(skl1 + 425);
    const auto *skl1_426 = buffer.data(skl1 + 426);
    const auto *skl1_428 = buffer.data(skl1 + 428);
    const auto *skl1_429 = buffer.data(skl1 + 429);
    const auto *skl1_430 = buffer.data(skl1 + 430);
    const auto *skl1_432 = buffer.data(skl1 + 432);
    const auto *skl1_449 = buffer.data(skl1 + 449);

    const auto *sli0_357 = buffer.data(sli0 + 357);
    const auto *sli0_359 = buffer.data(sli0 + 359);
    const auto *sli0_360 = buffer.data(sli0 + 360);
    const auto *sli0_361 = buffer.data(sli0 + 361);
    const auto *sli0_362 = buffer.data(sli0 + 362);
    const auto *sli0_363 = buffer.data(sli0 + 363);
    const auto *sli0_385 = buffer.data(sli0 + 385);
    const auto *sli0_387 = buffer.data(sli0 + 387);
    const auto *sli0_388 = buffer.data(sli0 + 388);
    const auto *sli0_389 = buffer.data(sli0 + 389);
    const auto *sli0_390 = buffer.data(sli0 + 390);
    const auto *sli0_391 = buffer.data(sli0 + 391);
    const auto *sli0_392 = buffer.data(sli0 + 392);
    const auto *sli0_395 = buffer.data(sli0 + 395);
    const auto *sli0_397 = buffer.data(sli0 + 397);
    const auto *sli0_398 = buffer.data(sli0 + 398);
    const auto *sli0_401 = buffer.data(sli0 + 401);
    const auto *sli0_402 = buffer.data(sli0 + 402);
    const auto *sli0_404 = buffer.data(sli0 + 404);
    const auto *sli0_406 = buffer.data(sli0 + 406);
    const auto *sli0_407 = buffer.data(sli0 + 407);
    const auto *sli0_409 = buffer.data(sli0 + 409);
    const auto *sli0_410 = buffer.data(sli0 + 410);
    const auto *sli0_412 = buffer.data(sli0 + 412);
    const auto *sli0_413 = buffer.data(sli0 + 413);
    const auto *sli0_415 = buffer.data(sli0 + 415);
    const auto *sli0_416 = buffer.data(sli0 + 416);
    const auto *sli0_417 = buffer.data(sli0 + 417);
    const auto *sli0_418 = buffer.data(sli0 + 418);
    const auto *sli0_419 = buffer.data(sli0 + 419);

    const auto *sli1_357 = buffer.data(sli1 + 357);
    const auto *sli1_359 = buffer.data(sli1 + 359);
    const auto *sli1_360 = buffer.data(sli1 + 360);
    const auto *sli1_361 = buffer.data(sli1 + 361);
    const auto *sli1_362 = buffer.data(sli1 + 362);
    const auto *sli1_363 = buffer.data(sli1 + 363);
    const auto *sli1_385 = buffer.data(sli1 + 385);
    const auto *sli1_387 = buffer.data(sli1 + 387);
    const auto *sli1_388 = buffer.data(sli1 + 388);
    const auto *sli1_389 = buffer.data(sli1 + 389);
    const auto *sli1_390 = buffer.data(sli1 + 390);
    const auto *sli1_391 = buffer.data(sli1 + 391);
    const auto *sli1_392 = buffer.data(sli1 + 392);
    const auto *sli1_395 = buffer.data(sli1 + 395);
    const auto *sli1_397 = buffer.data(sli1 + 397);
    const auto *sli1_398 = buffer.data(sli1 + 398);
    const auto *sli1_401 = buffer.data(sli1 + 401);
    const auto *sli1_402 = buffer.data(sli1 + 402);
    const auto *sli1_404 = buffer.data(sli1 + 404);
    const auto *sli1_406 = buffer.data(sli1 + 406);
    const auto *sli1_407 = buffer.data(sli1 + 407);
    const auto *sli1_409 = buffer.data(sli1 + 409);
    const auto *sli1_410 = buffer.data(sli1 + 410);
    const auto *sli1_412 = buffer.data(sli1 + 412);
    const auto *sli1_413 = buffer.data(sli1 + 413);
    const auto *sli1_415 = buffer.data(sli1 + 415);
    const auto *sli1_416 = buffer.data(sli1 + 416);
    const auto *sli1_417 = buffer.data(sli1 + 417);
    const auto *sli1_418 = buffer.data(sli1 + 418);
    const auto *sli1_419 = buffer.data(sli1 + 419);

    const auto *slk_452 = buffer.data(slk + 452);
    const auto *slk_455 = buffer.data(slk + 455);
    const auto *slk_456 = buffer.data(slk + 456);
    const auto *slk_457 = buffer.data(slk + 457);
    const auto *slk_459 = buffer.data(slk + 459);
    const auto *slk_460 = buffer.data(slk + 460);
    const auto *slk_461 = buffer.data(slk + 461);
    const auto *slk_462 = buffer.data(slk + 462);
    const auto *slk_463 = buffer.data(slk + 463);
    const auto *slk_464 = buffer.data(slk + 464);
    const auto *slk_465 = buffer.data(slk + 465);
    const auto *slk_466 = buffer.data(slk + 466);
    const auto *slk_467 = buffer.data(slk + 467);
    const auto *slk_468 = buffer.data(slk + 468);
    const auto *slk_470 = buffer.data(slk + 470);
    const auto *slk_471 = buffer.data(slk + 471);
    const auto *slk_473 = buffer.data(slk + 473);
    const auto *slk_474 = buffer.data(slk + 474);
    const auto *slk_477 = buffer.data(slk + 477);
    const auto *slk_478 = buffer.data(slk + 478);
    const auto *slk_482 = buffer.data(slk + 482);
    const auto *slk_483 = buffer.data(slk + 483);
    const auto *slk_488 = buffer.data(slk + 488);
    const auto *slk_496 = buffer.data(slk + 496);
    const auto *slk_497 = buffer.data(slk + 497);
    const auto *slk_498 = buffer.data(slk + 498);
    const auto *slk_499 = buffer.data(slk + 499);
    const auto *slk_500 = buffer.data(slk + 500);
    const auto *slk_501 = buffer.data(slk + 501);
    const auto *slk_502 = buffer.data(slk + 502);
    const auto *slk_503 = buffer.data(slk + 503);
    const auto *slk_504 = buffer.data(slk + 504);
    const auto *slk_506 = buffer.data(slk + 506);
    const auto *slk_507 = buffer.data(slk + 507);
    const auto *slk_509 = buffer.data(slk + 509);
    const auto *slk_510 = buffer.data(slk + 510);
    const auto *slk_513 = buffer.data(slk + 513);
    const auto *slk_514 = buffer.data(slk + 514);
    const auto *slk_516 = buffer.data(slk + 516);
    const auto *slk_518 = buffer.data(slk + 518);
    const auto *slk_519 = buffer.data(slk + 519);
    const auto *slk_521 = buffer.data(slk + 521);
    const auto *slk_522 = buffer.data(slk + 522);
    const auto *slk_524 = buffer.data(slk + 524);
    const auto *slk_525 = buffer.data(slk + 525);
    const auto *slk_527 = buffer.data(slk + 527);
    const auto *slk_528 = buffer.data(slk + 528);
    const auto *slk_529 = buffer.data(slk + 529);
    const auto *slk_531 = buffer.data(slk + 531);
    const auto *slk_532 = buffer.data(slk + 532);
    const auto *slk_533 = buffer.data(slk + 533);
    const auto *slk_534 = buffer.data(slk + 534);
    const auto *slk_535 = buffer.data(slk + 535);
    const auto *slk_536 = buffer.data(slk + 536);
    const auto *slk_537 = buffer.data(slk + 537);
    const auto *slk_538 = buffer.data(slk + 538);
    const auto *slk_539 = buffer.data(slk + 539);

#pragma omp simd aligned(t_563, t_564, t_565, pc_x, skk_455, skk_456, skk_457, sli0_359, \
                         sli0_360, sli0_361, sli1_359, sli1_360, sli1_361, slk_455, slk_456, \
                         slk_457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_18 * skk_455[k]
                   + f_12 * sli0_359[k]
                   - f_13 * sli1_359[k]
                   + f_3 * pc_x[k] * slk_455[k];

        t_564[k] = f_18 * skk_456[k]
                   + f_12 * sli0_360[k]
                   - f_13 * sli1_360[k]
                   + f_3 * pc_x[k] * slk_456[k];

        t_565[k] = f_18 * skk_457[k]
                   + f_12 * sli0_361[k]
                   - f_13 * sli1_361[k]
                   + f_3 * pc_x[k] * slk_457[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, t_569, pc_x, pc_y, skk_308, skk_459, skk_460, \
                         skk_461, sli0_363, sli1_363, slk_452, slk_459, slk_460, \
                         slk_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = f_16 * skk_308[k]
                   + f_3 * pc_y[k] * slk_452[k];

        t_567[k] = f_18 * skk_459[k]
                   + f_12 * sli0_363[k]
                   - f_13 * sli1_363[k]
                   + f_3 * pc_x[k] * slk_459[k];

        t_568[k] = f_18 * skk_460[k]
                   + f_3 * pc_x[k] * slk_460[k];

        t_569[k] = f_18 * skk_461[k]
                   + f_3 * pc_x[k] * slk_461[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, t_574, pc_x, skk_462, skk_463, skk_464, \
                         skk_465, skk_466, slk_462, slk_463, slk_464, slk_465, \
                         slk_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_18 * skk_462[k]
                   + f_3 * pc_x[k] * slk_462[k];

        t_571[k] = f_18 * skk_463[k]
                   + f_3 * pc_x[k] * slk_463[k];

        t_572[k] = f_18 * skk_464[k]
                   + f_3 * pc_x[k] * slk_464[k];

        t_573[k] = f_18 * skk_465[k]
                   + f_3 * pc_x[k] * slk_465[k];

        t_574[k] = f_18 * skk_466[k]
                   + f_3 * pc_x[k] * slk_466[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, pc_z, skk_280, skk_316, skk_467, \
                         sli0_357, sli1_357, slk_460, slk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_18 * skk_467[k]
                   + f_3 * pc_x[k] * slk_467[k];

        t_576[k] = f_16 * skk_316[k]
                   + f_1 * sli0_357[k]
                   - f_2 * sli1_357[k]
                   + f_3 * pc_y[k] * slk_460[k];

        t_577[k] = f_16 * skk_280[k]
                   + f_3 * pc_z[k] * slk_460[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, pc_y, skk_318, skk_319, skk_320, sli0_359, \
                         sli0_360, sli0_361, sli1_359, sli1_360, sli1_361, slk_462, slk_463, \
                         slk_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_16 * skk_318[k]
                   + f_4 * sli0_359[k]
                   - f_5 * sli1_359[k]
                   + f_3 * pc_y[k] * slk_462[k];

        t_579[k] = f_16 * skk_319[k]
                   + f_6 * sli0_360[k]
                   - f_7 * sli1_360[k]
                   + f_3 * pc_y[k] * slk_463[k];

        t_580[k] = f_16 * skk_320[k]
                   + f_8 * sli0_361[k]
                   - f_9 * sli1_361[k]
                   + f_3 * pc_y[k] * slk_464[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, skk_321, skk_322, skk_323, sli0_362, \
                         sli0_363, sli1_362, sli1_363, slk_465, slk_466, \
                         slk_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_16 * skk_321[k]
                   + f_10 * sli0_362[k]
                   - f_11 * sli1_362[k]
                   + f_3 * pc_y[k] * slk_465[k];

        t_582[k] = f_16 * skk_322[k]
                   + f_12 * sli0_363[k]
                   - f_13 * sli1_363[k]
                   + f_3 * pc_y[k] * slk_466[k];

        t_583[k] = f_16 * skk_323[k]
                   + f_3 * pc_y[k] * slk_467[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, t_587, pb_y, pc_y, pc_z, skl0_405, skk_287, \
                         skk_288, skk_324, skl1_405, sli0_363, sli1_363, slk_467, \
                         slk_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_16 * skk_287[k]
                   + f_1 * sli0_363[k]
                   - f_2 * sli1_363[k]
                   + f_3 * pc_z[k] * slk_467[k];

        t_585[k] = pb_y[k] * skl0_405[k]
                   - f_14 * pc_y[k] * skl1_405[k];

        t_586[k] = f_15 * skk_324[k]
                   + f_3 * pc_y[k] * slk_468[k];

        t_587[k] = f_17 * skk_288[k]
                   + f_3 * pc_z[k] * slk_468[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, t_591, pb_y, pc_y, skl0_408, skl0_410, skl0_411, \
                         skk_325, skk_326, skk_327, skl1_408, skl1_410, skl1_411, \
                         slk_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pb_y[k] * skl0_408[k]
                   + f_16 * skk_325[k]
                   - f_14 * pc_y[k] * skl1_408[k];

        t_589[k] = f_15 * skk_326[k]
                   + f_3 * pc_y[k] * slk_470[k];

        t_590[k] = pb_y[k] * skl0_410[k]
                   - f_14 * pc_y[k] * skl1_410[k];

        t_591[k] = pb_y[k] * skl0_411[k]
                   + f_17 * skk_327[k]
                   - f_14 * pc_y[k] * skl1_411[k];
    }

#pragma omp simd aligned(t_592, t_593, t_594, t_595, pb_y, pc_y, pc_z, skl0_414, skl0_415, \
                         skk_291, skk_329, skk_330, skl1_414, skl1_415, slk_471, \
                         slk_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_592[k] = f_17 * skk_291[k]
                   + f_3 * pc_z[k] * slk_471[k];

        t_593[k] = f_15 * skk_329[k]
                   + f_3 * pc_y[k] * slk_473[k];

        t_594[k] = pb_y[k] * skl0_414[k]
                   - f_14 * pc_y[k] * skl1_414[k];

        t_595[k] = pb_y[k] * skl0_415[k]
                   + f_18 * skk_330[k]
                   - f_14 * pc_y[k] * skl1_415[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_y, pc_y, pc_z, skl0_417, skl0_419, \
                         skk_294, skk_332, skk_333, skl1_417, skl1_419, slk_474, \
                         slk_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_17 * skk_294[k]
                   + f_3 * pc_z[k] * slk_474[k];

        t_597[k] = pb_y[k] * skl0_417[k]
                   + f_16 * skk_332[k]
                   - f_14 * pc_y[k] * skl1_417[k];

        t_598[k] = f_15 * skk_333[k]
                   + f_3 * pc_y[k] * slk_477[k];

        t_599[k] = pb_y[k] * skl0_419[k]
                   - f_14 * pc_y[k] * skl1_419[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pb_y, pc_y, pc_z, skl0_420, skl0_422, skk_298, \
                         skk_334, skk_336, skl1_420, skl1_422, \
                         slk_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pb_y[k] * skl0_420[k]
                   + f_19 * skk_334[k]
                   - f_14 * pc_y[k] * skl1_420[k];

        t_601[k] = f_17 * skk_298[k]
                   + f_3 * pc_z[k] * slk_478[k];

        t_602[k] = pb_y[k] * skl0_422[k]
                   + f_17 * skk_336[k]
                   - f_14 * pc_y[k] * skl1_422[k];
    }

#pragma omp simd aligned(t_603, t_604, t_605, t_606, pb_y, pc_y, skl0_423, skl0_425, skl0_426, \
                         skk_337, skk_338, skk_339, skl1_423, skl1_425, skl1_426, \
                         slk_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = pb_y[k] * skl0_423[k]
                   + f_16 * skk_337[k]
                   - f_14 * pc_y[k] * skl1_423[k];

        t_604[k] = f_15 * skk_338[k]
                   + f_3 * pc_y[k] * slk_482[k];

        t_605[k] = pb_y[k] * skl0_425[k]
                   - f_14 * pc_y[k] * skl1_425[k];

        t_606[k] = pb_y[k] * skl0_426[k]
                   + f_20 * skk_339[k]
                   - f_14 * pc_y[k] * skl1_426[k];
    }

#pragma omp simd aligned(t_607, t_608, t_609, pb_y, pc_y, pc_z, skl0_428, skl0_429, skk_303, \
                         skk_341, skk_342, skl1_428, skl1_429, \
                         slk_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_607[k] = f_17 * skk_303[k]
                   + f_3 * pc_z[k] * slk_483[k];

        t_608[k] = pb_y[k] * skl0_428[k]
                   + f_18 * skk_341[k]
                   - f_14 * pc_y[k] * skl1_428[k];

        t_609[k] = pb_y[k] * skl0_429[k]
                   + f_17 * skk_342[k]
                   - f_14 * pc_y[k] * skl1_429[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pb_y, pc_x, pc_y, skl0_430, skl0_432, \
                         skk_343, skk_344, skk_496, skl1_430, skl1_432, slk_488, \
                         slk_496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pb_y[k] * skl0_430[k]
                   + f_16 * skk_343[k]
                   - f_14 * pc_y[k] * skl1_430[k];

        t_611[k] = f_15 * skk_344[k]
                   + f_3 * pc_y[k] * slk_488[k];

        t_612[k] = pb_y[k] * skl0_432[k]
                   - f_14 * pc_y[k] * skl1_432[k];

        t_613[k] = f_18 * skk_496[k]
                   + f_3 * pc_x[k] * slk_496[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, t_618, pc_x, skk_497, skk_498, skk_499, \
                         skk_500, skk_501, slk_497, slk_498, slk_499, slk_500, \
                         slk_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_18 * skk_497[k]
                   + f_3 * pc_x[k] * slk_497[k];

        t_615[k] = f_18 * skk_498[k]
                   + f_3 * pc_x[k] * slk_498[k];

        t_616[k] = f_18 * skk_499[k]
                   + f_3 * pc_x[k] * slk_499[k];

        t_617[k] = f_18 * skk_500[k]
                   + f_3 * pc_x[k] * slk_500[k];

        t_618[k] = f_18 * skk_501[k]
                   + f_3 * pc_x[k] * slk_501[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, pc_x, pc_y, pc_z, skk_316, skk_352, \
                         skk_502, skk_503, sli0_385, sli1_385, slk_496, slk_502, \
                         slk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_18 * skk_502[k]
                   + f_3 * pc_x[k] * slk_502[k];

        t_620[k] = f_18 * skk_503[k]
                   + f_3 * pc_x[k] * slk_503[k];

        t_621[k] = f_15 * skk_352[k]
                   + f_1 * sli0_385[k]
                   - f_2 * sli1_385[k]
                   + f_3 * pc_y[k] * slk_496[k];

        t_622[k] = f_17 * skk_316[k]
                   + f_3 * pc_z[k] * slk_496[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_y, skk_354, skk_355, skk_356, sli0_387, \
                         sli0_388, sli0_389, sli1_387, sli1_388, sli1_389, slk_498, slk_499, \
                         slk_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_15 * skk_354[k]
                   + f_4 * sli0_387[k]
                   - f_5 * sli1_387[k]
                   + f_3 * pc_y[k] * slk_498[k];

        t_624[k] = f_15 * skk_355[k]
                   + f_6 * sli0_388[k]
                   - f_7 * sli1_388[k]
                   + f_3 * pc_y[k] * slk_499[k];

        t_625[k] = f_15 * skk_356[k]
                   + f_8 * sli0_389[k]
                   - f_9 * sli1_389[k]
                   + f_3 * pc_y[k] * slk_500[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pc_y, skk_357, skk_358, skk_359, sli0_390, \
                         sli0_391, sli1_390, sli1_391, slk_501, slk_502, \
                         slk_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = f_15 * skk_357[k]
                   + f_10 * sli0_390[k]
                   - f_11 * sli1_390[k]
                   + f_3 * pc_y[k] * slk_501[k];

        t_627[k] = f_15 * skk_358[k]
                   + f_12 * sli0_391[k]
                   - f_13 * sli1_391[k]
                   + f_3 * pc_y[k] * slk_502[k];

        t_628[k] = f_15 * skk_359[k]
                   + f_3 * pc_y[k] * slk_503[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, pb_y, pc_x, pc_y, pc_z, skl0_449, \
                         skk_324, skk_504, skl1_449, sli0_392, sli1_392, \
                         slk_504 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = pb_y[k] * skl0_449[k]
                   - f_14 * pc_y[k] * skl1_449[k];

        t_630[k] = f_18 * skk_504[k]
                   + f_1 * sli0_392[k]
                   - f_2 * sli1_392[k]
                   + f_3 * pc_x[k] * slk_504[k];

        t_631[k] = f_3 * pc_y[k] * slk_504[k];

        t_632[k] = f_18 * skk_324[k]
                   + f_3 * pc_z[k] * slk_504[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pc_x, pc_y, skk_507, skk_509, sli0_395, \
                         sli0_397, sli1_395, sli1_397, slk_506, slk_507, \
                         slk_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = f_18 * skk_507[k]
                   + f_4 * sli0_395[k]
                   - f_5 * sli1_395[k]
                   + f_3 * pc_x[k] * slk_507[k];

        t_634[k] = f_3 * pc_y[k] * slk_506[k];

        t_635[k] = f_18 * skk_509[k]
                   + f_4 * sli0_397[k]
                   - f_5 * sli1_397[k]
                   + f_3 * pc_x[k] * slk_509[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pc_x, pc_y, pc_z, skk_327, skk_510, sli0_398, \
                         sli1_398, slk_507, slk_509, slk_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = f_18 * skk_510[k]
                   + f_6 * sli0_398[k]
                   - f_7 * sli1_398[k]
                   + f_3 * pc_x[k] * slk_510[k];

        t_637[k] = f_18 * skk_327[k]
                   + f_3 * pc_z[k] * slk_507[k];

        t_638[k] = f_3 * pc_y[k] * slk_509[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_z, skk_330, skk_513, skk_514, sli0_401, \
                         sli0_402, sli1_401, sli1_402, slk_510, slk_513, \
                         slk_514 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_18 * skk_513[k]
                   + f_6 * sli0_401[k]
                   - f_7 * sli1_401[k]
                   + f_3 * pc_x[k] * slk_513[k];

        t_640[k] = f_18 * skk_514[k]
                   + f_8 * sli0_402[k]
                   - f_9 * sli1_402[k]
                   + f_3 * pc_x[k] * slk_514[k];

        t_641[k] = f_18 * skk_330[k]
                   + f_3 * pc_z[k] * slk_510[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, skk_516, skk_518, sli0_404, \
                         sli0_406, sli1_404, sli1_406, slk_513, slk_516, \
                         slk_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_18 * skk_516[k]
                   + f_8 * sli0_404[k]
                   - f_9 * sli1_404[k]
                   + f_3 * pc_x[k] * slk_516[k];

        t_643[k] = f_3 * pc_y[k] * slk_513[k];

        t_644[k] = f_18 * skk_518[k]
                   + f_8 * sli0_406[k]
                   - f_9 * sli1_406[k]
                   + f_3 * pc_x[k] * slk_518[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, pc_z, skk_334, skk_519, skk_521, sli0_407, \
                         sli0_409, sli1_407, sli1_409, slk_514, slk_519, \
                         slk_521 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_18 * skk_519[k]
                   + f_10 * sli0_407[k]
                   - f_11 * sli1_407[k]
                   + f_3 * pc_x[k] * slk_519[k];

        t_646[k] = f_18 * skk_334[k]
                   + f_3 * pc_z[k] * slk_514[k];

        t_647[k] = f_18 * skk_521[k]
                   + f_10 * sli0_409[k]
                   - f_11 * sli1_409[k]
                   + f_3 * pc_x[k] * slk_521[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, pc_y, skk_522, skk_524, sli0_410, \
                         sli0_412, sli1_410, sli1_412, slk_518, slk_522, \
                         slk_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_18 * skk_522[k]
                   + f_10 * sli0_410[k]
                   - f_11 * sli1_410[k]
                   + f_3 * pc_x[k] * slk_522[k];

        t_649[k] = f_3 * pc_y[k] * slk_518[k];

        t_650[k] = f_18 * skk_524[k]
                   + f_10 * sli0_412[k]
                   - f_11 * sli1_412[k]
                   + f_3 * pc_x[k] * slk_524[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, pc_z, skk_339, skk_525, skk_527, sli0_413, \
                         sli0_415, sli1_413, sli1_415, slk_519, slk_525, \
                         slk_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_18 * skk_525[k]
                   + f_12 * sli0_413[k]
                   - f_13 * sli1_413[k]
                   + f_3 * pc_x[k] * slk_525[k];

        t_652[k] = f_18 * skk_339[k]
                   + f_3 * pc_z[k] * slk_519[k];

        t_653[k] = f_18 * skk_527[k]
                   + f_12 * sli0_415[k]
                   - f_13 * sli1_415[k]
                   + f_3 * pc_x[k] * slk_527[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, skk_528, skk_529, sli0_416, \
                         sli0_417, sli1_416, sli1_417, slk_524, slk_528, \
                         slk_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_18 * skk_528[k]
                   + f_12 * sli0_416[k]
                   - f_13 * sli1_416[k]
                   + f_3 * pc_x[k] * slk_528[k];

        t_655[k] = f_18 * skk_529[k]
                   + f_12 * sli0_417[k]
                   - f_13 * sli1_417[k]
                   + f_3 * pc_x[k] * slk_529[k];

        t_656[k] = f_3 * pc_y[k] * slk_524[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, pc_x, skk_531, skk_532, skk_533, skk_534, \
                         sli0_419, sli1_419, slk_531, slk_532, slk_533, \
                         slk_534 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_18 * skk_531[k]
                   + f_12 * sli0_419[k]
                   - f_13 * sli1_419[k]
                   + f_3 * pc_x[k] * slk_531[k];

        t_658[k] = f_18 * skk_532[k]
                   + f_3 * pc_x[k] * slk_532[k];

        t_659[k] = f_18 * skk_533[k]
                   + f_3 * pc_x[k] * slk_533[k];

        t_660[k] = f_18 * skk_534[k]
                   + f_3 * pc_x[k] * slk_534[k];
    }

#pragma omp simd aligned(t_661, t_662, t_663, t_664, t_665, pc_x, skk_535, skk_536, skk_537, \
                         skk_538, skk_539, slk_535, slk_536, slk_537, slk_538, \
                         slk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_661[k] = f_18 * skk_535[k]
                   + f_3 * pc_x[k] * slk_535[k];

        t_662[k] = f_18 * skk_536[k]
                   + f_3 * pc_x[k] * slk_536[k];

        t_663[k] = f_18 * skk_537[k]
                   + f_3 * pc_x[k] * slk_537[k];

        t_664[k] = f_18 * skk_538[k]
                   + f_3 * pc_x[k] * slk_538[k];

        t_665[k] = f_18 * skk_539[k]
                   + f_3 * pc_x[k] * slk_539[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, skk_352, sli0_413, sli0_415, \
                         sli0_416, sli1_413, sli1_415, sli1_416, slk_532, slk_534, \
                         slk_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * sli0_413[k]
                   - f_2 * sli1_413[k]
                   + f_3 * pc_y[k] * slk_532[k];

        t_667[k] = f_18 * skk_352[k]
                   + f_3 * pc_z[k] * slk_532[k];

        t_668[k] = f_4 * sli0_415[k]
                   - f_5 * sli1_415[k]
                   + f_3 * pc_y[k] * slk_534[k];

        t_669[k] = f_6 * sli0_416[k]
                   - f_7 * sli1_416[k]
                   + f_3 * pc_y[k] * slk_535[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, sli0_417, sli0_418, sli0_419, \
                         sli1_417, sli1_418, sli1_419, slk_536, slk_537, slk_538, \
                         slk_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_8 * sli0_417[k]
                   - f_9 * sli1_417[k]
                   + f_3 * pc_y[k] * slk_536[k];

        t_671[k] = f_10 * sli0_418[k]
                   - f_11 * sli1_418[k]
                   + f_3 * pc_y[k] * slk_537[k];

        t_672[k] = f_12 * sli0_419[k]
                   - f_13 * sli1_419[k]
                   + f_3 * pc_y[k] * slk_538[k];

        t_673[k] = f_3 * pc_y[k] * slk_539[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_450 = buffer.data(skl0 + 450);
    const auto *skl0_453 = buffer.data(skl0 + 453);
    const auto *skl0_456 = buffer.data(skl0 + 456);
    const auto *skl0_460 = buffer.data(skl0 + 460);
    const auto *skl0_462 = buffer.data(skl0 + 462);
    const auto *skl0_465 = buffer.data(skl0 + 465);
    const auto *skl0_467 = buffer.data(skl0 + 467);
    const auto *skl0_468 = buffer.data(skl0 + 468);
    const auto *skl0_471 = buffer.data(skl0 + 471);
    const auto *skl0_473 = buffer.data(skl0 + 473);
    const auto *skl0_474 = buffer.data(skl0 + 474);
    const auto *skl0_475 = buffer.data(skl0 + 475);
    const auto *skl0_486 = buffer.data(skl0 + 486);

    const auto *skk_359 = buffer.data(skk + 359);
    const auto *skk_360 = buffer.data(skk + 360);
    const auto *skk_362 = buffer.data(skk + 362);
    const auto *skk_363 = buffer.data(skk + 363);
    const auto *skk_365 = buffer.data(skk + 365);
    const auto *skk_366 = buffer.data(skk + 366);
    const auto *skk_367 = buffer.data(skk + 367);
    const auto *skk_369 = buffer.data(skk + 369);
    const auto *skk_370 = buffer.data(skk + 370);
    const auto *skk_371 = buffer.data(skk + 371);
    const auto *skk_372 = buffer.data(skk + 372);
    const auto *skk_374 = buffer.data(skk + 374);
    const auto *skk_375 = buffer.data(skk + 375);
    const auto *skk_376 = buffer.data(skk + 376);
    const auto *skk_377 = buffer.data(skk + 377);
    const auto *skk_378 = buffer.data(skk + 378);
    const auto *skk_380 = buffer.data(skk + 380);
    const auto *skk_388 = buffer.data(skk + 388);
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
    const auto *skk_416 = buffer.data(skk + 416);
    const auto *skk_426 = buffer.data(skk + 426);
    const auto *skk_427 = buffer.data(skk + 427);
    const auto *skk_428 = buffer.data(skk + 428);
    const auto *skk_429 = buffer.data(skk + 429);
    const auto *skk_430 = buffer.data(skk + 430);
    const auto *skk_431 = buffer.data(skk + 431);
    const auto *skk_432 = buffer.data(skk + 432);
    const auto *skk_434 = buffer.data(skk + 434);
    const auto *skk_437 = buffer.data(skk + 437);
    const auto *skk_441 = buffer.data(skk + 441);
    const auto *skk_540 = buffer.data(skk + 540);
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
    const auto *skk_581 = buffer.data(skk + 581);
    const auto *skk_585 = buffer.data(skk + 585);
    const auto *skk_590 = buffer.data(skk + 590);
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
    const auto *skk_615 = buffer.data(skk + 615);
    const auto *skk_617 = buffer.data(skk + 617);
    const auto *skk_618 = buffer.data(skk + 618);
    const auto *skk_621 = buffer.data(skk + 621);
    const auto *skk_622 = buffer.data(skk + 622);
    const auto *skk_624 = buffer.data(skk + 624);
    const auto *skk_626 = buffer.data(skk + 626);
    const auto *skk_627 = buffer.data(skk + 627);

    const auto *skl1_450 = buffer.data(skl1 + 450);
    const auto *skl1_453 = buffer.data(skl1 + 453);
    const auto *skl1_456 = buffer.data(skl1 + 456);
    const auto *skl1_460 = buffer.data(skl1 + 460);
    const auto *skl1_462 = buffer.data(skl1 + 462);
    const auto *skl1_465 = buffer.data(skl1 + 465);
    const auto *skl1_467 = buffer.data(skl1 + 467);
    const auto *skl1_468 = buffer.data(skl1 + 468);
    const auto *skl1_471 = buffer.data(skl1 + 471);
    const auto *skl1_473 = buffer.data(skl1 + 473);
    const auto *skl1_474 = buffer.data(skl1 + 474);
    const auto *skl1_475 = buffer.data(skl1 + 475);
    const auto *skl1_486 = buffer.data(skl1 + 486);

    const auto *sli0_419 = buffer.data(sli0 + 419);
    const auto *sli0_420 = buffer.data(sli0 + 420);
    const auto *sli0_423 = buffer.data(sli0 + 423);
    const auto *sli0_425 = buffer.data(sli0 + 425);
    const auto *sli0_426 = buffer.data(sli0 + 426);
    const auto *sli0_429 = buffer.data(sli0 + 429);
    const auto *sli0_430 = buffer.data(sli0 + 430);
    const auto *sli0_432 = buffer.data(sli0 + 432);
    const auto *sli0_434 = buffer.data(sli0 + 434);
    const auto *sli0_435 = buffer.data(sli0 + 435);
    const auto *sli0_437 = buffer.data(sli0 + 437);
    const auto *sli0_438 = buffer.data(sli0 + 438);
    const auto *sli0_440 = buffer.data(sli0 + 440);
    const auto *sli0_441 = buffer.data(sli0 + 441);
    const auto *sli0_443 = buffer.data(sli0 + 443);
    const auto *sli0_444 = buffer.data(sli0 + 444);
    const auto *sli0_445 = buffer.data(sli0 + 445);
    const auto *sli0_446 = buffer.data(sli0 + 446);
    const auto *sli0_447 = buffer.data(sli0 + 447);
    const auto *sli0_453 = buffer.data(sli0 + 453);
    const auto *sli0_457 = buffer.data(sli0 + 457);
    const auto *sli0_462 = buffer.data(sli0 + 462);
    const auto *sli0_468 = buffer.data(sli0 + 468);
    const auto *sli0_471 = buffer.data(sli0 + 471);
    const auto *sli0_472 = buffer.data(sli0 + 472);
    const auto *sli0_473 = buffer.data(sli0 + 473);
    const auto *sli0_474 = buffer.data(sli0 + 474);
    const auto *sli0_475 = buffer.data(sli0 + 475);
    const auto *sli0_476 = buffer.data(sli0 + 476);
    const auto *sli0_479 = buffer.data(sli0 + 479);
    const auto *sli0_481 = buffer.data(sli0 + 481);
    const auto *sli0_482 = buffer.data(sli0 + 482);
    const auto *sli0_485 = buffer.data(sli0 + 485);
    const auto *sli0_486 = buffer.data(sli0 + 486);
    const auto *sli0_488 = buffer.data(sli0 + 488);
    const auto *sli0_490 = buffer.data(sli0 + 490);
    const auto *sli0_491 = buffer.data(sli0 + 491);

    const auto *sli1_419 = buffer.data(sli1 + 419);
    const auto *sli1_420 = buffer.data(sli1 + 420);
    const auto *sli1_423 = buffer.data(sli1 + 423);
    const auto *sli1_425 = buffer.data(sli1 + 425);
    const auto *sli1_426 = buffer.data(sli1 + 426);
    const auto *sli1_429 = buffer.data(sli1 + 429);
    const auto *sli1_430 = buffer.data(sli1 + 430);
    const auto *sli1_432 = buffer.data(sli1 + 432);
    const auto *sli1_434 = buffer.data(sli1 + 434);
    const auto *sli1_435 = buffer.data(sli1 + 435);
    const auto *sli1_437 = buffer.data(sli1 + 437);
    const auto *sli1_438 = buffer.data(sli1 + 438);
    const auto *sli1_440 = buffer.data(sli1 + 440);
    const auto *sli1_441 = buffer.data(sli1 + 441);
    const auto *sli1_443 = buffer.data(sli1 + 443);
    const auto *sli1_444 = buffer.data(sli1 + 444);
    const auto *sli1_445 = buffer.data(sli1 + 445);
    const auto *sli1_446 = buffer.data(sli1 + 446);
    const auto *sli1_447 = buffer.data(sli1 + 447);
    const auto *sli1_453 = buffer.data(sli1 + 453);
    const auto *sli1_457 = buffer.data(sli1 + 457);
    const auto *sli1_462 = buffer.data(sli1 + 462);
    const auto *sli1_468 = buffer.data(sli1 + 468);
    const auto *sli1_471 = buffer.data(sli1 + 471);
    const auto *sli1_472 = buffer.data(sli1 + 472);
    const auto *sli1_473 = buffer.data(sli1 + 473);
    const auto *sli1_474 = buffer.data(sli1 + 474);
    const auto *sli1_475 = buffer.data(sli1 + 475);
    const auto *sli1_476 = buffer.data(sli1 + 476);
    const auto *sli1_479 = buffer.data(sli1 + 479);
    const auto *sli1_481 = buffer.data(sli1 + 481);
    const auto *sli1_482 = buffer.data(sli1 + 482);
    const auto *sli1_485 = buffer.data(sli1 + 485);
    const auto *sli1_486 = buffer.data(sli1 + 486);
    const auto *sli1_488 = buffer.data(sli1 + 488);
    const auto *sli1_490 = buffer.data(sli1 + 490);
    const auto *sli1_491 = buffer.data(sli1 + 491);

    const auto *slk_539 = buffer.data(slk + 539);
    const auto *slk_540 = buffer.data(slk + 540);
    const auto *slk_542 = buffer.data(slk + 542);
    const auto *slk_543 = buffer.data(slk + 543);
    const auto *slk_545 = buffer.data(slk + 545);
    const auto *slk_546 = buffer.data(slk + 546);
    const auto *slk_549 = buffer.data(slk + 549);
    const auto *slk_550 = buffer.data(slk + 550);
    const auto *slk_552 = buffer.data(slk + 552);
    const auto *slk_554 = buffer.data(slk + 554);
    const auto *slk_555 = buffer.data(slk + 555);
    const auto *slk_557 = buffer.data(slk + 557);
    const auto *slk_558 = buffer.data(slk + 558);
    const auto *slk_560 = buffer.data(slk + 560);
    const auto *slk_561 = buffer.data(slk + 561);
    const auto *slk_563 = buffer.data(slk + 563);
    const auto *slk_564 = buffer.data(slk + 564);
    const auto *slk_565 = buffer.data(slk + 565);
    const auto *slk_567 = buffer.data(slk + 567);
    const auto *slk_568 = buffer.data(slk + 568);
    const auto *slk_569 = buffer.data(slk + 569);
    const auto *slk_570 = buffer.data(slk + 570);
    const auto *slk_571 = buffer.data(slk + 571);
    const auto *slk_572 = buffer.data(slk + 572);
    const auto *slk_573 = buffer.data(slk + 573);
    const auto *slk_574 = buffer.data(slk + 574);
    const auto *slk_575 = buffer.data(slk + 575);
    const auto *slk_576 = buffer.data(slk + 576);
    const auto *slk_578 = buffer.data(slk + 578);
    const auto *slk_579 = buffer.data(slk + 579);
    const auto *slk_581 = buffer.data(slk + 581);
    const auto *slk_582 = buffer.data(slk + 582);
    const auto *slk_585 = buffer.data(slk + 585);
    const auto *slk_586 = buffer.data(slk + 586);
    const auto *slk_590 = buffer.data(slk + 590);
    const auto *slk_591 = buffer.data(slk + 591);
    const auto *slk_596 = buffer.data(slk + 596);
    const auto *slk_603 = buffer.data(slk + 603);
    const auto *slk_604 = buffer.data(slk + 604);
    const auto *slk_605 = buffer.data(slk + 605);
    const auto *slk_606 = buffer.data(slk + 606);
    const auto *slk_607 = buffer.data(slk + 607);
    const auto *slk_608 = buffer.data(slk + 608);
    const auto *slk_609 = buffer.data(slk + 609);
    const auto *slk_610 = buffer.data(slk + 610);
    const auto *slk_611 = buffer.data(slk + 611);
    const auto *slk_612 = buffer.data(slk + 612);
    const auto *slk_614 = buffer.data(slk + 614);
    const auto *slk_615 = buffer.data(slk + 615);
    const auto *slk_617 = buffer.data(slk + 617);
    const auto *slk_618 = buffer.data(slk + 618);
    const auto *slk_621 = buffer.data(slk + 621);
    const auto *slk_622 = buffer.data(slk + 622);
    const auto *slk_624 = buffer.data(slk + 624);
    const auto *slk_626 = buffer.data(slk + 626);
    const auto *slk_627 = buffer.data(slk + 627);

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pc_x, pc_y, pc_z, skk_359, skk_360, \
                         skk_540, sli0_419, sli0_420, sli1_419, sli1_420, slk_539, \
                         slk_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_18 * skk_359[k]
                   + f_1 * sli0_419[k]
                   - f_2 * sli1_419[k]
                   + f_3 * pc_z[k] * slk_539[k];

        t_675[k] = f_17 * skk_540[k]
                   + f_1 * sli0_420[k]
                   - f_2 * sli1_420[k]
                   + f_3 * pc_x[k] * slk_540[k];

        t_676[k] = f_19 * skk_360[k]
                   + f_3 * pc_y[k] * slk_540[k];

        t_677[k] = f_3 * pc_z[k] * slk_540[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pc_x, pc_y, skk_362, skk_543, skk_545, sli0_423, \
                         sli0_425, sli1_423, sli1_425, slk_542, slk_543, \
                         slk_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = f_17 * skk_543[k]
                   + f_4 * sli0_423[k]
                   - f_5 * sli1_423[k]
                   + f_3 * pc_x[k] * slk_543[k];

        t_679[k] = f_19 * skk_362[k]
                   + f_3 * pc_y[k] * slk_542[k];

        t_680[k] = f_17 * skk_545[k]
                   + f_4 * sli0_425[k]
                   - f_5 * sli1_425[k]
                   + f_3 * pc_x[k] * slk_545[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, pc_x, pc_y, pc_z, skk_365, skk_546, sli0_426, \
                         sli1_426, slk_543, slk_545, slk_546 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_17 * skk_546[k]
                   + f_6 * sli0_426[k]
                   - f_7 * sli1_426[k]
                   + f_3 * pc_x[k] * slk_546[k];

        t_682[k] = f_3 * pc_z[k] * slk_543[k];

        t_683[k] = f_19 * skk_365[k]
                   + f_3 * pc_y[k] * slk_545[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, pc_x, pc_z, skk_549, skk_550, sli0_429, \
                         sli0_430, sli1_429, sli1_430, slk_546, slk_549, \
                         slk_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_17 * skk_549[k]
                   + f_6 * sli0_429[k]
                   - f_7 * sli1_429[k]
                   + f_3 * pc_x[k] * slk_549[k];

        t_685[k] = f_17 * skk_550[k]
                   + f_8 * sli0_430[k]
                   - f_9 * sli1_430[k]
                   + f_3 * pc_x[k] * slk_550[k];

        t_686[k] = f_3 * pc_z[k] * slk_546[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, pc_x, pc_y, skk_369, skk_552, skk_554, sli0_432, \
                         sli0_434, sli1_432, sli1_434, slk_549, slk_552, \
                         slk_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = f_17 * skk_552[k]
                   + f_8 * sli0_432[k]
                   - f_9 * sli1_432[k]
                   + f_3 * pc_x[k] * slk_552[k];

        t_688[k] = f_19 * skk_369[k]
                   + f_3 * pc_y[k] * slk_549[k];

        t_689[k] = f_17 * skk_554[k]
                   + f_8 * sli0_434[k]
                   - f_9 * sli1_434[k]
                   + f_3 * pc_x[k] * slk_554[k];
    }

#pragma omp simd aligned(t_690, t_691, t_692, pc_x, pc_z, skk_555, skk_557, sli0_435, \
                         sli0_437, sli1_435, sli1_437, slk_550, slk_555, \
                         slk_557 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_690[k] = f_17 * skk_555[k]
                   + f_10 * sli0_435[k]
                   - f_11 * sli1_435[k]
                   + f_3 * pc_x[k] * slk_555[k];

        t_691[k] = f_3 * pc_z[k] * slk_550[k];

        t_692[k] = f_17 * skk_557[k]
                   + f_10 * sli0_437[k]
                   - f_11 * sli1_437[k]
                   + f_3 * pc_x[k] * slk_557[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pc_x, pc_y, skk_374, skk_558, skk_560, sli0_438, \
                         sli0_440, sli1_438, sli1_440, slk_554, slk_558, \
                         slk_560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = f_17 * skk_558[k]
                   + f_10 * sli0_438[k]
                   - f_11 * sli1_438[k]
                   + f_3 * pc_x[k] * slk_558[k];

        t_694[k] = f_19 * skk_374[k]
                   + f_3 * pc_y[k] * slk_554[k];

        t_695[k] = f_17 * skk_560[k]
                   + f_10 * sli0_440[k]
                   - f_11 * sli1_440[k]
                   + f_3 * pc_x[k] * slk_560[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, pc_x, pc_z, skk_561, skk_563, sli0_441, \
                         sli0_443, sli1_441, sli1_443, slk_555, slk_561, \
                         slk_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_17 * skk_561[k]
                   + f_12 * sli0_441[k]
                   - f_13 * sli1_441[k]
                   + f_3 * pc_x[k] * slk_561[k];

        t_697[k] = f_3 * pc_z[k] * slk_555[k];

        t_698[k] = f_17 * skk_563[k]
                   + f_12 * sli0_443[k]
                   - f_13 * sli1_443[k]
                   + f_3 * pc_x[k] * slk_563[k];
    }

#pragma omp simd aligned(t_699, t_700, t_701, pc_x, pc_y, skk_380, skk_564, skk_565, sli0_444, \
                         sli0_445, sli1_444, sli1_445, slk_560, slk_564, \
                         slk_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_699[k] = f_17 * skk_564[k]
                   + f_12 * sli0_444[k]
                   - f_13 * sli1_444[k]
                   + f_3 * pc_x[k] * slk_564[k];

        t_700[k] = f_17 * skk_565[k]
                   + f_12 * sli0_445[k]
                   - f_13 * sli1_445[k]
                   + f_3 * pc_x[k] * slk_565[k];

        t_701[k] = f_19 * skk_380[k]
                   + f_3 * pc_y[k] * slk_560[k];
    }

#pragma omp simd aligned(t_702, t_703, t_704, t_705, pc_x, skk_567, skk_568, skk_569, skk_570, \
                         sli0_447, sli1_447, slk_567, slk_568, slk_569, \
                         slk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_702[k] = f_17 * skk_567[k]
                   + f_12 * sli0_447[k]
                   - f_13 * sli1_447[k]
                   + f_3 * pc_x[k] * slk_567[k];

        t_703[k] = f_17 * skk_568[k]
                   + f_3 * pc_x[k] * slk_568[k];

        t_704[k] = f_17 * skk_569[k]
                   + f_3 * pc_x[k] * slk_569[k];

        t_705[k] = f_17 * skk_570[k]
                   + f_3 * pc_x[k] * slk_570[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, t_710, pc_x, skk_571, skk_572, skk_573, \
                         skk_574, skk_575, slk_571, slk_572, slk_573, slk_574, \
                         slk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_17 * skk_571[k]
                   + f_3 * pc_x[k] * slk_571[k];

        t_707[k] = f_17 * skk_572[k]
                   + f_3 * pc_x[k] * slk_572[k];

        t_708[k] = f_17 * skk_573[k]
                   + f_3 * pc_x[k] * slk_573[k];

        t_709[k] = f_17 * skk_574[k]
                   + f_3 * pc_x[k] * slk_574[k];

        t_710[k] = f_17 * skk_575[k]
                   + f_3 * pc_x[k] * slk_575[k];
    }

#pragma omp simd aligned(t_711, t_712, t_713, pc_y, pc_z, skk_388, skk_390, sli0_441, \
                         sli0_443, sli1_441, sli1_443, slk_568, \
                         slk_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_711[k] = f_19 * skk_388[k]
                   + f_1 * sli0_441[k]
                   - f_2 * sli1_441[k]
                   + f_3 * pc_y[k] * slk_568[k];

        t_712[k] = f_3 * pc_z[k] * slk_568[k];

        t_713[k] = f_19 * skk_390[k]
                   + f_4 * sli0_443[k]
                   - f_5 * sli1_443[k]
                   + f_3 * pc_y[k] * slk_570[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_y, skk_391, skk_392, skk_393, sli0_444, \
                         sli0_445, sli0_446, sli1_444, sli1_445, sli1_446, slk_571, slk_572, \
                         slk_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_19 * skk_391[k]
                   + f_6 * sli0_444[k]
                   - f_7 * sli1_444[k]
                   + f_3 * pc_y[k] * slk_571[k];

        t_715[k] = f_19 * skk_392[k]
                   + f_8 * sli0_445[k]
                   - f_9 * sli1_445[k]
                   + f_3 * pc_y[k] * slk_572[k];

        t_716[k] = f_19 * skk_393[k]
                   + f_10 * sli0_446[k]
                   - f_11 * sli1_446[k]
                   + f_3 * pc_y[k] * slk_573[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, t_720, pb_z, pc_y, pc_z, skl0_450, skk_394, \
                         skk_395, skl1_450, sli0_447, sli1_447, slk_574, \
                         slk_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_19 * skk_394[k]
                   + f_12 * sli0_447[k]
                   - f_13 * sli1_447[k]
                   + f_3 * pc_y[k] * slk_574[k];

        t_718[k] = f_19 * skk_395[k]
                   + f_3 * pc_y[k] * slk_575[k];

        t_719[k] = f_1 * sli0_447[k]
                   - f_2 * sli1_447[k]
                   + f_3 * pc_z[k] * slk_575[k];

        t_720[k] = pb_z[k] * skl0_450[k]
                   - f_14 * pc_z[k] * skl1_450[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pb_z, pc_y, pc_z, skl0_453, skk_360, \
                         skk_396, skk_398, skl1_453, slk_576, slk_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_18 * skk_396[k]
                   + f_3 * pc_y[k] * slk_576[k];

        t_722[k] = f_15 * skk_360[k]
                   + f_3 * pc_z[k] * slk_576[k];

        t_723[k] = pb_z[k] * skl0_453[k]
                   - f_14 * pc_z[k] * skl1_453[k];

        t_724[k] = f_18 * skk_398[k]
                   + f_3 * pc_y[k] * slk_578[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, pb_z, pc_x, pc_z, skl0_456, skk_363, skk_581, \
                         skl1_456, sli0_453, sli1_453, slk_579, \
                         slk_581 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_17 * skk_581[k]
                   + f_4 * sli0_453[k]
                   - f_5 * sli1_453[k]
                   + f_3 * pc_x[k] * slk_581[k];

        t_726[k] = pb_z[k] * skl0_456[k]
                   - f_14 * pc_z[k] * skl1_456[k];

        t_727[k] = f_15 * skk_363[k]
                   + f_3 * pc_z[k] * slk_579[k];
    }

#pragma omp simd aligned(t_728, t_729, t_730, pb_z, pc_x, pc_y, pc_z, skl0_460, skk_401, \
                         skk_585, skl1_460, sli0_457, sli1_457, slk_581, \
                         slk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_728[k] = f_18 * skk_401[k]
                   + f_3 * pc_y[k] * slk_581[k];

        t_729[k] = f_17 * skk_585[k]
                   + f_6 * sli0_457[k]
                   - f_7 * sli1_457[k]
                   + f_3 * pc_x[k] * slk_585[k];

        t_730[k] = pb_z[k] * skl0_460[k]
                   - f_14 * pc_z[k] * skl1_460[k];
    }

#pragma omp simd aligned(t_731, t_732, t_733, pb_z, pc_y, pc_z, skl0_462, skk_366, skk_367, \
                         skk_405, skl1_462, slk_582, slk_585 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_731[k] = f_15 * skk_366[k]
                   + f_3 * pc_z[k] * slk_582[k];

        t_732[k] = pb_z[k] * skl0_462[k]
                   + f_16 * skk_367[k]
                   - f_14 * pc_z[k] * skl1_462[k];

        t_733[k] = f_18 * skk_405[k]
                   + f_3 * pc_y[k] * slk_585[k];
    }

#pragma omp simd aligned(t_734, t_735, t_736, pb_z, pc_x, pc_z, skl0_465, skk_370, skk_590, \
                         skl1_465, sli0_462, sli1_462, slk_586, \
                         slk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_734[k] = f_17 * skk_590[k]
                   + f_8 * sli0_462[k]
                   - f_9 * sli1_462[k]
                   + f_3 * pc_x[k] * slk_590[k];

        t_735[k] = pb_z[k] * skl0_465[k]
                   - f_14 * pc_z[k] * skl1_465[k];

        t_736[k] = f_15 * skk_370[k]
                   + f_3 * pc_z[k] * slk_586[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pb_z, pc_y, pc_z, skl0_467, skl0_468, skk_371, \
                         skk_372, skk_410, skl1_467, skl1_468, \
                         slk_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = pb_z[k] * skl0_467[k]
                   + f_16 * skk_371[k]
                   - f_14 * pc_z[k] * skl1_467[k];

        t_738[k] = pb_z[k] * skl0_468[k]
                   + f_17 * skk_372[k]
                   - f_14 * pc_z[k] * skl1_468[k];

        t_739[k] = f_18 * skk_410[k]
                   + f_3 * pc_y[k] * slk_590[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, pb_z, pc_x, pc_z, skl0_471, skk_375, skk_596, \
                         skl1_471, sli0_468, sli1_468, slk_591, \
                         slk_596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_17 * skk_596[k]
                   + f_10 * sli0_468[k]
                   - f_11 * sli1_468[k]
                   + f_3 * pc_x[k] * slk_596[k];

        t_741[k] = pb_z[k] * skl0_471[k]
                   - f_14 * pc_z[k] * skl1_471[k];

        t_742[k] = f_15 * skk_375[k]
                   + f_3 * pc_z[k] * slk_591[k];
    }

#pragma omp simd aligned(t_743, t_744, t_745, pb_z, pc_z, skl0_473, skl0_474, skl0_475, \
                         skk_376, skk_377, skk_378, skl1_473, skl1_474, \
                         skl1_475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_743[k] = pb_z[k] * skl0_473[k]
                   + f_16 * skk_376[k]
                   - f_14 * pc_z[k] * skl1_473[k];

        t_744[k] = pb_z[k] * skl0_474[k]
                   + f_17 * skk_377[k]
                   - f_14 * pc_z[k] * skl1_474[k];

        t_745[k] = pb_z[k] * skl0_475[k]
                   + f_18 * skk_378[k]
                   - f_14 * pc_z[k] * skl1_475[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pc_x, pc_y, skk_416, skk_603, skk_604, \
                         skk_605, sli0_475, sli1_475, slk_596, slk_603, slk_604, \
                         slk_605 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_18 * skk_416[k]
                   + f_3 * pc_y[k] * slk_596[k];

        t_747[k] = f_17 * skk_603[k]
                   + f_12 * sli0_475[k]
                   - f_13 * sli1_475[k]
                   + f_3 * pc_x[k] * slk_603[k];

        t_748[k] = f_17 * skk_604[k]
                   + f_3 * pc_x[k] * slk_604[k];

        t_749[k] = f_17 * skk_605[k]
                   + f_3 * pc_x[k] * slk_605[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, t_753, t_754, pc_x, skk_606, skk_607, skk_608, \
                         skk_609, skk_610, slk_606, slk_607, slk_608, slk_609, \
                         slk_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_17 * skk_606[k]
                   + f_3 * pc_x[k] * slk_606[k];

        t_751[k] = f_17 * skk_607[k]
                   + f_3 * pc_x[k] * slk_607[k];

        t_752[k] = f_17 * skk_608[k]
                   + f_3 * pc_x[k] * slk_608[k];

        t_753[k] = f_17 * skk_609[k]
                   + f_3 * pc_x[k] * slk_609[k];

        t_754[k] = f_17 * skk_610[k]
                   + f_3 * pc_x[k] * slk_610[k];
    }

#pragma omp simd aligned(t_755, t_756, t_757, pb_z, pc_x, pc_z, skl0_486, skk_388, skk_611, \
                         skl1_486, slk_604, slk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_755[k] = f_17 * skk_611[k]
                   + f_3 * pc_x[k] * slk_611[k];

        t_756[k] = pb_z[k] * skl0_486[k]
                   - f_14 * pc_z[k] * skl1_486[k];

        t_757[k] = f_15 * skk_388[k]
                   + f_3 * pc_z[k] * slk_604[k];
    }

#pragma omp simd aligned(t_758, t_759, t_760, pc_y, skk_426, skk_427, skk_428, sli0_471, \
                         sli0_472, sli0_473, sli1_471, sli1_472, sli1_473, slk_606, slk_607, \
                         slk_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_758[k] = f_18 * skk_426[k]
                   + f_4 * sli0_471[k]
                   - f_5 * sli1_471[k]
                   + f_3 * pc_y[k] * slk_606[k];

        t_759[k] = f_18 * skk_427[k]
                   + f_6 * sli0_472[k]
                   - f_7 * sli1_472[k]
                   + f_3 * pc_y[k] * slk_607[k];

        t_760[k] = f_18 * skk_428[k]
                   + f_8 * sli0_473[k]
                   - f_9 * sli1_473[k]
                   + f_3 * pc_y[k] * slk_608[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, pc_y, skk_429, skk_430, skk_431, sli0_474, \
                         sli0_475, sli1_474, sli1_475, slk_609, slk_610, \
                         slk_611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_18 * skk_429[k]
                   + f_10 * sli0_474[k]
                   - f_11 * sli1_474[k]
                   + f_3 * pc_y[k] * slk_609[k];

        t_762[k] = f_18 * skk_430[k]
                   + f_12 * sli0_475[k]
                   - f_13 * sli1_475[k]
                   + f_3 * pc_y[k] * slk_610[k];

        t_763[k] = f_18 * skk_431[k]
                   + f_3 * pc_y[k] * slk_611[k];
    }

#pragma omp simd aligned(t_764, t_765, t_766, pc_x, pc_y, pc_z, skk_395, skk_432, skk_612, \
                         sli0_475, sli0_476, sli1_475, sli1_476, slk_611, \
                         slk_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_764[k] = f_15 * skk_395[k]
                   + f_1 * sli0_475[k]
                   - f_2 * sli1_475[k]
                   + f_3 * pc_z[k] * slk_611[k];

        t_765[k] = f_17 * skk_612[k]
                   + f_1 * sli0_476[k]
                   - f_2 * sli1_476[k]
                   + f_3 * pc_x[k] * slk_612[k];

        t_766[k] = f_17 * skk_432[k]
                   + f_3 * pc_y[k] * slk_612[k];
    }

#pragma omp simd aligned(t_767, t_768, t_769, pc_x, pc_y, pc_z, skk_396, skk_434, skk_615, \
                         sli0_479, sli1_479, slk_612, slk_614, \
                         slk_615 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_767[k] = f_16 * skk_396[k]
                   + f_3 * pc_z[k] * slk_612[k];

        t_768[k] = f_17 * skk_615[k]
                   + f_4 * sli0_479[k]
                   - f_5 * sli1_479[k]
                   + f_3 * pc_x[k] * slk_615[k];

        t_769[k] = f_17 * skk_434[k]
                   + f_3 * pc_y[k] * slk_614[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, pc_x, pc_z, skk_399, skk_617, skk_618, sli0_481, \
                         sli0_482, sli1_481, sli1_482, slk_615, slk_617, \
                         slk_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_17 * skk_617[k]
                   + f_4 * sli0_481[k]
                   - f_5 * sli1_481[k]
                   + f_3 * pc_x[k] * slk_617[k];

        t_771[k] = f_17 * skk_618[k]
                   + f_6 * sli0_482[k]
                   - f_7 * sli1_482[k]
                   + f_3 * pc_x[k] * slk_618[k];

        t_772[k] = f_16 * skk_399[k]
                   + f_3 * pc_z[k] * slk_615[k];
    }

#pragma omp simd aligned(t_773, t_774, t_775, pc_x, pc_y, skk_437, skk_621, skk_622, sli0_485, \
                         sli0_486, sli1_485, sli1_486, slk_617, slk_621, \
                         slk_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_773[k] = f_17 * skk_437[k]
                   + f_3 * pc_y[k] * slk_617[k];

        t_774[k] = f_17 * skk_621[k]
                   + f_6 * sli0_485[k]
                   - f_7 * sli1_485[k]
                   + f_3 * pc_x[k] * slk_621[k];

        t_775[k] = f_17 * skk_622[k]
                   + f_8 * sli0_486[k]
                   - f_9 * sli1_486[k]
                   + f_3 * pc_x[k] * slk_622[k];
    }

#pragma omp simd aligned(t_776, t_777, t_778, pc_x, pc_y, pc_z, skk_402, skk_441, skk_624, \
                         sli0_488, sli1_488, slk_618, slk_621, \
                         slk_624 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_776[k] = f_16 * skk_402[k]
                   + f_3 * pc_z[k] * slk_618[k];

        t_777[k] = f_17 * skk_624[k]
                   + f_8 * sli0_488[k]
                   - f_9 * sli1_488[k]
                   + f_3 * pc_x[k] * slk_624[k];

        t_778[k] = f_17 * skk_441[k]
                   + f_3 * pc_y[k] * slk_621[k];
    }

#pragma omp simd aligned(t_779, t_780, t_781, pc_x, pc_z, skk_406, skk_626, skk_627, sli0_490, \
                         sli0_491, sli1_490, sli1_491, slk_622, slk_626, \
                         slk_627 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_779[k] = f_17 * skk_626[k]
                   + f_8 * sli0_490[k]
                   - f_9 * sli1_490[k]
                   + f_3 * pc_x[k] * slk_626[k];

        t_780[k] = f_17 * skk_627[k]
                   + f_10 * sli0_491[k]
                   - f_11 * sli1_491[k]
                   + f_3 * pc_x[k] * slk_627[k];

        t_781[k] = f_16 * skk_406[k]
                   + f_3 * pc_z[k] * slk_622[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *skl0_630 = buffer.data(skl0 + 630);
    const auto *skl0_633 = buffer.data(skl0 + 633);
    const auto *skl0_635 = buffer.data(skl0 + 635);
    const auto *skl0_636 = buffer.data(skl0 + 636);
    const auto *skl0_639 = buffer.data(skl0 + 639);
    const auto *skl0_640 = buffer.data(skl0 + 640);
    const auto *skl0_642 = buffer.data(skl0 + 642);
    const auto *skl0_644 = buffer.data(skl0 + 644);
    const auto *skl0_645 = buffer.data(skl0 + 645);
    const auto *skl0_647 = buffer.data(skl0 + 647);
    const auto *skl0_648 = buffer.data(skl0 + 648);
    const auto *skl0_650 = buffer.data(skl0 + 650);
    const auto *skl0_651 = buffer.data(skl0 + 651);
    const auto *skl0_653 = buffer.data(skl0 + 653);
    const auto *skl0_654 = buffer.data(skl0 + 654);
    const auto *skl0_655 = buffer.data(skl0 + 655);
    const auto *skl0_657 = buffer.data(skl0 + 657);

    const auto *skk_411 = buffer.data(skk + 411);
    const auto *skk_424 = buffer.data(skk + 424);
    const auto *skk_431 = buffer.data(skk + 431);
    const auto *skk_432 = buffer.data(skk + 432);
    const auto *skk_435 = buffer.data(skk + 435);
    const auto *skk_438 = buffer.data(skk + 438);
    const auto *skk_442 = buffer.data(skk + 442);
    const auto *skk_446 = buffer.data(skk + 446);
    const auto *skk_447 = buffer.data(skk + 447);
    const auto *skk_452 = buffer.data(skk + 452);
    const auto *skk_460 = buffer.data(skk + 460);
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
    const auto *skk_498 = buffer.data(skk + 498);
    const auto *skk_499 = buffer.data(skk + 499);
    const auto *skk_500 = buffer.data(skk + 500);
    const auto *skk_501 = buffer.data(skk + 501);
    const auto *skk_502 = buffer.data(skk + 502);
    const auto *skk_503 = buffer.data(skk + 503);
    const auto *skk_504 = buffer.data(skk + 504);
    const auto *skk_505 = buffer.data(skk + 505);
    const auto *skk_506 = buffer.data(skk + 506);
    const auto *skk_507 = buffer.data(skk + 507);
    const auto *skk_509 = buffer.data(skk + 509);
    const auto *skk_510 = buffer.data(skk + 510);
    const auto *skk_512 = buffer.data(skk + 512);
    const auto *skk_513 = buffer.data(skk + 513);
    const auto *skk_514 = buffer.data(skk + 514);
    const auto *skk_516 = buffer.data(skk + 516);
    const auto *skk_517 = buffer.data(skk + 517);
    const auto *skk_518 = buffer.data(skk + 518);
    const auto *skk_519 = buffer.data(skk + 519);
    const auto *skk_521 = buffer.data(skk + 521);
    const auto *skk_522 = buffer.data(skk + 522);
    const auto *skk_523 = buffer.data(skk + 523);
    const auto *skk_524 = buffer.data(skk + 524);
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
    const auto *skk_712 = buffer.data(skk + 712);
    const auto *skk_713 = buffer.data(skk + 713);
    const auto *skk_714 = buffer.data(skk + 714);
    const auto *skk_715 = buffer.data(skk + 715);
    const auto *skk_716 = buffer.data(skk + 716);
    const auto *skk_717 = buffer.data(skk + 717);

    const auto *skl1_630 = buffer.data(skl1 + 630);
    const auto *skl1_633 = buffer.data(skl1 + 633);
    const auto *skl1_635 = buffer.data(skl1 + 635);
    const auto *skl1_636 = buffer.data(skl1 + 636);
    const auto *skl1_639 = buffer.data(skl1 + 639);
    const auto *skl1_640 = buffer.data(skl1 + 640);
    const auto *skl1_642 = buffer.data(skl1 + 642);
    const auto *skl1_644 = buffer.data(skl1 + 644);
    const auto *skl1_645 = buffer.data(skl1 + 645);
    const auto *skl1_647 = buffer.data(skl1 + 647);
    const auto *skl1_648 = buffer.data(skl1 + 648);
    const auto *skl1_650 = buffer.data(skl1 + 650);
    const auto *skl1_651 = buffer.data(skl1 + 651);
    const auto *skl1_653 = buffer.data(skl1 + 653);
    const auto *skl1_654 = buffer.data(skl1 + 654);
    const auto *skl1_655 = buffer.data(skl1 + 655);
    const auto *skl1_657 = buffer.data(skl1 + 657);

    const auto *sli0_493 = buffer.data(sli0 + 493);
    const auto *sli0_494 = buffer.data(sli0 + 494);
    const auto *sli0_496 = buffer.data(sli0 + 496);
    const auto *sli0_497 = buffer.data(sli0 + 497);
    const auto *sli0_499 = buffer.data(sli0 + 499);
    const auto *sli0_500 = buffer.data(sli0 + 500);
    const auto *sli0_501 = buffer.data(sli0 + 501);
    const auto *sli0_502 = buffer.data(sli0 + 502);
    const auto *sli0_503 = buffer.data(sli0 + 503);
    const auto *sli0_504 = buffer.data(sli0 + 504);
    const auto *sli0_507 = buffer.data(sli0 + 507);
    const auto *sli0_509 = buffer.data(sli0 + 509);
    const auto *sli0_510 = buffer.data(sli0 + 510);
    const auto *sli0_513 = buffer.data(sli0 + 513);
    const auto *sli0_514 = buffer.data(sli0 + 514);
    const auto *sli0_516 = buffer.data(sli0 + 516);
    const auto *sli0_518 = buffer.data(sli0 + 518);
    const auto *sli0_519 = buffer.data(sli0 + 519);
    const auto *sli0_521 = buffer.data(sli0 + 521);
    const auto *sli0_522 = buffer.data(sli0 + 522);
    const auto *sli0_524 = buffer.data(sli0 + 524);
    const auto *sli0_525 = buffer.data(sli0 + 525);
    const auto *sli0_527 = buffer.data(sli0 + 527);
    const auto *sli0_528 = buffer.data(sli0 + 528);
    const auto *sli0_529 = buffer.data(sli0 + 529);
    const auto *sli0_530 = buffer.data(sli0 + 530);
    const auto *sli0_531 = buffer.data(sli0 + 531);

    const auto *sli1_493 = buffer.data(sli1 + 493);
    const auto *sli1_494 = buffer.data(sli1 + 494);
    const auto *sli1_496 = buffer.data(sli1 + 496);
    const auto *sli1_497 = buffer.data(sli1 + 497);
    const auto *sli1_499 = buffer.data(sli1 + 499);
    const auto *sli1_500 = buffer.data(sli1 + 500);
    const auto *sli1_501 = buffer.data(sli1 + 501);
    const auto *sli1_502 = buffer.data(sli1 + 502);
    const auto *sli1_503 = buffer.data(sli1 + 503);
    const auto *sli1_504 = buffer.data(sli1 + 504);
    const auto *sli1_507 = buffer.data(sli1 + 507);
    const auto *sli1_509 = buffer.data(sli1 + 509);
    const auto *sli1_510 = buffer.data(sli1 + 510);
    const auto *sli1_513 = buffer.data(sli1 + 513);
    const auto *sli1_514 = buffer.data(sli1 + 514);
    const auto *sli1_516 = buffer.data(sli1 + 516);
    const auto *sli1_518 = buffer.data(sli1 + 518);
    const auto *sli1_519 = buffer.data(sli1 + 519);
    const auto *sli1_521 = buffer.data(sli1 + 521);
    const auto *sli1_522 = buffer.data(sli1 + 522);
    const auto *sli1_524 = buffer.data(sli1 + 524);
    const auto *sli1_525 = buffer.data(sli1 + 525);
    const auto *sli1_527 = buffer.data(sli1 + 527);
    const auto *sli1_528 = buffer.data(sli1 + 528);
    const auto *sli1_529 = buffer.data(sli1 + 529);
    const auto *sli1_530 = buffer.data(sli1 + 530);
    const auto *sli1_531 = buffer.data(sli1 + 531);

    const auto *slk_626 = buffer.data(slk + 626);
    const auto *slk_627 = buffer.data(slk + 627);
    const auto *slk_629 = buffer.data(slk + 629);
    const auto *slk_630 = buffer.data(slk + 630);
    const auto *slk_632 = buffer.data(slk + 632);
    const auto *slk_633 = buffer.data(slk + 633);
    const auto *slk_635 = buffer.data(slk + 635);
    const auto *slk_636 = buffer.data(slk + 636);
    const auto *slk_637 = buffer.data(slk + 637);
    const auto *slk_639 = buffer.data(slk + 639);
    const auto *slk_640 = buffer.data(slk + 640);
    const auto *slk_641 = buffer.data(slk + 641);
    const auto *slk_642 = buffer.data(slk + 642);
    const auto *slk_643 = buffer.data(slk + 643);
    const auto *slk_644 = buffer.data(slk + 644);
    const auto *slk_645 = buffer.data(slk + 645);
    const auto *slk_646 = buffer.data(slk + 646);
    const auto *slk_647 = buffer.data(slk + 647);
    const auto *slk_648 = buffer.data(slk + 648);
    const auto *slk_650 = buffer.data(slk + 650);
    const auto *slk_651 = buffer.data(slk + 651);
    const auto *slk_653 = buffer.data(slk + 653);
    const auto *slk_654 = buffer.data(slk + 654);
    const auto *slk_657 = buffer.data(slk + 657);
    const auto *slk_658 = buffer.data(slk + 658);
    const auto *slk_660 = buffer.data(slk + 660);
    const auto *slk_662 = buffer.data(slk + 662);
    const auto *slk_663 = buffer.data(slk + 663);
    const auto *slk_665 = buffer.data(slk + 665);
    const auto *slk_666 = buffer.data(slk + 666);
    const auto *slk_668 = buffer.data(slk + 668);
    const auto *slk_669 = buffer.data(slk + 669);
    const auto *slk_671 = buffer.data(slk + 671);
    const auto *slk_672 = buffer.data(slk + 672);
    const auto *slk_673 = buffer.data(slk + 673);
    const auto *slk_675 = buffer.data(slk + 675);
    const auto *slk_676 = buffer.data(slk + 676);
    const auto *slk_677 = buffer.data(slk + 677);
    const auto *slk_678 = buffer.data(slk + 678);
    const auto *slk_679 = buffer.data(slk + 679);
    const auto *slk_680 = buffer.data(slk + 680);
    const auto *slk_681 = buffer.data(slk + 681);
    const auto *slk_682 = buffer.data(slk + 682);
    const auto *slk_683 = buffer.data(slk + 683);
    const auto *slk_684 = buffer.data(slk + 684);
    const auto *slk_686 = buffer.data(slk + 686);
    const auto *slk_687 = buffer.data(slk + 687);
    const auto *slk_689 = buffer.data(slk + 689);
    const auto *slk_690 = buffer.data(slk + 690);
    const auto *slk_693 = buffer.data(slk + 693);
    const auto *slk_694 = buffer.data(slk + 694);
    const auto *slk_698 = buffer.data(slk + 698);
    const auto *slk_699 = buffer.data(slk + 699);
    const auto *slk_704 = buffer.data(slk + 704);
    const auto *slk_712 = buffer.data(slk + 712);
    const auto *slk_713 = buffer.data(slk + 713);
    const auto *slk_714 = buffer.data(slk + 714);
    const auto *slk_715 = buffer.data(slk + 715);
    const auto *slk_716 = buffer.data(slk + 716);
    const auto *slk_717 = buffer.data(slk + 717);

#pragma omp simd aligned(t_782, t_783, t_784, pc_x, pc_y, skk_446, skk_629, skk_630, sli0_493, \
                         sli0_494, sli1_493, sli1_494, slk_626, slk_629, \
                         slk_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_782[k] = f_17 * skk_629[k]
                   + f_10 * sli0_493[k]
                   - f_11 * sli1_493[k]
                   + f_3 * pc_x[k] * slk_629[k];

        t_783[k] = f_17 * skk_630[k]
                   + f_10 * sli0_494[k]
                   - f_11 * sli1_494[k]
                   + f_3 * pc_x[k] * slk_630[k];

        t_784[k] = f_17 * skk_446[k]
                   + f_3 * pc_y[k] * slk_626[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, pc_x, pc_z, skk_411, skk_632, skk_633, sli0_496, \
                         sli0_497, sli1_496, sli1_497, slk_627, slk_632, \
                         slk_633 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_17 * skk_632[k]
                   + f_10 * sli0_496[k]
                   - f_11 * sli1_496[k]
                   + f_3 * pc_x[k] * slk_632[k];

        t_786[k] = f_17 * skk_633[k]
                   + f_12 * sli0_497[k]
                   - f_13 * sli1_497[k]
                   + f_3 * pc_x[k] * slk_633[k];

        t_787[k] = f_16 * skk_411[k]
                   + f_3 * pc_z[k] * slk_627[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, pc_x, skk_635, skk_636, skk_637, sli0_499, \
                         sli0_500, sli0_501, sli1_499, sli1_500, sli1_501, slk_635, slk_636, \
                         slk_637 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_17 * skk_635[k]
                   + f_12 * sli0_499[k]
                   - f_13 * sli1_499[k]
                   + f_3 * pc_x[k] * slk_635[k];

        t_789[k] = f_17 * skk_636[k]
                   + f_12 * sli0_500[k]
                   - f_13 * sli1_500[k]
                   + f_3 * pc_x[k] * slk_636[k];

        t_790[k] = f_17 * skk_637[k]
                   + f_12 * sli0_501[k]
                   - f_13 * sli1_501[k]
                   + f_3 * pc_x[k] * slk_637[k];
    }

#pragma omp simd aligned(t_791, t_792, t_793, t_794, pc_x, pc_y, skk_452, skk_639, skk_640, \
                         skk_641, sli0_503, sli1_503, slk_632, slk_639, slk_640, \
                         slk_641 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_791[k] = f_17 * skk_452[k]
                   + f_3 * pc_y[k] * slk_632[k];

        t_792[k] = f_17 * skk_639[k]
                   + f_12 * sli0_503[k]
                   - f_13 * sli1_503[k]
                   + f_3 * pc_x[k] * slk_639[k];

        t_793[k] = f_17 * skk_640[k]
                   + f_3 * pc_x[k] * slk_640[k];

        t_794[k] = f_17 * skk_641[k]
                   + f_3 * pc_x[k] * slk_641[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, t_798, t_799, pc_x, skk_642, skk_643, skk_644, \
                         skk_645, skk_646, slk_642, slk_643, slk_644, slk_645, \
                         slk_646 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_17 * skk_642[k]
                   + f_3 * pc_x[k] * slk_642[k];

        t_796[k] = f_17 * skk_643[k]
                   + f_3 * pc_x[k] * slk_643[k];

        t_797[k] = f_17 * skk_644[k]
                   + f_3 * pc_x[k] * slk_644[k];

        t_798[k] = f_17 * skk_645[k]
                   + f_3 * pc_x[k] * slk_645[k];

        t_799[k] = f_17 * skk_646[k]
                   + f_3 * pc_x[k] * slk_646[k];
    }

#pragma omp simd aligned(t_800, t_801, t_802, pc_x, pc_y, pc_z, skk_424, skk_460, skk_647, \
                         sli0_497, sli1_497, slk_640, slk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_800[k] = f_17 * skk_647[k]
                   + f_3 * pc_x[k] * slk_647[k];

        t_801[k] = f_17 * skk_460[k]
                   + f_1 * sli0_497[k]
                   - f_2 * sli1_497[k]
                   + f_3 * pc_y[k] * slk_640[k];

        t_802[k] = f_16 * skk_424[k]
                   + f_3 * pc_z[k] * slk_640[k];
    }

#pragma omp simd aligned(t_803, t_804, t_805, pc_y, skk_462, skk_463, skk_464, sli0_499, \
                         sli0_500, sli0_501, sli1_499, sli1_500, sli1_501, slk_642, slk_643, \
                         slk_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_803[k] = f_17 * skk_462[k]
                   + f_4 * sli0_499[k]
                   - f_5 * sli1_499[k]
                   + f_3 * pc_y[k] * slk_642[k];

        t_804[k] = f_17 * skk_463[k]
                   + f_6 * sli0_500[k]
                   - f_7 * sli1_500[k]
                   + f_3 * pc_y[k] * slk_643[k];

        t_805[k] = f_17 * skk_464[k]
                   + f_8 * sli0_501[k]
                   - f_9 * sli1_501[k]
                   + f_3 * pc_y[k] * slk_644[k];
    }

#pragma omp simd aligned(t_806, t_807, t_808, pc_y, skk_465, skk_466, skk_467, sli0_502, \
                         sli0_503, sli1_502, sli1_503, slk_645, slk_646, \
                         slk_647 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_806[k] = f_17 * skk_465[k]
                   + f_10 * sli0_502[k]
                   - f_11 * sli1_502[k]
                   + f_3 * pc_y[k] * slk_645[k];

        t_807[k] = f_17 * skk_466[k]
                   + f_12 * sli0_503[k]
                   - f_13 * sli1_503[k]
                   + f_3 * pc_y[k] * slk_646[k];

        t_808[k] = f_17 * skk_467[k]
                   + f_3 * pc_y[k] * slk_647[k];
    }

#pragma omp simd aligned(t_809, t_810, t_811, pc_x, pc_y, pc_z, skk_431, skk_468, skk_648, \
                         sli0_503, sli0_504, sli1_503, sli1_504, slk_647, \
                         slk_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_809[k] = f_16 * skk_431[k]
                   + f_1 * sli0_503[k]
                   - f_2 * sli1_503[k]
                   + f_3 * pc_z[k] * slk_647[k];

        t_810[k] = f_17 * skk_648[k]
                   + f_1 * sli0_504[k]
                   - f_2 * sli1_504[k]
                   + f_3 * pc_x[k] * slk_648[k];

        t_811[k] = f_16 * skk_468[k]
                   + f_3 * pc_y[k] * slk_648[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_x, pc_y, pc_z, skk_432, skk_470, skk_651, \
                         sli0_507, sli1_507, slk_648, slk_650, \
                         slk_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_17 * skk_432[k]
                   + f_3 * pc_z[k] * slk_648[k];

        t_813[k] = f_17 * skk_651[k]
                   + f_4 * sli0_507[k]
                   - f_5 * sli1_507[k]
                   + f_3 * pc_x[k] * slk_651[k];

        t_814[k] = f_16 * skk_470[k]
                   + f_3 * pc_y[k] * slk_650[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, pc_x, pc_z, skk_435, skk_653, skk_654, sli0_509, \
                         sli0_510, sli1_509, sli1_510, slk_651, slk_653, \
                         slk_654 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_17 * skk_653[k]
                   + f_4 * sli0_509[k]
                   - f_5 * sli1_509[k]
                   + f_3 * pc_x[k] * slk_653[k];

        t_816[k] = f_17 * skk_654[k]
                   + f_6 * sli0_510[k]
                   - f_7 * sli1_510[k]
                   + f_3 * pc_x[k] * slk_654[k];

        t_817[k] = f_17 * skk_435[k]
                   + f_3 * pc_z[k] * slk_651[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, pc_x, pc_y, skk_473, skk_657, skk_658, sli0_513, \
                         sli0_514, sli1_513, sli1_514, slk_653, slk_657, \
                         slk_658 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_16 * skk_473[k]
                   + f_3 * pc_y[k] * slk_653[k];

        t_819[k] = f_17 * skk_657[k]
                   + f_6 * sli0_513[k]
                   - f_7 * sli1_513[k]
                   + f_3 * pc_x[k] * slk_657[k];

        t_820[k] = f_17 * skk_658[k]
                   + f_8 * sli0_514[k]
                   - f_9 * sli1_514[k]
                   + f_3 * pc_x[k] * slk_658[k];
    }

#pragma omp simd aligned(t_821, t_822, t_823, pc_x, pc_y, pc_z, skk_438, skk_477, skk_660, \
                         sli0_516, sli1_516, slk_654, slk_657, \
                         slk_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_821[k] = f_17 * skk_438[k]
                   + f_3 * pc_z[k] * slk_654[k];

        t_822[k] = f_17 * skk_660[k]
                   + f_8 * sli0_516[k]
                   - f_9 * sli1_516[k]
                   + f_3 * pc_x[k] * slk_660[k];

        t_823[k] = f_16 * skk_477[k]
                   + f_3 * pc_y[k] * slk_657[k];
    }

#pragma omp simd aligned(t_824, t_825, t_826, pc_x, pc_z, skk_442, skk_662, skk_663, sli0_518, \
                         sli0_519, sli1_518, sli1_519, slk_658, slk_662, \
                         slk_663 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_824[k] = f_17 * skk_662[k]
                   + f_8 * sli0_518[k]
                   - f_9 * sli1_518[k]
                   + f_3 * pc_x[k] * slk_662[k];

        t_825[k] = f_17 * skk_663[k]
                   + f_10 * sli0_519[k]
                   - f_11 * sli1_519[k]
                   + f_3 * pc_x[k] * slk_663[k];

        t_826[k] = f_17 * skk_442[k]
                   + f_3 * pc_z[k] * slk_658[k];
    }

#pragma omp simd aligned(t_827, t_828, t_829, pc_x, pc_y, skk_482, skk_665, skk_666, sli0_521, \
                         sli0_522, sli1_521, sli1_522, slk_662, slk_665, \
                         slk_666 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_827[k] = f_17 * skk_665[k]
                   + f_10 * sli0_521[k]
                   - f_11 * sli1_521[k]
                   + f_3 * pc_x[k] * slk_665[k];

        t_828[k] = f_17 * skk_666[k]
                   + f_10 * sli0_522[k]
                   - f_11 * sli1_522[k]
                   + f_3 * pc_x[k] * slk_666[k];

        t_829[k] = f_16 * skk_482[k]
                   + f_3 * pc_y[k] * slk_662[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, pc_x, pc_z, skk_447, skk_668, skk_669, sli0_524, \
                         sli0_525, sli1_524, sli1_525, slk_663, slk_668, \
                         slk_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_17 * skk_668[k]
                   + f_10 * sli0_524[k]
                   - f_11 * sli1_524[k]
                   + f_3 * pc_x[k] * slk_668[k];

        t_831[k] = f_17 * skk_669[k]
                   + f_12 * sli0_525[k]
                   - f_13 * sli1_525[k]
                   + f_3 * pc_x[k] * slk_669[k];

        t_832[k] = f_17 * skk_447[k]
                   + f_3 * pc_z[k] * slk_663[k];
    }

#pragma omp simd aligned(t_833, t_834, t_835, pc_x, skk_671, skk_672, skk_673, sli0_527, \
                         sli0_528, sli0_529, sli1_527, sli1_528, sli1_529, slk_671, slk_672, \
                         slk_673 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_833[k] = f_17 * skk_671[k]
                   + f_12 * sli0_527[k]
                   - f_13 * sli1_527[k]
                   + f_3 * pc_x[k] * slk_671[k];

        t_834[k] = f_17 * skk_672[k]
                   + f_12 * sli0_528[k]
                   - f_13 * sli1_528[k]
                   + f_3 * pc_x[k] * slk_672[k];

        t_835[k] = f_17 * skk_673[k]
                   + f_12 * sli0_529[k]
                   - f_13 * sli1_529[k]
                   + f_3 * pc_x[k] * slk_673[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_x, pc_y, skk_488, skk_675, skk_676, \
                         skk_677, sli0_531, sli1_531, slk_668, slk_675, slk_676, \
                         slk_677 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_16 * skk_488[k]
                   + f_3 * pc_y[k] * slk_668[k];

        t_837[k] = f_17 * skk_675[k]
                   + f_12 * sli0_531[k]
                   - f_13 * sli1_531[k]
                   + f_3 * pc_x[k] * slk_675[k];

        t_838[k] = f_17 * skk_676[k]
                   + f_3 * pc_x[k] * slk_676[k];

        t_839[k] = f_17 * skk_677[k]
                   + f_3 * pc_x[k] * slk_677[k];
    }

#pragma omp simd aligned(t_840, t_841, t_842, t_843, t_844, pc_x, skk_678, skk_679, skk_680, \
                         skk_681, skk_682, slk_678, slk_679, slk_680, slk_681, \
                         slk_682 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_840[k] = f_17 * skk_678[k]
                   + f_3 * pc_x[k] * slk_678[k];

        t_841[k] = f_17 * skk_679[k]
                   + f_3 * pc_x[k] * slk_679[k];

        t_842[k] = f_17 * skk_680[k]
                   + f_3 * pc_x[k] * slk_680[k];

        t_843[k] = f_17 * skk_681[k]
                   + f_3 * pc_x[k] * slk_681[k];

        t_844[k] = f_17 * skk_682[k]
                   + f_3 * pc_x[k] * slk_682[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, pc_x, pc_y, pc_z, skk_460, skk_496, skk_683, \
                         sli0_525, sli1_525, slk_676, slk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = f_17 * skk_683[k]
                   + f_3 * pc_x[k] * slk_683[k];

        t_846[k] = f_16 * skk_496[k]
                   + f_1 * sli0_525[k]
                   - f_2 * sli1_525[k]
                   + f_3 * pc_y[k] * slk_676[k];

        t_847[k] = f_17 * skk_460[k]
                   + f_3 * pc_z[k] * slk_676[k];
    }

#pragma omp simd aligned(t_848, t_849, t_850, pc_y, skk_498, skk_499, skk_500, sli0_527, \
                         sli0_528, sli0_529, sli1_527, sli1_528, sli1_529, slk_678, slk_679, \
                         slk_680 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_848[k] = f_16 * skk_498[k]
                   + f_4 * sli0_527[k]
                   - f_5 * sli1_527[k]
                   + f_3 * pc_y[k] * slk_678[k];

        t_849[k] = f_16 * skk_499[k]
                   + f_6 * sli0_528[k]
                   - f_7 * sli1_528[k]
                   + f_3 * pc_y[k] * slk_679[k];

        t_850[k] = f_16 * skk_500[k]
                   + f_8 * sli0_529[k]
                   - f_9 * sli1_529[k]
                   + f_3 * pc_y[k] * slk_680[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, pc_y, skk_501, skk_502, skk_503, sli0_530, \
                         sli0_531, sli1_530, sli1_531, slk_681, slk_682, \
                         slk_683 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_16 * skk_501[k]
                   + f_10 * sli0_530[k]
                   - f_11 * sli1_530[k]
                   + f_3 * pc_y[k] * slk_681[k];

        t_852[k] = f_16 * skk_502[k]
                   + f_12 * sli0_531[k]
                   - f_13 * sli1_531[k]
                   + f_3 * pc_y[k] * slk_682[k];

        t_853[k] = f_16 * skk_503[k]
                   + f_3 * pc_y[k] * slk_683[k];
    }

#pragma omp simd aligned(t_854, t_855, t_856, t_857, pb_y, pc_y, pc_z, skl0_630, skk_467, \
                         skk_468, skk_504, skl1_630, sli0_531, sli1_531, slk_683, \
                         slk_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_854[k] = f_17 * skk_467[k]
                   + f_1 * sli0_531[k]
                   - f_2 * sli1_531[k]
                   + f_3 * pc_z[k] * slk_683[k];

        t_855[k] = pb_y[k] * skl0_630[k]
                   - f_14 * pc_y[k] * skl1_630[k];

        t_856[k] = f_15 * skk_504[k]
                   + f_3 * pc_y[k] * slk_684[k];

        t_857[k] = f_18 * skk_468[k]
                   + f_3 * pc_z[k] * slk_684[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, t_861, pb_y, pc_y, skl0_633, skl0_635, skl0_636, \
                         skk_505, skk_506, skk_507, skl1_633, skl1_635, skl1_636, \
                         slk_686 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = pb_y[k] * skl0_633[k]
                   + f_16 * skk_505[k]
                   - f_14 * pc_y[k] * skl1_633[k];

        t_859[k] = f_15 * skk_506[k]
                   + f_3 * pc_y[k] * slk_686[k];

        t_860[k] = pb_y[k] * skl0_635[k]
                   - f_14 * pc_y[k] * skl1_635[k];

        t_861[k] = pb_y[k] * skl0_636[k]
                   + f_17 * skk_507[k]
                   - f_14 * pc_y[k] * skl1_636[k];
    }

#pragma omp simd aligned(t_862, t_863, t_864, t_865, pb_y, pc_y, pc_z, skl0_639, skl0_640, \
                         skk_471, skk_509, skk_510, skl1_639, skl1_640, slk_687, \
                         slk_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_862[k] = f_18 * skk_471[k]
                   + f_3 * pc_z[k] * slk_687[k];

        t_863[k] = f_15 * skk_509[k]
                   + f_3 * pc_y[k] * slk_689[k];

        t_864[k] = pb_y[k] * skl0_639[k]
                   - f_14 * pc_y[k] * skl1_639[k];

        t_865[k] = pb_y[k] * skl0_640[k]
                   + f_18 * skk_510[k]
                   - f_14 * pc_y[k] * skl1_640[k];
    }

#pragma omp simd aligned(t_866, t_867, t_868, t_869, pb_y, pc_y, pc_z, skl0_642, skl0_644, \
                         skk_474, skk_512, skk_513, skl1_642, skl1_644, slk_690, \
                         slk_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_866[k] = f_18 * skk_474[k]
                   + f_3 * pc_z[k] * slk_690[k];

        t_867[k] = pb_y[k] * skl0_642[k]
                   + f_16 * skk_512[k]
                   - f_14 * pc_y[k] * skl1_642[k];

        t_868[k] = f_15 * skk_513[k]
                   + f_3 * pc_y[k] * slk_693[k];

        t_869[k] = pb_y[k] * skl0_644[k]
                   - f_14 * pc_y[k] * skl1_644[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, pb_y, pc_y, pc_z, skl0_645, skl0_647, skk_478, \
                         skk_514, skk_516, skl1_645, skl1_647, \
                         slk_694 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pb_y[k] * skl0_645[k]
                   + f_19 * skk_514[k]
                   - f_14 * pc_y[k] * skl1_645[k];

        t_871[k] = f_18 * skk_478[k]
                   + f_3 * pc_z[k] * slk_694[k];

        t_872[k] = pb_y[k] * skl0_647[k]
                   + f_17 * skk_516[k]
                   - f_14 * pc_y[k] * skl1_647[k];
    }

#pragma omp simd aligned(t_873, t_874, t_875, t_876, pb_y, pc_y, skl0_648, skl0_650, skl0_651, \
                         skk_517, skk_518, skk_519, skl1_648, skl1_650, skl1_651, \
                         slk_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_873[k] = pb_y[k] * skl0_648[k]
                   + f_16 * skk_517[k]
                   - f_14 * pc_y[k] * skl1_648[k];

        t_874[k] = f_15 * skk_518[k]
                   + f_3 * pc_y[k] * slk_698[k];

        t_875[k] = pb_y[k] * skl0_650[k]
                   - f_14 * pc_y[k] * skl1_650[k];

        t_876[k] = pb_y[k] * skl0_651[k]
                   + f_20 * skk_519[k]
                   - f_14 * pc_y[k] * skl1_651[k];
    }

#pragma omp simd aligned(t_877, t_878, t_879, pb_y, pc_y, pc_z, skl0_653, skl0_654, skk_483, \
                         skk_521, skk_522, skl1_653, skl1_654, \
                         slk_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_877[k] = f_18 * skk_483[k]
                   + f_3 * pc_z[k] * slk_699[k];

        t_878[k] = pb_y[k] * skl0_653[k]
                   + f_18 * skk_521[k]
                   - f_14 * pc_y[k] * skl1_653[k];

        t_879[k] = pb_y[k] * skl0_654[k]
                   + f_17 * skk_522[k]
                   - f_14 * pc_y[k] * skl1_654[k];
    }

#pragma omp simd aligned(t_880, t_881, t_882, t_883, pb_y, pc_x, pc_y, skl0_655, skl0_657, \
                         skk_523, skk_524, skk_712, skl1_655, skl1_657, slk_704, \
                         slk_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_880[k] = pb_y[k] * skl0_655[k]
                   + f_16 * skk_523[k]
                   - f_14 * pc_y[k] * skl1_655[k];

        t_881[k] = f_15 * skk_524[k]
                   + f_3 * pc_y[k] * slk_704[k];

        t_882[k] = pb_y[k] * skl0_657[k]
                   - f_14 * pc_y[k] * skl1_657[k];

        t_883[k] = f_17 * skk_712[k]
                   + f_3 * pc_x[k] * slk_712[k];
    }

#pragma omp simd aligned(t_884, t_885, t_886, t_887, t_888, pc_x, skk_713, skk_714, skk_715, \
                         skk_716, skk_717, slk_713, slk_714, slk_715, slk_716, \
                         slk_717 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_884[k] = f_17 * skk_713[k]
                   + f_3 * pc_x[k] * slk_713[k];

        t_885[k] = f_17 * skk_714[k]
                   + f_3 * pc_x[k] * slk_714[k];

        t_886[k] = f_17 * skk_715[k]
                   + f_3 * pc_x[k] * slk_715[k];

        t_887[k] = f_17 * skk_716[k]
                   + f_3 * pc_x[k] * slk_716[k];

        t_888[k] = f_17 * skk_717[k]
                   + f_3 * pc_x[k] * slk_717[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_674 = buffer.data(skl0 + 674);
    const auto *skl0_675 = buffer.data(skl0 + 675);
    const auto *skl0_678 = buffer.data(skl0 + 678);
    const auto *skl0_681 = buffer.data(skl0 + 681);

    const auto *skk_496 = buffer.data(skk + 496);
    const auto *skk_504 = buffer.data(skk + 504);
    const auto *skk_507 = buffer.data(skk + 507);
    const auto *skk_510 = buffer.data(skk + 510);
    const auto *skk_514 = buffer.data(skk + 514);
    const auto *skk_519 = buffer.data(skk + 519);
    const auto *skk_532 = buffer.data(skk + 532);
    const auto *skk_534 = buffer.data(skk + 534);
    const auto *skk_535 = buffer.data(skk + 535);
    const auto *skk_536 = buffer.data(skk + 536);
    const auto *skk_537 = buffer.data(skk + 537);
    const auto *skk_538 = buffer.data(skk + 538);
    const auto *skk_539 = buffer.data(skk + 539);
    const auto *skk_540 = buffer.data(skk + 540);
    const auto *skk_542 = buffer.data(skk + 542);
    const auto *skk_543 = buffer.data(skk + 543);
    const auto *skk_545 = buffer.data(skk + 545);
    const auto *skk_549 = buffer.data(skk + 549);
    const auto *skk_554 = buffer.data(skk + 554);
    const auto *skk_560 = buffer.data(skk + 560);
    const auto *skk_568 = buffer.data(skk + 568);
    const auto *skk_570 = buffer.data(skk + 570);
    const auto *skk_571 = buffer.data(skk + 571);
    const auto *skk_572 = buffer.data(skk + 572);
    const auto *skk_573 = buffer.data(skk + 573);
    const auto *skk_574 = buffer.data(skk + 574);
    const auto *skk_575 = buffer.data(skk + 575);
    const auto *skk_576 = buffer.data(skk + 576);
    const auto *skk_578 = buffer.data(skk + 578);
    const auto *skk_718 = buffer.data(skk + 718);
    const auto *skk_719 = buffer.data(skk + 719);
    const auto *skk_720 = buffer.data(skk + 720);
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
    const auto *skk_759 = buffer.data(skk + 759);
    const auto *skk_761 = buffer.data(skk + 761);
    const auto *skk_762 = buffer.data(skk + 762);
    const auto *skk_765 = buffer.data(skk + 765);
    const auto *skk_766 = buffer.data(skk + 766);
    const auto *skk_768 = buffer.data(skk + 768);
    const auto *skk_770 = buffer.data(skk + 770);
    const auto *skk_771 = buffer.data(skk + 771);
    const auto *skk_773 = buffer.data(skk + 773);
    const auto *skk_774 = buffer.data(skk + 774);
    const auto *skk_776 = buffer.data(skk + 776);
    const auto *skk_777 = buffer.data(skk + 777);
    const auto *skk_779 = buffer.data(skk + 779);
    const auto *skk_780 = buffer.data(skk + 780);
    const auto *skk_781 = buffer.data(skk + 781);
    const auto *skk_783 = buffer.data(skk + 783);
    const auto *skk_784 = buffer.data(skk + 784);
    const auto *skk_785 = buffer.data(skk + 785);
    const auto *skk_786 = buffer.data(skk + 786);
    const auto *skk_787 = buffer.data(skk + 787);
    const auto *skk_788 = buffer.data(skk + 788);
    const auto *skk_789 = buffer.data(skk + 789);
    const auto *skk_790 = buffer.data(skk + 790);
    const auto *skk_791 = buffer.data(skk + 791);
    const auto *skk_797 = buffer.data(skk + 797);

    const auto *skl1_674 = buffer.data(skl1 + 674);
    const auto *skl1_675 = buffer.data(skl1 + 675);
    const auto *skl1_678 = buffer.data(skl1 + 678);
    const auto *skl1_681 = buffer.data(skl1 + 681);

    const auto *sli0_553 = buffer.data(sli0 + 553);
    const auto *sli0_555 = buffer.data(sli0 + 555);
    const auto *sli0_556 = buffer.data(sli0 + 556);
    const auto *sli0_557 = buffer.data(sli0 + 557);
    const auto *sli0_558 = buffer.data(sli0 + 558);
    const auto *sli0_559 = buffer.data(sli0 + 559);
    const auto *sli0_560 = buffer.data(sli0 + 560);
    const auto *sli0_563 = buffer.data(sli0 + 563);
    const auto *sli0_565 = buffer.data(sli0 + 565);
    const auto *sli0_566 = buffer.data(sli0 + 566);
    const auto *sli0_569 = buffer.data(sli0 + 569);
    const auto *sli0_570 = buffer.data(sli0 + 570);
    const auto *sli0_572 = buffer.data(sli0 + 572);
    const auto *sli0_574 = buffer.data(sli0 + 574);
    const auto *sli0_575 = buffer.data(sli0 + 575);
    const auto *sli0_577 = buffer.data(sli0 + 577);
    const auto *sli0_578 = buffer.data(sli0 + 578);
    const auto *sli0_580 = buffer.data(sli0 + 580);
    const auto *sli0_581 = buffer.data(sli0 + 581);
    const auto *sli0_583 = buffer.data(sli0 + 583);
    const auto *sli0_584 = buffer.data(sli0 + 584);
    const auto *sli0_585 = buffer.data(sli0 + 585);
    const auto *sli0_586 = buffer.data(sli0 + 586);
    const auto *sli0_587 = buffer.data(sli0 + 587);
    const auto *sli0_588 = buffer.data(sli0 + 588);
    const auto *sli0_591 = buffer.data(sli0 + 591);
    const auto *sli0_593 = buffer.data(sli0 + 593);
    const auto *sli0_594 = buffer.data(sli0 + 594);
    const auto *sli0_597 = buffer.data(sli0 + 597);
    const auto *sli0_598 = buffer.data(sli0 + 598);
    const auto *sli0_600 = buffer.data(sli0 + 600);
    const auto *sli0_602 = buffer.data(sli0 + 602);
    const auto *sli0_603 = buffer.data(sli0 + 603);
    const auto *sli0_605 = buffer.data(sli0 + 605);
    const auto *sli0_606 = buffer.data(sli0 + 606);
    const auto *sli0_608 = buffer.data(sli0 + 608);
    const auto *sli0_609 = buffer.data(sli0 + 609);
    const auto *sli0_611 = buffer.data(sli0 + 611);
    const auto *sli0_612 = buffer.data(sli0 + 612);
    const auto *sli0_613 = buffer.data(sli0 + 613);
    const auto *sli0_614 = buffer.data(sli0 + 614);
    const auto *sli0_615 = buffer.data(sli0 + 615);
    const auto *sli0_621 = buffer.data(sli0 + 621);

    const auto *sli1_553 = buffer.data(sli1 + 553);
    const auto *sli1_555 = buffer.data(sli1 + 555);
    const auto *sli1_556 = buffer.data(sli1 + 556);
    const auto *sli1_557 = buffer.data(sli1 + 557);
    const auto *sli1_558 = buffer.data(sli1 + 558);
    const auto *sli1_559 = buffer.data(sli1 + 559);
    const auto *sli1_560 = buffer.data(sli1 + 560);
    const auto *sli1_563 = buffer.data(sli1 + 563);
    const auto *sli1_565 = buffer.data(sli1 + 565);
    const auto *sli1_566 = buffer.data(sli1 + 566);
    const auto *sli1_569 = buffer.data(sli1 + 569);
    const auto *sli1_570 = buffer.data(sli1 + 570);
    const auto *sli1_572 = buffer.data(sli1 + 572);
    const auto *sli1_574 = buffer.data(sli1 + 574);
    const auto *sli1_575 = buffer.data(sli1 + 575);
    const auto *sli1_577 = buffer.data(sli1 + 577);
    const auto *sli1_578 = buffer.data(sli1 + 578);
    const auto *sli1_580 = buffer.data(sli1 + 580);
    const auto *sli1_581 = buffer.data(sli1 + 581);
    const auto *sli1_583 = buffer.data(sli1 + 583);
    const auto *sli1_584 = buffer.data(sli1 + 584);
    const auto *sli1_585 = buffer.data(sli1 + 585);
    const auto *sli1_586 = buffer.data(sli1 + 586);
    const auto *sli1_587 = buffer.data(sli1 + 587);
    const auto *sli1_588 = buffer.data(sli1 + 588);
    const auto *sli1_591 = buffer.data(sli1 + 591);
    const auto *sli1_593 = buffer.data(sli1 + 593);
    const auto *sli1_594 = buffer.data(sli1 + 594);
    const auto *sli1_597 = buffer.data(sli1 + 597);
    const auto *sli1_598 = buffer.data(sli1 + 598);
    const auto *sli1_600 = buffer.data(sli1 + 600);
    const auto *sli1_602 = buffer.data(sli1 + 602);
    const auto *sli1_603 = buffer.data(sli1 + 603);
    const auto *sli1_605 = buffer.data(sli1 + 605);
    const auto *sli1_606 = buffer.data(sli1 + 606);
    const auto *sli1_608 = buffer.data(sli1 + 608);
    const auto *sli1_609 = buffer.data(sli1 + 609);
    const auto *sli1_611 = buffer.data(sli1 + 611);
    const auto *sli1_612 = buffer.data(sli1 + 612);
    const auto *sli1_613 = buffer.data(sli1 + 613);
    const auto *sli1_614 = buffer.data(sli1 + 614);
    const auto *sli1_615 = buffer.data(sli1 + 615);
    const auto *sli1_621 = buffer.data(sli1 + 621);

    const auto *slk_712 = buffer.data(slk + 712);
    const auto *slk_714 = buffer.data(slk + 714);
    const auto *slk_715 = buffer.data(slk + 715);
    const auto *slk_716 = buffer.data(slk + 716);
    const auto *slk_717 = buffer.data(slk + 717);
    const auto *slk_718 = buffer.data(slk + 718);
    const auto *slk_719 = buffer.data(slk + 719);
    const auto *slk_720 = buffer.data(slk + 720);
    const auto *slk_722 = buffer.data(slk + 722);
    const auto *slk_723 = buffer.data(slk + 723);
    const auto *slk_725 = buffer.data(slk + 725);
    const auto *slk_726 = buffer.data(slk + 726);
    const auto *slk_729 = buffer.data(slk + 729);
    const auto *slk_730 = buffer.data(slk + 730);
    const auto *slk_732 = buffer.data(slk + 732);
    const auto *slk_734 = buffer.data(slk + 734);
    const auto *slk_735 = buffer.data(slk + 735);
    const auto *slk_737 = buffer.data(slk + 737);
    const auto *slk_738 = buffer.data(slk + 738);
    const auto *slk_740 = buffer.data(slk + 740);
    const auto *slk_741 = buffer.data(slk + 741);
    const auto *slk_743 = buffer.data(slk + 743);
    const auto *slk_744 = buffer.data(slk + 744);
    const auto *slk_745 = buffer.data(slk + 745);
    const auto *slk_747 = buffer.data(slk + 747);
    const auto *slk_748 = buffer.data(slk + 748);
    const auto *slk_749 = buffer.data(slk + 749);
    const auto *slk_750 = buffer.data(slk + 750);
    const auto *slk_751 = buffer.data(slk + 751);
    const auto *slk_752 = buffer.data(slk + 752);
    const auto *slk_753 = buffer.data(slk + 753);
    const auto *slk_754 = buffer.data(slk + 754);
    const auto *slk_755 = buffer.data(slk + 755);
    const auto *slk_756 = buffer.data(slk + 756);
    const auto *slk_758 = buffer.data(slk + 758);
    const auto *slk_759 = buffer.data(slk + 759);
    const auto *slk_761 = buffer.data(slk + 761);
    const auto *slk_762 = buffer.data(slk + 762);
    const auto *slk_765 = buffer.data(slk + 765);
    const auto *slk_766 = buffer.data(slk + 766);
    const auto *slk_768 = buffer.data(slk + 768);
    const auto *slk_770 = buffer.data(slk + 770);
    const auto *slk_771 = buffer.data(slk + 771);
    const auto *slk_773 = buffer.data(slk + 773);
    const auto *slk_774 = buffer.data(slk + 774);
    const auto *slk_776 = buffer.data(slk + 776);
    const auto *slk_777 = buffer.data(slk + 777);
    const auto *slk_779 = buffer.data(slk + 779);
    const auto *slk_780 = buffer.data(slk + 780);
    const auto *slk_781 = buffer.data(slk + 781);
    const auto *slk_783 = buffer.data(slk + 783);
    const auto *slk_784 = buffer.data(slk + 784);
    const auto *slk_785 = buffer.data(slk + 785);
    const auto *slk_786 = buffer.data(slk + 786);
    const auto *slk_787 = buffer.data(slk + 787);
    const auto *slk_788 = buffer.data(slk + 788);
    const auto *slk_789 = buffer.data(slk + 789);
    const auto *slk_790 = buffer.data(slk + 790);
    const auto *slk_791 = buffer.data(slk + 791);
    const auto *slk_792 = buffer.data(slk + 792);
    const auto *slk_794 = buffer.data(slk + 794);
    const auto *slk_795 = buffer.data(slk + 795);
    const auto *slk_797 = buffer.data(slk + 797);

#pragma omp simd aligned(t_889, t_890, t_891, t_892, pc_x, pc_y, pc_z, skk_496, skk_532, \
                         skk_718, skk_719, sli0_553, sli1_553, slk_712, slk_718, \
                         slk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_889[k] = f_17 * skk_718[k]
                   + f_3 * pc_x[k] * slk_718[k];

        t_890[k] = f_17 * skk_719[k]
                   + f_3 * pc_x[k] * slk_719[k];

        t_891[k] = f_15 * skk_532[k]
                   + f_1 * sli0_553[k]
                   - f_2 * sli1_553[k]
                   + f_3 * pc_y[k] * slk_712[k];

        t_892[k] = f_18 * skk_496[k]
                   + f_3 * pc_z[k] * slk_712[k];
    }

#pragma omp simd aligned(t_893, t_894, t_895, pc_y, skk_534, skk_535, skk_536, sli0_555, \
                         sli0_556, sli0_557, sli1_555, sli1_556, sli1_557, slk_714, slk_715, \
                         slk_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_893[k] = f_15 * skk_534[k]
                   + f_4 * sli0_555[k]
                   - f_5 * sli1_555[k]
                   + f_3 * pc_y[k] * slk_714[k];

        t_894[k] = f_15 * skk_535[k]
                   + f_6 * sli0_556[k]
                   - f_7 * sli1_556[k]
                   + f_3 * pc_y[k] * slk_715[k];

        t_895[k] = f_15 * skk_536[k]
                   + f_8 * sli0_557[k]
                   - f_9 * sli1_557[k]
                   + f_3 * pc_y[k] * slk_716[k];
    }

#pragma omp simd aligned(t_896, t_897, t_898, pc_y, skk_537, skk_538, skk_539, sli0_558, \
                         sli0_559, sli1_558, sli1_559, slk_717, slk_718, \
                         slk_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_896[k] = f_15 * skk_537[k]
                   + f_10 * sli0_558[k]
                   - f_11 * sli1_558[k]
                   + f_3 * pc_y[k] * slk_717[k];

        t_897[k] = f_15 * skk_538[k]
                   + f_12 * sli0_559[k]
                   - f_13 * sli1_559[k]
                   + f_3 * pc_y[k] * slk_718[k];

        t_898[k] = f_15 * skk_539[k]
                   + f_3 * pc_y[k] * slk_719[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pb_y, pc_x, pc_y, pc_z, skl0_674, \
                         skk_504, skk_720, skl1_674, sli0_560, sli1_560, \
                         slk_720 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pb_y[k] * skl0_674[k]
                   - f_14 * pc_y[k] * skl1_674[k];

        t_900[k] = f_17 * skk_720[k]
                   + f_1 * sli0_560[k]
                   - f_2 * sli1_560[k]
                   + f_3 * pc_x[k] * slk_720[k];

        t_901[k] = f_3 * pc_y[k] * slk_720[k];

        t_902[k] = f_19 * skk_504[k]
                   + f_3 * pc_z[k] * slk_720[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, pc_x, pc_y, skk_723, skk_725, sli0_563, \
                         sli0_565, sli1_563, sli1_565, slk_722, slk_723, \
                         slk_725 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = f_17 * skk_723[k]
                   + f_4 * sli0_563[k]
                   - f_5 * sli1_563[k]
                   + f_3 * pc_x[k] * slk_723[k];

        t_904[k] = f_3 * pc_y[k] * slk_722[k];

        t_905[k] = f_17 * skk_725[k]
                   + f_4 * sli0_565[k]
                   - f_5 * sli1_565[k]
                   + f_3 * pc_x[k] * slk_725[k];
    }

#pragma omp simd aligned(t_906, t_907, t_908, pc_x, pc_y, pc_z, skk_507, skk_726, sli0_566, \
                         sli1_566, slk_723, slk_725, slk_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_906[k] = f_17 * skk_726[k]
                   + f_6 * sli0_566[k]
                   - f_7 * sli1_566[k]
                   + f_3 * pc_x[k] * slk_726[k];

        t_907[k] = f_19 * skk_507[k]
                   + f_3 * pc_z[k] * slk_723[k];

        t_908[k] = f_3 * pc_y[k] * slk_725[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, pc_x, pc_z, skk_510, skk_729, skk_730, sli0_569, \
                         sli0_570, sli1_569, sli1_570, slk_726, slk_729, \
                         slk_730 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_17 * skk_729[k]
                   + f_6 * sli0_569[k]
                   - f_7 * sli1_569[k]
                   + f_3 * pc_x[k] * slk_729[k];

        t_910[k] = f_17 * skk_730[k]
                   + f_8 * sli0_570[k]
                   - f_9 * sli1_570[k]
                   + f_3 * pc_x[k] * slk_730[k];

        t_911[k] = f_19 * skk_510[k]
                   + f_3 * pc_z[k] * slk_726[k];
    }

#pragma omp simd aligned(t_912, t_913, t_914, pc_x, pc_y, skk_732, skk_734, sli0_572, \
                         sli0_574, sli1_572, sli1_574, slk_729, slk_732, \
                         slk_734 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_912[k] = f_17 * skk_732[k]
                   + f_8 * sli0_572[k]
                   - f_9 * sli1_572[k]
                   + f_3 * pc_x[k] * slk_732[k];

        t_913[k] = f_3 * pc_y[k] * slk_729[k];

        t_914[k] = f_17 * skk_734[k]
                   + f_8 * sli0_574[k]
                   - f_9 * sli1_574[k]
                   + f_3 * pc_x[k] * slk_734[k];
    }

#pragma omp simd aligned(t_915, t_916, t_917, pc_x, pc_z, skk_514, skk_735, skk_737, sli0_575, \
                         sli0_577, sli1_575, sli1_577, slk_730, slk_735, \
                         slk_737 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_915[k] = f_17 * skk_735[k]
                   + f_10 * sli0_575[k]
                   - f_11 * sli1_575[k]
                   + f_3 * pc_x[k] * slk_735[k];

        t_916[k] = f_19 * skk_514[k]
                   + f_3 * pc_z[k] * slk_730[k];

        t_917[k] = f_17 * skk_737[k]
                   + f_10 * sli0_577[k]
                   - f_11 * sli1_577[k]
                   + f_3 * pc_x[k] * slk_737[k];
    }

#pragma omp simd aligned(t_918, t_919, t_920, pc_x, pc_y, skk_738, skk_740, sli0_578, \
                         sli0_580, sli1_578, sli1_580, slk_734, slk_738, \
                         slk_740 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_918[k] = f_17 * skk_738[k]
                   + f_10 * sli0_578[k]
                   - f_11 * sli1_578[k]
                   + f_3 * pc_x[k] * slk_738[k];

        t_919[k] = f_3 * pc_y[k] * slk_734[k];

        t_920[k] = f_17 * skk_740[k]
                   + f_10 * sli0_580[k]
                   - f_11 * sli1_580[k]
                   + f_3 * pc_x[k] * slk_740[k];
    }

#pragma omp simd aligned(t_921, t_922, t_923, pc_x, pc_z, skk_519, skk_741, skk_743, sli0_581, \
                         sli0_583, sli1_581, sli1_583, slk_735, slk_741, \
                         slk_743 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_921[k] = f_17 * skk_741[k]
                   + f_12 * sli0_581[k]
                   - f_13 * sli1_581[k]
                   + f_3 * pc_x[k] * slk_741[k];

        t_922[k] = f_19 * skk_519[k]
                   + f_3 * pc_z[k] * slk_735[k];

        t_923[k] = f_17 * skk_743[k]
                   + f_12 * sli0_583[k]
                   - f_13 * sli1_583[k]
                   + f_3 * pc_x[k] * slk_743[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, pc_x, pc_y, skk_744, skk_745, sli0_584, \
                         sli0_585, sli1_584, sli1_585, slk_740, slk_744, \
                         slk_745 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_17 * skk_744[k]
                   + f_12 * sli0_584[k]
                   - f_13 * sli1_584[k]
                   + f_3 * pc_x[k] * slk_744[k];

        t_925[k] = f_17 * skk_745[k]
                   + f_12 * sli0_585[k]
                   - f_13 * sli1_585[k]
                   + f_3 * pc_x[k] * slk_745[k];

        t_926[k] = f_3 * pc_y[k] * slk_740[k];
    }

#pragma omp simd aligned(t_927, t_928, t_929, t_930, pc_x, skk_747, skk_748, skk_749, skk_750, \
                         sli0_587, sli1_587, slk_747, slk_748, slk_749, \
                         slk_750 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_927[k] = f_17 * skk_747[k]
                   + f_12 * sli0_587[k]
                   - f_13 * sli1_587[k]
                   + f_3 * pc_x[k] * slk_747[k];

        t_928[k] = f_17 * skk_748[k]
                   + f_3 * pc_x[k] * slk_748[k];

        t_929[k] = f_17 * skk_749[k]
                   + f_3 * pc_x[k] * slk_749[k];

        t_930[k] = f_17 * skk_750[k]
                   + f_3 * pc_x[k] * slk_750[k];
    }

#pragma omp simd aligned(t_931, t_932, t_933, t_934, t_935, pc_x, skk_751, skk_752, skk_753, \
                         skk_754, skk_755, slk_751, slk_752, slk_753, slk_754, \
                         slk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_931[k] = f_17 * skk_751[k]
                   + f_3 * pc_x[k] * slk_751[k];

        t_932[k] = f_17 * skk_752[k]
                   + f_3 * pc_x[k] * slk_752[k];

        t_933[k] = f_17 * skk_753[k]
                   + f_3 * pc_x[k] * slk_753[k];

        t_934[k] = f_17 * skk_754[k]
                   + f_3 * pc_x[k] * slk_754[k];

        t_935[k] = f_17 * skk_755[k]
                   + f_3 * pc_x[k] * slk_755[k];
    }

#pragma omp simd aligned(t_936, t_937, t_938, t_939, pc_y, pc_z, skk_532, sli0_581, sli0_583, \
                         sli0_584, sli1_581, sli1_583, sli1_584, slk_748, slk_750, \
                         slk_751 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_936[k] = f_1 * sli0_581[k]
                   - f_2 * sli1_581[k]
                   + f_3 * pc_y[k] * slk_748[k];

        t_937[k] = f_19 * skk_532[k]
                   + f_3 * pc_z[k] * slk_748[k];

        t_938[k] = f_4 * sli0_583[k]
                   - f_5 * sli1_583[k]
                   + f_3 * pc_y[k] * slk_750[k];

        t_939[k] = f_6 * sli0_584[k]
                   - f_7 * sli1_584[k]
                   + f_3 * pc_y[k] * slk_751[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pc_y, sli0_585, sli0_586, sli0_587, \
                         sli1_585, sli1_586, sli1_587, slk_752, slk_753, slk_754, \
                         slk_755 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_8 * sli0_585[k]
                   - f_9 * sli1_585[k]
                   + f_3 * pc_y[k] * slk_752[k];

        t_941[k] = f_10 * sli0_586[k]
                   - f_11 * sli1_586[k]
                   + f_3 * pc_y[k] * slk_753[k];

        t_942[k] = f_12 * sli0_587[k]
                   - f_13 * sli1_587[k]
                   + f_3 * pc_y[k] * slk_754[k];

        t_943[k] = f_3 * pc_y[k] * slk_755[k];
    }

#pragma omp simd aligned(t_944, t_945, t_946, t_947, pc_x, pc_y, pc_z, skk_539, skk_540, \
                         skk_756, sli0_587, sli0_588, sli1_587, sli1_588, slk_755, \
                         slk_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_19 * skk_539[k]
                   + f_1 * sli0_587[k]
                   - f_2 * sli1_587[k]
                   + f_3 * pc_z[k] * slk_755[k];

        t_945[k] = f_16 * skk_756[k]
                   + f_1 * sli0_588[k]
                   - f_2 * sli1_588[k]
                   + f_3 * pc_x[k] * slk_756[k];

        t_946[k] = f_20 * skk_540[k]
                   + f_3 * pc_y[k] * slk_756[k];

        t_947[k] = f_3 * pc_z[k] * slk_756[k];
    }

#pragma omp simd aligned(t_948, t_949, t_950, pc_x, pc_y, skk_542, skk_759, skk_761, sli0_591, \
                         sli0_593, sli1_591, sli1_593, slk_758, slk_759, \
                         slk_761 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_948[k] = f_16 * skk_759[k]
                   + f_4 * sli0_591[k]
                   - f_5 * sli1_591[k]
                   + f_3 * pc_x[k] * slk_759[k];

        t_949[k] = f_20 * skk_542[k]
                   + f_3 * pc_y[k] * slk_758[k];

        t_950[k] = f_16 * skk_761[k]
                   + f_4 * sli0_593[k]
                   - f_5 * sli1_593[k]
                   + f_3 * pc_x[k] * slk_761[k];
    }

#pragma omp simd aligned(t_951, t_952, t_953, pc_x, pc_y, pc_z, skk_545, skk_762, sli0_594, \
                         sli1_594, slk_759, slk_761, slk_762 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_951[k] = f_16 * skk_762[k]
                   + f_6 * sli0_594[k]
                   - f_7 * sli1_594[k]
                   + f_3 * pc_x[k] * slk_762[k];

        t_952[k] = f_3 * pc_z[k] * slk_759[k];

        t_953[k] = f_20 * skk_545[k]
                   + f_3 * pc_y[k] * slk_761[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, pc_x, pc_z, skk_765, skk_766, sli0_597, \
                         sli0_598, sli1_597, sli1_598, slk_762, slk_765, \
                         slk_766 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_16 * skk_765[k]
                   + f_6 * sli0_597[k]
                   - f_7 * sli1_597[k]
                   + f_3 * pc_x[k] * slk_765[k];

        t_955[k] = f_16 * skk_766[k]
                   + f_8 * sli0_598[k]
                   - f_9 * sli1_598[k]
                   + f_3 * pc_x[k] * slk_766[k];

        t_956[k] = f_3 * pc_z[k] * slk_762[k];
    }

#pragma omp simd aligned(t_957, t_958, t_959, pc_x, pc_y, skk_549, skk_768, skk_770, sli0_600, \
                         sli0_602, sli1_600, sli1_602, slk_765, slk_768, \
                         slk_770 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_957[k] = f_16 * skk_768[k]
                   + f_8 * sli0_600[k]
                   - f_9 * sli1_600[k]
                   + f_3 * pc_x[k] * slk_768[k];

        t_958[k] = f_20 * skk_549[k]
                   + f_3 * pc_y[k] * slk_765[k];

        t_959[k] = f_16 * skk_770[k]
                   + f_8 * sli0_602[k]
                   - f_9 * sli1_602[k]
                   + f_3 * pc_x[k] * slk_770[k];
    }

#pragma omp simd aligned(t_960, t_961, t_962, pc_x, pc_z, skk_771, skk_773, sli0_603, \
                         sli0_605, sli1_603, sli1_605, slk_766, slk_771, \
                         slk_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_960[k] = f_16 * skk_771[k]
                   + f_10 * sli0_603[k]
                   - f_11 * sli1_603[k]
                   + f_3 * pc_x[k] * slk_771[k];

        t_961[k] = f_3 * pc_z[k] * slk_766[k];

        t_962[k] = f_16 * skk_773[k]
                   + f_10 * sli0_605[k]
                   - f_11 * sli1_605[k]
                   + f_3 * pc_x[k] * slk_773[k];
    }

#pragma omp simd aligned(t_963, t_964, t_965, pc_x, pc_y, skk_554, skk_774, skk_776, sli0_606, \
                         sli0_608, sli1_606, sli1_608, slk_770, slk_774, \
                         slk_776 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_963[k] = f_16 * skk_774[k]
                   + f_10 * sli0_606[k]
                   - f_11 * sli1_606[k]
                   + f_3 * pc_x[k] * slk_774[k];

        t_964[k] = f_20 * skk_554[k]
                   + f_3 * pc_y[k] * slk_770[k];

        t_965[k] = f_16 * skk_776[k]
                   + f_10 * sli0_608[k]
                   - f_11 * sli1_608[k]
                   + f_3 * pc_x[k] * slk_776[k];
    }

#pragma omp simd aligned(t_966, t_967, t_968, pc_x, pc_z, skk_777, skk_779, sli0_609, \
                         sli0_611, sli1_609, sli1_611, slk_771, slk_777, \
                         slk_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_966[k] = f_16 * skk_777[k]
                   + f_12 * sli0_609[k]
                   - f_13 * sli1_609[k]
                   + f_3 * pc_x[k] * slk_777[k];

        t_967[k] = f_3 * pc_z[k] * slk_771[k];

        t_968[k] = f_16 * skk_779[k]
                   + f_12 * sli0_611[k]
                   - f_13 * sli1_611[k]
                   + f_3 * pc_x[k] * slk_779[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, pc_x, pc_y, skk_560, skk_780, skk_781, sli0_612, \
                         sli0_613, sli1_612, sli1_613, slk_776, slk_780, \
                         slk_781 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_16 * skk_780[k]
                   + f_12 * sli0_612[k]
                   - f_13 * sli1_612[k]
                   + f_3 * pc_x[k] * slk_780[k];

        t_970[k] = f_16 * skk_781[k]
                   + f_12 * sli0_613[k]
                   - f_13 * sli1_613[k]
                   + f_3 * pc_x[k] * slk_781[k];

        t_971[k] = f_20 * skk_560[k]
                   + f_3 * pc_y[k] * slk_776[k];
    }

#pragma omp simd aligned(t_972, t_973, t_974, t_975, pc_x, skk_783, skk_784, skk_785, skk_786, \
                         sli0_615, sli1_615, slk_783, slk_784, slk_785, \
                         slk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_972[k] = f_16 * skk_783[k]
                   + f_12 * sli0_615[k]
                   - f_13 * sli1_615[k]
                   + f_3 * pc_x[k] * slk_783[k];

        t_973[k] = f_16 * skk_784[k]
                   + f_3 * pc_x[k] * slk_784[k];

        t_974[k] = f_16 * skk_785[k]
                   + f_3 * pc_x[k] * slk_785[k];

        t_975[k] = f_16 * skk_786[k]
                   + f_3 * pc_x[k] * slk_786[k];
    }

#pragma omp simd aligned(t_976, t_977, t_978, t_979, t_980, pc_x, skk_787, skk_788, skk_789, \
                         skk_790, skk_791, slk_787, slk_788, slk_789, slk_790, \
                         slk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_976[k] = f_16 * skk_787[k]
                   + f_3 * pc_x[k] * slk_787[k];

        t_977[k] = f_16 * skk_788[k]
                   + f_3 * pc_x[k] * slk_788[k];

        t_978[k] = f_16 * skk_789[k]
                   + f_3 * pc_x[k] * slk_789[k];

        t_979[k] = f_16 * skk_790[k]
                   + f_3 * pc_x[k] * slk_790[k];

        t_980[k] = f_16 * skk_791[k]
                   + f_3 * pc_x[k] * slk_791[k];
    }

#pragma omp simd aligned(t_981, t_982, t_983, pc_y, pc_z, skk_568, skk_570, sli0_609, \
                         sli0_611, sli1_609, sli1_611, slk_784, \
                         slk_786 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_981[k] = f_20 * skk_568[k]
                   + f_1 * sli0_609[k]
                   - f_2 * sli1_609[k]
                   + f_3 * pc_y[k] * slk_784[k];

        t_982[k] = f_3 * pc_z[k] * slk_784[k];

        t_983[k] = f_20 * skk_570[k]
                   + f_4 * sli0_611[k]
                   - f_5 * sli1_611[k]
                   + f_3 * pc_y[k] * slk_786[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, pc_y, skk_571, skk_572, skk_573, sli0_612, \
                         sli0_613, sli0_614, sli1_612, sli1_613, sli1_614, slk_787, slk_788, \
                         slk_789 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_20 * skk_571[k]
                   + f_6 * sli0_612[k]
                   - f_7 * sli1_612[k]
                   + f_3 * pc_y[k] * slk_787[k];

        t_985[k] = f_20 * skk_572[k]
                   + f_8 * sli0_613[k]
                   - f_9 * sli1_613[k]
                   + f_3 * pc_y[k] * slk_788[k];

        t_986[k] = f_20 * skk_573[k]
                   + f_10 * sli0_614[k]
                   - f_11 * sli1_614[k]
                   + f_3 * pc_y[k] * slk_789[k];
    }

#pragma omp simd aligned(t_987, t_988, t_989, t_990, pb_z, pc_y, pc_z, skl0_675, skk_574, \
                         skk_575, skl1_675, sli0_615, sli1_615, slk_790, \
                         slk_791 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_987[k] = f_20 * skk_574[k]
                   + f_12 * sli0_615[k]
                   - f_13 * sli1_615[k]
                   + f_3 * pc_y[k] * slk_790[k];

        t_988[k] = f_20 * skk_575[k]
                   + f_3 * pc_y[k] * slk_791[k];

        t_989[k] = f_1 * sli0_615[k]
                   - f_2 * sli1_615[k]
                   + f_3 * pc_z[k] * slk_791[k];

        t_990[k] = pb_z[k] * skl0_675[k]
                   - f_14 * pc_z[k] * skl1_675[k];
    }

#pragma omp simd aligned(t_991, t_992, t_993, t_994, pb_z, pc_y, pc_z, skl0_678, skk_540, \
                         skk_576, skk_578, skl1_678, slk_792, slk_794 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_991[k] = f_19 * skk_576[k]
                   + f_3 * pc_y[k] * slk_792[k];

        t_992[k] = f_15 * skk_540[k]
                   + f_3 * pc_z[k] * slk_792[k];

        t_993[k] = pb_z[k] * skl0_678[k]
                   - f_14 * pc_z[k] * skl1_678[k];

        t_994[k] = f_19 * skk_578[k]
                   + f_3 * pc_y[k] * slk_794[k];
    }

#pragma omp simd aligned(t_995, t_996, t_997, pb_z, pc_x, pc_z, skl0_681, skk_543, skk_797, \
                         skl1_681, sli0_621, sli1_621, slk_795, \
                         slk_797 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_995[k] = f_16 * skk_797[k]
                   + f_4 * sli0_621[k]
                   - f_5 * sli1_621[k]
                   + f_3 * pc_x[k] * slk_797[k];

        t_996[k] = pb_z[k] * skl0_681[k]
                   - f_14 * pc_z[k] * skl1_681[k];

        t_997[k] = f_15 * skk_543[k]
                   + f_3 * pc_z[k] * slk_795[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t skl0,
                                                          const size_t skk, const size_t skl1,
                                                          const size_t sli0, const size_t sli1,
                                                          const size_t slk, const size_t ncols,
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

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_685 = buffer.data(skl0 + 685);
    const auto *skl0_687 = buffer.data(skl0 + 687);
    const auto *skl0_690 = buffer.data(skl0 + 690);
    const auto *skl0_692 = buffer.data(skl0 + 692);
    const auto *skl0_693 = buffer.data(skl0 + 693);
    const auto *skl0_696 = buffer.data(skl0 + 696);
    const auto *skl0_698 = buffer.data(skl0 + 698);
    const auto *skl0_699 = buffer.data(skl0 + 699);
    const auto *skl0_700 = buffer.data(skl0 + 700);
    const auto *skl0_711 = buffer.data(skl0 + 711);

    const auto *skk_546 = buffer.data(skk + 546);
    const auto *skk_547 = buffer.data(skk + 547);
    const auto *skk_550 = buffer.data(skk + 550);
    const auto *skk_551 = buffer.data(skk + 551);
    const auto *skk_552 = buffer.data(skk + 552);
    const auto *skk_555 = buffer.data(skk + 555);
    const auto *skk_556 = buffer.data(skk + 556);
    const auto *skk_557 = buffer.data(skk + 557);
    const auto *skk_558 = buffer.data(skk + 558);
    const auto *skk_568 = buffer.data(skk + 568);
    const auto *skk_575 = buffer.data(skk + 575);
    const auto *skk_576 = buffer.data(skk + 576);
    const auto *skk_579 = buffer.data(skk + 579);
    const auto *skk_581 = buffer.data(skk + 581);
    const auto *skk_582 = buffer.data(skk + 582);
    const auto *skk_585 = buffer.data(skk + 585);
    const auto *skk_586 = buffer.data(skk + 586);
    const auto *skk_590 = buffer.data(skk + 590);
    const auto *skk_591 = buffer.data(skk + 591);
    const auto *skk_596 = buffer.data(skk + 596);
    const auto *skk_604 = buffer.data(skk + 604);
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
    const auto *skk_626 = buffer.data(skk + 626);
    const auto *skk_627 = buffer.data(skk + 627);
    const auto *skk_632 = buffer.data(skk + 632);
    const auto *skk_640 = buffer.data(skk + 640);
    const auto *skk_642 = buffer.data(skk + 642);
    const auto *skk_643 = buffer.data(skk + 643);
    const auto *skk_644 = buffer.data(skk + 644);
    const auto *skk_645 = buffer.data(skk + 645);
    const auto *skk_646 = buffer.data(skk + 646);
    const auto *skk_647 = buffer.data(skk + 647);
    const auto *skk_648 = buffer.data(skk + 648);
    const auto *skk_650 = buffer.data(skk + 650);
    const auto *skk_653 = buffer.data(skk + 653);
    const auto *skk_657 = buffer.data(skk + 657);
    const auto *skk_662 = buffer.data(skk + 662);
    const auto *skk_801 = buffer.data(skk + 801);
    const auto *skk_806 = buffer.data(skk + 806);
    const auto *skk_812 = buffer.data(skk + 812);
    const auto *skk_819 = buffer.data(skk + 819);
    const auto *skk_820 = buffer.data(skk + 820);
    const auto *skk_821 = buffer.data(skk + 821);
    const auto *skk_822 = buffer.data(skk + 822);
    const auto *skk_823 = buffer.data(skk + 823);
    const auto *skk_824 = buffer.data(skk + 824);
    const auto *skk_825 = buffer.data(skk + 825);
    const auto *skk_826 = buffer.data(skk + 826);
    const auto *skk_827 = buffer.data(skk + 827);
    const auto *skk_828 = buffer.data(skk + 828);
    const auto *skk_831 = buffer.data(skk + 831);
    const auto *skk_833 = buffer.data(skk + 833);
    const auto *skk_834 = buffer.data(skk + 834);
    const auto *skk_837 = buffer.data(skk + 837);
    const auto *skk_838 = buffer.data(skk + 838);
    const auto *skk_840 = buffer.data(skk + 840);
    const auto *skk_842 = buffer.data(skk + 842);
    const auto *skk_843 = buffer.data(skk + 843);
    const auto *skk_845 = buffer.data(skk + 845);
    const auto *skk_846 = buffer.data(skk + 846);
    const auto *skk_848 = buffer.data(skk + 848);
    const auto *skk_849 = buffer.data(skk + 849);
    const auto *skk_851 = buffer.data(skk + 851);
    const auto *skk_852 = buffer.data(skk + 852);
    const auto *skk_853 = buffer.data(skk + 853);
    const auto *skk_855 = buffer.data(skk + 855);
    const auto *skk_856 = buffer.data(skk + 856);
    const auto *skk_857 = buffer.data(skk + 857);
    const auto *skk_858 = buffer.data(skk + 858);
    const auto *skk_859 = buffer.data(skk + 859);
    const auto *skk_860 = buffer.data(skk + 860);
    const auto *skk_861 = buffer.data(skk + 861);
    const auto *skk_862 = buffer.data(skk + 862);
    const auto *skk_863 = buffer.data(skk + 863);
    const auto *skk_864 = buffer.data(skk + 864);
    const auto *skk_867 = buffer.data(skk + 867);
    const auto *skk_869 = buffer.data(skk + 869);
    const auto *skk_870 = buffer.data(skk + 870);
    const auto *skk_873 = buffer.data(skk + 873);
    const auto *skk_874 = buffer.data(skk + 874);
    const auto *skk_876 = buffer.data(skk + 876);
    const auto *skk_878 = buffer.data(skk + 878);
    const auto *skk_879 = buffer.data(skk + 879);
    const auto *skk_881 = buffer.data(skk + 881);
    const auto *skk_882 = buffer.data(skk + 882);
    const auto *skk_884 = buffer.data(skk + 884);
    const auto *skk_885 = buffer.data(skk + 885);

    const auto *skl1_685 = buffer.data(skl1 + 685);
    const auto *skl1_687 = buffer.data(skl1 + 687);
    const auto *skl1_690 = buffer.data(skl1 + 690);
    const auto *skl1_692 = buffer.data(skl1 + 692);
    const auto *skl1_693 = buffer.data(skl1 + 693);
    const auto *skl1_696 = buffer.data(skl1 + 696);
    const auto *skl1_698 = buffer.data(skl1 + 698);
    const auto *skl1_699 = buffer.data(skl1 + 699);
    const auto *skl1_700 = buffer.data(skl1 + 700);
    const auto *skl1_711 = buffer.data(skl1 + 711);

    const auto *sli0_625 = buffer.data(sli0 + 625);
    const auto *sli0_630 = buffer.data(sli0 + 630);
    const auto *sli0_636 = buffer.data(sli0 + 636);
    const auto *sli0_639 = buffer.data(sli0 + 639);
    const auto *sli0_640 = buffer.data(sli0 + 640);
    const auto *sli0_641 = buffer.data(sli0 + 641);
    const auto *sli0_642 = buffer.data(sli0 + 642);
    const auto *sli0_643 = buffer.data(sli0 + 643);
    const auto *sli0_644 = buffer.data(sli0 + 644);
    const auto *sli0_647 = buffer.data(sli0 + 647);
    const auto *sli0_649 = buffer.data(sli0 + 649);
    const auto *sli0_650 = buffer.data(sli0 + 650);
    const auto *sli0_653 = buffer.data(sli0 + 653);
    const auto *sli0_654 = buffer.data(sli0 + 654);
    const auto *sli0_656 = buffer.data(sli0 + 656);
    const auto *sli0_658 = buffer.data(sli0 + 658);
    const auto *sli0_659 = buffer.data(sli0 + 659);
    const auto *sli0_661 = buffer.data(sli0 + 661);
    const auto *sli0_662 = buffer.data(sli0 + 662);
    const auto *sli0_664 = buffer.data(sli0 + 664);
    const auto *sli0_665 = buffer.data(sli0 + 665);
    const auto *sli0_667 = buffer.data(sli0 + 667);
    const auto *sli0_668 = buffer.data(sli0 + 668);
    const auto *sli0_669 = buffer.data(sli0 + 669);
    const auto *sli0_670 = buffer.data(sli0 + 670);
    const auto *sli0_671 = buffer.data(sli0 + 671);
    const auto *sli0_672 = buffer.data(sli0 + 672);
    const auto *sli0_675 = buffer.data(sli0 + 675);
    const auto *sli0_677 = buffer.data(sli0 + 677);
    const auto *sli0_678 = buffer.data(sli0 + 678);
    const auto *sli0_681 = buffer.data(sli0 + 681);
    const auto *sli0_682 = buffer.data(sli0 + 682);
    const auto *sli0_684 = buffer.data(sli0 + 684);
    const auto *sli0_686 = buffer.data(sli0 + 686);
    const auto *sli0_687 = buffer.data(sli0 + 687);
    const auto *sli0_689 = buffer.data(sli0 + 689);
    const auto *sli0_690 = buffer.data(sli0 + 690);
    const auto *sli0_692 = buffer.data(sli0 + 692);
    const auto *sli0_693 = buffer.data(sli0 + 693);

    const auto *sli1_625 = buffer.data(sli1 + 625);
    const auto *sli1_630 = buffer.data(sli1 + 630);
    const auto *sli1_636 = buffer.data(sli1 + 636);
    const auto *sli1_639 = buffer.data(sli1 + 639);
    const auto *sli1_640 = buffer.data(sli1 + 640);
    const auto *sli1_641 = buffer.data(sli1 + 641);
    const auto *sli1_642 = buffer.data(sli1 + 642);
    const auto *sli1_643 = buffer.data(sli1 + 643);
    const auto *sli1_644 = buffer.data(sli1 + 644);
    const auto *sli1_647 = buffer.data(sli1 + 647);
    const auto *sli1_649 = buffer.data(sli1 + 649);
    const auto *sli1_650 = buffer.data(sli1 + 650);
    const auto *sli1_653 = buffer.data(sli1 + 653);
    const auto *sli1_654 = buffer.data(sli1 + 654);
    const auto *sli1_656 = buffer.data(sli1 + 656);
    const auto *sli1_658 = buffer.data(sli1 + 658);
    const auto *sli1_659 = buffer.data(sli1 + 659);
    const auto *sli1_661 = buffer.data(sli1 + 661);
    const auto *sli1_662 = buffer.data(sli1 + 662);
    const auto *sli1_664 = buffer.data(sli1 + 664);
    const auto *sli1_665 = buffer.data(sli1 + 665);
    const auto *sli1_667 = buffer.data(sli1 + 667);
    const auto *sli1_668 = buffer.data(sli1 + 668);
    const auto *sli1_669 = buffer.data(sli1 + 669);
    const auto *sli1_670 = buffer.data(sli1 + 670);
    const auto *sli1_671 = buffer.data(sli1 + 671);
    const auto *sli1_672 = buffer.data(sli1 + 672);
    const auto *sli1_675 = buffer.data(sli1 + 675);
    const auto *sli1_677 = buffer.data(sli1 + 677);
    const auto *sli1_678 = buffer.data(sli1 + 678);
    const auto *sli1_681 = buffer.data(sli1 + 681);
    const auto *sli1_682 = buffer.data(sli1 + 682);
    const auto *sli1_684 = buffer.data(sli1 + 684);
    const auto *sli1_686 = buffer.data(sli1 + 686);
    const auto *sli1_687 = buffer.data(sli1 + 687);
    const auto *sli1_689 = buffer.data(sli1 + 689);
    const auto *sli1_690 = buffer.data(sli1 + 690);
    const auto *sli1_692 = buffer.data(sli1 + 692);
    const auto *sli1_693 = buffer.data(sli1 + 693);

    const auto *slk_797 = buffer.data(slk + 797);
    const auto *slk_798 = buffer.data(slk + 798);
    const auto *slk_801 = buffer.data(slk + 801);
    const auto *slk_802 = buffer.data(slk + 802);
    const auto *slk_806 = buffer.data(slk + 806);
    const auto *slk_807 = buffer.data(slk + 807);
    const auto *slk_812 = buffer.data(slk + 812);
    const auto *slk_819 = buffer.data(slk + 819);
    const auto *slk_820 = buffer.data(slk + 820);
    const auto *slk_821 = buffer.data(slk + 821);
    const auto *slk_822 = buffer.data(slk + 822);
    const auto *slk_823 = buffer.data(slk + 823);
    const auto *slk_824 = buffer.data(slk + 824);
    const auto *slk_825 = buffer.data(slk + 825);
    const auto *slk_826 = buffer.data(slk + 826);
    const auto *slk_827 = buffer.data(slk + 827);
    const auto *slk_828 = buffer.data(slk + 828);
    const auto *slk_830 = buffer.data(slk + 830);
    const auto *slk_831 = buffer.data(slk + 831);
    const auto *slk_833 = buffer.data(slk + 833);
    const auto *slk_834 = buffer.data(slk + 834);
    const auto *slk_837 = buffer.data(slk + 837);
    const auto *slk_838 = buffer.data(slk + 838);
    const auto *slk_840 = buffer.data(slk + 840);
    const auto *slk_842 = buffer.data(slk + 842);
    const auto *slk_843 = buffer.data(slk + 843);
    const auto *slk_845 = buffer.data(slk + 845);
    const auto *slk_846 = buffer.data(slk + 846);
    const auto *slk_848 = buffer.data(slk + 848);
    const auto *slk_849 = buffer.data(slk + 849);
    const auto *slk_851 = buffer.data(slk + 851);
    const auto *slk_852 = buffer.data(slk + 852);
    const auto *slk_853 = buffer.data(slk + 853);
    const auto *slk_855 = buffer.data(slk + 855);
    const auto *slk_856 = buffer.data(slk + 856);
    const auto *slk_857 = buffer.data(slk + 857);
    const auto *slk_858 = buffer.data(slk + 858);
    const auto *slk_859 = buffer.data(slk + 859);
    const auto *slk_860 = buffer.data(slk + 860);
    const auto *slk_861 = buffer.data(slk + 861);
    const auto *slk_862 = buffer.data(slk + 862);
    const auto *slk_863 = buffer.data(slk + 863);
    const auto *slk_864 = buffer.data(slk + 864);
    const auto *slk_866 = buffer.data(slk + 866);
    const auto *slk_867 = buffer.data(slk + 867);
    const auto *slk_869 = buffer.data(slk + 869);
    const auto *slk_870 = buffer.data(slk + 870);
    const auto *slk_873 = buffer.data(slk + 873);
    const auto *slk_874 = buffer.data(slk + 874);
    const auto *slk_876 = buffer.data(slk + 876);
    const auto *slk_878 = buffer.data(slk + 878);
    const auto *slk_879 = buffer.data(slk + 879);
    const auto *slk_881 = buffer.data(slk + 881);
    const auto *slk_882 = buffer.data(slk + 882);
    const auto *slk_884 = buffer.data(slk + 884);
    const auto *slk_885 = buffer.data(slk + 885);

#pragma omp simd aligned(t_998, t_999, t_1000, pb_z, pc_x, pc_y, pc_z, skl0_685, skk_581, \
                         skk_801, skl1_685, sli0_625, sli1_625, slk_797, \
                         slk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_19 * skk_581[k]
                   + f_3 * pc_y[k] * slk_797[k];

        t_999[k] = f_16 * skk_801[k]
                   + f_6 * sli0_625[k]
                   - f_7 * sli1_625[k]
                   + f_3 * pc_x[k] * slk_801[k];

        t_1000[k] = pb_z[k] * skl0_685[k]
                    - f_14 * pc_z[k] * skl1_685[k];
    }

#pragma omp simd aligned(t_1001, t_1002, t_1003, pb_z, pc_y, pc_z, skl0_687, skk_546, skk_547, \
                         skk_585, skl1_687, slk_798, slk_801 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1001[k] = f_15 * skk_546[k]
                    + f_3 * pc_z[k] * slk_798[k];

        t_1002[k] = pb_z[k] * skl0_687[k]
                    + f_16 * skk_547[k]
                    - f_14 * pc_z[k] * skl1_687[k];

        t_1003[k] = f_19 * skk_585[k]
                    + f_3 * pc_y[k] * slk_801[k];
    }

#pragma omp simd aligned(t_1004, t_1005, t_1006, pb_z, pc_x, pc_z, skl0_690, skk_550, skk_806, \
                         skl1_690, sli0_630, sli1_630, slk_802, \
                         slk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1004[k] = f_16 * skk_806[k]
                    + f_8 * sli0_630[k]
                    - f_9 * sli1_630[k]
                    + f_3 * pc_x[k] * slk_806[k];

        t_1005[k] = pb_z[k] * skl0_690[k]
                    - f_14 * pc_z[k] * skl1_690[k];

        t_1006[k] = f_15 * skk_550[k]
                    + f_3 * pc_z[k] * slk_802[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pb_z, pc_y, pc_z, skl0_692, skl0_693, \
                         skk_551, skk_552, skk_590, skl1_692, skl1_693, \
                         slk_806 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = pb_z[k] * skl0_692[k]
                    + f_16 * skk_551[k]
                    - f_14 * pc_z[k] * skl1_692[k];

        t_1008[k] = pb_z[k] * skl0_693[k]
                    + f_17 * skk_552[k]
                    - f_14 * pc_z[k] * skl1_693[k];

        t_1009[k] = f_19 * skk_590[k]
                    + f_3 * pc_y[k] * slk_806[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, pb_z, pc_x, pc_z, skl0_696, skk_555, skk_812, \
                         skl1_696, sli0_636, sli1_636, slk_807, \
                         slk_812 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_16 * skk_812[k]
                    + f_10 * sli0_636[k]
                    - f_11 * sli1_636[k]
                    + f_3 * pc_x[k] * slk_812[k];

        t_1011[k] = pb_z[k] * skl0_696[k]
                    - f_14 * pc_z[k] * skl1_696[k];

        t_1012[k] = f_15 * skk_555[k]
                    + f_3 * pc_z[k] * slk_807[k];
    }

#pragma omp simd aligned(t_1013, t_1014, t_1015, pb_z, pc_z, skl0_698, skl0_699, skl0_700, \
                         skk_556, skk_557, skk_558, skl1_698, skl1_699, \
                         skl1_700 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1013[k] = pb_z[k] * skl0_698[k]
                    + f_16 * skk_556[k]
                    - f_14 * pc_z[k] * skl1_698[k];

        t_1014[k] = pb_z[k] * skl0_699[k]
                    + f_17 * skk_557[k]
                    - f_14 * pc_z[k] * skl1_699[k];

        t_1015[k] = pb_z[k] * skl0_700[k]
                    + f_18 * skk_558[k]
                    - f_14 * pc_z[k] * skl1_700[k];
    }

#pragma omp simd aligned(t_1016, t_1017, t_1018, t_1019, pc_x, pc_y, skk_596, skk_819, \
                         skk_820, skk_821, sli0_643, sli1_643, slk_812, slk_819, slk_820, \
                         slk_821 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1016[k] = f_19 * skk_596[k]
                    + f_3 * pc_y[k] * slk_812[k];

        t_1017[k] = f_16 * skk_819[k]
                    + f_12 * sli0_643[k]
                    - f_13 * sli1_643[k]
                    + f_3 * pc_x[k] * slk_819[k];

        t_1018[k] = f_16 * skk_820[k]
                    + f_3 * pc_x[k] * slk_820[k];

        t_1019[k] = f_16 * skk_821[k]
                    + f_3 * pc_x[k] * slk_821[k];
    }

#pragma omp simd aligned(t_1020, t_1021, t_1022, t_1023, t_1024, pc_x, skk_822, skk_823, \
                         skk_824, skk_825, skk_826, slk_822, slk_823, slk_824, slk_825, \
                         slk_826 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1020[k] = f_16 * skk_822[k]
                    + f_3 * pc_x[k] * slk_822[k];

        t_1021[k] = f_16 * skk_823[k]
                    + f_3 * pc_x[k] * slk_823[k];

        t_1022[k] = f_16 * skk_824[k]
                    + f_3 * pc_x[k] * slk_824[k];

        t_1023[k] = f_16 * skk_825[k]
                    + f_3 * pc_x[k] * slk_825[k];

        t_1024[k] = f_16 * skk_826[k]
                    + f_3 * pc_x[k] * slk_826[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, pb_z, pc_x, pc_z, skl0_711, skk_568, skk_827, \
                         skl1_711, slk_820, slk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_16 * skk_827[k]
                    + f_3 * pc_x[k] * slk_827[k];

        t_1026[k] = pb_z[k] * skl0_711[k]
                    - f_14 * pc_z[k] * skl1_711[k];

        t_1027[k] = f_15 * skk_568[k]
                    + f_3 * pc_z[k] * slk_820[k];
    }

#pragma omp simd aligned(t_1028, t_1029, t_1030, pc_y, skk_606, skk_607, skk_608, sli0_639, \
                         sli0_640, sli0_641, sli1_639, sli1_640, sli1_641, slk_822, slk_823, \
                         slk_824 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1028[k] = f_19 * skk_606[k]
                    + f_4 * sli0_639[k]
                    - f_5 * sli1_639[k]
                    + f_3 * pc_y[k] * slk_822[k];

        t_1029[k] = f_19 * skk_607[k]
                    + f_6 * sli0_640[k]
                    - f_7 * sli1_640[k]
                    + f_3 * pc_y[k] * slk_823[k];

        t_1030[k] = f_19 * skk_608[k]
                    + f_8 * sli0_641[k]
                    - f_9 * sli1_641[k]
                    + f_3 * pc_y[k] * slk_824[k];
    }

#pragma omp simd aligned(t_1031, t_1032, t_1033, pc_y, skk_609, skk_610, skk_611, sli0_642, \
                         sli0_643, sli1_642, sli1_643, slk_825, slk_826, \
                         slk_827 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1031[k] = f_19 * skk_609[k]
                    + f_10 * sli0_642[k]
                    - f_11 * sli1_642[k]
                    + f_3 * pc_y[k] * slk_825[k];

        t_1032[k] = f_19 * skk_610[k]
                    + f_12 * sli0_643[k]
                    - f_13 * sli1_643[k]
                    + f_3 * pc_y[k] * slk_826[k];

        t_1033[k] = f_19 * skk_611[k]
                    + f_3 * pc_y[k] * slk_827[k];
    }

#pragma omp simd aligned(t_1034, t_1035, t_1036, pc_x, pc_y, pc_z, skk_575, skk_612, skk_828, \
                         sli0_643, sli0_644, sli1_643, sli1_644, slk_827, \
                         slk_828 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1034[k] = f_15 * skk_575[k]
                    + f_1 * sli0_643[k]
                    - f_2 * sli1_643[k]
                    + f_3 * pc_z[k] * slk_827[k];

        t_1035[k] = f_16 * skk_828[k]
                    + f_1 * sli0_644[k]
                    - f_2 * sli1_644[k]
                    + f_3 * pc_x[k] * slk_828[k];

        t_1036[k] = f_18 * skk_612[k]
                    + f_3 * pc_y[k] * slk_828[k];
    }

#pragma omp simd aligned(t_1037, t_1038, t_1039, pc_x, pc_y, pc_z, skk_576, skk_614, skk_831, \
                         sli0_647, sli1_647, slk_828, slk_830, \
                         slk_831 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1037[k] = f_16 * skk_576[k]
                    + f_3 * pc_z[k] * slk_828[k];

        t_1038[k] = f_16 * skk_831[k]
                    + f_4 * sli0_647[k]
                    - f_5 * sli1_647[k]
                    + f_3 * pc_x[k] * slk_831[k];

        t_1039[k] = f_18 * skk_614[k]
                    + f_3 * pc_y[k] * slk_830[k];
    }

#pragma omp simd aligned(t_1040, t_1041, t_1042, pc_x, pc_z, skk_579, skk_833, skk_834, \
                         sli0_649, sli0_650, sli1_649, sli1_650, slk_831, slk_833, \
                         slk_834 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1040[k] = f_16 * skk_833[k]
                    + f_4 * sli0_649[k]
                    - f_5 * sli1_649[k]
                    + f_3 * pc_x[k] * slk_833[k];

        t_1041[k] = f_16 * skk_834[k]
                    + f_6 * sli0_650[k]
                    - f_7 * sli1_650[k]
                    + f_3 * pc_x[k] * slk_834[k];

        t_1042[k] = f_16 * skk_579[k]
                    + f_3 * pc_z[k] * slk_831[k];
    }

#pragma omp simd aligned(t_1043, t_1044, t_1045, pc_x, pc_y, skk_617, skk_837, skk_838, \
                         sli0_653, sli0_654, sli1_653, sli1_654, slk_833, slk_837, \
                         slk_838 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1043[k] = f_18 * skk_617[k]
                    + f_3 * pc_y[k] * slk_833[k];

        t_1044[k] = f_16 * skk_837[k]
                    + f_6 * sli0_653[k]
                    - f_7 * sli1_653[k]
                    + f_3 * pc_x[k] * slk_837[k];

        t_1045[k] = f_16 * skk_838[k]
                    + f_8 * sli0_654[k]
                    - f_9 * sli1_654[k]
                    + f_3 * pc_x[k] * slk_838[k];
    }

#pragma omp simd aligned(t_1046, t_1047, t_1048, pc_x, pc_y, pc_z, skk_582, skk_621, skk_840, \
                         sli0_656, sli1_656, slk_834, slk_837, \
                         slk_840 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1046[k] = f_16 * skk_582[k]
                    + f_3 * pc_z[k] * slk_834[k];

        t_1047[k] = f_16 * skk_840[k]
                    + f_8 * sli0_656[k]
                    - f_9 * sli1_656[k]
                    + f_3 * pc_x[k] * slk_840[k];

        t_1048[k] = f_18 * skk_621[k]
                    + f_3 * pc_y[k] * slk_837[k];
    }

#pragma omp simd aligned(t_1049, t_1050, t_1051, pc_x, pc_z, skk_586, skk_842, skk_843, \
                         sli0_658, sli0_659, sli1_658, sli1_659, slk_838, slk_842, \
                         slk_843 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1049[k] = f_16 * skk_842[k]
                    + f_8 * sli0_658[k]
                    - f_9 * sli1_658[k]
                    + f_3 * pc_x[k] * slk_842[k];

        t_1050[k] = f_16 * skk_843[k]
                    + f_10 * sli0_659[k]
                    - f_11 * sli1_659[k]
                    + f_3 * pc_x[k] * slk_843[k];

        t_1051[k] = f_16 * skk_586[k]
                    + f_3 * pc_z[k] * slk_838[k];
    }

#pragma omp simd aligned(t_1052, t_1053, t_1054, pc_x, pc_y, skk_626, skk_845, skk_846, \
                         sli0_661, sli0_662, sli1_661, sli1_662, slk_842, slk_845, \
                         slk_846 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1052[k] = f_16 * skk_845[k]
                    + f_10 * sli0_661[k]
                    - f_11 * sli1_661[k]
                    + f_3 * pc_x[k] * slk_845[k];

        t_1053[k] = f_16 * skk_846[k]
                    + f_10 * sli0_662[k]
                    - f_11 * sli1_662[k]
                    + f_3 * pc_x[k] * slk_846[k];

        t_1054[k] = f_18 * skk_626[k]
                    + f_3 * pc_y[k] * slk_842[k];
    }

#pragma omp simd aligned(t_1055, t_1056, t_1057, pc_x, pc_z, skk_591, skk_848, skk_849, \
                         sli0_664, sli0_665, sli1_664, sli1_665, slk_843, slk_848, \
                         slk_849 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1055[k] = f_16 * skk_848[k]
                    + f_10 * sli0_664[k]
                    - f_11 * sli1_664[k]
                    + f_3 * pc_x[k] * slk_848[k];

        t_1056[k] = f_16 * skk_849[k]
                    + f_12 * sli0_665[k]
                    - f_13 * sli1_665[k]
                    + f_3 * pc_x[k] * slk_849[k];

        t_1057[k] = f_16 * skk_591[k]
                    + f_3 * pc_z[k] * slk_843[k];
    }

#pragma omp simd aligned(t_1058, t_1059, t_1060, pc_x, skk_851, skk_852, skk_853, sli0_667, \
                         sli0_668, sli0_669, sli1_667, sli1_668, sli1_669, slk_851, slk_852, \
                         slk_853 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1058[k] = f_16 * skk_851[k]
                    + f_12 * sli0_667[k]
                    - f_13 * sli1_667[k]
                    + f_3 * pc_x[k] * slk_851[k];

        t_1059[k] = f_16 * skk_852[k]
                    + f_12 * sli0_668[k]
                    - f_13 * sli1_668[k]
                    + f_3 * pc_x[k] * slk_852[k];

        t_1060[k] = f_16 * skk_853[k]
                    + f_12 * sli0_669[k]
                    - f_13 * sli1_669[k]
                    + f_3 * pc_x[k] * slk_853[k];
    }

#pragma omp simd aligned(t_1061, t_1062, t_1063, t_1064, pc_x, pc_y, skk_632, skk_855, \
                         skk_856, skk_857, sli0_671, sli1_671, slk_848, slk_855, slk_856, \
                         slk_857 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1061[k] = f_18 * skk_632[k]
                    + f_3 * pc_y[k] * slk_848[k];

        t_1062[k] = f_16 * skk_855[k]
                    + f_12 * sli0_671[k]
                    - f_13 * sli1_671[k]
                    + f_3 * pc_x[k] * slk_855[k];

        t_1063[k] = f_16 * skk_856[k]
                    + f_3 * pc_x[k] * slk_856[k];

        t_1064[k] = f_16 * skk_857[k]
                    + f_3 * pc_x[k] * slk_857[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, t_1068, t_1069, pc_x, skk_858, skk_859, \
                         skk_860, skk_861, skk_862, slk_858, slk_859, slk_860, slk_861, \
                         slk_862 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_16 * skk_858[k]
                    + f_3 * pc_x[k] * slk_858[k];

        t_1066[k] = f_16 * skk_859[k]
                    + f_3 * pc_x[k] * slk_859[k];

        t_1067[k] = f_16 * skk_860[k]
                    + f_3 * pc_x[k] * slk_860[k];

        t_1068[k] = f_16 * skk_861[k]
                    + f_3 * pc_x[k] * slk_861[k];

        t_1069[k] = f_16 * skk_862[k]
                    + f_3 * pc_x[k] * slk_862[k];
    }

#pragma omp simd aligned(t_1070, t_1071, t_1072, pc_x, pc_y, pc_z, skk_604, skk_640, skk_863, \
                         sli0_665, sli1_665, slk_856, slk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1070[k] = f_16 * skk_863[k]
                    + f_3 * pc_x[k] * slk_863[k];

        t_1071[k] = f_18 * skk_640[k]
                    + f_1 * sli0_665[k]
                    - f_2 * sli1_665[k]
                    + f_3 * pc_y[k] * slk_856[k];

        t_1072[k] = f_16 * skk_604[k]
                    + f_3 * pc_z[k] * slk_856[k];
    }

#pragma omp simd aligned(t_1073, t_1074, t_1075, pc_y, skk_642, skk_643, skk_644, sli0_667, \
                         sli0_668, sli0_669, sli1_667, sli1_668, sli1_669, slk_858, slk_859, \
                         slk_860 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1073[k] = f_18 * skk_642[k]
                    + f_4 * sli0_667[k]
                    - f_5 * sli1_667[k]
                    + f_3 * pc_y[k] * slk_858[k];

        t_1074[k] = f_18 * skk_643[k]
                    + f_6 * sli0_668[k]
                    - f_7 * sli1_668[k]
                    + f_3 * pc_y[k] * slk_859[k];

        t_1075[k] = f_18 * skk_644[k]
                    + f_8 * sli0_669[k]
                    - f_9 * sli1_669[k]
                    + f_3 * pc_y[k] * slk_860[k];
    }

#pragma omp simd aligned(t_1076, t_1077, t_1078, pc_y, skk_645, skk_646, skk_647, sli0_670, \
                         sli0_671, sli1_670, sli1_671, slk_861, slk_862, \
                         slk_863 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1076[k] = f_18 * skk_645[k]
                    + f_10 * sli0_670[k]
                    - f_11 * sli1_670[k]
                    + f_3 * pc_y[k] * slk_861[k];

        t_1077[k] = f_18 * skk_646[k]
                    + f_12 * sli0_671[k]
                    - f_13 * sli1_671[k]
                    + f_3 * pc_y[k] * slk_862[k];

        t_1078[k] = f_18 * skk_647[k]
                    + f_3 * pc_y[k] * slk_863[k];
    }

#pragma omp simd aligned(t_1079, t_1080, t_1081, pc_x, pc_y, pc_z, skk_611, skk_648, skk_864, \
                         sli0_671, sli0_672, sli1_671, sli1_672, slk_863, \
                         slk_864 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1079[k] = f_16 * skk_611[k]
                    + f_1 * sli0_671[k]
                    - f_2 * sli1_671[k]
                    + f_3 * pc_z[k] * slk_863[k];

        t_1080[k] = f_16 * skk_864[k]
                    + f_1 * sli0_672[k]
                    - f_2 * sli1_672[k]
                    + f_3 * pc_x[k] * slk_864[k];

        t_1081[k] = f_17 * skk_648[k]
                    + f_3 * pc_y[k] * slk_864[k];
    }

#pragma omp simd aligned(t_1082, t_1083, t_1084, pc_x, pc_y, pc_z, skk_612, skk_650, skk_867, \
                         sli0_675, sli1_675, slk_864, slk_866, \
                         slk_867 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1082[k] = f_17 * skk_612[k]
                    + f_3 * pc_z[k] * slk_864[k];

        t_1083[k] = f_16 * skk_867[k]
                    + f_4 * sli0_675[k]
                    - f_5 * sli1_675[k]
                    + f_3 * pc_x[k] * slk_867[k];

        t_1084[k] = f_17 * skk_650[k]
                    + f_3 * pc_y[k] * slk_866[k];
    }

#pragma omp simd aligned(t_1085, t_1086, t_1087, pc_x, pc_z, skk_615, skk_869, skk_870, \
                         sli0_677, sli0_678, sli1_677, sli1_678, slk_867, slk_869, \
                         slk_870 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1085[k] = f_16 * skk_869[k]
                    + f_4 * sli0_677[k]
                    - f_5 * sli1_677[k]
                    + f_3 * pc_x[k] * slk_869[k];

        t_1086[k] = f_16 * skk_870[k]
                    + f_6 * sli0_678[k]
                    - f_7 * sli1_678[k]
                    + f_3 * pc_x[k] * slk_870[k];

        t_1087[k] = f_17 * skk_615[k]
                    + f_3 * pc_z[k] * slk_867[k];
    }

#pragma omp simd aligned(t_1088, t_1089, t_1090, pc_x, pc_y, skk_653, skk_873, skk_874, \
                         sli0_681, sli0_682, sli1_681, sli1_682, slk_869, slk_873, \
                         slk_874 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1088[k] = f_17 * skk_653[k]
                    + f_3 * pc_y[k] * slk_869[k];

        t_1089[k] = f_16 * skk_873[k]
                    + f_6 * sli0_681[k]
                    - f_7 * sli1_681[k]
                    + f_3 * pc_x[k] * slk_873[k];

        t_1090[k] = f_16 * skk_874[k]
                    + f_8 * sli0_682[k]
                    - f_9 * sli1_682[k]
                    + f_3 * pc_x[k] * slk_874[k];
    }

#pragma omp simd aligned(t_1091, t_1092, t_1093, pc_x, pc_y, pc_z, skk_618, skk_657, skk_876, \
                         sli0_684, sli1_684, slk_870, slk_873, \
                         slk_876 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1091[k] = f_17 * skk_618[k]
                    + f_3 * pc_z[k] * slk_870[k];

        t_1092[k] = f_16 * skk_876[k]
                    + f_8 * sli0_684[k]
                    - f_9 * sli1_684[k]
                    + f_3 * pc_x[k] * slk_876[k];

        t_1093[k] = f_17 * skk_657[k]
                    + f_3 * pc_y[k] * slk_873[k];
    }

#pragma omp simd aligned(t_1094, t_1095, t_1096, pc_x, pc_z, skk_622, skk_878, skk_879, \
                         sli0_686, sli0_687, sli1_686, sli1_687, slk_874, slk_878, \
                         slk_879 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1094[k] = f_16 * skk_878[k]
                    + f_8 * sli0_686[k]
                    - f_9 * sli1_686[k]
                    + f_3 * pc_x[k] * slk_878[k];

        t_1095[k] = f_16 * skk_879[k]
                    + f_10 * sli0_687[k]
                    - f_11 * sli1_687[k]
                    + f_3 * pc_x[k] * slk_879[k];

        t_1096[k] = f_17 * skk_622[k]
                    + f_3 * pc_z[k] * slk_874[k];
    }

#pragma omp simd aligned(t_1097, t_1098, t_1099, pc_x, pc_y, skk_662, skk_881, skk_882, \
                         sli0_689, sli0_690, sli1_689, sli1_690, slk_878, slk_881, \
                         slk_882 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1097[k] = f_16 * skk_881[k]
                    + f_10 * sli0_689[k]
                    - f_11 * sli1_689[k]
                    + f_3 * pc_x[k] * slk_881[k];

        t_1098[k] = f_16 * skk_882[k]
                    + f_10 * sli0_690[k]
                    - f_11 * sli1_690[k]
                    + f_3 * pc_x[k] * slk_882[k];

        t_1099[k] = f_17 * skk_662[k]
                    + f_3 * pc_y[k] * slk_878[k];
    }

#pragma omp simd aligned(t_1100, t_1101, t_1102, pc_x, pc_z, skk_627, skk_884, skk_885, \
                         sli0_692, sli0_693, sli1_692, sli1_693, slk_879, slk_884, \
                         slk_885 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1100[k] = f_16 * skk_884[k]
                    + f_10 * sli0_692[k]
                    - f_11 * sli1_692[k]
                    + f_3 * pc_x[k] * slk_884[k];

        t_1101[k] = f_16 * skk_885[k]
                    + f_12 * sli0_693[k]
                    - f_13 * sli1_693[k]
                    + f_3 * pc_x[k] * slk_885[k];

        t_1102[k] = f_17 * skk_627[k]
                    + f_3 * pc_z[k] * slk_879[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece10(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t sli0, const size_t sli1,
                                                           const size_t slk, const size_t ncols,
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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_900 = buffer.data(skl0 + 900);
    const auto *skl0_903 = buffer.data(skl0 + 903);
    const auto *skl0_905 = buffer.data(skl0 + 905);
    const auto *skl0_906 = buffer.data(skl0 + 906);
    const auto *skl0_909 = buffer.data(skl0 + 909);
    const auto *skl0_910 = buffer.data(skl0 + 910);
    const auto *skl0_912 = buffer.data(skl0 + 912);
    const auto *skl0_914 = buffer.data(skl0 + 914);
    const auto *skl0_915 = buffer.data(skl0 + 915);
    const auto *skl0_917 = buffer.data(skl0 + 917);
    const auto *skl0_918 = buffer.data(skl0 + 918);
    const auto *skl0_920 = buffer.data(skl0 + 920);
    const auto *skl0_921 = buffer.data(skl0 + 921);
    const auto *skl0_923 = buffer.data(skl0 + 923);
    const auto *skl0_924 = buffer.data(skl0 + 924);
    const auto *skl0_925 = buffer.data(skl0 + 925);
    const auto *skl0_927 = buffer.data(skl0 + 927);

    const auto *skk_640 = buffer.data(skk + 640);
    const auto *skk_647 = buffer.data(skk + 647);
    const auto *skk_648 = buffer.data(skk + 648);
    const auto *skk_651 = buffer.data(skk + 651);
    const auto *skk_654 = buffer.data(skk + 654);
    const auto *skk_658 = buffer.data(skk + 658);
    const auto *skk_663 = buffer.data(skk + 663);
    const auto *skk_668 = buffer.data(skk + 668);
    const auto *skk_676 = buffer.data(skk + 676);
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
    const auto *skk_714 = buffer.data(skk + 714);
    const auto *skk_715 = buffer.data(skk + 715);
    const auto *skk_716 = buffer.data(skk + 716);
    const auto *skk_717 = buffer.data(skk + 717);
    const auto *skk_718 = buffer.data(skk + 718);
    const auto *skk_719 = buffer.data(skk + 719);
    const auto *skk_720 = buffer.data(skk + 720);
    const auto *skk_721 = buffer.data(skk + 721);
    const auto *skk_722 = buffer.data(skk + 722);
    const auto *skk_723 = buffer.data(skk + 723);
    const auto *skk_725 = buffer.data(skk + 725);
    const auto *skk_726 = buffer.data(skk + 726);
    const auto *skk_728 = buffer.data(skk + 728);
    const auto *skk_729 = buffer.data(skk + 729);
    const auto *skk_730 = buffer.data(skk + 730);
    const auto *skk_732 = buffer.data(skk + 732);
    const auto *skk_733 = buffer.data(skk + 733);
    const auto *skk_734 = buffer.data(skk + 734);
    const auto *skk_735 = buffer.data(skk + 735);
    const auto *skk_737 = buffer.data(skk + 737);
    const auto *skk_738 = buffer.data(skk + 738);
    const auto *skk_739 = buffer.data(skk + 739);
    const auto *skk_740 = buffer.data(skk + 740);
    const auto *skk_748 = buffer.data(skk + 748);
    const auto *skk_750 = buffer.data(skk + 750);
    const auto *skk_751 = buffer.data(skk + 751);
    const auto *skk_752 = buffer.data(skk + 752);
    const auto *skk_887 = buffer.data(skk + 887);
    const auto *skk_888 = buffer.data(skk + 888);
    const auto *skk_889 = buffer.data(skk + 889);
    const auto *skk_891 = buffer.data(skk + 891);
    const auto *skk_892 = buffer.data(skk + 892);
    const auto *skk_893 = buffer.data(skk + 893);
    const auto *skk_894 = buffer.data(skk + 894);
    const auto *skk_895 = buffer.data(skk + 895);
    const auto *skk_896 = buffer.data(skk + 896);
    const auto *skk_897 = buffer.data(skk + 897);
    const auto *skk_898 = buffer.data(skk + 898);
    const auto *skk_899 = buffer.data(skk + 899);
    const auto *skk_900 = buffer.data(skk + 900);
    const auto *skk_903 = buffer.data(skk + 903);
    const auto *skk_905 = buffer.data(skk + 905);
    const auto *skk_906 = buffer.data(skk + 906);
    const auto *skk_909 = buffer.data(skk + 909);
    const auto *skk_910 = buffer.data(skk + 910);
    const auto *skk_912 = buffer.data(skk + 912);
    const auto *skk_914 = buffer.data(skk + 914);
    const auto *skk_915 = buffer.data(skk + 915);
    const auto *skk_917 = buffer.data(skk + 917);
    const auto *skk_918 = buffer.data(skk + 918);
    const auto *skk_920 = buffer.data(skk + 920);
    const auto *skk_921 = buffer.data(skk + 921);
    const auto *skk_923 = buffer.data(skk + 923);
    const auto *skk_924 = buffer.data(skk + 924);
    const auto *skk_925 = buffer.data(skk + 925);
    const auto *skk_927 = buffer.data(skk + 927);
    const auto *skk_928 = buffer.data(skk + 928);
    const auto *skk_929 = buffer.data(skk + 929);
    const auto *skk_930 = buffer.data(skk + 930);
    const auto *skk_931 = buffer.data(skk + 931);
    const auto *skk_932 = buffer.data(skk + 932);
    const auto *skk_933 = buffer.data(skk + 933);
    const auto *skk_934 = buffer.data(skk + 934);
    const auto *skk_935 = buffer.data(skk + 935);
    const auto *skk_964 = buffer.data(skk + 964);
    const auto *skk_965 = buffer.data(skk + 965);
    const auto *skk_966 = buffer.data(skk + 966);
    const auto *skk_967 = buffer.data(skk + 967);
    const auto *skk_968 = buffer.data(skk + 968);
    const auto *skk_969 = buffer.data(skk + 969);
    const auto *skk_970 = buffer.data(skk + 970);
    const auto *skk_971 = buffer.data(skk + 971);

    const auto *skl1_900 = buffer.data(skl1 + 900);
    const auto *skl1_903 = buffer.data(skl1 + 903);
    const auto *skl1_905 = buffer.data(skl1 + 905);
    const auto *skl1_906 = buffer.data(skl1 + 906);
    const auto *skl1_909 = buffer.data(skl1 + 909);
    const auto *skl1_910 = buffer.data(skl1 + 910);
    const auto *skl1_912 = buffer.data(skl1 + 912);
    const auto *skl1_914 = buffer.data(skl1 + 914);
    const auto *skl1_915 = buffer.data(skl1 + 915);
    const auto *skl1_917 = buffer.data(skl1 + 917);
    const auto *skl1_918 = buffer.data(skl1 + 918);
    const auto *skl1_920 = buffer.data(skl1 + 920);
    const auto *skl1_921 = buffer.data(skl1 + 921);
    const auto *skl1_923 = buffer.data(skl1 + 923);
    const auto *skl1_924 = buffer.data(skl1 + 924);
    const auto *skl1_925 = buffer.data(skl1 + 925);
    const auto *skl1_927 = buffer.data(skl1 + 927);

    const auto *sli0_693 = buffer.data(sli0 + 693);
    const auto *sli0_695 = buffer.data(sli0 + 695);
    const auto *sli0_696 = buffer.data(sli0 + 696);
    const auto *sli0_697 = buffer.data(sli0 + 697);
    const auto *sli0_698 = buffer.data(sli0 + 698);
    const auto *sli0_699 = buffer.data(sli0 + 699);
    const auto *sli0_700 = buffer.data(sli0 + 700);
    const auto *sli0_703 = buffer.data(sli0 + 703);
    const auto *sli0_705 = buffer.data(sli0 + 705);
    const auto *sli0_706 = buffer.data(sli0 + 706);
    const auto *sli0_709 = buffer.data(sli0 + 709);
    const auto *sli0_710 = buffer.data(sli0 + 710);
    const auto *sli0_712 = buffer.data(sli0 + 712);
    const auto *sli0_714 = buffer.data(sli0 + 714);
    const auto *sli0_715 = buffer.data(sli0 + 715);
    const auto *sli0_717 = buffer.data(sli0 + 717);
    const auto *sli0_718 = buffer.data(sli0 + 718);
    const auto *sli0_720 = buffer.data(sli0 + 720);
    const auto *sli0_721 = buffer.data(sli0 + 721);
    const auto *sli0_723 = buffer.data(sli0 + 723);
    const auto *sli0_724 = buffer.data(sli0 + 724);
    const auto *sli0_725 = buffer.data(sli0 + 725);
    const auto *sli0_726 = buffer.data(sli0 + 726);
    const auto *sli0_727 = buffer.data(sli0 + 727);
    const auto *sli0_749 = buffer.data(sli0 + 749);
    const auto *sli0_751 = buffer.data(sli0 + 751);
    const auto *sli0_752 = buffer.data(sli0 + 752);
    const auto *sli0_753 = buffer.data(sli0 + 753);

    const auto *sli1_693 = buffer.data(sli1 + 693);
    const auto *sli1_695 = buffer.data(sli1 + 695);
    const auto *sli1_696 = buffer.data(sli1 + 696);
    const auto *sli1_697 = buffer.data(sli1 + 697);
    const auto *sli1_698 = buffer.data(sli1 + 698);
    const auto *sli1_699 = buffer.data(sli1 + 699);
    const auto *sli1_700 = buffer.data(sli1 + 700);
    const auto *sli1_703 = buffer.data(sli1 + 703);
    const auto *sli1_705 = buffer.data(sli1 + 705);
    const auto *sli1_706 = buffer.data(sli1 + 706);
    const auto *sli1_709 = buffer.data(sli1 + 709);
    const auto *sli1_710 = buffer.data(sli1 + 710);
    const auto *sli1_712 = buffer.data(sli1 + 712);
    const auto *sli1_714 = buffer.data(sli1 + 714);
    const auto *sli1_715 = buffer.data(sli1 + 715);
    const auto *sli1_717 = buffer.data(sli1 + 717);
    const auto *sli1_718 = buffer.data(sli1 + 718);
    const auto *sli1_720 = buffer.data(sli1 + 720);
    const auto *sli1_721 = buffer.data(sli1 + 721);
    const auto *sli1_723 = buffer.data(sli1 + 723);
    const auto *sli1_724 = buffer.data(sli1 + 724);
    const auto *sli1_725 = buffer.data(sli1 + 725);
    const auto *sli1_726 = buffer.data(sli1 + 726);
    const auto *sli1_727 = buffer.data(sli1 + 727);
    const auto *sli1_749 = buffer.data(sli1 + 749);
    const auto *sli1_751 = buffer.data(sli1 + 751);
    const auto *sli1_752 = buffer.data(sli1 + 752);
    const auto *sli1_753 = buffer.data(sli1 + 753);

    const auto *slk_884 = buffer.data(slk + 884);
    const auto *slk_887 = buffer.data(slk + 887);
    const auto *slk_888 = buffer.data(slk + 888);
    const auto *slk_889 = buffer.data(slk + 889);
    const auto *slk_891 = buffer.data(slk + 891);
    const auto *slk_892 = buffer.data(slk + 892);
    const auto *slk_893 = buffer.data(slk + 893);
    const auto *slk_894 = buffer.data(slk + 894);
    const auto *slk_895 = buffer.data(slk + 895);
    const auto *slk_896 = buffer.data(slk + 896);
    const auto *slk_897 = buffer.data(slk + 897);
    const auto *slk_898 = buffer.data(slk + 898);
    const auto *slk_899 = buffer.data(slk + 899);
    const auto *slk_900 = buffer.data(slk + 900);
    const auto *slk_902 = buffer.data(slk + 902);
    const auto *slk_903 = buffer.data(slk + 903);
    const auto *slk_905 = buffer.data(slk + 905);
    const auto *slk_906 = buffer.data(slk + 906);
    const auto *slk_909 = buffer.data(slk + 909);
    const auto *slk_910 = buffer.data(slk + 910);
    const auto *slk_912 = buffer.data(slk + 912);
    const auto *slk_914 = buffer.data(slk + 914);
    const auto *slk_915 = buffer.data(slk + 915);
    const auto *slk_917 = buffer.data(slk + 917);
    const auto *slk_918 = buffer.data(slk + 918);
    const auto *slk_920 = buffer.data(slk + 920);
    const auto *slk_921 = buffer.data(slk + 921);
    const auto *slk_923 = buffer.data(slk + 923);
    const auto *slk_924 = buffer.data(slk + 924);
    const auto *slk_925 = buffer.data(slk + 925);
    const auto *slk_927 = buffer.data(slk + 927);
    const auto *slk_928 = buffer.data(slk + 928);
    const auto *slk_929 = buffer.data(slk + 929);
    const auto *slk_930 = buffer.data(slk + 930);
    const auto *slk_931 = buffer.data(slk + 931);
    const auto *slk_932 = buffer.data(slk + 932);
    const auto *slk_933 = buffer.data(slk + 933);
    const auto *slk_934 = buffer.data(slk + 934);
    const auto *slk_935 = buffer.data(slk + 935);
    const auto *slk_936 = buffer.data(slk + 936);
    const auto *slk_938 = buffer.data(slk + 938);
    const auto *slk_939 = buffer.data(slk + 939);
    const auto *slk_941 = buffer.data(slk + 941);
    const auto *slk_942 = buffer.data(slk + 942);
    const auto *slk_945 = buffer.data(slk + 945);
    const auto *slk_946 = buffer.data(slk + 946);
    const auto *slk_950 = buffer.data(slk + 950);
    const auto *slk_951 = buffer.data(slk + 951);
    const auto *slk_956 = buffer.data(slk + 956);
    const auto *slk_964 = buffer.data(slk + 964);
    const auto *slk_965 = buffer.data(slk + 965);
    const auto *slk_966 = buffer.data(slk + 966);
    const auto *slk_967 = buffer.data(slk + 967);
    const auto *slk_968 = buffer.data(slk + 968);
    const auto *slk_969 = buffer.data(slk + 969);
    const auto *slk_970 = buffer.data(slk + 970);
    const auto *slk_971 = buffer.data(slk + 971);

#pragma omp simd aligned(t_1103, t_1104, t_1105, pc_x, skk_887, skk_888, skk_889, sli0_695, \
                         sli0_696, sli0_697, sli1_695, sli1_696, sli1_697, slk_887, slk_888, \
                         slk_889 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1103[k] = f_16 * skk_887[k]
                    + f_12 * sli0_695[k]
                    - f_13 * sli1_695[k]
                    + f_3 * pc_x[k] * slk_887[k];

        t_1104[k] = f_16 * skk_888[k]
                    + f_12 * sli0_696[k]
                    - f_13 * sli1_696[k]
                    + f_3 * pc_x[k] * slk_888[k];

        t_1105[k] = f_16 * skk_889[k]
                    + f_12 * sli0_697[k]
                    - f_13 * sli1_697[k]
                    + f_3 * pc_x[k] * slk_889[k];
    }

#pragma omp simd aligned(t_1106, t_1107, t_1108, t_1109, pc_x, pc_y, skk_668, skk_891, \
                         skk_892, skk_893, sli0_699, sli1_699, slk_884, slk_891, slk_892, \
                         slk_893 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1106[k] = f_17 * skk_668[k]
                    + f_3 * pc_y[k] * slk_884[k];

        t_1107[k] = f_16 * skk_891[k]
                    + f_12 * sli0_699[k]
                    - f_13 * sli1_699[k]
                    + f_3 * pc_x[k] * slk_891[k];

        t_1108[k] = f_16 * skk_892[k]
                    + f_3 * pc_x[k] * slk_892[k];

        t_1109[k] = f_16 * skk_893[k]
                    + f_3 * pc_x[k] * slk_893[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, t_1113, t_1114, pc_x, skk_894, skk_895, \
                         skk_896, skk_897, skk_898, slk_894, slk_895, slk_896, slk_897, \
                         slk_898 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_16 * skk_894[k]
                    + f_3 * pc_x[k] * slk_894[k];

        t_1111[k] = f_16 * skk_895[k]
                    + f_3 * pc_x[k] * slk_895[k];

        t_1112[k] = f_16 * skk_896[k]
                    + f_3 * pc_x[k] * slk_896[k];

        t_1113[k] = f_16 * skk_897[k]
                    + f_3 * pc_x[k] * slk_897[k];

        t_1114[k] = f_16 * skk_898[k]
                    + f_3 * pc_x[k] * slk_898[k];
    }

#pragma omp simd aligned(t_1115, t_1116, t_1117, pc_x, pc_y, pc_z, skk_640, skk_676, skk_899, \
                         sli0_693, sli1_693, slk_892, slk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1115[k] = f_16 * skk_899[k]
                    + f_3 * pc_x[k] * slk_899[k];

        t_1116[k] = f_17 * skk_676[k]
                    + f_1 * sli0_693[k]
                    - f_2 * sli1_693[k]
                    + f_3 * pc_y[k] * slk_892[k];

        t_1117[k] = f_17 * skk_640[k]
                    + f_3 * pc_z[k] * slk_892[k];
    }

#pragma omp simd aligned(t_1118, t_1119, t_1120, pc_y, skk_678, skk_679, skk_680, sli0_695, \
                         sli0_696, sli0_697, sli1_695, sli1_696, sli1_697, slk_894, slk_895, \
                         slk_896 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1118[k] = f_17 * skk_678[k]
                    + f_4 * sli0_695[k]
                    - f_5 * sli1_695[k]
                    + f_3 * pc_y[k] * slk_894[k];

        t_1119[k] = f_17 * skk_679[k]
                    + f_6 * sli0_696[k]
                    - f_7 * sli1_696[k]
                    + f_3 * pc_y[k] * slk_895[k];

        t_1120[k] = f_17 * skk_680[k]
                    + f_8 * sli0_697[k]
                    - f_9 * sli1_697[k]
                    + f_3 * pc_y[k] * slk_896[k];
    }

#pragma omp simd aligned(t_1121, t_1122, t_1123, pc_y, skk_681, skk_682, skk_683, sli0_698, \
                         sli0_699, sli1_698, sli1_699, slk_897, slk_898, \
                         slk_899 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1121[k] = f_17 * skk_681[k]
                    + f_10 * sli0_698[k]
                    - f_11 * sli1_698[k]
                    + f_3 * pc_y[k] * slk_897[k];

        t_1122[k] = f_17 * skk_682[k]
                    + f_12 * sli0_699[k]
                    - f_13 * sli1_699[k]
                    + f_3 * pc_y[k] * slk_898[k];

        t_1123[k] = f_17 * skk_683[k]
                    + f_3 * pc_y[k] * slk_899[k];
    }

#pragma omp simd aligned(t_1124, t_1125, t_1126, pc_x, pc_y, pc_z, skk_647, skk_684, skk_900, \
                         sli0_699, sli0_700, sli1_699, sli1_700, slk_899, \
                         slk_900 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1124[k] = f_17 * skk_647[k]
                    + f_1 * sli0_699[k]
                    - f_2 * sli1_699[k]
                    + f_3 * pc_z[k] * slk_899[k];

        t_1125[k] = f_16 * skk_900[k]
                    + f_1 * sli0_700[k]
                    - f_2 * sli1_700[k]
                    + f_3 * pc_x[k] * slk_900[k];

        t_1126[k] = f_16 * skk_684[k]
                    + f_3 * pc_y[k] * slk_900[k];
    }

#pragma omp simd aligned(t_1127, t_1128, t_1129, pc_x, pc_y, pc_z, skk_648, skk_686, skk_903, \
                         sli0_703, sli1_703, slk_900, slk_902, \
                         slk_903 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1127[k] = f_18 * skk_648[k]
                    + f_3 * pc_z[k] * slk_900[k];

        t_1128[k] = f_16 * skk_903[k]
                    + f_4 * sli0_703[k]
                    - f_5 * sli1_703[k]
                    + f_3 * pc_x[k] * slk_903[k];

        t_1129[k] = f_16 * skk_686[k]
                    + f_3 * pc_y[k] * slk_902[k];
    }

#pragma omp simd aligned(t_1130, t_1131, t_1132, pc_x, pc_z, skk_651, skk_905, skk_906, \
                         sli0_705, sli0_706, sli1_705, sli1_706, slk_903, slk_905, \
                         slk_906 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1130[k] = f_16 * skk_905[k]
                    + f_4 * sli0_705[k]
                    - f_5 * sli1_705[k]
                    + f_3 * pc_x[k] * slk_905[k];

        t_1131[k] = f_16 * skk_906[k]
                    + f_6 * sli0_706[k]
                    - f_7 * sli1_706[k]
                    + f_3 * pc_x[k] * slk_906[k];

        t_1132[k] = f_18 * skk_651[k]
                    + f_3 * pc_z[k] * slk_903[k];
    }

#pragma omp simd aligned(t_1133, t_1134, t_1135, pc_x, pc_y, skk_689, skk_909, skk_910, \
                         sli0_709, sli0_710, sli1_709, sli1_710, slk_905, slk_909, \
                         slk_910 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1133[k] = f_16 * skk_689[k]
                    + f_3 * pc_y[k] * slk_905[k];

        t_1134[k] = f_16 * skk_909[k]
                    + f_6 * sli0_709[k]
                    - f_7 * sli1_709[k]
                    + f_3 * pc_x[k] * slk_909[k];

        t_1135[k] = f_16 * skk_910[k]
                    + f_8 * sli0_710[k]
                    - f_9 * sli1_710[k]
                    + f_3 * pc_x[k] * slk_910[k];
    }

#pragma omp simd aligned(t_1136, t_1137, t_1138, pc_x, pc_y, pc_z, skk_654, skk_693, skk_912, \
                         sli0_712, sli1_712, slk_906, slk_909, \
                         slk_912 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1136[k] = f_18 * skk_654[k]
                    + f_3 * pc_z[k] * slk_906[k];

        t_1137[k] = f_16 * skk_912[k]
                    + f_8 * sli0_712[k]
                    - f_9 * sli1_712[k]
                    + f_3 * pc_x[k] * slk_912[k];

        t_1138[k] = f_16 * skk_693[k]
                    + f_3 * pc_y[k] * slk_909[k];
    }

#pragma omp simd aligned(t_1139, t_1140, t_1141, pc_x, pc_z, skk_658, skk_914, skk_915, \
                         sli0_714, sli0_715, sli1_714, sli1_715, slk_910, slk_914, \
                         slk_915 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1139[k] = f_16 * skk_914[k]
                    + f_8 * sli0_714[k]
                    - f_9 * sli1_714[k]
                    + f_3 * pc_x[k] * slk_914[k];

        t_1140[k] = f_16 * skk_915[k]
                    + f_10 * sli0_715[k]
                    - f_11 * sli1_715[k]
                    + f_3 * pc_x[k] * slk_915[k];

        t_1141[k] = f_18 * skk_658[k]
                    + f_3 * pc_z[k] * slk_910[k];
    }

#pragma omp simd aligned(t_1142, t_1143, t_1144, pc_x, pc_y, skk_698, skk_917, skk_918, \
                         sli0_717, sli0_718, sli1_717, sli1_718, slk_914, slk_917, \
                         slk_918 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1142[k] = f_16 * skk_917[k]
                    + f_10 * sli0_717[k]
                    - f_11 * sli1_717[k]
                    + f_3 * pc_x[k] * slk_917[k];

        t_1143[k] = f_16 * skk_918[k]
                    + f_10 * sli0_718[k]
                    - f_11 * sli1_718[k]
                    + f_3 * pc_x[k] * slk_918[k];

        t_1144[k] = f_16 * skk_698[k]
                    + f_3 * pc_y[k] * slk_914[k];
    }

#pragma omp simd aligned(t_1145, t_1146, t_1147, pc_x, pc_z, skk_663, skk_920, skk_921, \
                         sli0_720, sli0_721, sli1_720, sli1_721, slk_915, slk_920, \
                         slk_921 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1145[k] = f_16 * skk_920[k]
                    + f_10 * sli0_720[k]
                    - f_11 * sli1_720[k]
                    + f_3 * pc_x[k] * slk_920[k];

        t_1146[k] = f_16 * skk_921[k]
                    + f_12 * sli0_721[k]
                    - f_13 * sli1_721[k]
                    + f_3 * pc_x[k] * slk_921[k];

        t_1147[k] = f_18 * skk_663[k]
                    + f_3 * pc_z[k] * slk_915[k];
    }

#pragma omp simd aligned(t_1148, t_1149, t_1150, pc_x, skk_923, skk_924, skk_925, sli0_723, \
                         sli0_724, sli0_725, sli1_723, sli1_724, sli1_725, slk_923, slk_924, \
                         slk_925 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1148[k] = f_16 * skk_923[k]
                    + f_12 * sli0_723[k]
                    - f_13 * sli1_723[k]
                    + f_3 * pc_x[k] * slk_923[k];

        t_1149[k] = f_16 * skk_924[k]
                    + f_12 * sli0_724[k]
                    - f_13 * sli1_724[k]
                    + f_3 * pc_x[k] * slk_924[k];

        t_1150[k] = f_16 * skk_925[k]
                    + f_12 * sli0_725[k]
                    - f_13 * sli1_725[k]
                    + f_3 * pc_x[k] * slk_925[k];
    }

#pragma omp simd aligned(t_1151, t_1152, t_1153, t_1154, pc_x, pc_y, skk_704, skk_927, \
                         skk_928, skk_929, sli0_727, sli1_727, slk_920, slk_927, slk_928, \
                         slk_929 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1151[k] = f_16 * skk_704[k]
                    + f_3 * pc_y[k] * slk_920[k];

        t_1152[k] = f_16 * skk_927[k]
                    + f_12 * sli0_727[k]
                    - f_13 * sli1_727[k]
                    + f_3 * pc_x[k] * slk_927[k];

        t_1153[k] = f_16 * skk_928[k]
                    + f_3 * pc_x[k] * slk_928[k];

        t_1154[k] = f_16 * skk_929[k]
                    + f_3 * pc_x[k] * slk_929[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, skk_930, skk_931, \
                         skk_932, skk_933, skk_934, slk_930, slk_931, slk_932, slk_933, \
                         slk_934 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_16 * skk_930[k]
                    + f_3 * pc_x[k] * slk_930[k];

        t_1156[k] = f_16 * skk_931[k]
                    + f_3 * pc_x[k] * slk_931[k];

        t_1157[k] = f_16 * skk_932[k]
                    + f_3 * pc_x[k] * slk_932[k];

        t_1158[k] = f_16 * skk_933[k]
                    + f_3 * pc_x[k] * slk_933[k];

        t_1159[k] = f_16 * skk_934[k]
                    + f_3 * pc_x[k] * slk_934[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, pc_x, pc_y, pc_z, skk_676, skk_712, skk_935, \
                         sli0_721, sli1_721, slk_928, slk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_16 * skk_935[k]
                    + f_3 * pc_x[k] * slk_935[k];

        t_1161[k] = f_16 * skk_712[k]
                    + f_1 * sli0_721[k]
                    - f_2 * sli1_721[k]
                    + f_3 * pc_y[k] * slk_928[k];

        t_1162[k] = f_18 * skk_676[k]
                    + f_3 * pc_z[k] * slk_928[k];
    }

#pragma omp simd aligned(t_1163, t_1164, t_1165, pc_y, skk_714, skk_715, skk_716, sli0_723, \
                         sli0_724, sli0_725, sli1_723, sli1_724, sli1_725, slk_930, slk_931, \
                         slk_932 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1163[k] = f_16 * skk_714[k]
                    + f_4 * sli0_723[k]
                    - f_5 * sli1_723[k]
                    + f_3 * pc_y[k] * slk_930[k];

        t_1164[k] = f_16 * skk_715[k]
                    + f_6 * sli0_724[k]
                    - f_7 * sli1_724[k]
                    + f_3 * pc_y[k] * slk_931[k];

        t_1165[k] = f_16 * skk_716[k]
                    + f_8 * sli0_725[k]
                    - f_9 * sli1_725[k]
                    + f_3 * pc_y[k] * slk_932[k];
    }

#pragma omp simd aligned(t_1166, t_1167, t_1168, pc_y, skk_717, skk_718, skk_719, sli0_726, \
                         sli0_727, sli1_726, sli1_727, slk_933, slk_934, \
                         slk_935 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1166[k] = f_16 * skk_717[k]
                    + f_10 * sli0_726[k]
                    - f_11 * sli1_726[k]
                    + f_3 * pc_y[k] * slk_933[k];

        t_1167[k] = f_16 * skk_718[k]
                    + f_12 * sli0_727[k]
                    - f_13 * sli1_727[k]
                    + f_3 * pc_y[k] * slk_934[k];

        t_1168[k] = f_16 * skk_719[k]
                    + f_3 * pc_y[k] * slk_935[k];
    }

#pragma omp simd aligned(t_1169, t_1170, t_1171, t_1172, pb_y, pc_y, pc_z, skl0_900, skk_683, \
                         skk_684, skk_720, skl1_900, sli0_727, sli1_727, slk_935, \
                         slk_936 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_18 * skk_683[k]
                    + f_1 * sli0_727[k]
                    - f_2 * sli1_727[k]
                    + f_3 * pc_z[k] * slk_935[k];

        t_1170[k] = pb_y[k] * skl0_900[k]
                    - f_14 * pc_y[k] * skl1_900[k];

        t_1171[k] = f_15 * skk_720[k]
                    + f_3 * pc_y[k] * slk_936[k];

        t_1172[k] = f_19 * skk_684[k]
                    + f_3 * pc_z[k] * slk_936[k];
    }

#pragma omp simd aligned(t_1173, t_1174, t_1175, t_1176, pb_y, pc_y, skl0_903, skl0_905, \
                         skl0_906, skk_721, skk_722, skk_723, skl1_903, skl1_905, skl1_906, \
                         slk_938 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1173[k] = pb_y[k] * skl0_903[k]
                    + f_16 * skk_721[k]
                    - f_14 * pc_y[k] * skl1_903[k];

        t_1174[k] = f_15 * skk_722[k]
                    + f_3 * pc_y[k] * slk_938[k];

        t_1175[k] = pb_y[k] * skl0_905[k]
                    - f_14 * pc_y[k] * skl1_905[k];

        t_1176[k] = pb_y[k] * skl0_906[k]
                    + f_17 * skk_723[k]
                    - f_14 * pc_y[k] * skl1_906[k];
    }

#pragma omp simd aligned(t_1177, t_1178, t_1179, t_1180, pb_y, pc_y, pc_z, skl0_909, skl0_910, \
                         skk_687, skk_725, skk_726, skl1_909, skl1_910, slk_939, \
                         slk_941 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1177[k] = f_19 * skk_687[k]
                    + f_3 * pc_z[k] * slk_939[k];

        t_1178[k] = f_15 * skk_725[k]
                    + f_3 * pc_y[k] * slk_941[k];

        t_1179[k] = pb_y[k] * skl0_909[k]
                    - f_14 * pc_y[k] * skl1_909[k];

        t_1180[k] = pb_y[k] * skl0_910[k]
                    + f_18 * skk_726[k]
                    - f_14 * pc_y[k] * skl1_910[k];
    }

#pragma omp simd aligned(t_1181, t_1182, t_1183, t_1184, pb_y, pc_y, pc_z, skl0_912, skl0_914, \
                         skk_690, skk_728, skk_729, skl1_912, skl1_914, slk_942, \
                         slk_945 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1181[k] = f_19 * skk_690[k]
                    + f_3 * pc_z[k] * slk_942[k];

        t_1182[k] = pb_y[k] * skl0_912[k]
                    + f_16 * skk_728[k]
                    - f_14 * pc_y[k] * skl1_912[k];

        t_1183[k] = f_15 * skk_729[k]
                    + f_3 * pc_y[k] * slk_945[k];

        t_1184[k] = pb_y[k] * skl0_914[k]
                    - f_14 * pc_y[k] * skl1_914[k];
    }

#pragma omp simd aligned(t_1185, t_1186, t_1187, pb_y, pc_y, pc_z, skl0_915, skl0_917, \
                         skk_694, skk_730, skk_732, skl1_915, skl1_917, \
                         slk_946 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1185[k] = pb_y[k] * skl0_915[k]
                    + f_19 * skk_730[k]
                    - f_14 * pc_y[k] * skl1_915[k];

        t_1186[k] = f_19 * skk_694[k]
                    + f_3 * pc_z[k] * slk_946[k];

        t_1187[k] = pb_y[k] * skl0_917[k]
                    + f_17 * skk_732[k]
                    - f_14 * pc_y[k] * skl1_917[k];
    }

#pragma omp simd aligned(t_1188, t_1189, t_1190, t_1191, pb_y, pc_y, skl0_918, skl0_920, \
                         skl0_921, skk_733, skk_734, skk_735, skl1_918, skl1_920, skl1_921, \
                         slk_950 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1188[k] = pb_y[k] * skl0_918[k]
                    + f_16 * skk_733[k]
                    - f_14 * pc_y[k] * skl1_918[k];

        t_1189[k] = f_15 * skk_734[k]
                    + f_3 * pc_y[k] * slk_950[k];

        t_1190[k] = pb_y[k] * skl0_920[k]
                    - f_14 * pc_y[k] * skl1_920[k];

        t_1191[k] = pb_y[k] * skl0_921[k]
                    + f_20 * skk_735[k]
                    - f_14 * pc_y[k] * skl1_921[k];
    }

#pragma omp simd aligned(t_1192, t_1193, t_1194, pb_y, pc_y, pc_z, skl0_923, skl0_924, \
                         skk_699, skk_737, skk_738, skl1_923, skl1_924, \
                         slk_951 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1192[k] = f_19 * skk_699[k]
                    + f_3 * pc_z[k] * slk_951[k];

        t_1193[k] = pb_y[k] * skl0_923[k]
                    + f_18 * skk_737[k]
                    - f_14 * pc_y[k] * skl1_923[k];

        t_1194[k] = pb_y[k] * skl0_924[k]
                    + f_17 * skk_738[k]
                    - f_14 * pc_y[k] * skl1_924[k];
    }

#pragma omp simd aligned(t_1195, t_1196, t_1197, t_1198, pb_y, pc_x, pc_y, skl0_925, skl0_927, \
                         skk_739, skk_740, skk_964, skl1_925, skl1_927, slk_956, \
                         slk_964 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1195[k] = pb_y[k] * skl0_925[k]
                    + f_16 * skk_739[k]
                    - f_14 * pc_y[k] * skl1_925[k];

        t_1196[k] = f_15 * skk_740[k]
                    + f_3 * pc_y[k] * slk_956[k];

        t_1197[k] = pb_y[k] * skl0_927[k]
                    - f_14 * pc_y[k] * skl1_927[k];

        t_1198[k] = f_16 * skk_964[k]
                    + f_3 * pc_x[k] * slk_964[k];
    }

#pragma omp simd aligned(t_1199, t_1200, t_1201, t_1202, t_1203, pc_x, skk_965, skk_966, \
                         skk_967, skk_968, skk_969, slk_965, slk_966, slk_967, slk_968, \
                         slk_969 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1199[k] = f_16 * skk_965[k]
                    + f_3 * pc_x[k] * slk_965[k];

        t_1200[k] = f_16 * skk_966[k]
                    + f_3 * pc_x[k] * slk_966[k];

        t_1201[k] = f_16 * skk_967[k]
                    + f_3 * pc_x[k] * slk_967[k];

        t_1202[k] = f_16 * skk_968[k]
                    + f_3 * pc_x[k] * slk_968[k];

        t_1203[k] = f_16 * skk_969[k]
                    + f_3 * pc_x[k] * slk_969[k];
    }

#pragma omp simd aligned(t_1204, t_1205, t_1206, t_1207, pc_x, pc_y, pc_z, skk_712, skk_748, \
                         skk_970, skk_971, sli0_749, sli1_749, slk_964, slk_970, \
                         slk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1204[k] = f_16 * skk_970[k]
                    + f_3 * pc_x[k] * slk_970[k];

        t_1205[k] = f_16 * skk_971[k]
                    + f_3 * pc_x[k] * slk_971[k];

        t_1206[k] = f_15 * skk_748[k]
                    + f_1 * sli0_749[k]
                    - f_2 * sli1_749[k]
                    + f_3 * pc_y[k] * slk_964[k];

        t_1207[k] = f_19 * skk_712[k]
                    + f_3 * pc_z[k] * slk_964[k];
    }

#pragma omp simd aligned(t_1208, t_1209, t_1210, pc_y, skk_750, skk_751, skk_752, sli0_751, \
                         sli0_752, sli0_753, sli1_751, sli1_752, sli1_753, slk_966, slk_967, \
                         slk_968 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1208[k] = f_15 * skk_750[k]
                    + f_4 * sli0_751[k]
                    - f_5 * sli1_751[k]
                    + f_3 * pc_y[k] * slk_966[k];

        t_1209[k] = f_15 * skk_751[k]
                    + f_6 * sli0_752[k]
                    - f_7 * sli1_752[k]
                    + f_3 * pc_y[k] * slk_967[k];

        t_1210[k] = f_15 * skk_752[k]
                    + f_8 * sli0_753[k]
                    - f_9 * sli1_753[k]
                    + f_3 * pc_y[k] * slk_968[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece11(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t sli0, const size_t sli1,
                                                           const size_t slk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_21 = 3.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_944 = buffer.data(skl0 + 944);
    const auto *skl0_945 = buffer.data(skl0 + 945);
    const auto *skl0_948 = buffer.data(skl0 + 948);
    const auto *skl0_951 = buffer.data(skl0 + 951);
    const auto *skl0_955 = buffer.data(skl0 + 955);
    const auto *skl0_960 = buffer.data(skl0 + 960);
    const auto *skl0_1260 = buffer.data(skl0 + 1260);
    const auto *skl0_1263 = buffer.data(skl0 + 1263);
    const auto *skl0_1265 = buffer.data(skl0 + 1265);
    const auto *skl0_1266 = buffer.data(skl0 + 1266);
    const auto *skl0_1269 = buffer.data(skl0 + 1269);
    const auto *skl0_1270 = buffer.data(skl0 + 1270);
    const auto *skl0_1272 = buffer.data(skl0 + 1272);
    const auto *skl0_1274 = buffer.data(skl0 + 1274);
    const auto *skl0_1275 = buffer.data(skl0 + 1275);
    const auto *skl0_1277 = buffer.data(skl0 + 1277);
    const auto *skl0_1278 = buffer.data(skl0 + 1278);
    const auto *skl0_1280 = buffer.data(skl0 + 1280);
    const auto *skl0_1281 = buffer.data(skl0 + 1281);
    const auto *skl0_1283 = buffer.data(skl0 + 1283);
    const auto *skl0_1284 = buffer.data(skl0 + 1284);
    const auto *skl0_1285 = buffer.data(skl0 + 1285);
    const auto *skl0_1287 = buffer.data(skl0 + 1287);
    const auto *skl0_1296 = buffer.data(skl0 + 1296);
    const auto *skl0_1298 = buffer.data(skl0 + 1298);
    const auto *skl0_1299 = buffer.data(skl0 + 1299);
    const auto *skl0_1300 = buffer.data(skl0 + 1300);
    const auto *skl0_1301 = buffer.data(skl0 + 1301);
    const auto *skl0_1302 = buffer.data(skl0 + 1302);
    const auto *skl0_1304 = buffer.data(skl0 + 1304);
    const auto *skl0_1310 = buffer.data(skl0 + 1310);
    const auto *skl0_1314 = buffer.data(skl0 + 1314);
    const auto *skl0_1317 = buffer.data(skl0 + 1317);
    const auto *skl0_1319 = buffer.data(skl0 + 1319);
    const auto *skl0_1322 = buffer.data(skl0 + 1322);
    const auto *skl0_1323 = buffer.data(skl0 + 1323);

    const auto *skk_720 = buffer.data(skk + 720);
    const auto *skk_723 = buffer.data(skk + 723);
    const auto *skk_726 = buffer.data(skk + 726);
    const auto *skk_730 = buffer.data(skk + 730);
    const auto *skk_735 = buffer.data(skk + 735);
    const auto *skk_748 = buffer.data(skk + 748);
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
    const auto *skk_776 = buffer.data(skk + 776);
    const auto *skk_791 = buffer.data(skk + 791);
    const auto *skk_792 = buffer.data(skk + 792);
    const auto *skk_794 = buffer.data(skk + 794);
    const auto *skk_797 = buffer.data(skk + 797);
    const auto *skk_801 = buffer.data(skk + 801);
    const auto *skk_806 = buffer.data(skk + 806);
    const auto *skk_972 = buffer.data(skk + 972);
    const auto *skk_975 = buffer.data(skk + 975);
    const auto *skk_977 = buffer.data(skk + 977);
    const auto *skk_978 = buffer.data(skk + 978);
    const auto *skk_981 = buffer.data(skk + 981);
    const auto *skk_982 = buffer.data(skk + 982);
    const auto *skk_984 = buffer.data(skk + 984);
    const auto *skk_986 = buffer.data(skk + 986);
    const auto *skk_987 = buffer.data(skk + 987);
    const auto *skk_989 = buffer.data(skk + 989);
    const auto *skk_990 = buffer.data(skk + 990);
    const auto *skk_992 = buffer.data(skk + 992);
    const auto *skk_993 = buffer.data(skk + 993);
    const auto *skk_995 = buffer.data(skk + 995);
    const auto *skk_996 = buffer.data(skk + 996);
    const auto *skk_997 = buffer.data(skk + 997);
    const auto *skk_999 = buffer.data(skk + 999);
    const auto *skk_1000 = buffer.data(skk + 1000);
    const auto *skk_1001 = buffer.data(skk + 1001);
    const auto *skk_1002 = buffer.data(skk + 1002);
    const auto *skk_1003 = buffer.data(skk + 1003);
    const auto *skk_1004 = buffer.data(skk + 1004);
    const auto *skk_1005 = buffer.data(skk + 1005);
    const auto *skk_1006 = buffer.data(skk + 1006);
    const auto *skk_1007 = buffer.data(skk + 1007);
    const auto *skk_1008 = buffer.data(skk + 1008);
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
    const auto *skk_1049 = buffer.data(skk + 1049);
    const auto *skk_1053 = buffer.data(skk + 1053);
    const auto *skk_1056 = buffer.data(skk + 1056);
    const auto *skk_1058 = buffer.data(skk + 1058);
    const auto *skk_1061 = buffer.data(skk + 1061);
    const auto *skk_1062 = buffer.data(skk + 1062);

    const auto *skl1_944 = buffer.data(skl1 + 944);
    const auto *skl1_945 = buffer.data(skl1 + 945);
    const auto *skl1_948 = buffer.data(skl1 + 948);
    const auto *skl1_951 = buffer.data(skl1 + 951);
    const auto *skl1_955 = buffer.data(skl1 + 955);
    const auto *skl1_960 = buffer.data(skl1 + 960);
    const auto *skl1_1260 = buffer.data(skl1 + 1260);
    const auto *skl1_1263 = buffer.data(skl1 + 1263);
    const auto *skl1_1265 = buffer.data(skl1 + 1265);
    const auto *skl1_1266 = buffer.data(skl1 + 1266);
    const auto *skl1_1269 = buffer.data(skl1 + 1269);
    const auto *skl1_1270 = buffer.data(skl1 + 1270);
    const auto *skl1_1272 = buffer.data(skl1 + 1272);
    const auto *skl1_1274 = buffer.data(skl1 + 1274);
    const auto *skl1_1275 = buffer.data(skl1 + 1275);
    const auto *skl1_1277 = buffer.data(skl1 + 1277);
    const auto *skl1_1278 = buffer.data(skl1 + 1278);
    const auto *skl1_1280 = buffer.data(skl1 + 1280);
    const auto *skl1_1281 = buffer.data(skl1 + 1281);
    const auto *skl1_1283 = buffer.data(skl1 + 1283);
    const auto *skl1_1284 = buffer.data(skl1 + 1284);
    const auto *skl1_1285 = buffer.data(skl1 + 1285);
    const auto *skl1_1287 = buffer.data(skl1 + 1287);
    const auto *skl1_1296 = buffer.data(skl1 + 1296);
    const auto *skl1_1298 = buffer.data(skl1 + 1298);
    const auto *skl1_1299 = buffer.data(skl1 + 1299);
    const auto *skl1_1300 = buffer.data(skl1 + 1300);
    const auto *skl1_1301 = buffer.data(skl1 + 1301);
    const auto *skl1_1302 = buffer.data(skl1 + 1302);
    const auto *skl1_1304 = buffer.data(skl1 + 1304);
    const auto *skl1_1310 = buffer.data(skl1 + 1310);
    const auto *skl1_1314 = buffer.data(skl1 + 1314);
    const auto *skl1_1317 = buffer.data(skl1 + 1317);
    const auto *skl1_1319 = buffer.data(skl1 + 1319);
    const auto *skl1_1322 = buffer.data(skl1 + 1322);
    const auto *skl1_1323 = buffer.data(skl1 + 1323);

    const auto *sli0_754 = buffer.data(sli0 + 754);
    const auto *sli0_755 = buffer.data(sli0 + 755);
    const auto *sli0_756 = buffer.data(sli0 + 756);
    const auto *sli0_759 = buffer.data(sli0 + 759);
    const auto *sli0_761 = buffer.data(sli0 + 761);
    const auto *sli0_762 = buffer.data(sli0 + 762);
    const auto *sli0_765 = buffer.data(sli0 + 765);
    const auto *sli0_766 = buffer.data(sli0 + 766);
    const auto *sli0_768 = buffer.data(sli0 + 768);
    const auto *sli0_770 = buffer.data(sli0 + 770);
    const auto *sli0_771 = buffer.data(sli0 + 771);
    const auto *sli0_773 = buffer.data(sli0 + 773);
    const auto *sli0_774 = buffer.data(sli0 + 774);
    const auto *sli0_776 = buffer.data(sli0 + 776);
    const auto *sli0_777 = buffer.data(sli0 + 777);
    const auto *sli0_779 = buffer.data(sli0 + 779);
    const auto *sli0_780 = buffer.data(sli0 + 780);
    const auto *sli0_781 = buffer.data(sli0 + 781);
    const auto *sli0_782 = buffer.data(sli0 + 782);
    const auto *sli0_783 = buffer.data(sli0 + 783);

    const auto *sli1_754 = buffer.data(sli1 + 754);
    const auto *sli1_755 = buffer.data(sli1 + 755);
    const auto *sli1_756 = buffer.data(sli1 + 756);
    const auto *sli1_759 = buffer.data(sli1 + 759);
    const auto *sli1_761 = buffer.data(sli1 + 761);
    const auto *sli1_762 = buffer.data(sli1 + 762);
    const auto *sli1_765 = buffer.data(sli1 + 765);
    const auto *sli1_766 = buffer.data(sli1 + 766);
    const auto *sli1_768 = buffer.data(sli1 + 768);
    const auto *sli1_770 = buffer.data(sli1 + 770);
    const auto *sli1_771 = buffer.data(sli1 + 771);
    const auto *sli1_773 = buffer.data(sli1 + 773);
    const auto *sli1_774 = buffer.data(sli1 + 774);
    const auto *sli1_776 = buffer.data(sli1 + 776);
    const auto *sli1_777 = buffer.data(sli1 + 777);
    const auto *sli1_779 = buffer.data(sli1 + 779);
    const auto *sli1_780 = buffer.data(sli1 + 780);
    const auto *sli1_781 = buffer.data(sli1 + 781);
    const auto *sli1_782 = buffer.data(sli1 + 782);
    const auto *sli1_783 = buffer.data(sli1 + 783);

    const auto *slk_969 = buffer.data(slk + 969);
    const auto *slk_970 = buffer.data(slk + 970);
    const auto *slk_971 = buffer.data(slk + 971);
    const auto *slk_972 = buffer.data(slk + 972);
    const auto *slk_974 = buffer.data(slk + 974);
    const auto *slk_975 = buffer.data(slk + 975);
    const auto *slk_977 = buffer.data(slk + 977);
    const auto *slk_978 = buffer.data(slk + 978);
    const auto *slk_981 = buffer.data(slk + 981);
    const auto *slk_982 = buffer.data(slk + 982);
    const auto *slk_984 = buffer.data(slk + 984);
    const auto *slk_986 = buffer.data(slk + 986);
    const auto *slk_987 = buffer.data(slk + 987);
    const auto *slk_989 = buffer.data(slk + 989);
    const auto *slk_990 = buffer.data(slk + 990);
    const auto *slk_992 = buffer.data(slk + 992);
    const auto *slk_993 = buffer.data(slk + 993);
    const auto *slk_995 = buffer.data(slk + 995);
    const auto *slk_996 = buffer.data(slk + 996);
    const auto *slk_997 = buffer.data(slk + 997);
    const auto *slk_999 = buffer.data(slk + 999);
    const auto *slk_1000 = buffer.data(slk + 1000);
    const auto *slk_1001 = buffer.data(slk + 1001);
    const auto *slk_1002 = buffer.data(slk + 1002);
    const auto *slk_1003 = buffer.data(slk + 1003);
    const auto *slk_1004 = buffer.data(slk + 1004);
    const auto *slk_1005 = buffer.data(slk + 1005);
    const auto *slk_1006 = buffer.data(slk + 1006);
    const auto *slk_1007 = buffer.data(slk + 1007);
    const auto *slk_1008 = buffer.data(slk + 1008);
    const auto *slk_1010 = buffer.data(slk + 1010);
    const auto *slk_1011 = buffer.data(slk + 1011);
    const auto *slk_1013 = buffer.data(slk + 1013);
    const auto *slk_1014 = buffer.data(slk + 1014);
    const auto *slk_1017 = buffer.data(slk + 1017);
    const auto *slk_1018 = buffer.data(slk + 1018);
    const auto *slk_1022 = buffer.data(slk + 1022);
    const auto *slk_1023 = buffer.data(slk + 1023);
    const auto *slk_1028 = buffer.data(slk + 1028);
    const auto *slk_1036 = buffer.data(slk + 1036);
    const auto *slk_1037 = buffer.data(slk + 1037);
    const auto *slk_1038 = buffer.data(slk + 1038);
    const auto *slk_1039 = buffer.data(slk + 1039);
    const auto *slk_1040 = buffer.data(slk + 1040);
    const auto *slk_1041 = buffer.data(slk + 1041);
    const auto *slk_1042 = buffer.data(slk + 1042);
    const auto *slk_1043 = buffer.data(slk + 1043);
    const auto *slk_1044 = buffer.data(slk + 1044);
    const auto *slk_1046 = buffer.data(slk + 1046);
    const auto *slk_1047 = buffer.data(slk + 1047);
    const auto *slk_1049 = buffer.data(slk + 1049);
    const auto *slk_1050 = buffer.data(slk + 1050);
    const auto *slk_1053 = buffer.data(slk + 1053);
    const auto *slk_1054 = buffer.data(slk + 1054);
    const auto *slk_1058 = buffer.data(slk + 1058);

#pragma omp simd aligned(t_1211, t_1212, t_1213, pc_y, skk_753, skk_754, skk_755, sli0_754, \
                         sli0_755, sli1_754, sli1_755, slk_969, slk_970, \
                         slk_971 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1211[k] = f_15 * skk_753[k]
                    + f_10 * sli0_754[k]
                    - f_11 * sli1_754[k]
                    + f_3 * pc_y[k] * slk_969[k];

        t_1212[k] = f_15 * skk_754[k]
                    + f_12 * sli0_755[k]
                    - f_13 * sli1_755[k]
                    + f_3 * pc_y[k] * slk_970[k];

        t_1213[k] = f_15 * skk_755[k]
                    + f_3 * pc_y[k] * slk_971[k];
    }

#pragma omp simd aligned(t_1214, t_1215, t_1216, t_1217, pb_y, pc_x, pc_y, pc_z, skl0_944, \
                         skk_720, skk_972, skl1_944, sli0_756, sli1_756, \
                         slk_972 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1214[k] = pb_y[k] * skl0_944[k]
                    - f_14 * pc_y[k] * skl1_944[k];

        t_1215[k] = f_16 * skk_972[k]
                    + f_1 * sli0_756[k]
                    - f_2 * sli1_756[k]
                    + f_3 * pc_x[k] * slk_972[k];

        t_1216[k] = f_3 * pc_y[k] * slk_972[k];

        t_1217[k] = f_20 * skk_720[k]
                    + f_3 * pc_z[k] * slk_972[k];
    }

#pragma omp simd aligned(t_1218, t_1219, t_1220, pc_x, pc_y, skk_975, skk_977, sli0_759, \
                         sli0_761, sli1_759, sli1_761, slk_974, slk_975, \
                         slk_977 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1218[k] = f_16 * skk_975[k]
                    + f_4 * sli0_759[k]
                    - f_5 * sli1_759[k]
                    + f_3 * pc_x[k] * slk_975[k];

        t_1219[k] = f_3 * pc_y[k] * slk_974[k];

        t_1220[k] = f_16 * skk_977[k]
                    + f_4 * sli0_761[k]
                    - f_5 * sli1_761[k]
                    + f_3 * pc_x[k] * slk_977[k];
    }

#pragma omp simd aligned(t_1221, t_1222, t_1223, pc_x, pc_y, pc_z, skk_723, skk_978, sli0_762, \
                         sli1_762, slk_975, slk_977, slk_978 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1221[k] = f_16 * skk_978[k]
                    + f_6 * sli0_762[k]
                    - f_7 * sli1_762[k]
                    + f_3 * pc_x[k] * slk_978[k];

        t_1222[k] = f_20 * skk_723[k]
                    + f_3 * pc_z[k] * slk_975[k];

        t_1223[k] = f_3 * pc_y[k] * slk_977[k];
    }

#pragma omp simd aligned(t_1224, t_1225, t_1226, pc_x, pc_z, skk_726, skk_981, skk_982, \
                         sli0_765, sli0_766, sli1_765, sli1_766, slk_978, slk_981, \
                         slk_982 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1224[k] = f_16 * skk_981[k]
                    + f_6 * sli0_765[k]
                    - f_7 * sli1_765[k]
                    + f_3 * pc_x[k] * slk_981[k];

        t_1225[k] = f_16 * skk_982[k]
                    + f_8 * sli0_766[k]
                    - f_9 * sli1_766[k]
                    + f_3 * pc_x[k] * slk_982[k];

        t_1226[k] = f_20 * skk_726[k]
                    + f_3 * pc_z[k] * slk_978[k];
    }

#pragma omp simd aligned(t_1227, t_1228, t_1229, pc_x, pc_y, skk_984, skk_986, sli0_768, \
                         sli0_770, sli1_768, sli1_770, slk_981, slk_984, \
                         slk_986 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1227[k] = f_16 * skk_984[k]
                    + f_8 * sli0_768[k]
                    - f_9 * sli1_768[k]
                    + f_3 * pc_x[k] * slk_984[k];

        t_1228[k] = f_3 * pc_y[k] * slk_981[k];

        t_1229[k] = f_16 * skk_986[k]
                    + f_8 * sli0_770[k]
                    - f_9 * sli1_770[k]
                    + f_3 * pc_x[k] * slk_986[k];
    }

#pragma omp simd aligned(t_1230, t_1231, t_1232, pc_x, pc_z, skk_730, skk_987, skk_989, \
                         sli0_771, sli0_773, sli1_771, sli1_773, slk_982, slk_987, \
                         slk_989 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1230[k] = f_16 * skk_987[k]
                    + f_10 * sli0_771[k]
                    - f_11 * sli1_771[k]
                    + f_3 * pc_x[k] * slk_987[k];

        t_1231[k] = f_20 * skk_730[k]
                    + f_3 * pc_z[k] * slk_982[k];

        t_1232[k] = f_16 * skk_989[k]
                    + f_10 * sli0_773[k]
                    - f_11 * sli1_773[k]
                    + f_3 * pc_x[k] * slk_989[k];
    }

#pragma omp simd aligned(t_1233, t_1234, t_1235, pc_x, pc_y, skk_990, skk_992, sli0_774, \
                         sli0_776, sli1_774, sli1_776, slk_986, slk_990, \
                         slk_992 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1233[k] = f_16 * skk_990[k]
                    + f_10 * sli0_774[k]
                    - f_11 * sli1_774[k]
                    + f_3 * pc_x[k] * slk_990[k];

        t_1234[k] = f_3 * pc_y[k] * slk_986[k];

        t_1235[k] = f_16 * skk_992[k]
                    + f_10 * sli0_776[k]
                    - f_11 * sli1_776[k]
                    + f_3 * pc_x[k] * slk_992[k];
    }

#pragma omp simd aligned(t_1236, t_1237, t_1238, pc_x, pc_z, skk_735, skk_993, skk_995, \
                         sli0_777, sli0_779, sli1_777, sli1_779, slk_987, slk_993, \
                         slk_995 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1236[k] = f_16 * skk_993[k]
                    + f_12 * sli0_777[k]
                    - f_13 * sli1_777[k]
                    + f_3 * pc_x[k] * slk_993[k];

        t_1237[k] = f_20 * skk_735[k]
                    + f_3 * pc_z[k] * slk_987[k];

        t_1238[k] = f_16 * skk_995[k]
                    + f_12 * sli0_779[k]
                    - f_13 * sli1_779[k]
                    + f_3 * pc_x[k] * slk_995[k];
    }

#pragma omp simd aligned(t_1239, t_1240, t_1241, pc_x, pc_y, skk_996, skk_997, sli0_780, \
                         sli0_781, sli1_780, sli1_781, slk_992, slk_996, \
                         slk_997 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1239[k] = f_16 * skk_996[k]
                    + f_12 * sli0_780[k]
                    - f_13 * sli1_780[k]
                    + f_3 * pc_x[k] * slk_996[k];

        t_1240[k] = f_16 * skk_997[k]
                    + f_12 * sli0_781[k]
                    - f_13 * sli1_781[k]
                    + f_3 * pc_x[k] * slk_997[k];

        t_1241[k] = f_3 * pc_y[k] * slk_992[k];
    }

#pragma omp simd aligned(t_1242, t_1243, t_1244, t_1245, pc_x, skk_999, skk_1000, skk_1001, \
                         skk_1002, sli0_783, sli1_783, slk_999, slk_1000, slk_1001, \
                         slk_1002 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1242[k] = f_16 * skk_999[k]
                    + f_12 * sli0_783[k]
                    - f_13 * sli1_783[k]
                    + f_3 * pc_x[k] * slk_999[k];

        t_1243[k] = f_16 * skk_1000[k]
                    + f_3 * pc_x[k] * slk_1000[k];

        t_1244[k] = f_16 * skk_1001[k]
                    + f_3 * pc_x[k] * slk_1001[k];

        t_1245[k] = f_16 * skk_1002[k]
                    + f_3 * pc_x[k] * slk_1002[k];
    }

#pragma omp simd aligned(t_1246, t_1247, t_1248, t_1249, t_1250, pc_x, skk_1003, skk_1004, \
                         skk_1005, skk_1006, skk_1007, slk_1003, slk_1004, slk_1005, slk_1006, \
                         slk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1246[k] = f_16 * skk_1003[k]
                    + f_3 * pc_x[k] * slk_1003[k];

        t_1247[k] = f_16 * skk_1004[k]
                    + f_3 * pc_x[k] * slk_1004[k];

        t_1248[k] = f_16 * skk_1005[k]
                    + f_3 * pc_x[k] * slk_1005[k];

        t_1249[k] = f_16 * skk_1006[k]
                    + f_3 * pc_x[k] * slk_1006[k];

        t_1250[k] = f_16 * skk_1007[k]
                    + f_3 * pc_x[k] * slk_1007[k];
    }

#pragma omp simd aligned(t_1251, t_1252, t_1253, t_1254, pc_y, pc_z, skk_748, sli0_777, \
                         sli0_779, sli0_780, sli1_777, sli1_779, sli1_780, slk_1000, slk_1002, \
                         slk_1003 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1251[k] = f_1 * sli0_777[k]
                    - f_2 * sli1_777[k]
                    + f_3 * pc_y[k] * slk_1000[k];

        t_1252[k] = f_20 * skk_748[k]
                    + f_3 * pc_z[k] * slk_1000[k];

        t_1253[k] = f_4 * sli0_779[k]
                    - f_5 * sli1_779[k]
                    + f_3 * pc_y[k] * slk_1002[k];

        t_1254[k] = f_6 * sli0_780[k]
                    - f_7 * sli1_780[k]
                    + f_3 * pc_y[k] * slk_1003[k];
    }

#pragma omp simd aligned(t_1255, t_1256, t_1257, t_1258, pc_y, sli0_781, sli0_782, sli0_783, \
                         sli1_781, sli1_782, sli1_783, slk_1004, slk_1005, slk_1006, \
                         slk_1007 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1255[k] = f_8 * sli0_781[k]
                    - f_9 * sli1_781[k]
                    + f_3 * pc_y[k] * slk_1004[k];

        t_1256[k] = f_10 * sli0_782[k]
                    - f_11 * sli1_782[k]
                    + f_3 * pc_y[k] * slk_1005[k];

        t_1257[k] = f_12 * sli0_783[k]
                    - f_13 * sli1_783[k]
                    + f_3 * pc_y[k] * slk_1006[k];

        t_1258[k] = f_3 * pc_y[k] * slk_1007[k];
    }

#pragma omp simd aligned(t_1259, t_1260, t_1261, pb_x, pc_x, pc_y, pc_z, skl0_1260, skk_755, \
                         skk_756, skk_1008, skl1_1260, sli0_783, sli1_783, slk_1007, \
                         slk_1008 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1259[k] = f_20 * skk_755[k]
                    + f_1 * sli0_783[k]
                    - f_2 * sli1_783[k]
                    + f_3 * pc_z[k] * slk_1007[k];

        t_1260[k] = pb_x[k] * skl0_1260[k]
                    + f_0 * skk_1008[k]
                    - f_14 * pc_x[k] * skl1_1260[k];

        t_1261[k] = f_21 * skk_756[k]
                    + f_3 * pc_y[k] * slk_1008[k];
    }

#pragma omp simd aligned(t_1262, t_1263, t_1264, pb_x, pc_x, pc_y, pc_z, skl0_1263, skk_758, \
                         skk_1011, skl1_1263, slk_1008, slk_1010 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1262[k] = f_3 * pc_z[k] * slk_1008[k];

        t_1263[k] = pb_x[k] * skl0_1263[k]
                    + f_20 * skk_1011[k]
                    - f_14 * pc_x[k] * skl1_1263[k];

        t_1264[k] = f_21 * skk_758[k]
                    + f_3 * pc_y[k] * slk_1010[k];
    }

#pragma omp simd aligned(t_1265, t_1266, t_1267, pb_x, pc_x, pc_z, skl0_1265, skl0_1266, \
                         skk_1013, skk_1014, skl1_1265, skl1_1266, \
                         slk_1011 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1265[k] = pb_x[k] * skl0_1265[k]
                    + f_20 * skk_1013[k]
                    - f_14 * pc_x[k] * skl1_1265[k];

        t_1266[k] = pb_x[k] * skl0_1266[k]
                    + f_19 * skk_1014[k]
                    - f_14 * pc_x[k] * skl1_1266[k];

        t_1267[k] = f_3 * pc_z[k] * slk_1011[k];
    }

#pragma omp simd aligned(t_1268, t_1269, t_1270, pb_x, pc_x, pc_y, skl0_1269, skl0_1270, \
                         skk_761, skk_1017, skk_1018, skl1_1269, skl1_1270, \
                         slk_1013 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1268[k] = f_21 * skk_761[k]
                    + f_3 * pc_y[k] * slk_1013[k];

        t_1269[k] = pb_x[k] * skl0_1269[k]
                    + f_19 * skk_1017[k]
                    - f_14 * pc_x[k] * skl1_1269[k];

        t_1270[k] = pb_x[k] * skl0_1270[k]
                    + f_18 * skk_1018[k]
                    - f_14 * pc_x[k] * skl1_1270[k];
    }

#pragma omp simd aligned(t_1271, t_1272, t_1273, pb_x, pc_x, pc_y, pc_z, skl0_1272, skk_765, \
                         skk_1020, skl1_1272, slk_1014, slk_1017 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1271[k] = f_3 * pc_z[k] * slk_1014[k];

        t_1272[k] = pb_x[k] * skl0_1272[k]
                    + f_18 * skk_1020[k]
                    - f_14 * pc_x[k] * skl1_1272[k];

        t_1273[k] = f_21 * skk_765[k]
                    + f_3 * pc_y[k] * slk_1017[k];
    }

#pragma omp simd aligned(t_1274, t_1275, t_1276, pb_x, pc_x, pc_z, skl0_1274, skl0_1275, \
                         skk_1022, skk_1023, skl1_1274, skl1_1275, \
                         slk_1018 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1274[k] = pb_x[k] * skl0_1274[k]
                    + f_18 * skk_1022[k]
                    - f_14 * pc_x[k] * skl1_1274[k];

        t_1275[k] = pb_x[k] * skl0_1275[k]
                    + f_17 * skk_1023[k]
                    - f_14 * pc_x[k] * skl1_1275[k];

        t_1276[k] = f_3 * pc_z[k] * slk_1018[k];
    }

#pragma omp simd aligned(t_1277, t_1278, t_1279, pb_x, pc_x, pc_y, skl0_1277, skl0_1278, \
                         skk_770, skk_1025, skk_1026, skl1_1277, skl1_1278, \
                         slk_1022 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1277[k] = pb_x[k] * skl0_1277[k]
                    + f_17 * skk_1025[k]
                    - f_14 * pc_x[k] * skl1_1277[k];

        t_1278[k] = pb_x[k] * skl0_1278[k]
                    + f_17 * skk_1026[k]
                    - f_14 * pc_x[k] * skl1_1278[k];

        t_1279[k] = f_21 * skk_770[k]
                    + f_3 * pc_y[k] * slk_1022[k];
    }

#pragma omp simd aligned(t_1280, t_1281, t_1282, pb_x, pc_x, pc_z, skl0_1280, skl0_1281, \
                         skk_1028, skk_1029, skl1_1280, skl1_1281, \
                         slk_1023 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1280[k] = pb_x[k] * skl0_1280[k]
                    + f_17 * skk_1028[k]
                    - f_14 * pc_x[k] * skl1_1280[k];

        t_1281[k] = pb_x[k] * skl0_1281[k]
                    + f_16 * skk_1029[k]
                    - f_14 * pc_x[k] * skl1_1281[k];

        t_1282[k] = f_3 * pc_z[k] * slk_1023[k];
    }

#pragma omp simd aligned(t_1283, t_1284, t_1285, pb_x, pc_x, skl0_1283, skl0_1284, skl0_1285, \
                         skk_1031, skk_1032, skk_1033, skl1_1283, skl1_1284, \
                         skl1_1285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1283[k] = pb_x[k] * skl0_1283[k]
                    + f_16 * skk_1031[k]
                    - f_14 * pc_x[k] * skl1_1283[k];

        t_1284[k] = pb_x[k] * skl0_1284[k]
                    + f_16 * skk_1032[k]
                    - f_14 * pc_x[k] * skl1_1284[k];

        t_1285[k] = pb_x[k] * skl0_1285[k]
                    + f_16 * skk_1033[k]
                    - f_14 * pc_x[k] * skl1_1285[k];
    }

#pragma omp simd aligned(t_1286, t_1287, t_1288, t_1289, pb_x, pc_x, pc_y, skl0_1287, skk_776, \
                         skk_1035, skk_1036, skk_1037, skl1_1287, slk_1028, slk_1036, \
                         slk_1037 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1286[k] = f_21 * skk_776[k]
                    + f_3 * pc_y[k] * slk_1028[k];

        t_1287[k] = pb_x[k] * skl0_1287[k]
                    + f_16 * skk_1035[k]
                    - f_14 * pc_x[k] * skl1_1287[k];

        t_1288[k] = f_15 * skk_1036[k]
                    + f_3 * pc_x[k] * slk_1036[k];

        t_1289[k] = f_15 * skk_1037[k]
                    + f_3 * pc_x[k] * slk_1037[k];
    }

#pragma omp simd aligned(t_1290, t_1291, t_1292, t_1293, t_1294, pc_x, skk_1038, skk_1039, \
                         skk_1040, skk_1041, skk_1042, slk_1038, slk_1039, slk_1040, slk_1041, \
                         slk_1042 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1290[k] = f_15 * skk_1038[k]
                    + f_3 * pc_x[k] * slk_1038[k];

        t_1291[k] = f_15 * skk_1039[k]
                    + f_3 * pc_x[k] * slk_1039[k];

        t_1292[k] = f_15 * skk_1040[k]
                    + f_3 * pc_x[k] * slk_1040[k];

        t_1293[k] = f_15 * skk_1041[k]
                    + f_3 * pc_x[k] * slk_1041[k];

        t_1294[k] = f_15 * skk_1042[k]
                    + f_3 * pc_x[k] * slk_1042[k];
    }

#pragma omp simd aligned(t_1295, t_1296, t_1297, t_1298, pb_x, pc_x, pc_z, skl0_1296, \
                         skl0_1298, skk_1043, skl1_1296, skl1_1298, slk_1036, \
                         slk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1295[k] = f_15 * skk_1043[k]
                    + f_3 * pc_x[k] * slk_1043[k];

        t_1296[k] = pb_x[k] * skl0_1296[k]
                    - f_14 * pc_x[k] * skl1_1296[k];

        t_1297[k] = f_3 * pc_z[k] * slk_1036[k];

        t_1298[k] = pb_x[k] * skl0_1298[k]
                    - f_14 * pc_x[k] * skl1_1298[k];
    }

#pragma omp simd aligned(t_1299, t_1300, t_1301, t_1302, pb_x, pc_x, skl0_1299, skl0_1300, \
                         skl0_1301, skl0_1302, skl1_1299, skl1_1300, skl1_1301, \
                         skl1_1302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1299[k] = pb_x[k] * skl0_1299[k]
                    - f_14 * pc_x[k] * skl1_1299[k];

        t_1300[k] = pb_x[k] * skl0_1300[k]
                    - f_14 * pc_x[k] * skl1_1300[k];

        t_1301[k] = pb_x[k] * skl0_1301[k]
                    - f_14 * pc_x[k] * skl1_1301[k];

        t_1302[k] = pb_x[k] * skl0_1302[k]
                    - f_14 * pc_x[k] * skl1_1302[k];
    }

#pragma omp simd aligned(t_1303, t_1304, t_1305, pb_x, pb_z, pc_x, pc_y, pc_z, skl0_945, \
                         skl0_1304, skk_791, skl1_945, skl1_1304, \
                         slk_1043 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1303[k] = f_21 * skk_791[k]
                    + f_3 * pc_y[k] * slk_1043[k];

        t_1304[k] = pb_x[k] * skl0_1304[k]
                    - f_14 * pc_x[k] * skl1_1304[k];

        t_1305[k] = pb_z[k] * skl0_945[k]
                    - f_14 * pc_z[k] * skl1_945[k];
    }

#pragma omp simd aligned(t_1306, t_1307, t_1308, t_1309, pb_z, pc_y, pc_z, skl0_948, skk_756, \
                         skk_792, skk_794, skl1_948, slk_1044, \
                         slk_1046 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1306[k] = f_20 * skk_792[k]
                    + f_3 * pc_y[k] * slk_1044[k];

        t_1307[k] = f_15 * skk_756[k]
                    + f_3 * pc_z[k] * slk_1044[k];

        t_1308[k] = pb_z[k] * skl0_948[k]
                    - f_14 * pc_z[k] * skl1_948[k];

        t_1309[k] = f_20 * skk_794[k]
                    + f_3 * pc_y[k] * slk_1046[k];
    }

#pragma omp simd aligned(t_1310, t_1311, t_1312, pb_x, pb_z, pc_x, pc_z, skl0_951, skl0_1310, \
                         skk_759, skk_1049, skl1_951, skl1_1310, \
                         slk_1047 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1310[k] = pb_x[k] * skl0_1310[k]
                    + f_20 * skk_1049[k]
                    - f_14 * pc_x[k] * skl1_1310[k];

        t_1311[k] = pb_z[k] * skl0_951[k]
                    - f_14 * pc_z[k] * skl1_951[k];

        t_1312[k] = f_15 * skk_759[k]
                    + f_3 * pc_z[k] * slk_1047[k];
    }

#pragma omp simd aligned(t_1313, t_1314, t_1315, pb_x, pb_z, pc_x, pc_y, pc_z, skl0_955, \
                         skl0_1314, skk_797, skk_1053, skl1_955, skl1_1314, \
                         slk_1049 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1313[k] = f_20 * skk_797[k]
                    + f_3 * pc_y[k] * slk_1049[k];

        t_1314[k] = pb_x[k] * skl0_1314[k]
                    + f_19 * skk_1053[k]
                    - f_14 * pc_x[k] * skl1_1314[k];

        t_1315[k] = pb_z[k] * skl0_955[k]
                    - f_14 * pc_z[k] * skl1_955[k];
    }

#pragma omp simd aligned(t_1316, t_1317, t_1318, pb_x, pc_x, pc_y, pc_z, skl0_1317, skk_762, \
                         skk_801, skk_1056, skl1_1317, slk_1050, \
                         slk_1053 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1316[k] = f_15 * skk_762[k]
                    + f_3 * pc_z[k] * slk_1050[k];

        t_1317[k] = pb_x[k] * skl0_1317[k]
                    + f_18 * skk_1056[k]
                    - f_14 * pc_x[k] * skl1_1317[k];

        t_1318[k] = f_20 * skk_801[k]
                    + f_3 * pc_y[k] * slk_1053[k];
    }

#pragma omp simd aligned(t_1319, t_1320, t_1321, pb_x, pb_z, pc_x, pc_z, skl0_960, skl0_1319, \
                         skk_766, skk_1058, skl1_960, skl1_1319, \
                         slk_1054 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1319[k] = pb_x[k] * skl0_1319[k]
                    + f_18 * skk_1058[k]
                    - f_14 * pc_x[k] * skl1_1319[k];

        t_1320[k] = pb_z[k] * skl0_960[k]
                    - f_14 * pc_z[k] * skl1_960[k];

        t_1321[k] = f_15 * skk_766[k]
                    + f_3 * pc_z[k] * slk_1054[k];
    }

#pragma omp simd aligned(t_1322, t_1323, t_1324, pb_x, pc_x, pc_y, skl0_1322, skl0_1323, \
                         skk_806, skk_1061, skk_1062, skl1_1322, skl1_1323, \
                         slk_1058 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1322[k] = pb_x[k] * skl0_1322[k]
                    + f_17 * skk_1061[k]
                    - f_14 * pc_x[k] * skl1_1322[k];

        t_1323[k] = pb_x[k] * skl0_1323[k]
                    + f_17 * skk_1062[k]
                    - f_14 * pc_x[k] * skl1_1323[k];

        t_1324[k] = f_20 * skk_806[k]
                    + f_3 * pc_y[k] * slk_1058[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece12(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t slk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_966 = buffer.data(skl0 + 966);
    const auto *skl0_1325 = buffer.data(skl0 + 1325);
    const auto *skl0_1328 = buffer.data(skl0 + 1328);
    const auto *skl0_1329 = buffer.data(skl0 + 1329);
    const auto *skl0_1330 = buffer.data(skl0 + 1330);
    const auto *skl0_1332 = buffer.data(skl0 + 1332);
    const auto *skl0_1341 = buffer.data(skl0 + 1341);
    const auto *skl0_1343 = buffer.data(skl0 + 1343);
    const auto *skl0_1344 = buffer.data(skl0 + 1344);
    const auto *skl0_1345 = buffer.data(skl0 + 1345);
    const auto *skl0_1346 = buffer.data(skl0 + 1346);
    const auto *skl0_1347 = buffer.data(skl0 + 1347);
    const auto *skl0_1349 = buffer.data(skl0 + 1349);
    const auto *skl0_1350 = buffer.data(skl0 + 1350);
    const auto *skl0_1353 = buffer.data(skl0 + 1353);
    const auto *skl0_1355 = buffer.data(skl0 + 1355);
    const auto *skl0_1356 = buffer.data(skl0 + 1356);
    const auto *skl0_1359 = buffer.data(skl0 + 1359);
    const auto *skl0_1360 = buffer.data(skl0 + 1360);
    const auto *skl0_1362 = buffer.data(skl0 + 1362);
    const auto *skl0_1364 = buffer.data(skl0 + 1364);
    const auto *skl0_1365 = buffer.data(skl0 + 1365);
    const auto *skl0_1367 = buffer.data(skl0 + 1367);
    const auto *skl0_1368 = buffer.data(skl0 + 1368);
    const auto *skl0_1370 = buffer.data(skl0 + 1370);
    const auto *skl0_1371 = buffer.data(skl0 + 1371);
    const auto *skl0_1373 = buffer.data(skl0 + 1373);
    const auto *skl0_1374 = buffer.data(skl0 + 1374);
    const auto *skl0_1375 = buffer.data(skl0 + 1375);
    const auto *skl0_1377 = buffer.data(skl0 + 1377);
    const auto *skl0_1386 = buffer.data(skl0 + 1386);
    const auto *skl0_1388 = buffer.data(skl0 + 1388);
    const auto *skl0_1389 = buffer.data(skl0 + 1389);
    const auto *skl0_1390 = buffer.data(skl0 + 1390);
    const auto *skl0_1391 = buffer.data(skl0 + 1391);
    const auto *skl0_1392 = buffer.data(skl0 + 1392);
    const auto *skl0_1394 = buffer.data(skl0 + 1394);
    const auto *skl0_1395 = buffer.data(skl0 + 1395);
    const auto *skl0_1398 = buffer.data(skl0 + 1398);
    const auto *skl0_1400 = buffer.data(skl0 + 1400);
    const auto *skl0_1401 = buffer.data(skl0 + 1401);
    const auto *skl0_1404 = buffer.data(skl0 + 1404);
    const auto *skl0_1405 = buffer.data(skl0 + 1405);
    const auto *skl0_1407 = buffer.data(skl0 + 1407);
    const auto *skl0_1409 = buffer.data(skl0 + 1409);
    const auto *skl0_1410 = buffer.data(skl0 + 1410);
    const auto *skl0_1412 = buffer.data(skl0 + 1412);
    const auto *skl0_1413 = buffer.data(skl0 + 1413);
    const auto *skl0_1415 = buffer.data(skl0 + 1415);
    const auto *skl0_1416 = buffer.data(skl0 + 1416);
    const auto *skl0_1418 = buffer.data(skl0 + 1418);
    const auto *skl0_1419 = buffer.data(skl0 + 1419);
    const auto *skl0_1420 = buffer.data(skl0 + 1420);
    const auto *skl0_1422 = buffer.data(skl0 + 1422);
    const auto *skl0_1431 = buffer.data(skl0 + 1431);
    const auto *skl0_1433 = buffer.data(skl0 + 1433);
    const auto *skl0_1434 = buffer.data(skl0 + 1434);
    const auto *skl0_1435 = buffer.data(skl0 + 1435);
    const auto *skl0_1436 = buffer.data(skl0 + 1436);
    const auto *skl0_1437 = buffer.data(skl0 + 1437);
    const auto *skl0_1439 = buffer.data(skl0 + 1439);
    const auto *skl0_1440 = buffer.data(skl0 + 1440);
    const auto *skl0_1443 = buffer.data(skl0 + 1443);

    const auto *skk_771 = buffer.data(skk + 771);
    const auto *skk_784 = buffer.data(skk + 784);
    const auto *skk_792 = buffer.data(skk + 792);
    const auto *skk_795 = buffer.data(skk + 795);
    const auto *skk_798 = buffer.data(skk + 798);
    const auto *skk_802 = buffer.data(skk + 802);
    const auto *skk_807 = buffer.data(skk + 807);
    const auto *skk_812 = buffer.data(skk + 812);
    const auto *skk_820 = buffer.data(skk + 820);
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
    const auto *skk_863 = buffer.data(skk + 863);
    const auto *skk_864 = buffer.data(skk + 864);
    const auto *skk_866 = buffer.data(skk + 866);
    const auto *skk_869 = buffer.data(skk + 869);
    const auto *skk_873 = buffer.data(skk + 873);
    const auto *skk_878 = buffer.data(skk + 878);
    const auto *skk_884 = buffer.data(skk + 884);
    const auto *skk_899 = buffer.data(skk + 899);
    const auto *skk_900 = buffer.data(skk + 900);
    const auto *skk_902 = buffer.data(skk + 902);
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
    const auto *skk_1083 = buffer.data(skk + 1083);
    const auto *skk_1085 = buffer.data(skk + 1085);
    const auto *skk_1086 = buffer.data(skk + 1086);
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
    const auto *skk_1155 = buffer.data(skk + 1155);

    const auto *skl1_966 = buffer.data(skl1 + 966);
    const auto *skl1_1325 = buffer.data(skl1 + 1325);
    const auto *skl1_1328 = buffer.data(skl1 + 1328);
    const auto *skl1_1329 = buffer.data(skl1 + 1329);
    const auto *skl1_1330 = buffer.data(skl1 + 1330);
    const auto *skl1_1332 = buffer.data(skl1 + 1332);
    const auto *skl1_1341 = buffer.data(skl1 + 1341);
    const auto *skl1_1343 = buffer.data(skl1 + 1343);
    const auto *skl1_1344 = buffer.data(skl1 + 1344);
    const auto *skl1_1345 = buffer.data(skl1 + 1345);
    const auto *skl1_1346 = buffer.data(skl1 + 1346);
    const auto *skl1_1347 = buffer.data(skl1 + 1347);
    const auto *skl1_1349 = buffer.data(skl1 + 1349);
    const auto *skl1_1350 = buffer.data(skl1 + 1350);
    const auto *skl1_1353 = buffer.data(skl1 + 1353);
    const auto *skl1_1355 = buffer.data(skl1 + 1355);
    const auto *skl1_1356 = buffer.data(skl1 + 1356);
    const auto *skl1_1359 = buffer.data(skl1 + 1359);
    const auto *skl1_1360 = buffer.data(skl1 + 1360);
    const auto *skl1_1362 = buffer.data(skl1 + 1362);
    const auto *skl1_1364 = buffer.data(skl1 + 1364);
    const auto *skl1_1365 = buffer.data(skl1 + 1365);
    const auto *skl1_1367 = buffer.data(skl1 + 1367);
    const auto *skl1_1368 = buffer.data(skl1 + 1368);
    const auto *skl1_1370 = buffer.data(skl1 + 1370);
    const auto *skl1_1371 = buffer.data(skl1 + 1371);
    const auto *skl1_1373 = buffer.data(skl1 + 1373);
    const auto *skl1_1374 = buffer.data(skl1 + 1374);
    const auto *skl1_1375 = buffer.data(skl1 + 1375);
    const auto *skl1_1377 = buffer.data(skl1 + 1377);
    const auto *skl1_1386 = buffer.data(skl1 + 1386);
    const auto *skl1_1388 = buffer.data(skl1 + 1388);
    const auto *skl1_1389 = buffer.data(skl1 + 1389);
    const auto *skl1_1390 = buffer.data(skl1 + 1390);
    const auto *skl1_1391 = buffer.data(skl1 + 1391);
    const auto *skl1_1392 = buffer.data(skl1 + 1392);
    const auto *skl1_1394 = buffer.data(skl1 + 1394);
    const auto *skl1_1395 = buffer.data(skl1 + 1395);
    const auto *skl1_1398 = buffer.data(skl1 + 1398);
    const auto *skl1_1400 = buffer.data(skl1 + 1400);
    const auto *skl1_1401 = buffer.data(skl1 + 1401);
    const auto *skl1_1404 = buffer.data(skl1 + 1404);
    const auto *skl1_1405 = buffer.data(skl1 + 1405);
    const auto *skl1_1407 = buffer.data(skl1 + 1407);
    const auto *skl1_1409 = buffer.data(skl1 + 1409);
    const auto *skl1_1410 = buffer.data(skl1 + 1410);
    const auto *skl1_1412 = buffer.data(skl1 + 1412);
    const auto *skl1_1413 = buffer.data(skl1 + 1413);
    const auto *skl1_1415 = buffer.data(skl1 + 1415);
    const auto *skl1_1416 = buffer.data(skl1 + 1416);
    const auto *skl1_1418 = buffer.data(skl1 + 1418);
    const auto *skl1_1419 = buffer.data(skl1 + 1419);
    const auto *skl1_1420 = buffer.data(skl1 + 1420);
    const auto *skl1_1422 = buffer.data(skl1 + 1422);
    const auto *skl1_1431 = buffer.data(skl1 + 1431);
    const auto *skl1_1433 = buffer.data(skl1 + 1433);
    const auto *skl1_1434 = buffer.data(skl1 + 1434);
    const auto *skl1_1435 = buffer.data(skl1 + 1435);
    const auto *skl1_1436 = buffer.data(skl1 + 1436);
    const auto *skl1_1437 = buffer.data(skl1 + 1437);
    const auto *skl1_1439 = buffer.data(skl1 + 1439);
    const auto *skl1_1440 = buffer.data(skl1 + 1440);
    const auto *skl1_1443 = buffer.data(skl1 + 1443);

    const auto *slk_1059 = buffer.data(slk + 1059);
    const auto *slk_1064 = buffer.data(slk + 1064);
    const auto *slk_1072 = buffer.data(slk + 1072);
    const auto *slk_1073 = buffer.data(slk + 1073);
    const auto *slk_1074 = buffer.data(slk + 1074);
    const auto *slk_1075 = buffer.data(slk + 1075);
    const auto *slk_1076 = buffer.data(slk + 1076);
    const auto *slk_1077 = buffer.data(slk + 1077);
    const auto *slk_1078 = buffer.data(slk + 1078);
    const auto *slk_1079 = buffer.data(slk + 1079);
    const auto *slk_1080 = buffer.data(slk + 1080);
    const auto *slk_1082 = buffer.data(slk + 1082);
    const auto *slk_1083 = buffer.data(slk + 1083);
    const auto *slk_1085 = buffer.data(slk + 1085);
    const auto *slk_1086 = buffer.data(slk + 1086);
    const auto *slk_1089 = buffer.data(slk + 1089);
    const auto *slk_1090 = buffer.data(slk + 1090);
    const auto *slk_1094 = buffer.data(slk + 1094);
    const auto *slk_1095 = buffer.data(slk + 1095);
    const auto *slk_1100 = buffer.data(slk + 1100);
    const auto *slk_1108 = buffer.data(slk + 1108);
    const auto *slk_1109 = buffer.data(slk + 1109);
    const auto *slk_1110 = buffer.data(slk + 1110);
    const auto *slk_1111 = buffer.data(slk + 1111);
    const auto *slk_1112 = buffer.data(slk + 1112);
    const auto *slk_1113 = buffer.data(slk + 1113);
    const auto *slk_1114 = buffer.data(slk + 1114);
    const auto *slk_1115 = buffer.data(slk + 1115);
    const auto *slk_1116 = buffer.data(slk + 1116);
    const auto *slk_1118 = buffer.data(slk + 1118);
    const auto *slk_1119 = buffer.data(slk + 1119);
    const auto *slk_1121 = buffer.data(slk + 1121);
    const auto *slk_1122 = buffer.data(slk + 1122);
    const auto *slk_1125 = buffer.data(slk + 1125);
    const auto *slk_1126 = buffer.data(slk + 1126);
    const auto *slk_1130 = buffer.data(slk + 1130);
    const auto *slk_1131 = buffer.data(slk + 1131);
    const auto *slk_1136 = buffer.data(slk + 1136);
    const auto *slk_1144 = buffer.data(slk + 1144);
    const auto *slk_1145 = buffer.data(slk + 1145);
    const auto *slk_1146 = buffer.data(slk + 1146);
    const auto *slk_1147 = buffer.data(slk + 1147);
    const auto *slk_1148 = buffer.data(slk + 1148);
    const auto *slk_1149 = buffer.data(slk + 1149);
    const auto *slk_1150 = buffer.data(slk + 1150);
    const auto *slk_1151 = buffer.data(slk + 1151);
    const auto *slk_1152 = buffer.data(slk + 1152);
    const auto *slk_1154 = buffer.data(slk + 1154);

#pragma omp simd aligned(t_1325, t_1326, t_1327, pb_x, pb_z, pc_x, pc_z, skl0_966, skl0_1325, \
                         skk_771, skk_1064, skl1_966, skl1_1325, \
                         slk_1059 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1325[k] = pb_x[k] * skl0_1325[k]
                    + f_17 * skk_1064[k]
                    - f_14 * pc_x[k] * skl1_1325[k];

        t_1326[k] = pb_z[k] * skl0_966[k]
                    - f_14 * pc_z[k] * skl1_966[k];

        t_1327[k] = f_15 * skk_771[k]
                    + f_3 * pc_z[k] * slk_1059[k];
    }

#pragma omp simd aligned(t_1328, t_1329, t_1330, pb_x, pc_x, skl0_1328, skl0_1329, skl0_1330, \
                         skk_1067, skk_1068, skk_1069, skl1_1328, skl1_1329, \
                         skl1_1330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1328[k] = pb_x[k] * skl0_1328[k]
                    + f_16 * skk_1067[k]
                    - f_14 * pc_x[k] * skl1_1328[k];

        t_1329[k] = pb_x[k] * skl0_1329[k]
                    + f_16 * skk_1068[k]
                    - f_14 * pc_x[k] * skl1_1329[k];

        t_1330[k] = pb_x[k] * skl0_1330[k]
                    + f_16 * skk_1069[k]
                    - f_14 * pc_x[k] * skl1_1330[k];
    }

#pragma omp simd aligned(t_1331, t_1332, t_1333, t_1334, pb_x, pc_x, pc_y, skl0_1332, skk_812, \
                         skk_1071, skk_1072, skk_1073, skl1_1332, slk_1064, slk_1072, \
                         slk_1073 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1331[k] = f_20 * skk_812[k]
                    + f_3 * pc_y[k] * slk_1064[k];

        t_1332[k] = pb_x[k] * skl0_1332[k]
                    + f_16 * skk_1071[k]
                    - f_14 * pc_x[k] * skl1_1332[k];

        t_1333[k] = f_15 * skk_1072[k]
                    + f_3 * pc_x[k] * slk_1072[k];

        t_1334[k] = f_15 * skk_1073[k]
                    + f_3 * pc_x[k] * slk_1073[k];
    }

#pragma omp simd aligned(t_1335, t_1336, t_1337, t_1338, t_1339, pc_x, skk_1074, skk_1075, \
                         skk_1076, skk_1077, skk_1078, slk_1074, slk_1075, slk_1076, slk_1077, \
                         slk_1078 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1335[k] = f_15 * skk_1074[k]
                    + f_3 * pc_x[k] * slk_1074[k];

        t_1336[k] = f_15 * skk_1075[k]
                    + f_3 * pc_x[k] * slk_1075[k];

        t_1337[k] = f_15 * skk_1076[k]
                    + f_3 * pc_x[k] * slk_1076[k];

        t_1338[k] = f_15 * skk_1077[k]
                    + f_3 * pc_x[k] * slk_1077[k];

        t_1339[k] = f_15 * skk_1078[k]
                    + f_3 * pc_x[k] * slk_1078[k];
    }

#pragma omp simd aligned(t_1340, t_1341, t_1342, t_1343, pb_x, pc_x, pc_z, skl0_1341, \
                         skl0_1343, skk_784, skk_1079, skl1_1341, skl1_1343, slk_1072, \
                         slk_1079 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1340[k] = f_15 * skk_1079[k]
                    + f_3 * pc_x[k] * slk_1079[k];

        t_1341[k] = pb_x[k] * skl0_1341[k]
                    - f_14 * pc_x[k] * skl1_1341[k];

        t_1342[k] = f_15 * skk_784[k]
                    + f_3 * pc_z[k] * slk_1072[k];

        t_1343[k] = pb_x[k] * skl0_1343[k]
                    - f_14 * pc_x[k] * skl1_1343[k];
    }

#pragma omp simd aligned(t_1344, t_1345, t_1346, t_1347, pb_x, pc_x, skl0_1344, skl0_1345, \
                         skl0_1346, skl0_1347, skl1_1344, skl1_1345, skl1_1346, \
                         skl1_1347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1344[k] = pb_x[k] * skl0_1344[k]
                    - f_14 * pc_x[k] * skl1_1344[k];

        t_1345[k] = pb_x[k] * skl0_1345[k]
                    - f_14 * pc_x[k] * skl1_1345[k];

        t_1346[k] = pb_x[k] * skl0_1346[k]
                    - f_14 * pc_x[k] * skl1_1346[k];

        t_1347[k] = pb_x[k] * skl0_1347[k]
                    - f_14 * pc_x[k] * skl1_1347[k];
    }

#pragma omp simd aligned(t_1348, t_1349, t_1350, t_1351, pb_x, pc_x, pc_y, skl0_1349, \
                         skl0_1350, skk_827, skk_828, skk_1080, skl1_1349, skl1_1350, \
                         slk_1079, slk_1080 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1348[k] = f_20 * skk_827[k]
                    + f_3 * pc_y[k] * slk_1079[k];

        t_1349[k] = pb_x[k] * skl0_1349[k]
                    - f_14 * pc_x[k] * skl1_1349[k];

        t_1350[k] = pb_x[k] * skl0_1350[k]
                    + f_0 * skk_1080[k]
                    - f_14 * pc_x[k] * skl1_1350[k];

        t_1351[k] = f_19 * skk_828[k]
                    + f_3 * pc_y[k] * slk_1080[k];
    }

#pragma omp simd aligned(t_1352, t_1353, t_1354, pb_x, pc_x, pc_y, pc_z, skl0_1353, skk_792, \
                         skk_830, skk_1083, skl1_1353, slk_1080, \
                         slk_1082 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1352[k] = f_16 * skk_792[k]
                    + f_3 * pc_z[k] * slk_1080[k];

        t_1353[k] = pb_x[k] * skl0_1353[k]
                    + f_20 * skk_1083[k]
                    - f_14 * pc_x[k] * skl1_1353[k];

        t_1354[k] = f_19 * skk_830[k]
                    + f_3 * pc_y[k] * slk_1082[k];
    }

#pragma omp simd aligned(t_1355, t_1356, t_1357, pb_x, pc_x, pc_z, skl0_1355, skl0_1356, \
                         skk_795, skk_1085, skk_1086, skl1_1355, skl1_1356, \
                         slk_1083 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1355[k] = pb_x[k] * skl0_1355[k]
                    + f_20 * skk_1085[k]
                    - f_14 * pc_x[k] * skl1_1355[k];

        t_1356[k] = pb_x[k] * skl0_1356[k]
                    + f_19 * skk_1086[k]
                    - f_14 * pc_x[k] * skl1_1356[k];

        t_1357[k] = f_16 * skk_795[k]
                    + f_3 * pc_z[k] * slk_1083[k];
    }

#pragma omp simd aligned(t_1358, t_1359, t_1360, pb_x, pc_x, pc_y, skl0_1359, skl0_1360, \
                         skk_833, skk_1089, skk_1090, skl1_1359, skl1_1360, \
                         slk_1085 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1358[k] = f_19 * skk_833[k]
                    + f_3 * pc_y[k] * slk_1085[k];

        t_1359[k] = pb_x[k] * skl0_1359[k]
                    + f_19 * skk_1089[k]
                    - f_14 * pc_x[k] * skl1_1359[k];

        t_1360[k] = pb_x[k] * skl0_1360[k]
                    + f_18 * skk_1090[k]
                    - f_14 * pc_x[k] * skl1_1360[k];
    }

#pragma omp simd aligned(t_1361, t_1362, t_1363, pb_x, pc_x, pc_y, pc_z, skl0_1362, skk_798, \
                         skk_837, skk_1092, skl1_1362, slk_1086, \
                         slk_1089 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1361[k] = f_16 * skk_798[k]
                    + f_3 * pc_z[k] * slk_1086[k];

        t_1362[k] = pb_x[k] * skl0_1362[k]
                    + f_18 * skk_1092[k]
                    - f_14 * pc_x[k] * skl1_1362[k];

        t_1363[k] = f_19 * skk_837[k]
                    + f_3 * pc_y[k] * slk_1089[k];
    }

#pragma omp simd aligned(t_1364, t_1365, t_1366, pb_x, pc_x, pc_z, skl0_1364, skl0_1365, \
                         skk_802, skk_1094, skk_1095, skl1_1364, skl1_1365, \
                         slk_1090 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1364[k] = pb_x[k] * skl0_1364[k]
                    + f_18 * skk_1094[k]
                    - f_14 * pc_x[k] * skl1_1364[k];

        t_1365[k] = pb_x[k] * skl0_1365[k]
                    + f_17 * skk_1095[k]
                    - f_14 * pc_x[k] * skl1_1365[k];

        t_1366[k] = f_16 * skk_802[k]
                    + f_3 * pc_z[k] * slk_1090[k];
    }

#pragma omp simd aligned(t_1367, t_1368, t_1369, pb_x, pc_x, pc_y, skl0_1367, skl0_1368, \
                         skk_842, skk_1097, skk_1098, skl1_1367, skl1_1368, \
                         slk_1094 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1367[k] = pb_x[k] * skl0_1367[k]
                    + f_17 * skk_1097[k]
                    - f_14 * pc_x[k] * skl1_1367[k];

        t_1368[k] = pb_x[k] * skl0_1368[k]
                    + f_17 * skk_1098[k]
                    - f_14 * pc_x[k] * skl1_1368[k];

        t_1369[k] = f_19 * skk_842[k]
                    + f_3 * pc_y[k] * slk_1094[k];
    }

#pragma omp simd aligned(t_1370, t_1371, t_1372, pb_x, pc_x, pc_z, skl0_1370, skl0_1371, \
                         skk_807, skk_1100, skk_1101, skl1_1370, skl1_1371, \
                         slk_1095 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1370[k] = pb_x[k] * skl0_1370[k]
                    + f_17 * skk_1100[k]
                    - f_14 * pc_x[k] * skl1_1370[k];

        t_1371[k] = pb_x[k] * skl0_1371[k]
                    + f_16 * skk_1101[k]
                    - f_14 * pc_x[k] * skl1_1371[k];

        t_1372[k] = f_16 * skk_807[k]
                    + f_3 * pc_z[k] * slk_1095[k];
    }

#pragma omp simd aligned(t_1373, t_1374, t_1375, pb_x, pc_x, skl0_1373, skl0_1374, skl0_1375, \
                         skk_1103, skk_1104, skk_1105, skl1_1373, skl1_1374, \
                         skl1_1375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1373[k] = pb_x[k] * skl0_1373[k]
                    + f_16 * skk_1103[k]
                    - f_14 * pc_x[k] * skl1_1373[k];

        t_1374[k] = pb_x[k] * skl0_1374[k]
                    + f_16 * skk_1104[k]
                    - f_14 * pc_x[k] * skl1_1374[k];

        t_1375[k] = pb_x[k] * skl0_1375[k]
                    + f_16 * skk_1105[k]
                    - f_14 * pc_x[k] * skl1_1375[k];
    }

#pragma omp simd aligned(t_1376, t_1377, t_1378, t_1379, pb_x, pc_x, pc_y, skl0_1377, skk_848, \
                         skk_1107, skk_1108, skk_1109, skl1_1377, slk_1100, slk_1108, \
                         slk_1109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1376[k] = f_19 * skk_848[k]
                    + f_3 * pc_y[k] * slk_1100[k];

        t_1377[k] = pb_x[k] * skl0_1377[k]
                    + f_16 * skk_1107[k]
                    - f_14 * pc_x[k] * skl1_1377[k];

        t_1378[k] = f_15 * skk_1108[k]
                    + f_3 * pc_x[k] * slk_1108[k];

        t_1379[k] = f_15 * skk_1109[k]
                    + f_3 * pc_x[k] * slk_1109[k];
    }

#pragma omp simd aligned(t_1380, t_1381, t_1382, t_1383, t_1384, pc_x, skk_1110, skk_1111, \
                         skk_1112, skk_1113, skk_1114, slk_1110, slk_1111, slk_1112, slk_1113, \
                         slk_1114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1380[k] = f_15 * skk_1110[k]
                    + f_3 * pc_x[k] * slk_1110[k];

        t_1381[k] = f_15 * skk_1111[k]
                    + f_3 * pc_x[k] * slk_1111[k];

        t_1382[k] = f_15 * skk_1112[k]
                    + f_3 * pc_x[k] * slk_1112[k];

        t_1383[k] = f_15 * skk_1113[k]
                    + f_3 * pc_x[k] * slk_1113[k];

        t_1384[k] = f_15 * skk_1114[k]
                    + f_3 * pc_x[k] * slk_1114[k];
    }

#pragma omp simd aligned(t_1385, t_1386, t_1387, t_1388, pb_x, pc_x, pc_z, skl0_1386, \
                         skl0_1388, skk_820, skk_1115, skl1_1386, skl1_1388, slk_1108, \
                         slk_1115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1385[k] = f_15 * skk_1115[k]
                    + f_3 * pc_x[k] * slk_1115[k];

        t_1386[k] = pb_x[k] * skl0_1386[k]
                    - f_14 * pc_x[k] * skl1_1386[k];

        t_1387[k] = f_16 * skk_820[k]
                    + f_3 * pc_z[k] * slk_1108[k];

        t_1388[k] = pb_x[k] * skl0_1388[k]
                    - f_14 * pc_x[k] * skl1_1388[k];
    }

#pragma omp simd aligned(t_1389, t_1390, t_1391, t_1392, pb_x, pc_x, skl0_1389, skl0_1390, \
                         skl0_1391, skl0_1392, skl1_1389, skl1_1390, skl1_1391, \
                         skl1_1392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1389[k] = pb_x[k] * skl0_1389[k]
                    - f_14 * pc_x[k] * skl1_1389[k];

        t_1390[k] = pb_x[k] * skl0_1390[k]
                    - f_14 * pc_x[k] * skl1_1390[k];

        t_1391[k] = pb_x[k] * skl0_1391[k]
                    - f_14 * pc_x[k] * skl1_1391[k];

        t_1392[k] = pb_x[k] * skl0_1392[k]
                    - f_14 * pc_x[k] * skl1_1392[k];
    }

#pragma omp simd aligned(t_1393, t_1394, t_1395, t_1396, pb_x, pc_x, pc_y, skl0_1394, \
                         skl0_1395, skk_863, skk_864, skk_1116, skl1_1394, skl1_1395, \
                         slk_1115, slk_1116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1393[k] = f_19 * skk_863[k]
                    + f_3 * pc_y[k] * slk_1115[k];

        t_1394[k] = pb_x[k] * skl0_1394[k]
                    - f_14 * pc_x[k] * skl1_1394[k];

        t_1395[k] = pb_x[k] * skl0_1395[k]
                    + f_0 * skk_1116[k]
                    - f_14 * pc_x[k] * skl1_1395[k];

        t_1396[k] = f_18 * skk_864[k]
                    + f_3 * pc_y[k] * slk_1116[k];
    }

#pragma omp simd aligned(t_1397, t_1398, t_1399, pb_x, pc_x, pc_y, pc_z, skl0_1398, skk_828, \
                         skk_866, skk_1119, skl1_1398, slk_1116, \
                         slk_1118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1397[k] = f_17 * skk_828[k]
                    + f_3 * pc_z[k] * slk_1116[k];

        t_1398[k] = pb_x[k] * skl0_1398[k]
                    + f_20 * skk_1119[k]
                    - f_14 * pc_x[k] * skl1_1398[k];

        t_1399[k] = f_18 * skk_866[k]
                    + f_3 * pc_y[k] * slk_1118[k];
    }

#pragma omp simd aligned(t_1400, t_1401, t_1402, pb_x, pc_x, pc_z, skl0_1400, skl0_1401, \
                         skk_831, skk_1121, skk_1122, skl1_1400, skl1_1401, \
                         slk_1119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1400[k] = pb_x[k] * skl0_1400[k]
                    + f_20 * skk_1121[k]
                    - f_14 * pc_x[k] * skl1_1400[k];

        t_1401[k] = pb_x[k] * skl0_1401[k]
                    + f_19 * skk_1122[k]
                    - f_14 * pc_x[k] * skl1_1401[k];

        t_1402[k] = f_17 * skk_831[k]
                    + f_3 * pc_z[k] * slk_1119[k];
    }

#pragma omp simd aligned(t_1403, t_1404, t_1405, pb_x, pc_x, pc_y, skl0_1404, skl0_1405, \
                         skk_869, skk_1125, skk_1126, skl1_1404, skl1_1405, \
                         slk_1121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1403[k] = f_18 * skk_869[k]
                    + f_3 * pc_y[k] * slk_1121[k];

        t_1404[k] = pb_x[k] * skl0_1404[k]
                    + f_19 * skk_1125[k]
                    - f_14 * pc_x[k] * skl1_1404[k];

        t_1405[k] = pb_x[k] * skl0_1405[k]
                    + f_18 * skk_1126[k]
                    - f_14 * pc_x[k] * skl1_1405[k];
    }

#pragma omp simd aligned(t_1406, t_1407, t_1408, pb_x, pc_x, pc_y, pc_z, skl0_1407, skk_834, \
                         skk_873, skk_1128, skl1_1407, slk_1122, \
                         slk_1125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1406[k] = f_17 * skk_834[k]
                    + f_3 * pc_z[k] * slk_1122[k];

        t_1407[k] = pb_x[k] * skl0_1407[k]
                    + f_18 * skk_1128[k]
                    - f_14 * pc_x[k] * skl1_1407[k];

        t_1408[k] = f_18 * skk_873[k]
                    + f_3 * pc_y[k] * slk_1125[k];
    }

#pragma omp simd aligned(t_1409, t_1410, t_1411, pb_x, pc_x, pc_z, skl0_1409, skl0_1410, \
                         skk_838, skk_1130, skk_1131, skl1_1409, skl1_1410, \
                         slk_1126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1409[k] = pb_x[k] * skl0_1409[k]
                    + f_18 * skk_1130[k]
                    - f_14 * pc_x[k] * skl1_1409[k];

        t_1410[k] = pb_x[k] * skl0_1410[k]
                    + f_17 * skk_1131[k]
                    - f_14 * pc_x[k] * skl1_1410[k];

        t_1411[k] = f_17 * skk_838[k]
                    + f_3 * pc_z[k] * slk_1126[k];
    }

#pragma omp simd aligned(t_1412, t_1413, t_1414, pb_x, pc_x, pc_y, skl0_1412, skl0_1413, \
                         skk_878, skk_1133, skk_1134, skl1_1412, skl1_1413, \
                         slk_1130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1412[k] = pb_x[k] * skl0_1412[k]
                    + f_17 * skk_1133[k]
                    - f_14 * pc_x[k] * skl1_1412[k];

        t_1413[k] = pb_x[k] * skl0_1413[k]
                    + f_17 * skk_1134[k]
                    - f_14 * pc_x[k] * skl1_1413[k];

        t_1414[k] = f_18 * skk_878[k]
                    + f_3 * pc_y[k] * slk_1130[k];
    }

#pragma omp simd aligned(t_1415, t_1416, t_1417, pb_x, pc_x, pc_z, skl0_1415, skl0_1416, \
                         skk_843, skk_1136, skk_1137, skl1_1415, skl1_1416, \
                         slk_1131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1415[k] = pb_x[k] * skl0_1415[k]
                    + f_17 * skk_1136[k]
                    - f_14 * pc_x[k] * skl1_1415[k];

        t_1416[k] = pb_x[k] * skl0_1416[k]
                    + f_16 * skk_1137[k]
                    - f_14 * pc_x[k] * skl1_1416[k];

        t_1417[k] = f_17 * skk_843[k]
                    + f_3 * pc_z[k] * slk_1131[k];
    }

#pragma omp simd aligned(t_1418, t_1419, t_1420, pb_x, pc_x, skl0_1418, skl0_1419, skl0_1420, \
                         skk_1139, skk_1140, skk_1141, skl1_1418, skl1_1419, \
                         skl1_1420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1418[k] = pb_x[k] * skl0_1418[k]
                    + f_16 * skk_1139[k]
                    - f_14 * pc_x[k] * skl1_1418[k];

        t_1419[k] = pb_x[k] * skl0_1419[k]
                    + f_16 * skk_1140[k]
                    - f_14 * pc_x[k] * skl1_1419[k];

        t_1420[k] = pb_x[k] * skl0_1420[k]
                    + f_16 * skk_1141[k]
                    - f_14 * pc_x[k] * skl1_1420[k];
    }

#pragma omp simd aligned(t_1421, t_1422, t_1423, t_1424, pb_x, pc_x, pc_y, skl0_1422, skk_884, \
                         skk_1143, skk_1144, skk_1145, skl1_1422, slk_1136, slk_1144, \
                         slk_1145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1421[k] = f_18 * skk_884[k]
                    + f_3 * pc_y[k] * slk_1136[k];

        t_1422[k] = pb_x[k] * skl0_1422[k]
                    + f_16 * skk_1143[k]
                    - f_14 * pc_x[k] * skl1_1422[k];

        t_1423[k] = f_15 * skk_1144[k]
                    + f_3 * pc_x[k] * slk_1144[k];

        t_1424[k] = f_15 * skk_1145[k]
                    + f_3 * pc_x[k] * slk_1145[k];
    }

#pragma omp simd aligned(t_1425, t_1426, t_1427, t_1428, t_1429, pc_x, skk_1146, skk_1147, \
                         skk_1148, skk_1149, skk_1150, slk_1146, slk_1147, slk_1148, slk_1149, \
                         slk_1150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1425[k] = f_15 * skk_1146[k]
                    + f_3 * pc_x[k] * slk_1146[k];

        t_1426[k] = f_15 * skk_1147[k]
                    + f_3 * pc_x[k] * slk_1147[k];

        t_1427[k] = f_15 * skk_1148[k]
                    + f_3 * pc_x[k] * slk_1148[k];

        t_1428[k] = f_15 * skk_1149[k]
                    + f_3 * pc_x[k] * slk_1149[k];

        t_1429[k] = f_15 * skk_1150[k]
                    + f_3 * pc_x[k] * slk_1150[k];
    }

#pragma omp simd aligned(t_1430, t_1431, t_1432, t_1433, pb_x, pc_x, pc_z, skl0_1431, \
                         skl0_1433, skk_856, skk_1151, skl1_1431, skl1_1433, slk_1144, \
                         slk_1151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1430[k] = f_15 * skk_1151[k]
                    + f_3 * pc_x[k] * slk_1151[k];

        t_1431[k] = pb_x[k] * skl0_1431[k]
                    - f_14 * pc_x[k] * skl1_1431[k];

        t_1432[k] = f_17 * skk_856[k]
                    + f_3 * pc_z[k] * slk_1144[k];

        t_1433[k] = pb_x[k] * skl0_1433[k]
                    - f_14 * pc_x[k] * skl1_1433[k];
    }

#pragma omp simd aligned(t_1434, t_1435, t_1436, t_1437, pb_x, pc_x, skl0_1434, skl0_1435, \
                         skl0_1436, skl0_1437, skl1_1434, skl1_1435, skl1_1436, \
                         skl1_1437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1434[k] = pb_x[k] * skl0_1434[k]
                    - f_14 * pc_x[k] * skl1_1434[k];

        t_1435[k] = pb_x[k] * skl0_1435[k]
                    - f_14 * pc_x[k] * skl1_1435[k];

        t_1436[k] = pb_x[k] * skl0_1436[k]
                    - f_14 * pc_x[k] * skl1_1436[k];

        t_1437[k] = pb_x[k] * skl0_1437[k]
                    - f_14 * pc_x[k] * skl1_1437[k];
    }

#pragma omp simd aligned(t_1438, t_1439, t_1440, t_1441, pb_x, pc_x, pc_y, skl0_1439, \
                         skl0_1440, skk_899, skk_900, skk_1152, skl1_1439, skl1_1440, \
                         slk_1151, slk_1152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1438[k] = f_18 * skk_899[k]
                    + f_3 * pc_y[k] * slk_1151[k];

        t_1439[k] = pb_x[k] * skl0_1439[k]
                    - f_14 * pc_x[k] * skl1_1439[k];

        t_1440[k] = pb_x[k] * skl0_1440[k]
                    + f_0 * skk_1152[k]
                    - f_14 * pc_x[k] * skl1_1440[k];

        t_1441[k] = f_17 * skk_900[k]
                    + f_3 * pc_y[k] * slk_1152[k];
    }

#pragma omp simd aligned(t_1442, t_1443, t_1444, pb_x, pc_x, pc_y, pc_z, skl0_1443, skk_864, \
                         skk_902, skk_1155, skl1_1443, slk_1152, \
                         slk_1154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1442[k] = f_18 * skk_864[k]
                    + f_3 * pc_z[k] * slk_1152[k];

        t_1443[k] = pb_x[k] * skl0_1443[k]
                    + f_20 * skk_1155[k]
                    - f_14 * pc_x[k] * skl1_1443[k];

        t_1444[k] = f_17 * skk_902[k]
                    + f_3 * pc_y[k] * slk_1154[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece13(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t slk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
    const auto f_3 = p / q;
    const auto f_14 = gamma / q;
    const auto f_15 = 0.5 / q;
    const auto f_16 = 1.0 / q;
    const auto f_17 = 1.5 / q;
    const auto f_18 = 2.0 / q;
    const auto f_19 = 2.5 / q;
    const auto f_20 = 3.0 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_1215 = buffer.data(skl0 + 1215);
    const auto *skl0_1220 = buffer.data(skl0 + 1220);
    const auto *skl0_1224 = buffer.data(skl0 + 1224);
    const auto *skl0_1229 = buffer.data(skl0 + 1229);
    const auto *skl0_1235 = buffer.data(skl0 + 1235);
    const auto *skl0_1242 = buffer.data(skl0 + 1242);
    const auto *skl0_1445 = buffer.data(skl0 + 1445);
    const auto *skl0_1446 = buffer.data(skl0 + 1446);
    const auto *skl0_1449 = buffer.data(skl0 + 1449);
    const auto *skl0_1450 = buffer.data(skl0 + 1450);
    const auto *skl0_1452 = buffer.data(skl0 + 1452);
    const auto *skl0_1454 = buffer.data(skl0 + 1454);
    const auto *skl0_1455 = buffer.data(skl0 + 1455);
    const auto *skl0_1457 = buffer.data(skl0 + 1457);
    const auto *skl0_1458 = buffer.data(skl0 + 1458);
    const auto *skl0_1460 = buffer.data(skl0 + 1460);
    const auto *skl0_1461 = buffer.data(skl0 + 1461);
    const auto *skl0_1463 = buffer.data(skl0 + 1463);
    const auto *skl0_1464 = buffer.data(skl0 + 1464);
    const auto *skl0_1465 = buffer.data(skl0 + 1465);
    const auto *skl0_1467 = buffer.data(skl0 + 1467);
    const auto *skl0_1476 = buffer.data(skl0 + 1476);
    const auto *skl0_1478 = buffer.data(skl0 + 1478);
    const auto *skl0_1479 = buffer.data(skl0 + 1479);
    const auto *skl0_1480 = buffer.data(skl0 + 1480);
    const auto *skl0_1481 = buffer.data(skl0 + 1481);
    const auto *skl0_1482 = buffer.data(skl0 + 1482);
    const auto *skl0_1484 = buffer.data(skl0 + 1484);
    const auto *skl0_1485 = buffer.data(skl0 + 1485);
    const auto *skl0_1488 = buffer.data(skl0 + 1488);
    const auto *skl0_1490 = buffer.data(skl0 + 1490);
    const auto *skl0_1491 = buffer.data(skl0 + 1491);
    const auto *skl0_1494 = buffer.data(skl0 + 1494);
    const auto *skl0_1495 = buffer.data(skl0 + 1495);
    const auto *skl0_1497 = buffer.data(skl0 + 1497);
    const auto *skl0_1499 = buffer.data(skl0 + 1499);
    const auto *skl0_1500 = buffer.data(skl0 + 1500);
    const auto *skl0_1502 = buffer.data(skl0 + 1502);
    const auto *skl0_1503 = buffer.data(skl0 + 1503);
    const auto *skl0_1505 = buffer.data(skl0 + 1505);
    const auto *skl0_1506 = buffer.data(skl0 + 1506);
    const auto *skl0_1508 = buffer.data(skl0 + 1508);
    const auto *skl0_1509 = buffer.data(skl0 + 1509);
    const auto *skl0_1510 = buffer.data(skl0 + 1510);
    const auto *skl0_1512 = buffer.data(skl0 + 1512);
    const auto *skl0_1521 = buffer.data(skl0 + 1521);
    const auto *skl0_1523 = buffer.data(skl0 + 1523);
    const auto *skl0_1524 = buffer.data(skl0 + 1524);
    const auto *skl0_1525 = buffer.data(skl0 + 1525);
    const auto *skl0_1526 = buffer.data(skl0 + 1526);
    const auto *skl0_1527 = buffer.data(skl0 + 1527);
    const auto *skl0_1529 = buffer.data(skl0 + 1529);
    const auto *skl0_1533 = buffer.data(skl0 + 1533);
    const auto *skl0_1536 = buffer.data(skl0 + 1536);
    const auto *skl0_1540 = buffer.data(skl0 + 1540);
    const auto *skl0_1542 = buffer.data(skl0 + 1542);
    const auto *skl0_1545 = buffer.data(skl0 + 1545);
    const auto *skl0_1547 = buffer.data(skl0 + 1547);
    const auto *skl0_1548 = buffer.data(skl0 + 1548);
    const auto *skl0_1551 = buffer.data(skl0 + 1551);
    const auto *skl0_1553 = buffer.data(skl0 + 1553);
    const auto *skl0_1554 = buffer.data(skl0 + 1554);
    const auto *skl0_1555 = buffer.data(skl0 + 1555);

    const auto *skk_867 = buffer.data(skk + 867);
    const auto *skk_870 = buffer.data(skk + 870);
    const auto *skk_874 = buffer.data(skk + 874);
    const auto *skk_879 = buffer.data(skk + 879);
    const auto *skk_892 = buffer.data(skk + 892);
    const auto *skk_900 = buffer.data(skk + 900);
    const auto *skk_903 = buffer.data(skk + 903);
    const auto *skk_905 = buffer.data(skk + 905);
    const auto *skk_906 = buffer.data(skk + 906);
    const auto *skk_909 = buffer.data(skk + 909);
    const auto *skk_910 = buffer.data(skk + 910);
    const auto *skk_914 = buffer.data(skk + 914);
    const auto *skk_915 = buffer.data(skk + 915);
    const auto *skk_920 = buffer.data(skk + 920);
    const auto *skk_928 = buffer.data(skk + 928);
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
    const auto *skk_971 = buffer.data(skk + 971);
    const auto *skk_972 = buffer.data(skk + 972);
    const auto *skk_974 = buffer.data(skk + 974);
    const auto *skk_977 = buffer.data(skk + 977);
    const auto *skk_981 = buffer.data(skk + 981);
    const auto *skk_986 = buffer.data(skk + 986);
    const auto *skk_992 = buffer.data(skk + 992);
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
    const auto *skk_1188 = buffer.data(skk + 1188);
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
    const auto *skk_1227 = buffer.data(skk + 1227);
    const auto *skk_1230 = buffer.data(skk + 1230);
    const auto *skk_1234 = buffer.data(skk + 1234);
    const auto *skk_1236 = buffer.data(skk + 1236);
    const auto *skk_1239 = buffer.data(skk + 1239);
    const auto *skk_1241 = buffer.data(skk + 1241);
    const auto *skk_1242 = buffer.data(skk + 1242);
    const auto *skk_1245 = buffer.data(skk + 1245);
    const auto *skk_1247 = buffer.data(skk + 1247);
    const auto *skk_1248 = buffer.data(skk + 1248);
    const auto *skk_1249 = buffer.data(skk + 1249);
    const auto *skk_1252 = buffer.data(skk + 1252);
    const auto *skk_1253 = buffer.data(skk + 1253);

    const auto *skl1_1215 = buffer.data(skl1 + 1215);
    const auto *skl1_1220 = buffer.data(skl1 + 1220);
    const auto *skl1_1224 = buffer.data(skl1 + 1224);
    const auto *skl1_1229 = buffer.data(skl1 + 1229);
    const auto *skl1_1235 = buffer.data(skl1 + 1235);
    const auto *skl1_1242 = buffer.data(skl1 + 1242);
    const auto *skl1_1445 = buffer.data(skl1 + 1445);
    const auto *skl1_1446 = buffer.data(skl1 + 1446);
    const auto *skl1_1449 = buffer.data(skl1 + 1449);
    const auto *skl1_1450 = buffer.data(skl1 + 1450);
    const auto *skl1_1452 = buffer.data(skl1 + 1452);
    const auto *skl1_1454 = buffer.data(skl1 + 1454);
    const auto *skl1_1455 = buffer.data(skl1 + 1455);
    const auto *skl1_1457 = buffer.data(skl1 + 1457);
    const auto *skl1_1458 = buffer.data(skl1 + 1458);
    const auto *skl1_1460 = buffer.data(skl1 + 1460);
    const auto *skl1_1461 = buffer.data(skl1 + 1461);
    const auto *skl1_1463 = buffer.data(skl1 + 1463);
    const auto *skl1_1464 = buffer.data(skl1 + 1464);
    const auto *skl1_1465 = buffer.data(skl1 + 1465);
    const auto *skl1_1467 = buffer.data(skl1 + 1467);
    const auto *skl1_1476 = buffer.data(skl1 + 1476);
    const auto *skl1_1478 = buffer.data(skl1 + 1478);
    const auto *skl1_1479 = buffer.data(skl1 + 1479);
    const auto *skl1_1480 = buffer.data(skl1 + 1480);
    const auto *skl1_1481 = buffer.data(skl1 + 1481);
    const auto *skl1_1482 = buffer.data(skl1 + 1482);
    const auto *skl1_1484 = buffer.data(skl1 + 1484);
    const auto *skl1_1485 = buffer.data(skl1 + 1485);
    const auto *skl1_1488 = buffer.data(skl1 + 1488);
    const auto *skl1_1490 = buffer.data(skl1 + 1490);
    const auto *skl1_1491 = buffer.data(skl1 + 1491);
    const auto *skl1_1494 = buffer.data(skl1 + 1494);
    const auto *skl1_1495 = buffer.data(skl1 + 1495);
    const auto *skl1_1497 = buffer.data(skl1 + 1497);
    const auto *skl1_1499 = buffer.data(skl1 + 1499);
    const auto *skl1_1500 = buffer.data(skl1 + 1500);
    const auto *skl1_1502 = buffer.data(skl1 + 1502);
    const auto *skl1_1503 = buffer.data(skl1 + 1503);
    const auto *skl1_1505 = buffer.data(skl1 + 1505);
    const auto *skl1_1506 = buffer.data(skl1 + 1506);
    const auto *skl1_1508 = buffer.data(skl1 + 1508);
    const auto *skl1_1509 = buffer.data(skl1 + 1509);
    const auto *skl1_1510 = buffer.data(skl1 + 1510);
    const auto *skl1_1512 = buffer.data(skl1 + 1512);
    const auto *skl1_1521 = buffer.data(skl1 + 1521);
    const auto *skl1_1523 = buffer.data(skl1 + 1523);
    const auto *skl1_1524 = buffer.data(skl1 + 1524);
    const auto *skl1_1525 = buffer.data(skl1 + 1525);
    const auto *skl1_1526 = buffer.data(skl1 + 1526);
    const auto *skl1_1527 = buffer.data(skl1 + 1527);
    const auto *skl1_1529 = buffer.data(skl1 + 1529);
    const auto *skl1_1533 = buffer.data(skl1 + 1533);
    const auto *skl1_1536 = buffer.data(skl1 + 1536);
    const auto *skl1_1540 = buffer.data(skl1 + 1540);
    const auto *skl1_1542 = buffer.data(skl1 + 1542);
    const auto *skl1_1545 = buffer.data(skl1 + 1545);
    const auto *skl1_1547 = buffer.data(skl1 + 1547);
    const auto *skl1_1548 = buffer.data(skl1 + 1548);
    const auto *skl1_1551 = buffer.data(skl1 + 1551);
    const auto *skl1_1553 = buffer.data(skl1 + 1553);
    const auto *skl1_1554 = buffer.data(skl1 + 1554);
    const auto *skl1_1555 = buffer.data(skl1 + 1555);

    const auto *slk_1155 = buffer.data(slk + 1155);
    const auto *slk_1157 = buffer.data(slk + 1157);
    const auto *slk_1158 = buffer.data(slk + 1158);
    const auto *slk_1161 = buffer.data(slk + 1161);
    const auto *slk_1162 = buffer.data(slk + 1162);
    const auto *slk_1166 = buffer.data(slk + 1166);
    const auto *slk_1167 = buffer.data(slk + 1167);
    const auto *slk_1172 = buffer.data(slk + 1172);
    const auto *slk_1180 = buffer.data(slk + 1180);
    const auto *slk_1181 = buffer.data(slk + 1181);
    const auto *slk_1182 = buffer.data(slk + 1182);
    const auto *slk_1183 = buffer.data(slk + 1183);
    const auto *slk_1184 = buffer.data(slk + 1184);
    const auto *slk_1185 = buffer.data(slk + 1185);
    const auto *slk_1186 = buffer.data(slk + 1186);
    const auto *slk_1187 = buffer.data(slk + 1187);
    const auto *slk_1188 = buffer.data(slk + 1188);
    const auto *slk_1190 = buffer.data(slk + 1190);
    const auto *slk_1191 = buffer.data(slk + 1191);
    const auto *slk_1193 = buffer.data(slk + 1193);
    const auto *slk_1194 = buffer.data(slk + 1194);
    const auto *slk_1197 = buffer.data(slk + 1197);
    const auto *slk_1198 = buffer.data(slk + 1198);
    const auto *slk_1202 = buffer.data(slk + 1202);
    const auto *slk_1203 = buffer.data(slk + 1203);
    const auto *slk_1208 = buffer.data(slk + 1208);
    const auto *slk_1216 = buffer.data(slk + 1216);
    const auto *slk_1217 = buffer.data(slk + 1217);
    const auto *slk_1218 = buffer.data(slk + 1218);
    const auto *slk_1219 = buffer.data(slk + 1219);
    const auto *slk_1220 = buffer.data(slk + 1220);
    const auto *slk_1221 = buffer.data(slk + 1221);
    const auto *slk_1222 = buffer.data(slk + 1222);
    const auto *slk_1223 = buffer.data(slk + 1223);
    const auto *slk_1224 = buffer.data(slk + 1224);
    const auto *slk_1226 = buffer.data(slk + 1226);
    const auto *slk_1227 = buffer.data(slk + 1227);
    const auto *slk_1229 = buffer.data(slk + 1229);
    const auto *slk_1230 = buffer.data(slk + 1230);
    const auto *slk_1233 = buffer.data(slk + 1233);
    const auto *slk_1234 = buffer.data(slk + 1234);
    const auto *slk_1238 = buffer.data(slk + 1238);
    const auto *slk_1239 = buffer.data(slk + 1239);
    const auto *slk_1244 = buffer.data(slk + 1244);
    const auto *slk_1252 = buffer.data(slk + 1252);
    const auto *slk_1253 = buffer.data(slk + 1253);

#pragma omp simd aligned(t_1445, t_1446, t_1447, pb_x, pc_x, pc_z, skl0_1445, skl0_1446, \
                         skk_867, skk_1157, skk_1158, skl1_1445, skl1_1446, \
                         slk_1155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1445[k] = pb_x[k] * skl0_1445[k]
                    + f_20 * skk_1157[k]
                    - f_14 * pc_x[k] * skl1_1445[k];

        t_1446[k] = pb_x[k] * skl0_1446[k]
                    + f_19 * skk_1158[k]
                    - f_14 * pc_x[k] * skl1_1446[k];

        t_1447[k] = f_18 * skk_867[k]
                    + f_3 * pc_z[k] * slk_1155[k];
    }

#pragma omp simd aligned(t_1448, t_1449, t_1450, pb_x, pc_x, pc_y, skl0_1449, skl0_1450, \
                         skk_905, skk_1161, skk_1162, skl1_1449, skl1_1450, \
                         slk_1157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1448[k] = f_17 * skk_905[k]
                    + f_3 * pc_y[k] * slk_1157[k];

        t_1449[k] = pb_x[k] * skl0_1449[k]
                    + f_19 * skk_1161[k]
                    - f_14 * pc_x[k] * skl1_1449[k];

        t_1450[k] = pb_x[k] * skl0_1450[k]
                    + f_18 * skk_1162[k]
                    - f_14 * pc_x[k] * skl1_1450[k];
    }

#pragma omp simd aligned(t_1451, t_1452, t_1453, pb_x, pc_x, pc_y, pc_z, skl0_1452, skk_870, \
                         skk_909, skk_1164, skl1_1452, slk_1158, \
                         slk_1161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1451[k] = f_18 * skk_870[k]
                    + f_3 * pc_z[k] * slk_1158[k];

        t_1452[k] = pb_x[k] * skl0_1452[k]
                    + f_18 * skk_1164[k]
                    - f_14 * pc_x[k] * skl1_1452[k];

        t_1453[k] = f_17 * skk_909[k]
                    + f_3 * pc_y[k] * slk_1161[k];
    }

#pragma omp simd aligned(t_1454, t_1455, t_1456, pb_x, pc_x, pc_z, skl0_1454, skl0_1455, \
                         skk_874, skk_1166, skk_1167, skl1_1454, skl1_1455, \
                         slk_1162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1454[k] = pb_x[k] * skl0_1454[k]
                    + f_18 * skk_1166[k]
                    - f_14 * pc_x[k] * skl1_1454[k];

        t_1455[k] = pb_x[k] * skl0_1455[k]
                    + f_17 * skk_1167[k]
                    - f_14 * pc_x[k] * skl1_1455[k];

        t_1456[k] = f_18 * skk_874[k]
                    + f_3 * pc_z[k] * slk_1162[k];
    }

#pragma omp simd aligned(t_1457, t_1458, t_1459, pb_x, pc_x, pc_y, skl0_1457, skl0_1458, \
                         skk_914, skk_1169, skk_1170, skl1_1457, skl1_1458, \
                         slk_1166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1457[k] = pb_x[k] * skl0_1457[k]
                    + f_17 * skk_1169[k]
                    - f_14 * pc_x[k] * skl1_1457[k];

        t_1458[k] = pb_x[k] * skl0_1458[k]
                    + f_17 * skk_1170[k]
                    - f_14 * pc_x[k] * skl1_1458[k];

        t_1459[k] = f_17 * skk_914[k]
                    + f_3 * pc_y[k] * slk_1166[k];
    }

#pragma omp simd aligned(t_1460, t_1461, t_1462, pb_x, pc_x, pc_z, skl0_1460, skl0_1461, \
                         skk_879, skk_1172, skk_1173, skl1_1460, skl1_1461, \
                         slk_1167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1460[k] = pb_x[k] * skl0_1460[k]
                    + f_17 * skk_1172[k]
                    - f_14 * pc_x[k] * skl1_1460[k];

        t_1461[k] = pb_x[k] * skl0_1461[k]
                    + f_16 * skk_1173[k]
                    - f_14 * pc_x[k] * skl1_1461[k];

        t_1462[k] = f_18 * skk_879[k]
                    + f_3 * pc_z[k] * slk_1167[k];
    }

#pragma omp simd aligned(t_1463, t_1464, t_1465, pb_x, pc_x, skl0_1463, skl0_1464, skl0_1465, \
                         skk_1175, skk_1176, skk_1177, skl1_1463, skl1_1464, \
                         skl1_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1463[k] = pb_x[k] * skl0_1463[k]
                    + f_16 * skk_1175[k]
                    - f_14 * pc_x[k] * skl1_1463[k];

        t_1464[k] = pb_x[k] * skl0_1464[k]
                    + f_16 * skk_1176[k]
                    - f_14 * pc_x[k] * skl1_1464[k];

        t_1465[k] = pb_x[k] * skl0_1465[k]
                    + f_16 * skk_1177[k]
                    - f_14 * pc_x[k] * skl1_1465[k];
    }

#pragma omp simd aligned(t_1466, t_1467, t_1468, t_1469, pb_x, pc_x, pc_y, skl0_1467, skk_920, \
                         skk_1179, skk_1180, skk_1181, skl1_1467, slk_1172, slk_1180, \
                         slk_1181 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1466[k] = f_17 * skk_920[k]
                    + f_3 * pc_y[k] * slk_1172[k];

        t_1467[k] = pb_x[k] * skl0_1467[k]
                    + f_16 * skk_1179[k]
                    - f_14 * pc_x[k] * skl1_1467[k];

        t_1468[k] = f_15 * skk_1180[k]
                    + f_3 * pc_x[k] * slk_1180[k];

        t_1469[k] = f_15 * skk_1181[k]
                    + f_3 * pc_x[k] * slk_1181[k];
    }

#pragma omp simd aligned(t_1470, t_1471, t_1472, t_1473, t_1474, pc_x, skk_1182, skk_1183, \
                         skk_1184, skk_1185, skk_1186, slk_1182, slk_1183, slk_1184, slk_1185, \
                         slk_1186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1470[k] = f_15 * skk_1182[k]
                    + f_3 * pc_x[k] * slk_1182[k];

        t_1471[k] = f_15 * skk_1183[k]
                    + f_3 * pc_x[k] * slk_1183[k];

        t_1472[k] = f_15 * skk_1184[k]
                    + f_3 * pc_x[k] * slk_1184[k];

        t_1473[k] = f_15 * skk_1185[k]
                    + f_3 * pc_x[k] * slk_1185[k];

        t_1474[k] = f_15 * skk_1186[k]
                    + f_3 * pc_x[k] * slk_1186[k];
    }

#pragma omp simd aligned(t_1475, t_1476, t_1477, t_1478, pb_x, pc_x, pc_z, skl0_1476, \
                         skl0_1478, skk_892, skk_1187, skl1_1476, skl1_1478, slk_1180, \
                         slk_1187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1475[k] = f_15 * skk_1187[k]
                    + f_3 * pc_x[k] * slk_1187[k];

        t_1476[k] = pb_x[k] * skl0_1476[k]
                    - f_14 * pc_x[k] * skl1_1476[k];

        t_1477[k] = f_18 * skk_892[k]
                    + f_3 * pc_z[k] * slk_1180[k];

        t_1478[k] = pb_x[k] * skl0_1478[k]
                    - f_14 * pc_x[k] * skl1_1478[k];
    }

#pragma omp simd aligned(t_1479, t_1480, t_1481, t_1482, pb_x, pc_x, skl0_1479, skl0_1480, \
                         skl0_1481, skl0_1482, skl1_1479, skl1_1480, skl1_1481, \
                         skl1_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1479[k] = pb_x[k] * skl0_1479[k]
                    - f_14 * pc_x[k] * skl1_1479[k];

        t_1480[k] = pb_x[k] * skl0_1480[k]
                    - f_14 * pc_x[k] * skl1_1480[k];

        t_1481[k] = pb_x[k] * skl0_1481[k]
                    - f_14 * pc_x[k] * skl1_1481[k];

        t_1482[k] = pb_x[k] * skl0_1482[k]
                    - f_14 * pc_x[k] * skl1_1482[k];
    }

#pragma omp simd aligned(t_1483, t_1484, t_1485, t_1486, pb_x, pc_x, pc_y, skl0_1484, \
                         skl0_1485, skk_935, skk_936, skk_1188, skl1_1484, skl1_1485, \
                         slk_1187, slk_1188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1483[k] = f_17 * skk_935[k]
                    + f_3 * pc_y[k] * slk_1187[k];

        t_1484[k] = pb_x[k] * skl0_1484[k]
                    - f_14 * pc_x[k] * skl1_1484[k];

        t_1485[k] = pb_x[k] * skl0_1485[k]
                    + f_0 * skk_1188[k]
                    - f_14 * pc_x[k] * skl1_1485[k];

        t_1486[k] = f_16 * skk_936[k]
                    + f_3 * pc_y[k] * slk_1188[k];
    }

#pragma omp simd aligned(t_1487, t_1488, t_1489, pb_x, pc_x, pc_y, pc_z, skl0_1488, skk_900, \
                         skk_938, skk_1191, skl1_1488, slk_1188, \
                         slk_1190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1487[k] = f_19 * skk_900[k]
                    + f_3 * pc_z[k] * slk_1188[k];

        t_1488[k] = pb_x[k] * skl0_1488[k]
                    + f_20 * skk_1191[k]
                    - f_14 * pc_x[k] * skl1_1488[k];

        t_1489[k] = f_16 * skk_938[k]
                    + f_3 * pc_y[k] * slk_1190[k];
    }

#pragma omp simd aligned(t_1490, t_1491, t_1492, pb_x, pc_x, pc_z, skl0_1490, skl0_1491, \
                         skk_903, skk_1193, skk_1194, skl1_1490, skl1_1491, \
                         slk_1191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1490[k] = pb_x[k] * skl0_1490[k]
                    + f_20 * skk_1193[k]
                    - f_14 * pc_x[k] * skl1_1490[k];

        t_1491[k] = pb_x[k] * skl0_1491[k]
                    + f_19 * skk_1194[k]
                    - f_14 * pc_x[k] * skl1_1491[k];

        t_1492[k] = f_19 * skk_903[k]
                    + f_3 * pc_z[k] * slk_1191[k];
    }

#pragma omp simd aligned(t_1493, t_1494, t_1495, pb_x, pc_x, pc_y, skl0_1494, skl0_1495, \
                         skk_941, skk_1197, skk_1198, skl1_1494, skl1_1495, \
                         slk_1193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1493[k] = f_16 * skk_941[k]
                    + f_3 * pc_y[k] * slk_1193[k];

        t_1494[k] = pb_x[k] * skl0_1494[k]
                    + f_19 * skk_1197[k]
                    - f_14 * pc_x[k] * skl1_1494[k];

        t_1495[k] = pb_x[k] * skl0_1495[k]
                    + f_18 * skk_1198[k]
                    - f_14 * pc_x[k] * skl1_1495[k];
    }

#pragma omp simd aligned(t_1496, t_1497, t_1498, pb_x, pc_x, pc_y, pc_z, skl0_1497, skk_906, \
                         skk_945, skk_1200, skl1_1497, slk_1194, \
                         slk_1197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1496[k] = f_19 * skk_906[k]
                    + f_3 * pc_z[k] * slk_1194[k];

        t_1497[k] = pb_x[k] * skl0_1497[k]
                    + f_18 * skk_1200[k]
                    - f_14 * pc_x[k] * skl1_1497[k];

        t_1498[k] = f_16 * skk_945[k]
                    + f_3 * pc_y[k] * slk_1197[k];
    }

#pragma omp simd aligned(t_1499, t_1500, t_1501, pb_x, pc_x, pc_z, skl0_1499, skl0_1500, \
                         skk_910, skk_1202, skk_1203, skl1_1499, skl1_1500, \
                         slk_1198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1499[k] = pb_x[k] * skl0_1499[k]
                    + f_18 * skk_1202[k]
                    - f_14 * pc_x[k] * skl1_1499[k];

        t_1500[k] = pb_x[k] * skl0_1500[k]
                    + f_17 * skk_1203[k]
                    - f_14 * pc_x[k] * skl1_1500[k];

        t_1501[k] = f_19 * skk_910[k]
                    + f_3 * pc_z[k] * slk_1198[k];
    }

#pragma omp simd aligned(t_1502, t_1503, t_1504, pb_x, pc_x, pc_y, skl0_1502, skl0_1503, \
                         skk_950, skk_1205, skk_1206, skl1_1502, skl1_1503, \
                         slk_1202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1502[k] = pb_x[k] * skl0_1502[k]
                    + f_17 * skk_1205[k]
                    - f_14 * pc_x[k] * skl1_1502[k];

        t_1503[k] = pb_x[k] * skl0_1503[k]
                    + f_17 * skk_1206[k]
                    - f_14 * pc_x[k] * skl1_1503[k];

        t_1504[k] = f_16 * skk_950[k]
                    + f_3 * pc_y[k] * slk_1202[k];
    }

#pragma omp simd aligned(t_1505, t_1506, t_1507, pb_x, pc_x, pc_z, skl0_1505, skl0_1506, \
                         skk_915, skk_1208, skk_1209, skl1_1505, skl1_1506, \
                         slk_1203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1505[k] = pb_x[k] * skl0_1505[k]
                    + f_17 * skk_1208[k]
                    - f_14 * pc_x[k] * skl1_1505[k];

        t_1506[k] = pb_x[k] * skl0_1506[k]
                    + f_16 * skk_1209[k]
                    - f_14 * pc_x[k] * skl1_1506[k];

        t_1507[k] = f_19 * skk_915[k]
                    + f_3 * pc_z[k] * slk_1203[k];
    }

#pragma omp simd aligned(t_1508, t_1509, t_1510, pb_x, pc_x, skl0_1508, skl0_1509, skl0_1510, \
                         skk_1211, skk_1212, skk_1213, skl1_1508, skl1_1509, \
                         skl1_1510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1508[k] = pb_x[k] * skl0_1508[k]
                    + f_16 * skk_1211[k]
                    - f_14 * pc_x[k] * skl1_1508[k];

        t_1509[k] = pb_x[k] * skl0_1509[k]
                    + f_16 * skk_1212[k]
                    - f_14 * pc_x[k] * skl1_1509[k];

        t_1510[k] = pb_x[k] * skl0_1510[k]
                    + f_16 * skk_1213[k]
                    - f_14 * pc_x[k] * skl1_1510[k];
    }

#pragma omp simd aligned(t_1511, t_1512, t_1513, t_1514, pb_x, pc_x, pc_y, skl0_1512, skk_956, \
                         skk_1215, skk_1216, skk_1217, skl1_1512, slk_1208, slk_1216, \
                         slk_1217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1511[k] = f_16 * skk_956[k]
                    + f_3 * pc_y[k] * slk_1208[k];

        t_1512[k] = pb_x[k] * skl0_1512[k]
                    + f_16 * skk_1215[k]
                    - f_14 * pc_x[k] * skl1_1512[k];

        t_1513[k] = f_15 * skk_1216[k]
                    + f_3 * pc_x[k] * slk_1216[k];

        t_1514[k] = f_15 * skk_1217[k]
                    + f_3 * pc_x[k] * slk_1217[k];
    }

#pragma omp simd aligned(t_1515, t_1516, t_1517, t_1518, t_1519, pc_x, skk_1218, skk_1219, \
                         skk_1220, skk_1221, skk_1222, slk_1218, slk_1219, slk_1220, slk_1221, \
                         slk_1222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1515[k] = f_15 * skk_1218[k]
                    + f_3 * pc_x[k] * slk_1218[k];

        t_1516[k] = f_15 * skk_1219[k]
                    + f_3 * pc_x[k] * slk_1219[k];

        t_1517[k] = f_15 * skk_1220[k]
                    + f_3 * pc_x[k] * slk_1220[k];

        t_1518[k] = f_15 * skk_1221[k]
                    + f_3 * pc_x[k] * slk_1221[k];

        t_1519[k] = f_15 * skk_1222[k]
                    + f_3 * pc_x[k] * slk_1222[k];
    }

#pragma omp simd aligned(t_1520, t_1521, t_1522, t_1523, pb_x, pc_x, pc_z, skl0_1521, \
                         skl0_1523, skk_928, skk_1223, skl1_1521, skl1_1523, slk_1216, \
                         slk_1223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1520[k] = f_15 * skk_1223[k]
                    + f_3 * pc_x[k] * slk_1223[k];

        t_1521[k] = pb_x[k] * skl0_1521[k]
                    - f_14 * pc_x[k] * skl1_1521[k];

        t_1522[k] = f_19 * skk_928[k]
                    + f_3 * pc_z[k] * slk_1216[k];

        t_1523[k] = pb_x[k] * skl0_1523[k]
                    - f_14 * pc_x[k] * skl1_1523[k];
    }

#pragma omp simd aligned(t_1524, t_1525, t_1526, t_1527, pb_x, pc_x, skl0_1524, skl0_1525, \
                         skl0_1526, skl0_1527, skl1_1524, skl1_1525, skl1_1526, \
                         skl1_1527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1524[k] = pb_x[k] * skl0_1524[k]
                    - f_14 * pc_x[k] * skl1_1524[k];

        t_1525[k] = pb_x[k] * skl0_1525[k]
                    - f_14 * pc_x[k] * skl1_1525[k];

        t_1526[k] = pb_x[k] * skl0_1526[k]
                    - f_14 * pc_x[k] * skl1_1526[k];

        t_1527[k] = pb_x[k] * skl0_1527[k]
                    - f_14 * pc_x[k] * skl1_1527[k];
    }

#pragma omp simd aligned(t_1528, t_1529, t_1530, t_1531, pb_x, pb_y, pc_x, pc_y, skl0_1215, \
                         skl0_1529, skk_971, skk_972, skl1_1215, skl1_1529, slk_1223, \
                         slk_1224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1528[k] = f_16 * skk_971[k]
                    + f_3 * pc_y[k] * slk_1223[k];

        t_1529[k] = pb_x[k] * skl0_1529[k]
                    - f_14 * pc_x[k] * skl1_1529[k];

        t_1530[k] = pb_y[k] * skl0_1215[k]
                    - f_14 * pc_y[k] * skl1_1215[k];

        t_1531[k] = f_15 * skk_972[k]
                    + f_3 * pc_y[k] * slk_1224[k];
    }

#pragma omp simd aligned(t_1532, t_1533, t_1534, pb_x, pc_x, pc_y, pc_z, skl0_1533, skk_936, \
                         skk_974, skk_1227, skl1_1533, slk_1224, \
                         slk_1226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1532[k] = f_20 * skk_936[k]
                    + f_3 * pc_z[k] * slk_1224[k];

        t_1533[k] = pb_x[k] * skl0_1533[k]
                    + f_20 * skk_1227[k]
                    - f_14 * pc_x[k] * skl1_1533[k];

        t_1534[k] = f_15 * skk_974[k]
                    + f_3 * pc_y[k] * slk_1226[k];
    }

#pragma omp simd aligned(t_1535, t_1536, t_1537, pb_x, pb_y, pc_x, pc_y, pc_z, skl0_1220, \
                         skl0_1536, skk_939, skk_1230, skl1_1220, skl1_1536, \
                         slk_1227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1535[k] = pb_y[k] * skl0_1220[k]
                    - f_14 * pc_y[k] * skl1_1220[k];

        t_1536[k] = pb_x[k] * skl0_1536[k]
                    + f_19 * skk_1230[k]
                    - f_14 * pc_x[k] * skl1_1536[k];

        t_1537[k] = f_20 * skk_939[k]
                    + f_3 * pc_z[k] * slk_1227[k];
    }

#pragma omp simd aligned(t_1538, t_1539, t_1540, pb_x, pb_y, pc_x, pc_y, skl0_1224, skl0_1540, \
                         skk_977, skk_1234, skl1_1224, skl1_1540, \
                         slk_1229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1538[k] = f_15 * skk_977[k]
                    + f_3 * pc_y[k] * slk_1229[k];

        t_1539[k] = pb_y[k] * skl0_1224[k]
                    - f_14 * pc_y[k] * skl1_1224[k];

        t_1540[k] = pb_x[k] * skl0_1540[k]
                    + f_18 * skk_1234[k]
                    - f_14 * pc_x[k] * skl1_1540[k];
    }

#pragma omp simd aligned(t_1541, t_1542, t_1543, pb_x, pc_x, pc_y, pc_z, skl0_1542, skk_942, \
                         skk_981, skk_1236, skl1_1542, slk_1230, \
                         slk_1233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1541[k] = f_20 * skk_942[k]
                    + f_3 * pc_z[k] * slk_1230[k];

        t_1542[k] = pb_x[k] * skl0_1542[k]
                    + f_18 * skk_1236[k]
                    - f_14 * pc_x[k] * skl1_1542[k];

        t_1543[k] = f_15 * skk_981[k]
                    + f_3 * pc_y[k] * slk_1233[k];
    }

#pragma omp simd aligned(t_1544, t_1545, t_1546, pb_x, pb_y, pc_x, pc_y, pc_z, skl0_1229, \
                         skl0_1545, skk_946, skk_1239, skl1_1229, skl1_1545, \
                         slk_1234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1544[k] = pb_y[k] * skl0_1229[k]
                    - f_14 * pc_y[k] * skl1_1229[k];

        t_1545[k] = pb_x[k] * skl0_1545[k]
                    + f_17 * skk_1239[k]
                    - f_14 * pc_x[k] * skl1_1545[k];

        t_1546[k] = f_20 * skk_946[k]
                    + f_3 * pc_z[k] * slk_1234[k];
    }

#pragma omp simd aligned(t_1547, t_1548, t_1549, pb_x, pc_x, pc_y, skl0_1547, skl0_1548, \
                         skk_986, skk_1241, skk_1242, skl1_1547, skl1_1548, \
                         slk_1238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1547[k] = pb_x[k] * skl0_1547[k]
                    + f_17 * skk_1241[k]
                    - f_14 * pc_x[k] * skl1_1547[k];

        t_1548[k] = pb_x[k] * skl0_1548[k]
                    + f_17 * skk_1242[k]
                    - f_14 * pc_x[k] * skl1_1548[k];

        t_1549[k] = f_15 * skk_986[k]
                    + f_3 * pc_y[k] * slk_1238[k];
    }

#pragma omp simd aligned(t_1550, t_1551, t_1552, pb_x, pb_y, pc_x, pc_y, pc_z, skl0_1235, \
                         skl0_1551, skk_951, skk_1245, skl1_1235, skl1_1551, \
                         slk_1239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1550[k] = pb_y[k] * skl0_1235[k]
                    - f_14 * pc_y[k] * skl1_1235[k];

        t_1551[k] = pb_x[k] * skl0_1551[k]
                    + f_16 * skk_1245[k]
                    - f_14 * pc_x[k] * skl1_1551[k];

        t_1552[k] = f_20 * skk_951[k]
                    + f_3 * pc_z[k] * slk_1239[k];
    }

#pragma omp simd aligned(t_1553, t_1554, t_1555, pb_x, pc_x, skl0_1553, skl0_1554, skl0_1555, \
                         skk_1247, skk_1248, skk_1249, skl1_1553, skl1_1554, \
                         skl1_1555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1553[k] = pb_x[k] * skl0_1553[k]
                    + f_16 * skk_1247[k]
                    - f_14 * pc_x[k] * skl1_1553[k];

        t_1554[k] = pb_x[k] * skl0_1554[k]
                    + f_16 * skk_1248[k]
                    - f_14 * pc_x[k] * skl1_1554[k];

        t_1555[k] = pb_x[k] * skl0_1555[k]
                    + f_16 * skk_1249[k]
                    - f_14 * pc_x[k] * skl1_1555[k];
    }

#pragma omp simd aligned(t_1556, t_1557, t_1558, t_1559, pb_y, pc_x, pc_y, skl0_1242, skk_992, \
                         skk_1252, skk_1253, skl1_1242, slk_1244, slk_1252, \
                         slk_1253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1556[k] = f_15 * skk_992[k]
                    + f_3 * pc_y[k] * slk_1244[k];

        t_1557[k] = pb_y[k] * skl0_1242[k]
                    - f_14 * pc_y[k] * skl1_1242[k];

        t_1558[k] = f_15 * skk_1252[k]
                    + f_3 * pc_x[k] * slk_1252[k];

        t_1559[k] = f_15 * skk_1253[k]
                    + f_3 * pc_x[k] * slk_1253[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece14(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t sli0, const size_t sli1,
                                                           const size_t slk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_21 = 3.5 / q;

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
    auto *t_1620 = buffer.data(target + 1620);
    auto *t_1621 = buffer.data(target + 1621);
    auto *t_1622 = buffer.data(target + 1622);
    auto *t_1623 = buffer.data(target + 1623);
    auto *t_1624 = buffer.data(target + 1624);
    auto *t_1625 = buffer.data(target + 1625);
    auto *t_1626 = buffer.data(target + 1626);
    auto *t_1627 = buffer.data(target + 1627);
    auto *t_1628 = buffer.data(target + 1628);
    auto *t_1629 = buffer.data(target + 1629);
    auto *t_1630 = buffer.data(target + 1630);
    auto *t_1631 = buffer.data(target + 1631);
    auto *t_1632 = buffer.data(target + 1632);
    auto *t_1633 = buffer.data(target + 1633);
    auto *t_1634 = buffer.data(target + 1634);
    auto *t_1635 = buffer.data(target + 1635);
    auto *t_1636 = buffer.data(target + 1636);
    auto *t_1637 = buffer.data(target + 1637);
    auto *t_1638 = buffer.data(target + 1638);
    auto *t_1639 = buffer.data(target + 1639);
    auto *t_1640 = buffer.data(target + 1640);
    auto *t_1641 = buffer.data(target + 1641);
    auto *t_1642 = buffer.data(target + 1642);
    auto *t_1643 = buffer.data(target + 1643);
    auto *t_1644 = buffer.data(target + 1644);
    auto *t_1645 = buffer.data(target + 1645);
    auto *t_1646 = buffer.data(target + 1646);
    auto *t_1647 = buffer.data(target + 1647);
    auto *t_1648 = buffer.data(target + 1648);
    auto *t_1649 = buffer.data(target + 1649);
    auto *t_1650 = buffer.data(target + 1650);
    auto *t_1651 = buffer.data(target + 1651);
    auto *t_1652 = buffer.data(target + 1652);
    auto *t_1653 = buffer.data(target + 1653);
    auto *t_1654 = buffer.data(target + 1654);
    auto *t_1655 = buffer.data(target + 1655);
    auto *t_1656 = buffer.data(target + 1656);
    auto *t_1657 = buffer.data(target + 1657);
    auto *t_1658 = buffer.data(target + 1658);
    auto *t_1659 = buffer.data(target + 1659);
    auto *t_1660 = buffer.data(target + 1660);
    auto *t_1661 = buffer.data(target + 1661);
    auto *t_1662 = buffer.data(target + 1662);
    auto *t_1663 = buffer.data(target + 1663);
    auto *t_1664 = buffer.data(target + 1664);
    auto *t_1665 = buffer.data(target + 1665);
    auto *t_1666 = buffer.data(target + 1666);
    auto *t_1667 = buffer.data(target + 1667);
    auto *t_1668 = buffer.data(target + 1668);
    auto *t_1669 = buffer.data(target + 1669);
    auto *t_1670 = buffer.data(target + 1670);
    auto *t_1671 = buffer.data(target + 1671);
    auto *t_1672 = buffer.data(target + 1672);
    auto *t_1673 = buffer.data(target + 1673);
    auto *t_1674 = buffer.data(target + 1674);
    auto *t_1675 = buffer.data(target + 1675);
    auto *t_1676 = buffer.data(target + 1676);
    auto *t_1677 = buffer.data(target + 1677);
    auto *t_1678 = buffer.data(target + 1678);
    auto *t_1679 = buffer.data(target + 1679);
    auto *t_1680 = buffer.data(target + 1680);
    auto *t_1681 = buffer.data(target + 1681);
    auto *t_1682 = buffer.data(target + 1682);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_1260 = buffer.data(skl0 + 1260);
    const auto *skl0_1263 = buffer.data(skl0 + 1263);
    const auto *skl0_1266 = buffer.data(skl0 + 1266);
    const auto *skl0_1270 = buffer.data(skl0 + 1270);
    const auto *skl0_1275 = buffer.data(skl0 + 1275);
    const auto *skl0_1566 = buffer.data(skl0 + 1566);
    const auto *skl0_1568 = buffer.data(skl0 + 1568);
    const auto *skl0_1569 = buffer.data(skl0 + 1569);
    const auto *skl0_1570 = buffer.data(skl0 + 1570);
    const auto *skl0_1571 = buffer.data(skl0 + 1571);
    const auto *skl0_1572 = buffer.data(skl0 + 1572);
    const auto *skl0_1574 = buffer.data(skl0 + 1574);
    const auto *skl0_1575 = buffer.data(skl0 + 1575);
    const auto *skl0_1578 = buffer.data(skl0 + 1578);
    const auto *skl0_1580 = buffer.data(skl0 + 1580);
    const auto *skl0_1581 = buffer.data(skl0 + 1581);
    const auto *skl0_1584 = buffer.data(skl0 + 1584);
    const auto *skl0_1585 = buffer.data(skl0 + 1585);
    const auto *skl0_1587 = buffer.data(skl0 + 1587);
    const auto *skl0_1589 = buffer.data(skl0 + 1589);
    const auto *skl0_1590 = buffer.data(skl0 + 1590);
    const auto *skl0_1592 = buffer.data(skl0 + 1592);
    const auto *skl0_1593 = buffer.data(skl0 + 1593);
    const auto *skl0_1595 = buffer.data(skl0 + 1595);
    const auto *skl0_1596 = buffer.data(skl0 + 1596);
    const auto *skl0_1598 = buffer.data(skl0 + 1598);
    const auto *skl0_1599 = buffer.data(skl0 + 1599);
    const auto *skl0_1600 = buffer.data(skl0 + 1600);
    const auto *skl0_1602 = buffer.data(skl0 + 1602);
    const auto *skl0_1611 = buffer.data(skl0 + 1611);
    const auto *skl0_1613 = buffer.data(skl0 + 1613);
    const auto *skl0_1614 = buffer.data(skl0 + 1614);
    const auto *skl0_1615 = buffer.data(skl0 + 1615);
    const auto *skl0_1616 = buffer.data(skl0 + 1616);
    const auto *skl0_1617 = buffer.data(skl0 + 1617);
    const auto *skl0_1619 = buffer.data(skl0 + 1619);

    const auto *skk_964 = buffer.data(skk + 964);
    const auto *skk_972 = buffer.data(skk + 972);
    const auto *skk_975 = buffer.data(skk + 975);
    const auto *skk_978 = buffer.data(skk + 978);
    const auto *skk_982 = buffer.data(skk + 982);
    const auto *skk_987 = buffer.data(skk + 987);
    const auto *skk_1000 = buffer.data(skk + 1000);
    const auto *skk_1007 = buffer.data(skk + 1007);
    const auto *skk_1008 = buffer.data(skk + 1008);
    const auto *skk_1010 = buffer.data(skk + 1010);
    const auto *skk_1011 = buffer.data(skk + 1011);
    const auto *skk_1013 = buffer.data(skk + 1013);
    const auto *skk_1014 = buffer.data(skk + 1014);
    const auto *skk_1017 = buffer.data(skk + 1017);
    const auto *skk_1018 = buffer.data(skk + 1018);
    const auto *skk_1022 = buffer.data(skk + 1022);
    const auto *skk_1028 = buffer.data(skk + 1028);
    const auto *skk_1036 = buffer.data(skk + 1036);
    const auto *skk_1038 = buffer.data(skk + 1038);
    const auto *skk_1039 = buffer.data(skk + 1039);
    const auto *skk_1040 = buffer.data(skk + 1040);
    const auto *skk_1041 = buffer.data(skk + 1041);
    const auto *skk_1042 = buffer.data(skk + 1042);
    const auto *skk_1043 = buffer.data(skk + 1043);
    const auto *skk_1044 = buffer.data(skk + 1044);
    const auto *skk_1046 = buffer.data(skk + 1046);
    const auto *skk_1049 = buffer.data(skk + 1049);
    const auto *skk_1053 = buffer.data(skk + 1053);
    const auto *skk_1254 = buffer.data(skk + 1254);
    const auto *skk_1255 = buffer.data(skk + 1255);
    const auto *skk_1256 = buffer.data(skk + 1256);
    const auto *skk_1257 = buffer.data(skk + 1257);
    const auto *skk_1258 = buffer.data(skk + 1258);
    const auto *skk_1259 = buffer.data(skk + 1259);
    const auto *skk_1260 = buffer.data(skk + 1260);
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

    const auto *skl1_1260 = buffer.data(skl1 + 1260);
    const auto *skl1_1263 = buffer.data(skl1 + 1263);
    const auto *skl1_1266 = buffer.data(skl1 + 1266);
    const auto *skl1_1270 = buffer.data(skl1 + 1270);
    const auto *skl1_1275 = buffer.data(skl1 + 1275);
    const auto *skl1_1566 = buffer.data(skl1 + 1566);
    const auto *skl1_1568 = buffer.data(skl1 + 1568);
    const auto *skl1_1569 = buffer.data(skl1 + 1569);
    const auto *skl1_1570 = buffer.data(skl1 + 1570);
    const auto *skl1_1571 = buffer.data(skl1 + 1571);
    const auto *skl1_1572 = buffer.data(skl1 + 1572);
    const auto *skl1_1574 = buffer.data(skl1 + 1574);
    const auto *skl1_1575 = buffer.data(skl1 + 1575);
    const auto *skl1_1578 = buffer.data(skl1 + 1578);
    const auto *skl1_1580 = buffer.data(skl1 + 1580);
    const auto *skl1_1581 = buffer.data(skl1 + 1581);
    const auto *skl1_1584 = buffer.data(skl1 + 1584);
    const auto *skl1_1585 = buffer.data(skl1 + 1585);
    const auto *skl1_1587 = buffer.data(skl1 + 1587);
    const auto *skl1_1589 = buffer.data(skl1 + 1589);
    const auto *skl1_1590 = buffer.data(skl1 + 1590);
    const auto *skl1_1592 = buffer.data(skl1 + 1592);
    const auto *skl1_1593 = buffer.data(skl1 + 1593);
    const auto *skl1_1595 = buffer.data(skl1 + 1595);
    const auto *skl1_1596 = buffer.data(skl1 + 1596);
    const auto *skl1_1598 = buffer.data(skl1 + 1598);
    const auto *skl1_1599 = buffer.data(skl1 + 1599);
    const auto *skl1_1600 = buffer.data(skl1 + 1600);
    const auto *skl1_1602 = buffer.data(skl1 + 1602);
    const auto *skl1_1611 = buffer.data(skl1 + 1611);
    const auto *skl1_1613 = buffer.data(skl1 + 1613);
    const auto *skl1_1614 = buffer.data(skl1 + 1614);
    const auto *skl1_1615 = buffer.data(skl1 + 1615);
    const auto *skl1_1616 = buffer.data(skl1 + 1616);
    const auto *skl1_1617 = buffer.data(skl1 + 1617);
    const auto *skl1_1619 = buffer.data(skl1 + 1619);

    const auto *sli0_1008 = buffer.data(sli0 + 1008);
    const auto *sli0_1011 = buffer.data(sli0 + 1011);
    const auto *sli0_1013 = buffer.data(sli0 + 1013);
    const auto *sli0_1014 = buffer.data(sli0 + 1014);
    const auto *sli0_1017 = buffer.data(sli0 + 1017);
    const auto *sli0_1018 = buffer.data(sli0 + 1018);
    const auto *sli0_1020 = buffer.data(sli0 + 1020);
    const auto *sli0_1022 = buffer.data(sli0 + 1022);
    const auto *sli0_1023 = buffer.data(sli0 + 1023);
    const auto *sli0_1025 = buffer.data(sli0 + 1025);
    const auto *sli0_1026 = buffer.data(sli0 + 1026);
    const auto *sli0_1028 = buffer.data(sli0 + 1028);
    const auto *sli0_1029 = buffer.data(sli0 + 1029);
    const auto *sli0_1031 = buffer.data(sli0 + 1031);
    const auto *sli0_1032 = buffer.data(sli0 + 1032);
    const auto *sli0_1033 = buffer.data(sli0 + 1033);
    const auto *sli0_1034 = buffer.data(sli0 + 1034);
    const auto *sli0_1035 = buffer.data(sli0 + 1035);
    const auto *sli0_1041 = buffer.data(sli0 + 1041);
    const auto *sli0_1045 = buffer.data(sli0 + 1045);
    const auto *sli0_1048 = buffer.data(sli0 + 1048);
    const auto *sli0_1050 = buffer.data(sli0 + 1050);
    const auto *sli0_1053 = buffer.data(sli0 + 1053);

    const auto *sli1_1008 = buffer.data(sli1 + 1008);
    const auto *sli1_1011 = buffer.data(sli1 + 1011);
    const auto *sli1_1013 = buffer.data(sli1 + 1013);
    const auto *sli1_1014 = buffer.data(sli1 + 1014);
    const auto *sli1_1017 = buffer.data(sli1 + 1017);
    const auto *sli1_1018 = buffer.data(sli1 + 1018);
    const auto *sli1_1020 = buffer.data(sli1 + 1020);
    const auto *sli1_1022 = buffer.data(sli1 + 1022);
    const auto *sli1_1023 = buffer.data(sli1 + 1023);
    const auto *sli1_1025 = buffer.data(sli1 + 1025);
    const auto *sli1_1026 = buffer.data(sli1 + 1026);
    const auto *sli1_1028 = buffer.data(sli1 + 1028);
    const auto *sli1_1029 = buffer.data(sli1 + 1029);
    const auto *sli1_1031 = buffer.data(sli1 + 1031);
    const auto *sli1_1032 = buffer.data(sli1 + 1032);
    const auto *sli1_1033 = buffer.data(sli1 + 1033);
    const auto *sli1_1034 = buffer.data(sli1 + 1034);
    const auto *sli1_1035 = buffer.data(sli1 + 1035);
    const auto *sli1_1041 = buffer.data(sli1 + 1041);
    const auto *sli1_1045 = buffer.data(sli1 + 1045);
    const auto *sli1_1048 = buffer.data(sli1 + 1048);
    const auto *sli1_1050 = buffer.data(sli1 + 1050);
    const auto *sli1_1053 = buffer.data(sli1 + 1053);

    const auto *slk_1252 = buffer.data(slk + 1252);
    const auto *slk_1254 = buffer.data(slk + 1254);
    const auto *slk_1255 = buffer.data(slk + 1255);
    const auto *slk_1256 = buffer.data(slk + 1256);
    const auto *slk_1257 = buffer.data(slk + 1257);
    const auto *slk_1258 = buffer.data(slk + 1258);
    const auto *slk_1259 = buffer.data(slk + 1259);
    const auto *slk_1260 = buffer.data(slk + 1260);
    const auto *slk_1262 = buffer.data(slk + 1262);
    const auto *slk_1263 = buffer.data(slk + 1263);
    const auto *slk_1265 = buffer.data(slk + 1265);
    const auto *slk_1266 = buffer.data(slk + 1266);
    const auto *slk_1269 = buffer.data(slk + 1269);
    const auto *slk_1270 = buffer.data(slk + 1270);
    const auto *slk_1274 = buffer.data(slk + 1274);
    const auto *slk_1275 = buffer.data(slk + 1275);
    const auto *slk_1280 = buffer.data(slk + 1280);
    const auto *slk_1288 = buffer.data(slk + 1288);
    const auto *slk_1289 = buffer.data(slk + 1289);
    const auto *slk_1290 = buffer.data(slk + 1290);
    const auto *slk_1291 = buffer.data(slk + 1291);
    const auto *slk_1292 = buffer.data(slk + 1292);
    const auto *slk_1293 = buffer.data(slk + 1293);
    const auto *slk_1294 = buffer.data(slk + 1294);
    const auto *slk_1295 = buffer.data(slk + 1295);
    const auto *slk_1296 = buffer.data(slk + 1296);
    const auto *slk_1298 = buffer.data(slk + 1298);
    const auto *slk_1299 = buffer.data(slk + 1299);
    const auto *slk_1301 = buffer.data(slk + 1301);
    const auto *slk_1302 = buffer.data(slk + 1302);
    const auto *slk_1305 = buffer.data(slk + 1305);
    const auto *slk_1306 = buffer.data(slk + 1306);
    const auto *slk_1308 = buffer.data(slk + 1308);
    const auto *slk_1310 = buffer.data(slk + 1310);
    const auto *slk_1311 = buffer.data(slk + 1311);
    const auto *slk_1313 = buffer.data(slk + 1313);
    const auto *slk_1314 = buffer.data(slk + 1314);
    const auto *slk_1316 = buffer.data(slk + 1316);
    const auto *slk_1317 = buffer.data(slk + 1317);
    const auto *slk_1319 = buffer.data(slk + 1319);
    const auto *slk_1320 = buffer.data(slk + 1320);
    const auto *slk_1321 = buffer.data(slk + 1321);
    const auto *slk_1323 = buffer.data(slk + 1323);
    const auto *slk_1324 = buffer.data(slk + 1324);
    const auto *slk_1325 = buffer.data(slk + 1325);
    const auto *slk_1326 = buffer.data(slk + 1326);
    const auto *slk_1327 = buffer.data(slk + 1327);
    const auto *slk_1328 = buffer.data(slk + 1328);
    const auto *slk_1329 = buffer.data(slk + 1329);
    const auto *slk_1330 = buffer.data(slk + 1330);
    const auto *slk_1331 = buffer.data(slk + 1331);
    const auto *slk_1332 = buffer.data(slk + 1332);
    const auto *slk_1334 = buffer.data(slk + 1334);
    const auto *slk_1335 = buffer.data(slk + 1335);
    const auto *slk_1337 = buffer.data(slk + 1337);
    const auto *slk_1338 = buffer.data(slk + 1338);
    const auto *slk_1341 = buffer.data(slk + 1341);
    const auto *slk_1342 = buffer.data(slk + 1342);
    const auto *slk_1344 = buffer.data(slk + 1344);
    const auto *slk_1346 = buffer.data(slk + 1346);
    const auto *slk_1349 = buffer.data(slk + 1349);

#pragma omp simd aligned(t_1560, t_1561, t_1562, t_1563, t_1564, pc_x, skk_1254, skk_1255, \
                         skk_1256, skk_1257, skk_1258, slk_1254, slk_1255, slk_1256, slk_1257, \
                         slk_1258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1560[k] = f_15 * skk_1254[k]
                    + f_3 * pc_x[k] * slk_1254[k];

        t_1561[k] = f_15 * skk_1255[k]
                    + f_3 * pc_x[k] * slk_1255[k];

        t_1562[k] = f_15 * skk_1256[k]
                    + f_3 * pc_x[k] * slk_1256[k];

        t_1563[k] = f_15 * skk_1257[k]
                    + f_3 * pc_x[k] * slk_1257[k];

        t_1564[k] = f_15 * skk_1258[k]
                    + f_3 * pc_x[k] * slk_1258[k];
    }

#pragma omp simd aligned(t_1565, t_1566, t_1567, t_1568, pb_x, pc_x, pc_z, skl0_1566, \
                         skl0_1568, skk_964, skk_1259, skl1_1566, skl1_1568, slk_1252, \
                         slk_1259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1565[k] = f_15 * skk_1259[k]
                    + f_3 * pc_x[k] * slk_1259[k];

        t_1566[k] = pb_x[k] * skl0_1566[k]
                    - f_14 * pc_x[k] * skl1_1566[k];

        t_1567[k] = f_20 * skk_964[k]
                    + f_3 * pc_z[k] * slk_1252[k];

        t_1568[k] = pb_x[k] * skl0_1568[k]
                    - f_14 * pc_x[k] * skl1_1568[k];
    }

#pragma omp simd aligned(t_1569, t_1570, t_1571, t_1572, pb_x, pc_x, skl0_1569, skl0_1570, \
                         skl0_1571, skl0_1572, skl1_1569, skl1_1570, skl1_1571, \
                         skl1_1572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1569[k] = pb_x[k] * skl0_1569[k]
                    - f_14 * pc_x[k] * skl1_1569[k];

        t_1570[k] = pb_x[k] * skl0_1570[k]
                    - f_14 * pc_x[k] * skl1_1570[k];

        t_1571[k] = pb_x[k] * skl0_1571[k]
                    - f_14 * pc_x[k] * skl1_1571[k];

        t_1572[k] = pb_x[k] * skl0_1572[k]
                    - f_14 * pc_x[k] * skl1_1572[k];
    }

#pragma omp simd aligned(t_1573, t_1574, t_1575, t_1576, pb_x, pc_x, pc_y, skl0_1574, \
                         skl0_1575, skk_1007, skk_1260, skl1_1574, skl1_1575, slk_1259, \
                         slk_1260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1573[k] = f_15 * skk_1007[k]
                    + f_3 * pc_y[k] * slk_1259[k];

        t_1574[k] = pb_x[k] * skl0_1574[k]
                    - f_14 * pc_x[k] * skl1_1574[k];

        t_1575[k] = pb_x[k] * skl0_1575[k]
                    + f_0 * skk_1260[k]
                    - f_14 * pc_x[k] * skl1_1575[k];

        t_1576[k] = f_3 * pc_y[k] * slk_1260[k];
    }

#pragma omp simd aligned(t_1577, t_1578, t_1579, pb_x, pc_x, pc_y, pc_z, skl0_1578, skk_972, \
                         skk_1263, skl1_1578, slk_1260, slk_1262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1577[k] = f_21 * skk_972[k]
                    + f_3 * pc_z[k] * slk_1260[k];

        t_1578[k] = pb_x[k] * skl0_1578[k]
                    + f_20 * skk_1263[k]
                    - f_14 * pc_x[k] * skl1_1578[k];

        t_1579[k] = f_3 * pc_y[k] * slk_1262[k];
    }

#pragma omp simd aligned(t_1580, t_1581, t_1582, pb_x, pc_x, pc_z, skl0_1580, skl0_1581, \
                         skk_975, skk_1265, skk_1266, skl1_1580, skl1_1581, \
                         slk_1263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1580[k] = pb_x[k] * skl0_1580[k]
                    + f_20 * skk_1265[k]
                    - f_14 * pc_x[k] * skl1_1580[k];

        t_1581[k] = pb_x[k] * skl0_1581[k]
                    + f_19 * skk_1266[k]
                    - f_14 * pc_x[k] * skl1_1581[k];

        t_1582[k] = f_21 * skk_975[k]
                    + f_3 * pc_z[k] * slk_1263[k];
    }

#pragma omp simd aligned(t_1583, t_1584, t_1585, pb_x, pc_x, pc_y, skl0_1584, skl0_1585, \
                         skk_1269, skk_1270, skl1_1584, skl1_1585, \
                         slk_1265 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1583[k] = f_3 * pc_y[k] * slk_1265[k];

        t_1584[k] = pb_x[k] * skl0_1584[k]
                    + f_19 * skk_1269[k]
                    - f_14 * pc_x[k] * skl1_1584[k];

        t_1585[k] = pb_x[k] * skl0_1585[k]
                    + f_18 * skk_1270[k]
                    - f_14 * pc_x[k] * skl1_1585[k];
    }

#pragma omp simd aligned(t_1586, t_1587, t_1588, pb_x, pc_x, pc_y, pc_z, skl0_1587, skk_978, \
                         skk_1272, skl1_1587, slk_1266, slk_1269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1586[k] = f_21 * skk_978[k]
                    + f_3 * pc_z[k] * slk_1266[k];

        t_1587[k] = pb_x[k] * skl0_1587[k]
                    + f_18 * skk_1272[k]
                    - f_14 * pc_x[k] * skl1_1587[k];

        t_1588[k] = f_3 * pc_y[k] * slk_1269[k];
    }

#pragma omp simd aligned(t_1589, t_1590, t_1591, pb_x, pc_x, pc_z, skl0_1589, skl0_1590, \
                         skk_982, skk_1274, skk_1275, skl1_1589, skl1_1590, \
                         slk_1270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1589[k] = pb_x[k] * skl0_1589[k]
                    + f_18 * skk_1274[k]
                    - f_14 * pc_x[k] * skl1_1589[k];

        t_1590[k] = pb_x[k] * skl0_1590[k]
                    + f_17 * skk_1275[k]
                    - f_14 * pc_x[k] * skl1_1590[k];

        t_1591[k] = f_21 * skk_982[k]
                    + f_3 * pc_z[k] * slk_1270[k];
    }

#pragma omp simd aligned(t_1592, t_1593, t_1594, pb_x, pc_x, pc_y, skl0_1592, skl0_1593, \
                         skk_1277, skk_1278, skl1_1592, skl1_1593, \
                         slk_1274 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1592[k] = pb_x[k] * skl0_1592[k]
                    + f_17 * skk_1277[k]
                    - f_14 * pc_x[k] * skl1_1592[k];

        t_1593[k] = pb_x[k] * skl0_1593[k]
                    + f_17 * skk_1278[k]
                    - f_14 * pc_x[k] * skl1_1593[k];

        t_1594[k] = f_3 * pc_y[k] * slk_1274[k];
    }

#pragma omp simd aligned(t_1595, t_1596, t_1597, pb_x, pc_x, pc_z, skl0_1595, skl0_1596, \
                         skk_987, skk_1280, skk_1281, skl1_1595, skl1_1596, \
                         slk_1275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1595[k] = pb_x[k] * skl0_1595[k]
                    + f_17 * skk_1280[k]
                    - f_14 * pc_x[k] * skl1_1595[k];

        t_1596[k] = pb_x[k] * skl0_1596[k]
                    + f_16 * skk_1281[k]
                    - f_14 * pc_x[k] * skl1_1596[k];

        t_1597[k] = f_21 * skk_987[k]
                    + f_3 * pc_z[k] * slk_1275[k];
    }

#pragma omp simd aligned(t_1598, t_1599, t_1600, pb_x, pc_x, skl0_1598, skl0_1599, skl0_1600, \
                         skk_1283, skk_1284, skk_1285, skl1_1598, skl1_1599, \
                         skl1_1600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1598[k] = pb_x[k] * skl0_1598[k]
                    + f_16 * skk_1283[k]
                    - f_14 * pc_x[k] * skl1_1598[k];

        t_1599[k] = pb_x[k] * skl0_1599[k]
                    + f_16 * skk_1284[k]
                    - f_14 * pc_x[k] * skl1_1599[k];

        t_1600[k] = pb_x[k] * skl0_1600[k]
                    + f_16 * skk_1285[k]
                    - f_14 * pc_x[k] * skl1_1600[k];
    }

#pragma omp simd aligned(t_1601, t_1602, t_1603, t_1604, pb_x, pc_x, pc_y, skl0_1602, \
                         skk_1287, skk_1288, skk_1289, skl1_1602, slk_1280, slk_1288, \
                         slk_1289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1601[k] = f_3 * pc_y[k] * slk_1280[k];

        t_1602[k] = pb_x[k] * skl0_1602[k]
                    + f_16 * skk_1287[k]
                    - f_14 * pc_x[k] * skl1_1602[k];

        t_1603[k] = f_15 * skk_1288[k]
                    + f_3 * pc_x[k] * slk_1288[k];

        t_1604[k] = f_15 * skk_1289[k]
                    + f_3 * pc_x[k] * slk_1289[k];
    }

#pragma omp simd aligned(t_1605, t_1606, t_1607, t_1608, t_1609, pc_x, skk_1290, skk_1291, \
                         skk_1292, skk_1293, skk_1294, slk_1290, slk_1291, slk_1292, slk_1293, \
                         slk_1294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1605[k] = f_15 * skk_1290[k]
                    + f_3 * pc_x[k] * slk_1290[k];

        t_1606[k] = f_15 * skk_1291[k]
                    + f_3 * pc_x[k] * slk_1291[k];

        t_1607[k] = f_15 * skk_1292[k]
                    + f_3 * pc_x[k] * slk_1292[k];

        t_1608[k] = f_15 * skk_1293[k]
                    + f_3 * pc_x[k] * slk_1293[k];

        t_1609[k] = f_15 * skk_1294[k]
                    + f_3 * pc_x[k] * slk_1294[k];
    }

#pragma omp simd aligned(t_1610, t_1611, t_1612, t_1613, pb_x, pc_x, pc_z, skl0_1611, \
                         skl0_1613, skk_1000, skk_1295, skl1_1611, skl1_1613, slk_1288, \
                         slk_1295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1610[k] = f_15 * skk_1295[k]
                    + f_3 * pc_x[k] * slk_1295[k];

        t_1611[k] = pb_x[k] * skl0_1611[k]
                    - f_14 * pc_x[k] * skl1_1611[k];

        t_1612[k] = f_21 * skk_1000[k]
                    + f_3 * pc_z[k] * slk_1288[k];

        t_1613[k] = pb_x[k] * skl0_1613[k]
                    - f_14 * pc_x[k] * skl1_1613[k];
    }

#pragma omp simd aligned(t_1614, t_1615, t_1616, t_1617, pb_x, pc_x, skl0_1614, skl0_1615, \
                         skl0_1616, skl0_1617, skl1_1614, skl1_1615, skl1_1616, \
                         skl1_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1614[k] = pb_x[k] * skl0_1614[k]
                    - f_14 * pc_x[k] * skl1_1614[k];

        t_1615[k] = pb_x[k] * skl0_1615[k]
                    - f_14 * pc_x[k] * skl1_1615[k];

        t_1616[k] = pb_x[k] * skl0_1616[k]
                    - f_14 * pc_x[k] * skl1_1616[k];

        t_1617[k] = pb_x[k] * skl0_1617[k]
                    - f_14 * pc_x[k] * skl1_1617[k];
    }

#pragma omp simd aligned(t_1618, t_1619, t_1620, t_1621, t_1622, pb_x, pc_x, pc_y, pc_z, \
                         skl0_1619, skk_1008, skl1_1619, sli0_1008, sli1_1008, slk_1295, \
                         slk_1296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1618[k] = f_3 * pc_y[k] * slk_1295[k];

        t_1619[k] = pb_x[k] * skl0_1619[k]
                    - f_14 * pc_x[k] * skl1_1619[k];

        t_1620[k] = f_1 * sli0_1008[k]
                    - f_2 * sli1_1008[k]
                    + f_3 * pc_x[k] * slk_1296[k];

        t_1621[k] = f_0 * skk_1008[k]
                    + f_3 * pc_y[k] * slk_1296[k];

        t_1622[k] = f_3 * pc_z[k] * slk_1296[k];
    }

#pragma omp simd aligned(t_1623, t_1624, t_1625, pc_x, pc_y, skk_1010, sli0_1011, sli0_1013, \
                         sli1_1011, sli1_1013, slk_1298, slk_1299, \
                         slk_1301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1623[k] = f_4 * sli0_1011[k]
                    - f_5 * sli1_1011[k]
                    + f_3 * pc_x[k] * slk_1299[k];

        t_1624[k] = f_0 * skk_1010[k]
                    + f_3 * pc_y[k] * slk_1298[k];

        t_1625[k] = f_4 * sli0_1013[k]
                    - f_5 * sli1_1013[k]
                    + f_3 * pc_x[k] * slk_1301[k];
    }

#pragma omp simd aligned(t_1626, t_1627, t_1628, t_1629, pc_x, pc_y, pc_z, skk_1013, \
                         sli0_1014, sli0_1017, sli1_1014, sli1_1017, slk_1299, slk_1301, \
                         slk_1302, slk_1305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1626[k] = f_6 * sli0_1014[k]
                    - f_7 * sli1_1014[k]
                    + f_3 * pc_x[k] * slk_1302[k];

        t_1627[k] = f_3 * pc_z[k] * slk_1299[k];

        t_1628[k] = f_0 * skk_1013[k]
                    + f_3 * pc_y[k] * slk_1301[k];

        t_1629[k] = f_6 * sli0_1017[k]
                    - f_7 * sli1_1017[k]
                    + f_3 * pc_x[k] * slk_1305[k];
    }

#pragma omp simd aligned(t_1630, t_1631, t_1632, t_1633, pc_x, pc_y, pc_z, skk_1017, \
                         sli0_1018, sli0_1020, sli1_1018, sli1_1020, slk_1302, slk_1305, \
                         slk_1306, slk_1308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1630[k] = f_8 * sli0_1018[k]
                    - f_9 * sli1_1018[k]
                    + f_3 * pc_x[k] * slk_1306[k];

        t_1631[k] = f_3 * pc_z[k] * slk_1302[k];

        t_1632[k] = f_8 * sli0_1020[k]
                    - f_9 * sli1_1020[k]
                    + f_3 * pc_x[k] * slk_1308[k];

        t_1633[k] = f_0 * skk_1017[k]
                    + f_3 * pc_y[k] * slk_1305[k];
    }

#pragma omp simd aligned(t_1634, t_1635, t_1636, t_1637, pc_x, pc_z, sli0_1022, sli0_1023, \
                         sli0_1025, sli1_1022, sli1_1023, sli1_1025, slk_1306, slk_1310, \
                         slk_1311, slk_1313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1634[k] = f_8 * sli0_1022[k]
                    - f_9 * sli1_1022[k]
                    + f_3 * pc_x[k] * slk_1310[k];

        t_1635[k] = f_10 * sli0_1023[k]
                    - f_11 * sli1_1023[k]
                    + f_3 * pc_x[k] * slk_1311[k];

        t_1636[k] = f_3 * pc_z[k] * slk_1306[k];

        t_1637[k] = f_10 * sli0_1025[k]
                    - f_11 * sli1_1025[k]
                    + f_3 * pc_x[k] * slk_1313[k];
    }

#pragma omp simd aligned(t_1638, t_1639, t_1640, pc_x, pc_y, skk_1022, sli0_1026, sli0_1028, \
                         sli1_1026, sli1_1028, slk_1310, slk_1314, \
                         slk_1316 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1638[k] = f_10 * sli0_1026[k]
                    - f_11 * sli1_1026[k]
                    + f_3 * pc_x[k] * slk_1314[k];

        t_1639[k] = f_0 * skk_1022[k]
                    + f_3 * pc_y[k] * slk_1310[k];

        t_1640[k] = f_10 * sli0_1028[k]
                    - f_11 * sli1_1028[k]
                    + f_3 * pc_x[k] * slk_1316[k];
    }

#pragma omp simd aligned(t_1641, t_1642, t_1643, t_1644, pc_x, pc_z, sli0_1029, sli0_1031, \
                         sli0_1032, sli1_1029, sli1_1031, sli1_1032, slk_1311, slk_1317, \
                         slk_1319, slk_1320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1641[k] = f_12 * sli0_1029[k]
                    - f_13 * sli1_1029[k]
                    + f_3 * pc_x[k] * slk_1317[k];

        t_1642[k] = f_3 * pc_z[k] * slk_1311[k];

        t_1643[k] = f_12 * sli0_1031[k]
                    - f_13 * sli1_1031[k]
                    + f_3 * pc_x[k] * slk_1319[k];

        t_1644[k] = f_12 * sli0_1032[k]
                    - f_13 * sli1_1032[k]
                    + f_3 * pc_x[k] * slk_1320[k];
    }

#pragma omp simd aligned(t_1645, t_1646, t_1647, t_1648, pc_x, pc_y, skk_1028, sli0_1033, \
                         sli0_1035, sli1_1033, sli1_1035, slk_1316, slk_1321, slk_1323, \
                         slk_1324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1645[k] = f_12 * sli0_1033[k]
                    - f_13 * sli1_1033[k]
                    + f_3 * pc_x[k] * slk_1321[k];

        t_1646[k] = f_0 * skk_1028[k]
                    + f_3 * pc_y[k] * slk_1316[k];

        t_1647[k] = f_12 * sli0_1035[k]
                    - f_13 * sli1_1035[k]
                    + f_3 * pc_x[k] * slk_1323[k];

        t_1648[k] = f_3 * pc_x[k] * slk_1324[k];
    }

#pragma omp simd aligned(t_1649, t_1650, t_1651, t_1652, t_1653, t_1654, t_1655, pc_x, \
                         slk_1325, slk_1326, slk_1327, slk_1328, slk_1329, slk_1330, \
                         slk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1649[k] = f_3 * pc_x[k] * slk_1325[k];

        t_1650[k] = f_3 * pc_x[k] * slk_1326[k];

        t_1651[k] = f_3 * pc_x[k] * slk_1327[k];

        t_1652[k] = f_3 * pc_x[k] * slk_1328[k];

        t_1653[k] = f_3 * pc_x[k] * slk_1329[k];

        t_1654[k] = f_3 * pc_x[k] * slk_1330[k];

        t_1655[k] = f_3 * pc_x[k] * slk_1331[k];
    }

#pragma omp simd aligned(t_1656, t_1657, t_1658, pc_y, pc_z, skk_1036, skk_1038, sli0_1029, \
                         sli0_1031, sli1_1029, sli1_1031, slk_1324, \
                         slk_1326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1656[k] = f_0 * skk_1036[k]
                    + f_1 * sli0_1029[k]
                    - f_2 * sli1_1029[k]
                    + f_3 * pc_y[k] * slk_1324[k];

        t_1657[k] = f_3 * pc_z[k] * slk_1324[k];

        t_1658[k] = f_0 * skk_1038[k]
                    + f_4 * sli0_1031[k]
                    - f_5 * sli1_1031[k]
                    + f_3 * pc_y[k] * slk_1326[k];
    }

#pragma omp simd aligned(t_1659, t_1660, t_1661, pc_y, skk_1039, skk_1040, skk_1041, \
                         sli0_1032, sli0_1033, sli0_1034, sli1_1032, sli1_1033, sli1_1034, \
                         slk_1327, slk_1328, slk_1329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1659[k] = f_0 * skk_1039[k]
                    + f_6 * sli0_1032[k]
                    - f_7 * sli1_1032[k]
                    + f_3 * pc_y[k] * slk_1327[k];

        t_1660[k] = f_0 * skk_1040[k]
                    + f_8 * sli0_1033[k]
                    - f_9 * sli1_1033[k]
                    + f_3 * pc_y[k] * slk_1328[k];

        t_1661[k] = f_0 * skk_1041[k]
                    + f_10 * sli0_1034[k]
                    - f_11 * sli1_1034[k]
                    + f_3 * pc_y[k] * slk_1329[k];
    }

#pragma omp simd aligned(t_1662, t_1663, t_1664, t_1665, pb_z, pc_y, pc_z, skl0_1260, \
                         skk_1042, skk_1043, skl1_1260, sli0_1035, sli1_1035, slk_1330, \
                         slk_1331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1662[k] = f_0 * skk_1042[k]
                    + f_12 * sli0_1035[k]
                    - f_13 * sli1_1035[k]
                    + f_3 * pc_y[k] * slk_1330[k];

        t_1663[k] = f_0 * skk_1043[k]
                    + f_3 * pc_y[k] * slk_1331[k];

        t_1664[k] = f_1 * sli0_1035[k]
                    - f_2 * sli1_1035[k]
                    + f_3 * pc_z[k] * slk_1331[k];

        t_1665[k] = pb_z[k] * skl0_1260[k]
                    - f_14 * pc_z[k] * skl1_1260[k];
    }

#pragma omp simd aligned(t_1666, t_1667, t_1668, t_1669, pb_z, pc_y, pc_z, skl0_1263, \
                         skk_1008, skk_1044, skk_1046, skl1_1263, slk_1332, \
                         slk_1334 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1666[k] = f_21 * skk_1044[k]
                    + f_3 * pc_y[k] * slk_1332[k];

        t_1667[k] = f_15 * skk_1008[k]
                    + f_3 * pc_z[k] * slk_1332[k];

        t_1668[k] = pb_z[k] * skl0_1263[k]
                    - f_14 * pc_z[k] * skl1_1263[k];

        t_1669[k] = f_21 * skk_1046[k]
                    + f_3 * pc_y[k] * slk_1334[k];
    }

#pragma omp simd aligned(t_1670, t_1671, t_1672, t_1673, pb_z, pc_x, pc_y, pc_z, skl0_1266, \
                         skk_1011, skk_1049, skl1_1266, sli0_1041, sli1_1041, slk_1335, \
                         slk_1337 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1670[k] = f_4 * sli0_1041[k]
                    - f_5 * sli1_1041[k]
                    + f_3 * pc_x[k] * slk_1337[k];

        t_1671[k] = pb_z[k] * skl0_1266[k]
                    - f_14 * pc_z[k] * skl1_1266[k];

        t_1672[k] = f_15 * skk_1011[k]
                    + f_3 * pc_z[k] * slk_1335[k];

        t_1673[k] = f_21 * skk_1049[k]
                    + f_3 * pc_y[k] * slk_1337[k];
    }

#pragma omp simd aligned(t_1674, t_1675, t_1676, pb_z, pc_x, pc_z, skl0_1270, skk_1014, \
                         skl1_1270, sli0_1045, sli1_1045, slk_1338, \
                         slk_1341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1674[k] = f_6 * sli0_1045[k]
                    - f_7 * sli1_1045[k]
                    + f_3 * pc_x[k] * slk_1341[k];

        t_1675[k] = pb_z[k] * skl0_1270[k]
                    - f_14 * pc_z[k] * skl1_1270[k];

        t_1676[k] = f_15 * skk_1014[k]
                    + f_3 * pc_z[k] * slk_1338[k];
    }

#pragma omp simd aligned(t_1677, t_1678, t_1679, pc_x, pc_y, skk_1053, sli0_1048, sli0_1050, \
                         sli1_1048, sli1_1050, slk_1341, slk_1344, \
                         slk_1346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1677[k] = f_8 * sli0_1048[k]
                    - f_9 * sli1_1048[k]
                    + f_3 * pc_x[k] * slk_1344[k];

        t_1678[k] = f_21 * skk_1053[k]
                    + f_3 * pc_y[k] * slk_1341[k];

        t_1679[k] = f_8 * sli0_1050[k]
                    - f_9 * sli1_1050[k]
                    + f_3 * pc_x[k] * slk_1346[k];
    }

#pragma omp simd aligned(t_1680, t_1681, t_1682, pb_z, pc_x, pc_z, skl0_1275, skk_1018, \
                         skl1_1275, sli0_1053, sli1_1053, slk_1342, \
                         slk_1349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1680[k] = pb_z[k] * skl0_1275[k]
                    - f_14 * pc_z[k] * skl1_1275[k];

        t_1681[k] = f_15 * skk_1018[k]
                    + f_3 * pc_z[k] * slk_1342[k];

        t_1682[k] = f_10 * sli0_1053[k]
                    - f_11 * sli1_1053[k]
                    + f_3 * pc_x[k] * slk_1349[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece15(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t sli0, const size_t sli1,
                                                           const size_t slk, const size_t ncols,
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
    const auto f_21 = 3.5 / q;

    auto *t_1683 = buffer.data(target + 1683);
    auto *t_1684 = buffer.data(target + 1684);
    auto *t_1685 = buffer.data(target + 1685);
    auto *t_1686 = buffer.data(target + 1686);
    auto *t_1687 = buffer.data(target + 1687);
    auto *t_1688 = buffer.data(target + 1688);
    auto *t_1689 = buffer.data(target + 1689);
    auto *t_1690 = buffer.data(target + 1690);
    auto *t_1691 = buffer.data(target + 1691);
    auto *t_1692 = buffer.data(target + 1692);
    auto *t_1693 = buffer.data(target + 1693);
    auto *t_1694 = buffer.data(target + 1694);
    auto *t_1695 = buffer.data(target + 1695);
    auto *t_1696 = buffer.data(target + 1696);
    auto *t_1697 = buffer.data(target + 1697);
    auto *t_1698 = buffer.data(target + 1698);
    auto *t_1699 = buffer.data(target + 1699);
    auto *t_1700 = buffer.data(target + 1700);
    auto *t_1701 = buffer.data(target + 1701);
    auto *t_1702 = buffer.data(target + 1702);
    auto *t_1703 = buffer.data(target + 1703);
    auto *t_1704 = buffer.data(target + 1704);
    auto *t_1705 = buffer.data(target + 1705);
    auto *t_1706 = buffer.data(target + 1706);
    auto *t_1707 = buffer.data(target + 1707);
    auto *t_1708 = buffer.data(target + 1708);
    auto *t_1709 = buffer.data(target + 1709);
    auto *t_1710 = buffer.data(target + 1710);
    auto *t_1711 = buffer.data(target + 1711);
    auto *t_1712 = buffer.data(target + 1712);
    auto *t_1713 = buffer.data(target + 1713);
    auto *t_1714 = buffer.data(target + 1714);
    auto *t_1715 = buffer.data(target + 1715);
    auto *t_1716 = buffer.data(target + 1716);
    auto *t_1717 = buffer.data(target + 1717);
    auto *t_1718 = buffer.data(target + 1718);
    auto *t_1719 = buffer.data(target + 1719);
    auto *t_1720 = buffer.data(target + 1720);
    auto *t_1721 = buffer.data(target + 1721);
    auto *t_1722 = buffer.data(target + 1722);
    auto *t_1723 = buffer.data(target + 1723);
    auto *t_1724 = buffer.data(target + 1724);
    auto *t_1725 = buffer.data(target + 1725);
    auto *t_1726 = buffer.data(target + 1726);
    auto *t_1727 = buffer.data(target + 1727);
    auto *t_1728 = buffer.data(target + 1728);
    auto *t_1729 = buffer.data(target + 1729);
    auto *t_1730 = buffer.data(target + 1730);
    auto *t_1731 = buffer.data(target + 1731);
    auto *t_1732 = buffer.data(target + 1732);
    auto *t_1733 = buffer.data(target + 1733);
    auto *t_1734 = buffer.data(target + 1734);
    auto *t_1735 = buffer.data(target + 1735);
    auto *t_1736 = buffer.data(target + 1736);
    auto *t_1737 = buffer.data(target + 1737);
    auto *t_1738 = buffer.data(target + 1738);
    auto *t_1739 = buffer.data(target + 1739);
    auto *t_1740 = buffer.data(target + 1740);
    auto *t_1741 = buffer.data(target + 1741);
    auto *t_1742 = buffer.data(target + 1742);
    auto *t_1743 = buffer.data(target + 1743);
    auto *t_1744 = buffer.data(target + 1744);
    auto *t_1745 = buffer.data(target + 1745);
    auto *t_1746 = buffer.data(target + 1746);
    auto *t_1747 = buffer.data(target + 1747);
    auto *t_1748 = buffer.data(target + 1748);
    auto *t_1749 = buffer.data(target + 1749);
    auto *t_1750 = buffer.data(target + 1750);
    auto *t_1751 = buffer.data(target + 1751);
    auto *t_1752 = buffer.data(target + 1752);
    auto *t_1753 = buffer.data(target + 1753);
    auto *t_1754 = buffer.data(target + 1754);
    auto *t_1755 = buffer.data(target + 1755);
    auto *t_1756 = buffer.data(target + 1756);
    auto *t_1757 = buffer.data(target + 1757);
    auto *t_1758 = buffer.data(target + 1758);
    auto *t_1759 = buffer.data(target + 1759);
    auto *t_1760 = buffer.data(target + 1760);
    auto *t_1761 = buffer.data(target + 1761);
    auto *t_1762 = buffer.data(target + 1762);
    auto *t_1763 = buffer.data(target + 1763);
    auto *t_1764 = buffer.data(target + 1764);
    auto *t_1765 = buffer.data(target + 1765);
    auto *t_1766 = buffer.data(target + 1766);
    auto *t_1767 = buffer.data(target + 1767);
    auto *t_1768 = buffer.data(target + 1768);
    auto *t_1769 = buffer.data(target + 1769);
    auto *t_1770 = buffer.data(target + 1770);
    auto *t_1771 = buffer.data(target + 1771);
    auto *t_1772 = buffer.data(target + 1772);
    auto *t_1773 = buffer.data(target + 1773);
    auto *t_1774 = buffer.data(target + 1774);
    auto *t_1775 = buffer.data(target + 1775);
    auto *t_1776 = buffer.data(target + 1776);
    auto *t_1777 = buffer.data(target + 1777);
    auto *t_1778 = buffer.data(target + 1778);
    auto *t_1779 = buffer.data(target + 1779);
    auto *t_1780 = buffer.data(target + 1780);
    auto *t_1781 = buffer.data(target + 1781);
    auto *t_1782 = buffer.data(target + 1782);
    auto *t_1783 = buffer.data(target + 1783);
    auto *t_1784 = buffer.data(target + 1784);
    auto *t_1785 = buffer.data(target + 1785);
    auto *t_1786 = buffer.data(target + 1786);
    auto *t_1787 = buffer.data(target + 1787);
    auto *t_1788 = buffer.data(target + 1788);
    auto *t_1789 = buffer.data(target + 1789);
    auto *t_1790 = buffer.data(target + 1790);
    auto *t_1791 = buffer.data(target + 1791);
    auto *t_1792 = buffer.data(target + 1792);
    auto *t_1793 = buffer.data(target + 1793);
    auto *t_1794 = buffer.data(target + 1794);
    auto *t_1795 = buffer.data(target + 1795);
    auto *t_1796 = buffer.data(target + 1796);
    auto *t_1797 = buffer.data(target + 1797);
    auto *t_1798 = buffer.data(target + 1798);
    auto *t_1799 = buffer.data(target + 1799);
    auto *t_1800 = buffer.data(target + 1800);
    auto *t_1801 = buffer.data(target + 1801);
    auto *t_1802 = buffer.data(target + 1802);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_1281 = buffer.data(skl0 + 1281);
    const auto *skl0_1296 = buffer.data(skl0 + 1296);
    const auto *skl0_1298 = buffer.data(skl0 + 1298);
    const auto *skl0_1299 = buffer.data(skl0 + 1299);
    const auto *skl0_1300 = buffer.data(skl0 + 1300);
    const auto *skl0_1301 = buffer.data(skl0 + 1301);
    const auto *skl0_1302 = buffer.data(skl0 + 1302);

    const auto *skk_1023 = buffer.data(skk + 1023);
    const auto *skk_1036 = buffer.data(skk + 1036);
    const auto *skk_1037 = buffer.data(skk + 1037);
    const auto *skk_1038 = buffer.data(skk + 1038);
    const auto *skk_1039 = buffer.data(skk + 1039);
    const auto *skk_1040 = buffer.data(skk + 1040);
    const auto *skk_1041 = buffer.data(skk + 1041);
    const auto *skk_1043 = buffer.data(skk + 1043);
    const auto *skk_1044 = buffer.data(skk + 1044);
    const auto *skk_1047 = buffer.data(skk + 1047);
    const auto *skk_1050 = buffer.data(skk + 1050);
    const auto *skk_1054 = buffer.data(skk + 1054);
    const auto *skk_1058 = buffer.data(skk + 1058);
    const auto *skk_1059 = buffer.data(skk + 1059);
    const auto *skk_1064 = buffer.data(skk + 1064);
    const auto *skk_1072 = buffer.data(skk + 1072);
    const auto *skk_1079 = buffer.data(skk + 1079);
    const auto *skk_1080 = buffer.data(skk + 1080);
    const auto *skk_1082 = buffer.data(skk + 1082);
    const auto *skk_1083 = buffer.data(skk + 1083);
    const auto *skk_1085 = buffer.data(skk + 1085);
    const auto *skk_1086 = buffer.data(skk + 1086);
    const auto *skk_1089 = buffer.data(skk + 1089);
    const auto *skk_1090 = buffer.data(skk + 1090);
    const auto *skk_1094 = buffer.data(skk + 1094);
    const auto *skk_1095 = buffer.data(skk + 1095);
    const auto *skk_1100 = buffer.data(skk + 1100);
    const auto *skk_1108 = buffer.data(skk + 1108);
    const auto *skk_1110 = buffer.data(skk + 1110);
    const auto *skk_1111 = buffer.data(skk + 1111);
    const auto *skk_1112 = buffer.data(skk + 1112);
    const auto *skk_1113 = buffer.data(skk + 1113);
    const auto *skk_1114 = buffer.data(skk + 1114);
    const auto *skk_1115 = buffer.data(skk + 1115);
    const auto *skk_1116 = buffer.data(skk + 1116);
    const auto *skk_1118 = buffer.data(skk + 1118);
    const auto *skk_1121 = buffer.data(skk + 1121);
    const auto *skk_1125 = buffer.data(skk + 1125);
    const auto *skk_1130 = buffer.data(skk + 1130);
    const auto *skk_1136 = buffer.data(skk + 1136);
    const auto *skk_1144 = buffer.data(skk + 1144);
    const auto *skk_1146 = buffer.data(skk + 1146);
    const auto *skk_1147 = buffer.data(skk + 1147);
    const auto *skk_1148 = buffer.data(skk + 1148);
    const auto *skk_1149 = buffer.data(skk + 1149);
    const auto *skk_1150 = buffer.data(skk + 1150);
    const auto *skk_1151 = buffer.data(skk + 1151);
    const auto *skk_1152 = buffer.data(skk + 1152);

    const auto *skl1_1281 = buffer.data(skl1 + 1281);
    const auto *skl1_1296 = buffer.data(skl1 + 1296);
    const auto *skl1_1298 = buffer.data(skl1 + 1298);
    const auto *skl1_1299 = buffer.data(skl1 + 1299);
    const auto *skl1_1300 = buffer.data(skl1 + 1300);
    const auto *skl1_1301 = buffer.data(skl1 + 1301);
    const auto *skl1_1302 = buffer.data(skl1 + 1302);

    const auto *sli0_1054 = buffer.data(sli0 + 1054);
    const auto *sli0_1056 = buffer.data(sli0 + 1056);
    const auto *sli0_1059 = buffer.data(sli0 + 1059);
    const auto *sli0_1060 = buffer.data(sli0 + 1060);
    const auto *sli0_1061 = buffer.data(sli0 + 1061);
    const auto *sli0_1063 = buffer.data(sli0 + 1063);
    const auto *sli0_1064 = buffer.data(sli0 + 1064);
    const auto *sli0_1067 = buffer.data(sli0 + 1067);
    const auto *sli0_1069 = buffer.data(sli0 + 1069);
    const auto *sli0_1070 = buffer.data(sli0 + 1070);
    const auto *sli0_1073 = buffer.data(sli0 + 1073);
    const auto *sli0_1074 = buffer.data(sli0 + 1074);
    const auto *sli0_1076 = buffer.data(sli0 + 1076);
    const auto *sli0_1078 = buffer.data(sli0 + 1078);
    const auto *sli0_1079 = buffer.data(sli0 + 1079);
    const auto *sli0_1081 = buffer.data(sli0 + 1081);
    const auto *sli0_1082 = buffer.data(sli0 + 1082);
    const auto *sli0_1084 = buffer.data(sli0 + 1084);
    const auto *sli0_1085 = buffer.data(sli0 + 1085);
    const auto *sli0_1087 = buffer.data(sli0 + 1087);
    const auto *sli0_1088 = buffer.data(sli0 + 1088);
    const auto *sli0_1089 = buffer.data(sli0 + 1089);
    const auto *sli0_1090 = buffer.data(sli0 + 1090);
    const auto *sli0_1091 = buffer.data(sli0 + 1091);
    const auto *sli0_1092 = buffer.data(sli0 + 1092);
    const auto *sli0_1095 = buffer.data(sli0 + 1095);
    const auto *sli0_1097 = buffer.data(sli0 + 1097);
    const auto *sli0_1098 = buffer.data(sli0 + 1098);
    const auto *sli0_1101 = buffer.data(sli0 + 1101);
    const auto *sli0_1102 = buffer.data(sli0 + 1102);
    const auto *sli0_1104 = buffer.data(sli0 + 1104);
    const auto *sli0_1106 = buffer.data(sli0 + 1106);
    const auto *sli0_1107 = buffer.data(sli0 + 1107);
    const auto *sli0_1109 = buffer.data(sli0 + 1109);
    const auto *sli0_1110 = buffer.data(sli0 + 1110);
    const auto *sli0_1112 = buffer.data(sli0 + 1112);
    const auto *sli0_1113 = buffer.data(sli0 + 1113);
    const auto *sli0_1115 = buffer.data(sli0 + 1115);
    const auto *sli0_1116 = buffer.data(sli0 + 1116);
    const auto *sli0_1117 = buffer.data(sli0 + 1117);
    const auto *sli0_1118 = buffer.data(sli0 + 1118);
    const auto *sli0_1119 = buffer.data(sli0 + 1119);
    const auto *sli0_1120 = buffer.data(sli0 + 1120);

    const auto *sli1_1054 = buffer.data(sli1 + 1054);
    const auto *sli1_1056 = buffer.data(sli1 + 1056);
    const auto *sli1_1059 = buffer.data(sli1 + 1059);
    const auto *sli1_1060 = buffer.data(sli1 + 1060);
    const auto *sli1_1061 = buffer.data(sli1 + 1061);
    const auto *sli1_1063 = buffer.data(sli1 + 1063);
    const auto *sli1_1064 = buffer.data(sli1 + 1064);
    const auto *sli1_1067 = buffer.data(sli1 + 1067);
    const auto *sli1_1069 = buffer.data(sli1 + 1069);
    const auto *sli1_1070 = buffer.data(sli1 + 1070);
    const auto *sli1_1073 = buffer.data(sli1 + 1073);
    const auto *sli1_1074 = buffer.data(sli1 + 1074);
    const auto *sli1_1076 = buffer.data(sli1 + 1076);
    const auto *sli1_1078 = buffer.data(sli1 + 1078);
    const auto *sli1_1079 = buffer.data(sli1 + 1079);
    const auto *sli1_1081 = buffer.data(sli1 + 1081);
    const auto *sli1_1082 = buffer.data(sli1 + 1082);
    const auto *sli1_1084 = buffer.data(sli1 + 1084);
    const auto *sli1_1085 = buffer.data(sli1 + 1085);
    const auto *sli1_1087 = buffer.data(sli1 + 1087);
    const auto *sli1_1088 = buffer.data(sli1 + 1088);
    const auto *sli1_1089 = buffer.data(sli1 + 1089);
    const auto *sli1_1090 = buffer.data(sli1 + 1090);
    const auto *sli1_1091 = buffer.data(sli1 + 1091);
    const auto *sli1_1092 = buffer.data(sli1 + 1092);
    const auto *sli1_1095 = buffer.data(sli1 + 1095);
    const auto *sli1_1097 = buffer.data(sli1 + 1097);
    const auto *sli1_1098 = buffer.data(sli1 + 1098);
    const auto *sli1_1101 = buffer.data(sli1 + 1101);
    const auto *sli1_1102 = buffer.data(sli1 + 1102);
    const auto *sli1_1104 = buffer.data(sli1 + 1104);
    const auto *sli1_1106 = buffer.data(sli1 + 1106);
    const auto *sli1_1107 = buffer.data(sli1 + 1107);
    const auto *sli1_1109 = buffer.data(sli1 + 1109);
    const auto *sli1_1110 = buffer.data(sli1 + 1110);
    const auto *sli1_1112 = buffer.data(sli1 + 1112);
    const auto *sli1_1113 = buffer.data(sli1 + 1113);
    const auto *sli1_1115 = buffer.data(sli1 + 1115);
    const auto *sli1_1116 = buffer.data(sli1 + 1116);
    const auto *sli1_1117 = buffer.data(sli1 + 1117);
    const auto *sli1_1118 = buffer.data(sli1 + 1118);
    const auto *sli1_1119 = buffer.data(sli1 + 1119);
    const auto *sli1_1120 = buffer.data(sli1 + 1120);

    const auto *slk_1346 = buffer.data(slk + 1346);
    const auto *slk_1347 = buffer.data(slk + 1347);
    const auto *slk_1350 = buffer.data(slk + 1350);
    const auto *slk_1352 = buffer.data(slk + 1352);
    const auto *slk_1355 = buffer.data(slk + 1355);
    const auto *slk_1356 = buffer.data(slk + 1356);
    const auto *slk_1357 = buffer.data(slk + 1357);
    const auto *slk_1359 = buffer.data(slk + 1359);
    const auto *slk_1360 = buffer.data(slk + 1360);
    const auto *slk_1361 = buffer.data(slk + 1361);
    const auto *slk_1362 = buffer.data(slk + 1362);
    const auto *slk_1363 = buffer.data(slk + 1363);
    const auto *slk_1364 = buffer.data(slk + 1364);
    const auto *slk_1365 = buffer.data(slk + 1365);
    const auto *slk_1366 = buffer.data(slk + 1366);
    const auto *slk_1367 = buffer.data(slk + 1367);
    const auto *slk_1368 = buffer.data(slk + 1368);
    const auto *slk_1370 = buffer.data(slk + 1370);
    const auto *slk_1371 = buffer.data(slk + 1371);
    const auto *slk_1373 = buffer.data(slk + 1373);
    const auto *slk_1374 = buffer.data(slk + 1374);
    const auto *slk_1377 = buffer.data(slk + 1377);
    const auto *slk_1378 = buffer.data(slk + 1378);
    const auto *slk_1380 = buffer.data(slk + 1380);
    const auto *slk_1382 = buffer.data(slk + 1382);
    const auto *slk_1383 = buffer.data(slk + 1383);
    const auto *slk_1385 = buffer.data(slk + 1385);
    const auto *slk_1386 = buffer.data(slk + 1386);
    const auto *slk_1388 = buffer.data(slk + 1388);
    const auto *slk_1389 = buffer.data(slk + 1389);
    const auto *slk_1391 = buffer.data(slk + 1391);
    const auto *slk_1392 = buffer.data(slk + 1392);
    const auto *slk_1393 = buffer.data(slk + 1393);
    const auto *slk_1395 = buffer.data(slk + 1395);
    const auto *slk_1396 = buffer.data(slk + 1396);
    const auto *slk_1397 = buffer.data(slk + 1397);
    const auto *slk_1398 = buffer.data(slk + 1398);
    const auto *slk_1399 = buffer.data(slk + 1399);
    const auto *slk_1400 = buffer.data(slk + 1400);
    const auto *slk_1401 = buffer.data(slk + 1401);
    const auto *slk_1402 = buffer.data(slk + 1402);
    const auto *slk_1403 = buffer.data(slk + 1403);
    const auto *slk_1404 = buffer.data(slk + 1404);
    const auto *slk_1406 = buffer.data(slk + 1406);
    const auto *slk_1407 = buffer.data(slk + 1407);
    const auto *slk_1409 = buffer.data(slk + 1409);
    const auto *slk_1410 = buffer.data(slk + 1410);
    const auto *slk_1413 = buffer.data(slk + 1413);
    const auto *slk_1414 = buffer.data(slk + 1414);
    const auto *slk_1416 = buffer.data(slk + 1416);
    const auto *slk_1418 = buffer.data(slk + 1418);
    const auto *slk_1419 = buffer.data(slk + 1419);
    const auto *slk_1421 = buffer.data(slk + 1421);
    const auto *slk_1422 = buffer.data(slk + 1422);
    const auto *slk_1424 = buffer.data(slk + 1424);
    const auto *slk_1425 = buffer.data(slk + 1425);
    const auto *slk_1427 = buffer.data(slk + 1427);
    const auto *slk_1428 = buffer.data(slk + 1428);
    const auto *slk_1429 = buffer.data(slk + 1429);
    const auto *slk_1431 = buffer.data(slk + 1431);
    const auto *slk_1432 = buffer.data(slk + 1432);
    const auto *slk_1433 = buffer.data(slk + 1433);
    const auto *slk_1434 = buffer.data(slk + 1434);
    const auto *slk_1435 = buffer.data(slk + 1435);
    const auto *slk_1436 = buffer.data(slk + 1436);
    const auto *slk_1437 = buffer.data(slk + 1437);
    const auto *slk_1438 = buffer.data(slk + 1438);
    const auto *slk_1439 = buffer.data(slk + 1439);
    const auto *slk_1440 = buffer.data(slk + 1440);

#pragma omp simd aligned(t_1683, t_1684, t_1685, pc_x, pc_y, skk_1058, sli0_1054, sli0_1056, \
                         sli1_1054, sli1_1056, slk_1346, slk_1350, \
                         slk_1352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1683[k] = f_10 * sli0_1054[k]
                    - f_11 * sli1_1054[k]
                    + f_3 * pc_x[k] * slk_1350[k];

        t_1684[k] = f_21 * skk_1058[k]
                    + f_3 * pc_y[k] * slk_1346[k];

        t_1685[k] = f_10 * sli0_1056[k]
                    - f_11 * sli1_1056[k]
                    + f_3 * pc_x[k] * slk_1352[k];
    }

#pragma omp simd aligned(t_1686, t_1687, t_1688, pb_z, pc_x, pc_z, skl0_1281, skk_1023, \
                         skl1_1281, sli0_1059, sli1_1059, slk_1347, \
                         slk_1355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1686[k] = pb_z[k] * skl0_1281[k]
                    - f_14 * pc_z[k] * skl1_1281[k];

        t_1687[k] = f_15 * skk_1023[k]
                    + f_3 * pc_z[k] * slk_1347[k];

        t_1688[k] = f_12 * sli0_1059[k]
                    - f_13 * sli1_1059[k]
                    + f_3 * pc_x[k] * slk_1355[k];
    }

#pragma omp simd aligned(t_1689, t_1690, t_1691, pc_x, pc_y, skk_1064, sli0_1060, sli0_1061, \
                         sli1_1060, sli1_1061, slk_1352, slk_1356, \
                         slk_1357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1689[k] = f_12 * sli0_1060[k]
                    - f_13 * sli1_1060[k]
                    + f_3 * pc_x[k] * slk_1356[k];

        t_1690[k] = f_12 * sli0_1061[k]
                    - f_13 * sli1_1061[k]
                    + f_3 * pc_x[k] * slk_1357[k];

        t_1691[k] = f_21 * skk_1064[k]
                    + f_3 * pc_y[k] * slk_1352[k];
    }

#pragma omp simd aligned(t_1692, t_1693, t_1694, t_1695, t_1696, t_1697, pc_x, sli0_1063, \
                         sli1_1063, slk_1359, slk_1360, slk_1361, slk_1362, slk_1363, \
                         slk_1364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1692[k] = f_12 * sli0_1063[k]
                    - f_13 * sli1_1063[k]
                    + f_3 * pc_x[k] * slk_1359[k];

        t_1693[k] = f_3 * pc_x[k] * slk_1360[k];

        t_1694[k] = f_3 * pc_x[k] * slk_1361[k];

        t_1695[k] = f_3 * pc_x[k] * slk_1362[k];

        t_1696[k] = f_3 * pc_x[k] * slk_1363[k];

        t_1697[k] = f_3 * pc_x[k] * slk_1364[k];
    }

#pragma omp simd aligned(t_1698, t_1699, t_1700, t_1701, t_1702, pb_z, pc_x, pc_z, skl0_1296, \
                         skk_1036, skl1_1296, slk_1360, slk_1365, slk_1366, \
                         slk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1698[k] = f_3 * pc_x[k] * slk_1365[k];

        t_1699[k] = f_3 * pc_x[k] * slk_1366[k];

        t_1700[k] = f_3 * pc_x[k] * slk_1367[k];

        t_1701[k] = pb_z[k] * skl0_1296[k]
                    - f_14 * pc_z[k] * skl1_1296[k];

        t_1702[k] = f_15 * skk_1036[k]
                    + f_3 * pc_z[k] * slk_1360[k];
    }

#pragma omp simd aligned(t_1703, t_1704, t_1705, pb_z, pc_z, skl0_1298, skl0_1299, skl0_1300, \
                         skk_1037, skk_1038, skk_1039, skl1_1298, skl1_1299, \
                         skl1_1300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1703[k] = pb_z[k] * skl0_1298[k]
                    + f_16 * skk_1037[k]
                    - f_14 * pc_z[k] * skl1_1298[k];

        t_1704[k] = pb_z[k] * skl0_1299[k]
                    + f_17 * skk_1038[k]
                    - f_14 * pc_z[k] * skl1_1299[k];

        t_1705[k] = pb_z[k] * skl0_1300[k]
                    + f_18 * skk_1039[k]
                    - f_14 * pc_z[k] * skl1_1300[k];
    }

#pragma omp simd aligned(t_1706, t_1707, t_1708, pb_z, pc_y, pc_z, skl0_1301, skl0_1302, \
                         skk_1040, skk_1041, skk_1079, skl1_1301, skl1_1302, \
                         slk_1367 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1706[k] = pb_z[k] * skl0_1301[k]
                    + f_19 * skk_1040[k]
                    - f_14 * pc_z[k] * skl1_1301[k];

        t_1707[k] = pb_z[k] * skl0_1302[k]
                    + f_20 * skk_1041[k]
                    - f_14 * pc_z[k] * skl1_1302[k];

        t_1708[k] = f_21 * skk_1079[k]
                    + f_3 * pc_y[k] * slk_1367[k];
    }

#pragma omp simd aligned(t_1709, t_1710, t_1711, t_1712, pc_x, pc_y, pc_z, skk_1043, skk_1044, \
                         skk_1080, sli0_1063, sli0_1064, sli1_1063, sli1_1064, slk_1367, \
                         slk_1368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1709[k] = f_15 * skk_1043[k]
                    + f_1 * sli0_1063[k]
                    - f_2 * sli1_1063[k]
                    + f_3 * pc_z[k] * slk_1367[k];

        t_1710[k] = f_1 * sli0_1064[k]
                    - f_2 * sli1_1064[k]
                    + f_3 * pc_x[k] * slk_1368[k];

        t_1711[k] = f_20 * skk_1080[k]
                    + f_3 * pc_y[k] * slk_1368[k];

        t_1712[k] = f_16 * skk_1044[k]
                    + f_3 * pc_z[k] * slk_1368[k];
    }

#pragma omp simd aligned(t_1713, t_1714, t_1715, pc_x, pc_y, skk_1082, sli0_1067, sli0_1069, \
                         sli1_1067, sli1_1069, slk_1370, slk_1371, \
                         slk_1373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1713[k] = f_4 * sli0_1067[k]
                    - f_5 * sli1_1067[k]
                    + f_3 * pc_x[k] * slk_1371[k];

        t_1714[k] = f_20 * skk_1082[k]
                    + f_3 * pc_y[k] * slk_1370[k];

        t_1715[k] = f_4 * sli0_1069[k]
                    - f_5 * sli1_1069[k]
                    + f_3 * pc_x[k] * slk_1373[k];
    }

#pragma omp simd aligned(t_1716, t_1717, t_1718, pc_x, pc_y, pc_z, skk_1047, skk_1085, \
                         sli0_1070, sli1_1070, slk_1371, slk_1373, \
                         slk_1374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1716[k] = f_6 * sli0_1070[k]
                    - f_7 * sli1_1070[k]
                    + f_3 * pc_x[k] * slk_1374[k];

        t_1717[k] = f_16 * skk_1047[k]
                    + f_3 * pc_z[k] * slk_1371[k];

        t_1718[k] = f_20 * skk_1085[k]
                    + f_3 * pc_y[k] * slk_1373[k];
    }

#pragma omp simd aligned(t_1719, t_1720, t_1721, pc_x, pc_z, skk_1050, sli0_1073, sli0_1074, \
                         sli1_1073, sli1_1074, slk_1374, slk_1377, \
                         slk_1378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1719[k] = f_6 * sli0_1073[k]
                    - f_7 * sli1_1073[k]
                    + f_3 * pc_x[k] * slk_1377[k];

        t_1720[k] = f_8 * sli0_1074[k]
                    - f_9 * sli1_1074[k]
                    + f_3 * pc_x[k] * slk_1378[k];

        t_1721[k] = f_16 * skk_1050[k]
                    + f_3 * pc_z[k] * slk_1374[k];
    }

#pragma omp simd aligned(t_1722, t_1723, t_1724, pc_x, pc_y, skk_1089, sli0_1076, sli0_1078, \
                         sli1_1076, sli1_1078, slk_1377, slk_1380, \
                         slk_1382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1722[k] = f_8 * sli0_1076[k]
                    - f_9 * sli1_1076[k]
                    + f_3 * pc_x[k] * slk_1380[k];

        t_1723[k] = f_20 * skk_1089[k]
                    + f_3 * pc_y[k] * slk_1377[k];

        t_1724[k] = f_8 * sli0_1078[k]
                    - f_9 * sli1_1078[k]
                    + f_3 * pc_x[k] * slk_1382[k];
    }

#pragma omp simd aligned(t_1725, t_1726, t_1727, pc_x, pc_z, skk_1054, sli0_1079, sli0_1081, \
                         sli1_1079, sli1_1081, slk_1378, slk_1383, \
                         slk_1385 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1725[k] = f_10 * sli0_1079[k]
                    - f_11 * sli1_1079[k]
                    + f_3 * pc_x[k] * slk_1383[k];

        t_1726[k] = f_16 * skk_1054[k]
                    + f_3 * pc_z[k] * slk_1378[k];

        t_1727[k] = f_10 * sli0_1081[k]
                    - f_11 * sli1_1081[k]
                    + f_3 * pc_x[k] * slk_1385[k];
    }

#pragma omp simd aligned(t_1728, t_1729, t_1730, pc_x, pc_y, skk_1094, sli0_1082, sli0_1084, \
                         sli1_1082, sli1_1084, slk_1382, slk_1386, \
                         slk_1388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1728[k] = f_10 * sli0_1082[k]
                    - f_11 * sli1_1082[k]
                    + f_3 * pc_x[k] * slk_1386[k];

        t_1729[k] = f_20 * skk_1094[k]
                    + f_3 * pc_y[k] * slk_1382[k];

        t_1730[k] = f_10 * sli0_1084[k]
                    - f_11 * sli1_1084[k]
                    + f_3 * pc_x[k] * slk_1388[k];
    }

#pragma omp simd aligned(t_1731, t_1732, t_1733, pc_x, pc_z, skk_1059, sli0_1085, sli0_1087, \
                         sli1_1085, sli1_1087, slk_1383, slk_1389, \
                         slk_1391 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1731[k] = f_12 * sli0_1085[k]
                    - f_13 * sli1_1085[k]
                    + f_3 * pc_x[k] * slk_1389[k];

        t_1732[k] = f_16 * skk_1059[k]
                    + f_3 * pc_z[k] * slk_1383[k];

        t_1733[k] = f_12 * sli0_1087[k]
                    - f_13 * sli1_1087[k]
                    + f_3 * pc_x[k] * slk_1391[k];
    }

#pragma omp simd aligned(t_1734, t_1735, t_1736, pc_x, pc_y, skk_1100, sli0_1088, sli0_1089, \
                         sli1_1088, sli1_1089, slk_1388, slk_1392, \
                         slk_1393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1734[k] = f_12 * sli0_1088[k]
                    - f_13 * sli1_1088[k]
                    + f_3 * pc_x[k] * slk_1392[k];

        t_1735[k] = f_12 * sli0_1089[k]
                    - f_13 * sli1_1089[k]
                    + f_3 * pc_x[k] * slk_1393[k];

        t_1736[k] = f_20 * skk_1100[k]
                    + f_3 * pc_y[k] * slk_1388[k];
    }

#pragma omp simd aligned(t_1737, t_1738, t_1739, t_1740, t_1741, t_1742, pc_x, sli0_1091, \
                         sli1_1091, slk_1395, slk_1396, slk_1397, slk_1398, slk_1399, \
                         slk_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1737[k] = f_12 * sli0_1091[k]
                    - f_13 * sli1_1091[k]
                    + f_3 * pc_x[k] * slk_1395[k];

        t_1738[k] = f_3 * pc_x[k] * slk_1396[k];

        t_1739[k] = f_3 * pc_x[k] * slk_1397[k];

        t_1740[k] = f_3 * pc_x[k] * slk_1398[k];

        t_1741[k] = f_3 * pc_x[k] * slk_1399[k];

        t_1742[k] = f_3 * pc_x[k] * slk_1400[k];
    }

#pragma omp simd aligned(t_1743, t_1744, t_1745, t_1746, t_1747, pc_x, pc_y, pc_z, skk_1072, \
                         skk_1108, sli0_1085, sli1_1085, slk_1396, slk_1401, slk_1402, \
                         slk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1743[k] = f_3 * pc_x[k] * slk_1401[k];

        t_1744[k] = f_3 * pc_x[k] * slk_1402[k];

        t_1745[k] = f_3 * pc_x[k] * slk_1403[k];

        t_1746[k] = f_20 * skk_1108[k]
                    + f_1 * sli0_1085[k]
                    - f_2 * sli1_1085[k]
                    + f_3 * pc_y[k] * slk_1396[k];

        t_1747[k] = f_16 * skk_1072[k]
                    + f_3 * pc_z[k] * slk_1396[k];
    }

#pragma omp simd aligned(t_1748, t_1749, t_1750, pc_y, skk_1110, skk_1111, skk_1112, \
                         sli0_1087, sli0_1088, sli0_1089, sli1_1087, sli1_1088, sli1_1089, \
                         slk_1398, slk_1399, slk_1400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1748[k] = f_20 * skk_1110[k]
                    + f_4 * sli0_1087[k]
                    - f_5 * sli1_1087[k]
                    + f_3 * pc_y[k] * slk_1398[k];

        t_1749[k] = f_20 * skk_1111[k]
                    + f_6 * sli0_1088[k]
                    - f_7 * sli1_1088[k]
                    + f_3 * pc_y[k] * slk_1399[k];

        t_1750[k] = f_20 * skk_1112[k]
                    + f_8 * sli0_1089[k]
                    - f_9 * sli1_1089[k]
                    + f_3 * pc_y[k] * slk_1400[k];
    }

#pragma omp simd aligned(t_1751, t_1752, t_1753, pc_y, skk_1113, skk_1114, skk_1115, \
                         sli0_1090, sli0_1091, sli1_1090, sli1_1091, slk_1401, slk_1402, \
                         slk_1403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1751[k] = f_20 * skk_1113[k]
                    + f_10 * sli0_1090[k]
                    - f_11 * sli1_1090[k]
                    + f_3 * pc_y[k] * slk_1401[k];

        t_1752[k] = f_20 * skk_1114[k]
                    + f_12 * sli0_1091[k]
                    - f_13 * sli1_1091[k]
                    + f_3 * pc_y[k] * slk_1402[k];

        t_1753[k] = f_20 * skk_1115[k]
                    + f_3 * pc_y[k] * slk_1403[k];
    }

#pragma omp simd aligned(t_1754, t_1755, t_1756, t_1757, pc_x, pc_y, pc_z, skk_1079, skk_1080, \
                         skk_1116, sli0_1091, sli0_1092, sli1_1091, sli1_1092, slk_1403, \
                         slk_1404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1754[k] = f_16 * skk_1079[k]
                    + f_1 * sli0_1091[k]
                    - f_2 * sli1_1091[k]
                    + f_3 * pc_z[k] * slk_1403[k];

        t_1755[k] = f_1 * sli0_1092[k]
                    - f_2 * sli1_1092[k]
                    + f_3 * pc_x[k] * slk_1404[k];

        t_1756[k] = f_19 * skk_1116[k]
                    + f_3 * pc_y[k] * slk_1404[k];

        t_1757[k] = f_17 * skk_1080[k]
                    + f_3 * pc_z[k] * slk_1404[k];
    }

#pragma omp simd aligned(t_1758, t_1759, t_1760, pc_x, pc_y, skk_1118, sli0_1095, sli0_1097, \
                         sli1_1095, sli1_1097, slk_1406, slk_1407, \
                         slk_1409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1758[k] = f_4 * sli0_1095[k]
                    - f_5 * sli1_1095[k]
                    + f_3 * pc_x[k] * slk_1407[k];

        t_1759[k] = f_19 * skk_1118[k]
                    + f_3 * pc_y[k] * slk_1406[k];

        t_1760[k] = f_4 * sli0_1097[k]
                    - f_5 * sli1_1097[k]
                    + f_3 * pc_x[k] * slk_1409[k];
    }

#pragma omp simd aligned(t_1761, t_1762, t_1763, pc_x, pc_y, pc_z, skk_1083, skk_1121, \
                         sli0_1098, sli1_1098, slk_1407, slk_1409, \
                         slk_1410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1761[k] = f_6 * sli0_1098[k]
                    - f_7 * sli1_1098[k]
                    + f_3 * pc_x[k] * slk_1410[k];

        t_1762[k] = f_17 * skk_1083[k]
                    + f_3 * pc_z[k] * slk_1407[k];

        t_1763[k] = f_19 * skk_1121[k]
                    + f_3 * pc_y[k] * slk_1409[k];
    }

#pragma omp simd aligned(t_1764, t_1765, t_1766, pc_x, pc_z, skk_1086, sli0_1101, sli0_1102, \
                         sli1_1101, sli1_1102, slk_1410, slk_1413, \
                         slk_1414 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1764[k] = f_6 * sli0_1101[k]
                    - f_7 * sli1_1101[k]
                    + f_3 * pc_x[k] * slk_1413[k];

        t_1765[k] = f_8 * sli0_1102[k]
                    - f_9 * sli1_1102[k]
                    + f_3 * pc_x[k] * slk_1414[k];

        t_1766[k] = f_17 * skk_1086[k]
                    + f_3 * pc_z[k] * slk_1410[k];
    }

#pragma omp simd aligned(t_1767, t_1768, t_1769, pc_x, pc_y, skk_1125, sli0_1104, sli0_1106, \
                         sli1_1104, sli1_1106, slk_1413, slk_1416, \
                         slk_1418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1767[k] = f_8 * sli0_1104[k]
                    - f_9 * sli1_1104[k]
                    + f_3 * pc_x[k] * slk_1416[k];

        t_1768[k] = f_19 * skk_1125[k]
                    + f_3 * pc_y[k] * slk_1413[k];

        t_1769[k] = f_8 * sli0_1106[k]
                    - f_9 * sli1_1106[k]
                    + f_3 * pc_x[k] * slk_1418[k];
    }

#pragma omp simd aligned(t_1770, t_1771, t_1772, pc_x, pc_z, skk_1090, sli0_1107, sli0_1109, \
                         sli1_1107, sli1_1109, slk_1414, slk_1419, \
                         slk_1421 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1770[k] = f_10 * sli0_1107[k]
                    - f_11 * sli1_1107[k]
                    + f_3 * pc_x[k] * slk_1419[k];

        t_1771[k] = f_17 * skk_1090[k]
                    + f_3 * pc_z[k] * slk_1414[k];

        t_1772[k] = f_10 * sli0_1109[k]
                    - f_11 * sli1_1109[k]
                    + f_3 * pc_x[k] * slk_1421[k];
    }

#pragma omp simd aligned(t_1773, t_1774, t_1775, pc_x, pc_y, skk_1130, sli0_1110, sli0_1112, \
                         sli1_1110, sli1_1112, slk_1418, slk_1422, \
                         slk_1424 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1773[k] = f_10 * sli0_1110[k]
                    - f_11 * sli1_1110[k]
                    + f_3 * pc_x[k] * slk_1422[k];

        t_1774[k] = f_19 * skk_1130[k]
                    + f_3 * pc_y[k] * slk_1418[k];

        t_1775[k] = f_10 * sli0_1112[k]
                    - f_11 * sli1_1112[k]
                    + f_3 * pc_x[k] * slk_1424[k];
    }

#pragma omp simd aligned(t_1776, t_1777, t_1778, pc_x, pc_z, skk_1095, sli0_1113, sli0_1115, \
                         sli1_1113, sli1_1115, slk_1419, slk_1425, \
                         slk_1427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1776[k] = f_12 * sli0_1113[k]
                    - f_13 * sli1_1113[k]
                    + f_3 * pc_x[k] * slk_1425[k];

        t_1777[k] = f_17 * skk_1095[k]
                    + f_3 * pc_z[k] * slk_1419[k];

        t_1778[k] = f_12 * sli0_1115[k]
                    - f_13 * sli1_1115[k]
                    + f_3 * pc_x[k] * slk_1427[k];
    }

#pragma omp simd aligned(t_1779, t_1780, t_1781, pc_x, pc_y, skk_1136, sli0_1116, sli0_1117, \
                         sli1_1116, sli1_1117, slk_1424, slk_1428, \
                         slk_1429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1779[k] = f_12 * sli0_1116[k]
                    - f_13 * sli1_1116[k]
                    + f_3 * pc_x[k] * slk_1428[k];

        t_1780[k] = f_12 * sli0_1117[k]
                    - f_13 * sli1_1117[k]
                    + f_3 * pc_x[k] * slk_1429[k];

        t_1781[k] = f_19 * skk_1136[k]
                    + f_3 * pc_y[k] * slk_1424[k];
    }

#pragma omp simd aligned(t_1782, t_1783, t_1784, t_1785, t_1786, t_1787, pc_x, sli0_1119, \
                         sli1_1119, slk_1431, slk_1432, slk_1433, slk_1434, slk_1435, \
                         slk_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1782[k] = f_12 * sli0_1119[k]
                    - f_13 * sli1_1119[k]
                    + f_3 * pc_x[k] * slk_1431[k];

        t_1783[k] = f_3 * pc_x[k] * slk_1432[k];

        t_1784[k] = f_3 * pc_x[k] * slk_1433[k];

        t_1785[k] = f_3 * pc_x[k] * slk_1434[k];

        t_1786[k] = f_3 * pc_x[k] * slk_1435[k];

        t_1787[k] = f_3 * pc_x[k] * slk_1436[k];
    }

#pragma omp simd aligned(t_1788, t_1789, t_1790, t_1791, t_1792, pc_x, pc_y, pc_z, skk_1108, \
                         skk_1144, sli0_1113, sli1_1113, slk_1432, slk_1437, slk_1438, \
                         slk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1788[k] = f_3 * pc_x[k] * slk_1437[k];

        t_1789[k] = f_3 * pc_x[k] * slk_1438[k];

        t_1790[k] = f_3 * pc_x[k] * slk_1439[k];

        t_1791[k] = f_19 * skk_1144[k]
                    + f_1 * sli0_1113[k]
                    - f_2 * sli1_1113[k]
                    + f_3 * pc_y[k] * slk_1432[k];

        t_1792[k] = f_17 * skk_1108[k]
                    + f_3 * pc_z[k] * slk_1432[k];
    }

#pragma omp simd aligned(t_1793, t_1794, t_1795, pc_y, skk_1146, skk_1147, skk_1148, \
                         sli0_1115, sli0_1116, sli0_1117, sli1_1115, sli1_1116, sli1_1117, \
                         slk_1434, slk_1435, slk_1436 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1793[k] = f_19 * skk_1146[k]
                    + f_4 * sli0_1115[k]
                    - f_5 * sli1_1115[k]
                    + f_3 * pc_y[k] * slk_1434[k];

        t_1794[k] = f_19 * skk_1147[k]
                    + f_6 * sli0_1116[k]
                    - f_7 * sli1_1116[k]
                    + f_3 * pc_y[k] * slk_1435[k];

        t_1795[k] = f_19 * skk_1148[k]
                    + f_8 * sli0_1117[k]
                    - f_9 * sli1_1117[k]
                    + f_3 * pc_y[k] * slk_1436[k];
    }

#pragma omp simd aligned(t_1796, t_1797, t_1798, pc_y, skk_1149, skk_1150, skk_1151, \
                         sli0_1118, sli0_1119, sli1_1118, sli1_1119, slk_1437, slk_1438, \
                         slk_1439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1796[k] = f_19 * skk_1149[k]
                    + f_10 * sli0_1118[k]
                    - f_11 * sli1_1118[k]
                    + f_3 * pc_y[k] * slk_1437[k];

        t_1797[k] = f_19 * skk_1150[k]
                    + f_12 * sli0_1119[k]
                    - f_13 * sli1_1119[k]
                    + f_3 * pc_y[k] * slk_1438[k];

        t_1798[k] = f_19 * skk_1151[k]
                    + f_3 * pc_y[k] * slk_1439[k];
    }

#pragma omp simd aligned(t_1799, t_1800, t_1801, t_1802, pc_x, pc_y, pc_z, skk_1115, skk_1116, \
                         skk_1152, sli0_1119, sli0_1120, sli1_1119, sli1_1120, slk_1439, \
                         slk_1440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1799[k] = f_17 * skk_1115[k]
                    + f_1 * sli0_1119[k]
                    - f_2 * sli1_1119[k]
                    + f_3 * pc_z[k] * slk_1439[k];

        t_1800[k] = f_1 * sli0_1120[k]
                    - f_2 * sli1_1120[k]
                    + f_3 * pc_x[k] * slk_1440[k];

        t_1801[k] = f_18 * skk_1152[k]
                    + f_3 * pc_y[k] * slk_1440[k];

        t_1802[k] = f_18 * skk_1116[k]
                    + f_3 * pc_z[k] * slk_1440[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece16(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pc,
                                                           const size_t skk, const size_t sli0,
                                                           const size_t sli1, const size_t slk,
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
    const auto f_20 = 3.0 / q;

    auto *t_1803 = buffer.data(target + 1803);
    auto *t_1804 = buffer.data(target + 1804);
    auto *t_1805 = buffer.data(target + 1805);
    auto *t_1806 = buffer.data(target + 1806);
    auto *t_1807 = buffer.data(target + 1807);
    auto *t_1808 = buffer.data(target + 1808);
    auto *t_1809 = buffer.data(target + 1809);
    auto *t_1810 = buffer.data(target + 1810);
    auto *t_1811 = buffer.data(target + 1811);
    auto *t_1812 = buffer.data(target + 1812);
    auto *t_1813 = buffer.data(target + 1813);
    auto *t_1814 = buffer.data(target + 1814);
    auto *t_1815 = buffer.data(target + 1815);
    auto *t_1816 = buffer.data(target + 1816);
    auto *t_1817 = buffer.data(target + 1817);
    auto *t_1818 = buffer.data(target + 1818);
    auto *t_1819 = buffer.data(target + 1819);
    auto *t_1820 = buffer.data(target + 1820);
    auto *t_1821 = buffer.data(target + 1821);
    auto *t_1822 = buffer.data(target + 1822);
    auto *t_1823 = buffer.data(target + 1823);
    auto *t_1824 = buffer.data(target + 1824);
    auto *t_1825 = buffer.data(target + 1825);
    auto *t_1826 = buffer.data(target + 1826);
    auto *t_1827 = buffer.data(target + 1827);
    auto *t_1828 = buffer.data(target + 1828);
    auto *t_1829 = buffer.data(target + 1829);
    auto *t_1830 = buffer.data(target + 1830);
    auto *t_1831 = buffer.data(target + 1831);
    auto *t_1832 = buffer.data(target + 1832);
    auto *t_1833 = buffer.data(target + 1833);
    auto *t_1834 = buffer.data(target + 1834);
    auto *t_1835 = buffer.data(target + 1835);
    auto *t_1836 = buffer.data(target + 1836);
    auto *t_1837 = buffer.data(target + 1837);
    auto *t_1838 = buffer.data(target + 1838);
    auto *t_1839 = buffer.data(target + 1839);
    auto *t_1840 = buffer.data(target + 1840);
    auto *t_1841 = buffer.data(target + 1841);
    auto *t_1842 = buffer.data(target + 1842);
    auto *t_1843 = buffer.data(target + 1843);
    auto *t_1844 = buffer.data(target + 1844);
    auto *t_1845 = buffer.data(target + 1845);
    auto *t_1846 = buffer.data(target + 1846);
    auto *t_1847 = buffer.data(target + 1847);
    auto *t_1848 = buffer.data(target + 1848);
    auto *t_1849 = buffer.data(target + 1849);
    auto *t_1850 = buffer.data(target + 1850);
    auto *t_1851 = buffer.data(target + 1851);
    auto *t_1852 = buffer.data(target + 1852);
    auto *t_1853 = buffer.data(target + 1853);
    auto *t_1854 = buffer.data(target + 1854);
    auto *t_1855 = buffer.data(target + 1855);
    auto *t_1856 = buffer.data(target + 1856);
    auto *t_1857 = buffer.data(target + 1857);
    auto *t_1858 = buffer.data(target + 1858);
    auto *t_1859 = buffer.data(target + 1859);
    auto *t_1860 = buffer.data(target + 1860);
    auto *t_1861 = buffer.data(target + 1861);
    auto *t_1862 = buffer.data(target + 1862);
    auto *t_1863 = buffer.data(target + 1863);
    auto *t_1864 = buffer.data(target + 1864);
    auto *t_1865 = buffer.data(target + 1865);
    auto *t_1866 = buffer.data(target + 1866);
    auto *t_1867 = buffer.data(target + 1867);
    auto *t_1868 = buffer.data(target + 1868);
    auto *t_1869 = buffer.data(target + 1869);
    auto *t_1870 = buffer.data(target + 1870);
    auto *t_1871 = buffer.data(target + 1871);
    auto *t_1872 = buffer.data(target + 1872);
    auto *t_1873 = buffer.data(target + 1873);
    auto *t_1874 = buffer.data(target + 1874);
    auto *t_1875 = buffer.data(target + 1875);
    auto *t_1876 = buffer.data(target + 1876);
    auto *t_1877 = buffer.data(target + 1877);
    auto *t_1878 = buffer.data(target + 1878);
    auto *t_1879 = buffer.data(target + 1879);
    auto *t_1880 = buffer.data(target + 1880);
    auto *t_1881 = buffer.data(target + 1881);
    auto *t_1882 = buffer.data(target + 1882);
    auto *t_1883 = buffer.data(target + 1883);
    auto *t_1884 = buffer.data(target + 1884);
    auto *t_1885 = buffer.data(target + 1885);
    auto *t_1886 = buffer.data(target + 1886);
    auto *t_1887 = buffer.data(target + 1887);
    auto *t_1888 = buffer.data(target + 1888);
    auto *t_1889 = buffer.data(target + 1889);
    auto *t_1890 = buffer.data(target + 1890);
    auto *t_1891 = buffer.data(target + 1891);
    auto *t_1892 = buffer.data(target + 1892);
    auto *t_1893 = buffer.data(target + 1893);
    auto *t_1894 = buffer.data(target + 1894);
    auto *t_1895 = buffer.data(target + 1895);
    auto *t_1896 = buffer.data(target + 1896);
    auto *t_1897 = buffer.data(target + 1897);
    auto *t_1898 = buffer.data(target + 1898);
    auto *t_1899 = buffer.data(target + 1899);
    auto *t_1900 = buffer.data(target + 1900);
    auto *t_1901 = buffer.data(target + 1901);
    auto *t_1902 = buffer.data(target + 1902);
    auto *t_1903 = buffer.data(target + 1903);
    auto *t_1904 = buffer.data(target + 1904);
    auto *t_1905 = buffer.data(target + 1905);
    auto *t_1906 = buffer.data(target + 1906);
    auto *t_1907 = buffer.data(target + 1907);
    auto *t_1908 = buffer.data(target + 1908);
    auto *t_1909 = buffer.data(target + 1909);
    auto *t_1910 = buffer.data(target + 1910);
    auto *t_1911 = buffer.data(target + 1911);
    auto *t_1912 = buffer.data(target + 1912);
    auto *t_1913 = buffer.data(target + 1913);
    auto *t_1914 = buffer.data(target + 1914);
    auto *t_1915 = buffer.data(target + 1915);
    auto *t_1916 = buffer.data(target + 1916);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skk_1119 = buffer.data(skk + 1119);
    const auto *skk_1122 = buffer.data(skk + 1122);
    const auto *skk_1126 = buffer.data(skk + 1126);
    const auto *skk_1131 = buffer.data(skk + 1131);
    const auto *skk_1144 = buffer.data(skk + 1144);
    const auto *skk_1151 = buffer.data(skk + 1151);
    const auto *skk_1152 = buffer.data(skk + 1152);
    const auto *skk_1154 = buffer.data(skk + 1154);
    const auto *skk_1155 = buffer.data(skk + 1155);
    const auto *skk_1157 = buffer.data(skk + 1157);
    const auto *skk_1158 = buffer.data(skk + 1158);
    const auto *skk_1161 = buffer.data(skk + 1161);
    const auto *skk_1162 = buffer.data(skk + 1162);
    const auto *skk_1166 = buffer.data(skk + 1166);
    const auto *skk_1167 = buffer.data(skk + 1167);
    const auto *skk_1172 = buffer.data(skk + 1172);
    const auto *skk_1180 = buffer.data(skk + 1180);
    const auto *skk_1182 = buffer.data(skk + 1182);
    const auto *skk_1183 = buffer.data(skk + 1183);
    const auto *skk_1184 = buffer.data(skk + 1184);
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
    const auto *skk_1202 = buffer.data(skk + 1202);
    const auto *skk_1203 = buffer.data(skk + 1203);
    const auto *skk_1208 = buffer.data(skk + 1208);
    const auto *skk_1216 = buffer.data(skk + 1216);
    const auto *skk_1218 = buffer.data(skk + 1218);
    const auto *skk_1219 = buffer.data(skk + 1219);
    const auto *skk_1220 = buffer.data(skk + 1220);
    const auto *skk_1221 = buffer.data(skk + 1221);
    const auto *skk_1222 = buffer.data(skk + 1222);
    const auto *skk_1223 = buffer.data(skk + 1223);
    const auto *skk_1224 = buffer.data(skk + 1224);
    const auto *skk_1226 = buffer.data(skk + 1226);
    const auto *skk_1229 = buffer.data(skk + 1229);
    const auto *skk_1233 = buffer.data(skk + 1233);
    const auto *skk_1238 = buffer.data(skk + 1238);
    const auto *skk_1244 = buffer.data(skk + 1244);

    const auto *sli0_1123 = buffer.data(sli0 + 1123);
    const auto *sli0_1125 = buffer.data(sli0 + 1125);
    const auto *sli0_1126 = buffer.data(sli0 + 1126);
    const auto *sli0_1129 = buffer.data(sli0 + 1129);
    const auto *sli0_1130 = buffer.data(sli0 + 1130);
    const auto *sli0_1132 = buffer.data(sli0 + 1132);
    const auto *sli0_1134 = buffer.data(sli0 + 1134);
    const auto *sli0_1135 = buffer.data(sli0 + 1135);
    const auto *sli0_1137 = buffer.data(sli0 + 1137);
    const auto *sli0_1138 = buffer.data(sli0 + 1138);
    const auto *sli0_1140 = buffer.data(sli0 + 1140);
    const auto *sli0_1141 = buffer.data(sli0 + 1141);
    const auto *sli0_1143 = buffer.data(sli0 + 1143);
    const auto *sli0_1144 = buffer.data(sli0 + 1144);
    const auto *sli0_1145 = buffer.data(sli0 + 1145);
    const auto *sli0_1146 = buffer.data(sli0 + 1146);
    const auto *sli0_1147 = buffer.data(sli0 + 1147);
    const auto *sli0_1148 = buffer.data(sli0 + 1148);
    const auto *sli0_1151 = buffer.data(sli0 + 1151);
    const auto *sli0_1153 = buffer.data(sli0 + 1153);
    const auto *sli0_1154 = buffer.data(sli0 + 1154);
    const auto *sli0_1157 = buffer.data(sli0 + 1157);
    const auto *sli0_1158 = buffer.data(sli0 + 1158);
    const auto *sli0_1160 = buffer.data(sli0 + 1160);
    const auto *sli0_1162 = buffer.data(sli0 + 1162);
    const auto *sli0_1163 = buffer.data(sli0 + 1163);
    const auto *sli0_1165 = buffer.data(sli0 + 1165);
    const auto *sli0_1166 = buffer.data(sli0 + 1166);
    const auto *sli0_1168 = buffer.data(sli0 + 1168);
    const auto *sli0_1169 = buffer.data(sli0 + 1169);
    const auto *sli0_1171 = buffer.data(sli0 + 1171);
    const auto *sli0_1172 = buffer.data(sli0 + 1172);
    const auto *sli0_1173 = buffer.data(sli0 + 1173);
    const auto *sli0_1174 = buffer.data(sli0 + 1174);
    const auto *sli0_1175 = buffer.data(sli0 + 1175);
    const auto *sli0_1176 = buffer.data(sli0 + 1176);
    const auto *sli0_1179 = buffer.data(sli0 + 1179);
    const auto *sli0_1181 = buffer.data(sli0 + 1181);
    const auto *sli0_1182 = buffer.data(sli0 + 1182);
    const auto *sli0_1185 = buffer.data(sli0 + 1185);
    const auto *sli0_1186 = buffer.data(sli0 + 1186);
    const auto *sli0_1188 = buffer.data(sli0 + 1188);
    const auto *sli0_1190 = buffer.data(sli0 + 1190);
    const auto *sli0_1191 = buffer.data(sli0 + 1191);
    const auto *sli0_1193 = buffer.data(sli0 + 1193);
    const auto *sli0_1194 = buffer.data(sli0 + 1194);
    const auto *sli0_1196 = buffer.data(sli0 + 1196);
    const auto *sli0_1197 = buffer.data(sli0 + 1197);
    const auto *sli0_1199 = buffer.data(sli0 + 1199);
    const auto *sli0_1200 = buffer.data(sli0 + 1200);
    const auto *sli0_1201 = buffer.data(sli0 + 1201);

    const auto *sli1_1123 = buffer.data(sli1 + 1123);
    const auto *sli1_1125 = buffer.data(sli1 + 1125);
    const auto *sli1_1126 = buffer.data(sli1 + 1126);
    const auto *sli1_1129 = buffer.data(sli1 + 1129);
    const auto *sli1_1130 = buffer.data(sli1 + 1130);
    const auto *sli1_1132 = buffer.data(sli1 + 1132);
    const auto *sli1_1134 = buffer.data(sli1 + 1134);
    const auto *sli1_1135 = buffer.data(sli1 + 1135);
    const auto *sli1_1137 = buffer.data(sli1 + 1137);
    const auto *sli1_1138 = buffer.data(sli1 + 1138);
    const auto *sli1_1140 = buffer.data(sli1 + 1140);
    const auto *sli1_1141 = buffer.data(sli1 + 1141);
    const auto *sli1_1143 = buffer.data(sli1 + 1143);
    const auto *sli1_1144 = buffer.data(sli1 + 1144);
    const auto *sli1_1145 = buffer.data(sli1 + 1145);
    const auto *sli1_1146 = buffer.data(sli1 + 1146);
    const auto *sli1_1147 = buffer.data(sli1 + 1147);
    const auto *sli1_1148 = buffer.data(sli1 + 1148);
    const auto *sli1_1151 = buffer.data(sli1 + 1151);
    const auto *sli1_1153 = buffer.data(sli1 + 1153);
    const auto *sli1_1154 = buffer.data(sli1 + 1154);
    const auto *sli1_1157 = buffer.data(sli1 + 1157);
    const auto *sli1_1158 = buffer.data(sli1 + 1158);
    const auto *sli1_1160 = buffer.data(sli1 + 1160);
    const auto *sli1_1162 = buffer.data(sli1 + 1162);
    const auto *sli1_1163 = buffer.data(sli1 + 1163);
    const auto *sli1_1165 = buffer.data(sli1 + 1165);
    const auto *sli1_1166 = buffer.data(sli1 + 1166);
    const auto *sli1_1168 = buffer.data(sli1 + 1168);
    const auto *sli1_1169 = buffer.data(sli1 + 1169);
    const auto *sli1_1171 = buffer.data(sli1 + 1171);
    const auto *sli1_1172 = buffer.data(sli1 + 1172);
    const auto *sli1_1173 = buffer.data(sli1 + 1173);
    const auto *sli1_1174 = buffer.data(sli1 + 1174);
    const auto *sli1_1175 = buffer.data(sli1 + 1175);
    const auto *sli1_1176 = buffer.data(sli1 + 1176);
    const auto *sli1_1179 = buffer.data(sli1 + 1179);
    const auto *sli1_1181 = buffer.data(sli1 + 1181);
    const auto *sli1_1182 = buffer.data(sli1 + 1182);
    const auto *sli1_1185 = buffer.data(sli1 + 1185);
    const auto *sli1_1186 = buffer.data(sli1 + 1186);
    const auto *sli1_1188 = buffer.data(sli1 + 1188);
    const auto *sli1_1190 = buffer.data(sli1 + 1190);
    const auto *sli1_1191 = buffer.data(sli1 + 1191);
    const auto *sli1_1193 = buffer.data(sli1 + 1193);
    const auto *sli1_1194 = buffer.data(sli1 + 1194);
    const auto *sli1_1196 = buffer.data(sli1 + 1196);
    const auto *sli1_1197 = buffer.data(sli1 + 1197);
    const auto *sli1_1199 = buffer.data(sli1 + 1199);
    const auto *sli1_1200 = buffer.data(sli1 + 1200);
    const auto *sli1_1201 = buffer.data(sli1 + 1201);

    const auto *slk_1442 = buffer.data(slk + 1442);
    const auto *slk_1443 = buffer.data(slk + 1443);
    const auto *slk_1445 = buffer.data(slk + 1445);
    const auto *slk_1446 = buffer.data(slk + 1446);
    const auto *slk_1449 = buffer.data(slk + 1449);
    const auto *slk_1450 = buffer.data(slk + 1450);
    const auto *slk_1452 = buffer.data(slk + 1452);
    const auto *slk_1454 = buffer.data(slk + 1454);
    const auto *slk_1455 = buffer.data(slk + 1455);
    const auto *slk_1457 = buffer.data(slk + 1457);
    const auto *slk_1458 = buffer.data(slk + 1458);
    const auto *slk_1460 = buffer.data(slk + 1460);
    const auto *slk_1461 = buffer.data(slk + 1461);
    const auto *slk_1463 = buffer.data(slk + 1463);
    const auto *slk_1464 = buffer.data(slk + 1464);
    const auto *slk_1465 = buffer.data(slk + 1465);
    const auto *slk_1467 = buffer.data(slk + 1467);
    const auto *slk_1468 = buffer.data(slk + 1468);
    const auto *slk_1469 = buffer.data(slk + 1469);
    const auto *slk_1470 = buffer.data(slk + 1470);
    const auto *slk_1471 = buffer.data(slk + 1471);
    const auto *slk_1472 = buffer.data(slk + 1472);
    const auto *slk_1473 = buffer.data(slk + 1473);
    const auto *slk_1474 = buffer.data(slk + 1474);
    const auto *slk_1475 = buffer.data(slk + 1475);
    const auto *slk_1476 = buffer.data(slk + 1476);
    const auto *slk_1478 = buffer.data(slk + 1478);
    const auto *slk_1479 = buffer.data(slk + 1479);
    const auto *slk_1481 = buffer.data(slk + 1481);
    const auto *slk_1482 = buffer.data(slk + 1482);
    const auto *slk_1485 = buffer.data(slk + 1485);
    const auto *slk_1486 = buffer.data(slk + 1486);
    const auto *slk_1488 = buffer.data(slk + 1488);
    const auto *slk_1490 = buffer.data(slk + 1490);
    const auto *slk_1491 = buffer.data(slk + 1491);
    const auto *slk_1493 = buffer.data(slk + 1493);
    const auto *slk_1494 = buffer.data(slk + 1494);
    const auto *slk_1496 = buffer.data(slk + 1496);
    const auto *slk_1497 = buffer.data(slk + 1497);
    const auto *slk_1499 = buffer.data(slk + 1499);
    const auto *slk_1500 = buffer.data(slk + 1500);
    const auto *slk_1501 = buffer.data(slk + 1501);
    const auto *slk_1503 = buffer.data(slk + 1503);
    const auto *slk_1504 = buffer.data(slk + 1504);
    const auto *slk_1505 = buffer.data(slk + 1505);
    const auto *slk_1506 = buffer.data(slk + 1506);
    const auto *slk_1507 = buffer.data(slk + 1507);
    const auto *slk_1508 = buffer.data(slk + 1508);
    const auto *slk_1509 = buffer.data(slk + 1509);
    const auto *slk_1510 = buffer.data(slk + 1510);
    const auto *slk_1511 = buffer.data(slk + 1511);
    const auto *slk_1512 = buffer.data(slk + 1512);
    const auto *slk_1514 = buffer.data(slk + 1514);
    const auto *slk_1515 = buffer.data(slk + 1515);
    const auto *slk_1517 = buffer.data(slk + 1517);
    const auto *slk_1518 = buffer.data(slk + 1518);
    const auto *slk_1521 = buffer.data(slk + 1521);
    const auto *slk_1522 = buffer.data(slk + 1522);
    const auto *slk_1524 = buffer.data(slk + 1524);
    const auto *slk_1526 = buffer.data(slk + 1526);
    const auto *slk_1527 = buffer.data(slk + 1527);
    const auto *slk_1529 = buffer.data(slk + 1529);
    const auto *slk_1530 = buffer.data(slk + 1530);
    const auto *slk_1532 = buffer.data(slk + 1532);
    const auto *slk_1533 = buffer.data(slk + 1533);
    const auto *slk_1535 = buffer.data(slk + 1535);
    const auto *slk_1536 = buffer.data(slk + 1536);
    const auto *slk_1537 = buffer.data(slk + 1537);

#pragma omp simd aligned(t_1803, t_1804, t_1805, pc_x, pc_y, skk_1154, sli0_1123, sli0_1125, \
                         sli1_1123, sli1_1125, slk_1442, slk_1443, \
                         slk_1445 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1803[k] = f_4 * sli0_1123[k]
                    - f_5 * sli1_1123[k]
                    + f_3 * pc_x[k] * slk_1443[k];

        t_1804[k] = f_18 * skk_1154[k]
                    + f_3 * pc_y[k] * slk_1442[k];

        t_1805[k] = f_4 * sli0_1125[k]
                    - f_5 * sli1_1125[k]
                    + f_3 * pc_x[k] * slk_1445[k];
    }

#pragma omp simd aligned(t_1806, t_1807, t_1808, pc_x, pc_y, pc_z, skk_1119, skk_1157, \
                         sli0_1126, sli1_1126, slk_1443, slk_1445, \
                         slk_1446 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1806[k] = f_6 * sli0_1126[k]
                    - f_7 * sli1_1126[k]
                    + f_3 * pc_x[k] * slk_1446[k];

        t_1807[k] = f_18 * skk_1119[k]
                    + f_3 * pc_z[k] * slk_1443[k];

        t_1808[k] = f_18 * skk_1157[k]
                    + f_3 * pc_y[k] * slk_1445[k];
    }

#pragma omp simd aligned(t_1809, t_1810, t_1811, pc_x, pc_z, skk_1122, sli0_1129, sli0_1130, \
                         sli1_1129, sli1_1130, slk_1446, slk_1449, \
                         slk_1450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1809[k] = f_6 * sli0_1129[k]
                    - f_7 * sli1_1129[k]
                    + f_3 * pc_x[k] * slk_1449[k];

        t_1810[k] = f_8 * sli0_1130[k]
                    - f_9 * sli1_1130[k]
                    + f_3 * pc_x[k] * slk_1450[k];

        t_1811[k] = f_18 * skk_1122[k]
                    + f_3 * pc_z[k] * slk_1446[k];
    }

#pragma omp simd aligned(t_1812, t_1813, t_1814, pc_x, pc_y, skk_1161, sli0_1132, sli0_1134, \
                         sli1_1132, sli1_1134, slk_1449, slk_1452, \
                         slk_1454 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1812[k] = f_8 * sli0_1132[k]
                    - f_9 * sli1_1132[k]
                    + f_3 * pc_x[k] * slk_1452[k];

        t_1813[k] = f_18 * skk_1161[k]
                    + f_3 * pc_y[k] * slk_1449[k];

        t_1814[k] = f_8 * sli0_1134[k]
                    - f_9 * sli1_1134[k]
                    + f_3 * pc_x[k] * slk_1454[k];
    }

#pragma omp simd aligned(t_1815, t_1816, t_1817, pc_x, pc_z, skk_1126, sli0_1135, sli0_1137, \
                         sli1_1135, sli1_1137, slk_1450, slk_1455, \
                         slk_1457 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1815[k] = f_10 * sli0_1135[k]
                    - f_11 * sli1_1135[k]
                    + f_3 * pc_x[k] * slk_1455[k];

        t_1816[k] = f_18 * skk_1126[k]
                    + f_3 * pc_z[k] * slk_1450[k];

        t_1817[k] = f_10 * sli0_1137[k]
                    - f_11 * sli1_1137[k]
                    + f_3 * pc_x[k] * slk_1457[k];
    }

#pragma omp simd aligned(t_1818, t_1819, t_1820, pc_x, pc_y, skk_1166, sli0_1138, sli0_1140, \
                         sli1_1138, sli1_1140, slk_1454, slk_1458, \
                         slk_1460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1818[k] = f_10 * sli0_1138[k]
                    - f_11 * sli1_1138[k]
                    + f_3 * pc_x[k] * slk_1458[k];

        t_1819[k] = f_18 * skk_1166[k]
                    + f_3 * pc_y[k] * slk_1454[k];

        t_1820[k] = f_10 * sli0_1140[k]
                    - f_11 * sli1_1140[k]
                    + f_3 * pc_x[k] * slk_1460[k];
    }

#pragma omp simd aligned(t_1821, t_1822, t_1823, pc_x, pc_z, skk_1131, sli0_1141, sli0_1143, \
                         sli1_1141, sli1_1143, slk_1455, slk_1461, \
                         slk_1463 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1821[k] = f_12 * sli0_1141[k]
                    - f_13 * sli1_1141[k]
                    + f_3 * pc_x[k] * slk_1461[k];

        t_1822[k] = f_18 * skk_1131[k]
                    + f_3 * pc_z[k] * slk_1455[k];

        t_1823[k] = f_12 * sli0_1143[k]
                    - f_13 * sli1_1143[k]
                    + f_3 * pc_x[k] * slk_1463[k];
    }

#pragma omp simd aligned(t_1824, t_1825, t_1826, pc_x, pc_y, skk_1172, sli0_1144, sli0_1145, \
                         sli1_1144, sli1_1145, slk_1460, slk_1464, \
                         slk_1465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1824[k] = f_12 * sli0_1144[k]
                    - f_13 * sli1_1144[k]
                    + f_3 * pc_x[k] * slk_1464[k];

        t_1825[k] = f_12 * sli0_1145[k]
                    - f_13 * sli1_1145[k]
                    + f_3 * pc_x[k] * slk_1465[k];

        t_1826[k] = f_18 * skk_1172[k]
                    + f_3 * pc_y[k] * slk_1460[k];
    }

#pragma omp simd aligned(t_1827, t_1828, t_1829, t_1830, t_1831, t_1832, pc_x, sli0_1147, \
                         sli1_1147, slk_1467, slk_1468, slk_1469, slk_1470, slk_1471, \
                         slk_1472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1827[k] = f_12 * sli0_1147[k]
                    - f_13 * sli1_1147[k]
                    + f_3 * pc_x[k] * slk_1467[k];

        t_1828[k] = f_3 * pc_x[k] * slk_1468[k];

        t_1829[k] = f_3 * pc_x[k] * slk_1469[k];

        t_1830[k] = f_3 * pc_x[k] * slk_1470[k];

        t_1831[k] = f_3 * pc_x[k] * slk_1471[k];

        t_1832[k] = f_3 * pc_x[k] * slk_1472[k];
    }

#pragma omp simd aligned(t_1833, t_1834, t_1835, t_1836, t_1837, pc_x, pc_y, pc_z, skk_1144, \
                         skk_1180, sli0_1141, sli1_1141, slk_1468, slk_1473, slk_1474, \
                         slk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1833[k] = f_3 * pc_x[k] * slk_1473[k];

        t_1834[k] = f_3 * pc_x[k] * slk_1474[k];

        t_1835[k] = f_3 * pc_x[k] * slk_1475[k];

        t_1836[k] = f_18 * skk_1180[k]
                    + f_1 * sli0_1141[k]
                    - f_2 * sli1_1141[k]
                    + f_3 * pc_y[k] * slk_1468[k];

        t_1837[k] = f_18 * skk_1144[k]
                    + f_3 * pc_z[k] * slk_1468[k];
    }

#pragma omp simd aligned(t_1838, t_1839, t_1840, pc_y, skk_1182, skk_1183, skk_1184, \
                         sli0_1143, sli0_1144, sli0_1145, sli1_1143, sli1_1144, sli1_1145, \
                         slk_1470, slk_1471, slk_1472 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1838[k] = f_18 * skk_1182[k]
                    + f_4 * sli0_1143[k]
                    - f_5 * sli1_1143[k]
                    + f_3 * pc_y[k] * slk_1470[k];

        t_1839[k] = f_18 * skk_1183[k]
                    + f_6 * sli0_1144[k]
                    - f_7 * sli1_1144[k]
                    + f_3 * pc_y[k] * slk_1471[k];

        t_1840[k] = f_18 * skk_1184[k]
                    + f_8 * sli0_1145[k]
                    - f_9 * sli1_1145[k]
                    + f_3 * pc_y[k] * slk_1472[k];
    }

#pragma omp simd aligned(t_1841, t_1842, t_1843, pc_y, skk_1185, skk_1186, skk_1187, \
                         sli0_1146, sli0_1147, sli1_1146, sli1_1147, slk_1473, slk_1474, \
                         slk_1475 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1841[k] = f_18 * skk_1185[k]
                    + f_10 * sli0_1146[k]
                    - f_11 * sli1_1146[k]
                    + f_3 * pc_y[k] * slk_1473[k];

        t_1842[k] = f_18 * skk_1186[k]
                    + f_12 * sli0_1147[k]
                    - f_13 * sli1_1147[k]
                    + f_3 * pc_y[k] * slk_1474[k];

        t_1843[k] = f_18 * skk_1187[k]
                    + f_3 * pc_y[k] * slk_1475[k];
    }

#pragma omp simd aligned(t_1844, t_1845, t_1846, t_1847, pc_x, pc_y, pc_z, skk_1151, skk_1152, \
                         skk_1188, sli0_1147, sli0_1148, sli1_1147, sli1_1148, slk_1475, \
                         slk_1476 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1844[k] = f_18 * skk_1151[k]
                    + f_1 * sli0_1147[k]
                    - f_2 * sli1_1147[k]
                    + f_3 * pc_z[k] * slk_1475[k];

        t_1845[k] = f_1 * sli0_1148[k]
                    - f_2 * sli1_1148[k]
                    + f_3 * pc_x[k] * slk_1476[k];

        t_1846[k] = f_17 * skk_1188[k]
                    + f_3 * pc_y[k] * slk_1476[k];

        t_1847[k] = f_19 * skk_1152[k]
                    + f_3 * pc_z[k] * slk_1476[k];
    }

#pragma omp simd aligned(t_1848, t_1849, t_1850, pc_x, pc_y, skk_1190, sli0_1151, sli0_1153, \
                         sli1_1151, sli1_1153, slk_1478, slk_1479, \
                         slk_1481 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1848[k] = f_4 * sli0_1151[k]
                    - f_5 * sli1_1151[k]
                    + f_3 * pc_x[k] * slk_1479[k];

        t_1849[k] = f_17 * skk_1190[k]
                    + f_3 * pc_y[k] * slk_1478[k];

        t_1850[k] = f_4 * sli0_1153[k]
                    - f_5 * sli1_1153[k]
                    + f_3 * pc_x[k] * slk_1481[k];
    }

#pragma omp simd aligned(t_1851, t_1852, t_1853, pc_x, pc_y, pc_z, skk_1155, skk_1193, \
                         sli0_1154, sli1_1154, slk_1479, slk_1481, \
                         slk_1482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1851[k] = f_6 * sli0_1154[k]
                    - f_7 * sli1_1154[k]
                    + f_3 * pc_x[k] * slk_1482[k];

        t_1852[k] = f_19 * skk_1155[k]
                    + f_3 * pc_z[k] * slk_1479[k];

        t_1853[k] = f_17 * skk_1193[k]
                    + f_3 * pc_y[k] * slk_1481[k];
    }

#pragma omp simd aligned(t_1854, t_1855, t_1856, pc_x, pc_z, skk_1158, sli0_1157, sli0_1158, \
                         sli1_1157, sli1_1158, slk_1482, slk_1485, \
                         slk_1486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1854[k] = f_6 * sli0_1157[k]
                    - f_7 * sli1_1157[k]
                    + f_3 * pc_x[k] * slk_1485[k];

        t_1855[k] = f_8 * sli0_1158[k]
                    - f_9 * sli1_1158[k]
                    + f_3 * pc_x[k] * slk_1486[k];

        t_1856[k] = f_19 * skk_1158[k]
                    + f_3 * pc_z[k] * slk_1482[k];
    }

#pragma omp simd aligned(t_1857, t_1858, t_1859, pc_x, pc_y, skk_1197, sli0_1160, sli0_1162, \
                         sli1_1160, sli1_1162, slk_1485, slk_1488, \
                         slk_1490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1857[k] = f_8 * sli0_1160[k]
                    - f_9 * sli1_1160[k]
                    + f_3 * pc_x[k] * slk_1488[k];

        t_1858[k] = f_17 * skk_1197[k]
                    + f_3 * pc_y[k] * slk_1485[k];

        t_1859[k] = f_8 * sli0_1162[k]
                    - f_9 * sli1_1162[k]
                    + f_3 * pc_x[k] * slk_1490[k];
    }

#pragma omp simd aligned(t_1860, t_1861, t_1862, pc_x, pc_z, skk_1162, sli0_1163, sli0_1165, \
                         sli1_1163, sli1_1165, slk_1486, slk_1491, \
                         slk_1493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1860[k] = f_10 * sli0_1163[k]
                    - f_11 * sli1_1163[k]
                    + f_3 * pc_x[k] * slk_1491[k];

        t_1861[k] = f_19 * skk_1162[k]
                    + f_3 * pc_z[k] * slk_1486[k];

        t_1862[k] = f_10 * sli0_1165[k]
                    - f_11 * sli1_1165[k]
                    + f_3 * pc_x[k] * slk_1493[k];
    }

#pragma omp simd aligned(t_1863, t_1864, t_1865, pc_x, pc_y, skk_1202, sli0_1166, sli0_1168, \
                         sli1_1166, sli1_1168, slk_1490, slk_1494, \
                         slk_1496 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1863[k] = f_10 * sli0_1166[k]
                    - f_11 * sli1_1166[k]
                    + f_3 * pc_x[k] * slk_1494[k];

        t_1864[k] = f_17 * skk_1202[k]
                    + f_3 * pc_y[k] * slk_1490[k];

        t_1865[k] = f_10 * sli0_1168[k]
                    - f_11 * sli1_1168[k]
                    + f_3 * pc_x[k] * slk_1496[k];
    }

#pragma omp simd aligned(t_1866, t_1867, t_1868, pc_x, pc_z, skk_1167, sli0_1169, sli0_1171, \
                         sli1_1169, sli1_1171, slk_1491, slk_1497, \
                         slk_1499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1866[k] = f_12 * sli0_1169[k]
                    - f_13 * sli1_1169[k]
                    + f_3 * pc_x[k] * slk_1497[k];

        t_1867[k] = f_19 * skk_1167[k]
                    + f_3 * pc_z[k] * slk_1491[k];

        t_1868[k] = f_12 * sli0_1171[k]
                    - f_13 * sli1_1171[k]
                    + f_3 * pc_x[k] * slk_1499[k];
    }

#pragma omp simd aligned(t_1869, t_1870, t_1871, pc_x, pc_y, skk_1208, sli0_1172, sli0_1173, \
                         sli1_1172, sli1_1173, slk_1496, slk_1500, \
                         slk_1501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1869[k] = f_12 * sli0_1172[k]
                    - f_13 * sli1_1172[k]
                    + f_3 * pc_x[k] * slk_1500[k];

        t_1870[k] = f_12 * sli0_1173[k]
                    - f_13 * sli1_1173[k]
                    + f_3 * pc_x[k] * slk_1501[k];

        t_1871[k] = f_17 * skk_1208[k]
                    + f_3 * pc_y[k] * slk_1496[k];
    }

#pragma omp simd aligned(t_1872, t_1873, t_1874, t_1875, t_1876, t_1877, pc_x, sli0_1175, \
                         sli1_1175, slk_1503, slk_1504, slk_1505, slk_1506, slk_1507, \
                         slk_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1872[k] = f_12 * sli0_1175[k]
                    - f_13 * sli1_1175[k]
                    + f_3 * pc_x[k] * slk_1503[k];

        t_1873[k] = f_3 * pc_x[k] * slk_1504[k];

        t_1874[k] = f_3 * pc_x[k] * slk_1505[k];

        t_1875[k] = f_3 * pc_x[k] * slk_1506[k];

        t_1876[k] = f_3 * pc_x[k] * slk_1507[k];

        t_1877[k] = f_3 * pc_x[k] * slk_1508[k];
    }

#pragma omp simd aligned(t_1878, t_1879, t_1880, t_1881, t_1882, pc_x, pc_y, pc_z, skk_1180, \
                         skk_1216, sli0_1169, sli1_1169, slk_1504, slk_1509, slk_1510, \
                         slk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1878[k] = f_3 * pc_x[k] * slk_1509[k];

        t_1879[k] = f_3 * pc_x[k] * slk_1510[k];

        t_1880[k] = f_3 * pc_x[k] * slk_1511[k];

        t_1881[k] = f_17 * skk_1216[k]
                    + f_1 * sli0_1169[k]
                    - f_2 * sli1_1169[k]
                    + f_3 * pc_y[k] * slk_1504[k];

        t_1882[k] = f_19 * skk_1180[k]
                    + f_3 * pc_z[k] * slk_1504[k];
    }

#pragma omp simd aligned(t_1883, t_1884, t_1885, pc_y, skk_1218, skk_1219, skk_1220, \
                         sli0_1171, sli0_1172, sli0_1173, sli1_1171, sli1_1172, sli1_1173, \
                         slk_1506, slk_1507, slk_1508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1883[k] = f_17 * skk_1218[k]
                    + f_4 * sli0_1171[k]
                    - f_5 * sli1_1171[k]
                    + f_3 * pc_y[k] * slk_1506[k];

        t_1884[k] = f_17 * skk_1219[k]
                    + f_6 * sli0_1172[k]
                    - f_7 * sli1_1172[k]
                    + f_3 * pc_y[k] * slk_1507[k];

        t_1885[k] = f_17 * skk_1220[k]
                    + f_8 * sli0_1173[k]
                    - f_9 * sli1_1173[k]
                    + f_3 * pc_y[k] * slk_1508[k];
    }

#pragma omp simd aligned(t_1886, t_1887, t_1888, pc_y, skk_1221, skk_1222, skk_1223, \
                         sli0_1174, sli0_1175, sli1_1174, sli1_1175, slk_1509, slk_1510, \
                         slk_1511 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1886[k] = f_17 * skk_1221[k]
                    + f_10 * sli0_1174[k]
                    - f_11 * sli1_1174[k]
                    + f_3 * pc_y[k] * slk_1509[k];

        t_1887[k] = f_17 * skk_1222[k]
                    + f_12 * sli0_1175[k]
                    - f_13 * sli1_1175[k]
                    + f_3 * pc_y[k] * slk_1510[k];

        t_1888[k] = f_17 * skk_1223[k]
                    + f_3 * pc_y[k] * slk_1511[k];
    }

#pragma omp simd aligned(t_1889, t_1890, t_1891, t_1892, pc_x, pc_y, pc_z, skk_1187, skk_1188, \
                         skk_1224, sli0_1175, sli0_1176, sli1_1175, sli1_1176, slk_1511, \
                         slk_1512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1889[k] = f_19 * skk_1187[k]
                    + f_1 * sli0_1175[k]
                    - f_2 * sli1_1175[k]
                    + f_3 * pc_z[k] * slk_1511[k];

        t_1890[k] = f_1 * sli0_1176[k]
                    - f_2 * sli1_1176[k]
                    + f_3 * pc_x[k] * slk_1512[k];

        t_1891[k] = f_16 * skk_1224[k]
                    + f_3 * pc_y[k] * slk_1512[k];

        t_1892[k] = f_20 * skk_1188[k]
                    + f_3 * pc_z[k] * slk_1512[k];
    }

#pragma omp simd aligned(t_1893, t_1894, t_1895, pc_x, pc_y, skk_1226, sli0_1179, sli0_1181, \
                         sli1_1179, sli1_1181, slk_1514, slk_1515, \
                         slk_1517 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1893[k] = f_4 * sli0_1179[k]
                    - f_5 * sli1_1179[k]
                    + f_3 * pc_x[k] * slk_1515[k];

        t_1894[k] = f_16 * skk_1226[k]
                    + f_3 * pc_y[k] * slk_1514[k];

        t_1895[k] = f_4 * sli0_1181[k]
                    - f_5 * sli1_1181[k]
                    + f_3 * pc_x[k] * slk_1517[k];
    }

#pragma omp simd aligned(t_1896, t_1897, t_1898, pc_x, pc_y, pc_z, skk_1191, skk_1229, \
                         sli0_1182, sli1_1182, slk_1515, slk_1517, \
                         slk_1518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1896[k] = f_6 * sli0_1182[k]
                    - f_7 * sli1_1182[k]
                    + f_3 * pc_x[k] * slk_1518[k];

        t_1897[k] = f_20 * skk_1191[k]
                    + f_3 * pc_z[k] * slk_1515[k];

        t_1898[k] = f_16 * skk_1229[k]
                    + f_3 * pc_y[k] * slk_1517[k];
    }

#pragma omp simd aligned(t_1899, t_1900, t_1901, pc_x, pc_z, skk_1194, sli0_1185, sli0_1186, \
                         sli1_1185, sli1_1186, slk_1518, slk_1521, \
                         slk_1522 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1899[k] = f_6 * sli0_1185[k]
                    - f_7 * sli1_1185[k]
                    + f_3 * pc_x[k] * slk_1521[k];

        t_1900[k] = f_8 * sli0_1186[k]
                    - f_9 * sli1_1186[k]
                    + f_3 * pc_x[k] * slk_1522[k];

        t_1901[k] = f_20 * skk_1194[k]
                    + f_3 * pc_z[k] * slk_1518[k];
    }

#pragma omp simd aligned(t_1902, t_1903, t_1904, pc_x, pc_y, skk_1233, sli0_1188, sli0_1190, \
                         sli1_1188, sli1_1190, slk_1521, slk_1524, \
                         slk_1526 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1902[k] = f_8 * sli0_1188[k]
                    - f_9 * sli1_1188[k]
                    + f_3 * pc_x[k] * slk_1524[k];

        t_1903[k] = f_16 * skk_1233[k]
                    + f_3 * pc_y[k] * slk_1521[k];

        t_1904[k] = f_8 * sli0_1190[k]
                    - f_9 * sli1_1190[k]
                    + f_3 * pc_x[k] * slk_1526[k];
    }

#pragma omp simd aligned(t_1905, t_1906, t_1907, pc_x, pc_z, skk_1198, sli0_1191, sli0_1193, \
                         sli1_1191, sli1_1193, slk_1522, slk_1527, \
                         slk_1529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1905[k] = f_10 * sli0_1191[k]
                    - f_11 * sli1_1191[k]
                    + f_3 * pc_x[k] * slk_1527[k];

        t_1906[k] = f_20 * skk_1198[k]
                    + f_3 * pc_z[k] * slk_1522[k];

        t_1907[k] = f_10 * sli0_1193[k]
                    - f_11 * sli1_1193[k]
                    + f_3 * pc_x[k] * slk_1529[k];
    }

#pragma omp simd aligned(t_1908, t_1909, t_1910, pc_x, pc_y, skk_1238, sli0_1194, sli0_1196, \
                         sli1_1194, sli1_1196, slk_1526, slk_1530, \
                         slk_1532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1908[k] = f_10 * sli0_1194[k]
                    - f_11 * sli1_1194[k]
                    + f_3 * pc_x[k] * slk_1530[k];

        t_1909[k] = f_16 * skk_1238[k]
                    + f_3 * pc_y[k] * slk_1526[k];

        t_1910[k] = f_10 * sli0_1196[k]
                    - f_11 * sli1_1196[k]
                    + f_3 * pc_x[k] * slk_1532[k];
    }

#pragma omp simd aligned(t_1911, t_1912, t_1913, pc_x, pc_z, skk_1203, sli0_1197, sli0_1199, \
                         sli1_1197, sli1_1199, slk_1527, slk_1533, \
                         slk_1535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1911[k] = f_12 * sli0_1197[k]
                    - f_13 * sli1_1197[k]
                    + f_3 * pc_x[k] * slk_1533[k];

        t_1912[k] = f_20 * skk_1203[k]
                    + f_3 * pc_z[k] * slk_1527[k];

        t_1913[k] = f_12 * sli0_1199[k]
                    - f_13 * sli1_1199[k]
                    + f_3 * pc_x[k] * slk_1535[k];
    }

#pragma omp simd aligned(t_1914, t_1915, t_1916, pc_x, pc_y, skk_1244, sli0_1200, sli0_1201, \
                         sli1_1200, sli1_1201, slk_1532, slk_1536, \
                         slk_1537 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1914[k] = f_12 * sli0_1200[k]
                    - f_13 * sli1_1200[k]
                    + f_3 * pc_x[k] * slk_1536[k];

        t_1915[k] = f_12 * sli0_1201[k]
                    - f_13 * sli1_1201[k]
                    + f_3 * pc_x[k] * slk_1537[k];

        t_1916[k] = f_16 * skk_1244[k]
                    + f_3 * pc_y[k] * slk_1532[k];
    }
}

static auto
compute_prim_sll_three_center_electron_repulsion_0_piece17(CSimdMatrix &buffer,
                                                           const size_t target, const size_t pb,
                                                           const size_t pc, const size_t skl0,
                                                           const size_t skk, const size_t skl1,
                                                           const size_t sli0, const size_t sli1,
                                                           const size_t slk, const size_t ncols,
                                                           const double gamma, const double p,
                                                           const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.0 / q;
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
    const auto f_21 = 3.5 / q;

    auto *t_1917 = buffer.data(target + 1917);
    auto *t_1918 = buffer.data(target + 1918);
    auto *t_1919 = buffer.data(target + 1919);
    auto *t_1920 = buffer.data(target + 1920);
    auto *t_1921 = buffer.data(target + 1921);
    auto *t_1922 = buffer.data(target + 1922);
    auto *t_1923 = buffer.data(target + 1923);
    auto *t_1924 = buffer.data(target + 1924);
    auto *t_1925 = buffer.data(target + 1925);
    auto *t_1926 = buffer.data(target + 1926);
    auto *t_1927 = buffer.data(target + 1927);
    auto *t_1928 = buffer.data(target + 1928);
    auto *t_1929 = buffer.data(target + 1929);
    auto *t_1930 = buffer.data(target + 1930);
    auto *t_1931 = buffer.data(target + 1931);
    auto *t_1932 = buffer.data(target + 1932);
    auto *t_1933 = buffer.data(target + 1933);
    auto *t_1934 = buffer.data(target + 1934);
    auto *t_1935 = buffer.data(target + 1935);
    auto *t_1936 = buffer.data(target + 1936);
    auto *t_1937 = buffer.data(target + 1937);
    auto *t_1938 = buffer.data(target + 1938);
    auto *t_1939 = buffer.data(target + 1939);
    auto *t_1940 = buffer.data(target + 1940);
    auto *t_1941 = buffer.data(target + 1941);
    auto *t_1942 = buffer.data(target + 1942);
    auto *t_1943 = buffer.data(target + 1943);
    auto *t_1944 = buffer.data(target + 1944);
    auto *t_1945 = buffer.data(target + 1945);
    auto *t_1946 = buffer.data(target + 1946);
    auto *t_1947 = buffer.data(target + 1947);
    auto *t_1948 = buffer.data(target + 1948);
    auto *t_1949 = buffer.data(target + 1949);
    auto *t_1950 = buffer.data(target + 1950);
    auto *t_1951 = buffer.data(target + 1951);
    auto *t_1952 = buffer.data(target + 1952);
    auto *t_1953 = buffer.data(target + 1953);
    auto *t_1954 = buffer.data(target + 1954);
    auto *t_1955 = buffer.data(target + 1955);
    auto *t_1956 = buffer.data(target + 1956);
    auto *t_1957 = buffer.data(target + 1957);
    auto *t_1958 = buffer.data(target + 1958);
    auto *t_1959 = buffer.data(target + 1959);
    auto *t_1960 = buffer.data(target + 1960);
    auto *t_1961 = buffer.data(target + 1961);
    auto *t_1962 = buffer.data(target + 1962);
    auto *t_1963 = buffer.data(target + 1963);
    auto *t_1964 = buffer.data(target + 1964);
    auto *t_1965 = buffer.data(target + 1965);
    auto *t_1966 = buffer.data(target + 1966);
    auto *t_1967 = buffer.data(target + 1967);
    auto *t_1968 = buffer.data(target + 1968);
    auto *t_1969 = buffer.data(target + 1969);
    auto *t_1970 = buffer.data(target + 1970);
    auto *t_1971 = buffer.data(target + 1971);
    auto *t_1972 = buffer.data(target + 1972);
    auto *t_1973 = buffer.data(target + 1973);
    auto *t_1974 = buffer.data(target + 1974);
    auto *t_1975 = buffer.data(target + 1975);
    auto *t_1976 = buffer.data(target + 1976);
    auto *t_1977 = buffer.data(target + 1977);
    auto *t_1978 = buffer.data(target + 1978);
    auto *t_1979 = buffer.data(target + 1979);
    auto *t_1980 = buffer.data(target + 1980);
    auto *t_1981 = buffer.data(target + 1981);
    auto *t_1982 = buffer.data(target + 1982);
    auto *t_1983 = buffer.data(target + 1983);
    auto *t_1984 = buffer.data(target + 1984);
    auto *t_1985 = buffer.data(target + 1985);
    auto *t_1986 = buffer.data(target + 1986);
    auto *t_1987 = buffer.data(target + 1987);
    auto *t_1988 = buffer.data(target + 1988);
    auto *t_1989 = buffer.data(target + 1989);
    auto *t_1990 = buffer.data(target + 1990);
    auto *t_1991 = buffer.data(target + 1991);
    auto *t_1992 = buffer.data(target + 1992);
    auto *t_1993 = buffer.data(target + 1993);
    auto *t_1994 = buffer.data(target + 1994);
    auto *t_1995 = buffer.data(target + 1995);
    auto *t_1996 = buffer.data(target + 1996);
    auto *t_1997 = buffer.data(target + 1997);
    auto *t_1998 = buffer.data(target + 1998);
    auto *t_1999 = buffer.data(target + 1999);
    auto *t_2000 = buffer.data(target + 2000);
    auto *t_2001 = buffer.data(target + 2001);
    auto *t_2002 = buffer.data(target + 2002);
    auto *t_2003 = buffer.data(target + 2003);
    auto *t_2004 = buffer.data(target + 2004);
    auto *t_2005 = buffer.data(target + 2005);
    auto *t_2006 = buffer.data(target + 2006);
    auto *t_2007 = buffer.data(target + 2007);
    auto *t_2008 = buffer.data(target + 2008);
    auto *t_2009 = buffer.data(target + 2009);
    auto *t_2010 = buffer.data(target + 2010);
    auto *t_2011 = buffer.data(target + 2011);
    auto *t_2012 = buffer.data(target + 2012);
    auto *t_2013 = buffer.data(target + 2013);
    auto *t_2014 = buffer.data(target + 2014);
    auto *t_2015 = buffer.data(target + 2015);
    auto *t_2016 = buffer.data(target + 2016);
    auto *t_2017 = buffer.data(target + 2017);
    auto *t_2018 = buffer.data(target + 2018);
    auto *t_2019 = buffer.data(target + 2019);
    auto *t_2020 = buffer.data(target + 2020);
    auto *t_2021 = buffer.data(target + 2021);
    auto *t_2022 = buffer.data(target + 2022);
    auto *t_2023 = buffer.data(target + 2023);
    auto *t_2024 = buffer.data(target + 2024);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *skl0_1575 = buffer.data(skl0 + 1575);
    const auto *skl0_1580 = buffer.data(skl0 + 1580);
    const auto *skl0_1584 = buffer.data(skl0 + 1584);
    const auto *skl0_1589 = buffer.data(skl0 + 1589);
    const auto *skl0_1595 = buffer.data(skl0 + 1595);
    const auto *skl0_1602 = buffer.data(skl0 + 1602);
    const auto *skl0_1611 = buffer.data(skl0 + 1611);
    const auto *skl0_1613 = buffer.data(skl0 + 1613);
    const auto *skl0_1614 = buffer.data(skl0 + 1614);
    const auto *skl0_1615 = buffer.data(skl0 + 1615);
    const auto *skl0_1616 = buffer.data(skl0 + 1616);
    const auto *skl0_1617 = buffer.data(skl0 + 1617);
    const auto *skl0_1619 = buffer.data(skl0 + 1619);

    const auto *skk_1216 = buffer.data(skk + 1216);
    const auto *skk_1223 = buffer.data(skk + 1223);
    const auto *skk_1224 = buffer.data(skk + 1224);
    const auto *skk_1227 = buffer.data(skk + 1227);
    const auto *skk_1230 = buffer.data(skk + 1230);
    const auto *skk_1234 = buffer.data(skk + 1234);
    const auto *skk_1239 = buffer.data(skk + 1239);
    const auto *skk_1252 = buffer.data(skk + 1252);
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
    const auto *skk_1274 = buffer.data(skk + 1274);
    const auto *skk_1275 = buffer.data(skk + 1275);
    const auto *skk_1280 = buffer.data(skk + 1280);
    const auto *skk_1288 = buffer.data(skk + 1288);
    const auto *skk_1290 = buffer.data(skk + 1290);
    const auto *skk_1291 = buffer.data(skk + 1291);
    const auto *skk_1292 = buffer.data(skk + 1292);
    const auto *skk_1293 = buffer.data(skk + 1293);
    const auto *skk_1294 = buffer.data(skk + 1294);
    const auto *skk_1295 = buffer.data(skk + 1295);

    const auto *skl1_1575 = buffer.data(skl1 + 1575);
    const auto *skl1_1580 = buffer.data(skl1 + 1580);
    const auto *skl1_1584 = buffer.data(skl1 + 1584);
    const auto *skl1_1589 = buffer.data(skl1 + 1589);
    const auto *skl1_1595 = buffer.data(skl1 + 1595);
    const auto *skl1_1602 = buffer.data(skl1 + 1602);
    const auto *skl1_1611 = buffer.data(skl1 + 1611);
    const auto *skl1_1613 = buffer.data(skl1 + 1613);
    const auto *skl1_1614 = buffer.data(skl1 + 1614);
    const auto *skl1_1615 = buffer.data(skl1 + 1615);
    const auto *skl1_1616 = buffer.data(skl1 + 1616);
    const auto *skl1_1617 = buffer.data(skl1 + 1617);
    const auto *skl1_1619 = buffer.data(skl1 + 1619);

    const auto *sli0_1197 = buffer.data(sli0 + 1197);
    const auto *sli0_1199 = buffer.data(sli0 + 1199);
    const auto *sli0_1200 = buffer.data(sli0 + 1200);
    const auto *sli0_1201 = buffer.data(sli0 + 1201);
    const auto *sli0_1202 = buffer.data(sli0 + 1202);
    const auto *sli0_1203 = buffer.data(sli0 + 1203);
    const auto *sli0_1207 = buffer.data(sli0 + 1207);
    const auto *sli0_1210 = buffer.data(sli0 + 1210);
    const auto *sli0_1214 = buffer.data(sli0 + 1214);
    const auto *sli0_1216 = buffer.data(sli0 + 1216);
    const auto *sli0_1219 = buffer.data(sli0 + 1219);
    const auto *sli0_1221 = buffer.data(sli0 + 1221);
    const auto *sli0_1222 = buffer.data(sli0 + 1222);
    const auto *sli0_1225 = buffer.data(sli0 + 1225);
    const auto *sli0_1227 = buffer.data(sli0 + 1227);
    const auto *sli0_1228 = buffer.data(sli0 + 1228);
    const auto *sli0_1229 = buffer.data(sli0 + 1229);
    const auto *sli0_1232 = buffer.data(sli0 + 1232);
    const auto *sli0_1235 = buffer.data(sli0 + 1235);
    const auto *sli0_1237 = buffer.data(sli0 + 1237);
    const auto *sli0_1238 = buffer.data(sli0 + 1238);
    const auto *sli0_1241 = buffer.data(sli0 + 1241);
    const auto *sli0_1242 = buffer.data(sli0 + 1242);
    const auto *sli0_1244 = buffer.data(sli0 + 1244);
    const auto *sli0_1246 = buffer.data(sli0 + 1246);
    const auto *sli0_1247 = buffer.data(sli0 + 1247);
    const auto *sli0_1249 = buffer.data(sli0 + 1249);
    const auto *sli0_1250 = buffer.data(sli0 + 1250);
    const auto *sli0_1252 = buffer.data(sli0 + 1252);
    const auto *sli0_1253 = buffer.data(sli0 + 1253);
    const auto *sli0_1255 = buffer.data(sli0 + 1255);
    const auto *sli0_1256 = buffer.data(sli0 + 1256);
    const auto *sli0_1257 = buffer.data(sli0 + 1257);
    const auto *sli0_1258 = buffer.data(sli0 + 1258);
    const auto *sli0_1259 = buffer.data(sli0 + 1259);

    const auto *sli1_1197 = buffer.data(sli1 + 1197);
    const auto *sli1_1199 = buffer.data(sli1 + 1199);
    const auto *sli1_1200 = buffer.data(sli1 + 1200);
    const auto *sli1_1201 = buffer.data(sli1 + 1201);
    const auto *sli1_1202 = buffer.data(sli1 + 1202);
    const auto *sli1_1203 = buffer.data(sli1 + 1203);
    const auto *sli1_1207 = buffer.data(sli1 + 1207);
    const auto *sli1_1210 = buffer.data(sli1 + 1210);
    const auto *sli1_1214 = buffer.data(sli1 + 1214);
    const auto *sli1_1216 = buffer.data(sli1 + 1216);
    const auto *sli1_1219 = buffer.data(sli1 + 1219);
    const auto *sli1_1221 = buffer.data(sli1 + 1221);
    const auto *sli1_1222 = buffer.data(sli1 + 1222);
    const auto *sli1_1225 = buffer.data(sli1 + 1225);
    const auto *sli1_1227 = buffer.data(sli1 + 1227);
    const auto *sli1_1228 = buffer.data(sli1 + 1228);
    const auto *sli1_1229 = buffer.data(sli1 + 1229);
    const auto *sli1_1232 = buffer.data(sli1 + 1232);
    const auto *sli1_1235 = buffer.data(sli1 + 1235);
    const auto *sli1_1237 = buffer.data(sli1 + 1237);
    const auto *sli1_1238 = buffer.data(sli1 + 1238);
    const auto *sli1_1241 = buffer.data(sli1 + 1241);
    const auto *sli1_1242 = buffer.data(sli1 + 1242);
    const auto *sli1_1244 = buffer.data(sli1 + 1244);
    const auto *sli1_1246 = buffer.data(sli1 + 1246);
    const auto *sli1_1247 = buffer.data(sli1 + 1247);
    const auto *sli1_1249 = buffer.data(sli1 + 1249);
    const auto *sli1_1250 = buffer.data(sli1 + 1250);
    const auto *sli1_1252 = buffer.data(sli1 + 1252);
    const auto *sli1_1253 = buffer.data(sli1 + 1253);
    const auto *sli1_1255 = buffer.data(sli1 + 1255);
    const auto *sli1_1256 = buffer.data(sli1 + 1256);
    const auto *sli1_1257 = buffer.data(sli1 + 1257);
    const auto *sli1_1258 = buffer.data(sli1 + 1258);
    const auto *sli1_1259 = buffer.data(sli1 + 1259);

    const auto *slk_1539 = buffer.data(slk + 1539);
    const auto *slk_1540 = buffer.data(slk + 1540);
    const auto *slk_1541 = buffer.data(slk + 1541);
    const auto *slk_1542 = buffer.data(slk + 1542);
    const auto *slk_1543 = buffer.data(slk + 1543);
    const auto *slk_1544 = buffer.data(slk + 1544);
    const auto *slk_1545 = buffer.data(slk + 1545);
    const auto *slk_1546 = buffer.data(slk + 1546);
    const auto *slk_1547 = buffer.data(slk + 1547);
    const auto *slk_1548 = buffer.data(slk + 1548);
    const auto *slk_1550 = buffer.data(slk + 1550);
    const auto *slk_1551 = buffer.data(slk + 1551);
    const auto *slk_1553 = buffer.data(slk + 1553);
    const auto *slk_1554 = buffer.data(slk + 1554);
    const auto *slk_1557 = buffer.data(slk + 1557);
    const auto *slk_1558 = buffer.data(slk + 1558);
    const auto *slk_1560 = buffer.data(slk + 1560);
    const auto *slk_1562 = buffer.data(slk + 1562);
    const auto *slk_1563 = buffer.data(slk + 1563);
    const auto *slk_1565 = buffer.data(slk + 1565);
    const auto *slk_1566 = buffer.data(slk + 1566);
    const auto *slk_1568 = buffer.data(slk + 1568);
    const auto *slk_1569 = buffer.data(slk + 1569);
    const auto *slk_1571 = buffer.data(slk + 1571);
    const auto *slk_1572 = buffer.data(slk + 1572);
    const auto *slk_1573 = buffer.data(slk + 1573);
    const auto *slk_1576 = buffer.data(slk + 1576);
    const auto *slk_1577 = buffer.data(slk + 1577);
    const auto *slk_1578 = buffer.data(slk + 1578);
    const auto *slk_1579 = buffer.data(slk + 1579);
    const auto *slk_1580 = buffer.data(slk + 1580);
    const auto *slk_1581 = buffer.data(slk + 1581);
    const auto *slk_1582 = buffer.data(slk + 1582);
    const auto *slk_1583 = buffer.data(slk + 1583);
    const auto *slk_1584 = buffer.data(slk + 1584);
    const auto *slk_1586 = buffer.data(slk + 1586);
    const auto *slk_1587 = buffer.data(slk + 1587);
    const auto *slk_1589 = buffer.data(slk + 1589);
    const auto *slk_1590 = buffer.data(slk + 1590);
    const auto *slk_1593 = buffer.data(slk + 1593);
    const auto *slk_1594 = buffer.data(slk + 1594);
    const auto *slk_1596 = buffer.data(slk + 1596);
    const auto *slk_1598 = buffer.data(slk + 1598);
    const auto *slk_1599 = buffer.data(slk + 1599);
    const auto *slk_1601 = buffer.data(slk + 1601);
    const auto *slk_1602 = buffer.data(slk + 1602);
    const auto *slk_1604 = buffer.data(slk + 1604);
    const auto *slk_1605 = buffer.data(slk + 1605);
    const auto *slk_1607 = buffer.data(slk + 1607);
    const auto *slk_1608 = buffer.data(slk + 1608);
    const auto *slk_1609 = buffer.data(slk + 1609);
    const auto *slk_1611 = buffer.data(slk + 1611);
    const auto *slk_1612 = buffer.data(slk + 1612);
    const auto *slk_1613 = buffer.data(slk + 1613);
    const auto *slk_1614 = buffer.data(slk + 1614);
    const auto *slk_1615 = buffer.data(slk + 1615);
    const auto *slk_1616 = buffer.data(slk + 1616);
    const auto *slk_1617 = buffer.data(slk + 1617);
    const auto *slk_1618 = buffer.data(slk + 1618);
    const auto *slk_1619 = buffer.data(slk + 1619);

#pragma omp simd aligned(t_1917, t_1918, t_1919, t_1920, t_1921, t_1922, pc_x, sli0_1203, \
                         sli1_1203, slk_1539, slk_1540, slk_1541, slk_1542, slk_1543, \
                         slk_1544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1917[k] = f_12 * sli0_1203[k]
                    - f_13 * sli1_1203[k]
                    + f_3 * pc_x[k] * slk_1539[k];

        t_1918[k] = f_3 * pc_x[k] * slk_1540[k];

        t_1919[k] = f_3 * pc_x[k] * slk_1541[k];

        t_1920[k] = f_3 * pc_x[k] * slk_1542[k];

        t_1921[k] = f_3 * pc_x[k] * slk_1543[k];

        t_1922[k] = f_3 * pc_x[k] * slk_1544[k];
    }

#pragma omp simd aligned(t_1923, t_1924, t_1925, t_1926, t_1927, pc_x, pc_y, pc_z, skk_1216, \
                         skk_1252, sli0_1197, sli1_1197, slk_1540, slk_1545, slk_1546, \
                         slk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1923[k] = f_3 * pc_x[k] * slk_1545[k];

        t_1924[k] = f_3 * pc_x[k] * slk_1546[k];

        t_1925[k] = f_3 * pc_x[k] * slk_1547[k];

        t_1926[k] = f_16 * skk_1252[k]
                    + f_1 * sli0_1197[k]
                    - f_2 * sli1_1197[k]
                    + f_3 * pc_y[k] * slk_1540[k];

        t_1927[k] = f_20 * skk_1216[k]
                    + f_3 * pc_z[k] * slk_1540[k];
    }

#pragma omp simd aligned(t_1928, t_1929, t_1930, pc_y, skk_1254, skk_1255, skk_1256, \
                         sli0_1199, sli0_1200, sli0_1201, sli1_1199, sli1_1200, sli1_1201, \
                         slk_1542, slk_1543, slk_1544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1928[k] = f_16 * skk_1254[k]
                    + f_4 * sli0_1199[k]
                    - f_5 * sli1_1199[k]
                    + f_3 * pc_y[k] * slk_1542[k];

        t_1929[k] = f_16 * skk_1255[k]
                    + f_6 * sli0_1200[k]
                    - f_7 * sli1_1200[k]
                    + f_3 * pc_y[k] * slk_1543[k];

        t_1930[k] = f_16 * skk_1256[k]
                    + f_8 * sli0_1201[k]
                    - f_9 * sli1_1201[k]
                    + f_3 * pc_y[k] * slk_1544[k];
    }

#pragma omp simd aligned(t_1931, t_1932, t_1933, pc_y, skk_1257, skk_1258, skk_1259, \
                         sli0_1202, sli0_1203, sli1_1202, sli1_1203, slk_1545, slk_1546, \
                         slk_1547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1931[k] = f_16 * skk_1257[k]
                    + f_10 * sli0_1202[k]
                    - f_11 * sli1_1202[k]
                    + f_3 * pc_y[k] * slk_1545[k];

        t_1932[k] = f_16 * skk_1258[k]
                    + f_12 * sli0_1203[k]
                    - f_13 * sli1_1203[k]
                    + f_3 * pc_y[k] * slk_1546[k];

        t_1933[k] = f_16 * skk_1259[k]
                    + f_3 * pc_y[k] * slk_1547[k];
    }

#pragma omp simd aligned(t_1934, t_1935, t_1936, t_1937, pb_y, pc_y, pc_z, skl0_1575, \
                         skk_1223, skk_1224, skk_1260, skl1_1575, sli0_1203, sli1_1203, \
                         slk_1547, slk_1548 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1934[k] = f_20 * skk_1223[k]
                    + f_1 * sli0_1203[k]
                    - f_2 * sli1_1203[k]
                    + f_3 * pc_z[k] * slk_1547[k];

        t_1935[k] = pb_y[k] * skl0_1575[k]
                    - f_14 * pc_y[k] * skl1_1575[k];

        t_1936[k] = f_15 * skk_1260[k]
                    + f_3 * pc_y[k] * slk_1548[k];

        t_1937[k] = f_21 * skk_1224[k]
                    + f_3 * pc_z[k] * slk_1548[k];
    }

#pragma omp simd aligned(t_1938, t_1939, t_1940, pb_y, pc_x, pc_y, skl0_1580, skk_1262, \
                         skl1_1580, sli0_1207, sli1_1207, slk_1550, \
                         slk_1551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1938[k] = f_4 * sli0_1207[k]
                    - f_5 * sli1_1207[k]
                    + f_3 * pc_x[k] * slk_1551[k];

        t_1939[k] = f_15 * skk_1262[k]
                    + f_3 * pc_y[k] * slk_1550[k];

        t_1940[k] = pb_y[k] * skl0_1580[k]
                    - f_14 * pc_y[k] * skl1_1580[k];
    }

#pragma omp simd aligned(t_1941, t_1942, t_1943, pc_x, pc_y, pc_z, skk_1227, skk_1265, \
                         sli0_1210, sli1_1210, slk_1551, slk_1553, \
                         slk_1554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1941[k] = f_6 * sli0_1210[k]
                    - f_7 * sli1_1210[k]
                    + f_3 * pc_x[k] * slk_1554[k];

        t_1942[k] = f_21 * skk_1227[k]
                    + f_3 * pc_z[k] * slk_1551[k];

        t_1943[k] = f_15 * skk_1265[k]
                    + f_3 * pc_y[k] * slk_1553[k];
    }

#pragma omp simd aligned(t_1944, t_1945, t_1946, pb_y, pc_x, pc_y, pc_z, skl0_1584, skk_1230, \
                         skl1_1584, sli0_1214, sli1_1214, slk_1554, \
                         slk_1558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1944[k] = pb_y[k] * skl0_1584[k]
                    - f_14 * pc_y[k] * skl1_1584[k];

        t_1945[k] = f_8 * sli0_1214[k]
                    - f_9 * sli1_1214[k]
                    + f_3 * pc_x[k] * slk_1558[k];

        t_1946[k] = f_21 * skk_1230[k]
                    + f_3 * pc_z[k] * slk_1554[k];
    }

#pragma omp simd aligned(t_1947, t_1948, t_1949, pb_y, pc_x, pc_y, skl0_1589, skk_1269, \
                         skl1_1589, sli0_1216, sli1_1216, slk_1557, \
                         slk_1560 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1947[k] = f_8 * sli0_1216[k]
                    - f_9 * sli1_1216[k]
                    + f_3 * pc_x[k] * slk_1560[k];

        t_1948[k] = f_15 * skk_1269[k]
                    + f_3 * pc_y[k] * slk_1557[k];

        t_1949[k] = pb_y[k] * skl0_1589[k]
                    - f_14 * pc_y[k] * skl1_1589[k];
    }

#pragma omp simd aligned(t_1950, t_1951, t_1952, pc_x, pc_z, skk_1234, sli0_1219, sli0_1221, \
                         sli1_1219, sli1_1221, slk_1558, slk_1563, \
                         slk_1565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1950[k] = f_10 * sli0_1219[k]
                    - f_11 * sli1_1219[k]
                    + f_3 * pc_x[k] * slk_1563[k];

        t_1951[k] = f_21 * skk_1234[k]
                    + f_3 * pc_z[k] * slk_1558[k];

        t_1952[k] = f_10 * sli0_1221[k]
                    - f_11 * sli1_1221[k]
                    + f_3 * pc_x[k] * slk_1565[k];
    }

#pragma omp simd aligned(t_1953, t_1954, t_1955, pb_y, pc_x, pc_y, skl0_1595, skk_1274, \
                         skl1_1595, sli0_1222, sli1_1222, slk_1562, \
                         slk_1566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1953[k] = f_10 * sli0_1222[k]
                    - f_11 * sli1_1222[k]
                    + f_3 * pc_x[k] * slk_1566[k];

        t_1954[k] = f_15 * skk_1274[k]
                    + f_3 * pc_y[k] * slk_1562[k];

        t_1955[k] = pb_y[k] * skl0_1595[k]
                    - f_14 * pc_y[k] * skl1_1595[k];
    }

#pragma omp simd aligned(t_1956, t_1957, t_1958, pc_x, pc_z, skk_1239, sli0_1225, sli0_1227, \
                         sli1_1225, sli1_1227, slk_1563, slk_1569, \
                         slk_1571 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1956[k] = f_12 * sli0_1225[k]
                    - f_13 * sli1_1225[k]
                    + f_3 * pc_x[k] * slk_1569[k];

        t_1957[k] = f_21 * skk_1239[k]
                    + f_3 * pc_z[k] * slk_1563[k];

        t_1958[k] = f_12 * sli0_1227[k]
                    - f_13 * sli1_1227[k]
                    + f_3 * pc_x[k] * slk_1571[k];
    }

#pragma omp simd aligned(t_1959, t_1960, t_1961, pc_x, pc_y, skk_1280, sli0_1228, sli0_1229, \
                         sli1_1228, sli1_1229, slk_1568, slk_1572, \
                         slk_1573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1959[k] = f_12 * sli0_1228[k]
                    - f_13 * sli1_1228[k]
                    + f_3 * pc_x[k] * slk_1572[k];

        t_1960[k] = f_12 * sli0_1229[k]
                    - f_13 * sli1_1229[k]
                    + f_3 * pc_x[k] * slk_1573[k];

        t_1961[k] = f_15 * skk_1280[k]
                    + f_3 * pc_y[k] * slk_1568[k];
    }

#pragma omp simd aligned(t_1962, t_1963, t_1964, t_1965, t_1966, t_1967, pb_y, pc_x, pc_y, \
                         skl0_1602, skl1_1602, slk_1576, slk_1577, slk_1578, slk_1579, \
                         slk_1580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1962[k] = pb_y[k] * skl0_1602[k]
                    - f_14 * pc_y[k] * skl1_1602[k];

        t_1963[k] = f_3 * pc_x[k] * slk_1576[k];

        t_1964[k] = f_3 * pc_x[k] * slk_1577[k];

        t_1965[k] = f_3 * pc_x[k] * slk_1578[k];

        t_1966[k] = f_3 * pc_x[k] * slk_1579[k];

        t_1967[k] = f_3 * pc_x[k] * slk_1580[k];
    }

#pragma omp simd aligned(t_1968, t_1969, t_1970, t_1971, pb_y, pc_x, pc_y, skl0_1611, \
                         skk_1288, skl1_1611, slk_1581, slk_1582, \
                         slk_1583 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1968[k] = f_3 * pc_x[k] * slk_1581[k];

        t_1969[k] = f_3 * pc_x[k] * slk_1582[k];

        t_1970[k] = f_3 * pc_x[k] * slk_1583[k];

        t_1971[k] = pb_y[k] * skl0_1611[k]
                    + f_0 * skk_1288[k]
                    - f_14 * pc_y[k] * skl1_1611[k];
    }

#pragma omp simd aligned(t_1972, t_1973, t_1974, pb_y, pc_y, pc_z, skl0_1613, skl0_1614, \
                         skk_1252, skk_1290, skk_1291, skl1_1613, skl1_1614, \
                         slk_1576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1972[k] = f_21 * skk_1252[k]
                    + f_3 * pc_z[k] * slk_1576[k];

        t_1973[k] = pb_y[k] * skl0_1613[k]
                    + f_20 * skk_1290[k]
                    - f_14 * pc_y[k] * skl1_1613[k];

        t_1974[k] = pb_y[k] * skl0_1614[k]
                    + f_19 * skk_1291[k]
                    - f_14 * pc_y[k] * skl1_1614[k];
    }

#pragma omp simd aligned(t_1975, t_1976, t_1977, pb_y, pc_y, skl0_1615, skl0_1616, skl0_1617, \
                         skk_1292, skk_1293, skk_1294, skl1_1615, skl1_1616, \
                         skl1_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1975[k] = pb_y[k] * skl0_1615[k]
                    + f_18 * skk_1292[k]
                    - f_14 * pc_y[k] * skl1_1615[k];

        t_1976[k] = pb_y[k] * skl0_1616[k]
                    + f_17 * skk_1293[k]
                    - f_14 * pc_y[k] * skl1_1616[k];

        t_1977[k] = pb_y[k] * skl0_1617[k]
                    + f_16 * skk_1294[k]
                    - f_14 * pc_y[k] * skl1_1617[k];
    }

#pragma omp simd aligned(t_1978, t_1979, t_1980, t_1981, pb_y, pc_x, pc_y, skl0_1619, \
                         skk_1295, skl1_1619, sli0_1232, sli1_1232, slk_1583, \
                         slk_1584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1978[k] = f_15 * skk_1295[k]
                    + f_3 * pc_y[k] * slk_1583[k];

        t_1979[k] = pb_y[k] * skl0_1619[k]
                    - f_14 * pc_y[k] * skl1_1619[k];

        t_1980[k] = f_1 * sli0_1232[k]
                    - f_2 * sli1_1232[k]
                    + f_3 * pc_x[k] * slk_1584[k];

        t_1981[k] = f_3 * pc_y[k] * slk_1584[k];
    }

#pragma omp simd aligned(t_1982, t_1983, t_1984, t_1985, pc_x, pc_y, pc_z, skk_1260, \
                         sli0_1235, sli0_1237, sli1_1235, sli1_1237, slk_1584, slk_1586, \
                         slk_1587, slk_1589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1982[k] = f_0 * skk_1260[k]
                    + f_3 * pc_z[k] * slk_1584[k];

        t_1983[k] = f_4 * sli0_1235[k]
                    - f_5 * sli1_1235[k]
                    + f_3 * pc_x[k] * slk_1587[k];

        t_1984[k] = f_3 * pc_y[k] * slk_1586[k];

        t_1985[k] = f_4 * sli0_1237[k]
                    - f_5 * sli1_1237[k]
                    + f_3 * pc_x[k] * slk_1589[k];
    }

#pragma omp simd aligned(t_1986, t_1987, t_1988, t_1989, pc_x, pc_y, pc_z, skk_1263, \
                         sli0_1238, sli0_1241, sli1_1238, sli1_1241, slk_1587, slk_1589, \
                         slk_1590, slk_1593 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1986[k] = f_6 * sli0_1238[k]
                    - f_7 * sli1_1238[k]
                    + f_3 * pc_x[k] * slk_1590[k];

        t_1987[k] = f_0 * skk_1263[k]
                    + f_3 * pc_z[k] * slk_1587[k];

        t_1988[k] = f_3 * pc_y[k] * slk_1589[k];

        t_1989[k] = f_6 * sli0_1241[k]
                    - f_7 * sli1_1241[k]
                    + f_3 * pc_x[k] * slk_1593[k];
    }

#pragma omp simd aligned(t_1990, t_1991, t_1992, t_1993, pc_x, pc_y, pc_z, skk_1266, \
                         sli0_1242, sli0_1244, sli1_1242, sli1_1244, slk_1590, slk_1593, \
                         slk_1594, slk_1596 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1990[k] = f_8 * sli0_1242[k]
                    - f_9 * sli1_1242[k]
                    + f_3 * pc_x[k] * slk_1594[k];

        t_1991[k] = f_0 * skk_1266[k]
                    + f_3 * pc_z[k] * slk_1590[k];

        t_1992[k] = f_8 * sli0_1244[k]
                    - f_9 * sli1_1244[k]
                    + f_3 * pc_x[k] * slk_1596[k];

        t_1993[k] = f_3 * pc_y[k] * slk_1593[k];
    }

#pragma omp simd aligned(t_1994, t_1995, t_1996, pc_x, pc_z, skk_1270, sli0_1246, sli0_1247, \
                         sli1_1246, sli1_1247, slk_1594, slk_1598, \
                         slk_1599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1994[k] = f_8 * sli0_1246[k]
                    - f_9 * sli1_1246[k]
                    + f_3 * pc_x[k] * slk_1598[k];

        t_1995[k] = f_10 * sli0_1247[k]
                    - f_11 * sli1_1247[k]
                    + f_3 * pc_x[k] * slk_1599[k];

        t_1996[k] = f_0 * skk_1270[k]
                    + f_3 * pc_z[k] * slk_1594[k];
    }

#pragma omp simd aligned(t_1997, t_1998, t_1999, t_2000, pc_x, pc_y, sli0_1249, sli0_1250, \
                         sli0_1252, sli1_1249, sli1_1250, sli1_1252, slk_1598, slk_1601, \
                         slk_1602, slk_1604 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1997[k] = f_10 * sli0_1249[k]
                    - f_11 * sli1_1249[k]
                    + f_3 * pc_x[k] * slk_1601[k];

        t_1998[k] = f_10 * sli0_1250[k]
                    - f_11 * sli1_1250[k]
                    + f_3 * pc_x[k] * slk_1602[k];

        t_1999[k] = f_3 * pc_y[k] * slk_1598[k];

        t_2000[k] = f_10 * sli0_1252[k]
                    - f_11 * sli1_1252[k]
                    + f_3 * pc_x[k] * slk_1604[k];
    }

#pragma omp simd aligned(t_2001, t_2002, t_2003, pc_x, pc_z, skk_1275, sli0_1253, sli0_1255, \
                         sli1_1253, sli1_1255, slk_1599, slk_1605, \
                         slk_1607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2001[k] = f_12 * sli0_1253[k]
                    - f_13 * sli1_1253[k]
                    + f_3 * pc_x[k] * slk_1605[k];

        t_2002[k] = f_0 * skk_1275[k]
                    + f_3 * pc_z[k] * slk_1599[k];

        t_2003[k] = f_12 * sli0_1255[k]
                    - f_13 * sli1_1255[k]
                    + f_3 * pc_x[k] * slk_1607[k];
    }

#pragma omp simd aligned(t_2004, t_2005, t_2006, t_2007, pc_x, pc_y, sli0_1256, sli0_1257, \
                         sli0_1259, sli1_1256, sli1_1257, sli1_1259, slk_1604, slk_1608, \
                         slk_1609, slk_1611 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2004[k] = f_12 * sli0_1256[k]
                    - f_13 * sli1_1256[k]
                    + f_3 * pc_x[k] * slk_1608[k];

        t_2005[k] = f_12 * sli0_1257[k]
                    - f_13 * sli1_1257[k]
                    + f_3 * pc_x[k] * slk_1609[k];

        t_2006[k] = f_3 * pc_y[k] * slk_1604[k];

        t_2007[k] = f_12 * sli0_1259[k]
                    - f_13 * sli1_1259[k]
                    + f_3 * pc_x[k] * slk_1611[k];
    }

#pragma omp simd aligned(t_2008, t_2009, t_2010, t_2011, t_2012, t_2013, t_2014, pc_x, \
                         slk_1612, slk_1613, slk_1614, slk_1615, slk_1616, slk_1617, \
                         slk_1618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2008[k] = f_3 * pc_x[k] * slk_1612[k];

        t_2009[k] = f_3 * pc_x[k] * slk_1613[k];

        t_2010[k] = f_3 * pc_x[k] * slk_1614[k];

        t_2011[k] = f_3 * pc_x[k] * slk_1615[k];

        t_2012[k] = f_3 * pc_x[k] * slk_1616[k];

        t_2013[k] = f_3 * pc_x[k] * slk_1617[k];

        t_2014[k] = f_3 * pc_x[k] * slk_1618[k];
    }

#pragma omp simd aligned(t_2015, t_2016, t_2017, t_2018, pc_x, pc_y, pc_z, skk_1288, \
                         sli0_1253, sli0_1255, sli1_1253, sli1_1255, slk_1612, slk_1614, \
                         slk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2015[k] = f_3 * pc_x[k] * slk_1619[k];

        t_2016[k] = f_1 * sli0_1253[k]
                    - f_2 * sli1_1253[k]
                    + f_3 * pc_y[k] * slk_1612[k];

        t_2017[k] = f_0 * skk_1288[k]
                    + f_3 * pc_z[k] * slk_1612[k];

        t_2018[k] = f_4 * sli0_1255[k]
                    - f_5 * sli1_1255[k]
                    + f_3 * pc_y[k] * slk_1614[k];
    }

#pragma omp simd aligned(t_2019, t_2020, t_2021, pc_y, sli0_1256, sli0_1257, sli0_1258, \
                         sli1_1256, sli1_1257, sli1_1258, slk_1615, slk_1616, \
                         slk_1617 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2019[k] = f_6 * sli0_1256[k]
                    - f_7 * sli1_1256[k]
                    + f_3 * pc_y[k] * slk_1615[k];

        t_2020[k] = f_8 * sli0_1257[k]
                    - f_9 * sli1_1257[k]
                    + f_3 * pc_y[k] * slk_1616[k];

        t_2021[k] = f_10 * sli0_1258[k]
                    - f_11 * sli1_1258[k]
                    + f_3 * pc_y[k] * slk_1617[k];
    }

#pragma omp simd aligned(t_2022, t_2023, t_2024, pc_y, pc_z, skk_1295, sli0_1259, sli1_1259, \
                         slk_1618, slk_1619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_2022[k] = f_12 * sli0_1259[k]
                    - f_13 * sli1_1259[k]
                    + f_3 * pc_y[k] * slk_1618[k];

        t_2023[k] = f_3 * pc_y[k] * slk_1619[k];

        t_2024[k] = f_0 * skk_1295[k]
                    + f_1 * sli0_1259[k]
                    - f_2 * sli1_1259[k]
                    + f_3 * pc_z[k] * slk_1619[k];
    }
}

auto
compute_prim_sll_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t skl0, const size_t skk,
                                                   const size_t skl1, const size_t sli0,
                                                   const size_t sli1, const size_t slk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sll_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece5(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece6(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece7(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece8(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece9(buffer, target, pb, pc, skl0, skk,
                                                              skl1, sli0, sli1, slk, ncols,
                                                              gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece10(buffer, target, pb, pc, skl0,
                                                               skk, skl1, sli0, sli1, slk,
                                                               ncols, gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece11(buffer, target, pb, pc, skl0,
                                                               skk, skl1, sli0, sli1, slk,
                                                               ncols, gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece12(buffer, target, pb, pc, skl0,
                                                               skk, skl1, slk, ncols, gamma, p,
                                                               q);

    compute_prim_sll_three_center_electron_repulsion_0_piece13(buffer, target, pb, pc, skl0,
                                                               skk, skl1, slk, ncols, gamma, p,
                                                               q);

    compute_prim_sll_three_center_electron_repulsion_0_piece14(buffer, target, pb, pc, skl0,
                                                               skk, skl1, sli0, sli1, slk,
                                                               ncols, gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece15(buffer, target, pb, pc, skl0,
                                                               skk, skl1, sli0, sli1, slk,
                                                               ncols, gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece16(buffer, target, pc, skk, sli0,
                                                               sli1, slk, ncols, gamma, p, q);

    compute_prim_sll_three_center_electron_repulsion_0_piece17(buffer, target, pb, pc, skl0,
                                                               skk, skl1, sli0, sli1, slk,
                                                               ncols, gamma, p, q);
}

}  // namespace simdt3ceri
